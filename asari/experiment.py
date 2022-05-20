import os
import sys
import json
import pickle

from mass2chem.search import *
# jms-metabolite-services
from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

# from .samples import SimpleSample
from .mass_functions import complete_mass_paired_mapping
from .constructors import CompositeMap
from .json_encoder import NpEncoder

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from . import db

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


class ext_Experiment:
    '''
    Similar to metDataModel.core.Experiment with preprocessing methods.
    This encapsulates a set of LC-MS files using the same experimental method 
    (chromatography and ionization) to be processed together.
    E.g., data from postive ESI and negative ESI should not in the same ext_Experiment instance.
    '''
    def __init__(self, sample_registry, parameters):
        '''
        This is the overall container for all data in an experiment/project.
        sample_registry: dictionary of sample and selected data after mass track extraction.
            The bulk mass track data are not kept in memory unless specified so.
        parameters: processing parameters passed from main.py.
        '''
        self.sample_registry = sample_registry
        self.valid_sample_ids = self.get_valid_sample_ids()
        self.number_of_samples = len(self.valid_sample_ids)
        self.number_scans = self.get_max_scan_number(sample_registry)
        self.all_samples = self.all_sample_instances = []

        self.parameters = parameters
        self.output_dir = parameters['outdir']
        self.mode = parameters['mode']
        self.mz_tolerance_ppm = self.parameters['mz_tolerance_ppm']
        self.check_isotope_ratio = self.parameters['check_isotope_ratio']
        self.database_mode = parameters['database_mode']
        self.reference_sample_id = self.get_reference_sample_id()
        
    def get_reference_sample_id(self):
        '''
        get_reference_sample_id either by user specification, or
        as sample of most number_anchor_mz_pairs, limited to first 100 samples to search.
        '''
        if self.parameters['reference']:
            # match file name; k is sm['sample_id']
            for k,v in self.sample_registry.items():
                if os.path.basename(self.parameters['reference']) == os.path.basename(v['input_file']):
                    return k
        elif self.sample_registry:
            L = [(v['number_anchor_mz_pairs'], v['sample_id']) for v in self.sample_registry.values()]
            L.sort(reverse=True)
            ref = self.sample_registry[L[0][1]]
            print("\n    The reference sample is:\n    ||* %s *||\n" %ref['name'])
            print("Max reference retention time is %4.2f at scan number %d.\n" %(
                max(ref['list_retention_time']), ref['max_scan_number']))
            return ref['sample_id']
        else:
            return None
        
    def get_valid_sample_ids(self):
        return [k for k,v in self.sample_registry.items() if v['status:eic'] == 'passed']

    def get_max_scan_number(self, sample_registry):
        if sample_registry:
            return max([sample_registry[k]['max_scan_number'] for k in self.valid_sample_ids]) + 1
        else:
            return None

    def process_all(self):
        '''
        This is default asari workflow.
        1. Build MassGrid, using either pairwise (small study) or clustering method. 
            Choose one reference from all samples for the largest number of landmark m/z tracks.
        2. RT alignment via a LOWESS function, using selective landmark peaks.
        3. Build composite elution profile (composite_mass_tracks),
            by cumulative sum of mass tracks from all samples after RT correction.
        4. Global peak detection is performed on each composite massTrack.
        5. Mapping global peaks (i.e. features) back to all samples and extract sample specific peak areas.
            This completes the FeatureTable.
        '''
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid()
        if not self.parameters['rt_align_on']:
            self.CMAP.mock_rentention_alignment()

        self.CMAP.build_composite_tracks()
        self.CMAP.global_peak_detection()

    def export_all(self):
        '''
        Export all files.
        Annotation of features to empirical compounds is done here.
        '''
        self.CMAP.MassGrid.to_csv(
            os.path.join(self.parameters['outdir'], 'export', self.parameters['mass_grid_mapping']) )
        self.export_CMAP_pickle()
        self.annotate()
        self.export_feature_tables()
        self.export_log()

    def annotate(self):
        '''
        Export Feature_annotation as tsv, Annotated_empricalCompounds in both JSON and pickle.
        Reference databases can be pre-loaded.
        With db_mass_calibrate to theoretical values, 
        '''
        self.load_annotation_db()
        self.db_mass_calibrate()

        EED = ExperimentalEcpdDatabase(mode=self.mode, mz_tolerance_ppm=self.mz_tolerance_ppm)
        EED.build_from_list_peaks(self.CMAP.FeatureList, mz_tolerance_ppm=self.mz_tolerance_ppm, check_isotope_ratio=self.check_isotope_ratio)
        EED.extend_empCpd_annotation(self.KCD)
        EED.annotate_singletons(self.KCD)       
        # EED.dict_empCpds misses some features 
        EED.dict_empCpds = self.append_orphans_to_epmCpds(EED.dict_empCpds)

        self.export_peak_annotation(EED.dict_empCpds, self.KCD, 'Feature_annotation')

        if self.sample_registry:                    # check because some subcommands may not have sample_registry
            self.select_unique_compound_features(EED.dict_empCpds)
        
        # export JSON
        outfile = os.path.join(self.parameters['outdir'], 'Annotated_empricalCompounds.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(EED.dict_empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)

        # export pickle
        outfile = os.path.join(self.parameters['outdir'], 'export', 'epd.pickle')
        with open(outfile, 'wb') as f:
            pickle.dump(EED.dict_empCpds, f, pickle.HIGHEST_PROTOCOL)


    def export_CMAP_pickle(self):
        '''
        Export main CMAP data and MassGrid to pickle, which can be used for visual data exploration.
        '''
        _export = {
            '_number_of_samples_': self.CMAP._number_of_samples_,
            'rt_length': self.CMAP.rt_length,
            'dict_scan_rtime': self.CMAP.dict_scan_rtime,
            'list_mass_tracks': self.CMAP.composite_mass_tracks,
            'MassGrid': dict(self.CMAP.MassGrid),
        }
        outfile = os.path.join(self.parameters['outdir'], 'export', 'cmap.pickle')
        with open(outfile, 'wb') as f:
            pickle.dump(_export, f, pickle.HIGHEST_PROTOCOL)


    def load_annotation_db(self, src='hmdb4'):
        '''
        Database of known compound using JMS. Will add more options later.
        '''
        self.KCD = knownCompoundDatabase()
        self.KCD.mass_indexed_compounds = pickle.load( pkg_resources.open_binary(db, 'mass_indexed_compounds.pickle') )
        self.KCD.emp_cpds_trees = pickle.load( pkg_resources.open_binary(db, 'emp_cpds_trees.pickle') )


    def db_mass_calibrate(self, required_calibrate_threshold=0.000002):
        '''
        Use KCD.evaluate_mass_accuracy_ratio to check systematic mass shift.
        If greater than required_calibrate_threshold (default 2 ppm), 
        calibrate m/z values for the whole experiment by updating self.CMAP.FeatureList.

        good_reference_landmark_peaks: [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...]
        ref_id_num -> index number of mass track in MassGrid.
        '''
        mz_landmarks = [self.CMAP.MassGrid['mz'][p['ref_id_num']] for p in self.CMAP.good_reference_landmark_peaks]
        mass_accuracy_ratio = self.KCD.evaluate_mass_accuracy_ratio(mz_landmarks, mode=self.mode, mz_tolerance_ppm=10)
        if mass_accuracy_ratio:
            if abs(mass_accuracy_ratio) > required_calibrate_threshold:
                print("Mass shift is greater than %2.1f ppm. Correction applied." %(required_calibrate_threshold*1000000))
                _correction = mass_accuracy_ratio + 1
                for F in self.CMAP.FeatureList:
                    F['mz'] = F['mz'] / _correction
                    F['mz_corrected_by_division'] = _correction
        else:
            print("Mass accuracy check is skipped, too few mz_landmarks (%d) matched." %len(mz_landmarks))


    def append_orphans_to_epmCpds(self, dict_empCpds):
        '''
        # EED.dict_empCpds does not include features without formula match, and we add them there.
        '''
        all_feature_ids = []
        for _, V in dict_empCpds.items():
            all_feature_ids += [peak['id_number'] for peak in V['MS1_pseudo_Spectra']]

        orphans = [peak for peak in self.CMAP.FeatureList if peak['id_number'] not in all_feature_ids]
        new_id_start = len(dict_empCpds)
        for peak in orphans:
            dict_empCpds[new_id_start] = {'interim_id': new_id_start,
                    'neutral_formula_mass': '', 'neutral_formula': '',
                    'MS1_pseudo_Spectra': [peak] }
            new_id_start += 1

        return dict_empCpds


    def export_peak_annotation(self, dict_empCpds, KCD, export_file_name_prefix):
        '''
        Export feature annotation.
        interim_id is empCpd id. dict_empCpds as seen in JMS.
        
        '''
        s = "[peak]id_number\tmz\trtime\tapex(scan number)\t[EmpCpd]interim_id\t[EmpCpd]ion_relation\tneutral_formula\tneutral_formula_mass\
        \tname_1st_guess\tmatched_DB_shorts\tmatched_DB_records\n"
        
        for _, V in dict_empCpds.items():
            name_1st_guess, matched_DB_shorts, matched_DB_records = '', '', ''
            if 'list_matches' in V:
                list_matches = V['list_matches']
                if list_matches:
                    name_1st_guess = KCD.mass_indexed_compounds[list_matches[0][0]]['compounds'][0]['name']
                    matched_DB_shorts = ", ".join([ "(" + KCD.short_report_emp_cpd(xx[0]) + ")"  for xx in list_matches])
                    matched_DB_records = ", ".join([str(xx) for xx  in list_matches])

            for peak in V['MS1_pseudo_Spectra']:
                s += '\t'.join([str(x) for x in [
                    peak['id_number'], peak['mz'], peak['rtime'], peak['apex'], V['interim_id'], peak.get('ion_relation', ''),
                    V['neutral_formula'], V['neutral_formula_mass'],
                    name_1st_guess, matched_DB_shorts, matched_DB_records]]) + "\n"

        outfile = os.path.join(self.parameters['outdir'], export_file_name_prefix + '.tsv')
        with open(outfile, encoding='utf-8', mode='w') as O:
            O.write(s)

        print("\nAnnotation of %d Empirical compounds was written to %s.\n\n" %(len(dict_empCpds), outfile))

    def select_unique_compound_features(self, dict_empCpds):
        '''
        Get unique feature by highest composite peak area per empirical compound. 
        '''
        self.selected_unique_features = {}
        for interim_id, V in dict_empCpds.items():
            if len(V['MS1_pseudo_Spectra']) == 1:
                self.selected_unique_features[V['MS1_pseudo_Spectra'][0]['id_number']] = (
                        # empCpd id, neutral_formula, ion_relation
                        interim_id, V['neutral_formula'], 'singleton'
                )
            else:
                try:
                    all_peaks = sorted([(peak['peak_area'], peak['goodness_fitting'], peak) for peak in V['MS1_pseudo_Spectra']],
                                        reverse=True)
                    best_peak = all_peaks[0][2]
                except TypeError:
                    best_peak = V['MS1_pseudo_Spectra'][0]
                self.selected_unique_features[best_peak['id_number']] = (
                    # empCpd id, neutral_formula, ion_relation (will change not always anchor)
                    interim_id, V['neutral_formula'], best_peak.get('ion_relation', '')
                )


    def export_feature_tables(self, outfile='cmap_feature_table.tsv',
                                    _snr=5, _peak_shape=0.7, _cSelectivity=0.7):
        '''
        To export features tables:
        1) preferred table under `outdir`, after quality filtering by SNR, peak shape and chromatographic selectivity.
        2) full table under `outdir/export/`
        3) unique compound table under `outdir/export/`
        4) dependent on `target` extract option, a targeted_extraction table under `outdir`.
        Filtering parameters (_snr, _peak_shape, _cSelectivity) only apply to preferred table and unique_compound_table.
        Full table is filtered by initial peak detection parameters of lower snr and gaussian_shape.
        '''
        good_samples = [sample.name for sample in self.all_samples] 
        filtered_FeatureTable = self.CMAP.FeatureTable[good_samples]                       
        # non zero counts
        count = filtered_FeatureTable[filtered_FeatureTable>1].count(axis='columns')
        self.CMAP.FeatureTable['detection_counts'] = count

        use_cols = [ 'id_number', 'mz', 'rtime', 'rtime_left_base', 'rtime_right_base', 'parent_masstrack_id', 
                    'peak_area', 'cSelectivity', 'goodness_fitting', 'snr', 'detection_counts' ] + good_samples
        filtered_FeatureTable = self.CMAP.FeatureTable[use_cols]
        filtered_FeatureTable['mz'] = filtered_FeatureTable['mz'].round(4)
        filtered_FeatureTable['rtime'] = filtered_FeatureTable['rtime'].round(2)
        filtered_FeatureTable['rtime_left_base'] = filtered_FeatureTable['rtime_left_base'].round(2)
        filtered_FeatureTable['rtime_right_base'] = filtered_FeatureTable['rtime_right_base'].round(2)
        filtered_FeatureTable['cSelectivity'] = filtered_FeatureTable['cSelectivity'].round(2)
        filtered_FeatureTable['goodness_fitting'] = filtered_FeatureTable['goodness_fitting'].round(2)

        outfile = os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table'])
        filtered_FeatureTable.to_csv(outfile, index=False, sep="\t")
        print("\nFeature table (%d x %d) was written to %s." %(
                                filtered_FeatureTable.shape[0], self.number_of_samples, outfile))

        # extract targeted m/z features
        if 'target' in self.parameters and self.parameters['target']:  
            matched_list, _, target_unmapped = complete_mass_paired_mapping(
                filtered_FeatureTable['mz'].to_list(), self.parameters['target'], self.parameters['mz_tolerance_ppm']
            )
            print("\nIn targeted extraction, %d target mz values are not found in this dataset: " %len(target_unmapped))
            print('    ', [self.parameters['target'][ii] for ii in target_unmapped])
            matched_targets = [self.parameters['target'][ii[1]] for ii in matched_list]
            targeted_table = filtered_FeatureTable.iloc[[x[0] for x in matched_list], :]
            targeted_table.insert(0, "query_target", matched_targets)
            outfile = os.path.join(self.parameters['outdir'], 'targeted_extraction__'+self.parameters['output_feature_table'])
            targeted_table.to_csv(outfile, index=False, sep="\t")
            print("Targeted extraction Feature table (%d x %d) was written to %s.\n" %(
                                targeted_table.shape[0], self.number_of_samples, outfile))

        outfile = os.path.join(self.parameters['outdir'], 'preferred_'+self.parameters['output_feature_table'])
        # Some features can have all 0s, filtered here
        filtered_FeatureTable = filtered_FeatureTable[  filtered_FeatureTable['detection_counts'] > 0 ]
        filtered_FeatureTable = filtered_FeatureTable[  filtered_FeatureTable['snr']>_snr][
                                                        filtered_FeatureTable['goodness_fitting']>_peak_shape][
                                                        filtered_FeatureTable['cSelectivity']>_cSelectivity ]
        filtered_FeatureTable.to_csv(outfile, index=False, sep="\t")
        print("Filtered Feature table (%d x %d) was written to %s.\n" %(
                                filtered_FeatureTable.shape[0], self.number_of_samples, outfile))
        
        # in self.selected_unique_features: (empCpd id, neutral_formula, ion_relation)
        sel = [ii for ii in filtered_FeatureTable.index if filtered_FeatureTable['id_number'][ii] in
                                    self.selected_unique_features.keys()]
        unique_compound_table = filtered_FeatureTable.loc[sel, :]
        unique_compound_table.insert(3, "empCpd", [self.selected_unique_features[ii][0] for ii in unique_compound_table['id_number']])
        unique_compound_table.insert(4, "neutral_formula", [self.selected_unique_features[ii][1] for ii in unique_compound_table['id_number']])
        unique_compound_table.insert(5, "ion_relation", [self.selected_unique_features[ii][2] for ii in unique_compound_table['id_number']])
        
        outfile = os.path.join(self.parameters['outdir'], 'export', 'unique_compound__'+self.parameters['output_feature_table'])
        unique_compound_table.to_csv(outfile, index=False, sep="\t")
        print("Unique compound table (%d x %d) was written to %s.\n" %(
                            unique_compound_table.shape[0], self.number_of_samples, outfile))

        
    def export_log(self):
        '''
        Export project parameters to project.json,
        which is also used by asari viz.
        '''
        self.parameters['number_of_samples'] = self.number_of_samples
        self.parameters['number_scans'] = self.number_scans
        outfile = os.path.join(self.parameters['outdir'], 'project.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(self.parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)
