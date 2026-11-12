"""Experiment class, managing samples, composite map, annotation and export in a project."""

import os
import sys
import json
import pickle
import numpy as np
import importlib.resources as pkg_resources  # Python 3.11+, importlib_resources for older versions

from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

# to move out of experiment to separate annotation module
from .default_parameters import adduct_search_patterns, \
    adduct_search_patterns_neg, isotope_search_patterns, extended_adducts, \
    readme_doc_str

from functools import partial

from .constructors import CompositeMap, _cluster_track_ms2_worker
from .utils import NpEncoder, bulk_process
from . import db
from .tools.ms2 import json_ms2_to_msp, group_ms2_by_precursor

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


class ext_Experiment:
    """
    Similar to metDataModel.core.Experiment with preprocessing methods.

    This encapsulates a set of LC-MS files using the same experimental method
    (chromatography and ionization) to be processed together.

    E.g., data from postive ESI and negative ESI should not in the same ext_Experiment instance.

    This class has annotation and export functions.

    Default asari work flow is in `ext_Experiment.process_all`.
    """
    def __init__(self, sample_registry, parameters):
        """
        This is the overall container for all data in an experiment/project.

        Parameters
        ----------
        sample_registry : dict
            dictionary of sample and selected data after mass track extraction.
            The bulk mass track data are not kept in memory unless specified so.
        parameters : dict
            processing parameters passed from main.py.

        Updates
        -------
        Major class attributes including self.number_of_samples, number_scans, reference_sample_id.
        """
        self.sample_registry = sample_registry
        self.valid_sample_ids = self.get_valid_sample_ids()
        self.number_of_samples = len(self.valid_sample_ids)
        self.number_scans = self.get_max_scan_number()
        self.all_samples = self.all_sample_instances = []

        self.parameters = parameters
        self.output_dir = parameters['outdir']
        self.mode = parameters['mode']
        self.mz_tolerance_ppm = self.parameters['mz_tolerance_ppm']
        self.check_isotope_ratio = self.parameters['check_isotope_ratio']
        self.database_mode = parameters['database_mode']
        self.reference_sample_id = self.get_reference_sample_id()

    def get_reference_sample_id(self):
        """
        get_reference_sample_id either by user specification, or
        using the sample of most number_anchor_mz_pairs.
        This assumes the sample of most good m/z values has a good coverage of features.
        """
        if self.parameters['reference']:
            # match file name; k is sm['sample_id']
            for k,v in self.sample_registry.items():
                if os.path.basename(self.parameters['reference']) == os.path.basename(v['input_file']):
                    return k
                elif os.path.basename(self.parameters['reference']) + ".mzML" == os.path.basename(v['input_file']):
                    return k

        elif self.sample_registry:
            L = [(v['number_anchor_mz_pairs'], v['sample_id']) for v in self.sample_registry.values() if v.get('number_anchor_mz_pairs') is not None]
            if not L:       # e.g. MSMS workflow, no mass tracks
                return None
            L.sort(reverse=True)
            ref = self.sample_registry[L[0][1]]
            self.parameters['reference'] = ref['input_file']
            print("\n    The reference sample is:\n    ||* %s *||\n" %ref['name'])
            print("Max reference retention time is %4.2f at scan number %d.\n" %(
                max(ref['list_retention_time']), ref['max_scan_number']))
            return ref['sample_id']
        else:
            return None

    def get_valid_sample_ids(self):
        """
        Get valid sample ids, as some samples may not be extracted successfully.
        """
        return [k for k,v in self.sample_registry.items() if v['status:eic'] == 'passed']

    def determine_acquisition_order(self):
        """Return (sample id, acquisition time) tuples, sorted by acquisition time."""
        sample_order_by_timestamp = [(k, v['acquisition_time']) for k,v in self.sample_registry.items()]
        try:
            # fix this better in future
            return sorted(sample_order_by_timestamp, key=lambda x: x[1])
        except:
            return sorted(sample_order_by_timestamp, key=lambda x: x[0])

    def get_max_scan_number(self):
        """
        Return max scan number among samples, or None if no valid sample.
        self.sample_registry: a dict that maps sample IDs to sample data
        """
        scan_numbers = [self.sample_registry[k].get('max_scan_number') for k in self.valid_sample_ids]
        scan_numbers = [x for x in scan_numbers if x is not None]
        if scan_numbers:
            return max(scan_numbers) + 1
        else:
            return None

    def process_all_LC(self):
        """
        This is the default asari workflow for LC-MS and GC-MS procssing; annotation is separate after processing.

        1. Build MassGrid, using either pairwise (small study) or clustering method.
           Choose one reference from all samples for the largest number of landmark m/z tracks.
        2. RT alignment via a LOWESS function, using selective landmark peaks.
        3. Build composite elution profile (composite_mass_tracks),
           by cumulative sum of mass tracks from all samples after RT correction.
        4. Global peak detection is performed on each composite massTrack.
        5. Mapping global peaks (i.e. features) back to all samples and extract sample specific peak areas.
           This completes the FeatureTable.

        Updates
        -------
        self.CMAP as instance of CompositeMap, and MassGrid, composite map and features within.
        """
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid()
        if not self.parameters['rt_align_on']:
            self.CMAP.mock_rentention_alignment()
        self.CMAP.build_composite_tracks()
        self.CMAP.global_peak_detection()

    def process_all_DIMS(self):
        """
        For DI-MS, no RT alignment or peak detection.
        Only MassGrid construction then report max intensity per mass track per sample.

        """
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid()
        self.CMAP.get_DIMS_feature_table()

    def process_all_LCMSMS(self):
        """
        The LC-MS/MS workflow constructs massGrid and CMAP as LC-MS processing,
        then assign all MS/MS spectra to mass tracks.
        The MS/MS spectra per mass track are clustered based on RT and similarity.

        The MS1 data may be sparse in LC-MS/MS, limiting elution peak detection and resolution of close isomers.
        Since MS/MS collection is not consistent cross samples,
        full integrated with MS1 data in the same experiment is not always feasible.
        Even if we use MS1 peak height not area, MS1 and MS/MS features are not consistently measured cross samples.
        Therefore, quantitative information cross samples is only usable for a subset of MS1 features.

        The MS1 feature table and MS/MS clusters are reported separately, sharing masstrack IDs.
        The MS/MS clusters are searched against spectral libraries seprately and result can be used to annotate
        the MS1 feature table here or different acquisitions.
        """
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid()
        if not self.parameters['rt_align_on']:
            self.CMAP.mock_rentention_alignment()
        # this includes MS2 in mass tracks
        self.CMAP.build_composite_tracks()
        # cluster and clean up MS2; else MS2 spectra are kept as is and exported per input file
        if self.parameters.get('cluster_ms2', True):
            self.CMAP.cluster_ms2_spectra(
                # similarity_function=self.parameters.get('ms2_similarity_function', cosine_similarity),
                rt_gap=self.parameters.get('ms2_rt_gap', 5),
                mz_tolerance=self.parameters.get('ms2_mz_tolerance', 0.01),
                distance_threshold=self.parameters.get('ms2_distance_threshold', 0.8),
                method=self.parameters.get('ms2_cluster_method', 'greedy'),
                similarity_threshold=self.parameters.get('ms2_similarity_threshold', 0.7),
                score=self.parameters.get('ms2_cluster_score', 'entropy')
            )
        # MS1 peaks
        self.CMAP.global_peak_detection()

    def process_all_MSMS(self):
        """
        The MS/MS workflow extracts MS2 spectra per sample. 
        If cluster_ms2 == True, cluster all MS2 spectra and write clustered spectra to one JSON file and one MSP file;
        else: save one MSP file per input sample.

        MS1 data are ignored; there is no alignment of m/z or retention time.
        MS2 spectra are extracted in workflow_setup_MSMS and kept in self.sample_registry[sid]['ms2_spectra'].
        For clustering, spectra of all samples are grouped by precursor m/z (within mz_tolerance_ppm),
        then each precursor group is clustered on RT and spectral similarity, same as in the LCMSMS workflow.

        Updates
        -------
        self.ms2_precursor_groups : {id: {'id_number', 'mz', 'ms2_spectra'}, ...}, if cluster_ms2.
        """
        self.dropped_sample_names = [v['name'] for k,v in self.sample_registry.items()
                                     if k not in self.valid_sample_ids]
        all_spectra = []
        for k in self.valid_sample_ids:
            all_spectra += self.sample_registry[k]['ms2_spectra']
        print("\nTotal %d MS2 spectra from %d samples.\n" %(len(all_spectra), self.number_of_samples))

        if self.parameters.get('cluster_ms2', True):
            groups = group_ms2_by_precursor(all_spectra, self.mz_tolerance_ppm)
            print("Grouped into %d precursor m/z groups. Clustering MS/MS spectra (%s) ...\n" %(
                len(groups), self.parameters.get('ms2_cluster_method', 'greedy')))
            self.ms2_precursor_groups = {
                ii: {'id_number': ii,
                     'mz': float(np.median([spec['precursor_mz'] for spec in G])),
                     'ms2_spectra': G}
                for ii, G in enumerate(groups)
            }
            worker = partial(_cluster_track_ms2_worker,
                             similarity_function=None,
                             rt_gap=self.parameters.get('ms2_rt_gap', 5),
                             mz_tolerance=self.parameters.get('ms2_mz_tolerance', 0.01),
                             distance_threshold=self.parameters.get('ms2_distance_threshold', 0.8),
                             method=self.parameters.get('ms2_cluster_method', 'greedy'),
                             similarity_threshold=self.parameters.get('ms2_similarity_threshold', 0.7),
                             score=self.parameters.get('ms2_cluster_score', 'entropy'))
            self.ms2_precursor_groups = dict(bulk_process(worker, list(self.ms2_precursor_groups.items()),
                                                          jobs_per_worker=self.parameters['multicores']))

    def export_all(self, diagnosis=False):
        """
        Export all files.
        Annotation of features to empirical compounds is done here.

        Parameters
        ----------
        diagnosis: bool, optional, default: True
            if true, export CMAP pickle for diagnosis and viz
        """
        if self.parameters['workflow'] in ["LC", "GC"]:
            self.CMAP.MassGrid.to_csv(
                os.path.join(self.parameters['outdir'], 'export', self.parameters['mass_grid_mapping']) )
            if diagnosis:
                self.export_CMAP_pickle()

            self.export_feature_tables()
            self.export_log()
            self.export_readme()

        elif self.parameters['workflow'] == "DIMS":
            self.CMAP.MassGrid.to_csv(
                os.path.join(self.parameters['outdir'], 'export', self.parameters['mass_grid_mapping']) )
            if diagnosis:
                # add later
                pass
            self.export_feature_tables()
            self.export_log()
            self.export_readme()

        elif self.parameters['workflow'] == "LCMSMS":
            # annotation will be separate step
            self.CMAP.MassGrid.to_csv(
                os.path.join(self.parameters['outdir'], 'export', self.parameters['mass_grid_mapping']) )

            self.export_feature_tables()
            if self.parameters.get('cluster_ms2', True):
                self.export_clustered_ms2()
            else:
                self.export_ms2_per_sample()

            self.export_log()
            self.export_readme()

        elif self.parameters['workflow'] == "MSMS":
            if self.parameters.get('cluster_ms2', True):
                self.export_clustered_ms2(self.ms2_precursor_groups,
                                          track_label='precursor_group', id_prefix='precursor')
            else:
                dict_sample_spectra = {
                    self.sample_registry[k]['name']: [(spec, None) for spec in self.sample_registry[k]['ms2_spectra']]
                    for k in self.valid_sample_ids
                }
                self.export_ms2_per_sample(dict_sample_spectra)

            self.export_log()
            self.export_readme()

    @staticmethod
    def _format_ms2_spectrum(spec, track, spec_id, track_label='parent_masstrack'):
        """
        Round values and add parent track info to an MS/MS spectrum (in place).
        track : composite mass track (LCMSMS) or precursor group (MSMS); None to add no track info.
        track_label : prefix of the track fields, e.g. parent_masstrack_id, parent_masstrack_mz.
        """
        if track is not None:
            spec[track_label + '_id'] = track['id_number']
            spec[track_label + '_mz'] = round(float(track['mz']), 4)
        spec['id'] = spec_id
        spec['precursor_mz'] = round(spec['precursor_mz'], 4)
        spec['rtime'] = round(spec['rtime'], 2)
        spec['peaks'] = [[round(mz, 4), intensity] for mz, intensity in spec['peaks']]
        return spec

    def export_clustered_ms2(self, tracks=None, track_label='parent_masstrack', id_prefix='track'):
        """
        Export clustered MS/MS spectra (cluster representatives per composite mass track,
        or per precursor group in MSMS workflow) to export/ms2_spectra.json and export/ms2_spectra.msp.

        tracks : dict of tracks with clustered 'ms2_spectra'; default self.CMAP.composite_mass_tracks.
        """
        if tracks is None:
            tracks = self.CMAP.composite_mass_tracks
        _list_ms2_spectra = []
        for mt in tracks.values():
            # after clustering, the ms2_spectra in each mass track are the cluster representatives
            for ii, spec in enumerate(mt['ms2_spectra']):
                _list_ms2_spectra.append(
                    self._format_ms2_spectrum(spec, mt, f"{id_prefix}_{mt['mz']:.4f}_clus{ii}", track_label))

        outfile = os.path.join(self.parameters['outdir'], 'export', 'ms2_spectra')
        with open(outfile + '.json', 'w', encoding='utf-8') as f:
            json.dump(_list_ms2_spectra, f, cls=NpEncoder, ensure_ascii=False)
        json_ms2_to_msp(_list_ms2_spectra, outfile + '.msp')
        print("Exported %d clustered MS/MS spectra." %len(_list_ms2_spectra))

    def export_ms2_per_sample(self, dict_sample_spectra=None):
        """
        Export unclustered MS/MS spectra per input file, to export/ms2_spectra_<sample>.msp
        (MSP only, to save disk space).
        Each spectrum carries precursor_mz, rtime (seconds, as in input file),
        and the composite mass track it is assigned to (parent_masstrack_id, parent_masstrack_mz).
        Only MS/MS spectra mapped to a composite mass track are included.

        dict_sample_spectra : {sample: [(spectrum, track or None), ...]};
            default collects spectra from self.CMAP.composite_mass_tracks (LCMSMS).
            MSMS workflow passes all spectra per sample with no track.
        """
        if dict_sample_spectra is None:
            dict_sample_spectra = {}
            for mt in self.CMAP.composite_mass_tracks.values():
                for spec in mt['ms2_spectra']:
                    dict_sample_spectra.setdefault(spec.get('sample', 'unknown'), []).append((spec, mt))

        for sample, specs in dict_sample_spectra.items():
            specs.sort(key=lambda x: (x[0]['rtime'], x[0]['precursor_mz']))
            _list_ms2_spectra = [
                self._format_ms2_spectrum(spec, mt, f"{sample}_ms2_{ii}")
                for ii, (spec, mt) in enumerate(specs)
            ]
            json_ms2_to_msp(_list_ms2_spectra,
                            os.path.join(self.parameters['outdir'], 'export', f"ms2_spectra_{sample}.msp"))
        print("Exported unclustered MS/MS spectra of %d samples." %len(dict_sample_spectra))

    def export_CMAP_pickle(self):
        """
        Export main CMAP data and MassGrid to pickle, which can be used for visual data exploration.

        Included in exported pickle:{
            '_number_of_samples_': self.CMAP._number_of_samples_,
            'rt_length': self.CMAP.rt_length,
            'rt_reference_landmarks': [p['apex']
                                       for p in self.CMAP.good_reference_landmark_peaks],
            'rt_records': [sample.get_rt_calibration_records()
                                        for sample in self.all_samples
                                        ],
            'dict_scan_rtime': self.CMAP.dict_scan_rtime,
            'list_mass_tracks': self.CMAP.composite_mass_tracks,
            'MassGrid': dict(self.CMAP.MassGrid),}

        rt_records includes for each sample: {
            'sample_id': self.sample_id,
            'name': self.name,
            'rt_landmarks': self.rt_landmarks,
            'reverse_rt_cal_dict': self.reverse_rt_cal_dict,
        }

        Note
        ----
            RT calibration is exported to include sample.reverse_rt_cal_dict,
            i.e. {key=reference scan number, value=sample specific scan number}.
            May add data throttle in the future. The file cmap.pickle can get big.
        """
        _export = {
            '_number_of_samples_': self.CMAP._number_of_samples_,
            'rt_length': self.CMAP.rt_length,
            'rt_reference_landmarks': [p['apex']
                                       for p in self.CMAP.good_reference_landmark_peaks],
            'rt_records': [sample.get_rt_calibration_records()
                                        for sample in self.all_samples
                                        ],
            'dict_scan_rtime': self.CMAP.dict_scan_rtime,
            'list_mass_tracks': self.CMAP.composite_mass_tracks,
            'MassGrid': dict(self.CMAP.MassGrid),
        }
        outfile = os.path.join(self.parameters['outdir'], 'export', 'cmap.pickle')
        with open(outfile, 'wb') as f:
            pickle.dump(_export, f, pickle.HIGHEST_PROTOCOL)

    def export_feature_tables(self, _snr=2, _peak_shape=0.7, _cSelectivity=0.7):
        """
        To export multiple features tables.
        Filtering parameters (_snr, _peak_shape, _cSelectivity) only apply to preferred table.
        Full table is filtered by initial peak detection parameters,
        which contains lower values of snr and gaussian_shape.

        Parameters
        ----------
        _snr: float, optional, default: 2
            signal to noise ratio, peaks must have SNR above this value to be a preferred feature
        _peak_shape: float, optional, default: 0.7
            the goodness fitting for a peak must be above this value to be a preferred feature
        _cSelectivity: float, optional, default: 0.7
            the cSelectivity of a peak must be above this value to be a preferred feature

        Outputs
        -------
        Multiple features tables under output directory:
        1. preferred table under `outdir`, after quality filtering
           by SNR, peak shape and chromatographic selectivity.
        2. full table under `outdir/export/`
        """
        if self.parameters['drop_unaligned_samples']:
            good_samples = [sample.name for sample in self.all_samples if sample.is_rt_aligned]
        else:
            good_samples = [sample.name for sample in self.all_samples]
        filtered_FeatureTable = self.CMAP.FeatureTable[good_samples]          # this fixes order of samples
        self.number_of_samples = number_of_samples = len(good_samples)
        self.dropped_sample_names = [sample.name for sample in self.all_samples
                                     if sample.name not in good_samples]     # samples dropped due to alignment or quality
        # non zero counts
        count = filtered_FeatureTable[filtered_FeatureTable>1].count(axis='columns')
        self.CMAP.FeatureTable['detection_counts'] = count

        use_cols = [ 'id_number', 'mz', 'rtime', 'rtime_left_base', 'rtime_right_base', 'parent_masstrack_id',
                    'peak_area', 'cSelectivity', 'goodness_fitting', 'snr', 'detection_counts' ] + good_samples
        filtered_FeatureTable = self.CMAP.FeatureTable[use_cols].copy()
        # sample intensities as integers, to save disk space in .tsv
        filtered_FeatureTable[good_samples] = filtered_FeatureTable[good_samples].fillna(0).round().astype(np.int64)
        filtered_FeatureTable['mz'] = filtered_FeatureTable['mz'].round(4)
        filtered_FeatureTable['rtime'] = filtered_FeatureTable['rtime'].round(2)
        filtered_FeatureTable['rtime_left_base'] = filtered_FeatureTable['rtime_left_base'].round(2)
        filtered_FeatureTable['rtime_right_base'] = filtered_FeatureTable['rtime_right_base'].round(2)
        filtered_FeatureTable['cSelectivity'] = filtered_FeatureTable['cSelectivity'].round(2)
        filtered_FeatureTable['goodness_fitting'] = filtered_FeatureTable['goodness_fitting'].round(2)
        if self.parameters['workflow'] == "GC" and self.parameters.get('kovats'):
            filtered_FeatureTable = self.add_retention_index(filtered_FeatureTable)

        outfile = os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table'])
        filtered_FeatureTable.to_csv(outfile, index=False, sep="\t")
        print("\nFeature table (%d x %d) was written to %s." %(
                                filtered_FeatureTable.shape[0], number_of_samples, outfile))

        outfile = os.path.join(self.parameters['outdir'], 'preferred_'+self.parameters['output_feature_table'])
        # Some features can have all 0s, filtered here
        filtered_FeatureTable = filtered_FeatureTable[  filtered_FeatureTable['detection_counts'] > 0 ]
        filtered_FeatureTable = filtered_FeatureTable[  filtered_FeatureTable['snr']>_snr][
                                                        filtered_FeatureTable['goodness_fitting']>_peak_shape][
                                                        filtered_FeatureTable['cSelectivity']>_cSelectivity ]
        filtered_FeatureTable.to_csv(outfile, index=False, sep="\t")
        print("\nFiltered Feature table (%d x %d) was written to %s.\n" %(
                                filtered_FeatureTable.shape[0], number_of_samples, outfile))


    def add_retention_index(self, feature_table):
        """
        GC workflow: add a column 'RI' (Kovats retention index) right after 'detection_counts',
        so that sample columns start one column later (data_start_col=12 in gc_functions).
        RI is interpolated from feature rtime (seconds) by a lowess fit of the Kovats standards
        file parameters['kovats'] (columns: Kovats index, measured time in minutes).
        Features outside the range of the standards get the RI at the nearest end.
        """
        from .gc_functions import read_fit_KovatsIndex_rtime
        ri_model = read_fit_KovatsIndex_rtime(self.parameters['kovats'], frac=0.3)
        ri = np.interp(self.CMAP.FeatureTable.loc[feature_table.index, 'rtime'].values / 60,
                       ri_model[:, 0], ri_model[:, 1])
        feature_table = feature_table.copy()
        feature_table.insert(feature_table.columns.get_loc('detection_counts') + 1, 'RI', ri.round(2))
        print(f"Added retention index (RI) using Kovats standards in {self.parameters['kovats']}.")
        return feature_table

    def export_log(self):
        """
        Export project parameters to project.json,
        which is also used by asari viz.
        """
        self.parameters['number_of_samples'] = self.number_of_samples
        self.parameters['number_scans'] = self.number_scans
        # Record of dropped samples. Right place in log, wrong name in parameters
        self.parameters['dropped_samples'] = self.dropped_sample_names
        outfile = os.path.join(self.parameters['outdir'], 'project.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(self.parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)

    def export_readme(self):
        """
        Export a REAME.txt file as simple instruction to end users.
        """
        outfile = os.path.join(self.parameters['outdir'], 'README.txt')
        with open(outfile, 'w', encoding='utf-8') as f:
            f.write(readme_doc_str)


    #
    # to move out core Asari after 2.0
    #

    def load_annotation_db(self, src='hmdb4'):
        """
        Load database of known compound using jms.dbStructures.knownCompoundDatabase.
        The compound tree is precomputed indexing.
        The `src` parameter is not used now, but placeholder to add more options later.

        Parameters
        ----------
        src: str, optional, default: hmdb4
            not used but can, in the future, dictate which database is used to generate annotations


        To move annotation out of core package.

        as for GC and LCMSMS


        """
        self.KCD = knownCompoundDatabase()
        self.KCD.mass_indexed_compounds = pickle.load(
            pkg_resources.open_binary(db, 'mass_indexed_compounds.pickle') )
        self.KCD.emp_cpds_trees = pickle.load(
            pkg_resources.open_binary(db, 'emp_cpds_trees.pickle') )
