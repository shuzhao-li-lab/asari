'''
Functions for subcommand `annotate`
To test and solidify.

GC-HRMS annotation functions are added here.
This does not change default --anno for LC-MS data processing; but we will unify in next version.

Example use for GC-HRMS annotation:
asari annotate -i /Users/lish/li.proj/IndiPHARM/LookAhead/GCMS/lookahead_gcms_asari13_asari_project_311234834/export/full_Feature_table.tsv -o /Users/lish/li.proj/IndiPHARM/LookAhead/GCMS/ -j v16 --kovats /Users/lish/li.proj/asari_project/GCHRMS/KovatsIndex_alkanestandards_cuimc.tsv --db /Users/lish/li.proj/asari_project/GCHRMS/Resources/GCHRMS_Database_251217.msp --denovo F --workflow GC

LC: 
asari annotate -i /Users/lish/li.play/test16_v16_599432/preferred_Feature_table.tsv -o /Users/lish/li.play/test16_v16_599432/ -j v16lc --workflow LC --mode pos
'''

import os
import json
import time
import pickle
import importlib.resources as pkg_resources
from pathvalidate import sanitize_filename

from . import db
from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from .default_parameters import (extended_adducts, adduct_search_patterns_pos, 
                         adduct_search_patterns_neg, 
                         isotope_search_patterns)
from .utils import NpEncoder
from .tools.file_io import read_features_from_asari_table
from .tools.plot import mirror_plot
from .gcms import *

def annotate_project(infile, parameters):
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    subdir = '_'.join(['annotation', parameters['project_name'], ''.join(time_stamp)])
    outdir = os.path.join(parameters['outdir'], subdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"Annotation directory: {outdir}.\n\n")
    
    if parameters['workflow'] == "LC":
        ANN = LCMS_Annotation(parameters)
        ANN.parameters['outdir'] = outdir
        ANN.annotate_user_featuretable(infile)
        
    elif parameters['workflow'] == "GC":
        annotate_gcms_full(
            infile=infile,
            outdir=outdir,
            KovatsIndex=parameters['kovats'],
            database_file=parameters['db'],
            project_name_handle=parameters['project_name'],
            denovo=parameters['denovo'],
            ms2_tolerance_in_ppm=parameters['ms2_tolerance_in_ppm'], 
            ms2_tolerance_in_da=parameters['ms2_tolerance_in_da'], 
            
            ri_tolerance=parameters['ri_tolerance'],  
            score_cutoff_cosine=parameters['score_cutoff_cosine'], 
            score_cutoff_entropy=parameters['score_cutoff_entropy'],
            
            corr_cutoff=parameters['corr_cutoff'], 
            max_ri_delta=parameters['max_ri_delta'], 
            do_mirror_plot=parameters['do_mirror_plot'],
            max_core_features=parameters['max_core_features'],
        )
        
    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)

#
# main GC-HRMS function
#

def annotate_gcms_full(
    infile,
    outdir,
    KovatsIndex,
    database_file,
    project_name_handle='result',
    low_peak_filter_factor=100,
    ms2_tolerance_in_ppm=5,     # not used now
    ms2_tolerance_in_da=0.005, 
    ri_tolerance=50,
    cosine_penalty=1,           # don't change
    score_cutoff_cosine=0.5, 
    score_cutoff_entropy=0.4,
    corr_cutoff=0.6, 
    min_ri_delta=1,
    max_ri_delta=100, 
    hcl_distance_cut = 1,
    clustering_step_size = 2000,
    feature_distance_filter=None,
    do_mirror_plot=True,
    max_core_features=20000,
    denovo=False 
    ):
    '''Workflow using batch_lib_search_score. 
    
    low_peak_filter_factor : to exclude peaks in a spectrum if low_peak_filter_factor*peak < base_peak.
    
    '''
    # infile in asari format
    num_samples, list_features = read_features_from_asari_table(open(infile).read())
    ri_model = read_fit_KovatsIndex_rtime(KovatsIndex, sep='\t', frac=0.3)
    list_features = append_kovats_index(list_features, ri_model)
    list_features.sort(key=lambda x: x['peak_area'], reverse=True)
    dict_features = {f['id']: f for f in list_features}
    
    # load GCHRMS Database in MSP or JSON format
    list_lib_entries = reformat_gcms_lib( load_gcms_dbfile(database_file), 
                                     filter_factor=low_peak_filter_factor
                                     )
    dict_lib_entries = {e.id: e for e in list_lib_entries}
    print(f"Imported {len(list_lib_entries)} compound library entries.")
    
    feature_dataframe = pd.read_csv(infile, sep="\t", index_col=0)
    feature_dataframe = feature_dataframe.iloc[:, 10:]
    
    matched_results = batch_lib_search_score(
        list_lib_entries, list_features, dict_features, feature_dataframe,
        mz_tolerance_da=ms2_tolerance_in_da, ri_tolerance=ri_tolerance,
        cosine_penalty = cosine_penalty, corr_cutoff=corr_cutoff
    )
    list_empCpds, feature_anno_list = curate_batch_lib_search_result(
        matched_results, 
        mz_tolerance_da=ms2_tolerance_in_da,
        score_cutoff_cosine=score_cutoff_cosine, 
        score_cutoff_entropy=score_cutoff_entropy,
    )
    write_tsv_feature_anno(feature_anno_list, dict_features, dict_lib_entries, 
                        os.path.join(outdir, "Features_" + project_name_handle + '.tsv'))
    write_tsv_empCpd_anno(list_empCpds, dict_features, dict_lib_entries, 
                        os.path.join(outdir, "empCpds_" + project_name_handle + '.tsv'))
    print("\nDone targeted annotation.")
    print(f"Exported tsv results for {len(list_empCpds)} annotated empCpds and {len(set([f['feature'] for f in feature_anno_list]))} unique features.\n")
    
    if do_mirror_plot:
        print("Exporting PDF mirror plots..\n")
        path_mirrorplots = os.path.join(outdir, "mirrorplots/")
        os.makedirs(path_mirrorplots, exist_ok=True)
        for tt in list_empCpds:
            # print(tt['name'] )
            _score = f"\n(cosine score: {tt['cosine_score']:.2f}, entropy score: {tt['entropy_score']:.2f})"
            _title = tt['id'] + '__' + tt['name']   
            mirror_plot(tt['peaks_as_features'], 
                        # tt['peaks_in_lib'],   # filtered by unit mz match
                        dict_lib_entries[tt['lib_entry_id']].peaks,     # peaks from original lib, not filtered by unit mz match, to show the full spectrum.
                        match_tol=ms2_tolerance_in_da,
                        colors=["blue", "tab:red", "black"],
                        title=_title + _score, 
                        outfile=os.path.join(path_mirrorplots, sanitize_filename(_title)+'.pdf'))
    
    # Export list_empCpds JSON and MSP
    list_empCpds = serialize_annotated_empCpds(list_empCpds)
    json_output_file = os.path.join(outdir, project_name_handle+'_annotated_pseudospectra.json')
    with open(json_output_file, 'w', encoding='utf-8') as O:
        json.dump(list_empCpds, O, indent=2, ensure_ascii=False)
    json_pseudospectra_to_msp(list_empCpds, json_output_file.replace('.json', '.msp'))
    print(project_name_handle + "_annotated_pseudospectra was written to JSON and MSP formats.\n")

    # Deconvolution de novo
    if denovo:
        print("De novo construction of pseudospectra (deconvolution)...")
        core_features = set([f['feature'] for f in feature_anno_list])
        list_pseudospectra, core_features = iterative_build_pseudospectra_by_hcl(
            list_features, feature_dataframe, 
            init_core_features=core_features, 
            step_size=clustering_step_size,
            hcl_distance_cut = hcl_distance_cut,
            ri_tolerance=ri_tolerance,
            correlation_cut = corr_cutoff,
            max_core_features=max_core_features 
        )
        json_output_file = os.path.join(outdir, project_name_handle+'_denovo_pseudospectra.json')
        list_pseudospectra = [port_pseudospectrum_to_json(x, normalize_intensity=True) for x in list_pseudospectra]
        with open(json_output_file, 'w', encoding='utf-8') as O:
            json.dump(list_pseudospectra, O, indent=2, ensure_ascii=False)
        json_pseudospectra_to_msp(list_pseudospectra, json_output_file.replace('.json', '.msp'))
        print(f"Exported {len(list_pseudospectra)} de novo pseudospectra to JSON and MSP formats.\n")


#
# -----------------------------------------------------------------------------
#

class LCMS_Annotation:
    def __init__(self, parameters):
        self.parameters = parameters

    def annotate_user_featuretable(self, infile):
        '''
            Annotate a user supplied feature table. Used by asari subcommand `annotate`.

        Parameters
        ----------
        infile : str
            input feature table filepath, tab delimited file with first row as header, 
            first column m/z and 2nd column rtime.
        parameters : dict
            parameter dictionary passed from main.py, 
            which imports from default_parameters and updates the dict by user arguments.
        rtime_tolerance : int, optional, default: 2
            retention time tolerance to group adducts etc.

        Outputs
        -------
        two files in current directory, Feature_annotation.tsv and Annotated_empricalCompounds.json
    

        Note
        ----
        Annotate features via JMS (jms.dbStructures) and khipu.
        The pre-annotation step is khipu based emprical compound construction, followed by 
        three steps of annotation:

        1) Search known compound database, via neutral mass inferred by khipu

        2) Search singletons for a formula match

        3) Encapsulate remaining features in empCpd format, so that all are exported consistently.

        Export Feature_annotation as tsv, Annotated_empricalCompounds in both JSON and pickle.
        Reference databases can be pre-loaded. 
        Measured m/z values are calibrated to database based values (db_mass_calibrate).
        '''
        _n, list_features = read_features_from_asari_table(open(infile).read())
        for f in list_features:
                f['representative_intensity'] = f['peak_area']
        list_features.sort(key=lambda x: x['peak_area'], reverse=True)
        self.list_features = list_features
                
        self.load_annotation_db()
        self.db_mass_calibrate()

        # asari uses seconds for rt
        EED = ExperimentalEcpdDatabase(mode=self.parameters['mode'], 
                                       mz_tolerance_ppm=self.parameters['mz_tolerance_ppm'], 
                                       rt_tolerance=self.parameters['khipu_rtime_tolerance']
                                       )
        # passing patterns from .default_parameters
        if self.parameters['mode'] == 'pos':
            EED.adduct_patterns = adduct_search_patterns_pos
        else:
            EED.adduct_patterns = adduct_search_patterns_neg
        EED.isotope_search_patterns = isotope_search_patterns
        EED.extended_adducts = extended_adducts

        EED.build_from_list_peaks(list_features)
        # It takes three steps to take care of all features. First khipu organized empCpds
        EED.extend_empCpd_annotation(self.KCD)
        # Second, singletons that get a formula match in KCD
        EED.annotate_singletons(self.KCD)       
        # Third, the remaining features unmatched to anything (orphans). Exported for potential downstream work.
        EED.dict_empCpds = self.append_orphans_to_epmCpds(EED.dict_empCpds)

        self.export_peak_annotation(EED.dict_empCpds, self.KCD, 'Feature_annotation')

        # export JSON
        outfile = os.path.join(self.parameters['outdir'], 'Annotated_empiricalCompounds.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(EED.dict_empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)
        print(f"JSON version was written to {outfile}.\n\n")


    def load_annotation_db(self, src='hmdb4'):
        '''
        Load database of known compound using jms.dbStructures.knownCompoundDatabase.
        The compound tree is precomputed indexing.
        The `src` parameter is not used now, but placeholder to add more options later.

        Parameters
        ----------
        src: str, optional, default: hmdb4
            not used but can, in the future, dictate which database is used to generate annotations
              
        
        '''
        self.KCD = knownCompoundDatabase()
        self.KCD.mass_indexed_compounds = pickle.load( 
            pkg_resources.open_binary(db, 'mass_indexed_compounds.pickle') )
        self.KCD.emp_cpds_trees = pickle.load( 
            pkg_resources.open_binary(db, 'emp_cpds_trees.pickle') )

    def db_mass_calibrate(self, max_features=1000, required_calibrate_threshold=0.000002):
        '''
        Use KCD.evaluate_mass_accuracy_ratio to check systematic mass shift,
        which is calculated as the average ppm difference between measured m/z and theoretical values.
        If greater than required_calibrate_threshold (default 2 ppm), 
        calibrate m/z values for the whole experiment by updating self.CMAP.FeatureList.

        Parameters
        ---------
        required_calibrate_treshold: float, optional, default: 0.000002
            if the mass shift exceeds this value, mass correction will be applied. 

        Note
        ----
        Data format in good_reference_landmark_peaks: 
        [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...],
        where ref_id_num is index number of mass track in MassGrid.
        '''
    def db_mass_calibrate(self, max_features=1000, required_calibrate_threshold=0.000002):
        mz_landmarks = [f['mz'] for f in self.list_features[:max_features]]
        mass_accuracy_ratio = self.KCD.evaluate_mass_accuracy_ratio(
            mz_landmarks, mode=self.parameters['mode'], mz_tolerance_ppm=10)
        if mass_accuracy_ratio:
            if abs(mass_accuracy_ratio) > required_calibrate_threshold:
                print("Mass shift is greater than %2.1f ppm. Correction applied." 
                      %(required_calibrate_threshold*1000000))
                _correction = mass_accuracy_ratio + 1
                for F in self.list_features:
                    F['mz'] = F['mz'] / _correction
                    F['mz_corrected_by_division'] = _correction
        else:
            print("Mass accuracy check is skipped, too few mz_landmarks (%d) matched." 
                  %len(mz_landmarks))

    def append_orphans_to_epmCpds(self, dict_empCpds):
        '''
        This is the third step of feature annotation in self.annotate,
        to encapsulate features without annotation in empCpd format.
        Input via dict_empCpds and returns updated dict_empCpds.
        See also: annotate

        Parameters
        ----------
        dict_empCpds: dict
            a dictionary of empirical compounds in empCpd format

        '''
        all_feature_ids = []
        for _, V in dict_empCpds.items():
            all_feature_ids += [peak['id_number'] for peak in V['MS1_pseudo_Spectra']]

        orphans = [peak for peak in self.list_features if peak['id_number'] not in all_feature_ids]
        # will need better tracking of empCpd ID; already used in jms.dbStructures new_id_start = len(self.dict_empCpds) + 10000
        new_id_start = len(dict_empCpds) + 100000
        for peak in orphans:
            dict_empCpds[new_id_start] = {'interim_id': new_id_start,
                    'neutral_formula_mass': '', 
                    'neutral_formula': '',
                    'MS1_pseudo_Spectra': [peak],
                    'ion_relation': None,
                    'modification': None}
            new_id_start += 1

        return dict_empCpds

    def export_peak_annotation(self, dict_empCpds, KCD, export_file_name_prefix):
        '''
        Export feature annotation to tab delimited tsv file, where interim_id is empCpd id.
        
        Parameters
        ----------
        dict_empCpds : dict
            dictionary of empirical compounds, using interim_id as key, as seen in JMS.
        KCD : KnownCompoundDatabase instance
            the known compound database that was used in annotating the empirical compounds.
        export_file_name_prefix : str
            to used in output file name.
        '''
        s = "[peak]id_number\tmz\trtime\tapex(scan number)\t[EmpCpd]interim_id\
            \t[EmpCpd]ion_relation\tneutral_formula\tneutral_formula_mass\
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
                    peak['id_number'], peak['mz'], peak['rtime'], peak['apex'], 
                    V['interim_id'], peak.get('ion_relation', ''),
                    V['neutral_formula'], V['neutral_formula_mass'],
                    name_1st_guess, matched_DB_shorts, matched_DB_records]]) + "\n"

        outfile = os.path.join(self.parameters['outdir'], export_file_name_prefix + '.tsv')
        with open(outfile, encoding='utf-8', mode='w') as O:
            O.write(s)

        print("\nAnnotation of %d Empirical compounds was written to %s.\n\n" %(len(dict_empCpds), outfile))

    def select_unique_compound_features(self, dict_empCpds):
        '''
        Get unique feature by highest composite peak area per empirical compound. 
        One may consider alternatives to select the peak representing an empirical compound,
        e.g. by SNR or M+H, M-H ions. This is can be done separately on the exported files.

        Parameters
        ----------
        dict_empCpds : dict
            dictionary of empirical compounds, using interim_id as key, as seen in JMS.
            
        Not used now.
        
        
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
                    all_peaks = sorted([(peak['peak_area'], peak['goodness_fitting'], peak) 
                                        for peak in V['MS1_pseudo_Spectra']],
                                        reverse=True)
                    best_peak = all_peaks[0][2]
                except TypeError:
                    best_peak = V['MS1_pseudo_Spectra'][0]
                self.selected_unique_features[best_peak['id_number']] = (
                    # empCpd id, neutral_formula, ion_relation (will change not always anchor)
                    interim_id, V['neutral_formula'], best_peak.get('ion_relation', '')
                )



# merging with LCMS_Annotation class
def annotate_user_featuretable(infile, parameters):
    '''
    
    Annotate a user supplied feature table. Used by asari subcommand `annotate`.

    Parameters
    ----------
    infile : str
        input feature table filepath, tab delimited file with first row as header, 
        first column m/z and 2nd column rtime.
    parameters : dict
        parameter dictionary passed from main.py, 
        which imports from default_parameters and updates the dict by user arguments.
    rtime_tolerance : int, optional, default: 2
        retention time tolerance to group adducts etc.

    Outputs
    -------
    _PreAnnotated_list_khipus.json
    
    # currently redundant w/ --anno, but to be unified in next version.
    '''
    from khipu.extended import peaklist_to_khipu_list 
    
    mode = parameters['mode']
    if mode == 'pos':
        adduct_search_patterns = adduct_search_patterns_pos
    else:
        adduct_search_patterns = adduct_search_patterns_neg
    _n, list_features = read_features_from_asari_table(open(infile).read())
    for f in list_features:
            f['representative_intensity'] = f['peak_area']
    list_khipus, all_assigned_fids = peaklist_to_khipu_list(
                            list_features, 
                            isotope_search_patterns=isotope_search_patterns, 
                            adduct_search_patterns=adduct_search_patterns,
                            extended_adducts=extended_adducts, 
                            mz_tolerance_ppm=parameters['mz_tolerance_ppm'],
                            rt_tolerance=parameters['khipu_rtime_tolerance'],
                            mode=mode,
                            charges=[1, 2, 3],
                            )
    #  exporting JSON
    json_output_file = os.path.join(parameters['outdir'], 
                                    parameters['project_name']+'_PreAnnotated_list_khipus.json')
    with open(json_output_file, 'w', encoding='utf-8') as f:
        json.dump(list_khipus, f, cls=NpEncoder, ensure_ascii=False, indent=2)
        
    
    


def annoate_by_standards(list_features, cpdLib):
    matched = match_features.list_match_lcms_features(
        list_features, cpdLib,  mz_ppm=5,  rt_tolerance=30
    )
    
    dict_f = {feature['id']: feature for feature in list_features}
    dict_lib = {feature['id']: feature for feature in cpdLib}
    
    
    return matched 



def get_concise_annotation(feature_id, matched_list, feature_dict, lib_dict):
    '''Remove duplicate annotations for the same feature, and return the concise list of fields.
    '''
    feature = feature_dict[feature_id]
    result = {
        'id': feature_id,
        'mz': feature['mz'],
        'rtime': feature['rtime'],
    }
    matched_auth_libs = [x for x in matched_list if x.startswith('v2r2024')]
    matched_csm_libs = [x for x in matched_list if x.startswith('r1_')]
    # for authentic cpd lib and csm matches, keep one each
    if matched_auth_libs:
        lib_id = matched_auth_libs[0]
        lib_info = lib_dict[lib_id]
        result.update({
            'lib_id': lib_id,
            'lib_name': lib_info['name'],
            'lib_mz': lib_info['mz'],
            'lib_rtime': lib_info['rtime'],
            'lib_identifier': lib_info['identifier'],
            'lib_ion': lib_info['ion'],
            'lib_isotope': lib_info['isotope'],
        })
    else:
        result.update({
            'lib_id': '',
            'lib_name': '',
            'lib_mz': '',
            'lib_rtime': '',
            'lib_identifier': '',
            'lib_ion': '',
            'lib_isotope': '',
        })
    if matched_csm_libs:
        csm_id = matched_csm_libs[0]
        result.update({
            'CSMF_ID': csm_id,
            'CSM_ion': lib_dict[csm_id]['ion_csm'],
            'CSM_top_recommendation_name': lib_dict[csm_id]['top_recommendation_name'],
            'CSM_top_recommendation_score': lib_dict[csm_id]['top_recommendation_score'],
            'CSM_HMDB': lib_dict[csm_id]['HMDB'],
        })
    else:
        result.update({
            'CSMF_ID': '',
            'CSM_ion': '',
            'CSM_top_recommendation_name': '',
            'CSM_top_recommendation_score': '',
            'CSM_HMDB': '',
        })
    return result


def export_combined_anno_table(matched, feature_dict, lib_dict, outfile):
    header = ['id', 'mz', 'rtime', 'lib_id', 'lib_name', 'lib_mz', 'lib_rtime', 'lib_identifier',  
            'lib_ion', 'lib_isotope', 'CSMF_ID', 'CSM_ion', 
            'CSM_top_recommendation_name', 'CSM_top_recommendation_score',
            'CSM_HMDB']
    s = '\t'.join(header) + '\n'
    for feature_id, matched_list in matched.items():
        concise_info = get_concise_annotation(feature_id, matched_list, feature_dict, lib_dict)
        s += '\t'.join([str(concise_info[h]) for h in header]) + '\n'
        
    with open(outfile, 'w') as f:
        f.write(s)

def export_combined_anno_json(matched, feature_dict, lib_dict, outfile):
    '''
    '''
    json_output = []
    for feature_id, matched_list in matched.items():
        json_output.append(
            {feature_id: {'Feature': feature_dict[feature_id], 
                          'Matched_libs': [lib_dict[lib_id] for lib_id in matched_list]}}
        )
    with open(outfile, 'w') as f:
        json.dump(json_output, f, indent=4)


