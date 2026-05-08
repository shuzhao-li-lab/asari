'''
Functions for subcommand `annotate`
To test and solidify.

GC-HRMS annotation functions are added here.
This does not change default --anno for LC-MS data processing; but we will unify in next version.

Example use:
asari annotate -i /Users/lish/li.proj/IndiPHARM/LookAhead/GCMS/lookahead_gcms_asari13_asari_project_311234834/export/full_Feature_table.tsv -o /Users/lish/li.proj/IndiPHARM/LookAhead/GCMS/ -j v16 --kovats /Users/lish/li.proj/asari_project/GCHRMS/KovatsIndex_alkanestandards_cuimc.tsv --db /Users/lish/li.proj/asari_project/GCHRMS/Resources/GCHRMS_Database_251217.msp --denovo F --workflow GC

'''
import os
import json
import time

from khipu.extended import peaklist_to_khipu_list 
#, export_empCpd_khipu_list
# When use relative import, sys.path.append errors out
from .default_parameters import (extended_adducts, adduct_search_patterns_pos, 
                         adduct_search_patterns_neg, 
                         isotope_search_patterns)
from .utils import NpEncoder

from .gcms import *

def annotate_project(infile, parameters):
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    subdir = 'annotation_' + ''.join(time_stamp)
    outdir = os.path.join(parameters['outdir'], subdir)
    os.mkdir(outdir)
    
    if parameters['workflow'] == "LC":
        annotate_user_featuretable(
            infile, 
            parameters, 
            
        )
    elif parameters['workflow'] == "GC":
        annotate_gcms_full(
            infile=infile,
            outdir=outdir,
            KovatsIndex=parameters['kovats'],
            database_file=parameters['db'],
            project_name_handle=parameters['project_name'],
            denovo=parameters['denovo'],
            # more paras
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
    ms2_tolerance_in_ppm=5, 
    ms2_tolerance_in_da=0.005, 
    ri_tolerance=50,
    cosine_penalty=1, 
    score_cutoff_cosine=0.5, 
    score_cutoff_entropy=0.4,
    corr_cutoff=0.7, 
    min_ri_delta=1,
    max_ri_delta=100, 
    hcl_distance_cut = 1,
    clustering_step_size = 2000,
    feature_distance_filter=None,
    do_mirror_plot=True,
    max_core_features=20000,
    denovo=False
    ):
    '''
    Full workflow of annotating GC-HRMS. 
    - Step 1a. Match base peak first; 
    - Step 1b. construct pseudospectra around base match.
    - Step 2a. Get top N features to do RI-coupled HCL, cut high to get seed clusters (de novo).
    - Step 2b. Add additional features to de novo clusters.
    - Step 2c. Iterate 2a and 2b until max features are accounted for.
    
    infile : input feature table, in asari format.
    KovatsIndex : .tsv file to map Kovats Index to retention time. 
        A regression function is derived to convert RT to RI for all features.
        This is not optional here. But users can use the notebook version to accommodate their own data formats.
    
    database_json : GCHRMS database converted to JSON, as MSP input may not be consistent. 
    
    ri_tolerance can be more stringent than max_ri_delta, which is used in penalty function.
    
    '''
    # infile in asari format
    num_samples, list_features = read_features_from_asari_table(open(infile).read())
    ri_model = read_fit_KovatsIndex_rtime(KovatsIndex, sep='\t', frac=0.3)
    list_features = append_kovats_index(list_features, ri_model)
    list_features.sort(key=lambda x: x['peak_area'], reverse=True)
    dict_features = {f['id']: f for f in list_features}
    
    # load GCHRMS Database in MSP or JSON format
    list_lib_entries = reformat_gcms_lib( load_gcms_dbfile(database_file), 
                                     filter_factor=low_peak_filter_factor)
    dict_lib_entries = {e.inchikey: e for e in list_lib_entries}
        
    feature_dataframe = pd.read_csv(infile, sep="\t", index_col=0)
    feature_dataframe = feature_dataframe.iloc[:, 10:]
    
    matched_results = batch_lib_search_by_basepeaks(
        list_lib_entries, list_features, feature_dataframe, 
            ms2_tolerance_in_ppm=ms2_tolerance_in_ppm, 
            ms2_tolerance_in_da=ms2_tolerance_in_da, 
            ri_tolerance=ri_tolerance,
            cosine_penalty = cosine_penalty,
            min_ri_delta=min_ri_delta,
            max_ri_delta=max_ri_delta,
            low_peak_filter_factor=low_peak_filter_factor,
            feature_distance_filter=feature_distance_filter
    )
    list_empCpds, feature_anno_list = export_feature_annotation_bybasepeaksearch(
        matched_results, 
        feature_dataframe,
        score_cutoff_cosine=score_cutoff_cosine, 
        score_cutoff_entropy=score_cutoff_entropy,
        corr_cutoff=corr_cutoff, 
    )
    write_tsv_feature_anno(feature_anno_list, dict_features, dict_lib_entries, 
                        os.path.join(outdir, "Features_" + project_name_handle + '.tsv'))
    write_tsv_empCpd_anno(list_empCpds, dict_features, dict_lib_entries, 
                        os.path.join(outdir, "empCpds_" + project_name_handle + '.tsv'))
    print("Exported tsv results for annotated empCpds and features.")
    
    if do_mirror_plot:
        print("Exportoing PDF mirror plots..")
        path_mirrorplots = os.path.join(outdir, "mirrorplots/")
        os.makedirs(path_mirrorplots, exist_ok=True)
        for tt in list_empCpds:
            # print(tt['name'] )
            _score = f" (cosine score: {tt['cosine_score']:.2f})"
            _title = tt['id'] + '__' + tt['name']   
            mirror_plot(tt['peaks_as_features'], tt['peaks_in_lib'], 
                        title=_title + _score, 
                        outfile=os.path.join(path_mirrorplots, _title+'.pdf'))
    
    # Deconvolution de novo
    if denovo:
        print("De novo construction of pseudospectra (deconvolution)...")
        core_features = set([f['feature'] for f in feature_anno_list if f['is_core']])
        list_pseudospectra, core_features = iterative_build_pseudospectra_by_hcl(
            list_features, feature_dataframe, 
            init_core_features=core_features, 
            step_size=clustering_step_size,
            hcl_distance_cut = hcl_distance_cut,
            ri_tolerance=ri_tolerance,
            correlation_cut = corr_cutoff,
            max_core_features=max_core_features 
        )
        json_output_file = os.path.join(outdir, project_name_handle+'_list_pseudospectra.json')
        with open(json_output_file, 'w') as O:
            json.dump(
                [port_pseudospectrum_to_json(x) for x in list_pseudospectra], O, indent=2
                )
        print(f"Exported {len(list_pseudospectra)} de novo empCpds to {json_output_file}.")


#
# -----------------------------------------------------------------------------
#

# to test 
def annotate_user_featuretable(infile, parameters, rtime_tolerance=2):
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
    '''
    parameters['outdir'] = ''
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
                            mz_tolerance_ppm=5,
                            rt_tolerance=1,
                            mode='pos',
                            charges=[1, 2, 3],
                            )
    # also exporting JSON
    with open('PreAnnotated_list_khipus.json', 'w', encoding='utf-8') as f:
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


