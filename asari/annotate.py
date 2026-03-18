'''
Functions for subcommand `annotate`
To test and solidify.

GC-HRMS annotation functions are added here.
This does not change default --anno for LC-MS data processing; but we will unify in next version.

'''
import json
from typing import NamedTuple
import numpy as np
import pandas as pd
import statsmodels.api as sm
import ms_entropy as ME

from khipu.extended import peaklist_to_khipu_list 
#, export_empCpd_khipu_list

# When use relative import, sys.path.append errors out
from .default_parameters import (extended_adducts, adduct_search_patterns_pos, 
                         adduct_search_patterns_neg, 
                         isotope_search_patterns)
from .utils import NpEncoder
from .mass_functions import complete_mass_paired_mapping

from .tools.file_io import read_features_from_asari_table
from .tools.cosine import cosine_similarity
from .tools import match_features 


class PseudoSpectrum(NamedTuple):
    """MS features in GC-MS as pseudo spectrum"""
    id: str
    rtime: float
    RI: float       # rentention_index
    rounded_mzs: list
    num_features: int
    members: list
    peaks: np.array

class GC_lib_entry(NamedTuple):
    inchikey: str
    name: str
    RI: float
    exact_mass: float
    compound_formula: str
    rounded_list: list
    peaks: np.array

def filter_against_libentry(query_spectrum, libentry):
    '''return candidate features and peaks'''
    bool_ = [x in libentry.rounded_list for x in query_spectrum.rounded_mzs]
    candidate_features = [query_spectrum.members[ii] for ii in range(query_spectrum.num_features) if bool_[ii]]
    return candidate_features, query_spectrum.peaks[bool_, :]   # peaks are in 2-D array

def find_entries_in_rtwindow(query, gclib, tol = 30):
    '''
    get entries in gclib database within tol rentention_index of query pseudo spectrum.
    '''
    return [x for x in gclib if abs(query-x.RI) < tol]

def reformat_gcms_lib(list_cpd_standards):
    '''list_cpd_standards is list of dicts, e.g. result from parse_msp_to_listdict
    '''
    list_lib_entries = []
    for entry in list_cpd_standards:
        if 'peaks' in entry and 'RETENTIONTIME' in entry and entry['RETENTIONTIME']:
            cpd = GC_lib_entry(
                entry.get('InChIKey', ''),
                entry.get('Name', ''),
                float(entry['RETENTIONTIME']),
                float(entry['ExactMass']),
                entry.get('Formula', ''),
                [round(x[0]) for x in entry['peaks']],
                np.array(entry['peaks'])
            )
            list_lib_entries.append(cpd)
    return list_lib_entries

def read_fit_KovatsIndex_rtime(KovatsIndex_file, sep='\t', frac=0.3):
    kovats = pd.read_csv(KovatsIndex_file, sep=sep)
    # loess regression
    model = sm.nonparametric.lowess(kovats.values[:, 0], kovats.values[:, 1], frac=frac)
    return model

def append_kovats_index(list_features, ri_model):
    '''Add RI to each featulre in list_features. 
    ri_model : np.array of regression result.
    '''
    for feature in list_features:
        feature['RI'] = float(np.interp(feature['rtime']/60, ri_model[:, 0], ri_model[:, 1])
                    )   # asari rtime is in seconds, convert to minutes
    return list_features

def group_pseudospectra_from_features(list_features, rtime_window_in_seconds=1, throttle=0.25):
    '''
    Group GC features into pseudo spectra by rtime bins. 
    Starting with features of highest peaks. 
    list_features: list of feature dicts, decending order of representative peak area. 
    throttle: fraction of top features as seeds in grouping pseudo spectra.
    returns list of PseudoSpectrum instances.
    '''
    list_pseudo_spectra, features_counted_for = [], set()
    N = int(len(list_features) * throttle)
    for feature in list_features[:N]:
        if feature['id_number'] not in features_counted_for:
            features_in_range = [f for f in list_features if abs(f['rtime'] - feature['rtime']) < rtime_window_in_seconds]
            pseudo_ = PseudoSpectrum(
                feature['id_number'], 
                feature['rtime'],
                feature['RI'],
                [round(feature['mz']) for feature in features_in_range],
                len(features_in_range),
                [f['id_number'] for f in features_in_range],
                np.array([(feature['mz'], feature['peak_area']) for feature in features_in_range])
            )
            features_counted_for.update([f['id_number'] for f in features_in_range])
            list_pseudo_spectra.append(pseudo_)
            
    print("From %d features, %d pseudo spectra were constructed, of which %d have 3 or more peaks." %(
        len(list_features), len(list_pseudo_spectra), len([x for x in list_pseudo_spectra if x.num_features>=3])
    ))
    return list_pseudo_spectra

def reverse_spec_searches(list_pseudo_spectra,
                          list_lib_entries,
                          ri_window=100,
                          mz_tolerance=0.005,
                          cosine_penalty=1, 
                          score_cutoff=0.5
                          ):
    '''Search list_pseudo_spectra by standards in list_lib_entries.
    Using both everse cosine and ms_entrophy. 
    cosine_penalty=1 means traditional reverse cosine.
    Returns [[score, entry, candidate_features, spec], 
            [_score, entry, candidate_features, spec, _n_matches]]
    
    '''
    matched_cosine, matched_entropy = [], []
    for spec in list_pseudo_spectra:
        candidate_lib_entries = find_entries_in_rtwindow( spec.RI, list_lib_entries, tol=ri_window
            )
        for entry in candidate_lib_entries:
            candidate_features, spec_peaks = filter_against_libentry(spec, entry)
            # calculate_entropy_similarity
            _score = ME.calculate_entropy_similarity(
                spec_peaks, entry.peaks, 
                ms2_tolerance_in_da = mz_tolerance,
            )
            if _score > score_cutoff:
                matched_entropy.append((_score, entry, candidate_features, spec))
            # cosine_similarity
            _score, _n_matches = cosine_similarity(
                spec_peaks, entry.peaks, 
                tolerance=mz_tolerance, sqrt_transform=True, penalty=cosine_penalty
            )
            if _score > score_cutoff:
                matched_cosine.append((_score, entry, candidate_features, spec, _n_matches))
                
    print("Found %d by entropy serach and %d by cosine search. " %(len(matched_entropy), len(matched_cosine)))
    return matched_entropy, matched_cosine
    
    
def export_feature_annotations(matched_list, mz_tolerance_ppm=5):
    '''
    Simple method to export results from reverse_spec_searches
    matched_list: [score, entry, candidate_features, spec, *]
    return list of [(feature ID, score, cpd name, inchikey, RI), ...]
    '''
    feature_anno_list = []
    for MM in matched_list:
        mapped, list1_unmapped, list2_unmapped = complete_mass_paired_mapping(
            MM[1].peaks[:, 0], MM[3].peaks[:, 0], std_ppm=mz_tolerance_ppm
        )
        for pair in mapped:
            feature_anno_list.append(
                # feature ID, score, cpd name, inchikey, RI
                (MM[3].members[pair[1]], MM[0], MM[1].name, MM[1].inchikey, MM[3].RI)
            )
    return feature_anno_list

def export_feature_annotation_details(matched_list, feature_dataframe, corr_cutoff=0.7, mz_tolerance_ppm=5):
    '''
    matched_list: [score, entry, candidate_features, spec, *]
    
    Quant feature is chosen by highest intensity rank combining lib and pseudospectrum
    my_list = matched_entropy[55][1].peaks[:, 1]
    my_list2 = matched_entropy[55][3].peaks[:, 1]
    top_feature_idx = np.argmax(np.argsort(my_list) + np.argsort(my_list2))
    
    Correlation of members to quant feature are calculated; 
    members with corr >= corr_cutoff are designated as core members. 
    
    Data for mirror plot are in (peaks_as_features, peaks_in_lib)

    return list_empCpds, feature_anno_list
    '''
    ii = 0
    list_empCpds = []       # a pseudospectrum matched to a DB entry is considered as an empirical compound
    feature_anno_list = []
    for MM in matched_list: # feature to lib peak match
        mapped, list1_unmapped, list2_unmapped = complete_mass_paired_mapping(
            MM[1].peaks[:, 0], MM[3].peaks[:, 0], std_ppm=mz_tolerance_ppm
        )
        if mapped:
            ii += 1
            epd_id = 'empCpd_' + "{:04d}".format(ii)
            matched_features = [MM[3].members[pair[1]] for pair in mapped]
            peaks_as_features = MM[3].peaks[[pair[1] for pair in mapped], :]
            peaks_in_lib = MM[1].peaks[[pair[0] for pair in mapped], :]
            top_feature_idx = np.argmax(np.argsort(peaks_as_features[:, 1]) + np.argsort(peaks_in_lib[:, 1]))
            quant_feature = matched_features[top_feature_idx]
            # record empCpd
            list_empCpds.append({
                'id': epd_id,
                'RI': MM[3].RI,
                'score': MM[0], 
                'name': MM[1].name, 
                'inchikey': MM[1].inchikey, 
                'features': matched_features, 
                'quant_ion': quant_feature,
                'peaks_as_features': peaks_as_features,
                'peaks_in_lib': peaks_in_lib
            })
            # record features
            for feature in matched_features:
                _corr = feature_dataframe.loc[quant_feature, :].corr(feature_dataframe.loc[feature, :]) # calculate corr
                feature_anno_list.append({
                    'feature': feature,         # ID only
                    'empCpd': epd_id,
                    'quant_ion': quant_feature,
                    # 'RI_delta': MM[3].RI - MM[1].RI,  # using empCpd RI as proxy, minor error introduced. Alternative calculation on export.
                    'score': MM[0], 
                    'name': MM[1].name, 
                    'inchikey': MM[1].inchikey, 
                    'correlation': _corr,
                    'is_core': _corr >= corr_cutoff
                })
    
    return list_empCpds, feature_anno_list

def iterative_reverse_annotation(list_features, 
                                 list_lib_entries, 
                                 feature_dataframe, 
                                 binning_rtime_window_in_seconds=1,
                                 search_ri_window=50, 
                                 search_mz_tolerance=0.005, 
                                 cosine_penalty=1, 
                                 score_cutoff=0.5, 
                                 corr_cutoff=0.7, 
                                 export_mz_tolerance_ppm=5,
                                 iterations=3
                                 ):
    
    '''
    list_features : list of features with RI 
    entropy search could use offset on score_cutoff
    
    '''
    core_features = set()
    list_empCpds_cosine, feature_anno_list_cosine = [], []
    list_empCpds_entropy, feature_anno_list_entropy = [], []
    for step in range(1, iterations+1):
        print("Iteration step ", step)
        _fraction_ = 0.25
        if step == iterations:
            _fraction_ = 1
        list_pseudo_spectra = group_pseudospectra_from_features(
            [feat for feat in list_features if feat['id'] not in core_features], 
            rtime_window_in_seconds=binning_rtime_window_in_seconds, throttle=_fraction_)

        matched_entropy, matched_cosine = reverse_spec_searches(
            list_pseudo_spectra, list_lib_entries, ri_window=search_ri_window, 
            mz_tolerance=search_mz_tolerance, 
            cosine_penalty=cosine_penalty, score_cutoff=score_cutoff
        )
        
        _list_empCpds_, _feature_anno_list_ = export_feature_annotation_details(
            matched_cosine, feature_dataframe, 
            corr_cutoff=corr_cutoff, mz_tolerance_ppm=export_mz_tolerance_ppm
        )
        list_empCpds_cosine += _list_empCpds_
        feature_anno_list_cosine += _feature_anno_list_
        _list_empCpds_, _feature_anno_list_ = export_feature_annotation_details(
            matched_entropy, feature_dataframe, 
            corr_cutoff=corr_cutoff, mz_tolerance_ppm=export_mz_tolerance_ppm
        )
        list_empCpds_entropy += _list_empCpds_
        feature_anno_list_entropy += _feature_anno_list_
        core_features = core_features.union(set(
            [f['feature'] for f in feature_anno_list_cosine + feature_anno_list_entropy if f['is_core']]))
        
    return {
        'list_empCpds_cosine': list_empCpds_cosine,
        'feature_anno_list_cosine': feature_anno_list_cosine,
        'list_empCpds_entropy': list_empCpds_entropy,
        'feature_anno_list_entropy': feature_anno_list_entropy,
        'core_features': core_features
    }

def is_same_match(epd1, epd2):
    '''Determine two empCpds as same if same inchikey, quant_ion and features
    '''
    if epd1['inchikey'] == epd2['inchikey'] and \
        epd1['quant_ion'] == epd2['quant_ion'] and epd1['features'] == epd2['features']:
        return True
    else:
        return False

def is_same_matched_feature(f1, f2):
    if f1['inchikey'] == f2['inchikey'] and f1['quant_ion'] == f2['quant_ion']:
        return True
    else:
        return False

def merge_epd_with_list(epd1, list2_epds):
    '''If epd1 already in list2, merge to matched epd in list2;
    else append
    '''
    _matched = False
    for epd in list2_epds:
        if is_same_match(epd, epd1):
            epd.update(epd1)
            _matched = True
    if not _matched:
        list2_epds.append(epd1)
    return list2_epds

def merge_feature_with_list(feat, list2_features):
    _matched = False
    for f in list2_features:
        if is_same_matched_feature(feat, f):
            f.update(feat)
            _matched = True
    if not _matched:
        list2_features.append(feat)
    return list2_features

def format_arrays_to_lists(epd):
    '''Reformat the np.arrays in epd to lists. 
    peaks_as_features, peaks_in_lib
    '''
    epd['peaks_as_features'] = epd['peaks_as_features'].tolist()
    epd['peaks_in_lib'] = epd['peaks_in_lib'].tolist()
    return epd

def cleanup_anno_empcpds_features(resultDict):
    '''
    Combine repeated matches from two search algoirthms. 
    Include both scores in list_empCpds and feature_anno_list. 
    
    resultDict = {
        'list_empCpds_cosine': list_empCpds_cosine,
        'feature_anno_list_cosine': feature_anno_list_cosine,
        'list_empCpds_entropy': list_empCpds_entropy,
        'feature_anno_list_entropy': feature_anno_list_entropy,
        'core_features': core_features
    }
    
    Returns dict_empCpds, dict_feature_anno
    '''
    dict_empCpds = {}
    dict_feature_anno = {}
    # do empCpds, using inchikey as key
    for epd in resultDict['list_empCpds_cosine']:
        epd = format_arrays_to_lists(epd)
        inchikey = epd['inchikey']
        quant_ion = epd['quant_ion']
        epd['score_cosine'] = epd['score']
        _ = epd.pop('score')
        if inchikey not in dict_empCpds:
            dict_empCpds[inchikey] = [epd]
        else:
            if True not in [is_same_match(epd, e) for e in dict_empCpds[inchikey]]:
                dict_empCpds[inchikey].append(epd)
    for epd in resultDict['list_empCpds_entropy']:
        epd = format_arrays_to_lists(epd)
        inchikey = epd['inchikey']
        quant_ion = epd['quant_ion']
        epd['score_entropy'] = epd['score']
        _ = epd.pop('score')
        if inchikey not in dict_empCpds:
            dict_empCpds[inchikey] = [epd]
        else:
            dict_empCpds[inchikey] = merge_epd_with_list(epd, dict_empCpds[inchikey])
            
    # do features
    for feat in resultDict['feature_anno_list_cosine']:
        feat['score_cosine'] = feat['score']
        _ = feat.pop('score')
        if feat['feature'] not in dict_feature_anno:
            dict_feature_anno[feat['feature']] = [feat]
        else:
            if True not in [is_same_matched_feature(feat, f) for f in dict_feature_anno[feat['feature']]]:
                dict_feature_anno[feat['feature']].append(feat)
    for feat in resultDict['feature_anno_list_entropy']:
        feat['score_entropy'] = feat['score']
        _ = feat.pop('score')
        if feat['feature'] not in dict_feature_anno:
            dict_feature_anno[feat['feature']] = [feat]
        else:
            dict_feature_anno[feat['feature']] = merge_feature_with_list(feat, 
                                                    dict_feature_anno[feat['feature']])
            
    return dict_empCpds, dict_feature_anno
    
def write_tsv_feature_anno(dict_feature_anno, dict_features, dict_lib_entries, outfile):
    '''
    The empCpd output is more concise (annotation by empCpd or cluster). 
    This exports a table per annotated feature.

    '''
    header = [
        'feature', 'mz', 'rtime', 'RI', 
        'empCpd', 'quant_ion', 'score_cosine', 'score_entropy', 'name', 'inchikey',
            'correlation', 'is_core', 'lib_RI', 'delta_RI', 
            'peak_area', 'cSelectivity', 'peak_shape', 'snr', 'detection_counts'
    ]
    s = '\t'.join(header) + '\n'
    for k,LL in dict_feature_anno.items():
        FF = dict_features[k]
        for feat in LL:
            if 'score_cosine' in feat:
                score_cosine = round(feat['score_cosine'], 3)
            else:
                score_cosine = ''
            if 'score_entropy' in feat:
                score_entropy = round(feat['score_entropy'], 3)
            else:
                score_entropy = ''
            lib_RI = dict_lib_entries[feat['inchikey']].RI
            s += '\t'.join([str(x) for x in [
                k, dict_features[k]['mz'], dict_features[k]['rtime'], round(dict_features[k]['RI'], 3), 
                feat['empCpd'], feat['quant_ion'], score_cosine, score_entropy, feat['name'], feat['inchikey'],
                round(feat['correlation'], 3), feat['is_core'], round(lib_RI, 3), round(FF['RI']-lib_RI, 3),
                FF['peak_area'], FF['cSelectivity'], FF['goodness_fitting'], FF['snr'], FF['detection_counts'] 
            ]]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)
    
def write_tsv_empCpd_anno(dict_empCpds, dict_features, dict_lib_entries, outfile):
    '''
    Table output for list of empCpds
    '''
    header = [
        'empCpd', 'name', 'inchikey', 'formula', 'lib_exact_mass', 'score_cosine', 'score_entropy', 
        'num_matched_peaks', 
        'quant_ion', 'mz_quant_ion', 'rtime_quant_ion(min)', 'RI_quant_ion', 'RI_empCpd', 'lib_RI', 'delta_RI', 
    ]
    s = '\t'.join(header) + '\n'
    for inchikey,LL in dict_empCpds.items():
        formula = dict_lib_entries[inchikey].compound_formula
        lib_exact_mass = dict_lib_entries[inchikey].exact_mass
        lib_RI = dict_lib_entries[inchikey].RI
        for epd in LL:
            num_matched_peaks = len(epd['features'])
            quant_ion = epd['quant_ion']
            if 'score_cosine' in epd:
                score_cosine = round(epd['score_cosine'], 3)
            else:
                score_cosine = ''
            if 'score_entropy' in epd:
                score_entropy = round(epd['score_entropy'], 3)
            else:
                score_entropy = ''
            s += '\t'.join([str(x) for x in [
                epd['id'], epd['name'], inchikey, formula, lib_exact_mass, score_cosine, score_entropy,
                num_matched_peaks, 
                quant_ion, dict_features[quant_ion]['mz'], 
                round(dict_features[quant_ion]['rtime']/60, 2), int(dict_features[quant_ion]['RI']),
                int(epd['RI']), int(lib_RI), int(epd['RI']-lib_RI)
            ]]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)
    

def append_anno_json(list_features, search_matched_features):
    '''
    search_matched_features, result from export_feature_annotations, 
    like [('F9296', 0.7512193986080801, 'Indole', 'SIKJAQJRHWYJAI-UHFFFAOYSA-N', RI), ...]
    '''
    _d = {}
    for entry in search_matched_features:
        # round score, RI
        _e = (round(entry[1], 3), entry[2], entry[3])
        if entry[0] in _d:
            _d[entry[0]].append(_e)
        else:
            _d[entry[0]] = [_e]
            
    for f in list_features:
        if f['id'] in _d:
            f['annotation'] = sorted(_d[f['id']], reverse=True)
        
    return list_features
    
def export_jsonanno_tsv(list_features, outfile):
    header = ['id_number', 'mz', 'RI',
              'top_match', 'all_matches', #('score', 'name', 'inchikey')
              'rtime', 'left_base', 'right_base', 
              'parent_masstrack_id', 'peak_area', 'cSelectivity', 'goodness_fitting', 
              'snr', 'detection_counts'
              ]
    s = '\t'.join(header) + '\n'
    for f in list_features:
        f['top_match'] = ''
        f['all_matches'] = []
        if 'annotation' in f:
            f['top_match'] = f['annotation'][0][1]
            f['all_matches'] = f['annotation']
        s += '\t'.join([str(x) for x in [
            f[ii] for ii in header
        ]]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)


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


