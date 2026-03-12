'''
Functions for subcommand `annotate`



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
from default_parameters import (extended_adducts, adduct_search_patterns_pos, 
                         adduct_search_patterns_neg, 
                         isotope_search_patterns)
from utils import NpEncoder
from mass_functions import complete_mass_paired_mapping

from tools.file_io import read_features_from_asari_table
from tools.cosine import cosine_similarity
from tools import match_features 

# from jms.io import read_table_to_peaks

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
    '''get entries with tol rentention_index'''
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

def group_pseudospectra_from_features(list_features, ri_model, rtime_window_in_seconds=2):
    '''
    Group GC features into pseudo spectra by rtime bins. 
    Starting with features of highest peaks. 
    list_features: list of feature dicts, decending order of representative peak area. 
    returns list of PseudoSpectrum instances.
    '''
    list_pseudo_spectra, features_counted_for = [], set()
    for feature in list_features:
        if feature['id_number'] not in features_counted_for:
            features_in_range = [f for f in list_features if abs(f['rtime'] - feature['rtime']) < rtime_window_in_seconds]
            pseudo_ = PseudoSpectrum(
                feature['id_number'], 
                feature['rtime'],
                float(np.interp(feature['rtime']/60, ri_model[:, 0], ri_model[:, 1])
                    ),   # asari rtime is in seconds, convert to minutes
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
        candidates = find_entries_in_rtwindow( spec.RI, list_lib_entries, tol=ri_window
            )
        for entry in candidates:
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
    [score, entry, candidate_features, spec, *]
    
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

def append_anno_json(list_features, search_matched_features, ri_model):
    '''
    Search result as [('F9296', 0.7512193986080801, 'Indole', 'SIKJAQJRHWYJAI-UHFFFAOYSA-N', RI), ...]
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
        # add RI
        f['RI'] = round(float(np.interp(f['rtime']/60, ri_model[:, 0], ri_model[:, 1])), 2)
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


