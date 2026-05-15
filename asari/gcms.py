import json
from typing import NamedTuple
import numpy as np
import pandas as pd
import statsmodels.api as sm

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

import ms_entropy as ME
from jms.search import build_centurion_tree
from .mass_functions import complete_mass_paired_mapping, all_mass_paired_mapping
from .tools.cosine import cosine_similarity
from .tools import match_features 
from .tools.msp_parser import MSP_dict, parse_msp_to_listdict, msp_standarize

# Nonunique features are allowed. I.e. a feature can belong to multiple PseudoSpectrum. 
class PseudoSpectrum(NamedTuple):
    """List of MS features as pseudo spectrum. 
    """
    id: str
    rtime: float
    RI: float       # rentention_index
    rounded_mzs: list
    num_features: int
    members: list
    peaks: np.array
    annotation: str

class GC_lib_entry(NamedTuple):
    """Reformatted library entry per compound, for standardized and fast access. 
    """
    id: str
    inchikey: str
    name: str
    RI: float
    exact_mass: float
    compound_formula: str
    rounded_list: list
    peaks: np.array
    base_peak: tuple
    meta_text: str      # additional text of meta data, records other than peaks


def load_gcms_dbfile(infile):
    if infile.endswith('.msp') or infile.endswith('.MSP'):
        return msp_standarize(parse_msp_to_listdict(infile), MSP_dict)
    elif infile.endswith('.json') or infile.endswith('.JSON'):
        return json.load(open(infile))
    else:
        print("Databse only supports MSP or JSON formats.")
        return None

def reformat_gcms_lib(list_cpd_standards, 
                      peaks_key='peaks', rt_key='RETENTIONTIME', 
                      inchi_key = 'InChIKey', name_key='Name',
                      mw_key = 'ExactMass', formula_key = 'Formula',
                      filter_factor=None):
    '''list_cpd_standards is list of dicts, e.g. result from parse_msp_to_listdict. 
    MSP files are often not well standardized. Attention is required in porting, but this is not a frequenct task.
    '''
    list_lib_entries = []
    counter = 0
    for entry in list_cpd_standards:
        if peaks_key in entry and rt_key in entry and entry[rt_key]:
            _peaks = entry[peaks_key]
            _base_peak = designate_base_peak(_peaks)
            if filter_factor:
                _peaks = filter_peaks_by_low_intensity_factor(_peaks, _base_peak[1], filter_factor)  
            counter += 1
            cpd = GC_lib_entry(
                'cpd' + f"{counter:08}",
                entry.get(inchi_key, ''),
                entry.get(name_key, ''),
                float(entry[rt_key]),
                float(entry[mw_key]),
                entry.get(formula_key, ''),
                [round(x[0]) for x in _peaks],
                np.array(_peaks),
                _base_peak,
                ''          # meta_text
            )
            list_lib_entries.append(cpd)
    return list_lib_entries


def designate_base_peak(peaks):
    ''' 
    peaks : [[71.0491,  1.97  ], [72.0443,  0.32  ], ...]
    '''
    idx = np.argmax([x[1] for x in peaks])
    return tuple(peaks[idx])

def filter_peaks_by_intensity_factor(peaks, base_peak_intensity, filter_factor=100):
    '''Only keep peaks in the intensity bracket.
    peaks : [[71.0491,  1.97  ], [72.0443,  0.32  ], ...]
    '''
    upper, lower = base_peak_intensity*filter_factor, base_peak_intensity/filter_factor
    return [x for x in peaks if upper > x[1] > lower]

def filter_peaks_by_low_intensity_factor(peaks, base_peak_intensity, filter_factor):
    lower = base_peak_intensity/filter_factor
    return [x for x in peaks if x[1] > lower]

def filter_features_by_low_intensity_factor(features, base_peak_intensity, filter_factor):
    '''features in standard JSON notion
    '''
    lower = base_peak_intensity/filter_factor
    return [x for x in features if x['peak_area'] > lower]

def port_pseudospectrum_to_json(PS):
    '''Port an instance of PseudoSpectrum to JSON dict.
    '''
    return {
        'id': PS.id,
        'rtime': PS.rtime,
        'RI': PS.RI,
        # 'rounded_mzs': PS.rounded_mzs,
        'num_features': PS.num_features,
        'members': PS.members,
        'peaks': PS.peaks.tolist(),
        'annotation': PS.annotation
    }

def serialize_annotated_empCpds(list_empCpds):
    '''list_empCpds from curate_batch_lib_search_result
    '''
    for x in list_empCpds:
        x['peaks_in_lib'] = x['peaks_in_lib'].tolist()
        x['peaks_as_features'] = x['peaks_as_features'].tolist()
    return list_empCpds


def ri_penalty_function(abs_ri_delta, min_delta=1, max_delta=100):
    '''Penalty by RI distance, using a 1-side sigmoid function
    '''
    if abs_ri_delta > max_delta:
        return 0
    elif abs_ri_delta < min_delta:
        return 1
    else:
        return 2 / (1 + np.exp(np.e * abs_ri_delta/max_delta))

def filter_peaks_by_penalized_distance(
    selected_features, 
    seed_feature, 
    feature_dataframe,
    min_ri_delta=1,
    max_ri_delta=100, 
    feature_distance_filter=0.5
    ):
    '''
    penalized_distance = Pearsonr * ri_penalty_function
    Slow, not good for large number of features.
    '''
    return [
        p for p in selected_features if 
        feature_dataframe.loc[p['id'], :].corr(feature_dataframe.loc[seed_feature['id'], :]
                                ) * ri_penalty_function(
                                    abs(p['RI'] - seed_feature['RI']), min_ri_delta, max_ri_delta
                                ) >= feature_distance_filter
    ]


#
# Multiple methods of constructing pseudospectra 
#
def find_all_matches_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=5):
    '''
    Return matched peaks in mz_centurion_tree.
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    results = []
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            if abs(peak['mz']-query_mz) < mz_tol:
                results.append(peak)
                
    return results


def get_matched_features_per_cpd(cpd, mz_centurion_tree,
                    mz_tolerance_da, ri_tolerance
                    ):
    '''Get all expt features that match to any peaks in a library compound.
    cpd : GC_lib_entry instance. 
    return list of matched feature IDs. 
    '''
    _matched_feature_ids = []
    for peak in cpd.peaks:      # (mz, intensity)
        q = int(peak[0] * 100)
        for ii in (q-1, q, q+1):
            _matched_feature_ids += [feat['id'] for feat in mz_centurion_tree.get(ii, [])
                               if abs(feat['mz'] - peak[0]) < mz_tolerance_da and 
                               abs(feat['RI'] - cpd.RI) < ri_tolerance
                               ]
    # redundant features are fetched above, thus set to uniques
    return list(set(_matched_feature_ids))

def distill_correlated_features(matched_feature_ids, feature_dataframe, corr_cutoff):
    '''Select core set of correlated features from matched_features.
    '''
    features_intensity = feature_dataframe.loc[matched_feature_ids, :].T
    features_intensity_corr = features_intensity.corr()
    quant_feature = features_intensity_corr.sum(axis=0).idxmax()
    selected = [feat for feat in matched_feature_ids if features_intensity_corr.loc[quant_feature, feat] > corr_cutoff]
    return quant_feature, selected


def get_seeded_pseudospectrum(
    seed_tag, 
    ref_ri,
    seed_feature, 
    list_features, 
    feature_dataframe,
    min_ri_delta=1,
    max_ri_delta=100, 
    low_peak_filter_factor=1000,
    feature_distance_filter=None    # 0.5 is a good default when in use
    ):
    '''
    Use base peaks in lib entries to initiate pseudospectra construction.   
    seed_tag : to track linked information, e.g. annotation, cpd library entry. 
    seed_feature : seed, e.g. the feature matched a lib query or one chosen from a group by top intensity.
    list_features : all features to work with here, in standard JSON notion.
    feature_dataframe : pandas dataframe of sample intensities, indexed by feature ID.
    
    ref_ri : the RI from lib entry
    
    if low_peak_filter_factor: remove features smaller than base peak intensity / filter_factor. 
                    E.g. features have to be above 1% of base peak intensity if filter_factor=100.
    if feature_distance_filter: filter features by penalized distance function. 
    
    '''
    selected_features = [f for f in list_features if abs(f['RI'] - ref_ri) < max_ri_delta]
    
    if low_peak_filter_factor:
        selected_features = filter_features_by_low_intensity_factor( 
            # filter_features_by_low_intensity_factor
            # filter_peaks_by_intensity_factor
            
            selected_features, seed_feature['peak_area'], low_peak_filter_factor
            
            )  
        
    if feature_distance_filter:
        selected_features = filter_peaks_by_penalized_distance(
            selected_features, seed_feature, feature_dataframe, min_ri_delta, max_ri_delta, feature_distance_filter
        )
    num_features = len(selected_features)
    if num_features > 1:
        return PseudoSpectrum(
                    seed_feature['id_number'], 
                    seed_feature['rtime'],
                    ref_ri,
                    [round(x['mz']) for x in selected_features],
                    num_features,
                    [x['id'] for x in selected_features],           # members
                    np.array([(x['mz'], x['peak_area']) for x in selected_features]), 
                    seed_tag
                )
    else:
        return None

def get_spaced_top_features(list_features_sorted, ri_gap=100):
    '''
    Batch mode to get top N features that are not overlapping in RI. 
    The pseudospectra built on the batch then don't have overlapping features (small chance of exception). 
    All core features can be collected before next iteration.
    '''
    selected = []
    for ii in range(1000, 4500, ri_gap): 
        # get first feature that is btw ii and gap
        selected.append(next((f for f in list_features_sorted if ii <= f['RI'] < ii+ri_gap), None))
    return [f for f in selected if f]

def iterative_build_pseudospectra_by_penalizeddistance(
    list_features_sorted,
    feature_dataframe,
    init_core_features = {},
    min_ri_delta=1,
    max_ri_delta=100, 
    low_peak_filter_factor=1000,
    feature_distance_filter=0.5,
    ri_gap=100,
    max_core_features = 20000,
    ):
    '''
    Build pseudospectra from list_features_sorted iteratively until reaching max_core_features 
    (no more than 90% of all features)
    
    returns list_pseudospectra, core_features
    '''
    core_features = init_core_features
    list_pseudospectra = []
    max_core_features = min(max_core_features, 0.9*len(list_features_sorted))
    while len(core_features) <= max_core_features:
        list_features = [f for f in list_features_sorted if f['id'] not in core_features]
        _tmp = []
        for seed_feature in get_spaced_top_features(list_features, ri_gap):
            new = get_seeded_pseudospectrum('',
                                            seed_feature['RI'],
                                            seed_feature, 
                                            list_features, 
                                            feature_dataframe,
                                            min_ri_delta=min_ri_delta,
                                            max_ri_delta=max_ri_delta, 
                                            low_peak_filter_factor=low_peak_filter_factor,
                                            feature_distance_filter=feature_distance_filter
            )
            if new:
                list_pseudospectra.append(new)
                _tmp += new.members
        core_features.union(set(_tmp))
        
    return list_pseudospectra, core_features
        
def have_basepeak_molecularion(base_mz, mole_mz, _peaks, mz_tolerance_da=0.005):
    has_basepeak, has_molecularion = False, False
    if abs(_peaks[:, 0] - base_mz).min() < mz_tolerance_da:
        has_basepeak = True
    if abs(_peaks[:, 0] - mole_mz).min() < mz_tolerance_da:
        has_molecularion = True
    return has_basepeak, has_molecularion


#
# Workflows 
#

def batch_lib_search_score(
    list_cpds, 
    list_features, 
    dict_features,
    feature_dataframe,
    mz_tolerance_da=0.005, 
    ri_tolerance=30,
    cosine_penalty=1,
    corr_cutoff=0.6, 
    ):
    '''Match all features from an experiment to compound library of spectra.
    
    The featuers within mz and RI window are grouped as a pseudospectrum, 
    which excludes features that do not correlated with the core set. 
    
    cosine_score and entropy_score are calculated on each pseudospectrum. 
    '''
    Results = []
    ctree = build_centurion_tree(list_features)
    for cpd in list_cpds:
        _matched_feature_ids = get_matched_features_per_cpd(
            cpd, ctree, mz_tolerance_da, ri_tolerance
            )
        if _matched_feature_ids:
            # Filter by feature corr
            quant_feature, _matched_feature_ids = distill_correlated_features(
                _matched_feature_ids, feature_dataframe, corr_cutoff
            )
            _peaks = np.array([[dict_features[feat]['mz'], dict_features[feat]['peak_area']] 
                            for feat in _matched_feature_ids])
            # calculate_entropy_similarity
            entropy_score = ME.calculate_entropy_similarity(
                _peaks, cpd.peaks, 
                ms2_tolerance_in_da = mz_tolerance_da,
            )
            # cosine_similarity
            cosine_score, num_matched_peaks = cosine_similarity(
                _peaks, cpd.peaks, 
                tolerance=mz_tolerance_da, 
                sqrt_transform=True, penalty=cosine_penalty
            )
            Results.append({
                'lib_entry': cpd,
                'quant_feature': quant_feature,
                'quant_feature_RI': dict_features[quant_feature]['RI'],
                'candidate_feature_ids': _matched_feature_ids,
                'entropy_score': entropy_score,
                'cosine_score': cosine_score,
                'pseudo_spec': _peaks,
                'num_matched_peaks': num_matched_peaks  # should = len(_matched_features)
            })
    print(f"{len(Results)} matched to library compounds.")
    return Results

def curate_batch_lib_search_result(
    matched_results, 
    mz_tolerance_da=0.005, 
    score_cutoff_cosine=0.5, 
    score_cutoff_entropy=0.4,
    ):
    '''Filter and organize results for export. 
    
    determine quant_ion: using the feature of highest sum(correlations to others) in all matched features.

    matched_results : a list of dicts of matched lib_entry and pseudo_spec, from batch_lib_search_score.
    feature_dataframe : dataframe from feature table, indexed on feaature IDs, sample intensity starts in first col.
    
    Data for mirror plot are in (peaks_as_features, peaks_in_lib).

    return list_empCpds, feature_anno_list
    '''
    ii = 0
    list_empCpds = []       # a pseudospectrum matched to a DB entry is considered as an empirical compound
    feature_anno_list = []
    for MM in matched_results:
        if MM['entropy_score'] >= score_cutoff_entropy or MM['cosine_score'] >= score_cutoff_cosine:
            ii += 1
            epd_id = 'empCpd_' + "{:04d}".format(ii)
            matched_features = MM['candidate_feature_ids']
            peaks_as_features = MM['pseudo_spec']
            peaks_in_lib = MM['lib_entry'].peaks
            has_basepeak, has_molecularion = have_basepeak_molecularion(
                MM['lib_entry'].base_peak[0], peaks_in_lib[:,0].max(), peaks_as_features, 
                mz_tolerance_da
            )
            # record empCpd
            list_empCpds.append({
                'id': epd_id,
                'RI': MM['quant_feature_RI'], # use expt data, not lib data MM['lib_entry'].RI,
                'entropy_score': MM['entropy_score'],
                'cosine_score': MM['cosine_score'],
                'lib_entry_id': MM['lib_entry'].id,
                'name': MM['lib_entry'].name, 
                'inchikey': MM['lib_entry'].inchikey, 
                'features': matched_features, 
                'quant_ion': MM['quant_feature'],
                'peaks_as_features': peaks_as_features,
                'peaks_in_lib': peaks_in_lib,
                'has_basepeak': has_basepeak,
                'has_molecularion': has_molecularion
            })
            # record features
            for feature in matched_features:
                feature_anno_list.append({
                    'feature': feature,         # ID only
                    'empCpd': epd_id,
                    'quant_ion': MM['quant_feature'],
                    'entropy_score': MM['entropy_score'],
                    'cosine_score': MM['cosine_score'],
                    'lib_entry_id': MM['lib_entry'].id,
                    'name': MM['lib_entry'].name, 
                    'inchikey': MM['lib_entry'].inchikey, 
                    # 'correlation': round(_corr, 2),
                    # 'is_core': _corr >= corr_cutoff
                })
    
    return list_empCpds, feature_anno_list




def batch_lib_search_by_basepeaks(
    list_lib_entries, 
    list_features, 
    feature_dataframe,
    ms2_tolerance_in_ppm=5, 
    ms2_tolerance_in_da=0.005, 
    ri_tolerance=50,
    cosine_penalty=1, 
    min_ri_delta=1,
    max_ri_delta=100, 
    low_peak_filter_factor=1000,
    feature_distance_filter=None,
    ):
    '''
    Search all PseudoSpectra constructed using base peak matching, and
    return results.
    
    feature_distance_filter : a filter combining feature intensity correlation and RI delta. Slow, not recommended here.
    
    Note: 
    this version of implementation gets poor Cosine scores. 
    Base peak is not necessarily the best entry as it varies by experimental methods/conditions. 
    
    The exact matched features are calcuated in each of the scoring functions. 
    We have yet to find them again for result inspection later.  
    Future improvements can rewrite scoring functions and have them output matched features. 
    
    feature matching needs to be consistent cross steps; base peak or molecular ion ----
    
    '''
    Results = []
    cpd_list = []
    cpdDict = {x.id: x for x in list_lib_entries}
    featureDict = {x['id']: x for x in list_features}
    for entry in list_lib_entries:
        cpd_list.append({
            'id': entry.id, 
            'mz': entry.base_peak[0],   # (mz, intensity)[0]
            'rtime': entry.RI,          # use RI in place of RT
        })
    tmp_list_features = []
    for feat in list_features:
        tmp_list_features.append({
            'id': feat['id'],
            'mz': feat['mz'],
            'rtime': feat['RI']
        })
    matched = match_features.list_match_lcms_features(
        cpd_list, tmp_list_features, 
        mz_ppm=ms2_tolerance_in_ppm, 
        rt_tolerance=ri_tolerance
    )
    for lib_entry_id, v in matched.items():
        for _feature_id in v:
            pseudo_spec = get_seeded_pseudospectrum(lib_entry_id, 
                                                    cpdDict[lib_entry_id].RI,
                                                    featureDict[_feature_id], 
                                                    list_features, 
                                                    feature_dataframe,
                                                    min_ri_delta=min_ri_delta,
                                                    max_ri_delta=max_ri_delta, 
                                                    low_peak_filter_factor=low_peak_filter_factor,
                                                    feature_distance_filter=feature_distance_filter
                                                    )
            if pseudo_spec:
                entry = cpdDict[lib_entry_id]
                candidate_features, spec_peaks = filter_against_libentry(
                    pseudo_spec, entry
                    )            
                # calculate_entropy_similarity
                entropy_score = ME.calculate_entropy_similarity(
                    spec_peaks, entry.peaks, 
                    ms2_tolerance_in_da = ms2_tolerance_in_da,
                )
                # cosine_similarity
                cosine_score, num_matched_peaks = cosine_similarity(
                    spec_peaks, entry.peaks, 
                    tolerance=ms2_tolerance_in_da, 
                    sqrt_transform=True, penalty=cosine_penalty
                )
                Results.append({
                    'lib_entry': entry,
                    'pseudo_spec': pseudo_spec,
                    'candidate_features': candidate_features,   # unit mz matched before precise m/z matching
                    'base_feature_id': _feature_id,
                    'entropy_score': entropy_score,
                    'cosine_score': cosine_score,
                    'num_matched_peaks': num_matched_peaks 
                })
    print(f"{len(Results)} matched results found by base peak search.")
    return Results


def export_feature_annotation_bybasepeaksearch(
    matched_results, 
    feature_dataframe, 
    score_cutoff_cosine=0.5, 
    score_cutoff_entropy=0.4,
    corr_cutoff=0.7, 
    mz_tolerance_ppm=5
    ):
    '''
    matched_results : a dict of matched lib_entry and pseudo_spec, from batch_lib_search_by_basepeaks
    feature_dataframe : dataframe from feature table, indexed on feaature IDs, sample intensity starts in first col.
    Quant feature is using highest intensity in lib.
    Data for mirror plot are in (peaks_as_features, peaks_in_lib)

    return list_empCpds, feature_anno_list
    '''
    ii = 0
    list_empCpds = []       # a pseudospectrum matched to a DB entry is considered as an empirical compound
    feature_anno_list = []
    for MM in matched_results: # feature to lib peak match
        # complete_mass_paired_mapping
        mapped, list1_unmapped, list2_unmapped = all_mass_paired_mapping(
            MM['lib_entry'].peaks[:, 0], MM['pseudo_spec'].peaks[:, 0], 
            std_ppm=mz_tolerance_ppm
        )
        if mapped:
            ii += 1
            epd_id = 'empCpd_' + "{:04d}".format(ii)
            matched_features = [MM['pseudo_spec'].members[pair[1]] for pair in mapped]
            peaks_as_features = MM['pseudo_spec'].peaks[[pair[1] for pair in mapped], :]
            peaks_in_lib = MM['lib_entry'].peaks[[pair[0] for pair in mapped], :]
            quant_feature = MM['base_feature_id']
            if MM['entropy_score'] >= score_cutoff_entropy or MM['cosine_score'] >= score_cutoff_cosine:
                # record empCpd
                list_empCpds.append({
                    'id': epd_id,
                    'RI': MM['pseudo_spec'].RI,
                    'entropy_score': MM['entropy_score'],
                    'cosine_score': MM['cosine_score'],
                    'lib_entry_id': MM['lib_entry'].id,
                    'name': MM['lib_entry'].name, 
                    'inchikey': MM['lib_entry'].inchikey, 
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
                        'entropy_score': MM['entropy_score'],
                        'cosine_score': MM['cosine_score'],
                        'lib_entry_id': MM['lib_entry'].id,
                        'name': MM['lib_entry'].name, 
                        'inchikey': MM['lib_entry'].inchikey, 
                        'correlation': _corr,
                        'is_core': _corr >= corr_cutoff
                    })
    
    return list_empCpds, feature_anno_list

# write out annotated results
def write_tsv_feature_anno(feature_anno_list, dict_features, dict_lib_entries, outfile):
    '''
    The empCpd output is more concise (annotation by empCpd or cluster). 
    This exports a table per annotated feature.
    '''
    header = [
        'feature', 'mz', 'rtime', 'RI', 
        'empCpd', 'quant_ion', 'score_cosine', 'score_entropy', 'name', 'inchikey',
            #'correlation', 'is_core', 
        'lib_RI', 'delta_RI', 'peak_area', 'cSelectivity', 'peak_shape', 'snr', 'detection_counts'
    ]
    s = '\t'.join(header) + '\n'
    for feat in feature_anno_list:
        k = feat['feature']
        FF = dict_features[k]
        score_cosine = round(feat['cosine_score'], 3)
        score_entropy = round(feat['entropy_score'], 3)
        lib_RI = dict_lib_entries[feat['lib_entry_id']].RI
        s += '\t'.join([str(x) for x in [
            k, dict_features[k]['mz'], dict_features[k]['rtime'], round(dict_features[k]['RI'], 3), 
            feat['empCpd'], feat['quant_ion'], score_cosine, score_entropy, feat['name'], feat['inchikey'],
            # round(feat['correlation'], 3), feat['is_core'], 
            round(lib_RI, 3), round(FF['RI']-lib_RI, 3),
            FF['peak_area'], FF['cSelectivity'], FF['goodness_fitting'], FF['snr'], FF['detection_counts'] 
        ]]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)
    
def write_tsv_empCpd_anno(list_empCpds, dict_features, dict_lib_entries, outfile):
    '''Table output for list of empCpds. 
    '''
    header = [
        'empCpd', 'name', 'inchikey', 'formula', 'lib_exact_mass', 'score_cosine', 'score_entropy', 
        'num_matched_peaks', 
        'quant_ion', 'mz_quant_ion', 'rtime_quant_ion(min)', 'RI_quant_ion', 'lib_RI', 'delta_RI', 
        'has_basepeak', 'has_molecularion', 'presence_ratio'
    ]
    s = '\t'.join(header) + '\n'
    for epd in list_empCpds:
        _formula = dict_lib_entries[epd['lib_entry_id']].compound_formula
        lib_exact_mass = dict_lib_entries[epd['lib_entry_id']].exact_mass
        lib_RI = dict_lib_entries[epd['lib_entry_id']].RI
        num_matched_peaks = len(epd['peaks_as_features'])
        quant_ion = epd['quant_ion']
        score_cosine = round(epd['cosine_score'], 3)
        score_entropy = round(epd['entropy_score'], 3)
        s += '\t'.join([str(x) for x in [
            epd['id'], epd['name'], epd['inchikey'], _formula, lib_exact_mass, score_cosine, score_entropy,
            num_matched_peaks, 
            quant_ion, dict_features[quant_ion]['mz'], 
            round(dict_features[quant_ion]['rtime']/60, 2), 
            int(epd['RI']), int(lib_RI), int(epd['RI']-lib_RI), 
            epd['has_basepeak'], epd['has_molecularion'], 
            round(num_matched_peaks/dict_lib_entries[epd['lib_entry_id']].peaks.shape[0], 2)
        ]]) + '\n' 
    with open(outfile, 'w') as O:
        O.write(s)
    
def get_clusters_by_ri_hcl(
    feature_dataframe_sorted, 
    dict_features,
    ri_penalty_function=ri_penalty_function, 
    min_delta=1, 
    max_delta=100, 
    hcl_distance_cut = 1,
    limit_fearture_num=None
    ):
    '''
    Using RI (retention index) coupled hierarchical clustering (hcl) to find feature clusters.
    
    feature_dataframe_sorted : dataframe from feature table, indexed on feaature IDs, 
                    sample intensity starts in first col.
                    Sorted by combined peak area, descending order.
    dict_features : dictionay of features, feature IDs as keys.
    limit_fearture_num : if specified, to use a subset of top-intensity features, 
                    as the low-intensity features are less relevant.
    
    '''
    if limit_fearture_num:
        feature_dataframe_sorted = feature_dataframe_sorted.iloc[:limit_fearture_num, :]
    # ri_penalties
    ri_penalties = []
    featunre_number = feature_dataframe_sorted.shape[0]
    print(f"Working with {featunre_number} features.")
    retention_index = [dict_features[x]['RI'] for x in feature_dataframe_sorted.index] 
    for ii in range(featunre_number):
        for jj in range(ii+1, featunre_number):
            ri_penalties.append(abs(retention_index[ii] - retention_index[jj]))
    ri_penalties = np.array(
        [ri_penalty_function(x,  min_delta, max_delta) for x in ri_penalties]
    )
    # distance function
    YM = pdist(feature_dataframe_sorted.values, 'correlation')
    # filter correlation to accelerate? 
    YM = 1 - ri_penalties * (1 - YM)
    # linkage and cluster
    ZM = linkage(YM, method='ward')
    metClus = fcluster(ZM, hcl_distance_cut, criterion='distance')
    # print(metClus[:3]). # which cluster the feature belongs to

    number_features, number_clusters = len(metClus), len(set(list(metClus)))
    print("number of features: ", number_features)
    print("number of clusters: ", number_clusters)

    # Compile clusters, members being feature IDs
    metClusDict = {}
    for ii in range(number_features):
        if metClus[ii] in metClusDict:
            metClusDict[ metClus[ii] ].append(feature_dataframe_sorted.index[ii])
        else:
            metClusDict[ metClus[ii] ] = [feature_dataframe_sorted.index[ii]]

    return metClus, metClusDict
    
def extend_cluster(fcluster, 
                   featureDict,
                   list_features, 
                   feature_dataframe,
                   ri_tolerance=50, 
                   correlation_cut=0.7
                   ):
    '''
    Extend a cluster of features by correlated features in RI window.
    list_features : features excluding core features already accounted for a cluster
    correlation_cut : filter feature pairs by correlation coefficient
    '''
    cluster_features = [featureDict[f] for f in fcluster]
    cluster_features.sort(key=lambda  x: x['peak_area'], reverse=True)
    seed_feature = cluster_features[0]
    selected_features = [f for f in list_features if abs(f['RI'] - seed_feature['RI']) < ri_tolerance]
    selected_features = [f['id'] for f in selected_features 
                         if feature_dataframe.loc[f['id'], :].corr(feature_dataframe.loc[seed_feature['id'], :]) >= correlation_cut
                         ]
    return (fcluster + selected_features, seed_feature)

def format_fcluster_to_pseudospectrum(selected_features, seed_feature):
    return PseudoSpectrum(
                    seed_feature['id_number'], 
                    seed_feature['rtime'],
                    seed_feature['RI'],
                    [round(x['mz']) for x in selected_features],
                    len(selected_features),
                    [x['id'] for x in selected_features],           # members
                    np.array([(x['mz'], x['peak_area']) for x in selected_features]), 
                    ''
                )

def iterative_build_pseudospectra_by_hcl(
    list_features_sorted,
    feature_dataframe,
    init_core_features = set(),
    step_size = 2000,   # batch size for HCL. Slow if batch too big. 
    hcl_distance_cut = 1, 
    ri_tolerance=50,
    correlation_cut=0.7,
    max_core_features = 20000,
    ):
    '''
    Build seudospectra from list_features_sorted iteratively until reaching max_core_features 
    (no more than 90% of all features)
    
    hcl_distance_cut : parameter to cut linkage tree in HCL. 
        Smaller values give tighter clusters, which are desired since they are extended after this step.
    returns list_pseudospectra, core_features
    '''
    core_features = init_core_features
    list_pseudospectra = []
    max_core_features = min(max_core_features, 0.9*len(list_features_sorted))  # safeguard
    featureDict = {x['id']: x for x in list_features_sorted}
    counter = 0
    while len(core_features) <= max_core_features:
        counter += 1
        remaining_features = [f for f in list_features_sorted if f['id'] not in core_features]
        print(f"Iteration {counter} of building pseudospectra from {len(remaining_features)} remaining features.")
        accounted_features = []
        _, metClusDict = get_clusters_by_ri_hcl(
            feature_dataframe.loc[[f['id'] for f in remaining_features], :], 
            featureDict,
            hcl_distance_cut = hcl_distance_cut,
            limit_fearture_num = step_size
        )
        new_clusters = [
            extend_cluster(clu, featureDict, remaining_features[step_size:],
                           feature_dataframe, ri_tolerance=ri_tolerance, correlation_cut=correlation_cut)
            for clu in metClusDict.values() if len(clu) > 1
        ]
        list_pseudospectra += [format_fcluster_to_pseudospectrum(
            [featureDict[f] for f in x[0]], x[1]) for x in new_clusters]
        for x in new_clusters: 
            accounted_features += x[0]
        core_features = core_features.union(set(accounted_features))
        print(f"Core features collected: {len(core_features)}.")
        
    return list_pseudospectra, core_features


#
# First implementation
#
def filter_against_libentry(query_spectrum, libentry):
    '''return candidate features and peaks, by unit m/z (integer) matching only
    '''
    bool_ = [x in libentry.rounded_list for x in query_spectrum.rounded_mzs]
    candidate_features = [query_spectrum.members[ii] for ii in range(query_spectrum.num_features) if bool_[ii]]
    return candidate_features, query_spectrum.peaks[bool_, :]   # peaks are in 2-D array

def find_entries_in_rtwindow(query, gclib, tol = 30):
    '''
    get entries in gclib database within tol rentention_index of query pseudo spectrum.
    '''
    return [x for x in gclib if abs(query-x.RI) < tol]

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

#
# Functions below may be based on the iterative pseudospectra construction,
# which is bypassed in the base-peak enforced search. 
# -----------------------------------------------------------------------------
# 
# some depend on matched_list format, 
# before implementation of base peak based searches
# 

def group_pseudospectra_from_features(list_features, rtime_window_in_seconds=1, bin_fraction=0.2):
    '''
    Group GC features into pseudo spectra by rtime bins. 
    Starting with features of highest peaks. 
    list_features: list of feature dicts, decending order of representative peak area. 
    throttle: fraction of top features as seeds in grouping pseudo spectra.
    returns list of PseudoSpectrum instances.
    '''
    list_pseudo_spectra, features_counted_for = [], set()
    N = int(len(list_features) * bin_fraction)
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
                np.array([(feature['mz'], feature['peak_area']) for feature in features_in_range]), 
                ''
            )
            features_counted_for.update([f['id_number'] for f in features_in_range])
            list_pseudo_spectra.append(pseudo_)
            
    print("From %d features, %d pseudo spectra were constructed, of which %d have 100 or more peaks." %(
        len(list_features), len(list_pseudo_spectra), len([x for x in list_pseudo_spectra if x.num_features>=100])
    ))
    return list_pseudo_spectra

def reverse_spec_searches(list_pseudo_spectra,
                          list_lib_entries,
                          ri_window=100,
                          mz_tolerance=0.005,
                          cosine_penalty=1, 
                          score_cutoff=0.5, 
                          score_cutoff_entropy=0.4,
                          ):
    '''Search list_pseudo_spectra by standards in list_lib_entries, without looking for base peaks.
    Using both everse cosine and ms_entrophy. 
    cosine_penalty=1 means traditional reverse cosine.
    
    score_cutoff : cutoff for reverse cosine
    score_cutoff_entropy : cutoff for MSentropy
    
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
            if _score > score_cutoff_entropy:
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

def iterative_reverse_annotation(
    list_features, 
    list_lib_entries, 
    feature_dataframe, 
    binning_rtime_window_in_seconds=1,
    search_ri_window=50, 
    search_mz_tolerance=0.005, 
    cosine_penalty=1, 
    score_cutoff=0.5,          # for reverse cosine
    score_cutoff_entropy=0.4,  # for MSentropy
    corr_cutoff=0.7, 
    export_mz_tolerance_ppm=5,
    bin_fraction=0.2,
    iterations=3
    ):
    '''
    list_features : list of features with RI in NamedTuples. 
    list_lib_entries : lib/database entries in NamedTuples. 
    feature_dataframe : dataframe of feature/sample intensities, to calculate feature correlations. 
    iterations : the pseudo spectra are constructed by seeding with top features. 
                Default 25% each iteration, increased to all in last iteration. 
    
    Returns dict of cosine and entropy results, and core features.
    '''
    core_features = set()
    list_empCpds_cosine, feature_anno_list_cosine = [], []
    list_empCpds_entropy, feature_anno_list_entropy = [], []
    for step in range(1, iterations+1):
        print("Iteration step ", step)
        if step == iterations:
            bin_fraction = 1
        list_pseudo_spectra = group_pseudospectra_from_features(
            [feat for feat in list_features if feat['id'] not in core_features], 
            rtime_window_in_seconds=binning_rtime_window_in_seconds, bin_fraction=bin_fraction)

        matched_entropy, matched_cosine = reverse_spec_searches(
            list_pseudo_spectra, list_lib_entries, ri_window=search_ri_window, 
            mz_tolerance=search_mz_tolerance, 
            cosine_penalty=cosine_penalty, score_cutoff=score_cutoff, score_cutoff_entropy=score_cutoff_entropy
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
