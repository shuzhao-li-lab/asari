"""
GC-MS feature processing is basically the same as LC-MS processing.
Retention index calculation can be either in the `process` step or in `deconvolute` or `annotate`,
based on a Kovats RI table to derive a regression model.

Default GC-MS annotation is using a reverse search against a high-resoluiton database, which requires no deconvolution.

Forward approaches are included too, first `deconvolute` then `annotate`.

Example use:

asari process -i mydir/projectx_dir -o test/ --workflow GC --kovats KovatsIndex_alkanestandards_cuimc.tsv

asari annotate -i project_311234834/export/full_Feature_table.tsv -o project_311234834 \
    --workflow GC  --search_algorithm reverse
    --db GCHRMS_Database_251217.msp

asari deconvolute -i project_311234834/export/full_Feature_table.tsv -o project_311234834 --deconvolution_algorithm hcl --workflow GC

asari annotate -i project_311234834/v19h_denovo_pseudospectra.json -o project_311234834 \
    --workflow GC  --search_algorithm forward --highres False
    --db lowres_database.msp
"""

import json
import re
from typing import NamedTuple
import numpy as np
import pandas as pd
import statsmodels.api as sm


from scipy import sparse

import ms_entropy as ME
from jms.search import build_centurion_tree
from .tools.cosine import cosine_similarity
from .tools.msp_parser import MSP_dict, parse_msp_to_listdict, msp_standarize, parse_ri_value

# from .mass_functions import complete_mass_paired_mapping, all_mass_paired_mapping
# from .tools import match_features

# ---------------------------------------------------------------------------
# Data structures and I/O
# ---------------------------------------------------------------------------

# Nonunique features are allowed. I.e. a feature can belong to multiple PseudoSpectrum.
class PseudoSpectrum(NamedTuple):
    """List of MS features as pseudo spectrum.
    """
    id: str
    rtime: float
    RI: float       # rentention_index; None if no ri_model
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
    """Load a GC-MS library or spectra from MSP or JSON file."""
    if infile.endswith('.msp') or infile.endswith('.MSP'):
        return msp_standarize(parse_msp_to_listdict(infile), MSP_dict)
    elif infile.endswith('.json') or infile.endswith('.JSON'):
        return json.load(open(infile))
    else:
        print("Databse only supports MSP or JSON formats.")
        return None

def reformat_gcms_lib(list_cpd_standards,
                      peaks_key='peaks', ri_keys=('RI', 'RETENTIONTIME'),
                      inchi_key = 'InChIKey', name_key='Name',
                      mw_key = 'ExactMass', formula_key = 'Formula',
                      filter_factor=None,
                      mz_range=(50, 750),
                      min_peak_num=2,
                      ):
    """list_cpd_standards is list of dicts, e.g. result from parse_msp_to_listdict.
    MSP files are often not well standardized. Attention is required in porting, but this is not a frequenct task.
    ri_keys : fields to read retention index from, first non-empty one wins.
        'RI' is standardized by msp_parser from variants (Retention_index, Kovats, ...);
        'RETENTIONTIME' is how some GC libraries (e.g. gcms_corelib) store RI.
    Entries without RI are skipped.
    """
    list_lib_entries = []
    counter = 0
    for entry in list_cpd_standards:
        ri = next((parse_ri_value(entry[k]) for k in ri_keys
                   if entry.get(k) not in (None, '') and parse_ri_value(entry[k]) is not None), None)
        if peaks_key in entry and ri is not None:
            _peaks = entry[peaks_key]
            if mz_range:    # must do before designating base peak
                _peaks = [x for x in _peaks if mz_range[0] < x[0] < mz_range[1]]
            if _peaks and len(_peaks) >= min_peak_num:
                _base_peak = designate_base_peak(_peaks)
                if filter_factor:
                    _peaks = filter_peaks_by_low_intensity_factor(_peaks, _base_peak[1], filter_factor)
                counter += 1
                cpd = GC_lib_entry(
                    'cpd' + f"{counter:08}",
                    entry.get(inchi_key, ''),
                    entry.get(name_key, ''),
                    ri,
                    float(entry[mw_key]) if entry.get(mw_key) not in (None, '') else None,
                    entry.get(formula_key, ''),
                    [round(x[0]) for x in _peaks],
                    np.array(_peaks),
                    _base_peak,
                    ''          # meta_text
                )
                list_lib_entries.append(cpd)
    return list_lib_entries

def json_pseudospectra_to_msp(records, outfile='test.msp'):
    """Convert a JSON list of pseudospectra (as produced by asari GC workflow) to MSP format.

    Each entry in the JSON list becomes one MSP record. The m/z values are written with
    4 decimal places and intensities as integers, matching the convention used for
    high-resolution GC-MS data.
    records = json.load(open(infile))

    Parameters
    ----------
    infile : str
        Path to the input JSON file (list of pseudospectrum dicts with keys
        ``id``, ``rtime``, ``RI``, ``peaks``, ``annotation``, etc.).
    outfile : str, optional
        Path for the output MSP file.  Defaults to ``infile`` with ``.json``
        replaced by ``.msp``.

    Returns
    -------
    str
        Path of the written MSP file.
    """
    lines = []
    for entry in records:
        peaks = entry.get('peaks', [])
        if not peaks:
            peaks = entry.get('peaks_as_features', [])

        lines.append(f"ID: {entry.get('id', '')}")
        lines.append(f"Name: {entry.get('name', '')}")
        ri = entry.get('RI')
        lines.append(f"RETENTIONTIME: {'' if ri is None else ri}")
        lines.append(f"InChIKey: {entry.get('inchikey', '')}")
        features = entry.get('features', []) or entry.get('members', [])
        if features:
            lines.append(f"Features: {','.join(features)}")
        annotation = entry.get('annotation', '')
        if annotation:
            lines.append(f"Comment: {annotation}")

        lines.append(f"Num Peaks: {len(peaks)}")
        for mz, intensity in peaks:
            lines.append(f"{mz:.4f} {intensity:.4f}")

        lines.append("")   # blank line between records

    with open(outfile, 'w') as fh:
        fh.write('\n'.join(lines))

def designate_base_peak(peaks):
    """
    peaks : [[71.0491,  1.97  ], [72.0443,  0.32  ], ...]
    """
    idx = np.argmax([x[1] for x in peaks])
    return tuple(peaks[idx])

def filter_peaks_by_intensity_factor(peaks, base_peak_intensity, filter_factor=100):
    """Only keep peaks in the intensity bracket.
    peaks : [[71.0491,  1.97  ], [72.0443,  0.32  ], ...]
    """
    upper, lower = base_peak_intensity*filter_factor, base_peak_intensity/filter_factor
    return [x for x in peaks if upper > x[1] > lower]

def filter_peaks_by_low_intensity_factor(peaks, base_peak_intensity, filter_factor):
    """Keep peaks with intensity above base_peak_intensity/filter_factor."""
    lower = base_peak_intensity/filter_factor
    return [x for x in peaks if x[1] > lower]

def filter_features_by_low_intensity_factor(features, base_peak_intensity, filter_factor):
    """features in standard JSON notion
    """
    lower = base_peak_intensity/filter_factor
    return [x for x in features if x['peak_area'] > lower]

def filter_against_libentry(query_spectrum, libentry):
    """return candidate features and peaks, by unit m/z (integer) matching only
    """
    bool_ = [x in libentry.rounded_list for x in query_spectrum.rounded_mzs]
    candidate_features = [query_spectrum.members[ii] for ii in range(query_spectrum.num_features) if bool_[ii]]
    return candidate_features, query_spectrum.peaks[bool_, :]   # peaks are in 2-D array

def find_entries_in_rtwindow(query, gclib, tol = 30):
    """
    get entries in gclib database within tol rentention_index of query pseudo spectrum.
    """
    return [x for x in gclib if abs(query-x.RI) < tol]

def read_fit_KovatsIndex_rtime(KovatsIndex_file, sep='\t', frac=0.3):
    """Fit a lowess model of Kovats index on retention time (min) from a Kovats standards file."""
    kovats = pd.read_csv(KovatsIndex_file, sep=sep)
    # loess regression
    model = sm.nonparametric.lowess(kovats.values[:, 0], kovats.values[:, 1], frac=frac)
    return model

def append_kovats_index(list_features, ri_model):
    """Add RI to each featulre in list_features.
    ri_model : np.array of regression result.
    """
    for feature in list_features:
        feature['RI'] = float(np.interp(feature['rtime']/60, ri_model[:, 0], ri_model[:, 1])
                    )   # asari rtime is in seconds, convert to minutes
    return list_features

def port_pseudospectrum_to_json(PS, normalize_intensity=False):
    """Port an instance of PseudoSpectrum to JSON dict.
    """
    peaks = PS.peaks.copy()
    if normalize_intensity:
        peaks[:, 1] /= peaks[:, 1].max()
    return {
        'id': PS.id,
        'rtime': PS.rtime,
        'RI': PS.RI,
        # 'rounded_mzs': PS.rounded_mzs,
        'num_features': PS.num_features,
        'members': PS.members,
        'peaks': peaks.tolist(),
        'annotation': PS.annotation
    }

def serialize_annotated_empCpds(list_empCpds):
    """list_empCpds from curate_batch_lib_reverse_search_result or curate_batch_lib_forward_search_result
    """
    for x in list_empCpds:
        x['peaks_in_lib'] = x['peaks_in_lib'].tolist()
        x['peaks_as_features'] = x['peaks_as_features'].tolist()
    return list_empCpds

def ri_penalty_function(abs_ri_delta, min_delta=1, max_delta=100):
    """Penalty by RI distance, using a 1-side sigmoid function
    """
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
    """
    penalized_distance = Pearsonr * ri_penalty_function
    Slow, not good for large number of features.
    """
    return [
        p for p in selected_features if
        feature_dataframe.loc[p['id'], :].corr(feature_dataframe.loc[seed_feature['id'], :]
                                ) * ri_penalty_function(
                                    abs(p['RI'] - seed_feature['RI']), min_ri_delta, max_ri_delta
                                ) >= feature_distance_filter
    ]

def get_matched_features_per_cpd(cpd, mz_centurion_tree,
                    mz_tolerance_da, ri_tolerance
                    ):
    """Get all expt features that match to any peaks in a library compound.
    cpd : GC_lib_entry instance.
    return list of matched feature IDs.
    """
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
    """Select core set of correlated features from matched_features.
    """
    features_intensity = feature_dataframe.loc[matched_feature_ids, :].T
    features_intensity_corr = features_intensity.corr()
    quant_feature = features_intensity_corr.sum(axis=0).idxmax()
    selected = [feat for feat in matched_feature_ids if features_intensity_corr.loc[quant_feature, feat] > corr_cutoff]
    return quant_feature, selected

def calculate_molecular_ion(exact_mass, mz_shift=-0.000549):
    """Assuming M[+] as molecular ion in EI GC-MS.
    Change mz_shift to 1.00728 for M+H[+]
    """
    return exact_mass + mz_shift

def have_basepeak_molecularion(base_mz, mole_mz, _peaks, mz_tolerance_da=0.005):
    """Check if base peak and molecular ion m/z are present in _peaks."""
    has_basepeak, has_molecularion = False, False
    if abs(_peaks[:, 0] - base_mz).min() < mz_tolerance_da:
        has_basepeak = True
    if abs(_peaks[:, 0] - mole_mz).min() < mz_tolerance_da:
        has_molecularion = True
    return has_basepeak, has_molecularion

def has_basepeak(base_mz, _peaks, mz_tolerance_da=0.005):
    """Check if base_mz is present in _peaks within mz_tolerance_da."""
    if abs(_peaks[:, 0] - base_mz).min() < mz_tolerance_da:
        return True
    else:
        return False

def has_molecular_ion(exact_mass, _peaks, mz_tolerance_da=0.005):
    """Check if the molecular ion of exact_mass is present in _peaks; 'NA' if no exact_mass."""
    if exact_mass:
        mole_mz = calculate_molecular_ion(exact_mass, mz_shift=-0.000549)
        if abs(_peaks[:, 0] - mole_mz).min() < mz_tolerance_da:
            return True
        else:
            return False
    else:
        return 'NA'


def unit_mz(mz):
    """Plain rounding (half up) to unit m/z, not Python's banker's rounding."""
    return int(np.floor(mz + 0.5))


def rtime_to_ri(rtime, ri_model):
    """rtime in seconds -> RI by interpolation on ri_model; None if no model."""
    if ri_model is None:
        return None
    return float(np.interp(rtime / 60, ri_model[:, 0], ri_model[:, 1]))

def append_kovats_index(list_features, ri_model):
    """Add RI to each feature in list_features (asari rtime in seconds)."""
    for feature in list_features:
        feature['RI'] = rtime_to_ri(feature['rtime'], ri_model)
    return list_features

def check_RI_presence(infile):
    """
    Check if a retention index column is present in an asari feature table (tsv).
    Matches 'RI' as a word token (e.g. 'RI', 'ri', 'R.I.', 'RI_kovats', 'kovats RI')
    or 'retention index' / 'kovats' in any case or separator style.
    Returns the first matching column name, or None if not found.
    """
    header = pd.read_csv(infile, sep='\t', nrows=0).columns
    for col in header:
        name = str(col).strip().lower()
        tokens = re.split(r'[^a-z0-9]+', name)
        compact = ''.join(tokens)
        if 'ri' in tokens or compact == 'ri' or 'retentionindex' in compact or 'kovats' in compact:
            return col
    return None


def read_feature_table_has_RI(infile, data_start_col=12):
    """
    Read an asari feature table (tsv) that already has a retention index column
    (detected by check_RI_presence), e.g. one column more than full_Feature_table.tsv.
    Returns list_features (sorted by peak_area, descending; with 'RI' from the table)
    and feature_dataframe (sample intensities, indexed by feature ID).
    Missing RI values are set to None.
    """
    df = pd.read_csv(infile, sep='\t', index_col=0)
    ri_col = check_RI_presence(infile)
    if ri_col is None:
        raise ValueError(f"No retention index column found in {infile}.")
    meta = df.iloc[:, :data_start_col - 1]
    if ri_col not in meta.columns:
        raise ValueError(f"RI column '{ri_col}' is not within the first {data_start_col} columns; check data_start_col.")
    list_features = []
    for fid, row in meta.iterrows():
        ri = row[ri_col]
        list_features.append({
            'id_number': fid,
            'id': fid,
            'mz': float(row['mz']),
            'rtime': float(row['rtime']),
            'apex': float(row['rtime']),
            'rtime_left_base': float(row['rtime_left_base']),
            'rtime_right_base': float(row['rtime_right_base']),
            'parent_masstrack_id': int(row['parent_masstrack_id']),
            'peak_area': float(row['peak_area']),
            'cSelectivity': float(row['cSelectivity']),
            'goodness_fitting': float(row['goodness_fitting']),
            'snr': float(row['snr']),
            'detection_counts': int(row['detection_counts']),
            'RI': None if pd.isna(ri) else float(ri),
        })
    list_features.sort(key=lambda x: x['peak_area'], reverse=True)

    feature_dataframe = df.iloc[:, data_start_col - 1:]

    return list_features, feature_dataframe


def read_feature_table_and_kovats(infile, KovatsIndex_file=None, frac=0.3, data_start_col=11):
    """
    Read an asari full_Feature_table.tsv.
    Returns list_features (sorted by peak_area, descending; with 'RI' if KovatsIndex_file),
    feature_dataframe (sample intensities, indexed by feature ID) and ri_model (or None).
    """
    df = pd.read_csv(infile, sep='\t', index_col=0)
    ri_model = read_fit_KovatsIndex_rtime(KovatsIndex_file, frac=frac) if KovatsIndex_file else None
    meta = df.iloc[:, :data_start_col - 1]
    list_features = []
    for fid, row in zip(meta.index, meta.itertuples(index=False)):
        list_features.append({
            'id_number': fid,
            'id': fid,
            'mz': float(row.mz),
            'rtime': float(row.rtime),
            'apex': float(row.rtime),
            'rtime_left_base': float(row.rtime_left_base),
            'rtime_right_base': float(row.rtime_right_base),
            'parent_masstrack_id': int(row.parent_masstrack_id),
            'peak_area': float(row.peak_area),
            'cSelectivity': float(row.cSelectivity),
            'goodness_fitting': float(row.goodness_fitting),
            'snr': float(row.snr),
            'detection_counts': int(row.detection_counts),
        })
    list_features = append_kovats_index(list_features, ri_model)
    list_features.sort(key=lambda x: x['peak_area'], reverse=True)

    feature_dataframe = df.iloc[:, data_start_col - 1:]

    return list_features, feature_dataframe, ri_model


def write_list_pseudospectra(
    list_pseudospectra,
    json_output_file
    ):
    """Export list of PseudoSpectrum to JSON and MSP files."""
    list_pseudospectra = [port_pseudospectrum_to_json(x, normalize_intensity=True) for x in list_pseudospectra]
    with open(json_output_file, 'w', encoding='utf-8') as O:
        json.dump(list_pseudospectra, O, indent=2, ensure_ascii=False)
    json_pseudospectra_to_msp(list_pseudospectra, json_output_file.replace('.json', '.msp'))
    print(f"Exported {len(list_pseudospectra)} de novo pseudospectra to JSON and MSP formats.\n")


#
# Workflows
#

def batch_lib_reverse_search_score(
    list_cpds,
    list_features,
    dict_features,
    feature_dataframe,
    mz_tolerance_da=0.005,
    ri_tolerance=30,
    cosine_penalty=1,   # not used
    corr_cutoff=0.6,
    ):
    """Match compound library of spectra to features from an experiment,
    using a reverse search algorithm.
    The featuers within mz and RI window are grouped as a pseudospectrum,
    which excludes features that do not correlated with the core set.
    cosine_score and entropy_score are calculated on each pseudospectrum.
    """
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
            # select single sample
            best_sample = feature_dataframe.loc[quant_feature, :].idxmax()
            _peaks = np.array([[dict_features[feat]['mz'],
                                feature_dataframe.loc[feat, best_sample]]
                                # dict_features[feat]['peak_area']] # if using composite peak area
                            for feat in _matched_feature_ids])
            if _peaks.any():
                # calculate_entropy_similarity
                entropy_score = ME.calculate_entropy_similarity(
                    _peaks, cpd.peaks,
                    ms2_tolerance_in_da = mz_tolerance_da,
                )
                # cosine_similarity
                cosine_score, num_matched_peaks = cosine_similarity(
                    _peaks, cpd.peaks,
                    mz_tolerance=mz_tolerance_da,
                )
                Results.append({
                    'lib_entry': cpd,
                    'quant_feature': quant_feature,
                    'sample': best_sample,
                    'quant_feature_RI': dict_features[quant_feature]['RI'],
                    'candidate_feature_ids': _matched_feature_ids,
                    'entropy_score': entropy_score,
                    'cosine_score': cosine_score,
                    'pseudo_spec': _peaks,
                    'num_matched_peaks': num_matched_peaks  # should = len(_matched_features)
                })
    print(f"{len(Results)} matched to library compounds.\n")
    return Results


def curate_batch_lib_reverse_search_result(
    matched_results,
    mz_tolerance_da=0.005,
    score_cutoff_cosine=0.5,
    score_cutoff_entropy=0.4,
    ):
    """Filter and organize results for export.

    determine quant_ion: using the feature of highest sum(correlations to others) in all matched features.

    matched_results : a list of dicts of matched lib_entry and pseudo_spec, from batch_lib_reverse_search_score.
    feature_dataframe : dataframe from feature table, indexed on feaature IDs, sample intensity starts in first col.

    Data for mirror plot are in (peaks_as_features, peaks_in_lib).

    return list_empCpds, feature_anno_list
    """
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
            _has_basepeak = has_basepeak(MM['lib_entry'].base_peak[0],
                peaks_as_features, mz_tolerance_da)
            _has_molecular_ion = has_molecular_ion(MM['lib_entry'].exact_mass,
                peaks_as_features, mz_tolerance_da)
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
                'sample': MM['sample'],
                'peaks_as_features': peaks_as_features,
                'peaks_in_lib': peaks_in_lib,
                'has_basepeak': _has_basepeak,
                'has_molecularion': _has_molecular_ion
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

def _lib_entry_view(lib_entry):
    """Uniform access to a library entry given as GC_lib_entry or dict (e.g. parsed MSP).
    Returns id, name, inchikey, exact_mass, peaks (np.array), base_peak_mz."""
    if isinstance(lib_entry, dict):
        peaks = np.asarray(lib_entry.get('peaks', []), dtype=float).reshape(-1, 2)
        base = lib_entry.get('base_peak')
        exact_mass = lib_entry.get('exact_mass', lib_entry.get('ExactMass'))
        return {
            'id': lib_entry.get('id', lib_entry.get('DB#')),
            'name': lib_entry.get('name', lib_entry.get('Name')),
            'inchikey': lib_entry.get('inchikey', lib_entry.get('InChIKey')),
            'exact_mass': float(exact_mass) if exact_mass not in (None, '') else None,
            'peaks': peaks,
            'base_peak_mz': base[0] if base else (peaks[peaks[:, 1].argmax(), 0] if len(peaks) else None),
        }
    return {
        'id': lib_entry.id, 'name': lib_entry.name, 'inchikey': lib_entry.inchikey,
        'exact_mass': lib_entry.exact_mass, 'peaks': np.asarray(lib_entry.peaks, dtype=float),
        'base_peak_mz': lib_entry.base_peak[0] if lib_entry.base_peak else None,
    }


def curate_batch_lib_forward_search_result(
    matched_results,
    mz_tolerance_da=0.005,
    score_cutoff_cosine=0.5,
    score_cutoff_entropy=0.4,
    lowres_mz_tolerance_da=0.5,
    ):
    """
    Filter and organize forward search results for export, same output format as
    curate_batch_lib_reverse_search_result.

    matched_results : list of dicts from batch_forward_unit_search_score (low res)
        or batch_forward_hires_search_score (high res). Each has 'query_entry'
        (a pseudospectrum in JSON format, with 'id', 'RI', 'members', 'peaks').
        Low-res results have no 'entropy_score'; they pass on cosine_score only,
        and base peak / molecular ion are checked with lowres_mz_tolerance_da.
        High-res results pass if entropy_score >= score_cutoff_entropy
        or cosine_score >= score_cutoff_cosine, as in the reverse search.

    quant_ion : the member feature with the highest intensity in the pseudospectrum
        (the base feature of deconvolution). 'sample' is None (not determined in forward search).

    Data for mirror plot are in (peaks_as_features, peaks_in_lib).

    return list_empCpds, feature_anno_list
    """
    ii = 0
    list_empCpds = []
    feature_anno_list = []
    for MM in matched_results:
        entropy_score = MM.get('entropy_score')
        cosine_score = MM.get('cosine_score', 0.0)
        highres = entropy_score is not None
        if highres:
            passed = entropy_score >= score_cutoff_entropy or cosine_score >= score_cutoff_cosine
            tol = mz_tolerance_da
        else:
            passed = cosine_score >= score_cutoff_cosine
            tol = lowres_mz_tolerance_da
        if not passed:
            continue

        query = MM['query_entry']
        lib = _lib_entry_view(MM['lib_entry'])
        peaks_as_features = np.asarray(query.get('peaks', []), dtype=float).reshape(-1, 2)
        matched_features = list(query.get('members') or [])
        if matched_features and len(matched_features) == len(peaks_as_features):
            quant_feature = matched_features[int(peaks_as_features[:, 1].argmax())]
        else:
            quant_feature = query.get('id')
        if len(peaks_as_features):
            _has_basepeak = has_basepeak(lib['base_peak_mz'], peaks_as_features, tol) \
                if lib['base_peak_mz'] is not None else 'NA'
            _has_molecular_ion = has_molecular_ion(lib['exact_mass'], peaks_as_features, tol)
        else:
            _has_basepeak, _has_molecular_ion = False, False

        ii += 1
        epd_id = 'empCpd_' + "{:04d}".format(ii)
        list_empCpds.append({
            'id': epd_id,
            'pseudospectrum_id': query.get('id'),
            'rtime': query.get('rtime'),
            'RI': query.get('RI'),          # use expt data, not lib data
            'quant_ion_mz': float(peaks_as_features[peaks_as_features[:, 1].argmax(), 0]) if len(peaks_as_features) else None,
            'entropy_score': entropy_score,
            'cosine_score': cosine_score,
            'num_matched_peaks': MM.get('num_matched_peaks'),
            'lib_entry_id': lib['id'],
            'name': lib['name'],
            'inchikey': lib['inchikey'],
            'features': matched_features,
            'quant_ion': quant_feature,
            'sample': None,
            'peaks_as_features': peaks_as_features,
            'peaks_in_lib': lib['peaks'],
            'has_basepeak': _has_basepeak,
            'has_molecularion': _has_molecular_ion,
        })
        for feature in matched_features:
            feature_anno_list.append({
                'feature': feature,
                'empCpd': epd_id,
                'quant_ion': quant_feature,
                'entropy_score': entropy_score,
                'cosine_score': cosine_score,
                'lib_entry_id': lib['id'],
                'name': lib['name'],
                'inchikey': lib['inchikey'],
            })

    return list_empCpds, feature_anno_list


#
# Low resolution functions
#

# ---------------------------------------------------------------------------
# Unit m/z conversion
# ---------------------------------------------------------------------------

def _to_unit(peaks):
    """
    peaks: iterable of [mz, intensity].
    Returns (unit_mz int32 sorted unique, unit_intensity float64, L2-normalized).
    Intensities falling into the same nominal mass are summed.
    """
    arr = np.asarray(peaks, dtype=np.float64).reshape(-1, 2)
    arr = arr[arr[:, 1] > 0]
    if arr.shape[0] == 0:
        return np.empty(0, np.int32), np.empty(0, np.float64)
    # plain rounding (half up), not Python's banker's rounding
    mz = np.floor(arr[:, 0] + 0.5).astype(np.int32)
    unit_mz, inv = np.unique(mz, return_inverse=True)
    inten = np.bincount(inv, weights=arr[:, 1])
    inten /= np.sqrt(np.dot(inten, inten))
    return unit_mz, inten


def convert_unit_mz(list_spectra):
    """
    Convert spectra to unit m/z. Exact m/z is converted to unit mass by plain rounding.

    Applies to both query spectra and library spectra.
    Input: list of spectra, each a dict with 'peaks', a GC_lib_entry, or a list of [mz, intensity].
    Returns a new list of dicts (input is not modified): the original fields plus
        'unit_mz': sorted unique int32 array of nominal masses,
        'unit_intensity': float64 array, summed per nominal mass and L2-normalized,
    so that cosine similarity reduces to a dot product over shared integers.
    """
    result = []
    for spec in list_spectra:
        if isinstance(spec, dict):
            new = dict(spec)
            peaks = spec.get('peaks', [])
        elif hasattr(spec, '_asdict'):          # GC_lib_entry or other NamedTuple
            new = spec._asdict()
            peaks = new.get('peaks', [])
        else:
            new, peaks = {}, spec
        new['unit_mz'], new['unit_intensity'] = _to_unit(peaks)
        result.append(new)
    return result


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _get_unit(spec):
    if isinstance(spec, dict):
        if 'unit_mz' not in spec:
            return _to_unit(spec.get('peaks', []))
        return spec['unit_mz'], spec['unit_intensity']
    return _to_unit(spec)


def unit_search_cosine(spec1, spec2):
    """
    The function accepts two spectra, preferably already converted by convert_unit_mz
    (raw dicts with 'peaks' or [mz, intensity] lists are converted on the fly, slower).

    Because unit m/z values are sorted unique integers and intensities are
    pre-normalized, cosine = sum of intensity products over shared m/z.

    The score is a float in [0, 1], where 1.0 is a perfect match.
    Returns score, num_matched_peaks
    """
    mz1, int1 = _get_unit(spec1)
    mz2, int2 = _get_unit(spec2)
    if mz1.size == 0 or mz2.size == 0:
        return 0.0, 0
    # search the shorter array in the longer one
    if mz1.size > mz2.size:
        mz1, int1, mz2, int2 = mz2, int2, mz1, int1
    idx = np.searchsorted(mz2, mz1)
    idx[idx == mz2.size] = 0
    hit = mz2[idx] == mz1
    n = int(hit.sum())
    if n == 0:
        return 0.0, 0
    score = float(np.dot(int1[hit], int2[idx[hit]]))
    return min(score, 1.0), n


def _build_matrix(unit_spectra, n_cols):
    """CSR matrix (n_spectra x n_cols) of normalized unit intensities."""
    lengths = np.fromiter((s['unit_mz'].size for s in unit_spectra), dtype=np.int64,
                          count=len(unit_spectra))
    indptr = np.concatenate([[0], np.cumsum(lengths)])
    if indptr[-1] == 0:
        return sparse.csr_matrix((len(unit_spectra), n_cols))
    indices = np.concatenate([s['unit_mz'] for s in unit_spectra])
    data = np.concatenate([s['unit_intensity'] for s in unit_spectra])
    return sparse.csr_matrix((data, indices, indptr), shape=(len(unit_spectra), n_cols))


def batch_forward_unit_search_score(
        query_list,
        list_cpds,  # library
        ri_tolerance=30,
        include_missing_ri=True,
        min_matched_peaks=3,
        top_n=1,
        chunk_size=2000,
    ):
    """
    Match EI spectra between query list and a library of known compounds.

    Limit spectral comparison to the window within ri_tolerance.
    - ri_tolerance=None disables the RI filter.
    - include_missing_ri: library entries (or queries) without RI cannot be ruled
      out by RI; keep them as candidates if True (MoNA: most entries lack RI).

    Both inputs are converted to unit m/z if not already done by convert_unit_mz.
    All scores are computed as sparse matrix products, in chunks of queries:
        cosine  = Q @ L.T            (rows L2-normalized)
        matched = (Q>0) @ (L>0).T

    Returns one dict per query that has a match (cosine > 0 and
    num_matched_peaks >= min_matched_peaks), with the best match fields at top
    level and up to top_n candidates in 'candidates'.
    """
    if query_list and 'unit_mz' not in query_list[0]:
        query_list = convert_unit_mz(query_list)
    if list_cpds and 'unit_mz' not in list_cpds[0]:
        list_cpds = convert_unit_mz(list_cpds)
    if not query_list or not list_cpds:
        return []

    max_mz = 1 + max(
        max((int(s['unit_mz'][-1]) for s in query_list if s['unit_mz'].size), default=0),
        max((int(s['unit_mz'][-1]) for s in list_cpds if s['unit_mz'].size), default=0),
    )
    L = _build_matrix(list_cpds, max_mz)
    LT = L.T.tocsc()
    LT_bin = LT.copy()
    LT_bin.data[:] = 1.0

    lib_ri = np.array([np.nan if s.get('RI') is None else s['RI'] for s in list_cpds])
    lib_ri_missing = np.isnan(lib_ri)

    Results = []
    for start in range(0, len(query_list), chunk_size):
        chunk = query_list[start:start + chunk_size]
        Q = _build_matrix(chunk, max_mz)
        Q_bin = Q.copy()
        Q_bin.data[:] = 1.0
        scores = (Q @ LT).toarray()
        matched = (Q_bin @ LT_bin).toarray()

        if ri_tolerance is not None:
            q_ri = np.array([np.nan if s.get('RI') is None else s['RI'] for s in chunk])
            with np.errstate(invalid='ignore'):
                in_window = np.abs(q_ri[:, None] - lib_ri[None, :]) <= ri_tolerance
            if include_missing_ri:
                in_window |= lib_ri_missing[None, :] | np.isnan(q_ri)[:, None]
            scores[~in_window] = 0.0
        scores[matched < min_matched_peaks] = 0.0

        for i, query in enumerate(chunk):
            row = scores[i]
            if top_n == 1:
                best = np.array([row.argmax()])
            else:
                k = min(top_n, row.size)
                best = np.argpartition(-row, k - 1)[:k]
                best = best[np.argsort(-row[best], kind='stable')]
            best = best[row[best] > 0]
            if best.size == 0:
                continue
            candidates = [{
                'lib_index': int(j),
                'lib_entry': list_cpds[j],
                'lib_name': _spec_field(list_cpds[j], 'Name'),
                'lib_RI': list_cpds[j].get('RI'),
                'cosine_score': float(min(row[j], 1.0)),
                'num_matched_peaks': int(matched[i, j]),
            } for j in best]
            top = candidates[0]
            Results.append({
                'query_id': query.get('id'),
                'query_entry': query,
                'query_RI': query.get('RI'),
                'lib_entry': top['lib_entry'],
                'lib_name': top['lib_name'],
                'lib_RI': top['lib_RI'],
                'cosine_score': top['cosine_score'],
                'num_matched_peaks': top['num_matched_peaks'],
                'candidates': candidates,
            })
    return Results


def _spec_field(spec, key):
    """Read a field from a spectrum given as dict or GC_lib_entry (NamedTuple)."""
    if isinstance(spec, dict):
        if key == 'Name':
            return spec.get('Name', spec.get('name'))
        return spec.get(key)
    return getattr(spec, {'Name': 'name'}.get(key, key), None)


def _spec_ri_array(specs):
    ri = (_spec_field(s, 'RI') for s in specs)
    return np.array([np.nan if (x is None or x == '') else float(x) for x in ri], dtype=float)


def batch_forward_hires_search_score(
        query_list,
        list_cpds,  # library
        mz_tolerance_da=0.005,
        ri_tolerance=30,
        include_missing_ri=True,
        min_matched_peaks=3,
        top_n=1,
        chunk_size=500,
        noise_threshold=0.01,
        n_threads=1,
    ):
    """
    Match high-resolution EI spectra between query list (e.g. pseudospectra) and
    a library of known compounds. High-res counterpart of batch_forward_unit_search_score.

    Spectra are dicts with 'peaks' ([[mz, intensity], ...]) and optional 'RI',
    or GC_lib_entry for the library.

    Scoring uses ms2batch open search (no precursor in EI): entropy similarity over
    peaks matched within mz_tolerance_da, computed in batch for each chunk of queries.
    Candidates are filtered by RI window and min_matched_peaks, ranked by entropy score;
    greedy cosine (tools.cosine.cosine_similarity) is then computed for the top_n kept.
    - ri_tolerance=None disables the RI filter.
    - include_missing_ri: entries without RI (library or query) are kept as candidates.
    - mz_tolerance_da must be < 0.025 (half of ms2batch min_peak_distance 0.05).

    Returns one dict per query that has a match, with the best match fields at top
    level and up to top_n candidates in 'candidates'. Same keys as the unit version,
    plus 'entropy_score'.
    """
    from .ms2batch import SpectralIndex, search

    if not query_list or not list_cpds:
        return []
    params = dict(noise_threshold=noise_threshold, precursor_removal_da=None)
    lib_peaks = [np.asarray(_spec_field(s, 'peaks'), dtype=float).reshape(-1, 2) for s in list_cpds]
    library = SpectralIndex.build([{'peaks': p} for p in lib_peaks],
                                  fragment_index=True, keep_metadata=False, **params)
    lib_ri = _spec_ri_array(list_cpds)
    lib_ri_missing = np.isnan(lib_ri)

    Results = []
    for start in range(0, len(query_list), chunk_size):
        chunk = query_list[start:start + chunk_size]
        q_peaks = [np.asarray(_spec_field(s, 'peaks') or [], dtype=float).reshape(-1, 2) for s in chunk]
        hits = search(library, [{'peaks': p} for p in q_peaks], method='open',
                      ms2_tolerance=mz_tolerance_da, top_k=None, n_threads=n_threads)
        qi, li, score, nm = hits.query_idx, hits.library_idx, hits.score, hits.n_matched

        keep = nm >= min_matched_peaks
        if ri_tolerance is not None:
            q_ri = _spec_ri_array(chunk)
            with np.errstate(invalid='ignore'):
                ok = np.abs(q_ri[qi] - lib_ri[li]) <= ri_tolerance
            if include_missing_ri:
                ok |= lib_ri_missing[li] | np.isnan(q_ri[qi])
            keep &= ok
        qi, li, score, nm = qi[keep], li[keep], score[keep], nm[keep]
        if qi.size == 0:
            continue

        # hits are grouped by query, best score first; take top_n per query
        first = np.flatnonzero(np.r_[True, qi[1:] != qi[:-1]])
        ends = np.r_[first[1:], qi.size]
        for a, b in zip(first, ends):
            i = int(qi[a])
            query = chunk[i]
            candidates = []
            for k in range(a, min(b, a + top_n)):
                j = int(li[k])
                cos, _ = cosine_similarity(q_peaks[i], lib_peaks[j], mz_tolerance=mz_tolerance_da)
                candidates.append({
                    'lib_index': j,
                    'lib_entry': list_cpds[j],
                    'lib_name': _spec_field(list_cpds[j], 'Name'),
                    'lib_RI': _spec_field(list_cpds[j], 'RI'),
                    'entropy_score': float(min(score[k], 1.0)),
                    'cosine_score': float(cos),
                    'num_matched_peaks': int(nm[k]),
                })
            top = candidates[0]
            Results.append({
                'query_id': _spec_field(query, 'id'),
                'query_entry': query,
                'query_RI': _spec_field(query, 'RI'),
                'lib_entry': top['lib_entry'],
                'lib_name': top['lib_name'],
                'lib_RI': top['lib_RI'],
                'entropy_score': top['entropy_score'],
                'cosine_score': top['cosine_score'],
                'num_matched_peaks': top['num_matched_peaks'],
                'candidates': candidates,
            })
    return Results
