"""
join as 1st class module, as it will be used to divide & conquer large studies.

The mass alignment is same as landmark_guided_mapping;
RT alignment is as default LOESS.
Merge feature tables and zero fill absent features in samples.
Track input feature IDs and parent_masstrack_id; keep mapping record.

Example use:
python3 -m asari.main join -i /Users/lish/li.play/join_MT_files.txt -o /Users/lish/li.play/ -j test_join_mt01
"""

import os
import time
import json
import numpy as np
import pandas as pd

from scipy import interpolate
from scipy.signal import find_peaks
from statsmodels.nonparametric.smoothers_lowess import lowess

from mass2chem.search import find_mzdiff_pairs_from_masstracks

from .mass_functions import (flatten_tuplelist,
                             landmark_guided_dict_mapping,
                             check_close_mzs)
from .experiment import ext_Experiment
from .constructors import (CompositeMap, MassGrid)
from .tools.file_io import read_features_from_asari_table


def join_tables(file_list_paths, parameters):
    """
    This joins data from multiple feature tables using the same LC/GC-MS method.
    Batch correction of intensities should be performed after joining.

    file_list_paths : a text file containing list of paths, one per line
    parameters : dict, containing parameters as typical Asari processing

    Treat each table as a virtual "sample".
    Build a 'MassGrid' first by aligning masstracks btw tables;
    RT alignment is performed on features that are high and unique per mass track.
    Elution peaks are reconstructed using Gaussian simulation, and peak detection is performed on the composite mass track.
    Features on each composite mass track are mapped back to the original tables.

    The IDs of mass tracks are not compatible with default Asari processing,
    because the latter uses positional indexes but this function reuses old IDs.
    Therefore, the sample-wise MassGrid building is standalone version here.

    The landmark_guided_dict_mapping function has mass calibration,
    triggered by correction_tolerance_ppm (1 default).

    Format of composite feature is like:
    {'id': 'F00000501m', 'mz': 284.2945, 'rtime': 1116.97, 'parent_masstrack_id': 107, 'height': 150503047, 'rtime_left_base': 1115.77, 'rtime_right_base': 1118.17,
    'members': {'3A5_H_full_Feature_table.tsv': [{'id': 'F388', 'mz': 284.2947, 'apex': 1116.94}],
    '2C19_H_full_Feature_table.tsv': [{'id': 'F366', 'mz': 284.2946, 'apex': 1116.55}],
    '2B6_H_full_Feature_table.tsv': [{'id': 'F360', 'mz': 284.2945, 'apex': 1116.81}],
    'WT_H_full_Feature_table.tsv': []}}

    TO-DO:
    options to use standards-calibrated RT alignment or turn off.
    E.g. environmental samples may not have consistent landmarks to align;
    GC RI can be applied per batch, no need to realign in join.

    """
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    subdir = '_'.join(['join', parameters['project_name'], ''.join(time_stamp)])
    outdir = os.path.join(parameters['outdir'], subdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"Output directory: {outdir}.\n\n")

    # Scan numbers are not recorded in the feature tables.
    # We use arbituary number max_scan_number here to mock raw data processing.
    max_scan_number = 10000         # mock number
    tables = [read_table_to_json(table) for table in verify_path_list(file_list_paths)]
    tables.sort(key=lambda x: x['number_features'], reverse=True)
    for ii, table in enumerate(tables):
        anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(
            table['sample_data']['list_mass_tracks'],
            mz_tolerance_ppm=parameters['mz_tolerance_ppm']
            )
        table['anchor_mz_pairs'] = anchor_mz_pairs
        table['_mz_landmarks_'] = flatten_tuplelist(anchor_mz_pairs)
        table['input_file'] = table['name']
        table['number_anchor_mz_pairs'] = len(anchor_mz_pairs)
        table['sample_id'] = ii
        table['data_location'] = None
        table['track_mzs'] = [(x['mz'], x['id_number']) for x in table['sample_data']['list_mass_tracks']]
        table['max_scan_number'] = max_scan_number
        table['list_scan_numbers'] = list(range(1, max_scan_number+1))
        max_rtime = max([f['rtime'] for f in table['sample_data']['list_features']])
        table['list_retention_time'] = np.linspace(0, max_rtime, num=max_scan_number).tolist()

    sample_registry = {d['sample_id'] : d for d in tables}
    EE = ext_Experiment(sample_registry, parameters)
    EE.database_mode = 'memory'
    EE.reference_sample_id = 0
    EE.valid_sample_ids = list(sample_registry.keys())
    EE.number_of_samples = len(tables)
    EE.number_scans = max_scan_number

    CMAP = CompositeMap(EE)
    MG = MassGrid( CMAP, CMAP.experiment )  # MG.build_grid_sample_wise()
    MG._initiate_mass_grid()
    # the _initiate_mass_grid step assigned reference_sample to MG.MassGrid
    # we now update REF_reference_mzdict using MG.MassGrid indexes
    REF_reference_mzdict = {ii: mz for ii,mz in enumerate(MG.MassGrid['mz'])}
    reverse_ref_dict = {x['id_number']:ii for ii,x in enumerate(MG.reference_sample_instance.list_mass_tracks)}
    # redo index
    REF_landmarks = [reverse_ref_dict[x] for x in MG._mz_landmarks_]
    print("Scan number is arbitrary in the `join` process. \nMassGrid alignment and composite peak detection are similar to raw data processing.\n")

    for table in tables[1:]:
        print(f"Adding to MassGrid {table['input_file']}.")
        _N1 = len(REF_reference_mzdict)
        new_reference_map2 = {}
        # similar flow as MG.add_sample() but using dict as input
        SM_mzdict = {x['id_number']: x['mz']
                     for x in table['sample_data']['list_mass_tracks']}
        mapped_pairs, dict2_unmapped, _r = landmark_guided_dict_mapping(
            REF_reference_mzdict, REF_landmarks, SM_mzdict, table['_mz_landmarks_']
        )
        for p in mapped_pairs:  # updating ref m/z here
            REF_reference_mzdict[p[0]] = 0.5*( REF_reference_mzdict[p[0]] + SM_mzdict[p[1]] )
            new_reference_map2[p[0]] = p[1]
        for ii,idx in enumerate(dict2_unmapped):
            REF_reference_mzdict[_N1 + ii] = SM_mzdict[idx]
            new_reference_map2[_N1 + ii] = idx
            # update landmark m/z values using the new index numbers as part of new_reference_mzlist
            # Not important to add more landmarks.
            # if idx in table['_mz_landmarks_']:
            #    REF_landmarks.append(_N1 + ii)

        # Update MassGrid
        new_num_rows = len(REF_reference_mzdict)
        NewGrid = pd.DataFrame(
            np.full((new_num_rows, 1+EE.number_of_samples), None),
            columns=['mz'] + MG.list_sample_names,
        )
        NewGrid[ :MG.MassGrid.shape[0]] = MG.MassGrid       # replicate last grid
        NewGrid['mz'] = [REF_reference_mzdict[ii] for ii in range(new_num_rows)]    # update ref m/z list
        NewGrid[ table['name'] ] = [new_reference_map2.get(ii, None) for ii in range(new_num_rows)]   # update col of this sample
        NewGrid[ table['name'] ] = NewGrid[ table['name'] ].astype('Int64')
        MG.MassGrid = NewGrid

    MG.MassGrid.to_csv(os.path.join(outdir, 'MassGrid.csv'))

    # Linked thru (Feature ID per table - parent_masstrack_id - MassGrid)
    # Do RT alignment as if table were sample; start by selecting good_landmark_peaks
    #
    MIN_PEAK_AREA = np.mean([f['peak_area'] for f in tables[0]['sample_data']['list_features']])
    reference_landmark_peaks = collect_rt_calibration_features(
        tables[0]['sample_data']['list_features'], min_area=MIN_PEAK_AREA, min_snr=10
    )           # MassGrid iloc can be looked up in reverse_ref_dict
    for table in tables[1:]:
        mass_grid_values = MG.MassGrid[table['name']].values        # contains NA
        candidate_landmark_peaks = collect_rt_calibration_features(
            table['sample_data']['list_features'], min_area=MIN_PEAK_AREA, min_snr=10
        )
        track_to_feature = {f['parent_masstrack_id']: f for f in candidate_landmark_peaks}
        ref_features, work_features = [], []
        for feat in reference_landmark_peaks:
            ii = mass_grid_values[reverse_ref_dict[feat['parent_masstrack_id']]]
            if not pd.isna(ii):
                ii = int(ii)
                if ii in track_to_feature:
                    ref_features.append(feat)
                    work_features.append(track_to_feature[ii])
        smoothed = lowess(
            [f['rtime'] for f in ref_features], [f['rtime'] for f in work_features], frac=0.2, it=2
        )
        interf = interpolate.interp1d(smoothed[:, 0], smoothed[:, 1], fill_value="extrapolate")
        sample_rt_values = [f['rtime'] for f in table['sample_data']['list_features']]
        new_rt_values = interf( sample_rt_values )
        for ii, feat in enumerate(table['sample_data']['list_features']):
            feat['rtime'] = new_rt_values[ii]

    #
    master_track2features = {}      # 2-level dict: {table_name: {parent_masstrack_id: [features]}}
    for table in tables:
        _d = {'': []}           # catch NA
        for feat in  table['sample_data']['list_features']:
            if feat['parent_masstrack_id'] in _d:
                _d[feat['parent_masstrack_id']].append(feat)
            else:
                _d[feat['parent_masstrack_id']] = [feat]
        master_track2features[table['name']] = _d

    print("\nCombining features and detecting elution peaks on the composite mass tracks.\n\n")

    #
    # def combined tracks, to merge tracks of close m/z values, by parameters['mz_tolerance_ppm'] (default 5 ppm)
    #
    operational_tracks, operational_tracks_mz = {}, {}
    mztracks = list(zip(MG.MassGrid['mz'].tolist(), MG.MassGrid.index))     # [(83.0603, 0), (84.0443, 1), ...]
    mztracks.sort()
    tracks_to_merge = check_close_mzs([x[0] for x in mztracks], parameters['mz_tolerance_ppm'])
    to_merge_clusters, _registered_ = {}, {}
    for ii, (a,b) in enumerate(tracks_to_merge):
        if a in _registered_:
            to_merge_clusters[_registered_[a]].append(mztracks[b])
        elif b in _registered_:
            to_merge_clusters[_registered_[b]].append(mztracks[a])
        else:
            to_merge_clusters[ii] = [mztracks[a], mztracks[b]]
            _registered_[a] = ii
            _registered_[b] = ii

    print(f"Number of composite mass tracks to merge due to close m/z values: {len(to_merge_clusters)}")
    for ii, track in enumerate(mztracks):           # track = (mz, index)
        if ii in _registered_:
            cluster_id = _registered_[ii]
            if cluster_id in operational_tracks:
                operational_tracks[cluster_id].append(track[1])
            else:
                operational_tracks[cluster_id] = [track[1]]
        else:
            operational_tracks[ii] = [track[1]]

    # Update feature list; decide feature correspondence on each composite mass track
    new_feature_list = []
    names = [table['name'] for table in tables]
    for jj, track_indices in operational_tracks.items():
        median_mz = np.median([MG.MassGrid['mz'][ii] for ii in track_indices])
        operational_tracks_mz[jj] = median_mz
        row_features = {name: [] for name in names}
        for ii in track_indices:
            row_dict = MG.MassGrid.iloc[ii, :].to_dict()        # mapping name to mass track index
            for t, track2feats in master_track2features.items():
                if pd.notna(row_dict[t]):
                    row_features[t] += track2feats[row_dict[t]]

        # if str(ii).endswith('1000'): print(row_features)
        new_feature_list += consolidate_track_features(
            row_features, parameters, jj, median_mz
            )

    # export operational_tracks, operational_tracks_mz to tsv file
    with open(os.path.join(outdir, 'operational_mass_tracks_mz.tsv'), 'w', encoding='utf-8') as f:
        for jj, mz in operational_tracks_mz.items():
            f.write(f"{jj}\t{mz}\t{operational_tracks[jj]}\n")

    # write mapping record; A few features would fail new peak matching, stored in empty_features
    new_feature_list, empty_features = reformat_composite_feature_list(new_feature_list, names)

    header = "id,mz,rtime,rtime_left_base,rtime_right_base,parent_masstrack_id_operational,height_combined".split(",")
    for name in names:
        header += [f'id({name})', f'mz({name})', f'apex({name})']
    s = '\t'.join(header) + '\n'
    for feat in new_feature_list:
        member_fields = []
        for name in names:
            members = feat['members'][name]
            if len(members) == 1:
                # id, mz, apex
                member_fields += [members[0]['id'], members[0]['mz'], members[0]['apex']]
            elif len(members) == 0:
                member_fields += ['', '', '']
            else:
                member_fields += [','.join([x['id'] for x in members]),
                                  members[0]['mz'], members[0]['apex']]

        s += '\t'.join([str(x) for x in [
            feat[ii] for ii in "id,mz,rtime,rtime_left_base,rtime_right_base,parent_masstrack_id,height".split(",")
        ] + member_fields
        ]) + '\n'
    with open( os.path.join(outdir, 'asari_join_record.tsv'), 'w') as O:
        O.write(s)

    #
    # Join all data, i.e. Merge feature tables with sample/intensity values
    # fill 0s if feature not present in a table
    # merge if more than 1 features are from same table

    DATADICT, sample_columns, sample_numbers = {}, [], {}
    for name in names:
        # {Feature ID: sample data} for all samples
        DATADICT[name] = extract_sample_data_dict(name)
        # check headers per table if all col names are unique?
        sample_columns += DATADICT[name]['id_number']   # 1st row actually has sample cols
        sample_numbers[name] = len(DATADICT[name]['id_number'])

    use_cols = [ 'id', 'mz', 'rtime', 'rtime_left_base', 'rtime_right_base', 'parent_masstrack_id',
                    'height', 'cSelectivity', 'goodness_fitting', 'snr', 'detection_counts' ]
    s = '\t'.join(use_cols  + sample_columns) + '\n'
    for feat in new_feature_list:
        member_fields = []
        for name in names:
            members = feat['members'][name]
            if len(members) == 1:
                member_fields += DATADICT[name][members[0]['id']]
            elif len(members) == 0:
                member_fields += [0] * sample_numbers[name]
            else:
                member_fields += join_split_features([DATADICT[name][row['id']] for row in members])

        s += '\t'.join([str(x) for x in [feat[ii] for ii in use_cols] + member_fields
        ]) + '\n'
    with open(os.path.join(outdir, 'asari_join_result.tsv'), 'w') as O:
        O.write(s)

    # export parameters
    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, ensure_ascii=False, indent=2)

    with open(os.path.join(outdir, 'failed_features.log'), 'w') as O:
        O.write('\n'.join([str(f) for f in empty_features]))


def join_split_features(list_str_intensities):
    """Sum intensities of split features, as lists of strings."""
    data = np.zeros(len(list_str_intensities[0])).astype(int)
    for L in list_str_intensities:
        data += np.array(L).astype(int)
    return [str(x) for x in data]

def extract_sample_data_dict(filepath, data_start_col=11):
    """
    filepath : Table in Asari format, which has no empty field.
    return data in str, as conversion is only needed occassionally.
    """
    d = {}
    for row in open(filepath).readlines():
        a = row.rstrip().split('\t')
        d[a[0]] = a[data_start_col: ]
    return d


def reformat_composite_feature_list(new_feature_list, table_names):
    """
    Example composite feature:
    {'mz': 284.2945, 'parent_masstrack_id': 107, 'rtime': 1116.97, 'height': 150503047, 'rtime_left_base': 1115.77, 'rtime_right_base': 1118.17,
    'members': [{'id_number': 'F388', 'id': 'F388', 'mz': 284.2947, 'rtime': 1116.94, 'apex': 1116.94, 'rtime_left_base': 1115.37, 'rtime_right_base': 1118.52, 'parent_masstrack_id': 139, 'peak_area': 75869861.0, 'cSelectivity': 0.64, 'goodness_fitting': 0.94, 'snr': 3.0, 'detection_counts': 9, 'source': '/Users/sli/play/join_test/DnCl_pellets_3A5_H_full_Feature_table.tsv'},
    {'id_number': 'F366', 'id': 'F366', 'mz': 284.2946, 'rtime': 1116.8590394349912, 'apex': 1116.55, 'rtime_left_base': 1114.97, 'rtime_right_base': 1118.12, 'parent_masstrack_id': 128, 'peak_area': 72230862.0, 'cSelectivity': 0.62, 'goodness_fitting': 0.95, 'snr': 3.0, 'detection_counts': 9, 'source': '/Users/sli/play/join_test/DnCl_pellets_2C19_H_full_Feature_table.tsv'},
    {'id_number': 'F360', 'id': 'F360', 'mz': 284.2945, 'rtime': 1116.9675010942542, 'apex': 1116.81, 'rtime_left_base': 1115.23, 'rtime_right_base': 1118.39, 'parent_masstrack_id': 127, 'peak_area': 67566816.0, 'cSelectivity': 0.54, 'goodness_fitting': 0.92, 'snr': 3.0, 'detection_counts': 9, 'source': '/Users/sli/play/join_test/DnCl_pellets_2B6_H_full_Feature_table.tsv'}]}
    """
    # parent_masstrack_id	peak_area	cSelectivity	goodness_fitting	snr	detection_counts
    def shorten(feature):
        return {
            'id': feature['id'], 'mz': feature['mz'], 'apex': feature['apex'],
        }
    new, empty_features = [], []
    for ii, feat in enumerate(new_feature_list):
        if feat['members']:
            tt = {name: [] for name in table_names}
            for f in feat['members']:
                tt[f['source']].append(shorten(f))
            new.append({
                'id': f"F{ii+1:08d}m",
                'mz': feat['mz'],
                'rtime': feat['rtime'],
                'parent_masstrack_id': feat['parent_masstrack_id'],
                'height': feat['height'],
                'rtime_left_base': feat['rtime_left_base'],
                'rtime_right_base': feat['rtime_right_base'],
                'cSelectivity': max([x['cSelectivity'] for x in feat['members']]),
                'goodness_fitting': max([x['goodness_fitting'] for x in feat['members']]),
                'snr': max([x['snr'] for x in feat['members']]),
                'detection_counts': sum([x['detection_counts'] for x in feat['members']]),
                'members': tt
            })
        else:
            # print("Empty composite feature - ", feat)
            empty_features.append(feat)

    return new, empty_features


def consolidate_track_features(dict_sample_features, parameters, id_number, mz, scan_time=0.1):
    """
    Combining signals from multiple samples on the same mass track, and detect elution peaks on the composite mass track.

    dict_sample_features : groups of features from each sample (table) on same mass track
    scan_time : interval used to resample peak shape (seconds)

    E.g.,
    {'3A5_H_full_Feature_table.tsv': [{'id_number': 'F25109', 'id': 'F25109', 'mz': 583.5966, 'rtime': 1221.2,
    'apex': 1221.2, 'rtime_left_base': 1198.36, 'rtime_right_base': 1244.05, 'parent_masstrack_id': 8030, 'peak_area': 1083540630.0, 'cSelectivity': 1.0, 'goodness_fitting': 0.98, 'snr': 106.0, 'detection_counts': 9},
                        {'id_number': 'F25110', 'id': 'F25110', 'mz': 583.5966, 'rtime': 1323.5,
                        'apex': 1323.5, 'rtime_left_base': 1308.92, 'rtime_right_base': 1323.89, 'parent_masstrack_id': 8030, 'peak_area': 15210424.0, 'cSelectivity': 0.68, 'goodness_fitting': 0.58, 'snr': 7.0, 'detection_counts': 9}],
    '2C19_H_full_Feature_table.tsv': [{'id_number': 'F23241', 'id': 'F23241', 'mz': 583.5964, 'rtime': 1220.358334439543,
    'apex': 1220.15, 'rtime_left_base': 1219.89, 'rtime_right_base': 1221.47, 'parent_masstrack_id': 7661, 'peak_area': 160745369.0, 'cSelectivity': 0.19, 'goodness_fitting': 0.66, 'snr': 8.0, 'detection_counts': 9},
                        {'id_number': 'F23242', 'id': 'F23242', 'mz': 583.5964, 'rtime': 1321.0629938933155,
                        'apex': 1320.87, 'rtime_left_base': 1320.48, 'rtime_right_base': 1321.66, 'parent_masstrack_id': 7661, 'peak_area': 1622804.0, 'cSelectivity': 0.03, 'goodness_fitting': 0.73, 'snr': 5.0, 'detection_counts': 9}],
    '2B6_H_full_Feature_table.tsv': [{'id_number': 'F23969', 'id': 'F23969', 'mz': 583.5963, 'rtime': 1219.4791091268726,
    'apex': 1219.23, 'rtime_left_base': 1218.84, 'rtime_right_base': 1219.89, 'parent_masstrack_id': 7579, 'peak_area': 189816475.0, 'cSelectivity': 0.15, 'goodness_fitting': 0.56, 'snr': 8.0, 'detection_counts': 9},
                        {'id_number': 'F23970', 'id': 'F23970', 'mz': 583.5963, 'rtime': 1220.9179464905783,
                        'apex': 1220.68, 'rtime_left_base': 1220.68, 'rtime_right_base': 1220.68, 'parent_masstrack_id': 7579, 'peak_area': 85019880.0, 'cSelectivity': 0.13, 'goodness_fitting': 0.89, 'snr': 7.0, 'detection_counts': 9}],
    'WT_H_full_Feature_table.tsv': [{'id_number': 'F22899', 'id': 'F22899', 'mz': 583.5964, 'rtime': 1223.5767717154906,
    'apex': 1223.57, 'rtime_left_base': 1223.04, 'rtime_right_base': 1225.14, 'parent_masstrack_id': 7830, 'peak_area': 157757246.0, 'cSelectivity': 0.23, 'goodness_fitting': 0.76, 'snr': 5.0, 'detection_counts': 9},
                        {'id_number': 'F22900', 'id': 'F22900', 'mz': 583.5964, 'rtime': 1320.4266479964485,
                        'apex': 1320.35, 'rtime_left_base': 1319.17, 'rtime_right_base': 1320.61, 'parent_masstrack_id': 7830, 'peak_area': 1890897.0, 'cSelectivity': 0.04, 'goodness_fitting': 0.64, 'snr': 6.0, 'detection_counts': 9}]}

    From the above examples:
    a) RT calibration did not shift much of original values, e.g. 2B6_H was  1219.23 & 1220.68
    b) Most tables have two distinct features around 1221 and 1321
    c) Split features exist, e.g. 2B6_H was  1219.23 & 1220.68 (if chromatography isn't good enough, not a problem of the software).

    We use the Composite Mass Track and elution peak detection as in LC-MS Asari to combine data.
    Elution peaks per feature are constructed using Gaussian simulation
    (e.g. parameters - 'apex': 1320.87, 'rtime_left_base': 1320.48, 'rtime_right_base': 1321.66, 'peak_area': 1622804.0).
    Composite Mass Track is sum of individual tracks. The method should take care of sporatic split features.
    Results should include like {consensus_feature:
    track_mz, rtime, {table_name: [{'id_number': 'F25109', 'id': 'F25109', 'mz': 583.5966, 'rtime': 1221.2, 'peak_area': 1083540630.0,}, ], }
    }.
    The returned composite features keep SNR as max value.

    returns list of new consolidated consensus features, with mapping to input features
    """
    # determine rtime range
    all_features = []
    for name, list_features in dict_sample_features.items():
        for feat in list_features:
            feat['source'] = name
            all_features.append(feat)

    if not all_features:
        return []
    else:
        left_base = min([feat['rtime_left_base'] for feat in all_features])
        right_base = max([feat['rtime_right_base'] for feat in all_features])
        padding = max(0.2 * (right_base - left_base), 10)           # avoid 0 or too small value
        # time points sampled on track
        list_rtime = np.arange(max(0, left_base-padding), right_base+padding, scan_time)        #.tolist()
        if not list_rtime.any():
            print(all_features)
            return []
        else:
            composite_mass_track = np.zeros(len(list_rtime))
            for feat in all_features:
                # This simulates a peak using Gaussian function
                sigma = 0.33 * max([feat['rtime']-feat['rtime_left_base'],  feat['rtime_right_base']-feat['rtime']])
                peak_height = feat['peak_area'] / np.sqrt(2 * np.pi * sigma**2)
                peak_data = peak_height * np.exp(-(list_rtime - feat['rtime'])**2/(2*sigma**2))
                composite_mass_track += peak_data

            # detection of elution peaks on composite_mass_track
            min_peak_height = parameters['min_peak_height']
            min_prominence_threshold = parameters['min_prominence_threshold']
            min_fwhm = round( 0.5 * parameters['min_timepoints'] )
            min_intensity_threshold = parameters['min_intensity_threshold']
            wlen = parameters['wlen']

            scaling_factor, LOW, HIGH = 1, min_intensity_threshold, 1E8
            scaling_factor = HIGH/max(composite_mass_track)
            composite_mass_track = composite_mass_track * scaling_factor
            composite_mass_track = np.where(composite_mass_track > LOW, composite_mass_track, 0)

            peaks, properties = find_peaks(composite_mass_track,
                                            height=min_peak_height,
                                            distance=min_fwhm,
                                            prominence=min_prominence_threshold,
                                            width=min_fwhm,
                                            wlen=wlen,
                                            )
            new_peaks = []
            for ii in range(peaks.size):
                new_peaks.append({
                        'mz': round(mz, 4),
                        'parent_masstrack_id': id_number,
                        'rtime': round(list_rtime[peaks[ii]], 2),
                        'height': int(properties['peak_heights'][ii] / scaling_factor),
                        'rtime_left_base': round(list_rtime[properties['left_bases'][ii]], 2),
                        'rtime_right_base': round(list_rtime[properties['right_bases'][ii]], 2)
                })

            # map input features to new_peaks
            for peak in new_peaks:
                peak['members'] = get_features_in_range(peak, all_features)

            return new_peaks


def get_features_in_range(peak, all_features):
    """Return features with rtime within the base range of peak."""
    return [feat for feat in all_features if peak['rtime_left_base'] < feat['rtime'] < peak['rtime_right_base']]


def collect_rt_calibration_features(list_features, min_area, min_snr=10):
    """
    Get features that are high and unique per mass track, used to align rtime.
    """
    selected = [f for f in list_features if f['peak_area']>min_area and f['snr']>min_snr]
    parent_masstrack_ids = [f['parent_masstrack_id'] for f in selected]
    track_counts = {x: parent_masstrack_ids.count(x) for x in parent_masstrack_ids}
    selected = [f for f in selected if track_counts[f['parent_masstrack_id']] == 1]
    return selected


def verify_path_list(file_list_paths, filename_contains='Feature_table.tsv'):
    """
    file_list_paths : a text file containing list of paths, one per line,
    each a feature table or an asari project directory.
    For a directory, Feature_table.tsv in it is used
    (or export/full_Feature_table.tsv, as in projects before asari 2.0).
    """
    new = []
    for line in open(file_list_paths).read().splitlines():
        line = line.strip()
        if line:
            if filename_contains in line:
                new.append(line)
            else:
                infile = os.path.join(line, 'Feature_table.tsv')
                if not os.path.exists(infile):
                    legacy = os.path.join(line, 'export', 'full_Feature_table.tsv')
                    if os.path.exists(legacy):
                        infile = legacy
                new.append(infile)
    return new

def read_table_to_json(feature_table_file):
    """
    Parse feature_table_file to internal format with list_features and list_mass_tracks.

    list_mass_tracks: [{ 'id_number': ii,  'mz': xic[0], ...]
    """
    num_samples, L = read_features_from_asari_table(open(feature_table_file).read())
    mass_tracks = {f['parent_masstrack_id']: f['mz'] for f in L}
    list_mass_tracks = [{'id_number': int(k), 'mz': v} for k,v in mass_tracks.items()]
    return {
        'name': feature_table_file,
        'status:eic': 'passed',
        'num_samples': num_samples,
        'number_features': len(L),
        'sample_data': {
            'list_features': L,
            'list_mass_tracks': list_mass_tracks
        }
    }


# -----------------------------------------------------------------------------
# join_algorithm dispatch
# -----------------------------------------------------------------------------

JOIN_ALGORITHMS = ['composite_table_merge', 'annotation_combine']

def join_project(file_list_paths, parameters):
    """
    Run subcommand `join` by parameters['join_algorithm']:
    'composite_table_merge' (default) joins LC/GC-MS feature tables (join_tables);
    'annotation_combine' combines annotations from different runs or methods, by workflow
    (GC: combine_gc_annotations; LC or LCMSMS: combine_lc_annotations).
    """
    algorithm = parameters.get('join_algorithm', 'composite_table_merge')
    if algorithm == 'composite_table_merge':
        join_tables(file_list_paths, parameters)
    elif algorithm == 'annotation_combine':
        if parameters['workflow'] == 'GC':
            combine_gc_annotations(file_list_paths, parameters)
        elif parameters['workflow'] in ('LC', 'LCMSMS'):
            combine_lc_annotations(file_list_paths, parameters)
        else:
            print(f"annotation_combine is not supported for workflow {parameters['workflow']}; GC, LC or LCMSMS.")
    else:
        print(f"Unknown join_algorithm {algorithm}, must be one of {', '.join(JOIN_ALGORITHMS)}.")


# -----------------------------------------------------------------------------
# annotation_combine for GC
# -----------------------------------------------------------------------------

def read_annotation_path_list(file_list_paths, pattern='Annotated_pseudospectra_*.json'):
    """
    file_list_paths : a text file containing list of paths, one per line,
    each an annotation JSON file or an annotation directory (containing pattern).
    """
    import glob
    paths = []
    for line in open(file_list_paths).read().splitlines():
        line = line.strip().strip('",')
        if not line or line in ('[', ']'):
            continue
        if os.path.isdir(line):
            found = sorted(glob.glob(os.path.join(line, pattern)))
            if not found:
                print(f"Warning: no {pattern} in {line}, skipped.")
            paths += found[:1]
        elif os.path.isfile(line):
            paths.append(line)
        else:
            print(f"Warning: {line} not found, skipped.")
    return paths

def _read_tsv_rows(infile, key):
    """Read a tsv into {row[key]: row}; empty dict if file not found."""
    if not os.path.exists(infile):
        return {}
    df = pd.read_csv(infile, sep='\t', dtype=str, keep_default_na=False)
    return {row[key]: row for row in df.to_dict('records')}

def load_gc_annotation_source(json_file, source_index):
    """
    Load an annotation JSON (from GC annotate) and its companion tsv files in the same directory
    (empCpds_<handle>.tsv, Features_<handle>.tsv), used for library and feature details.
    """
    data = json.load(open(json_file))
    handle = os.path.basename(json_file).replace('Annotated_pseudospectra_', '').replace('.json', '')
    indir = os.path.dirname(json_file)
    return {
        'source': source_index,
        'file': json_file,
        'database': data.get('database'),
        'algorithm': data.get('algorithm'),
        'time_stamp': data.get('time_stamp'),
        'list_of_empCpds': data.get('list_of_empCpds', []),
        'empCpd_rows': _read_tsv_rows(os.path.join(indir, 'empCpds_' + handle + '.tsv'), 'empCpd'),
        'feature_rows': _read_tsv_rows(os.path.join(indir, 'Features_' + handle + '.tsv'), 'feature'),
    }

def compound_identity_key(epd, inchikey_block=1):
    """
    Common identifier of the annotated compound across databases:
    first inchikey_block blocks of InChIKey (default 1, the skeleton, ignoring stereochemistry and charge);
    else normalized name; else lib_entry_id.
    """
    inchikey = (epd.get('inchikey') or '').strip()
    if inchikey and inchikey.lower() not in ('none', 'nan'):
        return 'inchikey:' + '-'.join(inchikey.split('-')[:inchikey_block])
    name = (epd.get('name') or '').strip().lower()
    if name:
        return 'name:' + name
    return 'lib:' + str(epd.get('lib_entry_id'))

def _representative_score(epd):
    """Rank to select representative entry: cosine score, then entropy score."""
    return (epd.get('cosine_score') or 0, epd.get('entropy_score') or 0)

def combine_gc_empCpds(sources, inchikey_block=1):
    """
    Combine annotated empCpds from multiple sources (runs or methods on the same feature table).
    Feature IDs are consistent across sources, empCpd IDs are not.

    Entries are grouped if they have the same compound identity (compound_identity_key)
    and share at least one feature (connected components, so A~B and B~C are grouped).
    The same compound on disjoint features (e.g. at different RI) is kept separate.
    Within a group, the entry of highest cosine score is the representative,
    for name, quant_ion, RI, scores and peaks (peaks_as_features, peaks_in_lib for mirror plots).

    Returns list of combined empCpds, sorted by RI, with new IDs, and extra fields:
    identity_key, number_of_sources, sources ([{source, empCpd, lib_entry_id, name, inchikey, scores}]),
    shared_features_with (combined empCpds of different identity sharing features, i.e. conflicting annotations).
    """
    entries = []        # (source_index, epd)
    for S in sources:
        for epd in S['list_of_empCpds']:
            entries.append((S['source'], epd))

    # union-find within each identity key, by shared features
    parent = list(range(len(entries)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen = {}       # (identity_key, feature) -> entry index
    for ii, (_, epd) in enumerate(entries):
        key = compound_identity_key(epd, inchikey_block)
        for f in set(epd.get('features') or [epd.get('quant_ion')]):
            if (key, f) in seen:
                a, b = find(ii), find(seen[(key, f)])
                if a != b:
                    parent[a] = b
            else:
                seen[(key, f)] = ii
    groups = {}
    for ii in range(len(entries)):
        groups.setdefault(find(ii), []).append(ii)

    combined = []
    for members in groups.values():
        rep_source, rep = max((entries[i] for i in members), key=lambda x: _representative_score(x[1]))
        features = list(rep.get('features') or [])
        _fset = set(features)
        for i in members:
            for f in entries[i][1].get('features') or []:
                if f not in _fset:
                    features.append(f)
                    _fset.add(f)
        new = {k: v for k, v in rep.items() if k not in ('id', 'features')}
        new.update({
            'identity_key': compound_identity_key(rep, inchikey_block),
            'features': features,
            'representative': {'source': rep_source, 'empCpd': rep['id']},
            'number_of_sources': len(set(entries[i][0] for i in members)),
            'sources': [{
                'source': entries[i][0],
                'empCpd': entries[i][1].get('id'),
                'lib_entry_id': entries[i][1].get('lib_entry_id'),
                'name': entries[i][1].get('name'),
                'inchikey': entries[i][1].get('inchikey'),
                'cosine_score': entries[i][1].get('cosine_score'),
                'entropy_score': entries[i][1].get('entropy_score'),
            } for i in members],
        })
        combined.append(new)

    combined.sort(key=lambda x: (x.get('RI') is None, x.get('RI') or 0, x['identity_key']))
    for ii, epd in enumerate(combined):
        epd['id'] = 'empCpd_' + "{:04d}".format(ii + 1)

    # conflicting annotations: different identity on shared features
    feature_to_cpds = {}
    for epd in combined:
        for f in epd['features']:
            feature_to_cpds.setdefault(f, []).append(epd['id'])
    for epd in combined:
        others = set()
        for f in epd['features']:
            others.update(feature_to_cpds[f])
        others.discard(epd['id'])
        epd['shared_features_with'] = sorted(others)
    return combined

def write_tsv_combined_empCpds(combined, sources, outfile):
    """
    Table of combined empCpds, same columns as ``empCpds_`` tsv from GC annotate
    (library and quant ion details from the representative source), plus combination columns.
    """
    base_header = ['empCpd', 'name', 'inchikey', 'formula', 'lib_exact_mass', 'score_cosine', 'score_entropy',
        'num_matched_peaks', 'quant_ion', 'mz_quant_ion', 'rtime_quant_ion(min)', 'RI_quant_ion', 'lib_RI', 'delta_RI',
        'has_basepeak', 'has_molecular_ion', 'presence_ratio']
    extra_header = ['number_of_sources', 'sources', 'number_of_features', 'number_shared_features_with', 'database']
    dict_sources = {S['source']: S for S in sources}
    s = '\t'.join(base_header + extra_header) + '\n'
    for epd in combined:
        rep = epd['representative']
        row = dict(dict_sources[rep['source']]['empCpd_rows'].get(rep['empCpd'], {}))
        row.update({'empCpd': epd['id'], 'name': epd.get('name'), 'inchikey': epd.get('inchikey'),
                    'quant_ion': epd.get('quant_ion')})
        for k, v in [('score_cosine', 'cosine_score'), ('score_entropy', 'entropy_score')]:
            if not row.get(k) and epd.get(v) is not None:
                row[k] = round(epd[v], 3)
        if not row.get('RI_quant_ion') and epd.get('RI') is not None:
            row['RI_quant_ion'] = int(epd['RI'])
        if not row.get('mz_quant_ion') and epd.get('quant_ion_mz') is not None:
            row['mz_quant_ion'] = epd['quant_ion_mz']
        row.update({
            'number_of_sources': epd['number_of_sources'],
            'sources': ';'.join(f"{x['source']}:{x['empCpd']}" for x in epd['sources']),
            'number_of_features': len(epd['features']),
            'number_shared_features_with': len(epd['shared_features_with']),
            'database': dict_sources[rep['source']]['database'],
        })
        s += '\t'.join(['' if row.get(k) is None else str(row.get(k)) for k in base_header + extra_header]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)

def write_tsv_combined_features(combined, sources, outfile):
    """
    Table per annotated feature, same columns as ``Features_`` tsv from GC annotate.
    Feature details (mz, rtime, RI, peak_area etc.) are taken from any source that has them,
    since pseudospectra-based annotations may lack feature details.
    """
    header = ['feature', 'mz', 'rtime', 'RI', 'empCpd', 'quant_ion', 'score_cosine', 'score_entropy',
              'name', 'inchikey', 'lib_RI', 'delta_RI', 'peak_area', 'cSelectivity', 'peak_shape', 'snr',
              'detection_counts', 'number_of_sources']
    feature_cols = ['mz', 'rtime', 'RI', 'peak_area', 'cSelectivity', 'peak_shape', 'snr', 'detection_counts']
    dict_features = {}
    for S in sources:
        for f, row in S['feature_rows'].items():
            FF = dict_features.setdefault(f, {})
            for k in feature_cols:
                if not FF.get(k) and row.get(k):
                    FF[k] = row[k]
    dict_sources = {S['source']: S for S in sources}
    s = '\t'.join(header) + '\n'
    for epd in combined:
        rep = epd['representative']
        lib_RI = dict_sources[rep['source']]['empCpd_rows'].get(rep['empCpd'], {}).get('lib_RI', '')
        for f in epd['features']:
            FF = dict_features.get(f, {})
            try:
                delta_RI = round(float(FF['RI']) - float(lib_RI), 2)
            except (KeyError, ValueError, TypeError):
                delta_RI = ''
            row = dict(FF, feature=f, empCpd=epd['id'], quant_ion=epd.get('quant_ion'),
                       score_cosine='' if epd.get('cosine_score') is None else round(epd['cosine_score'], 3),
                       score_entropy='' if epd.get('entropy_score') is None else round(epd['entropy_score'], 3),
                       name=epd.get('name'), inchikey=epd.get('inchikey'), lib_RI=lib_RI, delta_RI=delta_RI,
                       number_of_sources=epd['number_of_sources'])
            s += '\t'.join(['' if row.get(k) is None else str(row.get(k)) for k in header]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)

def combine_gc_annotations(file_list_paths, parameters):
    """
    join --join_algorithm annotation_combine --workflow GC.
    Combine multiple GC annotation JSON files (Annotated_pseudospectra_*.json, from runs or methods
    on the same feature table) to one JSON, and export ``empCpds_`` and ``Features_`` tsv files and MSP.
    See combine_gc_empCpds for the combination rules.

    file_list_paths : a text file containing list of paths to annotation JSON files or directories, one per line.
    """
    from .utils import NpEncoder
    from .gc_functions import json_pseudospectra_to_msp

    json_files = read_annotation_path_list(file_list_paths)
    if len(json_files) < 2:
        print("Need at least two annotation files to combine, found %d." %len(json_files))
        return
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    handle = parameters['project_name']
    outdir = os.path.join(parameters['outdir'], '_'.join([handle, 'annotation_combine', ''.join(time_stamp)]))
    os.makedirs(outdir, exist_ok=True)
    print(f"Output directory: {outdir}.\n")

    sources = [load_gc_annotation_source(f, ii) for ii, f in enumerate(json_files)]
    for S in sources:
        print("  source %d: %d empCpds, %s, %s" %(S['source'], len(S['list_of_empCpds']), S['algorithm'], S['file']))
    combined = combine_gc_empCpds(sources, inchikey_block=parameters.get('combine_inchikey_block', 1))

    write_tsv_combined_empCpds(combined, sources, os.path.join(outdir, 'empCpds_' + handle + '.tsv'))
    write_tsv_combined_features(combined, sources, os.path.join(outdir, 'Features_' + handle + '.tsv'))
    json_output_file = os.path.join(outdir, 'Annotated_pseudospectra_' + handle + '.json')
    json_pseudospectra_to_msp(combined, json_output_file.replace('.json', '.msp'))
    json_export = {
        'type': 'GCMS',
        'database': '; '.join(sorted(set(str(S['database']) for S in sources))),
        'algorithm': 'annotation_combine',
        'time_stamp': ':'.join(time_stamp),
        'sources': [{k: S[k] for k in ['source', 'file', 'database', 'algorithm', 'time_stamp']}
                    | {'number_of_empCpds': len(S['list_of_empCpds'])} for S in sources],
        'list_of_empCpds': combined,
    }
    with open(json_output_file, 'w', encoding='utf-8') as O:
        json.dump(json_export, O, ensure_ascii=False, cls=NpEncoder)
    parameters['join_sources'] = json_files
    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, ensure_ascii=False, indent=2, cls=NpEncoder)

    multi = sum(1 for x in combined if x['number_of_sources'] > 1)
    print("\nCombined %d empCpds into %d (%d supported by multiple sources, %d with features shared by other annotations)."
          %(sum(len(S['list_of_empCpds']) for S in sources), len(combined), multi,
            sum(1 for x in combined if x['shared_features_with'])))
    print(f"Results written to {outdir}.\n")


# -----------------------------------------------------------------------------
# annotation_combine for LC (LC-MS, LC-MS/MS)
# -----------------------------------------------------------------------------
#
# Sources are result directories, identified by project.json ('run', 'workflow'):
#   process (LC or LCMSMS)    : Feature_table.tsv, ms2_spectra.json (reference of all others)
#   deconvolute, LC           : <project>_list_khipus.json (preannotation, empCpds)
#   annotate, LC              : Annotated_features_<project>.json (MS1, mz or rtmz)
#   annotate, LCMSMS or MSMS  : <project>_MS2_annotation.json (MS/MS spectral search)
# All must come from the same process run (same feature table).

LC_SOURCE_FILES = {
    ('deconvolute', 'LC'): '*_list_khipus.json',
    ('annotate', 'LC'): 'Annotated_features_*.json',
    ('annotate', 'LCMSMS'): '*_MS2_annotation.json',
    ('annotate', 'MSMS'): '*_MS2_annotation.json',
}

def read_dir_path_list(file_list_paths):
    """
    file_list_paths : a text file containing list of result directories (or files in them), one per line.
    Returns list of directories, in input order, without duplicates.
    """
    dirs = []
    for line in open(file_list_paths).read().splitlines():
        line = line.strip().strip('",')
        if not line or line in ('[', ']') or line.startswith('#'):
            continue
        d = line if os.path.isdir(line) else os.path.dirname(line)
        if not os.path.isdir(d):
            print(f"Warning: {line} not found, skipped.")
        elif os.path.realpath(d) not in [os.path.realpath(x) for x in dirs]:
            dirs.append(d)
    return dirs

def _read_project_json(d):
    f = os.path.join(d, 'project.json')
    if not os.path.exists(f):
        return None
    with open(f, encoding='utf-8') as F:
        return json.load(F)

def _input_process_dir(desc):
    """The process directory of a run's input (a process directory, or a file in it)."""
    p = (desc.get('input') or '').rstrip('/')
    if os.path.isfile(p) or os.path.splitext(p)[1] in ('.tsv', '.json', '.msp', '.txt'):
        p = os.path.dirname(p)
    if os.path.basename(p) == 'export':         # projects before asari 2.0
        p = os.path.dirname(p)
    return p

def _same_dir(recorded, actual):
    """True if recorded path is actual directory; by name if recorded path no longer exists (moved)."""
    if os.path.exists(recorded):
        return os.path.realpath(recorded) == os.path.realpath(actual)
    return os.path.basename(recorded) == os.path.basename(os.path.realpath(actual))

def _find_feature_table(d, desc=None):
    for name in [(desc or {}).get('output_feature_table', 'Feature_table.tsv'), 'Feature_table.tsv',
                 os.path.join('export', 'full_Feature_table.tsv')]:
        f = os.path.join(d, name)
        if os.path.exists(f):
            return f
    return None

def collect_lc_sources(dirs, parameters):
    """
    Classify result directories by project.json, and verify they came from the same process run.

    The reference is the process directory in dirs, or else the common input of the other runs.
    A source is excluded (with a warning) if its input is not the reference, its mode differs,
    or its content does not match the reference feature table.

    Returns (reference, sources), reference = {'dir', 'desc', 'feature_table'},
    sources = [{'source', 'dir', 'run', 'workflow', 'file', 'desc'}, ...]; or (None, []) if failed.
    """
    import glob
    process_dirs, candidates = [], []
    for d in dirs:
        desc = _read_project_json(d)
        if desc is None:
            print(f"Warning: no project.json in {d}, skipped.")
            continue
        run, workflow = desc.get('run'), desc.get('workflow')
        if run == 'process':
            process_dirs.append((d, desc))
        elif (run, workflow) in LC_SOURCE_FILES:
            files = sorted(glob.glob(os.path.join(d, LC_SOURCE_FILES[(run, workflow)])))
            if files:
                candidates.append({'dir': d, 'run': run, 'workflow': workflow, 'file': files[0], 'desc': desc})
            else:
                print(f"Warning: no {LC_SOURCE_FILES[(run, workflow)]} in {d}, skipped.")
        else:
            print(f"Warning: {d} (run {run}, workflow {workflow}) is not used in LC annotation_combine, skipped.")

    if len(process_dirs) > 1:
        print("Error: more than one process directory given: " + ', '.join(x[0] for x in process_dirs))
        return None, []
    if process_dirs:
        ref_dir, ref_desc = process_dirs[0]
    else:
        inputs = sorted(set(_input_process_dir(S['desc']) for S in candidates))
        if len(inputs) != 1:
            print("Error: no process directory given, and the sources have different inputs: " + ', '.join(inputs))
            return None, []
        ref_dir, ref_desc = inputs[0], _read_project_json(inputs[0]) or {}
        print(f"Using process directory from the input of sources: {ref_dir}")
    feature_table = _find_feature_table(ref_dir, ref_desc)
    if not feature_table:
        print(f"Error: no feature table found in {ref_dir}.")
        return None, []
    reference = {'dir': ref_dir, 'desc': ref_desc, 'feature_table': feature_table}

    mode = ref_desc.get('mode') or parameters.get('mode')
    sources, labels = [], set()
    for S in candidates:
        if not _same_dir(_input_process_dir(S['desc']), ref_dir):
            print(f"Warning: {S['dir']} is from input {S['desc'].get('input')}, not {ref_dir}; excluded.")
            continue
        if mode and S['desc'].get('mode') and S['desc']['mode'] != mode:
            print(f"Warning: {S['dir']} is in mode {S['desc']['mode']}, not {mode}; excluded.")
            continue
        label = S['desc'].get('project_name') or os.path.basename(S['dir'])
        while label in labels:
            label += '_' + str(len(sources))
        labels.add(label)
        S['source'] = label
        sources.append(S)
    return reference, sources

def _valid_str(x):
    return bool(x) and str(x).strip().lower() not in ('none', 'nan', 'null')

def compound_keys(name='', inchikey='', inchikey_block=1):
    """Keys to compare compounds across databases: InChIKey (first inchikey_block blocks) and lower case name."""
    keys = set()
    if _valid_str(inchikey):
        keys.add('inchikey:' + '-'.join(str(inchikey).strip().split('-')[:inchikey_block]))
    if _valid_str(name):
        keys.add('name:' + str(name).strip().lower())
    return keys

def ms2_entry_fields(db_entry):
    """name, inchikey, precursor_type, formula, lib_id and precursor m/z from an MS/MS library entry (MSP fields, any case)."""
    low = {str(k).lower(): v for k, v in db_entry.items()}
    name = str(low.get('name') or '').split(';')[0].strip()       # MS-DIAL style 'Choline; CE30; InChIKey'
    return {
        'name': name,
        'inchikey': low.get('inchikey') or '',
        'ion': low.get('precursor_type') or low.get('precursortype') or low.get('adduct') or '',
        'formula': low.get('formula') or '',
        'lib_id': low.get('db#') or low.get('accession') or low.get('id') or '',
        'lib_mz': low.get('precursor_mz') or low.get('precursormz') or '',
    }

def map_ms2_to_features(list_spectra, list_features, mz_tolerance_ppm=5, rt_tolerance=10):
    """
    Map MS/MS spectra to MS1 features, by parent_masstrack_id (from asari LCMSMS process) if available,
    else precursor m/z within mz_tolerance_ppm. Among candidates, a feature whose peak range
    (rtime_left_base, rtime_right_base) contains the MS/MS rtime is preferred, then closest rtime
    within rt_tolerance.

    list_spectra : [{'id', 'precursor_mz', 'rtime', 'parent_masstrack_id' (optional)}, ...]
    Returns {spectrum id: feature id}.
    """
    by_track = {}
    for F in list_features:
        by_track.setdefault(str(F.get('parent_masstrack_id')), []).append(F)
    sorted_features = sorted(list_features, key=lambda x: x['mz'])
    mzs = np.array([F['mz'] for F in sorted_features])
    mapping = {}
    for sp in list_spectra:
        mz, rt = sp.get('precursor_mz'), sp.get('rtime')
        if rt is None:
            continue
        track = sp.get('parent_masstrack_id')
        if track is not None and str(track) in by_track:
            cands = by_track[str(track)]
        elif mz:
            tol = mz * mz_tolerance_ppm * 1e-6
            a, b = np.searchsorted(mzs, [mz - tol, mz + tol])
            cands = sorted_features[a: b + 1]
            cands = [F for F in cands if abs(F['mz'] - mz) <= tol]
        else:
            continue
        inside = [F for F in cands if F.get('rtime_left_base', F['rtime']) <= rt <= F.get('rtime_right_base', F['rtime'])]
        near = inside or [F for F in cands if abs(F['rtime'] - rt) <= rt_tolerance]
        if near:
            mapping[sp['id']] = min(near, key=lambda F: abs(F['rtime'] - rt))['id']
    return mapping

def combine_lc_annotations(file_list_paths, parameters):
    """
    join --join_algorithm annotation_combine --workflow LC (or LCMSMS).
    Combine LC-MS annotation results of the same process run:
    MS1 feature annotation (annotate LC, mz or rtmz, one or more databases),
    MS/MS annotation (annotate LCMSMS/MSMS, one or more databases),
    and preannotation (deconvolute LC, khipu empCpds).

    - MS1 direct matches from all sources are combined; with preannotation, matches are combined within
      empCpds by combine_preannotation (neutral mass consistency, then M0, then common ions).
    - MS/MS spectra are mapped to MS1 features by parent mass track (or precursor m/z) and rtime
      (map_ms2_to_features); their annotations are added to the feature, and to features in the same empCpd.
    - number_of_sources: sources supporting the same compound (InChIKey first block or name) on the feature.

    file_list_paths : a text file listing result directories, one per line (process, deconvolute, annotate).
    Sources are identified and verified by project.json (collect_lc_sources).

    Output, in <project>_annotation_combine_<time stamp>/:
    - <project>_MS2_annotation.json and .tsv, all MS/MS annotations with mapped features (viewable by `asari viz`);
    - Annotated_features_<project>.json and Features_<project>.tsv, MS1 features with MS1 and MS/MS annotations.
    """
    from .utils import NpEncoder
    from .annotate import combine_preannotation

    dirs = read_dir_path_list(file_list_paths)
    reference, sources = collect_lc_sources(dirs, parameters)
    if reference is None:
        return
    mode = reference['desc'].get('mode') or parameters['mode']
    ppm = parameters['mz_tolerance_ppm']
    inchikey_block = parameters.get('combine_inchikey_block', 1)

    with open(reference['feature_table']) as f:
        _n, list_features = read_features_from_asari_table(f.read())
    dict_features = {F['id']: F for F in list_features}
    print(f"Reference: {len(list_features)} features in {reference['feature_table']}.")

    # verify and load sources
    ms1_sources, ms2_sources, empCpd_sources = [], [], []
    for S in sources:
        with open(S['file'], encoding='utf-8') as f:
            data = json.load(f)
        if S['run'] == 'deconvolute':
            ids = [F['id'] for emp in data for F in emp.get('MS1_pseudo_Spectra', [])]
            missing = [x for x in ids if x not in dict_features]
            if missing:
                print(f"Warning: {len(missing)} features in {S['file']} not in reference feature table; excluded.")
                continue
            S['data'] = data
            empCpd_sources.append(S)
        elif S['workflow'] == 'LC':
            n = data.get('number_of_features')
            missing = [F['id'] for F in data.get('list_of_features', []) if F['id'] not in dict_features]
            if (n is not None and n != len(list_features)) or missing:
                print(f"Warning: {S['file']} was on {n} features, not matching the reference; excluded.")
                continue
            S['data'] = data
            ms1_sources.append(S)
        else:
            S['data'] = data
            ms2_sources.append(S)
        S['database'] = data.get('database', '') if isinstance(data, dict) else ''
        S['algorithm'] = data.get('algorithm', '') if isinstance(data, dict) else 'khipu'
        S['number_of_entries'] = len(data.get('list_of_features', [])) if isinstance(data, dict) else len(data)
        print(f"  source {S['source']}: {S['run']} {S['workflow']} {S['algorithm']}, "
              f"{S['number_of_entries']} entries, {S['file']}")
    if len(empCpd_sources) > 1:
        print(f"Warning: multiple preannotation sources; using {empCpd_sources[0]['file']}.")
    if not ms1_sources and not ms2_sources:
        print("No annotation source to combine.")
        return

    # MS1: direct matches from all sources, ids prefixed by source
    direct, dict_lib_entries = {}, {}
    for S in ms1_sources:
        for F in S['data'].get('list_of_features', []):
            for m in F.get('matches', []):
                if m.get('match', 'direct') != 'direct':        # recomputed from preannotation below
                    continue
                lib_id = S['source'] + ':' + str(m['id'])
                E = {k: v for k, v in m.items() if k not in ('match', 'ppm', 'delta_rtime')}
                E.update({'lib_id': m['id'], 'source': S['source'], 'database': S['database']})
                dict_lib_entries[lib_id] = E
                direct.setdefault(F['id'], [])
                if lib_id not in direct[F['id']]:
                    direct[F['id']].append(lib_id)
    combined_ms1, dict_feature_empCpd, empCpd_members = direct, {}, {}
    if empCpd_sources:
        list_empCpds = empCpd_sources[0]['data']
        combined_ms1, dict_feature_empCpd = combine_preannotation(
            direct, list_empCpds, dict_lib_entries, mode=mode, mz_tolerance_ppm=ppm)
        for emp in list_empCpds:
            empCpd_members[emp.get('interim_id')] = [F['id'] for F in emp.get('MS1_pseudo_Spectra', [])]

    feature_items = {}
    for fid, lib_ids in combined_ms1.items():
        F = dict_features[fid]
        for lib_id in lib_ids:
            E = dict_lib_entries[lib_id]
            item = {'evidence': 'MS1', 'match': 'direct' if lib_id in direct.get(fid, []) else 'empCpd'}
            item.update({k: E.get(k, '') for k in ('source', 'database', 'lib_id', 'name', 'identifier', 'ion',
                                                     'isotope', 'formula', 'neutral_mass', 'inchikey')})
            item['lib_mz'], item['lib_rtime'] = E.get('mz'), E.get('rtime')
            if item['match'] == 'direct':
                item['ppm'] = round((F['mz'] - E['mz']) / E['mz'] * 1e6, 2)
                if E.get('rtime') not in (None, ''):
                    item['delta_rtime'] = round(F['rtime'] - float(E['rtime']), 2)
            feature_items.setdefault(fid, []).append(item)

    # MS/MS: map spectra to features
    ms2_spectra = {}
    ms2_file = os.path.join(reference['dir'], 'ms2_spectra.json')
    if os.path.exists(ms2_file):
        with open(ms2_file, encoding='utf-8') as f:
            ms2_spectra = {sp['id']: sp for sp in json.load(f) if 'id' in sp}
    list_ms2_records = []
    for S in ms2_sources:
        for r in S['data'].get('list_of_features', []):
            rec = dict(r)
            rec.update({'source': S['source'], 'database': S['database']})
            rec.update({k: v for k, v in ms2_entry_fields(r.get('database_entry', {})).items()
                        if k in ('name', 'inchikey', 'ion', 'formula')})
            list_ms2_records.append(rec)
    query = {}
    for rec in list_ms2_records:
        sp = ms2_spectra.get(rec['spid'], {})
        query[rec['spid']] = {'id': rec['spid'], 'precursor_mz': sp.get('precursor_mz', rec.get('mz')),
                              'rtime': sp.get('rtime', rec.get('rtime')),
                              'parent_masstrack_id': sp.get('parent_masstrack_id')}
    ms2_to_feature = map_ms2_to_features(list(query.values()), list_features, ppm,
                                         parameters.get('lc_rt_tolerance', 10))
    print(f"MS/MS: {len(list_ms2_records)} annotations on {len(query)} spectra, "
          f"{sum(1 for x in query if x in ms2_to_feature)} spectra mapped to MS1 features.")

    for ii, rec in enumerate(list_ms2_records):
        fid = ms2_to_feature.get(rec['spid'])
        rec['feature'] = fid or ''
        rec['empCpd'] = dict_feature_empCpd.get(fid, {}).get('interim_id', '') if fid else ''
        if not fid:
            continue
        F = dict_features[fid]
        fields = ms2_entry_fields(rec.get('database_entry', {}))
        base = {'evidence': 'MS2', 'source': rec['source'], 'database': rec['database'],
                'lib_id': fields['lib_id'] or rec.get('database_key', ''), 'name': rec['name'],
                'identifier': rec.get('database_key', ''), 'ion': rec['ion'], 'formula': rec['formula'],
                'inchikey': rec['inchikey'], 'lib_mz': fields['lib_mz'], 'spid': rec['spid'],
                'cosine_score': rec.get('cosine_score'), 'entropy_score': rec.get('entropy_score'),
                'number_matched_peaks': rec.get('number_matched_peaks'), '_ms2_index': ii}
        targets = empCpd_members.get(rec['empCpd'], [fid]) if rec['empCpd'] else [fid]
        for t in targets:
            item = dict(base, match='direct' if t == fid else 'empCpd')
            if t == fid:
                item['ppm'] = round((F['mz'] - rec['mz']) / rec['mz'] * 1e6, 2) if rec.get('mz') else ''
                item['delta_rtime'] = round(F['rtime'] - rec['rtime'], 2) if rec.get('rtime') is not None else ''
            feature_items.setdefault(t, []).append(item)

    # number of sources supporting the same compound, per feature
    for fid, items in feature_items.items():
        keys = [compound_keys(x.get('name'), x.get('inchikey'), inchikey_block) for x in items]
        for x, k in zip(items, keys):
            x['number_of_sources'] = len(set(y['source'] for y, k2 in zip(items, keys) if k & k2) | {x['source']})
            if x['evidence'] == 'MS2' and x['match'] == 'direct':
                list_ms2_records[x['_ms2_index']]['number_of_sources'] = x['number_of_sources']
        items.sort(key=lambda x: (x['match'] != 'direct', -x['number_of_sources'], x['evidence'],
                                  -(x.get('cosine_score') or 0), abs(x.get('ppm') or 0) if x.get('ppm') != '' else 0))
        for x in items:
            x.pop('_ms2_index', None)
    for rec in list_ms2_records:        # unmapped spectra: count among MS/MS annotations of the same spectrum
        if 'number_of_sources' not in rec:
            k = compound_keys(rec['name'], rec['inchikey'], inchikey_block)
            rec['number_of_sources'] = len(set(r['source'] for r in list_ms2_records if r['spid'] == rec['spid']
                                    and k & compound_keys(r['name'], r['inchikey'], inchikey_block)) | {rec['source']})

    # export
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    handle = parameters['project_name']
    outdir = os.path.join(parameters['outdir'], '_'.join([handle, 'annotation_combine', ''.join(time_stamp)]))
    os.makedirs(outdir, exist_ok=True)
    meta_sources = [{k: S[k] for k in ('source', 'dir', 'run', 'workflow', 'file', 'database', 'algorithm',
                                        'number_of_entries')} for S in ms1_sources + ms2_sources + empCpd_sources[:1]]
    common = {
        'algorithm': 'annotation_combine',
        'mode': mode,
        'mz_tolerance_ppm': ppm,
        'lc_rt_tolerance': parameters.get('lc_rt_tolerance', 10),
        'time_stamp': ':'.join(time_stamp),
        'process_dir': reference['dir'],
        'feature_table': reference['feature_table'],
        'sources': meta_sources,
    }
    empCpd_fields = ['interim_id', 'neutral_formula_mass', 'isotope', 'modification']

    # MS/MS annotations
    ms2_keys = ['spid', 'mz', 'rtime', 'feature', 'empCpd', 'source', 'database', 'entropy_score', 'cosine_score',
                'number_matched_peaks', 'number_of_sources', 'database_key', 'name', 'inchikey', 'ion', 'formula']
    list_ms2_records.sort(key=lambda r: (r['spid'], -r['number_of_sources'], -(r.get('cosine_score') or 0)))
    if list_ms2_records:
        with open(os.path.join(outdir, handle + '_MS2_annotation.tsv'), 'w', encoding='utf-8') as O:
            O.write('\t'.join(ms2_keys) + '\n')
            for r in list_ms2_records:
                O.write('\t'.join(['' if r.get(k) is None else str(r.get(k)) for k in ms2_keys]) + '\n')
        with open(os.path.join(outdir, handle + '_MS2_annotation.json'), 'w', encoding='utf-8') as O:
            json.dump(dict({'type': 'MS/MS', 'database': '; '.join(S['database'] for S in ms2_sources)},
                           **common, list_of_features=list_ms2_records),
                      O, indent=2, ensure_ascii=False, cls=NpEncoder)

    # MS1 features
    list_annotated = []
    for F in list_features:
        if F['id'] in feature_items:
            F2 = dict(F, matches=feature_items[F['id']])
            if F['id'] in dict_feature_empCpd:
                F2['empCpd'] = dict_feature_empCpd[F['id']]
            list_annotated.append(F2)
    item_fields = ['evidence', 'match', 'source', 'number_of_sources', 'lib_id', 'name', 'identifier', 'ion',
                   'isotope', 'formula', 'inchikey', 'lib_mz', 'ppm', 'lib_rtime', 'delta_rtime', 'neutral_mass',
                   'spid', 'cosine_score', 'entropy_score', 'number_matched_peaks']
    header = ['feature', 'mz', 'rtime', 'empCpd', 'empCpd_neutral_mass', 'feature_isotope', 'feature_ion',
              'number_of_matches'] + item_fields
    with open(os.path.join(outdir, 'Features_' + handle + '.tsv'), 'w', encoding='utf-8') as O:
        O.write('\t'.join(header) + '\n')
        for F in list_annotated:
            emp = F.get('empCpd', {})
            for x in F['matches']:
                row = [F['id'], F['mz'], F['rtime']] + [emp.get(k, '') for k in empCpd_fields] + [
                       len(F['matches'])] + [x.get(k, '') for k in item_fields]
                O.write('\t'.join(['' if v is None else str(v) for v in row]) + '\n')
    with open(os.path.join(outdir, 'Annotated_features_' + handle + '.json'), 'w', encoding='utf-8') as O:
        json.dump(dict({'type': 'LCMS', 'database': '; '.join(S['database'] for S in ms1_sources + ms2_sources)},
                       **common, number_of_features=len(list_features),
                       number_of_annotated_features=len(list_annotated),
                       preannotation=empCpd_sources[0]['file'] if empCpd_sources else '',
                       list_of_features=list_annotated),
                  O, indent=2, ensure_ascii=False, cls=NpEncoder)

    parameters['join_sources'] = [S['dir'] for S in sources]
    parameters['process_dir'] = reference['dir']
    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, ensure_ascii=False, indent=2, cls=NpEncoder)

    n_ms1 = sum(1 for F in list_annotated if any(x['evidence'] == 'MS1' for x in F['matches']))
    n_ms2 = sum(1 for F in list_annotated if any(x['evidence'] == 'MS2' for x in F['matches']))
    n_both = sum(1 for F in list_annotated if any(x['number_of_sources'] > 1 for x in F['matches']))
    print(f"\n{len(list_annotated)} of {len(list_features)} features annotated: {n_ms1} by MS1, {n_ms2} by MS/MS, "
          f"{n_both} with a compound supported by multiple sources.")
    print(f"{len(list_ms2_records)} MS/MS annotations, {sum(1 for r in list_ms2_records if r['feature'])} mapped to features.")
    print(f"Results written to {outdir}.\n")
