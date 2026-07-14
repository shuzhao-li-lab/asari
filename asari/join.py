'''
join as 1st class module, as it will be used to divide & conquer large studies.

The mass alignment is same as landmark_guided_mapping; 
RT alignment is as default LOESS.
Merge feature tables and zero fill absent features in samples.
Track input feature IDs and parent_masstrack_id; keep mapping record.

Example use:
python3 -m asari.main join -i /Users/lish/li.play/join_MT_files.txt -o /Users/lish/li.play/ -j test_join_mt01
'''

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
                             )
from .experiment import ext_Experiment
from .constructors import (CompositeMap, MassGrid)
from .tools.file_io import read_features_from_asari_table


def join_tables(file_list_paths, parameters):
    '''
    This joins data from multiple feature tables using the same LC/GC-MS method.
    Batch correction should be performed after joining. 

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
    triggered by correction_tolerance_ppm (1 default)

    Format of composite feature is like:
    {'id': 'F00000501m', 'mz': 284.2945, 'rtime': 1116.97, 'parent_masstrack_id': 107, 'height': 150503047, 'rtime_left_base': 1115.77, 'rtime_right_base': 1118.17, 
    'members': {'3A5_H_full_Feature_table.tsv': [{'id': 'F388', 'mz': 284.2947, 'apex': 1116.94}], 
    '2C19_H_full_Feature_table.tsv': [{'id': 'F366', 'mz': 284.2946, 'apex': 1116.55}], 
    '2B6_H_full_Feature_table.tsv': [{'id': 'F360', 'mz': 284.2945, 'apex': 1116.81}], 
    'WT_H_full_Feature_table.tsv': []}}

    '''
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
            if idx in table['_mz_landmarks_']:
                REF_landmarks.append(_N1 + ii)

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
    master_track2features = {}
    for table in tables:
        _d = {'': []}           # catch NA
        for feat in  table['sample_data']['list_features']:
            if feat['parent_masstrack_id'] in _d:
                _d[feat['parent_masstrack_id']].append(feat)
            else:
                _d[feat['parent_masstrack_id']] = [feat]
        master_track2features[table['name']] = _d

    print("\nCombining features and detecting elution peaks on the composite mass tracks.\n\n")
    # Update feature list; decide feature correspondence on each composite mass track
    new_feature_list = []
    for ii in MG.MassGrid.index:
        row_dict = MG.MassGrid.iloc[ii, :].to_dict()        # mapping name to mass track index
        row_features = {}
        for t, track2feats in master_track2features.items():
            if row_dict[t]:
                row_features[t] = track2feats[row_dict[t]]
            else:
                row_features[t] = []
        # if str(ii).endswith('1000'): print(row_features)
        new_feature_list += consolidate_track_features(row_features, parameters, ii, MG.MassGrid['mz'][ii])

    # write mapping record 
    names = [table['name'] for table in tables]
    new_feature_list, empty_features = reformat_composite_feature_list(new_feature_list, names)
    # A few features would fail new peak matching, stored in empty_features

    header = "id,mz,rtime,rtime_left_base,rtime_right_base,parent_masstrack_id,height_combined".split(",")
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
    data = np.zeros(len(list_str_intensities[0])).astype(int)
    for L in list_str_intensities:
        data += np.array(L).astype(int)
    return [str(x) for x in data]

def extract_sample_data_dict(filepath, data_start_col=11):
    '''
    filepath : Table in Asari format, which has no empty field.
    return data in str, as conversion is only needed occassionally.
    '''
    d = {}
    for row in open(filepath).readlines():
        a = row.rstrip().split('\t')
        d[a[0]] = a[data_start_col: ]
    return d


def reformat_composite_feature_list(new_feature_list, table_names):
    '''
    Example composite feature:
    {'mz': 284.2945, 'parent_masstrack_id': 107, 'rtime': 1116.97, 'height': 150503047, 'rtime_left_base': 1115.77, 'rtime_right_base': 1118.17, 
    'members': [{'id_number': 'F388', 'id': 'F388', 'mz': 284.2947, 'rtime': 1116.94, 'apex': 1116.94, 'rtime_left_base': 1115.37, 'rtime_right_base': 1118.52, 'parent_masstrack_id': 139, 'peak_area': 75869861.0, 'cSelectivity': 0.64, 'goodness_fitting': 0.94, 'snr': 3.0, 'detection_counts': 9, 'source': '/Users/sli/play/join_test/DnCl_pellets_3A5_H_full_Feature_table.tsv'}, 
    {'id_number': 'F366', 'id': 'F366', 'mz': 284.2946, 'rtime': 1116.8590394349912, 'apex': 1116.55, 'rtime_left_base': 1114.97, 'rtime_right_base': 1118.12, 'parent_masstrack_id': 128, 'peak_area': 72230862.0, 'cSelectivity': 0.62, 'goodness_fitting': 0.95, 'snr': 3.0, 'detection_counts': 9, 'source': '/Users/sli/play/join_test/DnCl_pellets_2C19_H_full_Feature_table.tsv'}, 
    {'id_number': 'F360', 'id': 'F360', 'mz': 284.2945, 'rtime': 1116.9675010942542, 'apex': 1116.81, 'rtime_left_base': 1115.23, 'rtime_right_base': 1118.39, 'parent_masstrack_id': 127, 'peak_area': 67566816.0, 'cSelectivity': 0.54, 'goodness_fitting': 0.92, 'snr': 3.0, 'detection_counts': 9, 'source': '/Users/sli/play/join_test/DnCl_pellets_2B6_H_full_Feature_table.tsv'}]}
    '''
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
    '''
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
    '''
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
    return [feat for feat in all_features if peak['rtime_left_base'] < feat['rtime'] < peak['rtime_right_base']]


def collect_rt_calibration_features(list_features, min_area, min_snr=10):
    '''
    Get features that are high and unique per mass track, used to align rtime.
    '''
    selected = [f for f in list_features if f['peak_area']>min_area and f['snr']>min_snr]
    parent_masstrack_ids = [f['parent_masstrack_id'] for f in selected]
    track_counts = {x: parent_masstrack_ids.count(x) for x in parent_masstrack_ids}
    selected = [f for f in selected if track_counts[f['parent_masstrack_id']] == 1]
    return selected


def verify_path_list(file_list_paths, filename_contains='Feature_table.tsv', dir_contains='export/'):
    '''
    file_list_paths : a text file containing list of paths, one per line
    If not specified, assume full_Feature_table.tsv from Asari is to be used.
    '''
    new = []
    for line in open(file_list_paths).read().splitlines():
        if line.strip():
            if filename_contains in line:
                new.append(line)
            elif dir_contains in line:
                new.append(os.path.join(line, 'full_Feature_table.tsv'))
            else:
                new.append(os.path.join(line, 'export', 'full_Feature_table.tsv'))
    return new

def read_table_to_json(feature_table_file):
    '''
    Parse feature_table_file to internal format with list_features and list_mass_tracks.

    list_mass_tracks: [{ 'id_number': ii,  'mz': xic[0], ...]
    '''
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
