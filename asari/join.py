'''
join as 1st class module, as it will be used to divide & conquer large studies.

The mass alignment is same as landmark_guided_mapping; 
RT alignment is default LOESS plus a similar approach in CSM, rank prioritized.
Merge feature tables and zero fill absent features in samples.

Track input feature IDs and parent_masstrack_id; keep mapping record.
       
'''

import os
import numpy as np
import pandas as pd

from mass2chem.search import find_mzdiff_pairs_from_masstracks

from .mass_functions import (flatten_tuplelist,
                             landmark_guided_dict_mapping,
                             calculate_selectivity)
from .experiment import ext_Experiment
from .constructors import (CompositeMap, MassGrid)
from .tools.file_io import read_features_from_asari_table



def join_tables(file_list_paths, parameters):
    '''
    Treat each table as a virtual "sample".
    Build a 'MassGrid' first by aligning masstracks btw tables;
    then map features on each composite mass track. 

    The IDs of mass tracks are not compatible with default Asari processing, 
    because the latter uses positional indexes but this function reuses old IDs. 
    Therefore, the sample-wise MassGrid building is standalone version here.

    The landmark_guided_mapping function has mass calibration, 
    triggered by correction_tolerance_ppm (1 default)

    '''
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
    for table in tables[1:]:
        print(f"Adding {table['input_file']}.")
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
        MG.MassGrid = NewGrid

    print(MG.MassGrid.iloc[:25, :].values)
    # export 
    MG.MassGrid.to_csv('test_MassGrid.csv')

    # Next RT mapping









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


