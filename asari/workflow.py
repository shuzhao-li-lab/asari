'''
ext_Experiment is the container for whole project data.
sample_registry is the dict to track sample information and status,
outside class to facilitate multiprocessing.
Heavy lifting is in constructors.CompositeMap, 
    which contains MassGrid for correspondence, and FeatureList from feature/peak detection.
Annotation is facilitated by jms-metabolite-services, mass2chem. 
'''
import time
import multiprocessing as mp
import os
import pymzml
import pickle
import uuid
import numpy as np

from .mass_functions import flatten_tuplelist
from .experiment import ext_Experiment
from .chromatograms import extract_massTracks_ 

from mass2chem.search import find_mzdiff_pairs_from_masstracks


# -----------------------------------------------------------------------------
# main workflow for `process`
# -----------------------------------------------------------------------------
def process_project(list_input_files, parameters):
    '''
    This defines the main work flow in processing a list of LC-MS files, 
    creates the output folder with a time stamp, and uses sample registry to coordinate parallel processing.
    The whole project data are tracked in experiment.ext_Experiment class.

    Parameters
    ----------
    list_input_files : list[str]
        list of centroided mzML filepaths from LC-MS metabolomics. Usually found in a folder.
    parameters : dict
        parameter dictionary passed from main.py, 
        which imports from default_parameters and updates the dict by user arguments.

    Outputs
    -------
    A local folder with asari processing result, e.g::
    
        rsvstudy_asari_project_427105156
        ├── Annotated_empricalCompounds.json
        ├── Feature_annotation.tsv
        ├── export
        │   ├── _mass_grid_mapping.csv
        │   ├── cmap.pickle
        │   ├── full_Feature_table.tsv
        │   └── unique_compound__Feature_table.tsv
        ├── pickle
        │   ├── Blank_20210803_003.pickle
        │   ├── ...
        ├── preferred_Feature_table.tsv
        └── project.json
        

    The pickle folder is removed after the processing by default.
    '''
    sample_registry = register_samples(list_input_files)
    
    if parameters['database_mode'] == 'auto':
        if len(list_input_files) <= parameters['project_sample_number_small']:
            parameters['database_mode'] = 'memory'
        else:
            parameters['database_mode'] = 'ondisk'        # yet to implement mongo etc

    # time_stamp is `month daay hour minute second``
    # if parameters['database_mode'] == 'ondisk':
    create_export_folders(parameters)
        
    # samples are processed to mass tracks (EICs) here
    shared_dict = batch_EIC_from_samples_(sample_registry, parameters)
    for sid, sam in sample_registry.items():
        sam.update(shared_dict[sid])
        sam['name'] = os.path.basename(sam['input_file']).replace('.mzML', '')
    
    EE = ext_Experiment(sample_registry, parameters)
    EE.process_all()
    EE.export_all(anno=parameters["anno"])

    if not parameters['pickle'] and parameters['database_mode'] != 'memory':
        remove_intermediate_pickles(parameters)

def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files from directory.
    Returns a list of files that match file_pattern.

    Parameters
    ----------
    directory: str
        path to a directory containing mzML files
    file_pattern: str, optional, default: '.mzML'
        files with this substring will be ingested

    Return
    ------
    list of paths to files containing the file_pattern

    '''
    print("Working on ~~ %s ~~ \n\n" %directory)
    return [os.path.join(directory, f) for f in os.listdir(directory) if file_pattern in f]

def register_samples(list_input_files):
    '''
    Establish sample_id here, return sample_registry as a dictionary.

    Parameters
    ---------
    list_input_files: list[str]
        list of input filepaths, each representing a sample

    Return
    ------
    sample_registry, a dictionary of integer sample id's to filepaths
    '''

    return {ii: {'sample_id': ii, 'input_file': file} for ii, file in enumerate(list_input_files)}

def create_export_folders(parameters):
    '''
    Creates local directory for storing temporary files and output result.
    A time stamp is added to directory name to avoid overwriting existing projects.

    Parameters
    ----------
    paramaters: dict
        passed from main.py to get outdir and project_name
    time_stamp: str
        a time_stamp string to prevent overwriting existing projects
    '''

    parameters['outdir'] = os.path.abspath('_'.join([parameters['outdir'], parameters['project_name'], ''.join(parameters['time_stamp'])]))
    parameters['tmp_pickle_dir'] = os.path.join(parameters['outdir'], 'pickle')
    os.mkdir(parameters['outdir'])
    os.mkdir(parameters['tmp_pickle_dir'])
    os.mkdir(os.path.join(parameters['outdir'], 'export'))

def remove_intermediate_pickles(parameters):
    '''
    Remove all temporary files under pickle/ to free up disk space. 

    Parameters
    ----------
    paramaters: dict
        passed from main.py to get tmp_pickle_dir
    '''
    print("Removing temporary pickle files...")
    for f in os.listdir(parameters['tmp_pickle_dir']):
        os.remove( os.path.join(parameters['tmp_pickle_dir'], f) )
    try:
        os.rmdir(parameters['tmp_pickle_dir'])
    except:
        print("Failed to remove directory %s." %parameters['tmp_pickle_dir'])

# -----------------------------------------------------------------------------
# main workflow for `join`
# -----------------------------------------------------------------------------

def join_projects(list_input_projects, parameters):
    projects = []
    for project_path in list_input_projects:
        projects.append(pickle.load(open(os.path.join(project_path, "export/experiment.json"), 'rb')))
    create_export_folders(parameters)
    parameters.update(determine_meta_params(projects))
    paths = [os.path.join(project_path, "export/cmap.pickle") for project_path in list_input_projects]
    meta_sample_registry = {ii: generate_mock_sample_registry(project, id=ii, path=path) for ii, (project, path) in enumerate(zip(projects, paths))}
    meta_experiment = ext_Experiment(meta_sample_registry, parameters)
    for project in projects:
        project.meta_experiment = meta_experiment
    meta_experiment.process_all()
    meta_experiment.export_all()
    meta_experiment.export_log()

def determine_meta_params(list_projects, mode='auto'):
    if mode != "auto":
        meta_param_mode_map = {
            "min": min,
            "median": np.median
        }
        return {
            'min_peak_height': meta_param_mode_map[mode]([x.min_peak_height for x in list_projects]),
            'min_prominence_threshold': meta_param_mode_map[mode]([x.min_prominence_threshold for x in list_projects])
        }
    elif mode == 'auto':
        min_heights = []        
        for project in list_projects:
            composite_mass_tracks = project.CMAP.composite_mass_tracks
            anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(composite_mass_tracks.values(), list_mz_diff = [1.003355,], mz_tolerance_ppm=5)
            _mz_landmarks_ = flatten_tuplelist(anchor_mz_pairs)
            composite_mass_tracks = [composite_mass_tracks[x] for x in _mz_landmarks_]
            min_heights.append(int(min([x['intensity'].max() for x in composite_mass_tracks])))
        recommended = np.median(min_heights)
        return {
            'min_peak_height': recommended,
            'min_prominence_threshold': .33 * recommended
        }
    else:
        raise Exception()


def generate_mock_sample_registry(project, id=None, path=None):
    dict_scan_rtime = project.CMAP.dict_scan_rtime
    list_scan_times = list(dict_scan_rtime.keys())
    list_mass_tracks = list(project.CMAP.composite_mass_tracks.values())
    anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(list_mass_tracks)
    return {
            'sample_id': id,
            'name': path,
            'input_file': path,
            'status:mzml_parsing': 'passed',
            'status:eic': 'passed',
            'data_location': path,
            'max_scan_number': max(list_scan_times),
            'list_scan_numbers': list_scan_times,
            'list_retention_time': list(dict_scan_rtime.values()),
            'track_mzs': [(z['mz'], z['id_number']) for z in list_mass_tracks],
            'anchor_mz_pairs': anchor_mz_pairs,
            'number_anchor_mz_pairs': len(anchor_mz_pairs),
            'composite_of': [project],
            'sample_mapping': project.CMAP.MassGrid,
            'sample_data': {
                'sample_id': id,
                'input_file': None,
                'ion_mode': project.mode,
                'max_scan_number': max(list_scan_times),
                'list_mass_tracks': list_mass_tracks, 
                'anchor_mz_pairs': anchor_mz_pairs,
                'number_anchor_mz_pairs': len(anchor_mz_pairs)
                }
            }



# -----------------------------------------------------------------------------
# Mass track (EIC) extraction, multi-core parralization via multiprocessing

def make_iter_parameters(sample_registry, parameters):
    '''
    Generate iterables for multiprocess.starmap for getting sample mass tracks.

    Parameters
    ----------
    sample_registry : dict
        dictionary like {'sample_id': ii, 'input_file': file}, 
        generated by register_samples.
    parameters : dict
        parameter dictionary passed from main.py, 
        which imports from default_parameters and updates the dict by user arguments.
    shared_dict : dict
        dictionary object used to pass data btw multiple processing.

    Returns
    -------
    A list of iterative parameters, e.g.
    [('sample_id', input_file, mode, mz_tolerance_ppm, min_intensity, min_timepoints, 
    min_peak_height, output_file, shared_dict), ...]
    '''
    iters = []
    for sample in sample_registry.values():
        outfile = os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle')
        iters.append((sample['sample_id'], 
                      sample['input_file'], 
                      parameters['mode'], 
                      parameters['database_mode'],
                      parameters['mz_tolerance_ppm'], 
                      parameters['min_intensity_threshold'], 
                      parameters['min_timepoints'], 
                      parameters['min_peak_height'], 
                      outfile
            )
        )
    return iters

def batch_EIC_from_samples_(sample_registry, parameters):
    '''
    Batch extraction of mass tracks from samples via multiprocessing.
    
    Parameters
    ----------
    sample_registry : dict
        dictionary like {'sample_id': ii, 'input_file': file}, 
        generated by register_samples.
    parameters : dict
        parameter dictionary passed from main.py, 
        which imports from default_parameters and updates the dict by user arguments.

    Returns
    -------
    shared_dict : dict
        dictionary object used to pass data btw multiple processing.

    See also
    --------
    single_sample_EICs_
    '''

    with mp.Pool(parameters['multicores']) as pool:
        EICs = pool.starmap(single_sample_EICs_, make_iter_parameters(sample_registry, parameters))
        return {k: v for d in EICs for k, v in d.items()}

def single_sample_EICs_(sample_id, infile, ion_mode, database_mode,
                    mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, outfile):
    '''
    Extraction of mass tracks from a single sample. Used by multiprocessing in batch_EIC_from_samples_.
    `shared_dict` is used to pass back information, thus critical. Designed here as
    sample_id: 
        ('status:mzml_parsing', 'status:eic', outfile,
        max_scan_number, list_scan_numbers, list_retention_time,
        track_mzs,
        number_anchor_mz_pairs, anchor_mz_pairs, 
        dict({mass tracks}) )

    track_mzs or anchor_mz_pairs are used later for aligning m/z tracks.

    list of scans starts from 0. 
    
    `anchor_mz_pairs` are defined as m/z pairs that match to 13C/12C pattern.
    More anchors mean better coverage of features, helpful to select reference sample.

    Parameters
    ----------
    sample_id : int
        sample id passed from sample_registry by make_iter_parameters.
    infile : str
        input mzML filepath, passed from sample_registry by make_iter_parameters.
    ion_mode: str
        from parameter dictionary, marks if acquisition was in positive or negative mode 
    database_mode: str
        from parameter dictionary, marks if intermediatesare kept on disk or in memory
    mz_tolerance_ppm: float
        from parameter dictionray, the assumed mz resolution of the instrument
    min_intensity: float
        peaks below this value are ignored
    min_timepoints: int
        then number of time points a peak must span to be considered a peak
    min_peak_height: float
        peaks below this height are ignored
    outfile : str
        where the output will be written
        passed from parameter dictionary by make_iter_parameters.
    shared_dict : dict
        dictionary object used to pass data btw multiple processing.

    Updates
    -------
    shared_dict : dict
        dictionary object used to pass data btw multiple processing.

    See also
    --------
    batch_EIC_from_samples_, make_iter_parameters
    '''

    exp = pymzml.run.Reader(infile)
    xdict = extract_massTracks_(exp, 
                mz_tolerance_ppm=mz_tolerance_ppm, 
                min_intensity=min_intensity, 
                min_timepoints=min_timepoints, 
                min_peak_height=min_peak_height)

    # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
    list_mass_tracks = []
    track_mzs = []
    for ii, track in enumerate(xdict['tracks']):                         
        list_mass_tracks.append( {
            'id_number': ii, 
            'mz': track[0],
            'intensity': track[1], 
            } )
        track_mzs.append((track[0], ii))       # keep a reconrd in sample registry for fast MassGrid align
    anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(list_mass_tracks, mz_tolerance_ppm=mz_tolerance_ppm)
    

    database_mode_returns = {
        'memory': {
            'sample_id': sample_id,
            'input_file': infile,
            'ion_mode': ion_mode,
            'max_scan_number': max(xdict['rt_numbers']),
            'list_mass_tracks': list_mass_tracks,
            'anchor_mz_pairs': anchor_mz_pairs,
            # find_mzdiff_pairs_from_masstracks is not too sensitive to massTrack format
            'number_anchor_mz_pairs': len(anchor_mz_pairs)
            },
        'ondisk': {},
        'auto': {},
        'smart': {}
    }


    if database_mode != "memory":
        with open(outfile, 'wb') as f:
            pickle.dump(database_mode_returns["memory"], f, pickle.HIGHEST_PROTOCOL)


    print("Extracted %s to %d mass tracks." %(os.path.basename(infile), ii))
    return {
        sample_id: {
            'status:mzml_parsing': 'passed',
            "status:eic": 'passed',
            'data_location': outfile,
            'max_scan_number': max(xdict['rt_numbers']),
            'list_scan_numbers': xdict['rt_numbers'],
            'list_retention_time': xdict['rt_times'],
            'track_mzs': track_mzs,
            'number_anchor_mz_pairs': len(anchor_mz_pairs),
            'anchor_mz_pairs': anchor_mz_pairs,
            'sample_data': database_mode_returns[database_mode],
            'mem_footprint': os.stat(outfile).st_size * 2 if database_mode != "memory" else 0
            }
    }

# -----------------------------------------------------------------------------
# main workflow for `xic`
# -----------------------------------------------------------------------------

def process_xics(list_input_files, parameters):
    '''
    Get XICs (aka EICs or mass tracks) from a folder of centroid mzML files 
    and store in local pickle files.

    Parameters
    ----------
    list_input_files : list[str]
        list of centroided mzML filepaths from LC-MS metabolomics. 
        Usually found in a folder.
    parameters : dict
        parameter dictionary passed from main.py, 
        which imports from default_parameters and updates the dict by user arguments.

    Outputs
    -------
    A local folder with asari extracted XICs (aka EICs or mass tracks), 
    without full processing, in pickle files.
    '''
    sample_registry = register_samples(list_input_files)
    parameters['database_mode'] = 'ondisk'
    create_export_folders(parameters)
    # samples are processed to mass tracks (EICs) here
    _ = batch_EIC_from_samples_(sample_registry, parameters)
    print("XICs were stored as pickle objects under %s" 
          %os.path.join(parameters['outdir'], 'pickle'))

# -----------------------------------------------------------------------------
# main workflow for `extract`, short cut by passing mzlist to parameters['target']
# -----------------------------------------------------------------------------

def get_mz_list(infile):
    '''
    Get a list of m/z valuies from infile, to be used as targets to extract from LC-MS features.
    Currently, `extract` is a short cut function - asari still processes full feature tables, 
    then searchs for the given targets.

    Parameters
    ----------
    infile : str
        filepath to input table, which is tab or comma delimited and has m/z in the first column,
        header as first row.

    Returns
    -------
    A list of m/z values.
    '''
    lines = open(infile).read().splitlines()[1:]
    return [float(x.split('\t')[0].split(",")[0]) for x in lines]
