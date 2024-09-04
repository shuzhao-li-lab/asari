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
import gzip

from .experiment import ext_Experiment
from .chromatograms import extract_massTracks_ 

from mass2chem.search import find_mzdiff_pairs_from_masstracks
from .samples import SimpleSample


# -----------------------------------------------------------------------------
# main workflow for `process`
# -----------------------------------------------------------------------------

def workflow_setup(list_input_files, parameters):
    sample_registry = register_samples(list_input_files)
    if parameters['database_mode'] == 'auto':
        if len(list_input_files) <= parameters['project_sample_number_small']:
            parameters['database_mode'] = 'memory'
        else:
            parameters['database_mode'] = 'ondisk'
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    parameters['time_stamp'] = ':'.join(time_stamp)
    create_export_folders(parameters, ''.join(time_stamp))
    shared_dict = batch_EIC_from_samples_(sample_registry, parameters)
    for sid, sam in sample_registry.items():
        sam.update(shared_dict[sid])
    EE = ext_Experiment(sample_registry, parameters)
    return EE

def workflow_cleanup(EE, list_input_files, parameters):
    if not parameters['pickle'] and parameters['database_mode'] != 'memory':
        remove_intermediate_pickles(parameters)

def process_project(list_input_files, parameters):
    EE = workflow_setup(list_input_files, parameters)
    workflows = {
        "GC": process_GC_project,
        "LC": process_LC_project
    }
    workflows[parameters['workflow']](EE, list_input_files, parameters)
    workflow_cleanup(EE, list_input_files, parameters)

def process_GC_project(EE, list_input_files, parameters):
    print("Processing Project using GC Workflow")
    EE.process_all_GC()
    #EE.export_all(anno=False)

def process_LC_project(EE, list_input_files, paramaters):
    print("Processing Project using LC Workflow")
    EE.process_all_LC()
    EE.export_all(anno=False)

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

def read_project_file(project_file, file_pattern='.mzML'):
    print("Working on ~~ %s ~~ \n\n" %project_file)
    return [os.path.abspath(l.strip()) for l in open(project_file).readlines() if file_pattern in l]

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
    return {ii: {'sample_id': ii, "input_file": file} for ii, file in enumerate(sorted(list_input_files))}

def create_export_folders(parameters, time_stamp):
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
    parameters['outdir'] = '_'.join([parameters['outdir'], parameters['project_name'], time_stamp]) 
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
    for sample in sample_registry.values():
        yield (
            sample['sample_id'], 
            sample['input_file'], 
            parameters['mode'], 
            parameters['database_mode'],
            parameters['mz_tolerance_ppm'], 
            parameters['min_intensity_threshold'], 
            parameters['min_timepoints'], 
            parameters['min_peak_height'], 
            parameters['intensity_multiplier'], 
            os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle')
            )


def wrapped_EIC(args):
    return single_sample_EICs_(*args)

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
    shared_dict = {}
    iters = list(make_iter_parameters(sample_registry, parameters))
    with mp.Pool(parameters['multicores']) as workers:
        while iters:
            batch = []
            for _ in range(parameters['multicores']):
                try:
                    batch.append(iters.pop())
                except:
                    pass
            for r in workers.imap(wrapped_EIC, batch):
                sid, sam = r
                SimpleSample.save(sam, parameters)
                del sam['sample_data']
                shared_dict[sid] = sam
    return shared_dict

def single_sample_EICs_(sample_id, 
                        infile, 
                        ion_mode, 
                        database_mode, 
                        mz_tolerance_ppm, 
                        min_intensity, 
                        min_timepoints, 
                        min_peak_height, 
                        intensity_multiplier,
                        outfile):
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

    Returns
    -------
    EIC_dict: dict
        used to update the mapping dict in the main thread

    See also
    --------
    batch_EIC_from_samples_, make_iter_parameters
    '''
    try:
        exp = pymzml.run.Reader(infile)
        try:
            import datetime
            timestamp = int(datetime.datetime.strptime(exp.info['start_time'], "%Y-%m-%dT%H:%M:%SZ").timestamp())
        except:
            timestamp = None
        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    min_intensity=min_intensity, 
                    min_timepoints=min_timepoints, 
                    min_peak_height=min_peak_height,
                    intensity_multiplier=intensity_multiplier)
        
        list_mass_tracks = [{'id_number': ii, 'mz': t[0], 'intensity': t[1]} for ii, t in enumerate(xdict['tracks'])]
        anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(list_mass_tracks, mz_tolerance_ppm=mz_tolerance_ppm)

        new = {
            'sample_id': sample_id,
            'input_file': infile,
            'ion_mode': ion_mode,
            'timestamp': timestamp,
            'max_scan_number': max(xdict['rt_numbers']),
            'list_mass_tracks': [{'id_number': ii, 'mz': t[0], 'intensity': t[1]} for ii, t in enumerate(xdict['tracks'])],
            'anchor_mz_pairs': anchor_mz_pairs,
            'number_anchor_mz_pairs': len(anchor_mz_pairs),
            'outfile': outfile

        }

        if database_mode == 'memory':
            _, to_return, size = "memory", new, 0

        print("Extracted %s to %d mass tracks." %(os.path.basename(infile), len(new['list_mass_tracks'])))
        return sample_id,  {
                "status:mzml_parsing": 'passed',
                "status:eic": 'passed',
                "data_location": outfile,
                "max_scan_number": new['max_scan_number'],
                "list_scan_numbers": xdict['rt_numbers'],
                "list_retention_time": xdict["rt_times"],
                "track_mzs": [(t['mz'], t['id_number']) for t in new['list_mass_tracks']],
                "number_anchor_mz_pairs": new['number_anchor_mz_pairs'],
                "anchor_mz_pairs": new['anchor_mz_pairs'],
                "sample_data": to_return,
                "acquisition_time": timestamp,
                "name": os.path.basename(infile).replace('.mzML', ''),
                "size": size
            }
    except:
        # xml.etree.ElementTree.ParseError
        print("mzML processing error in sample %s, skipped." %infile)
        return sample_id, {
                "status:mzml_parsing": 'failed',
                "status:eic": '',
                "data_location": '',
                "max_scan_number": 0,
                "list_scan_numbers": [],
                "list_retention_time": [],
                "track_mzs": [],
                "number_anchor_mz_pairs": 0,
                "anchor_mz_pairs": [],
                "sample_data": {},
                "size": 0
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
    time_stamp = ''.join([str(x) for x in time.localtime()[1:6]])
    create_export_folders(parameters, time_stamp)
    # samples are processed to mass tracks (EICs) here
    shared_dict = batch_EIC_from_samples_(sample_registry, parameters)
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
