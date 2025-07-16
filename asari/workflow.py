'''
ext_Experiment is the container for whole project data.
sample_registry is the dict to track sample information and status,
outside class to facilitate multiprocessing.
Heavy lifting is in constructors.CompositeMap, 
    which contains MassGrid for correspondence, and FeatureList from feature/peak detection.
Annotation is facilitated by jms-metabolite-services, mass2chem. 
'''
import time
import os
import pickle
import zipfile

import json_tricks as json 
from mass2chem.search import find_mzdiff_pairs_from_masstracks

from .experiment import ext_Experiment
from .chromatograms import extract_massTracks_ 
from .peaks import audit_mass_track
from .utils import bulk_process
from .samples import SimpleSample

# -----------------------------------------------------------------------------
# main workflow for `process`
# -----------------------------------------------------------------------------

def workflow_setup(list_input_files, parameters):
    '''
    This defines the main work flow in processing a list of LC-MS files,
    '''
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
    if shared_dict:
        for sid, sam in sample_registry.items():
            # todo - this is a bit of a mess, should be a simpler return or something?
            sam['status:mzml_parsing'], sam['status:eic'], sam['data_location'], sam['max_scan_number'], sam['list_scan_numbers'], sam['list_retention_time'], sam['track_mzs'], sam['number_anchor_mz_pairs'], sam['anchor_mz_pairs'], sam['sample_data'], sam['sparsified'], sam['acquisition_time'] = shared_dict[sid]
            sam['name'] = os.path.basename(sam['input_file']).replace('.mzML', '')
        EE = ext_Experiment(sample_registry, parameters)
    else:
        raise Exception("No data was processed, check the input files.")
    EE = ext_Experiment(sample_registry, parameters)
    return EE

def workflow_cleanup(EE, list_input_files, parameters):
    print(parameters)
    if not parameters['keep_intermediates'] and parameters['database_mode'] != 'memory':
        remove_intermediate_pickles(parameters)

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
    EE = workflow_setup(list_input_files, parameters)
    workflow_export_mode = {
        'GC': (EE.process_all_GC, 'GC'),
        'LC': (EE.process_all_LC, 'LC'),
        'LC_start': (EE.process_all_LC_start, 'LC')
    }
    print(f'Processing Experiment Using {parameters["workflow"]} Workflow...')
    workflow_export_mode[parameters['workflow']][0]()
    print("Exporting...")
    EE.export_all(anno=parameters['anno'], mode=workflow_export_mode[parameters['workflow']][1])
    print("Done")
    exit()
    workflow_cleanup(EE, list_input_files, parameters)

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
    '''
    This reads centroided LC-MS files from a file list.
    Returns a list of files that match file_pattern.

    This is useful for running Asari across filesystems, 
    or when another tool needs to request Asari processing but 
    we need to read the files in place.

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
    print("Working on ~~ %s ~~ \n\n" %project_file)
    with open(project_file) as project_filehandle:
        return [os.path.abspath(l.strip()) for l in project_filehandle.readlines() if file_pattern in l]

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
    return {ii : {'sample_id': ii, 'input_file': file} for ii, file in enumerate(list_input_files)}

def create_export_folders(parameters, time_stamp=None):
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
    if parameters['project_name'] in parameters['outdir']:
        print("Export folders already exist, will not overwrite.")
        return None

    if time_stamp is None:
        parameters['outdir'] = '_'.join([parameters['outdir'], parameters['project_name'], parameters['time_stamp_for_dir']])
    else:
        parameters['outdir'] = '_'.join([parameters['outdir'], parameters['project_name'], time_stamp])


    os.makedirs(parameters['outdir'])
    # add additional subdirectories here.
    for subdir in ['export', 'qaqc_reports', 'ms2_spectra']:
        try:
            os.mkdir(os.path.join(parameters['outdir'], subdir))
        except FileExistsError:
            print(f"Warning, {subdir} subdirectory already exists")
        parameters[f'{subdir}_outdir'] = os.path.join(parameters['outdir'], subdir)


    if parameters['reuse_intermediates']:
        parameters['tmp_pickle_dir'] = os.path.abspath(parameters['reuse_intermediates'])
        assert os.path.exists(os.path.abspath(parameters['reuse_intermediates'])), "The reuse_intermediates directory does not exist."
    else:
        parameters['tmp_pickle_dir'] = os.path.join(parameters['outdir'], 'pickle')
        try:
            os.mkdir(parameters['tmp_pickle_dir'])
        except FileExistsError:
            print("Warning, pickle directory exists, this is normal in some circumstances")

def remove_intermediate_pickles(parameters):
    '''
    Remove all temporary files under pickle/ to free up disk space. 

    Parameters
    ----------
    paramaters: dict
        passed from main.py to get tmp_pickle_dir
    '''
    assert parameters['reuse_intermediates'] is None, "Cannot remove when reuse_intermediates is set."
    print("Removing temporary pickle files...")
    for f in os.listdir(parameters['tmp_pickle_dir']):
        os.remove( os.path.join(parameters['tmp_pickle_dir'], f) )
    try:
        os.rmdir(parameters['tmp_pickle_dir'])
    except Exception as _:
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

    # todo - this should just pass parameters maybe? 
    iters = []
    for sample in sample_registry.values():
        outfile = os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle')
        iters.append((
            sample['sample_id'], 
            sample['input_file'], 
            outfile,
            parameters
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
    sample_data = {}
    for sample_datum in bulk_process(single_sample_EICs_, 
                                     make_iter_parameters(sample_registry, parameters), 
                                     dask_ip=parameters['dask_ip'], 
                                     jobs_per_worker=parameters['multicores']):
        sample_data.update(sample_datum)
    return sample_data

def single_sample_EICs_(job):
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
    #todo, maybe job should be a dict or something, else, we can just pass the parameters right?

    sample_id, infile, outfile, parameters = job
    if True:
        if parameters['reuse_intermediates']:
            if os.path.exists(os.path.join(parameters['reuse_intermediates'], 'pickle')): 
                parameters['reuse_intermediates'] = os.path.join(parameters['reuse_intermediates'], 'pickle')
            for file in os.listdir(parameters['reuse_intermediates']):
                if os.path.basename(file).split(".")[0] == os.path.basename(outfile).split(".")[0]:
                    print("Reusing Intermediate: %s." %file)
                    new = SimpleSample.load_intermediate(os.path.join(parameters['reuse_intermediates'], file))
                    return {sample_id: ('passed', 
                                        'passed', 
                                        os.path.join(parameters['reuse_intermediates'], file),
                                        new['max_scan_number'], 
                                        new['xdict']['rt_numbers'], 
                                        new['xdict']['rt_times'],
                                        new['track_mzs'],
                                        new['number_anchor_mz_pairs'], 
                                        new['anchor_mz_pairs'],  
                                        new['acquisition_time'],
                                        {}, 
                                        zipfile.is_zipfile(file))} 
                
        new = {
            'sample_id': sample_id, 
            'input_file': infile, 
            'ion_mode': parameters['mode'], 
            'list_mass_tracks': []
            }        
        xdict = extract_massTracks_(infile, 
                    mz_tolerance_ppm=parameters['mz_tolerance_ppm'], 
                    min_intensity=parameters['min_intensity_threshold'], 
                    min_timepoints=parameters['min_timepoints'], 
                    min_peak_height=parameters['min_peak_height'])
        if xdict['tracks']:
            audit_fields = "baseline", "noise_level", "scaling_factor", "min_peak_height", "list_intensity"
            audits = [dict(zip(audit_fields, audit_mass_track(t[1], 
                                       round(0.5 * parameters['min_timepoints']), 
                                       parameters['min_intensity_threshold'], 
                                       parameters['min_peak_height'], 
                                       parameters['signal_noise_ratio']))) for t in xdict['tracks']]
            new['list_mass_tracks'] = [{'id_number': ii, 'mz': t[0], 'intensity': t[1], 'audit_results': audits[ii]} for ii, t in enumerate(xdict['tracks'])]
            print("Extracted %s to %d mass tracks." %(os.path.basename(infile), len(xdict['tracks'])))

        anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(new['list_mass_tracks'], mz_tolerance_ppm=parameters['mz_tolerance_ppm'])
        # find_mzdiff_pairs_from_masstracks is not too sensitive to massTrack format
        new.update({
            'anchor_mz_pairs': anchor_mz_pairs,
            'number_anchor_mz_pairs': len(anchor_mz_pairs),
            'xdict': xdict,
            'track_mzs': [(t[0], ii) for ii, t in enumerate(xdict['tracks'])],
            'ms2_spectra': xdict['ms2_spectra'],
            'max_scan_number': max(xdict['rt_numbers']),
            'acquisition_time': xdict['acquisition_time']
        })

        #todo - clean this up
        data_filepath = outfile
        if parameters['database_mode'] == 'ondisk':
            if not parameters['compress']:
                if parameters['storage_format'] == 'pickle':
                    data_filepath = outfile
                    with open(outfile, 'wb') as f:
                        pickle.dump(new, f, pickle.HIGHEST_PROTOCOL)
                elif parameters['storage_format'] == 'json':
                    data_filepath = outfile.replace(".pickle", ".json")
                    with open(outfile.replace(".pickle", ".json"), 'w') as f:
                        json.dump(new, f)
            else:
                data_filepath = outfile.replace(".pickle", ".zip")
                with zipfile.ZipFile(data_filepath, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    if parameters['storage_format'] == 'pickle':
                        with zipf.open(os.path.basename(outfile), 'w') as f:
                            pickle.dump(new, f, pickle.HIGHEST_PROTOCOL)
                    elif parameters['storage_format'] == 'json':
                        with zipf.open(os.path.basename(outfile).replace(".pickle", ".json"), 'w') as f:
                            f.write(json.dumps(new).encode('utf-8'))
            print(f"\tExtracted to {data_filepath}, {round(os.path.getsize(data_filepath)/1024/1024, 2)} MiB.")
        return {sample_id: ('passed', 
                            'passed', 
                            data_filepath,
                            new['max_scan_number'], 
                            xdict['rt_numbers'], 
                            xdict['rt_times'],
                            new['track_mzs'],
                            new['number_anchor_mz_pairs'], 
                            anchor_mz_pairs,
                            new['acquisition_time'], 
                            new if parameters['database_mode'] == 'memory' else {}, 
                            parameters['compress'])} 
    #except Exception as _:
    #    print("Failed to extract: %s." %os.path.basename(infile))
    #    return {sample_id: ('failed', # status:mzml_parsing
    #                        'failed', # status:eic
    #                        None, # outfile
    #                        None, # max_scan_number
    #                        None, # rt_numbers
    #                        None, # rt_times
    #                        None, # track_mzs
    #                        None, # number_anchor_mz_pairs
    #                        None, # anchor_mz_pairs
    #                        None, # acquisition_time
    #                        None, # sample_data
    #                        parameters['compress'] # compress
    #                        )}

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
    parameters['database_mode'] = 'ondisk'
    time_stamp = ''.join([str(x) for x in time.localtime()[1:6]])
    create_export_folders(parameters, time_stamp)
    # samples are processed to mass tracks (EICs) here
    _ = batch_EIC_from_samples_(register_samples(register_samples(list_input_files)), parameters)
    print("XICs were stored as pickle objects under %s" %os.path.join(parameters['outdir'], 'pickle'))

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
