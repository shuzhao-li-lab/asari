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
import json_tricks as json


from .experiment import ext_Experiment
from .chromatograms import extract_massTracks_ 
from .peaks import audit_mass_track
from .utils import bulk_process


from mass2chem.search import find_mzdiff_pairs_from_masstracks
import zipfile


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
            sam['status:mzml_parsing'], sam['status:eic'], sam['data_location'], sam['max_scan_number'], sam['list_scan_numbers'], sam['list_retention_time'], sam['track_mzs'], sam['number_anchor_mz_pairs'], sam['anchor_mz_pairs'], sam['sample_data'], sam['sparsified'] = shared_dict[sid]
            sam['name'] = os.path.basename(sam['input_file']).replace('.mzML', '')
        
        EE = ext_Experiment(sample_registry, parameters)
        #EE.process_all()
        #EE.export_all(anno=parameters["anno"])
    else:
        raise Exception("No data was processed, check the input files.")
    EE = ext_Experiment(sample_registry, parameters)
    return EE

def workflow_cleanup(EE, list_input_files, parameters):
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
    workflows = {
        "GC": process_GC_project,
        "LC": process_LC_project
    }
    workflows[parameters['workflow']](EE, list_input_files, parameters)
    workflow_cleanup(EE, list_input_files, parameters)

def process_GC_project(EE, list_input_files, parameters):
    print("Processing Project using GC Workflow")
    EE.process_all_GC()
    EE.export_all(anno=False)

def process_LC_project(EE, list_input_files, parameters):
    print("Processing Project using LC Workflow")
    EE.process_all_LC()
    EE.export_all(anno=True)

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

    os.mkdir(parameters['outdir'])
    try:
        os.mkdir(os.path.join(parameters['outdir'], 'export'))
    except FileExistsError:
        print("Warning, export directory exists, this is normal in some circumstances")
        
    parameters['export_outdir'] = os.path.join(parameters['outdir'], 'export')
    try:
        os.mkdir(os.path.join(parameters['outdir'], 'qaqc_reports'))
    except FileExistsError:
        print("Warning, qaqc_reports directory exists, this is normal in some circumstances")
        
    parameters['qaqc_reports_outdir'] = os.path.join(parameters['outdir'], 'qaqc_reports')

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
    assert parameters['reuse_intermediates'] is None, "Cannot remove intermediates when reuse_intermediates is set."
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
    iters = []
    mz_tolerance_ppm = parameters['mz_tolerance_ppm']
    min_intensity = parameters['min_intensity_threshold']
    min_timepoints = parameters['min_timepoints']
    min_peak_height = parameters['min_peak_height']
    for sample in sample_registry.values():
        outfile = os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle')
        iters.append((
            sample['sample_id'], 
            sample['input_file'], 
            parameters['mode'], 
            parameters['database_mode'],
            mz_tolerance_ppm, 
            min_intensity, 
            min_timepoints, 
            min_peak_height, 
            outfile,
            parameters['compress'],
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
    shared_dict = {}
    iters = make_iter_parameters(sample_registry, parameters)
    for result in bulk_process(single_sample_EICs_, iters, dask_ip=parameters['dask_ip']):
        # oversubscribe for dask
        shared_dict.update(result)
    return shared_dict

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
    
    sample_id, infile, ion_mode, database_mode, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, outfile, compress, parameters = job
    try:
        if parameters['reuse_intermediates']:
            for file in os.listdir(parameters['reuse_intermediates']):
                if os.path.basename(file).split(".")[0] == os.path.basename(outfile).split(".")[0]:
                    print("Reusing Intermediate: %s." %file)
                    from .samples import SimpleSample
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
                                        {}, 
                                        zipfile.is_zipfile(file))} 
                
        new = {'sample_id': sample_id, 'input_file': infile, 'ion_mode': ion_mode, 'list_mass_tracks': []}
        track_mzs = []
        
        xdict = extract_massTracks_(infile, 
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    min_intensity=min_intensity, 
                    min_timepoints=min_timepoints, 
                    min_peak_height=min_peak_height)
        new['max_scan_number'] = max(xdict['rt_numbers'])
        # already in asc ending order of m/z from extract_massTracks_, get_thousandth_regions
        
        for ii, track in enumerate(xdict['tracks']):
            audit_results = audit_mass_track(
                track[1],
                min_fwhm=round(0.5 * parameters['min_timepoints']),
                min_intensity_threshold=parameters['min_intensity_threshold'],
                min_peak_height=parameters['min_peak_height'],
                min_peak_ratio=parameters['signal_noise_ratio']
            )
            _baseline_, noise_level, scaling_factor, min_peak_height, list_intensity = audit_results
            new['list_mass_tracks'].append({
                'id_number': ii, 
                'mz': track[0],
                # 'rt_scan_numbers': track[1],       # format changed after v1.5
                'intensity': track[1], 
                'audit_results': {
                    "baseline": _baseline_,
                    "noise_level": noise_level,
                    "scaling_factor": scaling_factor,
                    "min_peak_height": min_peak_height,
                    "list_intensity": list_intensity
                    }
                })
            track_mzs.append( (track[0], ii) )       # keep a reconrd in sample registry for fast MassGrid align
        print("Extracted %s to %d mass tracks." %(os.path.basename(infile), ii))

        anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(new['list_mass_tracks'], mz_tolerance_ppm=mz_tolerance_ppm)
        # find_mzdiff_pairs_from_masstracks is not too sensitive to massTrack format
        new['anchor_mz_pairs'] = anchor_mz_pairs
        new['number_anchor_mz_pairs'] = len(anchor_mz_pairs)
        new['xdict'] = xdict
        new['track_mzs'] = track_mzs
        new['ms2_spectra'] = xdict['ms2_spectra']

        data_filepath = outfile
        if database_mode == 'ondisk':
            if not compress:
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
            print("\tExtracted to %s, %s MiB." % (data_filepath, round(os.path.getsize(data_filepath)/1024/1024, 2)))
            return {sample_id: ('passed', 
                                'passed', 
                                data_filepath,
                                new['max_scan_number'], 
                                xdict['rt_numbers'], 
                                xdict['rt_times'],
                                track_mzs,
                                new['number_anchor_mz_pairs'], 
                                anchor_mz_pairs,  
                                {}, 
                                compress)} 
            

        elif database_mode == 'memory':
            return {sample_id: ('passed', 
                                'passed', 
                                outfile,
                                new['max_scan_number'], 
                                xdict['rt_numbers'], 
                                xdict['rt_times'],
                                track_mzs,
                                new['number_anchor_mz_pairs'], 
                                anchor_mz_pairs,  
                                new, 
                                compress)}
    except Exception as e:
        print("Failed to extract: %s." %os.path.basename(infile))
        return {sample_id: ('failed', 'failed', None, None, None, None, None, None, None, None, compress)}



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
