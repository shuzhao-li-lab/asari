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
import sqlite3

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
    workflow_export_mode[parameters['workflow']][0]()
    EE.export_all(anno=parameters['anno'], mode=workflow_export_mode[parameters['workflow']][1])
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
        files with this substring ßbe ingested

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
        time_stamp = parameters['time_stamp_for_dir']
    
    parameters['outdir'] = os.path.join(os.path.abspath(parameters['outdir']), parameters['project_name'] + "_" + time_stamp)


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
    for sample_datum in bulk_process(single_sample_EICs, 
                                     make_iter_parameters(sample_registry, parameters), 
                                     dask_ip=parameters['dask_ip'],
                                     num_workers=parameters['multicores']):
        sample_data.update(sample_datum)
    return sample_data

def _save_to_database(sample_data, db_path):
    """
    Saves a sample's data into the appropriate tables in an SQLite database.

    Assumes the database has two tables: 'samples' for metadata and 
    'mass_tracks' for the detailed track data points.

    Parameters
    ----------
    sample_data : dict
        The sample data dictionary to save.
    db_path : str
        The file path of the SQLite database.
    """
    # Use a 'with' statement to ensure the connection is safely managed
    raise NotImplementedError()
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        
        # Use a transaction to ensure data integrity
        try:
            # 1. Insert the sample's metadata into the 'samples' table
            # Complex Python objects are serialized to JSON strings for storage
            cursor.execute('''
                INSERT OR REPLACE INTO samples (
                    sample_id, name, input_file, data_location, max_scan_number,
                    rt_numbers, list_retention_time, track_mzs, anchor_mz_pairs, ms2_spectra
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                sample_data.get('sample_id'),
                sample_data.get('name', os.path.basename(sample_data.get('input_file', ''))),
                sample_data.get('input_file'),
                # The data_location will be the URI pointing to this DB record
                f"sqlite://{os.path.abspath(db_path)}?sample_id={sample_data.get('sample_id')}",
                sample_data.get('max_scan_number'),
                json.dumps(sample_data.get('rt_numbers', [])),
                json.dumps(sample_data.get('list_retention_time', [])),
                json.dumps(sample_data.get('track_mzs', [])),
                json.dumps(sample_data.get('anchor_mz_pairs', [])),
                json.dumps(sample_data.get('ms2_spectra', []))
            ))

            # 2. Prepare the detailed mass track data for batch insertion
            track_data_to_insert = []
            for track in sample_data.get('list_mass_tracks', []):
                track_id = track.get('id_number')
                mz = track.get('mz')
                
                # Assume intensities correspond 1:1 with scan numbers in xdict
                # A more robust implementation might store scans directly with intensities
                scans = sample_data.get('xdict', {}).get('rt_numbers', [])
                
                for i, intensity in enumerate(track.get('intensity', [])):
                    # This assumes the intensity array matches the full scan list
                    if i < len(scans):
                        scan_number = scans[i]
                        track_data_to_insert.append(
                            (sample_data.get('sample_id'), track_id, mz, scan_number, intensity)
                        )
            
            # 3. Clear any old track data for this sample before inserting new data
            cursor.execute("DELETE FROM mass_tracks WHERE sample_id = ?", (sample_data.get('sample_id'),))
            
            # 4. Perform a batch insert into the 'mass_tracks' table
            if track_data_to_insert:
                cursor.executemany('''
                    INSERT INTO mass_tracks (sample_id, track_id, mz, scan_number, intensity)
                    VALUES (?, ?, ?, ?, ?)
                ''', track_data_to_insert)

            # Commit the transaction to save all changes
            conn.commit()
            
        except Exception as e:
            print(f"❌ Database Error: {e}. Rolling back transaction.")
            conn.rollback()
            raise # Re-raise the exception after rolling back

def _save_sample_data(data, outfile_base, storage_format, compress):
    """
    Handles saving sample data to disk in the specified format, with optional compression.
    Falls back to uncompressed saving if compression fails.
    """
    was_compressed = compress
    
    # Determine the file extension and writer function
    if storage_format == 'pickle':
        ext = '.pickle'
        writer = lambda f: pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
        binary_mode = True
    elif storage_format == 'json':
        ext = '.json'
        writer = lambda f: f.write(json.dumps(data).encode('utf-8'))
        binary_mode = True # Write bytes to zip
    else:
        raise ValueError(f"Unsupported storage format: {storage_format}")

    outfile = outfile_base.replace('.pickle', ext)

    if not compress:
        with open(outfile, 'wb' if binary_mode else 'w') as f:
            writer(f)
        return outfile, was_compressed

    # Attempt to save compressed
    zip_outfile = outfile_base.replace('.pickle', '.zip')
    try:
        with zipfile.ZipFile(zip_outfile, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # The name of the file inside the zip archive
            arcname = os.path.basename(outfile)
            with zipf.open(arcname, 'w', force_zip64=True) as f:
                writer(f)
        return zip_outfile, was_compressed
    except Exception as e:
        print(f"⚠️ Warning: Compression failed for {os.path.basename(outfile)} ({e}). Saving uncompressed.")
        # Fallback to uncompressed
        was_compressed = False
        with open(outfile, 'wb' if binary_mode else 'w') as f:
            writer(f)
        return outfile, was_compressed

def _get_failure_payload(sample_id):
    """Generates the standard return tuple for a failed job."""
    return {sample_id: ('failed', 'failed', None, None, None, None, None, None, None, None, {}, False)}

def single_sample_EICs(job):
    """
    Extracts mass tracks from a single sample mzML file.

    This function is designed to be used by a multiprocessing pool. It processes
    a single file, extracts ion chromatograms (mass tracks), identifies 
    13C/12C isotope pairs (anchor pairs), and saves the result.

    Parameters
    ----------
    job : tuple
        A tuple containing the necessary parameters for processing a single sample:
        - sample_id (int): A unique identifier for the sample.
        - infile (str): The full path to the input mzML file.
        - outfile (str): The base path for the output file (e.g., '.pickle').
        - parameters (dict): A dictionary of processing parameters, including:
            - 'reuse_intermediates' (str): Path to look for pre-computed results.
            - 'mode' (str): Ion mode ('positive' or 'negative').
            - 'database_mode' (str): 'ondisk' or 'memory'.
            - 'mz_tolerance_ppm' (float): Mass tolerance in PPM.
            - 'min_intensity_threshold' (float): Minimum intensity to consider.
            - 'min_timepoints' (int): Minimum length of a valid mass track.
            - 'min_peak_height' (float): Minimum peak height for auditing.
            - 'signal_noise_ratio' (float): S/N ratio for peak auditing.
            - 'storage_format' (str): 'pickle' or 'json'.
            - 'compress' (bool): Whether to compress the output file.

    Returns
    -------
    dict
        A dictionary with the sample_id as the key. The value is a tuple 
        containing the processing results and status:
        (
            str,  # 'passed' or 'failed' (mzML parsing status)
            str,  # 'passed' or 'failed' (EIC status)
            str,  # Output file path
            int,  # Maximum scan number
            list, # List of scan numbers (timepoints)
            list, # List of retention times
            list, # List of (m/z, track_index) tuples
            int,  # Number of anchor m/z pairs found
            list, # List of anchor m/z pairs
            float, # Acquisition time in minutes
            dict, # The extracted data dict (if database_mode=='memory')
            bool  # True if the output file was successfully compressed
        )
    """
    sample_id, infile, outfile, parameters = job

    #try:
    if True:
        # 1. Check for and reuse intermediate files if available
        if interm_dir := parameters.get('reuse_intermediates'):
            pickle_dir = os.path.join(interm_dir, 'pickle')
            search_dir = pickle_dir if os.path.exists(pickle_dir) else interm_dir
            
            base_name = os.path.basename(outfile).split('.')[0]
            for file in os.listdir(search_dir):
                if file.startswith(base_name):
                    print(f"Reusing intermediate: {file}")
                    reused_path = os.path.join(search_dir, file)
                    data = SimpleSample.load_intermediate(reused_path)
                    return {sample_id: ('passed', 'passed', reused_path,
                                        data['max_scan_number'], data['xdict']['rt_numbers'],
                                        data['xdict']['rt_times'], data['track_mzs'],
                                        data['number_anchor_mz_pairs'], data['anchor_mz_pairs'],
                                        data['acquisition_time'], {}, zipfile.is_zipfile(reused_path))}

        # 2. Extract mass tracks from the mzML file
        xdict = extract_massTracks_(
            infile,
            mz_tolerance_ppm=parameters['mz_tolerance_ppm'],
            min_intensity=parameters['min_intensity_threshold'],
            min_timepoints=parameters['min_timepoints'],
            min_peak_height=parameters['min_peak_height']
        )
        
        sample_data = {
            'sample_id': sample_id,
            'input_file': infile,
            'ion_mode': parameters['mode'],
            'list_mass_tracks': []
        }

        # 3. Audit extracted tracks and format results
        if xdict['tracks']:
            #print(f"Extracted {len(xdict['tracks'])} mass tracks from {os.path.basename(infile)}.")
            audit_fields = ("baseline", "noise_level", "scaling_factor", "min_peak_height", "list_intensity")
            audits = [
                dict(zip(audit_fields, audit_mass_track(
                    t[1], round(0.5 * parameters['min_timepoints']),
                    parameters['min_intensity_threshold'], parameters['min_peak_height'],
                    parameters['signal_noise_ratio']
                ))) for t in xdict['tracks']
            ]
            sample_data['list_mass_tracks'] = [
                {'id_number': i, 'mz': t[0], 'intensity': t[1], 'audit_results': audits[i]}
                for i, t in enumerate(xdict['tracks'])
            ]

        # 4. Find anchor pairs (e.g., 13C/12C isotopes)
        anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(
            sample_data['list_mass_tracks'], mz_tolerance_ppm=parameters['mz_tolerance_ppm']
        )
        
        # 5. Assemble final data dictionary
        sample_data.update({
            'anchor_mz_pairs': anchor_mz_pairs,
            'number_anchor_mz_pairs': len(anchor_mz_pairs),
            'xdict': xdict,
            'track_mzs': [(t[0], i) for i, t in enumerate(xdict['tracks'])],
            'ms2_spectra': xdict['ms2_spectra'],
            'max_scan_number': max(xdict['rt_numbers']) if xdict['rt_numbers'] else 0,
            'acquisition_time': xdict['acquisition_time']
        })

        # 6. Save results to the specified backend and generate the data location URI
        data_location_uri = None
        was_compressed = parameters.get('compress', False)
        
        # database_mode now determines the storage backend
        storage_backend = parameters.get('database_mode', 'ondisk')

        if storage_backend == 'ondisk':
            data_filepath, was_compressed = _save_sample_data(
                sample_data,
                outfile,
                parameters.get('storage_format', 'pickle'),
                was_compressed
            )
            data_location_uri = f"file://{os.path.abspath(data_filepath)}"
            #print(f"\tSaved to file: {data_filepath}")

        elif storage_backend == 'memory':
            # In memory mode, there is no URI, but the data is passed directly
            pass

        # 7. Construct and return the success payload with the new URI
        return {sample_id: (
            'passed',
            'passed',
            data_location_uri,  # CRITICAL: This is now a URI
            sample_data.get('max_scan_number', 0),
            xdict['rt_numbers'],
            xdict['rt_times'],
            sample_data.get('track_mzs', []),
            sample_data.get('number_anchor_mz_pairs', 0),
            anchor_mz_pairs,
            sample_data.get('acquisition_time', 0.0),
            sample_data if storage_backend == 'memory' else {},
            was_compressed
        )}

    #except Exception as e:
    #    print(f"Error processing {os.path.basename(infile)}: {e}")
    #    # Return a failure payload
    #    return {sample_id: ('failed', 'failed', None, None, None, None, None, None, None, None, {}, False)}



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
