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

import pymzml

from .experiment import *
from .chromatograms import extract_massTracks_ 
from .mass_functions import *

# -----------------------------------------------------------------------------
# main workflow for `process`
# -----------------------------------------------------------------------------
def process_project(list_input_files, parameters):
    '''
    list_input_files: mzML files.
    parameters: dictionary of parameters.
    Will adjust disk writing according to database_mode
    '''
    sample_registry = register_samples(list_input_files)
    if parameters['database_mode'] == 'auto':
        if len(list_input_files) <= parameters['project_sample_number_small']:
            parameters['database_mode'] = 'memory'
        else:
            parameters['database_mode'] = 'ondisk'        # yet to implement mongo etc

    # time_stamp is `month daay hour minute second``
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    parameters['time_stamp'] = ':'.join(time_stamp)
    time_stamp = ''.join(time_stamp)
    # if parameters['database_mode'] == 'ondisk':
    create_export_folders(parameters, time_stamp)
        
    # samples are processed to mass tracks (EICs) here
    shared_dict = batch_EIC_from_samples_(sample_registry, parameters)
    for sid, sam in sample_registry.items():
        sam['status:mzml_parsing'], sam['status:eic'], sam['data_location'
            ], sam['max_scan_number'], sam['list_scan_numbers'], sam['list_retention_time'
            ], sam['track_mzs'
            ], sam['number_anchor_mz_pairs'], sam['anchor_mz_pairs'
            ], sam['sample_data'] = shared_dict[sid]

        sam['name'] = os.path.basename(sam['input_file']).replace('.mzML', '')
    
    EE = ext_Experiment(sample_registry, parameters)
    EE.process_all()
    EE.export_all(anno=parameters["anno"])

    if not parameters['pickle'] and parameters['database_mode'] != 'memory':
        remove_intermediate_pickles(parameters)

def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files from directory.
    '''
    print("Working on ~~ %s ~~ \n\n" %directory)
    return [os.path.join(directory, f) for f in os.listdir(directory) if file_pattern in f]

def register_samples(list_input_files):
    '''
    Establish sample_id here, return sample_registry as a dictionary.
    '''
    sample_registry = {}
    ii = 0
    for file in sorted(list_input_files):
        sample_registry[ii] = {'sample_id': ii, 'input_file': file}
        ii += 1
    return sample_registry

def create_export_folders(parameters, time_stamp):
    parameters['outdir'] = '_'.join([parameters['outdir'], parameters['project_name'], time_stamp]) 
    parameters['tmp_pickle_dir'] = os.path.join(parameters['outdir'], 'pickle')
    os.mkdir(parameters['outdir'])
    os.mkdir(parameters['tmp_pickle_dir'])
    os.mkdir(os.path.join(parameters['outdir'], 'export'))

def remove_intermediate_pickles(parameters):
    '''
    Remove all temporary files under pickle/ to free up disk space. 
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

def make_iter_parameters(sample_registry, parameters, shared_dict):
    '''
    Generate iterables for multiprocess.starmap for getting sample mass tracks.
    return:
    [('sample_id', input_file, mode, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, output_file, shared_dict), ...]
    '''
    iters = []
    mz_tolerance_ppm = parameters['mz_tolerance_ppm']
    min_intensity = parameters['min_intensity_threshold']
    min_timepoints = parameters['min_timepoints']
    min_peak_height = parameters['min_peak_height']
    for sample in sample_registry.values():
        outfile = os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle')
        iters.append(
            (sample['sample_id'], sample['input_file'], parameters['mode'], parameters['database_mode'],
            mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, outfile,
            shared_dict
            )
        )
    return iters

def batch_EIC_from_samples_(sample_registry, parameters):
    '''
    multiprocessing of mass track extraction of samples, and return shared_dict.
    More anchors mean better coverage of features, helpful to select reference sample.
    '''
    with mp.Manager() as manager:
        shared_dict = manager.dict()
        iters = make_iter_parameters(sample_registry, parameters, shared_dict)
        # print("Number of processes ", number_processes)
        with mp.Pool( parameters['multicores'] ) as pool:
            pool.starmap( single_sample_EICs_, iters )

        _d = dict(shared_dict)
    return _d

def single_sample_EICs_(sample_id, infile, ion_mode, database_mode,
                    mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, outfile, 
                    shared_dict):
    '''
    Process infile. 
    `shared_dict` is used to pass back information, thus critical. Designed here as
    sample_id: ('status:mzml_parsing', 'status:eic', outfile,
                max_scan_number, list_scan_numbers, list_retention_time,
                track_mzs,
                number_anchor_mz_pairs, anchor_mz_pairs, 
                dict({mass tracks}) )
    track_mzs or anchor_mz_pairs are used later for aligning m/z tracks.
    list of scans starts from 0. 
    '''
    new = {'sample_id': sample_id, 'input_file': infile, 'ion_mode': ion_mode,}
    list_mass_tracks = []
    track_mzs = []
    try:
        exp = pymzml.run.Reader(infile)
        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    min_intensity=min_intensity, 
                    min_timepoints=min_timepoints, 
                    min_peak_height=min_peak_height)
        new['max_scan_number'] = max(xdict['rt_numbers'])
        ii = 0
        # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
        for track in xdict['tracks']:                         
            list_mass_tracks.append( {
                'id_number': ii, 
                'mz': track[0],
                # 'rt_scan_numbers': track[1],              # format changed after v1.5
                'intensity': track[1], 
                } )
            track_mzs.append( (track[0], ii) )                  # keep a reconrd in sample registry for fast MassGrid align
            ii += 1

        new['list_mass_tracks'] = list_mass_tracks
        anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(list_mass_tracks, mz_tolerance_ppm=mz_tolerance_ppm)
        # find_mzdiff_pairs_from_masstracks is not too sensitive to massTrack format
        new['anchor_mz_pairs'] = anchor_mz_pairs
        new['number_anchor_mz_pairs'] = len(anchor_mz_pairs)

        if database_mode == 'ondisk':
            shared_dict[new['sample_id']] = ('passed', 'passed', outfile,
                                            new['max_scan_number'], xdict['rt_numbers'], xdict['rt_times'],
                                            track_mzs,
                                            new['number_anchor_mz_pairs'], anchor_mz_pairs,  
                                            {} )  
            with open(outfile, 'wb') as f:
                pickle.dump(new, f, pickle.HIGHEST_PROTOCOL)

        elif database_mode == 'memory':
            shared_dict[new['sample_id']] = ('passed', 'passed', outfile,
                                            new['max_scan_number'], xdict['rt_numbers'], xdict['rt_times'],
                                            track_mzs,
                                            new['number_anchor_mz_pairs'], anchor_mz_pairs,  
                                            new )

        print("Extracted %s to %d mass tracks." %(os.path.basename(infile), ii))

    except:
        # xml.etree.ElementTree.ParseError
        shared_dict[new['sample_id']] = ('failed', '', '', 
                                            0, [], [], [], 0, [], 
                                            {})
        print("mzML processing error in sample %s, skipped." %infile)


# -----------------------------------------------------------------------------
# main workflow for `xic`
# -----------------------------------------------------------------------------

def process_xics(list_input_files, parameters):
    '''
    Get XICs (mass tracks) from a folder of centroid mzML files.
    list_input_files: mzML files.
    parameters: dictionary of most parameters.
    '''
    sample_registry = register_samples(list_input_files)
    parameters['database_mode'] = 'ondisk'
    time_stamp = ''.join([str(x) for x in time.localtime()[1:6]])
    create_export_folders(parameters, time_stamp)
    # samples are processed to mass tracks (EICs) here
    shared_dict = batch_EIC_from_samples_(sample_registry, parameters)
    print("XICs were stored as pickle objects under %s" %os.path.join(parameters['outdir'], 'pickle'))

# -----------------------------------------------------------------------------
# main workflow for `extract`
# -----------------------------------------------------------------------------

def get_mz_list(infile):
    '''
    Get a list of m/z valuies from infile, which is tab or comma delimited and has m/z in the first column.
    Assuming first line as header and skip it.
    '''
    lines = open(infile).read().splitlines()[1:]
    return [float(x.split('\t')[0].split(",")[0]) for x in lines]

