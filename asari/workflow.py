'''
ext_Experiment is the container for whole project data.
Heavy lifting is in constructors.CompositeMap, 
    which contains MassGrid for correspondence, and FeatureList from feature/peak detection.
Annotation is facilitated by jms-metabolite-services, mass2chem. 



        if sample_N <= 10:
            Start by matching anchor pairs, then work thru remaining traces.
            1. create a reference based on anchor pairs
            2. align each sample to the reference_anchor_pairs
        elif sample_N <= 1000:     # more samples use a different method, as peak density will be apparent in more samples.
            single batch build
        else:                    # even more samples should be split to batches for performance reasons
            multiple batch build


'''
import time
import random
import multiprocessing as mp

import pymzml

from .experiment import *
from .chromatograms import extract_massTracks_ 
from .mass_functions import *

# -----------------------------------------------------------------------------
# main workflow for `process`
# -----------------------------------------------------------------------------

def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files.
    For OpenMS based XIC workflow, file_pattern='chrom.mzML'.
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
    parameters['outdir'] = '_'.join([parameters['project_name'], parameters['outdir'], time_stamp]) 
    os.mkdir(parameters['outdir'])
    os.mkdir(os.path.join(parameters['outdir'], 'pickle'))
    os.mkdir(os.path.join(parameters['outdir'], 'export'))


# main workflow for `process`
def process_project(list_input_files, parameters):
    '''
    list_input_files: Extracted ion chromatogram files.
    parameters: dictionary of most parameters.
    sample_registry include: 'status:mzml_parsing', 'status:eic', 'number_anchor_mz_pairs', 
                'data_location', track_mzs: (mz, masstrack_id), 'sample_data': {}
    '''
    sample_registry = register_samples(list_input_files)
    if parameters['database_mode'] == 'auto':
        if len(list_input_files) <= parameters['project_sample_number_small']:
            parameters['database_mode'] = 'memory'
        elif len(list_input_files) <= parameters['project_sample_number_large']:
            parameters['database_mode'] = 'ondisk'
        else:
            parameters['database_mode'] = 'ondisk'       # to implement 'split' 

    # time_stamp is `month daay hour minute second``
    time_stamp = ''.join([str(x) for x in time.localtime()[1:6]])
    if parameters['database_mode'] == 'ondisk':
        create_export_folders(parameters, time_stamp)
        
    # samples are processed to mass tracks (EICs) here
    shared_dict = batch_EIC_from_samples_(sample_registry, parameters)
    for sid, sam in sample_registry.items():
        sam['status:mzml_parsing'], sam['status:eic'], sam['number_anchor_mz_pairs'
                ], sam['data_location'], sam['track_mzs'], sam['sample_data'] = shared_dict[sid]
        sam['name'] = os.path.basename(sam['input_file']).replace('.mzML', '')
    
    # print(sample_registry)
    EE = ext_Experiment(sample_registry, parameters)
    EE.process_all()

    # export is separated, as it may be limited in some environments
    if parameters['database_mode'] != 'ondisk':
        create_export_folders(parameters, time_stamp)
    EE.export_all()


# -----------------------------------------------------------------------------
# estimate_min_peak_height

def estimate_min_peak_height(list_input_files, 
            mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=500,
            num_files_to_use=3):
    '''
    return an estimated parameter for min peak_height as half of the min verified landmark peaks.
    '''
    estimated = []
    if len(list_input_files) <= num_files_to_use:
        selected = list_input_files
    else:
        selected = random.sample(list_input_files, num_files_to_use)
    print("Estimating parameter for min peak_height based on ", selected)
    for infile in selected:
        try:
            mz_landmarks, mode, min_peak_height_ = get_file_masstrack_stats(infile,
                        mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)
                        # not all above parameters are used or relevant
            estimated.append(min_peak_height_)
        except:
            print("Error in analyzing ", infile)
    recommended = int(0.5 * np.median(estimated))
    print("Estimated parameter for min peak_height is %d \n" %recommended)
    return recommended

def ext_estimate_min_peak_height(list_input_files, 
            mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=500,
            num_files_to_use=3):
    '''
    return dict of
    ion mode and
    an estimated parameter for min peak_height as half of the min verified landmark peaks.

    Extended estimate_min_peak_height for X-asari use.
    '''
    estimated, _ionmode = [], []
    if len(list_input_files) <= num_files_to_use:
        selected = list_input_files
    else:
        selected = random.sample(list_input_files, num_files_to_use)
    print("Estimating parameter for min peak_height based on ", selected)
    for infile in selected:
        try:
            mz_landmarks, mode, min_peak_height_ = get_file_masstrack_stats(infile,
                        mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)
                        # not all above parameters are used or relevant
            estimated.append(min_peak_height_)
            _ionmode.append(mode)
        except:
            print("Error in analyzing ", infile)
    recommended = int(0.5 * np.median(estimated))
    if len(set(_ionmode)) > 1:
        print("Error occured due to inconsistent ion mode." )
        print(selected, _ionmode)
        return None
    else:
        return {'mode': _ionmode[0], 'min_peak_height': recommended}



# -----------------------------------------------------------------------------
# Mass track (EIC) extraction, multi-core parralization via multiprocessing

def make_iter_parameters(sample_registry, parameters, shared_dict):
    '''
    Generate iterables for multiprocess.starmap for getting sample mass tracks.
    return:
    [('sample_id', input_file, mode, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, output_file, shared_dict), ...]
    '''
    iters = []
    mz_tolerance_ppm = parameters['mz_tolerance']
    min_intensity = parameters['min_intensity_threshold']
    min_timepoints = parameters['min_timepoints']
    min_peak_height = parameters['min_peak_height']
    for sample in sample_registry.values():
        iters.append(
            (sample['sample_id'], sample['input_file'], parameters['mode'], parameters['database_mode'],
            mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height,
            os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle'),
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
    Process infile. `shared_dict` is used to pass back information, 
    sample_id: ('status:mzml_parsing', 'status:eic', number_anchor_mz_pairs, outfile, track_mzs, list_mass_tracks)
    track_mzs are used later for aligning m/z tracks.
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
        new['list_scan_numbers'] = xdict['rt_numbers']            # list of scans, starting from 0
        new['list_retention_time'] = xdict['rt_times']        # full RT time points in sample
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

            shared_dict[new['sample_id']] = ('passed', 'passed', new['number_anchor_mz_pairs'], outfile, track_mzs, {})
            with open(outfile, 'wb') as f:
                pickle.dump(new, f, pickle.HIGHEST_PROTOCOL)

        elif database_mode == 'memory':
            shared_dict[new['sample_id']] = ('passed', 'passed', new['number_anchor_mz_pairs'], outfile, track_mzs, 
                                            new)

        print("Processed %s with %d mass tracks." %(os.path.basename(infile), ii))

    except:
        # xml.etree.ElementTree.ParseError
        shared_dict[new['sample_id']] = ('failed', '', 0, '', [], {})
        print("mzML processing error in sample %s, skipped." %infile)
