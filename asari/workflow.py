'''
ext_Experiment is the container for whole project data.
Heavy lifting is in constructors.CompositeMap, 
    which contains MassGrid for correspondence, and FeatureList from feature/peak detection.
Annotation is facilitated by jms-metabolite-services, mass2chem

        
'''
import multiprocessing as mp
import pymzml

from .experiment import *
from .chromatograms import extract_massTracks_        # extract_massTracks, 
from .mass_functions import *

#
# -----------------------------------------------------------------------------
#

def register_samples(list_input_files):
    '''
    Establish sample_id here, return sample_registry as a dictionary.
        samples.append( {
            'input_file': file,
            'ion_mode': parameters['mode'],
            # below will be populated after processing
            'name': '',
            'list_scan_numbers': [],
            'list_retention_time': [],
            'list_mass_tracks': [], 
            'anchor_mz_pairs': [],                    # mostly isotopic pairs to establish m/z anchors (landmarks)
            'number_anchor_mz_pairs': -1,
        } )
    , dict_meta_data, parameters
    dict_meta_data for future use.
    '''
    sample_registry = {}
    ii = 0
    for file in sorted(list_input_files):
        sample_registry[ii] = {'sample_id': ii, 'input_file': file}
        ii += 1
    return sample_registry

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
            (sample['sample_id'], sample['input_file'], parameters['mode'],
            mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height,
            os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle'),
            shared_dict
            )
        )
    return iters


def batch_EIC_from_samples_in_memory(sample_registry, parameters):
    pass

def batch_EIC_from_samples_to_mongo(sample_registry, parameters, cursor):
    pass



def batch_EIC_from_samples_ondisk(sample_registry, parameters):
    '''
    multiprocessing of mass track extraction of samples, and return shared_dict.
    More anchors mean better coverage of features, helpful to select reference sample.
    '''
    with mp.Manager() as manager:
        shared_dict = manager.dict()
        iters = make_iter_parameters(sample_registry, parameters, shared_dict)
        # print("Number of processes ", number_processes)
        with mp.Pool( parameters['multicores'] ) as pool:
            pool.starmap( single_sample_EICs_ondisk, iters )

        _d = dict(shared_dict)
    return _d

def single_sample_EICs_ondisk(sample_id, infile, ion_mode, 
                    mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, outfile, shared_dict):

    '''
    Process infile.
    shared_dict is used to pass back information, 
    sample_id: ('mzml_parsing', 'eic', number_anchor_mz_pairs, outfile)
    '''
    new = {'sample_id': sample_id, 'input_file': infile, 'ion_mode': ion_mode,}
    list_mass_tracks = []
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
                'rt_scan_numbers': track[1], 
                'intensity': track[2], 
                } )
            ii += 1

        new['list_mass_tracks'] = list_mass_tracks
        anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(list_mass_tracks, mz_tolerance_ppm=mz_tolerance_ppm)
        new['anchor_mz_pairs'] = anchor_mz_pairs
        new['number_anchor_mz_pairs'] = len(anchor_mz_pairs)
        shared_dict[new['sample_id']] = ('passed', 'passed', new['number_anchor_mz_pairs'], outfile)

        with open(outfile, 'wb') as f:
            pickle.dump(new, f, pickle.HIGHEST_PROTOCOL)

        print("Processed %s with %d mass tracks." %(os.path.basename(infile), ii))

    except:
        shared_dict[new['sample_id']] = ('failed', '', 0, '')
        print("mzML processing error in sample %s, skipped." %infile)
