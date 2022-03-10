'''
ext_Experiment is the container for whole project data.
Heavy lifting is in constructors.CompositeMap, 
    which contains MassGrid for correspondence, and FeatureList from feature/peak detection.
Annotation is facilitated by jms-metabolite-services, mass2chem

        pool = mp.Pool( min(mp.cpu_count(), self.experiment.parameters['multicores']) )
        pool.starmap(
            self.process_and_add_sample, 
            [f for f in self.list_input_files if f not in self.experiment.initiation_samples]
        )
        pool.close()
        

'''

import multiprocessing as mp
# from pyopenms import MSExperiment, MzMLFile

from .experiment import *
from .chromatograms import extract_massTracks_        # extract_massTracks, 
from .mass_functions import *


#
# -----------------------------------------------------------------------------
#

def register_samples(list_input_files, dict_meta_data, parameters):
    '''
    dict_meta_data for future use.
    '''
    samples = []
    for file in list_input_files:
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
    return samples

def make_iter_parameters(samples, parameters):
    '''
    Generate iterables for multiprocess.starmap for getting sample mass tracks.
    return:
    [(input_file, mode, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, output_file), ...]
    '''
    iters = []
    mz_tolerance_ppm = parameters['mz_tolerance']
    min_intensity = parameters['min_intensity_threshold']
    min_timepoints = parameters['min_timepoints']
    min_peak_height = parameters['min_peak_height']
    for sample in samples:
        iters.append(
            (sample['input_file'], sample['ion_mode'],
            mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height,
            os.path.join(parameters['outdir'], 'pickle', os.path.basename(sample['input_file']).replace('.mzML', '')+'.pickle')
            )
        )
    return iters


def batch_EIC_from_samples_in_memory(samples, parameters):
    pass

def batch_EIC_from_samples_to_mongo(samples, parameters, cursor):
    pass



def batch_EIC_from_samples_ondisk(samples, parameters):
    '''
    multiprocessing of mass track extraction
    '''
    number_processes = min(mp.cpu_count(), parameters['multicores'])
    iters = make_iter_parameters(samples, parameters)
    print("Number of processes ", number_processes)
    pool = mp.Pool( number_processes )
    pool.starmap( single_sample_EICs_ondisk, iters )
    pool.close()


def single_sample_EICs_ondisk(infile, ion_mode, 
                    mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height, outfile):
    name = os.path.basename(infile).replace('.mzML', '')
    new = { 'input_file': infile, 'ion_mode': ion_mode, 'name': name,}
    list_mass_tracks = []

    print(new)

    exp = MSExperiment()
    MzMLFile().load(infile, exp)
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

    with open(outfile, 'wb') as f:
        pickle.dump(new, f, pickle.HIGHEST_PROTOCOL)

    print("Processed %s with %d mass tracks." %(name, ii))
