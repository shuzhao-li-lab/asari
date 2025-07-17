'''
Functions for subcommand `analyze`
'''

import random
import pymzml
import numpy as np

from mass2chem.search import find_mzdiff_pairs_from_masstracks

from .chromatograms import extract_massTracks_ 
from .experiment import ext_Experiment
from .mass_functions import flatten_tuplelist
from .utils import bulk_process

def analyze_single_sample(infile, 
            mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000,
            parameters={}):
    '''
    Analyze single mzML file and print statistics. Used by asari subcommand `analyze`.
    This uses ext_Experiment and default HMDB data to estimate mass accuracy.

    Parameters
    ----------
    infile : str
        input mzML filepath.
    mz_tolerance_ppm : float, optional, default: 5
        m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    min_intensity : float, optional, default: 100
        minimal intentsity value to consider, also used to filter out 0s.
    min_timepoints : int, optional, default: 5
        minimal consecutive scans to be considered real signal.
    min_peak_height : float, optional, default: 1000
        a bin is not considered if the max intensity < min_peak_height.
    parameters : dict, optional, default: {}
        not used, just place holder to use ext_Experiment class.
    '''
    print("Analysis of %s\n" %infile)
    mz_landmarks, mode, min_peak_height_ = get_file_masstrack_stats(infile,mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)
    EE = ext_Experiment({}, parameters)
    EE.load_annotation_db()
    mass_accuracy_ratio = EE.KCD.evaluate_mass_accuracy_ratio(mz_landmarks, mode, mz_tolerance_ppm=10)
    print("  Mass accuracy is estimated as %2.1f ppm." %(mass_accuracy_ratio*1000000))
    print("\n")

def get_file_masstrack_stats(infile, parameters, return_sample=False):
                        #mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000,
                        #return_sample=False):
    '''
    Extract mass tracks from a file and get statistics.
    The ionization_mode is assumed on one scan, thus not supporting polarity switch in a single file.
    Landmark tracks are the mass tracks with paired 13C/12C pattern (based on m/z diff only).

    Parameters
    ----------
    infile : str
        input mzML filepath.
    mz_tolerance_ppm : float, optional, default: 5
        m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    min_intensity : float, optional, default: 100
        minimal intentsity value to consider, also used to filter out 0s.
    min_timepoints : int, optional, default: 5
        minimal consecutive scans to be considered real signal.
    min_peak_height : float, optional, default : 1000
        a bin is not considered if the max intensity < min_peak_height.
    return_sample : bool, optional, default: False
        if True, return full sample dictionary with mass tracks.
        Else, return _mz_landmarks_, ionization_mode, min_peak_height_.

        
    Note
    ----
    Example output::

        Total number of MS1 spectra: 741
        of which 0 are positive ionization mode.

        Assuming ionization mode is neg.

        Maxium retention time (sec): 299.818228
        m/z range: (min 80.011578, median 358.010062, max 997.616794)

        Found 14063 mass tracks.
        Found 4054 12C/13C isotopic pairs as landmarks.
        Max intensity in any landmark track:  687,714,048
        Minimal height of landmark tracks:  2,334 

        Mass accuracy was estimated on 203 matched values as -0.4 ppm.

    To-do: to add output info on instrumentation.
    '''
    mz_tolerance_ppm = parameters['mz_tolerance_ppm']
    min_intensity = parameters['min_intensity_threshold']
    min_timepoints = parameters['min_timepoints']
    min_peak_height = parameters['min_peak_height']
    new = {'sample_id': infile, 'input_file': infile, 'ion_mode': '',}
    list_mass_tracks = []
    
    jj = 0
    with pymzml.run.Reader(infile) as exp:
        for spec in exp:
            if spec.ms_level == 1:                          # MS Level 1 only
                if spec["positive scan"]:
                    ionization_mode = 'pos'
                    jj += 1
                else:
                    ionization_mode = 'neg'

    xdict = extract_massTracks_(infile, 
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
            'intensity': track[1], 
            } )
        ii += 1

    new['list_mass_tracks'] = list_mass_tracks
    _max_scan_number = len(new['list_scan_numbers'])
    anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(list_mass_tracks, 
                            list_mz_diff = [1.003355,], mz_tolerance_ppm=mz_tolerance_ppm)
    
    # Filter anchor_mz_pairs, to enforce apex within 5% of RT in a pair
    anchor_mz_pairs = match_mzdiff_pairs_by_rt(anchor_mz_pairs,
                            list_mass_tracks, _max_scan_number, 
                            rt_window_perc=0.05, min_scans_window=10)
    
    _mz_landmarks_ = flatten_tuplelist(anchor_mz_pairs)
    all_mz = [x['mz'] for x in list_mass_tracks]
    # down scale list_mass_tracks to verified by _mz_landmarks_
    list_mass_tracks = [list_mass_tracks[ii] for ii in _mz_landmarks_]
    peak_heights = [x['intensity'].max() for x in list_mass_tracks]
    max_peak_height = int(max(peak_heights))
    min_peak_height_ = int(min(peak_heights))
    # recommended_min_peak_height = int(0.5 * min(peak_heights))

    print("Total number of MS1 spectra: %d" %_max_scan_number)
    if jj > 1:
        print("of which %d are positive ionization mode." %jj)
    else:
        print("of which %d is positive ionization mode." %jj)
    print("Assuming ionization mode is %s.\n" %ionization_mode)

    print("Maxium retention time (sec): %f" %max(new['list_retention_time']))
    print("m/z range: (min %f, median %f, max %f)\n" %(np.min(all_mz), np.median(all_mz), np.max(all_mz)))

    print("Found %d mass tracks." %ii)
    print("Found %d 12C/13C isotopic pairs as landmarks." %len(anchor_mz_pairs))    
    print("Max intensity in any landmark track: ", f"{max_peak_height:,}")
    print("Minimal height of landmark tracks: ", f"{min_peak_height_:,}", "\n")

    if return_sample:
        return new
    else:
        return _mz_landmarks_, ionization_mode, min_peak_height_


def match_mzdiff_pairs_by_rt(matched_mz_pairs, 
            list_mass_tracks, max_scan_number, rt_window_perc=0.05, min_scans_window=10):
    '''
    Filter matched_mz_pairs by a window of retention time.
    Use max of (min_scans_window, max_scan_number*rt_window_perc) as match requirement on RT.
    This is used to enforce coelution of isotopes or adducts.

    Parameters
    ----------
    matched_mz_pairs : result form find_mzdiff_pairs_from_masstracks, [(id pair), ...]
    list_mass_tracks : [{'id_number': ii, 'mz': track[0], 'intensity': track[1]}, ...]
    
    Returns
    -------
    list of filtered matched_mz_pairs.
    '''
    _window = max(min_scans_window, max_scan_number*rt_window_perc)
    _dict_mass_tracks = {}
    for L in list_mass_tracks:
        _dict_mass_tracks[L['id_number']] = L['intensity']
    new = []
    for pair in matched_mz_pairs:
        apex1, apex2 = np.argmax(_dict_mass_tracks[pair[0]]), np.argmax(_dict_mass_tracks[pair[1]])
        if abs(apex1 - apex2) < _window:
            new.append(pair)
    return new


# -----------------------------------------------------------------------------
# estimate_min_peak_height

def __wrapped_get_file_masstrack_stats(job):
        try:
            infile, parameters = job
            return get_file_masstrack_stats(infile, parameters, return_sample=False)
        except:
            print("Error in analyzing ", infile)
            return None

def estimate_min_peak_height(list_input_files, parameters):
    '''
    Compute estimated min peak height. This gets min peak height from the andmark tracks in each file,
    which is the min of the mass tracks with paired 13C/12C pattern (based on m/z diff only).

    Parameters
    ----------
    list_input_files : list[str]
        input mzML filepaths.
    parameters : dict
        parameters for the analysis, including min_intensity_threshold, min_timepoints, min_peak_height,
        min_min_peak_height, mz_tolerance_ppm, dynamic_range, num_files_to_check.

    Returns
    -------
    int, an estimated parameter for min peak_height as half of the min verified landmark peaks.
    '''

    min_peak_height = parameters['min_peak_height']
    min_min_peak_height = parameters['min_min_peak_height']
    num_files_to_use = parameters['num_files_to_check']

    min_peak_height = max(min_min_peak_height, min_peak_height/parameters['dynamic_range'])

    estimated = []
    if num_files_to_use is None:
        selected = list_input_files
    elif len(list_input_files) <= num_files_to_use:
        selected = list_input_files
    elif isinstance(num_files_to_use, int):
        selected = random.sample(list_input_files, num_files_to_use)
    elif isinstance(num_files_to_use, float) and 0 <= num_files_to_use <= 1:
        to_select = max(int(num_files_to_use*len(list_input_files)), 10)
        selected = random.sample(list_input_files, min(to_select, len(list_input_files)))
    else:
        raise ValueError("num_files_to_use should be int or float between 0 and 1 or None")
    
    print("Estimating parameter for min peak_height based on ", selected)

    results = bulk_process(__wrapped_get_file_masstrack_stats, 
                           list(zip(selected, [parameters]*len(selected))),
                           jobs_per_worker=parameters['multicores'])
    for result, infile in zip(results, selected):
        if result is not None:
            _, _, min_peak_height_ = result
            estimated.append(min_peak_height_)
        else:
            print("Error in analyzing ", infile)
    recommended = int(0.5 * np.median(estimated))
    print("Estimated parameter for min peak_height is %d \n" %recommended)
    return max(recommended, min_min_peak_height)