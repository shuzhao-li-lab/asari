'''
Functions for subcommand `analyze`
'''

import random
import pymzml
import numpy as np

from .chromatograms import extract_massTracks_ 
from .experiment import ext_Experiment
from .mass_functions import flatten_tuplelist

from mass2chem.search import find_mzdiff_pairs_from_masstracks



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
    mz_landmarks, mode, min_peak_height_ = get_file_masstrack_stats(infile,
                        mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)

    EE = ext_Experiment({}, parameters)
    EE.load_annotation_db()
    mass_accuracy_ratio = EE.KCD.evaluate_mass_accuracy_ratio(mz_landmarks, mode, mz_tolerance_ppm=10)
    # print("  Mass accuracy is estimated as %2.1f ppm." %(mass_accuracy_ratio*1000000))
    print("\n")


def get_file_masstrack_stats(infile, 
                        mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000,
                        return_sample=False):
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

        Mass accuracy was estimated on 124 matched values as -1.8 ppm.

    To-do: to add output info on instrumentation.
    '''
    new = {'sample_id': infile, 'input_file': infile, 'ion_mode': '',}
    list_mass_tracks = []
    exp = pymzml.run.Reader(infile)
    jj = 0

    for spec in exp:
        if spec.ms_level == 1:                          # MS Level 1 only
            if spec["positive scan"]:
                ionization_mode = 'pos'
                jj += 1
            else:
                ionization_mode = 'neg'

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
            'intensity': track[1], 
            } )
        ii += 1

    new['list_mass_tracks'] = list_mass_tracks
    anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(list_mass_tracks, 
                            list_mz_diff = [1.003355,], mz_tolerance_ppm=mz_tolerance_ppm)
    _mz_landmarks_ = flatten_tuplelist(anchor_mz_pairs)
    all_mz = [x['mz'] for x in list_mass_tracks]
    # down scale list_mass_tracks to verified by _mz_landmarks_
    list_mass_tracks = [list_mass_tracks[ii] for ii in _mz_landmarks_]
    peak_heights = [x['intensity'].max() for x in list_mass_tracks]
    max_peak_height = int(max(peak_heights))
    min_peak_height_ = int(min(peak_heights))
    # recommended_min_peak_height = int(0.5 * min(peak_heights))

    print("Total number of MS1 spectra: %d" %len(new['list_scan_numbers']))
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


# -----------------------------------------------------------------------------
# estimate_min_peak_height

def estimate_min_peak_height(list_input_files, 
            mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=500,
            num_files_to_use=3):
    '''
    Compute estimated min peak height. This gets min peak height from the andmark tracks in each file,
    which is the min of the mass tracks with paired 13C/12C pattern (based on m/z diff only).

    Parameters
    ----------
    list_input_files : list[str]
        input mzML filepaths, but only using num_files_to_use.
    mz_tolerance_ppm : float, optional, default: 5
        m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    min_intensity : float, optional, default: 100
        minimal intentsity value to consider, also used to filter out 0s.
    min_timepoints : int, optional, default: 5
        minimal consecutive scans to be considered real signal.
    min_peak_height : float, optional, default: 500
        a bin is not considered if the max intensity < min_peak_height.
    num_files_to_use : int, optional, default: 3
        Use randomly chosen num_files_to_use from list_input_files.

    Returns
    -------
    int, an estimated parameter for min peak_height as half of the min verified landmark peaks.
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
    Extended estimate_min_peak_height for Xasari use.

    Parameters
    ----------
    list_input_files : list[str]
        input mzML filepaths, but only using num_files_to_use.
    mz_tolerance_ppm : float, optional, default: 5
        m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    min_intensity : float, optional, default: 100
        minimal intentsity value to consider, also used to filter out 0s.
    min_timepoints : int, optional, default: 5
        minimal consecutive scans to be considered real signal.
    min_peak_height : float, optional, default: 500
        a bin is not considered if the max intensity < min_peak_height.
    num_files_to_use : int, optional, default: 3
        Use randomly chosen num_files_to_use from list_input_files.

    Returns
    -------
    A dict of ion mode and recommended min_peak_height.
    The latter is an estimated parameter for min peak_height as half of the min verified landmark peaks.

    See also
    --------
    estimate_min_peak_height
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

