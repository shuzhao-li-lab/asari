'''
Functions for subcommand `analyze`
'''
import pymzml
from .chromatograms import extract_massTracks_ 
from .experiment import ext_Experiment
from .mass_functions import *

def analyze_single_sample(infile, 
            mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000,
            parameters={}):
    '''
    Analyze single mzML file and print statistics.
    parameters are not used, just place holder to use ext_Experiment class.
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
    Example output:
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

    ionization_mode is assumed on one scan, thus not supporting polarity switch in a single file.
        to add output info on instrumentation
    '''
    new = {'sample_id': infile, 'input_file': infile, 'ion_mode': '',}
    list_mass_tracks = []
    exp = pymzml.run.Reader(infile)
    jj = 0
    for spec in exp:
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




#-------------------------------------------------------
# Alternative using pandas dataframe

# import random
# import pandas as pd
# import pymzml

def get_file_stats_dataframes(infile, max_spectra=100, return_dataframe=False):
    '''
    Quick estimate statistics based on max_spectra number.
    Example output: 
        Total number of spectra: 780
                        mz     intensity
        count  106174.000000  1.061740e+05
        mean      316.125161  1.627240e+05
        std       202.411866  2.045669e+06
        min        80.024223  1.201407e+03
        25%       159.091602  1.091480e+04
        50%       239.175484  1.943417e+04
        75%       428.384628  4.698321e+04
        max       999.849426  2.154302e+08
    '''
    mzData, intensityData = np.array([]), np.array([])
    exp = pymzml.run.Reader(infile)
    N_all_spectra = exp.get_spectrum_count()
    if max_spectra < N_all_spectra:
        chosen = random.sample(range(N_all_spectra), max_spectra)
    else:
        chosen = range(N_all_spectra)
    ii = -1   # scan_number starts with 0
    for spec in exp:
        ii += 1
        if ii in chosen:
            if spec.ms_level == 1:
                mzData = np.concatenate((mzData, spec.mz))
                intensityData = np.concatenate((intensityData, spec.i))
                
    print("Total number of spectra: %d" %N_all_spectra)
    DF = pd.DataFrame({'mz': mzData, 'intensity': intensityData})
    print(DF.describe())
    if return_dataframe:
        return DF
    