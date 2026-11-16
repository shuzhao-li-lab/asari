"""
Functions for subcommand `peek` (renamed from `analyze` in asari 2.0)

Mass accuracy is estimated against reference m/z values in db/calibration_mz.json
(KEGG and CSM based, {"pos": [...], "neg": [...]}). No other database is used.
A plot of scan sequence (tools.plot_scan_seq.plot_msn) is exported if matplotlib is available.
"""

import os
import json
import time
import random
import bisect
import importlib.resources as pkg_resources
import pymzml
import numpy as np

from mass2chem.search import find_mzdiff_pairs_from_masstracks

from . import db
from .chromatograms import extract_massTracks_
from .mass_functions import flatten_tuplelist
from .utils import bulk_process

DEFAULT_CALIBRATION_DB = 'calibration_mz.json'
# scans kept with full peaks for plot_msn, matching its nspec_plot and offset_n
NSPEC_PLOT, OFFSET_N = 50, 100


def peek_single_sample(infile, parameters={}):
    """
    Peek into a single mzML file and print statistics. Used by asari subcommand `peek`.
    Mass accuracy is estimated by matching landmark m/z values (12C/13C pairs) to
    db/calibration_mz.json (or parameters['db'] if given, same format).
    A PDF plot of MS1/MS2 scan sequence is written to a time-stamped subdirectory under parameters['outdir'],
    <project_name>_peek_<timestamp>/, if matplotlib is available.

    Parameters
    ----------
    infile : str
        input mzML filepath.
    parameters : dict
        asari parameters; uses mz_tolerance_ppm, min_intensity_threshold, min_timepoints,
        min_peak_height, outdir and db.
    """
    print("Peek into %s\n" %infile)
    mz_landmarks, mode, min_peak_height_, plot_data = get_file_masstrack_stats(
        infile, parameters, return_sample=False, return_plot_data=True)

    calibration_file = parameters.get('db') or str(
        pkg_resources.files(db).joinpath(DEFAULT_CALIBRATION_DB))
    reference_mzs = load_calibration_mz(calibration_file, mode)
    mass_accuracy_ratio = evaluate_mass_accuracy_ratio(mz_landmarks, reference_mzs, mz_tolerance_ppm=10)
    if mass_accuracy_ratio is not None:
        print("  Mass accuracy is estimated as %2.1f ppm." %(mass_accuracy_ratio*1000000))
    print("\n")

    # time-stamped subdirectory, as in other subcommands; created only when the plot is written
    time_stamp = ''.join([str(x) for x in time.localtime()[1:6]])
    subdir = '_'.join([parameters.get('project_name', 'asari_project'), 'peek', time_stamp])
    outdir = os.path.join(parameters.get('outdir') or '', subdir)
    export_scan_plot(infile, plot_data, outdir)


def load_calibration_mz(calibration_file, mode='pos'):
    """Read reference m/z values of ionization mode from calibration JSON, {"pos": [...], "neg": [...]}."""
    with open(calibration_file, encoding='utf-8') as f:
        mzs = json.load(f).get(mode, [])
    return sorted(mzs)


def evaluate_mass_accuracy_ratio(query_mz_list, reference_mzs, mz_tolerance_ppm=10):
    """
    Match each query m/z to the nearest value in sorted reference_mzs within mz_tolerance_ppm,
    then estimate systematic shift as mean of (query - reference)/reference.

    Returns ratio, or None if no match.
    """
    results = []
    for mz in query_mz_list:
        ii = bisect.bisect_left(reference_mzs, mz)
        nearest = [reference_mzs[jj] for jj in (ii-1, ii) if 0 <= jj < len(reference_mzs)]
        if nearest:
            ref = min(nearest, key=lambda x: abs(x - mz))
            if abs(mz - ref) <= ref * mz_tolerance_ppm * 1e-6:
                results.append(mz/ref - 1)
    if results:
        ratio = float(np.mean(results))
        print("Mass accuracy was estimated on %d matched values as %2.1f ppm." %(len(results), ratio*1000000))
        return ratio
    else:
        print("No m/z match found in calibration database to estimate mass accuracy.")
        return None


def export_scan_plot(infile, plot_data, outdir=None):
    """
    Plot MS1/MS2 scan sequence via tools.plot_scan_seq.plot_msn, if matplotlib is available.
    Output is <outdir>/<file name>_peek.pdf.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        from .tools.plot_scan_seq import plot_msn
    except ImportError:
        print("matplotlib is not available, skipping scan plot.")
        return None
    ms1_spectra, ms2_spectra = plot_data
    if not ms1_spectra and not ms2_spectra:
        return None
    outdir = outdir or '.'
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, os.path.splitext(os.path.basename(infile))[0] + '_peek.pdf')
    title = os.path.basename(infile)
    plot_msn(ms1_spectra, ms2_spectra, nspec_plot=NSPEC_PLOT, offset_n=OFFSET_N, title=title, outfile=outfile)
    print("Scan plot was written to %s\n" %outfile)
    return outfile


def get_file_masstrack_stats(infile, parameters, return_sample=False, return_plot_data=False):
    """
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
        Else, return landmark_mzs (m/z values of landmark tracks), ionization_mode, ``min_peak_height_``.
    return_plot_data : bool, optional, default: False
        if True (and not return_sample), also return (ms1_spectra, ms2_spectra) for plot_msn.
        Full peaks are kept only for the first NSPEC_PLOT + OFFSET_N scans, which plot_msn may use;
        other scans keep only their max intensity.


    Note
    ----
    Example output::

        Total number of MS1 spectra: 741
        of which 0 are positive ionization mode.

        Assuming ionization mode is neg.

        Maximum retention time (sec): 299.818228
        m/z range: (min 80.011578, median 358.010062, max 997.616794)

        Found 14063 mass tracks.
        Found 4054 12C/13C isotopic pairs as landmarks.
        Max intensity in any landmark track:  687,714,048
        Minimal height of landmark tracks:  2,334

        Mass accuracy was estimated on 203 matched values as -0.4 ppm.


    To-do: to add output info on instrumentation? (Not sure how pymzml handles the info if at all.)
        To merge from ``chromatograms.extract_massTracks_``.

    """
    mz_tolerance_ppm = parameters['mz_tolerance_ppm']
    min_intensity = parameters['min_intensity_threshold']
    min_timepoints = parameters['min_timepoints']
    min_peak_height = parameters['min_peak_height']

    new = {'sample_id': infile, 'input_file': infile, 'ion_mode': '',}
    list_mass_tracks = []

    jj, L1, L2 = 0, 0, 0
    counter = 0
    ms1_spectra, ms2_spectra = [], []      # for plot_msn, ids as 'sp' + sequential number
    with pymzml.run.Reader(infile) as exp:
        for spec in exp:
            counter += 1
            if spec.ms_level == 1:
                L1 += 1
                if spec["positive scan"]:
                    ionization_mode = 'pos'
                    jj += 1
                else:
                    ionization_mode = 'neg'
            elif spec.ms_level == 2:
                L2 += 1
            if return_plot_data and spec.ms_level in (1, 2):
                peaks = [(mz, intensity) for mz, intensity in spec.peaks('centroided') if intensity > 0]
                if counter > NSPEC_PLOT + OFFSET_N:
                    peaks = [(0, max([x[1] for x in peaks]))] if peaks else []
                record = {'id': 'sp'+str(counter), 'rtime': spec.scan_time_in_minutes() * 60, 'peaks': peaks}
                if spec.ms_level == 1:
                    ms1_spectra.append(record)
                elif peaks:
                    ms2_spectra.append(record)

    xdict = extract_massTracks_(infile,                 # MS level 1
                mz_tolerance_ppm=mz_tolerance_ppm,
                min_intensity=min_intensity,
                min_timepoints=min_timepoints,
                min_peak_height=min_peak_height)
    new['list_scan_numbers'] = xdict['rt_numbers']        # list of scans, starting from 0
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

    print("Total number of MS2 spectra: %d" %L2)

    print("Maximum retention time (sec): %f" %max(new['list_retention_time']))
    print("m/z range: (min %f, median %f, max %f)\n" %(np.min(all_mz), np.median(all_mz), np.max(all_mz)))

    print("Found %d mass tracks." %ii)
    print("Found %d 12C/13C isotopic pairs as landmarks." %len(anchor_mz_pairs))
    print("Max intensity in any landmark track: ", f"{max_peak_height:,}")
    print("Minimal height of landmark tracks: ", f"{min_peak_height_:,}", "\n")

    landmark_mzs = [x['mz'] for x in list_mass_tracks]      # _mz_landmarks_ are index numbers
    if return_sample:
        return new
    elif return_plot_data:
        return landmark_mzs, ionization_mode, min_peak_height_, (ms1_spectra, ms2_spectra)
    else:
        return landmark_mzs, ionization_mode, min_peak_height_


def match_mzdiff_pairs_by_rt(matched_mz_pairs,
            list_mass_tracks, max_scan_number, rt_window_perc=0.05, min_scans_window=10):
    """
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
    """
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
    """
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
    """

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

