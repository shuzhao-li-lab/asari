'''
Functions related to chromatograms and mass tracks.
Use integers for RT scan numbers and intensities.
Flexible binning based on ppm accuracy.
Use mz tol (default 5 pmm) in XIC construction. 
XICs without neighbors within x ppm are considered specific (i.e. high selectivity). 
Low selectivity regions will be still inspected to determine the true number of XICs.
'''
from operator import itemgetter
import numpy as np

from scipy import interpolate
from scipy.ndimage import uniform_filter1d
from statsmodels.nonparametric.smoothers_lowess import lowess

from .mass_functions import (check_close_mzs, 
                             nn_cluster_by_mz_seeds)

INTENSITY_DATA_TYPE = np.int64
# int32 uses less memory - for large data one can check if int32 is safe, 
# i.e. under 2E9, or transform data

# -----------------------------------------------------------------------------
# mass Tracks
# -----------------------------------------------------------------------------

def extract_massTracks_(ms_expt, 
                        mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, 
                        min_peak_height=1000):
    '''
    Extract mass tracks from an object of parsed LC-MS data file.
    A mass track is an EIC for full RT range, without separating the mass traces of same m/z. 

    Parameters
    ----------
    ms_expt : pymzml.run.Reader(f)
        instance of pymzml.run.Reader(f), a parsed object of LC-MS data file
    mz_tolerance_ppm : float, optional, default: 5
        m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    min_intensity : float, optional, default: 100
        minimal intentsity value, needed because some instruments keep 0s 
    min_timepoints : int, optional, default: 5
        minimal consecutive scans to be considered real signal.
    min_peak_height : float, optional, default: 1000
        a bin is not considered if the max intensity < min_peak_height.

    Returns 
    -------
    Dict with keys 'rt_numbers', 'rt_times', and 'tracks', each with a value that is a list of
    co-indexed rt_numbers, rt_times, and mass tracks respectively. Mass tracks are represented
    as [( mz, np.array(intensities at full rt range) ), ...]
    '''
    alldata = []
    rt_times = []           # in seconds
    ii = 0
    for spec in ms_expt:
        if spec.ms_level == 1:                         # MS Level 1 only
            rt_times.append( spec.scan_time_in_minutes()*60 )
            intensities = spec.i.astype(int)
            good_positions = intensities > min_intensity
            intensities = intensities[good_positions]
            mzs = spec.mz[good_positions]
            alldata += [(mz, ii, inten) for mz, inten in zip(mzs, intensities)]
            ii += 1

    #print("extracted %d valide data points." %len(alldata))
    mzTree = {}
    for x in alldata:
        ii = int(x[0]*1000)
        if ii in mzTree:
            mzTree[ii].append(x)
        else:
            mzTree[ii] = [x]

    del alldata
    rt_numbers = list(range(len(rt_times)))
    rt_length = len(rt_numbers)

    good_bins = get_thousandth_bins(mzTree, mz_tolerance_ppm, min_timepoints, min_peak_height)
    tracks = []
    for bin in good_bins:
        tracks += bin_to_mass_tracks(bin, rt_length, mz_tolerance_ppm)

    # merge tracks if m/z overlap
    tracks.sort()
    merged, to_remove = [], []
    tracks_to_merge = check_close_mzs([x[0] for x in tracks], mz_tolerance_ppm)
    # this returns [(ii, ii-1), ...]
    for (a,b) in tracks_to_merge:
        merged.append( merge_two_mass_tracks(tracks[a], tracks[b]) )
        to_remove += [a, b]

    updated_tracks = [tracks[ii] for ii in range(len(tracks)) if ii not in to_remove] + merged

    return {
        'rt_numbers': rt_numbers,
        'rt_times': rt_times,
        'tracks': updated_tracks,
    }


def extract_single_track_fullrt_length(bin, rt_length, INTENSITY_DATA_TYPE=INTENSITY_DATA_TYPE):
    '''
    Build a mass track from a bin of data points already in limited m/z range.
    A mass track is an EIC for full RT range, without separating the mass traces. 
    Consensus m/z is taken as the mean of median m/z and the m/z of highest intensity, to be more robust.
    When multiple data points exist in the same scan (same RT), max intensity is used.

    Parameters
    ----------
    bin : list[tuple]
        data points, in format of [(mz_int, scan_num, intensity_int), ...].
    rt_length : int     
        full number of scans.
    INTENSITY_DATA_TYPE : any int numpy.dtype, optional, default: INTENSITY_DATA_TYPE (currently np.int64)
        default to np.int64. 
        Being future safe, but int32 may be adequate and more efficient.
     
    Returns
    -------
    a massTrack as ( mz, np.array(intensities at full rt range) ).
    '''
    mzs = [x[0] for x in bin]
    ints = [x[2] for x in bin]
    mz = 0.5 * (np.median(mzs) + mzs[np.argmax(ints)])  # safer option
    intensity_track = np.zeros(rt_length, dtype=INTENSITY_DATA_TYPE)
    for r in bin:                                       # this gets max intensity on the same RT scan
        intensity_track[r[1]] = max(r[2], intensity_track[r[1]])
    return ( mz, intensity_track )


def bin_to_mass_tracks(bin_data_tuples, rt_length, mz_tolerance_ppm=5):
    '''
    Construct mass tracks from data that may exceed proper m/z range (i.e. over mz_tolerance_ppm).
    A mass track is an EIC for full RT range, without separating the mass traces. 
    Imperfect ROIs will be examined in extract_massTracks_ and merged if needed.

    Parameters
    ----------
    bin_data_tuples : list[tuple]
        a flexible bin by units of 0.001 amu, in format of 
        [(mz, scan_num, intensity_int), ...]. This may or may not be within mz_tolerance_ppm.
    rt_length : int   
        full number of scans.
    mz_tolerance_ppm: float, optional, default: 5

    Returns
    -------
    a list of massTracks as [(mz, np.array(intensities at full rt range)), ...].

    See also
    --------
    extract_single_track_fullrt_length
    '''
    bin_data_tuples.sort()   # by m/z, ascending
    mz_range = bin_data_tuples[-1][0] - bin_data_tuples[0][0]
    mz_tolerance = bin_data_tuples[0][0] * 0.000001 * mz_tolerance_ppm   
    # important: double tol_ range here as mean_mz falls in single tol_
    if mz_range < mz_tolerance * 2:
        return [extract_single_track_fullrt_length(bin_data_tuples, rt_length)]
    else:
        ROIs = build_chromatogram_by_mz_clustering(bin_data_tuples, rt_length, mz_tolerance)
        return [extract_single_track_fullrt_length(R, rt_length) for R in ROIs]


def build_chromatogram_intensity_aware(bin_data_tuples, rt_length, mz_tolerance):
    '''
    Chromatogram builder as in ADAP for testing. 
    Start with highest intensity, going down by intensity and include data points within mz_tolerance.
    Repeat until no track is detected. Without requiring continuous RT, 
    which is handled in extract_single_track_fullrt_length.
    Default in asari is build_chromatogram_by_mz_clustering.

    Parameters
    ----------
    bin_data_tuples: list[tuple]
        a flexible bin in format of [(mz, scan_num, intensity_int), ...].
    rt_length: int
        number of scans (NOT USED, remove from call signature in future)
    mz_tolerance_ppm: float
        m/z tolerance in part-per-million. Used to seggregate m/z regions here.

    Returns
    -------
    assigned: 
        separated bins of [(mz, scan_num, intensity_int), ...], 
        prototype of extracted ion chromatograms to be used by extract_single_track_fullrt_length.
    '''
    bin_data_tuples.sort(key=itemgetter(2), reverse=True)
    assigned, remaining = [], bin_data_tuples
    while remaining:
        seed, tmp_remaining = [remaining[0]], []
        for T in remaining[1:]:
            if abs(T[0] - seed[0][0]) < mz_tolerance:
                seed.append(T)
            else:
                tmp_remaining.append(T)
        assigned.append( seed )
        remaining = tmp_remaining

    return assigned


def build_chromatogram_by_mz_clustering(bin_data_tuples, rt_length, mz_tolerance):
    '''
    Generates a list of prototype extracted ion chromatograms to be used 
    by extract_single_track_fullrt_length.

    Parameters
    ----------
    bin_data_tuples : list[tuple]
        a flexible bin in format of [(mz, scan_num, intensity_int), ...].
    rt_length: int
        number of scans (NOT USED, remove from call signature in future)
    mz_tolerance : float
        precomputed based on m/z and ppm, e.g. 5 ppm of 80 = 0.0004;  5 ppm of 800 = 0.0040.

    Returns
    -------
    A list of clusters as separated bins, each bin as [(mz, scan_num, intensity_int), ...]
    '''
    return nn_cluster_by_mz_seeds(bin_data_tuples, mz_tolerance, presorted=False)


def merge_two_mass_tracks(T1, T2):
    '''
    Merge two mass tracks, each massTrack as ( mz, np.array(intensities at full rt range) ).

    Parameters
    ----------
    T1 : list
        a massTrack rerpesented as ( mz, np.array(intensities at full rt range))
    T2: list
        a massTrack rerpesented as ( mz, np.array(intensities at full rt range))

    Returns
    -------
    The merged mass track. The mz value is averaged between tracks, intensities are summed pair-wise.
    '''
    return ( 0.5 * (T1[0] + T2[0]), T1[1] + T2[1] )


def get_thousandth_bins(mzTree, mz_tolerance_ppm=5, min_timepoints=5, min_peak_height=1000):
    '''
    Bin an mzTree into a list of data bins, to feed to `bin_to_mass_tracks`.
    These data bins can form a single mass track, or span larger m/z region 
    if the m/z values cannot be resolved into discrete tracks here.

    Parameters
    ----------
    mzTree: dict[list[tuples]]
        indexed data points, {thousandth_mz: [(mz, ii, intensity_int)...], ...}. 
        (all data points indexed by m/z to thousandth precision, i.e. 0.001 amu).
    mz_tolerance_ppm:  float, optional, default: 5
        m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    min_timepoints: int, optional, default: 5
        minimal consecutive scans to be considered real signal.
    min_peak_height: float, optional, default: 1000
        a bin is not considered if the max intensity < min_peak_height.
    
    Returns
    -------
    a list of flexible bins, [ [(mz, scan_num, intensity_int), ...], ... ]
    '''
    def __rough_check_consecutive_scans__(datatuples, gap_allowed=2, min_timepoints=min_timepoints):
        # check if the mass trace has at least min_timepoints consecutive scans
        # datatuples are a list of data points in format of (mz_int, scan_num, intensity_int)
        _checked = True
        check_max_len = 4 * min_timepoints                  # give longer traces a pass for now
        if len(datatuples) < check_max_len:
            min_check_val = gap_allowed + min_timepoints -1 # step is 4 in five consecutive values without gap
            rts = sorted([x[1] for x in datatuples])
            steps = [rts[ii]-rts[ii-min_timepoints+1] for ii in range(min_timepoints-1, len(rts))]
            if min(steps) > min_check_val:
                _checked = False
        return _checked

    def __check_min_peak_height__(datatuples, min_peak_height):
        return max([x[2] for x in datatuples]) >= min_peak_height

    tol_ = 0.000001 * mz_tolerance_ppm
    ks = sorted([k for k,v in mzTree.items() if len(v) >= min_timepoints]) # ascending order enforced
    bins_of_bins = []
    tmp = [ks[0]]
    for ii in range(1, len(ks)):
        _delta = ks[ii] - ks[ii-1]
        # merge adjacent bins if they are next to each other or within ppm tolerance
        if _delta==1 or _delta < tol_ * ks[ii]:
            tmp.append(ks[ii])
        else:
            bins_of_bins.append(tmp)
            tmp = [ks[ii]]

    bins_of_bins.append(tmp)
    good_bins = []
    for bin in bins_of_bins:
        datatuples = []
        for b in bin:
            datatuples += mzTree[b]
        # check the presence of min consecutive RT in small traces, to filter out more noises
        # e.g. in an example datafile, 5958 reduced to 4511 traces
        if __check_min_peak_height__(datatuples, min_peak_height) and __rough_check_consecutive_scans__(datatuples):
            good_bins.append(datatuples)
    
    # del mzTree
    return good_bins


# -----------------------------------------------------------------------------
# retention time alignment
# -----------------------------------------------------------------------------

def rt_lowess_calibration(good_landmark_peaks, 
                          selected_reference_landmark_peaks, 
                          sample_rt_numbers, reference_rt_numbers, num_iterations, sample_name, outdir):
    '''
    This is the alignment function of retention time between samples.
    Use LOWESS, Locally Weighted Scatterplot Smoothing, to create 
    correspondence between sample_rt_numbers, reference_rt_numbers.
    Predicted numbers are skipped when outside real sample boundaries.

    Parameters
    ----------
    good_landmark_peaks : list[peak]
        landmark peaks selected from this working sample.
        Landmark peaks are usually defined by 13C/12C patterns.
    selected_reference_landmark_peaks : list[peak]
        landmark peaks selected from the reference sample,
        matched and equal-length to good_landmark_peaks.
    sample_rt_numbers : list
        all scan numbers in this sample.
    reference_rt_numbers : list
        all scan numbers in the ref sample.
    sample_name : str
        sample's name

    Returns
    -------
    rt_cal_dict : dict
        dictionary converting scan number in sample_rt_numbers to calibrated integer values.
        Range matched. Only changed numbers are kept for efficiency.
    reverse_rt_cal_dict : dict
        from ref RT scan numbers to sample RT scan numbers. 
        Range matched. Only changed numbers are kept for efficiency.
    
        
    Note
    ----
        LOWESS method available in statsmodels.nonparametric.smoothers_lowess, v 0.12, 0.13+
        https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html    
        But xvals have to be forced as float array until the fix is in new release.
        See __hacked_lowess__.
    '''
    # force left and right ends, to prevent runaway curve functions
    reference_rt_bound = max(reference_rt_numbers)
    sample_rt_bound = max(sample_rt_numbers)
    rt_rightend_ = 1.1 * sample_rt_bound
    xx, yy = [-0.1 * sample_rt_bound,]*3, [-0.1 * sample_rt_bound,]*3
    rt_cal = clean_rt_calibration_points(
        [(x[0]['apex'], x[1]['apex']) for x in zip(good_landmark_peaks, selected_reference_landmark_peaks)]
    )
    xx += [L[0] for L in rt_cal] + [rt_rightend_]*3
    yy += [L[1] for L in rt_cal] + [rt_rightend_]*3
    # scale frac parameter like a sigmoid of number of data points when len(rt_cal) is in (50,150).
    FRAC = 0.6 - 0.004*(len(rt_cal)-50)
    FRAC = max(0.2, min(FRAC, 0.6))    # bound frac in (0.2, 0.6)
    
    # This requires statsmodels > v 0.12.
    # float conversion on xvals is to bypass a bug in statsmodels, which was fixed today 2022-01-27
    #lowess_predicted = lowess(yy, xx, frac= .2, it=1, xvals=np.array(sample_rt_numbers, dtype=float)) 
    #
    # downgrade now for compatibility to older statsmodels
    lowess_predicted = __hacked_lowess__(yy, xx, frac=FRAC, it=num_iterations, xvals=sample_rt_numbers)

    interf = interpolate.interp1d(lowess_predicted, sample_rt_numbers, fill_value="extrapolate")
    ref_interpolated = interf( reference_rt_numbers )
    lowess_predicted = [int(round(ii)) for ii in lowess_predicted]

    rt_cal_dict = dict( 
        [(x,y) for x,y in zip(sample_rt_numbers, lowess_predicted) if x!=y and 0<=y<=reference_rt_bound] )

    ref_interpolated = [int(round(ii)) for ii in ref_interpolated]
    reverse_rt_cal_dict = dict(
        [(x,y) for x,y in zip(reference_rt_numbers, ref_interpolated) if x!=y and 0<=y<=sample_rt_bound] )
        
    return rt_cal_dict, reverse_rt_cal_dict


def clean_rt_calibration_points(rt_cal_pairs):
    '''
    Remove redundant RT calibration data points and outliers (out of 3x stdev).
    This does not force even distribution of calibration data points.

    Parameters
    ----------
    rt_cal_pairs : list of paired scan numbers from this sample and from the reference sample.

    Returns
    -------
    rt_cal_pairs : clean and sorted version of rt_cal_pairs.
    '''
    rt_cal_pairs = [(x, x[0]-x[1]) for x in set(rt_cal_pairs)]
    _deltas = [x[1] for x in rt_cal_pairs]
    m, std3x = np.mean(_deltas), np.std(_deltas) * 3
    _low, _high = m-std3x, m+std3x

    return sorted([x[0] for x in rt_cal_pairs if _low < x[1] < _high])


def rt_lowess_calibration_debug(good_landmark_peaks, 
                          selected_reference_landmark_peaks, 
                          sample_rt_numbers, reference_rt_numbers, num_iterations, sample_name, outdir):
    '''
    This is the debug version of rt_lowess_calibration.

    Parameters
    ----------
    good_landmark_peaks : list[peak]
        landmark peaks selected from this working sample.
        Landmark peaks are usually defined by 13C/12C patterns.
    selected_reference_landmark_peaks : list[peak]
        landmark peaks selected from the reference sample,
        matched and equal-length to good_landmark_peaks.
    sample_rt_numbers : list
        all scan numbers in this sample.
    reference_rt_numbers : list
        all scan numbers in the ref sample.
    sample_name : str
        sample's name

    Returns
    -------
    rt_cal_dict : dict
        dictionary converting scan number in sample_rt_numbers to calibrated integer values.
        Range matched. Only changed numbers are kept for efficiency.
    reverse_rt_cal_dict : dict
        from ref RT scan numbers to sample RT scan numbers. 
        Range matched. Only changed numbers are kept for efficiency.
    
    Outputs
    -------
    A PDF file showing the LOWESS regression result.
        
    Note
    ----
        LOWESS method available in statsmodels.nonparametric.smoothers_lowess, v 0.12, 0.13+
        https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html    
        But xvals have to be forced as float array until the fix is in new release.
        See __hacked_lowess__.
    '''
    import matplotlib.pyplot as plt
    import os

    # force left and right ends, to prevent runaway curve functions
    reference_rt_bound = max(reference_rt_numbers)
    sample_rt_bound = max(sample_rt_numbers)
    rt_rightend_ = 1.1 * sample_rt_bound
    xx, yy = [-0.1 * sample_rt_bound,]*3, [-0.1 * sample_rt_bound,]*3
    rt_cal = clean_rt_calibration_points(
        [(x[0]['apex'], x[1]['apex']) for x in zip(good_landmark_peaks, selected_reference_landmark_peaks)]
    )
    xx += [L[0] for L in rt_cal] + [rt_rightend_]*3
    yy += [L[1] for L in rt_cal] + [rt_rightend_]*3
    # scale frac parameter like a sigmoid of number of data points when len(rt_cal) is in (50,150).
    FRAC = 0.6 - 0.004*(len(rt_cal)-50)
    FRAC = max(0.2, min(FRAC, 0.6))    # bound frac in (0.2, 0.6)

    # This requires statsmodels > v 0.12.
    # float conversion on xvals is to bypass a bug in statsmodels, which was fixed today 2022-01-27
    #lowess_predicted = lowess(yy, xx, frac= .2, it=1, xvals=np.array(sample_rt_numbers, dtype=float)) 
    #
    # downgrade now for compatibility to older statsmodels
    lowess_predicted = __hacked_lowess__(yy, xx, frac=FRAC, it=num_iterations, xvals=sample_rt_numbers)
    interf = interpolate.interp1d(lowess_predicted, sample_rt_numbers, fill_value="extrapolate")
    ref_interpolated = interf( reference_rt_numbers )
    lowess_predicted = [int(round(ii)) for ii in lowess_predicted]
    rt_cal_dict = dict( 
        [(x,y) for x,y in zip(sample_rt_numbers, lowess_predicted) if x!=y and 0<=y<=reference_rt_bound] )

    ref_interpolated = [int(round(ii)) for ii in ref_interpolated]
    reverse_rt_cal_dict = dict(
        [(x,y) for x,y in zip(reference_rt_numbers, ref_interpolated) if x!=y and 0<=y<=sample_rt_bound] )
    
    plt.figure()
    # export result figures for lowess regression
    # Create the scatter plot and line plot
    plt.scatter(xx[3:-3], yy[3:-3], label='Landmark peaks', color='black', marker='o')
    plt.plot(sample_rt_numbers, lowess_predicted, label='Alignment function', color='red', linestyle='-')

    # Customize your plot
    plt.xlabel('Sample Retention Time (sec)')
    plt.ylabel('Reference Retention Time (sec)')
    plt.title('Retention Time Alignment')
    plt.legend()
    # Save the figure as a PNG file
    plt.savefig(os.path.join(outdir, 'export', sample_name + '_rtime_alignment_result.png'))
    return rt_cal_dict, reverse_rt_cal_dict


def __hacked_lowess__(yy, xx, frac, it, xvals):
    '''
    This was workaround of a problem with Statsmodel 0.13, which was dealt by casting xvals as floats.
    Not checking possible range error.
    The bug was reported and fixed in newer Statsmodel, thus this can be phased out in future.
    '''
    lxy = lowess(yy, xx, frac, it)
    newx, newy = list(zip(*lxy))
    interf = interpolate.interp1d(newx, newy)
    return interf(xvals)


def savitzky_golay_spline(good_landmark_peaks, selected_reference_landmark_peaks, sample_rt_numbers, reference_rt_numbers):
    '''
    Placeholder.
    Modified Savitzky-Golay filter followed by spline fitting - pls follow format in rt_lowess.
    Because our data are not equally spaced, sav-gol method may produce undesired errors.
    # UnivariateSpline can't handle redundant values -
    spl = UnivariateSpline(xx, yy, )
    sample.rt_calibration_function = spl
    # rt_remap_dict will be used for index mapping to the reference sample; 
    for ii in sample.rt_numbers:
        sample.rt_remap_dict[ii] = round(spl(ii), None)
    '''
    pass

def dwt_rt_calibrate(good_landmark_peaks, selected_reference_landmark_peaks, sample_rt_numbers, reference_rt_numbers):
    '''
    Placeholder.
    Not implemented.
    '''
    pass


def remap_intensity_track(intensity_track, new, rt_cal_dict):
    '''
    Remap intensity_track based on rt_cal_dict, used by constructors.MassGrid.remap_intensity_track.

    Parameters
    ----------
    intensity_track : list[int]
        list of intensity values from a mass track.
    new : np.zeros
        new copy of np.zeros of RT length, possible longer than intensity_track,
        because samples may have different RT lengthes.
    rt_cal_dict : dict
        sample specific mapping dictionary of RT.
        dictionary converting scan number in sample_rt_numbers to calibrated integer values.

    Returns
    -------
    Updated list of intensity, using coordinates in composite mass track.
    '''
    new[ :intensity_track.shape[0]] = intensity_track.copy()
    for k,v in rt_cal_dict.items():
        new[v] = intensity_track[k]
    return new


# -----------------------------------------------------------------------------
# smoothing functions
# -----------------------------------------------------------------------------

def smooth_moving_average(list_intensity, size=9):
    '''
    Smooth data of a noisy mass track using simple moving average.

    Parameters
    ----------
    list_intensity : list[int]
        list of intensity values from a mass track.
    size : int, optional, default: 9
        window size for moving average.

    Returns
    -------
    New list of smoothed intensity values.

    
    Note
    ----
        For very noise data, one may use smooth_lowess.
    '''
    return uniform_filter1d(list_intensity, size, mode='nearest')

def smooth_lowess(list_intensity, frac=0.02):
    '''
    Smooth data of a very noisy mass track via LOWESS regression.

    Parameters
    ----------
    list_intensity : list[int]
        list of intensity values from a mass track.
    frac : float, optional, default: 0.02
        fraction of data used in LOWESS regression.

    Returns
    -------
    New list of smoothed intensity values.

        
    Note
    ----
        smooth_moving_average is preferred for most data. LOWESS is not good for small peaks.
    '''
    lxy = lowess(list_intensity, range(len(list_intensity)), frac=frac, it=1)
    _, newy = list(zip(*lxy))
    return newy
