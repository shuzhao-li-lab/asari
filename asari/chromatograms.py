'''
Use integers for RT scan numbers and intensities.
Flexible binning based on ppm accuracy.
Use mz tol (default 5 pmm) in XIC construction. 
XICs without neighbors within x ppm are considered specific (i.e. high selectivity). 
Low selectivity regions will be still inspected to determine the true number of XICs.
Leave calibration to Correspondence step.

import os
from itertools import combinations
# import matplotlib.pyplot as plt
# from sklearn.neighbors import KernelDensity
# from pyopenms import *
'''

from operator import itemgetter
import numpy as np

from scipy.signal import find_peaks 
from scipy import interpolate
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit 

from statsmodels.nonparametric.smoothers_lowess import lowess


def sum_dict(dict1, dict2):
    new = {}
    for k in dict2:
        if k in dict1:
            new[k] = dict1[k] + dict2[k]
        else:
            new[k] = dict2[k]
    return new


# -----------------------------------------------------------------------------
# peak evaluation
# -----------------------------------------------------------------------------

def gaussian_function__(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) 

def goodness_fitting__(y_orignal, y_fitted):                  # R^2 as goodness of fitting
    return 1 - (np.sum((y_fitted-y_orignal)**2) / np.sum((y_orignal-np.mean(y_orignal))**2))

def evaluate_gaussian_peak(mass_track, peak):
    '''
    Use Gaussian models to fit peaks, R^2 as goodness of fitting.
    mass_track: {'id_number': k, 'mz': mz, 'rt_scan_numbers': [..], 'intensity': [..]}
    Peak: {'parent_masstrace_id': 2812, 'mz': 359.9761889867791, 'apex': 91.0, 'height': 491241.0, 'left_base': 49.0, 'right_base': 267.0}
    return: goodness_fitting
    '''
    goodness_fitting = 0
    xx = mass_track['rt_scan_numbers'][peak['left_index']: peak['right_index']+1]
    yy = mass_track['intensity'][peak['left_index']: peak['right_index']+1]
    # set initial parameters
    a, mu, sigma =  peak['height'], peak['apex'], np.std(xx)
    try:
        popt, pcov = curve_fit(gaussian_function__, xx, yy, p0=[a, mu, sigma])
        goodness_fitting = goodness_fitting__( yy, gaussian_function__(xx, *popt))
    # failure to fit
    except (RuntimeError, ValueError):
        # about 50 occurancies on one dataset # print(peak['parent_masstrace_id'], peak['apex'])
        goodness_fitting = 0

    return goodness_fitting


# -----------------------------------------------------------------------------
# indexing function
# -----------------------------------------------------------------------------

def get_thousandth_regions(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5):
    '''
    Process all LC-MS spectral data into flexible bins by units of 0.001 amu.
    mz_tolerance_ppm: m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    ms_expt: pyopenms MSExperiment instance, loaded with LC-MS data.
    min_intensity: minimal intentsity value, needed because some instruments keep 0s 
    min_timepoints: minimal consecutive scans to be considered real signal.
    
    Return: a list of flexible bins
    '''
    def __rough_check_consecutive_scans__(datatuples, check_max_len=20, gap_allowed=2, min_timepoints=min_timepoints):
        # a list of data points in format of (mz_int, scan_num, intensity_int)
        if len(datatuples) < check_max_len:
            min_check_val = gap_allowed + min_timepoints -1 # step is 4 in five consecutive values without gap
            rts = sorted([x[1] for x in datatuples])
            steps = [rts[ii]-rts[ii-min_timepoints+1] for ii in range(min_timepoints-1, len(rts))]
            if min(steps) > min_check_val:
                return False
            else:
                return True
        else:
            return True
            
    tol_ = 0.000001 * mz_tolerance_ppm
    number_spectra = ms_expt.getNrSpectra()
    alldata = []
    for ii in range(number_spectra):
        if ms_expt[ii].getMSLevel() == 1:                               # MS Level 1 only
            for (mz, intensity) in zip(*ms_expt[ii].get_peaks()):
                if intensity > min_intensity:
                    alldata.append((mz, ii, int(intensity)))

    #print("extracted %d valide data points." %len(alldata))
    mzTree = {}
    for x in alldata:
        ii = int(x[0]*1000)
        if ii in mzTree:
            mzTree[ii].append(x)
        else:
            mzTree[ii] = [x]

    del alldata
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
        if __rough_check_consecutive_scans__(datatuples):
            good_bins.append(datatuples)
    
    # del mzTree
    return good_bins


# -----------------------------------------------------------------------------
# mass Traces
# -----------------------------------------------------------------------------

def extract_massTraces(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5):
    '''
    ms_expt: pyopenms MSExperiment instance, loaded with LC-MS data.
    return 
    rt_numbers, rt_times,
    XICs as [( mz, rtlist, intensities ), ...]
    '''
    rt_numbers = range(ms_expt.getNrSpectra())
    rt_times = [spec.getRT() for spec in ms_expt]
    good_bins = get_thousandth_regions(ms_expt, mz_tolerance_ppm, min_intensity, min_timepoints)
    xics = []
    for bin in good_bins:
        xics += bin_to_xics(bin, mz_tolerance_ppm, gap_allowed=2, min_timepoints=5)
        
    return {
        'rt_numbers': rt_numbers,
        'rt_times': rt_times,
        'xics': xics,
    }

def bin_to_xics(bin_data_tuples, mz_tolerance_ppm=5, gap_allowed=2, min_timepoints=5):
    '''
    input a flexible bin by units of 0.001 amu, in format of [(mz_int, scan_num, intensity_int), ...].
    return XICs as [( mz, rtlist, intensities ), ...]
    '''
    bin_data_tuples.sort()   # by m/z, ascending
    mz_range = bin_data_tuples[-1][0] - bin_data_tuples[0][0]
    mz_tolerance = bin_data_tuples[0][0] * 0.000001 * mz_tolerance_ppm   
    if mz_range < mz_tolerance:
        return extract_single_trace(bin_data_tuples, gap_allowed, min_timepoints)
    else:
        num_steps = int(5*mz_range/mz_tolerance)     # step in 1/5 mz_tolerance
        hist, bin_edges = np.histogram([x[0] for x in bin_data_tuples], num_steps)
        # example hist: array([  8,   33,  11,  24,  31,  50,  81, 164, 269,  28,   7])
        hist_starts = [0] + [hist[:ii].sum() for ii in range(1,num_steps+1)]
        # find_peaks returns e.g. array([ 1, 8]), {}. distance=3 because it's edge of tolerance
        peaks, _ = find_peaks(hist, distance=3)
        if peaks.any():
            XICs = []
            for p in peaks:
                left = max(0, p-2)
                right = min(p+3, num_steps)
                XICs += extract_single_trace(bin_data_tuples[hist_starts[left]: hist_starts[right]], gap_allowed, min_timepoints)
            return XICs
        else:
            peak = np.argmax(hist)  # if find_peaks fails, get position of max value
            left = max(0, peak-2)
            right = min(peak+3, num_steps)
            return extract_single_trace(bin_data_tuples[hist_starts[left]: hist_starts[right]], gap_allowed, min_timepoints)

def extract_single_trace(bin, gap_allowed=2, min_timepoints=5):
    mz = np.mean([x[0] for x in bin])
    bin.sort(key=itemgetter(1))   # sort by increasing RT (scan number)
    traces = []
    tmp = [bin[0]]
    for ii in range(1,len(bin)):
        _d = bin[ii][1] - bin[ii-1][1]
        if _d > gap_allowed:
            traces.append(tmp)
            tmp = [bin[ii]]
        else:
            tmp.append(bin[ii])
    traces.append(tmp)
    traces = [t for t in traces if len(t) >= min_timepoints]
    # filtered by consecutive RT scans. Now integrate intensity per RT
    XICs = []
    for tt in traces:
        _d = {}
        for r in tt: 
            _d[r[1]] = 0
        for r in tt:     # this gets max intensity on the same RT scan
            _d[r[1]] = max(r[2], _d[r[1]])
        rtlist = sorted(list(_d))
        intensities = [_d[x] for x in rtlist]
        XICs.append(( mz, rtlist, intensities ))
    
    return XICs


# -----------------------------------------------------------------------------
# mass Tracks
# -----------------------------------------------------------------------------

def extract_massTracks_(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5):
    '''
    A mass track is an EIC for full RT range, without separating the mass traces. 
    ms_expt: pyopenms MSExperiment instance, loaded with LC-MS data.
    return 
    rt_numbers, rt_times,
    tracks as [( mz, rtlist, intensities ), ...]

    To-do: check how rtlist becomes float numbers while they should be integers.
    '''
    rt_numbers = range(ms_expt.getNrSpectra())
    rt_times = [spec.getRT() for spec in ms_expt]
    good_bins = get_thousandth_regions(ms_expt, mz_tolerance_ppm, min_intensity, min_timepoints)
    tracks = []
    for bin in good_bins:
        tracks += bin_to_mass_tracks(bin, mz_tolerance_ppm)
        
    return {
        'rt_numbers': rt_numbers,
        'rt_times': rt_times,
        'tracks': tracks,
    }

def extract_single_track_(bin):
    '''
    A mass track is an EIC for full RT range, without separating the mass traces. 
    input bins in format of [(mz_int, scan_num, intensity_int), ...].
    RT gaps are filled by zeros to create one continuous trace, to be used for later peak detection.
    return a massTrack as ( mz, rtlist, intensities ), using continuous RT numbers.
    '''
    mz = np.mean([x[0] for x in bin])
    rtlist = [x[1] for x in bin]
    min_rt, max_rt = min(rtlist), max(rtlist)
    rtlist = range(min_rt, max_rt+1)    # filling gaps of RT this way if needed, and sorted
    # bin.sort(key=itemgetter(1))   # sort by increasing RT (scan number), not needed any more
    _d = {}                             # dict to hold rt to intensity mapping
    for ii in rtlist: 
        _d[ii] = 0
    for r in bin:                       # this gets max intensity on the same RT scan
        _d[r[1]] = max(r[2], _d[r[1]])
    intensities = [_d[x] for x in rtlist]
    return ( mz, rtlist, intensities )

def bin_to_mass_tracks(bin_data_tuples, mz_tolerance_ppm=5):
    '''
    input a flexible bin by units of 0.001 amu, in format of [(mz_int, scan_num, intensity_int), ...].
    A mass track is an EIC for full RT range, without separating the mass traces. 
    return massTracks as [( mz, rtlist, intensities ), ...]
    '''
    bin_data_tuples.sort()   # by m/z, ascending
    mz_range = bin_data_tuples[-1][0] - bin_data_tuples[0][0]
    mz_tolerance = bin_data_tuples[0][0] * 0.000001 * mz_tolerance_ppm   
    if mz_range < mz_tolerance:
        return [extract_single_track_(bin_data_tuples)]
    else:
        num_steps = int(5*mz_range/mz_tolerance)     # step in 1/5 mz_tolerance
        hist, bin_edges = np.histogram([x[0] for x in bin_data_tuples], num_steps)
        # example hist: array([  8,   33,  11,  24,  31,  50,  81, 164, 269,  28,   7])
        hist_starts = [0] + [hist[:ii].sum() for ii in range(1,num_steps+1)]
        # find_peaks returns e.g. array([ 1, 8]), {}. distance=3 because it's edge of tolerance
        peaks, _ = find_peaks(hist, distance=3)
        if peaks.any():
            tracks = []
            for p in peaks:
                left = max(0, p-2)
                right = min(p+3, num_steps)
                tracks.append( extract_single_track_(bin_data_tuples[hist_starts[left]: hist_starts[right]]) )
            return tracks
        else:
            peak = np.argmax(hist)  # if find_peaks fails, get position of max value
            left = max(0, peak-2)
            right = min(peak+3, num_steps)
            return [extract_single_track_(bin_data_tuples[hist_starts[left]: hist_starts[right]])]



# -----------------------------------------------------------------------------
# retention time alignment
# -----------------------------------------------------------------------------

def rt_lowess_calibration(good_landmark_peaks, selected_reference_landmark_peaks, full_rt_range):
    '''
    Use LOWESS, Locally Weighted Scatterplot Smoothing
    checked available in statsmodels.nonparametric.smoothers_lowess, v 0.12, 0.13+
        https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html    
        But xvals have to be forced as float array.

    Input
    =====
    good_landmark_peaks and selected_reference_landmark_peaks are equal-length lists of peaks,
    selected from working sample and reference sample as landmarks for RT alignment. 
    full_rt_range: all scan numbers in this sample.

    Return
    ======
    rt_remap_dict, dictionary converting scan number in full_rt_range to calibrated integer values.
                    so no need to recalculate every time for each track.
    '''
    # force left and right ends, to prevent runaway curve functions
    rt_rightend_ = 1.1 * max(full_rt_range)
    xx, yy = [0, 0, 0, ], [0, 0, 0, ]
    rt_cal = sorted([(x[0]['apex'], x[1]['apex']) for x in zip(good_landmark_peaks, selected_reference_landmark_peaks)])

    xx += [L[0] for L in rt_cal] + [rt_rightend_]*3
    yy += [L[1] for L in rt_cal] + [rt_rightend_]*3
    # float conversion on xvals is to bypass a bug in statsmodels, which was fixed today 2022-01-27
    lowess_predicted = lowess(yy, xx, frac= .2, it=1, xvals=np.array(full_rt_range, dtype=float)) 

    return dict(zip( full_rt_range, [round(ii,ndigits=None) for ii in lowess_predicted] ))

def __hacked_lowess__(yy, xx, frac, it, xvals):
    '''
    This was workaround of a problem with Statsmodel 0.13, which was dealt by casting xvals as floats.
    Not checking possible range error
    '''
    lxy = lowess(yy, xx, frac, it)
    newx, newy = list(zip(*lxy))
    interf = interpolate.interp1d(newx, newy)
    return interf(xvals)


def savitzky_golay_spline(good_landmark_peaks, selected_reference_landmark_peaks, full_rt_range):
    '''
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

def dwt_rt_calibrate(good_landmark_peaks, selected_reference_landmark_peaks, full_rt_range):
    '''Other RT alignment method to include - pls follow format in rt_lowess.
    not implemented.
    '''
    pass


# -----------------------------------------------------------------------------
# smoothing function
# -----------------------------------------------------------------------------

def smooth_moving_average(list_intensity, size=9):
    '''
    Smooth data for very noisy mass tracks.
    Using simple moving average here; LOWESS does not look good for this.
    '''
    return uniform_filter1d(list_intensity, size, mode='nearest')
