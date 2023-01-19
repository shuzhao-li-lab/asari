'''
Functions for elution peak detection and evaluation.
stats_detect_elution_peaks and deep_detect_elution_peaks use different supporting functions,
as the former deals with conversion of RT coordinates.
'''

import multiprocessing as mp
import numpy as np
from scipy.signal import detrend, find_peaks 
from scipy.optimize import curve_fit 

from .chromatograms import *

# -----------------------------------------------------------------------------
# multicore processing for peak detection
# -----------------------------------------------------------------------------

def batch_deep_detect_elution_peaks(list_mass_tracks, number_of_scans, parameters):
    with mp.Manager() as manager:
        shared_list = manager.list()
        iters = iter_peak_detection_parameters(list_mass_tracks, number_of_scans, parameters, shared_list)
        with mp.Pool( parameters['multicores'] ) as pool:
            # call peak detection function here
            pool.starmap( stats_detect_elution_peaks, iters )

        FeatureList = list(shared_list)
    return FeatureList

def iter_peak_detection_parameters(list_mass_tracks, number_of_scans, parameters, shared_list):
    '''
    Prominence requirement is critical to peak detection. 
    This can be dynamically determined based on parameters['min_prominence_threshold'] 
    and noise level in a region.
    In current default method, iteration is not used.
    '''
    iters = []
    min_peak_height = parameters['min_peak_height']
    min_prominence_threshold = parameters['min_prominence_threshold']
    min_intensity_threshold = parameters['min_intensity_threshold']
    min_fwhm = round( 0.5 * parameters['min_timepoints'] )
    snr = parameters['signal_noise_ratio']
    peakshape = parameters['gaussian_shape']
    wlen = parameters['wlen']            # divided window to the left and to right
    min_prominence_ratio = 0.05
    iteration = False                    # a 2nd round of peak detection if enough remaining datapoints
    for mass_track in list_mass_tracks:
        if mass_track['intensity'].max() > min_peak_height:
            iters.append(
                (mass_track, number_of_scans, min_peak_height, min_fwhm, min_prominence_threshold,
                wlen, snr, peakshape, min_prominence_ratio, iteration, min_intensity_threshold, 
                shared_list)
            )
    return iters


# -----------------------------------------------------------------------------
# Statistics guided peak detection
# -----------------------------------------------------------------------------

def stats_detect_elution_peaks(mass_track, number_of_scans, 
                min_peak_height, min_fwhm, min_prominence_threshold,
                wlen, snr, peakshape, min_prominence_ratio, iteration, min_intensity_threshold,
                shared_list):
    '''
    Stats guided peak detection. 
    The mass_track is first cleaned: rescale, detrend, smooth and subtract baseline if needed.
    Implicitly, track wide noise_level is estimated on stdev of bottom quartile values in cleaned mass track.
    ROIs are then separated by noise_level.
    peakshape is calculated on cleaned mass track.
    Final SNR is computed on local noise (average of up to 100 nonpeak data points on each side of a peak).

    The advantages of separating ROI are 1) better performance of long LC runs, and
    2) big peaks are less likely to shallow over small peaks. 
    iteration: if True, a 2nd round of peak detection if enough remaining datapoints. Not used now.
    Gap allowed for 2 scan in constructing ROIs. 
    Actually in composite mass track - gaps should not exist after combining many samples.
    Now position indices have to be converted internally in this function.

    Input
    =====
    mass_track: {'id_number': k, 'mz': mz, 'rt_scan_numbers': [..], 'intensity': [..]}
                Mass tracks are expected to be continuous per m/z value, with 0s for gaps.
    snr: minimal signal to noise ratio required for a peak. 
                Noise is defined by the mean of all non-peak data points in mass_track.
    min_peak_height, min_fwhm, min_prominence_threshold (default min_peak_height/3) as in main parameters.
    wlen: window size for evaluating prominence in peaks. Important to resolve clustered narrow peaks.
    
    min_prominence_ratio: require ratio of prominence relative to peak height. Not used now.
    reverse_detection: use reversed intensity array during peak detection. Not used now.
    iteration, not used in this function but possible to add iteration on high intensity ROIs 
                if peaks are not detected 1st round.

    Update
    ======
    shared_list: list of peaks in JSON format, to pool with batch_deep_detect_elution_peaks.
    '''
    list_json_peaks, list_peaks = [], []
    list_scans = np.arange(number_of_scans)
    _baseline_, noise_level, scaling_factor, min_prominence_threshold, list_intensity = audit_mass_track(
                mass_track['intensity'], min_fwhm, min_intensity_threshold, min_peak_height, min_prominence_threshold
                )       
    # # get ROIs by separation/filtering with noise_level, allowing 2 gap
    __selected_scans__ = list_scans[list_intensity > noise_level]
    if __selected_scans__.any():
        ROIs = []
        tmp = [__selected_scans__[0]]
        for ii in __selected_scans__[1:]:
            if ii - tmp[-1] < 3:
                tmp += range(tmp[-1]+1, ii+1)
            else:
                ROIs.append(tmp)
                tmp = [ii]

        ROIs.append(tmp)
        ROIs = [r for r in ROIs if len(r) >= min_fwhm + 2]
        if ROIs:
            for R in ROIs:
                if len(R) < 3 * min_fwhm:       # extend if too short - help narrow peaks
                    R = extend_ROI(R, number_of_scans)
                list_intensity_roi = list_intensity[R]
                peaks = detect_evaluate_peaks_on_roi( 
                            list_intensity_roi, R, 
                            min_peak_height, min_fwhm, min_prominence_threshold, wlen,
                            snr, peakshape, min_prominence_ratio, noise_level
                    )
                list_json_peaks += peaks

            # evaluate and format peaks
            list_cSelectivity, set_peak_datapoints = __peaks_cSelectivity_stats_(list_intensity, list_json_peaks)
            noise_data_points = [x for x in list_scans if x not in set_peak_datapoints]
            for ii in range(len(list_json_peaks)):
                list_json_peaks[ii]['cSelectivity'] = list_cSelectivity[ii]

            for peak in list_json_peaks:
                peak['parent_masstrack_id'] = mass_track['id_number']
                peak['mz'] = mass_track['mz']
                peak['height'] = scaling_factor * peak['height'] + _baseline_
                __height = min(peak['height'], 99999999)        # cap upper limit and avoid INF 
                peak['snr'] = int(__height/compute_noise_by_flanks(
                                    peak, mass_track['intensity'], noise_data_points, 
                                    min_intensity_threshold, noise_level)) 
                if peak['snr'] >= snr:
                    peak['height'] = int(peak['height'])
                    list_peaks.append(peak)

    shared_list += list_peaks

def compute_noise_by_flanks(peak, list_intensity, noise_data_points, min_intensity_threshold, old_noise_level):
    '''Compute SNR by the adjacent nonpeak data points.
    '''
    padding, N = 30, 100          # number of half data points to use
    noise_data_points = sorted(noise_data_points + [peak['apex']])
    center = noise_data_points.index(peak['apex'])
    use_index = noise_data_points[center-padding-N: center-padding] + \
                noise_data_points[center+padding: center+N+padding]
    if use_index:
        return max(min_intensity_threshold, (list_intensity[use_index]+1).mean())
    else:
        return old_noise_level

def audit_mass_track(list_intensity, min_fwhm, min_intensity_threshold, min_peak_height, min_prominence_threshold):
    '''Get stats on a mass track (list_intensity), rescale, detrend, smooth and subtract baseline if needed.
    If noise_level is higher than 1% of max value in cleaned track, or track is low intensity, 
        smooth_moving_average is applied.
    noise_level : modeled noise level, as stdev of bottom 25% values or arbitrary min as 100.
    If noise is high, min_prominence_threshold is matched to noise level.
    returns    _baseline_, noise_level, scaling_factor, new_prominence, clean_list_intensity 
    '''
    scaling_factor, LOW, HIGH = 1, min_intensity_threshold, 1E8
    _baseline_, noise_level = LOW, LOW                     # will not change on a clean track
    max_intensity, median_intensity = list_intensity.max(), np.median(list_intensity)
    if max_intensity > HIGH:
        scaling_factor = max_intensity/HIGH
        list_intensity = list_intensity/scaling_factor

    if median_intensity > LOW: 
        LL = list_intensity[list_intensity > LOW]
        if len(LL) > len(list_intensity) * 0.5 and median_intensity > 10 * min_peak_height:
            list_intensity = detrend(list_intensity)        # detrend
        bottom_x_perc = list_intensity[list_intensity < LOW + np.quantile(list_intensity, 0.25)]
        _baseline_, noise_level = bottom_x_perc.mean(), bottom_x_perc.std()
        
    _baseline_, noise_level = max(_baseline_, LOW), max(noise_level, LOW)
    new_prominence = max(noise_level, min_prominence_threshold)
    # decision on smoothing
    if 100 * noise_level > max_intensity or max_intensity < 10 * min_peak_height:        
        list_intensity = smooth_moving_average(list_intensity, size=min_fwhm + 2)

    list_intensity = list_intensity - _baseline_

    return _baseline_, noise_level, scaling_factor, new_prominence, list_intensity 

def detect_evaluate_peaks_on_roi(list_intensity_roi, rt_numbers_roi, 
                    min_peak_height, min_fwhm, min_prominence_threshold, wlen,
                    snr, peakshape, min_prominence_ratio,
                    noise_level,
                    ):
    '''
    Return list of peaks based on detection in ROI, 
        defined by (list_intensity_roi, rt_numbers_roi) after mass track is cleaned.
    Min peak height is check by (noise_level * snr).
    Find peaks first, then evaluate their SNR and peakshape.

    min_prominence_ratio : Prominence required as a ratio of peak height. 
        Top intensity is specific to each ROI, but there can be small peaks when ROI is not cleanly separated.
    find_peak is from scipy.signal.find_peaks;
        threshold is not meaningful here as it refers to amplitude difference to other peaks.
    '''
    list_peaks = []
    # raise prominence if high intensity and high noise
    max_intensity = list_intensity_roi.max()
    if max_intensity > 10 * min_peak_height and 10 * noise_level > min_peak_height:
        min_prominence_threshold = max(min_prominence_threshold, max_intensity*min_prominence_ratio)
    peaks, properties = find_peaks(list_intensity_roi, 
                                    height=min_peak_height, 
                                    distance=min_fwhm,
                                    prominence=min_prominence_threshold,
                                    width=min_fwhm, 
                                    wlen=wlen,
                                    ) 
    for ii in range(peaks.size):
        if properties['right_bases'][ii] - properties['left_bases'][ii] >= min_fwhm + 2:
            _jpeak = evaluate_roi_peak_json_(
                ii, list_intensity_roi, rt_numbers_roi, peaks, properties, peakshape, min_fwhm
                )
            if _jpeak:
                list_peaks.append(_jpeak)

    return check_overlap_peaks(list_peaks)

def evaluate_roi_peak_json_( ii, list_intensity_roi, rt_numbers_roi, peaks, properties, peakshape, min_fwhm):
    '''
    Return the ii-th peak in peaks with basic properties assigned in dict.
    peaks, properties as from scipy find_peaks; list_intensity_roi, rt_numbers_roi from mass_track.
    Peak apex position is peaks[ii].
    This handles the conversion btw indices of ROI and indices of mass_track.
    The peak shape is evluated here on a Gaussian model, 
    and the left, right bases are constrained by N*stdev in the fitted model.
    This avoids long tails from scipy.signal.find_peaks.

    '''
    left_index, right_index = properties['left_bases'][ii], properties['right_bases'][ii]   # index positions on ROI
    goodness_fitting, sigma = evaluate_gaussian_peak_on_intensity_list(list_intensity_roi,
                properties['peak_heights'][ii], peaks[ii], left_index, right_index)
    if goodness_fitting > peakshape:
        _halfwidth = int(max(abs(sigma), min_fwhm) * 6)     # sigma needs to be greater than min_fwhm
        _L, _R = peaks[ii] - _halfwidth, peaks[ii] + _halfwidth
        left_index, right_index = max(left_index, _L), min(right_index, _R)
        left_base, right_base = rt_numbers_roi[left_index], rt_numbers_roi[right_index]
        peak_area = int( list_intensity_roi[left_index: right_index+1].sum() ) 
        return {
                'apex': rt_numbers_roi[peaks[ii]], 
                'peak_area': peak_area,
                'height': int(properties['peak_heights'][ii]),
                'left_base': left_base,
                'right_base': right_base, 
                'goodness_fitting': goodness_fitting,
        }
    else:
        return None

def reverse_reverse_peak_detection(list_json_peaks, max_rt_index):
    '''
    Recover json peaks when reverse_peak_detection is used.
    '''
    LL = []
    for J in list_json_peaks:
        old_left, old_right = J['left_base'], J['right_base']
        J['right_base'] = max_rt_index - old_left
        J['left_base'] = max_rt_index - old_right
        J['apex'] = max_rt_index - J['apex']
        LL.append(J)
    return LL


def __peaks_cSelectivity_stats_(__list_intensity, _jpeaks):
    '''
    peaks, properties as from find_peaks; 
    __list_intensity is np.array of the mass_track.
    See also __peaks_cSelectivity__.
    return list_cSelectivity, set_peak_datapoints
    '''
    set_peak_datapoints = []
    for peak in _jpeaks:
            set_peak_datapoints += list(range(peak['left_base'], peak['right_base']+1))

    set_peak_datapoints = set(set_peak_datapoints)
    _peak_datapoints = __list_intensity[list(set_peak_datapoints)]    # peaks may overlap
    list_cSelectivity = []
    for ii in range(len(_jpeaks)):
        _threshold = 0.5 * _jpeaks[ii]['height']
        _peak_datapoints_level = _peak_datapoints[_peak_datapoints > _threshold].size
        _background_level = __list_intensity[__list_intensity > _threshold].size
        if _background_level >= _peak_datapoints_level > 0:                  # to rule out artifact from smooth_lowess
            list_cSelectivity.append(_peak_datapoints_level / _background_level)
        else:
            list_cSelectivity.append( 0 )
            
    return list_cSelectivity, set_peak_datapoints


# -----------------------------------------------------------------------------
# Generic functions for detection and evaluation
# -----------------------------------------------------------------------------

def gaussian_function__(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) 

def goodness_fitting__(y_orignal, y_fitted):                  # R^2 as goodness of fitting
    return 1 - (np.sum((y_fitted-y_orignal)**2) / np.sum((y_orignal-np.mean(y_orignal))**2))

def evaluate_gaussian_peak_on_intensity_list(intensity_list, height, apex, left, right):
    '''
    Use Gaussian models to fit peaks, R^2 as goodness of fitting.
    Peak shapes may be more Voigt or bigaussian, but the impact on fitting result is negligible.
    The parameters height, apex, left, right are relevant to intensity_list.
    When intensity_list is of ROI, these parameters are not referring to full mass track.
    return: goodness_fitting, fitted_sigma
    '''
    goodness_fitting = 0
    xx = range(left, right+1, 1)
    yy = intensity_list[xx]            # consider re-weighted by np.sqrt for very high data points?
    a, mu, sigma =  height, apex, np.std(xx)    # set initial parameters
    try:
        popt, pcov = curve_fit(gaussian_function__, xx, yy, p0=[a, mu, sigma])
        a, mu, sigma = popt
        goodness_fitting = goodness_fitting__( yy, gaussian_function__(xx, a, mu, sigma))
    # failure to fit
    except (RuntimeError, ValueError, TypeError):
        # fitting errors  # print(peak['mz'], peak['apex'], yy)
        goodness_fitting = 0

    return goodness_fitting, sigma


def lowess_smooth_track(list_intensity, number_of_scans):
    '''
    For noisy data, smooth_lowess is applied before peak detection. 
    The other method, smooth_moving_average, does not work well for very noisy data.
    Smoothing will reduce the height of narrow peaks in CMAP, but not on the reported values,  
    because peak area is extracted from each sample after. 
    The likely slight expansion of peak bases can add to robustness.
    '''
    _frac_ = 0.05
    if number_of_scans < 200:
        _frac_ = min(10.0/number_of_scans, 0.8)
    return np.array( smooth_lowess(list_intensity, frac=_frac_) )


def quick_detect_unique_elution_peak(intensity_track, 
                            min_peak_height=100000, min_fwhm=3, min_prominence_threshold_ratio=0.2):
    '''
    Quick peak detection, only looking for a high peak with high prominence.
    This can be used for quick check on good peaks, or selecting landmarks for alignment purposes.
    '''
    max_intensity = intensity_track.max()
    prominence = min_prominence_threshold_ratio * max_intensity
    unique_peak = None
    if max_intensity > min_peak_height:
        peaks, properties = find_peaks(intensity_track, height=min_peak_height, width=min_fwhm, 
                                                        prominence=prominence) 
        if peaks.size == 1:
            unique_peak = {
                'apex': peaks[0], 
                'height': properties['peak_heights'][0], # not used now
            }
    return unique_peak

def check_overlap_peaks(list_peaks):
    '''
    Check overlap btw a list of JSON peaks. 
    list_peaks are already ordered by RT from find_peaks. Overlap usually from splitting.
    Return unique peaks.
    '''
    if len(list_peaks) < 2:
        return list_peaks
    else:
        clusters = []
        tmp = [list_peaks[0]]
        for peak in list_peaks[1:]:
            if _check_overlap(peak, tmp[-1]):
                tmp.append(peak)
            else:
                clusters.append(tmp)
                tmp = [peak]
        clusters.append(tmp)
        new = []
        for C in clusters:
            new += cleanup_peak_cluster(C)
            
        return new

def _check_overlap(peak1, peak2):
    '''
    Check if overlap btw (peak1, peak2) exceeds 3 scans.
    '''
    tuple1, tuple2 = (peak1['left_base'], peak1['right_base']), (peak2['left_base'], peak2['right_base'])
    left_peak, right_peak = tuple1, tuple2
    if tuple1[1] > tuple2[1]:
        left_peak, right_peak = tuple2, tuple1
    overlap = max(0, left_peak[1] - right_peak[0])
    if overlap > 3: 
        return True
    else:
        return False

def _merge_peak_cluster(cluster_peaks):
    '''
    Merge overlap reported peaks that stem from unclear separation.
    If two true peaks overlap, find_peak should estabish their boundary at the valley,
    and the reported peak bases should not overlap by N scans.
    But good to have reevaluation of ROI if future developement desires so.

    Return: merged peaks by extending bases, inheritting other attributes from largest peak.
    '''
    if len(cluster_peaks) == 1:
        return cluster_peaks[0]
    else:
        peak_sizes = [peak['right_base'] - peak['left_base'] for peak in cluster_peaks]
        largest = np.argmax(peak_sizes)
        largest = cluster_peaks[largest]
        largest['left_base'] = min([peak['left_base'] for peak in cluster_peaks])
        largest['right_base'] = max([peak['right_base'] for peak in cluster_peaks])
        return largest

def cleanup_peak_cluster(cluster_peaks):
    '''
    scipy find_peaks sometimes report two overlap peak: one small and the other joined with the small peak. 
    Separate them here. 
    Return list of peaks.
    '''
    if len(cluster_peaks) == 1:
        return cluster_peaks
    elif len(cluster_peaks) == 2:
        [peak1, peak2] = cluster_peaks
        bases = list(set([peak1['left_base'], peak1['right_base'], peak2['left_base'], peak2['right_base']]))
        bases.sort()
        peak1['left_base'], peak1['right_base'] = bases[:2]
        peak2['left_base'], peak2['right_base'] = bases[-2:]
        return [peak1, peak2]
    else:
        return [_merge_peak_cluster(cluster_peaks)]

def extend_ROI(ROI, number_of_scans):
    '''
    Add 3 datapoints to each end if ROI is too short (< 3*min_fwhm)
    '''
    left = [x for x in [ROI[0]-3, ROI[0]-2, ROI[0]-1] if x >=0]
    right = [x for x in [ROI[-1]+1, ROI[-1]+2, ROI[-1]+3] if x < number_of_scans]
    return left + ROI + right
