'''
Functions for elution peak detection and evaluation.
stats_detect_elution_peaks and deep_detect_elution_peaks use different supporting functions,
as the former deals with conversion of RT coordinates.
'''

import multiprocessing as mp
import numpy as np
from scipy.signal import find_peaks 
from scipy.optimize import curve_fit 

from .chromatograms import *

# -----------------------------------------------------------------------------
# multicore processing for peak detection
# -----------------------------------------------------------------------------

def batch_deep_detect_elution_peaks(list_mass_tracks, max_rt_number, parameters):
    with mp.Manager() as manager:
        shared_list = manager.list()
        iters = iter_peak_detection_parameters(list_mass_tracks, max_rt_number, parameters, shared_list)
        with mp.Pool( parameters['multicores'] ) as pool:
            # switch peak detection algorithms
            # pool.starmap( deep_detect_elution_peaks, iters )
            pool.starmap( stats_detect_elution_peaks, iters )

        FeatureList = list(shared_list)
    return FeatureList

def iter_peak_detection_parameters(list_mass_tracks, max_rt_number, parameters, shared_list):
    '''
    Prominence requirement is critical to peak detection. 
    This is dynamically determined based on parameters['min_prominence_threshold'] 
    and baseline level in a region.
    '''
    iters = []
    min_peak_height = parameters['min_peak_height']
    min_prominence_threshold = parameters['min_prominence_threshold']
    min_fwhm = round( 0.5 * parameters['min_timepoints'] )
    snr = parameters['signal_noise_ratio']
    peakshape = parameters['gaussian_shape']
    wlen = 101                          # if used, this limits peak evaluation 50 scans to the left and 50 to right
    min_prominence_ratio = 0.1
    iteration = True                    # a 2nd round of peak detection if enough remaining datapoints
    for mass_track in list_mass_tracks:
        if mass_track['intensity'].max() > min_peak_height:
            iters.append(
                (mass_track, max_rt_number, min_peak_height, min_fwhm, min_prominence_threshold,
                wlen, snr, peakshape, min_prominence_ratio, iteration, shared_list)
            )
    return iters


# -----------------------------------------------------------------------------
# Statistics guided peak detection
# -----------------------------------------------------------------------------

def stats_detect_elution_peaks(mass_track, max_rt_number, 
                min_peak_height, min_fwhm, min_prominence_threshold,
                wlen, snr, peakshape, min_prominence_ratio, iteration,
                shared_list):
    '''
    Stats guided peak detection. Use statistical summary (stdev) to get ROI (region of interest), 
    apply smooth_moving_average and determine prominence per ROI.
    The advantages of separating ROI are 1) better performance of long LC runs, and
    2) big peaks are less likely to shallow over small peaks. 
    iteration: if True, a 2nd round of peak detection if enough remaining datapoints.

    Gap allowed for 1 scan. 
    Actually in composite mass track - gaps should not exist after combining many samples.
    Now position indices have to be converted internally in this function.

    Input
    =====
    mass_track: {'id_number': k, 'mz': mz, 'rt_scan_numbers': [..], 'intensity': [..]}
                Mass tracks are expected to be continuous per m/z value, with 0s for gaps.
    snr: minimal signal to noise ratio required for a peak. 
                Noise is defined by the mean of all non-peak data points in mass_track.
    min_peak_height, min_fwhm, min_prominence_threshold as in main parameters.
    min_prominence_ratio: require ratio of prominence relative to peak height.
    wlen, iteration, are not used in this function but kept for consistent format with others.

    Update
    ======
    shared_list: list of peaks in JSON format, to pool with batch_deep_detect_elution_peaks.
    '''
    list_json_peaks, list_peaks = [], []
    list_intensity = mass_track['intensity']
    list_scans = np.arange(max_rt_number)
    # This baseline method down weight the peak data points
    __baseline__ = 2**(np.log2(list_intensity + min_peak_height*0.01).mean())
    __selected_scans__ = list_scans[list_intensity > __baseline__]
    __noise_level__ = max(__baseline__,
                            np.quantile( list_intensity[__selected_scans__], 0.05 ))

    ROIs = []                           # get ROIs by separation/filtering with __baseline__, allowing 2 gap
    tmp = [__selected_scans__[0]]
    for ii in __selected_scans__[1:]:
        if ii - tmp[-1] < 3:
            tmp += range(tmp[-1]+1, ii+1)
        else:
            ROIs.append(tmp)
            tmp = [ii]

    ROIs.append(tmp)
    ROIs = [r for r in ROIs if len(r) >= min_fwhm * 2]
    _do_smoothing = True                # noisy mass track
    if list_intensity.max() > 100 * min_peak_height and len(__selected_scans__) < 0.5 * max_rt_number:
        _do_smoothing = False           # clean mass track 

    # min_prominence_ratio * list_intensity_roi.max() is not good because high noise will suppress peaks
    prominence = max(min_prominence_threshold, snr*__baseline__)
    for R in ROIs:
        if len(R) < 3 * min_fwhm:       # extend if too short - help narrow peaks
            R = extend_ROI(R, max_rt_number)
        list_json_peaks += _detect_regional_peaks( list_intensity[R], R, 
                    min_peak_height, min_fwhm, prominence, min_prominence_ratio, 
                    _do_smoothing
        )

    # evaluate and format peaks
    list_cSelectivity = __peaks_cSelectivity_stats_(list_intensity, list_json_peaks)
    for ii in range(len(list_json_peaks)):
        list_json_peaks[ii]['cSelectivity'] = list_cSelectivity[ii]
    #__noise_level__ = __baseline__
    for peak in list_json_peaks:
        peak['parent_masstrack_id'] = mass_track['id_number']
        peak['mz'] = mass_track['mz']
        peak['snr'] = int( min(peak['height'], 99999999) / __noise_level__+1)         # cap upper limit and avoid INF
        peak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, peak)
        if peak['snr'] >= snr and peak['goodness_fitting'] > peakshape:
            list_peaks.append(peak)

    shared_list += list_peaks


def _detect_regional_peaks(list_intensity_roi, rt_numbers_roi, 
                    min_peak_height, min_fwhm, min_prominence_threshold, min_prominence_ratio,
                    smoothing=True):
    '''
    Return list of peaks based on detection in ROI, defined by (list_intensity_roi, rt_numbers_roi).
    smooth_moving_average is applied before peak detection.
    Smoothing will reduce the height of narrow peaks in CMAP, but not on the reported values,  
    because peak area is extracted from each sample after. 
    The likely slight expansion of peak bases can add to robustness.
    '''
    list_peaks = []
    if smoothing:
        list_intensity_roi = smooth_moving_average(list_intensity_roi, size=min_fwhm)

    prominence = max(min_prominence_ratio * list_intensity_roi.max(), min_prominence_threshold)
    peaks, properties = find_peaks(list_intensity_roi, 
                                    height=min_peak_height, width=min_fwhm, prominence=prominence) 
    for ii in range(peaks.size):
        _jpeak = convert_roi_peak_json_(ii, list_intensity_roi, rt_numbers_roi, peaks, properties, )
        list_peaks.append(_jpeak)

    return check_overlap_peaks(list_peaks)
    

def convert_roi_peak_json_( ii, list_intensity_roi, rt_numbers_roi, peaks, properties):
    '''
    peaks, properties as from scipy find_peaks; list_intensity_roi, rt_numbers_roi from mass_track.
    This handles the conversion btw indices of ROI and indices of mass_track.
    '''
    left_index, right_index = properties['left_bases'][ii], properties['right_bases'][ii]   # index positions on ROI
    left_base, right_base = rt_numbers_roi[left_index], rt_numbers_roi[right_index]
    peak_area = int( list_intensity_roi[left_index: right_index+1].sum() ) 
    return {
            'apex': rt_numbers_roi[peaks[ii]], 
            'peak_area': peak_area,
            'height': int(properties['peak_heights'][ii]),
            'left_base': left_base,
            'right_base': right_base, 
    }

def __peaks_cSelectivity_stats_(__list_intensity, _jpeaks):
    '''
    peaks, properties as from find_peaks; 
    __list_intensity is np.array of the mass_track.
    See also __peaks_cSelectivity__.
    '''
    _peak_datapoints = []
    for peak in _jpeaks:
            _peak_datapoints += list(range(peak['left_base'], peak['right_base']))

    _peak_datapoints = __list_intensity[list(set(_peak_datapoints))]    # peaks may overlap
    list_cSelectivity = []
    for ii in range(len(_jpeaks)):
        _threshold = 0.5 * _jpeaks[ii]['height']
        _peak_datapoints_level = _peak_datapoints[_peak_datapoints > _threshold].size
        _background_level = __list_intensity[__list_intensity > _threshold].size
        if _background_level >= _peak_datapoints_level > 0:                  # to rule out artifact from smooth_lowess
            list_cSelectivity.append(_peak_datapoints_level / _background_level)
        else:
            list_cSelectivity.append( 0 )
            
    return list_cSelectivity


# -----------------------------------------------------------------------------
# Heuristics peak detection
# -----------------------------------------------------------------------------

def deep_detect_elution_peaks(mass_track, max_rt_number, 
                min_peak_height, min_fwhm, min_prominence_threshold,
                wlen, snr, peakshape, min_prominence_ratio, iteration,
                shared_list):
    '''
    Optimized peak detection on a mass track. 
    Local maxima is used in find_peak, but prominence and cSelectivity are used to control detection,
    in step-wise optimization. 
    Peak shape (Gaussian fitting goodness > 0.8) is required for noisy data or iterated detection.
    Noisy tracks are smoothed using moving average of N data points.
    peak area is integrated by summing up intensities of included scans.
    Reported left/right bases are not based on Gaussian shape or similar, just local maxima.
    Noise level is defined as mean of all non-peak, nonzero data points.

    Iterations should only be done when good peaks are detected in the first round, 
    to ensure no small peaks are skipped due to prominence control in round 1,
    by removing round 1 peaks from the track in the next peak detection.
    Iterations are not performed on noisy data tracks.

    Input
    =====
    mass_track: {'id_number': k, 'mz': mz, 'intensity': [..]}
                Mass tracks are expected to be continuous per m/z value, with 0s for gaps.
    iteration:  Second iteration of peak detection; done on remaining data points.
                The 2nd iteration may catch small peaks overshadowed by big ones. No clear need to go over 2.
    snr: signal to noise ratio. 
    min_prominence_ratio: require ratio of prominence relative to peak height.
    wlen: impacts the peak detection window. Peak boundaries should be re-examined in noisy tracks.

    Return list of peaks in JSON format, e.g. [{ 'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 
                                                'height', 'parent_masstrace_id', 'rtime', 'peak_area', 
                                                'goodness_fitting', 'snr', 'cSelectivity',}, ...]
    '''
    list_json_peaks, list_peaks = [], []
    list_intensity = __list_intensity = mass_track['intensity']   # np.array(list_intensity)
    
    __max_intensity = __list_intensity.max()
    __exam_window_size = 0.1*max_rt_number
    prominence = max(min_prominence_threshold, 
                    min_prominence_ratio * __max_intensity)               # larger prominence gets cleaner data
    peaks, properties = find_peaks(__list_intensity, height=min_peak_height, width=min_fwhm, 
                                                    prominence=prominence, wlen=wlen) 
    list_cSelectivity = __peaks_cSelectivity__(__list_intensity, peaks, properties)
    # cSelectivity indicates how much noise exists relative to a peak
    if peaks.size > 0:
        new_list_intensity = list_intensity.copy()
        if peaks.size <= 2:
            for ii in range(peaks.size):
                _jpeak = convert_peak_json__(ii, mass_track, peaks, properties, cSelectivity=list_cSelectivity[ii])
                # evaluate peak quality by goodness_fitting of Gaussian model
                _jpeak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, _jpeak)
                list_json_peaks.append(_jpeak)
                for jj in range(properties['left_bases'][ii], properties['right_bases'][ii]+1):
                    new_list_intensity[jj] = 0

            if iteration and __list_intensity.size > __exam_window_size and min(list_cSelectivity) > 0.9:  
                # very good peak(s) detected, do iterations
                # no need to iterate on short tracks, and requiring high quality peak on iteration 2
                # new iterations after removing the earlier peaks; adjusting parameters for min_int, prominence
                prominence = max(min_prominence_threshold, 0.3 * max(new_list_intensity))
                min_peak_height = max(min_peak_height, 2 * prominence)
                peaks, properties = find_peaks(new_list_intensity, height=min_peak_height, width=min_fwhm, 
                                                                prominence=prominence, wlen=wlen) 
                list_cSelectivity = __peaks_cSelectivity__(__list_intensity, peaks, properties)
                for ii in range(peaks.size):
                    _jpeak = convert_peak_json__(ii, mass_track, peaks, properties, cSelectivity=list_cSelectivity[ii])
                    # evaluate peak quality by goodness_fitting of Gaussian model
                    _jpeak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, _jpeak)
                    if _jpeak['goodness_fitting'] > 0.8:    # requiring good peaks after iteration or in noisy data
                        list_json_peaks.append(_jpeak)

        else:  
            # smoothing will suppress narrow peaks. So make sure good peaks are fetched before smoothihng
            for ii in range(peaks.size):
                # if list_cSelectivity[ii] > 0.9: # this is too dependent on baseline
                _jpeak = convert_peak_json__(ii, mass_track, peaks, properties, cSelectivity=list_cSelectivity[ii])
                _jpeak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, _jpeak)
                if _jpeak['goodness_fitting'] > 0.8:
                    list_json_peaks.append(_jpeak)
                    for jj in range(properties['left_bases'][ii], properties['right_bases'][ii]+1):
                        new_list_intensity[jj] = 0

            new_list_intensity = smooth_moving_average(new_list_intensity, size=9)
            # raise prominence potentially according to __max_intensity
            prominence = max(min_prominence_threshold, 0.3 * max(new_list_intensity)) 
            min_peak_height = max(min_peak_height, 2 * prominence)
            peaks, properties = find_peaks(new_list_intensity, height=min_peak_height, width=min_fwhm, 
                                                            prominence=prominence, wlen=wlen)
            list_cSelectivity = __peaks_cSelectivity__(__list_intensity, peaks, properties)
            for ii in range(peaks.size):
                _jpeak = convert_peak_json__(ii, mass_track, peaks, properties, cSelectivity=list_cSelectivity[ii])
                _jpeak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, _jpeak)
                if _jpeak['goodness_fitting'] > 0.6:    # requiring good peaks after iteration or in noisy data
                    list_json_peaks.append(_jpeak)

            # check too crowded peaks, usually from splitting; fix by more smoothing and open wlen
            __NN = len(list_json_peaks)
            if __NN > 4:
                list_apexes = sorted([_jpeak['apex'] for _jpeak in list_json_peaks])
                _crowded = False
                for ii in range(3, __NN):
                    if list_apexes[ii] - list_apexes[ii-3] > __exam_window_size:
                        _crowded = True
                if _crowded:
                    min_fwhm = 3 * min_fwhm
                    new_list_intensity = smooth_lowess(__list_intensity, frac=0.1)
                    # new_smooth_size = max(15, int(0.02*__list_intensity.size))
                    # new_list_intensity = smooth_moving_average(list_intensity, size=new_smooth_size)
                    list_json_peaks = []
                    peaks, properties = find_peaks(new_list_intensity, height=min_peak_height, width=min_fwhm, 
                                                            prominence=prominence) #, wlen=2*__exam_window_size)  # larger wlen
                    list_cSelectivity = __peaks_cSelectivity__(__list_intensity, peaks, properties)
                    for ii in range(peaks.size):
                        _jpeak = convert_peak_json__(ii, mass_track, peaks, properties, cSelectivity=list_cSelectivity[ii])
                        _jpeak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, _jpeak)
                        if _jpeak['goodness_fitting'] > 0.6:    # requiring good peaks after iteration or in noisy data
                            list_json_peaks.append(_jpeak)

        # do SNR, noise defined as mean of all non-peak data points
        _peak_datapoints = []
        for _jpeak in list_json_peaks:
            _peak_datapoints += list(range( _jpeak['left_base'], _jpeak['right_base'] ))
        _valid_datapoints = [x for x in range(len(list_intensity)) if x not in _peak_datapoints]
        __noise_level__ =  __list_intensity[_valid_datapoints]
        __noise_level__ = __noise_level__[__noise_level__>0]
        if __noise_level__.size > 0:
            __noise_level__ = __noise_level__.mean()
        else:
            __noise_level__ = max(__list_intensity[_peak_datapoints].min(), 1)
        for _jpeak in list_json_peaks:
            _jpeak['snr'] = int( min(_jpeak['height'], 99999999) / __noise_level__)         # cap upper limit and avoid INF
            if _jpeak['snr'] > snr:
                list_peaks.append(_jpeak)

    shared_list += check_overlap_peaks(list_peaks)

def convert_peak_json__( ii, mass_track, peaks, properties, cSelectivity=None):
    '''
    peaks, properties as from scipy find_peaks; rt_numbers from mass_track
    When asari switched from mass traces to mass tracks, 
    the left and right indices no longer needed converstion btw local indices and full RT indices.
    '''
    #  rt_numbers  = mass_track['rt_scan_numbers']
    left_index, right_index = properties['left_bases'][ii], properties['right_bases'][ii]   # index positions on mass track
    #  left_base, right_base = rt_numbers[left_index], rt_numbers[right_index]
    peak_area = int( mass_track['intensity'][left_index: right_index+1].sum() ) 
    return {
            'parent_masstrack_id': mass_track['id_number'],
            'mz': mass_track['mz'],
            'apex': peaks[ii],   #  rt_numbers[peaks[ii]], 
            # 'rtime': rt_numbers[peaks[ii]], 
            'peak_area': peak_area,
            'height': int(properties['peak_heights'][ii]),
            'left_base': left_index, # left_base,                                            # rt_numbers
            'right_base': right_index, # right_base,   
            #'left_index': left_index,                                                       # specific to the referred mass track
            #'right_index': right_index,
            'cSelectivity': cSelectivity,
    }

def __peaks_cSelectivity__(__list_intensity, peaks, properties):
    '''
    peaks, properties as from find_peaks; 
    __list_intensity is np.array of the mass_track.

    Chromatographic peak selectivity (cSelectivity) is defined by 
    the ratio of the data points in all peaks above 1/2 this peak height and all data points above 1/2 this peak height.
    This value will correlate with SNR in most cases. 
    It's not a measure how good the chromatograph is, but how good the data are in defining clean peaks.
    E.g. chromatography may fail to separate two isomers, 
    but the one mixture peak can have perfect cSelectivity as far as computational processing is concerned.
    '''
    _peak_datapoints = []
    for ii in range(peaks.size):
        _peak_datapoints += list(range(properties['left_bases'][ii], properties['right_bases'][ii]))
    _peak_datapoints = __list_intensity[list(set(_peak_datapoints))]    # peaks may overlap
    list_cSelectivity = []
    for ii in range(peaks.size):
        _threshold = 0.5 * properties['peak_heights'][ii]
        _peak_datapoints_level = _peak_datapoints[_peak_datapoints > _threshold].size
        _background_level = __list_intensity[__list_intensity > _threshold].size
        if _background_level >= _peak_datapoints_level > 0:                  # to rule out artifact from smooth_lowess
            list_cSelectivity.append(_peak_datapoints_level / _background_level)
        else:
            list_cSelectivity.append( 0 )
            
    return list_cSelectivity


# -----------------------------------------------------------------------------
# Generic functions for detection and evaluation
# -----------------------------------------------------------------------------

def gaussian_function__(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) 

def goodness_fitting__(y_orignal, y_fitted):                  # R^2 as goodness of fitting
    return 1 - (np.sum((y_fitted-y_orignal)**2) / np.sum((y_orignal-np.mean(y_orignal))**2))

def evaluate_gaussian_peak(mass_track, peak):
    '''
    Use Gaussian models to fit peaks, R^2 as goodness of fitting.
    Peak shapes may be more Voigt or bigaussian, but the impact on fitting result is negligible.
    mass_track: {'id_number': k, 'mz': mz, 'intensity': [..]}
    Peak: {'parent_masstrace_id': 2812, 'mz': 359.9761889867791, 'apex': 91.0, 'height': 491241.0, 'left_base': 49.0, 'right_base': 267.0}
    return: goodness_fitting
    '''
    goodness_fitting = 0
    xx = range(peak['left_base'], peak['right_base']+1, 1)
    yy = mass_track['intensity'][xx]
    # yy = mass_track['intensity'][peak['left_index']: peak['right_index']+1]
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

def detect_elution_peaks( mass_track, 
            min_peak_height=10000, min_fwhm=3, min_prominence_threshold=5000, wlen=50 ):
    '''
    Vanila peak detection, to be used for custom processes or exploration.
    Refer to stats_detect_elution_peaks or deep_detect_elution_peaks for regular functions.
    '''
    list_intensity = mass_track['intensity']
    peaks, properties = find_peaks(list_intensity, height=min_peak_height, width=min_fwhm, 
                                                    prominence=min_prominence_threshold, wlen=wlen) 
    list_peaks = []
    for ii in range(peaks.size):
        list_peaks.append(convert_peak_json__(ii, mass_track, peaks, properties))

    for peak in list_peaks:
        peak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, peak)

    return list_peaks

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

def extend_ROI(ROI, max_rt_number):
    '''
    Add 3 datapoints to each end if ROI is too short (< 3*min_fwhm)
    '''
    left = [x for x in [ROI[0]-3, ROI[0]-2, ROI[0]-1] if x >=0]
    right = [x for x in [ROI[-1]+1, ROI[-1]+2, ROI[-1]+3] if x < max_rt_number]
    return left + ROI + right
