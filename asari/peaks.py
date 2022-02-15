'''
While classes are provided in metDataModel, 
mass tracks (i.e. XIC, EIC or chromatogram) and peaks are not instanced as classes in asari
internal processing for efficiency.
XICs as [( mz, rtlist, intensities ), ...].
Peak format: 
{
    'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 'height', 'parent_masstrace_id', 
    'rtime', 'peak_area', 'goodness_fitting', 'snr', 'cSelectivity',
}
'''

import numpy as np
from scipy.signal import find_peaks 
from scipy.optimize import curve_fit 

from .chromatograms import *

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
# peak detection
# -----------------------------------------------------------------------------

def deep_detect_elution_peaks( mass_track, max_rt_number,
                min_peak_height=10000, min_fwhm=3, min_prominence_threshold=5000, wlen=50, 
                snr=2, min_prominence_ratio=0.1,
                iteration=True):
    '''
    Optimized peak detection on a mass track. 
    Local maxima is used in find_peak, but prominence and cSelectivity are used to control detection,
    in step-wise optimization. 
    Peak shape (Gaussian fitting goodness > 0.8) is required for noisy data or iterated detection.
    Noisy tracks are smoothed using moving average of 3 data points.
    peak area is integrated by summing up intensities of included scans.
    Reported left/right bases are not based on Gaussian shape or similar, just local maxima.
    Noise level is defined as mean of all non-peak, nonzero data points.

    Iterations are only done when good peaks are detected in the first round, 
    to ensure no small peaks are skipped due to prominence control in round 1,
    by removing round 1 peaks from the track in the next peak detection.
    Iterations are not performed on noisy data tracks.

    Input
    =====
    mass_track: {'id_number': k, 'mz': mz, 'rt_scan_numbers': [..], 'intensity': [..]}
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
    list_intensity = mass_track['intensity']
    __list_intensity = np.array(list_intensity)
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

        # do SNR, defined as mean of all non-peak data points
        _peak_datapoints = []
        for _jpeak in list_json_peaks:
            _peak_datapoints += list(range( _jpeak['left_index'], _jpeak['right_base'] ))
        _valid_datapoints = [x for x in range(len(list_intensity)) if x not in _peak_datapoints]
        __noise_level__ =  __list_intensity[_valid_datapoints]
        __noise_level__ = __noise_level__[__noise_level__>0]
        if __noise_level__.size > 0:
            __noise_level__ = __noise_level__.mean()
        else:
            __noise_level__ = 1
        for _jpeak in list_json_peaks:
            _jpeak['snr'] = int(_jpeak['height'] / __noise_level__)
            if _jpeak['snr'] > snr:
                list_peaks.append(_jpeak)

    return list_peaks


def convert_peak_json__( ii, mass_track, peaks, properties, cSelectivity=None):
    '''peaks, properties as from find_peaks; rt_numbers from mass_track
    '''
    rt_numbers  = mass_track['rt_scan_numbers']
    left_index, right_index = properties['left_bases'][ii], properties['right_bases'][ii]   # index positions on mass track
    left_base, right_base = rt_numbers[left_index], rt_numbers[right_index]
    peak_area = sum(mass_track['intensity'][left_index: right_index+1]) 
    return {
            'parent_masstrack_id': mass_track['id_number'],
            'mz': mass_track['mz'],
            'apex': rt_numbers[peaks[ii]], 
            # 'rtime': rt_numbers[peaks[ii]], 
            'peak_area': peak_area,
            'height': properties['peak_heights'][ii],
            'left_base': left_base,                                                         # rt_numbers
            'right_base': right_base,   
            'left_index': left_index,                                                       # specific to the referred mass track
            'right_index': right_index,
            'cSelectivity': cSelectivity,
    }

def __peaks_cSelectivity__(__list_intensity, peaks, properties):
    '''
    peaks, properties as from find_peaks; __list_intensity is np.array of the mass_track.

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
        list_cSelectivity.append(
            _peak_datapoints[_peak_datapoints > _threshold].size / __list_intensity[__list_intensity > _threshold].size
        )
    return list_cSelectivity



def detect_elution_peaks( mass_track, 
            min_peak_height=10000, min_fwhm=3, min_prominence_threshold=5000, wlen=50 ):
    '''Standard peak detection, to be used for exploration.
    Default in asari is the deep_detect_elution_peaks function.
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


def quick_detect_unique_elution_peak(rt_numbers, list_intensity, 
                            min_peak_height=100000, min_fwhm=3, min_prominence_threshold_ratio=0.2):
    '''Quick peak detection, only looking for highest peak with high prominence.
    This is used for quick check on good peaks, or selecting landmarks for alignment purposes.
    rt_numbers, list_intensity are matched vectors from a mass trace/track.
    '''
    max_intensity = max(list_intensity)
    prominence = min_prominence_threshold_ratio * max_intensity
    unique_peak = None
    if max_intensity > min_peak_height:
        peaks, properties = find_peaks(list_intensity, height=min_peak_height, width=min_fwhm, 
                                                        prominence=prominence) 
        if peaks.size == 1:
            unique_peak = {
                'apex': rt_numbers[peaks[0]], 
                'height': properties['peak_heights'][0], # not used now
            }
    return unique_peak

