'''
While classes are provided in metDataModel, 
mass tracks (i.e. XIC, EIC or chromatogram) and peaks are not instanced as classes in asari
internal processing for efficiency.
XICs as [( mz, rtlist, intensities ), ...].
Peak format: 
{
    'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 'height', 'parent_masstrace_id', 
    'rtime', 'peak_area', 'goodness_fitting'
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


def deep_detect_elution_peaks( mass_track, 
                min_intensity_threshold=10000, min_fwhm=3, min_prominence_threshold=5000, wlen=50, 
                snr=2, min_prominence_ratio=0.1, iterations=2):
    '''
    Peak detection on composite mass tracks; because this is at Experiment level, these are deemed as features.
    Mass tracks are expected to be continuous per m/z value, with 0s for gaps.

    Input
    =====
    mass_track: {'id_number': k, 'mz': mz, 'rt_scan_numbers': [..], 'intensity': [..]}
    iterations: default 2, must >=1. Number of iterations of peak detection; each done on remaining data points.
                The 2nd iteration may catch small peaks overshadowed by big ones. No clear need to go over 2.
    snr: signal to noise ratio. 
    min_prominence_ratio: require ratio of prominence relative to peak height.
    wlen: impacts the peak detection window. Peak boundaries should be re-examined in noisy tracks.

    Return list of peaks. 
    Reported left/right bases are not based on Gaussian shape or similar, just local maxima.
    Prominence has a key control of how peaks are considered, but peak shape evaluation later can filter out most bad peaks.

    peak area is integrated by summing up intensities of included scans.
    Peak area and height are cumulated from all samples. Not trying to average because some peaks are only few samples.

    Chromatographic peak selectivity (cSelectivity) is defined by 
    the ratio of the data points in all peaks above 1/2 this peak height and all data points above 1/2 this peak height.
    This value will correlate with SNR in most cases. 
    It's not a measure how good the chromatograph is, but how good the data are in defining clean peaks.
    E.g. chromatography may fail to separate two isomers, 
    but the one mixture peak can have perfect cSelectivity as far as computational processing is concerned.



    Peak format: {
        'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 'height', 'parent_masstrace_id', 
        'rtime', 'peak_area', 'goodness_fitting', 'snr',
    }

    
    To-do:
    step-wise peak detection; and cSelectivity calculation 

    change SNR back

    examine peak boundaries, maybe gaussian fit to guide?


    '''
    list_intensity = mass_track['intensity']
    __list_intensity = np.array(list_intensity)
    __list_intensity = __list_intensity[__list_intensity > 0]
    __max_intensity = __list_intensity.max()
    __noise_level_to_control__ = snr * np.percentile(__list_intensity, 10)
    # Noise level is defined as lower 10% of non-zero values - will change back to mean of non-peak data points

    prominence = max(min_prominence_threshold, 
                    min_prominence_ratio * __max_intensity)               # larger prominence gets cleaner data
    peaks, properties = find_peaks(list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                    prominence=prominence, wlen=wlen) 
    list_peaks = []
    # smooth noisy data, i.e. when >2 initial peaks are detected
    if peaks.size > 3:
        # noisy mass track
        new_list_intensity = smooth_moving_average(list_intensity)
        # raise prominence potentially according to __max_intensity
        prominence = max(min_prominence_threshold, 0.3 * __max_intensity) 
        peaks, properties = find_peaks(new_list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                        prominence=prominence, wlen=wlen)
        for ii in range(peaks.size):
            if properties['peak_heights'][ii] > __noise_level_to_control__:
                list_peaks.append(convert_peak_json__(ii, mass_track, peaks, properties))

    else:
        new_list_intensity = list_intensity.copy()                              # to substract peaks on the go
        for ii in range(peaks.size):
            if properties['peak_heights'][ii] > __noise_level_to_control__:
                list_peaks.append(convert_peak_json__(ii, mass_track, peaks, properties))

                for jj in range(properties['left_bases'][ii], properties['right_bases'][ii]+1):
                    new_list_intensity[jj] = 0

        # new iterations after removing the earlier peaks; adjusting parameters for min_int, prominence
        if peaks.size == 1:
            for _ in range(iterations-1): 
                #min_intensity_threshold = 0.1 * __max_intensity
                prominence = max(min_prominence_threshold, 0.3 * max(new_list_intensity))
                min_intensity_threshold = max(min_intensity_threshold, 2 * prominence)
                peaks, properties = find_peaks(new_list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                                prominence=prominence, wlen=wlen) 
                for ii in range(peaks.size):
                    if properties['peak_heights'][ii] > __noise_level_to_control__:
                        list_peaks.append(convert_peak_json__(ii, mass_track, peaks, properties))

    for peak in list_peaks:
        # evaluate peak quality by goodness_fitting of Gaussian model
        peak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, peak)

    return [peak for peak in list_peaks if peak['goodness_fitting'] > 0.01]


def convert_peak_json__( ii, mass_track, peaks, properties):
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
    }


def detect_elution_peaks( mass_track, 
            min_intensity_threshold=10000, min_fwhm=3, min_prominence_threshold=5000, wlen=50 ):
    '''Standard peak detection, to be used for exploration.
    Default in asari is the deep_detect_elution_peaks function.
    '''
    list_intensity = mass_track['intensity']
    peaks, properties = find_peaks(list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                    prominence=min_prominence_threshold, wlen=wlen) 
    list_peaks = []
    for ii in range(peaks.size):
        list_peaks.append(convert_peak_json__(ii, mass_track, peaks, properties))

    for peak in list_peaks:
        peak['goodness_fitting'] = evaluate_gaussian_peak(mass_track, peak)

    return list_peaks


def quick_detect_unique_elution_peak(rt_numbers, list_intensity, 
                            min_intensity_threshold=100000, min_fwhm=3, min_prominence_threshold_ratio=0.2):
    '''Quick peak detection, only looking for highest peak with high prominence.
    This is used for quick check on good peaks, or selecting landmarks for alignment purposes.
    rt_numbers, list_intensity are matched vectors from a mass trace/track.
    '''
    max_intensity = max(list_intensity)
    prominence = min_prominence_threshold_ratio * max_intensity
    unique_peak = None
    if max_intensity > min_intensity_threshold:
        peaks, properties = find_peaks(list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                        prominence=prominence) 
        if peaks.size == 1:
            unique_peak = {
                'apex': rt_numbers[peaks[0]], 
                'height': properties['peak_heights'][0], # not used now
            }
    return unique_peak

