'''

to try performance optimization w/
numba JITC


'''

import numpy as np

#from .lcms_experiment import Feature

from collections import namedtuple

#
# -----------------------------------------------------------------------------
# feature id will be assigned at the end; intensities is a list; mass_id links to MassTrace
Feature = namedtuple('Feature', ['feature_id', 'mass_id', 'mz', 'rtime', 'rt_min', 'rt_max', 
                                'peak_quality_max', 'peak_quality_median', 'number_peaks', 'perc_peaks',
                                'selectivity_combined', 'selectivity_mz', 'intensity_mean',
                                'intensities'])




def gaussian_function__(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) 

def goodness_fitting__(y_orignal, y_fitted):                  # R^2 as goodness of fitting
    return 1 - (np.sum((y_fitted-y_orignal)**2) / np.sum((y_orignal-np.mean(y_orignal))**2))

def calculate_selectivity(sorted_mz_list, std_ppm=5):
    '''
    To calculate selectivity for all valid mass traces (thus specific m/z values).

    The mass selectivity between two m/z values, between (0, 1), is defined as: 
    (1 - Probability(confusing two peaks)), further formalized as an exponential model:

    P = exp( -x/std_ppm ),
    whereas x is ppm distance between two peaks, 
    std_ppm standard deviation of ppm between true peaks and theoretical values, default at 5 pmm.

    Close to 1 means high selectivity.
    If multiple adjacent peaks are present, we multiply the selectivity scores.
    Considering 2 lower and 2 higher neighbors approximately here.

    Future direction: std_ppm can be dependent on m/z. 
    This can be taken into account by a higher order model.
    '''
    def __sel__(x, std_ppm=std_ppm): 
        if x > 100:         # too high, not bother
            return 1
        else:
            return 1 - np.exp(-x/std_ppm)

    mz_list = np.array(sorted_mz_list)
    ppm_distances = 1000000 * (mz_list[1:] - mz_list[:-1])/mz_list[:-1]
    # first two MassTraces
    selectivities = [
        __sel__(ppm_distances[0]) * __sel__(ppm_distances[0]+ppm_distances[1]),
        __sel__(ppm_distances[0]) * __sel__(ppm_distances[1])* __sel__(ppm_distances[1]+ppm_distances[2]),
        ]
    for ii in range(2, mz_list.size-2):
        selectivities.append(
            __sel__(ppm_distances[ii-2]+ppm_distances[ii-1]) * __sel__(ppm_distances[ii-1]) * __sel__(ppm_distances[ii]) * __sel__(ppm_distances[ii]+ppm_distances[ii+1])
        )
    # last two MassTraces
    selectivities += [
        __sel__(ppm_distances[-3]+ppm_distances[-2]) * __sel__(ppm_distances[-2]) * __sel__(ppm_distances[-1]),
        __sel__(ppm_distances[-2]+ppm_distances[-1]) * __sel__(ppm_distances[-1]),
        ]

    return selectivities


def bin_by_median(List_of_tuples, func_tolerance):
    '''
    Not perfect because left side may deviate out of tolerance, but LC-MS data always have enough gaps for separation.
    Will add kernel density method for grouping m/z features.
    List_of_tuples: [(value, object), (value, object), ...], to be separated into bins by values (either rt or mz).
                    objects have attribute of sample_name if to align elsewhere.
    return: [seprated bins], each as a list of objects as [X[1] for X in L]. Possible all falls in same bin.
    '''
    new = [[List_of_tuples[0], ], ]
    for X in List_of_tuples[1:]:
        if X[0]-np.median([ii[0] for ii in new[-1]]) < func_tolerance(X[0]):       # median moving with list change
            new[-1].append(X)
        else:
            new.append([X])
    PL = []
    for L in new:
        PL.append([X[1] for X in L])
    return PL

def peaks_to_features(peak_dict, rtime_tolerance, ordered_sample_names):
    '''
    peak_dict: {formula_mass or _M_id: list of Peaks}
    return List of Features (namedTuples, 'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'), following the input sample order.
    '''
    def __get_peaks_intensities__(peaks, ordered_sample_names):
        dict_intensities = {}
        for P in peaks: 
            if P.sample_name in dict_intensities:
                dict_intensities[P.sample_name] += P.peak_area      # redundant peaks in same m/z & rt are merged/summed here
            else:
                dict_intensities[P.sample_name] = P.peak_area
        return [dict_intensities.get(name, 0) for name in ordered_sample_names]

    def __bin_by_median_rt__(List_of_peaks, tolerance):
        List_of_tuples = [(P.cal_rtime, P) for P in List_of_peaks]
        List_of_tuples.sort()
        return bin_by_median(List_of_tuples, lambda x: max(tolerance, 0.1*x))

    FeatureList = []            # still possibly not distinguishing close peaks well, which may be treated as combined peaks. 
    for k,v in peak_dict.items():
        for F in __bin_by_median_rt__(v, rtime_tolerance):
            median_mz, median_rt = np.median([P.mz for P in F]), np.median([P.cal_rtime for P in F])
            feature_id = str(round(median_mz,4)) + '@' + str(round(median_rt,2))
            rt_min = float(min([P.left_rtime for P in F]))          # checking if this is np.array ???
            rt_max = float(max([P.right_rtime for P in F]))
            peak_quality_max = max([P.goodness_fitting for P in F])
            peak_quality_median = np.median([P.goodness_fitting for P in F])
            number_peaks = len(set([P.sample_name for P in F]))
            perc_peaks = 100.0 * number_peaks/len(ordered_sample_names)
            selectivity_mz = np.mean([P.selectivity for P in F])
            selectivity_combined = 9
            intensities = __get_peaks_intensities__(F, ordered_sample_names)
            intensity_mean = int(sum(intensities)/number_peaks)
            FeatureList += [Feature(feature_id, k, median_mz, median_rt, rt_min, rt_max, 
                            peak_quality_max, peak_quality_median, number_peaks, perc_peaks,
                            selectivity_combined, selectivity_mz, intensity_mean, intensities)]

    return FeatureList
