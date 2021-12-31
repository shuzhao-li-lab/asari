'''

A few common data structures here:

For two m/z list comparison:
[(mz, list_origin, index_origin), ...]



to try performance optimization w/
numba JITC


'''

import numpy as np


from scipy.interpolate import UnivariateSpline



from collections import namedtuple
# from scipy.stats import norm as normal_distribution

from .search import *

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

# not used now
def _get_downstair_neighbor_mixedlist_(all_sorted_list, ii, ii_limit=0):
    # return nearest m/z neighbor of same list_origin in sorted [(mz, list_origin, index_origin), ...]
    # -1 indicates none found
    list_origin = all_sorted_list[ii][1]        # from list 1 or 2
    ii -= 1
    while ii >= ii_limit and all_sorted_list[ii][1] != list_origin:
        ii -= 1
    if all_sorted_list[ii][1] == list_origin:
        return ii
    else:
        return -1

def _get_upstair_neighbor_mixedlist_(all_sorted_list, ii, ii_limit):
    # return nearest m/z neighbor of same list_origin in sorted [(mz, list_origin, index_origin), ...]
    # -1 indicates none found
    list_origin = all_sorted_list[ii][1]        # from list 1 or 2
    ii += 1
    while ii < ii_limit and all_sorted_list[ii][1] != list_origin:
        ii += 1
    if all_sorted_list[ii][1] == list_origin:
        return ii
    else:
        return -1


# can be numba optimized
def mass_paired_mapping(list1, list2, std_ppm=5):
    '''
    To find unambiguous matches of m/z values between two lists.
    This sorts all m/z values first, then compare their differences in sequential neighbors.
    To be considered as an unambiguous match, the m/z values from two lists 
    should have no overlap neighbors in either direction in either list other than their own pair.
    Not necessary to have full mapping. FeatureMap is done in multiple steps. This is step 1.

    This shares some similarity to the RANSAC algorithm but prioritizes selectivity.
    For illustration, one can use one-step Gaussian model for mass shift.
    Since only mean shift is used here, and stdev is implicitly enforced in matching, no need to do model fitting.

    Input
    =====
    Two lists of m/z values, not ncessarily same length.
    std_ppm: instrument accuracy to guide value matching. Low-selectiviy values are not considered in matching.

    Return
    ======
    mapped: mapping list [(index from list1, index from list2), ...]
    ratio_deltas: mean m/z ratio shift between two lists. This is ppm*10^6. No need to convert btw ppm here.

    Test
    ====
    list1 = [101.0596, 101.061, 101.0708, 101.0708, 101.1072, 101.1072, 101.1072, 102.0337, 102.0337, 102.0548, 102.0661, 102.0912, 102.0912, 102.1276, 102.1276, 103.0501, 103.0501, 103.0541, 103.0865, 103.0865, 103.9554, 104.0368, 104.0705, 104.0705, 104.1069, 104.1069, 104.9922, 105.0422, 105.0698, 105.0698, 105.0738, 105.1039, 105.1102, 105.9955, 106.0497, 106.065, 106.065, 106.0683, 106.0683, 106.0861, 106.0861, 106.0861, 106.1111, 106.9964, 107.0475, 107.0602, 107.0653, 107.0895, 107.9667, 108.0443, 108.0555, 108.0807, 109.0632, 109.0759]
    list2 = [101.0087, 101.035, 101.0601, 101.0601, 101.0601, 101.0601, 101.0713, 101.0714, 101.1077, 101.1077, 101.1077, 101.1077, 101.1077, 101.1158, 101.1158, 102.0286, 102.0376, 102.0468, 102.0539, 102.0554, 102.0554, 102.0554, 102.0554, 102.0666, 102.0917, 102.0917, 102.0917, 102.0918, 102.1281, 102.1281, 102.1282, 103.0394, 103.0505, 103.0507, 103.0547, 103.1233, 103.8162, 103.956, 103.956, 103.956, 104.0532, 104.0533, 104.0641, 104.0709, 104.071, 104.0831, 104.0878, 104.0895, 104.0953, 104.1073, 104.1073, 104.1074, 104.1074, 104.1182, 104.1199, 104.1265, 104.1318, 104.1354, 104.1725, 104.3998, 104.9927, 104.9927, 104.9927, 104.9927, 105.0654, 105.0703, 105.1043, 105.1133, 106.049, 106.0503, 106.0655, 106.0688, 106.0866, 106.0867, 106.0867, 106.0867, 106.114, 107.048, 107.0481, 107.0496, 107.0608, 107.0658, 108.0109, 108.0482, 108.0604, 108.0812, 108.0812, 108.9618, 109.0507, 109.0637, 109.0637, 109.0764, 109.1015]
    mass_paired_mapping(list1, list2) >>>
        ([(33, 34),
        (94, 95),
        (97, 98),
        (106, 107),
        (109, 110),
        (127, 128),
        (144, 145)],
        [4.898762180656323e-06,
        4.758718686464085e-06,
        3.805743437700149e-06,
        4.714068193732999e-06,
        4.713921530199148e-06,
        4.670025348919892e-06,
        4.583942997773922e-06])
    '''
    all = [(list1[ii], 1, ii) for ii in range(len(list1))] + [(list2[jj], 2, jj) for jj in range(len(list2))]
    # [(mz, list_origin, index_origin), ...]
    all.sort()
    NN = len(all)
    mapped, ratio_deltas = [], []
    for ii in range(1, NN-1):
        if all[ii][1] != all[ii-1][1]:          # from two diff list_origin
            _tolerance = all[ii][0] * std_ppm * 0.000001
            _d = all[ii][0]-all[ii-1][0]
            if _d < _tolerance and all[ii+1][0]-all[ii][0] > _tolerance:
                # not allowing ii to be matched to both ii-1 and ii+1
                if all[ii][1] > all[ii-1][1]:   # list1, list2
                    mapped.append( (ii-1, ii) )
                    ratio_deltas.append( _d/all[ii][0] )
                else:
                    mapped.append( (ii, ii-1) )
                    ratio_deltas.append( -_d/all[ii][0] )

    #check if (NN-2, NN-1) is valid
    if mapped[-1][1] < NN-2 and all[NN-1][1] != all[NN-2][1]:
        _tolerance = all[NN-2][0] * std_ppm * 0.000001
        _d = all[NN-1][0] - all[NN-2][0]
        if _d < _tolerance:
            mapped.append( (NN-2, NN-1) )
            ratio_deltas.append( _d/all[NN-2][0] )

    return mapped, ratio_deltas


def mass_paired_mapping_with_correction(list1, list2, std_ppm=5, correction_tolerance_ppm=1):
    '''
    To find unambiguous matches of m/z values between two lists, 
    with correciton on list2 if m/z shift exceeds correction_tolerance_ppm.
    See `mass_paired_mapping` for details.

    Input
    =====
    Two lists of m/z values, not ncessarily same length.
    std_ppm: instrument accuracy to guide value matching. Low-selectiviy values are not considered in matching.

    Return
    ======
    mapped: mapping list [(index from list1, index from list2), ...]
    _r: correction ration on list2
    # ratio_deltas, corrected_list2
    '''
    corrected_list2 = []
    mapped, ratio_deltas = mass_paired_mapping(list1, list2, std_ppm)
    _r = np.mean(ratio_deltas)
    if _r > correction_tolerance_ppm*0.000001:
        corrected_list2 = [x/(1+_r) for x in list2]
        mapped, ratio_deltas = mass_paired_mapping(list1, corrected_list2, std_ppm)

    # return mapped, ratio_deltas, corrected_list2
    return mapped, _r






def mass_epdtree_mapping(epdtrees1, epdtrees2):
    '''


    To-do


    '''
    mapped = []
    return mapped






def peaks_to_features(peak_dict, rtime_tolerance, ordered_sample_names):
    
    '''peak_dict: {formula_mass or _M_id: list of Peaks}
    return List of Features (namedTuples, 'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'), 
    following the input sample order.
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



