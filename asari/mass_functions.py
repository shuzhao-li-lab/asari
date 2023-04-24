'''
Functions related to mass operations, inlcuding mapping and clustering.

Functions here could potentially be sped up by JIT, but the code is currently weakly typed. 
In order to use Numba for JIT, they need to be rewritten with clear typing and likely compartmentalized.
Alternatively, some of the mass functions can be implemented in C and compiled to interface Python.
'''

import numpy as np
from scipy.signal import find_peaks 
from scipy.ndimage import uniform_filter1d

from mass2chem.search import *


def flatten_tuplelist(L):
    '''
    Given a list of tuples e.g., [(a,b), ...] flatten list to [a, b, ...] keeping unique entries only

    Parameters
    ----------
    L : list[tuple]
        list of tuples [(a,b), ...]

    Return
    ------
    list of the unique elements in L, flattened to one dimension: [a, b, ...]
    
    '''
    return list(set([x[0] for x in L] + [x[1] for x in L]))

def check_close_mzs(mzlist, ppm_tol=5):
    '''
    Given a list of mz_values and a mz tolerance in ppm, determine if any pair of values differ by less than the mz 
    tolerance. mzlist must be sorted in ascending order!

    Parameters
    ----------
    mzlist : list[float]
        list of floating point values represented m/z's in ascending order

    Return
    ------
    a list of the index pairs representing m/z values within the m/z tolerance
    
    '''
    warning = []
    for ii in range(1, len(mzlist)):
        _tolerance = mzlist[ii] * ppm_tol * 0.000001
        _d = mzlist[ii] - mzlist[ii-1]
        if _d < _tolerance:
            warning.append( (ii, ii-1) ) # (mzlist[ii], mzlist[ii-1]) )

    return warning


# -----------------------------------------------------------------------------
# Mass selectivity
# -----------------------------------------------------------------------------

def calculate_selectivity(sorted_mz_list, std_ppm=5):
    '''
    To calculate m or d-selectivity for a list of m/z values, 
    which can be all features in an experiment or a database.
    The mass selectivity between two m/z values is defined as: 
    (1 - Probability(confusing two peaks)), further formalized as an exponential model:
    P = exp( -x/std_ppm ),
    whereas x is ppm distance between two peaks, 
    std_ppm standard deviation of ppm between true peaks and theoretical values, default at 5 pmm.

    The selectivity value is between (0, 1), close to 1 meaning high selectivity.
    If multiple adjacent peaks are present, we multiply the selectivity scores.
    It is good approximation by considering 2 lower and 2 higher neighbors here.

    Parameters
    ----------
    sorted_mz_list: list 
        a list of m/z values, sorted from low to high, length > 3.
    std_ppm: float, optional, default: 5
        mass resolution in ppm (part per million).

    Returns
    -------
    A list of selectivity values, in matched order as the input m/z list.

    
    Note
    ----
    ppm is actually dependent on m/z, not an ideal method. 
    But it's in common practice and good enough approximation.
    '''
    def __sel__(x, std_ppm=std_ppm): 
        if x > 100:         # too high, not bother
            return 1
        elif x < 0.1:
            return 0
        else:
            return 1 - np.exp(-x/std_ppm)

    mz_list = np.array(sorted_mz_list)
    ppm_distances = 1000000 * (mz_list[1:] - mz_list[:-1])/mz_list[:-1]
    # first two MassTraces
    selectivities = [
        __sel__(ppm_distances[0]) * __sel__(ppm_distances[0]+ppm_distances[1]),
        __sel__(ppm_distances[0]) * __sel__(ppm_distances[1]
                                            )* __sel__(ppm_distances[1]+ppm_distances[2]),
        ]
    for ii in range(2, mz_list.size-2):
        selectivities.append(
            __sel__(ppm_distances[ii-2]+ppm_distances[ii-1]) * __sel__(
            ppm_distances[ii-1]) * __sel__(ppm_distances[ii]) * __sel__(
            ppm_distances[ii]+ppm_distances[ii+1])
        )
    # last two MassTraces
    selectivities += [
        __sel__(ppm_distances[-3]+ppm_distances[-2]) * __sel__(
        ppm_distances[-2]) * __sel__(ppm_distances[-1]),
        __sel__(ppm_distances[-2]+ppm_distances[-1]) * __sel__(ppm_distances[-1]),
        ]

    return selectivities


# -----------------------------------------------------------------------------
# Mapping and Alignment methods
# -----------------------------------------------------------------------------

def mass_paired_mapping(list1, list2, std_ppm=5):
    '''
    To find unambiguous matches of m/z values between two lists. 
    Nonunique matches are left out.

    This sorts all m/z values first, then compare their differences in sequential neighbors.
    To be considered as an unambiguous match, the m/z values from two lists 
    should have no overlap neighbors in either direction in either list other than their own pair.
    Thus, oow-selectiviy values are not considered in matching.
    This and related functions are for m/z alignment only, not used for general search. 
    See asari.tools.match_features for general search.

    This shares some similarity to the RANSAC algorithm but prioritizes selectivity.
    For illustration, one can use one-step Gaussian model for mass shift.
    Since only mean shift is used here, and stdev is implicitly enforced in matching, 
    no need to do model fitting.

    Parameters
    ----------
    list1: list[float]
        list of m/z values, not necessarially the same length as list2
    list2: list[float]
        list of m/z values, not necessarially the same length as list1
    std_ppm: float, optional, default: 5 
        limit of instrument accuracy in matching m/z values. 

    Returns
    -------
    mapped : 
        mapping list [(index from list1, index from list2), ...]
    ratio_deltas : 
        mean m/z ratio shift between two lists. This is ppm*10^-6. 
        No need to convert btw ppm here.

    Examples
    --------
        >>> list1 = [101.0596, 101.061, 101.0708, 101.0708, 101.1072, 101.1072, 101.1072, 102.0337, 
            102.0337, 102.0548, 102.0661, 102.0912, 102.0912, 102.1276, 102.1276, 103.0501, 
            103.0501, 103.0541, 103.0865, 103.0865, 103.9554, 104.0368, 104.0705, 104.0705, 
            104.1069, 104.1069, 104.9922, 105.0422, 105.0698, 105.0698, 105.0738, 105.1039, 
            105.1102, 105.9955, 106.0497, 106.065, 106.065, 106.0683, 106.0683, 106.0861, 106.0861, 
            106.0861, 106.1111, 106.9964, 107.0475, 107.0602, 107.0653, 107.0895, 107.9667, 108.0443, 
            108.0555, 108.0807, 109.0632, 109.0759]
        >>> list2 = [101.0087, 101.035, 101.0601, 101.0601, 101.0601, 101.0601, 101.0713, 101.0714, 
            101.1077, 101.1077, 101.1077, 101.1077, 101.1077, 101.1158, 101.1158, 102.0286, 102.0376, 
            102.0468, 102.0539, 102.0554, 102.0554, 102.0554, 102.0554, 102.0666, 102.0917, 102.0917, 
            102.0917, 102.0918, 102.1281, 102.1281, 102.1282, 103.0394, 103.0505, 103.0507, 103.0547, 
            103.1233, 103.8162, 103.956, 103.956, 103.956, 104.0532, 104.0533, 104.0641, 104.0709, 
            104.071, 104.0831, 104.0878, 104.0895, 104.0953, 104.1073, 104.1073, 104.1074, 104.1074, 
            104.1182, 104.1199, 104.1265, 104.1318, 104.1354, 104.1725, 104.3998, 104.9927, 104.9927, 
            104.9927, 104.9927, 105.0654, 105.0703, 105.1043, 105.1133, 106.049, 106.0503, 106.0655, 
            106.0688, 106.0866, 106.0867, 106.0867, 106.0867, 106.114, 107.048, 107.0481, 107.0496, 
            107.0608, 107.0658, 108.0109, 108.0482, 108.0604, 108.0812, 108.0812, 108.9618, 109.0507, 
            109.0637, 109.0637, 109.0764, 109.1015]
        >>> mass_paired_mapping(list1, list2) >>>
            ([(10, 23), (29, 65), (31, 66), (36, 70), (38, 71), (46, 81), (53, 91)],
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
    # Add a mock entry to allow loop goes through NN.
    all.append((999999, 2, None))
    mapped, ratio_deltas = [], []
    for ii in range(1, NN):
        if all[ii][1] != all[ii-1][1]:          # from two diff list_origin
            _tolerance = all[ii][0] * std_ppm * 0.000001
            _d = all[ii][0]-all[ii-1][0]
            if _d < _tolerance and all[ii+1][0]-all[ii][0] > _tolerance:
                # not allowing ii to be matched to both ii-1 and ii+1
                if all[ii][1] > all[ii-1][1]:   # always ordered as list1, list2
                    mapped.append( (all[ii-1][2], all[ii][2]) )
                    ratio_deltas.append( _d/all[ii][0] )
                else:
                    mapped.append( (all[ii][2], all[ii-1][2]) )
                    ratio_deltas.append( -_d/all[ii][0] )

    return mapped, ratio_deltas

def complete_mass_paired_mapping(list1, list2, std_ppm=5):
    '''
    To find best matched pairs of m/z values between two lists. 
    When multiple matches occur within std_ppm, this keeps the pair of closest m/z.

    This is different from mass_paired_mapping, which only keeps unique matches.
    Singletons are included in returned result if no matches are found.
    This does not calculate ratio_deltas.
    
    Parameters
    ----------
    list1: list[float]
        list of m/z values, not necessarially the same length as list2
    list2: list[float]
        list of m/z values, not necessarially the same length as list1
    std_ppm: float, optional, default: 5 
        limit of instrument accuracy in matching m/z values. 

    Returns
    -------
    mapped: 
        list of mapped index pairs. E.g. [ (3, 6), (6, 8), (33, 151), ...] 
    list1_unmapped : 
        list of unmapped items in list1.
    list2_unmapped : 
        list of unmapped items in list2.

    Note
    ----
    This and related functions are for m/z alignment only, not used for general search. 
    See asari.tools.match_features for general search.
    '''
    all = [(list1[ii], 1, ii) for ii in range(len(list1))] + [(list2[jj], 2, jj) for jj in range(len(list2))]
    # [(mz, list_origin, index_origin), ...]
    all.sort()
    NN = len(all)
    mapped = []
    for ii in range(1, NN):
        if all[ii][1] != all[ii-1][1]:          # from two diff list_origin
            _tolerance = all[ii][0] * std_ppm * 0.000001
            _d = all[ii][0]-all[ii-1][0]
            if _d < _tolerance:
                if all[ii][1] > all[ii-1][1]:   # always ordered as list1, list2
                    mapped.append( (all[ii-1][2], _d, all[ii][2]) )
                else:
                    mapped.append( (all[ii][2], _d, all[ii-1][2]) )

    # Now deal with multiple matches in either List1 or List2, smallest _d wins
    mapped2 = []
    mapped.append( (-1, -1, -1) )           # mock entry to allow loop below completes
    staged = mapped[0]
    for ii in range(1, len(mapped)):
        if mapped[ii][0] == staged[0] or mapped[ii][2] == staged[2]:
            if mapped[ii][1] < staged[1]:   # smaller _d 
                staged = mapped[ii]
        else:
            mapped2.append(staged)
            staged =  mapped[ii]         

    mapped = [(x[0], x[2]) for x in mapped2]
    # Now deal with singletons
    list1_unmapped = [x for x in range(len(list1)) if x not in [y[0] for y in mapped]]
    list2_unmapped = [x for x in range(len(list2)) if x not in [y[1] for y in mapped]]
    return mapped, list1_unmapped, list2_unmapped

def all_mass_paired_mapping(list1, list2, std_ppm=5):
    '''
    To find all matched pairs of m/z values between two lists. 
    When multiple matches occur within std_ppm, this keeps all pairs. 

    Singletons are included in returned result if no matches are found.
    This does not calculate ratio_deltas.

    Parameters
    ----------
    list1: list[float]
        list of m/z values, not necessarially the same length as list2
    list2: list[float]
        list of m/z values, not necessarially the same length as list1
    std_ppm: float, optional, default: 5 
        limit of instrument accuracy in matching m/z values. 
    Returns
    -------
    mapped: list of mapped index pairs. E.g. [ (3, 6), (6, 8), (33, 151), ...] 
    list1_unmapped : list of unmapped items in list1.
    list2_unmapped : list of unmapped items in list2.

    Note
    ----
    This function returns all matched pairs within std_ppm using 
    a tree search approach (mass2chem.search.build_centurion_tree_mzlist).
    The other functions, mass_paired_mapping and complete_mass_paired_mapping, 
    return only one pair per match and use a sorting based algorithm.
    '''
    mz_centurion_tree = build_centurion_tree_mzlist(list1)
    mapped = []
    for ii in range(len(list2)):
        result = _find_all_mzmatches_centurion_indexed_list(list2[ii], mz_centurion_tree, limit_ppm=std_ppm)
        mapped += [(x[1], ii) for x in result]

    # Now deal with singletons
    list1_unmapped = [x for x in range(len(list1)) if x not in [y[0] for y in mapped]]
    list2_unmapped = [x for x in range(len(list2)) if x not in [y[1] for y in mapped]]
    return mapped, list1_unmapped, list2_unmapped

def _find_all_mzmatches_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=5):
    '''
    Return matched peaks in mz_centurion_tree by m/z diff within limit_ppm.
    Internal use only.

    Parameters
    ----------
    query_mz: float
        one query m/z value
    mz_centurion_tree: dict
        indexed tree of m/z values, 
        mz_centurion_tree[cent] = [(mzList[ii], ii), ...]. 
        See mass2chem.search.build_centurion_tree_mzlist.
    std_ppm: float, optional, default: 5 
        limit of ppm accuracy in matching m/z values.

    Returns
    -------
    A list of matched m/z values and indices, [(mzList[ii], ii), ...].
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    results = []
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for x in L:         # (mzList[ii], ii)
            if abs(x[0] - query_mz) < mz_tol:
                results.append(x)
                
    return results

def mass_paired_mapping_with_correction(list1, list2, std_ppm=5, correction_tolerance_ppm=1):
    '''
    To find unambiguous matches of m/z values between two lists, 
    with correciton on list2 if m/z shift exceeds correction_tolerance_ppm.
    See `mass_paired_mapping` for details.

    Parameters
    ----------
    list1: list[float]
        list of m/z values, not necessarially the same length as list2
    list2: list[float]
        list of m/z values, not necessarially the same length as list1
    std_ppm: float, optional, default: 5 
        limit of instrument accuracy in matching m/z values. 
    correction_tolerance_ppm: float, optional, default: 1 
        threshold to trigger m/z recalibration of list2.

    Returns
    -------
    mapped: mapping list [(index from list1, index from list2), ...]
    _r: correction ratios on list2
    '''
    corrected_list2 = []
    mapped, ratio_deltas = mass_paired_mapping(list1, list2, std_ppm)
    _r = np.mean(ratio_deltas)
    if _r > correction_tolerance_ppm*0.000001:
        corrected_list2 = [x/(1+_r) for x in list2]
        mapped, ratio_deltas = mass_paired_mapping(list1, corrected_list2, std_ppm)

    return mapped, _r

def landmark_guided_mapping(REF_reference_mzlist, REF_landmarks, 
                            SM_mzlist, SM_landmarks, 
                            std_ppm=5, 
                            correction_tolerance_ppm=1
                            ):
    '''
    This is the main m/z alignment function in asari.
    A new sample is compared to the reference in a MassGrid, 
    then added to the MassGrid based on matched m/z values. 
    The MassGrid records the alignment information.

    This is a two-step process: aligning the anchors (landmarks) first, 
    mz correction if needed, then completing the remaining m/z values.
    Since the landmarks are of high confidence, this improves the quality 
    of m/z alignment.

    Parameters
    ----------
    REF_reference_mzlist: list
        the list of mz values from the REF sample
    REF_landmarks: list
        the list of landmarks in the REF sample
    SM_mzlist: list
        the list of mz values in the sample
    SM_landmarks: list
        the list of landmarks in the sample
    std_ppm: float, optional, default: 5
        the assumed mass resolution in ppm
    correction_tolerance_ppm: float, optional, default: 1
        a mass correction is applied if the mass shift is above this value

    Returns
    -------
    new_reference_mzlist : 
        combined list of all unique m/z values, 
        maintaining original order of REF_reference_mzlist but updating the values 
        as mean of the two lists.
    new_reference_map2 : 
        mapping index numbers from SM_malist, 
        to be used to update MassGrid[Sample.input_file]
    REF_landmarks : 
        updated landmark m/z values using the new index numbers 
        as part of new_reference_mzlist
    _r : 
        correction ratios on SM_mzlist, to be attached to Sample class instance.

    Note
    ----
    The m/z values are updated here because this is the best place to do it: 
    SM_mzlist is already corrected if needed; no need to look up irregular values in MassGrid.
    This mixes features from samples and they need to be consistent on how they are calibrated.

    The mzlists are already in ascending order when a Sample is processed,
    but the order of REF_reference_mzlist will be disrupted during building MassGrid.
    Do correciton on list2 if m/z shift exceeds correction_tolerance_ppm.
    See `MassGrid.add_sample`.

    '''
    _N1 = len(REF_reference_mzlist)
    # tracking how SM_mzlist is mapped to ref
    _d2, _r = {}, None   
    for ii in range(_N1): 
        _d2[ii] = None
    # first align to landmark mz values
    anchors_1 = [REF_reference_mzlist[x] for x in REF_landmarks]
    anchors_2 = [SM_mzlist[x] for x in SM_landmarks]

    mapped, ratio_deltas = mass_paired_mapping(anchors_1, anchors_2, std_ppm)
    # check number of mapped, to avoid forcing outlier sample into mass correction
    if len(mapped) > 0.2*_N1:
        _r = np.mean(ratio_deltas)
        if abs(_r) > correction_tolerance_ppm*0.000001:          # do m/z correction
            SM_mzlist = [x/(1+_r) for x in SM_mzlist]
            # rerun after mz correction
            anchors_2 = [SM_mzlist[x] for x in SM_landmarks]
            mapped, ratio_deltas = mass_paired_mapping(anchors_1, anchors_2, std_ppm)

    # convert back to index numbers in mzlists
    mapped = [( REF_landmarks[x[0]], SM_landmarks[x[1]] ) for x in mapped]
    # move onto remaining ions
    indices_remaining1 = [ii for ii in range(len(REF_reference_mzlist)) 
                          if ii not in [x[0] for x in mapped]]
    indices_remaining2 = [ii for ii in range(len(SM_mzlist)) 
                          if ii not in [x[1] for x in mapped]]
    mapped2, list1_unmapped, list2_unmapped = complete_mass_paired_mapping(
            [REF_reference_mzlist[ii] for ii in indices_remaining1], 
            [SM_mzlist[ii] for ii in indices_remaining2], 
            std_ppm)

    list2_unmapped = [indices_remaining2[ii] for ii in list2_unmapped]
    # Here we can have a few peaks that should have been merged during extraction of mass tracks
    mapped_pairs = mapped + [ ( indices_remaining1[x[0]], indices_remaining2[x[1]] ) 
                             for x in mapped2 ]
    print("    mapped pairs = %d / %d " %(len(mapped_pairs), len(SM_mzlist)))
    for p in mapped_pairs: 
        _d2[p[0]] = p[1]
        # updating ref m/z here
        REF_reference_mzlist[p[0]] = 0.5*( REF_reference_mzlist[p[0]] + SM_mzlist[p[1]] )
    for ii in range(len(list2_unmapped)): 
        _d2[_N1 + ii] = list2_unmapped[ii]
        # update landmark m/z values using the new index numbers as part of new_reference_mzlist
        if list2_unmapped[ii] in SM_landmarks:
            REF_landmarks.append(_N1 + ii)
    new_reference_mzlist = REF_reference_mzlist + [SM_mzlist[ii] for ii in list2_unmapped]
    new_reference_map2 = [_d2[x] for x in range(len(new_reference_mzlist))]

    return new_reference_mzlist, new_reference_map2, REF_landmarks, _r


# -----------------------------------------------------------------------------
# Binning and Clustering methods
# -----------------------------------------------------------------------------

def bin_by_median(List_of_tuples, func_tolerance):
    '''
    Separate List_of_tuples into bins within tolearance, by the distance of edge to median.
    
    Parameters
    ----------
    List_of_tuples : list[tuple]
        [(value, object), (value, object), ...], 
        to be separated into bins by values (either rt or mz).
        objects have attribute of sample_name if to align elsewhere.
    func_tolerance : function
        tolearance function to define bounary to separate bins.
        
    Returns
    -------
    A list of seprated bins, 
    each as a list of objects as [X[1] for X in L]. Possible all falls in same bin.

    Note
    ----
    Not perfect because left side may deviate out of tolerance, 
    even though LC-MS data always have enough gaps for separation.
    Use with caution.
    '''
    new = [[List_of_tuples[0], ], ]
    for X in List_of_tuples[1:]:
        if X[0]-np.median([ii[0] for ii in new[-1]]) < func_tolerance(X[0]):       
            # median moving with list change
            new[-1].append(X)
        else:
            new.append([X])
    PL = []
    for L in new:
        PL.append([X[1] for X in L])
    return PL

def gap_divide_mz_cluster(bin_data_tuples, mz_tolerance):
    '''
    Divides bin_data_tuples by the largest gap in  m/z values.
    mz_tolerance is not used now, assuming rarely needed after prior steps.
    This is a fallback method when `identify_mass_peaks` fails.
    See `nn_cluster_by_mz_seeds`.

    Parameters
    ----------
    bin_data_tuples: list[tuple]
        a flexible bin in format of [(mz, scan_num, intensity_int), ...], 
        or [(m/z, track_id, sample_id), ...].
    mz_tolerance: 
        the allowed tolerance in m/z values, NOT USED.

    Returns
    -------
    Two lists after dividing bin_data_tuples by the largest gap.
    '''
    def __divide_by_largest_gap__(L):
        gaps = [L[ii][0]-L[ii-1][0] for ii in range(1, len(L))]
        largest = np.argmax(gaps) + 1
        return L[:largest], L[largest:]
    return __divide_by_largest_gap__(bin_data_tuples)


def identify_mass_peaks(bin_data_tuples, mz_tolerance, presorted=True):
    '''
    Get the centroid m/z values as peaks in m/z values distribution, 
    at least mz_tolerance apart.
    This can be used to choose concensus m/z values in constructing 
    or after aligning mass tracks.
    See `nn_cluster_by_mz_seeds`.

    Parameters
    ----------
    bin_data_tuples : list[tuple]
        a flexible bin in format of [(mz, scan_num, intensity_int), ...], 
        or [(m/z, track_id, sample_id), ...].
    mz_tolerance : float
        precomputed based on m/z and ppm, 
        e.g. 5 ppm of 80 = 0.0004;  5 ppm of 800 = 0.0040.
    presorted : boolean, optional, default: True
        flag to determine if sorting is needed on bin_data_tuples.
    
    Returns
    -------
    A list of m/z values, peaks in m/z values distribution, at least mz_tolerance apart.
    '''

    #todo - should identify mass peaks have a default mz_tolerance of 5 ppm? like the other asari functions?
    tol4 = int(mz_tolerance * 10000)
    size = max(2, int(0.5 * tol4))
    mz4 = [int(x[0]*10000) for x in bin_data_tuples]
    if not presorted:
        mz4.sort()
    dict_count = {}
    positioned = range(mz4[0], mz4[-1]+1)
    for ii in positioned:
        dict_count[ii] = 0
    for x in mz4:
        dict_count[x] += 1
    values = uniform_filter1d(
        [dict_count[ii] for ii in positioned], size=size, mode='nearest'
        )
    peaks, _ = find_peaks( values, distance = tol4 )

    return [0.0001*positioned[ii] for ii in peaks]

def nn_cluster_by_mz_seeds(bin_data_tuples, mz_tolerance, presorted=True):
    '''
    Complete NN clustering, by assigning each data tuple to its closest m/z seed.
    This is used for both m/z alignment and mass track construction.
    See `chromatograms.build_chromatogram_by_mz_clustering` 
    and `constructors.MassGrid.bin_track_mzs`.

    Parameters
    ----------
    bin_data_tuples : list[tuple]
        a flexible bin in format of [(mz, scan_num, intensity_int), ...], 
        or [(m/z, track_id, sample_id), ...].
    mz_tolerance : float
        precomputed based on m/z and ppm, 
        e.g. 5 ppm of 80 = 0.0004;  5 ppm of 800 = 0.0040.
    presorted : boolean, optional, default: True
        flag to determine if sorting is needed on bin_data_tuples.

    Returns
    -------
    A list of clusters as separated bins, each bin as [(mz, scan_num, intensity_int), ...]

    Note
    ----
    Bug in Python compiler: when clusters = [[]] * _NN is used, 
    it causes occassional duplication of entries. 2022-05-21.

    Future consideration: np.argmin will be faster than sorted here.
    '''
    mz_seeds = identify_mass_peaks(bin_data_tuples, mz_tolerance, presorted)
    if mz_seeds:
        _NN = len(mz_seeds)
        # clusters = [[]] * _NN # see problem in Note
        clusters = []
        for x in mz_seeds:
            clusters.append([])
        # assign cluster number by nearest distance to a seed
        for x in bin_data_tuples:
            deltas = sorted([(abs(x[0] - mz_seeds[ii]), ii) for ii in range(_NN)])
            clusters[deltas[0][1]].append(x)
    else:
        clusters = [C for C in gap_divide_mz_cluster(bin_data_tuples, mz_tolerance)]

    return clusters
