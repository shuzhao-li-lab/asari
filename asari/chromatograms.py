'''
Use integers for RT scan numbers and intensities.

Flexible binning based on ppm accuracy.


Changing XIC procedure.
Use mz tol (default 5 pmm) in XIC construction. 
XICs without neighbors within x ppm are considered specific (i.e. high selectivity). 
Low selectivity regions will be still inspected to determine the true number of XICs.

Leave out calibration for now.


import os

from itertools import combinations
# import matplotlib.pyplot as plt
# from sklearn.neighbors import KernelDensity
# from pyopenms import *
'''

from operator import itemgetter
import numpy as np
from scipy.signal import find_peaks 


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
        for (mz, intensity) in zip(*ms_expt[ii].get_peaks()):
            if intensity > min_intensity:
                alldata.append((mz, ii, int(intensity)))

    print("extracted %d valide data points." %len(alldata))
    
    mzTree = {}
    for x in alldata:
        ii = int(x[0]*1000)
        if ii in mzTree:
            mzTree[ii].append(x)
        else:
            mzTree[ii] = [x]

    del alldata
    ks = sorted([k for k,v in mzTree.items() if len(v) >= min_timepoints]) # ascending order
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
        # e.g. in this datafile, 5958 reduced to 4511 traces
        if __rough_check_consecutive_scans__(datatuples):
            good_bins.append(datatuples)
    
    # del mzTree
    return good_bins

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
        # find_peaks returns array([ 1, 8]), {}. distance=3 because it's edge of tolerance
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
