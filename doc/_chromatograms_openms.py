'''
Old functions using pyOpenMS to be converted

in progress

'''

from operator import itemgetter
import numpy as np

from scipy.signal import find_peaks 
from scipy import interpolate
from scipy.ndimage import uniform_filter1d

from statsmodels.nonparametric.smoothers_lowess import lowess

from .mass_functions import check_close_mzs

from pyopenms import MSExperiment, MzMLFile


def sum_dict(dict1, dict2):
    new = {}
    for k in dict2:
        if k in dict1:
            new[k] = dict1[k] + dict2[k]
        else:
            new[k] = dict2[k]
    return new

# -----------------------------------------------------------------------------
# indexing function
# -----------------------------------------------------------------------------

def get_thousandth_regions(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
    '''
    Process all LC-MS spectral data into flexible bins by units of 0.001 amu.
    mz_tolerance_ppm: m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    ms_expt: pyopenms MSExperiment instance, loaded with LC-MS data.
    min_intensity: minimal intentsity value, needed because some instruments keep 0s 
    min_timepoints: minimal consecutive scans to be considered real signal.
    min_peak_height: a bin is not considered if the max intensity < min_peak_height.
    
    Return: a list of flexible bins, [ [(mz, scan_num, intensity_int), ...], ... ]
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

    def __check_min_peak_height__(datatuples, min_peak_height):
        if max([x[2] for x in datatuples]) < min_peak_height:
            return False
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
        if __check_min_peak_height__(datatuples, min_peak_height) and __rough_check_consecutive_scans__(datatuples):
            good_bins.append(datatuples)
    
    # del mzTree
    return good_bins


# -----------------------------------------------------------------------------
# mass Traces
# -----------------------------------------------------------------------------

def extract_massTracks_openms(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
    '''
    A mass track is an EIC for full RT range, without separating the mass traces. 
    ms_expt: pyopenms MSExperiment instance, loaded with LC-MS data.
    return 
    rt_numbers, rt_times,
    tracks as [( mz, rtlist, intensities ), ...]
    '''
    # rt_numbers = range(ms_expt.getNrSpectra())
    rt_times = [spec.getRT() for spec in ms_expt]
    rt_numbers = list(range(len(rt_times)))
    good_bins = get_thousandth_regions(ms_expt, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)
    tracks = []
    for bin in good_bins:
        tracks += bin_to_mass_tracks(bin, mz_tolerance_ppm)
    #
    # merge tracks if m/z overlap
    #
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

def extract_massTraces(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
    '''
    ms_expt: pyopenms MSExperiment instance, loaded with LC-MS data.
    return 
    rt_numbers, rt_times,
    XICs as [( mz, rtlist, intensities ), ...]
    '''
    rt_numbers = range(ms_expt.getNrSpectra())
    rt_times = [spec.getRT() for spec in ms_expt]
    good_bins = get_thousandth_regions(ms_expt, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)
    xics = []
    for bin in good_bins:
        xics += bin_to_xics(bin, mz_tolerance_ppm, min_timepoints)
        
    return {
        'rt_numbers': rt_numbers,
        'rt_times': rt_times,
        'xics': xics,
    }

def bin_to_xics(bin_data_tuples, mz_tolerance_ppm=5, min_timepoints=5, gap_allowed=2):
    '''
    input a flexible bin by units of 0.001 amu, in format of [(mz, scan_num, intensity_int), ...].
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
        peaks, _ = find_peaks(hist, distance=5)
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

def extract_massTracks_(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
    '''
    A mass track is an EIC for full RT range, without separating the mass traces. 
    ms_expt: pyopenms MSExperiment instance, loaded with LC-MS data.
    return 
    rt_numbers, rt_times,
    tracks as [( mz, rtlist, intensities ), ...]
    '''
    # rt_numbers = range(ms_expt.getNrSpectra())
    rt_times = [spec.getRT() for spec in ms_expt]
    rt_numbers = list(range(len(rt_times)))
    good_bins = get_thousandth_regions(ms_expt, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)
    tracks = []
    for bin in good_bins:
        tracks += bin_to_mass_tracks(bin, mz_tolerance_ppm)
    #
    # merge tracks if m/z overlap
    #
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


def extract_single_track_(bin):
    '''
    A mass track is an EIC for full RT range, without separating the mass traces. 
    input bins in format of [(mz_int, scan_num, intensity_int), ...].

    # To consider without zeros to optimize performance -

    For peak detection, we need RT as one continuous trace (RT gaps filled by zeros).
    Peak detection is performed at sample level when landmarks are sought, 
    but CMAP.MassGrid has its own more careful peak deteciton.
    
    return a massTrack as ( mz, rtlist, intensities ).
    '''
    mz = np.mean([x[0] for x in bin])
    # bin.sort(key=itemgetter(1))   # sort by increasing RT (scan number), not needed any more
    rtlist = [x[1] for x in bin]
    min_rt, max_rt = min(rtlist), max(rtlist)
    rtlist = range(min_rt, max_rt+1)    # filling gaps of RT this way if needed, and sorted
    _d = {}                             # dict to hold rt to intensity mapping
    for ii in rtlist: 
        _d[ii] = 0
    for r in bin:                       # this gets max intensity on the same RT scan
        _d[r[1]] = max(r[2], _d[r[1]])
    intensities = [_d[x] for x in rtlist]
    return ( mz, list(rtlist), intensities ) # range object is not desired - use list


def bin_to_mass_tracks(bin_data_tuples, mz_tolerance_ppm=5):
    '''
    input a flexible bin by units of 0.001 amu, in format of [(mz, scan_num, intensity_int), ...].
    A mass track is an EIC for full RT range, without separating the mass traces. 
    An optimization step is carried out in CMAP.MassGrid, to verify m/z separation.

    return massTracks as [( mz, rtlist, intensities ), ...]
    '''
    bin_data_tuples.sort()   # by m/z, ascending
    mz_range = bin_data_tuples[-1][0] - bin_data_tuples[0][0]
    mz_tolerance = bin_data_tuples[0][0] * 0.000001 * mz_tolerance_ppm   
    if mz_range < mz_tolerance:
        return [extract_single_track_(bin_data_tuples)]
    else:
        num_steps = int( 5* mz_range/ mz_tolerance )     # step in 1/5 mz_tolerance
        hist, bin_edges = np.histogram([x[0] for x in bin_data_tuples], num_steps)
        # example hist: array([  8,   33,  11,  24,  31,  50,  81, 164, 269,  28,   7])
        hist_starts = [0] + [hist[:ii].sum() for ii in range(1,num_steps+1)]
        # find_peaks returns e.g. array([ 1, 8]), {}. distance=5 because it's edge of tolerance
        peaks, _ = find_peaks( hist, distance = 5 )
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

def merge_two_mass_tracks(T1, T2):
    '''
    massTracks as [( mz, rtlist, intensities ), ...]
    '''
    mz = 0.5 * (T1[0] + T2[0])
    d_ = sum_dict( dict(zip(T1[1], T1[2])), dict(zip(T2[1], T2[2])) )
    return (mz, list(d_.keys()), list(d_.values()))

