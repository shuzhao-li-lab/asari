'''
Use integers for RT scan numbers and intensities.
Flexible binning based on ppm accuracy.
Use mz tol (default 5 pmm) in XIC construction. 
XICs without neighbors within x ppm are considered specific (i.e. high selectivity). 
Low selectivity regions will be still inspected to determine the true number of XICs.
Leave calibration to Correspondence step.
'''
import numpy as np
from scipy.signal import find_peaks 
from scipy import interpolate
from scipy.ndimage import uniform_filter1d

from statsmodels.nonparametric.smoothers_lowess import lowess

from .mass_functions import check_close_mzs


def sum_dict(dict1, dict2):
    new = {}
    for k in dict2:
        if k in dict1:
            new[k] = dict1[k] + dict2[k]
        else:
            new[k] = dict2[k]
    return new

def get_thousandth_bins(mzTree, mz_tolerance_ppm=5, min_timepoints=5, min_peak_height=1000):
    '''
    Process all LC-MS spectral data into flexible bins by units of 0.001 amu.
    mz_tolerance_ppm: m/z tolerance in part-per-million. Used to seggregate m/z regsions here.
    ms_expt: = pymzml.run.Reader(f)
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


def extract_massTracks_(ms_expt, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
    '''
    A mass track is an EIC for full RT range, without separating the mass traces of same m/z. 
    ms_expt = pymzml.run.Reader(f)
    return 
    rt_numbers, rt_times,
    tracks as [( mz, np.array(intensities at full rt range) ), ...]
    '''
    alldata = []
    rt_times = []           # in seconds
    ii = 0
    for spec in ms_expt:
        if spec.ms_level == 1:                         # MS Level 1 only
            rt_times.append( spec.scan_time_in_minutes()*60 )
            intensities = spec.i.astype(int)
            good_positions = intensities > min_intensity
            intensities = intensities[good_positions]
            mzs = spec.mz[good_positions]
            alldata += [(mz, ii, inten) for mz, inten in zip(mzs, intensities)]
            ii += 1

    #print("extracted %d valide data points." %len(alldata))
    mzTree = {}
    for x in alldata:
        ii = int(x[0]*1000)
        if ii in mzTree:
            mzTree[ii].append(x)
        else:
            mzTree[ii] = [x]

    del alldata
    rt_numbers = list(range(len(rt_times)))
    rt_length = len(rt_numbers)

    good_bins = get_thousandth_bins(mzTree, mz_tolerance_ppm, min_timepoints, min_peak_height)
    tracks = []
    for bin in good_bins:
        tracks += bin_to_mass_tracks(bin, rt_length, mz_tolerance_ppm)

    # merge tracks if m/z overlap
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
# mass Tracks
# -----------------------------------------------------------------------------

def extract_single_track_fullrt_length(bin, rt_length):
    '''
    New format after v1.5.
    A mass track is an EIC for full RT range, without separating the mass traces. 
    input bins in format of [(mz_int, scan_num, intensity_int), ...].
    return a massTrack as ( mz, np.array(intensities at full rt range) ).
    '''
    mz = np.mean([x[0] for x in bin])
    intensity_track = np.zeros(rt_length, dtype=np.int64)
    for r in bin:                       # this gets max intensity on the same RT scan
        intensity_track[r[1]] = max(r[2], intensity_track[r[1]])
    return ( mz, intensity_track )


def extract_single_track_old(bin):
    '''
    Phased out after v1.4.
    A mass track is an EIC for full RT range, without separating the mass traces. 
    input bins in format of [(mz_int, scan_num, intensity_int), ...].
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


def bin_to_mass_tracks(bin_data_tuples, rt_length, mz_tolerance_ppm=5):
    '''
    input a flexible bin by units of 0.001 amu, in format of 
                                    [(mz, scan_num, intensity_int), ...].
    A mass track is an EIC for full RT range, without separating the mass traces. 
    An optimization step is carried out in CMAP.MassGrid, to verify m/z separation.

    return massTracks as ( mz, np.array(intensities at full rt range) )
    '''
    bin_data_tuples.sort()   # by m/z, ascending
    mz_range = bin_data_tuples[-1][0] - bin_data_tuples[0][0]
    mz_tolerance = bin_data_tuples[0][0] * 0.000001 * mz_tolerance_ppm   
    if mz_range < mz_tolerance:
        return [extract_single_track_fullrt_length(bin_data_tuples, rt_length)]
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
                tracks.append( extract_single_track_fullrt_length(
                            bin_data_tuples[hist_starts[left]: hist_starts[right]], rt_length) )
            return tracks
        else:
            peak = np.argmax(hist)  # if find_peaks fails, get position of max value
            left = max(0, peak-2)
            right = min(peak+3, num_steps)
            return [extract_single_track_fullrt_length(bin_data_tuples[hist_starts[left]: hist_starts[right]], rt_length)]

def merge_two_mass_tracks(T1, T2):
    '''
    massTrack as ( mz, np.array(intensities at full rt range) )
    '''
    return ( 0.5 * (T1[0] + T2[0]), T1[1] + T2[1] )

def merge_two_mass_tracks_old(T1, T2):
    '''
    massTrack as ( mz, rtlist, intensities )
    '''
    mz = 0.5 * (T1[0] + T2[0])
    d_ = sum_dict( dict(zip(T1[1], T1[2])), dict(zip(T2[1], T2[2])) )
    return (mz, list(d_.keys()), list(d_.values()))


# -----------------------------------------------------------------------------
# retention time alignment
# -----------------------------------------------------------------------------

def rt_lowess_calibration(good_landmark_peaks, selected_reference_landmark_peaks, sample_rt_numbers, reference_rt_numbers):
    '''
    Use LOWESS, Locally Weighted Scatterplot Smoothing, to reate correspondence between sample_rt_numbers, reference_rt_numbers.
    Predicted numbers are skipped when outside real sample boundaries.
    checked available in statsmodels.nonparametric.smoothers_lowess, v 0.12, 0.13+
        https://www.statsmodels.org/stable/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html    
        But xvals have to be forced as float array until the fix is in new release.

    Input
    =====
    good_landmark_peaks and selected_reference_landmark_peaks are equal-length lists of peaks,
    selected from working sample and reference sample as landmarks for RT alignment. 
    sample_rt_numbers: all scan numbers in this sample.
    reference_rt_numbers: all scan numbers in the ref sample.

    Return
    ======
    rt_cal_dict, dictionary converting scan number in sample_rt_numbers to calibrated integer values.
                Range matched. Only changed numbers are kept for efficiency.
    reverse_rt_cal_dict, from ref RT scan numbers to sample RT scan numbers. 
                Range matched. Only changed numbers are kept for efficiency.
    '''
    # force left and right ends, to prevent runaway curve functions
    reference_rt_bound = max(reference_rt_numbers)
    sample_rt_bound = max(sample_rt_numbers)
    rt_rightend_ = 1.1 * sample_rt_bound
    xx, yy = [0, 0, 0, ], [0, 0, 0, ]
    rt_cal = sorted([(x[0]['apex'], x[1]['apex']) for x in zip(good_landmark_peaks, selected_reference_landmark_peaks)])

    xx += [L[0] for L in rt_cal] + [rt_rightend_]*3
    yy += [L[1] for L in rt_cal] + [rt_rightend_]*3
    
    # This requires statsmodels > v 0.12.
    # float conversion on xvals is to bypass a bug in statsmodels, which was fixed today 2022-01-27
    #lowess_predicted = lowess(yy, xx, frac= .2, it=1, xvals=np.array(sample_rt_numbers, dtype=float)) 
    #
    # downgrade now for compatibility to older statsmodels
    lowess_predicted = __hacked_lowess__(yy, xx, frac= .2, it=1, xvals=sample_rt_numbers)
    interf = interpolate.interp1d(lowess_predicted, sample_rt_numbers, fill_value="extrapolate")
    ref_interpolated = interf( reference_rt_numbers )

    lowess_predicted = [int(round(ii)) for ii in lowess_predicted]
    rt_cal_dict = dict( 
        [(x,y) for x,y in zip(sample_rt_numbers, lowess_predicted) if x!=y and 0<=y<=reference_rt_bound] )
    ref_interpolated = [int(round(ii)) for ii in ref_interpolated]
    reverse_rt_cal_dict = dict(
        [(x,y) for x,y in zip(reference_rt_numbers, ref_interpolated) if x!=y and 0<=y<=sample_rt_bound] )
        
    return rt_cal_dict, reverse_rt_cal_dict


def __hacked_lowess__(yy, xx, frac, it, xvals):
    '''
    This was workaround of a problem with Statsmodel 0.13, which was dealt by casting xvals as floats.
    Not checking possible range error
    '''
    lxy = lowess(yy, xx, frac, it)
    newx, newy = list(zip(*lxy))
    interf = interpolate.interp1d(newx, newy)
    return interf(xvals)


def mock_rt_calibration(sample_rt_numbers, reference_rt_numbers):
    _total = sorted(set(sample_rt_numbers + reference_rt_numbers))
    rt_cal_dict = reverse_rt_cal_dict = dict(zip( _total, _total ))
    return rt_cal_dict, reverse_rt_cal_dict
    

def savitzky_golay_spline(good_landmark_peaks, selected_reference_landmark_peaks, sample_rt_numbers, reference_rt_numbers):
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

def dwt_rt_calibrate(good_landmark_peaks, selected_reference_landmark_peaks, sample_rt_numbers, reference_rt_numbers):
    '''Other RT alignment method to include - pls follow format in rt_lowess.
    not implemented.
    '''
    pass


def remap_intensity_track(intensity_track, new, rt_cal_dict):
    '''
    new = basetrack.copy(), possible longer than intensity_track
    Remap intensity_track based on rt_cal_dict. 
    Can we make this faster?
    '''
    for k,v in rt_cal_dict.items():
        new[v] = intensity_track[k]
    return new



# -----------------------------------------------------------------------------
# smoothing function
# -----------------------------------------------------------------------------

def smooth_moving_average(list_intensity, size=9):
    '''
    Smooth data for very noisy mass tracks.
    Using simple moving average here; LOWESS does not look good for this.
    '''
    return uniform_filter1d(list_intensity, size, mode='nearest')

def smooth_rt_intensity_remap(L_rt_scan_numbers, L_intensity):
    '''
    After remapping/calibratoing RT, smooth intensity curve.
    '''
    _d = dict(zip( L_rt_scan_numbers, L_intensity )) # this overwrites repeated scan numbers by using last intensity value
    newx = range(min(L_rt_scan_numbers), max(L_rt_scan_numbers))
    for ii in newx:
        if ii not in _d:
            _d[ii] = _d[ii-1]
    newy = [_d[ii] for ii in newx]
    # smoothing by averaging 3 consecutive values
    for ii in newx[1: -1]:
        _d[ii] = (_d[ii-1]+_d[ii]+_d[ii+1])/3.0
    return _d

def smooth_lowess(list_intensity, frac=0.02):
    '''
    Smooth data for very noisy mass tracks.
    Using simple moving average here; LOWESS does not look good for this.
    '''
    lxy = lowess(list_intensity, range(len(list_intensity)), frac=frac, it=1)
    _, newy = list(zip(*lxy))
    return newy
