'''
A place to store some prototype or experimental code, 
which may be garbage or useful in someway elsewhere.
'''

# from metDataModel.core import MassTrace, Peak
# from .utils import *


def metafile_to_dict(infile):
    '''
    Optional.
    Meta data file, tab delimited, first two columns corresponding to [file_name, sample_type].
    '''
    meta = {}
    for line in open(infile).read().splitlines():
        a = line.split('\t')
        meta[a[0]] = a[1]
    return {}


def combine_mass_traces(list_mass_traces, mz_tolerance=0.0002):
    '''
    Combine mass traces of same m/z in the input list.
    We bypass this by using mass tracks instead of mass traces.
    
    Input
    =====
    list_mass_traces: likely to contain multiple entries per m/z, due to different retention time.
        They will be combined into unique m/z bins, but float m/z values may have very minor discrepency,
        thus they are verified by comparisons under mz_tolerance.
        mass_trace format: {'id_number': 99, 'mz': 87.0551970014688, 'rt_scan_numbers': [], 'intensity': []}

    Return
    ======
    dict_combined_mass_traces: {index: (mz, id1, id2), ...}
    '''
    mz_mass_traces = sorted([(mt['mz'], mt) for mt in list_mass_traces])
    new = {}  
    ii = 0
    new[ii] = [mz_mass_traces[0][0], mz_mass_traces[0][1]['id_number'], ]
    for jj in range(1, len(mz_mass_traces)):
        if mz_mass_traces[jj][0] - mz_mass_traces[jj-1][0] < mz_tolerance:
            new[ii].append(mz_mass_traces[jj][1]['id_number'])
        else:
            ii += 1
            new[ii] = [mz_mass_traces[jj][0], mz_mass_traces[jj][1]['id_number']]
                
    return new



def anchor_guided_mapping(REF_reference_mzlist, REF_reference_anchor_pairs, 
                            SM_mzlist, SM_anchor_mz_pairs, std_ppm=5, correction_tolerance_ppm=1):
    '''
    Align the mzlists btw CMAP (i.e. REF) and a new Sample,
    prioritizing paired anchors (from isotope/adduct patterns).
    The mzlists are already in ascending order when a Sample is processed.
    Do correciton on list2 if m/z shift exceeds correction_tolerance_ppm.
    The anchor_pairs use corresponding indices.

    Nested indices are: 
        mass_paired_mapping functions return positions of input lists -> 
        which refer to positions in anchor_pairs ->
        which refer to positions in list_mass_tracks or MassGridDict.

    Return
    ======
    new_reference_mzlist, new_reference_map2, _r

    mapped, list1_unmapped, list2_unmapped
    _r: correction ratios on SM_mzlist, to be attached to Sample class instance

    Updated reference list:  because it mixes features from samples 
    and they need to be consistent on how they are calibrated.

    '''
    _N1 = len(REF_reference_mzlist)
    _d2 = {}                                                # tracking how SM_mzlist is mapped to ref
    for ii in range(_N1): 
        _d2[ii] = None
    # first align to reference_anchor_pairs
    anchors_1 = [REF_reference_mzlist[x[0]] for x in REF_reference_anchor_pairs]
    anchors_2 = [SM_mzlist[x[0]] for x in SM_anchor_mz_pairs]

    mapped, ratio_deltas = mass_paired_mapping(anchors_1, anchors_2, std_ppm)
    _r = np.mean(ratio_deltas)
    if abs(_r) > correction_tolerance_ppm*0.000001:          # do m/z correction
        SM_mzlist = [x/(1+_r) for x in SM_mzlist]
        # rerun after mz correction
        anchors_2 = [SM_mzlist[x[0]] for x in SM_anchor_mz_pairs]
        mapped, ratio_deltas = mass_paired_mapping(anchors_1, anchors_2, std_ppm)

    # move onto paired ion in anchors
    anchors_1 = [REF_reference_mzlist[x[1]] for x in REF_reference_anchor_pairs]
    anchors_2 = [SM_mzlist[x[1]] for x in SM_anchor_mz_pairs]
    mapped2, ratio_deltas = mass_paired_mapping(anchors_1, anchors_2, std_ppm)
    # mapped & mapped2 refer to indices in anchor pairs, and have to be converted back to indices of mzlists
    mapped_pairs = [
        ( REF_reference_anchor_pairs[x[0]][0], SM_anchor_mz_pairs[x[1]][0] ) for x in mapped
    ] + [
        ( REF_reference_anchor_pairs[x[0]][1], SM_anchor_mz_pairs[x[1]][1] ) for x in mapped2
    ]
    # move onto remaining ions
    indices_remaining1 = [ii for ii in range(len(REF_reference_mzlist)) if ii not in [x[0] for x in mapped_pairs]]
    indices_remaining2 = [ii for ii in range(len(SM_mzlist)) if ii not in [x[1] for x in mapped_pairs]]
    mapped, list1_unmapped, list2_unmapped = complete_mass_paired_mapping(
            [REF_reference_mzlist[ii] for ii in indices_remaining1], [SM_mzlist[ii] for ii in indices_remaining2], 
            std_ppm)

    mapped_pairs += mapped
    for p in mapped_pairs: 
        _d2[p[0]] = p[1]
    for ii in range(len(list2_unmapped)): 
        _d2[_N1 + ii] = list2_unmapped[ii]
    new_reference_mzlist = REF_reference_mzlist + [SM_mzlist[ii] for ii in list2_unmapped]
    new_reference_map2 = [_d2[x] for x in range(len(new_reference_mzlist))]

    # return mapped_pairs, list1_unmapped, list2_unmapped, _r, REF_reference_mzlist
    return new_reference_mzlist, new_reference_map2, _r



def epd_paired_mapping_with_correction(empCpd_mzlist_1, empCpd_mzlist_2, std_ppm=5, correction_tolerance_ppm=1):
    '''
    Match empCpds as grouped m/z values.
    empCpds:
    [{'id': 358, 'list_peaks': [(4215, 'anchor'), (4231, '13C/12C'), (4339, 'anchor,+NH4')]},...]
    to be converted to 
    # [(id, anchor_peak_mz, other m/z ...), ...]

    Steps:
    1. check anchor m/z matches, not enforcing selectivity here

    2. verify at least another ion matches following the anchor ions
    3. mass correction on list 2 if 


    Return
    ======
    mapped: mapping list [(index from list1, index from list2), ...]

    list1, list2 = [x[1] for x in empCpd_flatlist_1], [x[1] for x in empCpd_flatlist_2]

    mapped, ratio_deltas = mass_paired_mapping(
        
        list1, list2, std_ppm)
    '''
    mapped = []
    # to do
    return mapped



#---------------------------------------------------------------------------------------------------------------

class ext_MassTrace_old(MassTrace):
    '''
    peak detection is in this class

    Extending metDataModel.core.MassTrace
    Peak detection using scipy.signal.find_peaks, a local maxima method with prominence control.
    Not keeping Peaks in this class; Peaks are attached to Sample, then to store in SQLite (next version).


        Useful to have a 2nd rule on prominence; can mark for review during correspondence.
        
    '''
    def __init2__(self, mz, RT, INTS):
        # np.array or list -
        self.mz, self.list_retention_time, self.list_intensity = [mz, RT, INTS]
        # self.cal_list_retention_time = []             # not using now
        self.formula_mass = None
        self.raw_mzs = []                               # multiple values possible when merging traces
        self.__mass_corrected_by_asari__ = False        # True if updated by calibration
        self.features_assigned = False                  # tracking if assigned to Experiment features
        self.sample_name = ''

    def detect_peaks(self, min_intensity_threshold, min_fwhm, min_prominence_threshold, prominence_window, gaussian_shape, snr):
        list_peaks = []
        peaks, properties = find_peaks(self.list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                        prominence=min_prominence_threshold, wlen=prominence_window) 
        # This is noise estimate, but doesn't have to correspond to final list of peaks
        _noise_level_ = self.__get_noise_level__(peaks, properties)
        if peaks.size > 1:
            # Rerun, raising min_prominence_threshold. Not relying on peak shape for small peaks, 
            # because chromatography is often not good so that small peaks can't be separated from noise.
            _tenth_height = 0.1* max(self.list_intensity)
            min_prominence_threshold = max(min_prominence_threshold, _tenth_height)
            peaks, properties = find_peaks(self.list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                        prominence=min_prominence_threshold, wlen=prominence_window)
        
        for ii in range(peaks.size):
            if properties['peak_heights'][ii] > snr*_noise_level_:
                P = ext_Peak()
                # scipy.signal.find_peaks works on 1-D array. RT coordinates need to be added back.
                P.__init2__(parent_mass_trace=self, mz = self.mz, apex = peaks[ii],
                                peak_height = properties['peak_heights'][ii], left_base = properties['left_bases'][ii],
                                right_base = properties['right_bases'][ii])
                P.evaluate_peak_model()
                if P.goodness_fitting > gaussian_shape:
                    list_peaks.append(P)

        return list_peaks

    def extract_targeted_peak(self, rt_range):
        pass

    def __get_noise_level__(self, peaks, properties):
        peak_data_points = []
        for ii in range(peaks.size):
            peak_data_points += range(properties['left_bases'][ii], properties['right_bases'][ii]+1)
        noise_data_points = [ii for ii in range(len(self.list_intensity)) if ii not in peak_data_points]
        if noise_data_points:
            return np.median([self.list_intensity[ii] for ii in noise_data_points])
        else:
            return 0


class ext_Peak_old(Peak):
    '''
    Extending metDataModel.core.Peak, 
    
    including pointer to parent MassTrace. peak evaluation is in this class
    
    Not storing RT or intensity arrays, but looking up in MassTrace if needed.
    '''
    def __init2__(self, parent_mass_trace, mz, apex, peak_height, left_base, right_base):
        [ self.parent_mass_trace, 
        self.mz, self.apex, self.peak_height, self.left_base, self.right_base
                ] = [ parent_mass_trace, mz, apex, peak_height, left_base, right_base ]

        self.left_rtime = float(self.parent_mass_trace.list_retention_time[left_base])
        self.right_rtime = float(self.parent_mass_trace.list_retention_time[right_base])

    def evaluate_peak_model(self):
        '''
        Use Gaussian models to fit peaks, R^2 as goodness of fitting.
        Peak area is defined by Gaussian model, the integral of a Gaussian function being a * c *sqrt(2*pi).
        Good: Left to right base may not be full represenation of the peak. The Gaussian function will propose a good boundary.
        Less good: peaks are not always in Gaussian shape. But we are comparing same features across samples, 
        same bias in peak shape is applied to all samples.
        '''
        xx = self.parent_mass_trace.list_retention_time[self.left_base: self.right_base+1]
        # set initial parameters
        a, mu, sigma =  self.peak_height, \
                        self.parent_mass_trace.list_retention_time[self.apex], \
                        np.std(xx)
        try:
            popt, pcov = curve_fit(gaussian_function__, 
                            xx, self.parent_mass_trace.list_intensity[self.left_base: self.right_base+1],
                            p0=[a, mu, sigma])
            self.gaussian_parameters = popt
            self.rtime = popt[1]
            self.peak_area = popt[0]*abs(popt[2])*2.506628274631        # abs because curve_fit may return negative sigma
            self.goodness_fitting = goodness_fitting__(
                            self.parent_mass_trace.list_intensity[self.left_base: self.right_base+1], 
                            gaussian_function__(xx, *popt))

        # failure to fit
        except RuntimeError:
            self.rtime = mu
            self.peak_area = a*sigma*2.506628274631
            self.goodness_fitting = 0

    def extend_model_range(self):
        # extend model xrange, as the initial peak definition may not be complete, for visualization
        _extended = self.right_base - self.left_base
        # negative index does not work here, thus max to 0
        self.rt_extended = self.parent_mass_trace.list_retention_time[
                                            max(0, self.apex-_extended): self.apex+_extended]
        self.y_fitted_extended = gaussian_function__(self.rt_extended, *self.gaussian_parameters)


#---------------------------------------------------------------------------------------------------------------

class Sample(SimpleSample):
    '''
    LC-MS Sample with fuller functions on mass traces and peak detections.
    This can be used externally for data mining, inferring m/z relationships (i.e. to construct epdTrees), etc.
    '''
        
    def get_masstraces(self, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5):
        '''
        Get chromatograms, this class is using asari.chromatograms algorithm.
        Use another class for a different algorithm (e.g. Sample_openms).
        XICs as [( mz, rtlist, intensities ), ...].
        mass_trace format: {'id_number':0, 'mz', 'rt_scan_numbers', 'intensity', 'number_peaks':0}
        
        Multiple mass tracess can have the same m/z values, i.e. they belong to the same `mass track`.
        '''
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xic_dict = extract_massTraces(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, min_intensity=min_intensity, min_timepoints=min_timepoints)
        self.rt_numbers = xic_dict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xic_dict['rt_times']     # full RT time points in sample
        ii = 0
        for xic in xic_dict['xics']:                        # xic = ( mz, rtlist, intensities )
            self.list_mass_traces.append( {
                'id_number': ii, 
                'mz': xic[0],
                'rt_scan_numbers': xic[1],                  # list
                'intensity': xic[2],                        # list
                # 'number_peaks': 0,
                } )
            ii += 1

        print("Processing %s, found %d mass traces." %(self.input_file, ii))

    def get_peaks(self, min_intensity_threshold=10000, min_fwhm=3, min_prominence_threshold=5000, snr=2):
        '''
        Initial elution peak detection; will do another round of peak detection on compositeMap.
        '''
        list_peaks = []
        for xic in self.list_mass_traces:
            list_peaks += self.detect_peaks(xic, min_intensity_threshold, min_fwhm, min_prominence_threshold, snr)
        ii = 0
        for p in list_peaks:
            p['id_number'] = ii
            self.list_peaks.append(p)
            ii += 1

    def detect_peaks(self, mass_trace, min_intensity_threshold=10000, min_fwhm=3, min_prominence_threshold=5000, snr=2):
        '''
        Return list of peaks. No stringent filtering of peaks. 
        Reported left/right bases are not based on Gaussian shape or similar, just local maxima.
        Prominence has a key control of how peaks are considered, but peak shape evaluation later can filter out most bad peaks.
        Peak format: {
            'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 'height', 'parent_masstrace_id', 
            'rtime', 'peak_area', 'goodness_fitting'
        }
        '''
        list_peaks = []
        rt_numbers, list_intensity = mass_trace['rt_scan_numbers'], mass_trace['intensity']
        # 
        prominence = max(min_prominence_threshold, 0.05*max(list_intensity))
        peaks, properties = find_peaks(list_intensity, height=min_intensity_threshold, width=min_fwhm, 
                                                        prominence=prominence) 
        _noise_level_ = self.get_noise_level__(list_intensity, peaks, properties)
        for ii in range(peaks.size):
            if properties['peak_heights'][ii] > snr*_noise_level_:
                list_peaks.append({
                    'parent_masstrace_id': mass_trace['id_number'],
                    'mz': mass_trace['mz'],
                    'apex': rt_numbers[peaks[ii]], 
                    'height': properties['peak_heights'][ii],
                    'left_base': rt_numbers[properties['left_bases'][ii]],
                    'right_base': rt_numbers[properties['right_bases'][ii]],
                })
        return list_peaks

    def get_noise_level__(self, list_intensity, peaks, properties):
        peak_data_points = []
        for ii in range(peaks.size):
            peak_data_points += range(properties['left_bases'][ii], properties['right_bases'][ii]+1)
        noise_data_points = [ii for ii in range(len(list_intensity)) if ii not in peak_data_points]
        if noise_data_points:
            return np.mean([list_intensity[ii] for ii in noise_data_points])        # mean more stringent than median
        else:
            return 0

    def get_empCompounds(self):
        '''
        update self.list_empCpds
        e.g. [{'id': 358, 'list_peaks': [(4215, 'anchor'), (4231, '13C/12C'), (4339, 'anchor,+NH4')]}, ...]
        '''
        ECCON = epdsConstructor(self.list_peaks)
        self.list_empCpds = ECCON.peaks_to_epds()

    def export_peaklist(self):                                  # for diagnosis etc.
        '''
        in progress -
        '''
        outfile = os.path.join(self.experiment.output_dir, self.name.replace('.mzML', '') + '.peaklist')
        header = ['m/z', 'retention time', 'area', 'shape_quality', 'gaussian_amplitude', 'gaussian_variance', 'mz_selectivity']
        peaklist = []
        for P in self.good_peaks:
            formatted = [str(round(P.mz, 6)), str(round(P.rtime, 2)), str(int(P.peak_area)), str(round(P.goodness_fitting, 2)), 
                            str(int(P.gaussian_parameters[0])), str(round(P.gaussian_parameters[2], 2)), str(round(P.selectivity, 2)),]
            if P.mzstr in self.mzstr_2_formula_mass:
                formatted.append(self.mzstr_2_formula_mass[P.mzstr])
            peaklist.append(formatted)

        with open(outfile, 'w') as O:
            O.write( '\t'.join(header) + '\n' + '\n'.join([ '\t'.join(L) for L in peaklist ]) + '\n' )

    def export_mass_traces(self):
        '''
         for diagnosis etc.
        '''
        outfile = os.path.join(self.experiment.output_dir, 
                            os.path.basename(self.input_file).replace('.mzML', '.peaklist') )
        header = ['m/z', 'retention time', 'area', 'shape_quality', 'gaussian_amplitude', 'gaussian_variance', 'mz_selectivity']
        peaklist = []
        for P in self.list_mass_tracks:
            peaklist.append(str(P))

        with open(outfile, 'w') as O:
            O.write( '\t'.join(header) + '\n' + '\n'.join( peaklist ) + '\n' )


    #---------------------------------------------------------------------------------------------------------------

    def convert_empCpd_mzlist(self, SM):
        '''
        Convert [{'id': 358, 'list_peaks': [(4215, 'anchor'), (4231, '13C/12C'), (4339, 'anchor,+NH4')]}, ...]
        to [{'id': 358, 'mz_peaks': [mz1, mz2, ...]}, ...]
        '''
        new = []
        for EC in SM.list_empCpds:
            new.append(
                {'id': EC['id'], 'mz_peaks': [SM.list_peaks[x[0]]['mz'] for x in EC['list_peaks']]}
            )
        return new



#---------------------------------------------------------------------------------------------------------------


class fullSample(SimpleSample):

    def process(self, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height):
        '''From input file to list_MassTraces with detected peaks and selectivity on peaks.
        '''
        self.generate_mass_tracks( mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)

    def generate_mass_tracks(self, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
        '''
        A mass track is an EIC for full RT range, without separating the mass traces,
        using asari.chromatograms algorithm.
        tracks as [( mz, rtlist, intensities ), ...].
        '''
        list_mass_tracks = []

        exp = pymzml.run.Reader(self.input_file)
        
        # exp = MSExperiment()                                                                                          
        # MzMLFile().load(self.input_file, exp)

        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    min_intensity=min_intensity, 
                    min_timepoints=min_timepoints, 
                    min_peak_height=min_peak_height)
        self.rt_numbers = xdict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xdict['rt_times']     # full RT time points in sample
        ii = 0
        # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
        for track in xdict['tracks']:                         
            list_mass_tracks.append( {
                'id_number': ii, 
                'mz': track[0],
                'rt_scan_numbers': track[1],                  # list
                'intensity': track[2],                        # list
                } )
            ii += 1

        print("Processing %s, found %d mass tracks." %(os.path.basename(self.input_file), ii))
        self.generate_anchor_mz_pairs(list_mass_tracks)
        print("    Number of anchor m/z pairs = %d" %self._number_anchor_mz_pairs_)

        if self.database_mode == 'memory':
            self.list_mass_tracks = list_mass_tracks
        elif self.database_mode == 'ondisk': 
            self.push_to_disk(list_mass_tracks)
        else:                         # requiring self.experiment and db connection
            self.push_to_db(list_mass_tracks
                # to implement
            )

    def generate_mass_tracks_openms(self, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
        '''
        A mass track is an EIC for full RT range, without separating the mass traces,
        using asari.chromatograms algorithm.
        tracks as [( mz, rtlist, intensities ), ...].
        '''
        list_mass_tracks = []
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    min_intensity=min_intensity, 
                    min_timepoints=min_timepoints, 
                    min_peak_height=min_peak_height)
        self.rt_numbers = xdict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xdict['rt_times']     # full RT time points in sample
        ii = 0
        # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
        for track in xdict['tracks']:                         
            list_mass_tracks.append( {
                'id_number': ii, 
                'mz': track[0],
                'rt_scan_numbers': track[1],                  # list
                'intensity': track[2],                        # list
                } )
            ii += 1

        print("Processing %s, found %d mass tracks." %(os.path.basename(self.input_file), ii))
        self.generate_anchor_mz_pairs(list_mass_tracks)
        print("    Number of anchor m/z pairs = %d" %self._number_anchor_mz_pairs_)

        if self.database_mode == 'memory':
            self.list_mass_tracks = list_mass_tracks
        elif self.database_mode == 'ondisk': 
            self.push_to_disk(list_mass_tracks)
        else:                         # requiring self.experiment and db connection
            self.push_to_db(list_mass_tracks
                # to implement
            )

    def push_to_disk(self, list_mass_tracks):
        '''Write pickle file for list_mass_tracks under project directory, pickle/samplename.pickle.
        '''
        with open(self.pickle_file, 'wb') as f:
            pickle.dump(list_mass_tracks, f, pickle.HIGHEST_PROTOCOL)


#
# -----------------------------------------------------------------------------
#
# not used now
# 


class Sample_old:
    '''
    A sample or injection, corresponding to one raw data file, either from positive or negative ionization.
    Each sample has a series of scans in LC-MS (retentiion time). 
    The RT is recorded as integers of scan number in asari, and converted to seconds on demand.
    The OpenMS workflow uses seconds directly, and users should class Sample_openms.

    RT and m/z values in a sample may be modified due to alignment to other samples.
    The correction functions are stored per sample.


    '''
    def __init__(self, experiment, mode='pos', input_file=''):
        self.input_file = input_file
        self.experiment = experiment                        # parent Experiment instance
        self.name = os.path.basename(input_file)            # must be unique
        self.mode = self.experiment.mode
        self.parameters = self.experiment.parameters
        self.dict_masstraces = {}                           # indexed by str(round(mz,6))
        self.mzstr_2_formula_mass = {}
        self.good_peaks = []

        self.__valid__ = True
        self.__mass_accuracy__, self.__mass_stdev__ = None, None

        # to apply correction/calibration 
        self.__mass_ppm_shift__ = None                      # _r in utils.mass_paired_mapping_with_correction
        self.__rt_calibration__ = None                      # This will be the calibration function

        self.rt_numbers = []                                # list of scans, starting from 0
        self.list_retention_time = []                       # full RT time points in sample

        self.list_mass_traces = []
        self.dict_mass_traces = {}
        self.list_peaks = []
        self.dict_peaks = {}

    def process_step_1(self):
        '''
        From input file to list_MassTraces with detected peaks and selectivity on peaks.
        '''
        self._get_masstraces_()
        self._detect_peaks_()
        self._assign_selectivity_()

    def process_step_2(self, DFDB):
        '''
        Annotate MassTraces with unique formula_mass. This is based on reference DB and may not cover all MassTraces.
        The remaining steps (correspondence, alignment) will be performed at Experiment level.
        '''
        self._match_mass_formula_(DFDB)
        for P in self.good_peaks:
            P.sample_name = self.name

    def export_peaklist(self):                                  # for diagnosis etc.
        outfile = os.path.join(self.experiment.output_dir, self.name.replace('.mzML', '') + '.peaklist')
        header = ['m/z', 'retention time', 'area', 'shape_quality', 'gaussian_amplitude', 'gaussian_variance', 'mz_selectivity']
        peaklist = []
        for P in self.good_peaks:
            formatted = [str(round(P.mz, 6)), str(round(P.rtime, 2)), str(int(P.peak_area)), str(round(P.goodness_fitting, 2)), 
                            str(int(P.gaussian_parameters[0])), str(round(P.gaussian_parameters[2], 2)), str(round(P.selectivity, 2)),]
            if P.mzstr in self.mzstr_2_formula_mass:
                formatted.append(self.mzstr_2_formula_mass[P.mzstr])
            peaklist.append(formatted)

        with open(outfile, 'w') as O:
            O.write( '\t'.join(header) + '\n' + '\n'.join([ '\t'.join(L) for L in peaklist ]) + '\n' )

    def create_peak_dict(self):
        dict_peaks = {}                                # index Peaks by by mzstr
        for P in self.good_peaks:
            if P.mzstr in dict_peaks:
                dict_peaks[P.mzstr].append(P)
            else:
                dict_peaks[P.mzstr] = P
        return dict_peaks

    def _detect_peaks_(self):
        for mzstr, ML in self.dict_masstraces.items():
            for M in ML:
                list_peaks = M.detect_peaks(self.parameters['min_intensity_threshold'], int(0.5 * self.parameters['min_timepoints']), 
                                self.parameters['min_prominence_threshold'], self.parameters['prominence_window'], self.parameters['gaussian_shape'],
                                self.parameters['signal_noise_ratio'],)
                if list_peaks:
                    for P in list_peaks:
                        P.mzstr = mzstr
                        self.good_peaks.append(P)

        print("Detected %d good peaks." %len(self.good_peaks))

    def _get_masstraces_(self):
        '''
        Get chromatograms, this class is using asari.chromatograms algorithm.
        Use another class for a different algorithm (e.g. Sample_openms).
        '''
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xic_dict = self._get_masstraces_from_centroided_rawdata_(exp)
        self.rt_numbers = xic_dict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xic_dict['rt_times']     # full RT time points in sample
        for xic in xic_dict['xics']:                        # ( mz, rtlist, intensities )
            _low, _high = self.parameters['mass_range']
            if _low < xic[0] < _high:
                mz_str = str(round(xic[0],4))
                [RT, INTS] = xic[1:] 
                # convert RT from scan number to seconds
                RT = [self.list_retention_time[ii] for ii in RT]
                if max(INTS) > self.experiment.parameters['min_intensity_threshold']:
                    M = ext_MassTrace()
                    M.__init2__(xic[0], RT, INTS)
                    if mz_str in self.dict_masstraces:
                        self.dict_masstraces[mz_str].append(M)
                    else:
                        self.dict_masstraces[mz_str] = [M]

        print("Processing %s, found %d mass traces." %(self.input_file, len(self.dict_masstraces)))

    def _get_masstraces_from_centroided_rawdata_(self, exp):
        '''
        Get chromatograms, using asari.chromatograms algorithm.
        XICs as [( mz, rtlist, intensities ), ...].
        Only using pyOpenMS to parse mzML centroided raw data. This can be replaced by any other good parser. 
        '''
        return extract_massTraces(exp,  mz_tolerance_ppm=self.experiment.parameters['mz_tolerance'], 
                                        min_intensity=100, 
                                        min_timepoints=self.experiment.parameters['min_timepoints'])

    def _assign_selectivity_(self, std_ppm=5):
        Q = [(P.mz, P.rtime, P) for P in self.good_peaks]                                       # need rtime to break ties in sorting
        Q.sort()
        selectivities = calculate_selectivity([x[0] for x in Q], std_ppm)
        for ii in range(len(Q)):
            Q[ii][2].selectivity = selectivities[ii]

    def _match_mass_formula_(self, DFDB, check_mass_accuracy=True, ppm=10):
        '''
        Match peaks to formula_mass database, including mass accuracy check and correction if needed.
        Only high-selectivity (> 0.9) peaks will be used downstream for alignment calibration.
        Mass accuracy has no hard limit at the moment, but should be restricted in future versions.
        E.g. if mass shift is larger than 20 ppm, user should do own mass calibrate first.

        DFDB: a reference DB in DataFrame.
        check_mass_accuracy should be on for all samples in asari processing.
        ppm: default to 10, because features (although in smaller number) will still match if the instrument has larger variations.
        '''
        list_ppm_errors, mzstr_2_formula_mass = [], {}
        for P in self.good_peaks:
            query = search_formula_mass_dataframe(P.mz, DFDB, ppm)
            if query:
                _formula_mass, delta_ppm = query
                mzstr_2_formula_mass[P.mzstr] = _formula_mass
                list_ppm_errors.append(delta_ppm)

        self.__mass_accuracy__, self.__mass_stdev__ = mass_accuracy, ppm_std = normal_distribution.fit(list_ppm_errors)
        print("    (mass_accuracy, stdev in matched results) = (%5.2f, %5.2f ) ppm." %(mass_accuracy, ppm_std))
        if check_mass_accuracy and abs(mass_accuracy) > 5:   # this is considered significant mass shift, requiring m/z correction for all
            list_ppm_errors, mzstr_2_formula_mass = [], {}
            for P in self.good_peaks:
                P.mz = P.mz - P.mz*0.000001*mass_accuracy
                # redo search because we may find different matches after mass correction; update search ppm too
                query = search_formula_mass_dataframe(P.mz, DFDB, 2*ppm_std)
                if query:
                    _formula_mass, delta_ppm = query
                    mzstr_2_formula_mass[P.mzstr] = _formula_mass
                    list_ppm_errors.append(delta_ppm)

            # update ppm_std because it may be used for downstream search parameters
            _mu, self.__mass_stdev__ = normal_distribution.fit(list_ppm_errors)
            for ML in self.dict_masstraces.items():
                for M in ML:
                    M.raw_mzs = [M.mz]      # may not need to be a list anymore
                    M.mz = M.mz - M.mz*0.000001*mass_accuracy
                    M.__mass_corrected_by_asari__ = True

        self.mzstr_2_formula_mass = mzstr_2_formula_mass



class Sample_openms(Sample):
    '''
    Reusing Sample class; only modifying mass trace functions as RT is used differently.

    in progress

    '''


    def _get_masstraces_(self):
        '''
        Get chromatograms, default using asari.chromatograms algorithm.

        # diff algorithms on chromatogram construction
        if algorithm == 'asari':
            xic_dict = self._get_masstraces_from_centroided_rawdata_(exp)
        elif algorithm == 'openms':
            xic_dict = self._get_masstraces_from_chromatogram_file_(exp)
        else:
            raise ValueError("New chromatogram algorithm? Supporting asari and openms now.")


        '''

        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xic_dict = self._get_masstraces_from_chromatogram_file_(exp)


        self.rt_numbers = xic_dict['rt_numbers']
        self.list_retention_time = xic_dict['rt_times']
        for xic in xic_dict['xics']:             # ( mz, rtlist, intensities )
            _low, _high = self.parameters['mass_range']
            if _low < xic[0] < _high:
                mz_str = str(round(xic[0],4))
                [RT, INTS] = xic[1:] 
                if max(INTS) > self.experiment.parameters['min_intensity_threshold']:
                    M = ext_MassTrace()
                    M.__init2__(xic[0], RT, INTS)
                    if mz_str in self.dict_masstraces:
                        self.dict_masstraces[mz_str].append(M)
                    else:
                        self.dict_masstraces[mz_str] = [M]

        print("\nProcessing %s, found %d mass traces." %(self.input_file, len(self.dict_masstraces)))


    def _get_masstraces_from_chromatogram_file_(self, exp):
        '''
        Get chromatograms, via pyOpenMS functions (only used in this function). 
        An mzML file of chromatograms (i.e. EIC or XIC, or mass trace). Currently using a MassTrace file produced by FeatureFinderMetabo (OpenMS). 
        min_intensity_threshold: minimal intensity requirement for a mass trace to be considered. Default value 10,000 is lenient in our Orbitrap data.
        Updates self.dict_masstraces, a dict of list of MassTrace objects, indexed by str(round(mz,4))

        list_retention_time from OpenMS is likely list of numbers in seconds.



        '''

        for chromatogram in exp.getChromatograms():
            mz = chromatogram.getPrecursor().getMZ()
            _low, _high = self.parameters['mass_range']
            if _low < mz < _high:
                mz_str = str(round(mz,4))
                RT, INTS = chromatogram.get_peaks() 
                if INTS.max() > self.experiment.parameters['min_intensity_threshold']:
                    # * 0.1:    # can lower if want more XICs for weak signal recovery, but higher threshold still for good peaks
                    M = ext_MassTrace()
                    M.__init2__(mz, RT, INTS)
                    if mz_str in self.dict_masstraces:
                        self.dict_masstraces[mz_str].append(M)
                    else:
                        self.dict_masstraces[mz_str] = [M]

        return xics




#
# -----------------------------------------------------------------------------
#
# not used now

class ext_Experiment_old(Experiment):

    def calibrate_retention_time(self, method='spline', smoothing_factor=0.5):
        '''
        Calibrate (align) RT using selected `good` peaks.
        Overlay all samples to take median values as reference RT for each peak.

        method:spline:
        Do Spline regression of each sample against the reference (composite),
        and apply the spline function to all RTs in the sample as calibrated RT.
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html

        method:dtw:
        https://dynamictimewarping.github.io/
        Will compare with spline later, and implement if desired (?).    

        '''
        rt_table = self.get_rt_calibration_ref()    # This is the pd.DataFrame containing peak data for RT calibration
        rt_table['median'] = rt_table.median(axis=1)
        # rt_table = rt_table.sort_values(by='median')
        # rtbins = np.linspace(0, self.max_rtime, 11)
        for SM in self.samples:
            rt_cal = rt_table[[SM.name, 'median']].dropna(axis=0, how='any').values.tolist() 
            # now this is listed converted from numpy.ndarray 
            if len(rt_cal) < self.parameters['peak_number_rt_calibration']:
                SM.__valid__ = False 
                print("\n\n*** Warning, RT regression using too few features (%d) ***" %len(rt_cal))
                print("*** Sample %s removed from processing. ***\n\n" %SM.name)
            else:
                rt_cal.sort()
                # to-do: need down sample to spread better over rt range
                xx, yy = [0, ], [0, ]                   # force left to 0
                for L in rt_cal:
                    if abs(L[0]-L[1])/L[1] < 0.2:       # control shift < 20%
                        xx.append(L[0])
                        yy.append(L[1])

                # force right match, to avoid erradic spline running out of sanity
                right_end = 1.1 * max( self.parameters['max_rtime'], L[0], L[1] )
                xx.append(right_end)
                yy.append(right_end)
                # leave out s=smoothing_factor
                spl = UnivariateSpline(xx, yy, )
                SM.__rt_calibration__ = spl
                SM.__rt_calibration__data__ = (xx, yy)

                # calibrate all detected RT for all peaks, and raw RT from MassTraces
                for P in SM.good_peaks:
                    P.cal_rtime = SM.__rt_calibration__(P.rtime)
                    P.left_rtime = SM.__rt_calibration__(P.left_rtime)
                    P.right_rtime = SM.__rt_calibration__(P.right_rtime)

    def get_rt_calibration_ref(self):
        '''
        Get N good peaks per sample, single peaks in mass_trace with R^2 > 0.9;
        consistently present in most samples (70%).
        RT variation less than 20% of RT range. !!! important !!! - because we don't fix if chromatography is too bad.
        No imputation solution for a peak missing in a sample. Because one can't assume the peak RT has no shift when the next peak may do.
        Therefore, it's better to skip the peak altogether for that sample.

        return
        ------
        A pd.DataFrame with RT values as reference for selected features.
        '''
        # get reference features for RT calibration/alignment
        d = {}
        for SM in self.samples:
            # SM.export_peaklist()  # test
            good_peaks_rtime, good_peaks_formula_mass = [], []
            for P in SM.good_peaks:                                    # selectivity rules out redundant formulae
                if P.mzstr in SM.mzstr_2_formula_mass and P.selectivity > 0.98 and P.goodness_fitting > 0.9:
                    good_peaks_rtime.append( P.rtime )
                    good_peaks_formula_mass.append(SM.mzstr_2_formula_mass[P.mzstr])

            d[SM.name] = pd.Series(good_peaks_rtime, index=good_peaks_formula_mass)
        
        rt_table = pd.DataFrame(d)      # merge into a table, each row as feature, col as sample
        # avoiding pd.DataFrame whenever possible, unpredictable behaivors
        # drop rows by min presence in > 50% of samples. Not significnat, but potentially tricky for very large studies. 
        # QC or blank samples may behave differently
        rt_table = rt_table.dropna(axis=0, thresh=min(int(0.5 * self.number_of_samples), 10))
        rt_table.to_csv("raw_rt_calibration_matrix.tsv", sep="\t")    # to export 
        return rt_table
        
    def correspondence(self):
        '''
        In each sample: Peak.mzstr links to Sample.mzstr_2_formula_mass, Sample.dict_masstraces
        Start feature table using good peaks (quality > 0.8), then fill weak peaks based on them. 
        Because no feature should be considered if no single good peak exists.
        To do: detailed peak info will be pushed in to SQLite DB.
        '''
        self.samples = [SM for SM in self.samples if SM.__valid__]      # remove sample!!!
        self.ordered_sample_names = [SM.name for SM in self.samples]    # used to order intensity values in Features 
        unassigned = []
        peak_dict = {}
        for SM in self.samples:
            self.name_to_Sample[SM.name] = SM
            for P in SM.good_peaks:
                #if P.goodness_fitting > 0.9:
                if P.mzstr not in SM.mzstr_2_formula_mass:              # next to get a label of consensus m/z
                    unassigned.append((P.mz, P.rtime, P))               # need rtime to break ties in sorting
                else:                                                   # those with formula_mass labelled
                    k = SM.mzstr_2_formula_mass[P.mzstr]
                    if k in peak_dict:
                        peak_dict[k].append(P)
                    else:
                        peak_dict[k] = [P]

        if unassigned:
            unassigned.sort()
            unassigned = [(x[0], x[2]) for x in unassigned]
            mz_peak_bins = bin_by_median(unassigned, lambda x: 2 * self.__mass_stdev__ * 0.000001 * x)
            for BIN in mz_peak_bins:
                peak_dict[ '_M_' + str(round(np.median([P.mz for P in BIN]),6)) ] = BIN
        FeatureList = peaks_to_features(peak_dict, self.parameters['rtime_tolerance'], self.ordered_sample_names)
        print("Additional features are assembled based on 2x stdev (%5.2f ppm) seen in this experiment, " % self.__mass_stdev__)
        self.FeatureTable = FeatureList
        # to update selectivity_combined



    #---------------------------------------------------------------------------------------------------------------

    def __obsolete__process_all(self):
        '''
        This will shift to a DB design in next version.
        '''
        self.init_hot_db( self._get_ref_db_() )                         # initial processing of 3 samples to set up HOT_DB
        for f in self.list_input_files:                                 # run remaining samples
            if f not in self.initiation_samples:
                SM = Sample(self, self.mode, f)
                SM.process_step_1()
                SM.process_step_2(self.HOT_DB)
                if not self.parameters['cache_mass_traces']:
                    del(SM.dict_masstraces)
                self.samples.append(SM)
        
        self.calibrate_retention_time()                                 # samples may be marked to drop
        self.correspondence()
        self.annotate_final()
        self.export_feature_table(self.FeatureTable, self.parameters['output_filename'])
        
    def _get_ref_db_(self):
        '''
        Dispatch for ref DB.
        Earlier version used INIT_DFDB = DB_to_DF( extend_DB1(DB_1) ), which was moved to mass2chem.
        '''
        if self.mode == 'pos':
            dbfile = os.path.join(os.path.dirname(__file__), 'ref_db_v0.2.tsv')
        elif self.mode == 'neg':
            dbfile = os.path.join(os.path.dirname(__file__), 'neg_ref_db_v0.2.tsv')
        else:
            print("Ionization mode is either `pos` or `neg`.")
        return tsv2refDB(dbfile)


    def init_hot_db(self, DFDB):
        '''
        Use three samples to initiate a hot DB to house feature annotation specific to this Experiment, and speed up subsequent search.
        The HOT_DB will be used during sample processing, and have another update after correspondence and additional annotation.
        The HOT_DB will then be exported as Expt annotation.
        '''
        chosen_Samples, found_formula_masses = [], []
        for f in self.initiation_samples:
            SM = Sample(self, self.mode, f)
            SM.process_step_1()
            SM.process_step_2(DFDB)
            chosen_Samples.append(SM)
            found_formula_masses += list(SM.mzstr_2_formula_mass.values())

        self.samples += chosen_Samples
        # Experiment wide parameters
        self.__mass_stdev__ = np.median([SM.__mass_stdev__ for SM in chosen_Samples])       # ppm stdev, used for later searches
        # ver 1 HOT_DB will be a subset of INIT_DFDB

        found_formula_masses = set(found_formula_masses)
        if None in found_formula_masses:
            found_formula_masses.remove(None)
        self.HOT_DB = DFDB.loc[found_formula_masses]   

        print("\n[@.@] Anchoring with %d initial formula matches." %len(found_formula_masses))
        print("[@.@] Initial estimation done on\n" + '\n'.join(self.initiation_samples))
        #export_hot_db, without last col
        self.HOT_DB.iloc[:, :-1].to_csv(os.path.join(self.output_dir, '__intermediary__' + self.parameters['annotation_filename']), sep="\t")
        

    def annotate_final(self):
        '''
        More formula annotation of _M_ features using HMDB+PubChemLite;
        Group into empCpds via mass2chem.


        Still to do remaining features
        '''
        s = u'\t'.join(['feature_id', 'formula_mass', 'mz_dbrecord',	'intensity_mean', 'charged_formula', 'selectivity',	'neutral_formula_mass',
                                    'ion_relation', 'id_HMDB', 'name']) + '\n'
        for F in self.FeatureTable:
            if "_M_" == F.mass_id[:3]:
                s += u'\t'.join([F.feature_id, F.mass_id, str(round(F.mz,4)), str(F.intensity_mean)]) + '\n'
            else:
                [mz, charged_formula, selectivity, neutral_formula_mass, ion_relation] = [str(x) for x in list(self.HOT_DB.loc[F.mass_id])[:5]]
                name = massDict_hmdb.get(neutral_formula_mass, '')
                if name:
                    name = u'\t'.join( [';'.join(x) for x in name] ).encode('utf-8', 'ignore').decode('utf-8')
                s += u'\t'.join([F.feature_id, F.mass_id, mz, str(F.intensity_mean),
                                charged_formula, selectivity, neutral_formula_mass, ion_relation, name]) + '\n'
                
        with open(os.path.join(self.output_dir, self.parameters['annotation_filename']), 'w', encoding='utf-8') as O:
            O.write(s.encode('utf-8', 'ignore').decode('utf-8'))

        
    def _reformat_epds_(self, list_empCpds, FeatureList):
        '''
        usage: list_empCpds = self._reformat_epds_(list_empCpds, self.CMAP.FeatureList)
        '''
        fDict = {}
        for F in FeatureList:
            fDict[F['id_number']] = F
        new = []
        for E in list_empCpds:
            features = []
            for peak in E['list_peaks']:
                features.append(
                    {'feature_id': peak[0], 
                    'mz': fDict[peak[0]]['mz'], 
                    'rtime': fDict[peak[0]]['apex'], 
                    'charged_formula': '', 
                    'ion_relation': peak[1]}
                )
            new.append(
                {
                'interim_id': E['id'], 
                'neutral_formula_mass': None,
                'neutral_formula': None,
                'Database_referred': [],
                'identity': [],
                'MS1_pseudo_Spectra': features,
                'MS2_Spectra': [],
                }
            )
        return new

    def export_feature_table_old__(self, FeatureList, outfile='feature_table.tsv'):
        '''
        FeatureList: a list of namedTuples, i.e. Features; Output two files, one main, another low quality features.
        '''
        def __write__(FeatureList, outfile):
            s = '\t'.join(['feature_id', 'formula_mass', 'mz', 'rtime', 'rt_min', 'rt_max', 'number_peaks',
                                    'peak_quality_max', 'peak_quality_median', 'intensity_mean', 'selectivity_mz',
                                    ] + self.ordered_sample_names) + '\n'
            for F in FeatureList:
                s += '\t'.join(
                    [F.feature_id, F.mass_id, str(round(F.mz,4)), str(round(F.rtime,2)), str(round(F.rt_min,2)), str(round(F.rt_max,2)), str(F.number_peaks),
                    str(round(F.peak_quality_max,2)), str(round(F.peak_quality_median,2)), str(F.intensity_mean), str(round(F.selectivity_mz,2)),
                    ] + [str(int(x)) for x in F.intensities]
                    ) + '\n'
            with open( outfile, 'w') as O:
                O.write(s)

        high_quality_features, low_quality_features = [], []
        for F in FeatureList: 
            if F.peak_quality_max > 0.8 and F.perc_peaks > 15:
                high_quality_features.append(F)
            else:
                low_quality_features.append(F)

        high_quality_features.sort(key=lambda F: F.peak_quality_max, reverse=True)
        #print(FeatureList[99])
        __write__(high_quality_features, os.path.join(self.output_dir, outfile))
        __write__(low_quality_features, os.path.join(self.output_dir, 'low_quality_features_' + outfile))
        print("Feature tables were written under %s." %self.output_dir)
        print("The main feature table (%s) has %d samples and %d features.\n\n\n" %(
            self.parameters['output_filename'], len(self.ordered_sample_names), len(high_quality_features)))


#---------------------------------------------------------------------------------------------------------------










