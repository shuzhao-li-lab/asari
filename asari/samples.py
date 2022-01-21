'''
asari, LC-MS metabolomics data preprocessing - trackable, scalable.

In the asari/mummichog packages, the data entities are presented in any of the four types: 
class, namedtuple, JSON style dictionary or implicit list. 
The implicit lists are used sparely as they have reduced clarity. 
Namedtuple is immutable, then limited applications. 
In-memory searches are conducted using indexed dictionaries or dataframes.

SimpleSample is used by default asari workflow. 
For different purposes, one can use Sample or other classes.

Data formats:
===============
mass tracks as [( mz, rtlist, intensities ), ...].
Peak format: 
{
    'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 'height', 'parent_masstrace_id', 
    'rtime', 'peak_area', 'goodness_fitting'
}
isotopic_patterns = [(1.003355, 'M(13C)', 0, 0.2), ...]

Mass Tracks
===========
They are used for full RT ranges, thus each mass track has a unique m/z. 
Some chromatogram builders separate the mass traces if there are gaps in RT scans, 
but that creates complexity in m/z alignment and searches. 

Peak detection
==============
The main step uses scipy.signal.find_peaks, a local maxima method with prominence control.
Prominence is important, but it should be more tailored to individual peaks. 
Here, 
prominence = max(min_prominence_threshold, 0.05*max(list_intensity)).


'''
import os
import random
import numpy as np
from scipy.signal import find_peaks 

from pyopenms import MSExperiment, MzMLFile

# from mass2chem.annotate import annotate_formula_mass, massDict_hmdb

from .chromatograms import extract_massTracks_        # extract_massTraces, 
from .functions import *
from .sql import *

from .constructors import epdsConstructor


class SimpleSample:
    '''
    LC-MS Sample to get mass tracks from a mzML file, and as a container for peaks and empCpds.

    Use scan numbers whereas possible. initial peak detection for assembling empCpds only.

    Use dictionary formats for mass_trace, peak and empCpd for clarity.

    Registry for mass traces, peaks & empCpds.

    RT and m/z values will have calibration functions after CompositeMap.

    '''
    def __init__(self, experiment=None, mode='pos', input_file=''):

        self.input_file = input_file
        if experiment:
            self.experiment = experiment
            self.id = self.experiment.list_input_files.index(input_file)
        else:
            self.id = input_file
        
        self.mode = mode
        self.parameters = {}
        self.rt_numbers = []                                # list of scans, starting from 0
        self.list_retention_time = []                       # full RT time points in sample

        self.mz_calibration_ratio = None
        self.rt_calibration_function = None
        
        # lists to store data
        self.list_mass_tracks = []                          # index number = id_number, in ascending order of m/z
        self.anchor_mz_pairs = []
        self._number_anchor_mz_pairs_ = 0
        
        # will populate after CMAP
        self.list_peaks = []   
        self.list_empCpds = []
        
        # dict for future use
        self.list_mass_traces = []
        self.dict_mass_traces = {}
        self.dict_peaks = {}


    def process(self):
        '''
        From input file to list_MassTraces with detected peaks and selectivity on peaks.
        '''
        self.get_mass_tracks_()
        self.get_anchor_mz_pairs()
        print("_number_anchor_mz_pairs_ = %d \n" %self._number_anchor_mz_pairs_)

        self._mz_landmarks_ = flatten_tuplelist(self.anchor_mz_pairs)
        # to send to SQL DB here


    def get_mass_tracks_(self, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5):
        '''
        A mass track is an EIC for full RT range, without separating the mass traces,
        using asari.chromatograms algorithm.
        tracks as [( mz, rtlist, intensities ), ...].

        '''
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, min_intensity=min_intensity, min_timepoints=min_timepoints)
        self.rt_numbers = xdict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xdict['rt_times']     # full RT time points in sample
        ii = 0
        # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
        for tk in xdict['tracks']:                         
            self.list_mass_tracks.append( {
                'id_number': ii, 
                'mz': tk[0],
                'rt_scan_numbers': tk[1],                  # list
                'intensity': tk[2],                        # list
                } )
            ii += 1

        print("Processing %s, found %d mass tracks." %(self.input_file, ii))
        # For diagnosis only - check m/z split
        # warnings = check_close_mzs([x['mz'] for x in self.list_mass_tracks], mz_tolerance_ppm)
        # print("Warning - some mass tracks are too close to each other: ", len(warnings), warnings[:5])

    def get_anchor_mz_pairs(self):
        '''
        update self.get_anchor_mz_pairs
        e.g. [(5, 8), (6, 13), (17, 25), (20, 27), ...]
        '''
        self.anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(self.list_mass_tracks, mz_tolerance_ppm=5)
        self._number_anchor_mz_pairs_ = len(self.anchor_mz_pairs)
        
        



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

