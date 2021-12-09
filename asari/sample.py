'''


'''

import os


from .peaks import ext_MassTrace


from scipy.stats import norm as normal_distribution

from pyopenms import MSExperiment, MzMLFile


# from metDataModel.core import MassTrace, Peak, Experiment       # Feature, Sample later

from .chromatograms import extract_massTraces
from .search import search_formula_mass_dataframe

from .utils import *
from .sql import *



class Sample:
    '''
    A sample or injection, corresponding to one raw data file, either from positive or negative ionization.
    Each sample has a series of scans in LC-MS (retentiion time). 
    The RT is recorded as integers of scan number in asari, and converted to seconds on demand.
    The OpenMS workflow uses seconds directly, and users should class Sample_openms.

    '''
    def __init__(self, experiment, mode, input_file=''):
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
        self.__rt_calibration__ = None                      # This will be the calibration function

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
        Use another Sample class for a different algorithm.

        '''
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xic_dict = self._get_masstraces_from_centroided_rawdata_(exp)

        self.rt_numbers = xic_dict['rt_numbers']    # list of scans, starting from 0
        self.list_retention_time = xic_dict['rt_times']
        
        for xic in xic_dict['xics']:             # ( mz, rtlist, intensities )
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

        print("\nProcessing %s, found %d mass traces." %(self.input_file, len(self.dict_masstraces)))

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


