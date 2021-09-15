'''
asari, a simple program for LC-MS metabolomics data preprocessing.
Last update by Shuzhao Li, 2021-09-14
'''
import os
import random
from collections import namedtuple

import numpy as np
from scipy.signal import find_peaks  
from scipy.optimize import curve_fit 
from scipy.stats import norm as normal_distribution
from scipy.interpolate import UnivariateSpline

from pyopenms import MSExperiment, MzMLFile
from metDataModel.core import MassTrace, Peak, Experiment       # Feature, Sample later
from mass2chem.annotate import annotate_formula_mass 

from .sql import *

# starting point of ref DB, 
#INIT_DFDB = DB_to_DF( extend_DB1(DB_1) )
INIT_DFDB = tsv2refDB('hot_db.tsv')

# feature id will be assigned at the end; intensities is list; mass_id links to MassTrace
Feature = namedtuple('Feature', 'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'.split(',') )

def __gaussian_function__(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) 

def __goodness_fitting__(y_orignal, y_fitted):                  # R^2 as goodness of fitting
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


def search_formula_mass_dataframe(query_mz, DFDB, limit_ppm=10):
    '''
    return best match formula_mass in DFDB if under ppm limit.
    DFDB is using a Pandas DataFrame to house reference database.
    ppm is signed to capture the direction of mass shift.

    Not using idxmin,   #ii = DFDB.tmp.idxmin()
    # in pd.DF, index not necessarily integer; can be sequence if more than one match, but they are trying to fix in pandas dev version
    '''
    DFDB['tmp'] = abs(DFDB.mz - query_mz)
    ii = DFDB.tmp.values.argmin()
    ppm = 1000000 * (query_mz - DFDB.iloc[ii].mz)/query_mz
    if  abs(ppm) < limit_ppm:
        return (DFDB.iloc[ii].name, ppm)            # name is formula_mass
    else:
        return None


def bin_by_median_rt(List_of_peaks, tolerance):
    '''
    List_of_tuples: [(value, object), (value, object), ...], to be separated into bins by values.
                    objects have attribute of sample_name if to align elsewhere.
    return: [seprated bins], each as a list in the same input format. If all falls in same bin, i.e. [List_of_tuples].
    '''
    List_of_tuples = [(P.cal_rtime, P) for P in List_of_peaks]
    List_of_tuples.sort()
    new = [[List_of_tuples[0], ], ]
    for X in List_of_tuples[1:]:
        if X[0]-np.median([ii[0] for ii in new[-1]]) < tolerance:       # median moving with list change
            new[-1].append(X)
        else:
            new.append([X])
    PL = []
    for L in new:
        PL.append([X[1] for X in L])
    return PL

#
# revisit weak peaks???
#


def peaks_to_features(peak_dict, rtime_tolerance, ordered_sample_names):
    '''
    This uese high-selectivity peaks to establish features;
    still possibly not distinguishing close peaks well, which may be treated as combined peaks. 
    
    return
    ------
    List of Features (namedTuples, 'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'), following the input sample order.

    '''
    def __get_peaks_intensities__(peaks, ordered_sample_names):
        dict_intensities = {}
        for P in peaks: 
            if P.sample_name in dict_intensities:
                dict_intensities[P.sample_name] += P.peak_area      # redundant peaks in same m/z & rt
            else:
                dict_intensities[P.sample_name] = P.peak_area
        return [dict_intensities.get(name, 0) for name in ordered_sample_names]

    FeatureList = []
    for k,v in peak_dict.items():
        for F in bin_by_median_rt(v, rtime_tolerance):
            median_mz, median_rt = np.median([P.mz for P in F]), np.median([P.cal_rtime for P in F])
            peak_quality = max([P.goodness_fitting for P in F])                                 # will include np.median
            intensities = __get_peaks_intensities__(F, ordered_sample_names)
            #'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'
            FeatureList += [Feature(k, median_mz, median_rt, peak_quality, 1, intensities)]

    return FeatureList


# General data processing steps are in this class
class ext_Experiment(Experiment):
    '''
    Extend metDataModel.core.Experiment with preprocessing methods.
    This encapsulates a set of LC-MS files using the same method to be processed together.
    '''
    def __init2__(self, list_input_files, dict_meta_data, parameters, output_dir):
        '''
        This is the overall container for all data in an experiment/project.
        We don't sort sample orders. One can sort after getting the feature table.
        Input
        -----
        list_input_files: list of inputfiles, including directory path, to read
        dict_meta_data: description of sample types for each file, e.g. 'QC', 'pooled', 'sample'.
        parameters: including 'ionization_mode', 'min_intensity_threshold', 'min_timepoints'. See main.py.

        '''
        self.list_input_files = list_input_files
        self.output_dir = output_dir
        self.samples = []                   # out of initial input order in samples
        self.ordered_sample_names = []      # will populate after sample processing, enforced in assembling featuretable
        self.HOT_DB = {}                    # will be pd.DataFrame, to be initiated by chosen 3 smaples
        self.feature_table = {}             # will be pd.DataFrame
        
        self.number_of_samples = len(list_input_files)
        self.files_meta_data = dict_meta_data
        self.parameters = parameters
        self.max_rtime = parameters['max_rtime']
        self.mode = parameters['mode']

        self.initiation_samples = self.__choose_initiation_samples__()
        # self.interpolate_rtime_list = np.linspace(0, parameters['max_rtime'], 3000)       # ?? as only XICs are input
        

    def process_all(self):
        '''
        This will shift to a DB design in next version, not to keep all samples in memory.
        '''
        self.init_hot_db( INIT_DFDB )   # initial processing of 3 samples to set up HOT_DB
        # run remaining samples
        for f in self.list_input_files:
            if f not in self.initiation_samples:
                SM = Sample(self, self.mode, f)
                SM.process_step_1()
                SM.process_step_2(self.HOT_DB)
                self.samples.append(SM)
        
        self.calibrate_retention_time()
        self.correspondence()

        #self.retrieve_weak_peaks()
        #self.export_feature_table()



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

        #export_hot_db
        self.HOT_DB.to_csv("dynamic_hot_db.tsv", sep="\t")
        

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
            if len(rt_cal) < 30:
                print("\n\n-- Warning, RT regression using too few features (%s, %d) --\n\n" %(SM.name, len(xx)))
            rt_cal.sort()
            # to-do: need down sample, each bin no more than 10 data points

            xx, yy = [], []
            for L in rt_cal:
                if abs(L[0]-L[1])/L[1] < 0.2:       # control shift < 20%
                    xx.append(L[0])
                    yy.append(L[1])

            spl = UnivariateSpline(xx, yy, s=smoothing_factor)
            SM.__rt_calibration__ = spl
            SM.__rt_calibration__data__ = (xx, yy)

            # calibrate all detected RT for all peaks, and raw RT from MassTraces
            for P in SM.good_peaks:
                P.cal_rtime = SM.__rt_calibration__(P.rtime)
            for ML in SM.dict_masstraces.values():
                for M in ML:
                    M.cal_list_retention_time = SM.__rt_calibration__(M.list_retention_time)


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
        # drop rows by min presence in > 70% of samples, as QC or blank samples may behave differently
        rt_table = rt_table.dropna(axis=0, thresh=int(0.7 * self.number_of_samples))
        
        rt_table.to_csv("raw_rt_calibration_matrix.tsv", sep="\t")    # to export 
        return rt_table
        

    def correspondence(self):
        '''
        
        In each sample: Peak.mzstr links to Sample.mzstr_2_formula_mass, Sample.dict_masstraces

        Start feature table; not including Peak instances yet.
        detailed peak info will be pushed in to SQLite DB.
        '''
        self.ordered_sample_names = [SM.name for SM in self.samples]    # this is used to order intensity values in Features and help group/split Features

        # -- next == >>
        #  First get unique formula_mass or tentative mass_id for all good peaks.

        peak_dict = {}
        for SM in self.samples:
            for P in SM.good_peaks:
                try:
                    k = SM.mzstr_2_formula_mass[P.mzstr]
                    if k in peak_dict:
                        peak_dict[k].append(P)
                    else:
                        peak_dict[k] = [P]
                except KeyError:
                    pass

        # organize features by peak_dict; if >1 peaks in same sample fall into same feature, merge.
        FeatureList = peaks_to_features(peak_dict, self.parameters['rtime_tolerance'], self.ordered_sample_names)

        #with open('tuple_feature_table.tsv', 'w') as O:
        #    O.write('\n'.join([str(x) for x in FeatureList]))

        self.export_feature_table(FeatureList)

        self.FeatureTable = FeatureList



    def retrieve_weak_peaks(self):
        pass


    def export_feature_table(self, FeatureList, outfile='feature_table.tsv'):
        '''
        ftable = pd.DataFrame(ftable)
        ftable.to_csv(outfile, sep='\t')
        'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'
        '''
        outfile = os.path.join(self.output_dir, outfile)
        s = '\t'.join(['formula_mass', 'm/z', 'retention time', 'peak_quality', 'selectivity_rt',] + self.ordered_sample_names) + '\n'
        for F in FeatureList:
            s += '\t'.join(
                [F.mass_id, str(round(F.mz,4)), str(round(F.rtime,2)), str(round(F.peak_quality,2)), 
                str(round(F.selectivity_rt,2))] + [str(int(x)) for x in F.intensities]
            ) + '\n'
        with open(outfile, 'w') as O:
            O.write(s)


    def __process__unassigned_masstraces__(self):

        '''
        build index for unassigned
        determine how many mass_ids
        return [(mass_id, List_of_MassTraces), ...]

        '''


        pass

    def __choose_initiation_samples__(self):
        if self.parameters['initiation_samples']:
            return self.parameters['initiation_samples'][:3]
        else:
            if self.number_of_samples < 4:
                return self.list_input_files
            elif not self.files_meta_data:
                return random.sample(self.list_input_files, 3)
            else:
                chosen = []
                POOLED = [f for f in self.list_input_files if self.files_meta_data[f] == 'POOLED']
                if POOLED:
                    chosen.append(random.choice(POOLED))
                QC = [f for f in self.list_input_files if self.files_meta_data[f] == 'QC']
                if QC:
                    chosen.append( random.choice(QC) )
                OTHERS = [f for f in self.list_input_files if self.files_meta_data[f] not in ['POOLED', 'QC', 'BLANK']]
                chosen += random.sample(OTHERS, 3)
                return chosen[:3]


class Sample:
    '''
    A sample or injection, corresponding to one raw data file, either from positive or negative ionization.
    '''
    def __init__(self, experiment, mode, input_file=''):
        self.input_file = input_file
        self.experiment = experiment    # parent Experiment instance
        self.name = os.path.basename(input_file)            # must be unique
        self.mode = self.experiment.mode
        self.parameters = self.experiment.parameters

        self.dict_masstraces = {}       # indexed by str(round(mz,6))
        self.mzstr_2_formula_mass = {}
        self.good_peaks = []

        self.number_MassTraces = 0
        self.__mass_accuracy__, self.__mass_stdev__ = 0, 0
        self.__rt_calibration__ = None          # This will be the calibration function

        # internal use only
        self.__mass_index__ = {}
        self.__unassigned_masstraces__ = []
        self.list_MassTraces = []       # fixed sequence

    def process_step_1(self):
        '''
        From input file to list_MassTraces with detected peaks and selectivity on peaks.
        '''
        self._get_masstraces_from_chromatogram_file_()
        self._detect_peaks_()
        self._assign_selectivity_()

    def process_step_2(self, DFDB):
        '''
        Annotate MassTraces with unique formula_mass.
        This is based on reference DB and may not cover all MassTraces.
        The remaining steps (correspondence, alignment) will be performed at Experiment level.
        DFDB: either from positive or negative ionization

        There may be repeated formula_mass - traces here, leaving to step_3 and masstraces_to_features to handle.
        '''
        self._match_mass_formula_(DFDB)
        for P in self.good_peaks:
            P.sample_name = self.name

    def export_peaklist(self):
        outfile = os.path.join(self.experiment.output_dir, self.name.replace('.mzML', '') + '.peaklist')
        header = ['m/z', 'retention time', 'area', 'shape_quality', 'gaussian_amplitude', 'gaussian_variance', 'mz_selectivity']
        peaklist = []
        for P in self.good_peaks:
            formatted = [str(round(P.mz, 6)), str(round(P.rtime, 2)), str(int(P.peak_area)), 
                            str(round(P.goodness_fitting, 2)), str(int(P.gaussian_parameters[0])), 
                                str(round(P.gaussian_parameters[2], 2)), str(round(P.selectivity, 2)),]
            if P.mzstr in self.mzstr_2_formula_mass:
                formatted.append(self.mzstr_2_formula_mass[P.mzstr])
            peaklist.append(formatted)

        with open(outfile, 'w') as O:
            O.write( '\t'.join(header) + '\n' + '\n'.join([ '\t'.join(L) for L in peaklist ]) + '\n' )

    def _detect_peaks_(self):
        for mzstr, ML in self.dict_masstraces.items():
            for M in ML:
                list_peaks = M.detect_peaks(self.parameters['min_intensity_threshold'], self.parameters['min_timepoints'], 
                                self.parameters['min_prominence_threshold'], self.parameters['prominence_window'], self.parameters['gaussian_shape'])
                if list_peaks:
                    for P in list_peaks:
                        P.mzstr = mzstr
                        self.good_peaks.append(P)

        print("Detected %d good peaks." %len(self.good_peaks))

    def _get_masstraces_from_chromatogram_file_(self):
        '''
        Get chromatograms, via pyOpenMS functions (only used in this function). 
        An mzML file of chromatograms (i.e. EIC or XIC, or mass trace). Currently using a MassTrace file produced by FeatureFinderMetabo (OpenMS). 
        min_intensity_threshold: minimal intensity requirement for a mass trace to be considered. Default value 10,000 is lenient in our Orbitrap data.
        Updates self.dict_masstraces, a dict of list of MassTrace objects, indexed by str(round(mz,6))
        '''
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        for chromatogram in exp.getChromatograms():
            mz = chromatogram.getPrecursor().getMZ()
            _low, _high = self.parameters['mass_range']
            if _low < mz < _high:
                mz_str = str(round(mz,6))
                RT, INTS = chromatogram.get_peaks() 
                if INTS.max() > self.experiment.parameters['min_intensity_threshold'] * 0.1:    # use threshold for good peaks; keep more
                    M = ext_MassTrace()
                    M.__init2__(mz, RT, INTS)
                    if mz_str in self.dict_masstraces:
                        self.dict_masstraces[mz_str].append(M)
                    else:
                        self.dict_masstraces[mz_str] = [M]

        print("\nProcessing %s, found %d mass traces." %(self.input_file, len(self.dict_masstraces)))

    def _assign_selectivity_(self, std_ppm=5):
        Q = [(P.mz, P.rtime, P) for P in self.good_peaks]
        Q.sort()
        selectivities = calculate_selectivity([x[0] for x in Q], std_ppm)
        for ii in range(len(Q)):
            Q[ii][2].selectivity = selectivities[ii]

    def _match_mass_formula_(self, DFDB, check_mass_accuracy=True, ppm=15):
        '''
        Match peaks to formula_mass database, including mass accuracy check and correction if needed.
        Only high-selectivity (> 0.9) peaks will be used downstream for alignment calibration.
        Mass accuracy has no hard limit at the moment, but should be restricted in future versions.
        E.g. if mass shift is larger than 20 ppm, user should do own mass calibrate first.

        DFDB: a reference DB in DataFrame.
        check_mass_accuracy should be on for all samples in asari processing.
        ppm: use a larger initial number, and the check step will reduce based on stdev of ppm.
        '''
        list_ppm_errors, mzstr_2_formula_mass = [], {}
        for P in self.good_peaks:
            query = search_formula_mass_dataframe(P.mz, DFDB, ppm)
            if query:
                _formula_mass, delta_ppm = query
                mzstr_2_formula_mass[P.mzstr] = _formula_mass
                list_ppm_errors.append(delta_ppm)

        self.__mass_accuracy__, self.__mass_stdev__ = mass_accuracy, ppm_std = normal_distribution.fit(list_ppm_errors)
        print("    (mass_accuracy, mass_stdev) = (%5.2f, %5.2f ) " %(mass_accuracy, ppm_std))
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

        for mzstr in self.dict_masstraces.keys():
            if mzstr not in mzstr_2_formula_mass:
                query = search_formula_mass_dataframe(self.dict_masstraces[mzstr][0].mz, DFDB, 2*self.__mass_stdev__)
                if query:
                    mzstr_2_formula_mass[mzstr], _ = query

        self.mzstr_2_formula_mass = mzstr_2_formula_mass


# peak detection is in this class
class ext_MassTrace(MassTrace):
    '''
    Extending metDataModel.core.MassTrace
    Peak detection using scipy.signal.find_peaks, a local maxima method with prominence control.
    Not keeping Peaks in this class; Peaks are attached to Sample.
    To-do: further refine peak detection in future, 
    e.g. extra iteration of peak detection in remaining region; better deconvolution of overlapping peaks.

    '''
    def __init2__(self, mz, RT, INTS):
        self.mz, self.list_retention_time, self.list_intensity = [mz, RT, INTS]
        self.cal_list_retention_time = []
        self.formula_mass = None
        self.raw_mzs = []                               # multiple values possible when merging traces
        self.__mass_corrected_by_asari__ = False        # True if updated by calibration
        self.features_assigned = False                  # tracking if assigned to Experiment features
        self.sample_name = ''

    def detect_peaks(self, min_intensity_threshold, min_timepoints, min_prominence_threshold, prominence_window, gaussian_shape):
        list_peaks = []
        peaks, properties = find_peaks(self.list_intensity, 
                                    height=min_intensity_threshold, width=min_timepoints, 
                                    prominence=min_prominence_threshold, wlen=prominence_window) 
        
        if peaks.size > 1:
            # Rerun, raising min_prominence_threshold. Not relying on peak shape for small peaks, 
            # because chromatography is often not good so that small peaks can't be separated from noise.
            _tenth_height = 0.1*self.list_intensity.max()
            min_prominence_threshold = max(min_prominence_threshold, _tenth_height)
            peaks, properties = find_peaks(self.list_intensity, 
                                    height=min_intensity_threshold, width=min_timepoints, 
                                    prominence=min_prominence_threshold, wlen=prominence_window)
        
        for ii in range(peaks.size):
            P = ext_Peak()
            # scipy.signal.find_peaks works on 1-D array. RT coordinates need to be added back.
            P.peak_initiate(parent_mass_trace=self, mz = self.mz, apex = peaks[ii],
                            peak_height = properties['peak_heights'][ii],
                            left_base = properties['left_bases'][ii],
                            right_base = properties['right_bases'][ii])
            P.evaluate_peak_model()
            if P.goodness_fitting > gaussian_shape:
                list_peaks.append(P)

        return list_peaks


# peak evaluation is in this class
class ext_Peak(Peak):
    '''
    Extending metDataModel.core.Peak.
    Include pointer to parent MassTrace.
    Not storing RT or intensity arrays, but looking up in MassTrace.
    '''
    def peak_initiate(self, parent_mass_trace, mz, apex, peak_height, left_base, right_base):
        [ self.parent_mass_trace, self.mz, self.apex, self.peak_height, self.left_base, self.right_base
                ] = [ parent_mass_trace, mz, apex, peak_height, left_base, right_base ]

    def evaluate_peak_model(self):
        '''
        Use Gaussian models to fit peaks, R^2 as goodness of fitting.
        Peak area is defined by Gaussian model,
        the integral of a Gaussian function being a * c *sqrt(2*pi).

        Left to right base may not be full represenation of the peak. 
        The Gaussian function will propose a good boundary.
        '''
        xx = self.parent_mass_trace.list_retention_time[self.left_base: self.right_base+1]
        # set initial parameters
        a, mu, sigma =  self.peak_height, \
                        self.parent_mass_trace.list_retention_time[self.apex], \
                        xx.std()
        try:
            popt, pcov = curve_fit(__gaussian_function__, 
                            xx,
                            self.parent_mass_trace.list_intensity[self.left_base: self.right_base+1],
                            p0=[a, mu, sigma])

            self.gaussian_parameters = popt
            self.rtime = popt[1]
            self.peak_area = popt[0]*popt[2]*2.506628274631
            self.goodness_fitting = __goodness_fitting__(
                            self.parent_mass_trace.list_intensity[self.left_base: self.right_base+1], 
                            __gaussian_function__(xx, *popt))
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
        self.y_fitted_extended = __gaussian_function__(self.rt_extended, *self.gaussian_parameters)
