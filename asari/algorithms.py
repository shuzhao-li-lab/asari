'''
Reusing class definitions in metDataModel.

The flow in asari:
    Get MassTraces
        Clean up redundant traces, calculate selectivity (mz)

    Peak detection per MassTrace, using simple local maxima (required height and prominence)
        If multiple peaks are found, redo peak dection by raising prominence to (max_intensity/10).
        Further refining can be added, e.g. redo peak detection in the remaining region of trace.

    Peak evaluation, calculate peak quality, RT, area.

    Correspondence, initial
        on peaks of high selectivity, by assigning to formula_mass

    Mass calibration check, calibrate and recalculate selectivity if needed
    
    Correspondence, RANSAC on peaks of low selectivity

    Output feature table

    (Further annotation is done via mass2chem)
    # future: Retention indexing

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

# starting point of ref DB, in DF format. Will update
# this will be THE ref DB.
# INIT_DFDB = DB_to_DF( extend_DB1(DB_1) )
INIT_DFDB = tsv2refDB('hot_db.tsv')

# feature id will be assigned at the end; intensities is list; mass_id links to MassTrace
Feature = namedtuple('Feature', 
                    'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'.split(',')    
                    )

def __gaussian_function__(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) 

def __goodness_fitting__(y_orignal, y_fitted):
    # R^2 as goodness of fitting
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


def binning_by_median(List_of_tuples, tolerance):
    '''
    List_of_tuples: [(value, object), (value, object), ...], to be separated into bins by values.
                    objects have attribute of sample_name if to align elsewhere.
    return: [seprated bins], each as a list in the same input format. If all falls in same bin, i.e. [List_of_tuples].
    '''
    List_of_tuples.sort()
    new = [[List_of_tuples[0], ], ]
    for X in List_of_tuples[1:]:
        if X[0]-np.median([ii[0] for ii in new[-1]]) < tolerance:       # median moving with list change
            new[-1].append(X)
        else:
            new.append([X])
    return new


def composite_mass_trace(List_of_MassTraces):
    
    pass



#
# revisit weak peaks???
#


def masstrace_to_features(List_of_MassTraces, common_MassTraces_id, parameters, ordered_sample_names):
    '''
    To define features from MassTraces, which have detected Peaks,
    and report simple feature attributes towards an Experiment-wide feature table.
    Full queries should be done via the real DB.
    1. If each MassTrace contains no more than one peak and with good RT agreement,
    report back as one feature.
    2. else multiple peaks per MassTrace, redo peak detection on pooled elution profiles.

    input
    -----
    List_of_MassTraces, nested list of MassTrace instances of the same m/z or formula, [SM.__mass_index__.get(fid, []) for SM in self.samples],
                        one list per sample; each sample may return 0 to N MassTraces.
                        
    common_MassTraces_id: formula_mass or internal id.
    ordered_sample_names: from Experiment

    return
    ------
    List of Features (namedTuples, 'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'), following the input sample order.

    '''



    def __get_peaks_intensities__(peaks, ordered_sample_names):
        dict_intensities = {}
        for P in peaks: 
            dict_intensities[P.sample_name] = P.peak_area
        return [dict_intensities.get(name, 0) for name in ordered_sample_names]

    def __make_a_feature__():
        pass


    # flaten_to_Peaks, nested list, one list of peaks per Sample, [[Sample1.Peak1, peak2, ...], [Sample2 empty], [Sample3Peak1, ....], ...]
    peakLevelList = []
    for T in List_of_MassTraces:
        peakLevelList.append([peak for M in T for peak in M.list_peaks])

    FeatureList = []
    if max([len(T) for T in peakLevelList]) == 1:
        tupleList = [(T[0].cal_rtime, T[0]) for T in peakLevelList if T]
        for LL in binning_by_median(tupleList, parameters['rtime_tolerance']):
            # separated tupleList if RT not in tolerance range
            peaks = [X[1] for X in LL]
            median_mz, median_rt = np.median([P.mz for P in peaks]), np.median([X[0] for X in LL])
            peak_quality = np.median([P.goodness_fitting for P in peaks])
            intensities = __get_peaks_intensities__(peaks, ordered_sample_names)
            #'mass_id,mz,rtime,peak_quality,selectivity_rt,intensities'
            FeatureList += [Feature(common_MassTraces_id, median_mz, median_rt, peak_quality, 1, intensities)]

    else:
        # combine MassTraces per sample into orphan MassTrace

        # M = __composite_mass_trace__(List_of_MassTraces)

        M = [x for x in List_of_MassTraces if x][0][0]     # test

        M.detect_peaks(parameters['min_intensity_threshold'], 
                        parameters['min_timepoints'], 
                        parameters['min_prominence_threshold'], 
                        parameters['prominence_window'],
                        parameters['gaussian_shape'],
                        )
        PEAKS = M.list_peaks
        FeatureList = []

    return FeatureList





class ext_Experiment(Experiment):
    '''
    Extend metDataModel.core.Experiment with preprocessing methods.
    This encapsulates a set of LC-MS files using the same method to be processed together.

    Steps:
    1) init_hot_db: use three samples to establish expt wide parameters like RT vector, and 
    a HOT_DB to hold all expt specific annotations.
    All other samples will be annotated to formula_mass using the HOT_DB.
    Leave annotation of remaining featuers to the end (No other ext search at sample level than initial effort).
    2) Initial correspondence of `good` peaks, and RT calibration
    3) Correspondency of remaining peaks. Leave empCpd organization and correction at the end.

    Not used now:
    # create new DB for common isotopes and adducts based on what's found earlier
    extended_DB = create_extended_DB([x[4][0][0] for x in result if x], mode)
    # extend search to common isotopes and adducts
    if ppm_std:
        limit_ppm = 2 * ppm_std
    for ii in range(N):
        if not result[ii]:
            result[ii] = search_formula_mass_db(list_query_mz[ii], extended_DB, limit_ppm)
    
    '''

    def __init2__(self, list_input_files, dict_meta_data, parameters, output_dir):
        '''
        This is the overall container for all data in an experiment/project.
        We don't sort sample orders. One can sort after getting the feature table.

        # will cumulate mass accuracy and stdev for all samples and report

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
        # initial processing of 3 samples to set up HOT_DB
        self.init_hot_db( INIT_DFDB )

        # run remaining samples
        for f in self.list_input_files:
            if f not in self.initiation_samples:
                SM = Sample(self, self.mode, f)
                SM.process_step_1()
                SM.process_step_2(self.HOT_DB)
                self.samples.append(SM)
        
        self.calibrate_retention_time()
        self.correspondence()
        #self.export_feature_table()




    def init_hot_db(self, DFDB):
        '''
        DFDB formatted as:
        {...}
        Watch out for redundancy there...???


        [['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]], 

         adducts: e.g. [(58.53894096677, 'M+2H[2+]', result_formula), ...,]

        Because this custom_database is updated since the 1st sample, repeated features will be fast to access in later samples.


        Use three samples to initiate a hot DB to house feature annotation specific to this Experiment.
        The HOT_DB will be used during sample processing, and have another update after correspondence and additional annotation.
        The HOT_DB will then be exported as Expt annotation.
        
        '''
        chosen_Samples, found_formula_masses = [], []
        for f in self.initiation_samples:
            SM = Sample(self, self.mode, f)
            SM.process_step_1()
            SM.process_step_2(DFDB)
            chosen_Samples.append(SM)
            for M in SM.list_MassTraces:
                found_formula_masses.append(M.formula_mass)

        self.samples += chosen_Samples
        # Experiment wide parameters
        self.__mass_stdev__ = np.median([SM.__mass_stdev__ for SM in chosen_Samples])       # ppm stdev, used for later searches
        # ver 1 HOT_DB will be a subset of INIT_DFDB

        found_formula_masses = set(found_formula_masses)
        if None in found_formula_masses:
            found_formula_masses.remove(None)
        self.HOT_DB = DFDB.loc[found_formula_masses]   

        print("\n[@.@] Anchoring with %d initial formula matches." %len(found_formula_masses))
        print("[@.@] Initial estimation done on ", self.initiation_samples)

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

        For the  mass traces with more than one peaks, after RT calibration,
        do a composite chromatogram by summing up all traces,
        then do peak detection to decide how many peaks are there.

        # to-do: issue warning for samples with too few features

        Not used:
        rt_table['bin_index'] = 10*(rt_table.median(axis=1)/10).apply(int)
        for B in set(rt_table['bin_index']):
            slice = rt_table[rt_table['bin_index']==B].values.tolist()       

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
            # to-do:
            # need down sample, each bin no more than 10 data points
            #
            xx, yy = [], []
            for L in rt_cal:
                if abs(L[0]-L[1])/L[1] < 0.2:       # control shift < 20%
                    xx.append(L[0])
                    yy.append(L[1])

            spl = UnivariateSpline(xx, yy, s=smoothing_factor)
            SM.__rt_calibration__ = spl
            SM.__rt_calibration__data__ = (xx, yy)

            # calibrate all detected RT for all peaks, not raw RT from MassTraces, 
            # which will be done only in masstrace_to_features if redoing peak detection for complex XICs
            for M in SM.list_MassTraces:
                M.sample_name = SM.name
                for P in M.list_peaks:
                    P.sample_name = SM.name
                    P.cal_rtime = SM.__rt_calibration__(P.rtime)


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
            for M in SM.list_MassTraces:
                if M.formula_mass and M.number_of_peaks == 1 and M.selectivity > 0.98:
                    if M.list_peaks[0].goodness_fitting > 0.9:
                        good_peaks_rtime.append( M.list_peaks[0].rtime)
                        good_peaks_formula_mass.append(M.formula_mass)

            d[SM.name] = pd.Series(good_peaks_rtime, index=good_peaks_formula_mass)

        #print([v[:200] for v in d.values()])
        # merge into a table, each row as feature, col as sample
        rt_table = pd.DataFrame(d)
        
        # avoiding pd.DataFrame whenever possible, unpredictable behaivors
        # drop rows by min presence in > 70% of samples, as QC or blank samples may behave differently
        rt_table = rt_table.dropna(axis=0, thresh=int(0.7 * self.number_of_samples))

        # to export 
        # rt_table.to_csv("raw_rt_calibration_matrix.tsv", sep="\t")
        return rt_table
        

    def correspondence(self):
        '''
        Start feature table; not including Peak instances yet.
        Each sample has a unique formula_mass for MassTraces with single peaks.
        The feature table now only counts what formula_masses qualify in each sample.

        detailed peak info will be pushed in to SQLite DB.

        keep min for feature table here

        Step 1. get unique mass traces with single peaks and formula_mass; calibrate RT; join all samples 
        Step 2. for unique mass traces with multiple peaks and formula_mass, calibrate RT, redefine peaks by patterns across samples

        Step 3. ransac remaining mass traces and peaks, join across samples



        There is benefit of empCpd grouping per sample before correspondence, 
        to trade confidence for performance.




            good_peaks_areas, good_peaks_formula_mass = [], []
            for M in SM.list_MassTraces:
                if M.formula_mass and M.number_of_peaks == 1 and M.selectivity > 0.9:
                    good_peaks_areas.append( M.list_peaks[0].peak_area)
                    good_peaks_formula_mass.append(M.formula_mass)
                    #
                    M.features_assigned = True

            d[SM.name] = pd.Series(good_peaks_areas, index=good_peaks_formula_mass)




            
        '''
        self.ordered_sample_names = [SM.name for SM in self.samples]    # this is used to order intensity values in Features and help group/split Features

        for SM in self.samples:
            SM.process_step_3()     # updating masstrace index

        featuretable0 = []
        for fid in self.HOT_DB.index:
            # Each SM may return 0 to N MassTraces per index
            List_of_MassTraces = [ SM.__mass_index__.get(fid, []) for SM in self.samples ]
            featuretable0 += masstrace_to_features( List_of_MassTraces, fid, self.parameters, self.ordered_sample_names )


        with open('tuple_feature_table.tsv', 'w') as O:
            O.write('\n'.join([str(x) for x in featuretable0]))

        # recalculate selectivity on both m/z and RT; get a composite quality score

        # self.export_feature_table(featuretable0)

        for fid, L in self.__process__unassigned_masstraces__():
            featuretable0 += masstrace_to_features( L, fid, self.parameters, self.ordered_sample_names )

        self.FeatureTable = featuretable0






    def export_feature_table(self, ftable, outfile='feature_table.tsv'):
        '''
        Input
        -----
        list_intensityxxxx

        Return
        ------
        list of Peaksxxxxxxxx

                    # add new varialbe to peak, cal_rtime
                    # M.list_peaks[0].cal_rtime = SM.__rt_calibration__(M.list_peaks[0].rtime)
        header = 'FID\t' + '\t'.join(self.ordered_sample_names) + '\n'
            O.write('\n'.join(
                ['\t'.join([str(x) for x in row]) for row in ftable]
            ))
        '''
        for s in self.samples:
            print(s.name)
        ftable = pd.DataFrame(ftable)
        ftable.to_csv(outfile, sep='\t')

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

    Example
    -------        
    SM = Sample(f)
    SM.process_step_1()
    SM.export_peaklist()    # to save intermediate files; one can stop here and resume later
    SM.process_step_2()

    Two workflows before SM.process_step_2():
        A. from a XIC file, do read_chromatogram_file/detect_peaks, then assign_selectivity;
        B. read from an intermediate peaklist file (generated by option 1).

    Peak detection is performed at mass trace level.
    Mass calibration and annotating formula_mass is at sample level.
    Retention time calibration and correspondency is done across samples at experiment level.
    '''
    def __init__(self, experiment, mode, input_file=''):
        self.input_file = input_file
        self.experiment = experiment    # parent Experiment instance
        self.name = os.path.basename(input_file)            # must be unique
        self.mode = self.experiment.mode
        # self.parameters = self.experiment.parameters
        
        self.list_MassTraces = []       # fixed sequence
        self.number_MassTraces = 0
        self.mz_list = []               # array
        self.peak_table = {}            # will be pd.DataFrame

        self.__process_stage__ = 0          # 0 = unprocessed, 1 = step_1, ..
        self.__mass_accuracy__, self.__mass_stdev__ = 0, 0
        self.__rt_calibration__ = None          # This will be the calibration function

        # internal use only
        self.__mass_index__ = {}
        self.__unassigned_masstraces__ = []

    def process_step_1(self, from_peaklist=False):
        '''
        From input file to list_MassTraces with detected peaks and selectivity.
        One can pause here by following with export_peaklist(), which saves a peaklist file.
        '''
        if from_peaklist:
            self._read_asari_peakfile_(self.input_file)
        else:
            self._detect_peaks_( self.experiment.parameters['min_intensity_threshold'], 
                                 self.experiment.parameters['min_timepoints'], 
                                 self.experiment.parameters['min_prominence_threshold'], 
                                 self.experiment.parameters['prominence_window'],
                                 self.experiment.parameters['gaussian_shape'],
                                )
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

    def process_step_3(self):
        '''
        Indexing MassTraces with formula_mass;
        If multiple MassTraces have same formula_mass, return the list.
        Do not merge them here, not to confuse tracing back to source data. Dealt with at masstrace_to_features.

        updating
        --------
        self.__mass_index__
        self.__unassigned_masstraces__ 
        '''
        for M in self.list_MassTraces:
            if M.formula_mass:
                if M.formula_mass in self.__mass_index__:
                    self.__mass_index__[M.formula_mass].append(M)
                else:
                    self.__mass_index__[M.formula_mass] = [M]
            else:
                self.__unassigned_masstraces__.append(M)


    def export_peaklist(self):
        '''
        Export tsv peak list for intermediate storage or diagnosis.
        Requiring selectivity assigned.
        '''
        outfile = os.path.join(self.experiment.output_dir, self.name.replace('.mzML', '') + '.peaklist')
        header = ['m/z', 'retention time', 'area', 'shape_quality', 'gaussian_amplitude', 'gaussian_variance', 'mz_selectivity']
        peaklist = []
        for M in self.list_MassTraces:
            for P in M.list_peaks:
                formatted = [str(round(M.mz, 6)), str(round(P.rtime, 2)), str(int(P.peak_area)), 
                                str(round(P.goodness_fitting, 2)), str(int(P.gaussian_parameters[0])), 
                                  str(round(P.gaussian_parameters[2], 2)), str(round(M.selectivity, 2)),
                                  ]
                if M.formula_mass:
                    formatted.append(M.formula_mass)
                peaklist.append(formatted)

        with open(outfile, 'w') as O:
            O.write(
                    '\t'.join(header) + '\n' + '\n'.join([ '\t'.join(L) for L in peaklist ]) + '\n'
            )

    #
    # not to call below functions from outside
    #
    def _detect_peaks_(self, min_intensity_threshold, min_timepoints, min_prominence_threshold, prominence_window, gaussian_shape):
        list_MassTraces = self._read_chromatogram_file_(self.input_file, min_intensity_threshold)
        for M in list_MassTraces:
            M.detect_peaks(min_intensity_threshold, min_timepoints, min_prominence_threshold, prominence_window, gaussian_shape) 
            if M.list_peaks:
                self.list_MassTraces.append(M)

        self.number_MassTraces = len(self.list_MassTraces)
        print("After peak detection, got %d valid mass traces." %self.number_MassTraces)

    def _read_chromatogram_file_(self, infile, min_intensity_threshold):
        '''
        Get chromatograms, via pyOpenMS functions (only used in this function). 
        Convert to list MassTraces for transparency.
        The rest of software uses metDataModel class extensions.

        A decision point is how to merge XICs of identical m/z.
        From OpenMS, redundant peaks may come from big peaks spilled over (shoulder), 
        and they should be merged.
        But if the peaks don't overlap, they should be kept and fillers are added to the trace.
        Because this software centers on formula mass, we will keep one XIC (MassTrace) per m/z.

        Input
        -----
        An mzML file of chromatograms (i.e. EIC or XIC, or mass trace).
        Currently using a MassTrace file produced by FeatureFinderMetabo (OpenMS). 
        min_intensity_threshold: minimal intensity requirement for a mass trace to be considered. 
        Default value 10,000 is lenient in our Orbitrap data.

        Return
        ------
        A list of MassTrace objects, sorted by increasing m/z values.

        '''
        def __concatenate_rt__(T1, T2, number_steps):
            # filling gap between T1 and T2 retention times, T2 must be greater than T1
            return np.concatenate((T1, np.linspace(T1[-1], T2[0], number_steps), T2))

        def __concatenate_ints__(T1, T2, number_steps):
            # filing zeros between T1 and T2 intensities
            return np.concatenate((T1, np.zeros(number_steps), T2))

        exp = MSExperiment()                                                                                          
        MzMLFile().load(infile, exp)
        list_MassTraces = []
        mzDict = {}
        for chromatogram in exp.getChromatograms():
            mz = chromatogram.getPrecursor().getMZ()
            mz_str = str(round(mz,6))
            RT, INTS = chromatogram.get_peaks() 
            retention_time_step = 0.5 * (RT[2] - RT[0])     # assuming XIC must have min 3 data points
                                                            # this could be more efficient if all RT uses scan number integers
            if INTS.max() > min_intensity_threshold:
                if mz_str not in mzDict:
                    mzDict[mz_str] = [mz, RT, INTS]

                elif RT[0] > mzDict[mz_str][1][-1]:
                    number_steps = int((RT[0] - mzDict[mz_str][1][-1])/retention_time_step)
                    mzDict[mz_str] = [mz, __concatenate_rt__(mzDict[mz_str][1], RT, number_steps), __concatenate_ints__(mzDict[mz_str][2], INTS, number_steps)]

                elif mzDict[mz_str][1][0] > RT[-1]:
                    number_steps = int(( mzDict[mz_str][1][0] - RT[-1] )/retention_time_step)
                    mzDict[mz_str] = [mz, __concatenate_rt__(RT, mzDict[mz_str][1], number_steps), __concatenate_ints__(INTS, mzDict[mz_str][2], number_steps)]

                elif INTS.sum() > mzDict[mz_str][2].sum():    
                    # overwrite if bigger trace found; watch out nuances in the future updates
                    mzDict[mz_str] = [mz, RT, INTS]

        for [mz, RT, INTS] in sorted(mzDict.values()):
            M = ext_MassTrace()
            M.__init2__(mz, RT, INTS)
            list_MassTraces.append( M )

        print("\nProcessing %s, found %d mass traces." %(self.input_file, len(list_MassTraces)))
        return list_MassTraces


    def _read_asari_peakfile_(self, infile):
        '''
        Read back peaklist generated by self.export_peaklist(), into MassTraces.
        The export function allows intermediate results to be stored, and easy parallel processing of peak detection.

        A MassTrace is defined by a unique m/z str, rounded to 6 digits in self.export_peaklist().
        Pre-sorted prior to export.
        Each Peak is from one row of ['m/z', 'retention time', 'area', 'shape_quality', 'gaussian_amplitude', 'gaussian_variance', 'mz_selectivity']
        '''
        PL = open(infile).read().splitlines()
        _mz_order = []
        mzdict = {}
        for L in PL:
            a = L.rstrip().split('\t')
            mz_str = a[0]
            _mz_order.append(mz_str)
            # this is not full ext_Peak, not supporting plot of XIC, 
            # which will require more exported information, and/or read the initial chromatogram file
            P = Peak()
            [P.mz, P.rtime, P.peak_area, P.goodness_fitting, P.gaussian_amplitude, P.gaussian_variance, P.selectivity
                    ] = [ float(x) for x in a ]
            if mz_str in mzdict:
                mzdict[mz_str].list_peaks.append(P)
            else:
                M = ext_MassTrace()
                M.mz = P.mz
                M.list_peaks = [P]
                mzdict[mz_str] = M

        self.list_MassTraces = [ mzdict[x] for x in _mz_order ] # preserving presorted m/z order
        self.number_MassTraces = len(self.list_MassTraces)


    def _match_mass_formula_(self, DFDB, check_mass_accuracy=True, ppm=10):
        '''
        Match peaks to formula_mass database, including mass accuracy check and correction if needed.
        Only high-selectivity (> 0.9) peaks will be used downstream for alignment calibration.
        Mass accuracy has no hard limit at the moment, but should be restricted in future versions.
        E.g. if mass shift is larger than 20 ppm, user should do own mass calibrate first.

        DFDB: a reference DB in DataFrame.
        check_mass_accuracy should be on for all samples in asari processing.
        ppm: use a larger initial number, and the check step will reduce based on stdev of ppm.

        This modifies
        -------------
        if matched to formula mass, per MassTrace:
            M.formula_mass
            M.ppm_from_formula_mass - only calculated for initial matches (high selectivity), 
                                      as self.__mass_stdev__ will be used in subsequent steps.

        self.__mass_accuracy__, self.__mass_stdev__
        if mass calibration/correction is needed:
            M.mz, M.raw_mz, 
            M.__mass_corrected_by_asari__ = True

        Diagnosis
        ---------
        outfile = os.path.join(self.experiment.output_dir, self.name.replace('.mzML', '') + '_errors.tsv')
        with open(outfile, 'w') as O:
            O.write(
                    '\n'.join([ str(x) for x in list_ppm_errors ]) + '\n'
            )

        '''
        list_ppm_errors = []
        for M in self.list_MassTraces:
            query = search_formula_mass_dataframe(M.mz, DFDB, ppm)
            if query:
                M.formula_mass, delta_ppm = query
                list_ppm_errors.append(delta_ppm)

        self.__mass_accuracy__, self.__mass_stdev__ = mass_accuracy, ppm_std = normal_distribution.fit(list_ppm_errors)
        print("    (mass_accuracy ppm mass_stdev) = ", self.__mass_accuracy__, self.__mass_stdev__)

        if check_mass_accuracy:
            if abs(mass_accuracy) > 5:   # this is considered significant mass shift, requiring m/z correction for all
                list_ppm_errors = []
                for M in self.list_MassTraces:
                    M.raw_mzs = [M.mz]
                    M.mz = M.mz - M.mz*0.000001*mass_accuracy
                    M.__mass_corrected_by_asari__ = True
                    # redo search because we may find different matches after mass correction; update search ppm too
                    query = search_formula_mass_dataframe(M.mz, DFDB, 2*ppm_std)
                    if query:
                        M.formula_mass, delta_ppm = query
                        list_ppm_errors.append( delta_ppm )
                # update ppm_std because it may be used for downstream search parameters
                _mu, self.__mass_stdev__ = normal_distribution.fit(list_ppm_errors)

    def _assign_selectivity_(self, std_ppm=5):
        selectivities = calculate_selectivity([M.mz for M in self.list_MassTraces], std_ppm)
        for ii in range(self.number_MassTraces):
            self.list_MassTraces[ii].selectivity = selectivities[ii]



#
# peak detection is in this class
#
class ext_MassTrace(MassTrace):
    '''
    Extending metDataModel.core.MassTrace

    Peak detection using scipy.signal.find_peaks, a local maxima method with prominence control.
    Example
    -------

    peaks:    array([10, 16]),
    properties:    {'peak_heights': array([3736445.25, 5558352.5 ]),
                    'prominences': array([1016473.75, 3619788.  ]),
                    'left_bases': array([0, 6]),
                    'right_bases': array([12, 22]),
                    'widths': array([3.21032149, 5.58696106]),
                    'width_heights': array([3228208.375, 3748458.5  ]),
                    'left_ips': array([ 8.12272812, 13.25125172]),
                    'right_ips': array([11.33304961, 18.83821278])},
    mz:    786.599123189169,
    list_retention_time:    array([66.26447883, 66.63433189, 67.01434098, 67.38420814, 67.7534237 ,
                    68.12477678, 68.49404578, 68.87324949, 69.25232779, 69.63137115,
                    70.01254658, 70.39144237, 70.77048133, 71.14956203, 71.52869926,
                    71.90967366, 72.27886773, 72.65788714, 73.02685195, 73.39578426,
                    73.76490146, 74.13592309, 74.51518654]),
    list_intensity:    array([  24545.924,   40612.44 ,  151511.27 ,  343555.5  ,  885063.8  ,
                    1232703.6  , 1938564.5  , 2306679.2  , 3195485.5  , 3462114.5  ,
                    3736445.2  , 3482002.5  , 2719971.5  , 3400839.8  , 4784387.5  ,
                    4500411.   , 5558352.5  , 4896509.   , 5039317.5  , 3499304.   ,
                    2457701.2  , 1694263.9  , 1377836.1  ], dtype=float32))

    min_prominence_threshold is unlikely one-size fits all. But an iterative detection will work better.
    prominence_window is not important, but is set to improve performance.
    min_timepoints is taken at mid-height, therefore 2 is set to allow very narrow peaks.
    gaussian_shape is minimal R^2 in Gaussian peak goodness_fitting.

    To-do
    -----
    further refine peak detection in future, 
    e.g. extra iteration of peak detection in remaining region; better deconvolution of overlapping peaks.

    '''
    def __init2__(self, mz, RT, INTS):
        self.mz, self.list_retention_time, self.list_intensity = [mz, RT, INTS]
        self.formula_mass = None
        self.raw_mzs = []                               # multiple values possible when merging traces
        self.__mass_corrected_by_asari__ = False        # True if updated by calibration
        self.features_assigned = False                  # tracking if assigned to Experiment features
        self.sample_name = ''

    def detect_peaks(self, min_intensity_threshold, min_timepoints, min_prominence_threshold, prominence_window, gaussian_shape):
        self.list_peaks = []
        peaks, properties = find_peaks(self.list_intensity, 
                                    height=min_intensity_threshold, width=min_timepoints, 
                                    prominence=min_prominence_threshold, wlen=prominence_window) 
        
        if peaks.size > 1:
            # Rerun, raising min_prominence_threshold. Not relying on peak shape for small peaks, 
            # because chromatography is often not good so that small peaks can't be separated from noise.
            #
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
                self.list_peaks.append(P)

        self.number_of_peaks = len(self.list_peaks)

    def serialize2(self):
        '''
        Placeholder
        '''
        d = {
            'number_of_peaks': self.number_of_peaks,
            'peaks': [P.id for P in self.list_peaks],
        }
        return d.update( self.serialize() )
        


#
# peak evaluation is in this class
#
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

    def serialize2(self):
        pass
