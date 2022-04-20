'''
MassGrid for correspondence, and FeatureList from feature/peak detection.
We use a similar concept of FeatureMap as in OpenMS here, 
but the correspondence algorithms take adavantage of high m/z resolution first, 
then utilizes MS1_pseudo spectra and cumulative elution profiles.

The use of CompositeMap facilitates data visualization and exploration.
'''

import pandas as pd
from mass2chem.search import *

from .mass_functions import *
from .chromatograms import *
from .peaks import *
from .samples import SimpleSample


class MassGrid:
    def __init__(self, cmap=None, experiment=None):
        self.experiment = experiment
        self.CMAP = cmap
        self.reference_sample_instance = self.CMAP.reference_sample_instance
        self.max_ref_rtime = self.CMAP.max_ref_rtime
        self.list_sample_names = self.CMAP.list_sample_names
        self._number_of_samples_ = self.CMAP._number_of_samples_
        
    def build_grid_sample_wise(self):
        '''
        Align one sample a time to reference m/z grid, based on their anchor m/z tracks.
        This is better for reliable assembly of small number of samples.
        '''
        self._initiate_mass_grid()
        sample_ids = self.experiment.valid_sample_ids
        sample_ids.pop(self.experiment.reference_sample_id)
        for sid in sample_ids:
            SM = SimpleSample(self.experiment.sample_registry[sid],
                experiment=self.experiment, database_mode=self.experiment.database_mode, mode=self.experiment.mode)
            self.add_sample(SM)


    def build_grid_by_centroiding(self):
        '''
        assemble mass grid by grouping m/z values to centroids.
        Each centroid can have no more than one mass track per sample.
        This is more efficient for large number of samples.
        
        '''
        # self.grid = np.array()

        self.MassGrid = pd.DataFrame(
            np.full((len(reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_sample_names,
        )




    def _initiate_mass_grid(self):
        '''
        Initiate MassGrid using reference sample
        '''
        reference_sample = self.reference_sample_instance
        _d = dict(zip(reference_sample.rt_numbers, reference_sample.rt_numbers))
        reference_sample.rt_cal_dict = reference_sample.reverse_rt_cal_dict = _d
        ref_list_mass_tracks = reference_sample.list_mass_tracks

        print("\nInitiating MassGrid, ...\n")
        self._mz_landmarks_ = reference_sample._mz_landmarks_
        reference_mzlist = [ x['mz'] for x in ref_list_mass_tracks ]
        # setting up DataFrame for MassGrid
        # not forcing dtype on DataFrame, to avoid unreported errors; convert to int when using MassGrid
        self.MassGrid = pd.DataFrame(
            np.full((len(reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_sample_names,
        )
        # Add ref mz as a column to MassGrid; ref mzlist will be dynamic updated in MassGrid["mz"]
        self.MassGrid['mz'] = reference_mzlist
        self.MassGrid[ reference_sample.name ] = [ x['id_number'] for x in ref_list_mass_tracks ]
        self.experiment.all_samples.append(reference_sample)


    def add_sample(self, sample, database_cursor=None):
        '''
        Add Sample instance to and update MassGrid; 
        add Sample to self.experiment.samples.

        recalculate_ref is not done here, because it's easier to include unmatched features from Sample.
        If needed, the reference m/z values should be updated by revisiting DB samples.
        But the recalculation should be based on calibrated m/z values so that they are consistent across samples.

        '''
        print("Adding sample to MassGrid,", sample.name)
        list_mass_tracks = sample.get_masstracks_and_anchors()
        mzlist = [x['mz'] for x in list_mass_tracks]
        new_reference_mzlist, new_reference_map2, updated_REF_landmarks, _r = landmark_guided_mapping(
                                    list(self.MassGrid['mz']), self._mz_landmarks_, mzlist, sample._mz_landmarks_)
        # print("_r = %f, new_reference_mzlist = %d" %(_r, len(new_reference_mzlist)))

        NewGrid = pd.DataFrame(
            np.full((len(new_reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_sample_names,
        )
        NewGrid[ :self.MassGrid.shape[0]] = self.MassGrid
        NewGrid['mz'] = new_reference_mzlist
        NewGrid[ sample.name ] = new_reference_map2
        self.MassGrid = NewGrid
        self._mz_landmarks_ = updated_REF_landmarks
        sample.mz_calibration_ratio = _r
        
        self.experiment.all_samples.append(sample)




    def join(self, M2):
        '''
        Join with another MassGridCluster.
        Using a common reference, which should be the 1st sample in both clusters.
        Return the merged MassGridCluster.
        '''

        pass










class CompositeMap:
    '''
    Each experiment is summarized into a CompositeMap (CMAP), as a master feature map.
    i) MassGrid: a matrix (DataFrame) for correspondence of mass tracks to each sample 
    ii) FeatureList: list of feature definitions, i.e. peaks defined on composite mass tracks.
    iii) FeatureTable: a matrix for feature intensities per sample

    Steps:
    1. Build MassGrid. Choose one reference from all samples for the largest number of landmark m/z tracks.
    2. Add each of remaining samples to MassGrid, 
    m/z alignment iguided by matched isotope/adduct pairs. Reference m/z is updated every time.
    RT alignment function is determined per sample, using selective landmark peaks.
    Default is a LOWESS function, but open to others to plugin.
    3. Build composite elution profile (composite_mass_tracks)
    by cumulative sum of mass tracks from all samples after RT correction.
    4. Global peak detection is performed on each composite massTrack.
    5. Mapping global peaks (i.e. features) back to all samples and extract sample specific peak areas.
    This completes the FeatureTable.
    6. Grouping features to empirical compounds (defined in metDataModel package).
    '''
    def __init__(self, experiment):
        '''
        Composite map of mass tracks and features, with pointers to individual samples.
        '''
        self.experiment = experiment
        self._number_of_samples_ = experiment.number_of_samples
        self.list_sample_names = [experiment.sample_registry[ii]['name'] for ii in experiment.valid_sample_ids]

        # designated reference sample; all RT is aligned to this sample
        self.reference_sample_instance = self.reference_sample = self.get_reference_sample_instance(experiment.reference_sample_id)
        self.rt_length = self.experiment.number_scans
        self.dict_scan_rtime = self.get_reference_rtimes(self.rt_length)
        self.max_ref_rtime = self.dict_scan_rtime[self.rt_length-1]

        self.MassGrid = None                        # will be DF
        self.FeatureTable = None
        self.FeatureList = []

        self._mz_landmarks_ = []                    # m/z landmarks as index numbers
        self.good_reference_landmark_peaks = []     # used for RT alignment and m/z calibration to DB
        # self.reference_mzdict = {}
        self.composite_mass_tracks = {}             # following MassGrid indices

    def get_reference_sample_instance(self, reference_sample_id):
        SM = SimpleSample(self.experiment.sample_registry[reference_sample_id],
                experiment=self.experiment, database_mode=self.experiment.database_mode, mode=self.experiment.mode)
        SM.list_mass_tracks = SM.get_masstracks_and_anchors()
        return SM

    def get_reference_rtimes(self, rt_length):
        '''
        Extrapolate retention time on self.reference_sample_instance to max scan number in the experiment.
        '''
        X, Y = self.reference_sample.rt_numbers, self.reference_sample.list_retention_time
        interf = interpolate.interp1d(X, Y, fill_value="extrapolate")
        newX = range(rt_length)
        newY = interf(newX)
        return dict(zip(newX, newY))


    def construct_mass_grid(self):
        '''
        MassGrid for whole experiment. Use sample name as column identifiers.
        All mass tracks are included at this stage, regardless if peaks are detected, because
        peak detection will be an improved process on the composite tracks.
        Number of samples dictate workflow: 
        build_grid_by_centroiding is fast, but build_grid_sample_wise is used for small studies 
        to compensate limited size for statistical distribution.
        '''
        MG = MassGrid( self, self.experiment )
        if self._number_of_samples_ <= self.experiment.parameters['project_sample_number_small']:
            MG.build_grid_sample_wise()
        else:
            MG.build_grid_by_centroiding()
           
        self.MassGrid = MG.MassGrid
        self._mz_landmarks_ = MG._mz_landmarks_


    def mock_rentention_alignment(self):
        for sample in self.experiment.all_samples[1:]:      # first sample is reference
            sample.rt_cal_dict, sample.reverse_rt_cal_dict = {}, {}


    def align_retention_time(self):
        '''
        Because RT will not match precisely btw samples, it's remapped to a common set of time coordinates.
        The default is a LOWESS algorithm.
        Also used in the field include dynamic time warp (DTW) and univariate spline.
        We saw no reason to use them, but people are welcome to implement alternative functions.

        Do alignment function using high-selectivity mass tracks.
        Step 1. get high-selectivity mass tracks among landmarks.
        2. for tracks of highest intensities, do quick peak detection to identify RT apexes.
        Only masstracks with single peaks will be used for RT alignment.
        3. use RT from 2, do LOWESS fit to reference RT values. 
        The fitted function will be recorded for each sample, 
        and applied to all RT scan numbers in the samples when used for CMAP construction.
        Important:
        sample.rt_cal_dict, sample.reverse_rt_cal_dict are kept for changed values only and set within sample RT boundaries.
        This is efficient by ignoring unnecessary tracking, and {} is consistent with samples without RT alignment.
        When samples fail in RT alignment,they are logged in warning and treated at the end as if no alignment is required.
        '''
        print("\nCalibrating retention time to reference, ...\n")
        cal_min_peak_height = self.experiment.parameters['cal_min_peak_height']
        MIN_PEAK_NUM = self.experiment.parameters['peak_number_rt_calibration']
        self.good_reference_landmark_peaks = self.set_RT_reference(cal_min_peak_height)
        for SM in self.experiment.all_samples[1:]:      # first sample is reference
            self.calibrate_sample_RT(SM, cal_min_peak_height=cal_min_peak_height, MIN_PEAK_NUM=MIN_PEAK_NUM)
        

    def calibrate_sample_RT(self, 
                                sample, 
                                calibration_fuction=rt_lowess_calibration, 
                                cal_min_peak_height=100000,
                                MIN_PEAK_NUM=15):
        '''
        Calibrate retention time per sample.
        This is based on a set of unambiguous peaks: quich peak detection on anchor mass trakcs, 
        and peaks that are unique to each track are used for RT alignment.

        When calibration_fuction fails, e.g. inf on lowess_predicted,
        it is assumed that this sample is not amendable to computational alignment,
        and the sample will be attached later without adjusting retention time.
        One can also use alternative workflow if so desired, 
        to do peak detection in all samples and correspond the peaks after.
        '''
        candidate_landmarks = [self.MassGrid[sample.name].values[
                                p['ref_id_num']] for p in self.good_reference_landmark_peaks] # contains NaN
        good_landmark_peaks, selected_reference_landmark_peaks = [], []
        list_mass_tracks = sample.get_masstracks_and_anchors()
        for jj in range(len(self.good_reference_landmark_peaks)):
            ii = candidate_landmarks[jj]
            if not pd.isna(ii):
                ii = int(ii)
                this_mass_track = list_mass_tracks[ii]
                Upeak = quick_detect_unique_elution_peak(this_mass_track['intensity'], 
                            min_peak_height=cal_min_peak_height, 
                            min_fwhm=3, min_prominence_threshold_ratio=0.2)
                if Upeak:
                    Upeak.update({'ref_id_num': ii})
                    good_landmark_peaks.append(Upeak)
                    selected_reference_landmark_peaks.append(self.good_reference_landmark_peaks[jj])

        _NN, _CALIBRATED = len(good_landmark_peaks), False
        # only do RT calibration if MIN_PEAK_NUM is met
        if _NN >  MIN_PEAK_NUM:
            try:
                sample.rt_cal_dict, sample.reverse_rt_cal_dict = calibration_fuction( 
                                            good_landmark_peaks, selected_reference_landmark_peaks, 
                                            sample.rt_numbers, self.reference_sample.rt_numbers, )
                _CALIBRATED = True
            except ValueError:
                pass
        if not _CALIBRATED:
                sample.rt_cal_dict, sample.reverse_rt_cal_dict =  {}, {}
                print("    ~warning~ Faluire in retention time alignment (%d); %s is used without alignment.\n" 
                                            %( _NN, sample.name))
                

    def set_RT_reference(self, cal_peak_intensity_threshold=100000):
        '''
        Start with the referecne samples, usually set for an initial sample of most peaks.
        Do a quick peak detection for good peaks; use high selectivity m/z to avoid ambiguity in peak definitions.

        Return 
        good_reference_landmark_peaks: [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...]
        '''
        selectivities = calculate_selectivity( self.MassGrid['mz'][self._mz_landmarks_], 
                                                self.experiment.parameters['mz_tolerance_ppm'])
        good_reference_landmark_peaks = []
        ref_list_mass_tracks = self.reference_sample.list_mass_tracks
        for ii in range(len(self._mz_landmarks_)):
            if selectivities[ii] > 0.99:
                ref_ii = self.MassGrid[self.reference_sample.name][self._mz_landmarks_[ii]]
                if ref_ii and not pd.isna(ref_ii):
                    this_mass_track = ref_list_mass_tracks[ int(ref_ii) ]
                    Upeak = quick_detect_unique_elution_peak(this_mass_track['intensity'], 
                                min_peak_height=cal_peak_intensity_threshold, 
                                min_fwhm=3, min_prominence_threshold_ratio=0.2)
                    if Upeak:
                        Upeak.update({'ref_id_num': self._mz_landmarks_[ii]}) # as in MassGrid index
                        good_reference_landmark_peaks.append(Upeak)

        return good_reference_landmark_peaks


    def global_peak_detection(self):
        '''
        Using peaks.deep_detect_elution_peaks on composite mass tracks.
        Results are deemed as features, because it's at Experiment level.
        Peak area and height are cumulated from all samples. 
        Not trying to average because some peaks are in only few samples.

        Performance can be improved further - quick filter to ROI will help long LC runs.

        '''
        self.composite_mass_tracks = self.make_composite_mass_tracks()
        print("Peak detection on %d composite mass tracks, ...\n" %len(self.composite_mass_tracks))
        self.FeatureList = batch_deep_detect_elution_peaks(
            self.composite_mass_tracks.values(), self.experiment.number_scans, self.experiment.parameters
        )

        ii = 0
        for peak in self.FeatureList:
            ii += 1
            peak['id_number'] = 'F'+str(ii)
            # convert scan numbers to rtime
            try:
                peak['rtime'] = self.dict_scan_rtime[peak['apex']]
            except KeyError:
                peak['rtime'] = self.max_ref_rtime                # imputed value set at max rtime
                print("Feature rtime out of bound - ", peak['id_number'], peak['apex'])
            try:
                peak['rtime_left_base'], peak['rtime_right_base'] = self.dict_scan_rtime[peak['left_base']
                                ], self.dict_scan_rtime[peak['right_base']]
            except KeyError:
                print("Feature rtime out of bound on", peak['id_number'], (peak['apex'], peak['left_base'], peak['right_base']))

        self.generate_feature_table()


    def make_composite_mass_tracks(self):
        '''
        Generate composite mass tracks by summing up signals from all samples after RT calibration.
        Each mass track from a sample is np.array of full RT range.
        The SM.rt_cal_dict only records changes to be made. 
        It's also empty if no RT alignment was done, in which case RT is assumed to match to reference.

        Return:
        Dict of dict {mz_index_number: {'id_number': k, 'mz': mzDict[k], 'intensity': v}, ...}
        '''
        mzDict = dict(self.MassGrid['mz'])
        mzlist = list(self.MassGrid.index)                          # this gets indices as keys, per mass track
        basetrack = np.zeros(self.rt_length, dtype=np.int64)        # self.rt_length defines max rt number
        _comp_dict = {}
        for k in mzlist: 
            _comp_dict[k] = basetrack.copy()

        for SM in self.experiment.all_samples:
            list_mass_tracks = SM.get_masstracks_and_anchors()
            # if SM.rt_cal_dict:
            for k in mzlist:
                ref_index = self.MassGrid[SM.name][k]
                if not pd.isna(ref_index): # ref_index can be NA 
                    _comp_dict[k] += remap_intensity_track(list_mass_tracks[int(ref_index)]['intensity'], 
                                                            basetrack.copy(), SM.rt_cal_dict)


        result = {}
        for k,v in _comp_dict.items():
            result[k] = { 'id_number': k, 'mz': mzDict[k], 'intensity': v }

        return result


    def generate_feature_table(self):
        '''
        '''
        FeatureTable = pd.DataFrame(self.FeatureList)
        for SM in self.experiment.all_samples:
            FeatureTable[SM.name] = self.extract_features_per_sample(SM)

        self.FeatureTable = FeatureTable


    def extract_features_per_sample(self, sample):
        '''
        watch for range due to calibration/conversion.

        '''
        fList = []
        mass_track_map = self.MassGrid[sample.name]
        list_mass_tracks = sample.get_masstracks_and_anchors()
        for peak in self.FeatureList:
            track_number = mass_track_map[peak['parent_masstrack_id']]
            peak_area = 0
            if not pd.isna(track_number):           # watch out dtypes
                mass_track = list_mass_tracks[ int(track_number) ]
                # ?? 
                left_base = sample.reverse_rt_cal_dict.get(peak['left_base'], peak['left_base'])
                right_base = sample.reverse_rt_cal_dict.get(peak['right_base'], peak['right_base'])
                peak_area = mass_track['intensity'][left_base: right_base+1].sum()

            fList.append( peak_area )

        return fList







class SampleAssembly(CompositeMap):
    '''
    
    
    Alternative workflow:
    peak detection in each sample, then align by distance threshold.



    '''


    def align_samples(self, ):
        pass



