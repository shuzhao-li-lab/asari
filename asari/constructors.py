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
        self._number_of_valid_samples_ = len(self.list_sample_names)

        # designated reference sample; all RT is aligned to this sample
        self.reference_sample_instance = self.reference_sample = self.get_reference_sample_instance(experiment.reference_sample_id)
        self.rt_length = len(self.reference_sample.rt_numbers)

        self.MassGrid = None                        # will be DF
        self.FeatureTable = None
        self.FeatureList = []

        self._mz_landmarks_ = []                    # m/z landmarks as index numbers
        self.good_reference_landmark_peaks = []     # used for RT alignment and m/z calibration to DB
        # 
        self.reference_mzdict = {}
        self.composite_mass_tracks = {}             # following MassGrid indices


    def get_reference_sample_instance(self, reference_sample_id):
        SM = SimpleSample(self.experiment.sample_registry[reference_sample_id],
                experiment=self.experiment, database_mode=self.experiment.database_mode, mode=self.experiment.mode)
        SM.list_mass_tracks = SM.get_masstracks_and_anchors()
        self.dict_scan_rtime = dict(zip(SM.rt_numbers, SM.list_retention_time))
        self.max_ref_rtime = max(SM.list_retention_time)
        return SM

    def construct_mass_grid(self):
        '''
        MassGrid for whole experiment. Use sample name as column identifiers.
        All mass tracks are included at this stage, regardless if peaks are detected, because
        peak detection will be an improved process on the composite tracks.
        '''
        if self._number_of_valid_samples_ <= self.experiment.parameters['project_sample_number_small']:
            self._initiate_mass_grid()
            sample_ids = self.experiment.valid_sample_ids
            sample_ids.pop(self.experiment.reference_sample_id)
            for sid in sample_ids:
                SM = SimpleSample(self.experiment.sample_registry[sid],
                    experiment=self.experiment, database_mode=self.experiment.database_mode, mode=self.experiment.mode)

                self.add_sample(SM)
                
        elif self._number_of_valid_samples_ <= self.experiment.parameters['project_sample_number_large']:

            MGC = MassGridCluster(  )
            self.MassGrid = MGC.grid()


        else:   # split and do batch build

            pass


    def _initiate_mass_grid(self):
        '''
        Initiate MassGrid using reference sample
        '''
        reference_sample = self.reference_sample_instance
        _d = dict(zip(reference_sample.rt_numbers, reference_sample.rt_numbers))
        reference_sample.rt_cal_dict = reference_sample.reverse_rt_cal_dict = _d
        ref_list_mass_tracks = reference_sample.list_mass_tracks
        self.experiment.number_scans = max(reference_sample.rt_numbers)

        print("\nInitiating MassGrid, ...\n    The reference sample is:\n    ||* %s *||\n" %reference_sample.name)
        print("Max _retention_time is %4.2f at scan number %d.\n" %(self.max_ref_rtime,
                                                                    max(reference_sample.rt_numbers)))

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
        
        self.experiment.number_scans = max(self.experiment.number_scans, max(sample.rt_numbers))
        self.experiment.all_samples.append(sample)


    def mock_rentention_alignment(self):
        for sample in self.experiment.all_samples[1:]:      # first sample is reference
            sample.rt_cal_dict, sample.reverse_rt_cal_dict = mock_rt_calibration(sample.rt_numbers, self.reference_sample.rt_numbers)

    def align_retention_time(self):
        '''
        Because RT will not match precisely btw samples, it's remapped to a common set of time coordinates.
        The default is a LOWESS algorith, while 
        dynamic time warp (DTW) and univariate spline are also used in the field.

        Do alignment function using high-selectivity mass tracks.
        Step 1. get high-selectivity mass tracks among landmarks.
        2. for tracks of highest intensities, do quick peak detection to identify RT apexes.
        Only masstracks with single peaks will be used for RT alignment.
        3. use RT from 2, do LOWESS fit to reference RT values. 
        The fitted function will be recorded for each sample, 
        and applied to all RT scan numbers in the samples when used for CMAP construction.

        sample.rt_cal_dict, sample.reverse_rt_cal_dict are kept for changed values only and set within sample RT boundaries.
        
        Marked samples that fail in RT alignment, and deal with them at the end.

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
        Calibrate retention time, via spline func, per sample
        Use anchor mass trakcs, do quick peak detection. Use tracks of single peaks for RT alignment.
        This produces a new set of RT coordinates (intensity values shifted along the RT, no need to change).

        Because RT is using scan numbers, landmarks can overlap, e.g. rt_cal:
        (55, 55), (56, 56), (56, 57), (56, 59), (57, 55), (57, 59), (58, 60), (60, 61), (61, 59), (61, 61), (62, 62), 
        (63, 63), (67, 67), (69, 69), (69, 70), (70, 70), (71, 71), (72, 72), (73, 72), (73, 74), (74, 75), (76, 75), (76, 78), (77, 75), (77, 77), ...,
        (190, 190), (190, 191), (190, 192), (191, 189), (191, 191), (191, 192), (192, 192), (192, 193),...
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

        _NN = len(good_landmark_peaks)
        # only do RT calibration if MIN_PEAK_NUM is met
        if _NN >  MIN_PEAK_NUM:
            try:
                sample.rt_cal_dict, sample.reverse_rt_cal_dict = calibration_fuction( 
                                            good_landmark_peaks, selected_reference_landmark_peaks, 
                                            sample.rt_numbers, self.reference_sample.rt_numbers, )
            except ValueError:
                print("    ~warning~ Faluire in retention time alignment, %s" %sample.name)
        else:
            print("    ~Warning~ Faluire in retention time alignment on %s, due to too few aligned features (%d)." 
                                %(sample.name, _NN))


    def set_RT_reference(self, cal_peak_intensity_threshold=100000):
        '''
        Start with the referecne samples, usually set for an initial sample of most peaks.
        Do a quick peak detection for good peaks; use high selectivity m/z to avoid ambiguity in peak definitions.

        Return 
        good_reference_landmark_peaks: [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...]
        '''
        selectivities = calculate_selectivity( self.MassGrid['mz'][self._mz_landmarks_], 
                                                self.experiment.parameters['mz_tolerance'])
        good_reference_landmark_peaks = []
        ref_list_mass_tracks = self.reference_sample.list_mass_tracks
        for ii in range(len(self._mz_landmarks_)):
            if selectivities[ii] > 0.99:
                ref_ii = self.MassGrid[self.reference_sample.name][self._mz_landmarks_[ii]]
                if ref_ii:
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
        Peak area and height are cumulated from all samples. Not trying to average because some peaks are in only few samples.

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
                peak['rtime_left_base'], peak['rtime_right_base'] = self.dict_scan_rtime[peak['left_base']], self.dict_scan_rtime[peak['right_base']]
            except KeyError:
                peak['rtime'] = self.max_ref_rtime                                  # imputed value set at max rtime
                print("Peak rtime out of bound on", ii)

        self.generate_feature_table()


    def make_composite_mass_tracks(self):
        '''
        Generate composite mass tracks by summing up signals from all samples after RT calibration.

        Start a RT_index: intensity dict, and cummulate intensity values. 
        Convert back to rtlist, intensities at the end.
        Because RT may be distorted, e.g. from sample to ref mapping: 
            {384: 387, 385: 387, 386: 388, 387: 388, 388: 389, 389: 389, 390: 390, 391: 390, },
        smoothing is needed in doing the composite.

        Return:
        Dict {mz_index_number: intensity_track, ...}
        Dict of dict {mz_index_number: {'rt_scan_numbers': xx, 'intensity':yy }, ...}
        '''
        mzDict = dict(self.MassGrid['mz'])
        mzlist = list(self.MassGrid.index)                          # this gets indices as keys, per mass track
        basetrack = np.zeros(self.rt_length, dtype=np.int64)        # self.rt_length defines max rt number
        _comp_dict = {}
        for k in mzlist: 
            _comp_dict[k] = basetrack.copy()

        for SM in self.experiment.all_samples:
            list_mass_tracks = SM.get_masstracks_and_anchors()
            if SM.rt_cal_dict:
                for k in mzlist:
                    ref_index = self.MassGrid[SM.name][k]
                    if not pd.isna(ref_index): # ref_index and 
                        _comp_dict[k] += remap_intensity_track(list_mass_tracks[int(ref_index)]['intensity'], 
                                                                basetrack.copy(), SM.rt_cal_dict)
        result = {}
        for k,v in _comp_dict.items():
            result[k] = { 'id_number': k, 'mz': mzDict[k], 'intensity': v }

        return result


    def generate_feature_table(self):
        '''
        cmap_header = ['feature_id', 'cmap_masstrack_id', 'mz', 'rtime', 'rt_min', 'rt_max', 
                'peak_quality', 'selectivity_mz', 'height_mean', 'peak_area_mean', ]
        header = cmap_header + self.list_sample_names
        FeatureTable = pd.DataFrame( 
                np.zeros((len(self.FeatureList), len(header))), columns=header, )
        '''
        FeatureTable = pd.DataFrame(self.FeatureList)
        for SM in self.experiment.all_samples:
            if SM.rt_cal_dict:
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


class MassGridCluster:
    def __init__(self):
        
        self.grid = np.array()

    def build_grid(self):
        '''
        from [(mz, track_id, sample_num), ...]
        assemble mass grid
        
        '''
        pass


    def join(self, M2):
        '''
        Join with another MassGridCluster.
        Using a common reference, which should be the 1st sample in both clusters.
        Return the merged MassGridCluster.
        '''

        pass
