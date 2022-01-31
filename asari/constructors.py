'''

Empirical Compound, list and tree
=================================
Constructing trees of emperical compounds (epdTree) from a list of peaks or features,
which follow a format:
list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...]

Steps:
1. find isotopic signatures from list_peaks
2. initiate epdTree classes using the above isotopic signatures; extend by searching common adducts
3. In remaining peaks (isotopes not seen due to low intensity etc), search pairs for primary adducts. 
   Initiate epdTree classes using those.
4. Extend to other adducts for empCpds from the above steps.
5. Reconsolidate overlap epdTrees. E.g. combinations of isotopes and adducts are found separately.

Notes:
a. Common anchor peaks are M+H or M-H etc, but will calculate at later round. M* may not show.
b. Adducts and isotopes are combinatorial, under restriction of chemical formulae.
c. We curated isotopic/adduct patterns in this package.
d. epdTree can be improved in future based on real data statistics and more structured cheminformatics.

from scipy.interpolate import UnivariateSpline
from .utils import *
from .sql import *
'''

import os
import pandas as pd

from .search import *
from .mass_functions import *
from .chromatograms import *
from .peaks import *


class CompositeMap:
    '''
    Each experiment is summarized into a CompositeMap (CMAP), as a master feature map.
    i) MassGrid: a matrix (DataFrame) for correspondence of mass tracks to each sample 
    ii) FeatureList: list of feature definitions, i.e. peaks defined on composite mass tracks.
    iii) FeatureTable: a matrix for feature intensities per sample

    Steps:

    1. Build Mass Grid first. 
    Choose one reference from the initial samples with largest number of landmark tracks.

    2. Add each of remaining samples to MassGrid.
    Alignment is guided by matched isotope/adduct pairs. Reference m/z is updated every time.

    3. Optimize mass alignment.

    4. Determine RT alignment function per sample, using selective landmark peaks.
    Default is a LOWESS function, but open to others to plugin.

    5. Build composite elution profile (composite_mass_tracks)
    by cumulative sum of mass tracks from all samples after RT correction.

    6. Global peak detection is performed on each composite massTrack, by two rounds -
    stepping down in prominence. Smoothing (moving average) is performed on tracks with more that 
    two peaks detected (considered noisy), followed by a fresh round of peak detection.

    7. Mapping global peaks (i.e. features) back to all samples and extract sample specific peak areas.
    This completes the FeatureTable.

    8. Annotation, grouping features to epds, and DB matching. Optional, m/z calibration to refDB.

    '''

    def __init__(self, experiment):
        '''
        Composite map of mass tracks and features, with pointers to individual samples.

        Workflow is in `workflow.ext_Experiment.process_all`.
        '''
        self.experiment = experiment
        self._number_of_samples_ = experiment.number_of_samples
        self.list_input_files = experiment.list_input_files
        self.reference_sample = None                # designated reference sample; all RT is aligned to this sample

        self.MassGrid = None                        # will be DF
        self.FeatureTable = None
        self.FeatureList = []

        self._mz_landmarks_ = []                    # keeping anchor pairs as landmarks
        #
        # self.ref_empCpds = []
        self.reference_mzdict = {}
        self.composite_mass_tracks = {}             # following MassGrid indices


    def construct_mass_grid(self, init_Samples):
        self.initiate_mass_grid(init_Samples)
        # print("Done initiation_Samples.\n")
        for f in self.list_input_files:                                 # run remaining samples, 
            if f not in self.experiment.initiation_samples:
                SM = self.experiment.process_single_sample(f)
                # not via DB
                self.add_sample(SM)


    def initiate_mass_grid(self, init_samples):
        '''
        Create MassGrid for whole experiment
        Use list_samples (list of Sample instances) to initiate MassGrid.
        Each Sample has list_mass_tracks and anchor_mz_pairs.
        The m/z of reference is taken as mean of these initiation samples (id_mz_str not updated). Afterwards, 
        the ref m/z isn't changed when single samples are added for performance/complexity reasons.

        All mass traces are included at this stage, regardless if peaks are detected, because
        peak detection will be an improved process on the composite traces.

        Start by matching anchor pairs, then work thru remaining traces.
        1. create a reference based on anchor pairs
        2. align each sample to the reference_anchor_pairs

        '''
        # initiation using the Sample of most _number_anchor_mz_pairs_
        tmp, other_list_samples = [init_samples[0]], []
        for SM in init_samples[1:]:
            if SM._number_anchor_mz_pairs_ > tmp[0]._number_anchor_mz_pairs_:
                other_list_samples.append(tmp[0])
                tmp = [SM]
            else:
                other_list_samples.append(SM)

        reference_sample = tmp[0]

        _d = dict(zip(reference_sample.rt_numbers, reference_sample.rt_numbers))
        reference_sample.rt_cal_dict = reference_sample.reverse_rt_cal_dict = _d
        self.reference_sample = self.experiment.reference_sample = reference_sample
        # note: other samples are aligned to this ref

        print("\nInitiating MassGrid, ...\n    The reference sample is:\n    ||* %s *||\n" %reference_sample.input_file)

        self.reference_anchor_pairs = reference_sample.anchor_mz_pairs
        self._mz_landmarks_ = flatten_tuplelist(reference_sample.anchor_mz_pairs)
        reference_mzlist = [ x['mz'] for x in reference_sample.list_mass_tracks ]
        # setting up DataFrame for MassGrid
        # not forcing dtype on DataFrame, to avoid unreported errors; convert to int when using MassGrid
        self.MassGrid = pd.DataFrame(
            np.full((len(reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_input_files,
        )
        # Add ref mz as a column to MassGrid; ref mzlist will be dynamic updated in MassGrid["mz"]
        self.MassGrid['mz'] = reference_mzlist
        self.MassGrid[ reference_sample.input_file ] = [ x['id_number'] for x in reference_sample.list_mass_tracks ]
        # self.reference_sample.export_mass_traces()
        
        for SM in other_list_samples:
            self.add_sample(SM, database_cursor=None)

        
    def add_sample(self, sample, database_cursor=None):
        '''
        Add Sample instance to and update MassGrid; 
        add Sample to self.experiment.samples.

        To add: push each sample to SQLDB, - database_cursor;

        recalculate_ref is not done here, because it's easier to include unmatched features from Sample.
        If needed, the reference m/z values should be updated by revisiting DB samples.
        But the recalculation should be based on calibrated m/z values so that they are consistent across samples.

        '''
        print("Adding sample to MassGrid,", os.path.basename(sample.input_file))

        mzlist = [x['mz'] for x in sample.list_mass_tracks]
        new_reference_mzlist, new_reference_map2, updated_REF_landmarks, _r = landmark_guided_mapping(
                                    list(self.MassGrid['mz']),   self._mz_landmarks_, mzlist, sample._mz_landmarks_)
        # print("_r = %f, new_reference_mzlist = %d" %(_r, len(new_reference_mzlist)))

        NewGrid = pd.DataFrame(
            np.full((len(new_reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_input_files,
        )
        NewGrid[ :self.MassGrid.shape[0]] = self.MassGrid
        NewGrid['mz'] = new_reference_mzlist
        NewGrid[ sample.input_file ] = new_reference_map2
        self.MassGrid = NewGrid
        self._mz_landmarks_ = updated_REF_landmarks
        sample.mz_calibration_ratio = _r
        self.experiment.samples_nonreference.append(sample)
        self.experiment.number_scans = max(self.experiment.number_scans, max(sample.rt_numbers))


    def optimize_mass_grid(self):
        '''
        This inspects split or misaligned regions of m/z tracks.
        Problematic definition will persist in each sample's mass tracks,
        e.g. about 50 mass tracks are too close to each other in our plasma data; 200 in Zurich E. coli data.

        merge close mass tracks unless MassGrid suggests multiple features; 
        for split tracks in same samples, replace them by merged tracks and new id_numbers.

        to implement -
        For low-selectivity mass tracks, could do 2-D deconvolution.

        Not really needed now, after mass_tracks in each sample are cleaned up for m/z overlap.

        '''
        # also check for empty tracks after composite?

        # watch out for indices used by other variables, e.g. _mz_landmarks_
        pass



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
        
        Marked samples that fail in RT alignment, and deal with them at the end.

        '''
        print("\nCalibrating retention time to reference, ...\n")
        self.good_reference_landmark_peaks = self.set_RT_reference()
        for SM in self.experiment.samples_nonreference:
            self.calibrate_sample_RT(SM)
        

    def calibrate_sample_RT(self, 
                                sample, 
                                calibration_fuction=rt_lowess_calibration, 
                                cal_peak_intensity_threshold=100000,
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
        candidate_landmarks = [self.MassGrid[sample.input_file].values[
                                p['ref_id_num']] for p in self.good_reference_landmark_peaks] # contains NaN
        good_landmark_peaks, selected_reference_landmark_peaks = [], []
        for jj in range(len(self.good_reference_landmark_peaks)):
            ii = candidate_landmarks[jj]
            if not pd.isna(ii):
                ii = int(ii)
                this_mass_track = sample.list_mass_tracks[ii]
                rt_numbers, list_intensity = this_mass_track['rt_scan_numbers'], this_mass_track['intensity']
                # continuity in rt_scan_numbers is implemented in chromatograms.extract_single_track_ 
                Upeak = quick_detect_unique_elution_peak(rt_numbers, list_intensity, 
                            min_intensity_threshold=cal_peak_intensity_threshold, 
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
                print("    ~warning~ Faluire in retention time alignment, %s" %sample.input_file)
        else:
            print("    ~Warning~ Faluire in retention time alignment on %s, due to too few aligned features (%d)." 
                                %(sample.input_file, _NN))


    def set_RT_reference(self, cal_peak_intensity_threshold=100000):
        '''
        Start with the referecne samples, usually set for an initial sample of most peaks.
        Do a quick peak detection for good peaks; use high selectivity m/z to avoid ambiguity in peak definitions.

        Return 
        good_reference_landmark_peaks: [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...]
        '''
        selectivities = calculate_selectivity( self.MassGrid['mz'][self._mz_landmarks_] )
        good_reference_landmark_peaks = []
        for ii in range(len(self._mz_landmarks_)):
            if selectivities[ii] > 0.99:
                ref_ii = self.MassGrid[self.reference_sample.input_file][self._mz_landmarks_[ii]]
                if ref_ii:
                    this_mass_track = self.reference_sample.list_mass_tracks[ int(ref_ii) ]
                    rt_numbers, list_intensity = this_mass_track['rt_scan_numbers'], this_mass_track['intensity']
                    # continuity in rt_scan_numbers is implemented in chromatograms.extract_single_track_ 
                    Upeak = quick_detect_unique_elution_peak(rt_numbers, list_intensity, 
                                min_intensity_threshold=cal_peak_intensity_threshold, 
                                min_fwhm=3, min_prominence_threshold_ratio=0.2)
                    if Upeak:
                        Upeak.update({'ref_id_num': self._mz_landmarks_[ii]}) # as in MassGrid index
                        good_reference_landmark_peaks.append(Upeak)

        return good_reference_landmark_peaks


    def match_ref_db(self):

        pass


    def global_peak_detection(self):
        '''
        # print(self.composite_mass_tracks[55])
        
        '''
        self.composite_mass_tracks = self.make_composite_mass_tracks()
        print("\nPeak detection on %d composite mass tracks, ...\n" %len(self.composite_mass_tracks))

        for _, mass_track in self.composite_mass_tracks.items():
            self.FeatureList +=  deep_detect_elution_peaks( mass_track 

                    )           # to specify parameters here according to Experiment parameters

        ii = 0
        for peak in self.FeatureList:
            ii += 1
            peak['id_number'] = 'F'+str(ii)

        self.generate_feature_table()

    def make_composite_mass_tracks(self):
        '''
        Generate composite mass tracks by summing up signals from all samples after RT calibration.

        Start a RT_index: intensity dict, and cummulate intensity values. 
        Convert back to rtlist, intensities at the end.
        Because RT may be distorted, e.g. from sample to ref mapping: 
            {384: 387,
                385: 387,
                386: 388,
                387: 388,
                388: 389,
                389: 389,
                390: 390,
                391: 390, },
        smoothing is needed in doing the composite.

        Return:
        Dict of dict {mz_index_number: {'rt_scan_numbers': xx, 'intensity':yy }, ...}
        '''
        mzDict = dict(self.MassGrid['mz'])
        mzlist = list(self.MassGrid.index)              # this gets indices as keys, per mass track
        _comp_dict = {}
        for k in mzlist: 
            _comp_dict[k] = {}       # will be {rt_index: intensity}
            # init using reference_sample
            ref_index = self.MassGrid[self.reference_sample.input_file][k]
            if not pd.isna(ref_index):
                this_mass_track = self.reference_sample.list_mass_tracks[int(ref_index)]
                # ref mass track should be continuous
                _comp_dict[k] = dict(zip(this_mass_track['rt_scan_numbers'], this_mass_track['intensity']))

        for SM in self.experiment.samples_nonreference:
            if SM.rt_cal_dict:
                for k in mzlist:
                    ref_index = self.MassGrid[SM.input_file][k]
                    if not pd.isna(ref_index): # ref_index and 
                        this_mass_track = SM.list_mass_tracks[int(ref_index)]
                        # after remapping RT indices, no guaranty rt/intensity vector is continuous; thus smoothing is needed
                        remapped_rt = [ SM.rt_cal_dict[ii] for ii in this_mass_track['rt_scan_numbers'] ] # convert to ref RT
                        _comp_dict[k] = sum_dict( _comp_dict[k], smooth_rt_intensity_remap(remapped_rt, this_mass_track['intensity']))

        # reformat result
        result = {}
        for k,vdict in _comp_dict.items():
            if vdict:
                rt_scan_numbers = sorted(list(vdict.keys()))
                result[k] = { 'id_number': k, 'mz': mzDict[k], 'rt_scan_numbers': rt_scan_numbers, 
                                                    'intensity': [vdict[x] for x in rt_scan_numbers], }
            else:
                # From samples that failed to align 
                print("    ... mass track in unaligned sample dropped ... %4.4f ... " %mzDict[k])
        return result



    def generate_feature_table(self):
        '''
        cmap_header = ['feature_id', 'cmap_masstrack_id', 'mz', 'rtime', 'rt_min', 'rt_max', 
                'peak_quality', 'selectivity_mz', 'height_mean', 'peak_area_mean', ]
        header = cmap_header + self.list_input_files
        FeatureTable = pd.DataFrame( 
                np.zeros((len(self.FeatureList), len(header))), columns=header, )

        '''
        FeatureTable = pd.DataFrame(self.FeatureList)
        for SM in [self.experiment.reference_sample] + self.experiment.samples_nonreference:
            if SM.rt_cal_dict:
                FeatureTable[SM.input_file] = self.extract_features_per_sample(SM)

        self.FeatureTable = FeatureTable


    def extract_features_per_sample(self, sample):
        '''
        rt_scan_numbers can be out of range due to calibration/conversion. Fill with max.
        '''
        fList = []
        mass_track_map = self.MassGrid[sample.input_file]
        max_rt_number = max(sample.rt_numbers)
        for peak in self.FeatureList:
            track_number = mass_track_map[peak['parent_masstrack_id']]
            peak_area = 0
            if not pd.isna(track_number):           # watch out dtypes
                mass_track = sample.list_mass_tracks[ int(track_number) ]
                left_base = sample.reverse_rt_cal_dict[peak['left_base']]
                try:
                    right_base = sample.reverse_rt_cal_dict[peak['right_base']]
                except KeyError:
                    right_base = max_rt_number
                    print("    ... in %s ... incomplete elution peak at ... %4.4f ..." %(os.path.basename(sample.input_file), mass_track['mz']))

                for ii in range(len(mass_track['rt_scan_numbers'])):
                    if left_base <= mass_track['rt_scan_numbers'][ii] <= right_base:
                        peak_area += mass_track['intensity'][ii]

            fList.append( peak_area )

        return fList






class epdsConstructor:
    '''
    To organize a table of peaks/features into a list of empirical compounds
    (https://github.com/shuzhao-li/metDataModel).
    empCpds (or epds) can be initiated by signatures on isotopic relationships or adduct relationships.
    For low-intensity peaks, their isotopic counterparts may not be detectable.
    Only common adducts are considered at this step.
    This list of epds are used in m/z alignment/correspondence.
    For assigned peaks/features, will calculate selectivity/score later too.

    Input
    =====
    isotopic_patterns: [(1.003355, '13C/12C', (0, 0.2)), ...],
    peak_list: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
    'left_base': 648, 'right_base': 655, 'id_number': 555},
                 ...]
    
    To get list of empCpds, run
    ECCON = epdsConstructor(list_peaks)
    list_empCpds = ECCON.peaks_to_epds()
    # e.g. {'id': 358, 'list_peaks': [(4215, 'anchor'), (4231, '13C/12C'), (4339, 'anchor,+NH4')]},
    '''

    def __init__(self, peak_list, mode='pos'):
        # will add parameters dict
        self.peak_list = peak_list
        self.peak_dict = build_peak_id_dict(self.peak_list)
        self.mode = mode

    def peaks_to_epds(self):
        '''
        Anchor peak is the most abundant isotopic peak.
        Modification peaks are from adducts, neutral loss and fragments.
        Requires at least one of isotopic_signatures and adduct_signatures. Use peak numbers only.
        If initiated by isotopic_signatures, do adduct search;
        elif initiated by adduct_signatures, assuming other primary isotopes are not seen.
        
        We build indexed centurion trees here to assist searches.

        Return
        ======
        list_empCpds, [{'id': ii, 'list_peaks': [(peak_id, ion), (), ...],}, ...]
            The first peak is anchor ion.
        '''
        print("\n\nAnnotating empirical compounds on %d features/peaks, ..." %len(self.peak_list))
        list_empCpds = []
        epds, found = [], []
        mztree = build_centurion_tree(self.peak_list)
        isosignatures = find_isotopic_signatures(self.peak_list, mztree, isotopic_patterns)
        # [[(182, 'anchor'), (191, '13C/12C'), (205, '18O/16O')], ...]
        for L in isosignatures:
            found += [x[0] for x in L]

        print("Round 1 - numbers of epds and included peaks: ", (len(isosignatures), len(found)))
        remaining_peaks = [P for P in self.peak_list if P['id_number'] not in found]
        mztree = build_centurion_tree(remaining_peaks)
        for G in isosignatures:
            epds.append(
                self.extend_isosignatures_by_adducts(G, mztree, common_adducts[self.mode])
            )
        # reset remaining_peaks and mztree
        found2 = []
        for L in epds:
            found2 += [x[0] for x in L]
        remaining_peaks = [P for P in remaining_peaks if P['id_number'] not in found2]
        mztree = build_centurion_tree(remaining_peaks)
        _NN2 = len(found2)

        print("Round 2 - numbers of epds and included peaks: ", (len(epds), _NN2))
        # do de novo adduct initiations
        adduct_signatures = find_adduct_signatures(remaining_peaks, mztree, common_adducts[self.mode])
        for G in adduct_signatures:
            epds.append(G)
            _NN2 += len(G)

        print("Round 3 - numbers of epds and included peaks: ", (len(epds), _NN2))
        # the now remaining are to be assigned or singletons
        for ii in range(len(epds)):
            list_empCpds.append(
                {'id': ii, 'list_peaks': epds[ii]}
            )

        return list_empCpds

    def extend_isosignatures_by_adducts(self, isotopic_peaks, mztree, adduct_patterns, mz_tolerance_ppm=5):
        '''
        isotopic_peaks: e.g. [(182, 'anchor'), (191, '13C/12C'), (205, '18O/16O')]
        mztree: reset for remaining peaks

        Return one epd (empirical compound) based on initiation isotopic peaks and extended by adduct_patterns,
        as a list of [(peak_id, relation), ...]
        '''
        matched = isotopic_peaks
        for pp in isotopic_peaks:
            # use self.peak_dict to retrieve peaks
            P1, relation = self.peak_dict[pp[0]], pp[1]
            for adduct in adduct_patterns:
                # (1.0078, 'H'), (21.9820, 'Na/H'), ...
                tmp = find_all_matches_centurion_indexed_list( P1['mz'] + adduct[0], mztree, mz_tolerance_ppm )
                for P2 in tmp:
                    if is_coeluted(P1, P2):
                        matched.append( (P2['id_number'], relation +','+ adduct[1]) )
        return matched

    def extend_empCpds_by_adducts(self, seed_list_empCpds, list_peaks, adduct_patterns, mz_tolerance_ppm=5):
        '''
        Search list_peaks for adducts that fit patterns relative to anchor ions (?? ) in existing empCpds.
        Co-elution required.

        in progress -

        '''
        
        mztree = build_centurion_tree(list_peaks)
        for EPD in seed_list_empCpds:
            P1, relation = EPD['list_peaks'][0]
            for adduct in adduct_patterns:
                # (1.0078, 'H'), (21.9820, 'Na/H'), ...
                tmp = find_all_matches_centurion_indexed_list( P1['mz'] + adduct[0], mztree, mz_tolerance_ppm )
                for P2 in tmp:
                    if is_coeluted(P1, P2):
                        EPD['list_peaks'].append( (P2['id_number'], relation +','+ adduct[1]) )

        return seed_list_empCpds

