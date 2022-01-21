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


from .utils import *
from .sql import *
'''

from scipy.stats import norm as normal_distribution

import pandas as pd

from .search import *
from .functions import *


class CompositeMap:
    '''
    Each experiment is summarized into a CompositeMap (CMAP), as a master feature map.

    Use DataFrames to hold 
    i) MassGrid: a matrix for correspondence of mass tracks to each sample 
    ii) FeatureGrid: a matrix for feature-peak correspondence


    Steps:

    1. Build Mass Grid first. Starting with the initial samples.
    mass correspondence is unique by using mass tracks. One grid can have maxium one mass track in one sample.
    Reference m/z will be concensus among initial samples.
    Alignment is guided by matched isotope/adduct pairs.
    2. Add remaining samples to MassGrid, same anchor guided alignment.
    Reference m/z is not recalculated for performance reason.

    3. Optional, m/z calibration to refDB.

    4. Determine RT DTW function per sample, by matching to selective epd pairs.
    5. Build composite elution profile by cumulative sum of mass traces from all samples after RT correction.
    6. New samples are added to CMap, both aligning to the CMap and augment the CMap by new features.
    7. After all samples are processed, peak detection is performed carefully on each massTrace in CMap, 
    to determine number and range of peaks.
    8. Report peaks per sample based on the peak definition in CMap.
    9. New round of Annotation; add remaining features to epds or group into new epds.


    '''

    def __init__(self, experiment):
        '''
        Composite map of mass tracks and features, with pointers to individual samples.
            self.reference_anchor_pairs = []            # (mz_str, mz_str), use mz str as identifiers 
        '''
        
        self.experiment = experiment
        self._number_of_samples_ = experiment.number_of_samples
        self.list_input_files = experiment.list_input_files

        self.MassGrid = None                        # will be DF
        self.FeatureGrid = None
        self._mz_landmarks_ = []                      # keeping anchor pairs as landmarks
        #
        self.reference_mzdict = {}
        self.ref_empCpds = []
        

    def initiate_mass_grid(self, list_samples, recalculate_ref=False):
        '''
        Use list_samples (list of Sample instances) to initiate MassGrid.
        Each Sample has list_mass_tracks and anchor_mz_pairs.
        The m/z of reference is taken as mean of these initiation samples (id_mz_str not updated). Afterwards, 
        the ref m/z isn't changed when single samples are added for performance/complexity reasons.

        All mass traces are included at this stage, regardless if peaks are detected, because
        peak detection will be an improved process on the composite traces.

        Start by matching anchor pairs, then work thru remaining traces.
        1. create a reference based on anchor pairs
        2. align each sample to the reference_anchor_pairs

        MassGrid = {id_mz_str: [massTrack id, ...], ...}


        keep both int index and str index; so that the ascending order is used in segmentations.

        Do str index at the end; need sorting while aligning with new samples.

        row_indices = [ (str(round(x['mz'], 4)), x['mz']) for x in list_samples[0].list_mass_tracks ] # index=row_indices,
        self.reference_mzdict = dict(row_indices)
        row_indices = [x[0] for x in row_indices]
        for pair in list_samples[0].anchor_mz_pairs:
            self.reference_anchor_pairs.append(( row_indices[pair[0]], row_indices[pair[1]] ))


        '''
        # initiation using the first Sample
        print("Initiating MassGrid, ...", list_samples[0].input_file)
        self.reference_anchor_pairs = list_samples[0].anchor_mz_pairs
        self._mz_landmarks_ = flatten_tuplelist(list_samples[0].anchor_mz_pairs)
        reference_mzlist = [ x['mz'] for x in list_samples[0].list_mass_tracks ]
        # setting up DataFrame for MassGrid
        self.MassGrid = pd.DataFrame(
            np.full((len(reference_mzlist), 1+self._number_of_samples_), None),
            columns=['mz'] + self.list_input_files,
        )
        # Add ref mz as a column to MassGrid; ref mzlist will be dynamic updated in MassGrid["mz"]
        self.MassGrid['mz'] = reference_mzlist
        self.MassGrid[ list_samples[0].input_file ] = [ x['id_number'] for x in list_samples[0].list_mass_tracks ]
        
        for SM in list_samples[1:]:
            self.add_sample(SM, database_cursor=None)
        
    def add_sample(self, sample, database_cursor=None):
        '''
        Add Sample instance to and update MassGrid. 
        To add: push each sample to SQLDB, - database_cursor;

        recalculate_ref is not done here, because it's easier to include unmatched features from Sample.
        If needed, the reference m/z values should be updated by revisiting DB samples.
        But the recalculation should be based on calibrated m/z values so that they are consistent across samples.

        push Sample to SQLDB

        '''
        print("Adding sample to MassGrid ", sample.input_file)
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


    def optimize_mass_grid(self):
        '''
        This inspects split or misaligned regions of m/z tracks.
        Problematic definition will persist in each sample's mass tracks,
        e.g. about 50 mass tracks are too close to each other in our plasma data; 200 in Zurich E. coli data.

        merge close mass tracks unless MassGrid suggests multiple features; 
        for split tracks in same samples, replace them by merged tracks and new id_numbers.

        to implement -

        '''


        pass



    def align_retention_time(self):
        '''
        Do alignment function using high-selectivity mass tracks.
        Step 1. get high-selectivity mass tracks among landmarks.
        2. for tracks of highest intensities, do quick peak detection to identify RT apexes.
        3. use RT from 2, do DWT conversion. Only masstracks with single peaks will be used for RT alignment,

        For low-selectivity mass tracks, could do 2-D deconvolution. Also see optimize_mass_grid

        '''

        





    def get_reference_RT_peaks(self):
        '''
        Use anchor mass trakcs, do quick peak detection. Use tracks of single peaks for RT alignment.
        '''
        pass


    def set_RT_reference(self):
        '''
        Because RT will not match precisely btw samples, 
        DTW should remap to a common set of time coordinates.

        Start with an initial sample of most peaks.

        Do a quick peak detection for good peaks;
                
        to avoid ambiguity in peak definitions.

        '''

        pass



    def match_ref_db(self):

        pass


    def global_peak_detection(self):
        pass


    def detect_peaks(self, mass_trace, min_intensity_threshold=10000, min_fwhm=3, min_prominence_threshold=5000, snr=2):
        '''

        to rewrite, applying to composite mass tracks



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

    def __update_mass_grid_by_init_samples__(self):
        pass






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

        print("Round 2 - numbers of epds and included peaks: ", (len(epds), len(found2)))
        # do de novo adduct initiations
        adduct_signatures = find_adduct_signatures(remaining_peaks, mztree, common_adducts[self.mode])
        for G in adduct_signatures:
            epds.append(G)

        print("Round 3 - numbers of epds: ", len(epds))
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

