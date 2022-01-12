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
import os

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

    Build Mass Grid first.
    mass correspondence is unique by using mass tracks. One grid can have maxium one mass track in one sample.
    MassGridDict = {id_mz_str: [m/z, {sample_id: massTrack id}], ...}


    1. Among the initial samples, the sample of most empCpds is used to seed the CMap.
    2. EmpCpds from each sample are aligned to the seed sample,
    then CMap is augmented by matched empCpds not in the seed sample. This step establishes m/z grid.
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
        composite mass traces, with study-wide peaks.

        masster_mz_registry has the pointers to individual samples
        '''
        
        self.experiment = experiment
        self.list_input_files = experiment.list_input_files

        self.MassGrid = None                    # will be DF
        self.FeatureGrid = None

        self.reference_mzdict = {}
        self.ref_empCpds = []
        

    def initiate_mass_grid(self, list_samples):
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

        MassGridDict = {id_mz_str: [m/z, {sample_id: massTrack id}], ...}
        '''
        # initiation using the first Sample
        row_indices = [ (str(round(x['mz'], 4)), x['mz']) for x in list_samples[0].list_mass_tracks ]
        self.reference_mzdict = dict(row_indices)
        row_indices = [x[0] for x in row_indices]
        self.reference_anchor_pairs = []             # (mz_str, mz_str), use mz str as identifiers 
        for pair in list_samples[0].anchor_mz_pairs:
            self.reference_anchor_pairs.append(( row_indices[pair[0]], row_indices[pair[1]] ))
        # DataFrame for MassGrid
        self.MassGrid = pd.DataFrame(
            np.zeros((len(row_indices), len(self.list_input_files)), dtype=int),
            index=row_indices,
            columns=self.list_input_files,
        )
        self.MassGrid[ self.list_input_files[0] ] = [ x['id_number'] for x in list_samples[0].list_mass_tracks ]

        print(self.MassGrid[:8])

        '''
        sample_id = list_samples[0].id
        for x in list_samples[0].list_mass_tracks:
            self.MassGridDict[ str(round(x['mz'], 4)) ] = [ x['mz'], { sample_id: x['id_number'], } ]

        for SM in list_samples[1:]:
            self.__update_mass_grid_by_init_sample__(SM)
        '''
        


    def __update_mass_grid_by_init_sample__(self, Sample):
        '''
        Update self.MassGrid, self.reference_anchor_pairs after adding another initiation sample.
        in MassGrid, mz is updated by taking mean values; 
        and mapping is updated by adding reference `sample_id: masstrack id`.

        Nested indices are: 
        mass_paired_mapping functions return positions of input lists -> 
        which refer to positions in anchor_pairs ->
        which refer to positions in list_mass_tracks or MassGridDict.
        '''
        # first align to reference_anchor_pairs
        ref_anchors_1 = [self.MassGridDict[x[0]][0] for x in self.reference_anchor_pairs]
        anchors_2 = [Sample.list_mass_tracks[pair[0]]['mz'] for pair in Sample.anchor_mz_pairs]
        mapped, correction_ = mass_paired_mapping_with_correction(
                        ref_anchors_1, anchors_2, std_ppm=5, correction_tolerance_ppm=1)
        updated_list2 = [x['mz']/(1+correction_) for x in Sample.list_mass_tracks]

        # mapped has () from lists 1 and 2; now check the paired values
        check_pattern_2 = [self.reference_anchor_pair[ii[0]][1] for ii in mapped]   # positions 
        check_pattern_2 = [self.MassGridDict[x][0] for x in check_pattern_2]        # m/z values to check, ordered as in mapped
        anchors_2_2 = [Sample.anchor_mz_pairs[ii[1]][1] for ii in mapped]        
        anchors_2_2 = [updated_list2[x] for x in anchors_2_2]
        mapped2, _ = mass_paired_mapping(check_pattern_2, anchors_2_2, std_ppm=5)

        # now tally pairs matched btw ref and samples - these two lists below are in right order of matches
        validated_ref_anchors = [self.reference_anchor_pairs[x] for x in [ mapped[ii[0]] for ii in mapped2 ]]
        validated_sample_anchors = [Sample.anchor_mz_pairs[x] for x in [ mapped[ii[1]] for ii in mapped2 ]]

        aligned_anchors = []
        for ii in range(len(validated_ref_anchors)):
            a,b = validated_ref_anchors[ii]
            c,d = validated_sample_anchors[ii]
            aligned_anchors += [(a, c), (b, d)]

        # now align everything between anchors






        LL = sorted([(len(SM.list_empCpds), SM) for SM in list_samples], reverse=True)



        init_dict = combine_mass_traces(LL[0][1].list_mass_traces)

        
        self.ref_mass_traces = []


        
        
        self.ref_empCpds = LL[0][1].list_empCpds
        
         
        
        self.master_masstrace_registry = {}




        
        empCpd_mzlist_1, empCpd_mzlist_2 = self.convert_empCpd_mzlist(LL[0][1]), self.convert_empCpd_mzlist(LL[1][1])
        aligned = init_cmap(empCpd_mzlist_1, empCpd_mzlist_2)

        for SM in LL[2:]:
            aligned.add_sample(SM[1])




        self.cmap = aligned





    def add_sample_to_map(self, sample, recalculate_ref=False):
        '''
        Use initial samples to start CMAP;
        push each sample to SQLDB, 
        Not keeping them in memory
        
        
        If recalculate_ref is True, update the ref mass trace m/z values.
        But this should not be used for regular samples after the initial ones.

        
        '''

        
        pass






    def set_RT_reference(self):
        '''
        Because RT will not match precisely btw samples, 
        DTW should remap to a common set of time coordinates.


        Do a quick peak detection for good peaks;
                but only masstracks with single peaks will be used for RT alignment,
        to avoid ambiguity in peak definitions.

        '''

        pass



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

