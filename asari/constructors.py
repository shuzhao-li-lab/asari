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

from .search import *


class FeatureMap:

    def __init__(self, experiment):

        self.experiment = experiment

        self.init_index = []

        self.db_ref = []

        self.featuremap = []

    def construct_mz_map(self):
        pass

    def set_RT_reference(self):
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

        Result
        ======
        This updates self.epds, [{'id': ii, 'list_peaks': [(peak_id, ion), (), ...],}, ...]
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

