'''
Two lines of search methods:
- A centurion_tree is an indexed dictionary of peaks/features.
- We also use DataFrame based vector operations.

Emperical compounds are constructed by co-eluting isotopic and adduct patterns.
# isotopic_signatures: example as [(182, 'anchor'), (191, 'M(13C)'), (205, 'M(18O)')]
# consider Cl as adduct, both 35Cl and 37Cl are considered
# 34S and 37Cl are close, but should have diff ratios and are formula dependent
'''

# Mass difference and Abundance of natually occuring isotopic elements.
# (mz difference, notion, ratio low limit, ratio high limit)
# Ratio is loose to leave room for 1) multiple atoms, e.g. 40 C atoms lead to ~ 40*1% in abundance
# 2) measurement errors; not all peaks are in linear dynamic range
# 3) quantification by peak height is not same as peak area

isotopic_patterns = [
    # mass diff, isotopes, (intensity ratio constraint)
    (1.003355, '13C/12C', (0, 0.8)),      # 13C-12C, 12C~99%, 13C ~ 1%
    (0.997035, '15N/14N', (0, 0.2)),     # 15N-14N, 14N ~ 99.64%, 15N ~ 0.36%
    (2.004245, '18O/16O', (0, 0.2)),      # 18O-16O, 16O ~ 99.76, 16O ~ 0.2%
    (1.995796, '34S/32S', (0, 0.4)),      # 32S (95.02%), 33S (0.75%), 34S (4.21%)
    (0.999388, '33S/32S', (0, 0.1)),
    # double isotopes
    (2.00039, 'M(13C),M(15N)', (0, 0.2)),
    (2.999151, 'M(13C),M(34S)', (0, 0.4)),
    # double charged
    (0.5017, '13C/12C, double charged', (0, 0.8)),
    (0.4985, '15N/14N, double charged', (0, 0.2)),
]

common_adducts = {
    # mass diff, modification
    # not using (intensity ratio constraint), but it can be documented or learned
    'pos': [
        (1.0078, 'H'),
        (21.9820, 'Na/H'), # Na replacing H
        (10.991, 'Na/H, double charged'),
        (18.0106, '+H2O'), 
        (18.033823, '+NH4'),
        (37.9559, '39K/H'),
        (39.9540, '41K/H'),
        (41.026549, 'Acetonitrile'),
    ],
    'neg': [
        (1.0078, 'H'),
        (22.9893, 'Na'),
        (20.97474706646, '+Na-2H'),
        (18.0106, 'H2O'), 
        (34.9689, '35Cl'),
        (36.9659, '37Cl'),
        (40.01926853323, '+ACN-H'),
        (44.998201, 'COOH'),
        (59.013295, 'CH3COO'),
    ],
}

extended_adducts = {
    'pos': [],
    'neg': [],
}



#
# -----------------------------------------------------------------------------
#

def build_centurion_tree(list_peaks):
    '''
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...]
    Return a dictionary, indexing mzList by 100*mz bins.
    Because most high-resolution mass spectrometers measure well under 0.01 amu, 
    one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).
    list_mass_tracks has similar format as list_peaks.
    '''
    d = {}
    for p in list_peaks:
        cent = int(100 * p['mz'])
        if cent in d:
            d[cent].append(p)
        else:
            d[cent] = [p]
    return d


def build_peak_id_dict(list_peaks):
    d = {}
    for p in list_peaks:
        d[p['id_number']] = p
    return d


def __build_centurion_tree_mzlist(mzList):
    '''
    Return a dictionary, indexing mzList by 100*mz bins.
    Because most high-resolution mass spectrometers measure well under 0.01 amu, 
    one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).
    '''
    d = {}
    for ii in range(len(mzList)):
        cent = int(100*mzList[ii])
        if cent in d:
            d[cent].append((mzList[ii], ii))
        else:
            d[cent] = [(mzList[ii], ii)]
    return d


def find_all_matches_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=5):
    '''
    Return matched peaks in mz_centurion_tree.
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    results = []
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            if abs(peak['mz']-query_mz) < mz_tol:
                results.append(peak)
                
    return results


def find_best_match_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=2):
    '''
    Return matched indices in mz_centurion_tree (based on peak list).
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    result = (None, 999)
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            _d = abs(peak['mz']-query_mz)
            if _d < min(result[1], mz_tol):     # enforce mz_tol here
                result = (peak, _d)
                
    return result[0]


def is_coeluted(P1, P2):
    '''
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 'left_base': 648, 'right_base': 655, 'id_number': 555},
                 ...]
    coelution: Overlap more than half of the smaller peak.
    '''
    _coeluted = False
    len1, len2 = P1['right_base'] - P1['left_base'], P2['right_base'] - P2['left_base']
    # overlap is the max L and min R
    overlap = min(P1['right_base'], P2['right_base']) - max(P1['left_base'], P2['left_base'])
    if overlap > 0.5 * min(len1, len2):
        _coeluted = True
    return _coeluted

   
def find_isotopic_signatures(list_peaks, mztree, isotopic_patterns, mz_tolerance_ppm=5, rt_tolerance_scans=5):
    '''
    See find_isotopic_pairs. This extends to all related isotopic signatures.
    Allows ambiguous matches, which are dealt with in empCpd statistics.

    Input
    =====
    isotopic_patterns = [(1.003355, '13C/12C', (0, 0.2)), ...], the third item is optional limits of abundance ratio.

    Return
    ======
    list of lists of peak numbers that match isotopic patterns, [ [], ... ]

    Example
    =======
    [ [(195, 'anchor'), (206, '13C/12C')], 
      [(182, 'anchor'), (191, '13C/12C'), (205, '18O/16O')],
      [(295, 'anchor'), (335, '13C/12C'), (368, 'M(13C),M(34S)')], ...]
    '''

    signatures = []
    for P1 in list_peaks:
        matched = [ (P1['id_number'], 'anchor'), ]          # Nature is nice to have lowest mass for the most abundant 
        for isotopic_pair in isotopic_patterns:
            (mass_difference, relation) = isotopic_pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if abs(P1['apex']-P2['apex']) <= rt_tolerance_scans and is_coeluted(P1, P2):
                    if len(isotopic_pair) == 2:             # not checking abundance ratio
                        matched.append( (P2['id_number'], relation) )
                    else:
                        (abundance_ratio_min, abundance_ratio_max) = isotopic_pair[2]
                        # checking abundance ratio
                        if abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                            matched.append( (P2['id_number'], relation) )
        if len(matched) > 1:
            signatures.append(matched)

    return signatures



def find_adduct_signatures(list_peaks, mztree, adduct_patterns, mz_tolerance_ppm=5):
    '''
    Search adduct mass_diff in ceelution peaks, de novo. 
    Not requiring matched apex as in isotopic pairs, nor abundance ratio restriction.
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [ (P1['id_number'], 'anchor'), ]
        for adduct in adduct_patterns:
            tmp = find_all_matches_centurion_indexed_list( P1['mz'] + adduct[0], mztree, mz_tolerance_ppm )
            for P2 in tmp:
                if is_coeluted(P1, P2):
                    matched.append( (P2['id_number'], adduct[1]) )

        if len(matched) > 1:
            signatures.append(matched)

    return signatures


def find_mzdiff_pairs_from_masstracks(list_mass_tracks, list_mz_diff=[1.003355, 21.9820], mz_tolerance_ppm=5):
    '''
    Find all pairs in list_mass_tracks that match a pattern in list_mz_diff, and return their id_numbers as pairs.
    This function does not use coeluction (rtime) rules. 

    Input
    =====
    list_mass_tracks: [{ 'id_number': ii,  'mz': xic[0], 'rt_scan_numbers': xic[1],  'intensity': xic[2],  }, ...]
    list_mz_diff: defaul 1.003355, 21.9820 are 13C/12C, Na/H. Dependent on ionization mode.
    
    Return
    ======
    list of pairs of mass tracks numbers.
    '''
    pairs = []
    # list_mass_tracks has similar format as list_peaks.
    mztree = build_centurion_tree(list_mass_tracks)
    for mzdiff in list_mz_diff:
        for P1 in list_mass_tracks:
            P2 = find_best_match_centurion_indexed_list(P1['mz'] + mzdiff, mztree, mz_tolerance_ppm)
            if P2:
                pairs.append((P1['id_number'], P2['id_number']))

    return pairs






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









#
# -----------------------------------------------------------------------------
#
# search using indexed trees, refactoring code from mass2chem
# 


def search_formula_mass_db(query_mz, indexed_DB, limit_ppm=10):
    '''
    Find best matched formula_mass in indexed_DB for query_mz within ppm limit.
    
    indexed_DB: # [formula_mass, m/z, charged_formula, selectivity, [relations to neutral formula]]
            example DB_1[99] = 
            [['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]], 
            ['C4H5NS_99.013721', 99.01372099999999, 'C4H5NS', 0.9967247109333564, [('C4H5NS_99.01427', 'M[1+]')]], 
            ['C3H8NaO2_99.041706', 99.04170646677, 'C3H8NaO2', 0.999999999998776, [('C3H8O2_76.05243', 'M+Na[1+]')]], ...]

    return
    ------
    Closest match as an entry,
             (match_ppm, [formula_mass, m/z, charged_formula, selectivity, [relations to neutral formula]])
    None if out of ppm limit.
    '''
    _delta = query_mz * limit_ppm * 0.000001
    _low_lim, _high_lim = query_mz - _delta, query_mz + _delta
    result = []
    for ii in range(int(_low_lim), int(_high_lim)+1):
        if ii in indexed_DB:
            for F in indexed_DB[ii]:
                # F[1] is m/z, fixed DB format
                result.append( (abs(query_mz-F[1]), F) )

    result.sort()
    if result and result[0][0] < _delta:
        return result[0][1]
    else:
        return None



def annotate_formula_mass(list_query_mz, indexed_DB, check_mass_accuracy=False, limit_ppm=10):
    '''
    Annotate input m/z list using formula based mass in the indexed_DB,
    which should be generic, not assuming separation of primary ions and their common isotopes and adducts.
    Search efficiency will be improved by the construction strategy of the indexed_DB,
    not dependent on the algorithm in this function. E.g., more frequent features should be in first DB to search.

    If check_mass_accuracy is True, the ppm deviations of each matched primary ion are fitted to a normal distribution,
    to estimate mass_accuracy and ppm_std.

    input
    -----
    input_list: m/z values.
    indexed_DB: per entry - [formula_mass, m/z, charged_formula, ion, selectivity, [relations to neutral formula]]
                This can be a pre-loaded DB, or hot DB in asari.

    # mode='pos', 

    return
    ------
    {mass_accuracy: estimated based on difference from formula value
    ppm_std:         estimated based on precision distribution
    annotated_list: [matched DB entry or None, ...]
    corrected_list_query_mz: m/z values will be corrected if they deviate from formula values by mass_accuracy > 5 ppm.
    }
    '''
    result = list_search_formula_mass_db(list_query_mz, indexed_DB, limit_ppm)
    mass_accuracy, ppm_std, corrected_list_query_mz = None, None, None
    N = len(list_query_mz)
    if check_mass_accuracy:
        list_ppm_errors = []
        for ii in range(N):
            if result[ii]:
                list_ppm_errors.append( 1000000 * (list_query_mz[ii] - result[ii][1])/result[ii][1] )

        mass_accuracy, ppm_std = normal_distribution.fit(list_ppm_errors)
        if abs(mass_accuracy) > 5:   # this is considered significant mass shift, requiring m/z correction for all 
            list_query_mz = [x-x*0.000001*mass_accuracy for x in list_query_mz]
            # redo search because we may find different matches after mass correction
            # Also change search ppm if self.__mass_stdev__ > 5
            result = list_search_formula_mass_db(list_query_mz, indexed_DB, max(10, 2*ppm_std))
            corrected_list_query_mz = list_query_mz
            # update ppm_std because it may be used for downstream search parameters
            list_ppm_errors = []
            for ii in range(N):
                if result[ii]:
                    list_ppm_errors.append( 1000000 * (list_query_mz[ii] - result[ii][1])/result[ii][1] )
            _mu, ppm_std= normal_distribution.fit(list_ppm_errors)

        elif ppm_std > 5:      # no mass correction but need redo search using higher ppm for unmatched mz
            result = list_search_formula_mass_db(list_query_mz, indexed_DB, 2*ppm_std)

    return {
        'mass_accuracy': mass_accuracy,
        'ppm_std': ppm_std,
        'annotated_list': result,
        'corrected_list_query_mz': corrected_list_query_mz,
        }






#
# -----------------------------------------------------------------------------
#
# not used now
# 


def find_isotopic_pairs(list_peaks, mztree, isotopic_pair, mz_tolerance_ppm=5, rt_tolerance_scans=5):
    '''
    - to update

    To find naturally occuring isotopic pairs in LC-MS peaks, 
    using strict coelution criteria and loose abundance ratio.
    This is done within a sample, so that all elution time is precisely matched by scan numbers.
    This avoids any issues in alignment across samples.
    Rules: elution apex should be within 5 scans; overlap of elution peaks minimal 50% of the smaller peak;
            abundance ratio loosely applied, leaving room for measuring errors.
    
    Abundance of natually occuring isotopic elements:
    12C~99%, 13C ~ 1%
    14N ~ 99.64%, 15N ~ 0.36%
    Sulfur:  32S (95.02%), 33S (0.75%), 34S (4.21%), 
    Chlorine has two stable isotopes, 35Cl (75.77%) and 37Cl (24.23%).
    1H and 31P are basically ~100% in nature.

    Input
    =====
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
    'left_base': 648, 'right_base': 655, 'id_number': 555},
                 ...]
    mztree is constructed on list_peaks.
    isotopic_pair: (mass_difference, notion, abundance_ratio). E.g. (1.0034, 'M(C13)', 0.2). - to update
    
    Return
    ======
    list of pairs of peaks.
    '''
    pairs = []
    (mass_difference, relation, abundance_ratio_min, abundance_ratio_max) = isotopic_pair
    for P1 in list_peaks:
        tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
        for P2 in tmp:
            if abs(P1['apex']-P2['apex']) <= rt_tolerance_scans and is_coeluted(P1, P2) \
                                             and abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                pairs.append((P1, P2))

    return pairs


def init_eTrees_by_isotopic_signatures(isotopic_signatures, iso_type, peak_list):
    '''
    An example of isotopic_signatures: 
    [(296, 'anchor'), (339, 'M(13C)'), (332, 'M(15N)'), (333, 'M(15N)')], using peak numbers. 
    

    '''
    eTrees = []
    (P1, P2) = iso_pairs
    for pp in iso_pairs:
        ET = epdTree(peak_list)
        ET.anchor_peak = P1['id_number']
        ET.isotope_peaks = [P1['id_number'], P2['id_number']]
        ET.peak_relationships.append(( P1['id_number'], P2['id_number'], iso_type ))
    return eTrees

def merge_isotopic_eTrees(list_eTrees):
    '''
    For each pair of eTrees after init_eTrees_by_isopairs, merge if at least one overlap peak is anchor. 
    Ignor if overlap is all at endpoints. I.e. trees can branch, but not converge.
    '''
    new = []
    peak2etree = {}
    for ET in list_eTrees:
        for P in ET.isotope_peaks:
            if P in peak2etree:
                peak2etree[P].append()

    pass


