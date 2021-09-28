'''
Utilities for merging feature tables from same study.

# F1 = (m/z, rt)
def match2(F1, F2):
    if abs(F1[0]-F2[0])/F1[0] < PPM_tolerance and abs(F1[1] - F2[1]) < RTime_tolerance:
        return True
    else:
        return False

def __bin_by_median_rt__(List_of_peaks, tolerance):
    List_of_tuples = [(P.cal_rtime, P) for P in List_of_peaks]
    List_of_tuples.sort()
    return bin_by_median(List_of_tuples, lambda x: max(tolerance, 0.1*x))

def _formula_selectivity_(L):
    d = {}
    for x in L: d[x] = L.count(x)
    return d

'''


import numpy as np

# copied from algorithms
def bin_by_median(List_of_tuples, func_tolerance):
    '''
    Not perfect because left side may deviate out of tolerance, but LC-MS data always have enough gaps for separation.
    Will add kernel density method for grouping m/z features.
    List_of_tuples: [(value, object), (value, object), ...], to be separated into bins by values (either rt or mz).
                    objects have attribute of sample_name if to align elsewhere.
    return: [seprated bins], each as a list of objects as [X[1] for X in L]. Possible all falls in same bin.
    '''
    new = [[List_of_tuples[0], ], ]
    for X in List_of_tuples[1:]:
        if X[0]-np.median([ii[0] for ii in new[-1]]) < func_tolerance(X[0]):       # median moving with list change
            new[-1].append(X)
        else:
            new.append([X])
    PL = []
    for L in new:
        PL.append([X[1] for X in L])
    return PL


def _read_asari_ftables_(infile, make_feature_id=False):
    '''
    return two dictionaries, {key: [formula_mass, mz, rtime]}, and {key: line}. Make key if needed.
    '''
    key_dict, data_dict = {}, {}
    lines = open(infile).read().splitlines()
    data_dict['header'] = lines[0]
    for line in lines[1:]:
        a = line.split('\t')
        if make_feature_id:     # old format, starting cols are [formula_mass, mz, rtime]
            formula_mass, mz, rtime = a[:3]
            feature_id = mz + '@' + rtime
        else:                   # new format, starting cols are [feature_id, formula_mass, mz, rtime]
            feature_id, formula_mass, mz, rtime = a[:4]
        key_dict[feature_id] = [formula_mass, float(mz), float(rtime)]
        data_dict[feature_id] = line.rstrip()
    return key_dict, data_dict


def _masstrace2features_(LL, RTime_tolerance):
    '''
    input on same mass [[], [], [], ...], one list of features per input table, in fids like `148.1105@228.15`
    '''
    def find_min_delta(L):
        ii = 0
        diffs = []
        for x in sorted(L):
            diffs.append(x-ii)
            ii = x
        return min(diffs[1:])

    NN = len(LL)
    RT_LL = []
    for ii in range(NN):
        for x in LL[ii]:
            # print(ii, x)
            RT_LL.append( (float(x.split("@")[1]), (x,ii)) )
            
    # determine RT tolerance by comparing RTime_tolerance to min distances within one feature table
    min_multiplets = [find_min_delta([float(x.split("@")[1]) for x in L]) for L in LL if len(L) > 1]
    if min_multiplets:
        RTime_tolerance = min(RTime_tolerance, min(min_multiplets)-1)           # tolerance should not merge features in the same original table

    RT_LL.sort()
    list_features = bin_by_median(RT_LL, lambda x: RTime_tolerance)
    # now fit grouped features into a grid corresponding to original input tables
    grouped = []
    for L in list_features:
        tmp = [ None ] * NN
        for x,ii in L:
            tmp[ii] = x
        grouped.append(tmp)

    return grouped


def merge_feature_tables(list_of_feature_tables, make_feature_id=False, PPM_tolerance = 4, RTime_tolerance = 15, outfile="merged_featuretable.tsv"):
    '''
    To merge asari tables from the same study on same platform. I.e. not expecting major differences in m/z and RT.
    Step 1. Identify massTraces
    Step 2. Determine features in each massTrace
    Step 3. out to outfile
    '''

    #short_names = [os.path.basename(f) for f in list_of_feature_tables]

    NN = len(list_of_feature_tables)
    
    formula_mass_dict, M_pointer = {}, {}
    all_keys, all_tables, unassigned = [], [], []
    for f in list_of_feature_tables:
        key_dict, data_dict = _read_asari_ftables_(f, make_feature_id)          # {key: [formula_mass, mz, rtime]}, and {key: line}
        all_keys.append(key_dict)
        all_tables.append(data_dict)
        for k,v in key_dict.items():
            if '_M_' == v[0][:3]:
                unassigned.append((v[1], v[2], v[0]))                    # (mz, rtime, fid)
            else:
                formula_mass_dict[v[0]] = [[]]*NN

    # Identify massTraces. formula_mass can do; deal with features without formula_mass
    unassigned.sort()
    unassigned = [(x[0], x[2]) for x in unassigned]
    for bin in bin_by_median(unassigned, lambda x: PPM_tolerance * 0.000001 * x):
        for x in bin:
            M_pointer[x] = bin[0]
    # now all formula_mass unified across samples, either by formula_mass or _M_ identifier mapping in M_pointer
    for x in M_pointer.values():
        formula_mass_dict[x] = [[]]*NN

    # index features to mass in all tables
    List_dict_mass2features = []
    for kdict in all_keys:                                      # {key: [formula_mass, mz, rtime]}
        d = {}
        for k,v in kdict.items():
            if v[0] in d:
                d[v[0]].append(k)
            else:
                d[v[0]] = [k]
        List_dict_mass2features.append(d)
    # organize each table into formula_mass_dict, {formula_mass: [[fid, fid, ...], [], [], ...] }
    for ii in range(NN):
        for k,v in List_dict_mass2features[ii].items():                 # {formula_mass: list of feature_ids}
            if '_M_' == k[:3]:
                k = M_pointer[k     ]
            formula_mass_dict[k][ii] = v

    new_features = []
    for k,v in formula_mass_dict.items():
        new_features.append( (k, _masstrace2features_(v , RTime_tolerance)) )

    numbers_cols = [ len(data_dict['header'].rstrip().split('\t')) for data_dict in all_tables ]

    s = 'new_feature_id\tnumber_src_tables\told_ids\tformula_mass\t' + '\t'.join([data_dict['header'].rstrip() for data_dict in all_tables]) + '\n'
    for k,LL in new_features:       # formula_mass, list of features per input table
        for L in LL:
            new_id = [x for x in L if x][0]
            filled_id = [x or '_' for x in L]
            number_src_tables = str(len([x for x in L if x]))
            line = '\t'.join([new_id, number_src_tables, ','.join(filled_id), k])
            for ii in range(NN):
                if L[ii]:           # None if not found in a position
                    line += '\t' + all_tables[ii][L[ii]]
                else:
                    line += '\t'*numbers_cols[ii]               # '\t' + 
            s += line + '\n'

    with open(outfile, 'w') as O:
        O.write(s)
    print("Merged result was written to ", outfile)
    
#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':

    list_of_feature_tables = [
        'feature_table_b1.tsv',  'feature_table_b2.tsv', 'feature_table_b3.tsv', 'feature_tableb4.tsv',
    ]

    # make_feature_id=False if feature_id in first col
    merge_feature_tables(list_of_feature_tables, make_feature_id=False, 
                        PPM_tolerance = 4, RTime_tolerance = 15, 
                        outfile="merged_featuretable_0928.tsv")

