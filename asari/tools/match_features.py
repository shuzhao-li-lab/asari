from jms.search import *

def get_featureList(infile, start_row=1, mz_col=0, rt_col=1, sep='\t'):
    '''Read features from infile and return a list of json peaks.
    '''
    FL = []
    ii = 0
    for line in open(infile).read().splitlines()[start_row:]:
        ii += 1
        a = line.split(sep)
        FL.append(
            {'id': 'row'+str(ii+start_row), 'mz': float(a[mz_col]), 'rtime': float(a[rt_col])}
        )
    return FL

def list_match_lcms_features(list1, list2, mz_ppm=5, rt_tolerance=5):
    '''Match all features from list1 to list2 by mz_ppm and rt_tolerance. rt_tolerance unit follows the unit in lists.
    Using jms.search.build_centurion_tree(), jms.search.find_all_matches_centurion_indexed_list()
    input format list1 = [{'mz': 133.0970, 'rtime': 654.24, 'id': 'row555'}, ...]
    
    return dict_mapped
        dict_mapped: {id from list1: [id from list2, ...], ...}. No match is not included.
    '''
    dict_mapped = {}
    ctree = build_centurion_tree(list2)
    for p1 in list1:
        tmp = []
        mz_matched = find_all_matches_centurion_indexed_list(p1['mz'], ctree, limit_ppm=mz_ppm)
        for p2 in mz_matched:
            if abs(p1['rtime'] - p2['rtime']) < rt_tolerance:
                tmp.append( p2['id'] )
        if tmp:        
            dict_mapped[p1['id']] = tmp
            
    print("Of %d list1 features, number of uni-direction matched features is %d." %(len(list1), len(dict_mapped)))    
    return dict_mapped

def bidirectional_match(list1, list2, mz_ppm=5, rt_tolerance=5):
    '''
    There is 1:N matches in each direction.
    We resolve here 1:1 in both directions.
    '''
    dict1 = list_match_lcms_features(list1, list2, mz_ppm, rt_tolerance)
    dict2 = list_match_lcms_features(list2, list1, mz_ppm, rt_tolerance)
    
    print('    ~~~ match_numbers ~~~     \n')
    unique1 = [x for x,v in dict1.items() if len(v)==1]
    print("Unique Number of matched features in table 1: ", len(unique1))
    unique2 = [x for x,v in dict2.items() if len(v)==1]
    print("Unique Number of matched features in table 2: ", len(unique2))

    bi_unique = [x for x in unique1 if dict1[x][0] in unique2]
    print("Biodirectional, unique Number of matched feature pairs: ", len(bi_unique))
    
    return dict1, dict2


# Because 1:N matches occur often, we choose the best match btw two tables here
# The best match can be based on RT or mz 
#
def best_mz_match_lcms_features(list1, list2, mz_ppm=5, rt_tolerance=5):
    '''Match all features from list1 to list2 by mz_ppm and rt_tolerance, getting best match by minimal mz shift. 
    rt_tolerance unit follows the unit in lists.
    Using jms.search.build_centurion_tree(), jms.search.find_all_matches_centurion_indexed_list()
    input format list1 = [{'mz': 133.0970, 'rtime': 654.24, 'id': 'row555'}, ...]
    
    return dict_mapped
        dict_mapped: {id from list1: [id from list2, ...], ...}. No match is not included.
    '''
    dict_mapped = {}
    ctree = build_centurion_tree(list2)
    for p1 in list1:
        tmp = []
        mz_matched = find_all_matches_centurion_indexed_list(p1['mz'], ctree, limit_ppm=mz_ppm)
        for p2 in mz_matched:
            if abs(p1['rtime'] - p2['rtime']) < rt_tolerance:
                delta_mz = abs(p1['mz'] - p2['mz'])
                tmp.append( (delta_mz, p2['id']) )
        if tmp:
            if len(tmp) > 1:
                dict_mapped[p1['id']] = sorted(tmp)[0][1]
            else:
                dict_mapped[p1['id']] = tmp[0][1]
            
    print("Of %d list1 features, number of uni-direction matched features is %d." %(len(list1), len(dict_mapped)))    
    return dict_mapped


def best_rt_match_lcms_features(list1, list2, mz_ppm=5, rt_tolerance=5):
    '''Match all features from list1 to list2 by mz_ppm and rt_tolerance, getting best match by minimal mz shift. 
    rt_tolerance unit follows the unit in lists.
    Using jms.search.build_centurion_tree(), jms.search.find_all_matches_centurion_indexed_list()
    input format list1 = [{'mz': 133.0970, 'rtime': 654.24, 'id': 'row555'}, ...]
    
    return dict_mapped
        dict_mapped: {id from list1: [id from list2, ...], ...}. No match is not included.
    '''
    dict_mapped = {}
    ctree = build_centurion_tree(list2)
    for p1 in list1:
        tmp = []
        mz_matched = find_all_matches_centurion_indexed_list(p1['mz'], ctree, limit_ppm=mz_ppm)
        for p2 in mz_matched:
            delta_rt = abs(p1['rtime'] - p2['rtime'])
            if delta_rt < rt_tolerance:
                tmp.append( (delta_rt, p2['id']) )
        if tmp:
            if len(tmp) > 1:
                dict_mapped[p1['id']] = sorted(tmp)[0][1]
            else:
                dict_mapped[p1['id']] = tmp[0][1]
            
    print("Of %d list1 features, number of uni-direction matched features is %d." %(len(list1), len(dict_mapped)))    
    return dict_mapped


def bidirectional_best_match(list1, list2, mz_ppm=5, rt_tolerance=5):
    '''
    There is 1:N matches in each direction.
    We resolve here 1:1 in both directions, by two methods considering best RT or m/z matches.
    '''
    print('\n    ~~~ By best rtime matches ~~~     \n')
    dict1 = best_rt_match_lcms_features(list1, list2, mz_ppm, rt_tolerance)
    dict2 = best_rt_match_lcms_features(list2, list1, mz_ppm, rt_tolerance)
    
    bi_unique = [(k,v) for k,v in dict1.items()] + [(v,k) for k,v in dict2.items()]
    print("~~~ Biodirectional, unique Number of matched feature pairs: ~~~\n", len(bi_unique)-len(set(bi_unique)))

    print('\n\n########################################################################')
    print('    ~~~ By best m/z matches ~~~     \n')
    dict1 = best_mz_match_lcms_features(list1, list2, mz_ppm, rt_tolerance)
    dict2 = best_mz_match_lcms_features(list2, list1, mz_ppm, rt_tolerance)

    pairs2 = [(v,k) for k,v in dict2.items()]
    valid_matches = [pair for pair in [(k,v) for k,v in dict1.items()] if pair in pairs2]
    print("~~~ Biodirectional, unique Number of matched feature pairs: ~~~\n    ", len(valid_matches))

    print('########################################################################\n\n')

    return valid_matches, dict1, dict2


def convert_min2secs(LL):
    for p in LL:
        p['rtime'] = p['rtime']*60
    return LL

def convert_sec2mins(LL):
    for p in LL:
        p['rtime'] = p['rtime']/60
    return LL



#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    # xcms_ = get_featureList('YeastNeg2021_NetID_XCMS_featureTable.csv', start_row=1, mz_col=1, rt_col=4, sep=',')
    

    from sys import argv

    file1 = argv[1]
    list1 = get_featureList(file1, start_row=1, mz_col=1, rt_col=2, sep='\t')

    # true_ = get_featureList('manual_certified.txt', start_row=1, mz_col=1, rt_col=2, sep='\t')
    file2 = argv[2]
    list2 = get_featureList(file2, start_row=1, mz_col=1, rt_col=2, sep='\t')
    list2 = convert_min2secs(list2)

    print("\n  List based inclusive comparisons:")
    dict1, dict2 = bidirectional_match(list1, list2, mz_ppm=5, rt_tolerance=6)
    # dict1, dict2 = bidirectional_match(list2, list1, mz_ppm=5, rt_tolerance=6)

    print("\n  Best match comparisons:")
    valid_matches, dict1, dict2 = bidirectional_best_match(list1, list2, mz_ppm=5, rt_tolerance=6)
    # valid_matches, dict1, dict2 = bidirectional_best_match(list2, list1, mz_ppm=5, rt_tolerance=6)

    print("Unmatched features: ")
    for p in list2:
        if p['id'] not in [x[1] for x in valid_matches]:
            print(p)

