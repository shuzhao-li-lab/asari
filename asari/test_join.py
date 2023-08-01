import os 

import pandas as pd
import numpy as np
import numpy as np
import statsmodels.api as sm
import multiprocessing as mp

def find_anchors(t, similar_chromatography=False, rt_similarity=0.05, mz_ppm=5):
    anchors = []
    t['mz_min'] = t['mz'] - t['mz'] / 1e6 * mz_ppm
    t['mz_max'] = t['mz'] + t['mz'] / 1e6 * mz_ppm
    t['rt_min'] = t['rtime'] - rt_similarity * t['rtime']
    t['rt_max'] = t['rtime'] + rt_similarity * t['rtime']
    for mz_min, mz_max, rt_min, rt_max in zip(t['mz_min'], t['mz_max'], t['rt_min'], t['rt_max']):
        t_sel = t[t['mz'].between(mz_min, mz_max)]
        if similar_chromatography:
            t_sel = t_sel[t_sel['rtime'].between(rt_min, rt_max)]
        anchors.append(t_sel.shape[0]==1)
    t['anchor'] = anchors
    return t

def join_ftables(ftables, similar_chromatography=False, rt_similarity=0.05, mz_ppm=5):
    workers = mp.Pool(mp.cpu_count())
    mp.freeze_support()

    iter = 0
    index = len(ftables)
    ftable_objects = []
    for ft in ftables:
        try:
            ftable_objects.append(pd.read_csv(ft, sep="\t"))
        except:
            pass

    ftables = {i: t for i, t in enumerate(workers.starmap(find_anchors, [(x, similar_chromatography, rt_similarity, mz_ppm) for x in ftable_objects]))}
    fail_list = []
    scores = {}
    while len(list(ftables.keys())) > 1:
        iter += 1
        list_keys = list(ftables.keys())
        to_compute = []
        for i, k1 in enumerate(list_keys):
            for k2 in list_keys[i + 1:]:
                if (k1, k2) not in scores:
                    to_compute.append((k1, k2))

        new_scores = workers.starmap(map_anchors, [(ftables[x[0]], ftables[x[1]], similar_chromatography) for x in to_compute])  
        for x, anchor_mapping in zip(to_compute, new_scores):
            scores[(x[0],x[1])] = (len(anchor_mapping)/min(ftables[x[0]].shape[0], ftables[x[1]].shape[0]), anchor_mapping)

        best_pair = [0, None, None, None]
        for (key1, key2), (score, mapping) in scores.items():
            if score > best_pair[0] and (key1, key2) not in fail_list:
                best_pair = [score, mapping, key1, key2]
        if best_pair[1] is None:
            for id, table in ftables.items():
                table.sort_values(by=['mz', 'rtime'], inplace=True)
                table.to_csv("joined_table_" + str(id) + ".tsv", sep="\t")
            return ftables

        joined_table = join_table_pair(ftables[best_pair[2]], ftables[best_pair[3]], best_pair[1])
        del scores[(best_pair[2], best_pair[3])]
        if joined_table[0]:
            ftables[index + iter] = joined_table[1]
            to_del = []
            for k1, k2 in scores.keys():
                if k1 in (best_pair[2], best_pair[3]) or k2 in (best_pair[2], best_pair[3]):
                    to_del.append((k1, k2))
            for td in to_del:
                del scores[td]
            del ftables[best_pair[2]]
            del ftables[best_pair[3]]
        else:
            fail_list.append((best_pair[1], best_pair[2]))
    print("OUT")
    for id, table in ftables.items():
        table.sort_values(by=['mz', 'rtime'], inplace=True)
        table.to_csv("joined_table_" + str(id) + ".tsv", sep="\t")
    list(ftables.values())[0].to_csv("final_paired.tsv", sep="\t")

def anchor_mapping_wrapper(input):
    return map_anchors(input['t1'], input['t2'], input['similar_chromatography'])

def map_anchors(t1, t2, similar_chromatography):
    t1_x = t1[t1["anchor"] == True]
    t2_x = t2[t2["anchor"] == True]
    common_ion_rts = []
    for mz1, mz1_min, mz1_max, rt1, rt1_min, rt1_max in zip(t1_x['mz'], t1_x['mz_min'], t1_x['mz_max'], t1_x['rtime'], t1_x['rt_min'], t1_x['rt_max']):
        if similar_chromatography:
            t2_x_sel = t2_x[(t2_x['rtime'].between(rt1_min, rt1_max)) & (t2_x['mz'].between(mz1_min, mz1_max))]
        else:
            t2_x_sel = t2_x[t2_x['mz'].between(mz1_min, mz1_max)]
        if t2_x_sel.shape[0] == 1:
            common_ion_rts.append((rt1, t2_x_sel['rtime'].values[0], mz1, t2_x_sel['mz'].values[0]))
    return np.array(common_ion_rts)

def find_best_frac_lowess(t1, anchor_mapping, it=3):
    lowess = sm.nonparametric.lowess
    frac = 1.0
    last_frac = None
    last_error = np.inf
    while frac > 0.2:
        found_solution = False
        frac -= 0.05
        predicted = lowess(anchor_mapping[:,1], anchor_mapping[:,0], it=it, xvals=anchor_mapping[:,0], frac=frac)
        new_error = np.sum(anchor_mapping[:,1] - lowess(anchor_mapping[:,1], anchor_mapping[:,0], it=it, xvals=anchor_mapping[:,0], frac=frac))
        if abs(new_error) < last_error and np.count_nonzero(np.isnan(predicted)) == 0:
            full_range = lowess(anchor_mapping[:,1], anchor_mapping[:,0], it=it, xvals=range(int(min((t1['rtime']))), int(max(t1['rtime']))), frac=frac)
            if np.count_nonzero(np.isnan(full_range)) == 0:
                found_solution = True
                last_frac = frac
        if not found_solution:
            break
    return last_frac

def join_table_pair(t1, t2, anchor_mapping):
    def __extract_row(row, columns):
        return {column: row[column] for column in columns}
    best_frac = find_best_frac_lowess(t1, anchor_mapping)
    lowess = sm.nonparametric.lowess

    if best_frac is None:
        return (False, t1, t2)
    else:
        t1['pred_t2_rtime'] = lowess(anchor_mapping[:,1], anchor_mapping[:,0], it=3, xvals=t1['rtime'], frac=best_frac)

    skip_columns = ['parent_masstrack_id', 'peak_area', 'cSelectivity', 'goodness_fitting', 'snr', 'detection_counts', 'pred_t2_rtime', 'anchor']
    common_columns = set(list(t1.columns)).intersection(set(list(t2.columns)))
    t1_unique_columns = set([x for x in t1.columns if x not in common_columns])
    t2_unique_columns = set([x for x in t2.columns if x not in common_columns])

    used_t2 = set()
    new_table = []

    t1['rtime_left_base'] = round(t1["rtime_left_base"] - (t1['rtime'] - t1['pred_t2_rtime']), 2)
    t1['rtime_right_base'] = round(t1['rtime_right_base'] - (t1['rtime'] - t1['pred_t2_rtime']), 2)
    t1['rtime'] = round(t1['pred_t2_rtime'], 2)
    t1.drop(labels='pred_t2_rtime', inplace=True, axis=1)

    t1_rows = t1.apply(__extract_row, axis=1, args=([c for c in t1.columns if c not in skip_columns],))
    t2_rows = t2.apply(__extract_row, axis=1, args=([c for c in t2.columns if c not in skip_columns],))
    t2_row_map = {x['id_number']: x for x in t2_rows}

    fluff_1 = {c: 0 for c in t2_unique_columns if c not in skip_columns}
    fluff_2 = {c: 0 for c in t1_unique_columns if c not in skip_columns}

    x,y,z = 0,0,0
    for row_1 in t1_rows:
        best_match = (np.inf, None)
        t2_sel = t2[~t2['id_number'].isin(used_t2)]
        t2_sel = t2_sel[t2_sel['mz'].between(row_1['mz_min'], row_1['mz_max'])]
        t2_sel = t2_sel[(t2_sel['rtime_left_base'] < row_1['rtime']) & (t2_sel['rtime_right_base'] > row_1['rtime'])]
        for row_2 in [t2_row_map[t2_id] for t2_id in t2_sel['id_number']]:
            score = (row_2['mz'] - row_1['mz'])**2 + (row_1['rtime'] - row_2['rtime'])**2
            if score < best_match[0]:
                best_match = (score, row_2)
        if best_match[1]:
            x += 1
            row_2 = best_match[1]
            used_t2.add(row_2['id_number'])
            new_row = dict(row_1)
            new_row['id_number'] = 'F' + str(x + y + z)
            new_row.update(row_2)
            new_row['rtime_left_base'] = min(row_1['rtime_left_base'], row_2['rtime_left_base'])
            new_row['rtime_right_base'] = max(row_1['rtime_right_base'], row_2['rtime_right_base'])
            new_row['rtime'] = (row_1['rtime'] + row_2['rtime'])/2
            new_table.append(new_row)
        else:
            y += 1
            row_1.update(fluff_1)
            row_1['id_number'] = 'F' + str(x + y + z)
            new_table.append(row_1)
            y += 1
    for row_2 in t2_rows:
        if row_2['id_number'] not in used_t2:
            z += 1
            row_2.update(fluff_2)
            row_2['id_number'] = 'F' + str(x + y + z)
            new_table.append(row_2)
    print(x,y,z)
    for i, d in enumerate(new_table):
        d['id_number'] = 'F' + str(i)
    new_table = pd.DataFrame(new_table)
    new_table = find_anchors(new_table)
    return (True, new_table)


if __name__ == '__main__':
    similar_chromatography = True
    paths = []
    for dir in os.listdir("."):
        if "_Pos" in dir and "mzML" in dir:
            full_path = os.path.join(os.path.abspath("."), dir, "preferred_Feature_table.tsv")
            paths.append(full_path)
    mp.freeze_support()
    join_ftables(paths, similar_chromatography=similar_chromatography)