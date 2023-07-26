import pandas as pd
import sys
import matplotlib.pyplot as plt
import intervaltree
import numpy as np
from statistics import mode
import numpy as np
import statsmodels.api as sm

lowess = sm.nonparametric.lowess
t1 = pd.read_csv(sys.argv[1], sep="\t")
t2 = pd.read_csv(sys.argv[2], sep="\t")

common_ion_rts = []
for i, (mz1, rt1, pa1, lrt1, rrt1, snr1) in enumerate(zip(t1['mz'], t1['rtime'], t1['peak_area'], t1['rtime_left_base'], t1['rtime_right_base'], t1['snr'])):
    matches_in_t1 = 0
    for mz1_a, rt1_a in zip(t1['mz'], t1['rtime']):
        if mz1 - mz1/1e6 * 5 < mz1_a < mz1 + mz1/1e6 * 5 and abs(rt1 - rt1_a) < .05*rt1:
            matches_in_t1 += 1
    if matches_in_t1 == 1:
        matches_in_t2 = 0
        these_commmon_ion_rts = []
        for j, (mz2, rt2, pa2, snr2) in enumerate(zip(t2['mz'], t2['rtime'], t2['peak_area'], t2['snr'])):
            if mz1 - mz1/1e6 * 5 < mz2 < mz1 + mz1/1e6 *5 and abs(rt1 - rt2) < .05*rt1:
                matches_in_t2 += 1
                these_commmon_ion_rts.append((rt1, rt2, mz1, mz2))
        if len(these_commmon_ion_rts) == 1:
            common_ion_rts.append(these_commmon_ion_rts[0])
common_ion_rts = np.array(common_ion_rts)
last_predicted = lowess(common_ion_rts[:,1], common_ion_rts[:,0], it=50, xvals=common_ion_rts[:,0], frac=1)
last_error = abs(np.sum(common_ion_rts[:,1] - last_predicted))
last_frac = 1.0
frac = 1.0

while frac > 0.02:
    frac -= 0.01
    predicted = lowess(common_ion_rts[:,1], common_ion_rts[:,0], it=50, xvals=common_ion_rts[:,0], frac=frac)
    new_error = np.sum(common_ion_rts[:,1] - predicted)
    has_nan = False
    for value in np.isnan(predicted):
        has_nan = has_nan or value
    if not has_nan and abs(new_error) < last_error:
        full_range = lowess(common_ion_rts[:,1], common_ion_rts[:,0], it=50, xvals=list(range(0,500)), frac=frac)
        full_range_has_nan = False
        for value in np.isnan(predicted):
            full_range_has_nan = full_range_has_nan or value
        if full_range_has_nan is False:
            last_predicted = predicted
            last_frac = frac
    #else:
    #    break

#plt.plot(list(range(0,500)), lowess(common_ion_rts[:,1], common_ion_rts[:,0], it=50, xvals=list(range(0,500)), frac=last_frac))    
#plt.scatter(common_ion_rts[:,0], last_predicted, c='g')
#plt.scatter(common_ion_rts[:,0], common_ion_rts[:,1], c='k')
#plt.show()
#plt.hist(common_ion_rts[:,1] - last_predicted)
#plt.show()

t1['pred_t2_rtime'] = lowess(common_ion_rts[:,1], common_ion_rts[:,0], it=50, xvals=t1['rtime'], frac=last_frac)
#plt.scatter(t1['rtime'], t1['pred_t2_rtime'])
#plt.show()

skip_columns = ['parent_masstrack_id', 'peak_area', 'cSelectivity', 'goodness_fitting', 'snr']
common_columns = set(list(t1.columns)).intersection(set(list(t2.columns)))
all_columns = set(list(t1.columns)).intersection(set(list(t2.columns)))
t1_unique_columns = set([x for x in t1.columns if x not in common_columns])
t2_unique_columns = set([x for x in t2.columns if x not in common_columns])

def extract_row(row, columns):
    return {column: row[column] for column in columns}

used_t1 = set()
used_t2 = set()
new_table = []
for row_1 in t1.apply(extract_row, axis=1, args=([c for c in t1.columns if c not in skip_columns],)):
    mz1 = row_1['mz']
    row_1['rtime_left_base'] = round(row_1['rtime_left_base'] - (row_1['rtime'] - row_1['pred_t2_rtime']),2)
    row_1['rtime_right_base'] = round(row_1['rtime_right_base'] - (row_1['rtime'] - row_1['pred_t2_rtime']), 2)
    row_1['rtime'] = round(row_1['pred_t2_rtime'], 2)
    for row_2 in t2.apply(extract_row, axis=1, args=([c for c in t2.columns if c not in skip_columns],)):
        mz2 = row_2['mz']
        if mz1 - mz1/1e6 * 5 < mz2 < mz1 + mz1/1e6 * 5:
            if row_2['rtime_left_base'] < row_1['rtime'] < row_2['rtime_right_base']:
                print(row_1['id_number'], row_2['id_number'])
                used_t1.add(row_1['id_number'])
                used_t2.add(row_2['id_number'])
                #todo - handle duplicate matches
                new_row = dict(row_1)
                new_row.update(row_2)
                new_row['detection_counts'] = row_1['detection_counts'] + row_2['detection_counts']
                del new_row['pred_t2_rtime']
                new_table.append(new_row)

for row_1 in t1.apply(extract_row, axis=1, args=([c for c in t1.columns if c not in skip_columns],)):
    row_1['rtime_left_base'] = round(row_1['rtime_left_base'] - (row_1['rtime'] - row_1['pred_t2_rtime']),2)
    row_1['rtime_right_base'] = round(row_1['rtime_right_base'] - (row_1['rtime'] - row_1['pred_t2_rtime']), 2)
    row_1['rtime'] = round(row_1['pred_t2_rtime'], 2)
    if row_1['id_number'] not in used_t1:
        for column in t2_unique_columns:
            row_1[column] = 0
        del row_1['pred_t2_rtime']
        new_table.append(row_1)

for row_2 in t2.apply(extract_row, axis=1, args=([c for c in t2.columns if c not in skip_columns],)):
    if row_2['id_number'] not in used_t2:
        for column in t1_unique_columns:
            row_2[column] = 0
        new_table.append(row_2)
    
for i, d in enumerate(new_table):
    d['id_number'] = 'F' + str(i)


new_table = pd.DataFrame(new_table)
new_table.to_csv("new_table.tsv", sep="\t")
            
