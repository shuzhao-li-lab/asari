#
# I/O and formatting functions
#

def read_features_from_asari_table(text_table, 
                        id_col=0, mz_col=1, rtime_col=2, 
                        left_base=3,
                        right_base=4,
                        parent_masstrack_id=5, 
                        peak_area=6, cSelectivity=7, goodness_fitting=8, snr=9, detection_counts=10,
                        delimiter="\t"):
    '''
    Read a feature table from asari result.
    Returns num_samples, List of peaks.
    '''
    featureLines = text_table.splitlines()
    header = featureLines[0].split(delimiter)
    num_features = len(featureLines)-1
    num_samples = len(header) - 11    # samples start at col 11 in asari output
    # sanity check
    print("table header looks like: \n  ", header[:20])
    print("Read %d feature lines" %num_features)
    L = []
    for ii in range(1, num_features+1):
        if featureLines[ii].strip():
            a = featureLines[ii].split(delimiter)
            L.append({
                'id_number': a[id_col], 
                'id': a[id_col],            # deliberate 
                'mz': float(a[mz_col]), 
                'rtime': float(a[rtime_col]),
                'apex': float(a[rtime_col]),
                'left_base': float(a[left_base]),
                'right_base': float(a[right_base]),
                'parent_masstrack_id': a[parent_masstrack_id],
                'peak_area': float(a[peak_area]),
                'cSelectivity': float(a[cSelectivity]),
                'goodness_fitting': float(a[goodness_fitting]),
                'snr': float(a[snr]),
                'detection_counts': int(a[detection_counts]),
            })
    return num_samples, L

def export_json_to_table(j, outfile, sep='\t'):
    '''
    Convert btw JSON list and SQL styles when compatible (not verifying here), 
    as tables are more readable and smaller in file size.
    The JSON fields are used as table columns. 
    '''
    fields = list(j[0].keys())
    if 'id' in fields:
        fields.remove('id')
        fields = ['id'] + fields
        
    s = sep.join(fields) + '\n'
    for line in j:
        s += sep.join([ str(line[ii]) for ii in fields ]) + '\n'
        
    with open(outfile, 'w') as O:
        O.write(s)

def read_table_to_json(file, sep='\t'):
    '''
    Read btw SQL style table into list_dicts without using pandas, 
    assuming 1st line as header.
    The table columns are used as dict fields.
    Input is treated as str as types are not detected.
    '''
    list_dicts = []
    lines = open(file).readlines()
    header = lines[0].rstrip().split(sep)
    num_fields = len(header)
    for line in lines[1:]:
        a = line.rstrip().split(sep)
        list_dicts.append(dict(zip(header, a)))
        
    return list_dicts
