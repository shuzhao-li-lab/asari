'''
Convert Asari LC-MS/MS result files to GNPS2 compatible formats. 
convert_asari_featuretable_gnps is based on Joshua Mitchell notebook. 
'''

import pandas as pd

def convert_asari_featuretable_gnps(infile, outfile):
    """Convert Asari feature file to GNPS2 compatible format.
    """
    new_df = pd.read_csv(infile, sep='\t')
    for_GNPS = pd.DataFrame()
    for_GNPS["row ID"]  = list(range(1, new_df.shape[0] + 1))
    for_GNPS["row m/z"] = new_df['mz']
    for_GNPS['row retention time'] = new_df['rtime']
    for z in new_df.columns[11:]:
        for_GNPS[z + ' Peak area'] = new_df[z]
    for_GNPS = for_GNPS[for_GNPS["row ID"] != 0]
    for_GNPS.to_csv(outfile, index=False)
    print(f"Wrote {outfile} with {for_GNPS.shape[0]} Features")


def msp_to_mgf(infile, outfile=None):
    """Convert an MSP file of MS/MS spectra to MGF format for GNPS.

    Parameters
    ----------
    infile : str
        Path to the input .msp file (format produced by json_ms2_to_msp).
    outfile : str, optional
        Path for the output .mgf file. Defaults to infile with .msp replaced by .mgf.

    Returns
    -------
    str
        Path of the written MGF file.
    """
    if outfile is None:
        outfile = infile.replace('.msp', '.mgf')

    with open(infile) as fh:
        text = fh.read()

    # Split into per-spectrum blocks on blank lines
    blocks = [b.strip() for b in text.split('\n\n') if b.strip()]

    lines = []
    for scan_num, block in enumerate(blocks, start=1):
        blines = block.splitlines()
        header, peaks = {}, []
        in_peaks = False
        for line in blines:
            if in_peaks:
                peaks.append(line)
            elif line.lower().startswith('num peaks'):
                in_peaks = True
            else:
                for sep in (':', '='):
                    if sep in line:
                        k, _, v = line.partition(sep)
                        header[k.strip().upper()] = v.strip()
                        break

        precursor_mz = header.get('PRECURSORMZ', '0.0')
        rt = header.get('RETENTIONTIME', '')
        feature_id = header.get('ID', '')

        lines.append('BEGIN IONS')
        lines.append(f'SCANS={scan_num}')
        lines.append(f'PEPMASS={precursor_mz}')
        lines.append('CHARGE=1')
        lines.append('COLLISION_ENERGY=0.0')
        if rt:
            lines.append(f'RTINSECONDS={rt}')
        if feature_id:
            lines.append(f'FEATURE_ID={feature_id}')
        lines.extend(peaks)
        lines.append('END IONS')
        lines.append('')

    with open(outfile, 'w') as fh:
        fh.write('\n'.join(lines))

    return outfile

if __name__ == "__main__":
    '''
    Example use: python3 gnps.py full_Feature_table.tsv ms2_spectra.msp --output testing_gnps2
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Convert Asari LC-MS/MS result files to GNPS2 compatible formats.')
    parser.add_argument('ms1_fulltable', help='Asari feature table file (tsv format).')
    parser.add_argument('ms2_msp', help='MSP file containing MS2 spectra.')
    parser.add_argument('--output', type=str, default="converted_gnps", help='Output file name.')
    args = parser.parse_args()
    
    if '.csv' in args.output:
        outfile_ms1 = args.output
    else:
        outfile_ms1 = args.output + '.csv'
    convert_asari_featuretable_gnps(args.ms1_fulltable, outfile_ms1)

    if '.mgf' in args.output or '.MGF' in args.output:
        outfile_ms2 =  args.output
    else:
        outfile_ms2 =  args.output + '.mgf'
    msp_to_mgf(args.ms2_msp, outfile_ms2)
