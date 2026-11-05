'''
Functions to prepare Asari GC input file on Kovats indices.  

Kovats index is 100 * (Alkane number). 
Example format is tab delimited .tsv file:
(test/KovatsIndex_alkanestandards_cuimc.tsv)

Kovats_index	Measured_time(minutes)
1000.0	2.69
1000.0	2.69
1000.0	2.7
'''

def get_index_from_tsv_batches(infile, 
                               batch_start_col=5, 
                               outfile='KovatsIndex_alkanestandards.tsv',
                               sep='\t'):
    '''
    The alkane standards are usually measured per batch. 
    This function converts a typical table to Asari format of Kovats index.
    All columns after batch_start_col are assumed as retention time (min).

    input format:
        Alkane	name	Molecular Formula	Parent m/z	m/z to find	B1a	B1b	B2a	B2b	B3a	B3b
        10	n-decane	C10H22	142.1715	142.1715	2.87	2.87	2.87	2.85	2.85	2.86
        11	n-undecane	C11H24	156.1872	156.1872	3.83	3.83	3.84	3.83	3.83	3.83
        12	n-dodecane	C12H26	170.2029	170.2029	4.57	4.57	4.58	4.58	4.57	4.57
        13	n-tridecane	C13H28	184.2186	184.2186	5.2	5.2	5.2	5.2	5.19	5.2
        14	n-tetradecane	C14H30	198.2342	198.2342	5.76	5.76	5.77	5.77	5.76	5.77
        15	n-pentadecane	C15H32	212.2499	212.2499	6.29	6.29	6.29	6.29	6.29	6.29
    output format:
        Kovats_index	Measured_time(minutes)
        1000.0	2.69
        1000.0	2.69
        1000.0	2.7

    The first column is alkane carbon number; Kovats index = 100 * carbon number.
    Empty or non-numeric retention times are skipped.
    Output is written to outfile, sorted by Kovats index then by batch order.
    Returns outfile path, which can be used as `--kovats` input to asari.
    '''
    import pandas as pd

    df = pd.read_csv(infile, sep=sep)
    alkanes = pd.to_numeric(df.iloc[:, 0], errors='coerce')
    rtimes = df.iloc[:, batch_start_col:].apply(pd.to_numeric, errors='coerce')

    result = []
    for ii in range(df.shape[0]):
        if pd.isna(alkanes.iloc[ii]):
            continue
        kovats_index = 100.0 * float(alkanes.iloc[ii])
        for t in rtimes.iloc[ii].values:
            if not pd.isna(t):
                result.append((kovats_index, float(t)))

    with open(outfile, 'w') as O:
        O.write("Kovats_index\tMeasured_time(minutes)\n")
        for kovats_index, t in result:
            O.write(f"{kovats_index}\t{t}\n")

    print(f"Wrote {len(result)} Kovats index entries to {outfile}.")
    return outfile


def get_tsv_from_xlsx(infile, outfile=None, sheet_name=0):
    '''
    Convert a sheet in an Excel file (.xlsx) to tab delimited .tsv file.
    Default is first sheet; outfile defaults to infile with .tsv extension.
    Requires openpyxl, which pandas uses to read .xlsx files.
    Returns outfile path.
    '''
    import os
    import pandas as pd

    try:
        df = pd.read_excel(infile, sheet_name=sheet_name)
    except ImportError as err:
        raise ImportError(
            "Reading .xlsx requires openpyxl. Install by `pip install openpyxl`."
            ) from err

    if outfile is None:
        outfile = os.path.splitext(infile)[0] + '.tsv'
    df.dropna(how='all').to_csv(outfile, sep='\t', index=False)
    return outfile


def get_index_from_xlsx_batches(infile, 
                                batch_start_col=5, 
                                outfile='KovatsIndex_alkanestandards.tsv'):
    return get_index_from_tsv_batches(
        get_tsv_from_xlsx(infile), 
        batch_start_col=batch_start_col,
        outfile=outfile,
        )

#
# -----------------------------------------------------------
#
def main():
    '''
    Command line usage:
        python -m asari.tools.kovats 240806_Test_Batches_Alkanes.xlsx
        python -m asari.tools.kovats alkanes_batches.tsv -o KovatsIndex.tsv --batch_start_col 5
    Function is chosen by input file extension, .xlsx or .tsv.
    '''
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert alkane standards measured per batch to Asari Kovats index file (tsv).'
        )
    parser.add_argument('infile', help='Input table of alkane retention times per batch, .xlsx or .tsv.')
    parser.add_argument('--batch_start_col', type=int, default=5, 
                        help='0-based index of first column of retention times (default: 5).')
    parser.add_argument('-o', '--outfile', type=str, default='KovatsIndex_alkanestandards.tsv', 
                        help='Output tsv file (default: KovatsIndex_alkanestandards.tsv).')
    args = parser.parse_args()

    suffix = args.infile.lower().rsplit('.', 1)[-1]
    if suffix == 'xlsx':
        get_index_from_xlsx_batches(args.infile, 
                                    batch_start_col=args.batch_start_col, 
                                    outfile=args.outfile)
    elif suffix == 'tsv':
        get_index_from_tsv_batches(args.infile, 
                                   batch_start_col=args.batch_start_col, 
                                   outfile=args.outfile)
    else:
        parser.error(f"Input file must be .xlsx or .tsv, got: {args.infile}")


if __name__ == '__main__':
    main()
