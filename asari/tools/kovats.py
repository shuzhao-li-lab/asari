"""
Prepare the asari GC Kovats index file: convert a table of alkane standards measured in multiple batches to the asari Kovats index format.

Input (.xlsx, .tsv or .csv), one row per alkane, one column of retention times (min) per batch:

    Alkane_N  name           Molecular Formula  Parent m/z  m/z to find  RT_mins_B1  RT_mins_B2  ...
    10        decane         C10H22             170.2029    170.2029    2.69        2.69        ...
    11        n-undecane     C11H24             156.1872    156.1872    3.01        3.01        ...

Output (tab delimited), all batches stacked, sorted by Kovats index then batch order,
for use as `asari process --workflow GC --kovats <file>`:

    Kovats_index	Measured_time(minutes)
    1000.0	2.69
    1000.0	2.69

Kovats index = 100 * carbon number. The carbon number is taken from the molecular formula column
if present (names in such tables are often copy-paste errors), else from the first numeric column
(e.g. 'Alkane' or 'Alkane_N'); a mismatch between the two is reported.
Batch columns are, by default, the numeric columns after the last descriptive column
(name, formula, m/z); use --batch_start_col or --batch_cols to set them explicitly.
Empty or non-numeric retention times are skipped. Retention times that do not increase
with carbon number within a batch are reported as warnings.

Usage:
    python -m asari.tools.kovats Alkane_LA_250529.xlsx -o KovatsIndex_alkanestandards.tsv
    python -m asari.tools.kovats alkanes.tsv --batch_cols 'B\\d'
"""

import os
import re
import sys
import pandas as pd

DESCRIPTIVE = re.compile(r'name|formula|m/?z|mass|smiles|inchi|cas', re.I)
FORMULA = re.compile(r'^\s*C(\d*)H\d+\s*$')


def read_table(infile, sheet_name=0):
    """Read .xlsx/.xls (requires openpyxl), .tsv/.txt or .csv into a DataFrame; empty rows/columns dropped."""
    ext = os.path.splitext(infile)[1].lower()
    if ext in ('.xlsx', '.xls'):
        try:
            df = pd.read_excel(infile, sheet_name=sheet_name)
        except ImportError as err:
            raise ImportError("Reading .xlsx requires openpyxl: `pip install openpyxl`.") from err
    elif ext == '.csv':
        df = pd.read_csv(infile)
    else:
        df = pd.read_csv(infile, sep='\t')
    df = df.dropna(how='all').dropna(axis=1, how='all')
    df.columns = [str(c).strip() for c in df.columns]
    return df


def carbon_numbers(df):
    """Carbon number per row: from a molecular formula column (CnH2n+2) if found, else the first numeric column."""
    from_formula, from_col = None, None
    for c in df.columns:
        vals = df[c].astype(str).str.extract(FORMULA)[0]
        if vals.notna().mean() > 0.5:
            from_formula = pd.to_numeric(vals.replace('', '1'), errors='coerce')
            break
    for c in df.columns:
        if DESCRIPTIVE.search(c):
            continue
        nums = pd.to_numeric(df[c], errors='coerce')
        if nums.notna().mean() > 0.5 and (nums.dropna() % 1 == 0).all():
            from_col = nums
            break
    if from_formula is not None and from_col is not None:
        diff = (from_formula != from_col) & from_formula.notna() & from_col.notna()
        for ii in df.index[diff]:
            print(f"Warning: row {ii + 2}, carbon number {from_col[ii]:g} does not match formula "
                  f"C{from_formula[ii]:g}; using formula.")
    if from_formula is not None:
        return from_formula.fillna(from_col) if from_col is not None else from_formula
    if from_col is not None:
        return from_col
    raise ValueError("No molecular formula or alkane carbon number column found.")


def batch_columns(df, batch_start_col=None, batch_cols=None):
    """Retention time columns: explicit start index, header regex, or numeric columns after the last descriptive one."""
    if batch_start_col is not None:
        return list(df.columns[batch_start_col:])
    if batch_cols:
        cols = [c for c in df.columns if re.search(batch_cols, c)]
    else:
        last = max([ii for ii, c in enumerate(df.columns) if DESCRIPTIVE.search(c)], default=0)
        cols = [c for c in df.columns[last + 1:]
                if pd.to_numeric(df[c], errors='coerce').notna().any()]
    if not cols:
        raise ValueError("No batch columns of retention times found; set --batch_start_col or --batch_cols.")
    return cols


def format_kovats(infile, outfile='KovatsIndex_alkanestandards.tsv',
                  batch_start_col=None, batch_cols=None, sheet_name=0):
    """Convert an alkane standards table with batches in columns to asari Kovats index tsv. Returns outfile."""
    df = read_table(infile, sheet_name=sheet_name)
    carbons = carbon_numbers(df)
    cols = batch_columns(df, batch_start_col, batch_cols)
    rtimes = df[cols].apply(pd.to_numeric, errors='coerce')
    keep = carbons.notna()
    carbons, rtimes = carbons[keep], rtimes[keep]
    print(f"{len(cols)} batches: {', '.join(cols)}")
    print(f"{len(carbons)} alkanes, C{carbons.min():g} to C{carbons.max():g}")

    order = carbons.sort_values(kind='stable').index
    for c in cols:
        rt = rtimes.loc[order, c].dropna()
        if (rt.diff().dropna() <= 0).any():
            print(f"Warning: batch {c}, retention times do not increase with carbon number.")

    result = []
    for ii in order:
        kovats_index = 100.0 * float(carbons[ii])
        for c in cols:
            t = rtimes.at[ii, c]
            if not pd.isna(t):
                result.append((kovats_index, round(float(t), 4)))   # drop float noise from spreadsheets

    with open(outfile, 'w') as O:
        O.write("Kovats_index\tMeasured_time(minutes)\n")
        for kovats_index, t in result:
            O.write(f"{kovats_index}\t{t}\n")
    print(f"Wrote {len(result)} Kovats index entries to {outfile}.")
    return outfile


# names kept from earlier versions; batch_start_col=None detects batch columns automatically
def get_index_from_tsv_batches(infile, batch_start_col=None,
                               outfile='KovatsIndex_alkanestandards.tsv', sep='\t'):
    """Kovats index file from a .tsv/.csv table of alkane standard batches. See format_kovats."""
    return format_kovats(infile, outfile=outfile, batch_start_col=batch_start_col)


def get_index_from_xlsx_batches(infile, batch_start_col=None,
                                outfile='KovatsIndex_alkanestandards.tsv'):
    """Kovats index file from an .xlsx table of alkane standard batches. See format_kovats."""
    return format_kovats(infile, outfile=outfile, batch_start_col=batch_start_col)


def get_tsv_from_xlsx(infile, outfile=None, sheet_name=0):
    """Convert a sheet in an Excel file (.xlsx) to tab delimited .tsv; outfile defaults to infile.tsv."""
    if outfile is None:
        outfile = os.path.splitext(infile)[0] + '.tsv'
    read_table(infile, sheet_name=sheet_name).to_csv(outfile, sep='\t', index=False)
    return outfile


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert alkane standards with batches in columns to asari Kovats index file (tsv).')
    parser.add_argument('infile', help='Alkane standards table, .xlsx, .tsv or .csv.')
    parser.add_argument('-o', '--outfile', default='KovatsIndex_alkanestandards.tsv',
                        help='Output tsv file (default: KovatsIndex_alkanestandards.tsv).')
    parser.add_argument('--batch_start_col', type=int, default=None,
                        help='0-based index of the first retention time column; all columns after are batches.')
    parser.add_argument('--batch_cols', default=None,
                        help="Regex on column headers to select batch columns, e.g. 'RT|B\\d'.")
    parser.add_argument('--sheet', default=0, help='Sheet name or 0-based index for .xlsx (default: 0).')
    args = parser.parse_args()
    sheet = int(args.sheet) if str(args.sheet).isdigit() else args.sheet
    try:
        format_kovats(args.infile, args.outfile, args.batch_start_col, args.batch_cols, sheet)
    except ValueError as err:
        sys.exit(f"Error: {err}")


if __name__ == '__main__':
    main()
