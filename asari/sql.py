'''

Use SQLite3 as backend database to store all sample-related data.

This relieves the memory requirement and provides persistence to stop/resume operations.

Multiple processes can visit the same database, thus enabling parallel processing.


Design
------

One Experiment (project) per database.

Each sample links to two tables: mass_traces and peaks.






Use pandas for DB operations (ORM).

Use JOIN etc for data operation.






For ref DBs, flat the relationship list and use the 1st only in SQL table.
[['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]],

There are cases of non-unique relationships, e.g.
 ['C4H9O5_137.044485',
  137.04448546677,
  'C4H9O5',
  0.832493371827353,
  [('C4H8O5_136.037173', 'M+H[1+]'), ('C4H6O4_118.026609', 'M+H2O+H[1+]')]],

As these are small minority, not dealt with in this version of asari.


'''


import sqlite3 as sql
import pandas as pd
from mass2chem.annotate import compute_adducts_formulae

#
# this is tier 1 pos only, will update the DB in later iterations
from pos_ref_DBs import DB_1

def make_formula_mass_id(formula, mz): return formula + '_' + str(round(mz,6))

def extend_DB1(DB1, mode='pos'):
    flat = []
    for k,v in DB1.items():
        if 40 < k < 2000:
            flat += v
    new = []
    for L in flat:
        neutral_formula, mw = L[4][0][0].split("_")
        adducts = compute_adducts_formulae(float(mw), neutral_formula, mode)    
                                                 # e.g. [(58.53894096677, 'M+2H[2+]', result_formula), ...,]
        for A in adducts:
            new.append( [make_formula_mass_id(A[2], A[0]), A[0], A[2], None, [(L[4][0][0], A[1]),] ] )

    return flat + new


def DB_to_DF(DB):
    '''
    Convenience function to convert indexed DB format to pandas dataframe.
    Input is de-indexed DB: 
    [['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]],...
    '''
    header = 'formula_mass, mz, charged_formula, selectivity, neutral_formula_mass, ion_relation'.split(', ')
    flat = [x[1:4] + list(x[4][0]) for x in DB]
    index = [x[0] for x in DB]
    return pd.DataFrame(flat, columns=header[1:], index=index)
    








