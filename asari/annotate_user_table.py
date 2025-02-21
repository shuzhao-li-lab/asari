'''
Functions for subcommand `annotate`
'''
import json

from .experiment import (ext_Experiment,
                         ExperimentalEcpdDatabase, 
                         adduct_search_patterns, 
                         adduct_search_patterns_neg, 
                         isotope_search_patterns)
from .default_parameters import extended_adducts
from .json_encoder import NpEncoder

from jms.io import read_table_to_peaks


def annotate_user_featuretable(infile, parameters, rtime_tolerance=2):
    '''
    Annotate a user supplied feature table. Used by asari subcommand `annotate`.

    Parameters
    ----------
    infile : str
        input feature table filepath, tab delimited file with first row as header, 
        first column m/z and 2nd column rtime.
    parameters : dict
        parameter dictionary passed from main.py, 
        which imports from default_parameters and updates the dict by user arguments.
    rtime_tolerance : int, optional, default: 2
        retention time tolerance to group adducts etc.

    Outputs
    -------
    two files in current directory, Feature_annotation.tsv and Annotated_empricalCompounds.json
    '''
    parameters['outdir'] = ''
    mode = parameters['mode']
    list_peaks = read_table_to_peaks(infile, 
                                has_header=True, mz_col=0, rtime_col=1, feature_id=None ,
                                )
    # print("Read %d features." %len(list_peaks))
    EE = ext_Experiment({}, parameters)
    EE.load_annotation_db()

    EED = ExperimentalEcpdDatabase(mode=mode, 
                                   mz_tolerance_ppm=parameters['mz_tolerance_ppm'],
                                   rt_tolerance=rtime_tolerance)
    
    # passing patterns from .default_parameters
    if mode == 'pos':
        EED.adduct_patterns = adduct_search_patterns
    else:
        EED.adduct_patterns = adduct_search_patterns_neg
    EED.isotope_search_patterns = isotope_search_patterns
    EED.extended_adducts = extended_adducts
    
    EED.build_from_list_peaks(list_peaks)
    EED.extend_empCpd_annotation(EE.KCD)
    EED.annotate_singletons(EE.KCD)
    # EED.dict_empCpds misses some features 
    EED.dict_empCpds = EED.append_orphans_to_epmCpds(EED.dict_empCpds)

    EE.export_peak_annotation(EED.dict_empCpds, EE.KCD, 'Feature_annotation')
    # also exporting JSON
    with open('Annotated_empricalCompounds.json', 'w', encoding='utf-8') as f:
        json.dump(EED.dict_empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)
