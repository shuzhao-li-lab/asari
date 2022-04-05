'''
Functions for subcommand `annotate`
'''
from .experiment import *
from .mass_functions import *

from jms.io import read_table_to_peaks
from mass2chem.epdsConstructor import epdsConstructor


def annotate_user_featuretable(infile, parameters):
                        # mode='pos', mz_tolerance_ppm=5):
    '''
    infile: tab delimited file with first row as header, first column m/z and 2nd column rtime.
    output: two files in current directory, Feature_annotation.tsv and Annotated_empricalCompounds.json

    def is_coeluted_by_distance(P1, P2, rt_tolerance=10):
        _coeluted = False
        if abs(P1['apex']-P2['apex']) <= rt_tolerance:
            _coeluted = True
        return _coeluted

    '''
    parameters['outdir'] = ''
    mode = parameters['mode']
    list_peaks = read_table_to_peaks(infile, 
                                has_header=True, mz_col=0, rtime_col=1, feature_id=None ,
                                )
    # print("Read %d features." %len(list_peaks))
    EE = ext_Experiment({}, parameters)
    EE.load_annotation_db()
    ECCON = epdsConstructor(list_peaks, mode=mode)
    EED = ExperimentalEcpdDatabase(mode=mode)
    EED.dict_empCpds = ECCON.peaks_to_epdDict(
                seed_search_patterns = ECCON.seed_search_patterns, 
                ext_search_patterns = ECCON.ext_search_patterns,
                mz_tolerance_ppm= parameters['mz_tolerance_ppm'], 
                coelution_function='distance',
                check_isotope_ratio = False
                ) 
    EED.index_empCpds()
    EED.extend_empCpd_annotation(EE.KCD)
    EED.annotate_singletons(EE.KCD)

    EE.export_peak_annotation(EED.dict_empCpds, EE.KCD, 'Feature_annotation')
    # also exporting JSON
    with open('Annotated_empricalCompounds.json', 'w', encoding='utf-8') as f:
        json.dump(EED.dict_empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)
