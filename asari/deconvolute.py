"""
Functions for subcommand `deconvolute`

Supports LC and GC after ver 2.0.
Will add MS/MS deconvolution in long future.
"""
import os
import json
import time
from khipu.extended import (peaklist_to_khipu_list, export_empCpd_khipu_list, export_khipu_table)
from khipu.utils import make_peak_dict
from .utils import NpEncoder
from .tools.file_io import read_features_from_asari_table
from .gc_deconvolution import *
from .default_parameters import (
        PROTON, adduct_search_patterns_pos, adduct_search_patterns_neg, isotope_search_patterns, extended_adducts
)


def deconvolute_project(infile, parameters):
    """Run subcommand `deconvolute` on infile, by parameters['workflow']."""
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    subdir = '_'.join([parameters['project_name'], 'deconvolution', ''.join(time_stamp)])
    outdir = os.path.join(parameters['outdir'], subdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"Result directory: {outdir}.\n\n")

    if parameters['workflow'] == "LC" or parameters['workflow'] == "LCMSMS":
        # deconvolute LC-MS using khipu
        json_output_file = os.path.join(outdir, parameters['project_name']+'_list_khipus.json')
        khipu_deconvolute_featuretable(infile, parameters, json_output_file)
        print(f"Results were written to {json_output_file} and .tsv.")

    elif parameters['workflow'] == "GC":
        list_features, feature_dataframe, ri_model = ingest_gc_featuretable(infile, parameters)
        list_pseudospectra, core_features = run_gc_deconvolution(
            list_features, feature_dataframe, ri_model, parameters['deconvolution_algorithm'],
            parameters=parameters
        )
        json_output_file = os.path.join(outdir, parameters['project_name'] + '_denovo_pseudospectra.json')
        write_list_pseudospectra(list_pseudospectra, json_output_file)

    else:
       print("deconvolute applies to LC, LCMSMS or GC workflow.")

    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)


# ---------------------------------------------------------------------------
# GC-MS deconvolution 
# ---------------------------------------------------------------------------

def ingest_gc_featuretable(infile, parameters):
    """Read a GC feature table with RI, from an RI column or a Kovats standards file.

    Returns list_features, feature_dataframe, ri_model (None if RI is from the table).
    """
    if check_RI_presence(infile):
        list_features, feature_dataframe = read_feature_table_has_RI(
            infile, data_start_col=12
            )
        ri_model = None     # RI already in table; adap/nmf take PseudoSpectrum RI from base feature
    elif parameters['kovats']:
        list_features, feature_dataframe, ri_model = read_feature_table_and_kovats(
            infile, KovatsIndex_file=parameters['kovats'],
            data_start_col=11
            )
    else:
        print("Need retention index values or a KovatsIndex_file to complete this function.")
        return None

    print(f"{len(list_features)} features, {feature_dataframe.shape[1]} samples.")
    return list_features, feature_dataframe, ri_model

def run_gc_deconvolution(list_features, feature_dataframe, ri_model, method, parameters=None):
    """Run GC deconvolution by method ('hcl', 'adap' or 'nmf'). Returns list_pseudospectra, core_features."""
    limit = (parameters or {}).get('max_core_features', 20000)
    # returns
    runs = {
        'hcl': lambda: iterative_build_pseudospectra_by_hcl(
            list_features, feature_dataframe, max_core_features=limit),
        'adap': lambda: adap_gc_deconvolution(list_features, ri_model),
        'nmf': lambda: nmf_deconvolution(list_features, feature_dataframe, ri_model),
    }
    if method not in runs: method = 'hcl'
    list_pseudospectra, core_features = runs[method]()
    print(f"Got {len(list_pseudospectra)} pseudospectra, {len(core_features)} core features.")
    return list_pseudospectra, core_features


# ---------------------------------------------------------------------------
# LC-MS deconvolution using Khipu
# ---------------------------------------------------------------------------


def khipu_feature_adaptor(list_features, intensity_key='peak_area'):
    """
    Add the intensity names khipu uses to asari features, from intensity_key (default peak_area):
    'representative_intensity' (khipu model, M0 selection) and 'height' (isotope ratio check,
    `check_isotope_ratio`). Existing values are kept, e.g. a feature that already has 'height'.
    """
    for f in list_features:
        intensity = f.get(intensity_key, f.get('representative_intensity', f.get('height', 0)))
        f.setdefault('representative_intensity', intensity)
        f.setdefault('height', intensity)
    return list_features


def khipu_deconvolute_featuretable(infile, parameters, json_output_file):
    """
    LC-MS deconvolution using Khipu, grouping isotopes and adducts of the same compound
    into empirical compounds (empCpds). No database is used.

    Parameters
    ----------
    infile : str
        input feature table filepath, asari format (Feature_table.tsv),
        tab delimited with header; columns id, mz, rtime, ..., parent_masstrack_id, peak_area, ...
    parameters : dict
        parameter dictionary passed from main.py,
        which imports from default_parameters and updates the dict by user arguments.
        Uses mode, mz_tolerance_ppm, khipu_rtime_tolerance (seconds) and check_isotope_ratio
        (isotope abundance ratio limits in isotope_search_patterns, on peak_area).
    json_output_file : str
        output JSON filepath. A tsv of the same name is also written.

    Outputs
    -------
    JSON, list of empCpds (khipu format, with MS1_pseudo_Spectra, neutral_formula_mass etc);
    tsv, one row per feature, with empCpd assignment, isotope, modification and ion_relation;
    features not in any khipu are listed as singletons.

    Returns
    -------
    list of empCpds.
    """
    mode = parameters['mode']
    if mode == 'pos':
        adduct_search_patterns = adduct_search_patterns_pos
        singleton_ion, delta_mz = 'M+H+', PROTON
    else:
        adduct_search_patterns = adduct_search_patterns_neg
        singleton_ion, delta_mz = 'M-H-', -PROTON

    with open(infile) as f:
        _n, list_features = read_features_from_asari_table(f.read())
    list_features = khipu_feature_adaptor(list_features)
    if parameters['check_isotope_ratio']:
        print("Checking isotope abundance ratios (13C etc.), using feature peak_area as intensity.")
    peak_dict = make_peak_dict(list_features)

    list_khipus, _ = peaklist_to_khipu_list(
                            list_features,
                            isotope_search_patterns=isotope_search_patterns,
                            adduct_search_patterns=adduct_search_patterns,
                            extended_adducts=extended_adducts,
                            mz_tolerance_ppm=parameters['mz_tolerance_ppm'],
                            rt_tolerance=parameters['khipu_rtime_tolerance'],
                            check_isotope_ratio=parameters['check_isotope_ratio'],
                            mode=mode,
                            charges=[1, 2, 3],
                            has_parent_masstrack=True,
                            )
    empCpds = export_empCpd_khipu_list(list_khipus)
    num_in_khipus = sum(len(E['MS1_pseudo_Spectra']) for E in empCpds)
    print(f"Got {len(empCpds)} khipus, with {num_in_khipus} of {len(list_features)} features; "
          f"the rest are singletons in the tsv.")

    with open(json_output_file, 'w', encoding='utf-8') as f:
        json.dump(empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)
    export_khipu_table(json_output_file.replace('.json', '.tsv'), peak_dict, empCpds, singleton_ion, delta_mz)
    return empCpds
