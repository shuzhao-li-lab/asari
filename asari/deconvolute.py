'''
Functions for subcommand `deconvolute`

Supports LC and GC after ver 2.0. 

Will add MS/MS deconvolution in long future.


'''
import os
import json
import time
from .utils import NpEncoder
from .gc_deconvolution import *   


def deconvolute_project(infile, parameters):
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    subdir = '_'.join(['annotation', parameters['project_name'], ''.join(time_stamp)])
    outdir = os.path.join(parameters['outdir'], subdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"Result directory: {outdir}.\n\n")
    
    if parameters['workflow'] == "LC":
        # to add deconvolute and algorithms in version 2
        print("To be implemented.")
        
    elif parameters['workflow'] == "GC":
        list_features, feature_dataframe, ri_model = ingest_gc_featuretable(infile, parameters)
        list_pseudospectra, core_features = run_gc_deconvolution(
            list_features, feature_dataframe, ri_model, parameters['deconvolution_algorithm']
        )
        json_output_file = os.path.join(outdir, parameters['project_name'] + '_denovo_pseudospectra.json')
        write_list_pseudospectra(list_pseudospectra, json_output_file)

    elif parameters['workflow'] == "LCMSMS":
       print("To be implemented.")

    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)


def ingest_gc_featuretable(infile, parameters):
    if check_RI_presence(infile):
        list_features, feature_dataframe = read_feature_table_has_RI(
            infile, data_start_col=12
            )
        ri_model = None     # RI already in table; adap/nmf take PseudoSpectrum RI from base feature
    elif parameters['kovats']:
        list_features, feature_dataframe, ri_model = read_feature_table_and_kovats(
            infile, KovatsIndex_file=parameters['kovats'], 
            frac=0.3, data_start_col=11
            )
    else:
        print("Need retention index values or a KovatsIndex_file to complete this function.")
        return None

    print(f"{len(list_features)} features, {feature_dataframe.shape[1]} samples.")
    return list_features, feature_dataframe, ri_model

def run_gc_deconvolution(list_features, feature_dataframe, ri_model, method): 
    # returns 
    runs = {
        'hcl': lambda: iterative_build_pseudospectra_by_hcl(list_features, feature_dataframe),
        'adap': lambda: adap_gc_deconvolution(list_features, ri_model),
        'nmf': lambda: nmf_deconvolution(list_features, feature_dataframe, ri_model),
    }
    if method not in runs: method = 'hcl'
    list_pseudospectra, core_features = runs[method]()
    print(f"Got {len(list_pseudospectra)} pseudospectra, {len(core_features)} core features.")
    return list_pseudospectra, core_features
