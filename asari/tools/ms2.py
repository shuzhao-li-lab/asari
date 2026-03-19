'''
Various functions

from extracting MS2 spectra from mzML files  

to matching the best to LC-MS1 features from asari.

Example use :
python3 get_ms2_data.py pcpfm_hilicpos_lookahead/asari_asari_project_9620525/export/full_Feature_table.tsv \
    /Volumes/BLUE2T25/LookAHEAD/LCMS/ddMS2_hilicpos/ --rt_tol 30 \ 
    --output_json ms2_hilicpos_asari_project_92125241.json > processing_hilicpos_ms2_log.txt

'''

import os
import json
import numpy as np
import pymzml
from asari.mass_functions import complete_mass_paired_mapping
from .file_io import read_features_from_asari_table

def extract_ms2_form_file(infile, min_intensity=1000):
    '''
    Extract ms2 spectra from mzML file, clean up and return as a list of spectra.

    '''
    ms2_spectra = []
    ms_expt = exp = pymzml.run.Reader(infile)
    for spec in ms_expt:
        if spec.ms_level == 2:
            precursor_mz = spec.selected_precursors[0]['mz']
            if precursor_mz is None:
                continue
            peaks = [(mz, intensity) for mz, intensity in spec.peaks('centroided')
                               if intensity >= min_intensity and mz < precursor_mz-1]    # excluding precursor ion
            if peaks:
                ms2_spectra.append(
                    {
                        'precursor_mz': precursor_mz,
                        'rtime': spec.scan_time_in_minutes() * 60,   # in seconds
                        # 'charge': spec.selected_precursors[0]['charge'],
                        'peaks': peaks
                    }
                )
    return ms2_spectra

def get_top_n_ms2_spectra(ms2_spectra, n=5):
    '''
    Get top n MS2 spectra based on the number of peaks.

    Parameters
    ----------
    ms2_spectra : list
        A list of MS2 spectra.
    n : int
        The number of top MS2 spectra to return.

    Returns
    -------
    top_n_spectra : list
        A list of top n MS2 spectra.
    '''
    if len(ms2_spectra) <= n:
        return ms2_spectra
    else:
        ms2_spectra = sorted(ms2_spectra, key=lambda x: len(x['peaks']), reverse=True)
        return ms2_spectra[:n]

def get_top_n_peaks(spectrum, n=20):
    '''
    Get top n peaks from a spectrum based on intensity.

    Parameters
    ----------
    spectrum : list
        A list of (mz, intensity) tuples.
    n : int
        The number of top peaks to return.

    Returns
    -------
    top_n_peaks : list
        A list of top n (mz, intensity) tuples.
    '''
    if len(spectrum) <= n:
        return spectrum
    else:
        spectrum = sorted(spectrum, key=lambda x: x[1], reverse=True)
        return spectrum[:n] 
    
def get_best_ms2_spectrum(ms2_spectra):
    '''
    Get the best MS2 spectrum based on the intensity of all peaks (better after filtering for top n).

    Parameters
    ----------
    ms2_spectra : list
        A list of MS2 spectra.

    Returns
    -------
    best_spectrum : dict
        The best MS2 spectrum.
    '''
    if not ms2_spectra:
        return None
    ms2_spectra = sorted(ms2_spectra, key=lambda x: sum(p[1] for p in x['peaks']), reverse=True)
    return ms2_spectra[0]

# ~~~~~~~~~~ 

def get_matched_ms2_ms1(LCMS1_features, list_ms2_spectra, rt_tol=30, ppm_tol=5):
    '''
    Get matched MS2 spectra for each LC-MS1 feature based on retention time and m/z tolerance. 
    Similar to asari.tools.match_features.list_match_lcms_features

    Parameters
    ----------
    features : list
        A list of AlignedFeature objects.
    rt_tol : float
        The retention time tolerance in seconds.
    ppm_tol : float
        The m/z tolerance in parts per million (ppm).

    Returns
    -------
    features : list
        A list of AlignedFeature objects with matched MS2 spectra.
    '''
    mz_mapped = complete_mass_paired_mapping([f['mz'] for f in LCMS1_features], 
                                             [f['precursor_mz'] for f in list_ms2_spectra], std_ppm=ppm_tol)[0]
    matched = [x for x in mz_mapped if abs(LCMS1_features[x[0]]['rtime'] - list_ms2_spectra[x[1]]['rtime']) <= rt_tol]
    # select the best MS2 spectrum for each feature
    feature_to_ms2 = {}
    for f_idx, ms2_idx in matched:
        if f_idx not in feature_to_ms2:
            feature_to_ms2[f_idx] = []
        feature_to_ms2[f_idx].append(list_ms2_spectra[ms2_idx])
    for f_idx in feature_to_ms2:
        feature_to_ms2[f_idx] = (get_best_ms2_spectrum(feature_to_ms2[f_idx]), len(feature_to_ms2[f_idx]))

    result = []
    for k, v in feature_to_ms2.items():
        result.append((LCMS1_features[k]['id'], v))

    return result

def match_ms2files_to_features(ms1_fulltable, list_ms2_files, 
                               rt_tol=30, ppm_tol=5, output_json="matched_ms2_spectra.json"  
                               ):
    '''
    '''
    LCMS_features = read_features_from_asari_table(open(ms1_fulltable).read())[1]
    # MS2 record is updated per file, to save memory
    master_dict = {}
    for ms2_file in list_ms2_files:
        print("Processing MS2 file: %s" %ms2_file)
        ms2_spectra = extract_ms2_form_file(ms2_file)
        print("  Extracted %d MS2 spectra" %len(ms2_spectra))
        matched = get_matched_ms2_ms1(LCMS_features, ms2_spectra, rt_tol=rt_tol, ppm_tol=ppm_tol)
        print("  Found %d features with matched MS2 spectra" %len(matched))
        for feature_id, (ms2, count) in matched:
            if feature_id not in master_dict:
                master_dict[feature_id] = (ms2, count)
            else:
                new_count = master_dict[feature_id][1] + count
                new_ms2 = get_best_ms2_spectrum([master_dict[feature_id][0], ms2])
                master_dict[feature_id] = (new_ms2, new_count)
    
    # trim
    for feature_id, pp in master_dict.items():
        master_dict[feature_id] = (get_top_n_peaks(pp[0], 20), pp[1])

    # export JSON
    with open(output_json, "w") as fout:
        json.dump(master_dict, fout, indent=2)
    
    return master_dict
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Match MS2 spectra from mzML files to LC-MS1 features from Asari result table.')
    parser.add_argument('ms1_fulltable', help='Asari feature table file (tsv format).')
    parser.add_argument('ms2_dir', help='Directory containing MS2 mzML files.')
    parser.add_argument('--rt_tol', type=float, default=30, help='Retention time tolerance in seconds (default: 30).')
    parser.add_argument('--ppm_tol', type=float, default=5, help='m/z tolerance in ppm (default: 5).')
    parser.add_argument('--output_json', type=str, default="matched_ms2_spectra.json", help='Output JSON file to save matched MS2 spectra (default: matched_ms2_spectra.json).')
    args = parser.parse_args()
    
    files = os.listdir(args.ms2_dir)
    ms2_files = [os.path.join(args.ms2_dir, f) for f in files if f.endswith(".mzML")]

    _d_ = match_ms2files_to_features(args.ms1_fulltable,
                                     ms2_files, rt_tol=args.rt_tol, ppm_tol=args.ppm_tol,
                                     output_json=args.output_json)
    
