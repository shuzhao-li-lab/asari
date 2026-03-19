'''

Example use :
python3 enp_search.py --ms2_json ms2_hilicpos_asari_project_92125241.json --db_pickle src_annotation/metabolomics_positive_for_masscube_v12.pkl --output_json matched_hilicpos_asari_project_92125241.json
python3 enp_search.py --ms2_json ms2_hilicneg_asari_project_92123515.json --db_pickle src_annotation/metabolomics_negative_for_masscube_v12.pkl --output_json matched_hilicneg_asari_project_92123515.json

'''

import json
import numpy as np
from tqdm import tqdm
import pickle
from ms_entropy import FlashEntropySearch


param = {
    "ion_mode": "positive",
    "mz_tol_ms1": 0.01,
    "mz_tol_ms2": 0.02,
    "ms2_sim_tol": 0.6,
    'precursor_mz_offset': 1.6
}


def get_cleaned_experimental_ms2(dict_ms2_spectra, entropy_search, params={}):
    '''
    Clean experimental MS2 spectra for searching.

    Parameters
    ----------
    dict_ms2_spectra : dict
        A dictionary of MS2 spectra. Keys are feature IDs, values are lists of [MS2 spectra, counts].
        Precursor peaks are already removed in my preprocessing step (offset 1).
    entropy_search : EntropySearch object
        The EntropySearch object containing the MS2 database.
    params : dict
        The parameters for cleaning.

    Returns
    -------
    cleaned_spectra : list
        A list of cleaned MS2 spectra.
    '''
    cleaned_spectra = []
    for k, sp in dict_ms2_spectra.items():
        precursor_mz = sp[0]['precursor_mz']
        peaks = np.array(sp[0]['peaks'])
        cleaned_peaks = entropy_search.clean_spectrum_for_search(precursor_mz, peaks, 
                                        precursor_ions_removal_da=params['precursor_mz_offset'])
        if len(cleaned_peaks) > 0:
            cleaned_spectra.append(
                (k, precursor_mz, cleaned_peaks)
            )
    return cleaned_spectra

def search_ms2_spectra(cleaned_spectra, entropy_search, params={}):
    '''
    Search MS2 spectra against the database.

    Parameters
    ----------
    cleaned_spectra : list
        A list of cleaned MS2 spectra. Each item is a tuple of (feature ID, precursor m/z, peaks).
    entropy_search : EntropySearch object
        The EntropySearch object containing the MS2 database.
    params : dict
        The parameters for searching.

    Returns
    -------
    results : dict
        A dictionary of search results. Keys are feature IDs, values are best matched database entries.
    '''
    results = {}
    for feature_id, precursor_mz, peaks in tqdm(cleaned_spectra, desc="Searching MS2 spectra"):
        similarity, matched_num = entropy_search.identity_search(precursor_mz=precursor_mz, peaks=peaks,
                                                                 ms1_tolerance_in_da=params['mz_tol_ms1'],
                                                                 ms2_tolerance_in_da=params['mz_tol_ms2'],
                                                                 output_matched_peak_number=True)
        if np.max(similarity) > params['ms2_sim_tol']:
            best_idx = np.argmax(similarity)
            matched_entry = entropy_search[best_idx]
            _ = matched_entry.pop('peaks') # np.array hard to JSON serialize
            results[feature_id] = (matched_entry, float(similarity[best_idx]), int(matched_num[best_idx]))
    return results

#
#  ---
# 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Search MS2 spectra against database.')
    parser.add_argument('--ms2_json', type=str, default='extracted_ms2_spectra.json', help='Input JSON file containing extracted MS2 spectra.')
    parser.add_argument('--db_pickle', type=str, default='MassCube_DB_2025-06-25.pkl', help='Pickle file containing the MS2 database.')
    parser.add_argument('--output_json', type=str, default='ms2_search_results.json', help='Output JSON file to save search results.')
    args = parser.parse_args()

    # load MS2 database
    entropy_search = pickle.load(open(args.db_pickle, 'rb'))
    # load MS2 spectra
    dict_ms2_spectra = json.load(open(args.ms2_json))
    print(f"Loaded {len(dict_ms2_spectra)} MS2 spectra for searching.")

    # clean MS2 spectra
    cleaned_spectra = get_cleaned_experimental_ms2(dict_ms2_spectra, entropy_search, param)

    # search MS2 spectra
    results = search_ms2_spectra(cleaned_spectra, entropy_search, param)

    # save results
    json.dump(results, open(args.output_json, 'w'), indent=2)

    print(f"Total {len(results)} MS2 spectra matched to database entries.")

