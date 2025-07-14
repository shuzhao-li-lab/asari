import multiprocessing as mp
from functools import cache, partial
from itertools import product
import os
import json

import matchms
from matchms.similarity import CosineGreedy
from matchms.filtering import default_filters, normalize_intensities
from matchms.importing import load_spectra


import tqdm
import importlib.resources as pkg_resources

try:
    from .utils import download_and_unzip_to_pkg_resources
    from .feature_graph import FeatureGraph
except ImportError:
    # for testing only, remove before release
    from utils import download_and_unzip_to_pkg_resources
    from feature_graph import FeatureGraph

class EI_MS_Library():
    """
    Handles loading and annotation of EI-MS spectra from a local or remote MoNA library.
    """
    MONA_MSP_URL = (
        "https://mona.fiehnlab.ucdavis.edu/rest/downloads/retrieve/a09f9652-c7bc-48b2-9dd9-0dc4343bb360"
    )

    def __init__(self, parameters) -> None:
        self.parameters = parameters
        # Allow overriding via CLI (--GC_Database)
        user_db = self.parameters.get('GC_Database')
        if user_db:
            self.library_path = user_db
        else:
            # Default location inside package resources under `db/`
            pkg_root = pkg_resources.files(__name__).parent
            db_dir = pkg_root / 'db'
            self.library_path = str(db_dir)
        self.multicores = parameters.get('multicores', mp.cpu_count())
        self._ensure_library()
        self.library = self.load_library()
    
    def _ensure_library(self):
        """
        Ensure the MSP library exists locally. If not, download and extract from MoNA.
        """
        # If path is a file or directory and exists, nothing to do
        if os.path.exists(self.library_path):
            return

        # Otherwise, attempt to download and unzip the MoNA export
        try:
            raise NotImplementedError
            # print(f"Downloading MoNA GC-MS library from {self.MONA_MSP_URL} ...")
            # download_and_unzip_to_pkg_resources will extract under <package>/data by default;
            # here we target subdir 'db'
            # download_and_unzip_to_pkg_resources(
            #     self.MONA_MSP_URL,
            #     package=__name__.split('.')[0],
            #     subdir='db'
            # )
            # After extraction, the .msp file should be in the same path
        except Exception as e:
            raise RuntimeError(
                f"Failed to download or extract MoNA library: {e}"
            )

    @staticmethod
    def extract_spectra(msp_paths):
        combined_spectra = []
        for file in msp_paths:
            for spectrum in load_spectra(file):
                combined_spectra.append(spectrum)
        return combined_spectra

    def load_library(self):
        """
        Load spectra from MSP files in the library_path.
        Supports both a single .msp file or a directory containing multiple .msp files.
        """
        if os.path.isdir(self.library_path):
            msp_paths = []
            for filename in os.listdir(self.library_path):
                if filename.lower().endswith('.msp'):
                    msp_paths.append(
                        os.path.join(self.library_path, filename)
                    )
            return self.extract_spectra(msp_paths)

        elif os.path.isfile(self.library_path):
            if self.library_path.lower().endswith('.msp'):
                return self.extract_spectra([self.library_path])

        raise FileNotFoundError(
            f"EI-MS library not found at {self.library_path}."
        )
                
    def annotate_gc_feature_table(self, feature_table_path, coelute_threshold=None):
        if coelute_threshold is None:
            coelute_threshold = self.parameters['coelute_threshold']
        min_peaks = self.parameters.get('min_peaks', 1)
        min_peaks_common = self.parameters.get('min_peaks_common', 1)
        min_score_threshold = self.parameters.get('min_score_threshold', 0.70)

        raw_ftgraph = FeatureGraph.ftgraph_from_ft(feature_table_path)
        coelute_ftgraph = raw_ftgraph.filter_graph(drt=coelute_threshold)
        extracted_spectra = coelute_ftgraph.extract_fragmentation_spectrum(MIN_PEAKS_EXTRACTION=min_peaks, find_clusters=True)

        print(f"Total Cluster Spectra: {len(extracted_spectra)}")
        print(f"Total Library Spectra: {len(self.library)}")
        print(f"Total Comparisons: {len(extracted_spectra) * len(self.library)}, this may take some time...")

        comparator = CosineGreedy()
        matches = []
        pbar = tqdm.tqdm(extracted_spectra)
        for extracted_spectrum in pbar:
            for lib in self.library:
                score = comparator.pair(extracted_spectrum, lib)
                score = score.tolist()
                if len(score) == 2:
                    similarity = score[0]
                    match_peaks = score[1]
                    if similarity >= min_score_threshold and match_peaks >= min_peaks_common:
                        matches.append({
                            "extract": extracted_spectrum,
                            "library": lib,
                            "similarity": similarity,
                            "match_peaks": match_peaks,
                            "method": self.parameters['similarity_metric']
                        })
                        pbar.set_description_str(f"Comparing Spectra, matches found {len(matches)}")
        coelute_ftgraph.map_annotations(matches)

def wrapped_comparison(method, job):
    return job, method.pair(job[0], job[1])