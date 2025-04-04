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
    def __init__(self, parameters) -> None:
        self.parameters = parameters
        self.library_path = parameters['GC_Database']
        self.library = self.load_library()
        self.multicores = parameters['multicores']    
    
    @staticmethod
    def extract_spectra(msp_paths):
        combined_spectra = []
        for file in msp_paths:
            for spectrum in load_spectra(file):
                combined_spectra.append(spectrum)
        return combined_spectra

    def load_library(self):
        if os.path.isdir(self.library_path):
            msp_paths = []
            abs_lib_path = os.path.abspath(self.library_path)
            for filename in os.listdir(self.library_path):
                if filename.lower().endswith(".msp"):
                    msp_paths.append(abs_lib_path, filename)
            return self.extract_spectra(msp_paths)
        elif os.path.isfile(self.library_path):
            if self.library_path.lower().endswith(".msp"):
                return self.extract_spectra([os.path.abspath(self.library_path)])
        else:
            raise Exception("Download Not Implemented")
                
    def annotate_gc_feature_table(self, feature_table_path):
        coelute_threshold = self.parameters['coelute_threshold']
        min_peaks = self.parameters['min_peaks']
        min_peaks_common = self.parameters['min_peaks_common']
        min_score_threshold = self.parameters['min_score_threshold']

        raw_ftgraph = FeatureGraph.ftgraph_from_ft(feature_table_path)
        coelute_ftgraph = raw_ftgraph.filter_graph(drt=coelute_threshold)
        extracted_spectra = coelute_ftgraph.extract_fragmentation_spectrum(MIN_PEAKS_EXTRACTION=min_peaks, find_clusters=True)

        print(f"Total Cluster Spectra: {len(extracted_spectra)}")
        print(f"Total Library Spectra: {len(self.library)}")
        print(f"Total Comparisons: {len(extracted_spectra) * len(self.library)}, this may take some time...")

        modes = {
            "cosine": CosineGreedy()
        }
        comparator = partial(modes[self.parameters['similarity_metric']])
        matches = []
        with mp.Pool(self.parameters['multicores']) as pool:
            scores = pool.imap(comparator, product(extracted_spectra, self.library))
            pbar = tqdm.tqdm(scores, desc=f"Comparing Spectra, matches found {len(matches)}", total=len(extracted_spectra) * len(self.library))
            for (extract, library), score in pbar:
                score = score.tolist()
                if len(score) == 2:
                    similarity = score[0]
                    match_peaks = score[1]
                    if similarity >= min_score_threshold and match_peaks >= min_peaks_common:
                        matches.append({
                            "extract": extract,
                            "library": library,
                            "similarity": similarity,
                            "match_peaks": match_peaks,
                            "method": self.parameters['similarity_metric']
                        })
                        pbar.set_description_str(f"Comparing Spectra, matches found {len(matches)}")
        coelute_ftgraph.map_annotations(matches)

def wrapped_comparison(method, job):
    return job, method.pair(job[0], job[1])