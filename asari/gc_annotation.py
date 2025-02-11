import multiprocessing as mp
from functools import cache
from itertools import product
import os
import json

import matchms
from matchms.similarity import CosineGreedy
from matchms.filtering import default_filters, normalize_intensities
import tqdm
import importlib.resources as pkg_resources


try:
    from .utils import download_and_unzip_to_pkg_resources
    from .feature_graph import FeatureGraph
except:
    # for testing only, remove before release
    from utils import download_and_unzip_to_pkg_resources
    from feature_graph import FeatureGraph


class EI_MS_Library():
    def __init__(self, library_ID) -> None:
        self.library_meta = self.retrieve_library_meta(library_ID)
        self.library = None
        self.load_library()
        
    def load_library(self, limit=None):
        extension = self.library_meta['Extension'][1:] if self.library_meta['Extension'][0] == "." else self.library_meta['Extension']
        if self.library_meta.get("Parser", "").lower() == "matchms" and hasattr(matchms.importing, f"load_from_{extension}"):
            self.loader = getattr(matchms.importing, f"load_from_{extension}")
            print(f"Using load_from_{extension} in MatchMS to load library")
        else:
            raise NotImplementedError
        
        if limit:
            library = []
            for x in self.loader(self.library_meta['LIB_PATH']):
                library.append(x)
                if len(library) > 1000:
                    break
        else:
            library = list(self.loader(self.library_meta['LIB_PATH']))

        processed_library = []
        for spectrum in library:
            spectrum = default_filters(spectrum)
            spectrum = normalize_intensities(spectrum)
            processed_library.append(spectrum)

        self.library = processed_library

    def retrieve_library_meta(self, library_ID):
        EI_MS_Library_manifest = EI_MS_Library.load_library_manifest()
        if os.path.exists(library_ID) and os.path.isfile(library_ID):
            print("Assuming provided file is the library")
            print("Assuming MatchMS is okay for parsing if in JSON, MSP, or MGF format")
            library_to_load = {
                "Name": library_ID,
                "Description": "User Provided Library",
                "URL": None,
                "OnDiskName": os.path.basename(library_ID),
                "Extension": os.path.splitext(library_ID)[1],
                "Parser": "matchms",
                "License": "Unknown",
                "Required Citations": "Unknown",
                "LIB_PATH": library_ID}
            return library_to_load
        elif library_ID in EI_MS_Library_manifest["EI_MS"]:
            library_to_load = EI_MS_Library_manifest["EI_MS"][library_ID]
            on_disk_name = library_to_load['OnDiskName']
            on_disk_path = os.path.join(os.path.dirname(pkg_resources.files('asari')), 'db', on_disk_name)
            if not os.path.exists(on_disk_path):
                if library_to_load["URL"].endswith('zip'):
                    download_and_unzip_to_pkg_resources(library_to_load['URL'], 'asari', 'db')
                elif library_to_load.get('Compression', None) == 'zip':
                    download_and_unzip_to_pkg_resources(library_to_load['URL'], 'asari', 'db')
                assert os.path.exists(on_disk_path), f"Library not found at {on_disk_path} after Download"
                assert "LIB_PATH" not in library_to_load, "Library already loaded!"
                library_to_load["LIB_PATH"] = on_disk_path
                return library_to_load
            else:
                print("Exists")
                library_to_load["LIB_PATH"] = on_disk_path
                return library_to_load
        else:
            print(f"Library ID {library_ID} not found")
            print("Valid Selections:\n\t" + "\n\t".join(list(EI_MS_Library_manifest['EI_MS'].keys())))
            raise ValueError

    # to prevent multiple hits to disk, structure is immutable so this is safe
    @cache 
    @staticmethod
    def load_library_manifest():
        return json.load(open(os.path.join(pkg_resources.files('asari'), 'db', 'gcms_libraries.json')))
    
    def annotate_gc_feature_table(self, feature_table_path, drt=0.5, min_peaks=3, min_shared_peaks=1, min_score_threshold=0.7):
        
        raw_ftgraph = FeatureGraph.ftgraph_from_ft(feature_table_path)
        coelute_ftgraph = raw_ftgraph.filter_graph(drt=drt)
        extracted_spectra = coelute_ftgraph.extract_fragmentation_spectrum(MIN_PEAKS_EXTRACTION=min_peaks, find_clusters=True)
        total = len(extracted_spectra) * len(self.library)
        print(f"Total Cluster Spectra: {len(extracted_spectra)}")
        print(f"Total Library Spectra: {len(self.library)}")
        print(f"Total Comparisons: {total}, this may take some time...")
        matches = []
        with mp.Pool(mp.cpu_count()) as pool:
            scores = pool.imap(wrapped_cosine, product(extracted_spectra, self.library))
            pbar = tqdm.tqdm(scores, desc=f"Comparing Spectra, matches found {len(matches)}", total=total)
            for (extract, library), score in pbar:
                score = score.tolist()
                if len(score) == 2:
                    similarity = score[0]
                    match_peaks = score[1]
                    if similarity > min_score_threshold and match_peaks >= min_shared_peaks:
                        matches.append({
                            "extract": extract,
                            "library": library,
                            "similarity": similarity,
                            "match_peaks": match_peaks
                        })
                        pbar.set_description_str(f"Comparing Spectra, matches found {len(matches)}")
        coelute_ftgraph.map_annotations(matches)


    @staticmethod
    def annotate_gc_feature_table_with_library(feature_table_path, library_ID):
        library = EI_MS_Library(library_ID)
        library.annotate_gc_feature_table(feature_table_path)

def wrapped_cosine(job):
    return job, CosineGreedy().pair(job[0], job[1])

if __name__ == '__main__':
    EI_MS_Library.annotate_gc_feature_table_with_library('/Users/mitchjo/ATLAS/01172025_Oxygen_12C13C_isotope_labeling/Cellpellets/HILICneg/output_asari_project_12721586/preferred_Feature_table.tsv', "MoNA_GCMS")
