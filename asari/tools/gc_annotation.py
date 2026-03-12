#
# not using now
#

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
except ImportError:
    # for testing only, remove before release
    from utils import download_and_unzip_to_pkg_resources
    from feature_graph import FeatureGraph


class GC_Annotation:
    '''
    Do annotation here instead of during processing
    
    '''



    def populate_RI_lookup(self, sample_map):
        RI_maps = {}
        RI_models = {}
        reverse_RI_models = {}
        RI_list = pd.read_csv(self.parameters['retention_index_standards'])
        for reference_id in tqdm.tqdm(list(dict.fromkeys(list(sample_map.values())))):
            print(reference_id)
            RI_maps[reference_id] = {}
            reference_instance = SimpleSample(self.sample_registry[reference_id], experiment=self)
            prev_index, next_index = None, None
            prev_rt, next_rt = None, None
            RTs, indexes, scan_nos = [], [], []
            for rt, scan_no in zip(reference_instance.list_retention_time, reference_instance.list_scan_numbers):
                RTs.append(rt)
                scan_nos.append(scan_no)
                print(rt, scan_no)
                for index, index_rt in zip(RI_list['Index'], RI_list[reference_instance.name]):
                    index, index_rt = int(index), float(index_rt)
                    print("\t", index, index_rt)
                    if rt > index_rt:
                        prev_index, prev_rt = index, index_rt
                    elif rt <= index_rt:
                        _, next_rt = index, index_rt
                        break
                if next_rt is None:
                    next_rt = max(reference_instance.list_retention_time) * 1.1
                    _ = max(RI_list['Index']) + 1
                RI_value = 100 * (prev_index + ((rt - prev_rt)/(next_rt - rt)))
                indexes.append(RI_value)
                RI_maps[reference_id][rt] = RI_value
            model = lowess(indexes, RTs)
            model2 = lowess(indexes, scan_nos)
            newx, newy = list(zip(*model))
            interf = interpolate.interp1d(newx, newy, fill_value="extrapolate", bounds_error=False)
            RI_models[reference_id] = interf
            newx, newy = list(zip(*model2))
            reverse_RI_models[reference_id] = interpolate.interp1d(newy, newx, fill_value="extrapolate", bounds_error=False)
            
        self.RI_models = RI_models
        self.reverse_RI_models = reverse_RI_models

    def convert_to_RI(self, sample_map):
        if not self.RI_map:
            self.populate_RI_lookup(sample_map)
        for k, v in sample_map.items():
            sam = self.sample_registry[k]
            sam['list_retention_index'] = self.RI_models[v](sam['list_retention_time'])


    def annotate_GC(self):
        pref_ft = os.path.join(self.parameters['outdir'], 'preferred_'+self.parameters['output_feature_table'])
        full_ft = os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table'])
        EI_MS_Library.annotate_gc_feature_table_with_library(pref_ft, self.parameters['GC_Database'])
        EI_MS_Library.annotate_gc_feature_table_with_library(full_ft, self.parameters['GC_Database'])


class EI_MS_Library():
    def __init__(self, library_ID, multicores=None) -> None:
        self.library_meta = self.retrieve_library_meta(library_ID)
        self.library = None
        if multicores is None:
            self.multicores = mp.cpu_count()
        else:
            self.multicores = min(multicores, mp.cpu_count())
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
                if len(library) > limit:
                    break
        else:
            library = list(self.loader(self.library_meta['LIB_PATH']))

        processed_library = []
        for spectrum in library:
            # todo - replace this with matchms pipeline
            # todo - make the pipeline generic for all MS2 in asari ecosystem.
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
                print("Downloading Library...")
                if library_to_load["URL"].endswith('zip'):
                    download_and_unzip_to_pkg_resources(library_to_load['URL'], 'asari', 'db')
                elif library_to_load.get('Compression', None) == 'zip':
                    download_and_unzip_to_pkg_resources(library_to_load['URL'], 'asari', 'db')
                assert os.path.exists(on_disk_path), f"Library not found at {on_disk_path} after Download"
                assert "LIB_PATH" not in library_to_load, "Library already loaded!"
                library_to_load["LIB_PATH"] = on_disk_path
                return library_to_load
            else:
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
        with open(os.path.join(pkg_resources.files('asari'), 'db', 'gcms_libraries.json')) as f:
            return json.load(f)
    
    def annotate_gc_feature_table(self, feature_table_path, drt=0.5, min_peaks=3, min_shared_peaks=1, min_score_threshold=0.7):
        raw_ftgraph = FeatureGraph.ftgraph_from_ft(feature_table_path)
        coelute_ftgraph = raw_ftgraph.filter_graph(drt=drt)
        extracted_spectra = coelute_ftgraph.extract_fragmentation_spectrum(MIN_PEAKS_EXTRACTION=min_peaks, find_clusters=True)
        print(f"Total Cluster Spectra: {len(extracted_spectra)}")
        print(f"Total Library Spectra: {len(self.library)}")
        print(f"Total Comparisons: {len(extracted_spectra) * len(self.library)}, this may take some time...")
        matches = []
        with mp.Pool(self.multicores) as pool:
            scores = pool.imap(wrapped_cosine, product(extracted_spectra, self.library))
            pbar = tqdm.tqdm(scores, desc=f"Comparing Spectra, matches found {len(matches)}", total=len(extracted_spectra) * len(self.library))
            for (extract, library), score in pbar:
                score = score.tolist()
                if len(score) == 2:
                    similarity = score[0]
                    match_peaks = score[1]
                    if similarity >= min_score_threshold and match_peaks >= min_shared_peaks:
                        matches.append({
                            "extract": extract,
                            "library": library,
                            "similarity": similarity,
                            "match_peaks": match_peaks
                        })
                        pbar.set_description_str(f"Comparing Spectra, matches found {len(matches)}")
        coelute_ftgraph.map_annotations(matches)

    @staticmethod
    def annotate_gc_feature_table_with_library(feature_table_path, library_ID, multicores=None):
        library = EI_MS_Library(library_ID, multicores=multicores)
        library.annotate_gc_feature_table(feature_table_path)

def wrapped_cosine(job):
    return job, CosineGreedy().pair(job[0], job[1])


