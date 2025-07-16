import os
import tqdm
import networkx as nx
import pandas as pd 
import numpy as np
from itertools import combinations

import matchms
from matchms import Spectrum
from matchms.similarity import CosineGreedy
from matchms.filtering import default_filters, normalize_intensities
from matchms.importing import load_spectra

# This should be flipped, extractor should have a library or at the very least they should be decoupled.
class EI_MS_Library():
    """
    Handles loading and annotation of EI-MS spectra from a local or remote MoNA library.
    """
    MONA_MSP_URL = (
        "https://mona.fiehnlab.ucdavis.edu/rest/downloads/retrieve/a09f9652-c7bc-48b2-9dd9-0dc4343bb360"
    )

    def __init__(self, parameters) -> None:
        # kind of hackey
        database_path = parameters.get('GC_Database')
        if database_path is None:
            assert os.path.exists(database_path) is True
            raise Exception("Need to provide path to database MSP or directory of MSPs")

        # load required params
        self.library_path = str(parameters.get('GC_Database'))
        self.multicores = parameters.get('multicores', 4)
        self.min_peaks = parameters.get('min_peaks', 1)
        self.min_peaks_common = parameters.get('min_peaks_common', 1)
        self.min_score_threshold = parameters.get('min_score_threshold', 0.70)
        self.similarity_metric = parameters.get('similarity_metric', 'Default (Cosine)')
        self.coelute_threshold = parameters.get('coelute_threshold', 5)
        self.ppm = parameters.get('ppm', 5)
        
        # load the actual library
        self.library = self.load_library()

    
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
    
    def annotate_spectra(self, extracted_spectra):
        print(f"Total Cluster Spectra: {len(extracted_spectra)}")
        print(f"Total Library Spectra: {len(self.library)}")
        print(f"Total Comparisons: {len(extracted_spectra) * len(self.library)}, this may take some time...")

        comparator = CosineGreedy()
        annotations = []
        pbar = tqdm.tqdm(extracted_spectra)
        for extracted_spectrum in pbar:
            for lib in self.library:
                score = comparator.pair(extracted_spectrum, lib).tolist()
                if len(score) == 2:
                    if score[0] >= self.min_score_threshold and score[1] >= self.min_peaks_common:
                        annotations.append({
                            "extract": extracted_spectrum,
                            "library": lib,
                            "similarity": score[0],
                            "match_peaks": score[1],
                            "method": self.similarity_metric
                        })
                        pbar.set_description_str(f"Comparing Spectra, matches found {len(annotations)}")
        return annotations

    def annotate_gc_feature_table(self, feature_table_path):
        extractor = SpectraExtractor(self.coelute_threshold, 
                                     self.ppm, 
                                     feature_table_path, 
                                     min_peaks=self.min_peaks)
        extractor.extract_spectra()
        annotations = self.annotate_spectra(extractor.spectra)
        extractor.apply_annotations(annotations)

class SpectraExtractor():
    # this method extracts the EI_MS spectrum on the basis of co-eluting features
    # eventually extractor should have a library not the other way around.
    def __init__(self, drt, ppm, ft_path, min_peaks=None): 
        self.ft_path = ft_path
        self.drt = drt
        self.ppm = ppm
        if min_peaks is None:
            self.min_peaks = 2
        else:
            self.min_peaks = min_peaks

        self.clusters = None
        self.reverse_clusters = None
        self.spectra = None
    
    def __ft_to_graph(self, ft):
        ft_dict = ft.to_dict(orient='records')
        G = nx.Graph()
        G.add_nodes_from([(x['id_number'], x) for x in tqdm.tqdm(ft_dict, desc="Building Feature Graph Nodes")])
        for x, y in combinations(ft_dict, 2):
            d_rt = abs(x['rtime'] - y['rtime'])
            if d_rt < self.drt:
                d_mz = abs(x['mz'] - y['mz'])
                frag_type = ",".join(sorted(self.classify_edges(x['mz'], y['mz'], self.ppm)))
                G.add_edge(x['id_number'], 
                           y['id_number'], data={"dmz": d_mz, "drt": d_rt, "type": frag_type})
        return G

    def extract_clusters(self):
        df = pd.read_csv(self.ft_path, sep='\t')
        graph = self.__ft_to_graph(df)
        components = nx.connected_components(graph)

        feature_to_cluster = {}
        cluster_to_feature = {}
        for clique_id, component_graph in enumerate(components):
            cluster_to_feature[clique_id] = component_graph
            for feature in component_graph:
                graph.nodes[feature]['clique_id'] = clique_id
                feature_to_cluster[feature] = clique_id
        print(f"Clustering finds: {len(cluster_to_feature)} spectral clusters!")
        self.clusters = cluster_to_feature
        self.reverse_clusters = feature_to_cluster
                
    def extract_spectra(self):
        ft = pd.read_csv(self.ft_path, sep="\t")
        ft_dict = {x['id_number']: x for x in ft.to_dict(orient='records')}
        self.extract_clusters()
        
        processed_spectra = []
        for cluster, feature_list in tqdm.tqdm(self.clusters.items(), desc=f"Extracting Fragmentation Spectra for {len(self.clusters)} Clusters"):
            cluster_spectra = []
            for sample in ft.columns[11:]:
                cluster_spectrum_data = []
                for feature in feature_list:
                    cluster_spectrum_data.append((ft_dict[feature]['mz'], float(ft_dict[feature][sample])))
                sorted_cluster_data = sorted(cluster_spectrum_data)
                cluster_mzs = np.array([x[0] for x in sorted_cluster_data])
                cluster_intensities = np.array([x[1] for x in sorted_cluster_data])
                if np.sum(cluster_intensities) > 0:
                    cluster_spectrum = Spectrum(mz=np.array(cluster_mzs), intensities=np.array(cluster_intensities), metadata={"cluster_id": cluster, "sample": sample})
                    cluster_spectra.append(cluster_spectrum)
            if cluster_spectra:
                spectral_intensities = sorted([(np.sum(x.intensities), x) for x in cluster_spectra], key=lambda x: x[0], reverse=True)
                for _, selected_spectrum in spectral_intensities:
                    selected_spectrum = default_filters(selected_spectrum)
                    selected_spectrum = normalize_intensities(selected_spectrum)
                    if len(selected_spectrum.peaks) >= self.min_peaks:
                        processed_spectra.append(selected_spectrum)
        self.spectra = processed_spectra
    
    def apply_annotations(self, annotations, to_extract=['compound_name', 'inchikey', 'formula']):
        ft_dict = {x['id_number']: x for x in pd.read_csv(self.ft_path, sep="\t").to_dict(orient='records')}
        for match in annotations:
            cluster_id = match['extract'].metadata['cluster_id']
            for feature in self.clusters[cluster_id]:
                annotation = []
                for field in to_extract:
                    annotation.append(match['library'].metadata.get(field, 'Not Found'))
                    annotation.append(match['similarity'])
                    annotation.append(match['match_peaks'])
                    if "annotations" not in ft_dict[feature]:
                        ft_dict[feature]["annotations"] = [annotation]
                    else:
                        ft_dict[feature]["annotations"].append(annotation)
        self.df = pd.DataFrame([d for n, d in ft_dict.items()])
        self.df.to_csv(self.ft_path.replace(".tsv", "_annotated_gc_beta.tsv"), sep='\t', index=False)

    @staticmethod
    def classify_edges(m1, m2, ppm) -> list[str]:
        """
        Return all edge‐type names whose expected mass delta lies within ±tol of the observed delta_mz.

        Parameters
        ----------
        delta_mz : float
            Observed mass difference.
        edge_dict : dict, optional
            Mapping of edge‐type names to their theoretical mass deltas.
        tol : float, optional
            Mass tolerance in Daltons for a match.

        Returns
        -------
        list[str]
            List of matching edge‐type names. Empty if none match.
        """

        # this should be moved to somewhere else like mass2chem
        EDGE_TYPES = {
            # single atoms
            "H":      1.007825,    # hydrogen radical
            "H2":     2.015650,    # dihydrogen
            "C":     12.000000,    # carbon atom
            "N":     14.003074,    # nitrogen atom
            "O":     15.994915,    # oxygen atom

            # small radicals/groups
            "CH":    13.018860,    # methylidyne
            "CH2":   14.015650,    # methylene
            "CH3":   15.023475,    # methyl
            "NH":    15.010899,    # imidogen
            "NH2":   16.018724,    # amidogen
            "NH3":   17.026549,    # ammonia
            "OH":    17.002740,    # hydroxyl

            # common neutral losses
            "H2O":   18.010565,    # water
            "CO":    27.994915,    # carbon monoxide
            "CO2":   43.989830,    # carbon dioxide
        }

        delta_mz = m1 - m2
        allowed_err_m1 = m1 / 1e6 * ppm
        allowed_err_m2 = m2 / 1e6 * ppm
        total_err = allowed_err_m1 + allowed_err_m2
        matches = [
            name
            for name, expected_delta in EDGE_TYPES.items()
            if abs(delta_mz - expected_delta) < total_err
        ]

        return matches
