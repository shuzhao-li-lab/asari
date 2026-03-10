import numpy as np
import tqdm
import networkx as nx
import pandas as pd

from matchms import Spectrum
from matchms.filtering import default_filters, normalize_intensities

class FeatureGraph():
    def __init__(self, ft_path, graph=None):
        if ft_path:
            self.ft_path = ft_path
            self.df = pd.read_csv(self.ft_path, sep='\t')
            self.graph = None
        if graph:
            self.graph = graph
        else:
            self.graph = self.ft_to_graph()
        self.clusters = None
        self.reverse_clusters = None

        assert self.graph is not None, "Either a graph or a feature table path must be provided"
        assert self.df is not None, "Either a graph or a feature table path must be provided"
    
    @staticmethod
    def ftgraph_from_ft(ft_path):
        return FeatureGraph(ft_path)

    @staticmethod
    def ftgraph_from_graph(ft_path, graph):
        return FeatureGraph(ft_path, graph=graph)

    @staticmethod
    def metric(x, y):
        mz = abs(x['mz'] - y['mz'])
        rt = abs(x['rtime'] - y['rtime'])
        return mz, rt

    def ft_to_graph(self):
        G = nx.Graph()
        df_dict = self.df.to_dict(orient='records')
        G.add_nodes_from([(x['id_number'], x) for x in tqdm.tqdm(df_dict, desc="Building Feature Graph Nodes")])
        for x in tqdm.tqdm(df_dict, desc="Building Feature Graph Edges"):
            for y in df_dict:
                if x['id_number'] != y['id_number']:
                    m = FeatureGraph.metric(x, y)
                    G.add_edges_from([(x['id_number'], y['id_number'], {"dmz": m[0], "drt": m[1]})])
        return G
    
    def graph_to_ft(self):
        return pd.DataFrame([d for n, d in self.graph.nodes(data=True)])
    
    def filter_graph(self, drt=.5):
        H = nx.Graph()
        selected_edges = [(u, v) for u, v, d in tqdm.tqdm(self.graph.edges(data=True), desc="Filtering Edges") if d['drt'] < drt]
        H.add_edges_from(selected_edges)
        return FeatureGraph.ftgraph_from_graph(self.ft_path, graph=H)

    def find_spectral_clusters(self):
        assert self.graph is not None, "Must have a graph to find spectral clusters"
        components = nx.connected_components(self.graph)
        feature_to_cluster = {}
        cluster_to_feature = {}
        for clique_id, component_graph in enumerate(components):
            cluster_to_feature[clique_id] = component_graph
            for feature in component_graph:
                self.graph.nodes[feature]['clique_id'] = clique_id
                feature_to_cluster[feature] = clique_id
        print(f"Clustering finds: {len(cluster_to_feature)} spectral clusters!")
        self.clusters = cluster_to_feature
        self.reverse_clusters = feature_to_cluster
    
    def extract_fragmentation_spectrum(self, find_clusters=False, MIN_PEAKS_EXTRACTION=3):
        if find_clusters:
            self.find_spectral_clusters()
        else:
            assert self.clusters is not None, "Must find clusters before extracting fragmentation spectra"
        feature_table = {x['id_number']: x for x in self.df.to_dict(orient='records')}
        processed_spectra = []
        for cluster, feature_list in tqdm.tqdm(self.clusters.items(), desc=f"Extracting Fragmentation Spectra for {len(self.clusters)} Clusters"):
            cluster_spectra = []
            for sample in self.df.columns[11:]:

                cluster_spectrum_data = []
                for feature in feature_list:
                    cluster_spectrum_data.append((feature_table[feature]['mz'], float(feature_table[feature][sample])))
                cluster_mzs = np.array([x[0] for x in sorted(cluster_spectrum_data)])
                cluster_intensities = np.array([x[1] for x in sorted(cluster_spectrum_data)])
                if np.sum(cluster_intensities) > 0:
                    cluster_spectrum = Spectrum(mz=np.array(cluster_mzs), intensities=np.array(cluster_intensities), metadata={"cluster_id": cluster, "sample": sample})
                    cluster_spectra.append(cluster_spectrum)
            if cluster_spectra:
                if len(cluster_spectra) == 1:
                    selected_spectrum = cluster_spectra[0]
                else:
                    spectral_intensities = sorted([(np.sum(x.intensities), x) for x in cluster_spectra], key=lambda x: x[0], reverse=True)
                    for _, selected_spectrum in spectral_intensities:
                        selected_spectrum = default_filters(selected_spectrum)
                        selected_spectrum = normalize_intensities(selected_spectrum)
                        if len(selected_spectrum.peaks) >= MIN_PEAKS_EXTRACTION:
                            processed_spectra.append(selected_spectrum)
                            break
        return processed_spectra

    def map_annotations(self, matches, to_extract=['compound_name', 'inchikey', 'formula']):
        df_dict = {x['id_number']: x for x in self.df.to_dict(orient='records')}

        for match in matches:
            cluster_id = match['extract'].metadata['cluster_id']
            for feature in self.clusters[cluster_id]:
                annotation = []
                for field in to_extract:
                    annotation.append(match['library'].metadata.get(field, 'Not Found'))
                    annotation.append(match['similarity'])
                    annotation.append(match['match_peaks'])
                    if "annotations" not in df_dict[feature]:
                        df_dict[feature]["annotations"] = [annotation]
                    else:
                        df_dict[feature]["annotations"].append(annotation)
                        
        self.df = pd.DataFrame([d for n, d in df_dict.items()])
        self.df.to_csv(self.ft_path.replace(".tsv", "_annotated_gc_beta.tsv"), sep='\t', index=False)
