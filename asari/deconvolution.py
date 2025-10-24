# -*- coding: utf-8 -*-
"""
Asari Processing Module

This script consolidates the functionality from AutoKhipu (feature deconvolution),
MSMS Annotate (GC-MS spectral matching), and L4 Annotate (LC-MS formula inference)
into a single, configurable module for metabolomics data processing.

The main entry point is the `AsariProcessor` class, which orchestrates the
workflow based on the specified mode ('gc' or 'lc').

Required Libraries:
- pandas
- numpy
- tqdm
- intervaltree (pip install intervaltree)
- mass2chem (pip install mass2chem)
- matchms (pip install matchms)
- matplotlib (pip install matplotlib)
"""
import json
import os
import pandas as pd
import heapq
import numpy as np
import random
import argparse
import multiprocessing as mp
import tqdm
import shlex
import subprocess
from math import comb
from itertools import combinations, product
from collections import defaultdict

# --- Third-party library imports with installation notes ---
try:
    from mass2chem.formula import calculate_formula_mass
except ImportError:
    print("Warning: 'mass2chem' not found. Please install via 'pip install mass2chem'.")
    # Define a dummy function to avoid crashing if mass2chem is not installed
    def calculate_formula_mass(formula):
        print(f"Warning: mass2chem not installed. Cannot calculate mass for {formula}.")
        return 0.0

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages  # <-- ADD THIS
except ImportError:
    print("Warning: 'matplotlib' not found. Please install via 'pip install matplotlib'.")
    plt = None
    PdfPages = None  # <-- ADD THIS

try:
    from intervaltree import Interval, IntervalTree
except ImportError:
    print("Warning: 'intervaltree' not found. Please install via 'pip install intervaltree'.")
    Interval, IntervalTree = None, None

try:
    from matchms.importing import load_from_msp
    from matchms import Spectrum
    from matchms.filtering import default_filters, normalize_intensities, add_retention_index
    from matchms.similarity import CosineGreedy
except ImportError:
    print("Warning: 'matchms' not found. Please install via 'pip install matchms'.")
    load_from_msp, Spectrum, default_filters, normalize_intensities, CosineGreedy = [None] * 5

try:
    from pypdf import PdfReader, PdfWriter
except ImportError:
    print("Warning: 'pypdf' not found. Please install via 'pip install pypdf'. PDF bookmarks will not be generated.")
    PdfReader, PdfWriter = None, None


# ==============================================================================
# ==                         STATISTICAL & UTILITY FUNCTIONS                  ==
# ==============================================================================

def _hg_tail_obs(k, N, K, n):
    """Calculates the tail probability of the hypergeometric distribution."""
    if comb(N, n) == 0: return 1.0
    top = 0
    for i in range(k, min(n, K) + 1):
        if n - i > N - K:
            continue
        top += comb(K, i) * comb(N - K, n - i)
    return top / comb(N, n)

def _pearson(x, y):
    """Calculates Pearson correlation, handling small inputs."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 2: return np.nan
    x -= x.mean()
    y -= y.mean()
    den = np.linalg.norm(x) * np.linalg.norm(y)
    return float(np.dot(x, y) / den) if den > 0 else np.nan

def _sigma_eff(apex, left, right):
    """Calculates the effective standard deviation of an asymmetric peak."""
    sL = max((apex - left) / 3.0, 1e-12)
    sR = max((right - apex) / 3.0, 1e-12)
    return ((sL * sL + sR * sR) / 2.0) ** 0.5

def bc_overlap(apex1, left1, right1, apex2, left2, right2):
    """Calculates the Bhattacharyya coefficient for two Gaussian peaks."""
    s1 = _sigma_eff(apex1, left1, right1)
    s2 = _sigma_eff(apex2, left2, right2)
    ssum = s1 * s1 + s2 * s2
    if ssum == 0: return 0.0
    pref = np.sqrt((2.0 * s1 * s2) / ssum)
    expo = np.exp(-((apex1 - apex2) ** 2) / (2.0 * ssum))
    return float(pref * expo)


# ==============================================================================
# ==                       DECONVOLUTION FRAMEWORK (AUTOKHIPU)                ==
# ==============================================================================

class DeconvolutionFramework:
    """
    Performs hierarchical deconvolution of metabolomics features to group them
    into isotopic clusters and then into putative compounds based on neutral
    losses and statistical correlations.
    """
    def __init__(self, feature_dict, samples, params, default_abs_mz_tolerance=0.0002):
        self.feature_dict = feature_dict
        self.samples = samples
        self.params = params
        self.default_abs_mz_tolerance = default_abs_mz_tolerance

        self.neutral_loss_data = self._load_neutral_losses()
        self.neutral_loss_tree = self._build_neutral_loss_tree()

        self.comparison_cache = {}
        self.statistical_baselines = {}
        self.isotopic_clusters = {}
        self.final_compounds = {}
        self.ALPHA = 0.05

        self.ISOTOPE_DELTAS = {
            "13C": 1.00335483507, "15N": 0.99703489445, "37Cl": 1.99704992,
            "34S": 1.995795826, "2H": 1.00627674589, "18O": 2.004245,
            "29Si": 0.99956813, "81Br": 1.9979521, "78Se": -1.99921256
        }

    # --- Configuration Properties ---
    @property
    def neutral_loss_file(self):
        return self.params.get("neutral_loss_file")

    @property
    def neutral_loss_tolerance(self):
        # Renamed from nl_tolerance in params for clarity
        return self.params.get("nl_tolerance", self.default_abs_mz_tolerance)

    @property
    def isotope_delta_tolerance(self):
        return self.params.get("isotope_mz_tolerance", self.default_abs_mz_tolerance)

    @property
    def rt_window(self):
        return self.params.get("rt_window", 1.0)

    @property
    def correlation_threshold_bc(self):
        return self.params.get("correlation_threshold_bc", 0.8)

    def _load_neutral_losses(self):
        """Loads neutral loss data from the file specified in parameters."""
        if self.neutral_loss_file and os.path.exists(self.neutral_loss_file):
            print(f"--- Loading neutral loss file: {self.neutral_loss_file} ---")
            try:
                with open(self.neutral_loss_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error reading neutral loss file: {e}")
        return None

    def _build_neutral_loss_tree(self):
        """Builds an interval tree for efficient neutral loss lookup."""
        if not self.neutral_loss_data or not IntervalTree:
            return None
        print("Building neutral loss interval tree...")
        tree = IntervalTree()
        for loss in self.neutral_loss_data:
            mass = calculate_formula_mass(loss['neutral_loss_formula'])
            loss['mass'] = mass
            tree.add(Interval(mass - self.neutral_loss_tolerance, mass + self.neutral_loss_tolerance, loss))
        print(f"Neutral loss tree built with {len(tree)} intervals.")
        return tree

    def _compare_features(self, fid_a, fid_b):
        """Compares two features statistically and chromatographically."""
        key = tuple(sorted((fid_a, fid_b)))
        if key in self.comparison_cache: return self.comparison_cache[key]

        fa, fb = self.feature_dict[fid_a], self.feature_dict[fid_b]
        va = np.array([fa.get(s, np.nan) for s in self.samples], float)
        vb = np.array([fb.get(s, np.nan) for s in self.samples], float)

        ok = np.isfinite(va) & np.isfinite(vb)
        posA, posB = ok & (va > 0), ok & (vb > 0)
        setA = {s for s, m in zip(self.samples, posA) if m}
        setB = {s for s, m in zip(self.samples, posB) if m}
        N, n, m, k = len(self.samples), len(setA), len(setB), len(setA & setB)
        m_lin = posA & posB

        ratio_mean, ratio_std = np.nan, np.nan
        if m_lin.sum() > 1:
            ratios = vb[m_lin] / (va[m_lin] + 1e-9)
            ratio_mean = np.mean(ratios)
            ratio_std = np.std(ratios)
        
        result = {
            "p_BinA": 1.0 if k < m or N * n * m == 0 else _hg_tail_obs(k, N, n, m),
            "r_pearson": _pearson(va[m_lin], vb[m_lin]) if m_lin.sum() > 1 else np.nan,
            "bc": bc_overlap(
                fa['rtime'], fa['rtime_left_base'], fa['rtime_right_base'],
                fb['rtime'], fb['rtime_left_base'], fb['rtime_right_base']
            ),
            "intensity_ratio_mean": ratio_mean,
            "intensity_ratio_std": ratio_std,
        }
        self.comparison_cache[key] = result
        return result

    def _calculate_statistical_baselines(self, num_samples=1000):
        """Calculates null distribution for Pearson correlation."""
        print(f"Calculating statistical baselines from {num_samples} random feature pairs...")
        feature_ids = list(self.feature_dict.keys())
        if len(feature_ids) < 2:
            self.statistical_baselines['r_pearson'] = {'median': 0.1, 'std': 0.1}
            return

        null_dist = [self._compare_features(*random.sample(feature_ids, 2)) for _ in tqdm.tqdm(range(num_samples))]
        values = [res['r_pearson'] for res in null_dist if res and not np.isnan(res['r_pearson'])]
        self.statistical_baselines['r_pearson'] = {
            'median': np.nanmedian(values) if values else 0.1,
            'std': np.nanstd(values) if values else 0.1
        }
        print("Baselines calculated:", self.statistical_baselines)

    def _find_isotopic_envelope(self, anchor_fid, rt_window_fids, used_fids):
        """Finds isotopic peaks related to an anchor feature."""
        anchor_data = self.feature_dict[anchor_fid]
        envelope = [{'id': anchor_fid, 'type': 'anchor', 'mz': anchor_data['mz'], 'rtime': anchor_data['rtime'], 'evidence': None}]
        found_isotopes = {anchor_fid}
        base_isotopes_found = []

        # --- Forward and Reverse Search (M+n, M-n) ---
        for name, mass_delta in self.ISOTOPE_DELTAS.items():
            for direction in [1, -1]: # 1 for forward, -1 for reverse
                last_found_fid = anchor_fid
                for num_iso in range(1, 6):
                    expected_mz = anchor_data['mz'] + (direction * mass_delta * num_iso)
                    best_candidate = None
                    for fid in rt_window_fids:
                        if fid not in used_fids and fid not in found_isotopes:
                            candidate_mz = self.feature_dict[fid]['mz']
                            if abs(candidate_mz - expected_mz) < self.isotope_delta_tolerance:
                                stats = self._compare_features(last_found_fid, fid)
                                if stats['bc'] > 0.9 and stats.get('r_pearson', 0) > self.statistical_baselines['r_pearson']['median']:
                                    best_candidate = fid
                                    break
                    if best_candidate:
                        stats_vs_anchor = self._compare_features(anchor_fid, best_candidate)
                        candidate_data = self.feature_dict[best_candidate]
                        iso_type = f"{'+' if direction == 1 else '-'}{name}_{num_iso}"
                        envelope.append({
                            'id': best_candidate, 'type': 'isotope', 'mz': candidate_data['mz'],
                            'rtime': candidate_data['rtime'], 'evidence': {'linked_to': anchor_fid, 'isotope_type': iso_type}
                        })
                        if direction == 1 and num_iso == 1:
                            base_isotopes_found.append((name, mass_delta, best_candidate))
                        found_isotopes.add(best_candidate)
                        last_found_fid = best_candidate
                    else:
                        break
        return envelope

    def _build_isotopic_clusters(self):
        """Groups features into clusters based on isotopic patterns."""
        print("\n--- Phase 1: Building Isotopic Clusters ---")
        feature_heap = [(-data.get('detection_counts', 0), fid) for fid, data in self.feature_dict.items()]
        heapq.heapify(feature_heap)
        used_fids = set()
        rt_sorted_fids = sorted(self.feature_dict.keys(), key=lambda k: self.feature_dict[k]['rtime'])

        with tqdm.tqdm(total=len(feature_heap)) as pbar:
            while feature_heap:
                _, anchor_fid = heapq.heappop(feature_heap)
                if anchor_fid in used_fids:
                    pbar.update(1)
                    continue

                anchor_rt = self.feature_dict[anchor_fid]['rtime']
                rt_window_fids = [fid for fid in rt_sorted_fids if abs(self.feature_dict[fid]['rtime'] - anchor_rt) < self.rt_window]
                
                annotated_envelope = self._find_isotopic_envelope(anchor_fid, rt_window_fids, used_fids)
                envelope_fids = {f['id'] for f in annotated_envelope}
                
                pbar.update(len(envelope_fids - used_fids))
                
                mono_isotopic_feature = min(annotated_envelope, key=lambda f: f['mz'])
                self.isotopic_clusters[mono_isotopic_feature['id']] = annotated_envelope
                used_fids.update(envelope_fids)

        print(f"Found {len(self.isotopic_clusters)} isotopic clusters.")

    def _find_neutral_loss(self, mz1, mz2):
        """Finds matching neutral losses using the pre-built interval tree."""
        if not self.neutral_loss_tree: return []
        mass_diff = abs(mz1 - mz2)
        matches = []
        for interval in self.neutral_loss_tree.at(mass_diff):
            loss = interval.data
            matches.append({
                "neutral_loss_formula": loss['neutral_loss_formula'],
                "mass_diff": mass_diff,
                "neutral_loss_mass": loss['mass'],
                "mass_error_da": mass_diff - loss['mass']
            })
        return matches

    def _merge_fragment_groups(self):
        """Merges isotopic clusters into compounds."""
        cluster_ids = list(self.isotopic_clusters.keys())
        parent = {cid: cid for cid in cluster_ids}
        
        def find(v):
            if v == parent[v]: return v
            parent[v] = find(parent[v]); return parent[v]
        
        def unite(a, b):
            a, b = find(a), find(b)
            if a != b: parent[b] = a

        # --- Step 1: Targeted merging via neutral losses ---
        if self.neutral_loss_tree:
            print("\n--- Phase 2a: Targeted merging via neutral losses ---")
            for i in tqdm.tqdm(range(len(cluster_ids)), desc="Targeted Merge"):
                for j in range(i + 1, len(cluster_ids)):
                    c1_id, c2_id = cluster_ids[i], cluster_ids[j]
                    if abs(self.feature_dict[c1_id]['rtime'] - self.feature_dict[c2_id]['rtime']) > self.rt_window:
                        continue
                    if self._find_neutral_loss(self.feature_dict[c1_id]['mz'], self.feature_dict[c2_id]['mz']):
                        stats = self._compare_features(c1_id, c2_id)
                        if stats['bc'] > self.correlation_threshold_bc and stats.get('r_pearson', 0) > self.statistical_baselines['r_pearson']['median']:
                            unite(c1_id, c2_id)

        # --- Step 2: Salvage merging via statistical correlation ---
        print("\n--- Phase 2b: Salvage merging via statistical correlation ---")
        for i in tqdm.tqdm(range(len(cluster_ids)), desc="Salvage Merge"):
            for j in range(i + 1, len(cluster_ids)):
                c1_id, c2_id = cluster_ids[i], cluster_ids[j]
                if find(c1_id) == find(c2_id): continue
                if abs(self.feature_dict[c1_id]['rtime'] - self.feature_dict[c2_id]['rtime']) < self.rt_window:
                    stats = self._compare_features(c1_id, c2_id)
                    if stats['bc'] > self.correlation_threshold_bc and stats.get('r_pearson', 0) > self.statistical_baselines['r_pearson']['median']:
                        unite(c1_id, c2_id)
        
        # --- Assemble final compounds ---
        final_compounds_map = defaultdict(list)
        for cid in cluster_ids:
            final_compounds_map[find(cid)].append(cid)
        
        print(f"\n--- Assembling {len(final_compounds_map)} final compounds ---")
        for root_id, member_ids in final_compounds_map.items():
            compound = {"anchor_fragment": {}, "fragments": [], "annotations": {}}
            for member_id in sorted(member_ids, key=lambda mid: self.feature_dict[mid]['mz']):
                cluster_data = self.isotopic_clusters[member_id]
                anchor = next(f for f in cluster_data if f['type'] == 'anchor')
                isotopes = [f for f in cluster_data if f['type'] == 'isotope']
                
                fragment_data = {
                    "id": anchor['id'], "mz": anchor['mz'], "rtime": anchor['rtime'],
                    "isotopes": isotopes, "neutral_loss_annotations": []
                }

                if member_id == root_id:
                    compound['anchor_fragment'] = fragment_data
                else:
                    nl_annotations = self._find_neutral_loss(self.feature_dict[root_id]['mz'], anchor['mz'])
                    fragment_data['neutral_loss_annotations'] = nl_annotations
                    compound['fragments'].append(fragment_data)
            self.final_compounds[root_id] = compound
        print(f"Merged into {len(self.final_compounds)} final compounds.")

    def run(self):
        """Executes the full deconvolution workflow."""
        self._calculate_statistical_baselines()
        self._build_isotopic_clusters()
        self._merge_fragment_groups()
        return self.final_compounds

# ==============================================================================
# ==                      GC-MS ANNOTATOR (MSMS ANNOTATE)                     ==
# ==============================================================================

class MSMSAnnotator:
    """
    Annotates compounds by matching their reconstructed spectra against
    user-provided spectral libraries.
    """
    def __init__(self, feature_table_dict, sample_cols, params):
        if not all([load_from_msp, Spectrum, CosineGreedy]):
            raise ImportError("MatchMS library is not installed. Cannot perform GC-MS annotation.")
        
        self.ft_dict = feature_table_dict
        self.sample_cols = sample_cols
        self.params = params
        self.lib_spectra = self._load_libraries()
        self.cosine_greedy = CosineGreedy(tolerance=self.msms_tolerance)

    # --- Configuration Properties ---
    @property
    def msms_tolerance(self):
        return self.params.get("msms_tolerance", 0.005)

    @property
    def spectral_libraries(self):
        libs = self.params.get("GC_Database", [])
        if isinstance(libs, str):
            path = os.path.abspath(libs)
            if os.path.isdir(path):
                return [
                    os.path.join(path, f)
                    for f in os.listdir(path)
                    if f.lower().endswith(".msp")
                ]
            return [path]
        return libs

    @property
    def min_matched_peaks(self):
        return self.params.get("min_matched_peaks", 1)

    @property
    def min_cosine_score(self):
        return self.params.get("min_cosine_score", 0.7)

    @property
    def mirror_plot_pdf(self):
        """Defines the output path for the PDF plot report."""
        prefix = self.params.get("output_prefix", "asari_results")
        outdir = self.params.get("outdir", ".") 
        pdf_path = os.path.join(outdir, f"{prefix}_mirror_plots.pdf")
        return os.path.abspath(pdf_path)

    @property
    def generate_plots(self):
        """Flag to enable/disable plot generation."""
        # We need both matplotlib and PdfPages to be available
        return plt is not None and PdfPages is not None and self.mirror_plot_pdf


    def _load_libraries(self):
        """Loads and preprocesses spectral libraries."""
        lib_spectra = []
        print("\n--- Loading spectral libraries ---")
        for lib_path in self.spectral_libraries:
            if os.path.exists(lib_path):
                print(f"Loading {lib_path}...")
                try:
                    for s in load_from_msp(lib_path):
                        lib_spectra.append((os.path.basename(lib_path), self._prep_spectrum(s)))
                except Exception as e:
                    print(f"Could not load library {lib_path}: {e}")
            else:
                print(f"Warning: Library file not found: {lib_path}")
        print(f"Loaded {len(lib_spectra)} total library spectra.")
        return lib_spectra

    def _prep_spectrum(self, spec):
        """Applies standard filters to a spectrum."""
        spec = default_filters(spec)
        spec = normalize_intensities(spec)
        add_retention_index(spec)
        return spec
    
    def _plot_mirror(self, query, lib, pdf_object, title=None):
        """Generates a mirror plot and saves it directly to the PDF object."""
        if not plt or not pdf_object: return
        
        # Create a new figure and axes for this plot
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)

        q_mz, q_ints = query.peaks.mz, query.peaks.intensities
        l_mz, l_ints = lib.peaks.mz, lib.peaks.intensities
        q_ints = q_ints / q_ints.max() if q_ints.max() > 0 else q_ints
        l_ints = l_ints / l_ints.max() if l_ints.max() > 0 else l_ints

        ax.vlines(q_mz, 0, q_ints, color="blue", lw=1, label="Query")
        ax.vlines(l_mz, -l_ints, 0, color="red", lw=1, label="Library")
        ax.axhline(0, color='black', lw=0.5)
        ax.set_xlabel("m/z")
        ax.set_ylabel("Relative Intensity")
        ax.set_title(title or "Mirror Plot")
        ax.legend()
        plt.tight_layout()
        
        # Save the current figure to the PDF and close it to free memory
        pdf_object.savefig(fig)
        plt.close(fig)

    def _get_feature_row(self, fid):
        """Safely retrieves a feature row from the dictionary."""
        fid_str = str(fid)
        fid_int_str = str(int(float(fid))) if '.' in fid_str else fid_str
        return self.ft_dict.get(fid_str) or self.ft_dict.get(fid_int_str)

    def _spectrum_for_compound(self, cpd):
        """Constructs a pseudo-MS/MS spectrum for a compound."""
        members = [cpd["anchor_fragment"]["id"]]
        members.extend(iso["id"] for iso in cpd["anchor_fragment"].get("isotopes", []))
        for frag in cpd.get("fragments", []):
            members.append(frag["id"])
            members.extend(iso["id"] for iso in frag.get("isotopes", []))

        anchor_row = self._get_feature_row(cpd["anchor_fragment"]["id"])
        if not anchor_row: return None, None

        svals = pd.Series({c: pd.to_numeric(anchor_row.get(c, 0), errors="coerce") for c in self.sample_cols}).fillna(0)
        rep_sample = svals.idxmax() if not svals.empty else self.sample_cols[0]

        mzs, ints = [], []
        for mid in members:
            row = self._get_feature_row(mid)
            if row:
                mz = pd.to_numeric(row.get("mz", np.nan), errors="coerce")
                intensity = pd.to_numeric(row.get(rep_sample, 0), errors="coerce")
                if np.isfinite(mz) and intensity > 0:
                    mzs.append(float(mz))
                    ints.append(float(intensity))

        if not mzs: return None, None

        order = np.argsort(mzs)
        mzs, ints = np.array(mzs)[order], np.array(ints)[order]
        precursor = float(anchor_row.get("mz", np.nan))
        
        spec = Spectrum(mz=mzs, intensities=ints, metadata={"precursor_mz": precursor})
        return self._prep_spectrum(spec), rep_sample

    def _score_compound(self, fid_cpd_tuple):
        """
        Scores a compound against all library spectra.
        Returns annotation data and data required for plotting.
        """
        fid, cpd = fid_cpd_tuple
        qspec, rep_sample = self._spectrum_for_compound(cpd)
        if qspec is None: return fid, None, []  # Return empty list for plots

        scores = []
        plot_data = [] # List to hold plotting info

        for libname, lspec in self.lib_spectra:
            score_tuple = self.cosine_greedy.pair(qspec, lspec)
            score, n_matches = score_tuple['score'], score_tuple['matches']
            retention_index = lspec.get("retention_index", "Not Reported")

            if n_matches >= self.min_matched_peaks and score > self.min_cosine_score:
                lib_compound_name = lspec.get("compound_name") or lspec.get("name") or "unknown"
                
                scores.append({
                    "library": libname, 
                    "name": lib_compound_name,
                    "precursor_mz": lspec.get("precursor_mz"),
                    "score": float(score), 
                    "matched_peaks": int(n_matches),
                    "retention_index": retention_index
                })
                
                # If plotting is enabled, store the data needed to make the plot
                if self.generate_plots:
                    title = (
                        f"Compound {fid} (Query) vs {lib_compound_name} (Library)\n"
                        f"(Score: {score:.2f}, Matched Peaks: {n_matches}, "
                        f"Rep Sample: {rep_sample})"
                    )
                    plot_data.append({
                        "fid": fid,
                        "query_spec": qspec,
                        "lib_spec": lspec,
                        "title": title
                    })

        scores.sort(key=lambda x: (-x["score"], -x["matched_peaks"]))
        
        return fid, {
            "representative_sample": rep_sample,
            "top_matches": scores
        }, plot_data  # Return plot data

    def annotate_compounds(self, compounds, num_cores=4):
        """
        Runs the annotation process for all compounds, separating
        scoring (parallel) from PDF plotting (serial).
        """
        print("\n--- Phase 3 (GC): Annotating compounds via MS/MS search ---")
        if not self.lib_spectra:
            print("No libraries loaded, skipping MS/MS annotation.")
            return compounds

        pdf_object = None
        all_plot_data = []
        
        # 1. Initialize PDF object in the main process
        if self.generate_plots:
            pdf_path = self.mirror_plot_pdf
            # Ensure output directory exists
            os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
            print(f"--- Initializing PDF report at: {pdf_path} ---")
            pdf_object = PdfPages(pdf_path)
        else:
            print("Skipping PDF plot generation (matplotlib/PdfPages not found or no output path set).")

        # 2. Run multiprocessing pool to get scores and plot *data*
        with mp.Pool(num_cores) as pool:
            results = list(tqdm.tqdm(
                pool.imap(self._score_compound, compounds.items()),
                total=len(compounds),
                desc="GC-MS Annotation"
            ))

        # 3. Process results in the main process
        print("Processing results and annotations...")
        for fid, annotation, plot_data in results:
            if annotation:
                compounds[fid]['annotations']['gc_msms_search'] = annotation
            if plot_data:
                # Store plot data
                all_plot_data.extend(plot_data)
        
        # 4. Generate plots (in the main process)
        if pdf_object:
            print(f"--- Generating {len(all_plot_data)} mirror plots for PDF ---")
            
            # Sort plots by feature ID (fid) to "index" the report
            all_plot_data.sort(key=lambda x: (
                # Try to sort numerically by extracting digits from the ID
                int(''.join(filter(str.isdigit, str(x["fid"]))) or 0), 
                str(x["fid"]) # Fallback to string sort
            ))
            
            for plot_item in tqdm.tqdm(all_plot_data, desc="Saving plots to PDF"):
                self._plot_mirror(
                    plot_item["query_spec"],
                    plot_item["lib_spec"],
                    pdf_object,
                    plot_item["title"]
                )
            
            # 5. Close the PDF object
            pdf_object.close()
            print(f"--- PDF report generation complete. ---")
        
        return compounds

# ==============================================================================
# ==                    LC-MS ANNOTATOR (L4 FORMULA ANNOTATE)                 ==
# ==============================================================================

class FormulaAnnotator:
    """
    Infers possible chemical formulas for compounds using an external tool.
    NOTE: This class requires a pre-configured external command-line tool.
    """
    def __init__(self, params):
        self.params = params

    # --- Configuration Properties ---
    @property
    def formula_predictor_cmd(self):
        return self.params.get("formula_predictor_cmd")

    @property
    def temp_csv_dir(self):
        return self.params.get("temp_csv_dir", "./temp_l4")

    @property
    def formula_mass_tolerance(self):
        return self.params.get("formula_mass_tolerance", 0.001)

    def _run_external_predictor(self, input_csv_path, output_csv_path):
        """Executes the external formula prediction command using subprocess."""
        if not self.formula_predictor_cmd:
            print("Warning: 'formula_predictor_cmd' not defined in parameters. Skipping formula prediction.")
            return False
        
        command = self.formula_predictor_cmd.format(
            input_csv=shlex.quote(input_csv_path),
            output_csv=shlex.quote(output_csv_path)
        )
        print(f"Running external command: {command}")
        try:
            # Execute the command, raise exception on error, capture output
            result = subprocess.run(
                command,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            if result.stdout:
                print("Go program stdout:", result.stdout)
            if result.stderr:
                print("Go program stderr:", result.stderr)
            return True
        except FileNotFoundError:
            print(f"Error: Command not found. Make sure the Go executable is in your system's PATH or provide a full path.")
            return False
        except subprocess.CalledProcessError as e:
            # This block runs if the Go program returns a non-zero exit code
            print(f"An error occurred while running the external command.")
            print(f"Exit Code: {e.returncode}")
            print(f"Stdout: {e.stdout}")
            print(f"Stderr: {e.stderr}")
            return False

    def annotate_compounds(self, compounds):
        """Runs the formula annotation workflow."""
        print("\n--- Phase 3 (LC): Annotating compounds with possible formulas ---")
        masses_to_predict = []
        for fid, cpd in compounds.items():
            # Assuming an 'inferred_mass' key is added by a previous (unprovided) step.
            # For now, we'll use the anchor fragment's mz as a proxy.
            if 'inferred_mass' not in cpd:
                cpd['inferred_mass'] = {'calculated_monoisotopic_mz': cpd['anchor_fragment']['mz']}
            
            masses_to_predict.append({
                "Mass": cpd['inferred_mass']['calculated_monoisotopic_mz']
            })

        if not masses_to_predict:
            print("No masses found for formula prediction.")
            return compounds

        # Prepare CSV files for the external tool
        os.makedirs(self.temp_csv_dir, exist_ok=True)
        input_csv = os.path.join(self.temp_csv_dir, "l4_input_masses.csv")
        output_csv = os.path.join(self.temp_csv_dir, "l4_output_formulas.csv")
        pd.DataFrame(masses_to_predict).drop_duplicates().to_csv(input_csv, index=False)
        
        # Run predictor and process results
        if self._run_external_predictor(input_csv, output_csv) and os.path.exists(output_csv):
            possible_formulas = {}
            results_df = pd.read_csv(output_csv)
            
            for _, row in results_df.iterrows():
                mass = row['Mass']
                if mass not in possible_formulas:
                    possible_formulas[mass] = []
                if isinstance(row.get('FormulaList'), str):
                    for f_prob in row['FormulaList'].split(';'):
                        formula = f_prob.split('_')[0]
                        if abs(calculate_formula_mass(formula) - float(mass)) < self.formula_mass_tolerance:
                            possible_formulas[mass].append(formula)
            
            # Add predictions back to the compounds
            for cpd in compounds.values():
                mass = cpd['inferred_mass']['calculated_monoisotopic_mz']
                if mass in possible_formulas:
                    cpd['annotations']['lc_formula_predictions'] = possible_formulas[mass]
        
        return compounds

# ==============================================================================
# ==                          MAIN WORKFLOW ORCHESTRATOR                      ==
# ==============================================================================

class AsariProcessor:
    """
    Orchestrates the entire metabolomics data processing workflow, from
    feature loading to deconvolution and mode-specific annotation.
    """
    def __init__(self, params):
        self.params = params
        self.feature_table = None
        self.feature_dict = None
        self.samples = []
        self.deconvoluted_compounds = {}
        self.output_path = None

    # --- Configuration Properties ---
    @property
    def mode(self):
        return self.params.get("workflow", "gc").lower()

    @property
    def id_column(self):
        return self.params.get("id_column", "id_number")

    @property
    def sample_start_col(self):
        return self.params.get("sample_start_col", 11)
        
    @property
    def output_prefix(self):
        return self.params.get("output_prefix", "asari_results")

    @property
    def num_cores(self):
        return self.params.get("num_cores", 4)

    # Properties for sub-dictionaries to pass to other classes
    @property
    def deconvolution_params(self):
        return self.params.get("deconvolution", {})

    @property
    def gc_annotation_params(self):
        return self.params.get("gc_annotation", {})

    @property
    def lc_annotation_params(self):
        return self.params.get("lc_annotation", {})

    def default_workflow(self, filepath):
        """Runs the full, standard workflow on a given feature table file."""
        final_data = {}
        print(filepath)
        if self.output_path is None:
            self.output_path = os.path.join(self.parameters['outdir'], '/export/')
        if self.load_feature_table(filepath):
            if self.run_deconvolution():
                self.run_annotation()
                self.save_results()
                final_data = self.get_results()
                print(f"\nProcessing complete. Found {len(final_data)} compounds.")
        return final_data

    def load_feature_table(self, filepath):
        """Loads and prepares the feature table from a CSV or TSV file."""
        print(f"--- Loading feature table: {filepath} ---")
        try:
            sep = '\t' if filepath.endswith('.tsv') else ','
            self.feature_table = pd.read_csv(filepath, sep=sep)
            print(f"Loaded table with shape: {self.feature_table.shape}")

            if self.id_column not in self.feature_table.columns:
                raise ValueError(f"ID column '{self.id_column}' not found in the feature table.")

            self.samples = list(self.feature_table.columns[self.sample_start_col:])
            print(f"Found {len(self.samples)} samples.")
            
            self.feature_table[self.id_column] = self.feature_table[self.id_column].astype(str)
            self.feature_table.fillna(0.0, inplace=True)
            self.feature_dict = {row[self.id_column]: row for row in self.feature_table.to_dict(orient='records')}
            print(f"Loaded {len(self.feature_dict)} features.")
            return True
        
        except Exception as e:
            print(f"Error loading feature table: {e}")
            return False

    def run_deconvolution(self):
        """Runs the feature deconvolution part of the workflow."""
        if not self.feature_dict:
            print("Feature table not loaded. Cannot run deconvolution.")
            return False
        
        framework = DeconvolutionFramework(self.feature_dict, self.samples, self.deconvolution_params)
        self.deconvoluted_compounds = framework.run()
        return True

    def run_annotation(self):
        """Runs the annotation part of the workflow based on the selected mode."""
        if not self.deconvoluted_compounds:
            print("Deconvolution not performed. Cannot run annotation.")
            return False
        
        if self.mode == 'gc':
            annotator = MSMSAnnotator(self.feature_dict, self.samples, self.params)
            self.deconvoluted_compounds = annotator.annotate_compounds(self.deconvoluted_compounds, self.num_cores)
        elif self.mode == 'lc':
            annotator = FormulaAnnotator(self.lc_annotation_params)
            self.deconvoluted_compounds = annotator.annotate_compounds(self.deconvoluted_compounds)
        else:
            print(f"Warning: Unknown mode '{self.mode}'. No annotation will be performed.")
            return False
        return True

    def save_results(self):
        """Saves the final annotated compound data to a JSON file."""
        #output_path = os.path.join(self.output_path, f"{self.output_prefix}_final.json")
        print(f"\n--- Saving final results to: {self.output_path} ---")
        try:
            with open(self.output_path, 'w') as f:
                json.dump(self.deconvoluted_compounds, f, indent=2)
            print("Save complete.")
        except Exception as e:
            print(f"Error saving results: {e}")
            
    def get_results(self):
        """Returns the final annotated compound data."""
        return self.deconvoluted_compounds