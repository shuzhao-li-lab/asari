"""
Functions for spectral deconvolution of GC-MS data.

1) Clustering method, default in Asari. Use `--deconvolution_algorithm hcl`.
This uses intensity correlation between samples and penalize retention time difference, similar to RamClustR.

2) `--deconvolution_algorithm adap`: Port of ADAP-GC 4.0 in Java, described in https://pmc.ncbi.nlm.nih.gov/articles/PMC6705124/,
source code in mzMine version 2 and 3, as in https://github.com/du-lab/mzmine3/.

3) `--deconvolution_algorithm nmf`: NMF approach as implemented in MShub, Python source code in
https://bitbucket.org/iAnalytica/mshub_process/src/master/,
in function `peakdetect.deconvolve_peaks`.

Entry points, all returning (list_pseudospectra, core_features):
    iterative_build_pseudospectra_by_hcl(list_features_sorted, feature_dataframe, ...)   # 1
    adap_gc_deconvolution(list_features, ri_model=None, ...)                             # 2
    nmf_deconvolution(list_features, feature_dataframe, ri_model=None, ...)              # 3
list_pseudospectra is a list of PseudoSpectrum (same fields as asari.gcms.PseudoSpectrum);
core_features is a set of feature IDs that are members of any pseudospectrum.
Pseudospectrum peaks are [[mz, intensity], ...] on the peak_area scale of the feature table.

Units: asari rtime is in seconds; all time parameters here are in seconds.
The ADAP defaults in mzMine (window 0.2 min, RT tolerance 0.05 min) are converted to seconds.

Inputs: list_features in asari JSON notion (from read_feature_table or
asari.tools.file_io.read_features_from_asari_table, plus 'RI'), and
feature_dataframe of sample intensities indexed by feature ID (sample columns only).
"""

from typing import NamedTuple
import heapq

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.optimize import nnls
from scipy.signal import find_peaks
from scipy.spatial.distance import pdist

from .gc_functions import *

# requires sklearn for nmf method


# ---------------------------------------------------------------------------
# Simulated EICs (Gaussian models as in asari.join.consolidate_track_features)
# ---------------------------------------------------------------------------

def feature_sigmas(features, shape='gaussian', min_sigma=0.1):
    """
    Gaussian standard deviations from feature boundaries.
    gaussian (asari.join):  sigma = 0.33 * max(rtime - left_base, right_base - rtime)
    bigaussian:             left/right sigmas = 0.33 * each half width
    Returns (sigma_left, sigma_right) arrays; equal for 'gaussian'.
    min_sigma guards zero-width features (seen in asari tables).
    """
    rt = np.array([f['rtime'] for f in features], float)
    left = rt - np.array([f['rtime_left_base'] for f in features], float)
    right = np.array([f['rtime_right_base'] for f in features], float) - rt
    if shape == 'gaussian':
        s = np.maximum(0.33 * np.maximum(left, right), min_sigma)
        return s, s
    elif shape == 'bigaussian':
        return np.maximum(0.33 * left, min_sigma), np.maximum(0.33 * right, min_sigma)
    raise ValueError("shape must be 'gaussian' or 'bigaussian'")


def simulate_eics(rt_grid, apex, sigma_left, sigma_right, areas, truncate=4.0):
    """
    Simulate EICs on rt_grid as Gaussian (or bi-Gaussian) peaks of given areas.
    peak_height = area / (sqrt(2*pi) * mean(sigma)), as asari.join for the symmetric case.
    Signal beyond `truncate` sigmas is set to 0, so simulated EICs have finite support like real peaks.
    Returns array (len(rt_grid), n_features).
    """
    t = np.asarray(rt_grid, float)[:, None]
    apex = np.asarray(apex, float)[None, :]
    sl = np.asarray(sigma_left, float)[None, :]
    sr = np.asarray(sigma_right, float)[None, :]
    height = np.asarray(areas, float)[None, :] / (np.sqrt(2 * np.pi) * 0.5 * (sl + sr))
    s = np.where(t < apex, sl, sr)
    z = (t - apex) / s
    out = height * np.exp(-0.5 * z * z)
    out[np.abs(z) > truncate] = 0.0
    return out


# ---------------------------------------------------------------------------
# 1) RI-coupled HCL on sample intensity correlation (asari)
# ---------------------------------------------------------------------------


def ri_penalty_function(abs_ri_delta, min_delta=1, max_delta=100):
    """Penalty by RI distance, using a 1-side sigmoid function. Vectorized; accepts scalars or arrays.
    """
    d = np.asarray(abs_ri_delta, dtype=float)
    out = 2 / (1 + np.exp(np.e * d / max_delta))
    out = np.where(d < min_delta, 1.0, out)
    out = np.where(d > max_delta, 0.0, out)
    return out if out.ndim else float(out)


def _standardize_rows(M):
    """Center and L2-normalize rows so that Pearson r = dot product. Constant rows become 0."""
    M = np.asarray(M, dtype=float)
    M = M - M.mean(axis=1, keepdims=True)
    norm = np.sqrt((M * M).sum(axis=1, keepdims=True))
    norm[norm == 0] = np.inf
    return M / norm


def get_clusters_by_ri_hcl(
    feature_dataframe_sorted,
    dict_features,
    ri_penalty_function=ri_penalty_function,
    min_delta=1,
    max_delta=100,
    hcl_distance_cut=1,
    limit_fearture_num=None,
    linkage_method='ward',
    verbose=True,
    ):
    """
    Using RI (retention index) coupled hierarchical clustering (hcl) to find feature clusters.
    - the RI penalty matrix is computed vectorized (was an O(n^2) Python loop);
    - zero-variance features (NaN correlation) get distance 1 instead of breaking linkage;
    - linkage_method is exposed (default 'ward' as in asari).

    feature_dataframe_sorted : dataframe indexed on feature IDs, sample intensities only,
                    sorted by combined peak area, descending order.
    dict_features : dictionary of features, feature IDs as keys.
    Returns metClus, metClusDict.
    """
    if limit_fearture_num:
        feature_dataframe_sorted = feature_dataframe_sorted.iloc[:limit_fearture_num, :]
    n = feature_dataframe_sorted.shape[0]
    if verbose:
        print(f"Working with {n} features.")
    if n < 2:
        return np.ones(n, dtype=int), ({1: list(feature_dataframe_sorted.index)} if n else {})
    retention_index = np.array([dict_features[x]['RI'] for x in feature_dataframe_sorted.index], float)
    ri_penalties = ri_penalty_function(pdist(retention_index[:, None], 'cityblock'), min_delta, max_delta)
    with np.errstate(invalid='ignore', divide='ignore'):
        YM = pdist(feature_dataframe_sorted.values.astype(float), 'correlation')
    YM = np.nan_to_num(YM, nan=1.0)
    YM = 1 - ri_penalties * (1 - YM)
    ZM = linkage(YM, method=linkage_method)
    metClus = fcluster(ZM, hcl_distance_cut, criterion='distance')

    metClusDict = {}
    for ii, label in enumerate(metClus):
        metClusDict.setdefault(label, []).append(feature_dataframe_sorted.index[ii])
    if verbose:
        print("number of features: ", n)
        print("number of clusters: ", len(metClusDict))
    return metClus, metClusDict


def _split_by_ri(cluster, all_ri, pos, max_span):
    """Split a cluster (list of feature IDs) so that each piece spans <= max_span in RI,
    by complete linkage on RI. Pieces with fewer than 2 features are dropped."""
    ri = all_ri[[pos[f] for f in cluster]]
    if np.ptp(ri) <= max_span:
        return [cluster]
    labels = fcluster(linkage(ri[:, None], method='complete'), max_span, criterion='distance')
    pieces = {}
    for f, label in zip(cluster, labels):
        pieces.setdefault(label, []).append(f)
    return [p for p in pieces.values() if len(p) > 1]


def format_fcluster_to_pseudospectrum(selected_features, seed_feature, annotation=''):
    """Format a cluster of features as a PseudoSpectrum, using seed_feature for id, rtime and RI."""
    return PseudoSpectrum(
        seed_feature['id_number'],
        seed_feature['rtime'],
        seed_feature['RI'],
        [unit_mz(x['mz']) for x in selected_features],
        len(selected_features),
        [x['id'] for x in selected_features],           # members
        np.array([(x['mz'], x['peak_area']) for x in selected_features]),
        annotation
    )


def iterative_build_pseudospectra_by_hcl(
    list_features_sorted,
    feature_dataframe,
    init_core_features=None,
    step_size=2000,   # batch size for HCL. Slow if batch too big.
    hcl_distance_cut=1,
    ri_tolerance='std',
    ri_window_min=3,
    ri_window_max=50,
    split_wide_clusters=True,
    correlation_cut=0.7,
    max_core_features=20000,
    min_ri_delta=1,
    max_ri_delta=100,
    linkage_method='ward',
    verbose=True,
    ):
    """
    Build pseudospectra from list_features_sorted iteratively until reaching max_core_features
    (no more than 90% of all features):
    - take the top step_size remaining features, RI-penalized correlation HCL (get_clusters_by_ri_hcl);
    - each cluster of >1 features is extended by correlated features (Pearson >= correlation_cut,
      RI within the RI window) among the remaining features beyond the batch;
    - cluster members become core features and are excluded from the next iteration.

    - HCL clusters whose member RIs span more than max_ri_delta are split
      (split_wide_clusters=True): features that far apart get no RI-penalty credit, yet Ward
      linkage at hcl_distance_cut=1 can still merge them, giving clusters across the whole run.
      Clusters are cut by complete linkage on RI at max_ri_delta; pieces of 1 feature are dropped;
    - the RI window for cluster extension is adaptive: ri_tolerance='std' (default) is
      median(cluster RI) +/- 1 robust standard deviation (1.4826 * MAD) of the cluster's RIs,
      clamped to [ri_window_min, ri_window_max]. A number gives a fixed window around the seed
      RI, as in asari (ri_tolerance=50, split_wide_clusters=False reproduces asari);
    - cluster extension uses one precomputed standardized intensity matrix (dot products),
      instead of a pandas correlation per candidate feature;
    - stops when an iteration adds no new core features (asari could loop forever);
    - init_core_features is not a shared mutable default;
    - rounded_mzs use plain rounding.

    list_features_sorted : features in asari JSON notion with 'RI', sorted by peak_area descending.
    feature_dataframe : sample intensities, indexed by feature ID.
    hcl_distance_cut : parameter to cut linkage tree in HCL.
        Smaller values give tighter clusters, which are desired since they are extended after this step.
    returns list_pseudospectra, core_features
    """
    core_features = set(init_core_features or ())
    list_pseudospectra = []
    max_core_features = min(max_core_features, 0.9 * len(list_features_sorted))  # safeguard
    featureDict = {x['id']: x for x in list_features_sorted}
    all_ids = [f['id'] for f in list_features_sorted]
    pos = {fid: i for i, fid in enumerate(all_ids)}
    Z = _standardize_rows(feature_dataframe.loc[all_ids, :].values)
    all_ri = np.array([f['RI'] for f in list_features_sorted], float)

    counter = 0
    while len(core_features) <= max_core_features:
        counter += 1
        remaining = [fid for fid in all_ids if fid not in core_features]
        if len(remaining) < 2:
            break
        if verbose:
            print(f"Iteration {counter} of building pseudospectra from {len(remaining)} remaining features.")
        _, metClusDict = get_clusters_by_ri_hcl(
            feature_dataframe.loc[remaining[:step_size], :],
            featureDict,
            min_delta=min_ri_delta,
            max_delta=max_ri_delta,
            hcl_distance_cut=hcl_distance_cut,
            linkage_method=linkage_method,
            verbose=verbose,
        )
        pool = np.array([pos[fid] for fid in remaining[step_size:]], dtype=int)
        pool_ri = all_ri[pool] if pool.size else np.empty(0)
        accounted = set()
        clusters = [c for c in metClusDict.values() if len(c) > 1]
        if split_wide_clusters:
            clusters = [p for c in clusters for p in _split_by_ri(c, all_ri, pos, max_ri_delta)]
        for clu in clusters:
            seed = max((featureDict[f] for f in clu), key=lambda x: x['peak_area'])
            if ri_tolerance == 'std':
                clu_ri = all_ri[[pos[f] for f in clu]]
                ri_center = float(np.median(clu_ri))
                ri_window = float(np.clip(1.4826 * np.median(np.abs(clu_ri - ri_center)),
                                          ri_window_min, ri_window_max))
            else:
                ri_center, ri_window = seed['RI'], ri_tolerance
            cand = pool[np.abs(pool_ri - ri_center) < ri_window] if pool.size else pool
            extended = []
            if cand.size:
                corr = Z[cand] @ Z[pos[seed['id']]]
                extended = [all_ids[i] for i in cand[corr >= correlation_cut]]
            members = list(clu) + extended
            list_pseudospectra.append(
                format_fcluster_to_pseudospectrum([featureDict[f] for f in members], seed))
            accounted.update(members)

        new_core = core_features | accounted
        if len(new_core) == len(core_features):
            if verbose:
                print("No new clusters in this iteration; stopping.")
            break
        core_features = new_core
        if verbose:
            print(f"Core features collected: {len(core_features)}.")

    return list_pseudospectra, core_features


# ---------------------------------------------------------------------------
# 2) ADAP-GC 4.0: clustering-assisted multivariate curve resolution
# ---------------------------------------------------------------------------
#
# Follows dulab.adap.workflow.decomposition (ADAP 4.1.10, as bundled in mzMine 3):
#   RetTimeClusterer   -> _adap_rt_windows
#   BiGaussianDetector -> _adap_bigaussians
#   ComponentSelector  -> _adap_decompose_window (ALS with ModifiedUpdateRule constraint, then NNLS)
# Each feature is treated as one EIC peak, simulated from its apex, boundaries and peak_area.

_ADAP_MIN_WINDOW_PEAKS = 4     # RetTimeClusterer.MIN_CLUSTER_SIZE; windows need > 4 peaks
_HWHM_OVER_STD = np.sqrt(2 * np.log(2.0))


def _merge_intervals_complete_linkage(starts, ends, max_width):
    """
    Complete-linkage clustering of peak intervals with ADAP distance
        d(p1, p2) = max(p1.end - p2.start, p2.end - p1.start),
    cut at max_width. Between clusters, complete linkage gives
        D(A, B) = max(max_end(A) - min_start(B), max_end(B) - min_start(A)).
    For 1-D intervals, merges are taken between neighbours in midpoint order
    (a heap of adjacent pairs), which avoids the O(n^2) distance matrix.
    Returns list of index arrays.
    """
    n = len(starts)
    if n == 0:
        return []
    order = np.argsort((starts + ends) / 2, kind='stable')
    lo = starts[order].astype(float).copy()
    hi = ends[order].astype(float).copy()
    members = [[int(i)] for i in order]
    nxt = list(range(1, n)) + [-1]
    prv = [-1] + list(range(n - 1))
    alive = [True] * n
    version = [0] * n

    def dist(a, b):
        return max(hi[a] - lo[b], hi[b] - lo[a])

    heap = [(dist(i, i + 1), i, 0, 0) for i in range(n - 1)]
    heapq.heapify(heap)
    while heap:
        d, a, va, vb = heapq.heappop(heap)
        if d >= max_width:
            break
        b = nxt[a]
        if not alive[a] or b < 0 or version[a] != va or version[b] != vb:
            continue
        # merge b into a
        lo[a], hi[a] = min(lo[a], lo[b]), max(hi[a], hi[b])
        members[a] += members[b]
        alive[b] = False
        version[a] += 1
        nxt[a] = nxt[b]
        if nxt[b] >= 0:
            prv[nxt[b]] = a
        if prv[a] >= 0:
            p = prv[a]
            heapq.heappush(heap, (dist(p, a), p, version[p], version[a]))
        if nxt[a] >= 0:
            c = nxt[a]
            heapq.heappush(heap, (dist(a, c), a, version[a], version[c]))
    return [np.array(members[i]) for i in range(n) if alive[i]]


def _adap_rt_windows(starts, ends, apex, areas, sig_l, sig_r, max_width, scan_time,
                     min_window_peaks=_ADAP_MIN_WINDOW_PEAKS + 1, truncate=4.0):
    """
    RetTimeClusterer.execute: pre-cluster peaks by gaps; split clusters wider than max_width
    by complete linkage; move boundaries between overlapping adjacent windows to the TIC minimum;
    reassign peaks by apex and keep windows with >= min_window_peaks peaks.
    Returns list of (start, end, peak_indices).
    """
    order = np.argsort(starts, kind='stable')
    clusters = []
    current, cur_max = [], -np.inf
    for i in order:
        if starts[i] > cur_max and current:
            clusters.append(np.array(current))
            current = []
        cur_max = max(cur_max, ends[i])
        current.append(i)
    if current:
        clusters.append(np.array(current))

    windows = []
    for c in clusters:
        if ends[c].max() - starts[c].min() > max_width:
            for sub in _merge_intervals_complete_linkage(starts[c], ends[c], max_width):
                if sub.size >= min_window_peaks:
                    windows.append(c[sub])
        else:
            windows.append(c)
    if not windows:
        return []

    # boundaries: [min start, max end], ordered by mean apex
    windows.sort(key=lambda w: apex[w].mean())
    bounds = [[starts[w].min(), ends[w].max()] for w in windows]

    # TIC of all window peaks on a global grid (simulated EICs)
    grid_lo = min(b[0] for b in bounds) - scan_time
    grid = np.arange(grid_lo, max(b[1] for b in bounds) + 2 * scan_time, scan_time)
    tic = np.zeros(grid.size)
    used = np.unique(np.concatenate(windows))
    for i in used:
        a = max(0, int((apex[i] - truncate * sig_l[i] - grid_lo) / scan_time))
        b = min(grid.size, int((apex[i] + truncate * sig_r[i] - grid_lo) / scan_time) + 2)
        if b > a:
            tic[a:b] += simulate_eics(grid[a:b], [apex[i]], [sig_l[i]], [sig_r[i]], [areas[i]],
                                      truncate)[:, 0]
    for k in range(1, len(bounds)):
        w_start, w_end = bounds[k][0], bounds[k - 1][1]
        if w_end <= w_start:
            continue
        inside = (grid > w_start) & (grid < w_end)
        boundary = grid[inside][np.argmin(tic[inside])] if inside.any() else (w_start + w_end) / 2
        bounds[k - 1][1] = boundary
        bounds[k][0] = boundary

    # reassign peaks by apex
    all_peaks = np.unique(np.concatenate(windows))
    all_apex = apex[all_peaks]
    result, seen = [], set()
    for lo, hi in bounds:
        peaks = all_peaks[(all_apex > lo) & (all_apex < hi)]
        key = (peaks.size, float(apex[peaks].mean()) if peaks.size else 0.0)
        if peaks.size >= min_window_peaks and key not in seen:
            seen.add(key)
            result.append((lo, hi, peaks))
    return result


def _upper_median(x):
    """Median as in BiGaussianDetector.median: sorted()[n // 2]."""
    x = np.sort(np.asarray(x))
    return x[x.size // 2]


def _adap_bigaussians(apex, hwhm_left, hwhm_right, heights, ret_time_tolerance,
                      intensity_factor=0.0, min_cluster_size=1):
    """
    BiGaussianDetector: filter peaks by height > intensity_factor * max height;
    complete-linkage clustering of apex retention times cut at ret_time_tolerance;
    each cluster (>= min_cluster_size peaks) -> bi-Gaussian with median apex and median
    left/right half widths at half maximum, sigma = HWHM / sqrt(2 ln 2).
    Returns list of (center, sigma_left, sigma_right).
    """
    keep = heights > intensity_factor * heights.max()
    apex, hwhm_left, hwhm_right = apex[keep], hwhm_left[keep], hwhm_right[keep]
    if apex.size == 0:
        return []
    if apex.size == 1:
        labels = np.zeros(1, int)
    else:
        Z = linkage(apex[:, None], method='complete')
        labels = (np.zeros(apex.size, int) if Z[-1, 2] < ret_time_tolerance
                  else fcluster(Z, ret_time_tolerance, criterion='distance'))
    out = []
    for label in pd.unique(labels):
        idx = labels == label
        if idx.sum() < min_cluster_size:
            continue
        out.append((
            _upper_median(apex[idx]),
            max(_upper_median(hwhm_left[idx]), 1e-6) / _HWHM_OVER_STD,
            max(_upper_median(hwhm_right[idx]), 1e-6) / _HWHM_OVER_STD,
        ))
    return out


def _bigaussian(t, center, sl, sr):
    s = np.where(t < center, sl, sr)
    d = t - center
    return np.exp(-d * d / (2 * np.maximum(s * s, 1e-12)))


def _constrain_components(W):
    """
    ModifiedUpdateRule: make each component (column) unimodal, set the baseline (last column)
    to 1, and scale each column to max 1.
    """
    W = np.maximum(W, 0.0)
    for j in range(W.shape[1] - 1):
        w = W[:, j]
        m = int(np.argmax(w))
        w[:m + 1] = np.minimum.accumulate(w[:m + 1][::-1])[::-1]
        w[m:] = np.minimum.accumulate(w[m:])
        mx = w[m]
        if mx > 0:
            w /= mx
    W[:, -1] = 1.0
    return W


def _als_factorization(X, W, H, max_iter=500, tol=1e-6):
    """
    Alternating least squares X ~ W H with non-negativity, W constrained by _constrain_components
    (javanmf AlternatingLeastSquaresMatrixFactorization with ComponentConstraint in ADAP).
    """
    prev = np.inf
    norm_x = np.linalg.norm(X) or 1.0
    for _ in range(max_iter):
        H = np.maximum(np.linalg.lstsq(W, X, rcond=None)[0], 0.0)
        W = np.maximum(np.linalg.lstsq(H.T, X.T, rcond=None)[0].T, 0.0)
        W = _constrain_components(W)
        err = np.linalg.norm(X - W @ H) / norm_x
        if abs(prev - err) < tol:
            break
        prev = err
    return W, H


def _adap_decompose_window(t, X_peaks, peak_params, X_all, ret_time_tolerance,
                           intensity_factor, min_cluster_size, max_iter, tol):
    """
    ComponentSelector.execute for one window.
    X_peaks : simulated EICs (time x peaks) of peaks with apex in the window.
    peak_params : (apex, hwhm_left, hwhm_right) arrays for those peaks.
    X_all : simulated EICs (time x chromatograms) of all features with signal in the window.
    Returns (components, spectra): components (time x k) scaled to max 1, and
    spectra (k x chromatograms) of NNLS coefficients (heights).
    """
    apex, hl, hr = peak_params
    heights = X_peaks.max(axis=0)
    bgs = _adap_bigaussians(apex, hl, hr, heights, ret_time_tolerance, intensity_factor, min_cluster_size)
    if not bgs:
        return None, None

    # NMF initialized by bi-Gaussians, plus a flat baseline component
    colmax = X_peaks.max(axis=0)
    colmax[colmax <= 0] = 1.0
    Xp = X_peaks / np.sqrt(colmax)
    W0 = np.column_stack([_bigaussian(t, c, sl, sr) for c, sl, sr in bgs] + [np.ones(t.size)])
    H0 = np.vstack([Xp[np.argmax(W0[:, j]), :] for j in range(len(bgs))]
                   + [np.full(Xp.shape[1], Xp.max() * 1e-12)])
    W, _ = _als_factorization(Xp, W0, H0, max_iter=max_iter, tol=tol)

    # least-squares decomposition of all chromatograms on the components + baseline
    comps = W[:, :-1]
    comps = comps[:, comps.max(axis=0) > 0]
    if comps.shape[1] == 0:
        return None, None
    comps = comps / comps.max(axis=0)
    limit = X_all[np.argmax(comps, axis=0), :]          # intensity at each component apex
    A = np.column_stack([comps, np.ones(t.size)])
    coef = np.empty((A.shape[1], X_all.shape[1]))
    for j in range(X_all.shape[1]):
        coef[:, j] = nnls(A, X_all[:, j])[0]
    spectra = np.minimum(coef[:-1, :], limit)
    return comps, spectra


def adap_gc_deconvolution(
    list_features,
    ri_model=None,
    feature_dataframe=None,
    intensity_column=None,
    pref_window_width=12.0,     # mzMine default 0.2 min
    ret_time_tolerance=3.0,     # mzMine default 0.05 min
    min_cluster_size=1,
    intensity_factor=0.0,
    adjust_apex_ret_time=False, # not applicable: simulated EIC apex is the feature rtime
    eic_shape='gaussian',
    scan_time=0.1,
    min_sigma=0.1,
    truncate=4.0,
    min_feature_fraction=0.1,
    min_members=2,
    max_iter=500,
    tol=1e-6,
    verbose=True,
    ):
    """
    ADAP-GC 4.0 deconvolution (clustering-assisted multivariate curve resolution) on simulated EICs.

    Steps (dulab.adap.workflow.decomposition):
    1. Deconvolution windows from peak intervals [rtime_left_base, rtime_right_base]:
       split at gaps, split wider than pref_window_width by complete linkage,
       move boundaries of overlapping windows to the TIC minimum. Windows need > 4 peaks.
    2. Per window, cluster apex RTs by complete linkage at ret_time_tolerance; each cluster with
       >= min_cluster_size peaks gives a bi-Gaussian initial component (median apex and widths).
    3. Constrained ALS factorization of the window peaks (unimodal components + baseline).
    4. NNLS decomposition of every feature EIC with signal in the window on the components,
       coefficients capped at the EIC intensity at the component apex; the coefficients are the
       component spectrum.

    Input intensities: composite peak_area by default; or a sample column of feature_dataframe
    (intensity_column) to deconvolve one sample, as ADAP does per raw file.

    Output: a feature is a member of a component if its NNLS contribution (coefficient times
    component area) is >= min_feature_fraction of the feature's simulated area. Peak intensity
    = that contribution, on the peak_area scale. Components need >= min_members members.
    Pseudospectrum id is the id_number of the feature with the largest coefficient (the ADAP
    "best peak"); rtime is the component apex.

    returns list_pseudospectra, core_features
    """
    feats = list(list_features)
    if intensity_column is not None:
        col = feature_dataframe[intensity_column]
        areas = np.array([col.get(f['id'], 0.0) for f in feats], float)
    else:
        areas = np.array([f['peak_area'] for f in feats], float)
    keep = areas > 0
    feats = [f for f, k in zip(feats, keep) if k]
    areas = areas[keep]
    if not feats:
        return [], set()

    apex = np.array([f['rtime'] for f in feats], float)
    sig_l, sig_r = feature_sigmas(feats, eic_shape, min_sigma)
    starts = np.minimum(np.array([f['rtime_left_base'] for f in feats], float), apex - sig_l)
    ends = np.maximum(np.array([f['rtime_right_base'] for f in feats], float), apex + sig_r)
    sup_lo, sup_hi = apex - truncate * sig_l, apex + truncate * sig_r
    hwhm_l, hwhm_r = sig_l * _HWHM_OVER_STD, sig_r * _HWHM_OVER_STD

    windows = _adap_rt_windows(starts, ends, apex, areas, sig_l, sig_r,
                               pref_window_width, scan_time, truncate=truncate)
    if verbose:
        print(f"ADAP: {len(feats)} features in {len(windows)} deconvolution windows.")

    order_sup = np.argsort(sup_lo, kind='stable')
    sup_lo_sorted = sup_lo[order_sup]
    max_support = (sup_hi - sup_lo).max()
    list_pseudospectra, core_features, used_ids = [], set(), {}
    for w_start, w_end, peaks in windows:
        t = np.arange(w_start, w_end + scan_time / 2, scan_time)
        if t.size < 3:
            continue
        # all features with signal in [w_start, w_end]
        a = np.searchsorted(sup_lo_sorted, w_start - max_support)
        b = np.searchsorted(sup_lo_sorted, w_end, side='right')
        cand = order_sup[a:b]
        cand = cand[(sup_hi[cand] >= w_start) & (sup_lo[cand] <= w_end)]
        X_all = simulate_eics(t, apex[cand], sig_l[cand], sig_r[cand], areas[cand], truncate)
        has_signal = X_all.sum(axis=0) > 0
        cand, X_all = cand[has_signal], X_all[:, has_signal]
        local = {g: i for i, g in enumerate(cand)}
        peaks = np.array([p for p in peaks if p in local], dtype=int)
        if peaks.size == 0:
            continue
        X_peaks = X_all[:, [local[p] for p in peaks]]
        comps, spectra = _adap_decompose_window(
            t, X_peaks, (apex[peaks], hwhm_l[peaks], hwhm_r[peaks]), X_all,
            ret_time_tolerance, intensity_factor, min_cluster_size, max_iter, tol)
        if comps is None:
            continue

        comp_area = comps.sum(axis=0) * scan_time                 # area of unit-height profile
        contrib = spectra * comp_area[:, None]                    # k x chromatograms, area scale
        feat_area = X_all.sum(axis=0) * scan_time
        frac = contrib / np.where(feat_area > 0, feat_area, np.inf)
        for j in range(comps.shape[1]):
            sel = np.where((spectra[j] > 0) & (frac[j] >= min_feature_fraction))[0]
            if sel.size < min_members:
                continue
            best = cand[sel[np.argmax(spectra[j, sel])]]
            rtime = float(t[np.argmax(comps[:, j])])
            members = [feats[cand[i]] for i in sel]
            list_pseudospectra.append(_make_pseudospectrum(
                feats[best]['id_number'], rtime, _pseudospectrum_ri(rtime, ri_model, feats[best]),
                members, contrib[j, sel], used_ids))
            core_features.update(f['id'] for f in members)

    return list_pseudospectra, core_features


def _pseudospectrum_ri(rtime, ri_model, base_feature):
    """RI for a pseudospectrum: from ri_model at rtime if given,
    else the base feature's own RI (e.g. from a feature table with an RI column), else None."""
    if ri_model is not None:
        return rtime_to_ri(rtime, ri_model)
    ri = base_feature.get('RI')
    return None if ri is None or ri != ri else float(ri)     # ri != ri catches NaN


def _make_pseudospectrum(base_id, rtime, ri, members, intensities, used_ids, annotation=''):
    """Build PseudoSpectrum; ids are made unique by suffix if a base feature is reused."""
    n = used_ids.get(base_id, 0)
    used_ids[base_id] = n + 1
    ps_id = base_id if n == 0 else f"{base_id}_{n + 1}"
    return PseudoSpectrum(
        ps_id, rtime, ri,
        [unit_mz(f['mz']) for f in members],
        len(members),
        [f['id'] for f in members],
        np.array([(f['mz'], float(v)) for f, v in zip(members, intensities)]),
        annotation
    )


# ---------------------------------------------------------------------------
# 3) NMF deconvolution (MSHub, peakdetect.deconvolve_peaks)
# ---------------------------------------------------------------------------

def median_threshold(X):
    """As in MSHub: median + 3 * MAD."""
    md = np.median(X)
    MAD = np.median(np.abs(np.subtract(X, md))) * 1.4826
    return md + 3 * MAD


def _mshub_rt_groups(features, intensity_matrix, rt_tol='auto', bin_width=0.1,
                     weighted_density=False, smooth=False):
    """
    MSHub cluster_rts_kernel: histogram of peak apex RTs over all samples, local maxima with
    minimal separation rt_tol, maxima above max(median_threshold * 0.05, 5) define common RTs;
    each feature is assigned to its nearest common RT if within rt_tol, otherwise left out.

    A feature present in k samples counts k times (MSHub counts peaks per sample);
    weighted_density uses per-sample intensity fractions instead.
    rt_tol='auto' is median peak width / 5, as in MSHub.
    smooth: optional Gaussian smoothing of the histogram (sigma = rt_tol / 2), not in MSHub.
    Returns (crt array, list of feature index arrays).
    """
    rt = np.array([f['rtime'] for f in features], float)
    if rt_tol == 'auto':
        widths = np.array([f['rtime_right_base'] - f['rtime_left_base'] for f in features], float)
        rt_tol = float(np.median(widths)) / 5
    present = intensity_matrix > 0
    if weighted_density:
        col = intensity_matrix.sum(axis=0)
        weights = (intensity_matrix / np.where(col > 0, col, 1)).sum(axis=1)
    else:
        weights = present.sum(axis=1).astype(float)
    edges = np.arange(rt.min() - bin_width, rt.max() + 2 * bin_width, bin_width)
    histc = np.histogram(rt, edges, weights=weights)[0].astype(float)
    if smooth:
        from scipy.ndimage import gaussian_filter1d
        histc = gaussian_filter1d(histc, max(rt_tol / 2 / bin_width, 0.5))
    gap = max(int(rt_tol / bin_width), 1)
    idx, _ = find_peaks(histc, distance=gap)
    threshold = max(median_threshold(histc) * 0.05, 5) if not weighted_density \
        else median_threshold(histc) * 0.05
    idx = idx[histc[idx] > threshold]
    if idx.size == 0:
        return np.empty(0), []
    crt = (edges[idx] + edges[idx + 1]) / 2
    nearest = np.round(np.interp(rt, crt, np.arange(crt.size))).astype(int)
    ok = np.abs(crt[nearest] - rt) <= rt_tol
    groups = [np.where(ok & (nearest == k))[0] for k in range(crt.size)]
    return crt, groups


def _mshub_adaptive_nmf(X, max_components=10, min_var_gain=5.0, max_iter=500, random_state=0):
    """
    MSHub deconvolve_peaks: increase NMF components (init 'nndsvd') while the cumulative
    explained variance, round(100 * (1 - ||X - WH||^2 / ||X||^2)), gains >= min_var_gain;
    stop at the component limit or if the gain stops decreasing (equal consecutive gains).
    Returns W, H, n_components.
    """
    from sklearn.decomposition import NMF
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    denom = np.sum(X * X)
    max_comps = min(max_components, *X.shape)
    cum_var = []
    best = prev_best = None
    i_comp = 0
    while True:
        i_comp += 1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', ConvergenceWarning)
            model = NMF(n_components=i_comp, init='nndsvd', max_iter=max_iter, random_state=random_state)
            Wcur = model.fit_transform(X)
        Hcur = model.components_
        cum_var.append(np.round(100 * (1 - np.sum((X - Wcur @ Hcur) ** 2) / denom)))
        if i_comp == 1:
            best = prev_best = (Wcur, Hcur)
            if max_comps <= 1:
                break
            continue
        if i_comp >= max_comps:
            # MSHub breaks before accepting the component at the limit; keep it if it helps
            if cum_var[-1] - cum_var[-2] >= min_var_gain:
                prev_best, best = best, (Wcur, Hcur)
            else:
                i_comp -= 1
            break
        if cum_var[i_comp - 1] - cum_var[i_comp - 2] < min_var_gain:
            i_comp -= 1
            break
        if i_comp > 2:
            # linear decrease in component variance is likely due to chance
            if cum_var[i_comp - 1] - cum_var[i_comp - 2] == cum_var[i_comp - 2] - cum_var[i_comp - 3]:
                i_comp -= 2
                best = prev_best
                break
        prev_best, best = best, (Wcur, Hcur)
    W, H = best
    return W, H, W.shape[1]


def nmf_deconvolution(
    list_features,
    feature_dataframe,
    ri_model=None,
    rt_tol='auto',
    bin_width=0.1,
    weighted_density=False,
    smooth_histogram=False,
    use_time_profiles=True,
    eic_shape='gaussian',
    scan_time=0.2,
    min_sigma=0.1,
    truncate=4.0,
    max_samples=None,
    max_components=10,
    min_var_gain=5.0,
    min_feature_fraction=0.1,
    min_members=2,
    max_iter=500,
    random_state=0,
    verbose=True,
    ):
    """
    NMF deconvolution after MSHub (peakdetect.deconvolve_peaks), on simulated EICs.

    1. Features are grouped by apex RT (MSHub kernel/histogram grouping, _mshub_rt_groups).
    2. For each group, a data matrix X (rows x features) is built:
       - use_time_profiles=True: for every sample s, the simulated EIC of feature f is
         intensity[s, f] * unit-area Gaussian(t); rows are (sample, time) pairs stacked.
         Summing X over time gives MSHub's samples x m/z integral matrix, so this adds peak-shape
         information to MSHub's cross-sample factorization.
       - use_time_profiles=False: X = samples x features intensities (MSHub as is).
    3. NMF with an adaptive number of components (MSHub rule: stop when cumulative explained
       variance gains < min_var_gain percent). H rows are component spectra.
    4. Feature f is a member of component c if its share of X[:, f] reconstructed by c is
       >= min_feature_fraction; peak intensity = share * peak_area (peak_area scale).
       Components need >= min_members members.

    max_samples : use only the top samples by group total intensity (limits matrix size).
    Pseudospectrum id is the id_number of the feature with the largest normalized loading;
    rtime is the apex of the sample-summed component profile (or the loading-weighted mean apex
    of members if use_time_profiles=False).

    returns list_pseudospectra, core_features
    """
    feats = list(list_features)
    ids = [f['id'] for f in feats]
    I_all = feature_dataframe.loc[ids, :].values.astype(float)   # features x samples
    crt, groups = _mshub_rt_groups(feats, I_all, rt_tol, bin_width, weighted_density, smooth_histogram)
    if verbose:
        print(f"NMF: {len(feats)} features grouped into {len(groups)} RT groups; "
              f"{sum(len(g) for g in groups)} features assigned.")

    apex = np.array([f['rtime'] for f in feats], float)
    sig_l, sig_r = feature_sigmas(feats, eic_shape, min_sigma)
    areas = np.array([f['peak_area'] for f in feats], float)

    list_pseudospectra, core_features, used_ids = [], set(), {}
    for g in groups:
        if g.size < min_members:
            continue
        I = I_all[g, :].T                                  # samples x features
        if max_samples and I.shape[0] > max_samples:
            I = I[np.argsort(-I.sum(axis=1))[:max_samples]]
        I = I[I.sum(axis=1) > 0]
        if I.shape[0] == 0:
            continue
        if use_time_profiles:
            t = np.arange((apex[g] - truncate * sig_l[g]).min(),
                          (apex[g] + truncate * sig_r[g]).max() + scan_time / 2, scan_time)
            shapes = simulate_eics(t, apex[g], sig_l[g], sig_r[g], np.ones(g.size), truncate)
            shapes /= np.maximum(shapes.sum(axis=0) * scan_time, 1e-300)   # unit area
            X = (I[:, None, :] * shapes[None, :, :]).reshape(-1, g.size)
        else:
            X = I
        keep_cols = X.sum(axis=0) > 0
        if keep_cols.sum() < min_members:
            continue
        g, X = g[keep_cols], X[:, keep_cols]
        if use_time_profiles:
            shapes = shapes[:, keep_cols]

        W, H, k = _mshub_adaptive_nmf(X, max_components, min_var_gain, max_iter, random_state)
        contrib = H * W.sum(axis=0)[:, None]               # k x features, reconstructed totals
        frac = contrib / X.sum(axis=0)[None, :]
        explained = contrib.sum(axis=1)
        for c in np.argsort(-explained, kind='stable'):
            sel = np.where(frac[c] >= min_feature_fraction)[0]
            if sel.size < min_members:
                continue
            pattern = H[c] / H[c].max()
            best = g[sel[np.argmax(pattern[sel])]]
            if use_time_profiles:
                profile = W[:, c].reshape(-1, t.size).sum(axis=0)
                rtime = float(t[np.argmax(profile)])
            else:
                rtime = float(np.average(apex[g[sel]], weights=H[c, sel]))
            members = [feats[i] for i in g[sel]]
            list_pseudospectra.append(_make_pseudospectrum(
                feats[best]['id_number'], rtime, _pseudospectrum_ri(rtime, ri_model, feats[best]),
                members, np.minimum(frac[c, sel], 1.0) * areas[g[sel]], used_ids))
            core_features.update(f['id'] for f in members)

    return list_pseudospectra, core_features
