"""
Batch entropy-similarity search, vectorized across queries.

For a chunk of queries, all query peaks are matched against the library at
once, producing a flat list of matched (query peak, library peak) pairs.
Each pair contributes

    g(a + b) - g(a) - g(b),   g(x) = x log2 x,   a, b = half the peak intensities

to its (query, library spectrum) score; contributions are summed per pair of
spectra. Spectrum pairs that share no peak are never touched, so the cost is
proportional to the number of matched peaks, not to n_query * n_library.

Two candidate generators:
  open      every library peak within ms2_tolerance, from the fragment-sorted
            index (binary search on all query peaks at once)
  identity  only library spectra with precursor within ms1_tolerance; their
            peaks are contiguous in the precursor-sorted layout, so candidate
            pairs are generated directly and then filtered by m/z
"""

import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np

from .index import SpectralIndex
from .preprocess import expand_ranges, spectrum_ids


DENSE_FACTOR = 16            # use dense accumulation if n_query * n_library <= this * pairs
DENSE_MAX_CELLS = 1 << 20    # cells per dense chunk (score + count arrays: 16 MB)


@dataclass
class SearchResult:
    """
    Flat table of hits, sorted by query then descending score.
    query_idx / library_idx are positions in the original input lists.
    """
    query_idx: np.ndarray
    library_idx: np.ndarray
    score: np.ndarray
    n_matched: np.ndarray

    def __len__(self):
        return len(self.score)

    def to_dense(self, n_query, n_library):
        """Dense (n_query, n_library) score matrix; zeros for unscored pairs."""
        m = np.zeros((n_query, n_library), dtype=np.float32)
        m[self.query_idx, self.library_idx] = self.score
        return m

    def to_dataframe(self):
        """Return hits as a pandas DataFrame."""
        import pandas as pd
        return pd.DataFrame({'query_idx': self.query_idx, 'library_idx': self.library_idx,
                             'score': self.score, 'n_matched': self.n_matched})


def search(library, queries, method='identity', ms1_tolerance=0.01, ms2_tolerance=0.02,
           top_k=None, min_score=0.0, exclude_self=False, max_pairs=250_000, n_threads=1):
    """
    Search a batch of query spectra against a library index.

    library     SpectralIndex (needs the fragment index for open search).
    queries     SpectralIndex, or a list of spectrum dicts, which are cleaned
                with the library's parameters.
    method      'identity' (precursor within ms1_tolerance Da) or 'open'.
    top_k       keep at most this many hits per query (None: all).
    min_score   keep hits with score >= min_score (hits always share >= 1 peak).
    exclude_self  drop hits where query_idx == library_idx, for all-vs-all
                search of a dataset against itself.
    max_pairs   matched peak pairs per chunk; small chunks stay in CPU cache.
    n_threads   process chunks in parallel threads.

    Returns SearchResult.
    """
    queries = _as_queries(queries, library)
    _check_compatible(library, queries, ms2_tolerance)
    if method == 'open':
        if library.frag_mz is None:
            raise ValueError('open search needs a library built with fragment_index=True')
        start, per_peak = _open_candidates(library, queries, ms2_tolerance)
        pairs_fn = _open_pairs
    elif method == 'identity':
        start, per_peak = _identity_candidates(library, queries, ms1_tolerance)
        pairs_fn = _identity_pairs
    else:
        raise ValueError("method must be 'identity' or 'open'")

    # derived arrays, computed once here rather than inside worker threads
    qspec = queries.cached('spec_ids', lambda: spectrum_ids(queries.offsets))
    _half_and_g(queries, 'intensity')
    _half_and_g(library, 'frag_intensity' if method == 'open' else 'intensity')
    library.cached('spec_ids', lambda: spectrum_ids(library.offsets))
    _mz_bounds(queries, ms2_tolerance)

    per_query = np.bincount(qspec, weights=per_peak, minlength=len(queries))
    total_pairs = per_query.sum()
    if n_threads > 1:
        max_pairs = int(min(max_pairs, max(total_pairs / (4 * n_threads), 50_000)))
    # Dense accumulation (one cell per query x library spectrum) is linear time
    # and wins unless the library is much larger than the matches per query.
    dense = len(queries) * len(library) <= DENSE_FACTOR * max(total_pairs, 1)
    max_rows = max(1, DENSE_MAX_CELLS // max(len(library), 1)) if dense else None
    chunks = _chunk_queries(per_query, max_pairs, max_rows)
    self_map = (queries.source_idx, library.source_idx) if exclude_self else None

    def run(chunk):
        q0, q1 = chunk
        q, lib, vals = pairs_fn(library, queries, qspec, start, per_peak, q0, q1, ms2_tolerance)
        return _accumulate(q, lib, vals, q0, q1 - q0, len(library),
                           top_k, min_score, dense, self_map)

    if n_threads > 1:
        with ThreadPoolExecutor(n_threads) as ex:
            parts = list(ex.map(run, chunks))
    else:
        parts = [run(c) for c in chunks]

    q = np.concatenate([p[0] for p in parts]) if parts else np.zeros(0, np.int64)
    lib = np.concatenate([p[1] for p in parts]) if parts else np.zeros(0, np.int64)
    score = np.concatenate([p[2] for p in parts]) if parts else np.zeros(0, np.float32)
    n = np.concatenate([p[3] for p in parts]) if parts else np.zeros(0, np.int32)

    # Hits come grouped by query position, best first; reorder whole groups
    # into input order of the queries.
    first = np.flatnonzero(np.r_[True, q[1:] != q[:-1]]) if len(q) else np.zeros(0, np.int64)
    size = np.diff(np.r_[first, len(q)])
    groups = np.argsort(queries.source_idx[q[first]], kind='stable')
    order = expand_ranges(first[groups], size[groups])
    return SearchResult(queries.source_idx[q[order]], library.source_idx[lib[order]],
                        score[order], n[order])


# ----------------------------------------------------------------- candidates

def _open_candidates(library, queries, tol):
    """Library fragment range for every query peak: (start, count)."""
    lo32, hi32 = _mz_bounds(queries, tol)
    lo = np.searchsorted(library.frag_mz, lo32, side='left')
    hi = np.searchsorted(library.frag_mz, hi32, side='right')
    return lo, (hi - lo).astype(np.int64)


def _identity_candidates(library, queries, ms1_tol):
    """Candidate library peak range for every query peak, from precursor m/z: (start, count)."""
    prec = queries.precursor_mz
    s_lo = np.searchsorted(library.precursor_mz, prec - ms1_tol, side='left')
    s_hi = np.searchsorted(library.precursor_mz, prec + ms1_tol, side='right')
    s_hi[np.isnan(prec)] = s_lo[np.isnan(prec)]
    start = library.offsets[s_lo]
    n_cand = library.offsets[s_hi] - start
    qspec = queries.cached('spec_ids', lambda: spectrum_ids(queries.offsets))
    return start[qspec], n_cand[qspec]


def _open_pairs(library, queries, qspec, start, per_peak, q0, q1, tol):
    """All (query, library spectrum, contribution) for queries q0..q1, open search."""
    p0, p1 = queries.offsets[q0], queries.offsets[q1]
    counts = per_peak[p0:p1]
    lp = expand_ranges(start[p0:p1], counts)       # library windows are contiguous
    lb, lg = _half_and_g(library, 'frag_intensity')
    qa, qg = _half_and_g(queries, 'intensity')
    vals = _entropy_terms(np.repeat(qa[p0:p1], counts), np.repeat(qg[p0:p1], counts), lb[lp], lg[lp])
    return np.repeat(qspec[p0:p1], counts), library.frag_spec[lp], vals


def _identity_pairs(library, queries, qspec, start, per_peak, q0, q1, tol):
    """Same for identity search: candidate peaks of precursor-matched spectra, filtered by m/z."""
    p0, p1 = queries.offsets[q0], queries.offsets[q1]
    counts = per_peak[p0:p1]
    lp = expand_ranges(start[p0:p1], counts)
    lo32, hi32 = _mz_bounds(queries, tol)
    lmz = library.mz[lp]
    ok = (lmz >= np.repeat(lo32[p0:p1], counts)) & (lmz <= np.repeat(hi32[p0:p1], counts))
    qp = np.repeat(np.arange(p0, p1), counts)[ok]
    lp = lp[ok]
    lib = library.cached('spec_ids', lambda: spectrum_ids(library.offsets))[lp]
    lb, lg = _half_and_g(library, 'intensity')
    qa, qg = _half_and_g(queries, 'intensity')
    return qspec[qp], lib, _entropy_terms(qa[qp], qg[qp], lb[lp], lg[lp])


# ------------------------------------------------------------------- scoring

def _half_and_g(index, name):
    """Half intensities and g(x) = x log2 x of them, cached per index."""
    def make():
        x = getattr(index, name) / np.float32(2)
        return x, x * np.log2(x)
    return index.cached('half_' + name, make)


def _entropy_terms(a, ga, b, gb):
    ab = a + b
    return ab * np.log2(ab) - ga - gb


def _accumulate(q, lib, vals, q0, n_q, n_lib, top_k, min_score, dense, self_map):
    """
    Sum per-peak contributions into (query, library) scores for one chunk of
    queries q0 .. q0 + n_q, then select hits. Returns (query position,
    library position, score, n_matched), grouped by query, best first.
    """
    if len(vals) == 0:
        return (np.zeros(0, np.int64),) * 2 + (np.zeros(0, np.float32), np.zeros(0, np.int32))
    keys = (q - q0).astype(np.int64) * n_lib + lib
    if self_map is not None:
        q_src, lib_src = self_map
        keep = q_src[q] != lib_src[lib]
        keys, vals = keys[keep], vals[keep]

    if dense:
        n_cells = n_q * n_lib
        count = np.bincount(keys, minlength=n_cells)
        score = np.bincount(keys, weights=vals, minlength=n_cells)
        if top_k is not None and top_k < n_lib:
            hit = _topk_candidates(score.reshape(n_q, n_lib), top_k, min_score)
        else:
            hit = np.flatnonzero(count)
        score, count = score[hit], count[hit]
    else:
        hit, inv = np.unique(keys, return_inverse=True)
        score = np.bincount(inv, weights=vals)
        count = np.bincount(inv)

    keep = score >= min_score
    hit, score, count = hit[keep], score[keep].astype(np.float32), count[keep].astype(np.int32)
    q, lib = hit // n_lib, hit % n_lib
    # group by query, descending score (single-key sort; scores are in [0, 1])
    order = np.argsort(q * 4.0 - score, kind='stable')
    q, lib, score, count = q[order] + q0, lib[order], score[order], count[order]
    if top_k is not None:
        first = np.flatnonzero(np.r_[True, q[1:] != q[:-1]])
        rank = np.arange(len(q)) - np.repeat(first, np.diff(np.r_[first, len(q)]))
        keep = rank < top_k
        q, lib, score, count = q[keep], lib[keep], score[keep], count[keep]
    return q, lib, score, count


def _topk_candidates(S, k, min_score):
    """
    Flat indices of cells of the dense score block S (queries x library, 0 for
    no match) that include every row's k best scores >= min_score, plus few
    others. A per-row threshold, lowered from the row maximum until at least k
    cells pass, replaces a full per-row partition: each step is one fast
    comparison pass, and most rows are settled at the first step.
    """
    rowmax = S.max(axis=1)
    thr = np.full(len(S), np.inf)
    todo = np.flatnonzero(rowmax >= max(min_score, np.finfo(np.float64).tiny))
    for f in (0.5, 0.1, 0.01, 0.001, 0.0):
        if len(todo) == 0:
            break
        t = np.maximum(rowmax[todo] * f, max(min_score, np.finfo(np.float64).tiny))
        rows = S if len(todo) == len(S) else S[todo]
        done = (np.count_nonzero(rows >= t[:, None], axis=1) >= k) | (t <= min_score) | (f == 0.0)
        thr[todo[done]] = t[done]
        todo = todo[~done]
    return np.flatnonzero(S >= thr[:, None])


# -------------------------------------------------------------------- helpers

def _as_queries(queries, library):
    if isinstance(queries, SpectralIndex):
        return queries
    return SpectralIndex.build(queries, fragment_index=False, keep_metadata=False, **library.params)


def _check_compatible(library, queries, ms2_tolerance):
    lp, qp = library.params, queries.params
    if lp['intensity_weight'] != qp['intensity_weight']:
        raise ValueError('library and queries use different intensity_weight')
    for name in ('noise_threshold', 'precursor_removal_da', 'min_peak_distance', 'max_peak_num'):
        if lp[name] != qp[name]:
            warnings.warn(f'library and queries were cleaned with different {name}: {lp[name]} vs {qp[name]}')
    if 2 * ms2_tolerance >= min(lp['min_peak_distance'], qp['min_peak_distance']):
        raise ValueError('ms2_tolerance must be < min_peak_distance / 2 used in cleaning, '
                         'so that each peak matches at most one peak per spectrum')


def _mz_bounds(queries, tol):
    """
    float32 bounds [lo, hi] such that, for float32 library m/z x,
    lo <= x <= hi  exactly when  mz - tol <= x <= mz + tol  in float64.
    Keeps the library arrays float32 without changing tolerance edges.
    """
    def make():
        mz = queries.mz.astype(np.float64)
        lo, hi = mz - tol, mz + tol
        lo32, hi32 = lo.astype(np.float32), hi.astype(np.float32)
        lo32 = np.where(lo32 < lo, np.nextafter(lo32, np.float32(np.inf)), lo32)
        hi32 = np.where(hi32 > hi, np.nextafter(hi32, np.float32(-np.inf)), hi32)
        return lo32, hi32
    return queries.cached(('bounds', tol), make)


def _chunk_queries(per_query, max_pairs, max_rows=None):
    """
    Split queries into contiguous chunks of at most about max_pairs candidate
    peak pairs and at most max_rows queries.
    """
    n = len(per_query)
    if n == 0:
        return []
    cum = np.cumsum(per_query)
    chunks, q0 = [], 0
    while q0 < n:
        base = cum[q0 - 1] if q0 else 0
        q1 = max(int(np.searchsorted(cum, base + max_pairs, side='right')), q0 + 1)
        if max_rows is not None:
            q1 = min(q1, q0 + max_rows)
        q1 = min(q1, n)
        chunks.append((q0, q1))
        q0 = q1
    return chunks
