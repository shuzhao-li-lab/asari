"""
Batch preprocessing of MS/MS spectra, vectorized over all spectra at once.

Spectra are held in a packed (CSR-like) layout:
    mz, intensity : float32 arrays of all peaks, sorted by m/z within each spectrum
    offsets       : int64 array of length n_spectra + 1; spectrum i owns
                    peaks offsets[i]:offsets[i+1]

The cleaning steps follow ms_entropy.clean_spectrum / FlashEntropySearch
(Li et al. 2021 Nat Methods 18:1524; Li & Fiehn 2023 Nat Methods 20:1475),
so that similarity scores are directly comparable.
"""

import numpy as np


def spectrum_ids(offsets):
    """Spectrum index of every peak in a packed layout."""
    return np.repeat(np.arange(len(offsets) - 1), np.diff(offsets))


def offsets_from_ids(spec, n_spectra):
    """Offsets (n_spectra + 1) from sorted spectrum ids of peaks."""
    offsets = np.zeros(n_spectra + 1, dtype=np.int64)
    np.cumsum(np.bincount(spec, minlength=n_spectra), out=offsets[1:])
    return offsets


def expand_ranges(starts, lengths):
    """Concatenate arange(s, s + n) for all (s, n), vectorized."""
    lengths = np.asarray(lengths, dtype=np.int64)
    ends = np.cumsum(lengths)
    return np.repeat(np.asarray(starts, dtype=np.int64) - (ends - lengths), lengths) + np.arange(ends[-1] if len(ends) else 0)


def pack_peaks(peaks_list):
    """
    Pack a list of peak arrays, each of shape (n, 2) as [m/z, intensity],
    into flat float32 arrays. Returns mz, intensity, spectrum id per peak.
    """
    arrays = [np.asarray(p, dtype=np.float32).reshape(-1, 2) for p in peaks_list]
    lengths = np.array([len(a) for a in arrays], dtype=np.int64)
    flat = np.concatenate(arrays) if arrays else np.zeros((0, 2), dtype=np.float32)
    spec = np.repeat(np.arange(len(arrays)), lengths)
    return flat[:, 0].copy(), flat[:, 1].copy(), spec


def clean_spectra(peaks_list,
                  precursor_mz=None,
                  precursor_removal_da=1.6,
                  noise_threshold=0.01,
                  min_peak_distance=0.05,
                  max_peak_num=None):
    """
    Clean a batch of spectra:

        1. remove peaks with m/z <= 0 or intensity <= 0;
        2. remove peaks with m/z >= precursor_mz - precursor_removal_da
           (skipped if precursor_mz or precursor_removal_da is None);
        3. centroid: merge peaks closer than min_peak_distance, greedily from
           the most intense peak (only for spectra that need it);
        4. remove peaks below noise_threshold * base peak intensity;
        5. keep the max_peak_num most intense peaks (if given);
        6. normalize intensities to sum to 1.

    min_peak_distance must be > 2 * the MS2 tolerance used in search, so that
    each peak matches at most one peak of another spectrum.

    Returns (mz, intensity, offsets) in packed layout. Spectra left with no
    peaks keep an empty slot, so indices still refer to the input list.
    """
    n = len(peaks_list)
    mz, intensity, spec = pack_peaks(peaks_list)

    keep = (mz > 0) & (intensity > 0)
    if precursor_mz is not None and precursor_removal_da is not None:
        max_mz = np.asarray(precursor_mz, dtype=np.float64) - precursor_removal_da
        limit = max_mz[spec]
        keep &= ~(limit > 0) | (mz < limit)          # no limit if precursor is NaN
    mz, intensity, spec = mz[keep], intensity[keep], spec[keep]

    if np.any((spec[1:] == spec[:-1]) & (np.diff(mz) < 0)):
        order = np.lexsort((mz, spec))
        mz, intensity, spec = mz[order], intensity[order], spec[order]

    mz, intensity, spec = _centroid_where_needed(mz, intensity, spec, min_peak_distance)

    if noise_threshold is not None and len(mz):
        base = _segment_max(intensity, spec, n)
        # float32 arithmetic, as in ms_entropy, so that peaks stored at exactly
        # the threshold (e.g. 0.01 -> 0.0099999998) are kept
        keep = intensity >= np.float32(noise_threshold) * base[spec]
        mz, intensity, spec = mz[keep], intensity[keep], spec[keep]

    if max_peak_num is not None and max_peak_num > 0 and len(mz):
        mz, intensity, spec = _top_peaks(mz, intensity, spec, n, max_peak_num)

    totals = np.bincount(spec, weights=intensity, minlength=n)
    intensity = (intensity / totals[spec]).astype(np.float32)
    return mz, intensity, offsets_from_ids(spec, n)


def entropy_weight(intensity, offsets):
    """
    Entropy-based intensity weighting (Li et al. 2021), vectorized over spectra.
    Spectra with entropy S < 3 have intensities raised to 0.25 + 0.25 * S and
    renormalized; others are unchanged. Input intensities must sum to 1 per
    spectrum. Returns a new float32 array.
    """
    spec = spectrum_ids(offsets)
    n = len(offsets) - 1
    x = intensity.astype(np.float64)
    entropy = -np.bincount(spec, weights=x * np.log(x), minlength=n)
    low = entropy < 3
    power = np.where(low, 0.25 + 0.25 * entropy, 1.0)
    x = x ** power[spec]
    totals = np.bincount(spec, weights=x, minlength=n)
    x = np.where(low[spec], x / totals[spec], x)
    return x.astype(np.float32)


def _segment_max(values, spec, n):
    out = np.zeros(n, dtype=values.dtype)
    if len(values):
        starts = np.flatnonzero(np.r_[True, spec[1:] != spec[:-1]])
        out[spec[starts]] = np.maximum.reduceat(values, starts)
    return out


def _top_peaks(mz, intensity, spec, n, k):
    order = np.lexsort((-intensity, spec))
    starts = offsets_from_ids(spec, n)[:-1]
    rank = np.arange(len(order)) - starts[spec[order]]
    kept = order[rank < k]
    kept = kept[np.lexsort((mz[kept], spec[kept]))]
    return mz[kept], intensity[kept], spec[kept]


def _centroid_where_needed(mz, intensity, spec, min_distance):
    """
    Centroid, vectorized over all spectra, the spectra that have peaks closer
    than min_distance. Gives the same result as centroid_spectrum() per spectrum.
    Input must be sorted by (spec, mz).

    Only peaks with a neighbor within min_distance take part; isolated peaks
    stay as they are. Merged peaks are written in place of their center, so
    the arrays stay sorted without a new sort in the usual case.
    """
    mz, intensity = mz.astype(np.float64), intensity.astype(np.float64)
    while len(mz) > 1:
        same = spec[1:] == spec[:-1]
        gap = np.diff(mz)
        too_close = same & (gap < min_distance)
        if not too_close.any():
            break
        needs = np.zeros(spec[-1] + 1, dtype=bool)
        needs[spec[1:][too_close]] = True
        linked = same & (gap <= min_distance) & needs[spec[1:]]
        member = np.zeros(len(mz), dtype=bool)
        member[:-1] |= linked
        member[1:] |= linked
        sub = np.flatnonzero(member)

        owner = _greedy_owners(mz[sub], intensity[sub], spec[sub], min_distance)
        w = intensity[sub]
        total = np.bincount(owner, weights=w, minlength=len(sub))
        weighted_mz = np.bincount(owner, weights=w * mz[sub], minlength=len(sub))
        is_center = owner == np.arange(len(sub))
        c = sub[is_center]
        mz[c] = weighted_mz[is_center] / total[is_center]
        intensity[c] = total[is_center]
        keep = np.ones(len(mz), dtype=bool)
        keep[sub[~is_center]] = False
        mz, intensity, spec = mz[keep], intensity[keep], spec[keep]
        if np.any((spec[1:] == spec[:-1]) & (np.diff(mz) < 0)):
            order = np.lexsort((mz, spec))
            mz, intensity, spec = mz[order], intensity[order], spec[order]
    return mz.astype(np.float32), intensity.astype(np.float32), spec


def _greedy_owners(mz, intensity, spec, min_distance):
    """
    One centroiding pass over peaks sorted by (spec, mz), all of which have a
    neighbor within min_distance. Returns, for each peak, the index of the
    center that absorbs it (itself for centers).

    A pass of centroid_spectrum() visits peaks from most to least intense; a
    peak still unmerged at its turn becomes a center and absorbs all unmerged
    peaks within min_distance. The centers are thus the greedy maximal
    independent set, by intensity, of the graph linking peaks within
    min_distance, found exactly in parallel rounds: every undecided peak more
    intense than all its undecided neighbors becomes a center, and its
    neighbors are removed. Each other peak is absorbed by the most intense
    center within min_distance, the first one visited.
    """
    n = len(mz)
    pairs = []                    # (k, near): peaks i and i + k within min_distance
    k = 1
    while k < n:
        near = (spec[k:] == spec[:-k]) & (mz[k:] - mz[:-k] <= min_distance)
        if not near.any():
            break
        pairs.append((k, near))
        k += 1

    rank = np.empty(n, dtype=np.int64)         # higher rank = visited earlier
    rank[np.argsort(intensity, kind='stable')] = np.arange(n)
    undecided = np.ones(n, dtype=bool)
    center = np.zeros(n, dtype=bool)
    while undecided.any():
        best = np.full(n, -1, dtype=np.int64)  # best undecided neighbor rank
        for k, near in pairs:
            np.maximum(best[:-k], np.where(near & undecided[k:], rank[k:], -1), out=best[:-k])
            np.maximum(best[k:], np.where(near & undecided[:-k], rank[:-k], -1), out=best[k:])
        new = undecided & (rank > best)
        center |= new
        undecided &= ~new
        for k, near in pairs:
            undecided[k:] &= ~(near & new[:-k])
            undecided[:-k] &= ~(near & new[k:])

    owner = np.where(center, np.arange(n), -1)
    owner_rank = np.where(center, rank, -1)
    for k, near in pairs:
        take = near & center[k:] & (rank[k:] > owner_rank[:-k])
        owner[:-k][take] = np.flatnonzero(take) + k
        owner_rank[:-k][take] = rank[k:][take]
        take = near & center[:-k] & (rank[:-k] > owner_rank[k:])
        owner[k:][take] = np.flatnonzero(take)
        owner_rank[k:][take] = rank[:-k][take]
    return owner


def centroid_spectrum(mz, intensity, min_distance):
    """
    Centroid one m/z-sorted spectrum as ms_entropy does: from the most intense
    peak down, merge all peaks within min_distance into one peak (summed
    intensity, intensity-weighted m/z); repeat passes until all peaks are at
    least min_distance apart.
    """
    mz = mz.astype(np.float64)
    intensity = intensity.astype(np.float64)
    while len(mz) > 1 and np.any(np.diff(mz) < min_distance):
        remaining = intensity.copy()
        out_mz, out_int = [], []
        for idx in np.argsort(intensity)[::-1]:
            if remaining[idx] <= 0:
                continue
            left = np.searchsorted(mz, mz[idx] - min_distance, side='left')
            right = np.searchsorted(mz, mz[idx] + min_distance, side='right')
            w = remaining[left:right]
            total = w.sum()
            out_mz.append(np.dot(mz[left:right], w) / total)
            out_int.append(total)
            remaining[left:right] = 0
        out_mz, out_int = np.array(out_mz), np.array(out_int)
        order = np.argsort(out_mz)
        mz, intensity = out_mz[order], out_int[order]
    return mz.astype(np.float32), intensity.astype(np.float32)
