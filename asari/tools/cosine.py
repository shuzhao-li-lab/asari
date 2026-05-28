"""
Four spectral similarity algorithms for MS2 spectra comparison.

See evaluation notebook for benchmarking and discussion:
`notebooks/evaluate_cosine_scoring_functions.ipynb` 
"""

import numpy as np
from scipy.optimize import linear_sum_assignment


def cosine_similarity(spec1, spec2, mz_tolerance=0.05):
    """
    The function accepts two spectra as numpy arrays of shape (n, 2),
    where column 0 is m/z and column 1 is intensity.
    The score is a float in [0, 1], where 1.0 is a perfect match.
    Returns score, num_matched_peaks
    
    Greedy cosine similarity: peaks are matched in descending order of
    intensity product, each peak used at most once.
    This is simple_cosine, the most common cosine variant in MS/MS database search
    (equivalent to CosineGreedy in matchms).
    """
    s1 = np.array(spec1, dtype=float)
    s2 = np.array(spec2, dtype=float)

    norm1 = np.linalg.norm(s1[:, 1])
    norm2 = np.linalg.norm(s2[:, 1])
    if norm1 == 0 or norm2 == 0:
        return 0.0
    s1[:, 1] /= norm1
    s2[:, 1] /= norm2

    # Collect all candidate pairs within mz_tolerance, rank by intensity product
    candidates = []
    for i in range(len(s1)):
        for j in range(len(s2)):
            if abs(s1[i, 0] - s2[j, 0]) <= mz_tolerance:
                candidates.append((s1[i, 1] * s2[j, 1], i, j))
    candidates.sort(reverse=True)

    used1, used2 = set(), set()
    score = 0.0
    for prod, i, j in candidates:
        if i not in used1 and j not in used2:
            score += prod
            used1.add(i)
            used2.add(j)

    return min(score, 1.0), len(used1)


def hungarian_cosine(spec1, spec2, mz_tolerance=0.05):
    """
    Cosine similarity with optimal peak assignment via the Hungarian algorithm
    (linear_sum_assignment). Finds the globally optimal one-to-one peak
    matching rather than a greedy approximation.

    Equivalent to CosineHungarian in matchms; slower than simple_cosine
    but theoretically optimal for the bipartite matching problem.
    """
    s1 = np.array(spec1, dtype=float)
    s2 = np.array(spec2, dtype=float)

    norm1 = np.linalg.norm(s1[:, 1])
    norm2 = np.linalg.norm(s2[:, 1])
    if norm1 == 0 or norm2 == 0:
        return 0.0
    s1[:, 1] /= norm1
    s2[:, 1] /= norm2

    n1, n2 = len(s1), len(s2)
    # Cost matrix: 0 where peaks are outside tolerance (cannot be matched),
    # negative intensity product elsewhere (linear_sum_assignment minimizes).
    cost = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            if abs(s1[i, 0] - s2[j, 0]) <= mz_tolerance:
                cost[i, j] = -s1[i, 1] * s2[j, 1]

    row_idx, col_idx = linear_sum_assignment(cost)
    score = -cost[row_idx, col_idx].sum()

    return min(score, 1.0)


def dot_product(spec1, spec2, mz_tolerance=0.05, mz_power=2.0, int_power=0.5):
    """
    Weighted dot product similarity (NIST / AMDIS convention).

    Each peak is weighted by  mz^mz_power * intensity^int_power
    before normalization and matching.  The defaults (mz^2 * sqrt(intensity))
    reduce the influence of chemical noise and up-weight high-mass fragments,
    which improves discrimination for GC-HRMS library search.

    Greedy matching (same strategy as simple_cosine) is used after weighting.
    Set mz_power=0, int_power=1 to recover an intensity-only cosine.
    """
    s1 = np.array(spec1, dtype=float)
    s2 = np.array(spec2, dtype=float)

    w1 = (s1[:, 0] ** mz_power) * (s1[:, 1] ** int_power)
    w2 = (s2[:, 0] ** mz_power) * (s2[:, 1] ** int_power)

    norm1 = np.linalg.norm(w1)
    norm2 = np.linalg.norm(w2)
    if norm1 == 0 or norm2 == 0:
        return 0.0
    w1 /= norm1
    w2 /= norm2

    # Greedy matching on weighted intensities
    candidates = []
    for i in range(len(s1)):
        for j in range(len(s2)):
            if abs(s1[i, 0] - s2[j, 0]) <= mz_tolerance:
                candidates.append((w1[i] * w2[j], i, j))
    candidates.sort(reverse=True)

    used1, used2 = set(), set()
    score = 0.0
    for prod, i, j in candidates:
        if i not in used1 and j not in used2:
            score += prod
            used1.add(i)
            used2.add(j)

    return min(score, 1.0)


def linear_cosine(spec1, spec2, mz_tolerance=0.05):
    """
    Cosine similarity with O(n+m) two-pointer peak matching on mz-sorted spectra.

    Both spectra are sorted by m/z and scanned with two pointers.  At each
    step the pointer on the lower m/z side advances, and a match is recorded
    when both peaks fall within mz_tolerance.  This is the fastest matching
    strategy and works well when spectra are already sorted (typical for
    centroided data from asari / ms_entropy pipelines).

    Note: unlike simple_cosine, the first valid match per peak is taken rather
    than the highest-scoring one, so scores can differ slightly for spectra
    with overlapping m/z clusters.
    """
    s1 = np.array(spec1, dtype=float)
    s2 = np.array(spec2, dtype=float)

    # Sort by m/z
    s1 = s1[np.argsort(s1[:, 0])]
    s2 = s2[np.argsort(s2[:, 0])]

    norm1 = np.linalg.norm(s1[:, 1])
    norm2 = np.linalg.norm(s2[:, 1])
    if norm1 == 0 or norm2 == 0:
        return 0.0
    s1[:, 1] /= norm1
    s2[:, 1] /= norm2

    score = 0.0
    i, j = 0, 0
    while i < len(s1) and j < len(s2):
        delta = s1[i, 0] - s2[j, 0]
        if abs(delta) <= mz_tolerance:
            score += s1[i, 1] * s2[j, 1]
            i += 1
            j += 1
        elif delta < 0:
            i += 1
        else:
            j += 1

    return min(score, 1.0)
