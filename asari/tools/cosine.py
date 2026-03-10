"""
From https://github.com/Philipbear/reverse_search/blob/main/reverse_spectral_search/cosine.py

This file contains code modified from the matchms project
(https://github.com/matchms/matchms)
Copyright matchms Team 2020

Modified by Shipei Xing in 2024

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from typing import Tuple

import numba as nb
import numpy as np


@nb.njit
def find_matches(ref_spec_mz: np.ndarray, qry_spec_mz: np.ndarray,
                 tolerance: float, shift: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """Find matching peaks between two spectra."""
    matches_idx1 = np.empty(len(ref_spec_mz) * len(qry_spec_mz), dtype=np.int64)
    matches_idx2 = np.empty_like(matches_idx1)
    match_count = 0
    lowest_idx = 0

    for peak1_idx in range(len(ref_spec_mz)):
        mz = ref_spec_mz[peak1_idx]
        low_bound = mz - tolerance
        high_bound = mz + tolerance

        for peak2_idx in range(lowest_idx, len(qry_spec_mz)):
            mz2 = qry_spec_mz[peak2_idx] - shift
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                lowest_idx = peak2_idx
            else:
                matches_idx1[match_count] = peak1_idx
                matches_idx2[match_count] = peak2_idx
                match_count += 1

    return matches_idx1[:match_count], matches_idx2[:match_count]


@nb.njit
def collect_peak_pairs(ref_spec: np.ndarray, qry_spec: np.ndarray, min_matched_peak: int, sqrt_transform: bool,
                       tolerance: float, shift: float = 0.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Find and score matching peak pairs between spectra."""

    if len(ref_spec) == 0 or len(qry_spec) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float32)

    # Exact matching
    matches_idx1, matches_idx2 = find_matches(ref_spec[:, 0], qry_spec[:, 0], tolerance, 0.0)

    # If shift is not 0, perform hybrid search
    if abs(shift) > 1e-6:
        matches_idx1_shift, matches_idx2_shift = find_matches(ref_spec[:, 0], qry_spec[:, 0], tolerance, shift)
        matches_idx1 = np.concatenate((matches_idx1, matches_idx1_shift))
        matches_idx2 = np.concatenate((matches_idx2, matches_idx2_shift))

    if len(matches_idx1) < min_matched_peak:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float32)

    # Calculate scores for matches
    if sqrt_transform:
        scores = np.sqrt(ref_spec[matches_idx1, 1] * qry_spec[matches_idx2, 1]).astype(np.float32)
    else:
        scores = (ref_spec[matches_idx1, 1] * qry_spec[matches_idx2, 1]).astype(np.float32)

    # Sort by score descending
    sort_idx = np.argsort(-scores)
    return matches_idx1[sort_idx], matches_idx2[sort_idx], scores[sort_idx]


@nb.njit
def score_matches(matches_idx1: np.ndarray, matches_idx2: np.ndarray,
                  scores: np.ndarray, ref_spec: np.ndarray, qry_spec: np.ndarray,
                  sqrt_transform: bool, penalty: float):
    """Calculate final similarity score from matching peaks."""

    # Use boolean arrays for tracking used peaks - initialized to False
    used1 = np.zeros(len(ref_spec), dtype=nb.boolean)
    used2 = np.zeros(len(qry_spec), dtype=nb.boolean)

    total_score = 0.0
    used_matches = 0

    # Find best non-overlapping matches
    for i in range(len(matches_idx1)):
        idx1 = matches_idx1[i]
        idx2 = matches_idx2[i]
        if not used1[idx1] and not used2[idx2]:
            total_score += scores[i]
            used1[idx1] = True
            used2[idx2] = True
            used_matches += 1

    if used_matches == 0:
        return 0.0, 0

    # Calculate normalization factors
    if sqrt_transform:
        norm1 = np.sqrt(np.sum(np.sqrt(ref_spec[:, 1] * ref_spec[:, 1])))
    else:
        norm1 = np.sqrt(np.sum(ref_spec[:, 1] * ref_spec[:, 1]))

    # Sum intensities of matched peaks
    matched_intensities = np.zeros(used_matches, dtype=np.float32)
    # new intensities of qry peaks, matched peaks are the same, others are penalized
    new_qry_intensities = np.zeros(len(qry_spec), dtype=np.float32)
    match_idx = 0
    for i in range(len(qry_spec)):
        if used2[i]:
            matched_intensities[match_idx] = qry_spec[i, 1]
            new_qry_intensities[i] = qry_spec[i, 1]
            match_idx += 1
        else:
            new_qry_intensities[i] = qry_spec[i, 1] * (1 - penalty)

    if sqrt_transform:
        norm2 = np.sqrt(np.sum(np.sqrt(new_qry_intensities * new_qry_intensities)))
    else:
        norm2 = np.sqrt(np.sum(new_qry_intensities * new_qry_intensities))

    if norm1 == 0.0 or norm2 == 0.0:
        return 0.0, used_matches

    score = total_score / (norm1 * norm2)

    return min(float(score), 1.0), used_matches


def cosine_similarity(qry_spec: np.ndarray, ref_spec: np.ndarray,
                      tolerance: float = 0.1,
                      min_matched_peak: int = 1,
                      sqrt_transform: bool = True,
                      penalty: float = 0.,
                      shift: float = 0.0):
    """
    Calculate similarity between two spectra.

    Parameters
    ----------
    qry_spec: np.ndarray
        Query spectrum.
    ref_spec: np.ndarray
        Reference spectrum.
    tolerance: float
        Tolerance for m/z matching.
    min_matched_peak: int
        Minimum number of matched peaks.
    sqrt_transform: bool
        If True, use square root transformation.
    penalty: float
        Penalty for unmatched peaks. If set to 0, traditional cosine score; if set to 1, traditional reverse cosine score.
    shift: float
        Shift for m/z values. If not 0, hybrid search is performed. shift = prec_mz(qry) - prec_mz(ref)
    """
    tolerance = np.float32(tolerance)
    penalty = np.float32(penalty)
    shift = np.float32(shift)

    if qry_spec.size == 0 or ref_spec.size == 0:
        return 0.0, 0

    # normalize the intensity
    ref_spec[:, 1] /= np.max(ref_spec[:, 1])
    qry_spec[:, 1] /= np.max(qry_spec[:, 1])

    matches_idx1, matches_idx2, scores = collect_peak_pairs(
        ref_spec, qry_spec, min_matched_peak, sqrt_transform,
        tolerance, shift
    )

    if len(matches_idx1) == 0:
        return 0.0, 0

    return score_matches(
        matches_idx1, matches_idx2, scores,
        ref_spec, qry_spec, sqrt_transform, penalty
    )


if __name__ == "__main__":

    # Example usage
    peaks1 = np.array([[50, 8.0], [70, 100.0], [80, 50.0], [100, 50.0]], dtype=np.float32)

    peaks2 = np.array([[55, 38.0], [80, 66.0], [90, 999.0]], dtype=np.float32)

    # Example with standard cosine
    score, n_matches = cosine_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=0)
    print(f"Standard Score: {score:.3f}, Matches: {n_matches}")

    # Example with enhanced reverse cosine
    score, n_matches = cosine_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=0.6)
    print(f"Reverse Score: {score:.3f}, Matches: {n_matches}")

    # Example with traditional reverse cosine
    score, n_matches = cosine_similarity(peaks1, peaks2, tolerance=0.05, sqrt_transform=True, penalty=1)
    print(f"Reverse Score: {score:.3f}, Matches: {n_matches}")
