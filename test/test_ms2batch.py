"""
Tests for asari.ms2batch. Run with pytest, unittest, or directly:  PYTHONPATH=. python3 test/test_ms2batch.py
"""

import io
import unittest
import os
import tempfile
import zipfile

import numpy as np

from asari.ms2batch import SpectralIndex, search, read_msp
from asari.ms2batch import engine
from asari.ms2batch.preprocess import clean_spectra, centroid_spectrum, entropy_weight

rng = np.random.default_rng(1)


def random_spectrum(n_clusters=None, center=None):
    n_clusters = n_clusters or rng.integers(2, 25)
    centers = rng.uniform(50, 400, n_clusters) if center is None else center
    mz = np.concatenate([c + rng.normal(0, 0.02, rng.integers(1, 5)) for c in centers])
    return np.column_stack([mz, rng.exponential(1, len(mz)) ** 3]).astype(np.float32)


def random_library(n, shared=None):
    out = []
    for i in range(n):
        p = random_spectrum()
        if shared is not None and i % 3 == 0:     # some spectra share fragments
            p = np.vstack([p, shared + [rng.normal(0, 0.005), 0]])
        out.append({'precursor_mz': float(rng.uniform(300, 320)), 'peaks': p, 'name': f's{i}'})
    return out


def reference_similarity(a, b, tol):
    """Brute-force entropy similarity of two cleaned, weighted spectra (n, 2)."""
    score = 0.0
    used = set()
    for mz_a, ia in a:
        for j, (mz_b, ib) in enumerate(b):
            if j not in used and abs(float(mz_a) - float(mz_b)) <= tol:
                x, y = ia / 2, ib / 2
                score += (x + y) * np.log2(x + y) - x * np.log2(x) - y * np.log2(y)
                used.add(j)
                break
    return score


def dense_reference(L, Q, method, ms1_tol=0.01, ms2_tol=0.02):
    out = np.zeros((Q.header['n_input_spectra'], L.header['n_input_spectra']))
    for qi in range(len(Q)):
        for li in range(len(L)):
            if method == 'identity' and abs(Q.precursor_mz[qi] - L.precursor_mz[li]) > ms1_tol:
                continue
            out[Q.source_idx[qi], L.source_idx[li]] = reference_similarity(Q.peaks(qi), L.peaks(li), ms2_tol)
    return out


# ------------------------------------------------------------------ tests

def test_centroid_vectorized_matches_loop():
    spectra = [random_spectrum() for _ in range(300)]
    mz, it, off = clean_spectra(spectra, noise_threshold=None)
    for i, p in enumerate(spectra):
        p = p[np.argsort(p[:, 0])]
        m, x = centroid_spectrum(p[:, 0], p[:, 1], 0.05)
        a, b = off[i], off[i + 1]
        assert len(m) == b - a
        np.testing.assert_allclose(m, mz[a:b], atol=1e-4)
        np.testing.assert_allclose(x / x.sum(), it[a:b], rtol=1e-4)
        assert np.all(np.diff(mz[a:b]) >= 0.05)


def test_cleaning_steps():
    peaks = np.array([[50, 1.0], [80, 0.005], [100, 10.0], [100.02, 5.0], [199, 3.0], [0, 5], [120, 0]], np.float32)
    mz, it, off = clean_spectra([peaks], precursor_mz=[200.0])
    # 80 below 1% noise, 100/100.02 merged, 199 >= 200 - 1.6 removed, zero peaks removed
    np.testing.assert_allclose(mz, [50, 100.00667], atol=1e-4)
    np.testing.assert_allclose(it, [1 / 16, 15 / 16], rtol=1e-6)
    np.testing.assert_array_equal(off, [0, 2])


def test_entropy_weight():
    intensity = np.array([0.7, 0.2, 0.1, 1.0], np.float32)
    w = entropy_weight(intensity, np.array([0, 3, 4]))
    s = -np.sum(intensity[:3] * np.log(intensity[:3]))
    expected = intensity[:3] ** (0.25 + 0.25 * s)
    np.testing.assert_allclose(w[:3], expected / expected.sum(), rtol=1e-6)
    assert w[3] == 1.0


def test_scores_match_bruteforce():
    shared = random_spectrum(8)
    lib = random_library(60, shared)
    queries = random_library(20, shared)
    L = SpectralIndex.build(lib)
    Q = SpectralIndex.build(queries, fragment_index=False)
    for method in ('identity', 'open'):
        ref = dense_reference(L, Q, method)
        got = search(L, Q, method=method).to_dense(*ref.shape)
        np.testing.assert_allclose(got, ref, atol=1e-6)
        assert np.array_equal(got > 0, ref > 0)


def test_topk_min_score_and_order():
    shared = random_spectrum(8)
    L = SpectralIndex.build(random_library(80, shared))
    Q = SpectralIndex.build(random_library(15, shared), fragment_index=False)
    full = search(L, Q, method='open')
    top = search(L, Q, method='open', top_k=3, min_score=0.05)
    for q in range(15):
        s_full = np.sort(full.score[(full.query_idx == q) & (full.score >= 0.05)])[::-1][:3]
        s_top = top.score[top.query_idx == q]
        np.testing.assert_allclose(s_top, s_full, atol=1e-7)
    assert np.all(np.diff(top.query_idx) >= 0)
    for q in np.unique(top.query_idx):
        assert np.all(np.diff(top.score[top.query_idx == q]) <= 0)


def test_dense_sparse_and_threads_agree():
    shared = random_spectrum(8)
    L = SpectralIndex.build(random_library(100, shared))
    Q = SpectralIndex.build(random_library(40, shared), fragment_index=False)
    old = engine.DENSE_FACTOR
    try:
        results = []
        for factor, threads, max_pairs in ((16, 1, 10_000_000), (0, 1, 10_000_000), (16, 3, 500), (0, 2, 500)):
            engine.DENSE_FACTOR = factor
            h = search(L, Q, method='open', top_k=4, n_threads=threads, max_pairs=max_pairs)
            results.append(h)
    finally:
        engine.DENSE_FACTOR = old
    for h in results[1:]:
        np.testing.assert_array_equal(h.query_idx, results[0].query_idx)
        np.testing.assert_allclose(h.score, results[0].score, atol=1e-7)
        np.testing.assert_array_equal(h.n_matched, results[0].n_matched)


def test_identity_is_subset_of_open():
    lib = random_library(50)
    L = SpectralIndex.build(lib)
    ident = search(L, lib[:10], method='identity', ms1_tolerance=2.0)
    op = search(L, lib[:10], method='open')
    o = {(q, l): s for q, l, s in zip(op.query_idx, op.library_idx, op.score)}
    for q, l, s in zip(ident.query_idx, ident.library_idx, ident.score):
        assert abs(o[(q, l)] - s) < 1e-7


def test_self_search_excludes_self():
    lib = random_library(30)
    L = SpectralIndex.build(lib)
    h = search(L, L, method='open', top_k=3)
    top1 = h.library_idx[np.r_[True, h.query_idx[1:] != h.query_idx[:-1]]]
    np.testing.assert_array_equal(top1, np.arange(30))       # each spectrum finds itself
    np.testing.assert_allclose(h.score[np.r_[True, h.query_idx[1:] != h.query_idx[:-1]]], 1.0, atol=1e-5)
    h = search(L, L, method='open', top_k=3, exclude_self=True)
    assert not np.any(h.query_idx == h.library_idx)
    assert np.all(np.bincount(h.query_idx) <= 3)


def test_save_load_roundtrip():
    lib = random_library(40)
    L = SpectralIndex.build(lib)
    Q = SpectralIndex.build(lib[:10], fragment_index=False)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'lib.npz')
        L.save(path)
        with np.load(path, allow_pickle=False) as z:            # plain numpy can read it
            assert set(z.files) >= {'header', 'mz', 'frag_mz', 'precursor_mz'}
        for mmap in (False, True):
            L2 = SpectralIndex.load(path, mmap=mmap)
            for name in ('precursor_mz', 'offsets', 'mz', 'intensity', 'frag_mz', 'frag_spec', 'source_idx'):
                np.testing.assert_array_equal(getattr(L2, name), getattr(L, name))
            assert L2.metadata(7)['name'] == 's7'
            a, b = search(L, Q, method='open'), search(L2, Q, method='open')
            np.testing.assert_array_equal(a.score, b.score)


def test_load_rejects_pickle_and_foreign_files():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'bad.npz')
        np.savez(path, header=np.frombuffer(b'{"format": "other"}', np.uint8), x=np.zeros(3))
        try:
            SpectralIndex.load(path)
            assert False, 'foreign file accepted'
        except ValueError:
            pass
        # an object array inside the archive must not be unpickled
        L = SpectralIndex.build(random_library(5))
        good = os.path.join(d, 'good.npz')
        L.save(good)
        evil = os.path.join(d, 'evil.npz')
        with zipfile.ZipFile(good) as zin, zipfile.ZipFile(evil, 'w') as zout:
            for item in zin.infolist():
                data = zin.read(item.filename)
                if item.filename == 'mz.npy':
                    buf = io.BytesIO()
                    np.save(buf, np.array([{'a': 1}], dtype=object), allow_pickle=True)
                    data = buf.getvalue()
                zout.writestr(item, data)
        for mmap in (False, True):
            try:
                SpectralIndex.load(evil, mmap=mmap)
                assert False, 'object array accepted'
            except ValueError:
                pass


def test_empty_and_missing_precursor():
    lib = [{'precursor_mz': None, 'peaks': random_spectrum()} for _ in range(5)]
    lib.append({'precursor_mz': 100.0, 'peaks': np.zeros((0, 2))})
    L = SpectralIndex.build(lib)
    assert len(L) == 5 and L.header['n_input_spectra'] == 6
    h = search(L, lib, method='identity')          # no precursor: no identity hits
    assert len(h) == 0
    h = search(L, lib, method='open', top_k=1)
    np.testing.assert_array_equal(h.query_idx, h.library_idx)


def test_read_msp():
    text = ('NAME: A\nPRECURSORMZ: 150.5\nNum Peaks: 2\n50.1\t10\n60.2 20\n\n'
            'NAME: B\nPrecursorMZ: 200\nNum Peaks: 1\n70 5\n')
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'x.msp')
        with open(path, 'w') as f:
            f.write(text)
        s = read_msp(path)
    assert [x['name'] for x in s] == ['A', 'B']
    assert s[0]['precursor_mz'] == 150.5 and s[1]['precursor_mz'] == 200.0
    np.testing.assert_allclose(s[0]['peaks'], [[50.1, 10], [60.2, 20]], rtol=1e-6)


def load_tests(loader, tests, pattern):
    '''Run the test_ functions under unittest discovery as well.'''
    return unittest.TestSuite(unittest.FunctionTestCase(v) for k, v in sorted(globals().items())
                              if k.startswith('test_') and callable(v))


if __name__ == '__main__':
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith('test_') and callable(v)]
    for name, fn in tests:
        fn()
        print('ok  ', name)
    print(f'{len(tests)} tests passed')
