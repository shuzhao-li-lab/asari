"""
SpectralIndex: a cleaned, packed and indexed set of MS/MS spectra.

The same structure serves as a search library or as a batch of queries, so
dataset-vs-dataset search only indexes each side once.

On-disk format (.npz, uncompressed, no pickle)
----------------------------------------------
A standard NumPy .npz archive, readable with np.load(path, allow_pickle=False).
Every member is a plain numeric array; nothing is executed on load.

    header             uint8   UTF-8 JSON: format name/version, parameters, counts
    precursor_mz       float64 (n,)    sorted ascending, NaN last
    source_idx         int64   (n,)    position of each spectrum in the input
    offsets            int64   (n+1,)  spectrum i owns peaks offsets[i]:offsets[i+1]
    mz                 float32 (p,)    peaks, spectrum-major, m/z-sorted per spectrum
    intensity          float32 (p,)    search-ready intensities (weighted, sum 1)
    frag_mz            float32 (p,)    same peaks sorted by m/z across the library
    frag_intensity     float32 (p,)
    frag_spec          uint32  (p,)    spectrum of each fragment-sorted peak
    metadata_json      uint8   (b,)    optional, concatenated UTF-8 JSON objects
    metadata_offsets   int64   (n+1,)  optional

Members are stored uncompressed so they can be memory-mapped (mmap=True).
"""

import datetime
import json
import struct
import zipfile

import numpy as np

from . import __version__
from .preprocess import (clean_spectra, entropy_weight, expand_ranges, offsets_from_ids,
                         pack_peaks, spectrum_ids)

FORMAT_NAME = 'ms2batch.spectral_index'
FORMAT_VERSION = 1

DEFAULT_PARAMS = {
    'precursor_removal_da': 1.6,
    'noise_threshold': 0.01,
    'min_peak_distance': 0.05,
    'max_peak_num': None,
    'intensity_weight': 'entropy',   # or None for unweighted entropy similarity
}

ARRAY_DTYPES = {
    'precursor_mz': np.float64,
    'source_idx': np.int64,
    'offsets': np.int64,
    'mz': np.float32,
    'intensity': np.float32,
    'frag_mz': np.float32,
    'frag_intensity': np.float32,
    'frag_spec': np.uint32,
    'metadata_json': np.uint8,
    'metadata_offsets': np.int64,
}
REQUIRED = ('precursor_mz', 'source_idx', 'offsets', 'mz', 'intensity')
FRAGMENT_ARRAYS = ('frag_mz', 'frag_intensity', 'frag_spec')


class SpectralIndex:
    """
    Build with SpectralIndex.build(spectra), save with .save(path),
    load with SpectralIndex.load(path).
    """

    def __init__(self, arrays, header):
        self.header = header
        self.params = header['params']
        for name in ARRAY_DTYPES:
            setattr(self, name, arrays.get(name))
        self._cache = {}
        self._validate()

    # ------------------------------------------------------------------ build

    @classmethod
    def build(cls, spectra, fragment_index=True, keep_metadata=True, clean=True, **params):
        """
        Build an index from a list of spectra, each a dict with 'peaks'
        (array-like of shape (n, 2), [m/z, intensity]) and optionally
        'precursor_mz'. Other keys are kept as metadata if keep_metadata.

        params override DEFAULT_PARAMS. fragment_index=False skips the
        fragment-sorted arrays, which are only needed on the library side.
        clean=False takes peaks as already cleaned (e.g. by ms_entropy): they
        are only sorted and normalized, and must be >= min_peak_distance apart.
        """
        precursor = [s.get('precursor_mz') for s in spectra]
        precursor = np.array([np.nan if p is None else float(p) for p in precursor], dtype=np.float64)
        metadata = None
        if keep_metadata:
            metadata = [{k: v for k, v in s.items() if k != 'peaks'} for s in spectra]
        return cls.from_arrays([s['peaks'] for s in spectra], precursor, metadata,
                               fragment_index=fragment_index, clean=clean, **params)

    @classmethod
    def from_arrays(cls, peaks_list, precursor_mz=None, metadata=None, fragment_index=True,
                    clean=True, **params):
        """Build an index from a list of peak arrays, with optional precursor m/z and metadata."""
        unknown = set(params) - set(DEFAULT_PARAMS)
        if unknown:
            raise TypeError(f'unknown parameters: {sorted(unknown)}')
        p = {**DEFAULT_PARAMS, **params}
        n = len(peaks_list)
        precursor_mz = (np.full(n, np.nan) if precursor_mz is None
                        else np.asarray(precursor_mz, dtype=np.float64))

        if clean:
            mz, intensity, offsets = clean_spectra(
                peaks_list, precursor_mz,
                precursor_removal_da=p['precursor_removal_da'],
                noise_threshold=p['noise_threshold'],
                min_peak_distance=p['min_peak_distance'],
                max_peak_num=p['max_peak_num'])
        else:
            mz, intensity, offsets = _precleaned(peaks_list, p['min_peak_distance'])
        if p['intensity_weight'] == 'entropy':
            intensity = entropy_weight(intensity, offsets)
        elif p['intensity_weight'] is not None:
            raise ValueError("intensity_weight must be 'entropy' or None")

        # Drop empty spectra; order the rest by precursor m/z (NaN last).
        counts = np.diff(offsets)
        source_idx = np.flatnonzero(counts > 0)
        source_idx = source_idx[np.argsort(precursor_mz[source_idx], kind='stable')]
        new_offsets = np.zeros(len(source_idx) + 1, dtype=np.int64)
        np.cumsum(counts[source_idx], out=new_offsets[1:])
        peak_order = expand_ranges(offsets[source_idx], counts[source_idx])

        arrays = {
            'precursor_mz': precursor_mz[source_idx],
            'source_idx': source_idx.astype(np.int64),
            'offsets': new_offsets,
            'mz': mz[peak_order],
            'intensity': intensity[peak_order],
        }
        if fragment_index:
            arrays.update(_fragment_arrays(arrays['mz'], arrays['intensity'], new_offsets))
        if metadata is not None:
            arrays.update(_pack_metadata([metadata[i] for i in source_idx]))

        header = {
            'format': FORMAT_NAME,
            'format_version': FORMAT_VERSION,
            'software': f'ms2batch {__version__}',
            'created': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
            'n_input_spectra': n,
            'n_spectra': len(source_idx),
            'n_peaks': int(new_offsets[-1]),
            'params': p,
        }
        return cls(arrays, header)

    def with_fragment_index(self):
        """Add fragment-sorted arrays if missing (e.g. to use a query set as library)."""
        if self.frag_mz is None:
            for k, v in _fragment_arrays(self.mz, self.intensity, self.offsets).items():
                setattr(self, k, v)
        return self

    # --------------------------------------------------------------- file i/o

    def save(self, path):
        """Save index to an uncompressed .npz file."""
        arrays = {name: getattr(self, name) for name in ARRAY_DTYPES if getattr(self, name) is not None}
        header = np.frombuffer(json.dumps(self.header).encode('utf-8'), dtype=np.uint8)
        with open(path, 'wb') as f:
            np.savez(f, header=header, **arrays)

    @classmethod
    def load(cls, path, mmap=False):
        """
        Load an index. With mmap=True, arrays are memory-mapped read-only
        instead of read into memory. Never unpickles.
        """
        if mmap:
            arrays = _npz_memmap(path)
        else:
            with np.load(path, allow_pickle=False) as z:
                arrays = {k: z[k] for k in z.files}
        header = json.loads(bytes(arrays.pop('header')).decode('utf-8'))
        if header.get('format') != FORMAT_NAME:
            raise ValueError(f'{path} is not a {FORMAT_NAME} file')
        if header.get('format_version', 0) > FORMAT_VERSION:
            raise ValueError(f'{path} has format version {header["format_version"]}; '
                             f'this ms2batch reads up to {FORMAT_VERSION}')
        return cls({k: v for k, v in arrays.items() if k in ARRAY_DTYPES}, header)

    # ---------------------------------------------------------------- access

    def __len__(self):
        return len(self.precursor_mz)

    @property
    def n_peaks(self):
        """Total number of peaks."""
        return int(self.offsets[-1])

    def peaks(self, i):
        """Search-ready peaks of the i-th stored spectrum, shape (n, 2)."""
        a, b = self.offsets[i], self.offsets[i + 1]
        return np.column_stack([self.mz[a:b], self.intensity[a:b]])

    def metadata(self, source_idx):
        """Metadata dict for a spectrum, by its position in the input list."""
        if self.metadata_json is None:
            return None
        pos = self._source_to_pos()[source_idx]
        if pos < 0:
            return None
        a, b = self.metadata_offsets[pos], self.metadata_offsets[pos + 1]
        return json.loads(bytes(self.metadata_json[a:b]).decode('utf-8'))

    def _source_to_pos(self):
        if 'source_to_pos' not in self._cache:
            m = np.full(self.header['n_input_spectra'], -1, dtype=np.int64)
            m[self.source_idx] = np.arange(len(self))
            self._cache['source_to_pos'] = m
        return self._cache['source_to_pos']

    def cached(self, name, fn):
        """Per-index cache for derived arrays used in search."""
        if name not in self._cache:
            self._cache[name] = fn()
        return self._cache[name]

    def __repr__(self):
        p = self.params
        return (f'SpectralIndex({len(self)} spectra, {self.n_peaks} peaks, '
                f'weight={p["intensity_weight"]}, min_peak_distance={p["min_peak_distance"]}, '
                f'fragment_index={self.frag_mz is not None})')

    def _validate(self):
        for name in REQUIRED:
            if getattr(self, name) is None:
                raise ValueError(f'index is missing array {name!r}')
        for name, dtype in ARRAY_DTYPES.items():
            a = getattr(self, name)
            if a is not None and (a.dtype != dtype or a.ndim != 1):
                raise ValueError(f'array {name!r} must be 1-D {np.dtype(dtype)}, got {a.dtype} {a.shape}')
        n, p = len(self.precursor_mz), len(self.mz)
        if len(self.source_idx) != n or len(self.offsets) != n + 1 or self.offsets[-1] != p \
                or len(self.intensity) != p:
            raise ValueError('inconsistent array lengths in index')
        if self.frag_mz is not None and not (len(self.frag_mz) == len(self.frag_intensity) == len(self.frag_spec) == p):
            raise ValueError('inconsistent fragment index lengths')


def _precleaned(peaks_list, min_peak_distance):
    mz, intensity, spec = pack_peaks(peaks_list)
    order = np.lexsort((mz, spec))
    mz, intensity, spec = mz[order], intensity[order], spec[order]
    if np.any((spec[1:] == spec[:-1]) & (np.diff(mz).astype(np.float64) < min_peak_distance)):
        raise ValueError('pre-cleaned peaks must be at least min_peak_distance apart')
    if np.any(intensity <= 0):
        raise ValueError('pre-cleaned peaks must have positive intensity')
    n = len(peaks_list)
    totals = np.bincount(spec, weights=intensity, minlength=n)
    return mz, (intensity / totals[spec]).astype(np.float32), offsets_from_ids(spec, n)


def _fragment_arrays(mz, intensity, offsets):
    order = np.argsort(mz, kind='stable')
    return {
        'frag_mz': mz[order],
        'frag_intensity': intensity[order],
        'frag_spec': spectrum_ids(offsets)[order].astype(np.uint32),
    }


def _pack_metadata(records):
    blobs = [json.dumps(r, default=str).encode('utf-8') for r in records]
    offsets = np.zeros(len(blobs) + 1, dtype=np.int64)
    np.cumsum([len(b) for b in blobs], out=offsets[1:])
    return {'metadata_json': np.frombuffer(b''.join(blobs), dtype=np.uint8),
            'metadata_offsets': offsets}


def _npz_memmap(path):
    """Memory-map every member of an uncompressed .npz archive, read-only."""
    out = {}
    with zipfile.ZipFile(path) as zf, open(path, 'rb') as f:
        for info in zf.infolist():
            if not info.filename.endswith('.npy'):
                continue
            if info.compress_type != zipfile.ZIP_STORED:
                raise ValueError(f'{info.filename} is compressed; load with mmap=False')
            f.seek(info.header_offset)
            local = f.read(30)
            if local[:4] != b'PK\x03\x04':
                raise ValueError(f'bad zip entry {info.filename}')
            name_len, extra_len = struct.unpack('<HH', local[26:30])
            f.seek(info.header_offset + 30 + name_len + extra_len)
            version = np.lib.format.read_magic(f)
            if version == (1, 0):
                shape, fortran, dtype = np.lib.format.read_array_header_1_0(f)
            else:
                shape, fortran, dtype = np.lib.format.read_array_header_2_0(f)
            if dtype.hasobject:
                raise ValueError(f'{info.filename} holds Python objects; refusing to load')
            name = info.filename[:-4]
            if np.prod(shape) == 0:
                out[name] = np.zeros(shape, dtype=dtype)
            else:
                out[name] = np.memmap(path, dtype=dtype, mode='r', offset=f.tell(),
                                      shape=shape, order='F' if fortran else 'C')
    return out
