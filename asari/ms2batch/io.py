"""
Minimal readers for MSP and MGF spectral files.

Each spectrum is returned as a dict with 'peaks' (float32 array of shape
(n, 2)), 'precursor_mz' (float or None), and the other fields as strings,
keys lowercased.
"""

import json

import numpy as np

PRECURSOR_KEYS = ('precursor_mz', 'precursormz', 'pepmass', 'precursor mz', 'precursor_m/z')


def read_spectra(path):
    """Read .msp, .mgf or .json (list of dicts with 'peaks') by file extension."""
    p = str(path).lower()
    if p.endswith('.msp'):
        return read_msp(path)
    if p.endswith('.mgf'):
        return read_mgf(path)
    if p.endswith('.json'):
        with open(path) as f:
            spectra = json.load(f)
        for s in spectra:
            _set_precursor({k.lower(): v for k, v in s.items()}, s)
        return spectra
    raise ValueError(f'unsupported file type: {path}')


def read_msp(path):
    spectra = []
    with open(path, encoding='utf-8', errors='replace') as f:
        meta, peak_lines = {}, []
        for line in f:
            line = line.strip()
            if not line:
                if peak_lines or meta:
                    spectra.append(_make(meta, peak_lines))
                meta, peak_lines = {}, []
            elif line[0].isdigit() or line[0] == '.':
                peak_lines.append(line)
            elif ':' in line:
                k, v = line.split(':', 1)
                meta[k.strip().lower()] = v.strip()
        if peak_lines or meta:
            spectra.append(_make(meta, peak_lines))
    return spectra


def read_mgf(path):
    spectra = []
    with open(path, encoding='utf-8', errors='replace') as f:
        meta, peak_lines, inside = {}, [], False
        for line in f:
            line = line.strip()
            if line == 'BEGIN IONS':
                meta, peak_lines, inside = {}, [], True
            elif line == 'END IONS':
                spectra.append(_make(meta, peak_lines))
                inside = False
            elif inside and line:
                if line[0].isdigit():
                    peak_lines.append(line)
                elif '=' in line:
                    k, v = line.split('=', 1)
                    meta[k.strip().lower()] = v.strip()
    return spectra


def _make(meta, peak_lines):
    values = []
    for line in peak_lines:
        parts = line.replace(',', ' ').split()
        if len(parts) >= 2:
            values.append((parts[0], parts[1]))
    peaks = np.array(values, dtype=np.float32).reshape(-1, 2)
    spectrum = dict(meta)
    spectrum['peaks'] = peaks
    _set_precursor(meta, spectrum)
    return spectrum


def _set_precursor(meta, spectrum):
    for k in PRECURSOR_KEYS:
        if meta.get(k) not in (None, ''):
            try:
                spectrum['precursor_mz'] = float(str(meta[k]).split()[0])
                return
            except ValueError:
                pass
    spectrum.setdefault('precursor_mz', None)
