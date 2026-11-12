"""
ms2batch: fast batch MS/MS spectral search by entropy similarity.

    from ms2batch import SpectralIndex, search, read_spectra

    library = SpectralIndex.build(read_spectra('library.msp'))
    library.save('library.npz')

    library = SpectralIndex.load('library.npz')
    hits = search(library, read_spectra('queries.mgf'), method='identity', top_k=5)
"""

__version__ = '0.1.0'

from .index import SpectralIndex
from .engine import search, SearchResult
from .io import read_spectra, read_msp, read_mgf
