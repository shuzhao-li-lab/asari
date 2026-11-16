"""
Command line interface.

    python -m asari.ms2batch index library.msp library.npz
    python -m asari.ms2batch index flash_library.pkl library.npz     (convert legacy pickle, trusted files only)
    python -m asari.ms2batch search library.npz queries.mgf -o hits.tsv --method identity --top-k 5
"""

import argparse
import sys
import time

from . import SpectralIndex, read_spectra, search


def main(argv=None):
    """Command line interface of ms2batch."""
    parser = argparse.ArgumentParser(prog='ms2batch', description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('index', help='build an index from .msp/.mgf/.json, or convert a FlashEntropySearch .pkl')
    p.add_argument('spectra')
    p.add_argument('output', help='index file (.npz)')
    p.add_argument('--min-peak-distance', type=float, default=0.05)
    p.add_argument('--noise-threshold', type=float, default=0.01)
    p.add_argument('--precursor-removal-da', type=float, default=1.6,
                   help='remove peaks >= precursor m/z minus this; negative to disable')
    p.add_argument('--unweighted', action='store_true', help='no entropy intensity weighting')
    p.add_argument('--no-metadata', action='store_true')

    p = sub.add_parser('search', help='search queries against an index')
    p.add_argument('library', help='index file (.npz)')
    p.add_argument('queries', help='.msp/.mgf/.json, or an index file (.npz)')
    p.add_argument('-o', '--output', default='-', help='TSV output (default stdout)')
    p.add_argument('--method', choices=['identity', 'open'], default='identity')
    p.add_argument('--ms1-tolerance', type=float, default=0.01)
    p.add_argument('--ms2-tolerance', type=float, default=0.02)
    p.add_argument('--top-k', type=int, default=5)
    p.add_argument('--min-score', type=float, default=0.5)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--mmap', action='store_true', help='memory-map the library index')

    args = parser.parse_args(argv)
    t0 = time.perf_counter()
    if args.command == 'index':
        if args.spectra.lower().endswith(('.pkl', '.pickle')):
            index = _index_from_flash_pickle(args.spectra, keep_metadata=not args.no_metadata,
                                             intensity_weight=None if args.unweighted else 'entropy')
            index.save(args.output)
            print(f'{index} -> {args.output} ({time.perf_counter() - t0:.1f} s)', file=sys.stderr)
            return
        spectra = read_spectra(args.spectra)
        index = SpectralIndex.build(
            spectra, keep_metadata=not args.no_metadata,
            min_peak_distance=args.min_peak_distance, noise_threshold=args.noise_threshold,
            precursor_removal_da=None if args.precursor_removal_da < 0 else args.precursor_removal_da,
            intensity_weight=None if args.unweighted else 'entropy')
        index.save(args.output)
        print(f'{index} -> {args.output} ({time.perf_counter() - t0:.1f} s)', file=sys.stderr)

    elif args.command == 'search':
        library = SpectralIndex.load(args.library, mmap=args.mmap)
        if args.queries.endswith('.npz'):
            queries = SpectralIndex.load(args.queries)
        else:
            queries = read_spectra(args.queries)
        hits = search(library, queries, method=args.method, ms1_tolerance=args.ms1_tolerance,
                      ms2_tolerance=args.ms2_tolerance, top_k=args.top_k, min_score=args.min_score,
                      n_threads=args.threads)
        out = sys.stdout if args.output == '-' else open(args.output, 'w')
        name = _name_lookup(library)
        out.write('query_idx\tlibrary_idx\tscore\tn_matched\tlibrary_name\n')
        for q, l, s, n in zip(hits.query_idx, hits.library_idx, hits.score, hits.n_matched):
            out.write(f'{q}\t{l}\t{s:.4f}\t{n}\t{name(l)}\n')
        if out is not sys.stdout:
            out.close()
        print(f'{len(hits)} hits ({time.perf_counter() - t0:.1f} s)', file=sys.stderr)


def _index_from_flash_pickle(path, keep_metadata=True, intensity_weight='entropy'):
    """
    Convert a FlashEntropySearch pickle (ms_entropy) to SpectralIndex, for migration only.
    Unpickling can execute code: use only on trusted files. Requires ms_entropy.
    Spectra in the pickle are already cleaned by ms_entropy, thus clean=False.
    """
    import pickle
    with open(path, 'rb') as f:
        entropy_search = pickle.load(f)
    spectra = [entropy_search[i] for i in range(len(entropy_search.precursor_mz_array))]
    return SpectralIndex.build(spectra, keep_metadata=keep_metadata, clean=False,
                               intensity_weight=intensity_weight)


def _name_lookup(index):
    def name(i):
        m = index.metadata(i) or {}
        return m.get('name') or m.get('id') or m.get('title') or ''
    return name


if __name__ == '__main__':
    main()
