"""
Functions for subcommand `annotate`

Annotation links compound identity to features measured in an experiment.

- For LC-MS, annotation follows [L1.lib -> features].
- For LC-MS/MS, annotation follows [MS2 -> L2.lib -> features].
Where lib is a compound library with m/z and retention time (and CCS etc.), to be matched to experimentally measured parameters.

MS2 involves clustering of experimental spectra and database search. The result can be in compatible `lib` format.

Mulitple databases and libs can be used. Additional annotation methods (e.g. CSM) are possible. We use a `combine` step to combine results from multiple libs.

As multiple features can derive from the same compound, pre-annotation assigns such relationships (isotopes, adducts, fragments and conjugates) in MS and LC-MS.
Annotation should prioritize primary ions or neutral mass.
Mass match alone (low level of confidence) should center on neutral mass.

GC-MS data are dominated by fragments, requiring a different workflow.
"""

import os
import json
import time
import pickle
import numpy as np
import importlib.resources as pkg_resources
from pathvalidate import sanitize_filename

# from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase
from . import db
from .default_parameters import (PROTON, adduct_search_patterns_pos, adduct_search_patterns_neg)
from .utils import NpEncoder
from .tools.file_io import read_features_from_asari_table
from .tools.match_features import list_match_lcms_features
from .mass_functions import all_mass_paired_mapping
from .tools.plot import mirror_plot
from .tools.entropy_search import clean_list_ms2spectra
from .tools.msp_parser import parse_msp_to_listdict
from .tools.cosine import cosine_similarity
from .ms2batch import SpectralIndex, read_spectra
from .ms2batch import search as ms2batch_search
from .gc_functions import *
from .deconvolute import (ingest_gc_featuretable, run_gc_deconvolution)

def annotate_project(infile, parameters):
    """Run subcommand `annotate` on infile, by parameters['workflow'] (LC, GC or LCMSMS)."""
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    subdir = '_'.join([parameters['project_name'], 'annotation', ''.join(time_stamp)])
    outdir = os.path.join(parameters['outdir'], subdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"Annotation directory: {outdir}.\n\n")

    if parameters['workflow'] == "LCMSMS" and parameters['search_algorithm'] in ('mz', 'rtmz'):
        # shortcut: MS1 search of the LC-MS/MS feature table is the same as LC annotation
        print(f"Shortcut: --workflow LCMSMS with --search_algorithm {parameters['search_algorithm']} "
              "is an MS1 search, run as --workflow LC on the feature table.\n")
        parameters['workflow'] = 'LC'
        if os.path.isdir(infile):
            infile = os.path.join(infile, 'Feature_table.tsv')
            print(f"Using feature table {infile}.\n")

    if parameters['workflow'] == "LC":
        # to add deconvolute and algorithms in version 2
        if parameters['search_algorithm'] not in ('mz', 'rtmz'):
            print(f"search_algorithm {parameters['search_algorithm']} is not for LC; using 'mz'.")
            parameters['search_algorithm'] = 'mz'
        ANN = LCMS_Annotation(dict(parameters, outdir=outdir))
        list_empCpds = []
        if parameters.get('preannotation'):
            list_empCpds = load_khipu_empCpds(parameters['preannotation'])
        ANN.annotate_user_featuretable(infile, list_empCpds=list_empCpds)

    elif parameters['workflow'] == "GC":
        ANN = GCMS_Annotation(dict(parameters, outdir=outdir))
        if infile.endswith('.tsv'):
            ANN.table_workflow(infile)
        else:
            list_query_spectra = load_gcms_dbfile(infile)   # auto determine msp or json
            if list_query_spectra:
                ANN.pseudospectra_workflow(list_query_spectra)
            else:
                print(f"No pseudospectra found in {infile}.")

    elif parameters['workflow'] in ("LCMSMS", "MSMS"):
        # to handle diff anno methods, including pre-annotation in version 2
        ANN = LCMSMS_Annotation(dict(parameters, outdir=outdir))
        ANN.annotate_spectra(infile)

    # recording project parameters as meta data
    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)


#
# -----------------------------------------------------------------------------
#

class LCMSMS_Annotation:
    """
    Annotation of MS/MS spectra (LCMSMS or MSMS workflow) against a spectral library (--db).

    Library formats, by file extension:
    - .npz (default, recommended): ms2batch SpectralIndex, searched in batch mode (vectorized, NumPy only).
      Build once by `python -m asari.ms2batch index library.msp library.npz`.
    - .msp, .mgf, .json: plain spectra, indexed in memory as SpectralIndex, then searched in batch mode.
    - .pkl, .pickle: FlashEntropySearch pickle, searched per spectrum (legacy).
      Pickle is not secure for sharing databases; convert by `python -m asari.ms2batch index library.pkl library.npz`.
    Both search modes give the same entropy and cosine scores.

    Input spectra are from asari `process` (LCMSMS or MSMS workflow), in JSON (ms2_spectra.json)
    or MSP (ms2_spectra.msp, ms2_spectra_<sample>.msp), or a process directory with ms2_spectra.json/.msp.

    Steps: entropy identity search (low prefilter), then cosine similarity on the best entropy match;
    matches pass by parameters['score_algorithm'] and score_cutoff_cosine or score_cutoff_entropy;
    redundant library compounds (by database_key, e.g. INCHIKEY) are collapsed to the best scoring match.

    Outputs in outdir: <project>_MS2_annotation.tsv and <project>_MS2_annotation.json
    (with experimental 'peaks' and library 'peaks_in_lib' for mirror plots in viz Dashboard).
    """
    entropy_prefilter = 0.2         # intentionally low on entropy identity search, to filter later
    precursor_mz_offset = 1.6       # removal of precursor ions in cleaning query spectra, in Da
    database_keys = ('INCHIKEY', 'InChIKey', 'inchikey', 'NAME', 'name')    # to collapse redundant entries

    def __init__(self, parameters):
        self.parameters = parameters
        self.outdir = parameters['outdir']
        self.project_name = parameters['project_name']
        self.database_file = parameters.get('db')
        self.ms1_tolerance_in_da = parameters['lcmsms_ms1_tolerance_in_da']
        self.ms2_tolerance_in_da = parameters['lcmsms_ms2_tolerance_in_da']
        self.score_algorithm = parameters['score_algorithm']
        if self.score_algorithm == 'entropy':
            self.score_cutoff = parameters['score_cutoff_entropy']
        else:
            self.score_algorithm, self.score_cutoff = 'cosine', parameters['score_cutoff_cosine']
        self.entropy_search = None      # FlashEntropySearch, legacy pickle
        self.library = None             # ms2batch SpectralIndex, batch mode

    def annotate_spectra(self, infile):
        """
        Workflow: load database and query spectra, search, collapse and export.
        Returns list of annotation results (JSON style), or None if no database.
        """
        if not self.database_file:
            print("LCMSMS/MSMS annotation requires an MS/MS library via --db: "
                  ".npz index (preferred, see `python -m asari.ms2batch index -h`), or .msp, .mgf, .json.")
            return None
        self.load_database(self.database_file)
        list_query_spectra = self.load_query_spectra(infile)
        results = self.search(list_query_spectra)
        concise_results = self.collapse_by_database_key(results)
        json_results = self.format_results(concise_results, list_query_spectra)
        self.export(json_results)
        return json_results

    def load_database(self, database_file):
        """
        Load MS/MS library by file extension: .npz (SpectralIndex, default), .msp/.mgf/.json 
        (indexed in memory) or .pkl/.pickle (FlashEntropySearch, legacy).
        Sets self.library (batch mode) or self.entropy_search.
        """
        ext = os.path.splitext(database_file)[1].lower()
        if ext == '.npz':
            self.library = SpectralIndex.load(database_file)
        elif ext in ('.msp', '.mgf', '.json'):
            self.library = SpectralIndex.build(read_spectra(database_file))
            print(f"Indexed {database_file}. To skip indexing in future runs, save as .npz by "
                  f"`python -m asari.ms2batch index {database_file} <library.npz>`.")
        elif ext in ('.pkl', '.pickle'):
            print("Warning: pickle files are not secure to share; please convert the database to .npz by "
                  f"`python -m asari.ms2batch index {database_file} <library.npz>`.")
            with open(database_file, 'rb') as f:
                self.entropy_search = pickle.load(f)
            print(f"Loaded database of {len(self.entropy_search.precursor_mz_array)} spectra from {database_file}.")
            return self.entropy_search
        else:
            raise ValueError(f"Unsupported MS/MS database format: {database_file}. Use .npz, .msp, .mgf, .json or .pkl.")
        print(f"Loaded database of {len(self.library)} spectra from {database_file}.")
        return self.library

    def load_query_spectra(self, infile):
        """
        Load experimental MS/MS spectra from JSON or MSP, or from a process directory 
        (ms2_spectra.json, else ms2_spectra.msp).
        Returns list of spectra as dicts with 'id', 'precursor_mz', 'rtime', 'peaks'.
        """
        if os.path.isdir(infile):
            for f in ('ms2_spectra.json', 'ms2_spectra.msp'):
                if os.path.exists(os.path.join(infile, f)):
                    infile = os.path.join(infile, f)
                    break
            else:
                raise FileNotFoundError(f"No ms2_spectra.json or ms2_spectra.msp in {infile}.")

        if infile.lower().endswith('.msp'):
            list_query_spectra = []
            for ii, e in enumerate(parse_msp_to_listdict(infile)):
                precursor_mz = e.get('PRECURSORMZ') or e.get('PrecursorMZ')
                if precursor_mz in (None, ''):
                    continue
                rtime = e.get('RETENTIONTIME')
                list_query_spectra.append({
                    'id': e.get('ID') or e.get('NAME') or f"sp{ii}",
                    'precursor_mz': float(precursor_mz),
                    'rtime': float(rtime) if rtime not in (None, '') else None,
                    'peaks': [list(x) for x in e['peaks']],
                })
        else:
            with open(infile, encoding='utf-8') as f:
                list_query_spectra = json.load(f)
        print(f"Loaded {len(list_query_spectra)} experimental spectra from {infile}.")
        return list_query_spectra

    def search(self, list_query_spectra):
        """
        Entropy identity search of query spectra; cosine similarity on the best entropy match.
        Batch mode if library is a SpectralIndex (search_batch), else per spectrum by FlashEntropySearch.
        Returns list of (library entry, cosine_score, entropy_score, number_matched_peaks, spid, precursor_mz, rtime)
        that pass score cutoff by self.score_algorithm.
        """
        if self.library is not None:
            return self.search_batch(list_query_spectra)
        return self.search_flash(list_query_spectra)

    def search_batch(self, list_query_spectra):
        """
        Batch identity search (ms2batch) of all query spectra against self.library.
        Queries are cleaned with the library parameters (incl. precursor removal);
        cosine similarity uses cleaned, unweighted peaks, as in search_flash.
        """
        L = self.library
        queries = SpectralIndex.build(
            [{'peaks': sp['peaks'], 'precursor_mz': sp['precursor_mz']} for sp in list_query_spectra],
            fragment_index=False, keep_metadata=False, **L.params)
        print(f"Batch searching {len(queries)} cleaned experimental spectra, "
              f"{self.score_algorithm} score cutoff {self.score_cutoff}.\n")
        prefilter = min(self.entropy_prefilter, self.score_cutoff)
        hits = ms2batch_search(L, queries, method='identity',
                               ms1_tolerance=self.ms1_tolerance_in_da,
                               ms2_tolerance=self.ms2_tolerance_in_da,
                               top_k=1, min_score=prefilter,
                               n_threads=self.parameters.get('multicores', 1))
        lib_pos, query_pos = L._source_to_pos(), queries._source_to_pos()
        results = []
        for q, lib_idx, entropy_score, n_matched in zip(hits.query_idx, hits.library_idx, hits.score, hits.n_matched):
            if entropy_score <= prefilter:
                continue
            lib_peaks = L.peaks(lib_pos[lib_idx], unweighted=True)
            peaks = queries.peaks(query_pos[q], unweighted=True)
            cosine_score, _ = cosine_similarity(peaks, lib_peaks, mz_tolerance=self.ms2_tolerance_in_da)
            entropy_score = float(entropy_score)
            score = entropy_score if self.score_algorithm == 'entropy' else cosine_score
            if score > self.score_cutoff:
                sp = list_query_spectra[q]
                entry = dict(L.metadata(lib_idx) or {'id': int(lib_idx)})
                entry['peaks'] = lib_peaks
                results.append((entry, float(cosine_score), entropy_score, int(n_matched),
                                sp['id'], sp['precursor_mz'], sp.get('rtime')))
        print(f"{len(results)} spectra matched to database.")
        return results

    def search_flash(self, list_query_spectra):
        """Per spectrum identity search by FlashEntropySearch (legacy pickle database)."""
        cleaned_spectra = clean_list_ms2spectra(list_query_spectra, self.entropy_search,
                                                {'precursor_mz_offset': self.precursor_mz_offset})
        print(f"Searching {len(cleaned_spectra)} cleaned experimental spectra, "
              f"{self.score_algorithm} score cutoff {self.score_cutoff}.\n")
        prefilter = min(self.entropy_prefilter, self.score_cutoff)
        results = []
        for (precursor_mz, peaks, rtime, spid) in cleaned_spectra:
            similarity, matched_num = self.entropy_search.identity_search(
                precursor_mz=precursor_mz,
                peaks=peaks,
                ms1_tolerance_in_da=self.ms1_tolerance_in_da,
                ms2_tolerance_in_da=self.ms2_tolerance_in_da,
                output_matched_peak_number=True)
            if len(similarity) == 0 or np.max(similarity) <= prefilter:
                continue
            best_idx = int(np.argmax(similarity))
            entry = dict(self.entropy_search[best_idx])
            entropy_score = float(similarity[best_idx])
            cosine_score, _ = cosine_similarity(peaks, entry['peaks'], mz_tolerance=self.ms2_tolerance_in_da)
            score = entropy_score if self.score_algorithm == 'entropy' else cosine_score
            if score > self.score_cutoff:
                results.append((entry, float(cosine_score), entropy_score, int(matched_num[best_idx]),
                                spid, precursor_mz, rtime))
        print(f"{len(results)} spectra matched to database.")
        return results

    def get_database_key(self, entry):
        """Compound key of a library entry, first available of self.database_keys."""
        for k in self.database_keys:
            if entry.get(k):
                return entry[k]
        return None

    def collapse_by_database_key(self, results):
        """
        Collapse redundant library compounds (same database key), keeping the best match 
        by self.score_algorithm. Consider to keep top N later.
        Returns {database key: result}.
        """
        score_idx = 2 if self.score_algorithm == 'entropy' else 1
        concise_results = {}
        for res in results:
            k = self.get_database_key(res[0]) or res[4]
            if k not in concise_results or res[score_idx] > concise_results[k][score_idx]:
                concise_results[k] = res
        return concise_results

    def format_results(self, concise_results, list_query_spectra):
        """
        Format results as JSON style list, keeping experimental 'peaks' and library 'peaks_in_lib'
        (JSON only, not in tsv) for mirror plots.
        """
        dict_query_peaks = {sp['id']: sp['peaks'] for sp in list_query_spectra}
        json_results = []
        for k, res in concise_results.items():
            entry = dict(res[0])
            peaks_in_lib = [[round(float(mz), 4), float(intensity)] for mz, intensity in entry.pop('peaks')]
            json_results.append({
                'spid': res[4],
                'mz': round(float(res[5]), 4),          # expt_precursor_mz
                'rtime': res[6],                        # expt_rtime
                'entropy_score': round(res[2], 3),
                'cosine_score': round(res[1], 3),
                'number_matched_peaks': res[3],
                'database_key': k,
                'database_entry': entry,
                'peaks': dict_query_peaks.get(res[4], []),     # experimental MS/MS spectrum
                'peaks_in_lib': peaks_in_lib,                   # library spectrum
            })
        return json_results

    def export(self, json_results):
        """Export <project>_MS2_annotation.tsv and .json to self.outdir."""
        tsv_output_file = os.path.join(self.outdir, self.project_name + '_MS2_annotation.tsv')
        keys1 = ['spid', 'mz', 'rtime', 'entropy_score', 'cosine_score', 'number_matched_peaks']
        keys2 = []
        for row in json_results:        # union of database fields, in order of appearance
            keys2 += [k for k in row['database_entry'] if k not in keys2]
        with open(tsv_output_file, 'w', encoding='utf-8') as O:
            O.write('\t'.join(keys1 + keys2) + '\n')
            for row in json_results:
                O.write('\t'.join([str(row.get(k, '')) for k in keys1]
                                  + [str(row['database_entry'].get(k, '')) for k in keys2]) + '\n')

        json_export = {
            'type': 'MS/MS',
            'database': self.database_file,
            'algorithm': 'entropy,cosine',
            'score_algorithm': self.score_algorithm,
            'score_cutoff': self.score_cutoff,
            'ms1_tolerance_in_da': self.ms1_tolerance_in_da,
            'ms2_tolerance_in_da': self.ms2_tolerance_in_da,
            'time_stamp': ':'.join([str(x) for x in time.localtime()[1:6]]),
            'list_of_features': json_results,
        }
        json_output_file = os.path.join(self.outdir, self.project_name + '_MS2_annotation.json')
        with open(json_output_file, 'w', encoding='utf-8') as O:
            json.dump(json_export, O, indent=2, ensure_ascii=False, cls=NpEncoder)
        print(f"MS/MS annotation results ({len(json_results)}) in JSON and tsv were written to {self.outdir}.\n\n")


#
# -----------------------------------------------------------------------------
#

class LCMS_Annotation:
    """
    Annotation of LC-MS feature tables by matching to known compound attributes. 
    Most common are by m/z and retention time. 
    When MS/MS spectra are used to in the MS/MS workflow, 
    database matching is done there and the precursors can be matched here. 
    When m/z is used alone, e.g. MSI level 4, error rate may be high. 
    """
    def __init__(self, parameters):
        self.parameters = parameters
        self.outdir = parameters['outdir']
        self.project_name = parameters['project_name']

    def load_annotation_db(self):
        """
        Load database. Two flavors based on parameters['search_algorithm']. 
        The 'rtmz' format is like list of features with 'mz' and 'rtime'. 
        The 'mz' approach is usually on compound library with neutral mass, 
        where we will consider primary ions.  
        """
        
        if self.parameters['search_algorithm'] == 'mz':
            _default_database_ = 'core_blood_metabolites_20250212.json'
            _database_file = self.parameters.get('db') or str(
                pkg_resources.files(db).joinpath(_default_database_))
            list_lib_entries = self.read_cpd_entries(_database_file)
            
        elif self.parameters['search_algorithm'] == 'rtmz':
            _default_database_ = 'v3_authlib2024_hilicpos.json'
            _database_file = self.parameters.get('db') or str(
                pkg_resources.files(db).joinpath(_default_database_))
            list_lib_entries = json.load(open(_database_file))
            print(f"Read {len(list_lib_entries)} compounds from {_database_file}.")

        else:
            print(f"search_algorithm {self.parameters['search_algorithm']} is not supported in LC.")
            return None
        
        self.database_file = _database_file
        return list_lib_entries

    def read_cpd_entries(self, _database_file):
        """
        Input format example: [
                {
                "inchikey": "DDRJAANPRJIHGJ-UHFFFAOYSA-N",
                "name": "Creatinine",
                "chemical_formula": "C4H7N3O",
                "id": "MDB0000003r1",
                "monoisotopic_molecular_weight": 113.058911861,
                "RTI_RP_NORMAN": 107.23, ...
                }, ...
        ]

        Neutral mass is taken from `monoisotopic_molecular_weight` or `neutral_formula_mass`
        (or `neutral_mass`, `monoisotopic_mass`); entries without any are skipped.

        returns list of features, considering primary adducts.
        Each adduct of a compound is a "feature", in format similar to v3_authlib2024_hilicpos.json:
            {'id': 'MDB0000003r1_M+H+', 'mz': 114.066188, 'rtime': None, 'ion': 'M+H+', 'isotope': 'M0',
             'name': 'Creatinine', 'identifier': 'MDB0000003r1', 'formula': 'C4H7N3O',
             'neutral_mass': 113.058911861, 'inchikey': 'DDRJAANPRJIHGJ-UHFFFAOYSA-N'}
        Ions are M+H+ (pos) or M-H- (neg), plus adducts from adduct_search_patterns_pos/neg
        (relative to M+H+ or M-H-), named as khipu modifications (e.g. 'Na/H', 'ACN'),
        so that they are comparable to khipu results. Only singly charged M0 ions are generated.
        rtime is taken from 'rtime' if present in the entry, else None.
        """
        if self.parameters['mode'] == 'pos':
            adduct_search_patterns = adduct_search_patterns_pos
            base_ion, delta_mz = 'M+H+', PROTON
        else:
            adduct_search_patterns = adduct_search_patterns_neg
            base_ion, delta_mz = 'M-H-', -PROTON

        # Add M+H+ or M-H-
        ions = [(0, base_ion)] + [(x[0], x[1]) for x in adduct_search_patterns]

        with open(_database_file, encoding='utf-8') as f:
            entries = json.load(f)
        if isinstance(entries, dict):
            entries = list(entries.values())

        list_lib_features, num_skipped = [], 0
        for cpd in entries:
            neutral_mass = None
            for k in ('monoisotopic_molecular_weight', 'neutral_formula_mass', 'neutral_mass', 'monoisotopic_mass'):
                if cpd.get(k) not in (None, ''):
                    neutral_mass = float(cpd[k])
                    break
            if neutral_mass is None:
                num_skipped += 1
                continue
            cpd_id = cpd.get('id') or cpd.get('accession') or cpd.get('name')
            for delta, ion in ions:
                list_lib_features.append({
                    'id': f"{cpd_id}_{ion}",
                    'mz': round(neutral_mass + delta_mz + delta, 6),
                    'rtime': cpd.get('rtime'),
                    'ion': ion,
                    'isotope': 'M0',
                    'name': cpd.get('name', ''),
                    'identifier': cpd_id,
                    'formula': cpd.get('chemical_formula') or cpd.get('formula', ''),
                    'neutral_mass': neutral_mass,
                    'inchikey': cpd.get('inchikey', ''),
                })

        print(f"Read {len(entries) - num_skipped} compounds from {_database_file} ({num_skipped} without mass skipped), "
              f"as {len(list_lib_features)} library features of {len(ions)} ions.")
        return list_lib_features


    def load_tissue_specific_lib(self):
        """
        To add 'csm' etc. in the future.

        """
        pass

        
    def annotate_user_featuretable(self, infile, list_empCpds=None):
        """
        Annotate a user supplied feature table. Used by asari subcommand `annotate`.

        Parameters
        ----------
        infile : str
            asari feature table (tsv).
        list_empCpds : list, optional
            pre-annotation from `deconvolute` (khipu_deconvolute_featuretable, <project>_list_khipus.json).
            If given, library matches are combined within each empCpd and assigned to all its features,
            by combine_preannotation (more common ions prioritized in conflicts).

        Returns
        -------
        matched_results : {feature id: [library entry id, ...], ...}
        """
        with open(infile) as f:
            _n, list_features = read_features_from_asari_table(f.read())
        self.list_features = list_features
        list_lib_entries = self.load_annotation_db()
        if not list_lib_entries:
            return None
        matched_results = self.match(list_features, list_lib_entries)
        self.direct_matches = matched_results
        self.dict_feature_empCpd = {}

        if list_empCpds:
            # combine ion annotations in list_empCpds
            matched_results, self.dict_feature_empCpd = combine_preannotation(
                matched_results, list_empCpds, self.dict_lib_entries,
                mode=self.parameters['mode'], mz_tolerance_ppm=self.parameters['mz_tolerance_ppm'])

        self.export_feature_annotation(matched_results)
        return matched_results


    def match(self, list_features, list_lib_entries):
        """
        Use --algorithm 'rtmz' to match autentic library or MS/MS search results;
        --algorithm 'mz' to do mz only match.

        Parameters
        ----------
        list_features : list
            features as dicts with 'id', 'mz', 'rtime', e.g. from read_features_from_asari_table.
        list_lib_entries : list
            library entries as dicts with 'id', 'mz' (and 'rtime' for 'rtmz'),
            from load_annotation_db.

        Returns
        -------
        dict_matched : {feature id: [library entry id, ...], ...}. Unmatched features are not included.
            m/z tolerance is parameters['mz_tolerance_ppm'];
            rtime tolerance (rtmz) is parameters['lc_rt_tolerance'], in seconds.
            Library entries are kept in self.dict_lib_entries, indexed by id.
        """
        self.dict_lib_entries = {e['id']: e for e in list_lib_entries}
        mz_ppm = self.parameters['mz_tolerance_ppm']

        if self.parameters['search_algorithm'] == 'rtmz':
            lib_with_rt = [e for e in list_lib_entries if e.get('rtime') not in (None, '')]
            if len(lib_with_rt) < len(list_lib_entries):
                print(f"{len(list_lib_entries) - len(lib_with_rt)} library entries without rtime are excluded.")
            dict_matched = list_match_lcms_features(
                list_features, lib_with_rt,
                mz_ppm=mz_ppm, rt_tolerance=self.parameters['lc_rt_tolerance'])

        elif self.parameters['search_algorithm'] == 'mz':
            mapped, _, _ = all_mass_paired_mapping(
                [f['mz'] for f in list_features], [e['mz'] for e in list_lib_entries],
                std_ppm=mz_ppm)
            dict_matched = {}
            for ii, jj in mapped:
                dict_matched.setdefault(list_features[ii]['id'], []).append(list_lib_entries[jj]['id'])
            print(f"Of {len(list_features)} features, {len(dict_matched)} are matched by m/z.")

        else:
            print(f"search_algorithm {self.parameters['search_algorithm']} is not supported in LC.")
            return {}

        return dict_matched


    def export_feature_annotation(self, matched_results):
        """
        Export feature annotation to JSON and tsv files, in self.outdir:
        - Features_<project>.tsv, one row per matched pair of feature and library entry,
          in the order of the input feature table; features without match are not included.
        - Annotated_features_<project>.json, meta data and list of matched features,
          each with all feature fields and 'matches' (library entries plus ppm and delta_rtime),
          matches sorted by abs(ppm).

        Parameters
        ----------
        matched_results : dict
            {feature id: [library entry id, ...], ...}, from self.match.
            Feature details are from self.list_features, library entries from self.dict_lib_entries.
        """
        algorithm = self.parameters['search_algorithm']
        list_annotated = []
        for F in self.list_features:
            if F['id'] not in matched_results:
                continue
            direct = set(getattr(self, 'direct_matches', matched_results).get(F['id'], []))
            matches = []
            for lib_id in matched_results[F['id']]:
                E = dict(self.dict_lib_entries[lib_id])
                if lib_id in direct:
                    # ppm and delta_rtime only apply to direct matches, not those via empCpd
                    E['match'] = 'direct'
                    E['ppm'] = round((F['mz'] - E['mz']) / E['mz'] * 1e6, 2)
                    if E.get('rtime') not in (None, ''):
                        E['delta_rtime'] = round(F['rtime'] - float(E['rtime']), 2)
                else:
                    E['match'] = 'empCpd'
                matches.append(E)
            matches.sort(key=lambda x: (x['match'] != 'direct', abs(x.get('ppm', 0))))
            F2 = dict(F, matches=matches)
            if F['id'] in getattr(self, 'dict_feature_empCpd', {}):
                F2['empCpd'] = getattr(self, 'dict_feature_empCpd')[F['id']]
            list_annotated.append(F2)

        # tsv
        lib_fields = ['id', 'name', 'identifier', 'ion', 'isotope', 'formula', 'mz', 'ppm',
                      'rtime', 'delta_rtime', 'neutral_mass', 'inchikey']
        empCpd_fields = ['interim_id', 'neutral_formula_mass', 'isotope', 'modification']
        header = ['feature', 'mz', 'rtime', 'empCpd', 'empCpd_neutral_mass', 'feature_isotope', 'feature_ion',
                  'number_of_matches', 'match'] + [
                  'lib_' + x if x in ('id', 'mz', 'rtime') else x for x in lib_fields]
        tsv_file = os.path.join(self.outdir, 'Features_' + self.project_name + '.tsv')
        with open(tsv_file, 'w', encoding='utf-8') as O:
            O.write('\t'.join(header) + '\n')
            for F in list_annotated:
                for E in F['matches']:
                    emp = F.get('empCpd', {})
                    row = [F['id'], F['mz'], F['rtime']] + [emp.get(x, '') for x in empCpd_fields] + [
                           len(F['matches']), E['match']] + [E.get(x, '') for x in lib_fields]
                    O.write('\t'.join(['' if x is None else str(x) for x in row]) + '\n')

        # JSON
        json_file = os.path.join(self.outdir, 'Annotated_features_' + self.project_name + '.json')
        json_export = {
            'type': 'LCMS',
            'database': getattr(self, 'database_file', ''),
            'algorithm': algorithm,
            'mode': self.parameters['mode'],
            'mz_tolerance_ppm': self.parameters['mz_tolerance_ppm'],
            'time_stamp': ':'.join([str(x) for x in time.localtime()[1:6]]),
            'number_of_features': len(self.list_features),
            'number_of_annotated_features': len(list_annotated),
            'list_of_features': list_annotated,
        }
        if algorithm == 'rtmz':
            json_export['lc_rt_tolerance'] = self.parameters['lc_rt_tolerance']
        if getattr(self, 'dict_feature_empCpd', {}):
            json_export['preannotation'] = self.parameters.get('preannotation', '')
        with open(json_file, 'w', encoding='utf-8') as O:
            json.dump(json_export, O, indent=2, ensure_ascii=False, cls=NpEncoder)

        print("\nDone LC-MS annotation.")
        print(f"{len(list_annotated)} of {len(self.list_features)} features annotated, "
              f"{sum(len(F['matches']) for F in list_annotated)} matches, written to\n  {tsv_file}\n  {json_file}\n")


def load_khipu_empCpds(preannotation):
    """
    Load list of empCpds from `deconvolute` (LC, khipu), a JSON file <project>_list_khipus.json,
    or a deconvolution directory containing one.
    """
    if os.path.isdir(preannotation):
        files = [f for f in os.listdir(preannotation) if f.endswith('_list_khipus.json')]
        if not files:
            raise FileNotFoundError(f"No *_list_khipus.json found in {preannotation}.")
        preannotation = os.path.join(preannotation, files[0])
    with open(preannotation, encoding='utf-8') as f:
        list_empCpds = json.load(f)
    print(f"Read {len(list_empCpds)} empCpds as pre-annotation from {preannotation}.")
    return list_empCpds


def get_ion_priority(ion, mode='pos'):
    """
    Rank of an ion label (lower is more common): M+H+ (or M-H-) first, 
    then adducts in the order of adduct_search_patterns_pos/neg; others last.
    Accepts khipu labels (e.g. 'Na/H', as in read_cpd_entries) and 
    library labels (e.g. 'M+Na[1+]', as in v3_authlib2024_hilicpos.json).
    """
    if mode == 'pos':
        labels = ['M+H+'] + [x[1] for x in adduct_search_patterns_pos]
        aliases = {'M+H[1+]': 'M+H+', 'M+Na[1+]': 'Na/H', 'M+K[1+]': 'K/H',
                   'M+ACN+H[1+]': 'ACN', 'M+HCOONa+H[1+]': 'NaCOOH'}
    else:
        labels = ['M-H-'] + [x[1] for x in adduct_search_patterns_neg]
        aliases = {'M-H[1-]': 'M-H-', 'M-H[-]': 'M-H-', 'M+Cl[1-]': 'M+Cl-', 'M+Cl[-]': 'M+Cl-',
                   'M+HCOO[1-]': 'M+HCOO-', 'M+FA-H[1-]': 'M+HCOO-', 'M-2H+Na[1-]': 'Na/H'}
    ion = aliases.get(ion, ion)
    return labels.index(ion) if ion in labels else len(labels)


def combine_preannotation(matched_results, list_empCpds, dict_lib_entries, mode='pos', mz_tolerance_ppm=5):
    """
    Combine library matches with pre-annotation (empCpds from khipu, `deconvolute`), 
    so that all features in an empCpd get the same annotation.

    For each empCpd, candidate library entries are the direct matches of its features,
    whose neutral mass is consistent with empCpd neutral_formula_mass (within mz_tolerance_ppm).
    In conflicts, the best candidates are selected by, in order:
        1) M0 over isotopes in library entry;
        2) more common ions (get_ion_priority).
    All candidates tied for best (e.g. isomers) are assigned to all features in the empCpd,
    replacing their direct matches. 
    If no candidate is consistent with the empCpd neutral mass (i.e. library match and khipu disagree),
    the features keep their direct matches. Features not in any empCpd keep their direct matches.

    Parameters
    ----------
    matched_results : dict
        {feature id: [library entry id, ...], ...}, direct matches.
    list_empCpds : list
        empCpds as from khipu_deconvolute_featuretable, with 'interim_id', 'neutral_formula_mass',
        'MS1_pseudo_Spectra' (features with 'id', 'isotope', 'modification').
    dict_lib_entries : dict
        library entries indexed by id, with 'neutral_mass', 'ion', 'isotope'.
    mode : str
        ionization mode, 'pos' or 'neg'.
    mz_tolerance_ppm : float
        tolerance in checking neutral mass consistency.

    Returns
    -------
    combined_results : {feature id: [library entry id, ...], ...}
    dict_feature_empCpd : {feature id: {'interim_id', 'neutral_formula_mass', 'isotope', 'modification'}, ...},
        for all features in empCpds.
    """
    combined_results = dict(matched_results)
    dict_feature_empCpd = {}
    num_annotated_empCpds, num_inconsistent = 0, 0
    for emp in list_empCpds:
        features = emp.get('MS1_pseudo_Spectra', [])
        neutral_mass = emp.get('neutral_formula_mass')
        for F in features:
            dict_feature_empCpd[F['id']] = {
                'interim_id': emp.get('interim_id', ''),
                'neutral_formula_mass': neutral_mass,
                'isotope': F.get('isotope', ''),
                'modification': F.get('modification', ''),
            }
        candidates = []
        for F in features:
            for lib_id in matched_results.get(F['id'], []):
                if lib_id not in candidates:
                    candidates.append(lib_id)
        if not candidates:
            continue

        if not neutral_mass:
            continue
        candidates = [x for x in candidates if dict_lib_entries[x].get('neutral_mass') and 
                      abs(float(dict_lib_entries[x]['neutral_mass']) - neutral_mass) <= neutral_mass * mz_tolerance_ppm * 1e-6]
        if not candidates:
            num_inconsistent += 1
            continue

        def _rank(lib_id):
            E = dict_lib_entries[lib_id]
            return (E.get('isotope', 'M0') != 'M0', get_ion_priority(E.get('ion', ''), mode))

        ranks = {x: _rank(x) for x in candidates}
        best = min(ranks.values())
        selected = [x for x in candidates if ranks[x] == best]
        for F in features:
            combined_results[F['id']] = selected
        num_annotated_empCpds += 1

    print(f"Pre-annotation: {num_annotated_empCpds} of {len(list_empCpds)} empCpds annotated "
          f"({num_inconsistent} with matches inconsistent to empCpd neutral mass kept as direct matches); "
          f"{len(combined_results)} features annotated, from {len(matched_results)} by direct match.")
    return combined_results, dict_feature_empCpd


#
# ---------------------------------------------------------------------------
#

class GCMS_Annotation:
    """
    GC-MS annotation against a spectral library (MSP or JSON).

    Workflows, by input and parameters['search_algorithm']:
    - feature table (.tsv), 'reverse': library spectra are matched to features within
      m/z and RI windows (batch_lib_reverse_search_score). Requires high resolution data.
    - feature table (.tsv), 'forward': features are deconvoluted into pseudospectra
      (parameters['deconvolution_algorithm']), then searched forward (pseudospectra_workflow).
    - pseudospectra (.json or .msp, e.g. from `asari deconvolute`): forward search only.
    Forward search uses high resolution m/z if parameters['hires'], else unit m/z.

    RI comes from an RI column in the feature table, or from parameters['kovats'].
    """
    default_database = 'gcms_corelib_20260930.msp'

    def __init__(self, parameters):
        self.parameters = parameters
        self.outdir = parameters['outdir']
        self.project_name = parameters['project_name']
        self.database_file = parameters.get('db') or str(
            pkg_resources.files(db).joinpath(self.default_database))
        self.list_lib_entries, self.dict_lib_entries = self.load_gcms_database(
                    database_file=self.database_file,
                    low_peak_filter_factor=100,
                )

    def load_gcms_database(self, database_file, low_peak_filter_factor=100):
        """Load GC-MS Database in MSP or JSON format, as list and dict of GC_lib_entry."""
        list_lib_entries = reformat_gcms_lib(load_gcms_dbfile(database_file),
                                             filter_factor=low_peak_filter_factor)
        dict_lib_entries = {e.id: e for e in list_lib_entries}
        print(f"Imported {len(list_lib_entries)} compound library entries from {database_file}.")
        return list_lib_entries, dict_lib_entries

    def _mz_tolerance(self, highres=None):
        highres = self.parameters.get('hires', True) if highres is None else highres
        return self.parameters['ms2_tolerance_in_da'] if highres else 0.5

    # ------------------------------------------------------------------ workflows

    def table_workflow(self, infile):
        """Annotate an asari feature table (tsv). Returns list_empCpds, feature_anno_list."""
        parameters = self.parameters
        ingested = ingest_gc_featuretable(infile, parameters)
        if ingested is None:
            return None
        list_features, feature_dataframe, ri_model = ingested
        dict_features = {f['id']: f for f in list_features}

        if parameters['search_algorithm'] == 'reverse':
            if not parameters.get('hires', True):
                print("Reverse search requires high resolution data; use --search_algorithm forward.")
                return None
            list_empCpds, feature_anno_list = self.annotate_gcms_reverse(
                list_features,
                feature_dataframe,
                self.list_lib_entries,
                dict_features,
                ms2_tolerance_in_da=parameters['ms2_tolerance_in_da'],
                ri_tolerance=parameters['ri_tolerance'],
                score_cutoff_cosine=parameters['score_cutoff_cosine'],
                score_cutoff_entropy=parameters['score_cutoff_entropy'],
                corr_cutoff=parameters['corr_cutoff'],
            )
            self.export(list_empCpds, feature_anno_list, dict_features,
                        algorithm='asari-gc-reverse')

        elif parameters['search_algorithm'] == 'forward':
            # construct pseudospectra (feature groups) from scratch
            list_pseudospectra, core_features = run_gc_deconvolution(
                list_features, feature_dataframe, ri_model, parameters['deconvolution_algorithm']
                )
            json_output_file = os.path.join(
                self.outdir, self.project_name + '_denovo_pseudospectra.json')
            write_list_pseudospectra(list_pseudospectra, json_output_file)
            list_query_spectra = [port_pseudospectrum_to_json(x, normalize_intensity=True)
                                  for x in list_pseudospectra]
            list_empCpds, feature_anno_list = self.pseudospectra_workflow(      # contains export function
                list_query_spectra, dict_features=dict_features)

        else:
            print(f"search_algorithm {parameters['search_algorithm']} is not supported in GC.")
            return None

        return list_empCpds, feature_anno_list

    def pseudospectra_workflow(self, list_query_spectra, dict_features=None):
        """
        Forward search of pseudospectra against the library, then export.
        Only in parameters['search_algorithm'] == 'forward', or when input is pseudospectra.
        High resolution data are dealt differently from low resolution data (parameters['hires']).

        list_query_spectra : pseudospectra as dicts, from port_pseudospectrum_to_json
            or loaded from JSON/MSP files (normalized by standardize_query_spectra).
        dict_features : features by ID, if available (from a feature table), for export.

        Returns list_empCpds, feature_anno_list.
        """
        parameters = self.parameters
        list_query_spectra = standardize_query_spectra(list_query_spectra)
        print(f"Forward search of {len(list_query_spectra)} pseudospectra, "
              f"{'high' if parameters.get('hires', True) else 'unit'} resolution.")
        list_empCpds, feature_anno_list = self.annotate_gcms_forward(
            list_query_spectra,
            self.list_lib_entries,
            highres=parameters.get('hires', True),
            ms2_tolerance_in_da=parameters['ms2_tolerance_in_da'],
            ri_tolerance=parameters['ri_tolerance'],
            score_cutoff_cosine=parameters['score_cutoff_cosine'],
            score_cutoff_entropy=parameters['score_cutoff_entropy'],
        )
        self.export(list_empCpds, feature_anno_list, dict_features,
                    algorithm='asari-gc-forward-' + ('hires' if parameters.get('hires', True) else 'unit'))
        return list_empCpds, feature_anno_list

    # ------------------------------------------------------------------ search

    def annotate_gcms_reverse(self,
        list_features,
        feature_dataframe,
        list_lib_entries,
        dict_features,
        ms2_tolerance_in_da=0.005,
        ri_tolerance=50,
        cosine_penalty=1,           # don't change
        score_cutoff_cosine=0.5,
        score_cutoff_entropy=0.4,
        corr_cutoff=0.6,
        ):
        """
        Reverse search starting with HRMS database using def batch_lib_reverse_search_score

        returns list_empCpds, feature_anno_list
        """
        matched_results = batch_lib_reverse_search_score(
            list_lib_entries, list_features, dict_features, feature_dataframe,
            mz_tolerance_da=ms2_tolerance_in_da, ri_tolerance=ri_tolerance,
            cosine_penalty = cosine_penalty, corr_cutoff=corr_cutoff
        )
        return curate_batch_lib_reverse_search_result(
            matched_results,
            mz_tolerance_da=ms2_tolerance_in_da,
            score_cutoff_cosine=score_cutoff_cosine,
            score_cutoff_entropy=score_cutoff_entropy,
        )

    def annotate_gcms_forward(self,
        list_query_spectra,
        list_lib_entries,
        highres=True,
        ms2_tolerance_in_da=0.005,
        ri_tolerance=30,
        score_cutoff_cosine=0.5,
        score_cutoff_entropy=0.4,
        include_missing_ri=True,
        min_matched_peaks=3,
        top_n=1,
        ):
        """
        Forward search of pseudospectra against GC-MS database,
        separate algorithms for low-res and high-resolution:
        - highres: batch_forward_hires_search_score (entropy + cosine, m/z within ms2_tolerance_in_da);
        - low res: batch_forward_unit_search_score (cosine on unit m/z).

        list_query_spectra : list of pseudospectra dicts with 'id', 'RI', 'members', 'peaks'.

        returns list_empCpds, feature_anno_list
        """
        if highres:
            matched_results = batch_forward_hires_search_score(
                list_query_spectra, list_lib_entries,
                mz_tolerance_da=ms2_tolerance_in_da,
                ri_tolerance=ri_tolerance,
                include_missing_ri=include_missing_ri,
                min_matched_peaks=min_matched_peaks,
                top_n=top_n,
            )
        else:       # low resolution, using unit mz search
            matched_results = batch_forward_unit_search_score(
                convert_unit_mz(list_query_spectra),
                convert_unit_mz(list_lib_entries),
                ri_tolerance=ri_tolerance,
                include_missing_ri=include_missing_ri,
                min_matched_peaks=min_matched_peaks,
                top_n=top_n,
            )
        print(f"Forward search: {len(matched_results)} pseudospectra matched to library.")

        return curate_batch_lib_forward_search_result(
            matched_results,
            mz_tolerance_da=ms2_tolerance_in_da,
            score_cutoff_cosine=score_cutoff_cosine,
            score_cutoff_entropy=score_cutoff_entropy,
        )

    # ------------------------------------------------------------------ export

    def export(self, list_empCpds, feature_anno_list, dict_features=None, algorithm=''):
        """
        Export tsv tables (Features_, empCpds_), JSON and MSP of annotated pseudospectra,
        and mirror plots if parameters['do_mirror_plot'].
        dict_features can be None (pseudospectra input): feature details are then
        taken from the empCpds or left blank.
        """
        outdir, handle = self.outdir, self.project_name
        write_tsv_feature_anno(feature_anno_list, dict_features, self.dict_lib_entries,
                               os.path.join(outdir, "Features_" + handle + '.tsv'))
        write_tsv_empCpd_anno(list_empCpds, dict_features, self.dict_lib_entries,
                              os.path.join(outdir, "empCpds_" + handle + '.tsv'))
        print("\nDone GC-MS annotation.")
        print(f"Exported tsv results for {len(list_empCpds)} annotated empCpds and "
              f"{len(set([f['feature'] for f in feature_anno_list]))} unique features.\n")

        if self.parameters.get('do_mirror_plot', False):
            self.make_mirror_plots(list_empCpds)

        # Export list_empCpds JSON and MSP
        list_empCpds = serialize_annotated_empCpds(list_empCpds)
        json_output_file = os.path.join(outdir, 'Annotated_pseudospectra_' + handle + '.json')
        json_pseudospectra_to_msp(list_empCpds, json_output_file.replace('.json', '.msp'))
        json_export = {
            'type': 'GCMS',
            'database': self.database_file,
            'algorithm': algorithm,
            'time_stamp': ':'.join([str(x) for x in time.localtime()[1:6]]),
            'list_of_empCpds': list_empCpds
        }
        with open(json_output_file, 'w', encoding='utf-8') as O:
            json.dump(json_export, O, indent=2, ensure_ascii=False, cls=NpEncoder)
        print(f"{json_output_file} was written to JSON and MSP formats.\n")

    def make_mirror_plots(self, list_empCpds):
        """PDF mirror plot per empCpd, experimental peaks vs full library spectrum."""
        print("Exporting PDF mirror plots..\n")
        path_mirrorplots = os.path.join(self.outdir, "mirrorplots/")
        os.makedirs(path_mirrorplots, exist_ok=True)
        match_tol = self._mz_tolerance()
        for tt in list_empCpds:
            _score = f"\n(cosine score: {_fmt_score(tt['cosine_score'])}, entropy score: {_fmt_score(tt['entropy_score'])})"
            _title = tt['id'] + '__' + str(tt['name'])
            mirror_plot(tt['peaks_as_features'],
                        # peaks from original lib, not filtered by match, to show the full spectrum.
                        self.dict_lib_entries[tt['lib_entry_id']].peaks,
                        match_tol=match_tol,
                        colors=["blue", "tab:red", "black"],
                        title=_title + _score,
                        outfile=os.path.join(path_mirrorplots, sanitize_filename(_title)+'.pdf'))


def standardize_query_spectra(list_query_spectra):
    """
    Normalize pseudospectra loaded from different sources to dicts with
    'id', 'rtime', 'RI', 'members', 'peaks' (other fields are kept):
    - port_pseudospectrum_to_json / denovo JSON (already in this format);
    - JSON export of annotation ({'list_of_empCpds': [...]}, 'features' as members);
    - MSP written by json_pseudospectra_to_msp ('ID', 'RETENTIONTIME' as RI, 'Features').
    Spectra without peaks are dropped.
    """
    if isinstance(list_query_spectra, dict):
        list_query_spectra = list_query_spectra.get('list_of_empCpds', [])
    result = []
    for i, e in enumerate(list_query_spectra):
        q = dict(e)
        peaks = q.get('peaks') or q.get('peaks_as_features') or []
        if not len(peaks):
            continue
        q['peaks'] = [tuple(p[:2]) for p in peaks]
        q['id'] = q.get('id') or q.get('ID') or q.get('Name') or f'pseudospectrum_{i}'
        if q.get('RI') is None:
            q['RI'] = parse_ri_value(q.get('RETENTIONTIME'))
        members = q.get('members') or q.get('features') or q.get('Features') or []
        if isinstance(members, str):
            members = [x.strip() for x in members.split(',') if x.strip()]
        q['members'] = list(members)
        result.append(q)
    return result


# ---------------------------------------------------------------------------
#
# GC export
#

def _fmt_score(x, digits=3):
    """Round a score, or '' if missing (e.g. no entropy score in unit m/z search)."""
    return '' if x is None else round(float(x), digits)

def _ri_delta(ri, lib_ri, digits=3):
    if ri is None or lib_ri is None:
        return ''
    return round(ri - lib_ri, digits)

# write out annotated results
def write_tsv_feature_anno(feature_anno_list, dict_features, dict_lib_entries, outfile):
    """
    The empCpd output is more concise (annotation by empCpd or cluster).
    This exports a table per annotated feature.
    dict_features can be None (pseudospectra input); feature details are then blank.
    """
    header = [
        'feature', 'mz', 'rtime', 'RI',
        'empCpd', 'quant_ion', 'score_cosine', 'score_entropy', 'name', 'inchikey',
        'lib_RI', 'delta_RI', 'peak_area', 'cSelectivity', 'peak_shape', 'snr', 'detection_counts'
    ]
    dict_features = dict_features or {}
    s = '\t'.join(header) + '\n'
    for feat in feature_anno_list:
        k = feat['feature']
        FF = dict_features.get(k, {})
        ri = FF.get('RI')
        lib_RI = dict_lib_entries[feat['lib_entry_id']].RI
        s += '\t'.join([str(x) for x in [
            k, FF.get('mz', ''), FF.get('rtime', ''), _fmt_score(ri),
            feat['empCpd'], feat['quant_ion'], _fmt_score(feat['cosine_score']),
            _fmt_score(feat['entropy_score']), feat['name'], feat['inchikey'],
            _fmt_score(lib_RI), _ri_delta(ri, lib_RI),
            FF.get('peak_area', ''), FF.get('cSelectivity', ''), FF.get('goodness_fitting', ''),
            FF.get('snr', ''), FF.get('detection_counts', '')
        ]]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)

def write_tsv_empCpd_anno(list_empCpds, dict_features, dict_lib_entries, outfile):
    """Table output for list of empCpds.
    dict_features can be None (pseudospectra input); quant ion m/z and rtime are then
    taken from the empCpd ('quant_ion_mz', 'rtime') if available.
    """
    header = [
        'empCpd', 'name', 'inchikey', 'formula', 'lib_exact_mass', 'score_cosine', 'score_entropy',
        'num_matched_peaks',
        'quant_ion', 'mz_quant_ion', 'rtime_quant_ion(min)', 'RI_quant_ion', 'lib_RI', 'delta_RI',
        'has_basepeak', 'has_molecular_ion', 'presence_ratio'
    ]
    dict_features = dict_features or {}
    s = '\t'.join(header) + '\n'
    for epd in list_empCpds:
        lib = dict_lib_entries[epd['lib_entry_id']]
        # forward search reports matched peaks; in reverse search all peaks_as_features are matched
        num_matched_peaks = epd.get('num_matched_peaks') or len(epd['peaks_as_features'])
        quant_ion = epd['quant_ion']
        FF = dict_features.get(quant_ion, {})
        mz = FF.get('mz', epd.get('quant_ion_mz', ''))
        rtime = FF.get('rtime', epd.get('rtime'))
        ri = epd.get('RI')
        s += '\t'.join([str(x) for x in [
            epd['id'], epd['name'], epd['inchikey'], lib.compound_formula, lib.exact_mass,
            _fmt_score(epd['cosine_score']), _fmt_score(epd['entropy_score']),
            num_matched_peaks,
            quant_ion, mz,
            '' if rtime in (None, '') else round(rtime/60, 2),
            '' if ri is None else int(ri), '' if lib.RI is None else int(lib.RI), _ri_delta(ri, lib.RI, 0),
            epd['has_basepeak'], epd['has_molecularion'],
            round(num_matched_peaks/lib.peaks.shape[0], 2) if lib.peaks.shape[0] else ''
        ]]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)


# ---------------------------------------------------------------------------
#
# To be merged with LCMS_Annotation class
#



def annoate_by_standards(list_features, cpdLib):
    """Match features to a library of compound standards by m/z and rtime."""
    matched = match_features.list_match_lcms_features(
        list_features, cpdLib,  mz_ppm=5,  rt_tolerance=30
    )

    dict_f = {feature['id']: feature for feature in list_features}
    dict_lib = {feature['id']: feature for feature in cpdLib}


    return matched


def get_concise_annotation(feature_id, matched_list, feature_dict, lib_dict):
    """Remove duplicate annotations for the same feature, and return the concise list of fields.
    """
    feature = feature_dict[feature_id]
    result = {
        'id': feature_id,
        'mz': feature['mz'],
        'rtime': feature['rtime'],
    }
    matched_auth_libs = [x for x in matched_list if x.startswith('v2r2024')]
    matched_csm_libs = [x for x in matched_list if x.startswith('r1_')]
    # for authentic cpd lib and csm matches, keep one each
    if matched_auth_libs:
        lib_id = matched_auth_libs[0]
        lib_info = lib_dict[lib_id]
        result.update({
            'lib_id': lib_id,
            'lib_name': lib_info['name'],
            'lib_mz': lib_info['mz'],
            'lib_rtime': lib_info['rtime'],
            'lib_identifier': lib_info['identifier'],
            'lib_ion': lib_info['ion'],
            'lib_isotope': lib_info['isotope'],
        })
    else:
        result.update({
            'lib_id': '',
            'lib_name': '',
            'lib_mz': '',
            'lib_rtime': '',
            'lib_identifier': '',
            'lib_ion': '',
            'lib_isotope': '',
        })
    if matched_csm_libs:
        csm_id = matched_csm_libs[0]
        result.update({
            'CSMF_ID': csm_id,
            'CSM_ion': lib_dict[csm_id]['ion_csm'],
            'CSM_top_recommendation_name': lib_dict[csm_id]['top_recommendation_name'],
            'CSM_top_recommendation_score': lib_dict[csm_id]['top_recommendation_score'],
            'CSM_HMDB': lib_dict[csm_id]['HMDB'],
        })
    else:
        result.update({
            'CSMF_ID': '',
            'CSM_ion': '',
            'CSM_top_recommendation_name': '',
            'CSM_top_recommendation_score': '',
            'CSM_HMDB': '',
        })
    return result


def export_combined_anno_table(matched, feature_dict, lib_dict, outfile):
    """Export matched features with concise annotation to a tsv file."""
    header = ['id', 'mz', 'rtime', 'lib_id', 'lib_name', 'lib_mz', 'lib_rtime', 'lib_identifier',
            'lib_ion', 'lib_isotope', 'CSMF_ID', 'CSM_ion',
            'CSM_top_recommendation_name', 'CSM_top_recommendation_score',
            'CSM_HMDB']
    s = '\t'.join(header) + '\n'
    for feature_id, matched_list in matched.items():
        concise_info = get_concise_annotation(feature_id, matched_list, feature_dict, lib_dict)
        s += '\t'.join([str(concise_info[h]) for h in header]) + '\n'

    with open(outfile, 'w') as f:
        f.write(s)

def export_combined_anno_json(matched, feature_dict, lib_dict, outfile):
    """Export matched features and library entries to a JSON file."""
    json_output = []
    for feature_id, matched_list in matched.items():
        json_output.append(
            {feature_id: {'Feature': feature_dict[feature_id],
                          'Matched_libs': [lib_dict[lib_id] for lib_id in matched_list]}}
        )
    with open(outfile, 'w') as f:
        json.dump(json_output, f, indent=4)
