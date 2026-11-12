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

*******

GC-HRMS annotation functions are added here.
This does not change default --anno for LC-MS data processing; but we will unify in next version.

Example use for GC-HRMS annotation:
asari annotate -i /Users/lish/li.proj/IndiPHARM/LookAhead/GCMS/lookahead_gcms_asari13_asari_project_311234834/export/full_Feature_table.tsv \
    -o /Users/lish/li.proj/IndiPHARM/LookAhead/GCMS/ \
        -j v16 \
            --kovats /Users/lish/li.proj/asari_project/GCHRMS/KovatsIndex_alkanestandards_cuimc.tsv \
            --db /Users/lish/li.proj/asari_project/GCHRMS/Resources/GCHRMS_Database_251217.msp \
                --denovo F --workflow GC

LC:
asari annotate -i /Users/lish/li.play/test16_v16_599432/preferred_Feature_table.tsv -o /Users/lish/li.play/test16_v16_599432/ -j v16lc --workflow LC --mode pos

LC-MS/MS
asari annotate --db /Users/lish/li.proj/Resources/MSMS-Public_experimentspectra-pos-VS19.pkl -i /Users/lish/li.play/miao_u2c_IDX028_asari_project_6512215/export/ms2_spectra.json -o /Users/lish/li.play/miao_u2c_IDX02/ --workflow LCMSMS

"""

import os
import json
import time
import pickle
import importlib.resources as pkg_resources
from pathvalidate import sanitize_filename

from . import db
from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from .default_parameters import (extended_adducts, adduct_search_patterns_pos,
                         adduct_search_patterns_neg,
                         isotope_search_patterns)
from .utils import NpEncoder
from .tools.file_io import read_features_from_asari_table
from .tools.plot import mirror_plot
from .tools.entropy_search import clean_list_ms2spectra
from .gc_functions import *
from .deconvolute import (ingest_gc_featuretable, run_gc_deconvolution)

def annotate_project(infile, parameters):
    """Run subcommand `annotate` on infile, by parameters['workflow'] (LC, GC or LCMSMS)."""
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    subdir = '_'.join(['annotation', parameters['project_name'], ''.join(time_stamp)])
    outdir = os.path.join(parameters['outdir'], subdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"Annotation directory: {outdir}.\n\n")

    if parameters['workflow'] == "LC":
        # to add deconvolute and algorithms in version 2
        ANN = LCMS_Annotation(parameters)
        ANN.parameters['outdir'] = outdir
        ANN.annotate_user_featuretable(infile)

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

    elif parameters['workflow'] == "LCMSMS":
        # to handle diff anno methods, including pre-annotation in version 2
        annotate_lcmsms(
            infile=infile,
            outdir=outdir,
            database_file=parameters['db'],
            project_name_handle=parameters['project_name'],
            ms2_tolerance_in_da=parameters.get('lcmsms_ms2_tolerance_in_da', 0.01),

        )

    # recording project parameters as meta data
    with open(os.path.join(outdir, 'project.json'), 'w', encoding='utf-8') as f:
        json.dump(parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)


#
# -----------------------------------------------------------------------------
#

def annotate_lcmsms(
        infile,
        outdir,
        database_file,
        project_name_handle='lcmsmsresult',
        ms2_tolerance_in_da=0.01,
        cosine_similarity_filter=0.5,
        database_key='INCHIKEY'
    ):
    """
    LC-MS/MS annotation of the input JSON file on a pre-indexed database_file.

    infile : MS/MS spectra or clusters (as from Asari LCMSMS process workflow) in JSON
    database_file : pre-indexed database_file compatible with Flash Entropy
    """
    entropy_screen_params = {'mz_tol_ms1': 0.01,
        'mz_tol_ms2': ms2_tolerance_in_da,
        'ms2_sim_tol': 0.2,     # intentionally low, to filter later. ENTROPY_FILTER = 0.2
        'precursor_mz_offset': 1.6,
        }

    entropy_search = pickle.load(open(database_file, 'rb'))
    print(f"Loaded database of {len(entropy_search.precursor_mz_array)} spectra.")
    ms2_spectra = json.load(open(infile))
    cleaned_spectra = clean_list_ms2spectra(ms2_spectra, entropy_search, entropy_screen_params)
    print(f"Loadded {len(cleaned_spectra)} cleaned experimental spectra.\n")

    results = []
    for (precursor_mz, peaks, rtime, spid) in cleaned_spectra:
        similarity, matched_num = entropy_search.identity_search(precursor_mz=precursor_mz,
                                                                peaks=peaks,
                                                                ms1_tolerance_in_da=entropy_screen_params['mz_tol_ms1'],
                                                                ms2_tolerance_in_da=entropy_screen_params['mz_tol_ms2'],
                                                                output_matched_peak_number=True)
        if np.max(similarity) > entropy_screen_params['ms2_sim_tol']:
            # get best entropy match
            best_idx = np.argmax(similarity)
            entropy_matched_entry = entropy_search[best_idx]
            entropy_score = similarity[best_idx]
            number_matched_peaks = matched_num[best_idx]
            # get cosine score
            cosine_score, _ = cosine_similarity(peaks, entropy_matched_entry['peaks'],
                                            mz_tolerance=ms2_tolerance_in_da)
            if cosine_score > cosine_similarity_filter:
                results.append( (entropy_matched_entry, cosine_score, entropy_score, number_matched_peaks, spid, precursor_mz, rtime) )

    # Collapse redundant DB entries for export, selecting by best cosine_score. Consider to keep top N later.
    concise_results = {}
    for res in results:
            k = res[0][database_key]
            if k in concise_results:
                if res[1] > concise_results[k][1]:
                    concise_results[k] = res
            else:
                concise_results[k] = res

    #
    json_concise_results = []
    for k, res in concise_results.items():
        # if not isinstance(db_entry['peaks'], list): db_entry['peaks'] = db_entry['peaks'].tolist()
        _ = res[0].pop('peaks')         # exclude peaks
        json_concise_results.append({
            'spid': res[4],
            'mz': round(res[5], 4),     # expt_precursor_mz
            'rtime': res[6],            # expt_rtime
            'entropy_score': round(float(res[2]), 3),
            'cosine_score': round(float(res[1]), 3),
            'number_matched_peaks': int(res[3]),
            'database_key': k,
            'database_entry': res[0],
        })

    tsv_output_file = os.path.join(outdir, project_name_handle+'_MS2_annotation.tsv')
    keys1 = ['spid', 'mz', 'rtime', 'entropy_score', 'cosine_score', 'number_matched_peaks']
    keys2 = list(json_concise_results[0]['database_entry'].keys())
    s = '\t'.join(keys1 + keys2) + '\n'
    for row in json_concise_results:
        s += '\t'.join([str(row.get(k, '')) for k in keys1]
                                  + [str(row['database_entry'].get(k, '')) for k in keys2]
                                  ) + '\n'
    with open(tsv_output_file, 'w') as fh:
        fh.write(s)

    json_export = {
        'type': 'MS/MS',
        'database': database_file,
        'algorithm': 'entropy,cosine',
        'time_stamp': ':'.join([str(x) for x in time.localtime()[1:6]]),
        'list_of_features': json_concise_results
    }
    json_output_file = os.path.join(outdir, project_name_handle+'_MS2_annotation.json')
    with open(json_output_file, 'w', encoding='utf-8') as O:
        json.dump(json_export, O, indent=2, ensure_ascii=False)

    print(f"MS/MS annotation results ({len(json_concise_results)}) in JSON and tsv were written to {outdir}.\n\n")


#
# -----------------------------------------------------------------------------
#

class LCMS_Annotation:
    """Annotation of LC-MS feature tables."""
    def __init__(self, parameters):
        self.parameters = parameters

    def annotate_user_featuretable(self, infile):
        """
            Annotate a user supplied feature table. Used by asari subcommand `annotate`.

        Parameters
        ----------
        infile : str
            input feature table filepath, tab delimited file with first row as header,
            first column m/z and 2nd column rtime.
        parameters : dict
            parameter dictionary passed from main.py,
            which imports from default_parameters and updates the dict by user arguments.
        rtime_tolerance : int, optional, default: 2
            retention time tolerance to group adducts etc.

        Outputs
        -------
        two files in current directory, Feature_annotation.tsv and Annotated_empricalCompounds.json


        Note
        ----
        Annotate features via JMS (jms.dbStructures) and khipu.
        The pre-annotation step is khipu based emprical compound construction, followed by
        three steps of annotation:

        1) Search known compound database, via neutral mass inferred by khipu

        2) Search singletons for a formula match

        3) Encapsulate remaining features in empCpd format, so that all are exported consistently.

        Export Feature_annotation as tsv, Annotated_empricalCompounds in both JSON and pickle.
        Reference databases can be pre-loaded.
        Measured m/z values are calibrated to database based values (db_mass_calibrate).
        """
        _n, list_features = read_features_from_asari_table(open(infile).read())
        for f in list_features:
                f['representative_intensity'] = f['peak_area']
        list_features.sort(key=lambda x: x['peak_area'], reverse=True)
        self.list_features = list_features

        self.load_annotation_db()
        self.db_mass_calibrate()

        # asari uses seconds for rt
        EED = ExperimentalEcpdDatabase(mode=self.parameters['mode'],
                                       mz_tolerance_ppm=self.parameters['mz_tolerance_ppm'],
                                       rt_tolerance=self.parameters['khipu_rtime_tolerance']
                                       )
        # passing patterns from .default_parameters
        if self.parameters['mode'] == 'pos':
            EED.adduct_patterns = adduct_search_patterns_pos
        else:
            EED.adduct_patterns = adduct_search_patterns_neg
        EED.isotope_search_patterns = isotope_search_patterns
        EED.extended_adducts = extended_adducts

        EED.build_from_list_peaks(list_features)
        # It takes three steps to take care of all features. First khipu organized empCpds
        EED.extend_empCpd_annotation(self.KCD)
        # Second, singletons that get a formula match in KCD
        EED.annotate_singletons(self.KCD)
        # Third, the remaining features unmatched to anything (orphans). Exported for potential downstream work.
        EED.dict_empCpds = self.append_orphans_to_epmCpds(EED.dict_empCpds)

        self.export_peak_annotation(EED.dict_empCpds, self.KCD, 'Feature_annotation')

        # export JSON
        outfile = os.path.join(self.parameters['outdir'], 'Annotated_empiricalCompounds.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(EED.dict_empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)
        print(f"JSON version was written to {outfile}.\n\n")


    def load_annotation_db(self, src='hmdb4'):
        """
        Load database of known compound using jms.dbStructures.knownCompoundDatabase.
        The compound tree is precomputed indexing.
        The `src` parameter is not used now, but placeholder to add more options later.

        Parameters
        ----------
        src: str, optional, default: hmdb4
            not used but can, in the future, dictate which database is used to generate annotations


        """
        self.KCD = knownCompoundDatabase()
        self.KCD.mass_indexed_compounds = pickle.load(
            pkg_resources.open_binary(db, 'mass_indexed_compounds.pickle') )
        self.KCD.emp_cpds_trees = pickle.load(
            pkg_resources.open_binary(db, 'emp_cpds_trees.pickle') )

    def db_mass_calibrate(self, max_features=1000, required_calibrate_threshold=0.000002):
        """
        Use KCD.evaluate_mass_accuracy_ratio to check systematic mass shift,
        which is calculated as the average ppm difference between measured m/z and theoretical values.
        If greater than required_calibrate_threshold (default 2 ppm),
        calibrate m/z values for the whole experiment by updating self.CMAP.FeatureList.

        Parameters
        ---------
        required_calibrate_treshold: float, optional, default: 0.000002
            if the mass shift exceeds this value, mass correction will be applied.

        Note
        ----
        Data format in good_reference_landmark_peaks:
        [{'ref_id_num': 99, 'apex': 211, 'height': 999999}, ...],
        where ref_id_num is index number of mass track in MassGrid.
        """
        mz_landmarks = [f['mz'] for f in self.list_features[:max_features]]
        mass_accuracy_ratio = self.KCD.evaluate_mass_accuracy_ratio(
            mz_landmarks, mode=self.parameters['mode'], mz_tolerance_ppm=10)
        if mass_accuracy_ratio:
            if abs(mass_accuracy_ratio) > required_calibrate_threshold:
                print("Mass shift is greater than %2.1f ppm. Correction applied."
                      %(required_calibrate_threshold*1000000))
                _correction = mass_accuracy_ratio + 1
                for F in self.list_features:
                    F['mz'] = F['mz'] / _correction
                    F['mz_corrected_by_division'] = _correction
        else:
            print("Mass accuracy check is skipped, too few mz_landmarks (%d) matched."
                  %len(mz_landmarks))


    def append_orphans_to_epmCpds(self, dict_empCpds):
        """
        This is the third step of feature annotation in self.annotate,
        to encapsulate features without annotation in empCpd format.
        Input via dict_empCpds and returns updated dict_empCpds.
        See also: annotate

        Parameters
        ----------
        dict_empCpds: dict
            a dictionary of empirical compounds in empCpd format

        """
        all_feature_ids = []
        for _, V in dict_empCpds.items():
            all_feature_ids += [peak['id_number'] for peak in V['MS1_pseudo_Spectra']]

        orphans = [peak for peak in self.list_features if peak['id_number'] not in all_feature_ids]
        # will need better tracking of empCpd ID; already used in jms.dbStructures new_id_start = len(self.dict_empCpds) + 10000
        new_id_start = len(dict_empCpds) + 100000
        for peak in orphans:
            dict_empCpds[new_id_start] = {'interim_id': new_id_start,
                    'neutral_formula_mass': '',
                    'neutral_formula': '',
                    'MS1_pseudo_Spectra': [peak],
                    'ion_relation': None,
                    'modification': None}
            new_id_start += 1

        return dict_empCpds

    def export_peak_annotation(self, dict_empCpds, KCD, export_file_name_prefix):
        """
        Export feature annotation to tab delimited tsv file, where interim_id is empCpd id.

        Parameters
        ----------
        dict_empCpds : dict
            dictionary of empirical compounds, using interim_id as key, as seen in JMS.
        KCD : KnownCompoundDatabase instance
            the known compound database that was used in annotating the empirical compounds.
        export_file_name_prefix : str
            to used in output file name.
        """
        s = "[peak]id_number\tmz\trtime\tapex(scan number)\t[EmpCpd]interim_id\
            \t[EmpCpd]ion_relation\tneutral_formula\tneutral_formula_mass\
            \tname_1st_guess\tmatched_DB_shorts\tmatched_DB_records\n"

        for _, V in dict_empCpds.items():
            name_1st_guess, matched_DB_shorts, matched_DB_records = '', '', ''
            if 'list_matches' in V:
                list_matches = V['list_matches']
                if list_matches:
                    name_1st_guess = KCD.mass_indexed_compounds[list_matches[0][0]]['compounds'][0]['name']
                    matched_DB_shorts = ", ".join([ "(" + KCD.short_report_emp_cpd(xx[0]) + ")"  for xx in list_matches])
                    matched_DB_records = ", ".join([str(xx) for xx  in list_matches])

            for peak in V['MS1_pseudo_Spectra']:
                s += '\t'.join([str(x) for x in [
                    peak['id_number'], peak['mz'], peak['rtime'], peak['apex'],
                    V['interim_id'], peak.get('ion_relation', ''),
                    V['neutral_formula'], V['neutral_formula_mass'],
                    name_1st_guess, matched_DB_shorts, matched_DB_records]]) + "\n"

        outfile = os.path.join(self.parameters['outdir'], export_file_name_prefix + '.tsv')
        with open(outfile, encoding='utf-8', mode='w') as O:
            O.write(s)

        print("\nAnnotation of %d Empirical compounds was written to %s.\n\n" %(len(dict_empCpds), outfile))

    def select_unique_compound_features(self, dict_empCpds):
        """
        Get unique feature by highest composite peak area per empirical compound.
        One may consider alternatives to select the peak representing an empirical compound,
        e.g. by SNR or M+H, M-H ions. This is can be done separately on the exported files.

        Parameters
        ----------
        dict_empCpds : dict
            dictionary of empirical compounds, using interim_id as key, as seen in JMS.

        Not used now.


        """
        self.selected_unique_features = {}
        for interim_id, V in dict_empCpds.items():
            if len(V['MS1_pseudo_Spectra']) == 1:
                self.selected_unique_features[V['MS1_pseudo_Spectra'][0]['id_number']] = (
                        # empCpd id, neutral_formula, ion_relation
                        interim_id, V['neutral_formula'], 'singleton'
                )
            else:
                try:
                    all_peaks = sorted([(peak['peak_area'], peak['goodness_fitting'], peak)
                                        for peak in V['MS1_pseudo_Spectra']],
                                        reverse=True)
                    best_peak = all_peaks[0][2]
                except TypeError:
                    best_peak = V['MS1_pseudo_Spectra'][0]
                self.selected_unique_features[best_peak['id_number']] = (
                    # empCpd id, neutral_formula, ion_relation (will change not always anchor)
                    interim_id, V['neutral_formula'], best_peak.get('ion_relation', '')
                )


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
# merging with LCMS_Annotation class
#


def annotate_user_featuretable(infile, parameters):
    """

    Annotate a user supplied feature table. Used by asari subcommand `annotate`.

    Parameters
    ----------
    infile : str
        input feature table filepath, tab delimited file with first row as header,
        first column m/z and 2nd column rtime.
    parameters : dict
        parameter dictionary passed from main.py,
        which imports from default_parameters and updates the dict by user arguments.
    rtime_tolerance : int, optional, default: 2
        retention time tolerance to group adducts etc.

    Outputs
    -------
    _PreAnnotated_list_khipus.json

    # currently redundant w/ --anno, but to be unified in next version.
    """
    from khipu.extended import peaklist_to_khipu_list

    mode = parameters['mode']
    if mode == 'pos':
        adduct_search_patterns = adduct_search_patterns_pos
    else:
        adduct_search_patterns = adduct_search_patterns_neg
    _n, list_features = read_features_from_asari_table(open(infile).read())
    for f in list_features:
            f['representative_intensity'] = f['peak_area']
    list_khipus, all_assigned_fids = peaklist_to_khipu_list(
                            list_features,
                            isotope_search_patterns=isotope_search_patterns,
                            adduct_search_patterns=adduct_search_patterns,
                            extended_adducts=extended_adducts,
                            mz_tolerance_ppm=parameters['mz_tolerance_ppm'],
                            rt_tolerance=parameters['khipu_rtime_tolerance'],
                            mode=mode,
                            charges=[1, 2, 3],
                            )
    #  exporting JSON
    json_output_file = os.path.join(parameters['outdir'],
                                    parameters['project_name']+'_PreAnnotated_list_khipus.json')
    with open(json_output_file, 'w', encoding='utf-8') as f:
        json.dump(list_khipus, f, cls=NpEncoder, ensure_ascii=False, indent=2)


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
