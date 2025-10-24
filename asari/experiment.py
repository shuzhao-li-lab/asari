import os
import sys
import json
import pickle
import numpy as np
import pandas as pd
from .deconvolution import AsariProcessor
from scipy import interpolate
from .chromatograms import __hacked_lowess__
from statsmodels.nonparametric.smoothers_lowess import lowess
import tqdm

from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from .default_parameters import adduct_search_patterns, \
    adduct_search_patterns_neg, isotope_search_patterns, extended_adducts, \
    readme_doc_str
from .mass_functions import all_mass_paired_mapping
from .constructors import CompositeMap
from .json_encoder import NpEncoder

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from . import db

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


class ext_Experiment:
    '''
    Similar to metDataModel.core.Experiment with preprocessing methods.

    This encapsulates a set of LC-MS files using the same experimental method 
    (chromatography and ionization) to be processed together.

    E.g., data from postive ESI and negative ESI should not in the same ext_Experiment instance.

    This class has annotation and export functions.
    
    Default asari work flow is in `ext_Experiment.process_all`.
    '''
    def __init__(self, sample_registry, parameters):
        '''
        This is the overall container for all data in an experiment/project.

        Parameters
        ----------
        sample_registry : dict
            dictionary of sample and selected data after mass track extraction.
            The bulk mass track data are not kept in memory unless specified so.
        parameters : dict
            processing parameters passed from main.py.

        Updates
        -------
        Major class attributes including self.number_of_samples, number_scans, reference_sample_id.
        '''
        self.sample_registry = sample_registry
        self.valid_sample_ids = self.get_valid_sample_ids()
        self.number_of_samples = len(self.valid_sample_ids)
        self.number_scans = self.get_max_scan_number(sample_registry)
        self.all_samples = self.all_sample_instances = []

        self.parameters = parameters
        self.output_dir = parameters['outdir']
        self.mode = parameters['mode']
        self.mz_tolerance_ppm = self.parameters['mz_tolerance_ppm']
        self.check_isotope_ratio = self.parameters['check_isotope_ratio']
        self.database_mode = parameters['database_mode']
        self.reference_sample_id = self.get_reference_sample_id()
        self.list_sample_indices = {}
        
    def get_reference_sample_instance(self):
        return self.__find_reference()[0]
    
    def get_reference_sample_id(self):
        return self.__find_reference()[1]

    def __find_reference(self):
        if self.parameters['reference']:
            # match file name; k is sm['sample_id']
            for k, v in self.sample_registry.items():
                if os.path.basename(self.parameters['reference']) == os.path.basename(v['input_file']):
                    return v, k
                elif os.path.basename(self.parameters['reference']) + ".mzML" == os.path.basename(v['input_file']):
                    return v, k

        elif self.sample_registry:
            L = [(v['number_anchor_mz_pairs'], v['sample_id']) for v in self.sample_registry.values() if v['number_anchor_mz_pairs'] is not None]
            L.sort(reverse=True)
            ref = self.sample_registry[L[0][1]]
            self.parameters['reference'] = ref['input_file']
            print("\n    The reference sample is:\n    ||* %s *||\n" %ref['name'])
            print("Max reference retention time is %4.2f at scan number %d.\n" %(
                max(ref['list_retention_time']), ref['max_scan_number']))
            return ref, ref['sample_id']
        else:
            return None, None

    def get_valid_sample_ids(self):
        '''
        Get valid sample ids, as some samples may not be extracted successfully.
        '''
        return [k for k, v in self.sample_registry.items() if v['status:eic'] == 'passed']

    def determine_acquisition_order(self):
        sample_order_by_timestamp = [(k, v['acquisition_time']) for k,v in self.sample_registry.items()]
        try:
            # fix this better in future
            return sorted(sample_order_by_timestamp, key=lambda x: x[1])
        except:
            return sorted(sample_order_by_timestamp, key=lambda x: x[0])
        
    def associate_stds_batches(self, batch_map):
        found = set()
        association = {}
        for sample_name, batch_id in batch_map.items():
            for sample_id , sample_info in self.sample_registry.items():
                if sample_id not in found and sample_info['name'] == sample_name:
                    sample_info['batch_id'] = batch_id
                    association[sample_id] = batch_id
                    found.add(sample_id)
        return association
    
    def associate_stds_samples(self, sample_run_order):
        association = {}
        current_non_RI_samples = []
        last_reference = None
        for sample_id, runtime in sample_run_order:
            sample_name = self.sample_registry[sample_id]['name']
            RI_list = pd.read_csv(self.parameters['retention_index_standards'])
            if sample_name in RI_list.columns:
                last_reference = sample_id
                for non_RI_sample in current_non_RI_samples:
                    association[non_RI_sample] = sample_id
                current_non_RI_samples = []
            else:
                current_non_RI_samples.append(sample_id)
        for non_RI_sample in current_non_RI_samples:
            association[non_RI_sample] = last_reference
        return association
    
    def get_max_scan_number(self, sample_registry):
        '''
        Return max scan number among samples, or None if no valid sample.

        Parameters
        ----------
        sample_registry: dict
            a dict that maps sample IDs to sample data
        '''

        # todo - why does this function need to take sample_registry as an external argument vs. self.sample_registry?
        if sample_registry:
            return max([sample_registry[k]['max_scan_number'] for k in self.valid_sample_ids]) + 1
        else:
            return None

    def process_all_LC(self):
        '''
        This is the default asari workflow.
        
        1. Build MassGrid, using either pairwise (small study) or clustering method. 
           Choose one reference from all samples for the largest number of landmark m/z tracks.
        2. RT alignment via a LOWESS function, using selective landmark peaks.
        3. Build composite elution profile (composite_mass_tracks),
           by cumulative sum of mass tracks from all samples after RT correction.
        4. Global peak detection is performed on each composite massTrack.
        5. Mapping global peaks (i.e. features) back to all samples and extract sample specific peak areas.
           This completes the FeatureTable.

        Updates
        -------
        self.CMAP as instance of CompositeMap, and MassGrid, composite map and features within.
        '''
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid()
        if not self.parameters['rt_align_on']:
            self.CMAP.mock_rentention_alignment()
        self.CMAP.build_composite_tracks()
        self.CMAP.global_peak_detection()

    def process_all_LC_start(self):
        '''
        This is the default asari workflow.
        
        1. Build MassGrid, using either pairwise (small study) or clustering method. 
           Choose one reference from all samples for the largest number of landmark m/z tracks.
        2. RT alignment via a LOWESS function, using selective landmark peaks.
        3. Build composite elution profile (composite_mass_tracks),
           by cumulative sum of mass tracks from all samples after RT correction.
        4. Global peak detection is performed on each composite massTrack.
        5. Mapping global peaks (i.e. features) back to all samples and extract sample specific peak areas.
           This completes the FeatureTable.

        Updates
        -------
        self.CMAP as instance of CompositeMap, and MassGrid, composite map and features within.
        '''
        raise NotImplementedError()
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid()
        self.CMAP.START()
        self.CMAP.global_peak_detection()


    def populate_RI_lookup(self, sample_map):
        RI_maps, RI_models, reverse_RI_models = {}, {}, {}
        RI_list = pd.read_csv(self.parameters['retention_index_standards'])
        RI_list = RI_list.sort_values(by=list(sample_map.values())[0])
        standard_indexes = RI_list['Index'].astype(int).tolist()

        for reference_id in tqdm.tqdm(list(sample_map.keys())):
            reference_instance = self.sample_registry[reference_id]
            RI_maps[reference_id] = {}
            RTs, indexes, scan_nos = [], [], []

            # Convert both feature and standard RTs to minutes
            standard_RTs = (RI_list[sample_map[reference_id]].astype(float)).tolist()
            num_standards = len(standard_RTs)
            list_rt_minutes = [rt / 60.0 for rt in reference_instance['list_retention_time']]

            standard_idx = 0
            for rt, scan_no in zip(list_rt_minutes, reference_instance['list_scan_numbers']):
                RTs.append(rt)
                scan_nos.append(scan_no)
                while standard_idx < num_standards and rt > standard_RTs[standard_idx]:
                    standard_idx += 1
                if standard_idx == 0:
                    prev_index, prev_rt = 0, 0.0
                    next_index, next_rt = standard_indexes[0], standard_RTs[0]
                elif standard_idx == num_standards:
                    prev_index, prev_rt = standard_indexes[-1], standard_RTs[-1]
                    next_index = standard_indexes[-1]
                    next_rt = max(list_rt_minutes) * 1.1
                else:
                    prev_index, prev_rt = standard_indexes[standard_idx - 1], standard_RTs[standard_idx - 1]
                    next_index, next_rt = standard_indexes[standard_idx], standard_RTs[standard_idx]

                if next_rt == prev_rt:
                    RI_value = next_index
                else:
                    RI_value = 100 * (prev_index + (rt - prev_rt) / (next_rt - prev_rt))
                indexes.append(RI_value)
                RI_maps[reference_id][rt] = RI_value

            model = lowess(indexes, RTs)
            model2 = lowess(RTs, indexes)
            x1, y1 = zip(*model)
            RI_models[reference_id] = interpolate.interp1d(x1, y1, fill_value="extrapolate", bounds_error=False)
            x2, y2 = zip(*model2)
            reverse_RI_models[reference_id] = interpolate.interp1d(y2, x2, fill_value="extrapolate", bounds_error=False)

        self.RI_models = RI_models
        self.reverse_RI_models = reverse_RI_models
        for i, SM in enumerate(self.all_samples):
            SM.list_retention_index = self.RI_models[i]([rt/60 for rt in SM.list_retention_time])

    def determine_batches(self):
        sample_to_batch = {}
        sample_names = set([v['name'] for _, v in self.sample_registry.items()])
        retention_index_information = pd.read_csv(self.parameters['retention_index_standards'])
        for column in retention_index_information:
            column = column[-3:].rstrip()
            for sample in sample_names:
                #print(column, sample, column in sample)
                if column in sample and sample not in sample_to_batch:
                    sample_to_batch[sample] = column
                elif column in sample and sample in sample_to_batch:
                    print(f"Warning: Sample {sample} is present in multiple batches.  This may cause issues with alignment.")
                    print("Ignoring duplicate, assuming first encountered is first batch")
        for sample in sample_names:
            if sample not in sample_to_batch:
                print(f"Warning: {sample} unmapped!")
        return sample_to_batch

    def process_all_GC(self):
        print("Processing using GC Workflow...")
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid()
        retention_index_information = pd.read_csv(self.parameters['retention_index_standards'])
        sample_names = set([(k, v['name']) for k,v in self.sample_registry.items()])
        tripped = False
        for column in retention_index_information:
            if column in sample_names:
                print("Assuming RI Samples are Present - Will Align with Sample Method")
                print("RI Samples Method is Not As Tested!!!")
                tripped = True
                break
        if not tripped:
            print("Assuming Columns are Batch Substrings - Will Align with Batch Method")

        sample_run_order = self.determine_batches()
        mapping = self.associate_stds_batches(sample_run_order)
        self.mapping = mapping
        self.populate_RI_lookup(mapping)
        self.CMAP.build_composite_tracks_GC()
        self.CMAP.global_peak_detection()

    def export_all(self, anno=True, mode="LC"):
        '''
        Export all files.
        Annotation of features to empirical compounds is done here.

        Parameters
        ----------
        anno: bool, optional, default: True
            if true, generate annotation files, export CMAP pickle and do QC plot;
            else skip annotating.
        '''
        if self.parameters['workflow'] == "LC" or mode.rstrip().upper() == "LC":
            self.CMAP.MassGrid.to_csv(
                os.path.join(self.parameters['outdir'], 'export', self.parameters['mass_grid_mapping']) )
            if anno:
                for peak in self.CMAP.FeatureList:
                    peak['id'] = str(peak['id_number'])
                self.export_CMAP_pickle()
                self.annotate()
                self.generate_qc_plot_pdf()
            self.export_feature_tables()
            self.export_log()
            self.export_readme()
        elif self.parameters['workflow'] == "GC" or mode.rstrip().upper() == "GC":
            self.export_feature_tables()
            self.CMAP.MassGrid.to_csv(
                os.path.join(self.parameters['outdir'], 'export', self.parameters['mass_grid_mapping']))
            if self.parameters['anno']:
                deconvolution = AsariProcessor(self.parameters)
                deconvolution.output_path = os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table']).replace(".tsv", "_annotated_compounds.json")
                results = deconvolution.default_workflow(os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table']))
            self.export_gc_annotations(results, os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table']))
            self.export_gc_annotations(results, os.path.join(self.parameters['outdir'], 'preferred_'+self.parameters['output_feature_table']))
            self.export_all_gc_annotations(results, os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table']))
            self.export_all_gc_annotations(results, os.path.join(self.parameters['outdir'], 'preferred_'+self.parameters['output_feature_table']))

            self.export_log()
            self.export_readme()
        else:
            raise Exception("Invalid workflow provided, specify LC or GC workflow")


    def export_all_gc_annotations(self, results, table):
        return self.export_gc_annotations(results, table, mode='all', max_ri_error=0)

    def export_gc_annotations(self, results, table, mode='best', max_ri_error=None):
        """
        Exports GC annotations with selectable modes and retention index filtering.

        Modes:
            - 'best': Exports one row per feature, annotated with the single best
                    library match (based on score * matched_peaks).
            - 'all':  Exports a separate row for *every* library match found
                    for a compound.

        RI Filter:
            - If max_ri_error (float or int) is provided, library matches will be
            filtered. Only matches where:
                abs(feature_RI - library_RI) <= max_ri_error
            will be kept. Matches without a valid numeric library_RI are kept.
        """
        
        # --- 1. Common Setup: Load features and calculate RI ---
        try:
            features = pd.read_csv(table, sep="\t")
        except FileNotFoundError:
            print(f"Error: Feature table file not found at {table}")
            return
            
        master_reference_sample = self.get_reference_sample_id()
        if master_reference_sample not in self.RI_models:
            print(f"Error: RI model for reference sample '{master_reference_sample}' not found.")
            return
            
        # Calculate RI for all features
        features['RI'] = self.RI_models[master_reference_sample](features['rtime'] / 60)
        # Create a lookup dictionary from the feature table
        ft_dict = {x['id_number']: x for x in features.to_dict(orient='records')}

        # --- 2. Initialize Output Containers ---
        all_annotation_records = []  # Used for 'all' mode
        feature_final_annotations = {} # Used for 'best' mode
        
        # Check if RI filtering is active
        ri_filter_active = isinstance(max_ri_error, (float, int))
        if not ri_filter_active:
            exp_max_ri_error = self.parameters['max_retention_index_error']
            if isinstance(exp_max_ri_error, (float, int)):
                max_ri_error = exp_max_ri_error
                ri_filter_active = True

        # --- 3. Iterate Through Compounds and Process Annotations ---
        for _ii, (id, cpd) in enumerate(tqdm.tqdm(results.items(), desc=f"Processing Compounds (Mode: {mode})")):
            cpd_id = f"CPD_{_ii}"
            
            # --- 3a. Collect all feature IDs for this compound ---
            all_features_in_compound = []
            if "anchor_fragment" not in cpd: continue # Skip malformed compound
            
            all_features_in_compound.append(cpd["anchor_fragment"]["id"])
            for isotope in cpd["anchor_fragment"].get("isotopes", []):
                all_features_in_compound.append(isotope["id"])
            for fragment in cpd["fragments"]:
                all_features_in_compound.append(fragment["id"])
                for isotope in fragment.get("isotopes", []):
                    all_features_in_compound.append(isotope['id'])

            # --- 3b. Prepare base annotations (Compound_ID, Isotopologue, etc.) ---
            base_feature_annotations = {}
            anchor_mz = cpd["anchor_fragment"].get("mz", 0)

            # Anchor fragment
            iso_id = "M0" if cpd["anchor_fragment"].get("isotopes") else "M0?"
            base_feature_annotations[cpd["anchor_fragment"]["id"]] = {
                "Compound_ID": cpd_id, "Isotopologue": iso_id, "Fragment": "quant",
            }
            for isotope in cpd["anchor_fragment"].get("isotopes", []):
                base_feature_annotations[isotope["id"]] = {
                    "Compound_ID": cpd_id, "Isotopologue": isotope["evidence"]["isotope_type"], "Fragment": "quant"
                }
            
            # Other fragments
            for fragment in cpd["fragments"]:
                frag_mz = fragment["mz"]
                frag_tag = f"fragment_{round(anchor_mz - frag_mz, 3)}"
                iso_id = "M0" if fragment.get("isotopes") or fragment.get('fragments') else "M0?"
                base_feature_annotations[fragment["id"]] = {
                    "Compound_ID": cpd_id, "Isotopologue": iso_id, "Fragment": frag_tag,
                }
                for isotope in fragment.get("isotopes", []):
                    base_feature_annotations[isotope['id']] = {
                        "Compound_ID": cpd_id, "Isotopologue": isotope["evidence"]["isotope_type"], "Fragment": frag_tag
                    }
            
            # --- 3c. Get and Filter Library Annotations ---
            annotations_list = cpd.get("annotations", {}).get("gc_msms_search", {}).get("top_matches", [])
            
            # Get the compound's RI from its anchor feature
            anchor_id = cpd["anchor_fragment"]["id"]
            feature_ri = ft_dict.get(anchor_id, {}).get('RI')
            
            filtered_annotations = []
            if ri_filter_active and feature_ri is not None:
                for annot in annotations_list:
                    library_ri_str = annot.get("retention_index")
                    try:
                        # Try to convert library RI to a float
                        lib_ri_float = float(library_ri_str)
                        # If successful, apply the filter
                        if abs(feature_ri - lib_ri_float) <= max_ri_error:
                            filtered_annotations.append(annot)
                    except (ValueError, TypeError):
                        # If library_RI is "not reported", None, or other non-numeric,
                        # we cannot filter it, so we keep it.
                        filtered_annotations.append(annot)
            else:
                # RI filter is not active, or feature has no RI, so keep all
                filtered_annotations = annotations_list

            # --- 3d. Apply Mode-Specific Logic ---
            
            # Define a placeholder for compounds with no (or no passing) matches
            placeholder_annot = {
                "name": "?", "library": "", "score": 0, "matched_peaks": 0, "retention_index": np.nan
            }
            
            if mode == 'all':
                # Use placeholders if list is empty, otherwise use filtered list
                annots_to_process = filtered_annotations if filtered_annotations else [placeholder_annot]
                
                for feature_id in all_features_in_compound:
                    if feature_id in ft_dict:
                        feature_base_data = ft_dict[feature_id].copy()
                        feature_base_data.update(base_feature_annotations.get(feature_id, {}))
                        
                        for annot in annots_to_process:
                            annot_record = feature_base_data.copy()
                            annot_record.update({
                                "Annotation_Name": annot["name"],
                                "Annotation_Source": annot["library"],
                                "Annotation_Score": annot["score"],
                                "Annotation_Matched_Peaks": annot["matched_peaks"],
                                "Annotation_Library_Index": annot.get("retention_index", np.nan)
                            })
                            all_annotation_records.append(annot_record)

            elif mode == 'best':
                best_score = -1
                best_annot_data = {
                    "Annotation_Name": "?",
                    "Annotation_Source": "",
                    "Annotation_Score": 0,
                    "Annotation_Matched_Peaks": 0,
                    "Annotation_Library_Index": np.nan
                }

                for annot in filtered_annotations:
                    current_score = annot.get('score', 0) * annot.get('matched_peaks', 0)
                    if current_score > best_score:
                        best_score = current_score
                        best_annot_data = {
                            "Annotation_Name": annot["name"],
                            "Annotation_Source": annot["library"],
                            "Annotation_Score": annot["score"],
                            "Annotation_Matched_Peaks": annot["matched_peaks"],
                            "Annotation_Library_Index": annot.get("retention_index", np.nan)
                        }

                # Assign this best annotation to all features in the compound
                for feature_id in all_features_in_compound:
                    if feature_id in ft_dict:
                        # Combine base info (M0, quant) + best annotation
                        feature_data = base_feature_annotations.get(feature_id, {}).copy()
                        feature_data.update(best_annot_data)
                        feature_final_annotations[feature_id] = feature_data

        # --- 4. Final DataFrame Creation and Export ---
        
        # Get base columns and sample columns from original feature table
        base_cols = features.columns.tolist()
        # Assuming first 11 columns are metadata
        meta_cols = base_cols[:11] + ['RI']
        sample_cols = [col for col in base_cols[11:] if col not in meta_cols]
        
        new_annot_cols = [
            "Compound_ID", "Isotopologue", "Fragment",
            "Annotation_Name", "Annotation_Source", "Annotation_Score",
            "Annotation_Matched_Peaks", "Annotation_Library_Index"
        ]

        if mode == 'all':
            df_out = pd.DataFrame(all_annotation_records)
            # Define column order
            export_order = meta_cols + new_annot_cols + sample_cols
            # Reorder and filter columns
            df_out = df_out.reindex(columns=[col for col in export_order if col in df_out.columns])
            output_filename = table.replace(".tsv", "_all_annotations.tsv")

        elif mode == 'best':
            # Update the original feature data (in ft_dict) with annotations
            for feature_id, annotations in feature_final_annotations.items():
                if feature_id in ft_dict:
                    ft_dict[feature_id].update(annotations)
            
            df_out = pd.DataFrame(list(ft_dict.values()))
            # Define column order
            export_order = meta_cols + new_annot_cols + sample_cols
            # Reorder and filter columns
            df_out = df_out.reindex(columns=[col for col in export_order if col in df_out.columns])
            output_filename = table.replace(".tsv", "_annotated.tsv")

        else:
            print(f"Error: Invalid mode '{mode}'. Use 'best' or 'all'.")
            return

        # Export the final file
        df_out.to_csv(output_filename, sep="\t", index=False)
        print(f"Successfully exported annotations to: {output_filename}")

    def annotate(self):
        '''
        Annotate features via JMS (jms.dbStructures) and khipu.
        The pre-annotation step is khipu based emprical compound construction, followed by 
        three steps of annotation:

        1) Search known compound database, via neutral mass inferred by khipu

        2) Search singletons for a formula match

        3) Encapsulate remaining features in empCpd format, so that all are exported consistently.

        Export Feature_annotation as tsv, Annotated_empricalCompounds in both JSON and pickle.
        Reference databases can be pre-loaded. 
        Measured m/z values are calibrated to database based values (db_mass_calibrate).

        Note
        ----
        This produces default annotation with asari, 
        but one can redo annotation on the features afterwards, using a method of choice.
        With JMS/khipu, one can also pass custom adduct/isotopes to EED.adduct_patterns etc. 
        See ExperimentalEcpdDatabase.get_isotope_adduct_patterns().
        '''
        self.load_annotation_db()
        self.db_mass_calibrate()

        # asari uses seconds for rt
        EED = ExperimentalEcpdDatabase(mode=self.mode, 
                                       mz_tolerance_ppm=self.mz_tolerance_ppm, rt_tolerance=2)
        # passing patterns from .default_parameters
        if self.mode == 'pos':
            EED.adduct_patterns = adduct_search_patterns
        else:
            EED.adduct_patterns = adduct_search_patterns_neg
        EED.isotope_search_patterns = isotope_search_patterns
        EED.extended_adducts = extended_adducts

        EED.build_from_list_peaks(self.CMAP.FeatureList)
        # It takes three steps to take care of all features. First khipu organized empCpds
        EED.extend_empCpd_annotation(self.KCD)
        # Second, singletons that get a formula match in KCD
        EED.annotate_singletons(self.KCD)       
        # Third, the remaining features unmatched to anything (orphans). Exported for potential downstream work.
        EED.dict_empCpds = self.append_orphans_to_epmCpds(EED.dict_empCpds)

        self.export_peak_annotation(EED.dict_empCpds, self.KCD, 'Feature_annotation')

        if self.sample_registry:                    # check because some subcommands may not have sample_registry
            self.select_unique_compound_features(EED.dict_empCpds)
        
        # export JSON
        outfile = os.path.join(self.parameters['outdir'], 'Annotated_empiricalCompounds.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(EED.dict_empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)

        # export pickle
        outfile = os.path.join(self.parameters['outdir'], 'export', 'epd.pickle')
        with open(outfile, 'wb') as f:
            pickle.dump(EED.dict_empCpds, f, pickle.HIGHEST_PROTOCOL)

    def generate_qc_plot_pdf(self, outfile="qc_plot.pdf"):
        '''
        Generates a PDF figure of a combined plot feature quality metrics.
        Used only when --anno True (default).
        Skip if matplotlib is missing.
        '''
        try:
            from .qc import asari_qc_plot
            # using prefilter CMAP.FeatureTable
            asari_qc_plot(self.CMAP.FeatureTable,
                          outfile=os.path.join(self.parameters['outdir'], 'export', outfile))

        except ImportError:
            print("[QC plot] cannot import matplotlib, skipping.")

    def export_CMAP_pickle(self):
        '''
        Export main CMAP data and MassGrid to pickle, which can be used for visual data exploration.

        Included in exported pickle:{
            '_number_of_samples_': self.CMAP._number_of_samples_,
            'rt_length': self.CMAP.rt_length,
            'rt_reference_landmarks': [p['apex'] 
                                       for p in self.CMAP.good_reference_landmark_peaks],
            'rt_records': [sample.get_rt_calibration_records()
                                        for sample in self.all_samples
                                        ],
            'dict_scan_rtime': self.CMAP.dict_scan_rtime,
            'list_mass_tracks': self.CMAP.composite_mass_tracks,
            'MassGrid': dict(self.CMAP.MassGrid),}

        rt_records includes for each sample: {
            'sample_id': self.sample_id,
            'name': self.name,
            'rt_landmarks': self.rt_landmarks,
            'reverse_rt_cal_dict': self.reverse_rt_cal_dict,
        }

        Note
        ----
            RT calibration is exported to include sample.reverse_rt_cal_dict, 
            i.e. {key=reference scan number, value=sample specific scan number}.
            May add data throttle in the future. The file cmap.pickle can get big.
        '''
        _export = {
            '_number_of_samples_': self.CMAP._number_of_samples_,
            'rt_length': self.CMAP.rt_length,
            'rt_reference_landmarks': [p['apex'] 
                                       for p in self.CMAP.good_reference_landmark_peaks],
            'rt_records': [sample.get_rt_calibration_records()
                                        for sample in self.all_samples
                                        ],
            'dict_scan_rtime': self.CMAP.dict_scan_rtime,
            'list_mass_tracks': self.CMAP.composite_mass_tracks,
            'MassGrid': dict(self.CMAP.MassGrid),
        }
        outfile = os.path.join(self.parameters['outdir'], 'export', 'cmap.pickle')
        with open(outfile, 'wb') as f:
            pickle.dump(_export, f, pickle.HIGHEST_PROTOCOL)

    def load_annotation_db(self, src='hmdb4'):
        '''
        Load database of known compound using jms.dbStructures.knownCompoundDatabase.
        The compound tree is precomputed indexing.
        The `src` parameter is not used now, but placeholder to add more options later.

        Parameters
        ----------
        src: str, optional, default: hmdb4
            not used but can, in the future, dictate which database is used to generate annotations
        '''
        self.KCD = knownCompoundDatabase()
        self.KCD.mass_indexed_compounds = pickle.load( 
            pkg_resources.open_binary(db, 'mass_indexed_compounds.pickle') )
        self.KCD.emp_cpds_trees = pickle.load( 
            pkg_resources.open_binary(db, 'emp_cpds_trees.pickle') )

    def db_mass_calibrate(self, required_calibrate_threshold=0.000002):
        '''
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
        '''
        mz_landmarks = [self.CMAP.MassGrid['mz'][p['ref_id_num']] 
                        for p in self.CMAP.good_reference_landmark_peaks]
        mass_accuracy_ratio = self.KCD.evaluate_mass_accuracy_ratio(
            mz_landmarks, mode=self.mode, mz_tolerance_ppm=10)
        if mass_accuracy_ratio:
            if abs(mass_accuracy_ratio) > required_calibrate_threshold:
                print("Mass shift is greater than %2.1f ppm. Correction applied." 
                      %(required_calibrate_threshold*1000000))
                _correction = mass_accuracy_ratio + 1
                for F in self.CMAP.FeatureList:
                    F['mz'] = F['mz'] / _correction
                    F['mz_corrected_by_division'] = _correction
        else:
            print("Mass accuracy check is skipped, too few mz_landmarks (%d) matched." 
                  %len(mz_landmarks))

    def append_orphans_to_epmCpds(self, dict_empCpds):
        '''
        This is the third step of feature annotation in self.annotate,
        to encapsulate features without annotation in empCpd format.
        Input via dict_empCpds and returns updated dict_empCpds.
        See also: annotate

        Parameters
        ----------
        dict_empCpds: dict
            a dictionary of empirical compounds in empCpd format

        '''
        all_feature_ids = []
        for _, V in dict_empCpds.items():
            all_feature_ids += [peak['id_number'] for peak in V['MS1_pseudo_Spectra']]

        orphans = [peak for peak in self.CMAP.FeatureList if peak['id_number'] not in all_feature_ids]
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
        '''
        Export feature annotation to tab delimited tsv file, where interim_id is empCpd id.
        
        Parameters
        ----------
        dict_empCpds : dict
            dictionary of empirical compounds, using interim_id as key, as seen in JMS.
        KCD : KnownCompoundDatabase instance
            the known compound database that was used in annotating the empirical compounds.
        export_file_name_prefix : str
            to used in output file name.
        '''
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
        '''
        Get unique feature by highest composite peak area per empirical compound. 
        One may consider alternatives to select the peak representing an empirical compound,
        e.g. by SNR or M+H, M-H ions. This is can be done separately on the exported files.

        Parameters
        ----------
        dict_empCpds : dict
            dictionary of empirical compounds, using interim_id as key, as seen in JMS.
        '''
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

    def export_feature_tables(self, _snr=2, _peak_shape=0.7, _cSelectivity=0.7):
        '''
        To export multiple features tables. 
        Filtering parameters (_snr, _peak_shape, _cSelectivity) only 
        apply to preferred table and unique_compound_table.
        Full table is filtered by initial peak detection parameters, 
        which contains lower values of snr and gaussian_shape.

        Parameters
        ----------
        _snr: float, optional, default: 2
            signal to noise ratio, peaks must have SNR above this value to be a preferred feature
        _peak_shape: float, optional, default: 0.7
            the goodness fitting for a peak must be above this value to be a preferred feature
        _cSelectivity: float, optional, default: 0.7
            the cSelectivity of a peak must be above this value to be a preferred feature

        Outputs
        -------
        Multiple features tables under output directory:
        1. preferred table under `outdir`, after quality filtering 
           by SNR, peak shape and chromatographic selectivity.
        2. full table under `outdir/export/`
        3. unique compound table under `outdir/export/`
        4. dependent on `target` extract option, a targeted_extraction table under `outdir`. 
        '''
        if self.parameters['drop_unaligned_samples']:
            good_samples = [sample.name for sample in self.all_samples if sample.is_rt_aligned]
        else:
            good_samples = [sample.name for sample in self.all_samples]
        filtered_FeatureTable = self.CMAP.FeatureTable[good_samples]          # this fixes order of samples       
        self.number_of_samples = number_of_samples = len(good_samples)
        self.dropped_sample_names = [sample.name for sample in self.all_samples 
                                     if sample.name not in good_samples]     # samples dropped due to alignment or quality
        # non zero counts
        count = filtered_FeatureTable[filtered_FeatureTable>1].count(axis='columns')
        self.CMAP.FeatureTable['detection_counts'] = count

        use_cols = [ 'id_number', 'mz', 'rtime', 'rtime_left_base', 'rtime_right_base', 'parent_masstrack_id', 
                    'peak_area', 'cSelectivity', 'goodness_fitting', 'snr', 'detection_counts' ] + good_samples
        filtered_FeatureTable = self.CMAP.FeatureTable[use_cols]
        filtered_FeatureTable['mz'] = filtered_FeatureTable['mz'].round(4)
        filtered_FeatureTable['rtime'] = filtered_FeatureTable['rtime'].round(2)
        filtered_FeatureTable['rtime_left_base'] = filtered_FeatureTable['rtime_left_base'].round(2)
        filtered_FeatureTable['rtime_right_base'] = filtered_FeatureTable['rtime_right_base'].round(2)
        filtered_FeatureTable['cSelectivity'] = filtered_FeatureTable['cSelectivity'].round(2)
        filtered_FeatureTable['goodness_fitting'] = filtered_FeatureTable['goodness_fitting'].round(2)

        outfile = os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table'])
        filtered_FeatureTable.to_csv(outfile, index=False, sep="\t")
        print("\nFeature table (%d x %d) was written to %s." %(
                                filtered_FeatureTable.shape[0], number_of_samples, outfile))

        # extract targeted m/z features
        if 'target' in self.parameters and self.parameters['target']:  
            matched_list, _, target_unmapped = all_mass_paired_mapping(
                filtered_FeatureTable['mz'].to_list(), self.parameters['target'], self.parameters['mz_tolerance_ppm']
            )
            print("\nIn targeted extraction, %d target mz values are not found in this dataset: " %len(target_unmapped))
            print('    ', [self.parameters['target'][ii] for ii in target_unmapped])
            matched_targets = [self.parameters['target'][ii[1]] for ii in matched_list]
            targeted_table = filtered_FeatureTable.iloc[[x[0] for x in matched_list], :]
            targeted_table.insert(0, "query_target", matched_targets)
            outfile = os.path.join(self.parameters['outdir'], 'targeted_extraction__'+self.parameters['output_feature_table'])
            targeted_table.to_csv(outfile, index=False, sep="\t")
            print("Targeted extraction Feature table (%d x %d) was written to %s.\n" %(
                                targeted_table.shape[0], number_of_samples, outfile))

        outfile = os.path.join(self.parameters['outdir'], 'preferred_'+self.parameters['output_feature_table'])
        # Some features can have all 0s, filtered here
        filtered_FeatureTable = filtered_FeatureTable[  filtered_FeatureTable['detection_counts'] > 0 ]
        filtered_FeatureTable = filtered_FeatureTable[  filtered_FeatureTable['snr']>_snr][
                                                        filtered_FeatureTable['goodness_fitting']>_peak_shape][
                                                        filtered_FeatureTable['cSelectivity']>_cSelectivity ]
        filtered_FeatureTable.to_csv(outfile, index=False, sep="\t")
        print("\nFiltered Feature table (%d x %d) was written to %s.\n" %(
                                filtered_FeatureTable.shape[0], number_of_samples, outfile))
        
        # todo - this should be made more robust in the future
        if self.parameters['anno'] and self.parameters['workflow'] != "GC":
            # in self.selected_unique_features: (empCpd id, neutral_formula, ion_relation)
            sel = [ii for ii in filtered_FeatureTable.index if filtered_FeatureTable['id_number'][ii] in self.selected_unique_features.keys()]
            unique_compound_table = filtered_FeatureTable.loc[sel, :]
            unique_compound_table.insert(3, "empCpd", [self.selected_unique_features[ii][0] for ii in unique_compound_table['id_number']])
            unique_compound_table.insert(4, "neutral_formula", [self.selected_unique_features[ii][1] for ii in unique_compound_table['id_number']])
            unique_compound_table.insert(5, "ion_relation", [self.selected_unique_features[ii][2] for ii in unique_compound_table['id_number']])
            
            outfile = os.path.join(self.parameters['outdir'], 'export', 'unique_compound__'+self.parameters['output_feature_table'])
            unique_compound_table.to_csv(outfile, index=False, sep="\t")
            print("Unique compound table (%d x %d) was written to %s.\n" %(
                                unique_compound_table.shape[0], number_of_samples, outfile))

        
    def export_log(self):
        '''
        Export project parameters to project.json,
        which is also used by asari viz.
        '''
        self.parameters['number_of_samples'] = self.number_of_samples
        self.parameters['number_scans'] = self.number_scans
        # Record of dropped samples. Right place in log, wrong name in parameters
        self.parameters['dropped_samples'] = self.dropped_sample_names
        outfile = os.path.join(self.parameters['outdir'], 'project.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(self.parameters, f, cls=NpEncoder, ensure_ascii=False, indent=2)

    def export_readme(self):
        '''
        Export a REAME.txt file as simple instruction to end users.
        '''
        outfile = os.path.join(self.parameters['outdir'], 'README.txt')
        with open(outfile, 'w', encoding='utf-8') as f:
            f.write(readme_doc_str)

