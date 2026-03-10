'''
Functions for subcommand `annotate`
'''
import json

from jms.io import read_table_to_peaks

from .experiment import (ext_Experiment,
                         ExperimentalEcpdDatabase, 
                         adduct_search_patterns, 
                         adduct_search_patterns_neg, 
                         isotope_search_patterns)
from .default_parameters import extended_adducts
from .utils import NpEncoder


class GC_Annotation:
    '''
    Do annotation here instead of during processing
    
    '''



    def populate_RI_lookup(self, sample_map):
        RI_maps = {}
        RI_models = {}
        reverse_RI_models = {}
        RI_list = pd.read_csv(self.parameters['retention_index_standards'])
        for reference_id in tqdm.tqdm(list(dict.fromkeys(list(sample_map.values())))):
            print(reference_id)
            RI_maps[reference_id] = {}
            reference_instance = SimpleSample(self.sample_registry[reference_id], experiment=self)
            prev_index, next_index = None, None
            prev_rt, next_rt = None, None
            RTs, indexes, scan_nos = [], [], []
            for rt, scan_no in zip(reference_instance.list_retention_time, reference_instance.list_scan_numbers):
                RTs.append(rt)
                scan_nos.append(scan_no)
                print(rt, scan_no)
                for index, index_rt in zip(RI_list['Index'], RI_list[reference_instance.name]):
                    index, index_rt = int(index), float(index_rt)
                    print("\t", index, index_rt)
                    if rt > index_rt:
                        prev_index, prev_rt = index, index_rt
                    elif rt <= index_rt:
                        _, next_rt = index, index_rt
                        break
                if next_rt is None:
                    next_rt = max(reference_instance.list_retention_time) * 1.1
                    _ = max(RI_list['Index']) + 1
                RI_value = 100 * (prev_index + ((rt - prev_rt)/(next_rt - rt)))
                indexes.append(RI_value)
                RI_maps[reference_id][rt] = RI_value
            model = lowess(indexes, RTs)
            model2 = lowess(indexes, scan_nos)
            newx, newy = list(zip(*model))
            interf = interpolate.interp1d(newx, newy, fill_value="extrapolate", bounds_error=False)
            RI_models[reference_id] = interf
            newx, newy = list(zip(*model2))
            reverse_RI_models[reference_id] = interpolate.interp1d(newy, newx, fill_value="extrapolate", bounds_error=False)
            
        self.RI_models = RI_models
        self.reverse_RI_models = reverse_RI_models

    def convert_to_RI(self, sample_map):
        if not self.RI_map:
            self.populate_RI_lookup(sample_map)
        for k, v in sample_map.items():
            sam = self.sample_registry[k]
            sam['list_retention_index'] = self.RI_models[v](sam['list_retention_time'])


    def annotate_GC(self):
        pref_ft = os.path.join(self.parameters['outdir'], 'preferred_'+self.parameters['output_feature_table'])
        full_ft = os.path.join(self.parameters['outdir'], 'export', 'full_'+self.parameters['output_feature_table'])
        EI_MS_Library.annotate_gc_feature_table_with_library(pref_ft, self.parameters['GC_Database'])
        EI_MS_Library.annotate_gc_feature_table_with_library(full_ft, self.parameters['GC_Database'])






# todo - does this still work? complaining about missing append_orphans
def annotate_user_featuretable(infile, parameters, rtime_tolerance=2):
    '''
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
    '''
    parameters['outdir'] = ''
    mode = parameters['mode']
    list_peaks = read_table_to_peaks(infile, 
                                has_header=True, mz_col=0, rtime_col=1, feature_id=None ,
                                )
    # print("Read %d features." %len(list_peaks))
    EE = ext_Experiment({}, parameters)
    EE.load_annotation_db()

    EED = ExperimentalEcpdDatabase(mode=mode, 
                                   mz_tolerance_ppm=parameters['mz_tolerance_ppm'],
                                   rt_tolerance=rtime_tolerance)
    
    # passing patterns from .default_parameters
    if mode == 'pos':
        EED.adduct_patterns = adduct_search_patterns
    else:
        EED.adduct_patterns = adduct_search_patterns_neg
    EED.isotope_search_patterns = isotope_search_patterns
    EED.extended_adducts = extended_adducts
    
    EED.build_from_list_peaks(list_peaks)
    EED.extend_empCpd_annotation(EE.KCD)
    EED.annotate_singletons(EE.KCD)
    # EED.dict_empCpds misses some features 
    EED.dict_empCpds = EED.append_orphans_to_epmCpds(EED.dict_empCpds)

    EE.export_peak_annotation(EED.dict_empCpds, EE.KCD, 'Feature_annotation')
    # also exporting JSON
    with open('Annotated_empricalCompounds.json', 'w', encoding='utf-8') as f:
        json.dump(EED.dict_empCpds, f, cls=NpEncoder, ensure_ascii=False, indent=2)
