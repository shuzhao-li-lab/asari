'''
ext_Experiment is the container for whole project data.

Heavy lifting is in constructors.CompositeMap, 
which contains MassGrid for correspondence, and FeatureList from feature/peak detection.

'''
import os
import random
import json
import pickle

from metDataModel.core import Experiment
from mass2chem.search import *
from mass2chem.epdsConstructor import epdsConstructor
# jms-metabolite-services
from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from .samples import SimpleSample
from .constructors import CompositeMap
from .sql import *

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from . import db


# General data processing steps are in this class
class ext_Experiment(Experiment):
    '''
    Extend metDataModel.core.Experiment with preprocessing methods.
    This encapsulates a set of LC-MS files using the same method to be processed together.
    '''
    def __init2__(self, list_input_files, dict_meta_data, parameters, output_dir):
        '''
        This is the overall container for all data in an experiment/project.
        Samples are sorted by name and assigned unique IDs.

        Input
        -----
        list_input_files: list of inputfiles, including directory path, to read
        dict_meta_data: description of sample types for each file, e.g. 'QC', 'pooled', 'sample'.
        parameters: including 'ionization_mode', 'min_intensity_threshold', 'min_timepoints'. See main.py.

        '''
        self.list_input_files = sorted(list_input_files)        # ordered by name
        self.output_dir = output_dir
        self.number_of_samples = len(list_input_files)

        self.files_meta_data = dict_meta_data
        
        self.parameters = parameters
        # self.max_rtime = parameters['max_rtime'] # to update from each sample
        self.number_scans = 0                                   # max scan number, to update when samples are processed
        self.mode = parameters['mode']

        self.initiation_samples = self.__choose_initiation_samples__(N=6)

        # SAMPLE_REGISTRY
        self.all_samples = []               # list of Sample instances
        self.samples_nonreference = []
        self.samples_by_id = {}             # sample ID: Sample instance
        self.samples_by_name = {}           # input file name: Sample instance

        
    def process_all(self):
        '''
        initiation_Samples are used to select one most representative sample to seed MassGrid and RT alignment.
        If refDB is used, it's better to be used after all samples are processed, 
        because the m/z values of samples are closer to other samples than refDB.
        
        if refDB:
            self.CMAP.align_to_refdb(refDB)

        '''
        # start SQLite database
        # self.cursor = connect_sqlite_db(self.parameters['project_name'])
        # if self.number_of_samples > NNN:
        #      finish in memory first
        # else:
        #      start DB after init
        
        self.process_all_without_export()

        # 
        self.CMAP.MassGrid.to_csv(
            os.path.join(self.parameters['outdir'], self.parameters['mass_grid_mapping']) )
        self.annotate()
        self.export_feature_table()


    def process_all_without_export(self):
        self.CMAP = CompositeMap(self)
        self.CMAP.construct_mass_grid( self.process_initiation_samples() )
        self.CMAP.align_retention_time()
        # some samples could fail alignment; can be processed and aligned at the end
        self.CMAP.global_peak_detection()


    def process_initiation_samples(self):
        return [self.process_single_sample(f) for f in self.initiation_samples]


    def process_single_sample(self, input_file):
        '''
        Some parameters can be automatically determined here.

        To add DB function HERE ??
        '''
        mz_tolerance_ppm = self.parameters['mz_tolerance']
        min_intensity = self.parameters['min_intensity_threshold']
        min_timepoints = self.parameters['min_timepoints']
        try:
            SM = SimpleSample(self, self.mode, input_file)
            SM.process( mz_tolerance_ppm, min_intensity, min_timepoints)
            # sample id, assigned by index in self.list_input_files.
            # DB commit
            return SM
        except IndexError:
            print("Input error in sample %s, dropped from processing." %input_file)
            return None
        

    def __choose_initiation_samples__(self, N=3):
        '''
        N initial samples are chosen to be analyzed first.
        One best sample among them is chosen as the reference, especially for retention time alignment.
        '''
        if self.parameters['initiation_samples']:
            return self.parameters['initiation_samples']
        else:
            if self.number_of_samples < N+1:
                return self.list_input_files
            else:
                return random.sample(self.list_input_files, N)


    def annotate(self):
        '''
        
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(list_empCpds, f, ensure_ascii=False, indent=2)

        print("\nEmpirical compound annotaion (%d) was written to %s." %(len(list_empCpds), outfile))

        '''
        outfile = os.path.join(self.parameters['outdir'], 
                                                            'annotated_empricalCompounds_')
        ECCON = epdsConstructor(self.CMAP.FeatureList, mode=self.mode)
        list_empCpds = ECCON.peaks_to_epds()
        #list_empCpds = self._reformat_epds_(list_empCpds, self.CMAP.FeatureList)

        # use JMS
        KCD = knownCompoundDatabase()
        # will move pickle files to db/
        KCD.mass_indexed_compounds = pickle.load( pkg_resources.open_binary(db, 'mass_indexed_compounds.pickle') )
            # open('mass_indexed_compounds.pickle', 'rb') )
        KCD.emp_cpds_trees = pickle.load( pkg_resources.open_binary(db, 'emp_cpds_trees.pickle') )
            # open('emp_cpds_trees.pickle', 'rb') )

        EED = ExperimentalEcpdDatabase(mode=self.mode)
        EED.list_peaks = self.CMAP.FeatureList
        EED.dict_empCpds = EED.index_reformat_epds(list_empCpds, self.CMAP.FeatureList)
        EED.index_empCpds()
        search_result = EED.annotate_against_KCD(KCD)
        EED.export_annotations(search_result, KCD, outfile)


    def _reformat_epds_(self, list_empCpds, FeatureList):
        fDict = {}
        for F in FeatureList:
            fDict[F['id_number']] = F
        new = []
        for E in list_empCpds:
            features = []
            for peak in E['list_peaks']:
                features.append(
                    {'feature_id': peak[0], 
                    'mz': fDict[peak[0]]['mz'], 
                    'rtime': fDict[peak[0]]['apex'], 
                    'charged_formula': '', 
                    'ion_relation': peak[1]}
                )
            new.append(
                {
                'interim_id': E['id'], 
                'neutral_formula_mass': None,
                'neutral_formula': None,
                'Database_referred': [],
                'identity': [],
                'MS1_pseudo_Spectra': features,
                'MS2_Spectra': [],
                }
            )
        return new


    def export_feature_table(self, full=True, outfile='cmap_feature_table.csv'):
        '''
        Will need real RT time;
        Selectivity in m/z, RT and overall
        
        '''
        outfile = os.path.join(self.parameters['outdir'], self.parameters['output_feature_table'])
        if full:
            self.CMAP.FeatureTable.to_csv(outfile)
        else:
            # select columns to export
            pass

        print("\n\nFeature table (%d) was written to %s.\n\n" %(self.CMAP.FeatureTable.shape[0], outfile))







    def export_feature_table_old__(self, FeatureList, outfile='feature_table.tsv'):
        '''
        FeatureList: a list of namedTuples, i.e. Features; Output two files, one main, another low quality features.
        '''
        def __write__(FeatureList, outfile):
            s = '\t'.join(['feature_id', 'formula_mass', 'mz', 'rtime', 'rt_min', 'rt_max', 'number_peaks',
                                    'peak_quality_max', 'peak_quality_median', 'intensity_mean', 'selectivity_mz',
                                    ] + self.ordered_sample_names) + '\n'
            for F in FeatureList:
                s += '\t'.join(
                    [F.feature_id, F.mass_id, str(round(F.mz,4)), str(round(F.rtime,2)), str(round(F.rt_min,2)), str(round(F.rt_max,2)), str(F.number_peaks),
                    str(round(F.peak_quality_max,2)), str(round(F.peak_quality_median,2)), str(F.intensity_mean), str(round(F.selectivity_mz,2)),
                    ] + [str(int(x)) for x in F.intensities]
                    ) + '\n'
            with open( outfile, 'w') as O:
                O.write(s)

        high_quality_features, low_quality_features = [], []
        for F in FeatureList: 
            if F.peak_quality_max > 0.8 and F.perc_peaks > 15:
                high_quality_features.append(F)
            else:
                low_quality_features.append(F)

        high_quality_features.sort(key=lambda F: F.peak_quality_max, reverse=True)
        #print(FeatureList[99])
        __write__(high_quality_features, os.path.join(self.output_dir, outfile))
        __write__(low_quality_features, os.path.join(self.output_dir, 'low_quality_features_' + outfile))
        print("Feature tables were written under %s." %self.output_dir)
        print("The main feature table (%s) has %d samples and %d features.\n\n\n" %(
            self.parameters['output_filename'], len(self.ordered_sample_names), len(high_quality_features)))




    #---------------------------------------------------------------------------------------------------------------

    def __obsolete__process_all(self):
        '''
        This will shift to a DB design in next version.
        '''
        self.init_hot_db( self._get_ref_db_() )                         # initial processing of 3 samples to set up HOT_DB
        for f in self.list_input_files:                                 # run remaining samples
            if f not in self.initiation_samples:
                SM = Sample(self, self.mode, f)
                SM.process_step_1()
                SM.process_step_2(self.HOT_DB)
                if not self.parameters['cache_mass_traces']:
                    del(SM.dict_masstraces)
                self.samples.append(SM)
        
        self.calibrate_retention_time()                                 # samples may be marked to drop
        self.correspondence()
        self.annotate_final()
        self.export_feature_table(self.FeatureTable, self.parameters['output_filename'])
        
    def _get_ref_db_(self):
        '''
        Dispatch for ref DB.
        Earlier version used INIT_DFDB = DB_to_DF( extend_DB1(DB_1) ), which was moved to mass2chem.
        '''
        if self.mode == 'pos':
            dbfile = os.path.join(os.path.dirname(__file__), 'ref_db_v0.2.tsv')
        elif self.mode == 'neg':
            dbfile = os.path.join(os.path.dirname(__file__), 'neg_ref_db_v0.2.tsv')
        else:
            print("Ionization mode is either `pos` or `neg`.")
        return tsv2refDB(dbfile)


    def init_hot_db(self, DFDB):
        '''
        Use three samples to initiate a hot DB to house feature annotation specific to this Experiment, and speed up subsequent search.
        The HOT_DB will be used during sample processing, and have another update after correspondence and additional annotation.
        The HOT_DB will then be exported as Expt annotation.
        '''
        chosen_Samples, found_formula_masses = [], []
        for f in self.initiation_samples:
            SM = Sample(self, self.mode, f)
            SM.process_step_1()
            SM.process_step_2(DFDB)
            chosen_Samples.append(SM)
            found_formula_masses += list(SM.mzstr_2_formula_mass.values())

        self.samples += chosen_Samples
        # Experiment wide parameters
        self.__mass_stdev__ = np.median([SM.__mass_stdev__ for SM in chosen_Samples])       # ppm stdev, used for later searches
        # ver 1 HOT_DB will be a subset of INIT_DFDB

        found_formula_masses = set(found_formula_masses)
        if None in found_formula_masses:
            found_formula_masses.remove(None)
        self.HOT_DB = DFDB.loc[found_formula_masses]   

        print("\n[@.@] Anchoring with %d initial formula matches." %len(found_formula_masses))
        print("[@.@] Initial estimation done on\n" + '\n'.join(self.initiation_samples))
        #export_hot_db, without last col
        self.HOT_DB.iloc[:, :-1].to_csv(os.path.join(self.output_dir, '__intermediary__' + self.parameters['annotation_filename']), sep="\t")
        

    def annotate_final(self):
        '''
        More formula annotation of _M_ features using HMDB+PubChemLite;
        Group into empCpds via mass2chem.


        Still to do remaining features
        '''
        s = u'\t'.join(['feature_id', 'formula_mass', 'mz_dbrecord',	'intensity_mean', 'charged_formula', 'selectivity',	'neutral_formula_mass',
                                    'ion_relation', 'id_HMDB', 'name']) + '\n'
        for F in self.FeatureTable:
            if "_M_" == F.mass_id[:3]:
                s += u'\t'.join([F.feature_id, F.mass_id, str(round(F.mz,4)), str(F.intensity_mean)]) + '\n'
            else:
                [mz, charged_formula, selectivity, neutral_formula_mass, ion_relation] = [str(x) for x in list(self.HOT_DB.loc[F.mass_id])[:5]]
                name = massDict_hmdb.get(neutral_formula_mass, '')
                if name:
                    name = u'\t'.join( [';'.join(x) for x in name] ).encode('utf-8', 'ignore').decode('utf-8')
                s += u'\t'.join([F.feature_id, F.mass_id, mz, str(F.intensity_mean),
                                charged_formula, selectivity, neutral_formula_mass, ion_relation, name]) + '\n'
                
        with open(os.path.join(self.output_dir, self.parameters['annotation_filename']), 'w', encoding='utf-8') as O:
            O.write(s.encode('utf-8', 'ignore').decode('utf-8'))

        



