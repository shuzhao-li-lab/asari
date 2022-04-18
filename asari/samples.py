import pickle

from .mass_functions import flatten_tuplelist
from .workflow import *


class SimpleSample:
    '''
    Lightweight class of an experimental sample to facilitate workflow.
    Depending on database_mode, sample list_mass_tracks are stored in memory, or on disk, or in MongoDB.
    Function to get mass tracks from a mzML file is in workflow.process_project and batch_EIC_from_samples_.
    Peaks and empCpds are determined in constructors.CompositeMap.
    '''
    def __init__(self, registry={}, experiment=None, database_mode='ondisk', mode='pos'):
        '''
        experiment: mostly required pointer to ext_Experiment instance, in order to get parameters.
        database_mode: 'ondisk', 'mongo', or 'memory' (run in memory, only small studies).
                This is not necessary to be same as experiment wide setting, to have flexibility in handling individual samples.
                E.g. two samples can be one in memory and the other on disk.
        mode: ionization mode, 'pos' or 'neg'

        '''
        self.registry = registry
        self.input_file = registry['input_file']
        self.name = registry['name']
        self.sample_id = registry['sample_id']
        self.data_location = registry['data_location']
        self.sample_data = registry['sample_data']

        self.experiment = experiment
        self.mode = mode
        self.database_mode = database_mode 

        self.rt_numbers = []                                # list of scans, starting from 0
        self.list_retention_time = []                       # full RT time points in sample

        # These calibration functions are critical, though mz_calibration is not required if m/z accuracy is good.
        self.mz_calibration_ratio = None
        self.rt_cal_dict = None                             # index mapping to the reference sample
        self.reverse_rt_cal_dict = None
        
        # lists to store data; empty if not stored in memory and need retrieval from disk or DB
        self.list_mass_tracks = []                          # index number = id_number, in ascending order of m/z
        self.anchor_mz_pairs = []
        self._mz_landmarks_ = []
        self._number_anchor_mz_pairs_ = 0

    def get_masstracks_and_anchors(self):
        '''
        Retrieve stored data, including list_mass_tracks, anchors and supporting stats. 
        '''
        sample_data = self._get_sample_data()
        list_mass_tracks = sample_data['list_mass_tracks']
        if not self.anchor_mz_pairs:
            self.anchor_mz_pairs = sample_data['anchor_mz_pairs']
            self._number_anchor_mz_pairs_ = sample_data['number_anchor_mz_pairs']
            self.rt_numbers = sample_data['list_scan_numbers']
            self.list_retention_time = sample_data['list_retention_time']
            self._mz_landmarks_ = flatten_tuplelist(self.anchor_mz_pairs)
        return list_mass_tracks

    def _get_sample_data(self):
        '''
        Retrieve mass tracks, which are the bulk of data per sample, depending on where they are stored.
        '''
        if self.database_mode == 'memory':
            return self.sample_data
        elif self.database_mode == 'ondisk': 
            return self._retrieve_from_disk()
        else:                         # requiring self.experiment and db connection
            return self.retrieve_from_db(
                # to implement
            )

    def _retrieve_from_disk(self):
        with open(self.data_location, 'rb') as f:
            sample_data = pickle.load(f)
        return sample_data


    def push_to_db(self, cursor):
        '''
        Each sample -> TABLE mass tracks (pickle objects)
        '''

        pass


    def retrieve_from_db(self, cursor):
        pass

