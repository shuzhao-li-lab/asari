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
    def __init__(self, registry={}, experiment=None, database_mode='ondisk', mode='pos', is_reference=False):
        '''
        experiment: mostly required pointer to ext_Experiment instance, in order to get parameters.
        database_mode: 'ondisk', 'mongo', or 'memory' (run in memory, only small studies).
        mode: ionization mode, 'pos' or 'neg'. This should be consistent with experiment and registry if given.
        '''
        self.experiment = experiment
        self.mode = mode
        self.database_mode = database_mode 
        self.is_reference = is_reference 

        self.input_file = registry['input_file']
        self.name = registry['name']
        self.sample_id = registry['sample_id']
        self.data_location = registry['data_location']
        self.track_mzs = registry['track_mzs']
        self.max_scan_number = registry['max_scan_number']
        self.anchor_mz_pairs = registry['anchor_mz_pairs']
        self.rt_numbers = registry['list_scan_numbers']
        self.list_retention_time = registry['list_retention_time']

        if self.database_mode == 'memory':
            self.list_mass_tracks = registry['sample_data']['list_mass_tracks']
        else:
            self.list_mass_tracks = []
            
        self._mz_landmarks_ = flatten_tuplelist(self.anchor_mz_pairs)

        # These are critical RT calibration functions, index mapping with the reference sample
        self.rt_cal_dict = None
        self.reverse_rt_cal_dict = None

        # function mz_calibration is not used yet, but will be.
        # For "small studies", mass tracks are assembled to MassGrid via landmark peaks with m/z calibration
        # For larger studies, m/z accuracy is not checked during MassGrid construction. Thus,
        # m/z calibration can be supplied to this function by checking spike-ins etc.
        self.mz_calibration_function = None
                                   

    def get_masstracks_and_anchors(self):
        '''
        Retrieve stored data, including list_mass_tracks, anchors and supporting stats. 
        Except the reference sample, mass tracks are not kept in memory.

        # index number = id_number, in ascending order of m/z

        list_mass_tracks is accessed twice in this version of asari:
        1) RT calibration and building composite map
        2) extraction of peak areas for features
        '''
        if self.list_mass_tracks:                           # important, this is used to check if in memory
            return self.list_mass_tracks
        else:
            sample_data = self._get_sample_data()
            list_mass_tracks = sample_data['list_mass_tracks']
            return list_mass_tracks

    def _get_sample_data(self):
        '''
        Retrieve mass tracks, which are the bulk of data per sample, depending on where they are stored.
        if self.database_mode == 'memory':
            return self.sample_data
        elif self.database_mode == 'ondisk': 
            return self._retrieve_from_disk()
        else:                         # requiring self.experiment and db connection
            return self.retrieve_from_db(
                # to implement)
        '''
        return self._retrieve_from_disk()


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

