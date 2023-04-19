import pickle
from .mass_functions import flatten_tuplelist
from .workflow import *

class SimpleSample:
    '''
    Lightweight class of an experimental sample to facilitate workflow.
    The primary use of this class is managing data stroage and retreival at sample level.

    Depending on database_mode, sample list_mass_tracks are stored in memory, or on disk, or in MongoDB.
    Function to get mass tracks from a mzML file is in workflow.process_project and batch_EIC_from_samples_.
    Peaks and empCpds are determined in constructors.CompositeMap.
    '''
    def __init__(self, registry={}, experiment=None, database_mode='ondisk', mode='pos', is_reference=False):
        '''
        Build a lightweight sample class.

        Parameters
        ----------
        registry : sample_registry, a dictionary like {'sample_id': ii, 'input_file': file}, 
            generated by workflow.register_samples.
        experiment : mostly required pointer to ext_Experiment instance, in order to get parameters.
        database_mode : 'ondisk' or 'memory' (run in memory, only small studies).
        mode: ionization mode, 'pos' or 'neg'. This should be consistent with experiment and registry if given.

        Note:
            m/z calibration is performed in m/z alignment for "small studies", 
            where mass tracks are assembled to MassGrid via landmark peaks with m/z calibration per sample.
            For larger studies, m/z alignment is done via NN clustering, where m/z accuracy 
            is not checked during MassGrid construction. But it is checked during DB annotation.
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
        # placeholder
        self.mz_calibration_function = None
                                   

    def get_masstracks_and_anchors(self):
        '''
        Retrieve list_mass_tracks for this sample if not alrady in memory.

        Returns
        ------- 
        A list of all mass tracks in this sample.

        Note:
            Mass tracks are the bulk of data per sample, stored dependent on database_mode.
            list_mass_tracks is accessed twice in this version of asari:
            1) RT calibration and building composite map
            2) extraction of peak areas for features
        '''
        if self.list_mass_tracks:     # important, this is used to check if in memory
            return self.list_mass_tracks
        else:
            sample_data = self._get_sample_data()
            list_mass_tracks = sample_data['list_mass_tracks']
            return list_mass_tracks

    def _get_sample_data(self):
        '''
        Wrapper of _retrieve_from_disk function.
        Old function kept to leave room for additional logics.

        Note:
            Potential use of database_mode, e.g.
            if self.database_mode == 'ondisk': 
                return self._retrieve_from_disk()
            elif: self.database_mode == 'firebase': 
                return self.retrieve_from_db()
        '''
        return self._retrieve_from_disk()

    def _retrieve_from_disk(self):
        '''
        Retrieve sample data from local pickle file.
        '''
        with open(self.data_location, 'rb') as f:
            sample_data = pickle.load(f)
        return sample_data

    def push_to_db(self, cursor):
        '''
        Placeholder.
        '''
        pass

    def retrieve_from_db(self, cursor):
        '''
        Placeholder.
        '''
        pass
