import pickle
import pandas as pd
from .LRU_cache import lru_cache
from .mass_functions import flatten_tuplelist
from .peaks import get_gaussian_peakarea_on_intensity_list, peak_area_sum, peak_area_auc


class SimpleSample:
    '''
    Lightweight class of an experimental sample to facilitate workflow.
    The primary use of this class is managing data stroage and retreival at sample level.

    Depending on database_mode, sample list_mass_tracks are stored in memory, or on disk, or in MongoDB.
    Function to get mass tracks from a mzML file is in workflow.process_project and batch_EIC_from_samples_.
    Peaks and empCpds are determined in constructors.CompositeMap.
    '''
    peak_area_mode_map = {
        "auc": peak_area_auc,
        "gauss": get_gaussian_peakarea_on_intensity_list,
        "sum": peak_area_sum
    }

    def __init__(self, registry, experiment=None, is_reference=False):
        '''
        Build a lightweight sample class.

        Parameters
        ----------
        registry : dict
            sample_registry, a dictionary like {'sample_id': ii, 'input_file': file}, 
            generated by workflow.register_samples.
        experiment : ext_Experiment instance
            mostly required pointer to ext_Experiment instance, in order to get parameters.
        database_mode : str, optional, default: 'ondisk
            'ondisk' or 'memory' (run in memory, only small studies).
        mode: str, optional, default: 'pos'
            ionization mode, 'pos' or 'neg'. This should be consistent with experiment and registry if given.

        Note
        ----
        m/z calibration is performed in m/z alignment for "small studies", 
        where mass tracks are assembled to MassGrid via landmark peaks with m/z calibration per sample.
        For larger studies, m/z alignment is done via NN clustering, where m/z accuracy 
        is not checked during MassGrid construction. But it is checked during DB annotation.
        '''
        self.experiment = experiment
        self.is_reference = is_reference 
        self.registry = registry
        self.rt_landmarks = []  # to populate at CMAP.calibrate_sample_RT

        self.__list_mass_tracks = None
        self.__mz_landmarks = None
        self.__retrieve_mode = None

        # These are critical RT calibration functions, index mapping with the reference sample
        self.rt_cal_dict = None
        self.reverse_rt_cal_dict = None
        # placeholder
        self.mz_calibration_function = None
        self.rt_mapping = {}

    @property
    def input_file(self):
        return self.registry['input_file']
        
    @property
    def name(self):
        return self.registry['name']
    
    @property
    def sample_id(self):
        return self.registry['sample_id']
    
    @property
    def data_location(self):
        return self.registry['data_location']
    
    @property
    def track_mzs(self):
        return self.registry['track_mzs']
    
    @property
    def max_scan_number(self):
        return self.registry['max_scan_number']
    
    @property
    def anchor_mz_pairs(self):
        return self.registry['anchor_mz_pairs']

    @property
    def rt_numbers(self):
        return self.registry['list_scan_numbers']
    
    @property
    def list_retention_time(self):
        return self.registry['list_retention_time']
    
    @property
    def mem_footprint(self):
        return self.registry['mem_footprint']
    
    @property
    def sample_mapping(self):
        if 'sample_mapping' in self.registry:
            return self.registry['sample_mapping']
        else:
            return None

    @property
    def database_mode(self):
        return self.experiment.database_mode

    @property
    def composite_of(self):
        if 'composite_of' in self.registry:
            return self.registry['composite_of']
        else:
            return None
    
    @property
    def reserved_memory(self):
        return self.experiment.parameters['reserve_memory'] * 1024**3

    @property
    def list_mass_tracks(self):
        mass_tracks = None
        if self.database_mode == "smart":
            mass_tracks = self.list_mass_tracks_smart()
        elif self.database_mode == "ondisk":
            mass_tracks = self._get_sample_data()['list_mass_tracks']
        elif self.database_mode == "memory":
            mass_tracks = self.registry['sample_data']['list_mass_tracks']
        if type(mass_tracks) is dict:
            return list(mass_tracks.values())   
        return mass_tracks

    @lru_cache(use_memory_up_to=8 * 1024**3)
    def list_mass_tracks_smart(self):
        mass_tracks = self._get_sample_data()['list_mass_tracks']
        if type(mass_tracks) is dict:
            return list(mass_tracks.values())
        return mass_tracks
    
    @property
    def _mz_landmarks_(self):
        if self.__mz_landmarks is None:
            self.__mz_landmarks = flatten_tuplelist(self.anchor_mz_pairs)
        return self.__mz_landmarks
    
    @property
    def mode(self):
        return self.experiment.mode
    
    @property
    def peak_area_mode(self):
        return self.experiment.peak_area_mode

    def calc_peak_area(self, mass_track_id, left_base, right_base):
        if self.composite_of is None:
            if pd.isna(mass_track_id):
                return {self.name: 0}
            else:
                return {self.name: self.peak_area_mode_map[self.peak_area_mode](self.list_mass_tracks[mass_track_id]['intensity'], left_base, right_base)}
        else:
            return {k: v for d in [x.calc_peak_area(mass_track_id, left_base, right_base) for x in self.composite_of] for k, v in d.items()}
        
    def calc_peak_areas(self, peaks):
        if self.composite_of is None:
            areas = []
            for p in peaks:
                if p[0] is None or pd.isna(p[0]):
                    areas.append(0)
                else:
                    mass_tracks = self.list_mass_tracks
                    areas.append(self.peak_area_mode_map[self.peak_area_mode](mass_tracks[p[0]]['intensity'], p[1], p[2]))
            return {self.name: areas}
        else:
            return {k: v for d in [x.calc_peak_areas(peaks) for x in self.composite_of] for k,v in d.items()}

    @property
    def all_nested_samples(self):
        if self.composite_of is None:
            return [self.name]
        else:
            return [j for i in [x.all_nested_samples for x in self.composite_of] for j in i]
    
    @property
    def all_nested_sample_instances(self):
        if self.composite_of is None:
            return [self]
        else:
            return [j for i in [x.all_nested_sample_instances for x in self.composite_of] for j in i]

    def get_rt_calibration_records(self):
        '''
        Returns a dictionary of sample_id, name, 
        rt_landmarks (list of apex scan numbers for the peaks used in RT calibration), 
        reverse_rt_cal_dict (key=reference scan number, value=sample specific scan number).
        '''
        return {
            'sample_id': self.sample_id,
            'name': self.name,
            'rt_landmarks': self.rt_landmarks,
            'reverse_rt_cal_dict': self.reverse_rt_cal_dict,
        }
    
    def map_rtime_to_scan_number(self, rtime):
        best_match = (999999, None)
        for scan_rtime, scan_no in zip(self.list_retention_time, self.rt_numbers):
            if abs(rtime - scan_rtime) < best_match[0]:
                best_match = (abs(rtime - scan_rtime), scan_no)
        return best_match[1]

    def scan_time_to_rtime(self, scan_no):
        return dict(zip(self.rt_numbers, self.list_retention_time))[scan_no]

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
        mode_method_map = {
            "ondisk": self._retrieve_from_disk,
            }
        mode_method_map['smart'] = mode_method_map['ondisk']
        return mode_method_map[self.database_mode]()

    def _retrieve_from_disk(self):
        '''
        Retrieve sample data from local pickle file.
        '''
        with open(self.data_location, 'rb') as f:
            sample_data = pickle.load(f)
        return sample_data