import pickle
import multiprocessing as mp
from .mass_functions import flatten_tuplelist
import gzip
import os
import functools
from intervaltree import IntervalTree
import gc
import sys
import numpy as np
import scipy
import time
from scipy import interpolate

class MassTrackCache():
    def __init__(self, parameters) -> None:
        self.cache = {}
        self.in_memory = set()
        self.on_disk = set()
        self.compressed = set()
        self.pages = [[]]
        self.page_size = parameters['multicores']
        self.sparsify = parameters['sparsify']
        self.compress = parameters.get('compress', False)

    @property
    def fast_evictable(self):
        return list(self.in_memory.intersection(self.on_disk))
    
    @property
    def not_written(self):
        return list(self.in_memory.difference(self.on_disk))
    
    @property
    def compressable(self):
        return list(self.on_disk.difference(self.compressed))

    @property
    def page_full(self):
        return len(self.in_memory) >= self.page_size

    def __setitem__(self, key, value): 
        #print(round(cls.memory_use / cls.max_memory  * 100), round(cls.disk_use/(1024*2)))
        if key not in self.cache:
            if len(self.pages[-1]) == self.page_size:
                self.pages.append([])
            self.pages[-1].append(key)

            if self.sparsify:
                for mass_track in value['list_mass_tracks']:
                    mass_track['intensity'] = scipy.sparse.coo_array(mass_track['intensity'])

            self.cache[key] = {
                "value": value,
                "disk": None,
                "sparse": self.sparsify,
                "page": len(self.pages) - 1,
                "sizes": {
                    "memory": self.get_obj_size(value),
                    "disk": None,
                    "compressed": None,
                }
            }
            self.in_memory.add(key)
            self.evict()
        else:
            self.cache[key]["value"] = value
            self.in_memory.add(key)
            self.evict()

    def __getitem__(self, key):
        # Check if key is in cache
        if key in self.in_memory:
            return self.cache[key]["value"]
        else:
            if key in self.on_disk:
                if self.not_written:
                    to_write = list(self.not_written)
                    with mp.Pool(self.page_size) as workers:
                        results = workers.starmap(self.save_to_disk,[(self.cache[k]['value']['outfile'], self.cache[k]['value']) for k in to_write])
                    for k, (save_path, save_size) in zip(to_write, results):
                        self.cache[key]["value"] = None
                        self.in_memory.remove(k)
                        self.cache[k]["disk"] = save_path
                        self.on_disk.add(k)

                page_keys = [k for k in self.pages[self.cache[key]['page']] if k not in list(self.in_memory)]
            try:
                with mp.Pool(self.page_size) as workers:
                    results = workers.map(self.load_from_disk, [self.cache[k]["disk"] for k in page_keys])
                    for k, value in zip(page_keys, results):
                        if k == key:
                            content = value
                        self.cache[k]["value"] = value
                        self.in_memory.add(k)
            except BrokenPipeError:
                value = self.load_from_disk(self.cache[key]['disk'])
                self.cache[key]['value'] = value
                self.in_memory.add(key)
                content = value
        return content
    
    def evict(self):
        if self.page_full:
            # delete those on disk already
            if self.in_memory:
                for k in list(self.in_memory):
                    if k in self.on_disk:
                        self.cache[k]["value"] = None
                        self.in_memory.remove(k)
            #write unwritten to disk
            if self.in_memory:
                try:
                    with mp.Pool(self.page_size) as workers:
                        print("evicting: ")
                        results = workers.starmap(self.save_to_disk,[(self.cache[k]['value']['outfile'], self.cache[k]['value']) for k in list(self.in_memory)])
                        for k, (save_path, save_size) in zip(list(self.in_memory), results):
                            self.cache[k]["value"] = None
                            self.in_memory.remove(k)
                            self.cache[k]["disk"] = save_path
                            self.on_disk.add(k)
                except BrokenPipeError:
                    for to_evict in [(self.cache[k]['value']['outfile'], self.cache[k]['value']) for k in list(self.not_written)]:
                        save_path, _ = self.save_to_disk(to_evict)
                        self.cache[k]['value'] = None
                        self.in_memory.remove(k)
                        self.cache[k]["disk"] = save_path
                        self.on_disk.add(k)
            if self.compress:
                if self.compressable:
                    try:
                        with mp.Pool(self.page_size) as workers:
                            results = workers.map(self.compress_on_disk, [self.cache[k]['disk'] for k in list(self.compressable)])
                            print("compressing (very slow): ")
                            for k, (gz_path, _) in zip(self.compressable, results): 
                                print("\t", k)               
                                self.cache[k]['disk'] = gz_path
                                self.compressed.add(k)
                            workers.terminate()
                    except BrokenPipeError:
                        for to_compress in [self.cache[k]['disk'] for k in list(self.compressable)]:
                            gz_path, _ = self.compress_on_disk(to_compress)
                            self.cache[k]['disk'] = gz_path
                            self.compressed.add(k)

        gc.collect()


    @staticmethod
    def compress_on_disk(file_path):
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")
        output_file_path = f"{file_path}.gz"

        try:
            f_in = open(file_path, 'rb')
            f_out = gzip.open(output_file_path, 'wb', compresslevel=1)
            pickle.dump(pickle.load(f_in), f_out)
            f_in.close()
            f_out.close()
            time.sleep(1)
            os.remove(file_path)
            return output_file_path, os.path.getsize(output_file_path)
        except:
            print("failure_compressing")
            return file_path, os.path.getsize(file_path)
            

    @staticmethod
    def save_to_disk(path, data):
        if path.endswith(".pickle"):
            with open(path, 'wb') as fh:
                pickle.dump(data, fh, pickle.HIGHEST_PROTOCOL)
        elif path.endswith(".pickle.gz"):
            with gzip.GzipFile(path, 'wb', compresslevel=1) as fh:
                pickle.dump(data, fh, pickle.HIGHEST_PROTOCOL)
        return path, os.path.getsize(path)
    
    @staticmethod
    def load_from_disk(path):
        print("loading: ", path)
        if path.endswith(".pickle"):
            with open(path, 'rb') as fh:
                data = pickle.load(fh)
        elif path.endswith(".pickle.gz"):
            with gzip.GzipFile(path, 'rb') as fh:
                data = pickle.load(fh)
        return data

    @staticmethod
    def get_obj_size(obj):
        marked = {id(obj)}
        obj_q = [obj]
        sz = 0
        while obj_q:
            sz += sum(map(sys.getsizeof, obj_q))
            all_refr = ((id(o), o) for o in gc.get_referents(*obj_q))
            new_refr = {o_id: o for o_id, o in all_refr if o_id not in marked and not isinstance(o, type)}
            obj_q = new_refr.values()
            marked.update(new_refr.keys())
        return sz


class SimpleSample:
    '''
    Lightweight class of an experimental sample to facilitate workflow.
    The primary use of this class is managing data stroage and retreival at sample level.

    Depending on database_mode, sample list_mass_tracks are stored in memory, or on disk, or in MongoDB.
    Function to get mass tracks from a mzML file is in workflow.process_project and batch_EIC_from_samples_.
    Peaks and empCpds are determined in constructors.CompositeMap.
    '''
    _cached_mass_tracks = None

    def __init__(self, registry={}, experiment=None, database_mode='ondisk', mode='pos', is_reference=False):
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
        if SimpleSample._cached_mass_tracks is None:
            SimpleSample._cached_mass_tracks = MassTrackCache(self.experiment.parameters)


        self.__dict__.update(registry)

        self.experiment = experiment
        self.is_reference = is_reference 
        self.is_rt_aligned = is_reference
            
        self.rt_landmarks = []  # to populate at CMAP.calibrate_sample_RT

        # These are critical RT calibration functions, index mapping with the reference sample
        self.rt_cal_dict = None
        self.reverse_rt_cal_dict = None
        
        # placeholder
        self.mz_calibration_function = None
        self._reference_tracks = None

        self.calibrations = []
        self.reverse_calibrations = []

    @functools.cached_property
    def mz_tree(self):
        __mz_tree = IntervalTree()
        mz_tol = self.experiment.mz_tolerance_ppm
        for t in self.list_mass_tracks:
            t_mz = t['mz']
            t_mz_err = t_mz / 1e6 * mz_tol * 4
            __mz_tree.addi(t_mz - t_mz_err, t_mz + t_mz_err, t['id_number'])
        return __mz_tree

    @property
    def database_mode(self):
        return self.experiment.database_mode

    @property
    def mode(self):
        return self.experiment.mode

    @property
    def rt_numbers(self):
        return self.list_scan_numbers

    @property
    def _mz_landmarks_(self):
        return flatten_tuplelist(self.anchor_mz_pairs)  

    @property
    def rt_calibration_records(self):
        return {
            'sample_id': self.sample_id,
            'name': self.name,
            'rt_landmarks': self.rt_landmarks,
            'reverse_rt_cal_dict': self.reverse_rt_cal_dict,
        }

    @property
    def list_mass_tracks(self):
        if self.is_reference:
            if self._reference_tracks is None: 
                self._reference_tracks = SimpleSample._cached_mass_tracks[self.data_location]['list_mass_tracks']
            return self._reference_tracks
        return SimpleSample._cached_mass_tracks[self.data_location]['list_mass_tracks']

    @staticmethod
    def save(data, parameters):
        if SimpleSample._cached_mass_tracks is None:
            SimpleSample._cached_mass_tracks = MassTrackCache(parameters)
        SimpleSample._cached_mass_tracks[data['data_location']] = data['sample_data']

