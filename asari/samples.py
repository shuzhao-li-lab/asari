import pickle
import multiprocessing as mp
from .mass_functions import flatten_tuplelist
import gzip
import os
import pandas as pd
import functools
from intervaltree import IntervalTree
from functools import lru_cache, partial
from collections.abc import Iterable
import atexit
import threading
import gc
import sys
import asyncio
import aiofiles

def save_to_disk(path, data):
    if path.endswith(".pickle"):
        with open(path, 'wb') as fh:
            pickle.dump(data, fh, pickle.HIGHEST_PROTOCOL)
    elif path.endswith(".pickle.gz"):
        with gzip.GzipFile(path, 'wb', compresslevel=1) as fh:
            pickle.dump(data, fh, pickle.HIGHEST_PROTOCOL)
    return path, {}, get_obj_size(data)

def load_from_disk(path, direct=True):
    if path.endswith(".pickle"):
        data = pickle.load(open(path, 'rb'))
    elif path.endswith(".pickle.gz"):
        data = pickle.load(gzip.GzipFile(path, 'rb'))
    return path, data if not direct else data

def get_obj_size(obj):
    marked = {id(obj)}
    obj_q = [obj]
    sz = 0

    while obj_q:
        sz += sum(map(sys.getsizeof, obj_q))

        # Lookup all the object referred to by the object in obj_q.
        # See: https://docs.python.org/3.7/library/gc.html#gc.get_referents
        all_refr = ((id(o), o) for o in gc.get_referents(*obj_q))

        # Filter object that are already marked.
        # Using dict notation will prevent repeated objects.
        new_refr = {o_id: o for o_id, o in all_refr if o_id not in marked and not isinstance(o, type)}

        # The new obj_q will be the ones that were not marked,
        # and we will update marked with their ids so we will
        # not traverse them again.
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

    mass_track_cache = {}
    sizes = {}
    sample_order = []
    order_map = {}
    memory_limit = 50 * 1024 ** 3
    memory_use = 0
    enable_async_preload = True
    loop = asyncio.new_event_loop()
    _loop_thread = None
    _cache_lock = threading.Lock()
    memory_mode_achieved = False
    

    @classmethod
    async def __load_from_disk_async(cls, path):
        """
        Asynchronously load data from disk.
        
        :param path: The file path to load.
        :return: The loaded object.
        """
        if path.endswith(".pickle"):
            return await asyncio.to_thread(pickle.load, open(path, 'rb'))
        elif path.endswith(".pickle.gz"):
            return await asyncio.to_thread(pickle.load, gzip.GzipFile(path, 'rb'))

    @classmethod
    def start_preloading(cls):
            """
            Start the background preloading process in a separate thread.
            """
            print("Starting preloading...")
            atexit.register(cls.stop_event_loop)
            import threading
            # Run the event loop in a separate thread
            if cls._loop_thread is None:
                cls._loop_thread = threading.Thread(target=cls._run_event_loop, daemon=True)
                cls._loop_thread.start()

            asyncio.run_coroutine_threadsafe(cls._preload_files_async(), cls.loop)

    @classmethod
    def _run_event_loop(cls):
        """
        Run the event loop.
        """
        asyncio.set_event_loop(cls.loop)
        cls.loop.run_forever()

    @classmethod
    def stop_event_loop(cls):
        if cls.loop.is_running():
            cls.loop.call_soon_threadsafe(cls.loop.stop)
            cls._loop_thread.join()

    @classmethod
    async def _preload_files_async(cls):
        """
        Asynchronously preload files in the background.
        """
        while SimpleSample.enable_async_preload and cls.memory_use < cls.memory_limit:
            for data_location in SimpleSample.sample_order:
                if data_location not in SimpleSample.mass_track_cache:
                    if cls.memory_use + cls.sizes[data_location] < cls.memory_limit:
                        print("preloading: ", data_location)
                        with cls._cache_lock:
                            SimpleSample.mass_track_cache[data_location] = 'async_load'
                        # Load the file asynchronously
                        content = await cls.__load_from_disk_async(data_location)
                        with cls._cache_lock:
                            if SimpleSample.mass_track_cache[data_location] == 'async_load':
                                print("preloaded: ", data_location)
                                SimpleSample.mass_track_cache[data_location] = content['list_mass_tracks']
                                cls.memory_use += cls.sizes[data_location]
            await asyncio.sleep(0)  # Yield control to the event loop

    #@classmethod
    #def populate_order(cls, directory):
        #cls.sample_order = sorted([directory + x for x in os.listdir(directory)])
        #cls.order_map = {x: i for i, x in enumerate(cls.sample_order)}

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

        self.__dict__.update(registry)
        SimpleSample.sizes[self.data_location] = self.size

        self.experiment = experiment
        self.is_reference = is_reference 
        self.is_rt_aligned = is_reference
            
        self.rt_landmarks = []  # to populate at CMAP.calibrate_sample_RT

        # These are critical RT calibration functions, index mapping with the reference sample
        self.rt_cal_dict = None
        self.reverse_rt_cal_dict = None
        
        # placeholder
        self.mz_calibration_function = None
        self._cached_mass_tracks = None
        self.memory_mode_achieved = False

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
        if SimpleSample.memory_mode_achieved:
            #for k, v in SimpleSample.mass_track_cache.items():
            #    print(k, type(v))
            #    if isinstance(v, str):
            #        print("\t", v)
            #exit()
            return SimpleSample.mass_track_cache[self.data_location]
        elif self.sample_data:
            return self.sample_data['list_mass_tracks']
        elif self._cached_mass_tracks:
            return self._cached_mass_tracks
        else:
            SimpleSample.stop_event_loop()
            with SimpleSample._cache_lock:
                to_del = []
                for x, v in SimpleSample.mass_track_cache.items():
                    if isinstance(v, str):
                        to_del.append(x)
                for x in to_del:
                    del SimpleSample.mass_track_cache[x]
                    
                if self.data_location in SimpleSample.mass_track_cache:
                    current_content = SimpleSample.mass_track_cache[self.data_location]
                    if not isinstance(current_content, str):
                        return current_content
                SimpleSample.mass_track_cache[self.data_location] = 'sync_load'
            start = SimpleSample.order_map[self.data_location]
            to_load = [self.data_location]
            if self.is_reference:
                with SimpleSample._cache_lock:
                    results = [load_from_disk(self.data_location, direct=True)]
                    SimpleSample.mass_track_cache.update({x: y['list_mass_tracks'] for x,y in results})
                    mass_tracks = SimpleSample.mass_track_cache[self.data_location]
                    self._cached_mass_tracks = mass_tracks
                    return mass_tracks
            else:
                with SimpleSample._cache_lock:
                    # skip this if the reference sample as the reference sample can be requested out of order
                    if not self.is_reference:
                        ii = 0
                        preloadable = set()
                        for x in SimpleSample.sizes:
                            if x != self.data_location:
                                if x in SimpleSample.mass_track_cache:
                                    if isinstance(SimpleSample.mass_track_cache[x], str):
                                        preloadable.add(x)
                                else:
                                    preloadable.add(x)
                        while preloadable:
                            ii += 1
                            to_preload = SimpleSample.sample_order[(start + ii) % len(SimpleSample.sample_order)]
                            if to_preload in preloadable:
                                preloadable.remove(to_preload)
                                preload_size = SimpleSample.sizes[to_preload]
                                chunk_mem_size = 0
                                if SimpleSample.memory_use + preload_size < SimpleSample.memory_limit:
                                    if to_preload not in SimpleSample.mass_track_cache:
                                        to_load.append(to_preload)
                                        chunk_mem_size += preload_size
                                        SimpleSample.memory_use += preload_size
                                        SimpleSample.mass_track_cache[to_preload] = 'sync_load'
                                else:
                                    break
                with mp.Pool(min(self.experiment.parameters['multicores'], len(to_load))) as workers:
                    results = workers.map(load_from_disk, to_load)
                with SimpleSample._cache_lock:
                    if SimpleSample.memory_use + min(SimpleSample.sizes.values()) > SimpleSample.memory_limit:    
                        print("Dropping Cache")
                        SimpleSample.mass_track_cache = {x: y['list_mass_tracks'] for x,y in results}
                        SimpleSample.memory_use = chunk_mem_size
                    else:
                        print("Updating Cache")
                        SimpleSample.mass_track_cache.update({x: y['list_mass_tracks'] for x,y in results})
                    SimpleSample.stop_event_loop() # probably not needed
                    if not [x for x in SimpleSample.sizes if x not in SimpleSample.mass_track_cache]:
                        SimpleSample.memory_mode_achieved = True
                    else:
                        SimpleSample.start_preloading()
                    return SimpleSample.mass_track_cache[self.data_location]
    
    def tracks_by_mz(self, query_mass):
        track_set = set()
        for x in self.mz_tree.at(float(query_mass) - 1.0072764665789):
            track_set.add(x.data)
        return track_set
    
    def retrieve_tracks_id(self, id):
        if isinstance(id, Iterable):
            return [self.retrieve_tracks_id(x) for x in id]
        return self.retrieve_track_id(id)

    def retrieve_track_id(self, id):
        import matplotlib.pyplot as plt
        for x in self.list_mass_tracks:
            if x['id_number'] == id:
                return x
    
            
    def find_kovats(self, kovats_csv="/Users/mitchjo/asari/asari/db/kovats.csv"):
        import matplotlib.pyplot as plt
        from .peaks import stats_detect_elution_peaks

        kovats = pd.read_csv(kovats_csv)
        for kovat in kovats.to_dict(orient='records'):
            for t in self.retrieve_tracks_id(self.tracks_by_mz(kovat['mass'])):
                EP = stats_detect_elution_peaks(t,                                           
                                           len(t['intensity']), 
                                           self.experiment.parameters['min_peak_height'],
                                           self.experiment.parameters['min_peak_ratio'],
                                           round(0.5 * self.experiment.parameters['min_timepoints']),
                                           self.experiment.parameters['min_prominence_threshold'],
                                           self.experiment.parameters['wlen'],
                                           self.experiment.parameters['signal_noise_ratio'] * 100,
                                           self.experiment.parameters['gaussian_shape'],
                                           .02,
                                           False,
                                           self.experiment.parameters['min_intensity_threshold'])
                if EP:
                    plt.scatter(range(len(t['intensity'])), t['intensity'], c='k')
                    for p in EP:
                        print(p['apex'], p['height'])
                        plt.scatter(p['apex'], p['height'], c='r')
                    plt.show()
                print(kovat, len(EP))
                for p in EP:
                    print("\t", p)

        exit()

    def find_kovats2(self, kovats_csv="/Users/mitchjo/asari/asari/db/kovats.csv"):
        kovats = pd.read_csv(kovats_csv)
        kovats_hits = []
        for kovat in kovats.to_dict(orient='records'):
            kovat_result = dict(kovat)
            matching_track_ids = self.tracks_by_mz(kovat['mass'])
            matching_tracks = self.retrieve_tracks_id(matching_track_ids)
            likely_track = {'mass_track_id': None,
                            'max_intensity': 0,
                            'apex': None}
            for m_t in matching_tracks:
                apex_tracker = {
                    'intensity': 0,
                    'apex': None
                }
                for i, intensity in enumerate(m_t['intensity']):
                    if intensity > apex_tracker['intensity']:
                        apex_tracker['intensity'] = intensity 
                        apex_tracker['apex'] = i
                if apex_tracker['intensity'] > likely_track['max_intensity']:
                    likely_track['mass_track_id'] = m_t['id_number']
                    likely_track['max_intensity'] = apex_tracker['intensity']
                    likely_track['apex'] = apex_tracker['apex']
            if likely_track['apex']:
                kovat_result.update(likely_track)
                kovats_hits.append(kovat_result)
        return kovats_hits
                





