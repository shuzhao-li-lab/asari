'''
SimpleSample as a wrapper for a sample, lightweight to facilitate workflow.
Depending on database_mode, sample list_mass_tracks are stored in memory, or on disk, or in MongoDB.
'''
import pickle
import os

import pymzml

#from pyopenms import MSExperiment, MzMLFile

from .chromatograms import extract_massTracks_        # extract_massTracks, 
from .mass_functions import *


class SimpleSample:
    '''
    LC-MS Sample to get mass tracks from a mzML file, and as a container for peaks and empCpds.
    Use scan numbers whereas possible. 
    Use dictionary formats for mass_track, etc for clarity.
    {
    input_file: '',
    name: '',                               # usually short form of input_file
    ion_mode: '',
    # status
    'status:mzml_parsing': '',
    'status:eic': '',
    'mass_aligned': '',
    'rtime_aligned': '',
    'peaks_extracted': '',
    #
    list_scan_numbers = [],
    list_retention_time: [],
    list_mass_tracks: [
        {'id_number': ii, 
                'mz': float,
                'rt_scan_numbers': [],
                'intensity': [],
                }, ...
    ], 
    anchor_mz_pairs: [],                    # mostly isotopic pairs to establish m/z anchors (landmarks)
    number_anchor_mz_pairs: int,
}

    Registry for mass traces, peaks & empCpds.
    RT and m/z values will have calibration functions after CompositeMap.

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
            return self.list_mass_tracks
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



#---------------------------------------------------------------------------------------------------------------


class fullSample(SimpleSample):

    def process(self, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height):
        '''From input file to list_MassTraces with detected peaks and selectivity on peaks.
        '''
        self.generate_mass_tracks( mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height)

    def generate_mass_tracks(self, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
        '''
        A mass track is an EIC for full RT range, without separating the mass traces,
        using asari.chromatograms algorithm.
        tracks as [( mz, rtlist, intensities ), ...].
        '''
        list_mass_tracks = []

        exp = pymzml.run.Reader(self.input_file)
        
        # exp = MSExperiment()                                                                                          
        # MzMLFile().load(self.input_file, exp)

        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    min_intensity=min_intensity, 
                    min_timepoints=min_timepoints, 
                    min_peak_height=min_peak_height)
        self.rt_numbers = xdict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xdict['rt_times']     # full RT time points in sample
        ii = 0
        # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
        for track in xdict['tracks']:                         
            list_mass_tracks.append( {
                'id_number': ii, 
                'mz': track[0],
                'rt_scan_numbers': track[1],                  # list
                'intensity': track[2],                        # list
                } )
            ii += 1

        print("Processing %s, found %d mass tracks." %(os.path.basename(self.input_file), ii))
        self.generate_anchor_mz_pairs(list_mass_tracks)
        print("    Number of anchor m/z pairs = %d" %self._number_anchor_mz_pairs_)

        if self.database_mode == 'memory':
            self.list_mass_tracks = list_mass_tracks
        elif self.database_mode == 'ondisk': 
            self.push_to_disk(list_mass_tracks)
        else:                         # requiring self.experiment and db connection
            self.push_to_db(list_mass_tracks
                # to implement
            )

    def generate_mass_tracks_openms(self, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5, min_peak_height=1000):
        '''
        A mass track is an EIC for full RT range, without separating the mass traces,
        using asari.chromatograms algorithm.
        tracks as [( mz, rtlist, intensities ), ...].
        '''
        list_mass_tracks = []
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    min_intensity=min_intensity, 
                    min_timepoints=min_timepoints, 
                    min_peak_height=min_peak_height)
        self.rt_numbers = xdict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xdict['rt_times']     # full RT time points in sample
        ii = 0
        # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
        for track in xdict['tracks']:                         
            list_mass_tracks.append( {
                'id_number': ii, 
                'mz': track[0],
                'rt_scan_numbers': track[1],                  # list
                'intensity': track[2],                        # list
                } )
            ii += 1

        print("Processing %s, found %d mass tracks." %(os.path.basename(self.input_file), ii))
        self.generate_anchor_mz_pairs(list_mass_tracks)
        print("    Number of anchor m/z pairs = %d" %self._number_anchor_mz_pairs_)

        if self.database_mode == 'memory':
            self.list_mass_tracks = list_mass_tracks
        elif self.database_mode == 'ondisk': 
            self.push_to_disk(list_mass_tracks)
        else:                         # requiring self.experiment and db connection
            self.push_to_db(list_mass_tracks
                # to implement
            )

    def push_to_disk(self, list_mass_tracks):
        '''Write pickle file for list_mass_tracks under project directory, pickle/samplename.pickle.
        '''
        with open(self.pickle_file, 'wb') as f:
            pickle.dump(list_mass_tracks, f, pickle.HIGHEST_PROTOCOL)
