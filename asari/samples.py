'''
SimpleSample is used by default asari workflow. 
For different purposes, one can use Sample or other classes.

'''
import os

from pyopenms import MSExperiment, MzMLFile

# from mass2chem.annotate import annotate_formula_mass, massDict_hmdb

from .chromatograms import extract_massTracks_        # extract_massTraces, 
from .mass_functions import *
from .sql import *

# from .constructors import epdsConstructor


class SimpleSample:
    '''
    LC-MS Sample to get mass tracks from a mzML file, and as a container for peaks and empCpds.

    Use scan numbers whereas possible. initial peak detection for assembling empCpds only.

    Use dictionary formats for mass_trace, peak and empCpd for clarity.

    Registry for mass traces, peaks & empCpds.

    RT and m/z values will have calibration functions after CompositeMap.

    '''
    def __init__(self, experiment=None, mode='pos', input_file=''):

        self.input_file = input_file
        if experiment:
            self.experiment = experiment
            self.id = self.experiment.list_input_files.index(input_file)
        else:
            self.id = input_file
        
        self.mode = mode
        self.parameters = {}
        self.rt_numbers = []                                # list of scans, starting from 0
        self.list_retention_time = []                       # full RT time points in sample

        # These calibration functions are critical, though mz_calibration is not required if m/z accuracy is good.
        self.mz_calibration_ratio = None
        # 
        self.rt_cal_dict = None                             # index mapping to the reference sample
        self.reverse_rt_cal_dict = None
        
        # lists to store data
        self.list_mass_tracks = []                          # index number = id_number, in ascending order of m/z
        self.anchor_mz_pairs = []
        self._number_anchor_mz_pairs_ = 0
        
        # will populate after CMAP
        self.list_peaks = []   
        self.list_empCpds = []
        
        # dict for future use
        self.list_mass_traces = []
        self.dict_mass_traces = {}
        self.dict_peaks = {}


    def process(self, mz_tolerance_ppm, min_intensity, min_timepoints):
        '''
        From input file to list_MassTraces with detected peaks and selectivity on peaks.
        '''
        self.get_mass_tracks_( mz_tolerance_ppm, min_intensity, min_timepoints)
        self.get_anchor_mz_pairs()
        print("    Number of anchor m/z pairs = %d" %self._number_anchor_mz_pairs_)

        self._mz_landmarks_ = flatten_tuplelist(self.anchor_mz_pairs)
        # to send to SQL DB here


    def get_mass_tracks_(self, mz_tolerance_ppm=5, min_intensity=100, min_timepoints=5):
        '''
        A mass track is an EIC for full RT range, without separating the mass traces,
        using asari.chromatograms algorithm.
        tracks as [( mz, rtlist, intensities ), ...].
        '''
        exp = MSExperiment()                                                                                          
        MzMLFile().load(self.input_file, exp)
        xdict = extract_massTracks_(exp, 
                    mz_tolerance_ppm=mz_tolerance_ppm, min_intensity=min_intensity, min_timepoints=min_timepoints)
        self.rt_numbers = xdict['rt_numbers']            # list of scans, starting from 0
        self.list_retention_time = xdict['rt_times']     # full RT time points in sample
        ii = 0
        # already in ascending order of m/z from extract_massTracks_, get_thousandth_regions
        for tk in xdict['tracks']:                         
            self.list_mass_tracks.append( {
                'id_number': ii, 
                'mz': tk[0],
                'rt_scan_numbers': tk[1],                  # list
                'intensity': tk[2],                        # list
                } )
            ii += 1

        print("Processing %s, found %d mass tracks." %(os.path.basename(self.input_file), ii))


    def get_anchor_mz_pairs(self):
        '''
        This will be dependent on ion mode

        update self.get_anchor_mz_pairs
        e.g. [(5, 8), (6, 13), (17, 25), (20, 27), ...]
        '''
        self.anchor_mz_pairs = find_mzdiff_pairs_from_masstracks(self.list_mass_tracks, mz_tolerance_ppm=5)
        self._number_anchor_mz_pairs_ = len(self.anchor_mz_pairs)
        

    def export_mass_traces(self):
        '''
         for diagnosis etc.
        '''
        outfile = os.path.join(self.experiment.output_dir, 
                            os.path.basename(self.input_file).replace('.mzML', '.peaklist') )
        header = ['m/z', 'retention time', 'area', 'shape_quality', 'gaussian_amplitude', 'gaussian_variance', 'mz_selectivity']
        peaklist = []
        for P in self.list_mass_tracks:
            peaklist.append(str(P))

        with open(outfile, 'w') as O:
            O.write( '\t'.join(header) + '\n' + '\n'.join( peaklist ) + '\n' )
