'''
LC-MS metabolomics data pre-processing

Example use
-----------
python3 -m asari.main neg /Users/shuzhao/li.projects/asari/T03

asari, LC-MS metabolomics data preprocessing - trackable, scalable.

In the asari/mummichog packages, the data entities are presented in any of the four types: 
class, namedtuple, JSON style dictionary or implicit list. 
The implicit lists are used sparely as they have reduced clarity. 
Namedtuple is immutable, then limited applications. 
In-memory searches are conducted using indexed dictionaries or dataframes.


Data formats:
===============
mass tracks as [( mz, rtlist, intensities ), ...].
Peak format: 
{
    'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 'height', 'parent_masstrace_id', 
    'rtime', 'peak_area', 'goodness_fitting'
}
isotopic_patterns = [(1.003355, 'M(13C)', 0, 0.2), ...]

Mass Tracks
===========
They are used for full RT ranges, thus each mass track has a unique m/z. 
Some chromatogram builders separate the mass traces if there are gaps in RT scans, 
but that creates complexity in m/z alignment and searches. 

Peak detection
==============
The main step uses scipy.signal.find_peaks, a local maxima method with prominence control.
Prominence is important, but it should be more tailored to individual peaks. 
Here, 
prominence = max(min_prominence_threshold, 0.05*max(list_intensity)).

Empirical Compound, list and tree
=================================
Constructing trees of emperical compounds (epdTree) from a list of peaks or features,
which follow a format:
list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...]

Steps:
1. find isotopic signatures from list_peaks
2. initiate epdTree classes using the above isotopic signatures; extend by searching common adducts
3. In remaining peaks (isotopes not seen due to low intensity etc), search pairs for primary adducts. 
   Initiate epdTree classes using those.
4. Extend to other adducts for empCpds from the above steps.
5. Reconsolidate overlap epdTrees. E.g. combinations of isotopes and adducts are found separately.

Notes:
a. Common anchor peaks are M+H or M-H etc, but will calculate at later round. M* may not show.
b. Adducts and isotopes are combinatorial, under restriction of chemical formulae.
c. We curated isotopic/adduct patterns in this package.
d. epdTree can be improved in future based on real data statistics and more structured cheminformatics.

'''

import sys
import time
from .workflow import *


PARAMETERS = {
    'project_name': 'test_asari',
    'outdir': 'asari_output_',
    
    'min_intensity_threshold': 100,     # minimal intensity for mass track extraction, filtering baseline
    'min_peak_height': 5000,            # minimal peak height
    'cal_min_peak_height': 100000,      # minimal peak height required for peaks used for RT calibration
    'min_timepoints': 5,                # minimal number of data points in elution profile. scipy find_peaks treat `width` as FWHM, thus half of this value.
    'signal_noise_ratio': 2,            # peak height at least x fold over noise, which is median of non-peak data points.
    #
    'peak_number_rt_calibration': 15,   # minimal number of selected high-quality peaks required for RT calibration. 
                                        # Samples with fewer selected peaks are dropped out.

    # tweaking this section now -
    'cache_mass_traces': False,         # to save memory if not using DB; turn on if need to plot and diagnose
    'output_feature_table': '_cmap_feature_table.csv',
    'mass_grid_mapping': "_mass_grid_mapping.csv",
    'annotation_filename': "annotation_table.tsv",
    'json_empricalCompounds': "_empCpd_json.json",

    #
    'mode': 'pos',                      # ionization mode
    'mass_range': (50, 2000),
    #'max_rtime': 300,                   # retention time range (chromatography) 0-300 seconds

    'mz_tolerance': 5,                  # ppm, high selectivity meaning no overlap neighbors to confuse; 
                                        # Low selectivity regions will be still inspected to determine the true number of features
    'rtime_tolerance': 10,              # feature rtime shift threshold under 10 seconds; or 10% of rtime
                                        # not used now?
    #
    'initiation_samples': [],           # if user to specify 3 samples to initiate data processing, to init HOT_DB; 
                                        # otherwise they are chosen automatically

    # no need to modify below unless you know what you are doing
    'prominence_window': 30,
    'gaussian_shape': 0.3,              # min cutoff
    }
    
PARAMETERS['min_prominence_threshold'] = PARAMETERS['min_peak_height']/3.0


def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files.
    For OpenMS based XIC workflow, file_pattern='chrom.mzML'.
    '''
    print("Working on ~~ %s ~~ \n\n" %directory)
    return [os.path.join(directory, f) for f in os.listdir(directory) if file_pattern in f]

def metafile_to_dict(infile):
    '''
    Optional.
    Meta data file, tab delimited, first two columns corresponding to [file_name, sample_type].
    '''
    meta = {}
    for line in open(infile).read().splitlines():
        a = line.split('\t')
        meta[a[0]] = a[1]
    return {}

def process_project(list_input_files, dict_meta_data={}, parameters=PARAMETERS, output_dir=''):
    '''
    list_input_files: Extracted ion chromatogram files.
    parameters: dictionary of most parameters.
    '''
    if dict_meta_data:
        for k in dict_meta_data:
            dict_meta_data[k] = dict_meta_data[k].upper()       # upper characters to standardize
    if not list_input_files:
        print("No input file found. Please verify your pathway to files.")

    EE = ext_Experiment()
    EE.__init2__(list_input_files, dict_meta_data, parameters, output_dir)

    EE.process_all()


def main(directory):
    time_stamp = str(time.time())
    PARAMETERS['outdir'] += time_stamp
    os.mkdir(PARAMETERS['outdir'])

    print("\n\n~~~~~~~ Hello from Asari! ~~~~~~~~~\n")
    process_project(
            read_project_dir(directory), 
            {},         # not used now
            PARAMETERS, 
            PARAMETERS['outdir']   #setting output_dir as input dir
    )

#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    PARAMETERS['mode'] = sys.argv[1]
    directory = sys.argv[2]
    main(directory)
