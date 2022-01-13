'''
LC-MS metabolomics data pre-processing

Example use
-----------
python3 -m asari.main neg /Users/shuzhao/li.projects/asari/T03

'''

import sys
from .workflow import *


PARAMETERS = {
    'project_name': 'test_asari',
    
    'min_intensity_threshold': 5000,    # minimal peak height
    'min_timepoints': 5,                # minimal number of data points in elution profile. scipy find_peaks treat `width` as FWHM, thus half of this value.
    'signal_noise_ratio': 2,            # peak height at least x fold over noise, which is median of non-peak data points.
    #
    'peak_number_rt_calibration': 5,    # minimal number of selected high-quality peaks required for RT calibration. 
                                        # Samples with fewer selected peaks are dropped out.
    'cache_mass_traces': False,         # to save memory if not using DB; turn on if need to plot and diagnose
    'output_filename': 'feature_table.tsv',
    'annotation_filename': "annotation_table.tsv",
    #
    'mode': 'pos',                      # ionization mode
    'mass_range': (50, 2000),
    'max_rtime': 300,                   # retention time range (chromatography) 0-300 seconds

    'mz_tolerance': 5,                  # ppm, high selectivity meaning no overlap neighbors to confuse; 
                                        # Low selectivity regions will be still inspected to determine the true number of features

    'rtime_tolerance': 10,              # feature rtime shift threshold under 10 seconds; or 10% of rtime
                                        # will change to automated parameter using stdev??
    #
    'initiation_samples': [],           # if user to specify 3 samples to initiate data processing, to init HOT_DB; 
                                        # otherwise they are chosen automatically

    # no need to modify below unless you know what you are doing
    'prominence_window': 30,
    'gaussian_shape': 0.3,              # min cutoff
    }
PARAMETERS['min_prominence_threshold'] = PARAMETERS['min_intensity_threshold']/3.0


def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files.
    For OpenMS based XIC workflow, file_pattern='chrom.mzML'.
    '''
    print("\nWorking on ", directory)
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
    print("\n\n~~~~~~~ Hello from Asari! ~~~~~~~~~\n")
    process_project(
            read_project_dir(directory), {}, PARAMETERS, directory   #setting output_dir as input dir
    )

#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    PARAMETERS['mode'] = sys.argv[1]
    directory = sys.argv[2]
    main(directory)
