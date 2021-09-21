'''
LC-MS metabolomics data pre-processing

- only for high resolution data
- formula mass centric
- local maxima peak detection with prominence control


Typical chromatogram (XIC) extraction from mzML raw files:

> for file in ../mzML_IS_HILICposRPneg_05072021/*.mzML
>   do FeatureFinderMetabo -in $file -out ${file/.mzML/.featureXML} -out_chrom ${file/.mzML/_chrom.mzML} \
>   -algorithm:common:chrom_fwhm 2 -algorithm:mtd:mass_error_ppm 2 -algorithm:mtd:min_trace_length 5
> done

Parameters above: 2, 2, 5 seconds
Ref: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html



'''

import os, sys
from .algorithms import Sample, ext_MassTrace, ext_Experiment
# from .plot import plot_sample_rt_calibration

PARAMETERS = {
    'min_intensity_threshold': 5000,   # minimal peak height
    'min_timepoints': 5,                # minimal number of data points in elution profile
    #
    'peak_number_rt_calibration': 20,   # minimal number of selected high-quality peaks required for RT calibration. Samples with fewer selected peaks are dropped out.
    'cache_mass_traces': False,         # to save memory if not using DB; turn on if need to plot and diagnose
    'output_filename': 'feature_table.tsv',
    'annotation_filename': "annotation_db.tsv",
    #
    'mode': 'pos',                      # ionization mode
    'mass_range': (50, 2000),
    'max_rtime': 300,                   # retention time range (chromatography) 0-300 seconds
    # 'interpolate_factor': 10,           # per second. Can increase to 100 for very fast scan rate.
    #
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


def read_project_dir(directory, file_pattern='chrom.mzML'):
    print("\nWorking on ", directory)
    return [os.path.join(directory, f) for f in os.listdir(directory) if file_pattern in f]

def metafile_to_dict(infile):
    '''
    Meta data file, tab delimited, first two columns corresponding to [file_name, sample_type].
    '''
    meta = {}
    for line in open(infile).read().splitlines():
        a = line.split('\t')
        meta[a[0]] = a[1]
    return {}

def process_project(list_input_files, dict_meta_data={}, parameters=PARAMETERS, output_dir=''):
    '''
    Use ext_Experiment as a containing class to hold processed data.

    ionization_mode='positive'

    EE.assign_formula_masses()
    EE.calibrate_retention_time()
    EE.correspondency()
    EE.export_feature_table('out.tsv')

    if not output_dir:
        output_dir = './'
        
    for SM in EE.samples:
        # plot_sample_rt_calibration(SM)
        SM.export_peaklist()

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
    directory = sys.argv[1]
    main(directory)
