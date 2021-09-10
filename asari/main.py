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


PARAMETERS = {
    'mode': 'pos',                      # ionization mode
    'min_intensity_threshold': 30000,   # minimal peak intensity
    'min_timepoints': 5,                # minimal number of data points in elution profile
    #
    'initiation_samples': [],           # if user to specify 3 samples to initiate data processing, to init HOT_DB; 
                                        # otherwise they are chosen automatically
    # no need to modify below
    'prominence_window': 30,
    'gaussian_shape': 0.8,
}

PARAMETERS['min_prominence_threshold'] = PARAMETERS['min_intensity_threshold']/3.0


def read_projec_dir(directory, file_pattern='chrom.mzML'):
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

def process_project(list_input_files, dict_meta_data={}, parameters=PARAMETERS):
    '''
    Use ext_Experiment as a containing class to hold processed data.

    ionization_mode='positive'

    EE.assign_formula_masses()
    EE.calibrate_retention_time()
    EE.correspondency()
    EE.export_feature_table('out.tsv')


    '''
    if dict_meta_data:
        for k in dict_meta_data:
            dict_meta_data[k] = dict_meta_data[k].upper()       # upper characters to standardize
            
    EE = ext_Experiment()
    EE.__init2__(list_input_files, dict_meta_data, parameters)
    
    #EE.set_sample_order(list_input_files)
    
    for f in list_input_files:
        SM = Sample(f)
        SM.process_step_1()
        SM.export_peaklist()
        #
        #EE.samples.append(SM)




def main(directory):
    process_project(
            read_projec_dir(directory), {}, 
    )
    


#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':
    directory = sys.argv[1]
    main(directory)
