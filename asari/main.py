'''
LC-MS metabolomics data pre-processing

- only for high resolution data
- formula mass centric
- local maxima peak detection with prominence control

'''

import os, sys
from .algorithms import Sample, ext_MassTrace, ext_Experiment


def read_projec_dir(directory, file_pattern='chrom.mzML'):
    return [os.path.join(directory, f) for f in os.listdir(directory) if file_pattern in f]

def metafile_to_dict(infile):
    # placeholder
    return {}

def process_project(list_input_files, dict_meta_data, ionization_mode='positive'):
    '''
    Use ext_Experiment as a containing class to hold processed data.



    EE.assign_formula_masses()
    EE.calibrate_retention_time()
    EE.correspondency()
    EE.export_feature_table('out.tsv')


    '''
    EE = ext_Experiment()
    EE.set_sample_order(list_input_files)
    for f in list_input_files:
        SM = Sample(f)
        SM.detect_peaks()
        SM.assign_selectivity()
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
