'''
asari, LC-MS metabolomics data preprocessing - trackable, scalable.

subcommands:
    process: LC-MS data preprocessing
    xic: construct mass trakcs (chromatogram) from mzML files
    annotate: annotate a list of features
    join: merge multiple processed projects (possibly split a large dataset)
    viz: start interactive data visualization and exploration.


--params: allow passing paramters via a JSON file



Example use
-----------
python3 -m asari.main neg /Users/shuzhao/li.projects/asari/T03



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
# import sys
import time
import argparse

from asari import __version__
from .workflow import *
from .defaul_parameters import PARAMETERS

def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files.
    For OpenMS based XIC workflow, file_pattern='chrom.mzML'.
    '''
    print("Working on ~~ %s ~~ \n\n" %directory)
    return [os.path.join(directory, f) for f in os.listdir(directory) if file_pattern in f]

def process_project(list_input_files, parameters=PARAMETERS):
    '''
    list_input_files: Extracted ion chromatogram files.
    parameters: dictionary of most parameters.
    '''
    sample_registry = register_samples(list_input_files) #, dict_meta_data, parameters)
    shared_dict = batch_EIC_from_samples_ondisk(sample_registry, parameters)
    for sid, sam in sample_registry.items():
        sam['status:mzml_parsing'], sam['status:eic'], sam['number_anchor_mz_pairs'
                ], sam['data_location'] = shared_dict[sid]
        sam['name'] = os.path.basename(sam['input_file']).replace('.mzML', '')
    
    # print(sample_registry)
    EE = ext_Experiment(sample_registry, parameters)
    EE.process_all()
    EE.export_all()


def main(parameters=PARAMETERS):
    '''
    analyze: analyze a single mzML file to print summary of statistics and recommended parameters.
    process: LC-MS data preprocessing
    xic: construct mass trakcs (chromatogram) from mzML files
    annotate: annotate a list of features
    join: merge multiple processed projects (possibly split a large dataset)
    viz: start interactive data visualization and exploration.
            action='store_const', 
    '''
    print("\n\n~~~~~~~ Hello from Asari! ~~~~~~~~~\n")

    parameters['min_prominence_threshold'] = parameters['min_peak_height']/3.0
    parameters['multicores'] = min(mp.cpu_count(), parameters['multicores'])

    parser = argparse.ArgumentParser(description='asari, LC-MS metabolomics data preprocessing')

    parser.add_argument('-v', '--version', action='version', version=__version__, 
            help='print version and exit')
    parser.add_argument('run', metavar='subcommand', 
            help='one of the subcommands: analyze, process, xic, annotate, join, viz')
    parser.add_argument('-m', '--mode', default='pos', 
            help='mode of ionization, pos or neg')
    parser.add_argument('--ppm', default=5, type=int, 
            help='mass precision in ppm (part per million)')
    parser.add_argument('-i', '--input', 
            help='input directory of mzML files to process, or a single file to analyze')
    parser.add_argument('-o', '--output', 
            help='output directory')
    parser.add_argument('-j', '--project', 
            help='project name')
    parser.add_argument('-p', '--parameters', 
            help='user supplied paramter file in JSON')
    parser.add_argument('-c', '--cores', 
            help='nunmber of CPU cores intented to use')

    args = parser.parse_args()

    # update parameters
    if args.mode:
        parameters['mode'] = args.mode
    if args.ppm:
        parameters['mz_tolerance_ppm'] = args.ppm

    if args.run == 'process':
        # time_stamp is `month daay hour minute second``
        time_stamp = ''.join([str(x) for x in time.localtime()[1:6]])
        if parameters['database_mode'] == 'ondisk':
            parameters['outdir'] = '_'.join([parameters['project_name'], parameters['outdir'], time_stamp]) 
            os.mkdir(parameters['outdir'])
            os.mkdir(os.path.join(parameters['outdir'], 'pickle'))
            os.mkdir(os.path.join(parameters['outdir'], 'export'))

        process_project( read_project_dir(args.input),  parameters )        #directory = args.input

    elif args.run == 'analyze':
        # use a single sample file to analyze statistics
        analyze_single_sample(args.input, parameters=parameters)

    elif args.run == 'xic':
        # 
        pass

    elif args.run == 'annotate':
        # 
        annotate_user_featuretable(args.input, parameters=parameters)

    elif args.run == 'join':
        # input a list of directoreis, each a result of asari process
        pass

    elif args.run == 'viz':
        # launch data dashboard

        pass

    else:
        print("Expecting one of the subcommands: analyze, process, xic, annotate, join, viz.")



#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    # PARAMETERS['mode'] = sys.argv[1]
    # directory = sys.argv[2]
    # main(directory)

    main(PARAMETERS)
