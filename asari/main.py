import argparse
import multiprocessing as mp
import os
import yaml
import time
import sys

from asari import __version__
from .workflow import (get_mz_list, process_project, process_xics, read_project_dir, create_export_folders)
from .default_parameters import PARAMETERS
from .dashboard import read_project, dashboard
from .analyze import estimate_min_peak_height, analyze_single_sample
from .annotate_user_table import annotate_user_featuretable
from .utils import build_boolean_dict, bulk_process
from .qc import generate_qc_report

booleandict = build_boolean_dict()

PARAMETERS['asari_version'] = __version__

def __run_process__(parameters, args):
    # main process function
    list_input_files = read_project_dir(args.input)
    if not list_input_files:
        print("No valid mzML files are found in the input directory :(")
    else:
        process_project(list_input_files, parameters)

def convert(parameters, args):
    needs_conversion = []
    for x in os.listdir(args.input):
        if x.endswith('.raw'):
            needs_conversion.append(x)
    if needs_conversion:
        from .mzml_converter import mzMLconverter
        converter = mzMLconverter()
        converter.bulk_convert(needs_conversion)

def process(parameters, args):
    __run_process__(parameters, args)

def analyze(parameters, args):
    analyze_single_sample(args.input, parameters=parameters)

def xic(parameters, args):
    list_input_files = read_project_dir(args.input)
    process_xics(list_input_files, parameters)

def extract(parameters, args):
    mzlist = get_mz_list(args.target)
    print("Retrieved %d target mz values from %s.\n" %(len(mzlist), args.target))
    parameters['target'] = mzlist
    __run_process__(parameters, args)

def annotate(parameters, args):
    annotate_user_featuretable(args.input, parameters=parameters, rtime_tolerance=2)

def join(parameters, args):
    print("NOT IMPLEMENTED")

def viz(parameters, args):
    datadir = args.input
    project_desc, cmap, epd, Ftable, Ptable = read_project(datadir)
    if args.table_for_viz == 'full':
        dashboard(project_desc, cmap, epd, Ftable)
    elif args.table_for_viz == 'preferred':
        dashboard(project_desc, cmap, epd, Ptable)
    else:
        raise ValueError("Table for visualization must be either 'preferred' or 'full'.")

def qc_report(parameters, args):
    list_input_files = read_project_dir(args.input)
    create_export_folders(parameters)
    jobs = [(f, os.path.join(parameters['qaqc_reports_outdir'], os.path.basename(f).replace(".mzML", "_report.html")), parameters['spikeins']) for f in list_input_files]
    bulk_process(generate_qc_report, jobs)

def update_peak_detection_params(parameters, args):
    if parameters['autoheight']:
        try:
            parameters['min_peak_height'] = estimate_min_peak_height(read_project_dir(args.input))
            parameters['min_intensity_threshold'] = parameters['min_peak_height'] / 10
        except ValueError as err:
            print("Problems with input files: {0}. Back to default min_peak_height.".format(err))
            
    elif args and args.min_peak_height:
        try:
            parameters['min_peak_height'] = float(args.min_peak_height)
        except:
            print("Problems with specified min_height. Back to default min_peak_height.")

    parameters['min_prominence_threshold'] = int( 0.33 * parameters['min_peak_height'] )
    parameters['cal_min_peak_height'] = 10 * parameters['min_peak_height']

    if args and args.min_prominence_threshold:
        try:
            parameters['min_prominence_threshold'] = float(args.min_prominence_threshold)
        except ValueError as err:
            print("Problems with specified min_prominence_threshold. Back to default min_prominence_threshold.")

    if args and args.cal_min_peak_height:
        try:
            parameters['cal_min_peak_height'] = float(args.cal_min_peak_height)
        except ValueError as err:
            print("Problems with specified cal_min_peak_height. Back to default cal_min_peak_height.")

    if args and args.min_intensity_threshold:
        try:
            parameters['min_intensity_threshold'] = float(args.min_intensity_threshold)
        except ValueError as err:
            print("Problems with specified min_intensity_threshold. Back to default min_intensity_threshold.")                    
    return parameters


def main(parameters=PARAMETERS):
    '''
    asari, Trackable and scalable Python program for high-resolution LC-MS metabolomics data preprocessing.

        * analyze: analyze a single mzML file to print summary of statistics and recommended parameters.
        * process: LC-MS data preprocessing
        * xic: construct mass tracks (chromatogram) from mzML files
        * extract: targeted extraction of given m/z list
        * annotate: annotate a list of features
        * join: merge multiple processed projects (possibly split a large dataset)
        * viz: start interactive data visualization and exploration.

    all mzML files should be centroided.

    Parameters
    ----------
    parameters : dict
        This dictionary contains a number of key value pairs that determine the behavior of various aspects
        of the asari processing. The parameters can be seen in default_parameters.py. Command line arguments
        will override any defaults and any values provided in the parameters.json file. 
    '''


    print("\n\n~~~~~~~ Hello from Asari (%s) ~~~~~~~~~\n" %__version__)

    parser = argparse.ArgumentParser(description='asari, LC-MS metabolomics data preprocessing')

    parser.add_argument('-v', '--version', action='version', version=__version__, 
            help='print version and exit')
    parser.add_argument('run', metavar='subcommand', 
            help='one of the subcommands: analyze, process, xic, extract, annotate, join, viz')
    parser.add_argument('-m', '--mode', default='pos', 
            help='mode of ionization, pos or neg')
    parser.add_argument('--ppm', default=5, type=int, 
            help='mass precision in ppm (part per million), same as mz_tolerance_ppm')
    parser.add_argument('-i', '--input', 
            help='input directory of mzML files to process, or a single file to analyze')
    parser.add_argument('-o', '--output', 
            help='output directory')
    parser.add_argument('-j', '--project', 
            help='project name')
    parser.add_argument('-p', '--parameters', 
            help='Custom paramter file in YAML. Use parameters.yaml as template.')
    parser.add_argument('-c', '--multicores', type=int, default=0, 
            help='nunmber of CPU cores intented to use')
    parser.add_argument('-f', '--reference', 
            help='designated reference file for alignments')
    parser.add_argument('--target', 
            help='file of m/z list for targeted extraction')
    parser.add_argument('--database_mode', default='auto',
            help='determines how intermediates are stored, can be "ondisk" or "memory"')
    parser.add_argument('--wlen', default=25, type=int,
            help='determines the number of rt points used when calculating peak prominence')
    parser.add_argument('--max_retention_shift', default=None,
            help='alignment is attempted only using peak pairs differing by this value in seconds or fewer')
    parser.add_argument('--num_lowess_iterations', type=int, default=3,
            help='number of lowess iterations attempted during alignment')
    parser.add_argument('--autoheight', default=False,
            help='automatic determining min peak height')
    parser.add_argument('--min_peak_height', default=False,
            help='minimum height for peaks')
    parser.add_argument('--min_prominence_threshold', default=None,
            help='minimum prominence threshold for peak detection')
    parser.add_argument('--cal_min_peak_height', default=None,
            help='peaks with an intensity below this value are not used for rt calibration')
    parser.add_argument('--min_intensity_threshold', default=None,
            help='signal below this value is removed before peak picking')
    parser.add_argument('--peak_area', default='sum',
            help='peak area calculation, sum, auc or gauss for area under the curve')
    parser.add_argument('--keep_intermediates', default=False, 
            help='keep all intermediate files, ondisk mode only.')
    parser.add_argument('--anno', default=True, 
            help='perform default annotation after processing data')
    parser.add_argument('--debug_rtime_align', default=False, 
            help='Toggle on debug mode for retention alignment: output align figures and reference features.')
    parser.add_argument('--compress', default=False, 
            help='Compress mass tracks to reduce disk usage, default is False')
    parser.add_argument('--drop_unaligned_samples', default=False, 
            help='Drop samples that fail RT alignment from composite map., recommend true for data mining')
    parser.add_argument('--reuse_intermediates', default=None,
            help='Import pickle files for faster processing')
    parser.add_argument('--storage_format', default='pickle',
            help='Storage format for intermediate files, pickle or json')
    parser.add_argument('--single_file_qc_reports', default=False,
            help='Generate a QC report for mzML files during processing')
    parser.add_argument('--spikeins', default=None,
            help='Spike-in standards for QC report - JSON formatted list of lists (name, mz, rt - not checked currently)')
    parser.add_argument('--convert_raw', default=False,
            help='Convert found .raw files to mzML format before processing')
    parser.add_argument('--table_for_viz', default='preferred',
            help='Table to use for visualization, preferred or full')
    parser.add_argument('--vizualization_max_samples', default=20,
            help='Maximum number of samples to display in visualization')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        raise SystemExit

    # update parameters from user specified yaml file
    if args.parameters:
        parameters.update(
            yaml.load(open(args.parameters).read(), Loader=yaml.Loader)
        )
        
    if parameters['multicores'] == 0:
        parameters['multicores'] = mp.cpu_count()
    else:
        parameters['multicores'] = min(mp.cpu_count(), parameters['multicores'])

    parameters['input'] = args.input
    parameters['debug_rtime_align'] = booleandict[args.debug_rtime_align]
    parameters['drop_unaligned_samples'] = booleandict[args.drop_unaligned_samples]
    parameters['autoheight'] = booleandict[args.autoheight]
    
    if args.mode:
        parameters['mode'] = args.mode
    if args.ppm:
        parameters['mz_tolerance_ppm'] = args.ppm
    if args.multicores >= 0:
        if args.multicores == 0:
            print("Using all available cores.")
            parameters['multicores'] = mp.cpu_count()
        else:
            parameters['multicores'] = min(mp.cpu_count(), args.multicores)
    if args.project:
        parameters['project_name'] = args.project
    if args.output:
        parameters['outdir'] = os.path.abspath(args.output)
    if args.peak_area:
        parameters['peak_area'] = args.peak_area
    if args.keep_intermediates:
        parameters['keep_intermediates'] = booleandict[args.keep_intermediates]
        parameters['database_mode'] = 'ondisk'
    if args.anno:
        parameters['anno'] = booleandict[args.anno]
    if args.reference:
        parameters['reference'] = args.reference
    if args.database_mode:
        assert args.database_mode in ['auto', 'ondisk', 'memory'], "Database mode must be either auto, ondisk, or memory."
        parameters['database_mode'] = args.database_mode
    if args.wlen:
        parameters['wlen'] = args.wlen
    if args.compress:
        assert args.compress in ['True', 'False'], "Compress must be either True or False."
        parameters['compress'] = booleandict[args.compress]
        parameters['database_mode'] = 'ondisk'
    if args.storage_format:
        assert args.storage_format in ['pickle', 'json'], "Storage format must be either pickle or json."
        parameters['storage_format'] = args.storage_format
    if args.reuse_intermediates:
        parameters['reuse_intermediates'] = args.reuse_intermediates
        parameters['keep_intermediates'] = True
        parameters['database_mode'] = 'ondisk'
    if args.spikeins:
        parameters['spikeins'] = args.spikeins
    if args.num_lowess_iterations:
        parameters['num_lowess_iterations'] = args.num_lowess_iterations
    if args.table_for_viz:
        parameters['table_for_viz'] = args.table_for_viz

    # update peak detection parameters by autoheight then CLI args
    # min_peak_height, min_prominence_threshold, cal_min_peak_height, min_intensity_threshold
    parameters = update_peak_detection_params(parameters, args)


    # set timestamp here so it can be used in multiple places
    time_stamp = [str(x) for x in time.localtime()[1:6]]
    parameters['time_stamp_for_dir'] = ''.join(time_stamp)
    parameters['time_stamp'] = ':'.join(time_stamp)

    if args.convert_raw:
        convert(parameters, args)
    elif args.run == 'convert':
        convert(parameters, args)
        sys.exit()

    if args.single_file_qc_reports:
        qc_report(parameters, args)
    elif args.run == "qc_report":
        qc_report(parameters, args)
        sys.exit()

    if args.run == 'process':
        process(parameters, args)
    elif args.run == 'analyze':
        # analyze a single sample file to get descriptions
        analyze(parameters, args)
    elif args.run == 'xic':
        # Get XICs (mass tracks) from a folder of centroid mzML files.
        xic(parameters, args)
    elif args.run == 'extract':
        # targeted extraction from a file designated by --target
        extract(parameters, args)
    elif args.run == 'annotate':
        # Annotate a user supplied feature table
        annotate(parameters, args)
    elif args.run == 'join':
        # input a list of directories, each a result of asari process
        join(parameters, args)
    elif args.run == 'qc_report':
        qc_report(parameters, args)
    elif args.run == 'viz':
        # launch data dashboard
        viz(parameters, args)
    elif args.run == 'list_workflows':
        print("Available Worfklows:")
        print("\t1. LC - default option")
        print("\t2. GC, pass `--workflow GC` to enable")
        print("\t3. Lipidomics LC, pass `--workflow Lipidomics` NOT IMPLEMENTED")
    else:
        print("Expecting one of the subcommands: analyze, process, xic, annotate, join, viz, list_workflows.")
#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    main(PARAMETERS)
