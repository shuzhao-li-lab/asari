import argparse
import multiprocessing as mp
import os
import json 
import time
from functools import partial

import yaml

from asari import __version__
from .workflow import (get_mz_list, 
                       process_project, 
                       process_xics, 
                       read_project_dir, 
                       create_export_folders)
from .default_parameters import PARAMETERS
from .dashboard import read_project, dashboard
from .analyze import estimate_min_peak_height, analyze_single_sample
from .annotate_user_table import annotate_user_featuretable
from .utils import build_boolean_dict, bulk_process
from .qc import generate_qc_report
from .gc_annotation import EI_MS_Library

booleandict = build_boolean_dict()
SUBCOMMANDS = ["analyze", "process", "xic", "extract", "annotate", "join", "viz", "list_workflows"]


def convert(parameters, args):
    needs_conversion = []
    for x in os.listdir(args.input):
        if x.endswith('.raw'):
            needs_conversion.append(os.path.abspath(os.path.join(args.input, x)))
    if needs_conversion:
        from .mzml_converter import mzMLconverter
        converter = mzMLconverter(multicores=parameters['multicores'])
        converter.bulk_convert(needs_conversion)

def process(parameters):
    list_input_files = read_project_dir(parameters['input'])
    if not list_input_files:
        print("No valid mzML files are found in the input directory :(")
    else:
        process_project(list_input_files, parameters)

def analyze(parameters, args):
    analyze_single_sample(args.input, parameters=parameters)

def xic(parameters, args):
    list_input_files = read_project_dir(args.input)
    process_xics(list_input_files, parameters)

def extract(parameters, args):
    mzlist = get_mz_list(args.target)
    print("Retrieved %d target mz values from %s.\n" %(len(mzlist), args.target))
    parameters['target'] = mzlist
    process(parameters)

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
        print(f"Table for viz not recognized! Using Full table.")
        dashboard(project_desc, cmap, epd, Ftable)

def qc_report(parameters, args):
    list_input_files = read_project_dir(args.input)
    create_export_folders(parameters)
    jobs = [(f, os.path.join(parameters['qaqc_reports_outdir'], os.path.basename(f).replace(".mzML", "_report.html")), parameters['spikeins']) for f in list_input_files]
    bulk_process(generate_qc_report, jobs)

def update_peak_detection_params(parameters, args=None):
    if parameters['autoheight']:
        try:
            parameters['min_peak_height'] = estimate_min_peak_height(read_project_dir(args.input), parameters)
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

def update_params_from_CLI(parameters, args, debug_print=True):
    def __debug_print(debug_print, to_print):
        if debug_print:
            print(to_print)

    debug_print = partial(__debug_print, debug_print)
    # check if args and parameters are provided
    if args is None:
        raise Exception("No arguments provided.")
    if parameters is None:
        raise Exception("No parameters provided.")
    
    # update parameters from user specified yaml file
    if args.parameters:
        # can be yaml
        try:
            parameters.update(
                yaml.load(open(args.parameters).read(), Loader=yaml.Loader)
            )
        except:
            raise Exception("Failure parsing provided yaml parameters file.")
        # can be json
        try:
            parameters.update(
                json.load(open(args.parameters))
            )
        except:
            raise Exception("Failure parsing provided JSON parameters file.")
        debug_print(to_print=f"Updating default parameters from {args.parameters}")
    else:
        debug_print(to_print=f"Using default parameters")

    # update parameters from command line arguments
        
    # set the ionization mode
    if args.mode:
        assert args.mode in {'pos', 'neg'}, "Mode must be either pos or neg."
        parameters['mode'] = args.mode
        debug_print(to_print=f"Setting mode to {parameters['mode']}")
    else:
        debug_print(to_print=f"Using default mode {parameters['mode']}")
    
    # set the ppm
    if args.ppm:
        assert args.ppm > 0, "PPM must be greater than 0."
        parameters['ppm'] = args.ppm
        debug_print(to_print=f"Setting ppm to {parameters['ppm']}")
    else:
        debug_print(to_print=f"Using default ppm {parameters['mz_tolerance_ppm']}")

    # set the input directory
    if args.input:
        if os.path.isdir(args.input):
            parameters['input'] = args.input
            debug_print(to_print=f"Input determined to be directory")
            debug_print(to_print=f"Setting input to {parameters['input']}")
        elif os.path.isfile(args.input):
            parameters['input'] = os.path.dirname(args.input)
            debug_print(to_print=f"Input determined to be file")
            debug_print(to_print=f"Setting input to {parameters['input']}")
    
    # set the output directory
    if args.output:
        parameters['outdir'] = os.path.abspath(args.output)
        debug_print(to_print=f"Setting outdir to {parameters['outdir']}")
    else:
        parameters['outdir'] = os.path.abspath(os.path.join("./", parameters["outdir"]))
        debug_print(to_print=f"Using default outdir: {parameters['outdir']}")

    # set the project name
    if args.project:
        parameters['project_name'] = args.project
        debug_print(to_print=f"Setting project_name to {parameters['project_name']}")
    else:
        debug_print(to_print=f"Using default project_name: {parameters['project_name']}")

    # set the number of cores to use
    if args.multicores is not None:
        if args.multicores == 0:
            parameters['multicores'] = mp.cpu_count()
            debug_print(to_print=f"Using all available cores to {parameters['multicores']}")
        else:
            parameters['multicores'] = min(mp.cpu_count(), args.multicores)
            debug_print(to_print=f"Setting multicores to {parameters['multicores']}")
    else:
        debug_print(to_print=f"Using default multicores: {parameters['multicores']}")

    # set the reference file
    if args.reference:
        parameters['reference'] = args.reference
        debug_print(to_print=f"Setting reference to {parameters['reference']}")
    else:
        debug_print(to_print=f"Using default reference: {parameters['reference']}")

    # set the target mode
    if args.target:
        assert os.path.isfile(args.target), "Target file must be a valid file."
        parameters['target'] = args.target
        debug_print(to_print=f"Setting target to {parameters['target']}")
    else:
        debug_print(to_print=f"Using default target: {parameters['target']}")

    # set the database mode
    if args.database_mode:
        assert args.database_mode in {'auto', 'ondisk', 'memory'}, "Database mode must be either auto, ondisk, or memory."
        parameters['database_mode'] = args.database_mode
        debug_print(to_print=f"Setting database_mode to {parameters['database_mode']}")
    else:
        debug_print(to_print=f"Using default database_mode: {parameters['database_mode']}")

    # set wlen
    if args.wlen:
        assert args.wlen > 0, "Wlen must be greater than 0."
        parameters['wlen'] = args.wlen
        debug_print(to_print=f"Setting wlen to {parameters['wlen']}")
    else:
        debug_print(to_print=f"Using default wlen: {parameters['wlen']}")

    # set max retention shift
    if args.max_retention_shift:
        assert args.max_retention_shift > 0, "Max retention shift must be greater than 0."
        parameters['max_retention_shift'] = args.max_retention_shift
        debug_print(to_print=f"Setting max_retention_shift to {parameters['max_retention_shift']}")
    else:
        debug_print(to_print=f"Using default max_retention_shift: {parameters['max_retention_shift']}")

    # set num lowess iterations
    if args.num_lowess_iterations:
        assert args.num_lowess_iterations > 0, "Num lowess iterations must be greater than 0."
        parameters['num_lowess_iterations'] = args.num_lowess_iterations
        debug_print(to_print=f"Setting num_lowess_iterations to {parameters['num_lowess_iterations']}")
    else:
        debug_print(to_print=f"Using default num_lowess_iterations: {parameters['num_lowess_iterations']}")
    
    # set autoheight
    if args.autoheight:
        parameters['autoheight'] = booleandict[args.autoheight]
        assert parameters['autoheight'] in {True, False}, "Autoheight must be either True or False."
        debug_print(to_print=f"Setting autoheight to {parameters['autoheight']}")
    else:
        debug_print(to_print=f"Using default autoheight: {parameters['autoheight']}")

    # set min peak height
    if args.min_peak_height:
        parameters['min_peak_height'] = args.min_peak_height
        assert parameters['min_peak_height'] > 0, "Min peak height must be greater than 0."
        debug_print(to_print=f"Setting min_peak_height to {parameters['min_peak_height']}")
    else:
        debug_print(to_print=f"Using default min_peak_height: {parameters['min_peak_height']}")

    # set min prominence threshold
    if args.min_prominence_threshold:
        parameters['min_prominence_threshold'] = args.min_prominence_threshold
        assert parameters['min_prominence_threshold'] > 0, "Min prominence threshold must be greater than 0."
        debug_print(to_print=f"Setting min_prominence_threshold to {parameters['min_prominence_threshold']}")
    else:
        debug_print(to_print=f"Using default min_prominence_threshold of .33 * min_peak_height")

    # set cal min peak height
    if args.cal_min_peak_height:
        parameters['cal_min_peak_height'] = args.cal_min_peak_height
        assert parameters['cal_min_peak_height'] > 0, "Cal min peak height must be greater than 0."
        debug_print(to_print=f"Setting cal_min_peak_height to {parameters['cal_min_peak_height']}")
    else:
        debug_print(to_print=f"Using default cal_min_peak_height: {parameters['cal_min_peak_height']}")

    # set min intensity threshold
    if args.min_intensity_threshold:
        parameters['min_intensity_threshold'] = args.min_intensity_threshold
        assert parameters['min_intensity_threshold'] > 0, "Min intensity threshold must be greater than 0."
        debug_print(to_print=f"Setting min_intensity_threshold to {parameters['min_intensity_threshold']}")
    else:
        debug_print(to_print=f"Using default min_intensity_threshold: {parameters['min_intensity_threshold']}")

    # set peak area
    if args.peak_area:
        assert args.peak_area in {'sum', 'auc', 'gauss'}, "Peak area must be either sum, auc, or gauss."
        parameters['peak_area'] = args.peak_area
        debug_print(to_print=f"Setting peak_area to {parameters['peak_area']}")
    else:
        debug_print(to_print=f"Using default peak_area: {parameters['peak_area']}")

    # set keep intermediates
    if args.keep_intermediates:
        parameters['keep_intermediates'] = booleandict[args.keep_intermediates]
        assert parameters['keep_intermediates'] in {True, False}, "Keep intermediates must be either True or False."
        debug_print(to_print=f"Setting keep_intermediates to {parameters['keep_intermediates']}")
    else:
        debug_print(to_print=f"Using default keep_intermediates: {parameters['keep_intermediates']}")

    # set anno
    if args.anno:
        parameters['anno'] = booleandict[args.anno]
        assert parameters['anno'] in {True, False}, "Anno must be either True or False."
        debug_print(to_print=f"Setting anno to {parameters['anno']}")
    else:
        debug_print(to_print=f"Using default anno: {parameters['anno']}")

    # set debug rtime align
    if args.debug_rtime_align:
        parameters['debug_rtime_align'] = booleandict[args.debug_rtime_align]
        assert parameters['debug_rtime_align'] in {True, False}, "Debug rtime align must be either True or False."
        debug_print(to_print=f"Setting debug_rtime_align to {parameters['debug_rtime_align']}")
    else:
        debug_print(to_print=f"Using default debug_rtime_align: {parameters['debug_rtime_align']}")

    # set compress
    if args.compress:
        parameters['compress'] = booleandict[args.compress]
        assert parameters['compress'] in {True, False}, "Compress must be either True or False."
        debug_print(to_print=f"Setting compress to {parameters['compress']}")
    else:
        debug_print(to_print=f"Using default compress: {parameters['compress']}")

    # set drop unaligned samples
    if args.drop_unaligned_samples:
        parameters['drop_unaligned_samples'] = booleandict[args.drop_unaligned_samples]
        assert parameters['drop_unaligned_samples'] in {True, False}, "Drop unaligned samples must be either True or False."
        debug_print(to_print=f"Setting drop_unaligned_samples to {parameters['drop_unaligned_samples']}")
    else:
        debug_print(to_print=f"Using default drop_unaligned_samples: {parameters['drop_unaligned_samples']}")

    # set reuse intermediates
    if args.reuse_intermediates:
        assert os.path.isdir(args.reuse_intermediates), "Reuse intermediates must be a valid directory."
        parameters['reuse_intermediates'] = args.reuse_intermediates
        debug_print(to_print=f"Setting reuse_intermediates to {parameters['reuse_intermediates']}")
    else:
        debug_print(to_print=f"Using default reuse_intermediates: {parameters['reuse_intermediates']}")

    # set storage format
    if args.storage_format:
        assert args.storage_format in {'pickle', 'json'}, "Storage format must be either pickle or json."
        parameters['storage_format'] = args.storage_format
        debug_print(to_print=f"Setting storage_format to {parameters['storage_format']}")
    else:
        debug_print(to_print=f"Using default storage_format: {parameters['storage_format']}")

    # set single file qc reports
    if args.single_file_qc_reports:
        parameters['single_file_qc_reports'] = booleandict[args.single_file_qc_reports]
        assert parameters['single_file_qc_reports'] in {True, False}, "Single file qc reports must be either True or False."
        debug_print(to_print=f"Setting single_file_qc_reports to {parameters['single_file_qc_reports']}")
    else:
        debug_print(to_print=f"Using default single_file_qc_reports: {parameters['single_file_qc_reports']}")

    # set spikeins
    if args.spikeins:
        assert os.path.isfile(args.spikeins), "Spikeins must be a valid file."
        parameters['spikeins'] = args.spikeins
        debug_print(to_print=f"Setting spikeins to {parameters['spikeins']}")
    else:
        debug_print(to_print=f"Using default spikeins: {parameters['spikeins']}")

    # set convert raw
    if args.convert_raw:
        parameters['convert_raw'] = booleandict[args.convert_raw]
        assert parameters['convert_raw'] in {True, False}, "Convert raw must be either True or False."
        debug_print(to_print=f"Setting convert_raw to {parameters['convert_raw']}")
    else:
        debug_print(to_print=f"Using default convert_raw: {parameters['convert_raw']}")

    # set table for viz
    if args.table_for_viz:
        assert args.table_for_viz in {'preferred', 'full'}, "Table for viz must be either preferred or full."
        parameters['table_for_viz'] = args.table_for_viz
        debug_print(to_print=f"Setting table_for_viz to {parameters['table_for_viz']}")
    else:
        debug_print(to_print=f"Using default table_for_viz: {parameters['table_for_viz']}")
    
    # set visualization max samples
    if args.vizualization_max_samples:
        parameters['vizualization_max_samples'] = args.vizualization_max_samples
        assert parameters['vizualization_max_samples'] > 0, "Visualization max samples must be greater than 0."
        debug_print(to_print=f"Setting vizualization_max_samples to {parameters['vizualization_max_samples']}")
    else:
        debug_print(to_print=f"Using default vizualization_max_samples: {parameters['vizualization_max_samples']}")

    if args.workflow:
        parameters['workflow'] = args.workflow
        assert parameters['workflow'] in {'LC', 'GC', 'LC_start'}, "Workflow must be either LC, GC, or Lipidomics."
        debug_print(to_print=f"Setting workflow to {parameters['workflow']}")
    else:
        debug_print(to_print=f"Using default workflow: {parameters['workflow']}")

    # set retention index standards
    if args.retention_index_standards:
        assert os.path.isfile(args.retention_index_standards), "Retention index standards must be a valid file."
        parameters['retention_index_standards'] = args.retention_index_standards
        debug_print(to_print=f"Setting retention_index_standards to {parameters['retention_index_standards']}")
    else:
        debug_print(to_print=f"Using default retention_index_standards: {parameters['retention_index_standards']}")

    if args.GC_Database_Manifest:
        parameters['GC_Database_Manifest'] = args.GC_Database_Manifest
        debug_print(to_print=f"Setting GC_Database_Manifest to {parameters['GC_Database_Manifest']}")
    else:
        parameters['GC_Database_Manifest'] = None
        debug_print(to_print=f"Using default GC_Database_Manifest.")

    if args.GC_Database:
        parameters['GC_Database'] = args.GC_Database
        EI_MS_Library.load_library_manifest()
        debug_print(to_print=f"Setting GC_Database to {parameters['GC_Database']}")
    else:
        debug_print(to_print=f"Using default GC_Database: {parameters['GC_Database']}")

    print(args.run)
    if args.run:
        parameters['run'] = args.run.rstrip()
        debug_print(to_print=f"Setting run to {parameters['run']}")


def initialize_parameters(parameters, args):
    parameters['asari_version'] = __version__
    parameters['timestamp'] = time.strftime("%Y%m%d-%H%M%S")

def build_parser():
    # Since the default_parameters is used by default, we should not set default 
    # values in the CLI. This is good practice to avoid confusion or misconfiguration.
    # The CLI takes priority over an optional parameter file which takes priority over
    # any of the default_parameters.

    parser = argparse.ArgumentParser(description='asari, LC-MS metabolomics data preprocessing')

    parser.add_argument('-v', '--version', action='version', version=__version__, 
            help='print version and exit')
    parser.add_argument('run', metavar='subcommand', type=str,
            help="one of the subcommands: " + ", ".join(SUBCOMMANDS))
    parser.add_argument('-m', '--mode', type=str,
            help='mode of ionization, pos or neg')
    parser.add_argument('--ppm', type=int, 
            help='mass precision in ppm (part per million), same as mz_tolerance_ppm')
    parser.add_argument('-i', '--input', type=str,
            help='input directory of mzML files to process, or a single file to analyze')
    parser.add_argument('-o', '--output', type=str,
            help='output directory')
    parser.add_argument('-j', '--project', type=str,
            help='project name')
    parser.add_argument('-p', '--parameters', type=str, 
            help='Custom paramter file in YAML. Use parameters.yaml as template.')
    parser.add_argument('-c', '--multicores', type=int,
            help='nunmber of CPU cores intented to use')
    parser.add_argument('-f', '--reference', type=str,
            help='designated reference file for alignments')
    parser.add_argument('--target', type=str,
            help='file of m/z list for targeted extraction')
    parser.add_argument('--database_mode', type=str,
            help='determines how intermediates are stored, can be "ondisk" or "memory"')
    parser.add_argument('--wlen', type=int,
            help='determines the number of rt points used when calculating peak prominence')
    parser.add_argument('--max_retention_shift', default=None,
            help='alignment is attempted only using peak pairs differing by this value in seconds or fewer')
    parser.add_argument('--num_lowess_iterations', type=int,
            help='number of lowess iterations attempted during alignment')
    parser.add_argument('--autoheight',
            help='automatic determining min peak height')
    parser.add_argument('--min_peak_height', type=int,
            help='minimum height for peaks')
    parser.add_argument('--min_prominence_threshold', type=int,
            help='minimum prominence threshold for peak detection')
    parser.add_argument('--cal_min_peak_height', type=int,
            help='peaks with an intensity below this value are not used for rt calibration')
    parser.add_argument('--min_intensity_threshold', type=int,
            help='signal below this value is removed before peak picking')
    parser.add_argument('--peak_area', type=str,
            help='peak area calculation, sum, auc or gauss for area under the curve')
    parser.add_argument('--keep_intermediates', 
            help='keep all intermediate files, ondisk mode only.')
    parser.add_argument('--anno',
            help='perform default annotation after processing data')
    parser.add_argument('--debug_rtime_align', 
            help='Toggle on debug mode for retention alignment: output align figures and reference features.')
    parser.add_argument('--compress', 
            help='Compress mass tracks to reduce disk usage, default is False')
    parser.add_argument('--drop_unaligned_samples', 
            help='Drop samples that fail RT alignment from composite map., recommend true for data mining')
    parser.add_argument('--reuse_intermediates', 
            help='Import pickle files for faster processing')
    parser.add_argument('--storage_format', type=str,
            help='Storage format for intermediate files, pickle or json')
    parser.add_argument('--single_file_qc_reports',
            help='Generate a QC report for mzML files during processing')
    parser.add_argument('--spikeins', type=str,
            help='Spike-in standards for QC report - JSON formatted list of lists (name, mz, rt - not checked currently)')
    parser.add_argument('--convert_raw',
            help='Convert found .raw files to mzML format before processing')
    parser.add_argument('--table_for_viz', type=str,
            help='Table to use for visualization, preferred or full')
    parser.add_argument('--vizualization_max_samples', type=int,
            help='Maximum number of samples to display in visualization')
    parser.add_argument('--workflow', type=str,
            help='Workflow to use, LC by default')
    parser.add_argument('--retention_index_standards', type=str,
            help='Path to retention index standards, needed for GC workflow')
    parser.add_argument('--GC_Database', type=str, 
            help='Path to GC database, or GCMS database name for retrieval')
    parser.add_argument('--GC_Database_Manifest', type=str,
            help='Path to GC database manifest file')
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        raise SystemExit
    return args

def run_asari(parameters, args=None):
    if 'run_gui' in parameters:
        if parameters['run_gui'] == 'process':
            process(parameters)
            exit()

    if parameters['run'] == 'process':
        # these can be done before processing
        if args.convert_raw:
            convert(parameters, args)
        if args.single_file_qc_reports:
            qc_report(parameters, args)
        process(parameters)
    elif parameters['run'] == 'convert':
        convert(parameters, args)
    elif parameters['run']== "qc_report":
        qc_report(parameters, args)
    elif parameters['run'] == 'analyze':
        # analyze a single sample file to get descriptions
        analyze(parameters, args)
    elif parameters['run'] == 'xic':
        # Get XICs (mass tracks) from a folder of centroid mzML files.
        xic(parameters, args)
    elif parameters['run'] == 'extract':
        # targeted extraction from a file designated by --target
        extract(parameters, args)
    elif parameters['run'] == 'annotate':
        # Annotate a user supplied feature table
        annotate(parameters, args)
    elif parameters['run'] == 'join':
        # input a list of directories, each a result of asari process
        join(parameters, args)
    elif parameters['run'] == 'qc_report':
        qc_report(parameters, args)
    elif parameters['run'] == 'viz':
        # launch data dashboard
        viz(parameters, args)
    elif parameters['run'] == 'list_workflows':
        print("Available Worfklows:")
        print("\t1. LC - default option")
        print("\t2. GC, pass `--workflow GC` to enable")
        print("\t3. Lipidomics LC, pass `--workflow Lipidomics` NOT IMPLEMENTED")
    else:
        print("Expecting one of the subcommands: analyze, process, xic, annotate, join, viz, list_workflows.")

def main():
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

    # make a copy of the default parameters
    parameters = PARAMETERS

    # build CLI parser
    args = build_parser()

    # set time stamp, version, etc. in parameters
    initialize_parameters(parameters, args)

    # overwrite parameters with user provided arguments from CLI
    update_params_from_CLI(parameters, args)

    # update peak detection parameters
    # min_peak_height, min_prominence_threshold, cal_min_peak_height, min_intensity_threshold
    update_peak_detection_params(parameters, args)
    run_asari(parameters, args)


#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    main()
