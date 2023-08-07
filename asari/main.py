import argparse
from yaml import load, Loader
import multiprocessing as mp
import time 

from asari import __version__
from .workflow import (get_mz_list, 
                       process_project, 
                       process_xics, 
                       read_project_dir)
from .default_parameters import PARAMETERS
from .dashboard import read_project, dashboard
from .analyze import estimate_min_peak_height, analyze_single_sample
from .annotate_user_table import annotate_user_featuretable

booleandict = {
    'T': True, 
    'F': False, 
    1: True, 
    0: False, 
    'True': True, 
    'False': False, 
    'TRUE': True, 
    'FALSE': False, 
    'true': True, 
    'false': False
    }

PARAMETERS['asari_version'] = __version__

def __run_process__(parameters, args):
    
    # main process function
    list_input_files = read_project_dir(args.input)
    if not list_input_files:
        print("No valid mzML files are found in the input directory :(")
    else:
        if args.autoheight:
            try:
                parameters['min_peak_height'] = estimate_min_peak_height(list_input_files)
                if parameters['cal_autoheight_mode'] == 'multiplier':
                    parameters['cal_min_peak_height'] = parameters['cal_autoheight_multiplier'] * parameters['min_peak_height']
            except ValueError as err:
                print("Problems with input files: {0}. Back to default min_peak_height.".format(err))
        parameters['min_prominence_threshold'] = int( 0.33 * parameters['min_peak_height'] )
        process_project(list_input_files, parameters )

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
    import os
    import pickle
    import pandas as pd
    import uuid
    from . import workflow
    from . import constructors
    from . import experiment

    
    projects = []
    for x in parameters['to_join']:
        p_desc, cmap, epd, ftable, p_map = read_project(x)
        exp = pickle.load(open(os.path.join(x, "export/experiment.json"), 'rb'))
        projects.append({"experiment": exp, 
                         "cmap": cmap, 
                         "project": p_desc,
                         "emp_cpd": epd,
                         "ftable": ftable,
                         "mapping": p_map,
                         "input_file": x
                        })
    workflow.create_export_folders(parameters)
    parameters['min_peak_height'] = min([x["experiment"].parameters["min_peak_height"] for x in projects])
    parameters['min_prominence_threshold'] = min([x["experiment"].parameters["min_prominence_threshold"] for x in projects])
    parameters['database_mode'] = 'memory'

    
    sample_registry = {}
    mock_extracts = {}
    for ii, x in enumerate(projects):
        sample_registry[ii] = {
            "sample_id": ii,
            "input_file": x['input_file'],
            "name": uuid.uuid1()
        }
        dict_scan_rtime = x["cmap"]["dict_scan_rtime"]
        list_scan_times = list(dict_scan_rtime.keys())
        list_mass_tracks = list(x["cmap"]["list_mass_tracks"].values())
        anchor_mz_pairs = workflow.find_mzdiff_pairs_from_masstracks(list_mass_tracks)
        print("EXP: ", x['experiment'].all_samples)
        mock_extracts[ii] = {
            'status:mzml_parsing': 'passed',
            'status:eic': 'passed',
            'data_location': None,
            'max_scan_number': max(list_scan_times),
            'list_scan_numbers': list_scan_times,
            'list_retention_time': list(dict_scan_rtime.values()),
            'track_mzs': [(z['mz'], z['id_number']) for z in list_mass_tracks],
            'anchor_mz_pairs': anchor_mz_pairs,
            'number_anchor_mz_pairs': len(anchor_mz_pairs),
            'composite_of': x['experiment'],
            'sample_mapping': x['mapping'],
            'sample_data': {
                'sample_id': ii,
                'input_file': None,
                'ion_mode': x['experiment'].mode,
                'max_scan_number': max(list_scan_times),
                'list_mass_tracks': list_mass_tracks, 
                'anchor_mz_pairs': anchor_mz_pairs,
                'number_anchor_mz_pairs': len(anchor_mz_pairs)
            }
        }
        sample_registry[ii].update(mock_extracts[ii])

    meta_exp = experiment.ext_Experiment(sample_registry, parameters)
    meta_exp.process_all()
    #meta_exp.export_all()
    meta_ft = []
    for ii, peak in enumerate(meta_exp.CMAP.FeatureList):
        row = dict(peak)
        for meta_sample in meta_exp.all_samples:
            meta_l_base = meta_sample.reverse_rt_cal_dict.get(peak['left_base'], peak['left_base'])
            meta_r_base = meta_sample.reverse_rt_cal_dict.get(peak['right_base'], peak['right_base'])
            meta_track_number = meta_exp.CMAP.MassGrid[meta_sample.name][peak['parent_masstrack_id']]
            if not pd.isna(meta_track_number):
                meta_track_number = int(meta_track_number)
                for sample in meta_sample.composite_of.all_samples:
                    sample_track_number = meta_sample.sample_mapping[sample.name][meta_track_number]
                    if not pd.isna(sample_track_number):
                        sample_mass_tracks = sample.get_masstracks_and_anchors()
                        sample_mass_track = sample_mass_tracks[int(sample_track_number)]
                        sample_left_base = sample.reverse_rt_cal_dict.get(meta_l_base, meta_l_base)
                        sample_right_base = sample.reverse_rt_cal_dict.get(meta_r_base, meta_r_base)
                        sample_peak_area = sample_mass_track['intensity'][sample_left_base: sample_right_base+1].sum()
                        row[sample.name] = sample_peak_area
                    else:
                        row[sample.name] = None
        meta_ft.append(row)
    meta_ft = pd.DataFrame(meta_ft)
    meta_ft.to_csv('./new_ftable.tsv', sep="\t")

def viz(parameters, args):
    datadir = args.input
    project_desc, cmap, epd, Ftable = read_project(datadir)
    dashboard(project_desc, cmap, epd, Ftable)

def main(parameters=PARAMETERS):
    '''
    asari, Trackable and scalable Python program for high-resolution LC-MS metabolomics data preprocessing.

        * analyze: analyze a single mzML file to print summary of statistics and recommended parameters.
        * process: LC-MS data preprocessing
        * xic: construct mass trakcs (chromatogram) from mzML files
        * extract: targeted extraction of given m/z list
        * annotate: annotate a list of features
        * join: merge multiple processed projects (possibly split a large dataset)
        * viz: start interactive data visualization and exploration.

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
    parser.add_argument('--to_join', '--list', nargs='+',
            help='experiments to join')
    parser.add_argument('-i', '--input', 
            help='input directory of mzML files to process, or a single file to analyze')
    parser.add_argument('-i2', '--input2', 
                        help='input_file_2')
    parser.add_argument('-o', '--output', 
            help='output directory')
    parser.add_argument('-j', '--project', 
            help='project name')
    parser.add_argument('-p', '--parameters', 
            help='Custom paramter file in YAML. Use parameters.yaml as template.')
    parser.add_argument('-c', '--cores', type=int, 
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
    parser.add_argument('--num_lowess_iterations', type=int, default=1,
            help='number of lowess iterations attempted during alignment')
    parser.add_argument('--autoheight', default=False,
            help='automatic determining min peak height')
    parser.add_argument('--peak_area', default='sum',
            help='peak area culculation, sum, auc or gauss for area under the curve')
    parser.add_argument('--pickle', default=False, 
            help='keep all intermediate pickle files, ondisk mode only.')
    parser.add_argument('--anno', default=True, 
            help='perform default annotation after processing data')
    parser.add_argument('--debug_rtime_align', default=False, 
            help='Toggle on debug mode for retention alignment: output align figures and reference features.')

    args = parser.parse_args()

    # update parameters
    if args.parameters:
        parameters.update(
            load(open(args.parameters).read(), Loader=Loader)
        )
    parameters['multicores'] = min(mp.cpu_count(), parameters['multicores'])
    parameters['input'] = args.input
    parameters['input2'] = args.input2
    parameters['to_join'] = args.to_join
    parameters['debug_rtime_align'] = booleandict[args.debug_rtime_align]
    parameters['time_stamp'] = '_'.join([str(x) for x in time.localtime()[1:6]])


    if args.mode:
        parameters['mode'] = args.mode
    if args.ppm:
        parameters['mz_tolerance_ppm'] = args.ppm
    if args.cores:
        parameters['multicores'] = min(mp.cpu_count(), args.cores)
    if args.project:
        parameters['project_name'] = args.project
    if args.output:
        parameters['outdir'] = args.output
    if args.peak_area:
        parameters['peak_area'] = args.peak_area
    if args.pickle:
        parameters['pickle'] = booleandict[args.pickle]
    if args.anno:
        parameters['anno'] = booleandict[args.anno]
    if args.reference:
        parameters['reference'] = args.reference
    if args.database_mode:
        parameters['database_mode'] = args.database_mode
    if args.wlen:
        parameters['wlen'] = args.wlen
    if args.max_retention_shift:
        parameters['max_retention_shift'] = float(args.max_retention_shift)
    if args.num_lowess_iterations:
        parameters['num_lowess_iterations'] = args.num_lowess_iterations

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
        # input a list of directoreis, each a result of asari process
        join(parameters, args)
    elif args.run == 'viz':
        # launch data dashboard
        viz(parameters, args)
    else:
        print("Expecting one of the subcommands: analyze, process, xic, annotate, join, viz.")

#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    main(PARAMETERS)
