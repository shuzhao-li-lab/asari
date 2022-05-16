import argparse
from yaml import load, Loader

from asari import __version__
from .workflow import *
from .defaul_parameters import PARAMETERS


def main(parameters=PARAMETERS):
    '''
    asari, Trackable and scalable Python program for high-resolution LC-MS metabolomics data preprocessing.

        analyze: analyze a single mzML file to print summary of statistics and recommended parameters.
        process: LC-MS data preprocessing
        xic: construct mass trakcs (chromatogram) from mzML files
        extract: targeted extraction of given m/z list
        annotate: annotate a list of features
        join: merge multiple processed projects (possibly split a large dataset)
        viz: start interactive data visualization and exploration.
    '''
    def __run_process__(parameters, args):
        # main process function
        list_input_files = read_project_dir(args.input)
        if args.autoheight:
            from .analyze import estimate_min_peak_height
            try:
                parameters['min_peak_height'] = estimate_min_peak_height(list_input_files)
            except ValueError as err:
                print("Problems with input files: {0}. Back to default min_peak_height.".format(err))

        parameters['min_prominence_threshold'] = int( 0.33 * parameters['min_peak_height'] )
        process_project( list_input_files,  parameters )

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
    parser.add_argument('-c', '--cores', type=int, 
            help='nunmber of CPU cores intented to use')
    parser.add_argument('-f', '--reference', 
            help='designated reference file for alignments')
    parser.add_argument('--target', 
            help='file of m/z list for targeted extraction')
    parser.add_argument('--autoheight', default=False,
            help='automatic determining min peak height')

    args = parser.parse_args()

    # update parameters
    if args.parameters:
        parameters.update(
            load(open(args.parameters).read(), Loader=Loader)
        )
    parameters['multicores'] = min(mp.cpu_count(), parameters['multicores'])
    parameters['input'] = args.input
    
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
    if args.reference:
        parameters['reference'] = args.reference

    if args.run == 'process':
        __run_process__(parameters, args)

    elif args.run == 'analyze':
        # analyze a single sample file to get descriptions
        from .analyze import analyze_single_sample
        analyze_single_sample(args.input, parameters=parameters)

    elif args.run == 'xic':
        # Get XICs (mass tracks) from a folder of centroid mzML files.
        list_input_files = read_project_dir(args.input)
        process_xics(list_input_files, parameters)

    elif args.run == 'extract':
        # targeted extraction from a file designated by --target
        mzlist = get_mz_list(args.target)
        print("Retrieved %d target mz values from %s.\n" %(len(mzlist), args.target))
        parameters['target'] = mzlist
        __run_process__(parameters, args)

    elif args.run == 'annotate':
        # Annotate a user supplied feature table
        from .annotate_user_table import annotate_user_featuretable
        annotate_user_featuretable(args.input, parameters=parameters)

    elif args.run == 'join':
        # input a list of directoreis, each a result of asari process
        pass

    elif args.run == 'viz':
        # launch data dashboard
        datadir = args.input
        from .dashboard import read_project, dashboard
        project_desc, cmap, epd, Ftable = read_project(datadir)
        dashboard(project_desc, cmap, epd, Ftable)

    else:
        print("Expecting one of the subcommands: analyze, process, xic, annotate, join, viz.")

#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    PARAMETERS['asari_version'] = __version__
    main(PARAMETERS)
