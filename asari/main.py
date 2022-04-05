import argparse

from asari import __version__
from .workflow import *
from .defaul_parameters import PARAMETERS


def main(parameters=PARAMETERS):
    '''
    asari, Trackable and scalable Python program for high-resolution LC-MS metabolomics data preprocessing.
        analyze: analyze a single mzML file to print summary of statistics and recommended parameters.
        process: LC-MS data preprocessing
        xic: construct mass trakcs (chromatogram) from mzML files
        annotate: annotate a list of features
        join: merge multiple processed projects (possibly split a large dataset)
        viz: start interactive data visualization and exploration.
    '''
    print("\n\n~~~~~~~ Hello from Asari (%s) ~~~~~~~~~\n" %__version__)
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
    parser.add_argument('-c', '--cores', type=int, 
            help='nunmber of CPU cores intented to use')

    parser.add_argument('--autoheight', default=False,
            help='automatic determining min peak height')

    args = parser.parse_args()

    # update parameters
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
    
    if args.run == 'process':
        list_input_files = read_project_dir(args.input)
        if args.autoheight:
            parameters['min_peak_height'] = estimate_min_peak_height(list_input_files)
        parameters['min_prominence_threshold'] = int( 0.33 * parameters['min_peak_height'] )

        process_project( list_input_files,  parameters )

    elif args.run == 'analyze':
        # analyze a single sample file to get descriptions
        from .analyze import analyze_single_sample
        analyze_single_sample(args.input, parameters=parameters)

    elif args.run == 'xic':
        # 
        pass


    elif args.run == 'target':
        # targeted extraction
        pass


    elif args.run == 'annotate':
        # 
        from .annotate_user_table import annotate_user_featuretable
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

    main(PARAMETERS)
