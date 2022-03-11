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
import sys
import time
from .workflow import *

#
# important parameters to estimate: ppm of chrom construction; min_peak_height
#

PARAMETERS = {
    'project_name': 'test_asari',
    'outdir': 'output',
    'database_mode': 'ondisk',          # 'ondisk', 'memory' (run in memory, only small studies), 
                                        # 'mongo' (MongoDB, requiring installation of DB server)

    'multicores': 4,                    # number of threads allowed in parallel processing
    'init_samples_number': 3,           # initiation samples, kept in memory
    'initiation_samples': [],           # if user to specify N samples to initiate data processing, 
                                        # otherwise they are chosen automatically

    'rt_align': True,                   # False to bypass retention time alignment
    
    'min_intensity_threshold': 100,     # minimal intensity for mass track extraction, filtering baseline
    'min_peak_height': 10000,           # minimal peak height

    'cal_min_peak_height': 100000,      # minimal peak height required for peaks used for RT calibration
    'min_timepoints': 6,                # minimal number of data points in elution profile. scipy find_peaks treat `width` as FWHM, thus half of this value.
    'signal_noise_ratio': 5,            # peak height at least x fold over noise, which is median of non-peak data points.
    #
    'peak_number_rt_calibration': 15,   # minimal number of selected high-quality peaks required for RT calibration. 
                                        # Samples with fewer selected peaks are dropped out.

    # tweaking this section now -
    'cache_mass_traces': False,         # to save memory if not using DB; turn on if need to plot and diagnose
    'output_feature_table': 'Feature_table.tsv',
    'mass_grid_mapping': "_mass_grid_mapping.csv",
    'annotation_filename': "Annotation_table.tsv",
    'json_empricalCompounds': "_empCpd_json.json",

    #
    'mode': 'pos',                      # ionization mode
    'mass_range': (50, 2000),
    #'max_rtime': 300,                   # retention time range (chromatography) in seconds, to auto populate

    'mz_tolerance': 5,                  # ppm, high selectivity meaning no overlap neighbors to confuse; 
                                        # Low selectivity regions will be still inspected to determine the true number of features
    'rtime_tolerance': 10,              # feature rtime shift threshold under 10 seconds; or 10% of rtime
                                        # not used now?
    #
    
    # no need to modify below unless you know what you are doing
    'prominence_window': 30,
    'gaussian_shape': 0.3,              # min cutoff
    }
    
PARAMETERS['min_prominence_threshold'] = PARAMETERS['min_peak_height']/3.0
PARAMETERS['multicores'] = min(mp.cpu_count(), PARAMETERS['multicores'])


def read_project_dir(directory, file_pattern='.mzML'):
    '''
    This reads centroided LC-MS files.
    For OpenMS based XIC workflow, file_pattern='chrom.mzML'.
    '''
    print("Working on ~~ %s ~~ \n\n" %directory)
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

def process_project(list_input_files, dict_meta_data={}, parameters=PARAMETERS):
    '''
    list_input_files: Extracted ion chromatogram files.
    parameters: dictionary of most parameters.

    '''
    if dict_meta_data:
        for k in dict_meta_data:
            dict_meta_data[k] = dict_meta_data[k].upper()       # upper characters to standardize
    if not list_input_files:
        print("No input file found. Please verify your pathway to files.")

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


def main(directory, parameters=PARAMETERS):
    # time_stamp is `month daay hour minute second``
    time_stamp = ''.join([str(x) for x in time.localtime()[1:6]])
    if parameters['database_mode'] == 'ondisk':
        parameters['outdir'] = '_'.join([parameters['project_name'], parameters['outdir'], time_stamp]) 
        os.mkdir(parameters['outdir'])
        os.mkdir(os.path.join(parameters['outdir'], 'pickle'))
        os.mkdir(os.path.join(parameters['outdir'], 'export'))


    print("\n\n~~~~~~~ Hello from Asari! ~~~~~~~~~\n")
    process_project(
            read_project_dir(directory), 
            {},         # not used now
            parameters, 
               #setting output_dir as input dir
    )

#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    PARAMETERS['mode'] = sys.argv[1]
    directory = sys.argv[2]
    main(directory)
