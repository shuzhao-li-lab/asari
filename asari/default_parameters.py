# Here are default values of parameters. Please do NOT modify unless necesssary.
# 
# One can specify parameters via 
# 1) custom parameters.yaml
# 2) commandline arguments
# Priority is: commandline overwriting parameters.yaml overwriting this file.
#
# Default values of others should be acceptable for most studies. Really, 
# only ionization mode (default pos) and ppm of mz_tolerance (defulat 5) may need attention, 
# but they should be given via commandline arguments.
# 

PARAMETERS = {
    'project_name': 'asari_project',
    'outdir': 'output',
    'pickle': False,
    'database_mode': 'memory',            # `auto` determined on sample number 
                                        # 'ondisk', 'memory' (run in memory, only small studies), 
                                        # 'mongo' (MongoDB, requiring installation of DB server, to implement)
    
    'multicores': 8,                    # number of cores allowed in parallel processing
    'write_cache': 4 * 1024 ** 3,       # amount of memory in bytes used for mass track cache
    'disk_cache': 0 * 1024 ** 3,      # amount of disk space allowed to be used before compressing stored mass tracks
    'sparsify': False,                  # use sparse representations of mass tracks in memory for space saving                

    # mass parameters
    'mode': 'pos',                      # ionization mode
    'mass_range': (50, 2000),
    'mz_tolerance_ppm': 5,              # ppm, high selectivity meaning no overlap neighbors to confuse; 
                                        # Low selectivity regions will be still inspected to determine the true number of features

    # chromatogram and peak parameters
    'min_timepoints': 6,                # minimal number of data points in elution profile. scipy find_peaks treat `width` as FWHM, thus half of this value.
    'signal_noise_ratio': 2,            # peak height at least x fold over local noise
    'min_intensity_threshold': 1000,    # minimal intensity for mass track extraction, filtering baseline
    'min_peak_height': 100000,           # minimal peak height.
    'min_peak_ratio': 0.001,             # minimal ratio of a peak of the max height of its ROI, relevant to small peaks next to big ones.
    'wlen': 25,                         # window size for evaluating prominence in peaks. Important to resolve clustered narrow peaks.
    'autoheight': False,                # min_peak_height can be estimated automatically by setting autoheight on in CLI  
    'gaussian_shape': 0.5,              # min cutoff of goodness of fitting to Gauss model
    'peak_area': 'sum',                 # `sum` for simple sum, `auc` for area under the curve, `gauss` for gaussian
    
    # retention time alignment
    'reference': None,
    'rt_align_method': 'lowess',        # 'lowess', 'tolerance', or to implement           
    'rt_align_on': True,                # False to bypass retention time alignment
    'drop_unaligned_samples': False,    # Drop samples that fail RT alignment from composite map
    'rtime_tolerance': 50,              # feature rtime shift threshold under 10 seconds; or 10% of rtime   
    'cal_min_peak_height': 100000,      # minimal peak height required for peaks used for RT calibration
    'peak_number_rt_calibration': 20,   # minimal number of selected high-quality peaks required for RT calibration. 
                                        # Samples with fewer selected peaks are dropped out.
    'max_retention_shift': None,        # landmark peak pairs with a scan number delta greater than this are not used for RT calibration
    'num_lowess_iterations': 3,         # number of lowess iterations to perform for RT calibration, higher values take longer but less sensitive to outliers
    # Number of samples dictates workflow 
    'project_sample_number_small': 0,  # 10
    
    # default output names
    'output_feature_table': 'Feature_table.tsv',
    'mass_grid_mapping': "_mass_grid_mapping.csv",
    'annotation_filename': "Annotation_table.tsv",
    'json_empricalCompounds': "_empCpd_json.json",
    # for annotation
    'check_isotope_ratio': False,
    'workflow': "LC",


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 'max_rtime': 300,                   # retention time range (chromatography) in seconds, to auto populate
    # 'prominence_window': 101,           # not used now; no need to modify  
    # 'cache_mass_traces': False,         # to save memory if not using DB; turn on if need to plot and diagnose
    # 'init_samples_number': 3,           # initiation samples,
    # 'initiation_samples': [],           # if user to specify N samples to initiate data processing, 
                                          # otherwise they are chosen automatically
    # 'project_sample_number_large': 1000, # not used now
    # 
    }
    

#
# From khipu.utils. Placed here so that people can customize within asari.
#
PROTON = 1.00727646677
electron = 0.000549

# avoid confusing adducts in initial search, e.g. H, H2O
adduct_search_patterns = [  # initial patterns are relative to M+H+
                            (21.9820, 'Na/H'),
                            (41.026549, 'ACN'),     # Acetonitrile
                            (35.9767, 'HCl'),
                            (37.955882, 'K/H'),
                            ]

adduct_search_patterns_neg = [ (35.9767, 'HCl'), 
                            (46.00548, 'HCOOH'),
                            (21.9820, 'Na/H'), 
                            (41.026549, 'ACN'),
                            (37.955882, 'K/H'),
                            ]

isotope_search_patterns = [ (1.003355, '13C/12C', (0, 0.8)),
                            (2.00671, '13C/12C*2', (0, 0.8)),
                            # (3.010065, '13C/12C*3', (0, 0.8)),
                            # (4.01342, '13C/12C*4', (0, 0.8)),
                            # (5.016775, '13C/12C*5', (0, 0.8)),
                            # (6.02013, '13C/12C*6', (0, 0.8)),
                            # (7.023485, '13C/12C*7', (0, 0.8)),
                            # (8.02684, '13C/12C*8', (0, 0.8)),
                            # (9.030195, '13C/12C*9', (0, 0.8)),
                            # (10.03355, '13C/12C*10', (0, 0.8)),
                            # (11.036905, '13C/12C*11', (0, 0.8)),
                            # (12.04026, '13C/12C*12', (0, 0.8)),
                            ]

extended_adducts = [(1.0078, 'H'),
                            (-1.0078, '-H'),
                            (10.991, 'Na/H, double charged'),
                            (0.5017, '13C/12C, double charged'),
                            (117.02655, '-NH3'),
                            (17.02655, 'NH3'),
                            (-18.0106, '-H2O'),
                            (18.0106, 'H2O'),      # easy to confuse with bio reactions
                            (18.033823, 'NH4'),
                            (27.01089904, 'HCN'),
                            (37.94694, 'Ca/H2'),
                            (32.026215, 'MeOH'),
                            (43.96389, 'Na2/H2'),
                            (67.987424, 'NaCOOH'),
                            (83.961361, 'KCOOH'),
                            (97.96737927, 'H2SO4'),
                            (97.97689507, 'H3PO4'),
]


# to include in export for end users
readme_doc_str = """
The recommended feature table is `preferred_Feature_table.tsv`. 

All peaks are kept in `export/full_Feature_table.tsv` 
if they meet signal (snr) and shape standards 
(part of input parameters but default values are fine for most people). 
The filtering decisions are left to end users.

Annotation is in JSON (`Annotated_empricalCompounds.json`) 
and in tab delimited text (`Feature_annotation.tsv`).

The processing parameters and history are in `project.json`.

Please refer to https://github.com/shuzhao-li/asari for details, 
report bugs or request features.
"""
