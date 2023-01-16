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
    'database_mode': 'auto',            # `auto` determined on sample number 
                                        # 'ondisk', 'memory' (run in memory, only small studies), 
                                        # 'mongo' (MongoDB, requiring installation of DB server, to implement)
    'multicores': 4,                    # number of cores allowed in parallel processing

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
    'wlen': 25,                         # window size for evaluating prominence in peaks. Important to resolve clustered narrow peaks.
    'reverse_detection': False,         # use the reverse peak detection algorithm
    'autoheight': False,                # min_peak_height can be estimated automatically by setting autoheight on in CLI  
    'gaussian_shape': 0.5,              # min cutoff of goodness of fitting to Gauss model
    
    # retention time alignment
    'reference': None,
    'rt_align_method': 'lowess',        # 'lowess', 'tolerance', or to implement           
    'rt_align_on': True,                # False to bypass retention time alignment
    'rtime_tolerance': 50,              # feature rtime shift threshold under 10 seconds; or 10% of rtime   
    'cal_min_peak_height': 100000,      # minimal peak height required for peaks used for RT calibration
    'peak_number_rt_calibration': 15,   # minimal number of selected high-quality peaks required for RT calibration. 
                                        # Samples with fewer selected peaks are dropped out.

    # Number of samples dictates workflow 
    'project_sample_number_small': 10,  # 10
    
    # default output names
    'output_feature_table': 'Feature_table.tsv',
    'mass_grid_mapping': "_mass_grid_mapping.csv",
    'annotation_filename': "Annotation_table.tsv",
    'json_empricalCompounds': "_empCpd_json.json",
    # for annotation
    'check_isotope_ratio': False,


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
    