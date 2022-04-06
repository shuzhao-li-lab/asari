# Attention is required for ionization mode and ppm of mz_tolerance.
# Default values of others should be acceptable for most studies.

PARAMETERS = {
    'project_name': 'asari_project',
    'outdir': 'output',
    'database_mode': 'auto',            # `auto` determined on sample number 
                                        # 'ondisk', 'memory' (run in memory, only small studies), 
                                        # 'mongo' (MongoDB, requiring installation of DB server, to implement)
    'multicores': 4,                    # number of cores allowed in parallel processing

    # mass parameters
    'mode': 'pos',                      # ionization mode
    'mass_range': (50, 2000),
    'mz_tolerance': 5,                  # ppm, high selectivity meaning no overlap neighbors to confuse; 
                                        # Low selectivity regions will be still inspected to determine the true number of features

    # chromatogram and peak parameters
    'min_timepoints': 6,                # minimal number of data points in elution profile. scipy find_peaks treat `width` as FWHM, thus half of this value.
    'signal_noise_ratio': 5,            # peak height at least x fold over noise, which is median of non-peak data points.
    'min_intensity_threshold': 1000,    # minimal intensity for mass track extraction, filtering baseline
    'min_peak_height': 10000,           # minimal peak height.
                                        # min_peak_height can be estimated automatically by setting autoheight on in CLI
    'prominence_window': 30,            # no need to modify  
    'gaussian_shape': 0.3,              # min cutoff
    
    # retention time alignment
    'rt_align_method': 'lowess',        # 'lowess', 'tolerance', or to implement           
    'rt_align_on': True,                # False to bypass retention time alignment
    'rtime_tolerance': 50,              # feature rtime shift threshold under 10 seconds; or 10% of rtime   
    'cal_min_peak_height': 100000,      # minimal peak height required for peaks used for RT calibration
    'peak_number_rt_calibration': 15,   # minimal number of selected high-quality peaks required for RT calibration. 
                                        # Samples with fewer selected peaks are dropped out.

    # Number of samples dictate workflow - will set to 10, 1000
    'project_sample_number_small': 1000,
    'project_sample_number_large': 10000,

    # default output names
    'output_feature_table': 'Feature_table.tsv',
    'mass_grid_mapping': "_mass_grid_mapping.csv",
    'annotation_filename': "Annotation_table.tsv",
    'json_empricalCompounds': "_empCpd_json.json",


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 'max_rtime': 300,                   # retention time range (chromatography) in seconds, to auto populate
    # 'cache_mass_traces': False,         # to save memory if not using DB; turn on if need to plot and diagnose
    # 'init_samples_number': 3,           # initiation samples,
    # 'initiation_samples': [],           # if user to specify N samples to initiate data processing, 
                                          # otherwise they are chosen automatically

    }
    