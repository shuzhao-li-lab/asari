#
# important parameters to estimate: ppm of chrom construction; min_peak_height
#

PARAMETERS = {
    'project_name': 'test_asari',
    'outdir': 'output',
    'database_mode': 'ondisk',          # 'ondisk', 'memory' (run in memory, only small studies), 
                                        # 'mongo' (MongoDB, requiring installation of DB server)
    'mode': 'pos',                      # ionization mode
    'mass_range': (50, 2000),
    #'max_rtime': 300,                   # retention time range (chromatography) in seconds, to auto populate

    'mz_tolerance': 5,                  # ppm, high selectivity meaning no overlap neighbors to confuse; 
                                        # Low selectivity regions will be still inspected to determine the true number of features
    'rtime_tolerance': 10,              # feature rtime shift threshold under 10 seconds; or 10% of rtime
                                        # not used now?
    #

    'multicores': 4,                    # number of threads allowed in parallel processing
    'rt_align': True,                   # False to bypass retention time alignment
    
    'min_intensity_threshold': 100,     # minimal intensity for mass track extraction, filtering baseline
    'min_peak_height': 10000,           # minimal peak height

    'cal_min_peak_height': 100000,      # minimal peak height required for peaks used for RT calibration
    'min_timepoints': 6,                # minimal number of data points in elution profile. scipy find_peaks treat `width` as FWHM, thus half of this value.
    'signal_noise_ratio': 5,            # peak height at least x fold over noise, which is median of non-peak data points.
    #
    'peak_number_rt_calibration': 15,   # minimal number of selected high-quality peaks required for RT calibration. 
                                        # Samples with fewer selected peaks are dropped out.

    # no need to modify below 
    'prominence_window': 30,
    'gaussian_shape': 0.3,              # min cutoff

    # tweaking this section now -
    'cache_mass_traces': False,         # to save memory if not using DB; turn on if need to plot and diagnose
    'output_feature_table': 'Feature_table.tsv',
    'mass_grid_mapping': "_mass_grid_mapping.csv",
    'annotation_filename': "Annotation_table.tsv",
    'json_empricalCompounds': "_empCpd_json.json",

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    'init_samples_number': 3,           # initiation samples, kept in memory
    'initiation_samples': [],           # if user to specify N samples to initiate data processing, 
                                        # otherwise they are chosen automatically

    }
    