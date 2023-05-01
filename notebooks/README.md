# README for Developers

Please see README at upper directory for general information of asari.

## Overview of algorithms and work flows

- From each mzML data file, mass tracks of unique m/z values are extracted, then aligned into a MassGrid.
- Retention time is calibrated for each sample to a common reference sample.
- For each m/z value, corresponding mass tracks from all sample files are summarized into one composite mass track.
- Peak detection is performed on the composite mass tracks to generate a feature list for the experiment.
- The features are mapped back to each sample to extract peak area as intensity values.

The approach of MassGrid and composite mass tracks is highly scalable. 
This is different from the traditional approach, where peak detection is performed per sample upfront, 
followed by retention time alignment and correspondence.

Terminology is mostly defined in metDataModel (https://github.com/shuzhao-li/metDataModel).
Major difference to proteomics is that a feature in LC-MS metabolomics is defined by m/z and retention time across one experiment.
We use `empirical compound` to group degenerate features into a tentative compound/metabolite.
Mass track is used in asari, to cover full range of retention time, because alignment of m/z values is fixed in an early step.

## Data processing algorithms step by step

1. Build sample registry; see `workflow.process_project`, `workflow.register_samples`.

2. Mass track construction for each sample; see `chromatograms.extract_massTracks_`.

   a. Get all MS1 spectra from a data file, as a list of [(m/z, scan number, intensity), ...]. Index to a dictionary by int(mz * 1000) for efficient retrieval. This is an mzTree.

   b. Create a list of data bins from mzTree. Each bin starts with a value in the mzTree, which has a m/z range around 0.001 and is filtered by minimal required scan number. If two bins are adjacent by 0.001 or within tolerance ppm, they are merged. See `chromatograms.get_thousandth_bins`.

   c. Build mass tracks per data bin. If the m/z range in a data bin is within 2 x tolerance ppm, the bin leads to a single mass track.

   Else, a nearest neighbor (NN) clustering is performed to establish the number of mass tracks. The NN clustering assigns each data point to its nearest "peak mz value". The "peak mz values" are determined by finding peaks in the m/z value distribution, with the requirement that the two peaks need to be separated by the mz tolerance minimally. The m/z value distribution is approximated by a histogram for computing efficiency. See `mass_functions.nn_cluster_by_mz_seeds`.

    A mass track is defined by a consensus m/z value and a list of intensity values. The consensus m/z value is taken as the mean of (median m/z and the m/z at highest intensity). This avoids instable values caused by outliers. When multiple data points exist in the same scan (same RT), max intensity is used. This cleans up duplicate mass peaks from the centroiding process. The mass track is of full RT range, with zeros inserted where intensity is missing. See `chromatograms.extract_single_track_fullrt_length`.

   d. Establish anchor mass tracks by finding m/z differences that match to either 13C/12C isotopes or Na/H adducts. These anchors are considered of higher confidence and prioritized in m/z alignment.

   e. For all samples, use parallel processing; See `workflow.batch_EIC_from_samples_`.

3. Aignment of mass tracks across samples, resulting in the MassGrid; See `CompositeMap.construct_mass_grid`.
    
   a. The sample with the highest number of anchor mass tracks is designated as the reference sample, unless a user specifies a reference sample via input parameters. See `ext_Experiment.get_reference_sample_id`.

   b. If the sample number is no more than a predefined parameter ('project_sample_number_small', default 10), this is considered a small study and a pairwise alignment is performed. When a new sample is aligned to the reference list of m/z values (initially from the reference sample), the anchor mass tracks are compared and mapped first. If there is a systematic difference greater than a preset value (default 1 ppm), all m/z values in the sample are recalibrated against the reference. Because the anchors are of higher confidence and mostly have well resolved m/z values, their alignment has less chance for errors. By completing the alignment of anchors first, the remaining m/z values (significantly fewer than the total features) do not compete against the anchors during their alignment, leading to cleaner results. See `MassGrid.build_grid_sample_wise`, `MassGrid.add_sample`.

   c. Else, for a larger study, the mass alignment is performed by the same NN clustering method that is used in initial mass track construction. See `MassGrid.build_grid_by_centroiding`, `MassGrid.bin_track_mzs`

4. Retention time alignment 

   a. The reference sample was established prior to mass alignment. A set of landmark elution peaks are determined in this reference sample by the criteria: mSelectivity > 0.99, min_peak_height is satisfied (default 1e5), prominence > 20% of peak height and the peak is the only peak on its mass track. These are our "selected_reference_landmark_peaks". See `set_RT_reference`, `peaks.quick_detect_unique_elution_peak`.

   b. For each of the remaining samples, a set of good peaks are selected using the same criteria as above, but limited to from the mass tracks already aligned to the selected_reference_landmark_peaks. They constitute the "good_landmark_peaks" specific to a sample. See `CompositeMap.calibrate_sample_RT`.

   c. Perform a LOWESS (Locally Weighted Scatterplot Smoothing) regression to obtain a function to describe the relationship of the RT values between good_landmark_peaks and selected_reference_landmark_peaks. To prevent spurious results from the regression, we add 10% extension out of both ends as boundaries that the regression must converge to. In practice, the function extrapolates to a dictionary to map all scan numbers between the current sample and the reference sample. The dictionary skips numbers that are identical between the two samples to save memory and computing. The resulting "rt_cal_dict" is a correspondence between sample_rt_numbers and reference_rt_numbers. See `chromatograms.rt_lowess_calibration`

   d. Loop through all remaining samples to obtain the RT calibration function (rt_cal_dict) per sample.

5. Building the composite map

   a. Now that we have established both the m/z alignment in MassGrid and RT alignment between samples, they are linked in the `CompositeMap` class. 
   
   b. Each LC-MS feature should be on a specific location on the composite map (with a consensus m/z and RT). A feature can be conceptualized as a pattern that is repeated in multiple samples. Thus, the pattern will persist and be even enhanced when the signals from multiple samples are superimposed. This concept is implemented into "composite mass tracks", where the intensity values are summed on corresponding mass tracks across all samples after RT calibration. The mass tracks are intensity vectors of the same length, based on scan numbers, after retention time alignment. See `CompositeMap.build_composite_tracks`.
   
6. Detection of elution peaks (features)

   Instead of detecting the same peaks in every sample, asari detects an elution peak only once on the composite mass tracks. Because a composite mass track represents all samples, an elution peak here is equivalent to a feature in the experiment. The shape and prominence of a peak matter in relation to its surroundings. Some instruments generate higher intensity numbers than others arbitrarily. The absolute intensity number is not the most important; a peak is distinctive as long as its shape is good and its height is well above the noise level. We may state that all good peaks look similar. It is the "scaling" factor that differs between metabolites and between platforms. The peak detection in asari is largely dependent on the noise levels. The default parameters are functional over a large range and updated based on statistical analysis of the mass track. See `peaks.audit_mass_track`, `peaks.stats_detect_elution_peaks`.
   
   a. If the max intensity of a mass track is higher that a preset ceiling (1E8), the mass track is rescaled under the preset ceiling for the purpose of peak detection. This serves a normalization purpose to make peak detection parameters robust. Because the composite mass tracks sum up the cumulative intensity of all samples, the intensity could be exceedingly high in large studies. But the detection of elution peaks on composite mass tracks only needs to determine the apex and peak boundaries. The useful quantitative value is the peak area reported on individual samples, not on a composite mass track. After peak detection, the peak height is scaled back using the same scaling factor.
   
   b. If the median intensity on a mass track is below the preset min_intensity_threshold (default 1e3 for Orbitrap data), this is a low-intensity track. Both baseline level and noise level are set to min_intensity_threshold.
   
   c. Else, if over half the data points are above min_intensity_threshold and the median intensity is higher than 10 times of preset min_peak_height (default 1e5 for Orbitrap data), detrend (scipy.signal.detrend) is performed on the mass track. Detrend is a regression method to ensure the baseline is not significantly rising or decreasing over the chromatography. It is a computationally expensive operation and only required for high-signal and high-noise data.

   d. If a track is not low-intensity (see b), the bottom signals are taken as intensity values below the lower quartile plus min_intensity_threshold. The constant of min_intensity_threshold makes this method stable, even when zeros dominate the track. Here, the baseline level and noise level are assigned as the mean and standard deviation of the bottom signals, respectively. 

   e. Smoothing (chromatograms.smooth_moving_average) is applied when the noise level is higher than 1% of max intensity and max intensity is lower than 10 times of the preset min_peak_height. For low noise or high intensity tracks, smoothing is not needed.
   
   f. The mass track is subtracted by a filter (i.e. baseline + noise level). This creates multiple segments of positive intensity values, because some data points are below the filter. Peak detection is performed on the separate regions, because i) it is faster to skip very low signals, ii) this reduces the chance of dealing with multiple peaks simultaneous, and iii) this allows further determination of proper peak prominence based on local signal levels.

   g. Peak detection on a segment of signals. Asari uses a simple local maxima method (scipy.signal.find_peaks), with prominence control that is dynamically determined on each mass track then each segment. Prominence is the vertical distance from a peak top to its lowest contour line. Prominence is calculated in a series of steps: the initial value is 1/3 of min_peak_height; then the greater of prominence and the noise level of a track is used as the new prominence; if the segment is high-intensity and high-noise, the greater of prominence and 5% of max intensity is used as the new prominence. Another parameter, min_peak_height, applies to this method. The parameter min_timepoints is used for controlling both distance between peaks and peak width. The prominence is computed on a sliding window size (default 25 scans). The parameters other than prominence are fairly robust at their default values. See `peaks.stats_detect_elution_peaks`, `peaks.detect_evaluate_peaks_on_roi`.
   
   h. The detected elution peaks are evaluated for peakshape, cSelectivity and SNR. The default filters are set low for these values, so that users can decide on their filtering on the feature table. See `peaks.evaluate_gaussian_peak_on_intensity_list`, `peaks.__peaks_cSelectivity_stats_`, `peaks.compute_noise_by_flanks`. The peaks passing the thresholds (default SNR > 2 and peakshape > 0.5) are reported in a JSON format, with link to the composite mass track identifier. 

7. Feature asignment and annotation

   a. The elution peaks detected from composite mass tracks are stored in CompositeMap.FeatureList (they are experiment-wide features). The retention time is converted from scan numbers to seconds (see `CompositeMap.global_peak_detection`). The features have peak positions and boundaries on the composite mass tracks, which are translated to positions on the individual samples, via the RT alignment dictionaries (see `CompositeMap.extract_features_per_sample`). The sample specific peak areas (default method as the sum of all intensity values within peak boundaries) are recorded into `CompositeMap.FeatureTable`.

   b. The preannotaion is done via another package khipu (https://github.com/shuzhao-li-lab/khipu), where isotopes and adducts are grouped into empirical compounds. 

   c. The empirical compounds are searched against known compound database (default HMDB 4) via another package JMS (https://github.com/shuzhao-li/JMS). The matched isomers (not distinguished without additional information) and formulas are included in asari output. See `experiment.ext_Experiment.annotate`.

8. Data export

   The output directory by asari bears a time stamp, not to overwrite existing data. The recommended feature table is `preferred_Feature_table.tsv`. All peaks are kept in `export/full_Feature_table.tsv`. Annotation is exported into both JSON and tsv formats. See `experiment.ext_Experiment.export_peak_annotation`. MassGrid is exported as a csv file. The composite map is exported as a pickle file, which is used by the visual dashboard. 

9. The asari dashboard for data inpsection

   The dashboard is built using panel (https://panel.holoviz.org/). The dashboard uses these files under the result folder: 'project.json', 'export/cmap.pickle', 'export/epd.pickle' and 'export/full_Feature_table.tsv'. Multiple summary and QC metrics are plotted in one panel with multiple tabs.

   The Feature browser displays a feature on a single composite mass track as a scatter plot. It is important that asari does not use smoothed data for visual inspection, which may be misleading. The mass track viewer is similar to Feature browser, but displays information on all features that are found on a composite mass track.


## Data formats
Basic data concepts follow https://github.com/shuzhao-li/metDataModel.

In the asari/mummichog packages, the data entities are presented in any of the four types: 
class, namedtuple, JSON style dictionary or implicit list. 
The implicit lists are used sparely as they have reduced clarity. 
Namedtuple is immutable, then limited applications. 
In-memory searches are conducted using indexed dictionaries or dataframes.
JSON is used more common in asari for transparency. E.g.
```
isotopic_patterns = [(1.003355, 'M(13C)', 0, 0.2), ...]
```

### Sample
Sample Registry example format:
```
    {
    input_file: '',
    name: '',                               # usually short form of input_file
    ion_mode: '',
    # status
    'status:mzml_parsing': '',
    'status:eic': '',
    'mass_aligned': '',
    'rtime_aligned': '',
    'peaks_extracted': '',
    #
    list_scan_numbers = [],                 # Use scan numbers whereas possible. 
    list_retention_time: [],
    list_mass_tracks: [
        {'id_number': ii, 
                'mz': float,
                'intensity': [],
                }, ...
    ], 
    anchor_mz_pairs: [],                    # mostly isotopic pairs to establish m/z anchors (landmarks)
    number_anchor_mz_pairs: int,
    }
```
A SimpleSample class is later used to faciliate tracking.

### Mass Tracks
They are used for full RT ranges, thus each mass track has a unique m/z. 
Some chromatogram builders in the field separate the mass traces if there are gaps in RT scans, 
but that creates complexity in m/z alignment and searches, and we avoid that in asari.

**MassTrack using full RT range np.array**
Since version 1.5, the massTrack format was changed from ( mz, rtlist, intensities ) to ( mz, intensity_track ).
intensity_track is np.array(full RT length).
This increases storage for processed samples, but simplifies
i) CMAP construction, and
ii) RT index switching btw peak detection functions and others.

### Peaks
Example peak format:
```
    {
    "apex": 500,
    "peak_area": 418918,
    "height": 19686,
    "left_base": 480,
    "right_base": 509,
    "cSelectivity": 0.7652027027027027,
    "parent_masstrack_id": 3738,
    "mz": 332.3312801863069,
    "snr": 3,
    "goodness_fitting": 0.048468150246283814,
    "id_number": "F5474",
    "rtime": 192.6176975050002,
    "rtime_left_base": 185.00743377799978,
    "rtime_right_base": 196.091130897
    }
```
Peaks are organized as features.
A peak is specific to a sample, but a feature is defined at experiment level.

### Empirical Compound
Organizing features into empirical compounds is based on the khipu package:
    Li and Zheng (2023), Generalized tree structure to annotate untargeted metabolomics and stable isotope tracing data. Analytical Chemistry, https://pubs.acs.org/doi/full/10.1021/acs.analchem.2c05810 .
    doi: https://doi.org/10.1101/2023.01.04.522722 .
    Code repo - https://github.com/shuzhao-li-lab/khipu
    
The khipu pre-annotation on experimental features and the search against known compound databases are wrapped in JMS:
    https://github.com/shuzhao-li/JMS).

Briefly:
1. Extract all feature pairs that match to m/z difference patterns in selected isotopes and adducts.
2. Connect the feature pairs to networks via shared features.
3. Recompute each network to a 2-tier tree structure (i.e. khipu), with adducts as tier 1 and isotopes as tier 2.
4. Extend each tree by searching additional patterns of adducts and neutral loss.

Notes:
1. Common anchor peaks are M+H or M-H etc, but will calculate at later round. M* may not show.
2. Adducts and isotopes are combinatorial, under restriction of chemical formulae.
3. We curated isotopic/adduct patterns in this package.
4. epdTree can be improved in future based on real data statistics and more structured cheminformatics.

Empirical Compound example in asari:
```
    {
        "interim_id": "kp12_298.2507",
        "neutral_formula_mass": 298.2507398942377,
        "neutral_formula": null,
        "Database_referred": [],
        "identity": [],
        "MS1_pseudo_Spectra": [
        {
            "apex": 60,
            "peak_area": 28114867,
            "height": 3945532,
            "left_base": 57,
            "right_base": 71,
            "goodness_fitting": 0.8755371723144527,
            "cSelectivity": 1.0,
            "parent_masstrack_id": 1600,
            "mz": 300.26135993003845,
            "snr": 367,
            "id_number": "F122",
            "rtime": 23.439265712999998,
            "rtime_left_base": 22.30214305699998,
            "rtime_right_base": 27.59372140999998,
            "representative_intensity": 28114867,
            "id": "F122",
            "isotope": "13C/12C",
            "modification": "M+H+",
            "ion_relation": "13C/12C,M+H+",
            "parent_epd_id": "kp12_298.2507"
        },
        {
            "apex": 60,
            "peak_area": 157122979,
            "height": 21243341,
            "left_base": 54,
            "right_base": 71,
            "goodness_fitting": 0.9040879272465252,
            "cSelectivity": 1.0,
            "parent_masstrack_id": 1590,
            "mz": 299.25802779197693,
            "snr": 564,
            "id_number": "F46",
            "rtime": 23.439265712999998,
            "rtime_left_base": 21.159260257000017,
            "rtime_right_base": 27.59372140999998,
            "representative_intensity": 157122979,
            "id": "F46",
            "isotope": "M0",
            "modification": "M+H+",
            "ion_relation": "M0,M+H+",
            "parent_epd_id": "kp12_298.2507"
        }
        ],
        "MS2_Spectra": [],
        "list_matches": [
        [
            "C18H34O3_298.250795",
            "",
            1
        ]
        ]
    }
```


## Algorithmic topics

### Use scan numbers whereas possible
Because scan numbers are integers, they are efficient as indices and should be used for most low-level operations.
When real retention time is used, they are float numbers and not good for indices, 
requiring many more comparison operations and decreasing performance.

### Chromatogram building
This is in the chromatogram module. The default method uses `pymzml` to parse mzML files. 
All data points are processed as [(m/z, scan_number, intensity), ...].
The m/z values are binned to 0.001 amu, and assembled to mass tracks based on the expected mass resolution.
The bins can be merged or split in this process.
If multiple data points exist in the same scan on the same mass track, the highest intensity is used.
Mass tracks are superimposed to composite mass tracks, which are input to later peak detection.

### m/z alignment
Because asari mass tracks respect mass resolution, the alignment between samples is expected to be 1:1.
The m/z alignment functions (in asari.mass_functions) use a sort based approach, with marks of originating tracks (e.g. 1 or 2). When the delta between two adjacent m/z values from different samples is minimal and within tolerance, the two values are considered matched.

### Retention time calibration
Calibration or alignment is needed as slight chromatographic shift is common.
A set of set of unambiguous peaks are used for calibration, via LOWESS regression to a reference sample.
Also used in the field include dynamic time warp (DTW) and univariate spline.
We saw no reason to use them, but people are welcome to implement alternative functions.

1. get high-selectivity mass tracks among landmark mass tracks.
2. do quick peak detection in these tracks to identify RT apexes. Only mass tracks with single peaks are used for RT alignment.
3. do LOWESS fit of peaks from step 2 to reference sample. The RT alignment is recorded in two dictionaries: sample.rt_cal_dict, sample.reverse_rt_cal_dict.

Notes:

The RT alignment dictionaries only keep differing values and set within sample RT boundaries. This is efficient by ignoring unnecessary tracking, and {} is consistent with samples without RT alignment. When samples fail in RT alignment (common for blank controls,they are logged in warning and treated at the end as if no alignment is required. 
Example mapped values in rt_cal could look like:

        (55, 55), (56, 56), (56, 57), (56, 59), (57, 55), (57, 59), (58, 60), (60, 61), (61, 59), (61, 61), (62, 62), 
        (63, 63), (67, 67), (69, 69), (69, 70), (70, 70), (71, 71), (72, 72), (73, 72), (73, 74), (74, 75), (76, 75), (76, 78), (77, 75), (77, 77), ...,
        (190, 190), (190, 191), (190, 192), (191, 189), (191, 191), (191, 192), (192, 192), (192, 193),...

### Composite mass tracks
A feature is expected to correspond to peaks in multiple samples on the same mass track and similar elution time.
The conventional approach detects elution peaks in each sample, then match them across samples via correspondence.
Asari assembles composite mass tracks by summing up the signals from each sample on the same mass track.
Then peak detection is only performed on the composite mass tracks, and peaks are mapped back to individual samples based on the result in the composite mass tracks.
This approach offers
1) significant improvement in computational efficiency, as the expensive peak detection is only run once on the composite data not N times in N samples.
2) Peak profiles are often enhanced after combining signals from multiple samples.
3) The composite mass tracks facilitate data tracking and exploration.

### Elution peak detection
The main function is peaks.stats_detect_elution_peaks, whereas a mass_track is separated to ROIs 
(region of interests) based on background noise level. 
The peak detection is based on scipy.signal.find_peaks, a local maxima method with prominence control.
Prominence is important, and it is determined per RIO (e.g. 10% of max intensity). 
For tracks of lower intensity or high noise level, smoothing is applied. 

### Search and annotation functions
They are mostly in the JMS package. It uses many indexed dictionaries to facilitate fast search of m/z values.

Note the difference between m/z alignment functions (in asari.mass_functions) and feature match functions (asari.tools.match_features).
The latter are proper for searching and matching features. 
The former is efficient pair-wise matching for asari mass tracks when multiple matches are not really concerned, 
because asari does not allow ambiguous mass tracks within a data file.

### Multi-core processing
Python classes are not easy for parallal computing. Therefore, we use JSON and individual functions to facilitate multiprocessing.
Currently (version 1.10), two steps are parallel computed: extraction of mass tracks and peak detection.
The m/z alignment (in MassGrid) and RT alignment are using single CPU core.


## Notebooks, current and wishlist
Planned notebooks explaining the algorithms, and examples of how to use the library for advanced functions.

- Single sample processing, inspection, and determine ppm precision.
    Quick peak detection; plot m/z peaks and LC peaks; Figure 1 in paper.

- Annotate and search

- chain with annotation, and MS2 search

- batched processing; separate CMAP per processing batch and assemble afterwards. 
    Or, how to combine multiple asari projects (based on empCpds and massTracks).

- targeted extraction of ROI; allowing lowering m/z boundaries under ppm tolerance
