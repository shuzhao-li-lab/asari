# Technical details

This is more for Developers. Please see README at upper directory for general information of asari.

## Algorithms and work flows

- From each mzML data file, mass tracks of unique m/z values are extracted, then aligned into a MassGrid.
- Retention time is calibrated for each sample to a common reference sample.
- For each m/z value, corresponding mass tracks from all sample files are summarized into one composite mass track.
- Peak detection is performed on the composite mass tracks to generate a feature list for the experiment.
- The features and mapped back to each sample to extract peak area as intensity values.

The approach of MassGrid and composite mass tracks is highly scalable.
Should one wish to take the traditional way, an alternative workflow allows peak detection upfront, 
followed by retention time alignment and correspondence.

Terminology is mostly defined in metDataModel (https://github.com/shuzhao-li/metDataModel).
Major difference to proteomics is that a feature in LC-MS metabolomics is defined by m/z and retention time across one experiment.
We use `empirical compound` to group degenerate features into a tentative compound/metabolite.
Mass track is used in asari, to cover full range of retention time, because alignment of m/z values is fixed in an early step.

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
In version 1.5, the massTrack format was changed from ( mz, rtlist, intensities ) to ( mz, intensity_track ).
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
Organizing features into empirical compounds is done in the JMS package (https://github.com/shuzhao-li/JMS).

Briefly:
1. find isotopic signatures from list_peaks
2. initiate epdTree classes using the above isotopic signatures; extend by searching common adducts
3. In remaining peaks (isotopes not seen due to low intensity etc), search pairs for primary adducts. 
   Initiate epdTree classes using those.
4. Extend to other adducts for empCpds from the above steps.
5. Reconsolidate overlap epdTrees. E.g. combinations of isotopes and adducts are found separately.

Notes:
1. Common anchor peaks are M+H or M-H etc, but will calculate at later round. M* may not show.
2. Adducts and isotopes are combinatorial, under restriction of chemical formulae.
3. We curated isotopic/adduct patterns in this package.
4. epdTree can be improved in future based on real data statistics and more structured cheminformatics.

Empirical Compound example:
```
    {'interim_id': 15,
        'neutral_formula_mass': 100.112624453,
        'neutral_formula': 'C6H14N',
        'Database_referred': [],
        'identity': [],
        'MS1_pseudo_Spectra': [{'id_number': 'F117',
                'mz': 100.11207049661286,
                'apex': 221.0,
                'ion_relation': 'anchor',
                'parent_epd_id': 15},
            {'id_number': 'F132',
                'mz': 101.11543204162328,
                'apex': 221.0,
                'ion_relation': '13C/12C',
                'parent_epd_id': 15}],
        'MS2_Spectra': [],
        'list_matches': [('C6H14N_100.112624', 'M[1+]', 2),
        ('C6H13N_99.104799', 'M+H[1+]', 2)]},...
```


## Algorithms

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
Mass tracks are input to later peak detection.

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

### Search functions
They are mostly in the JMS package. It uses many indexed dictionaries to facilitate fast search of m/z values.

### Multi-core processing
Python classes are not easy for parallal computing. Therefore, we use JSON and individual functions to facilitate multiprocessing.
Currently (version 1.7), two steps are parallel computed: extraction of mass tracks and peak detection.
The m/z alignment (in MassGrid) and RT alignment are using single CPU core.


## Notebooks (wishlist) 
Planned notebooks explaining the algorithms, and examples of how to use the library for advanced functions.

- Single sample processing, inspection, and determine ppm precision.
    Quick peak detection; plot m/z peaks and LC peaks; Figure 1 in paper.

- Process data without upfront LC alignment

- regression match of datasets of different LC time

- Annotate and search

- OpenMS based workflow

- chain with annotation, and MS2 search

- batched processing; separate CMAP per processing batch and assemble afterwards. 
    Or, how to combine multiple asari projects (based on empCpds and massTracks).

- targeted extraction of ROI; allowing lowering m/z boundaries under ppm tolerance

- alternative workflow as XCMS, peak detection first followed by correspondence

