Notes for Developers

Algorithms and work flows
=========================
- From each mzML data file, mass tracks of unique m/z values are extracted, then aligned into a MassGrid.
- Retention time is calibrated for each sample to a common reference sample.
- For each m/z value, corresponding mass tracks from all sample files are summarized into one composite mass track.
- Peak detection is performed on the composite mass trackes to generate a feature list for the experiment.
- The features and mapped back to each sample to extract peak area as intensity values.

The approach of MassGrid and composite mass trackes is highly scalable.
Should one wish to take the traditional way, an alternative workflow allows peak detection upfront, 
followed by retention time alignment and correspondence.

Terminology is mostly defined in metDataModel (https://github.com/shuzhao-li/metDataModel).
Major difference to proteomics is that a feature in LC-MS metabolomics is defined by m/z and retention time across one experiment.
We use `empirical compound` to group degenerate features into a tentative compound/metabolite.
Mass track is used in asari, to cover full range of retention time, because alignment of m/z values is fixed in an early step.


Data formats
============
In the asari/mummichog packages, the data entities are presented in any of the four types: 
class, namedtuple, JSON style dictionary or implicit list. 
The implicit lists are used sparely as they have reduced clarity. 
Namedtuple is immutable, then limited applications. 
In-memory searches are conducted using indexed dictionaries or dataframes.
JSON is used more common in asari for transparency.

```
mass tracks as [( mz, rtlist, intensities ), ...].
Peak format: 
{
    'id_number': 0, 'mz', 'apex', 'left_base', 'right_base', 'height', 'parent_masstrace_id', 
    'rtime', 'peak_area', 'goodness_fitting'
}
isotopic_patterns = [(1.003355, 'M(13C)', 0, 0.2), ...]
```

Mass Tracks
===========
They are used for full RT ranges, thus each mass track has a unique m/z. 
Some chromatogram builders in the field separate the mass traces if there are gaps in RT scans, 
but that creates complexity in m/z alignment and searches, and we avoid that in asari.

Peak detection
==============
The main function is peaks.stats_detect_elution_peaks, whereas a mass_track is separated to RIOs 
(region of interests) based on background noise level. 
The peak detection is based on scipy.signal.find_peaks, a local maxima method with prominence control.
Prominence is important, and it is determined per RIO (e.g. 10% of max intensity). 
For tracks of lower intensity or high noise level, smoothing is applied. 

Retention time calibration
==========================
Calibration or alignment is needed as slight chromatographic shift is common.
A set of set of unambiguous peaks are used for calibration, via LOESS regression to a reference sample.
If a sample fails to align (common for blank controls), the method falls back to distance threshold based alignment.

RT is using scan numbers and they can overlap after calibration, e.g. rt_cal:
        (55, 55), (56, 56), (56, 57), (56, 59), (57, 55), (57, 59), (58, 60), (60, 61), (61, 59), (61, 61), (62, 62), 
        (63, 63), (67, 67), (69, 69), (69, 70), (70, 70), (71, 71), (72, 72), (73, 72), (73, 74), (74, 75), (76, 75), (76, 78), (77, 75), (77, 77), ...,
        (190, 190), (190, 191), (190, 192), (191, 189), (191, 191), (191, 192), (192, 192), (192, 193),...

The default is a LOWESS algorithm.
        Also used in the field include dynamic time warp (DTW) and univariate spline.
        We saw no reason to use them, but people are welcome to implement alternative functions.

        Do alignment function using high-selectivity mass tracks.
        Step 1. get high-selectivity mass tracks among landmarks.
        2. for tracks of highest intensities, do quick peak detection to identify RT apexes.
        Only masstracks with single peaks will be used for RT alignment.
        3. use RT from 2, do LOWESS fit to reference RT values. 
        The fitted function will be recorded for each sample, 
        and applied to all RT scan numbers in the samples when used for CMAP construction.
        Important:
        sample.rt_cal_dict, sample.reverse_rt_cal_dict are kept for changed values only and set within sample RT boundaries.
        This is efficient by ignoring unnecessary tracking, and {} is consistent with samples without RT alignment.
        When samples fail in RT alignment,they are logged in warning and treated at the end as if no alignment is required.



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


Sample
======
Python classes are not easy for parallal computing. 
JSON and individual functions are used in multiprocessing.
A SimpleSample class is later used to faciliate tracking.
Registry example format:

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


CenturionTree
=============
a dictionary, indexing mzList by 100*mz bins.
Because most high-resolution mass spectrometers measure well under 0.01 amu, 
one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).

Use scan numbers 
================
Because scan numbers are integers, they are efficient as indices and should be used for most low-level operations.
When real retention time is used, they are float numbers and not good for indices, 
requiring many more comparison operations and decreasing performance.


MassTrack as full RT range np.array
===================================
In version 1.5, the massTrack format was changed from
( mz, rtlist, intensities ) to ( mz, intensity_track ).
intensity_track is np.array(full RT length).
This increases storage for processed samples, but simplifies
i) CMAP construction
ii) RT index switching btw peak detection functions and others.


Notebooks (wishlist) 
====================

notebooks explaining the algorithms, and examples of how to use the library for advanced functions.

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


