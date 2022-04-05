
In the asari/mummichog packages, the data entities are presented in any of the four types: 
class, namedtuple, JSON style dictionary or implicit list. 
The implicit lists are used sparely as they have reduced clarity. 
Namedtuple is immutable, then limited applications. 
In-memory searches are conducted using indexed dictionaries or dataframes.
JSON is used more common in asari for transparency.

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




### to add to doc/

notebook explaining the algorithms.

notebook example of how to use the library for advanced functions.

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



    
Possible to improve performance in next version:

Use C to rewrite chromatogram constructor.
After initial samples, the peak detection of most features can start from building EIC in C, 
to reduce workload in scipy.signal.find_peak.



Notebooks 
=========

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






