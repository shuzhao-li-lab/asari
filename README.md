Asari
=====

Version 0.8

Simple Python program for LC-MS metabolomics data preprocessing.
This focuses on high-resolution data, where most features are specific to a formula based mass, 
and baseline is rarely a problem in peak detection.

- only for high resolution data. Prioritized leverage of high mass resolution.
- Simple peak dection based on local maxima and prominence.
- Tracking peak quality, selectiviy (on m/z, database, elution), and sequential selectivity.
- reproducible, trackable from features to XICs
- Peaks of high quality and selectivity are aligned via formula mass and epdTrees.
- Fast assembly and annotation of serum/plasma metabolomes based on a reference database.
- Use integers for RT scan numbers and intensities for computing efficiency. 
- Avoid mathematical curves whereas possible for computing efficiency. 
- Performance conscious, memory and CPU uses scalable.

Basic concepts follow https://github.com/shuzhao-li/metDataModel, as

    ├── Experiment
       ├── Sample
           ├── MassTrace
               ├── Peak
           ├── MassTrace 
               ├── Peak
               ├── Peak
        ...
       ├── Sample 
        ...
       ├── Sample 

A sample here is an injection in LC-MS experiments. A MassTrace is an extracted chromatogram (EIC or XIC).
Peak is specific to a sample, but a feature is defined per experiment.

Algorithms
==========

- Chromatogram construction is based on m/z values via flexible bins and frequency counts (in lieu histograms). 
- Peak dection based on local maxima and prominence.
- Align (correspondence) peaks of high selectivity in both measured data and in reference database via formula mass.
- Use information on observed features and epdTrees in the first few samples to guide data extraction and assembly in the remaining data.


Each sample is checked for mass accuracy. Mass calibration is done if systematic shift > 5 pm.

Because peak detection uses lenient default parameters (high quality peaks used to guide major steps),
weak signal reovery is not implemented nor desired. The missing peaks are considered under limit of detection.
One can impute with minimal detected values in each feature (row) in downstream analysis. 

Use
===
Currently 
`python3 -m asari.main neg /Users/shuzhao/projectx/chromdir`

The two arguments are ionization_mode and data_directory.

Next to-do
==========

Update FeatureMap algorithm (borrowed name from OpenMS).

The reference DB is not finalized. 

Group features into empicical compounds via mass2chem.

Add SQLite DB for storage.
Before then, we may consider to process data in limited sizes. But without caching mass traces, the memory use is small.

More work is needed for estimation of mass stdev. The current method is too dependent on the limit window. 
Also, if users have spike-in controls, they should be allowed to use first.


Repository
==========
https://github.com/shuzhao-li/asari
