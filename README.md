Asari
=====

Simple Python program for LC-MS metabolomics data preprocessing.
This focuses on high-resolution data, where most features are specific to a formula based mass, 
and baseline is rarely a problem in peak detection.

- Simple peak dection based on local maxima and prominence.
- Fully incorporating peak quality, selectiviy (on m/z, database, elution).
- Peaks of high quality and selectivity are aligned via formula mass.
- Fast assembly and annotation of serum/plasma metabolomes.
- Prioritized leverage of high mass resolution.

Basic concepts follow https://github.com/shuzhao-li/metDataModel, as

    ├── Experiment
       ├── Sample (injection) 
           ├── MassTrace (EIC)
               ├── Peak
           ├── MassTrace 
               ├── Peak
               ├── Peak
        ...
       ├── Sample 
        ...
       ├── Sample 

Algorithms
==========

Align (correspondence) peaks of high selectivity in both measured data and in reference database via formula mass.

Each sample is checked for mass accuracy. Mass calibration is done if systematic shift > 5 pm.

RT calibration is performed based on limited number of selected high quality peaks. A sample is dropped if too few such peaks are found. 

Because peak detection uses lenient default parameters (high quality peaks used to guide major steps),
weak signal reovery is not implemented. The missing peaks are considered under limit of detection.
One can impute with minimal detected values in each feature (row) in downstream analysis. 


Next to-do
==========

Group features into empicical compounds via mass2chem.

To add kernel density method for grouping m/z features (and RT). Right now, formula_mass anchors correspondence. The remaining peaks are groupped by binning median.

Add SQLite DB for storage.
Before then, we may consider to process data in limited sizes. But without caching mass traces, the memory use is small.

It will be significant to replace pyOpenMS (and OpenMS) and update the chromatogram file format.
Right now, each file takes ~10 seconds to process and half of the time is spent on reading the chromatogram file.

More work is needed for estimation of mass stdev. The current method is too dependent on the limit window. 
Also, if users have spike-in controls, they should be allowed to use first.


Repository
==========
https://github.com/shuzhao-li/asari
