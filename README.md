Asari
=====

Simple Python program for LC-MS metabolomics data preprocessing.
This focuses on high-resolution data, where most features are specific to a formula based mass.

- Simple peak dection based on local maxima and prominence.
- Fully incorporating peak quality, selectiviy (on m/z, database, elution).
- Peaks of high quality and selectivity are aligned via formula mass.
- Fast assembly and annotation of serum/plasma metabolomes.

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

Group degenerate peaks into empicical compounds (based on primary ions in previous step) in each sample.

Mass calibration if needed. RT calibration.

RANSAC alignment of the remaining peaks.





https://github.com/shuzhao-li/asari
