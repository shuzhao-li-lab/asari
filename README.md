Asari
=====

Version 0.8

Python program for high-resolution LC-MS metabolomics data preprocessing, 
with a design focus to be trackable and scalable.

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
           ├── MassTrack
               ├── Peak
           ├── MassTrack 
               ├── Peak
               ├── Peak
        ...
       ├── Sample 
        ...
       ├── Sample 

A sample here is an injection in LC-MS experiments. A MassTrace is an extracted chromatogram (EIC or XIC).
Peak is specific to a sample, but a feature is defined per experiment.

This uses `mass2chem` for mass search and annotation functions.

Algorithms
==========

- Chromatogram construction is based on m/z values via flexible bins and frequency counts (in lieu histograms). 
- Peak dection based on local maxima and prominence.
- Align (correspondence) peaks of high selectivity in both measured data and in reference database via formula mass.
- Use information on observed features and epdTrees in the first few samples to guide data extraction and assembly in the remaining data.


Each sample is checked for mass accuracy. 
Each sample has a recorded function of mass calibration and a function of RT calibration.

Use
===
Currently 
`python3 -m asari.main pos mydir/projectx_dir`

The two arguments are ionization_mode and data_directory.

Next to-do
==========

The reference DB is not finalized. 
Add SQLite DB for storage.


Repository
==========
https://github.com/shuzhao-li/asari
