Asari
=====


Python program for high-resolution LC-MS metabolomics data preprocessing, 
with a design focus to be trackable and scalable.

- only for high resolution data. Prioritized leverage of high mass resolution.
- Simple peak dection based on local maxima and prominence.
- Tracking peak quality, selectiviy (on m/z, database, elution).
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

This uses `mass2chem` and `JMS` for mass search and annotation functions.

Algorithms
==========

- Chromatogram construction is based on m/z values via flexible bins and frequency counts (in lieu histograms). 
- Peak dection based on local maxima and prominence.
- Align (correspondence) peaks of high selectivity in both measured data and in reference database via formula mass.
- Use information on observed features and epdTrees in the first few samples to guide data extraction and assembly in the remaining data.


Each sample is checked for mass accuracy. 
Each sample has a recorded function of mass calibration and a function of RT calibration.

Selectivity is tracked for
- mSelectivity, how distinct are m/z measurements 
- cSelectivity, how distinct are chromatograhic elution peaks
- dSelectivity, how distinct are database records 


Use 
===
Help information:

`python3 -m asari.main`

To process all mzML files under directory mydir/projectx_dir:

`python3 -m asari.main process --mode pos --input mydir/projectx_dir`

To get statistical description on a single file (useful to understand data and parameters):

`python3 -m asari.main analyze --input mydir/projectx_dir/file_to_analyze.mzML`

To get annotation on a tab delimited feature table:

`python3 -m asari.main annotate --mode pos --ppm 10 --input mydir/projectx_dir/feature_table_file.tsv`


Parameters
==========
Only one parameter in asari requires attention, i.e., precision ppm is set at 5 by default. 
Most modern instruments are fine with 5 ppm, but one may want to change if needed.

For the adventurous:
Default parameters are set in `defaul_parameters.py`. 
They can be (work in progress) overwritten by user supplied parameter file in JSON.
Lastly, parameters specified in command line overwrite all the above.


Next to-do
==========

The reference DB is not finalized. 

Repository
==========
https://github.com/shuzhao-li/asari
