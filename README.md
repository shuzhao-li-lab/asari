# Asari

[![Documentation Status](https://readthedocs.org/projects/asari/badge/?version=latest)](https://asari.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://img.shields.io/badge/DOI-doi.org%2F10.1038%2Fs41467--023--39889--1-blue)](https://doi.org/10.1038/s41467-023-39889-1)

Trackable and scalable Python program for high-resolution metabolomics and exposomics data processing.
Asari supports LC-MS, LC-MS/MS, MS/MS, GC-MS and direct infusion (DI-MS) workflows. 
- Fast enough to process 1000s samples on a laptop computer. 
- Precise enough to find a compound present in a single sample. 
- Annotation is separate from feature extraction, allowing flexible and modular workflows. 
- Many modules and functions for reuse by developers and analysts.
- Cloud friendly and AI ready. 

This is the quick start. For details on each workflow, parameters, output, algorithms and data formats, see [GUIDE.md](GUIDE.md).

The LC-MS method is described in [Li et al. Nature Communications 14.1 (2023): 4113](https://www.nature.com/articles/s41467-023-39889-1).
A web server (https://asari.app), a [full pipeline](https://pypi.org/project/pcpfm/) and [tutorials](https://github.com/shuzhao-li-lab/asari_pcpfm_tutorials/) are available.

## Install

```
pip3 install asari-metabolomics
```

- Add `--upgrade` to update. The GitHub repo (https://github.com/shuzhao-li/asari) is often ahead of PyPI; from source, run `python3 -m asari.main ...`.
- Requires Python 3.9+ (see `pyproject.toml`).
- The core package is compact and cloud friendly. Optional extras: `pip install "asari-metabolomics[viz]"` for the dashboard (`asari viz`), `[gc]` for NMF deconvolution, `[plot]` for plotting tools, or `[all]`.

Input data are centroided mzML files. Thermo .RAW files can be converted by [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser); msconvert in [ProteoWizard](https://proteowizard.sourceforge.io/tools.shtml) handles most vendor formats.

## Quick start

```
asari -h                                     # help
asari list_workflows                         # workflows and descriptions
asari peek -i mydir/sample1.mzML             # inspect a single file
asari process -i mydir/project_dir           # LC-MS feature extraction (default workflow)
```

Common options: `-i/--input` (directory of mzML files, or a file; `--files` for a text list of mzML paths), `-o/--output` (default: current directory), `-j/--project` (project name), `--mode pos|neg`, `--workflow`, `-p` (parameter YAML file; template at `test/parameters.yaml`).
Every run writes a time-stamped result directory, e.g. `<project>_<time stamp>/` for `process`, `<project>_annotation_<time stamp>/` for `annotate`, with a `project.json` recording all parameters.

## Workflows and subcommands

Subcommands are `peek, process, deconvolute, annotate, join, viz`. Each takes `--workflow` (`LC` default, `GC`, `LCMSMS`, `MSMS`, `DIMS`) and, where applicable, an algorithm option.

| Subcommand | Workflow | Algorithm option | Input | Output |
|---|---|---|---|---|
| `peek` | any | | one mzML file | statistics, mass accuracy (vs `db/calibration_mz.json`), scan plot PDF if matplotlib is installed |
| `process` | `LC` (default) | `--peak_area sum` (default), `auc`, `gauss` | mzML files, MS1 only | `Feature_table.tsv` |
| `process` | `GC` | optional `--kovats` | mzML files, optional Kovats RI table | `Feature_table.tsv`, with `RI` column if `--kovats` |
| `process` | `LCMSMS` | `--cluster_ms2 T` (default) or `F` | LC-MS/MS mzML files | `Feature_table.tsv`, MS/MS clusters `ms2_spectra.json`/`.msp` (or per file with `F`) |
| `process` | `MSMS` | `--cluster_ms2 T` (default) or `F` | mzML files, MS/MS only | MS/MS spectra `ms2_spectra.json`/`.msp` (or per file with `F`); no feature table |
| `process` | `DIMS` | | DI-MS mzML files, MS1 only | `Feature_table.tsv`, max intensity per mass track; no RT alignment or peak detection |
| `deconvolute` | `LC` | khipu | feature table | empirical compounds `<project>_list_khipus.json` and `.tsv` |
| `deconvolute` | `GC` | `--deconvolution_algorithm hcl` (default), `adap`, `nmf` | feature table with RI (or `--kovats`) | pseudospectra `<project>_denovo_pseudospectra.json` and `.msp` |
| `annotate` | `LC` | `--search_algorithm mz` (default) | feature table, compound database (`--db`, JSON; default built-in); optional `--preannotation` from `deconvolute` | `Features_<project>.tsv`, `Annotated_features_<project>.json` |
| `annotate` | `LC` | `--search_algorithm rtmz` | feature table, library with m/z and RT (`--db`, JSON; default built-in HILIC pos standards), `--lc_rt_tolerance`; optional `--preannotation` | same as above |
| `annotate` | `LCMSMS`, `MSMS` | entropy search, then `--score_algorithm cosine` (default) or `entropy` | MS/MS spectra (JSON, MSP or process directory), MS/MS library (`--db`, `.npz` recommended; `.msp`, `.mgf`, `.json`, legacy `.pkl`) | `<project>_MS2_annotation.tsv` and `.json` |
| `annotate` | `LCMSMS` | `--search_algorithm mz` or `rtmz` | feature table or `process` directory | shortcut, run as `LC` above (MS1 search) |
| `annotate` | `GC` | `--search_algorithm reverse` (default) | feature table with RI (or `--kovats`), GC-MS library (`--db`, MSP/JSON; default built-in) | `empCpds_`, `Features_<project>.tsv`, `Annotated_pseudospectra_<project>.json`/`.msp` |
| `annotate` | `GC` | `--search_algorithm forward`, `--hires T` (default) or `F` (unit m/z) | pseudospectra (JSON/MSP), or feature table (deconvoluted first by `--deconvolution_algorithm`) | same as above |
| `join` | `LC`, `GC` | `--join_algorithm composite_table_merge` (default) | list of feature tables or process directories | `join_<project>_<time stamp>/asari_join_result.tsv` |
| `join` | `GC` | `--join_algorithm annotation_combine` | list of GC annotation results | combined annotation, same files as GC `annotate` |
| `join` | `LC`, `LCMSMS` | `--join_algorithm annotation_combine` | list of result directories of one `process` run: LC/LCMSMS `annotate`, optional `deconvolute` (preannotation) and the `process` directory | `Features_<project>.tsv`, `Annotated_features_<project>.json` (MS1 features with MS1 and MS/MS annotations); `<project>_MS2_annotation.tsv`/`.json` (all MS/MS annotations) |
| `viz` | `LC`, `GC` | | `process` result directory (run with `--diagnosis T`) | dashboard in web browser |
| `viz` | `GC`, `LCMSMS`, `MSMS` | | `annotate` or `annotation_combine` result directory (LC `annotation_combine` with MS/MS sources) | dashboard with interactive mirror plots |

Notes:
- Scores for spectral matching: `--score_cutoff_cosine` (default 0.4); `--mirror_plot T` exports PDF mirror plots in GC annotation.
- MS/MS libraries are best indexed once to `.npz`: `python -m asari.ms2batch index library.msp library.npz`. See [GUIDE.md](GUIDE.md#indexing-an-msms-library-to-npz).
- `annotate` can be run on different databases and libraries; results of one `process` run (MS1, MS/MS and `deconvolute` preannotation) are combined by `join --join_algorithm annotation_combine`. See [GUIDE.md](GUIDE.md#joining-results).

## Examples

LC-MS, from feature extraction to annotation:

```
asari process -i mydir/project_dir -o mydir -j study1 --mode pos
asari deconvolute -i mydir/study1_<ts>/Feature_table.tsv -o mydir -j study1 --workflow LC
asari annotate -i mydir/study1_<ts>/Feature_table.tsv -o mydir -j study1 --workflow LC \
    --preannotation mydir/study1_deconvolution_<ts>
```

LC-MS/MS:

```
asari process -i mydir/project_dir -o mydir -j dda1 --workflow LCMSMS
asari annotate -i mydir/dda1_<ts> -o mydir -j dda1 --workflow LCMSMS --db MassBank.npz
asari annotate -i mydir/dda1_<ts> -o mydir -j dda1 --workflow LCMSMS --search_algorithm mz   # MS1, same as LC
asari deconvolute -i mydir/dda1_<ts>/Feature_table.tsv -o mydir -j dda1 --workflow LC
# list the result directories above, one per line, in dirs.txt
asari join -i dirs.txt -o mydir -j dda1 --workflow LC --join_algorithm annotation_combine
```

GC-MS:

```
asari process -i mydir/gc_dir -o mydir -j gc1 --workflow GC --kovats KovatsIndex_alkanestandards.tsv
asari annotate -i mydir/gc1_<ts>/Feature_table.tsv -o mydir -j gc1 --workflow GC --db my_GCMS_library.msp
```

Dashboard:

```
asari process -i mydir/project_dir -o mydir -j study1 --diagnosis T
asari viz -i mydir/study1_<ts>
```

## Output

`process` writes one feature table, `Feature_table.tsv`, under the result directory. All features passing signal and shape thresholds are kept, including those present in a single sample; filtering is left to users (e.g. by `snr`, `goodness_fitting`, `cSelectivity`, `detection_counts` or annotation). See [GUIDE.md](GUIDE.md#output) for the full layout.

## Asari suite and related projects

- asari (https://github.com/shuzhao-li/asari, https://pypi.org/project/asari-metabolomics/)
- metDataModel: data models for metabolomics (https://github.com/shuzhao-li-lab/metDataModel)
- mass2chem: utilities for mass spectrometry data and annotation (https://github.com/shuzhao-li-lab/mass2chem)
- khipu: generalized, low-level annotation of MS metabolomics (https://github.com/shuzhao-li-lab/khipu)
- JMS: Json's Metabolite Services (https://github.com/shuzhao-li-lab/JMS)
- pcpfm: full pipeline including asari processing, QC and annotation (https://github.com/shuzhao-li-lab/PythonCentricPipelineForMetabolomics)

Citation: Li et al. Nature Communications 14.1 (2023): 4113, https://doi.org/10.1038/s41467-023-39889-1.
