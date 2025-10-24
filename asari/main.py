import argparse
import multiprocessing as mp
import os
import json
import time
from functools import partial

import yaml

from asari import __version__
from .workflow import (
    get_mz_list,
    process_project,
    process_xics,
    read_project_dir,
    create_export_folders,
)
from .default_parameters import PARAMETERS
from .dashboard import read_project, dashboard
from .analyze import estimate_min_peak_height, analyze_single_sample
from .annotate_user_table import annotate_user_featuretable
from .utils import build_boolean_dict, bulk_process
from .qc import generate_qc_report

# -----------------------------------------------------------------------------
# GLOBALS
# -----------------------------------------------------------------------------

booleandict = build_boolean_dict()
SUBCOMMANDS = {
    "analyze": lambda p, a: analyze_single_sample(a.input, parameters=p),
    "process": lambda p, a: _process(p),
    "xic": lambda p, a: _xic(p, a),
    "extract": lambda p, a: _extract(p, a),
    "annotate": lambda p, a: annotate_user_featuretable(a.input, parameters=p),
    "join": lambda *_: print("NOT IMPLEMENTED"),
    "viz": lambda p, a: _viz(p, a),
    "list_workflows": lambda *_: _list_workflows(),
    "qc_report": lambda p, a: _qc_report(p, a),
    "convert": lambda p, a: _convert_raw_if_needed(p, a),
}

# -----------------------------------------------------------------------------
# UTILS
# -----------------------------------------------------------------------------

def _debug(msg: str, enabled: bool):
    if enabled:
        print(msg)

def _convert_raw_if_needed(params, args):
    """Convert Thermo .raw files to mzML when --convert_raw is set."""
    needs_conversion = [
        os.path.abspath(os.path.join(args.input, f))
        for f in os.listdir(args.input)
        if f.endswith(".raw")
    ]
    if needs_conversion:
        from .mzml_converter import mzMLconverter

        mzMLconverter(multicores=params["multicores"]).bulk_convert(needs_conversion)


def _process(params):
    files = read_project_dir(params["input"])
    if not files:
        print("No valid mzML files found in the input directory :(")
        return
    process_project(files, params)


def _xic(params, args):
    process_xics(read_project_dir(args.input), params)


def _extract(params, args):
    mzlist = get_mz_list(args.target)
    print(f"Retrieved {len(mzlist)} target m/z values from {args.target}.\n")
    params["target"] = mzlist
    _process(params)


def _viz(_, args):
    project_desc, cmap, epd, Ftable, Ptable = read_project(args.input)
    table = {
        "full": Ftable,
        "preferred": Ptable,
    }.get(args.table_for_viz, Ftable)
    dashboard(project_desc, cmap, epd, table)


def _qc_report(params, args):
    files = read_project_dir(args.input)
    create_export_folders(params)
    jobs = [
        (
            f,
            os.path.join(
                params["qaqc_reports_outdir"],
                os.path.basename(f).replace(".mzML", "_report.html"),
            ),
            params["spikeins"],
        )
        for f in files
    ]
    bulk_process(generate_qc_report, jobs)


def _list_workflows():
    print(
        "Available Workflows:\n"
        "\t1. LC - default option\n"
        "\t2. GC (pass `--workflow GC`)\n"
        "\t3. Lipidomics LC (NOT IMPLEMENTED)"
    )


# -----------------------------------------------------------------------------
# PARAMETER HANDLING
# -----------------------------------------------------------------------------

# Mapping of CLI argument names to post‑processing lambdas (validator/transformer)
SPECIAL_RULES = {
    # boolean-like strings
    "autoheight": lambda v: booleandict[v],
    "keep_intermediates": lambda v: booleandict[v],
    "anno": lambda v: booleandict[v],
    "debug_rtime_align": lambda v: booleandict[v],
    "compress": lambda v: booleandict[v],
    "drop_unaligned_samples": lambda v: booleandict[v],
    "single_file_qc_reports": lambda v: booleandict[v],
    "convert_raw": lambda v: booleandict[v],
    # numeric fields that must be > 0
    **{k: float for k in (
        "ppm",
        "wlen",
        "max_retention_shift",
        "num_lowess_iterations",
        "min_peak_height",
        "min_prominence_threshold",
        "cal_min_peak_height",
        "min_intensity_threshold",
        "vizualization_max_samples",
        "coelute_threshold",
        "min_peaks",
        "min_peaks_common",
        "min_score_threshold",
    )},
}

PATH_VARS = {
    "input": (os.path.isdir, "Input must be a directory."),
    "target": (os.path.isfile, "Target file not found."),
    "outdir": (lambda p: True, ""),  # allow auto-creation
    "reference": (os.path.exists, "Reference file not found."),
    "reuse_intermediates": (os.path.isdir, "Reuse intermediates must be a directory."),
    "spikeins": (os.path.isfile, "Spikein file not found."),
    "retention_index_standards": (os.path.isfile, "RI standards file not found."),
    "GC_Database": (os.path.exists, "GC database not found."),
}


def _apply_cli_overrides(params: dict, args, verbose: bool = True) -> dict:
    """Overlay CLI args (non‑None) onto parameter dict with minimal boilerplate."""
    pr = partial(_debug, enabled=verbose)

    for key, val in vars(args).items():
        if val is None or key == "run":
            continue
        if key in PATH_VARS and val is not None:
            check, msg = PATH_VARS[key]
            assert check(val), msg
            if key == "output":
                val = os.path.abspath(val)
        # cast / validate
        if key in SPECIAL_RULES:
            caster = SPECIAL_RULES[key]
            val = caster(val) if callable(caster) else caster(val)
            if isinstance(caster, type):
                assert val > 0, f"{key} must be > 0"
        if key == "mode":
            assert val in {"pos", "neg"}, "Mode must be pos or neg."
        if key == "database_mode":
            assert val in {"auto", "ondisk", "memory", "database"}, "Invalid database_mode"
        if key == "storage_format":
            assert val in {"pickle", "json"}, "Invalid storage_format"
        if key == "table_for_viz":
            assert val in {"preferred", "full"}, "Invalid table_for_viz"
        if key == "workflow":
            assert val in {"LC", "GC", "LC_start"}, "Invalid workflow"
        if key == "peak_area":
            assert val in {"sum", "auc", "gauss"}, "Invalid peak_area"
        pr(f"Setting {key} -> {val}")
        params[key] = val

    # derived values -----------------------------------------------
    if params.get("autoheight"):
        try:
            params["min_peak_height"] = estimate_min_peak_height(
                read_project_dir(params["input"]), params
            )
            params["min_intensity_threshold"] = params["min_peak_height"] / 10
        except ValueError as e:
            print(f"Autoheight failed ({e}); using defaults.")

    # ensure dependent thresholds are consistent
    mph = params.get("min_peak_height", 1)
    params.setdefault("min_prominence_threshold", 0.33 * mph)
    params.setdefault("cal_min_peak_height", 10 * mph)
    return params


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def _build_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="asari – GC and LC‑MS metabolomics preprocessing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("run", metavar="subcommand", choices=SUBCOMMANDS.keys())

    # simple one‑liners – defaults taken from default_parameters
    simple = {
        "-m --mode": dict(),
        "--ppm": dict(type=int),
        "-i --input": dict(),
        "-o --outdir": dict(),
        "-j --project": dict(dest="project_name"),
        "-c --multicores": dict(type=int),
        "--reference": dict(),
        "--target": dict(),
        "--database_mode": dict(),
        "--wlen": dict(type=int),
        "--max_retention_shift": dict(type=float),
        "--num_lowess_iterations": dict(type=int),
        "--autoheight": dict(),
        "--min_peak_height": dict(type=float),
        "--min_prominence_threshold": dict(type=float),
        "--cal_min_peak_height": dict(type=float),
        "--min_intensity_threshold": dict(type=float),
        "--peak_area": dict(),
        "--keep_intermediates": dict(),
        "--anno": dict(),
        "--debug_rtime_align": dict(),
        "--compress": dict(),
        "--drop_unaligned_samples": dict(),
        "--reuse_intermediates": dict(),
        "--storage_format": dict(),
        "--single_file_qc_reports": dict(),
        "--spikeins": dict(),
        "--convert_raw": dict(),
        "--table_for_viz": dict(),
        "--vizualization_max_samples": dict(type=int),
        "--workflow": dict(),
        "--retention_index_standards": dict(),
        "--GC_Database": dict(),
        "--coelute_threshold": dict(type=float),
        "--similarity_metric": dict(),
        "--min_peaks": dict(type=int),
        "--min_peaks_common": dict(type=int),
        "--min_score_threshold": dict(type=float),
        "--max_retention_index_error": dict(type=float)
    }

    for flags, kw in simple.items():
        flag_list = flags.split()
        parser.add_argument(*flag_list, **kw)

    parser.add_argument("-p", "--parameters", help="YAML/JSON parameter file")
    parser.add_argument("-v", "--version", action="version", version=__version__)

    return parser.parse_args()


# -----------------------------------------------------------------------------
# MAIN CONTROL FLOW
# -----------------------------------------------------------------------------

def main():
    print(f"\n\n~~~~~~~ Hello from Asari ({__version__}) ~~~~~~~~~\n")
    args = _build_parser()

    # clone default parameters so we don't mutate the import
    params = PARAMETERS.copy()

    # stamp
    params.update({"asari_version": __version__, "timestamp": time.strftime("%Y%m%d-%H%M%S")})

    # parameter file overlay (takes priority over defaults)
    if args.parameters:
        with open(args.parameters) as fh:
            content = fh.read()
        try:
            params.update(yaml.safe_load(content))
        except yaml.YAMLError:
            params.update(json.loads(content))

    # CLI overrides (highest priority)
    params = _apply_cli_overrides(params, args)

    # launch requested sub‑command
    SUBCOMMANDS[args.run](params, args)


# -----------------------------------------------------------------------------
# ENTRY POINT
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    main()
