# Organoid FLIM analysis

This repository accompanies the figure-linked intestinal organoid FLIM analyses in Results 1–6 and their supplementary figures. Its companion data package is `01_Zenodo_data`; the Zenodo DOI is pending. Imaging data are distributed separately from this code repository. The main entry point uses only deposited files and Python dependencies. Figure configurations select the original statistical branches, pairings and correction families.

## Install

Open a terminal in this folder. Use Python 3.10 (the validated version is recorded under `environments/`). Create a virtual environment, activate it, then install:

```sh
python -m pip install -r environments/requirements-lock.txt
```

## Run

```sh
python reproduce.py --data-root ../01_Zenodo_data --verify-only
python reproduce.py --data-root ../01_Zenodo_data --panel Fig05b --output results
python reproduce.py --data-root ../01_Zenodo_data --all --output results
python -m unittest discover -s tests
```

Input paths are resolved from `--data-root`, not the current workstation. Outputs are written to the requested output directory. The default verifies deposited file checksums before calculation; `--skip-hashes` is for repeated local debugging only. `--no-plots` performs numerical checks without rendering plots.

Each workflow exports recalculated summaries/statistics, comparison checks against the archived reference results, and PNG/SVG plots. `reproduction_report.csv` records successes and failures. Nonzero exit status indicates a computation error, mismatch or unmapped reference comparison. Plots reproduce the numerical analyses; they do not automatically reconstruct manually assembled manuscript composite layouts.

## Data flow

Figure configuration -> deposited pixel / mask / object inputs -> per-image or per-object observations -> summaries -> statistical comparisons -> figures. Calibration phasor centroids remain intensity-weighted; Native acquisition-to-FLIM conversion and historical manual ROIs are done in LAS X  and GSLab.

`workflows/reporter_masks.py` contains selected original segmentation functions with explicit parameters. Final measurement masks and their coordinate origins are indexed in the data package. `manifests/processing_parameters.json` records the common settings; available per-image Cellpose settings are deposited with their masks;  Original model weights/software versions were not inferred from a present-day installation.

## Reproducibility boundaries

Read the data package's `metadata/missing_inputs.csv` and this folder's `docs/REPRODUCIBILITY.md`. Native raw files, final representative-image assignments and upstream FLIM data export and export were done in LAS X and GSLab.
See `metadata.json`; the Zenodo DOI, authors and licences require confirmation before public release. The code repository URL is recorded; the repository is currently private. Upstream code hashes are recorded in `docs/upstream_code_provenance.csv`. No cloud credentials or local environment directories are required.

## Navigation

- [Methods-to-code map](docs/METHODS_CODE_MAP.md) and [CSV](docs/METHODS_CODE_MAP.csv)
- [Panel guide](docs/PANEL_GUIDE.md)
- [Figure 1j-l spatial analysis](docs/FIG01J_L.md)
- [Validation and remaining gaps](docs/VALIDATION.md)

Figure 2d is one combined condition/sample collection. Its code reads the exact source object rows, rebuilds the recorded exclusion-mask combinations and checks cell areas and centroids. The new acquisition crosswalk is `metadata/new_acquisition_index.csv` in the data package.

The companion data folder separates `raw_for_GSLab`, `GSLab_image_exports` and `masks` across every figure. See [the folder guide](docs/DATA_FOLDER_STRUCTURE.md). Figure 1 uses all eight acquisition settings (56 acquisitions), with the author grouping recorded in `parameters/Fig01b-d_acquisition_settings.csv`.

Raw acquisition channel sets and original XML are now linked to the newly resolved exports in `metadata/pixel_to_raw_acquisition.csv`. Both normal runs and `--verify-only` check raw channel groups, XML and these crosswalk links.
