# Organoid FLIM analysis

This public repository contains the analysis code for Figures 1–6 and Supplementary Figures S1–S15 of the intestinal organoid FLIM manuscript. Start with the [manuscript panel guide](docs/MANUSCRIPT_PANELS.md) to find the command and input files for a published panel. The [downloadable crosswalk](docs/manuscript_panel_crosswalk.csv) maps the 43 original workflow keys to the manuscript’s final numbering.

On 20 November 2026, the Zenodo record was a draft with a reserved DOI. Extract `01_Zenodo_data.rar` and `02_manuscript_figure_data.zip` when you have access to the files. The first contains the selected acquisition and measurement exports. The second contains the final artwork, fitting inputs and parameters, S9 transcriptomic plot data, tables and videos. [Reserved Zenodo DOI](https://doi.org/10.5281/zenodo.23103396)

## Install and run

Use Python 3.10. The original imaging workflows were validated with Python 3.10.9. In a virtual environment, run these commands from the repository root.

```sh
python -m pip install -r environments/requirements-lock.txt
python reproduce.py --data-root ../01_Zenodo_data --verify-only
python reproduce.py --data-root ../01_Zenodo_data --panel Fig05b --output results
python reproduce.py --data-root ../01_Zenodo_data --all --output results
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S7 --output results
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S9C --output results
python -m unittest discover -s tests
```

`--data-root` accepts the path to the extracted data folder. The original entry point checks file hashes before calculation. It writes results, numerical comparisons with the archived reference tables, plots and a run report to the output folder. `--no-plots` runs its numerical checks alone. `--skip-hashes` skips bulk file hashes during repeat runs. Existing figure configurations retain the original statistical branches, pairings and multiple-comparison families.

The companion entry point reconstructs the published fit summaries from the saved per-cell or per-ROI fits. S7 retains the distinction between summaries across cells and summaries that weight each image equally. S9A renders the stored embedding, and S9C renders lineage-level mean expression and detection percentages. The complete, manually assembled figure artwork is supplied in the companion ZIP. [Companion inputs and commands](docs/MANUSCRIPT_PANELS.md)

## Imaging provenance

The original analyses start from the deposited pixel, mask, object or summary exports identified in each configuration. The acquisition channels and original XML files are linked to exported measurements through `metadata/pixel_to_raw_acquisition.csv`. `metadata/missing_inputs.csv` identifies unresolved historical inputs. LAS X and GSLab were used for the upstream acquisition and FLIM exports. [Processing scope](docs/REPRODUCIBILITY.md)

Figure 1 uses eight acquisition settings across 56 acquisitions. Figure 2D combines the selected conditions and samples in one collection. Its workflow rebuilds the recorded exclusion-mask combinations and checks cell areas and centroids against the deposited object rows. [Methods and code](docs/METHODS_CODE_MAP.md) · [Original panel guide](docs/PANEL_GUIDE.md) · [Figure 1J–L](docs/FIG01J_L.md) · [Data folder guide](docs/DATA_FOLDER_STRUCTURE.md) · [Validation](docs/VALIDATION.md)
