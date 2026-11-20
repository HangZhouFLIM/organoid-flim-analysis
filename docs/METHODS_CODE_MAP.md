# Methods-to-code map

Paths below are relative to the repository root or to the extracted original data folder, `01_Zenodo_data/`. Current manuscript numbers and companion workflows are listed in the [manuscript panel guide](MANUSCRIPT_PANELS.md). The public code URL and reserved data DOI are recorded in `metadata.json`.

The [panel input table](PANEL_INPUTS.csv) and [input catalog](../manifests/workflow_inputs.json) link the updated Zenodo assets to all 43 workflows. `python reproduce.py --list-panels` lists valid commands without requiring the data download.

## How to find and run an analysis

1. Find the assay below and open its linked JSON configuration.
2. Run `reproduce.py` with that configuration's panel identifier; the filename stem is the `--panel` value.
3. Inspect the generated summaries, statistics and verification tables in `results/<panel>/`.

For example, from the analysis-folder root:

```sh
python reproduce.py --data-root ../01_Zenodo_data --panel Fig04d --output results
```

This selects ATP Red 1 calibration. Installation is described in [README.md](../README.md). The [machine-readable map](METHODS_CODE_MAP.csv) supplies a complete command for every configured analysis. Output PNG/SVG files reproduce analysis displays, not the manually assembled manuscript layout.

## Shared implementation

| File | Role |
|---|---|
| [reproduce.py](../reproduce.py) | Common entry point; configuration selection, numerical summaries, endpoint/paired/independent tests, correction families, descriptive correlations, specialized summaries and reference checks. |
| [workflows/raw_inputs.py](../workflows/raw_inputs.py) | Validates raw channel triplets, original XML and pixel-to-acquisition links before numerical analysis. |
| [workflows/pixel_measurements.py](../workflows/pixel_measurements.py) | `from_pixels` recalculates configured image means; `from_measurement_masks` extracts measurements using deposited final masks and explicit coordinate origins. |
| [workflows/object_measurements.py](../workflows/object_measurements.py) | `from_object_exports` retrieves exact retained object rows and verifies their IDs and metric values; `verify_object_masks` also rebuilds the recorded Figure 2d exclusion-mask recipes and validates cell geometry; the workflow retains the original observation selection. |
| [workflows/reporter_masks.py](../workflows/reporter_masks.py) | Selected original coordinate reconstruction, GFP percentile grouping and TMRM dot/filament segmentation functions. These functions are provided for inspection; the main reproduction entry point uses deposited final masks. |
| [workflows/specialized_outputs.py](../workflows/specialized_outputs.py) | ATP calibration phasor display, Nile Red equal-region emission comparison, and verification of additional single-cell sample-score summaries. |
| [workflows/plotting.py](../workflows/plotting.py) | Endpoint and paired analysis displays, including significance labels. |
| [manifests/processing_parameters.json](../manifests/processing_parameters.json) | Recorded reporter-mask parameters, coordinate conventions, normalized-lifetime threshold and common statistical rules. |
| [environments/requirements-lock.txt](../environments/requirements-lock.txt) | Pinned dependencies for the portable reproduction environment, not a record of every historical acquisition/segmentation installation. |

## Analysis-to-code map

Every row uses `reproduce.py`. The configuration specifies the data input and analysis branch; implementation modules appear in the CSV map above. Input levels describe what the reader workflow actually recalculates.

| Methods subsection | Analysis | Configuration | Reproduction starts at |
|---|---|---|---|
| NAD(P)H photon-count optimization | Acquisition-setting comparison | [Fig01b-d](../configs/Fig01b-d.json) | Deposited pixel exports linked to selected measurements |
| Spatial analysis along the crypt-villus axis | Cell-based NAD(P)H spatial profiles | [Fig01j-l](../configs/Fig01j-l.json) | Selected object tables and recovered anchor IDs; original MATLAB projections and SVG curves verified |
| NAD(P)H photon-count optimization | Low/intermediate/high photon comparison | [Fig01f-h](../configs/Fig01f-h.json) | Exact retained rows in deposited object exports |
| NAD(P)H photon-count optimization | Zoom comparison | [FigS01b-c](../configs/FigS01b-c.json) | Deposited pixel exports linked to selected measurements |
| NAD(P)H photon-count optimization | Photon-regime valid-ROI percentages | [FigS01d](../configs/FigS01d.json) | Deposited ROI counts |
| NAD(P)H photon-count optimization | Median-filter/downsampling comparison | [FigS02](../configs/FigS02.json) | Exact retained rows in deposited object exports |
| NAD(P)H photon-count optimization | Median-filter/downsampling valid-ROI percentages | [FigS02_ROI](../configs/FigS02_ROI.json) | Deposited ROI counts |
| Autofluorescence excitation scans | Two-photon excitation scans | [FigS03](../configs/FigS03.json) | Deposited excitation-scan values |
| Cell and nuclear segmentation | Differentiation: cellular NAD(P)H distributions | [Fig02d_cells](../configs/Fig02d_cells.json) | Exact retained rows in deposited object exports |
| Cell and nuclear segmentation | Differentiation: nuclear measurements | [Fig02d_nuclei](../configs/Fig02d_nuclei.json) | Retained exported nuclear measurements |
| Live labelling with fluorescent tracers | Culture-medium reporter comparison | [FigS05c-d](../configs/FigS05c-d.json) | Deposited pixel exports linked to selected measurements |
| Spatial analysis along the crypt–villus axis | Grouped TMRM/NAD(P)H spatial profiles | [Fig03b](../configs/Fig03b.json) | Deposited pooled-bin summary |
| Rhodamine 123 and TMRM perturbation assays | Initial mitochondrial-perturbation screen | [Fig03c](../configs/Fig03c.json) | Deposited pixel exports |
| Rhodamine 123 and TMRM perturbation assays | Oligomycin concentration series | [Fig03d-h](../configs/Fig03d-h.json) | Deposited pixel exports |
| Rhodamine 123 and TMRM perturbation assays | Tariquidar epithelial regional analyses | [Fig03i-m](../configs/Fig03i-m.json) | Deposited pixel exports |
| Rhodamine 123 and TMRM perturbation assays | Tariquidar luminal intensity | [Fig03m_lumen](../configs/Fig03m_lumen.json) | Deposited pixel exports |
| Rhodamine 123 and TMRM perturbation assays | Elacridar concentration series | [FigS08a-e](../configs/FigS08a-e.json) | Deposited pixel exports |
| Rhodamine 123 and TMRM perturbation assays | Rhodamine 123 paired luminal regions | [FigS07j-k](../configs/FigS07j-k.json) | Deposited pixel exports and recorded ROI bounds |
| Rhodamine 123 and TMRM perturbation assays | Rhodamine 123 treatment comparisons | [FigS07l-n](../configs/FigS07l-n.json) | Deposited pixel exports |
| Rhodamine 123 and TMRM perturbation assays | Tariquidar time course | [FigS08g-i](../configs/FigS08g-i.json) | Deposited time-resolved image summaries |
| Rhodamine 123 and TMRM perturbation assays | Oligomycin time course | [FigS08j-l](../configs/FigS08j-l.json) | Deposited time-resolved image summaries |
| Mitochondrial morphology analysis | Crypt/villus morphology comparisons | [Fig03n-p](../configs/Fig03n-p.json) | Deposited field-level morphology measurements |
| ATP Red 1 calibration and regional analysis | Control crypt–villus pairing | [Fig04c](../configs/Fig04c.json) | Deposited pixel exports |
| ATP Red 1 calibration and regional analysis | ATP solution calibration and phasor displays | [Fig04d](../configs/Fig04d.json) | Deposited pixel exports |
| ATP Red 1 calibration and regional analysis | Regional drug-treatment comparisons | [Fig04g-h](../configs/Fig04g-h.json) | Deposited pixel exports |
| Oxidation-sensitive probe assays and regional analysis | Control CellROX paired regions | [FigS10b](../configs/FigS10b.json) | Deposited pixel exports |
| Oxidation-sensitive probe assays and regional analysis | Inhibitor-associated CellROX/TMRM/NAD(P)H paired regions | [FigS10c-h](../configs/FigS10c-h.json) | Deposited pixel exports |
| Oxidation-sensitive probe assays and regional analysis | Control/tariquidar ROSTAR paired regions | [FigS10i-l](../configs/FigS10i-l.json) | Deposited pixel exports |
| Lgr5-GFP-defined compartments and mitochondrial NAD(P)H | First Lgr5-GFP-enriched mitochondrial dataset | [Fig05b](../configs/Fig05b.json) | NAD(P)H pixels and final measurement masks |
| Lgr5-GFP-defined compartments and mitochondrial NAD(P)H | Second Lgr5-GFP-enriched mitochondrial dataset | [Fig05d](../configs/Fig05d.json) | NAD(P)H pixels and final measurement masks |
| Lgr5-GFP-defined compartments and mitochondrial NAD(P)H | Wild-type villus mitochondrial regions | [Fig05f](../configs/Fig05f.json) | NAD(P)H pixels and final measurement masks |
| Lgr5-GFP-defined compartments and mitochondrial NAD(P)H | Wild-type crypt mitochondrial regions | [FigS11b](../configs/FigS11b.json) | Final masks/pixels where recovered; archived means for one missing mask |
| Muc2 reporter segmentation and exclusion of long-lifetime fluorescence | Muc2-positive crypt cellular regions | [Fig05h](../configs/Fig05h.json) | NAD(P)H pixels and final measurement masks |
| Muc2 reporter segmentation and exclusion of long-lifetime fluorescence | Muc2-positive villus cellular regions | [Fig05j](../configs/Fig05j.json) | NAD(P)H pixels and final measurement masks |
| Muc2 reporter segmentation and exclusion of long-lifetime fluorescence | Muc2-negative crypt mitochondrial regions | [FigS11c](../configs/FigS11c.json) | NAD(P)H pixels and final measurement masks |
| Muc2 reporter segmentation and exclusion of long-lifetime fluorescence | Muc2-negative villus mitochondrial regions | [FigS11d](../configs/FigS11d.json) | NAD(P)H pixels and final measurement masks |
| Fuel-pathway perturbations and luminal fluorescence | Individual UK-5099/BPTES/etomoxir treatments | [Fig06b](../configs/Fig06b.json) | Deposited pixel exports |
| Fuel-pathway perturbations and luminal fluorescence | Pairwise and triple fuel-inhibitor combinations | [FigS12b](../configs/FigS12b.json) | Deposited pixel exports |
| Fuel-pathway perturbations and luminal fluorescence | UK-5099 luminal comparison | [Fig06d](../configs/Fig06d.json) | Deposited pixel exports |
| Nile Red and BDP-DA1 regional analysis | Nile Red regional and green/red emission comparisons | [Fig06g](../configs/Fig06g.json) | Deposited regional measurements |
| Nile Red and BDP-DA1 regional analysis | BDP-DA1 paired regions | [FigS12d](../configs/FigS12d.json) | Deposited pixel exports |
| Public single-cell RNA-sequencing analysis | Equal-source-sample pathway summaries | [FigS07b](../configs/FigS07b.json) | Deposited source-sample module scores |
| Public single-cell RNA-sequencing analysis | Nuclear ATP-synthase scores and enterocyte contrasts | [FigS07c](../configs/FigS07c.json) | Deposited source-sample module scores |

Figure 1j-l uses [workflows/spatial_cells.py](../workflows/spatial_cells.py), [recorded parameters](../parameters/Fig01j-l_spatial.json) and [archived MATLAB sources](../matlab/Fig01j-l/README.md). This is distinct from the still-summary-based Figure 3 workflow.

## Steps outside the packaged reproduction entry point

- **Native FLIM and manual segmentation:** LAS X/GSLab conversion, historical GUI thresholds, manual ROI correction and original NR12S/Hoechst segmentation settings are not fully implemented or documented by `reproduce.py`. Cell/nuclear comparisons start from retained exported objects; they do not rerun the upstream 3×IQR selection.
- **Reporter segmentation:** `reporter_masks.py` contains GFP and TMRM functions, not a complete Cellpose-SAM driver. Accepted Muc2 masks and final compartment masks are the reproduction inputs. Normalized-lifetime exclusions are already embodied in final measurement masks. There is no separate configured workflow for all historical long-lifetime object analyses.
- **Wild-type crypt exception:** FigS11b rotenone sample 4 retains archived compartment measurements because its original accepted measurement mask remains unrecovered. Do not describe all observations in that configuration as newly extracted from masks.
- **Spatial analysis:** Fig03b starts at a pooled-bin summary; the original MATLAB coordinate-to-axis workflow is not supplied as an executable entry point here.
- **Time courses:** FigS08g-i and FigS08j-l summarize deposited time-resolved measurements. They do not rerun the original XML parsing, pixel-validity filters or baseline-normalization calculation.
- **Morphology:** Fig03n-p compares deposited field-level values. Its entry point does not rerun Nellie, tissue outlining, skeleton extraction or apparent-width measurements.
- **Single-cell RNA sequencing:** FigS07b and FigS07c aggregate deposited source-sample scores. They do not run Scanpy/scDblFinder/Harmony or reconstruct the expression-matrix-to-score/UMAP pipeline.
- **Upstream names and hashes:** [upstream_code_provenance.csv](upstream_code_provenance.csv) identifies original scripts where recorded; a provenance row does not mean that script is included in this repository.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md), [VALIDATION.md](VALIDATION.md), and the companion data package's `metadata/missing_inputs.csv` for the existing coverage record. These statements document the current package, not a change to the experiments or their results.

## Interpretation of configuration labels

Use configurations together with the companion data package's `metadata/file_manifest.csv` to identify biological provenance. In FigS07j-k, the original source is the Rhodamine 123 paired-lumen analysis, although the portable observations retain the channel label `NADPH`. In Fig05j, the original source is the Muc2 villus dataset, although a source-table region label still says Paneth-cell area. This map follows recorded source provenance; it does not silently rename the data fields.

Figure aliases and detailed data locations are listed in [PANEL_GUIDE.md](PANEL_GUIDE.md) and the companion `metadata/panel_index.csv`. Test assignments are indexed in the companion `metadata/statistical_methods_by_panel.csv`; raw and adjusted P values remain in the panel-specific results.
