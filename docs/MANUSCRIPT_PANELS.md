# Manuscript panels and their sources

The two archives have complementary entry points. Extract each archive into its named folder and pass that folder to the corresponding command. The original analysis keeps its historical filenames and configuration keys so that published commands and data manifests continue to work.

| Extracted folder | Contents | Entry point |
|---|---|---|
| `01_Zenodo_data/` | Acquisition channels and metadata, selected measurement exports, masks and original reference tables | `python reproduce.py --data-root ../01_Zenodo_data --panel KEY --output results` |
| `02_manuscript_figure_data/` | Native ROI decay and response-function samples, fitted parameters, S9 plotting tables, final figures, numbered workbooks and two time-lapse videos | `python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel PANEL --output results` |

The [manuscript-to-workflow CSV](manuscript_panel_crosswalk.csv) identifies the exact configuration, data tables, statistics and command for each of the 43 original workflows. The companion ZIP contains the same CSV at its root.

| Manuscript | Original figure or companion input |
|---|---|
| Figures 1–6 | Original `Fig01`–`Fig06`. The specific command for 3J is `Fig03m_lumen`. |
| S1, S2 | Original S1, S2 |
| S3 | Companion `fitting/S3/` |
| S4–S6 | Original S3–S5 |
| S7 | Companion `fitting/S7/` |
| S8–S10 | Original S6–S8, supplemented by `transcriptomics/S9A/` and `transcriptomics/S9C/` |
| S11–S13 | Original S9–S11. S11 uses `Fig04d` for the ATP-Red 1 calibration. |
| S14 | Companion `fitting/S14/` |
| S15 | Original S12 |

## Panel-specific mappings

- **3J** measures tariquidar-associated luminal TMRM intensity. Use `Fig03m_lumen`. The tissue measurements in 3K–M use `Fig03i-m`. Current 3M shows villus NAD(P)H.
- **S9K** shows total luminal rhodamine 123 intensity for control (8 images), tariquidar (9) and oligomycin (10). Use `FigS07l-n`. S9J shows representative images.
- **S9L–N** show rhodamine 123 intensity, phase lifetime and modulation lifetime in two luminal regions of the same six control images. Use `FigS07j-k`. The original observations carry the generic channel label `NADPH`. The paired input, acquisition assignment and reference metric names identify rhodamine 123. The plots use that measured probe name and retain the original observations and statistical matches.
- **S9B** uses `FigS07b` for source-sample pathway summaries. The five source groups receive equal weight. The companion also contains the archived S9B summary and per-sample tables. The historical `FigS07c` command calculates nuclear complex-V module scores and enterocyte contrasts. That analysis has no plotted panel in the final manuscript.
- **S10G and S10J** have accompanying tariquidar and oligomycin time-lapse videos in `videos/`. Their group time courses are reproduced separately by the original `FigS08g-i` and `FigS08j-l` commands.

## Fitting summaries and curves

Run any single panel or all five companion workflows. Each command compares reconstructed summaries and predictions with the saved checkpoints and returns a nonzero exit code if they disagree.

```sh
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S3 --output results
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S7 --output results
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S14 --output results
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel all --output results
```

| Panel | Fitted observation and source | Reconstructed results |
|---|---|---|
| S3 / Table 9 | 1,715 segmented cells from 12 images in three photon-count groups. Fits support 1F–H. | Six bi-exponential fit quantities, cell-level means and SEMs by photon-count group. |
| S7 / Table 11 | 2,298 segmented cells from 27 images across the four unsplit culture groups in 2D. | Six fit quantities. The primary image-level summary first averages cells within each image, then weights the images equally. A separate table summarizes pooled cells. The archived intensity assignments support the secondary low/high summaries. |
| S14 / Table 15 | 73 image-specific ROIs. S14A–D map to 5B, 5D, 5F and 5H. Each ROI has one bi-exponential and one tri-exponential fit, giving 365 component records. | All 140 summaries of component lifetimes and amplitude fractions, by cohort, condition and model order. The 100 existing comparison rows are supplied in `existing_statistics.csv`. |

Within each fitting folder, `native_decays.csv` gives native counts, `irf_samples.csv` gives the recorded response-function samples and normalized weights, and `acquisition_settings.csv` gives the bin width and pulse timing. `fit_parameters.csv` contains the coefficients. `predicted_decays.csv.gz` is the numerical checkpoint. Each curve retains 129 exported bins. Bins 0–127 contributed to fitting, and bin 128 is the overflow bin. Replaying the archived coefficients reconstructs the full saved curve without optimizing a new fit or testing a new hypothesis. The implementation uses the original S3/S7 reconvolution kernel. Its forward calculation also reproduces both archived S14 model predictions. [Fitting workflow](../workflows/manuscript_fitting.py) · [Reconvolution kernel](../workflows/reconvolution.py)

The `curve_comparisons.csv` output gives the maximum absolute difference in counts per cell or ROI/model. To export observed, replayed and archived counts for one observation, add `--curve-id` and an ID from that panel’s `fit_parameters.csv`. For example:

```sh
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S14 --curve-id Fig05b_control_S01 --output results
```

## Transcriptomic panels

```sh
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S9A --output results
python reproduce_companion.py --data-root ../02_manuscript_figure_data --panel S9C --output results
```

- **S9A** reads `transcriptomics/S9A/umap_coordinates.tsv` and plots the stored embedding of 18,806 cells from five source groups. Four highlighted annotations account for 12,455 cells. Other cell types are plotted in grey. The workflow exports cell counts by source group and annotation.
- **S9C** reads `transcriptomics/S9C/gene_dotplot_values.tsv`. It plots 21 genes, including *Atp5if1*, across four lineages. Dot area depends on the percentage of cells with detectable expression. Colour scales the mean log-normalized expression within each gene across the four lineages. The exported plotting table includes the original values, relative expression and point sizes.

These commands reconstruct the panels from derived plotting data. The originating count data and study metadata are available under [GSE164638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164638) and [GSE165092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165092). Embedding estimation, cell annotation and expression preprocessing precede these plotting steps.

## Artwork and source coverage

The companion folders `figures/main/` and `figures/supplementary/` contain the final assembled figures. `tables/` contains the workbooks for Supplementary Tables 1–7 and 8–16. `fitting/S14/cohorts.csv` maps each fitted ROI to the original main-figure mask. The native decay exports start after ROI selection and histogram extraction. Dataset-specific assignments and missing original image/crop records are documented in the first archive’s `metadata/missing_inputs.csv`.
