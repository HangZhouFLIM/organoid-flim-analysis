# Figure 1j-l spatial NAD(P)H profiles

Five selected acquisitions are indexed in `Fig01/acquisition_metadata/Fig01j-l_input_index.csv`. They contain 789 retained object rows, of which four have missing intensity/lifetime values. These four rows are retained for traceability and excluded metric by metric from bin summaries. The five acquisitions are technical image/sample identifiers; no additional biological replication hierarchy is inferred.

The data include 15 original FLIM channel TIFFs with 10 original XML files, five NR12S segmentation intensity images with their original XML, five source cell-label masks, five organoid masks, 15 GSLab image exports, five selected object tables and five archived per-sample projection tables. Raw inputs, image exports and masks occupy separate folders. The original MATLAB plots are under `Fig01/figures/Fig01j-l/`; regenerated source tables and a clearly named reproduction display are also deposited.

## Exact analysis replay

From this GitHub folder run:

```sh
python reproduce.py --data-root ../01_Zenodo_data --panel Fig01j-l --output results
```

The entry point reads the per-cell Labels tables. Crypt-base and villus-tip IDs were recovered as the unique pair of source cell centroids reproducing each archived projection; all 789 cell positions match to within 1e-10. The IDs and recovery evidence are in the input index. They are recovered parameters, not newly chosen biological landmarks.

The MATLAB calculation projects centroids onto the crypt-to-villus vector, clamps positions to [0,1], and assigns 50 equal-width bins (centres 0.01-0.99). The right endpoint belongs to bin 50. Four empty object rows have archived position zero and missing metrics. Per-bin means and SEM pool finite cell values; singleton SEM is zero. Percentiles use the original MATLAB-compatible interpolation. The plotted profiles are unsmoothed means with SEM. They are descriptive, with no cell-level hypothesis tests.

The supplied folder contains outputs from several runs using the same generic filenames. The deposited per-sample intensity summaries and phase-bin assignments match all five samples. The generic overall summary instead matches only the final sample's modulation lifetime; it is not used as the pooled result. The three final archived SVGs independently confirm every pooled mean and SEM band (150 bins, within 0.0003 SVG pixel rounding). Unused combined/unselected tables, LOWESS outputs and unrelated MATLAB analyses are excluded.

## Source code and remaining limits

GitHub `workflows/spatial_cells.py` provides the verified noninteractive replay; `parameters/Fig01j-l_spatial.json` records parameters. Original MATLAB intensity, phase-lifetime and modulation-lifetime variants are preserved byte-for-byte in `matlab/Fig01j-l/original/`.

No `pxData.csv` exports were supplied for these acquisitions. The deposited masks are source masks; the workflow does not infer the missing pixel-validity steps or recompute object lifetime means from display TIFFs. A native microscopy container/photon-arrival file and the historical segmentation software/model versions remain unavailable. Quantitative spatial reproduction starts at the selected object exports.
