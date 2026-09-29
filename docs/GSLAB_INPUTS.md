# Raw GSLab inputs and exported measurements

The data package keeps three imaging stages separate:

| Data folder | Contents | Use in the reader workflow |
|---|---|---|
| `raw_for_GSLab` | Original microscope-exported channel TIFFs and matching `MetaData` XML | Acquisition provenance and input for optional GSLab reprocessing |
| `GSLab_image_exports` | Processed intensity, phasor/lifetime-colour and label-reference images | Image inspection; display intensities are not substituted for pixel measurements |
| `masks` | Source labels, accepted measurement masks and documented reconstructed masks | Segmentation provenance and configured compartment measurements |
| `pixel_exports` | Exported per-pixel intensity, phasor and lifetime tables | Quantitative input where the panel configuration enables pixel reconstruction |
| `object_exports` | Exported per-object measurements or linked object input tables | Exact selected object-row and ID reproduction |

## Trace an export to its acquisition

1. Find the panel in `docs/PANEL_INPUTS.csv` or `manifests/workflow_inputs.json`.
2. Use the data package's `metadata/file_manifest.csv` to find its pixel export. Assets shared with supplementary panels can be stored under the companion figure.
3. Find that `pixel_file` in `metadata/pixel_to_raw_acquisition.csv`. The row lists the raw channel files, original XML and public acquisition prefix. Earlier Figure 1/2 additions use `metadata/new_acquisition_index.csv` instead. Every raw channel group is also listed in `metadata/GSLab_input_index.csv`.
4. Use the unchanged XML and the manifest's `original_reference` to interpret channel identities and acquisition settings. The file manifest resolves original source names referenced inside XML to the public figure-based filenames.

For the confirmed naming example, `control_1_06_pxData.csv` comes from acquisition `1` in the matched control/assay folder; `06` is the GSLab export sequence. Acquisition numbers are local to each assay and condition. Region and probe exports can share one raw channel set. A time-lapse name retains its frame identifier, and selected frames share the original series XML and timestamps. Frame counts are not biological replicate counts.

Channel numbers are retained as `_ch0`, `_ch1`, etc. Complete selected triplets are checked; probe identity must be read from the acquisition records rather than inferred from a universal channel number. Some original XML requires Windows-1252 decoding; the index records this and the deposited bytes remain unchanged.

## Reproduction boundary

The added raw files do not change the configured numerical starting points or sample selection. Running `reproduce.py` validates their integrity and linkage, then computes from the selected processed exports. GSLab installation, native photon-arrival fitting, missing historical vendor parameters and manual selections are not automated by this package. Available settings are retained in original XML, `parameters/`, per-panel `configs/` and `manifests/processing_parameters.json`. Missing settings and image assignments remain in the Zenodo gap register.

## Shared raw input sets

The input catalog explicitly links FigS11c to Fig05h and FigS11d to Fig05j: each pair uses the exact same deposited pixel population but different reporter compartments. FigS01d links to Fig01f-h and FigS02_ROI links to FigS02 through their identical deposited source-mask sets. These companion links identify existing raw inputs without duplicating or adding data. Counts refer to indexed channel groups/frame sets, which can include separate probe inputs for the same acquisition; they are not sample counts.
