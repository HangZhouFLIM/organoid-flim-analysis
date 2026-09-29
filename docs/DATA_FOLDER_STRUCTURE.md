# Imaging folder structure

The companion data package: every figure uses separate locations for raw acquisition inputs, GSLab-generated images and masks. A category is present when verified files are deposited; shared files remain stored once and are linked from their companion panels.

```text
FigXX/
  raw_for_GSLab/<panel>/<condition>/
    <figure-sample-prefix>_ch0.tif
    <figure-sample-prefix>_ch1.tif
    <figure-sample-prefix>_ch2.tif
    MetaData/<matching-prefix>.xml
  GSLab_image_exports/<panel>/
    <figure-sample-prefix>_Int.tif
    <figure-sample-prefix>_PP.tif
    <figure-sample-prefix>_PPI2.tif
  masks/
  pixel_exports/
  object_exports/
  source_data/
  statistics/
  acquisition_metadata/
```

The example shows one FLIM channel triplet. Multichannel acquisitions can contain additional triplets; original channel numbers are retained. `raw_for_GSLab` contains the microscope-exported TIFF channels used as GSLab inputs and their original XML files. The XML remains byte-for-byte unchanged, including original acquisition metadata. Some vendor XML uses Windows-1252 without an encoding declaration; its decoding is noted in `GSLab_input_index.csv`.

`GSLab_image_exports` contains processed intensity, phasor, lifetime-colour and label-reference images exported from GSLab. A file under a source folder named `raw_image` is still classified as a GSLab output when its provenance and export format identify it as one.

`masks` contains binary or integer-label masks; processing descriptions distinguish source masks, accepted masks and verified reconstructed masks. Quantitative pixel CSVs and per-object Labels CSVs retain their separate folders.

Raw excitation-scan and STED images use `raw_microscopy`, because those workflows are not GSLab processing. Original NR12S intensity images used for segmentation use `segmentation_inputs`. They are kept out of GSLab-rendered images.

Use `imaging_roles_by_figure.csv` for counts, `imaging_processing_levels.csv` for processing-stage assignments, `GSLab_input_index.csv` for channel/metadata groups, and `file_manifest.csv` for every original source and checksum. Folder names and file prefixes follow figure/condition/sample identities. Figure 2d remains one combined collection.

## Pixel-to-raw acquisition lookup

`metadata/pixel_to_raw_acquisition.csv` links the newly resolved pixel exports to their raw channel sets, original XML, acquisition numbers and GSLab export sequence numbers. For example, the acquisition number in `control_1_06_pxData.csv` is `1`; `06` is the export sequence. The acquisition number is interpreted within its matched assay and condition folder. The crosswalk includes only pixel exports already selected for the Results. Regions and probes from the same acquisition share one raw input set.

Every selected acquisition retains its complete available channel triplets. These are the original GSLab input sets, not a claim that each channel is a separate plotted measurement. Time-lapse frames retain their frame identifiers and share the acquisition XML containing the original timestamps. Original names referenced inside XML are resolved through the file manifest; XML contents were not rewritten when public files were renamed.

This lookup complements `metadata/new_acquisition_index.csv` for the earlier Figure 1/2 additions and `metadata/GSLab_input_index.csv` for all deposited GSLab channel groups. `reproduce.py` checks the complete raw input index and the new pixel-to-raw links before running analyses. Computation still starts at the documented pixel, mask, object or summary exports; native FLIM fitting and historical manual selections are not newly implemented.
