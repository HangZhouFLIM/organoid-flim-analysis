# Scope of verified reproduction

The scope below describes the original imaging archive and the `reproduce.py` entry point. The companion archive adds the submitted artwork, S9 plotting values and fit checkpoints. [Manuscript panel guide](MANUSCRIPT_PANELS.md)

The workflow starts at the deposited processed exports and accepted masks. Raw microscope-to-FLIM processing requires original native files and dataset-specific vendor calibration/processing settings that are not fully available here.

The statistical unit is preserved per dataset. Differentiation distributions remain pooled cell/nuclear analyses where organoid identifiers were not retained. Other endpoint inputs are image/FOV summaries unless additional biological hierarchy is documented. Phasor trajectories, pooled spatial profiles and time courses remain descriptive unless an explicit test is implemented.

The legacy label "Dunnett T3" in some original exports refers to a studentized-range calculation. The portable workflow names this Games-Howell-type while preserving the formula and original reference P values. Pairwise Holm families include supporting endpoints from the same samples, including normalized lifetime where used.

Measurements and coordinates are validated before statistics. Lgr5 coordinates use the one-based 246 x 246 grid. Muc2 masks use direct exported coordinate indices and per-image extents. Muc2-negative mitochondrial segmentation is performed within the negative compartment. Accepted masks are the authoritative starting point for the reader workflow; regenerating a segmentation with a new model is a new analysis.

Native acquisition metadata are retained for provenance. Public filenames and organization use figure/sample identifiers, not acquisition dates. Time-lapse calculations retain real elapsed times and channel-specific acquisition timing.

Figure rendering is a reproducible analysis display. Final artwork, image-only examples, crop/display assignments, missing original files and unresolved upstream scripts are itemized in the data coverage register. These are not silently filled with an alternative dataset.

The added Figure 1 and supplementary zoom/medium workflows read their linked pixel CSVs and retain explicitly configured source-table rounding. Figure 1 cell/photon, supplementary processing and Figure 2d cell workflows read exact deposited object rows and IDs. Figure 2d also rebuilds the five required exclusion-mask combinations and checks every exported area and centroid. The source cell/nuclear distribution unit and archived selection are preserved.

Figure 2d inputs form one collection across their source folders; public organization uses condition and sample IDs. No date-specific analyses were introduced. The medium comparison keeps the exact labels mOGM, mENR and IM, separately for Lgr5 and Muc2 reporter datasets.

## Pixel-to-raw acquisition lookup

`metadata/pixel_to_raw_acquisition.csv` links the newly resolved pixel exports to their raw channel sets, original XML, acquisition numbers and GSLab export sequence numbers. For example, the acquisition number in `control_1_06_pxData.csv` is `1`; `06` is the export sequence. The acquisition number is interpreted within its matched assay and condition folder. The crosswalk includes only pixel exports already selected for the Results. Regions and probes from the same acquisition share one raw input set.

Every selected acquisition retains its complete available channel triplets. These are the original GSLab input sets, not a claim that each channel is a separate plotted measurement. Time-lapse frames retain their frame identifiers and share the acquisition XML containing the original timestamps. Original names referenced inside XML are resolved through the file manifest; XML contents were not rewritten when public files were renamed.

This lookup complements `metadata/new_acquisition_index.csv` for the earlier Figure 1/2 additions and `metadata/GSLab_input_index.csv` for all deposited GSLab channel groups. `reproduce.py` checks the complete raw input index and the new pixel-to-raw links before running analyses. Computation still starts at the documented pixel, mask, object or summary exports; native FLIM fitting and historical manual selections are not newly implemented.

## Matching code and data releases

`manifests/data_release.json` binds this code candidate to the selected Zenodo file manifest, panel index, raw-channel index and pixel-to-acquisition lookup. These metadata hashes are checked on every analysis run, including runs using `--skip-hashes`. Bulk experimental-file hashes are checked by default. The portable workflow does not depend on acquisition-date folders or original workstation paths. `manifests/workflow_inputs.json` and `docs/PANEL_INPUTS.csv` provide the per-panel input map.

Figure 1j-l additionally replays the selected MATLAB object-centroid projection and 50-bin summaries from uniquely recovered anchor IDs. It verifies archived positions, sample summaries and original SVG curves. Source MATLAB code is preserved; the tested noninteractive implementation is Python. The supplied folder contained no pixel CSVs, so source-mask-to-object FLIM extraction is still outside that workflow.
