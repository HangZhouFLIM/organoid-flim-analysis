# Release coverage and acceptance

This is a figure-linked analysis release candidate. The code repository is private; the accompanying data package remains a local candidate pending its Zenodo deposit.

## Verified contents

- 9,133 experimental files; 3,856,634,699 bytes (approximately 3.86 GB), plus administrative documentation.
- 43 configured quantitative/descriptive workflows covering selected analyses in Results 1–6 and their supplementary panels.
- All 495 archived primary reference P-value comparisons matched after exact matching of channel, condition, reference comparison and analysis unit.
- 5,025 pixel-to-measurement values and 676 mask-to-measurement values matched the archived measurements. These are verification counts, not independent sample counts.
- 169 accepted or reconstructed final measurement masks passed all four compartment-metric checks; original accepted segmentation masks are also deposited where available.
- Additional sample-score summaries and the Nile Red emission comparison were checked separately; see `additional_summary_checks.csv`.
- Nineteen regression tests passed, including modified-input rejection, duplicate-coordinate rejection, coordinate-origin checks and rejection of ambiguous statistical matches.
- The pinned Python 3.10.9 environment passed dependency consistency checks. All 43 workflows were run from a copied code folder with full hashes and numerical checks. The new Figure 1j-l PNG/SVG display was separately generated and visually inspected; the other plotting implementations retain their earlier validation.

- 19,252 object-to-observation values match exact source rows and IDs. Figure 2d additionally verifies 2,298 cell areas and centroids across 27 acquisitions, including five reconstructed exclusion-mask combinations.

## Scope and selection

the data package’s `metadata/panel_index.csv` maps figures to data inputs and executable configurations. `file_manifest.csv` maps every experimental asset to its panel(s), processing role, original relative reference and SHA256 checksum. `selection_audit.csv` records table filtering. A shared asset is stored once and linked to every dependent panel.

Only acquisitions/observations supporting the selected Results were included. This includes all measurements behind quantitative comparisons, not only the illustrated examples, and supporting endpoints required by the original multiple-testing families. Whole experiment directories, unused ETC treatments, alternate calibrations, unused differentiation groups and segmentation parameter sweeps are excluded.

Public names use figure and sample identifiers. Source acquisition metadata remain unchanged; elapsed-time analysis uses the original timestamps. Replicate IDs remain technical image/FOV/observation IDs unless a source explicitly supplies a biological hierarchy.

## Remaining gaps before a complete public release

1. Original native microscopy containers and photon-arrival files have not been verified in the available workspace. TIFF/image channels and pixel exports are explicitly labelled as exports, not native raw acquisitions.
2. Final image-only panel identities, crops and display settings remain unresolved for the panels listed in `missing_inputs.csv`. Available quantitative acquisition channels are included, but exact representative assignments are not invented. In particular, full coverage is not claimed for Fig. 1i, Fig. 2 image-only panels, Fig. S4, Fig. S6 or the single-cell UMAP.
3. Fig. S11b, rotenone sample 4: a newly reconstructed mitochondrial mask did not reproduce the archived means. Its regenerated mask is withheld; the archived compartment measurements reproduce the reported statistics. Recovering the original accepted mask/preprocessing record is still necessary for that one upstream measurement.
4. Fig. 3b starts from the deposited pooled spatial-bin summary. The original coordinate-to-axis MATLAB processing is incomplete. Single-cell panels start from deposited sample-level scores; the complete expression-matrix-to-score/UMAP provenance is not yet reconstructed. Time courses start from the deposited time-resolved summaries.
5. Authors/ORCIDs, Zenodo and archived-release identifiers, and licences remain unset in `metadata.json`. They must be supplied before public publication. The private GitHub repository URL is recorded; the Zenodo data have not been uploaded.

The numerical workflows do not reconstruct manually assembled manuscript artwork. Reproduced plots retain their numerical data, pairings, statistical comparisons and explicit error-bar definitions; composite placement and representative-image assembly remain separate.

## Newly linked imaging inputs

113 selected pixel exports (56 Figure 1 acquisitions, 18 supplementary zoom acquisitions and 39 medium-comparison acquisitions) and 63 per-image object exports now link directly to the executable workflows. Figure 2d is one combined collection organized by condition and sample, with provenance retained. See the data package `metadata/added_imaging_coverage.csv` and `metadata/new_acquisition_index.csv`.

Cellular pixel tables for the photon/processing comparisons and Figure 2d, and upstream nuclear imaging inputs, remain unavailable. Source cell masks for Figure 1/its processing supplement are upstream segmentation masks; final pixel-validity/resampling steps cannot be reconstructed without those pixel tables.

## Imaging processing-stage separation

All figure folders distinguish `raw_for_GSLab`, `GSLab_image_exports` and `masks`. There are no mixed `image_exports` directories. The deposited raw GSLab inputs contain 3,552 channel TIFFs and 1,217 original XML files across 699 acquisition/frame groups. All groups have complete selected channel triplets and linked original XML. The pixel-to-raw acquisition crosswalk links 1,216 selected pixel exports to 388 acquisition/frame sets, including previously deposited ATP inputs. The earlier Figure 1/2 and reporter datasets retain their existing acquisition and measurement indexes. Source bytes and metadata are unchanged. Time-lapse frames link to their shared series XML; frame counts are not biological replicate counts. See the data package indexes `metadata/imaging_roles_by_figure.csv`, `metadata/GSLab_input_index.csv` and `metadata/pixel_to_raw_acquisition.csv`.

Figure 1 now contains all eight settings and all 56 acquisitions. Its 168 measurements are recomputed from pixel exports. An additional 21 posthoc implementation checks agree with the original analysis function; these are distinct from the 495 archived reference comparisons. The author grouping table keeps reported approximate photon counts separate from calculated means.

## GitHub package refresh

The reader now lists panels without a data download and checks that eight Zenodo selection/acquisition indexes match the code release before computation. `manifests/workflow_inputs.json` and `docs/PANEL_INPUTS.csv` map all 43 workflows to their deposited inputs and raw asset coverage, including four explicitly verified shared-input relationships. Raw channel sets, GSLab-rendered images and masks retain distinct roles. The 19 tests include altered-release rejection, incorrect frame-link rejection and incomplete-triplet rejection. The previous GitHub-only refresh preserved the data package; the current update adds the five selected Figure 1j-l acquisitions and their verified spatial outputs. See [validation_summary.json](validation_summary.json) for the machine-readable result.

## Figure 1j-l addition

The five author-identified acquisitions and three original MATLAB script variants are now linked. All 789 retained object rows reproduce their archived normalized positions; 250 per-sample intensity bins, 789 phase-bin assignments and the three 50-bin mean/SEM SVG curves verify the selected spatial calculation (5,806 checks in total). Four rows with missing metrics remain in the source tables and do not contribute to bin means. These pooled cell profiles are descriptive. The generic single-sample overwritten summary is not substituted for pooled results. Missing per-pixel exports remain explicitly listed in `metadata/missing_inputs.csv`. See the Figure 1j-l guide and input index.
