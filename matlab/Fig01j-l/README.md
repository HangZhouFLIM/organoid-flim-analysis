# Original Figure 1j-l MATLAB source

The `original/` folder preserves the author-supplied scripts byte-for-byte:

| File | Relevant metric |
|---|---|
| `crypt_to_villus_gradients_manual_2.m` | MeanInt |
| `crypt_to_villus_gradients_manual_3.m` | MeanTauPhase |
| `crypt_to_villus_gradients_manual_4.m` | MeanTauModulation |

These original variants request CSV files and anchor-cell IDs interactively. Their shared primary function name, `crypt_to_villus_gradients_manual_1`, is retained for provenance. For interactive MATLAB use, work on a copy and align its filename with the primary function name; do not add all three copies under that name simultaneously. The original intensity script also contains auxiliary LOWESS/heatmap/boxplot exports, which are outside the three archived unsmoothed profiles. The code calls functions such as `smooth`, `prctile` and `boxplot`; the original MATLAB/toolbox versions were not recorded.

The validated automated route is `python reproduce.py --data-root ../01_Zenodo_data --panel Fig01j-l --output results`. It uses the uniquely recovered anchor IDs in the Zenodo input index and reproduces the selected MATLAB calculations without dialogs. All archived cell projections, available per-sample summaries and the final SVG mean/SEM curves are checked. MATLAB itself was not run during this validation. Original script hashes are in `docs/upstream_code_provenance.csv` relative to the repository root.
