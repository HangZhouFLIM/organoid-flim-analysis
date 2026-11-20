"""Regression checks for manuscript numbering, saved-fit summaries and S9 input identity."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflows.manuscript_fitting import compare_summary
from workflows.manuscript_transcriptomics import render_gene_dots


class ManuscriptPanelTests(unittest.TestCase):
    def test_crosswalk_covers_all_original_keys_and_separates_rhodamine_panels(self):
        with (ROOT / "docs/manuscript_panel_crosswalk.csv").open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        mapped = {r["Executable legacy panel key"]: r["Current figure or panel"] for r in rows}
        self.assertEqual(set(mapped), {p.stem for p in (ROOT / "configs").glob("*.json")})
        self.assertEqual(mapped["Fig03m_lumen"], "3J")
        self.assertEqual(mapped["FigS07l-n"], "S9K")
        self.assertEqual(mapped["FigS07j-k"], "S9L–N")
        self.assertIn("No current equivalent", mapped["FigS07c"])
        for key in ["FigS07j-k", "FigS07l-n"]:
            cfg = json.loads((ROOT / "configs" / f"{key}.json").read_text())
            self.assertEqual(cfg["display_channel"], "Rhodamine 123")

    def test_saved_summary_requires_matching_group_keys_and_values(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            results = folder / "output"
            results.mkdir()
            original = pd.DataFrame({"condition": ["control"], "mean": [1.2],
                                     "unit": ["ns"], "n_images": [3]})
            original.to_csv(folder / "image_summary.csv", index=False)
            trial = original.set_index("condition")
            compare_summary(trial, folder, results, "image_summary.csv", ["condition"])
            with self.assertRaisesRegex(ValueError, "mean"):
                compare_summary(trial.assign(mean=1.1), folder, results, "image_summary.csv", ["condition"])
            with self.assertRaisesRegex(ValueError, "Group"):
                compare_summary(trial.rename(index={"control": "other"}), folder, results,
                                "image_summary.csv", ["condition"])

    def test_s9c_rejects_incomplete_gene_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            pd.DataFrame([{"annotated_cell_type": "absorptive_enterocyte", "canonical_gene": "Atp5if1",
                           "cell_type_label": "Enterocyte", "cells": 7881, "pct_cells_detected": 10,
                           "mean_log_normalized_expression": 0.1}]).to_csv(
                               folder / "gene_dotplot_values.tsv", index=False, sep="\t")
            with self.assertRaisesRegex(ValueError, "all 21 genes"):
                render_gene_dots(folder, folder)
