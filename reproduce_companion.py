"""Reconstruct the five manuscript panels added to the original analysis release."""
import argparse
from pathlib import Path

from workflows.manuscript_fitting import reconstruct
from workflows.manuscript_transcriptomics import render


PANELS = ("S3", "S7", "S9A", "S9C", "S14")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True, help="Extracted 02_manuscript_figure_data folder")
    parser.add_argument("--panel", choices=(*PANELS, "all"), default="all")
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--curve-id", help="Also write the binwise curve for a specified cell or ROI in one fitting panel")
    args = parser.parse_args()
    if args.curve_id and args.panel not in ("S3", "S7", "S14"):
        parser.error("--curve-id requires one fitting panel.")
    for panel in PANELS if args.panel == "all" else (args.panel,):
        output = args.output / panel
        output.mkdir(parents=True, exist_ok=True)
        if panel in ("S9A", "S9C"):
            render(panel, args.data_root, output)
            print(f"{panel}: rendered figure and plotted-data summary.")
        else:
            count = reconstruct(panel, args.data_root, output, args.curve_id)
            print(f"{panel}: matched descriptive summaries and all {count} saved curves.")


if __name__ == "__main__":
    main()
