"""Render S9A and S9C from the published coordinates and expression summaries."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


CELL_TYPES = {
    "absorptive_enterocyte": ("Enterocyte", "#0F4D92", 7881),
    "intestinal_stem_cell": ("Stem cell", "#42949E", 2607),
    "goblet_cell": ("Goblet cell", "#B86B3E", 1313),
    "paneth_cell": ("Paneth cell", "#9A4D8E", 654),
}
GENE_GROUPS = {
    "F1 catalytic": "Atp5f1a Atp5f1b Atp5f1c Atp5f1d Atp5f1e".split(),
    "Fo/stalk": (
        "Atp5mc1 Atp5mc2 Atp5mc3 Atp5me Atp5mf Atp5mg Atp5mj "
        "Atp5mk Atp5pb Atp5pd Atp5pf Atp5po Atp5s"
    ).split(),
    "mt Fo": ["mt-Atp6", "mt-Atp8"],
    "Atp5if1 regulator": ["Atp5if1"],
}


def save(fig, output, panel):
    for ext in ("png", "svg"):
        fig.savefig(output / f"{panel}.{ext}", dpi=450, bbox_inches="tight")
    plt.close(fig)


def render_umap(folder, output):
    cells = pd.read_csv(folder / "umap_coordinates.tsv", sep="\t")
    if (len(cells) != 18806 or cells.cell_id.duplicated().any() or
            cells[["cell_id", "cell_type", "sample"]].isna().any().any() or
            cells["sample"].nunique() != 5 or
            not np.isfinite(cells[["UMAP1", "UMAP2"]]).all().all()):
        raise ValueError("S9A requires the 18,806 unique cells from the five source samples.")
    counts = cells.groupby(["sample", "cell_type"], sort=True).size().rename("cells")
    counts.to_csv(output / "S9A_cell_counts.csv")
    fig, ax = plt.subplots(figsize=(4.9, 4.1))
    background = cells[~cells.cell_type.isin(CELL_TYPES)]
    ax.scatter(background.UMAP1, background.UMAP2, s=0.6, color="#D9D9D9",
               alpha=0.25, linewidths=0)
    for key, (label, color, expected) in CELL_TYPES.items():
        subset = cells[cells.cell_type.eq(key)]
        if len(subset) != expected:
            raise ValueError(f"S9A has {len(subset)} cells annotated as {label}, expected {expected}.")
        ax.scatter(subset.UMAP1, subset.UMAP2, s=1.8, color=color, alpha=0.78,
                   linewidths=0, label=f"{label} ({len(subset):,})")
        ax.text(subset.UMAP1.median(), subset.UMAP2.median(), label, color=color,
                fontsize=6.3, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.75))
    ax.set(xlabel="UMAP1", ylabel="UMAP2", xticks=[], yticks=[])
    ax.set_title("Strict-normal organoid map", loc="left", fontweight="bold", fontsize=8)
    ax.legend(loc="upper right", fontsize=5.8, frameon=False, markerscale=2.5)
    fig.tight_layout(pad=0.9)
    save(fig, output, "S9A")


def render_gene_dots(folder, output):
    dots = pd.read_csv(folder / "gene_dotplot_values.tsv", sep="\t")
    genes = [gene for group in GENE_GROUPS.values() for gene in group]
    expected_index = pd.MultiIndex.from_product([list(CELL_TYPES), genes])
    keys = pd.MultiIndex.from_frame(dots[["annotated_cell_type", "canonical_gene"]])
    observed_labels = dots.annotated_cell_type.map({k: v[0] for k, v in CELL_TYPES.items()})
    expected_counts = dots.annotated_cell_type.map({k: v[2] for k, v in CELL_TYPES.items()})
    values = dots.mean_log_normalized_expression
    if (len(keys) != len(expected_index) or keys.duplicated().any() or
            not expected_index.difference(keys).empty or
            not observed_labels.equals(dots.cell_type_label) or
            not expected_counts.equals(dots.cells) or not np.isfinite(values).all() or
            not dots.pct_cells_detected.between(0, 100).all()):
        raise ValueError("S9C requires detection and expression for all 21 genes in four cell types.")
    groups = values.groupby(dots.canonical_gene)
    low, high = groups.transform("min"), groups.transform("max")
    dots["relative_expression_within_gene"] = (values - low) / (high - low).replace(0, 1)
    dots["point_area"] = 8 + 1.15 * dots.pct_cells_detected
    dots.to_csv(output / "S9C_plotted_values.csv", index=False)
    fig, ax = plt.subplots(figsize=(4.9, 6.3))
    y = {}
    offset = 0
    for label, group in GENE_GROUPS.items():
        for j, gene in enumerate(group):
            y[gene] = offset + j
        lo, hi = y[group[0]], y[group[-1]]
        ax.axhspan(lo - 0.45, hi + 0.45, color="#F7F7F7", zorder=0)
        ax.text(-0.34, (lo + hi) / 2, label, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=5.8, clip_on=False)
        offset = hi + 1.9
    x = dots.annotated_cell_type.map({key: i for i, key in enumerate(CELL_TYPES)})
    cmap = LinearSegmentedColormap.from_list("gene_relative", ["#F4F4F4", "#80B4C5", "#0F4D92"])
    points = ax.scatter(x, dots.canonical_gene.map(y), s=dots.point_area,
                        c=dots.relative_expression_within_gene, cmap=cmap, vmin=0, vmax=1,
                        edgecolors="#272727", linewidths=0.25)
    ax.set_xticks(range(4), [v[0] for v in CELL_TYPES.values()], rotation=25, ha="right")
    ax.set_yticks(list(y.values()), list(y))
    ax.invert_yaxis()
    ax.set_xlim(-0.65, 5.55)
    ax.tick_params(length=0)
    ax.set_title("Complex V genes", loc="left", fontweight="bold", fontsize=8)
    for index, pct in enumerate((25, 50, 75)):
        row = 1.0 + index * 1.1
        ax.scatter([4.28], [row], s=8 + 1.15 * pct, facecolor="white", edgecolor="#272727", linewidth=0.35)
        ax.text(4.5, row, f"{pct}%", va="center", fontsize=6.2)
    ax.text(4.28, 0.35, "detected", fontsize=6.2)
    bar = fig.colorbar(points, ax=ax, fraction=0.04, pad=0.03, ticks=[0, 1])
    bar.ax.set_yticklabels(["lower", "higher"])
    bar.set_label("Relative expression per gene")
    fig.subplots_adjust(left=0.38, right=0.86, bottom=0.08, top=0.94)
    save(fig, output, "S9C")


def render(panel: str, root: Path, output: Path):
    with plt.rc_context({"font.size": 7, "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.7}):
        func = render_umap if panel == "S9A" else render_gene_dots
        func(root / "transcriptomics" / panel, output)
