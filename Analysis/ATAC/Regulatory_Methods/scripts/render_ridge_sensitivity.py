#!/usr/bin/env python3
"""Individual candidate panels for the separate wider-ridge RNA sensitivity."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import network_compare as core


LABELS = {"distance": "Distance", "cicero": "Donor-weighted Cicero",
          "directnet_adapted": "Adapted DIRECT-NET", "promoter_only": "Promoter only"}
METHODS = ("distance", "cicero", "directnet_adapted")
GRIDS = ("original_replay", "wide")


def read(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def style():
    plt.rcParams.update({"font.size": 6, "axes.titlesize": 6, "axes.labelsize": 6,
                         "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "legend.fontsize": 6, "font.weight": "normal",
                         "axes.titleweight": "normal", "axes.labelweight": "normal",
                         "font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "ps.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": .4,
                         "xtick.major.width": .4, "ytick.major.width": .4})


def forest(lineage, original, wide, limits, path):
    # KEY MESSAGE: The wider-grid distal models must beat a promoter-only baseline tuned on that same grid.
    data = {(row["grid"], row["method"]): row for row in original + wide if row["lineage"] == lineage}
    if set(data) != {(grid, method) for grid in GRIDS for method in METHODS}:
        raise ValueError("Forest comparison coverage is incomplete")
    fig, ax = plt.subplots(figsize=(4.8, 2.6))
    for grid, color, marker, offset, label in [("original_replay", "#9E9E9E", "o", .13, "Original grid"),
                                               ("wide", "#C9265E", "s", -.13, "Wider grid")]:
        for index, method in enumerate(METHODS):
            row = data[grid, method]
            mean, low, high = (float(row[field]) for field in ("improvement", "ci_low", "ci_high"))
            ax.errorbar(mean, index + offset, xerr=[[mean-low], [high-mean]],
                        color=color, marker=marker, markersize=3, linewidth=.7,
                        capsize=1.5, label=label if index == 0 else None)
    ax.axvline(0, color="#9E9E9E", linewidth=.6, linestyle="--")
    ax.set_xlim(*limits)
    ax.set_ylim(2.55, -.55)
    ax.set_yticks(range(3), [LABELS[method] for method in METHODS])
    ax.set_xlabel("RNA MSE gain vs promoter only (training-variance units)")
    title = "Hepatocytes, 39 donors" if lineage == "hepatocyte" else "Source Mesenchymal, 18 donors"
    ax.set_title(title, loc="left", pad=7)
    ax.legend(loc="upper left", bbox_to_anchor=(0, -.31), ncol=2, frameon=False,
              handlelength=1.5, borderaxespad=0)
    fig.subplots_adjust(left=.28, right=.80, top=.83, bottom=.32)
    for index, method in enumerate(METHODS):
        row = data["wide", method]
        ax.text(1.02, index, f"Holm p={float(row['holm_pvalue']):.3g}\n{int(row['n_common_genes']):,} genes",
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6)
    fig.text(.28, .035, "Unadjusted 95% donor-bootstrap intervals; wider-grid Holm family: six.", fontsize=6)
    fig.savefig(path)
    plt.close(fig)


def occupancy(lineage, rows, path):
    # KEY MESSAGE: Alpha occupancy is measured per fitted gene/fold model, rather than per held donor row.
    selected = [row for row in rows if row["lineage"] == lineage]
    methods = ("promoter_only", *METHODS)
    alphas = sorted({float(row["ridge_alpha"]) for row in selected})
    lookup = {(row["grid"], row["method"], float(row["ridge_alpha"])): float(row["fraction_models"]) for row in selected}
    display = [(method, grid) for method in methods for grid in GRIDS]
    matrix = np.array([[lookup.get((grid, method, alpha), np.nan) for alpha in alphas]
                       for method, grid in display])
    if np.any(np.abs(np.nansum(matrix, axis=1)-1) > 1e-10):
        raise ValueError("Alpha frequencies do not sum to one per method/grid")
    fig, ax = plt.subplots(figsize=(4.8, 2.8))
    palette = plt.get_cmap("Greys").copy()
    palette.set_bad("#f1edf5")
    image = ax.imshow(matrix, vmin=0, vmax=1, aspect="auto", cmap=palette)
    ax.set_yticks(range(len(display)), [f"{LABELS[method]}: {'original' if grid == 'original_replay' else 'wide'}"
                                      for method, grid in display])
    ax.set_xticks(range(len(alphas)), [f"{alpha:g}" for alpha in alphas], rotation=45, ha="right")
    ax.set_xlabel("Selected ridge alpha")
    for i in range(len(display)):
        for j in range(len(alphas)):
            value = matrix[i, j]
            if np.isfinite(value) and value >= .01:
                ax.text(j, i, f"{100*value:.0f}", ha="center", va="center", fontsize=6,
                        color="white" if value > .55 else "black")
    ax.set_title("Hepatocytes, 39 donors" if lineage == "hepatocyte" else "Source Mesenchymal, 18 donors", loc="left", pad=7)
    colorbar = fig.colorbar(image, ax=ax, fraction=.045, pad=.025, ticks=[0, .5, 1])
    colorbar.set_ticklabels(["0%", "50%", "100%"])
    colorbar.set_label("Fraction of gene/fold fits")
    fig.subplots_adjust(left=.35, right=.89, top=.84, bottom=.25)
    fig.text(.35, .015, "Numbers are percent of fits; tinted cells: alpha outside original grid.", fontsize=6)
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    core.require_compute()
    receipt = json.loads((args.summary / "completion.json").read_text())
    if receipt.get("status") != "completed" or receipt.get("n_primary_comparisons") != 6 or receipt.get("n_folds") != 10:
        raise ValueError("Completed full wider-grid sensitivity required")
    if args.output.exists():
        raise FileExistsError(args.output)
    original, wide = (read(args.summary / filename) for filename in ("original_grid_reference.tsv", "comparisons.tsv"))
    alpha = read(args.summary / "ridge_alpha_distribution.tsv")
    lower = min(float(row["ci_low"]) for row in original + wide)
    upper = max(float(row["ci_high"]) for row in original + wide)
    width = max(upper-lower, .001)
    limits = min(0, lower)-.08*width, max(0, upper)+.08*width
    style()
    args.output.mkdir(parents=True, exist_ok=False)
    for lineage in ("hepatocyte", "stellate"):
        forest(lineage, original, wide, limits, args.output / f"{lineage}_rna_gain.pdf")
        occupancy(lineage, alpha, args.output / f"{lineage}_alpha_occupancy.pdf")
    (args.output / "completion.json").write_text(json.dumps(dict(status="completed", n_panels=4,
        source_summary=str(args.summary.resolve()), forest_shared_x_limits=limits,
        bootstrap_intervals="unadjusted_95_percent_conditional_on_fitted_cv_predictions",
        primary_holm_family_size=6, scientific_adoption=False), indent=2) + "\n")


if __name__ == "__main__":
    main()
