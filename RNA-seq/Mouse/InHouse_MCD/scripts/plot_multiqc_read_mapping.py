#!/usr/bin/env python3
"""Generate a MultiQC-style plot of unique vs multi-mapped read percentages."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot per-sample uniquely and multi-mapped read percentages from MultiQC STAR output.",
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("qc/multiqc/multiqc_data/multiqc_star.txt"),
        help="Path to multiqc_star.txt (tab-delimited).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("qc/multiqc/plots/multiqc_read_mapping_percent.png"),
        help="Path for the output figure (PNG).",
    )
    parser.add_argument(
        "--include-pass1",
        action="store_true",
        help="Keep STAR pass1 rows (by default they are removed).",
    )
    return parser.parse_args()


def load_mapping_percentages(data_path: Path, include_pass1: bool) -> pd.DataFrame:
    df = pd.read_csv(data_path, sep="\t")
    if not include_pass1:
        df = df[~df["Sample"].str.contains("_STARpass1", na=False)]
    if df.empty:
        raise ValueError("No samples available to plot after applying filters.")

    columns = {
        "uniquely_mapped_percent": "Uniquely mapped",
        "multimapped_percent": "Multi mapped",
    }
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise KeyError(f"Missing expected columns in {data_path}: {', '.join(missing)}")

    plot_df = df[["Sample", *columns.keys()]].rename(columns=columns)
    plot_df["Sample"] = pd.Categorical(
        plot_df["Sample"], categories=list(dict.fromkeys(plot_df["Sample"])), ordered=True
    )
    return plot_df.sort_values("Sample")


def make_plot(plot_df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.size": 16,
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#333333",
            "axes.labelcolor": "#333333",
            "axes.spines.right": False,
            "axes.spines.top": False,
        }
    )

    fig_width = max(10.0, len(plot_df) * 0.5)
    fig, ax = plt.subplots(figsize=(fig_width, 8))

    palette = {
        "Uniquely mapped": "#F2A45E",  # peach
        "Multi mapped": "#C23B75",  # magenta
    }

    x = range(len(plot_df))
    bottom = None
    for column in ["Uniquely mapped", "Multi mapped"]:
        values = plot_df[column].to_numpy()
        bottom_values = bottom if bottom is not None else None
        ax.bar(x, values, label=column, bottom=bottom_values, color=palette[column])
        bottom = values if bottom is None else bottom + values

    ax.set_ylim(0, 100)
    ax.set_ylabel("Reads (%)", fontsize=18, fontweight="bold")
    ax.set_xlabel("Sample", fontsize=18, fontweight="bold")
    ax.set_title(
        "STAR alignment read mapping percentages", fontsize=22, fontweight="bold", pad=16
    )
    legend = ax.legend(
        title="Category", fontsize=14, title_fontsize=16, frameon=True, loc="upper right"
    )
    legend.get_title().set_fontweight("bold")
    legend.get_frame().set_edgecolor("#d0d0d0")
    legend.get_frame().set_linewidth(1)
    ax.set_xticks(list(x))
    ax.set_xticklabels(plot_df["Sample"], rotation=45, ha="right", fontsize=14)
    ax.tick_params(axis="y", labelsize=14, colors="#333333")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=100))
    ax.yaxis.grid(True, color="#e6e6e6", linewidth=1, linestyle="-")
    ax.xaxis.grid(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color("#333333")
        ax.spines[spine].set_linewidth(1.2)

    fig.tight_layout()

    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    plot_df = load_mapping_percentages(args.data, args.include_pass1)
    make_plot(plot_df, args.output)


if __name__ == "__main__":
    main()
