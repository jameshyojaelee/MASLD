#!/usr/bin/env python3
"""
Cross-dataset LFC overlap scatter plots for MCD RNA-seq.

Uses minimal cutoffs (no filtering); highlights genes with |LFC| >= 0.5 in both datasets.
Also writes a quadrant view (positive LFCs) colored by min(padj).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import seaborn as sns

PRIMARY_MAGENTA = "#C23B75"
ACCENT_PEACH = "#F2A45E"
BASE_GRAY = "#9c9c9c"
LABEL_COLOR = "#1d3557"

# Set context and style
sns.set_context("paper") 
sns.set_style("whitegrid", {'axes.grid': False})

# CRITICAL: Configure fonts for Illustrator (Type 42 = TrueType)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Strict fallback to DejaVu Sans to prevent segfaults
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = ['DejaVu Sans']

# Strict Font Sizes (Journal Standards)
plt.rcParams['font.size'] = 7           # Default text size
plt.rcParams['axes.titlesize'] = 8      # Title size
plt.rcParams['axes.labelsize'] = 8      # X and Y label size
plt.rcParams['xtick.labelsize'] = 6     # X tick size
plt.rcParams['ytick.labelsize'] = 6     # Y tick size
plt.rcParams['legend.fontsize'] = 6     # Legend size

TITLE_SIZE = 8
AXIS_LABEL_SIZE = 8
TICK_SIZE = 6
LEGEND_SIZE = 6
LABEL_TEXT_SIZE = 6

LFC_HIGHLIGHT = 0.5
PADJ_SIG_CUTOFF = 0.1


def load_de_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        raise ValueError(f"{path} missing gene column")
    if "log2FoldChange" not in df.columns:
        raise ValueError(f"{path} missing log2FoldChange column")
    if "padj" not in df.columns:
        raise ValueError(f"{path} missing padj column")

    cols = ["gene", "log2FoldChange", "padj"]
    if "gene_name" in df.columns:
        cols.append("gene_name")
    out = df[cols].copy()
    out = out.rename(columns={"log2FoldChange": "lfc"})
    out["padj"] = out["padj"].fillna(1)
    out["lfc"] = out["lfc"].fillna(0)
    out = out.dropna(subset=["gene", "lfc"])

    # Deduplicate by gene to ensure clean merges.
    agg = {"lfc": "mean", "padj": "min"}
    if "gene_name" in out.columns:
        agg["gene_name"] = "first"
    out = out.groupby("gene", as_index=False).agg(agg)
    return out


def merge_pair(df_x: pd.DataFrame, df_y: pd.DataFrame) -> pd.DataFrame:
    merged = pd.merge(df_x, df_y, on="gene", suffixes=("_x", "_y"))
    return merged


def choose_labels(df: pd.DataFrame) -> pd.Series:
    labels = df["gene"]
    if "gene_name_x" in df.columns or "gene_name_y" in df.columns:
        gx = df["gene_name_x"] if "gene_name_x" in df.columns else pd.Series("", index=df.index)
        gy = df["gene_name_y"] if "gene_name_y" in df.columns else pd.Series("", index=df.index)
        labels = gx.where(gx.notna() & (gx != ""), gy.where(gy.notna() & (gy != ""), labels))
    return labels


def plot_pair(
    merged: pd.DataFrame,
    title: str,
    out_file: Path,
    highlight_threshold: float = LFC_HIGHLIGHT,
    label_genes: Iterable[str] | None = None,
    quadrant_out_file: Path | None = None,
):
    if merged.empty:
        print(f"Skipping {title}: no overlapping genes")
        return

    merged = merged.copy()
    merged["label"] = choose_labels(merged)
    merged["both_up"] = (merged["lfc_x"] >= highlight_threshold) & (merged["lfc_y"] >= highlight_threshold)
    merged["both_down"] = (merged["lfc_x"] <= -highlight_threshold) & (merged["lfc_y"] <= -highlight_threshold)

    plt.figure(figsize=(6.8, 6.0))
    plt.scatter(
        merged["lfc_x"],
        merged["lfc_y"],
        c=BASE_GRAY,
        s=6,
        alpha=0.6,
        edgecolors="none",
    )
    up = merged[merged["both_up"]]
    down = merged[merged["both_down"]]
    if not up.empty:
        plt.scatter(up["lfc_x"], up["lfc_y"], c=PRIMARY_MAGENTA, s=8, alpha=0.85, edgecolors="none", label=f"Both LFC≥{highlight_threshold}")
    if not down.empty:
        plt.scatter(down["lfc_x"], down["lfc_y"], c=ACCENT_PEACH, s=8, alpha=0.85, edgecolors="none", label=f"Both LFC≤-{highlight_threshold}")

    # Correlation line
    r = np.corrcoef(merged["lfc_x"], merged["lfc_y"])[0, 1]
    m, b = np.polyfit(merged["lfc_x"], merged["lfc_y"], 1)
    xline = np.linspace(merged["lfc_x"].min(), merged["lfc_x"].max(), 100)
    plt.plot(xline, m * xline + b, color=LABEL_COLOR, lw=1.1, ls="--", alpha=0.7, label=f"r={r:.2f}")

    plt.axvline(0, color="#cccccc", lw=0.8, ls="--")
    plt.axhline(0, color="#cccccc", lw=0.8, ls="--")
    plt.axvline(highlight_threshold, color="#bbbbbb", lw=0.8, ls=":")
    plt.axhline(highlight_threshold, color="#bbbbbb", lw=0.8, ls=":")
    plt.axvline(-highlight_threshold, color="#bbbbbb", lw=0.8, ls=":")
    plt.axhline(-highlight_threshold, color="#bbbbbb", lw=0.8, ls=":")

    if label_genes:
        lab_df = merged[merged["label"].isin(set(label_genes))].copy()
        lab_df = lab_df.drop_duplicates(subset="label")
        offsets = np.linspace(-3.2, 3.2, num=max(1, len(lab_df)))
        for (offset, (_, row)) in zip(offsets, lab_df.iterrows()):
            plt.annotate(
                row["label"],
                (row["lfc_x"], row["lfc_y"]),
                xytext=(row["lfc_x"] + 1.1, row["lfc_y"] + offset),
                textcoords="data",
                fontsize=LABEL_TEXT_SIZE,
                color=LABEL_COLOR,
                arrowprops=dict(arrowstyle="-", color=LABEL_COLOR, lw=1.1),
            )

    plt.xlabel(title.split("|")[0].strip(), fontsize=AXIS_LABEL_SIZE)
    plt.ylabel(title.split("|")[1].strip(), fontsize=AXIS_LABEL_SIZE)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    plt.legend(loc="upper left", fontsize=LEGEND_SIZE, frameon=False)
    plt.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_file}")

    if quadrant_out_file is None:
        return

    quad = merged[(merged["lfc_x"] >= 0) & (merged["lfc_y"] >= 0)].copy()
    if quad.empty:
        return
    quad["min_padj"] = quad[["padj_x", "padj_y"]].min(axis=1)
    quad["neglog10_padj"] = -np.log10(quad["min_padj"].clip(lower=1e-300))
    sig = quad[quad["min_padj"] < PADJ_SIG_CUTOFF]

    plt.figure(figsize=(6.5, 5.8))
    plt.scatter(
        quad["lfc_x"],
        quad["lfc_y"],
        c="#b9b9b9",
        s=10,
        alpha=0.45,
        edgecolors="none",
    )
    if not sig.empty:
        sig = sig.sort_values("neglog10_padj")
        sc = plt.scatter(
            sig["lfc_x"],
            sig["lfc_y"],
            c=sig["neglog10_padj"],
            cmap="magma",
            s=18,
            alpha=0.95,
            edgecolors="none",
        )
        cb = plt.colorbar(sc)
        cb.set_label("-log10 min(padj_x, padj_y)", fontsize=AXIS_LABEL_SIZE)
        cb.ax.tick_params(labelsize=TICK_SIZE)

    plt.axvline(highlight_threshold, color="#bbbbbb", lw=0.8, ls=":")
    plt.axhline(highlight_threshold, color="#bbbbbb", lw=0.8, ls=":")
    plt.xlabel(title.split("|")[0].strip(), fontsize=AXIS_LABEL_SIZE)
    plt.ylabel(title.split("|")[1].strip(), fontsize=AXIS_LABEL_SIZE)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    xlim = max(1.0, quad["lfc_x"].max())
    ylim = max(1.0, quad["lfc_y"].max())
    plt.xlim(0, xlim * 1.05)
    plt.ylim(0, ylim * 1.05)
    plt.tight_layout()
    quadrant_out_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(quadrant_out_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Wrote {quadrant_out_file}")


def main():
    datasets: Dict[str, Tuple[str, Path]] = {
        "inhouse_week_pooled": (
            "In-house week-pooled LFC",
            Path("RNA-seq/in-house_MCD_RNAseq/analysis_v1_weekPooled_combined/deseq2_results.tsv"),
        ),
        "gse156918": (
            "Paquette et al. LFC",
            Path("RNA-seq/other_MCD_RNAseq/GSE156918/analysis_mcd_vs_control/deseq2_mcd_vs_control.tsv"),
        ),
        "gse205974": (
            "Yue et al. LFC",
            Path("RNA-seq/other_MCD_RNAseq/GSE205974/analysis_mcd_vs_control/deseq2_mcd_vs_control.tsv"),
        ),
    }

    pairs = [
        ("inhouse_week_pooled", "gse156918"),
        ("inhouse_week_pooled", "gse205974"),
        ("gse156918", "gse205974"),
    ]

    label_genes = {
        "Cidec",
        "Mmp12",
        "Fgf21",
        "Lcn2",
        "Nupr1",
        "Atf3",
        "Clec4e",
        "Vcam1",
        "Spp1",
        "Ccl2",
        "Col1a1",
        "Timp1",
    }

    de_tables: Dict[str, pd.DataFrame] = {}
    for key, (_, path) in datasets.items():
        de_tables[key] = load_de_tsv(path)

    out_dir = Path("RNA-seq/other_MCD_RNAseq/cross_dataset/lfc_overlap")
    summary_rows = []
    for xk, yk in pairs:
        x_label, _ = datasets[xk]
        y_label, _ = datasets[yk]
        title = f"{x_label} | {y_label}"
        merged = merge_pair(de_tables[xk], de_tables[yk])
        if merged.empty:
            continue
        out_file = out_dir / f"lfc_overlap_{xk}_vs_{yk}.png"
        quadrant_out = out_dir / f"lfc_overlap_{xk}_vs_{yk}_quadrant.png"
        plot_pair(
            merged,
            title,
            out_file,
            highlight_threshold=LFC_HIGHLIGHT,
            label_genes=label_genes,
            quadrant_out_file=quadrant_out,
        )
        same_sign = (np.sign(merged["lfc_x"]) == np.sign(merged["lfc_y"])).sum()
        both_up = ((merged["lfc_x"] >= LFC_HIGHLIGHT) & (merged["lfc_y"] >= LFC_HIGHLIGHT)).sum()
        both_down = ((merged["lfc_x"] <= -LFC_HIGHLIGHT) & (merged["lfc_y"] <= -LFC_HIGHLIGHT)).sum()
        r = np.corrcoef(merged["lfc_x"], merged["lfc_y"])[0, 1]
        summary_rows.append(
            {
                "pair": f"{xk}_vs_{yk}",
                "n_overlap": len(merged),
                "n_same_sign": int(same_sign),
                "n_both_up_lfc_ge_0.5": int(both_up),
                "n_both_down_lfc_le_-0.5": int(both_down),
                "pearson_r": float(r),
            }
        )

    if summary_rows:
        out_dir.mkdir(parents=True, exist_ok=True)
        summary_path = out_dir / "lfc_overlap_summary.csv"
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
        print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
