#!/usr/bin/env python3
"""Generate TPM distribution ridgeline plots for bundled RNA-seq datasets."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import gzip


ROOT = Path(__file__).resolve().parents[2]

ORDERED_DATASETS: list[tuple[str, str]] = [
    ("MCD Week 1", "mcd_week1.tsv.gz"),
    ("MCD Week 2", "mcd_week2.tsv.gz"),
    ("MCD Week 3", "mcd_week3.tsv.gz"),
    ("MCD Week pooled (combined)", "mcd_week_pooled_combined.tsv.gz"),
    ("GSE156918 (external MCD)", "other_mcd_gse156918.tsv.gz"),
    ("GSE205974 (external MCD)", "other_mcd_gse205974.tsv.gz"),
    ("GSE130970 NAS high", "gse130970_nas_high.csv.gz"),
    ("GSE130970 NAS low", "gse130970_nas_low.csv.gz"),
    ("GSE130970 Fibrosis", "gse130970_fibrosis.csv.gz"),
    ("GSE135251 NAS high", "gse135251_nas_high.csv.gz"),
    ("GSE135251 NAS low", "gse135251_nas_low.csv.gz"),
    ("GSE135251 Fibrosis", "gse135251_fibrosis.csv.gz"),
]

ALIAS_LABELS = {
    "MCD Week pooled": "MCD Week pooled (combined)",
    "MCD Week pooled (combined)": "MCD Week pooled (combined)",
    "GSE156918": "GSE156918 (external MCD)",
    "GSE156918 (external MCD)": "GSE156918 (external MCD)",
    "GSE205974": "GSE205974 (external MCD)",
    "GSE205974 (external MCD)": "GSE205974 (external MCD)",
    "GSE130970 NAS high": "GSE130970 NAS high",
    "GSE135251 NAS high": "GSE135251 NAS high",
}

RAW_MCD_COUNTS = ROOT / "RNA-seq" / "in-house_MCD_RNAseq" / "counts" / "featurecounts" / "gene_counts.txt"
RAW_MCD_METADATA = ROOT / "RNA-seq" / "in-house_MCD_RNAseq" / "metadata" / "samples.tsv"
BIOTYPE_MAP_PATH = ROOT / "streamlit_deg_explorer" / "data" / "ensembl_gene_biotypes.tsv.gz"

RAW_EXTERNAL_MCD = {
    "GSE156918 (external MCD)": {
        "counts": ROOT
        / "RNA-seq"
        / "other_MCD_RNAseq"
        / "GSE156918"
        / "counts"
        / "featurecounts"
        / "gene_counts.txt",
        "metadata": ROOT / "RNA-seq" / "other_MCD_RNAseq" / "GSE156918" / "metadata" / "samples.tsv",
    },
    "GSE205974 (external MCD)": {
        "counts": ROOT
        / "RNA-seq"
        / "other_MCD_RNAseq"
        / "GSE205974"
        / "counts"
        / "featurecounts"
        / "gene_counts.txt",
        "metadata": ROOT / "RNA-seq" / "other_MCD_RNAseq" / "GSE205974" / "metadata" / "samples.tsv",
    },
}

RAW_PATIENT_DATASETS = {
    "GSE130970 NAS high": "GSE130970",
    "GSE130970 NAS low": "GSE130970",
    "GSE130970 Fibrosis": "GSE130970",
    "GSE135251 NAS high": "GSE135251",
    "GSE135251 NAS low": "GSE135251",
    "GSE135251 Fibrosis": "GSE135251",
}


def load_biotype_map(path: Path) -> dict[str, str]:
    if not path.exists():
        print(f"Warning: Biotype map not found at {path}")
        return {}
    mapping = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("gene_id") or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                gid = parts[0].split(".")[0]
                bio = parts[1]
                mapping[gid] = bio
    return mapping


def _normalize_sample_ids(sample_ids: Iterable[str] | None) -> tuple[str, ...] | None:
    if not sample_ids:
        return None
    cleaned = [str(s).strip() for s in sample_ids if str(s).strip()]
    if not cleaned:
        return None
    return tuple(sorted(set(cleaned)))


def _featurecounts_sample_id(col: str) -> str:
    path = Path(str(col))
    name = path.name
    if ".Aligned" in name:
        return name.split(".Aligned", 1)[0]
    parts = path.parts
    if len(parts) >= 2:
        return parts[-2]
    return name.split(".", 1)[0]


def _subset_featurecounts_counts(
    counts: pd.DataFrame, sample_ids: Iterable[str] | None
) -> pd.DataFrame:
    sample_ids = _normalize_sample_ids(sample_ids)
    if not sample_ids:
        return counts
    wanted = set(sample_ids)
    selected = [col for col in counts.columns if _featurecounts_sample_id(col) in wanted]
    if not selected:
        return counts
    return counts[selected]


def _subset_counts_matrix(counts: pd.DataFrame, sample_ids: Iterable[str] | None) -> pd.DataFrame:
    sample_ids = _normalize_sample_ids(sample_ids)
    if not sample_ids:
        return counts
    wanted = set(sample_ids)
    selected = [col for col in counts.columns if col in wanted]
    if not selected:
        return counts
    return counts[selected]


def load_mcd_sample_ids(metadata_path: Path) -> set[str]:
    if not metadata_path.exists():
        return set()
    df = pd.read_csv(metadata_path, sep="\t")
    if "sample_id" not in df.columns:
        return set()
    if "diet" not in df.columns:
        return set(df["sample_id"].astype(str))
    mask = df["diet"].astype(str).str.contains("mcd", case=False, na=False)
    if not mask.any():
        return set(df["sample_id"].astype(str))
    return set(df.loc[mask, "sample_id"].astype(str))


def load_masld_sample_ids(dataset: str) -> set[str]:
    samplesheet = (
        ROOT
        / "RNA-seq"
        / "patient_RNAseq"
        / "data"
        / "samplesheets"
        / f"{dataset}_samplesheet.csv"
    )
    if not samplesheet.exists():
        return set()
    df = pd.read_csv(samplesheet)
    if "sample" not in df.columns:
        return set()
    if dataset == "GSE135251":
        if "disease" in df.columns:
            mask = df["disease"].astype(str).str.lower().ne("control")
        elif "group_in_paper" in df.columns:
            mask = df["group_in_paper"].astype(str).str.lower().ne("control")
        else:
            mask = pd.Series(True, index=df.index)
    else:
        if "nafld_activity_score" in df.columns:
            score = pd.to_numeric(df["nafld_activity_score"], errors="coerce")
            mask = score > 0
            if not mask.any():
                mask = pd.Series(True, index=df.index)
        elif "steatosis_grade" in df.columns:
            score = pd.to_numeric(df["steatosis_grade"], errors="coerce")
            mask = score > 0
        else:
            mask = pd.Series(True, index=df.index)
    return set(df.loc[mask, "sample"].astype(str))


def read_featurecounts_counts(path: Path) -> tuple[pd.Series, pd.DataFrame]:
    df = pd.read_csv(path, sep="\t", comment="#")
    cols = {c.lower(): c for c in df.columns}
    gene_col = cols.get("geneid") or cols.get("gene_id") or cols.get("gene")
    length_col = cols.get("length")
    if gene_col is None or length_col is None:
        raise ValueError(f"Missing Geneid/Length columns in {path}")
    meta_cols = {gene_col, length_col}
    for key in ("chr", "start", "end", "strand"):
        col = cols.get(key)
        if col is not None:
            meta_cols.add(col)
    sample_cols = [c for c in df.columns if c not in meta_cols]
    gene_ids = df[gene_col].astype(str)
    lengths = pd.to_numeric(df[length_col], errors="coerce")
    counts = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    counts.index = gene_ids
    lengths.index = gene_ids
    return lengths, counts


def compute_tpm_mean(counts: pd.DataFrame, lengths: pd.Series) -> pd.Series:
    lengths = lengths.dropna()
    lengths = lengths[lengths > 0]
    shared = counts.index.intersection(lengths.index)
    counts = counts.loc[shared]
    lengths = lengths.loc[shared]
    length_kb = lengths / 1000.0
    rpk = counts.div(length_kb, axis=0)
    scale = rpk.sum(axis=0) / 1e6
    tpm = rpk.div(scale, axis=1)
    return tpm.mean(axis=1, skipna=True)


def load_tpm_from_featurecounts(path: Path, sample_ids: Iterable[str] | None) -> pd.Series:
    lengths, counts = read_featurecounts_counts(path)
    counts = _subset_featurecounts_counts(counts, sample_ids)
    return compute_tpm_mean(counts, lengths)


def load_tpm_from_counts_matrix(
    counts_path: Path, length_source: Path, sample_ids: Iterable[str] | None
) -> pd.Series:
    counts_df = pd.read_csv(counts_path, sep="\t")
    gene_col = counts_df.columns[0]
    counts = counts_df.set_index(gene_col)
    counts = counts.apply(pd.to_numeric, errors="coerce")
    counts = _subset_counts_matrix(counts, sample_ids)

    lengths_df = pd.read_csv(length_source, sep="\t", comment="#")
    cols = {c.lower(): c for c in lengths_df.columns}
    gene_col = cols.get("geneid") or cols.get("gene_id") or cols.get("gene")
    length_col = cols.get("length")
    if gene_col is None or length_col is None:
        raise ValueError(f"Missing Geneid/Length columns in {length_source}")
    lengths = pd.to_numeric(lengths_df[length_col], errors="coerce")
    lengths.index = lengths_df[gene_col].astype(str)
    return compute_tpm_mean(counts, lengths)


def load_raw_tpm(label: str) -> pd.Series:
    if label in RAW_PATIENT_DATASETS:
        dataset = RAW_PATIENT_DATASETS[label]
        counts_dir = ROOT / "RNA-seq" / "patient_RNAseq" / "results" / dataset / "counts"
        counts_path = counts_dir / "gene_counts_matrix.txt"
        length_source = next(counts_dir.glob("individual/*_counts.txt"), None)
        if length_source is None:
            raise FileNotFoundError(f"No length source found for {dataset} in {counts_dir}/individual")
        sample_ids = load_masld_sample_ids(dataset)
        return load_tpm_from_counts_matrix(counts_path, length_source, sample_ids)

    if label in RAW_EXTERNAL_MCD:
        source = RAW_EXTERNAL_MCD[label]
        sample_ids = load_mcd_sample_ids(source["metadata"])
        return load_tpm_from_featurecounts(source["counts"], sample_ids)

    if label.startswith("MCD Week"):
        sample_ids = load_mcd_sample_ids(RAW_MCD_METADATA)
        return load_tpm_from_featurecounts(RAW_MCD_COUNTS, sample_ids)

    raise ValueError(f"No raw TPM source configured for label: {label}")


def read_tpm_filtered(path: Path, filter_mode: str = "standard", biotype_map: dict[str, str] = None) -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(path)
    sep = "\t" if path.suffixes[-2] == ".tsv" else ","
    df = pd.read_csv(path, sep=sep)
    
    needed = {"tpm_mean", "padj", "log2FoldChange"}
    if not needed.issubset(df.columns):
        print(f"Warning: Missing filtering columns in {path.name}, returning unfiltered.")
        return pd.to_numeric(df["tpm_mean"], errors="coerce")
    
    mask_sig = (df["padj"] < 0.1) & (df["log2FoldChange"] > 0.8)
    
    if filter_mode == "dual_threshold" and biotype_map:
        cols = {c.lower(): c for c in df.columns}
        gene_col = cols.get("gene_id") or cols.get("gene") or cols.get("gene_symbol")
        
        if not gene_col:
             print(f"Warning: No ID column to map biotypes in {path.name}. Returning standard filter.")
             filtered = df.loc[mask_sig, "tpm_mean"]
        else:
            clean_ids = df[gene_col].astype(str).apply(lambda x: x.split(".")[0])
            biotypes = clean_ids.map(biotype_map).fillna("unknown")
            
            mask_pcg = (biotypes == "protein_coding") & (df["tpm_mean"] >= 1.0)
            mask_lnc = (biotypes == "lncRNA") & (df["tpm_mean"] >= 0.5)
            
            final_mask = mask_sig & (mask_pcg | mask_lnc)
            
            filtered = df.loc[final_mask, "tpm_mean"]
    else:
        filtered = df.loc[mask_sig, "tpm_mean"]

    return pd.to_numeric(filtered, errors="coerce")


def make_ridgeline(
    datasets: list[tuple[str, pd.Series]],
    xmax: float,
    output_path: Path,
) -> None:
    # ... (existing code) ...
    # ...
    # NOTE: Function body continues... I am replacing just the changed parts if possible, but structure prevents small edits easily.
    # I will replace the read function and the main block where default output is defined.
    # But here I am replacing a chunk. Let's do the read function first.
    pass

# ... (skip to main)

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate TPM distribution ridgeline plots.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "streamlit_deg_explorer" / "data",
        help="Directory containing bundled DEG tables with tpm_mean.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "TPM_analysis" / "filtered" / "tpm_distribution_ridgeline_xcap5.png",
        help="Output PNG path.",
    )
    # ...
    
    datasets: list[tuple[str, pd.Series]] = []
    missing = []
    for label, filename in ORDERED_DATASETS:
        if selected_labels is not None and label not in selected_labels:
            continue
        try:
            if args.source == "raw":
                vals = load_raw_tpm(label)
            else:
                path = args.data_dir / filename
                vals = read_tpm_filtered(path, args.filter_mode, biotype_map)
        except Exception as exc:
            missing.append((label, str(exc)))
            continue
        datasets.append((label, vals))
    
    # ...



def make_ridgeline(
    datasets: list[tuple[str, pd.Series]],
    xmax: float,
    output_path: Path,
) -> None:
    raw_labels = [label for label, _ in datasets]
    display_labels = [DISPLAY_LABELS.get(lbl, lbl) for lbl in raw_labels]
    values = [vals.dropna().astype(float) for _, vals in datasets]
    
    # Use log2(TPM+1) for better visualization of skewed data
    # CAP AT 5 per user request
    log_values = [np.minimum(np.log2(vals + 1.0), 5.0) for vals in values]
    
    # Define range for log plot (0 to 5)
    log_xmax = 5.2
    bins = np.linspace(0.0, 5.0, 100)
    
    histograms = []
    max_density = 0.0
    for vals in log_values:
        if vals.empty:
            hist = np.zeros(len(bins) - 1)
        else:
            hist, _ = np.histogram(vals, bins=bins, density=True)
        max_density = max(max_density, float(hist.max()) if hist.size else 0.0)
        histograms.append(hist)

    if max_density <= 0:
        raise RuntimeError("No TPM variation available for plotting.")

    # Increased scale for overlap (Joyplot style)
    scale = 2.0 / max_density
    centers = (bins[:-1] + bins[1:]) / 2.0

    n = len(display_labels)
    # Adjust height
    height = max(6.0, 0.8 * n + 2.0)
    fig, ax = plt.subplots(figsize=(10, height))

    species = ["Mouse" if lbl in MOUSE_DATASETS else "Human" for lbl in raw_labels]
    mouse_color = "#D81B60"  # magenta
    human_color = "#FFB300"  # amber / orange-yellow

    offsets = list(range(n - 1, -1, -1))

    for idx, (label, hist, sp) in enumerate(zip(display_labels, histograms, species)):
        offset = offsets[idx]
        color = mouse_color if sp == "Mouse" else human_color
        ridge = hist * scale
        
        # Fill with alpha to see overlap
        ax.fill_between(centers, offset, ridge + offset, color=color, alpha=0.8, zorder=n-idx)
        # White outline for separation
        ax.plot(centers, ridge + offset, color="white", linewidth=1.0, zorder=n-idx+0.1)

    # Add visual indication of truncation
    ax.axvline(x=5.0, color="grey", linestyle="--", alpha=0.7)
    ax.text(5.05, n, "Capped at 5", ha="left", va="bottom", fontstyle="italic", color="grey")

    ax.set_xlim(0, log_xmax)
    ax.set_ylim(-0.2, n + 1.0)
    ax.set_yticks(offsets)
    ax.set_yticklabels(display_labels)
    ax.set_xlabel("log2(TPM + 1)")
    ax.set_ylabel("")
    ax.set_title("TPM distribution (ridgeline) - Capped at 5")
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=mouse_color, edgecolor="white", label="Mouse"),
        Patch(facecolor=human_color, edgecolor="white", label="Human"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", framealpha=0.9)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    # Save summary stats
    summary_path = output_path.with_suffix(".stats.tsv")
    summary_df = pd.DataFrame(
        {
            "dataset": display_labels,
            "species": species,
            "num_genes": [int(v.size) for v in values],
            "median_tpm": [float(v.median()) for v in values],
            "mean_tpm": [float(v.mean()) for v in values],
        }
    )
    summary_df.to_csv(summary_path, sep="\t", index=False)


def _default_violin_output(ridge_output: Path) -> Path:
    stem = ridge_output.stem
    if "ridgeline" in stem:
        stem = stem.replace("ridgeline", "violin_log2")
    else:
        stem = f"{stem}_violin_log2"
    return ridge_output.with_name(stem + ridge_output.suffix)


DISPLAY_LABELS = {
    "GSE130970 NAS high": "Human (Hoang)",
    "GSE130970 NAS low": "Human (Hoang) Low",
    "GSE130970 Fibrosis": "Human (Hoang) Fibrosis",
    "GSE135251 NAS high": "Human (Govaere)",
    "GSE135251 NAS low": "Human (Govaere) Low",
    "GSE135251 Fibrosis": "Human (Govaere) Fibrosis",
    "MCD Week pooled (combined)": "Cas13 mouse MCD",
    "GSE156918 (external MCD)": "Mouse MCD (Paquette)",
    "GSE205974 (external MCD)": "Mouse MCD (Yue)",
}

MOUSE_DATASETS = {
    "MCD Week 1",
    "MCD Week 2",
    "MCD Week 3",
    "MCD Week pooled (combined)",
    "GSE156918 (external MCD)",
    "GSE205974 (external MCD)",
}


def make_violin_plot(
    datasets: list[tuple[str, pd.Series]],
    output_path: Path,
) -> None:
    raw_labels = [label for label, _ in datasets]
    display_labels = [DISPLAY_LABELS.get(lbl, lbl) for lbl in raw_labels]
    values = [vals.dropna().astype(float) for _, vals in datasets]
    log2_values = [np.log2(vals.clip(lower=0) + 1.0) for vals in values]

    species = ["Mouse" if lbl in MOUSE_DATASETS else "Human" for lbl in raw_labels]

    n = len(display_labels)
    height = max(7.5, 0.9 * n + 4.0)
    fig, ax = plt.subplots(figsize=(10, height))

    positions = np.arange(1, n + 1)
    width = 0.9

    mouse_color = "#D81B60"  # magenta
    human_color = "#FFB300"  # amber / orange-yellow

    if log2_values:
        flat = np.concatenate([vals.values if hasattr(vals, "values") else np.array(vals) for vals in log2_values])
        if flat.size:
            y_min = float(np.nanmin(flat)) - 0.5
        else:
            y_min = 0.0
    else:
        y_min = 0.0

    # Draw violins one at a time to color by species
    for idx, (pos, data, sp) in enumerate(zip(positions, log2_values, species)):
        color = mouse_color if sp == "Mouse" else human_color
        vp = ax.violinplot(
            [data],
            positions=[pos],
            showmeans=False,
            showmedians=False,
            showextrema=False,
            widths=width,
        )
        for body in vp["bodies"]:
            body.set_facecolor(color)
            body.set_edgecolor("white")
            body.set_alpha(0.7)

    # Overlay boxplots
    box = ax.boxplot(
        log2_values,
        positions=positions,
        widths=0.22,
        patch_artist=True,
        showfliers=False,
    )
    for idx, patch in enumerate(box["boxes"]):
        color = mouse_color if species[idx] == "Mouse" else human_color
        patch.set_facecolor(color)
        patch.set_alpha(0.9)
        patch.set_edgecolor("black")
    for element in ("medians", "whiskers", "caps"):
        for item in box[element]:
            item.set_color("black")

    ax.set_ylim(y_min, 10)

    ax.set_ylabel("log2(TPM + 1)")
    ax.set_title("TPM distribution (log2 TPM + 1)")
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.grid(axis="y", alpha=0.25)
    ax.set_xticks(positions)
    ax.set_xticklabels(display_labels, rotation=30, ha="right")
    ax.set_xlabel("Dataset")

    # Legend for species
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=mouse_color, edgecolor="black", label="Mouse"),
        Patch(facecolor=human_color, edgecolor="black", label="Human"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", framealpha=0.9)

    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.18)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate TPM distribution ridgeline plots.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "streamlit_deg_explorer" / "data",
        help="Directory containing bundled DEG tables with tpm_mean.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "TPM_analysis" / "tpm_distribution_ridgeline_xcap5.png",
        help="Output PNG path.",
    )
    parser.add_argument("--xmax", type=float, default=5.0, help="X-axis cap for TPM values.")
    parser.add_argument(
        "--source",
        choices=["bundled", "raw"],
        default="bundled",
        help="TPM source: bundled tpm_mean columns or raw counts (disease-only samples).",
    )
    parser.add_argument(
        "--violin-output",
        type=Path,
        default=None,
        help="Output PNG path for log2(TPM+1) violin plot (defaults to a name based on --output).",
    )
    parser.add_argument(
        "--no-violin",
        action="store_true",
        help="Skip generating the log2(TPM+1) violin plot.",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        help="Optional list of dataset labels or aliases to include.",
    )
    parser.add_argument(
        "--filter-mode",
        choices=["standard", "dual_threshold"],
        default="standard",
        help="Filtering mode. dual_threshold: PCG>=1.0, lncRNA>=0.5.",
    )
    args = parser.parse_args()

    biotype_map = {}
    if args.filter_mode == "dual_threshold":
        print(f"Loading biotype map from {BIOTYPE_MAP_PATH}...")
        biotype_map = load_biotype_map(BIOTYPE_MAP_PATH)

    allowed_labels = {label for label, _ in ORDERED_DATASETS}
    selected_labels = None
    if args.only:
        normalized = []
        for label in args.only:
            key = ALIAS_LABELS.get(label, label)
            if key not in allowed_labels:
                raise SystemExit(f"Unknown dataset label: {label}")
            normalized.append(key)
        selected_labels = set(normalized)

    datasets: list[tuple[str, pd.Series]] = []
    missing = []
    for label, filename in ORDERED_DATASETS:
        if selected_labels is not None and label not in selected_labels:
            continue
        try:
            if args.source == "raw":
                vals = load_raw_tpm(label)
            else:
                path = args.data_dir / filename
                vals = read_tpm_filtered(path, args.filter_mode, biotype_map)
        except Exception as exc:
            missing.append((label, str(exc)))
            continue
        datasets.append((label, vals))

    if not datasets:
        raise SystemExit("No datasets with TPM values found.")

    make_ridgeline(datasets, args.xmax, args.output)
    if not args.no_violin:
        violin_output = args.violin_output or _default_violin_output(args.output)
        make_violin_plot(datasets, violin_output)

    if missing:
        print("Skipped datasets:")
        for label, reason in missing:
            print(f"- {label}: {reason}")


if __name__ == "__main__":
    main()
