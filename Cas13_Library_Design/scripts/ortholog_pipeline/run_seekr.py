#!/usr/bin/env python
"""
run_seekr.py — L6 Seekr k-mer profile layer for cross-species lncRNA orthology.

Seekr (Kirk et al. 2018 Mol Cell) compares lncRNAs by k-mer frequency profiles,
capturing "compositional orthology" for alignment-dark lncRNAs where BLAST/MMseqs2
fail. At 90 Myr human-mouse divergence, many orthologous lncRNAs have lost
contiguous sequence but retained short motif composition.

Strategy:
  1. Collapse transcripts to gene-level (longest transcript per gene)
  2. Write temporary gene-level FASTAs
  3. Run Seekr k-mer counting at k=5 and k=6
  4. Compute all-vs-all Pearson correlation (batched for memory safety)
  5. Extract top-N pairs per gene, identify RBH (reciprocal best hit) pairs
  6. Emit L6_seekr_kmer.tsv in the common ortholog-table schema
  7. Generate diagnostic plots

Inputs:
  Mouse FASTA: blast_work/mouse_vM37_lncRNA.fa
  Human FASTA: blast_work/human_v47_lncRNA.fa
  Header TSVs: blast_work/{mouse,human}_*_headers.tsv

Outputs:
  data/external/orthologs/layers/L6_seekr_kmer.tsv
  results/ortholog_validation/figs/seekr_*.pdf

Environment: .envs/seekr_env (seekr 2.0.2, numpy<2, pandas<2.1)
"""
from __future__ import annotations

import argparse
import gc
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
BLAST_WORK = (PROJECT / "Cas13_Library_Design/scripts"
              "/ortholog_pipeline/blast_work")
LAYERS_DIR = PROJECT / "data/external/orthologs/layers"
FIG_DIR = (PROJECT / "Cas13_Library_Design/results"
           "/ortholog_validation/figs")

MOUSE_FA = BLAST_WORK / "mouse_vM37_lncRNA.fa"
HUMAN_FA = BLAST_WORK / "human_v47_lncRNA.fa"
MOUSE_HDR = BLAST_WORK / "mouse_vM37_lncRNA_headers.tsv"
HUMAN_HDR = BLAST_WORK / "human_v47_lncRNA_headers.tsv"

OUT_LAYER = LAYERS_DIR / "L6_seekr_kmer.tsv"

# Canonical lncRNA pairs for spot-check
CANONICAL_PAIRS = {
    "MEG3":    ("Meg3",    "MEG3"),
    "MALAT1":  ("Malat1",  "MALAT1"),
    "NEAT1":   ("Neat1",   "NEAT1"),
    "HOTAIR":  ("Hotair",  "HOTAIR"),
    "XIST":    ("Xist",    "XIST"),
    "H19":     ("H19",     "H19"),
}

# Thresholds
TOP_N = 5           # Top-N per gene direction
RBH_TIER_M_MIN = 0.8   # Tier M: RBH + r > 0.8
TIER_L_MIN = 0.5       # Tier L: r > 0.5


def strip_version(s: str) -> str:
    """ENSMUSG00000089699.3 -> ENSMUSG00000089699"""
    return s.split(".")[0]


# ---------------------------------------------------------------------------
# Step 1: Gene-level representative selection
# ---------------------------------------------------------------------------
def select_gene_representatives(hdr_path: Path) -> pd.DataFrame:
    """Pick the longest transcript per gene from the header TSV.

    Returns DataFrame with columns: tx_id, gene_id, gene_id_base, gene_name,
    length, sorted by gene_id_base.
    """
    df = pd.read_csv(hdr_path, sep="\t", dtype=str)
    df["length"] = df["length"].astype(int)
    df["gene_id_base"] = df["gene_id"].apply(strip_version)
    # Pick longest transcript per gene
    df = (df.sort_values("length", ascending=False)
            .drop_duplicates(subset="gene_id_base", keep="first")
            .sort_values("gene_id_base")
            .reset_index(drop=True))
    return df


def write_gene_fasta(reps: pd.DataFrame, full_fa: Path, out_fa: Path) -> None:
    """Extract gene-level representative sequences from the full FASTA.

    Reads the full FASTA, keeps only sequences whose header matches
    the tx_id in reps.  Writes a new FASTA with header = gene_id_base.
    """
    wanted = set(reps["tx_id"].tolist())
    # Build tx_id -> gene_id_base map
    tx2gene = dict(zip(reps["tx_id"], reps["gene_id_base"]))

    keep = False
    current_tx = None
    with open(full_fa) as fin, open(out_fa, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                tx_id = line[1:].strip().split()[0]  # Just the ID
                if tx_id in wanted:
                    keep = True
                    current_tx = tx_id
                    fout.write(f">{tx2gene[tx_id]}\n")
                else:
                    keep = False
            else:
                if keep:
                    fout.write(line)
    print(f"  Wrote {len(wanted)} gene-level seqs -> {out_fa}")


# ---------------------------------------------------------------------------
# Step 2: Seekr k-mer counting
# ---------------------------------------------------------------------------
def count_kmers(fa_path: str, k: int) -> np.ndarray:
    """Run Seekr BasicCounter on a FASTA file. Returns counts array."""
    from seekr.kmer_counts import BasicCounter
    counter = BasicCounter(
        infasta=fa_path,
        k=k,
        binary=True,
        mean=True,
        std=True,
        log2="Log2.post",
        silent=False,
    )
    counter.get_counts()
    return counter.counts


def get_headers_from_fasta(fa_path: str) -> list[str]:
    """Extract headers (gene IDs) from the gene-level FASTA, in order."""
    from seekr.kmer_counts import Reader
    reader = Reader(infasta=fa_path)
    headers = reader.get_headers()
    # Headers look like ">ENSMUSG00000089699"
    return [h.lstrip(">").strip() for h in headers]


# ---------------------------------------------------------------------------
# Step 3: Batched Pearson correlation
# ---------------------------------------------------------------------------
def batched_pearson(counts1: np.ndarray, counts2: np.ndarray,
                    batch_size: int = 5000) -> np.ndarray:
    """Compute Pearson correlation between counts1 (n1, d) and counts2 (n2, d).

    Uses seekr's row-standardization then batched inner product to manage
    memory. Returns (n1, n2) float32 matrix.
    """
    n1, d = counts1.shape
    n2 = counts2.shape[0]

    # Row-standardize (seekr convention)
    c1 = counts1.copy()
    c2 = counts2.copy()
    c1 = (c1.T - np.mean(c1, axis=1)).T
    stds1 = np.std(c1, axis=1)
    stds1[stds1 == 0] = 1.0  # avoid division by zero
    c1 = (c1.T / stds1).T

    c2 = (c2.T - np.mean(c2, axis=1)).T
    stds2 = np.std(c2, axis=1)
    stds2[stds2 == 0] = 1.0
    c2 = (c2.T / stds2).T

    # Batched inner product
    result = np.zeros((n1, n2), dtype=np.float32)
    for i_start in range(0, n1, batch_size):
        i_end = min(i_start + batch_size, n1)
        batch = np.dot(c1[i_start:i_end], c2.T) / d
        result[i_start:i_end] = batch.astype(np.float32)
        if (i_start // batch_size) % 5 == 0:
            print(f"    Batch {i_start}:{i_end} / {n1}", flush=True)

    return result


# ---------------------------------------------------------------------------
# Step 4: Extract top pairs + RBH
# ---------------------------------------------------------------------------
def extract_top_pairs(sim: np.ndarray,
                      mouse_ids: list[str],
                      human_ids: list[str],
                      top_n: int = 5) -> pd.DataFrame:
    """For each mouse gene, find top-N human genes; for each human, top-N mouse.

    Returns deduplicated DataFrame of (mouse_id, human_id, pearson_r,
    mouse_rank, human_rank, is_rbh).
    """
    n_mouse, n_human = sim.shape
    print(f"  Extracting top-{top_n} pairs from {n_mouse}x{n_human} matrix...")

    # Mouse -> Human: top-N per mouse gene
    rows = []
    # Use argpartition for efficiency (partial sort)
    for i in range(n_mouse):
        row = sim[i]
        if top_n < n_human:
            top_idx = np.argpartition(row, -top_n)[-top_n:]
            top_idx = top_idx[np.argsort(row[top_idx])[::-1]]
        else:
            top_idx = np.argsort(row)[::-1][:top_n]
        for rank, j in enumerate(top_idx, 1):
            rows.append((mouse_ids[i], human_ids[j], float(row[j]), rank))

    m2h = pd.DataFrame(rows, columns=["mouse_id", "human_id", "r", "mouse_rank"])

    # Human -> Mouse: top-N per human gene
    rows = []
    for j in range(n_human):
        col = sim[:, j]
        if top_n < n_mouse:
            top_idx = np.argpartition(col, -top_n)[-top_n:]
            top_idx = top_idx[np.argsort(col[top_idx])[::-1]]
        else:
            top_idx = np.argsort(col)[::-1][:top_n]
        for rank, i in enumerate(top_idx, 1):
            rows.append((mouse_ids[i], human_ids[j], float(col[i]), rank))

    h2m = pd.DataFrame(rows, columns=["mouse_id", "human_id", "r", "human_rank"])

    # Merge to get both ranks
    merged = m2h.merge(h2m, on=["mouse_id", "human_id", "r"], how="outer")
    merged["mouse_rank"] = merged["mouse_rank"].fillna(999).astype(int)
    merged["human_rank"] = merged["human_rank"].fillna(999).astype(int)

    # RBH: mutual top-1
    merged["is_rbh"] = (merged["mouse_rank"] == 1) & (merged["human_rank"] == 1)

    # Deduplicate
    merged = (merged.sort_values(["mouse_id", "r"], ascending=[True, False])
                    .drop_duplicates(subset=["mouse_id", "human_id"])
                    .reset_index(drop=True))

    print(f"  Total pairs (union of top-{top_n} both directions): {len(merged)}")
    print(f"  RBH pairs (mutual top-1): {merged['is_rbh'].sum()}")

    return merged


# ---------------------------------------------------------------------------
# Step 5: Build L6 output
# ---------------------------------------------------------------------------
def build_layer(pairs_k5: pd.DataFrame, pairs_k6: pd.DataFrame,
                mouse_reps: pd.DataFrame, human_reps: pd.DataFrame,
                ) -> pd.DataFrame:
    """Merge k=5 and k=6 results, annotate with gene metadata, assign tiers."""

    # Rename r columns
    k5 = pairs_k5.rename(columns={
        "r": "seekr_pearson_k5",
        "is_rbh": "is_rbh_k5",
        "mouse_rank": "mouse_rank_k5",
        "human_rank": "human_rank_k5",
    })
    k6 = pairs_k6.rename(columns={
        "r": "seekr_pearson_k6",
        "is_rbh": "is_rbh_k6",
        "mouse_rank": "mouse_rank_k6",
        "human_rank": "human_rank_k6",
    })

    # Merge on (mouse_id, human_id)
    merged = k5.merge(k6, on=["mouse_id", "human_id"], how="outer")

    # Composite Seekr score: max of k5, k6
    merged["seekr_pearson_max"] = merged[
        ["seekr_pearson_k5", "seekr_pearson_k6"]
    ].max(axis=1)

    # RBH: mutual top-1 in EITHER k
    merged["seekr_is_rbh"] = (
        merged["is_rbh_k5"].fillna(False) | merged["is_rbh_k6"].fillna(False)
    )

    # Filter: keep pairs with max r > TIER_L_MIN
    merged = merged[merged["seekr_pearson_max"] >= TIER_L_MIN].copy()
    print(f"\n  After r >= {TIER_L_MIN} filter: {len(merged)} pairs")

    # Provenance
    has_k5 = merged["seekr_pearson_k5"].notna()
    has_k6 = merged["seekr_pearson_k6"].notna()
    merged["provenance_sources"] = "seekr_k5+k6"
    merged.loc[has_k5 & ~has_k6, "provenance_sources"] = "seekr_k5"
    merged.loc[~has_k5 & has_k6, "provenance_sources"] = "seekr_k6"

    # Tier assignment
    #   M: r > 0.8 AND is_rbh
    #   L: 0.5 <= r < 0.8, or r > 0.8 without RBH
    merged["tier_M_seekr"] = (
        (merged["seekr_pearson_max"] >= RBH_TIER_M_MIN)
        & merged["seekr_is_rbh"]
    ).astype(int)
    merged["tier_L_seekr"] = (
        (merged["tier_M_seekr"] == 0)
        & (merged["seekr_pearson_max"] >= TIER_L_MIN)
    ).astype(int)
    merged["confidence_tier"] = np.where(merged["tier_M_seekr"] == 1, "M", "L")

    # Add gene metadata
    mouse_meta = mouse_reps[["gene_id_base", "gene_name"]].rename(columns={
        "gene_id_base": "mouse_id", "gene_name": "mouse_symbol",
    })
    human_meta = human_reps[["gene_id_base", "gene_name"]].rename(columns={
        "gene_id_base": "human_id", "gene_name": "human_symbol",
    })
    merged = merged.merge(mouse_meta, on="mouse_id", how="left")
    merged = merged.merge(human_meta, on="human_id", how="left")

    # Rename to common schema
    merged = merged.rename(columns={
        "mouse_id": "mouse_ensembl",
        "human_id": "human_ensembl",
    })

    # Add biotype (all lncRNA from the FASTA source)
    merged["mouse_biotype"] = "lncRNA"
    merged["human_biotype"] = "lncRNA"

    # Select output columns
    out_cols = [
        "mouse_ensembl", "mouse_symbol", "mouse_biotype",
        "human_ensembl", "human_symbol", "human_biotype",
        "tier_M_seekr", "tier_L_seekr",
        "seekr_pearson_k5", "seekr_pearson_k6", "seekr_pearson_max",
        "seekr_is_rbh", "confidence_tier", "provenance_sources",
    ]
    out = merged[out_cols].sort_values(
        ["seekr_pearson_max", "mouse_ensembl"],
        ascending=[False, True],
    ).reset_index(drop=True)

    return out


# ---------------------------------------------------------------------------
# Step 6: Diagnostic plots
# ---------------------------------------------------------------------------
def make_plots(pairs_k5: pd.DataFrame, pairs_k6: pd.DataFrame,
               layer: pd.DataFrame, blast_rbh_path: Path) -> None:
    """Generate diagnostic PDFs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # --- Plot 1: Distribution of Pearson r values ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, (label, df) in zip(axes, [("k=5", pairs_k5), ("k=6", pairs_k6)]):
        # Only plot top-1 pairs (one per mouse gene) to avoid inflating counts
        top1 = df[df["mouse_rank"] == 1]["r"].values
        ax.hist(top1, bins=100, color="#4A90D9", alpha=0.8, edgecolor="none")
        ax.set_xlabel("Seekr Pearson r (top-1 match per mouse gene)")
        ax.set_ylabel("Count")
        ax.set_title(f"Seekr {label}: top-1 Pearson r distribution")
        ax.axvline(0.5, color="orange", ls="--", lw=1.2, label="r=0.5")
        ax.axvline(0.8, color="red", ls="--", lw=1.2, label="r=0.8")
        ax.legend()
        # Annotation: how many above thresholds
        n_05 = (top1 >= 0.5).sum()
        n_08 = (top1 >= 0.8).sum()
        ax.text(0.95, 0.95,
                f"r >= 0.5: {n_05:,}\nr >= 0.8: {n_08:,}\ntotal: {len(top1):,}",
                transform=ax.transAxes, va="top", ha="right",
                fontsize=9, bbox=dict(facecolor="white", alpha=0.8))

    plt.tight_layout()
    plt.savefig(FIG_DIR / "seekr_pearson_distribution.pdf", dpi=150)
    plt.close()
    print(f"  Saved {FIG_DIR / 'seekr_pearson_distribution.pdf'}")

    # --- Plot 2: Seekr r vs BLAST pident for overlapping pairs ---
    if blast_rbh_path.exists():
        blast = pd.read_csv(blast_rbh_path, sep="\t", dtype=str)
        blast["blast_pident"] = pd.to_numeric(blast["blast_pident"],
                                              errors="coerce")
        blast_pairs = blast[["mouse_ensembl", "human_ensembl",
                             "blast_pident"]].dropna()

        seekr_pairs = layer[["mouse_ensembl", "human_ensembl",
                             "seekr_pearson_k5", "seekr_pearson_k6",
                             "seekr_pearson_max"]].copy()

        overlap = blast_pairs.merge(seekr_pairs,
                                    on=["mouse_ensembl", "human_ensembl"])

        if len(overlap) > 0:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            for ax, (k_label, col) in zip(axes, [
                ("k=5", "seekr_pearson_k5"), ("k=6", "seekr_pearson_k6")
            ]):
                valid = overlap.dropna(subset=[col])
                ax.scatter(valid["blast_pident"], valid[col],
                           alpha=0.3, s=8, c="#4A90D9", edgecolors="none")
                ax.set_xlabel("BLAST percent identity")
                ax.set_ylabel(f"Seekr Pearson r ({k_label})")
                ax.set_title(f"BLAST pident vs Seekr r ({k_label})\n"
                             f"n={len(valid):,} overlapping pairs")
                # Pearson correlation annotation
                from scipy.stats import pearsonr, spearmanr
                if len(valid) > 2:
                    rp, pp = pearsonr(valid["blast_pident"], valid[col])
                    rs, ps = spearmanr(valid["blast_pident"], valid[col])
                    ax.text(0.05, 0.95,
                            f"Pearson r={rp:.3f}\nSpearman rho={rs:.3f}",
                            transform=ax.transAxes, va="top", ha="left",
                            fontsize=9,
                            bbox=dict(facecolor="white", alpha=0.8))

            plt.tight_layout()
            plt.savefig(FIG_DIR / "seekr_vs_blast_scatter.pdf", dpi=150)
            plt.close()
            print(f"  Saved {FIG_DIR / 'seekr_vs_blast_scatter.pdf'} "
                  f"({len(overlap)} overlapping pairs)")
        else:
            print("  No overlapping pairs between Seekr and BLAST RBH.")
    else:
        print(f"  BLAST RBH file not found at {blast_rbh_path}; "
              f"skipping BLAST comparison plot.")

    # --- Plot 3: k=5 vs k=6 agreement scatter ---
    both = layer.dropna(subset=["seekr_pearson_k5", "seekr_pearson_k6"])
    if len(both) > 100:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(both["seekr_pearson_k5"], both["seekr_pearson_k6"],
                   alpha=0.2, s=6, c="#4A90D9", edgecolors="none")
        ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5)
        ax.set_xlabel("Seekr Pearson r (k=5)")
        ax.set_ylabel("Seekr Pearson r (k=6)")
        ax.set_title(f"k=5 vs k=6 agreement (n={len(both):,} pairs)")
        from scipy.stats import pearsonr
        rp, _ = pearsonr(both["seekr_pearson_k5"], both["seekr_pearson_k6"])
        ax.text(0.05, 0.95, f"Pearson r={rp:.3f}",
                transform=ax.transAxes, va="top", ha="left",
                fontsize=10, bbox=dict(facecolor="white", alpha=0.8))
        plt.tight_layout()
        plt.savefig(FIG_DIR / "seekr_k5_vs_k6.pdf", dpi=150)
        plt.close()
        print(f"  Saved {FIG_DIR / 'seekr_k5_vs_k6.pdf'}")


# ---------------------------------------------------------------------------
# Step 7: Canonical pair spot-check
# ---------------------------------------------------------------------------
def spot_check_canonical(layer: pd.DataFrame) -> str:
    """Report Seekr r for canonical conserved lncRNA pairs."""
    lines = ["\n=== Canonical lncRNA Pair Spot-Check ==="]
    lines.append(f"{'Pair':<20} {'k5':>8} {'k6':>8} {'max':>8} {'RBH':>5}")
    lines.append("-" * 55)

    for label, (mouse_sym, human_sym) in CANONICAL_PAIRS.items():
        match = layer[
            (layer["mouse_symbol"].str.lower() == mouse_sym.lower())
            & (layer["human_symbol"].str.lower() == human_sym.lower())
        ]
        if len(match) == 0:
            # Try relaxed: just check if both symbols appear anywhere
            mouse_hit = layer[
                layer["mouse_symbol"].str.lower() == mouse_sym.lower()
            ]
            human_hit = layer[
                layer["human_symbol"].str.lower() == human_sym.lower()
            ]
            if len(mouse_hit) > 0 and len(human_hit) > 0:
                lines.append(
                    f"{mouse_sym}/{human_sym:<14} "
                    f"  NOT IN TOP-5 (mouse has {len(mouse_hit)} hits, "
                    f"human has {len(human_hit)} hits)")
            else:
                lines.append(
                    f"{mouse_sym}/{human_sym:<14}   NOT FOUND "
                    f"(m={len(mouse_hit)}, h={len(human_hit)})")
        else:
            row = match.iloc[0]
            k5 = f"{row['seekr_pearson_k5']:.4f}" if pd.notna(
                row["seekr_pearson_k5"]) else "N/A"
            k6 = f"{row['seekr_pearson_k6']:.4f}" if pd.notna(
                row["seekr_pearson_k6"]) else "N/A"
            mx = f"{row['seekr_pearson_max']:.4f}"
            rbh = "Yes" if row["seekr_is_rbh"] else "No"
            lines.append(f"{mouse_sym}/{human_sym:<14} {k5:>8} {k6:>8} "
                         f"{mx:>8} {rbh:>5}")

    report = "\n".join(lines)
    print(report)
    return report


# ---------------------------------------------------------------------------
# Step 8: Summary statistics
# ---------------------------------------------------------------------------
def summary_stats(layer: pd.DataFrame, canonical_report: str) -> str:
    """Print and return summary statistics."""
    lines = ["\n=== L6 Seekr K-mer Layer Summary ==="]
    lines.append(f"Total pairs (r >= {TIER_L_MIN}): {len(layer):,}")

    for thresh in [0.5, 0.6, 0.7, 0.8, 0.9]:
        n = (layer["seekr_pearson_max"] >= thresh).sum()
        lines.append(f"  Pairs with max r >= {thresh}: {n:,}")

    n_rbh = layer["seekr_is_rbh"].sum()
    lines.append(f"  RBH pairs (mutual top-1): {n_rbh:,}")

    n_tier_m = (layer["tier_M_seekr"] == 1).sum()
    n_tier_l = (layer["tier_L_seekr"] == 1).sum()
    lines.append(f"  Tier M (r>=0.8 + RBH): {n_tier_m:,}")
    lines.append(f"  Tier L (r>=0.5, non-M): {n_tier_l:,}")

    n_mouse = layer["mouse_ensembl"].nunique()
    n_human = layer["human_ensembl"].nunique()
    lines.append(f"  Unique mouse genes: {n_mouse:,}")
    lines.append(f"  Unique human genes: {n_human:,}")

    lines.append(canonical_report)

    report = "\n".join(lines)
    print(report)
    return report


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="L6 Seekr k-mer orthology")
    parser.add_argument("--k-values", nargs="+", type=int, default=[5, 6],
                        help="k-mer sizes to run (default: 5 6)")
    parser.add_argument("--top-n", type=int, default=TOP_N,
                        help="Top-N matches per gene (default: 5)")
    parser.add_argument("--batch-size", type=int, default=5000,
                        help="Batch size for Pearson computation")
    parser.add_argument("--skip-plots", action="store_true",
                        help="Skip diagnostic plots")
    args = parser.parse_args()

    t0 = time.time()

    # -----------------------------------------------------------------------
    # Step 1: Gene-level representatives
    # -----------------------------------------------------------------------
    print("=" * 70)
    print("Step 1: Selecting gene-level representatives (longest tx per gene)")
    print("=" * 70)

    mouse_reps = select_gene_representatives(MOUSE_HDR)
    human_reps = select_gene_representatives(HUMAN_HDR)
    print(f"  Mouse: {len(mouse_reps):,} genes")
    print(f"  Human: {len(human_reps):,} genes")

    # Write temporary gene-level FASTAs
    tmpdir = tempfile.mkdtemp(prefix="seekr_")
    mouse_gene_fa = os.path.join(tmpdir, "mouse_genes.fa")
    human_gene_fa = os.path.join(tmpdir, "human_genes.fa")

    print("\n  Writing gene-level FASTAs...")
    write_gene_fasta(mouse_reps, MOUSE_FA, Path(mouse_gene_fa))
    write_gene_fasta(human_reps, HUMAN_FA, Path(human_gene_fa))

    # Verify
    n_mouse_written = sum(1 for line in open(mouse_gene_fa) if line.startswith(">"))
    n_human_written = sum(1 for line in open(human_gene_fa) if line.startswith(">"))
    print(f"  Verified: {n_mouse_written} mouse, {n_human_written} human gene seqs")

    # -----------------------------------------------------------------------
    # Step 2 + 3: K-mer counting and correlation for each k
    # -----------------------------------------------------------------------
    all_pairs = {}

    for k in args.k_values:
        print(f"\n{'=' * 70}")
        print(f"Step 2-3: Seekr k={k} (counting + correlation)")
        print(f"{'=' * 70}")

        print(f"\n  Counting {k}-mers for mouse ({len(mouse_reps):,} genes)...")
        t1 = time.time()
        mouse_counts = count_kmers(mouse_gene_fa, k)
        print(f"  Mouse counts shape: {mouse_counts.shape} "
              f"({time.time()-t1:.1f}s)")

        print(f"\n  Counting {k}-mers for human ({len(human_reps):,} genes)...")
        t1 = time.time()
        human_counts = count_kmers(human_gene_fa, k)
        print(f"  Human counts shape: {human_counts.shape} "
              f"({time.time()-t1:.1f}s)")

        # Get gene IDs in FASTA order
        mouse_ids = get_headers_from_fasta(mouse_gene_fa)
        human_ids = get_headers_from_fasta(human_gene_fa)
        assert len(mouse_ids) == mouse_counts.shape[0], \
            f"Mouse ID/count mismatch: {len(mouse_ids)} vs {mouse_counts.shape[0]}"
        assert len(human_ids) == human_counts.shape[0], \
            f"Human ID/count mismatch: {len(human_ids)} vs {human_counts.shape[0]}"

        print(f"\n  Computing Pearson correlation "
              f"({len(mouse_ids)} x {len(human_ids)})...")
        t1 = time.time()
        sim = batched_pearson(mouse_counts, human_counts,
                              batch_size=args.batch_size)
        print(f"  Correlation matrix shape: {sim.shape} "
              f"({time.time()-t1:.1f}s)")
        print(f"  Correlation range: [{sim.min():.4f}, {sim.max():.4f}], "
              f"mean={sim.mean():.4f}, median={np.median(sim):.4f}")

        # Free counts
        del mouse_counts, human_counts
        gc.collect()

        # Extract top pairs
        print(f"\n  Extracting top-{args.top_n} pairs...")
        pairs = extract_top_pairs(sim, mouse_ids, human_ids, top_n=args.top_n)
        all_pairs[k] = pairs

        # Free similarity matrix
        del sim
        gc.collect()

    # -----------------------------------------------------------------------
    # Step 4: Build L6 layer
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("Step 4: Building L6 layer")
    print(f"{'=' * 70}")

    layer = build_layer(
        all_pairs.get(5, pd.DataFrame()),
        all_pairs.get(6, pd.DataFrame()),
        mouse_reps, human_reps,
    )

    # Save
    LAYERS_DIR.mkdir(parents=True, exist_ok=True)
    layer.to_csv(OUT_LAYER, sep="\t", index=False)
    print(f"\n  Saved L6 layer: {OUT_LAYER} ({len(layer):,} rows)")

    # -----------------------------------------------------------------------
    # Step 5: Diagnostic plots
    # -----------------------------------------------------------------------
    if not args.skip_plots:
        print(f"\n{'=' * 70}")
        print("Step 5: Diagnostic plots")
        print(f"{'=' * 70}")
        blast_rbh_path = LAYERS_DIR / "L4_blast_rbh.tsv"
        make_plots(all_pairs.get(5, pd.DataFrame()),
                   all_pairs.get(6, pd.DataFrame()),
                   layer, blast_rbh_path)

    # -----------------------------------------------------------------------
    # Step 6: Canonical pair spot-check
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("Step 6: Canonical pair spot-check")
    print(f"{'=' * 70}")
    canonical_report = spot_check_canonical(layer)

    # -----------------------------------------------------------------------
    # Step 7: Summary
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("Step 7: Summary")
    print(f"{'=' * 70}")
    report = summary_stats(layer, canonical_report)

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed/60:.1f} minutes")

    # Cleanup temp files
    for f in [mouse_gene_fa, human_gene_fa]:
        if os.path.exists(f):
            os.remove(f)
    if os.path.exists(tmpdir):
        os.rmdir(tmpdir)

    return report


if __name__ == "__main__":
    main()
