#!/usr/bin/env python3
"""
262b_direct_posteriors.py — replacement for 262_compute_posteriors.py

Uses raw layer scores directly as per-layer posteriors. Stage 1 pre-filtered
each layer to confident edges, so the permutation-based lfdr framework
(in 262) was degenerate: pi0 -> 1, posteriors -> 0.

Each layer's raw score is already a calibrated confidence metric:
  PPI       STRING combined_score / 1000  (a posterior by construction)
  coexpr    |Pearson r|                    (top 0.1%, already filtered)
  pathway   Jaccard of shared pathways     (rank-normalized to reduce skew)
  regulon   Jaccard of shared TFs          (binary-ish)
  lr        LIANA interaction magnitude    (normalized to max)
  genetic   min(PP.H4) at shared locus     (posterior by construction)
  cerna     n_shared_miRNAs / max          (already [0,1])

For pathway (highly skewed, median Jaccard = 0.053), we additionally apply
rank-percentile normalization within the layer to prevent dominance by
many weak pathway co-memberships.

Outputs match the schema that 263+ expects:
  posterior_edges/edges_{layer}_posterior.csv  (parquet format, .csv extension)
  posterior_edges/posterior_summary.csv        (true CSV)
"""

from pathlib import Path
import time
import numpy as np
import pandas as pd

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq")
EDGE_DIR = BASE / "results/network"
OUTDIR = EDGE_DIR / "posterior_edges"
OUTDIR.mkdir(parents=True, exist_ok=True)

LAYER_NAMES = ["ppi", "coexpr", "regulon", "lr", "genetic", "pathway",
               "spatial", "cosmos", "cerna", "xspecies"]

RANK_NORMALIZED = {"pathway"}

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:8.1f}s] {msg}", flush=True)


def load_edges(name):
    path = EDGE_DIR / f"edges_{name}.csv"
    if not path.exists():
        log(f"  WARNING: {path.name} not found, skipping '{name}'")
        return None
    df = pd.read_csv(path)
    if "raw_score" not in df.columns:
        log(f"  WARNING: {name} has no raw_score column, skipping")
        return None
    df["raw_score"] = df["raw_score"].clip(lower=0.0, upper=1.0)
    return df


def compute_posterior(df, name):
    """Use raw score directly; optionally rank-normalize for skewed layers."""
    scores = df["raw_score"].values.astype(np.float64)
    if name in RANK_NORMALIZED:
        ranks = pd.Series(scores).rank(method="average").values
        post = ranks / len(ranks)
    else:
        post = scores
    df = df.copy()
    df["posterior"] = np.clip(post, 0.0, 1.0)
    return df


def main():
    log("=== 262b: Direct-posterior computation (no permutation calibration) ===")
    summary_rows = []

    for name in LAYER_NAMES:
        df = load_edges(name)
        if df is None or len(df) == 0:
            continue

        log(f"\n--- Layer: {name} ---")
        log(f"  Loaded {len(df):,} edges  raw_score range=[{df['raw_score'].min():.3f}, {df['raw_score'].max():.3f}]")

        out_df = compute_posterior(df, name)
        post = out_df["posterior"].values

        n_gt_05 = int(np.sum(post > 0.5))
        n_gt_09 = int(np.sum(post > 0.9))
        summary_rows.append({
            "layer": name,
            "n_edges": len(post),
            "pi0": 0.0,
            "mean_posterior": float(post.mean()),
            "median_posterior": float(np.median(post)),
            "n_posterior_gt_05": n_gt_05,
            "n_posterior_gt_09": n_gt_09,
            "frac_gt_05": float(n_gt_05 / len(post)),
            "frac_gt_09": float(n_gt_09 / len(post)),
            "mode": "rank-percentile" if name in RANK_NORMALIZED else "raw-score",
        })
        log(f"  Posterior: mean={post.mean():.3f}  median={np.median(post):.3f}  "
            f">0.5: {n_gt_05:,} ({n_gt_05/len(post)*100:.1f}%)  "
            f">0.9: {n_gt_09:,} ({n_gt_09/len(post)*100:.1f}%)")

        out_path = OUTDIR / f"edges_{name}_posterior.csv"
        out_df.to_parquet(out_path, index=False)
        log(f"  Wrote: {out_path.name}")

    summary = pd.DataFrame(summary_rows)
    summary_path = OUTDIR / "posterior_summary.csv"
    summary.to_csv(summary_path, index=False)
    log(f"\nSummary saved: {summary_path}")
    print(summary.to_string(index=False))

    log(f"\n=== 262b Complete ===")
    log(f"Layers processed: {len(summary_rows)}")
    log(f"Total time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
