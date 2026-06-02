#!/usr/bin/env python3
"""
262c_confidence_scores.py — final replacement for 262/262b

=============================================================================
METHODS REVISION NOTE (2026-04-16)
=============================================================================
This script REPLACES both 262_compute_posteriors.py and 262b_direct_posteriors.py
in response to adversarial review (docs/network_analysis_adversarial_review.md)
and authored response (docs/network_analysis_response_to_review.md).

WHAT CHANGED vs. the original design spec (Section 1.3):
1. The Efron local-FDR framework (262) was abandoned because Stage 1 pre-filters
   each layer to confident edges, leaving no null component (pi0 -> 1, posteriors
   collapse). Fix: use direct confidence scores, not estimated Bayesian posteriors.
2. Rank-percentile transform (262b) was pathological because it guaranteed median
   pathway score = 0.5 by construction, letting pathway dominate 70-92% of
   composite edges. Fix: tail transform (bottom 90% -> 0, top 10% -> [0,1]).
3. The pathway layer is additionally pre-filtered at Jaccard > 0.15 (top ~5%)
   to prevent MSigDB co-membership from drowning out direct evidence.
4. Antisense-sense lncRNA pairs (GATA6/GATA6-AS1 etc.) are excluded from all
   layers — they score trivially high in regulon/coexpr and don't reflect
   functional biology.
5. Layers producing zero meaningful edges (spatial, cosmos, xspecies) are
   dropped. Reduced from "10-layer multiplex" to "7-layer multiplex."
6. Output column is `confidence_score`, not `posterior`. The summary CSV no
   longer contains a `pi0` column (it was hardcoded to 0.0 in 262b).

RESULT: the composite posterior is now a CONFIDENCE SCORE, not a probability.
Noisy-OR composition still yields values in [0,1], but these should be
interpreted as "aggregated cross-modality confidence", not "P(edge is real)".

HOW TO READ THE OUTPUT:
- `confidence_score` column = normalized score in [0,1] within each layer
- `K_multiplicity` (in composite) = count of layers supporting an edge — the
   primary ranking metric per design Section 1.5
- `p_composite` = 1 - product(1 - c_m) across layers — use AS TIEBREAKER only
=============================================================================
"""

from pathlib import Path
import re
import time
import numpy as np
import pandas as pd

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq")
EDGE_DIR = BASE / "results/network"
OUTDIR = EDGE_DIR / "posterior_edges"
OUTDIR.mkdir(parents=True, exist_ok=True)

ACTIVE_LAYERS = ["ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna"]
DROPPED_LAYERS = ["spatial", "cosmos", "xspecies"]

PATHWAY_JACCARD_MIN = 0.15
TAIL_CUTOFF_QUANTILE = 0.90

REGULON_SCORE_CAP = 0.9

t0 = time.time()
def log(msg):
    print(f"[{time.time()-t0:8.1f}s] {msg}", flush=True)


def is_antisense_pair(a, b):
    """Detect sense-antisense lncRNA pairs like HNF4A/HNF4A-AS1, GATA6/GATA6-AS1."""
    if pd.isna(a) or pd.isna(b):
        return False
    as_pat = re.compile(r"^(.+)-AS\d+$")
    for x, y in [(a, b), (b, a)]:
        m = as_pat.match(x)
        if m and m.group(1) == y:
            return True
    return False


def drop_antisense_pairs(df):
    if "gene_a" not in df.columns or "gene_b" not in df.columns:
        return df
    mask = df.apply(lambda r: is_antisense_pair(r["gene_a"], r["gene_b"]), axis=1)
    n_drop = int(mask.sum())
    if n_drop:
        log(f"  Dropped {n_drop:,} antisense-sense pairs")
    return df[~mask].reset_index(drop=True)


def tail_transform(scores, cutoff_q=TAIL_CUTOFF_QUANTILE):
    """Map bottom cutoff_q fraction -> 0, top (1-cutoff_q) -> [0,1] linearly."""
    scores = np.asarray(scores, dtype=np.float64)
    if len(scores) == 0:
        return scores
    q = float(np.quantile(scores, cutoff_q))
    mx = float(np.max(scores))
    if mx - q <= 1e-12:
        out = (scores >= q).astype(np.float64)
    else:
        out = np.clip((scores - q) / (mx - q), 0.0, 1.0)
    return out


def load_edges(name):
    path = EDGE_DIR / f"edges_{name}.csv"
    if not path.exists():
        log(f"  WARNING: {path.name} not found, skipping '{name}'")
        return None
    df = pd.read_csv(path)
    return df


def compute_confidence(df, name):
    """Layer-specific confidence transform."""
    df = df.copy()
    df["raw_score"] = pd.to_numeric(df["raw_score"], errors="coerce").fillna(0.0).clip(0.0, 1.0)

    if name == "ppi":
        df["confidence_score"] = df["raw_score"].values
        mode = "raw-score"
    elif name == "coexpr":
        df["confidence_score"] = df["raw_score"].values
        mode = "raw-score"
    elif name == "genetic":
        df["confidence_score"] = df["raw_score"].values
        mode = "raw-score"
    elif name == "cerna":
        df["confidence_score"] = df["raw_score"].values
        mode = "raw-score"
    elif name == "lr":
        df["confidence_score"] = df["raw_score"].values
        mode = "raw-score"
    elif name == "regulon":
        df["confidence_score"] = df["raw_score"].clip(upper=REGULON_SCORE_CAP).values
        mode = "raw-score-capped"
    elif name == "pathway":
        mask = df["raw_score"] >= PATHWAY_JACCARD_MIN
        n_pre = len(df)
        df = df[mask].reset_index(drop=True)
        log(f"  Pathway pre-filter at Jaccard >= {PATHWAY_JACCARD_MIN}: {n_pre:,} -> {len(df):,} edges")
        if len(df) > 0:
            df["confidence_score"] = tail_transform(df["raw_score"].values)
        else:
            df["confidence_score"] = df["raw_score"].values
        mode = f"tail-q{TAIL_CUTOFF_QUANTILE}-pre-filtered-J{PATHWAY_JACCARD_MIN}"
    else:
        df["confidence_score"] = df["raw_score"].values
        mode = "raw-score"

    return df, mode


def main():
    log("=== 262c: Confidence Score Computation (post-review) ===")
    log(f"Active layers: {ACTIVE_LAYERS}")
    log(f"Dropped layers (zero meaningful edges): {DROPPED_LAYERS}")

    summary_rows = []

    for name in ACTIVE_LAYERS:
        log(f"\n--- Layer: {name} ---")
        df = load_edges(name)
        if df is None or len(df) == 0:
            log(f"  Skipping (no edges)")
            continue

        df = drop_antisense_pairs(df)

        out_df, mode = compute_confidence(df, name)
        if len(out_df) == 0:
            log(f"  Skipping (empty after filters)")
            continue

        c = out_df["confidence_score"].values
        n_gt_05 = int(np.sum(c > 0.5))
        n_gt_09 = int(np.sum(c > 0.9))

        summary_rows.append({
            "layer": name,
            "mode": mode,
            "n_edges": len(c),
            "min_confidence": float(c.min()),
            "mean_confidence": float(c.mean()),
            "median_confidence": float(np.median(c)),
            "max_confidence": float(c.max()),
            "n_gt_05": n_gt_05,
            "n_gt_09": n_gt_09,
            "frac_gt_05": float(n_gt_05 / len(c)),
            "frac_gt_09": float(n_gt_09 / len(c)),
        })
        log(f"  Confidence: mean={c.mean():.3f}  median={np.median(c):.3f}  "
            f">0.5: {n_gt_05:,} ({n_gt_05/len(c)*100:.1f}%)  "
            f">0.9: {n_gt_09:,} ({n_gt_09/len(c)*100:.1f}%)")

        out_df["posterior"] = out_df["confidence_score"]
        out_path = OUTDIR / f"edges_{name}_posterior.csv"
        out_df.to_parquet(out_path, index=False)
        log(f"  Wrote: {out_path.name}  ({len(out_df):,} edges)")

    summary = pd.DataFrame(summary_rows)
    summary_path = OUTDIR / "posterior_summary.csv"
    summary.to_csv(summary_path, index=False)
    log(f"\nSummary saved: {summary_path}")
    print(summary.to_string(index=False))

    log(f"\n=== 262c Complete ===")
    log(f"Active layers processed: {len(summary_rows)}")
    log(f"Dropped layers (no meaningful edges): {len(DROPPED_LAYERS)}")
    log(f"Total time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
