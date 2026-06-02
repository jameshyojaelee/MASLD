#!/usr/bin/env python
"""
490d_mouse_null_aggregate.py — combine the 10 chunks written by 490c into the
banner-required TSV: validation/mouse/mouse_conservation_null_k16.tsv with
columns
    program, n_strict, observed_score, null_mean, null_sd, null_95p,
    perm_p, perm_q_bh

The observed score used here is the direct-mean version computed by
`490c --mode observed` (`mouse_observed_directmean_k16.tsv`).  We also keep a
diagnostic column with the legacy score_genes value (`obs_score_legacy`) read
from `mouse_projection_strict_k16.tsv` so reviewers can see both.

The null distribution per program is the slice of the merged null with
`program_size == n_strict` for that program.  Two-sided permutation p uses the
+1 smoothing convention (more_extreme + 1) / (n_finite + 1).
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
OUT_DIR = MCP / "validation/mouse"
NULL_CHUNK_DIR = OUT_DIR / "null_chunks"


def main(k: int = 16, score_kind: str = "score_donor_mean"):
    chunks = sorted(NULL_CHUNK_DIR.glob("null_chunk_*.tsv"))
    if not chunks:
        sys.exit(f"[490d] no chunks found in {NULL_CHUNK_DIR}")
    print(f"[490d] merging {len(chunks)} chunks")
    null_df = pd.concat([pd.read_csv(c, sep="\t") for c in chunks], ignore_index=True)
    print(f"[490d] total null rows: {len(null_df)}")

    # Observed direct-mean (matches null metric exactly).
    obs_path = OUT_DIR / f"mouse_observed_directmean_k{k}.tsv"
    if not obs_path.exists():
        sys.exit(f"[490d] missing observed table: {obs_path}; "
                 "run `python 490c_mouse_null_vectorized.py --mode observed`")
    obs_df = pd.read_csv(obs_path, sep="\t").rename(columns={
        "obs_score_cell_mean": "observed_score_cell_mean",
        "obs_score_donor_mean": "observed_score_donor_mean",
        "obs_score_donor_median": "observed_score_donor_median",
    })

    # Legacy (score_genes) observed for cross-check.
    legacy_path = OUT_DIR / f"mouse_projection_strict_k{k}.tsv"
    legacy_df = pd.read_csv(legacy_path, sep="\t")
    legacy_df["program"] = legacy_df["program"].astype(str)

    # Match score-kind to observed column.
    obs_col_map = {
        "score_cell_mean": "observed_score_cell_mean",
        "score_donor_mean": "observed_score_donor_mean",
        "score_donor_median": "observed_score_donor_median",
    }
    obs_col = obs_col_map[score_kind]

    rows = []
    for _, r in obs_df.iterrows():
        p = str(r["program"])
        sz = int(r["n_strict"])
        observed = float(r[obs_col])
        if sz < 5 or not np.isfinite(observed):
            rows.append({
                "program": p,
                "n_strict": sz,
                "observed_score": observed,
                "observed_score_cell_mean": float(r["observed_score_cell_mean"]),
                "observed_score_donor_mean": float(r["observed_score_donor_mean"]),
                "observed_score_donor_median": float(r["observed_score_donor_median"]),
                "null_mean": np.nan, "null_sd": np.nan, "null_95p": np.nan,
                "perm_p": np.nan, "n_perm_finite": 0, "score_kind": score_kind,
            })
            continue
        sub = null_df[(null_df["program_size"] == sz) & null_df[score_kind].notna()]
        nulls = sub[score_kind].astype(float).values
        n_finite = len(nulls)
        if n_finite < 10:
            rows.append({
                "program": p, "n_strict": sz,
                "observed_score": observed,
                "observed_score_cell_mean": float(r["observed_score_cell_mean"]),
                "observed_score_donor_mean": float(r["observed_score_donor_mean"]),
                "observed_score_donor_median": float(r["observed_score_donor_median"]),
                "null_mean": np.nan, "null_sd": np.nan, "null_95p": np.nan,
                "perm_p": np.nan, "n_perm_finite": n_finite, "score_kind": score_kind,
            })
            continue
        nm = float(nulls.mean())
        ns = float(nulls.std(ddof=1)) if n_finite > 1 else np.nan
        n95 = float(np.quantile(nulls, 0.95))
        more_extreme = (np.abs(nulls - nm) >= abs(observed - nm)).sum()
        p_perm = (more_extreme + 1) / (n_finite + 1)
        rows.append({
            "program": p, "n_strict": sz,
            "observed_score": observed,
            "observed_score_cell_mean": float(r["observed_score_cell_mean"]),
            "observed_score_donor_mean": float(r["observed_score_donor_mean"]),
            "observed_score_donor_median": float(r["observed_score_donor_median"]),
            "null_mean": nm, "null_sd": ns, "null_95p": n95,
            "perm_p": p_perm, "n_perm_finite": n_finite, "score_kind": score_kind,
        })

    df = pd.DataFrame(rows).sort_values("program", key=lambda s: s.astype(float).astype(int))

    # Add legacy score_genes obs (diagnostic).
    legacy_obs_lookup = {}
    # mouse_projection_strict_k16.tsv only has size info, not actual obs score.  We
    # pull the legacy obs from mouse_projection_strict_k16_mean_by_sample.tsv:
    # column means there are per-sample scores → take overall mean across samples
    # to match the cell-averaged metric in 490b's null_mode (which used the
    # full-cells mean).  Cell-averaged is better captured by our direct-mean
    # numbers above; this column is purely diagnostic.
    by_sample_path = OUT_DIR / f"mouse_projection_strict_k{k}_mean_by_sample.tsv"
    if by_sample_path.exists():
        bs = pd.read_csv(by_sample_path, sep="\t")
        for col in bs.columns:
            if col.startswith("prog_") and col.endswith("_strict"):
                pid = col.replace("prog_", "").replace("_strict", "")
                legacy_obs_lookup[pid] = float(bs[col].mean())
    df["obs_score_legacy_score_genes"] = df["program"].map(legacy_obs_lookup)

    # BH FDR.
    pmask = df["perm_p"].notna()
    df["perm_q_bh"] = np.nan
    if pmask.any():
        from scipy.stats import false_discovery_control
        try:
            df.loc[pmask, "perm_q_bh"] = false_discovery_control(
                df.loc[pmask, "perm_p"].values, method="bh"
            )
        except Exception:
            pv = df.loc[pmask, "perm_p"].values
            order = np.argsort(pv)
            ranks = np.empty_like(order); ranks[order] = np.arange(1, len(pv) + 1)
            qv = pv * len(pv) / ranks
            qv = np.minimum.accumulate(qv[order[::-1]])[::-1]
            qfull = np.empty_like(pv); qfull[order] = qv
            df.loc[pmask, "perm_q_bh"] = qfull

    out_f = OUT_DIR / f"mouse_conservation_null_k{k}.tsv"
    df.to_csv(out_f, sep="\t", index=False)
    print(f"[490d] wrote {out_f}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument(
        "--score-kind",
        choices=["score_cell_mean", "score_donor_mean", "score_donor_median"],
        default="score_donor_mean",
        help="which observed/null aggregation to compare. donor_mean is the "
             "level reviewers usually expect (same level as the existing "
             "mouse_projection_strict_k16_mean_by_sample.tsv).",
    )
    args = ap.parse_args()
    main(k=args.k, score_kind=args.score_kind)
