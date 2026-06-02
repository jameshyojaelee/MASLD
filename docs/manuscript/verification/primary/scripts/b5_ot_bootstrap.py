#!/usr/bin/env python
"""
B5 Optimal Transport bootstrap CI (V5 verification).

For each of 4 fibrosis transitions (F0→F1, F1→F2, F2→F3, F3→F4),
bootstrap 1000x (resample donors with replacement), recompute Sinkhorn OT
cost on NAS-VAE embeddings, and report mean ± 95% CI.

Pass criterion: CIs non-overlapping between adjacent stages, especially F2→F3 jump.
Output: verification/controls/b5_ot_bootstrap_ci.csv
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from pathlib import Path

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
STAG = INT / "results/staging_classifier"
OUT_CSV = BASE / "docs/manuscript/verification/controls/b5_ot_bootstrap_ci.csv"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

N_BOOT = 1000
RNG_SEED = 42
REG = 0.05  # Sinkhorn regularization (matches Script 116 line 537)

# ---- Sinkhorn OT (matches Script 116 lines 193-210) ----
def sinkhorn_ot(C, a, b, reg=0.1, max_iter=1000, tol=1e-8):
    K = np.exp(-C / reg)
    u = np.ones_like(a)
    for _ in range(max_iter):
        v = b / (K.T @ u + 1e-10)
        u_new = a / (K @ v + 1e-10)
        if np.max(np.abs(u_new - u)) < tol:
            break
        u = u_new
    return np.diag(u) @ K @ np.diag(v)


def ot_cost(emb_from, emb_to, reg=REG):
    C = cdist(emb_from, emb_to, metric="euclidean")
    n, m = len(emb_from), len(emb_to)
    a = np.ones(n) / n
    b = np.ones(m) / m
    plan = sinkhorn_ot(C, a, b, reg=reg)
    return float(np.sum(plan * C))


def main():
    # ---- Load NAS embeddings + metadata ----
    emb_path = STAG / "nas_embeddings_all_samples.csv"
    meta_path = STAG / "modeling_metadata.csv"
    print(f"Loading NAS embeddings: {emb_path}")
    emb = pd.read_csv(emb_path, index_col=0)
    print(f"  Embeddings: {emb.shape}")

    meta = pd.read_csv(meta_path)
    meta = meta[meta["fibrosis_stage"].isin([0, 1, 2, 3, 4])].copy()
    meta["fibrosis_stage"] = meta["fibrosis_stage"].astype(int)
    print(f"  Metadata (staged): {len(meta)}")

    # intersect with embedding samples
    shared = sorted(set(emb.index) & set(meta["sample_id"]))
    emb = emb.loc[shared]
    meta = meta[meta["sample_id"].isin(shared)].set_index("sample_id").loc[shared]
    print(f"  Shared: {len(shared)}")
    print("Per-stage N:")
    print(meta["fibrosis_stage"].value_counts().sort_index())

    stage_idx = {s: np.where(meta["fibrosis_stage"].values == s)[0] for s in range(5)}

    # ---- Observed OT costs ----
    transitions = [(0, 1), (1, 2), (2, 3), (3, 4)]
    obs_costs = {}
    for a, b in transitions:
        cost = ot_cost(emb.values[stage_idx[a]], emb.values[stage_idx[b]], reg=REG)
        obs_costs[(a, b)] = cost
        print(f"  Observed F{a}->F{b}: n_from={len(stage_idx[a])}, n_to={len(stage_idx[b])}, cost={cost:.4f}")

    # ---- Bootstrap ----
    rng = np.random.default_rng(RNG_SEED)
    boot_costs = {t: np.zeros(N_BOOT) for t in transitions}

    import time
    t0 = time.time()
    for b_iter in range(N_BOOT):
        for a, b in transitions:
            ia = rng.choice(stage_idx[a], size=len(stage_idx[a]), replace=True)
            ib = rng.choice(stage_idx[b], size=len(stage_idx[b]), replace=True)
            cost = ot_cost(emb.values[ia], emb.values[ib], reg=REG)
            boot_costs[(a, b)][b_iter] = cost
        if (b_iter + 1) % 50 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (b_iter + 1) * (N_BOOT - b_iter - 1)
            print(f"  Bootstrap {b_iter+1}/{N_BOOT}  ({elapsed:.0f}s elapsed, ~{eta:.0f}s ETA)")

    # ---- Summarize ----
    rows = []
    for (a, b) in transitions:
        arr = boot_costs[(a, b)]
        rows.append({
            "transition": f"F{a}_to_F{b}",
            "n_from": int(len(stage_idx[a])),
            "n_to": int(len(stage_idx[b])),
            "observed_ot_cost": obs_costs[(a, b)],
            "boot_mean": float(np.mean(arr)),
            "boot_sd": float(np.std(arr, ddof=1)),
            "ci_lo_2.5": float(np.quantile(arr, 0.025)),
            "ci_hi_97.5": float(np.quantile(arr, 0.975)),
            "n_boot": N_BOOT,
        })
    df = pd.DataFrame(rows)

    # Non-overlap check for adjacent transitions
    # Does transition i+1 ci_lo exceed transition i ci_hi?
    df["adjacent_ci_nonoverlap_with_prev"] = [np.nan] + [
        bool(df.loc[i, "ci_lo_2.5"] > df.loc[i-1, "ci_hi_97.5"])
        for i in range(1, len(df))
    ]
    print("\n=== RESULT ===")
    print(df.to_string(index=False))

    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved: {OUT_CSV}")


if __name__ == "__main__":
    main()
