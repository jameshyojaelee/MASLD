#!/usr/bin/env python
"""
507_scvi_subsample_cnmf.py — Echo Task 20 / A17.

Subsample atlas to 200K cells stratified by (dataset, cell_type, stage), train scVI
(n_latent=30, batch_key=dataset), then export scVI-denoised counts as a cNMF input h5ad.
The companion sbatch wrapper submits the standard cNMF prepare/factorize/consensus chain
on the denoised h5ad with name `global_scvi_denoised`.

After the chain finishes, run with --stage compare to compute top-100 Jaccard between
cNMF-on-integrated (denoised) and cNMF-on-uncorrected (current production global k=16).

Outputs:
  inputs/atlas_cnmf_global_scvi.h5ad           (denoised, ready for cNMF)
  reviewer_defense/integrated_vs_uncorrected_cnmf.tsv   (after compare stage)
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
INPUTS = MCP / "inputs"
CNMF_ROOT = MCP / "cnmf_runs"
RD = MCP / "reviewer_defense"
RD.mkdir(parents=True, exist_ok=True)

ATLAS = INPUTS / "atlas_cnmf_global.h5ad"


def stratified_subsample(adata, n_target: int, seed: int = 42) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    obs = adata.obs.copy()
    strat_cols = []
    for c in ["dataset", "cell_type", "disease_stage_coarse"]:
        if c in obs.columns:
            strat_cols.append(c)
    print(f"[507] strat cols: {strat_cols}")
    obs["_strat"] = obs[strat_cols].astype(str).agg("|".join, axis=1) if strat_cols else "all"
    counts = obs["_strat"].value_counts()
    n_groups = len(counts)
    # Proportional allocation, min 5 per group
    base = np.maximum(5, np.floor(counts.values / counts.sum() * n_target).astype(int))
    # Clip to actual group size
    base = np.minimum(base, counts.values)
    # Adjust to hit n_target (down/upweight)
    deficit = n_target - base.sum()
    if deficit > 0:
        # Distribute deficit to groups with remaining capacity, weighted by remaining
        rem = counts.values - base
        for _ in range(deficit):
            idx = np.argmax(rem)
            if rem[idx] <= 0:
                break
            base[idx] += 1
            rem[idx] -= 1
    elif deficit < 0:
        # Trim from largest allocations
        for _ in range(-deficit):
            idx = np.argmax(base)
            base[idx] -= 1

    keep_idx = []
    by_strat = {s: i for i, s in enumerate(counts.index)}
    for strat, idxs in obs.groupby("_strat", observed=True).groups.items():
        i = by_strat[strat]
        n = int(base[i])
        if n >= len(idxs):
            keep_idx.extend(list(idxs))
        else:
            keep_idx.extend(rng.choice(list(idxs), size=n, replace=False).tolist())
    sub = adata[keep_idx].copy()
    print(f"[507] subsampled: {sub.shape[0]:,} cells across {n_groups} strata")
    return sub


def train_scvi_and_denoise(adata, out_h5ad: Path, n_latent: int = 30, max_epochs: int = 400):
    import torch
    import scvi
    print(f"[507] torch={torch.__version__} cuda={torch.cuda.is_available()}")
    print(f"[507] scvi={scvi.__version__}")

    # scVI expects raw counts in adata.X (or .layers["counts"]). The cNMF norm_counts
    # h5ad we read holds normalized values; instead we use the raw counts that the atlas
    # carries in .layers if available, else assume .X is raw.
    layer_key = None
    for cand in ["counts", "raw_counts", "X_raw"]:
        if cand in adata.layers:
            layer_key = cand
            break
    print(f"[507] using layer={layer_key} for scVI training")

    if "dataset" not in adata.obs.columns:
        raise RuntimeError("dataset column required for batch_key")
    scvi.model.SCVI.setup_anndata(adata, batch_key="dataset", layer=layer_key)
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(max_epochs=max_epochs, early_stopping=True, batch_size=512)
    print("[507] scVI training complete")

    denoised = model.get_normalized_expression(return_mean=True, n_samples=1)
    # denoised is a DataFrame with cells × genes; cNMF expects raw count-like ints
    # We scale to library_size 1e4 and round/cast to integer counts to feed cNMF prep.
    arr = denoised.values
    # Per-cell scale: target library 1e4
    libsize = arr.sum(axis=1, keepdims=True) + 1e-12
    arr = (arr / libsize) * 1e4
    # NOTE: cnmf.prepare internally normalizes by library size and assigns back into
    # adata.X; if X is integer it raises a ufunc-cast error. Keep as float32.
    arr_f32 = np.asarray(arr, dtype=np.float32)

    out = ad.AnnData(
        X=arr_f32,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    out.var_names = adata.var_names
    out.obs_names = adata.obs_names
    out.uns["scvi_denoised"] = {"n_latent": int(n_latent), "max_epochs": int(max_epochs)}
    out.write_h5ad(out_h5ad, compression="gzip")
    print(f"[507] wrote scVI-denoised h5ad: {out_h5ad}")


def stage_subsample_train(args):
    print(f"[507] reading atlas: {ATLAS}")
    a = ad.read_h5ad(ATLAS)
    print(f"[507] atlas: {a.shape}")
    sub = stratified_subsample(a, n_target=args.n_target, seed=args.seed)
    out_h5ad = INPUTS / "atlas_cnmf_global_scvi.h5ad"
    train_scvi_and_denoise(
        sub,
        out_h5ad,
        n_latent=args.n_latent,
        max_epochs=args.max_epochs,
    )


def stage_compare(args):
    """Top-100 Jaccard: cNMF-on-scVI-denoised vs current production cNMF-on-uncorrected."""
    k = args.k
    dt_tag = f"{args.dthresh:.2f}".replace(".", "_")

    full_f = CNMF_ROOT / "global" / f"global.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
    if not full_f.exists():
        full_f = CNMF_ROOT / "global" / "global" / f"global.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"

    name = "global_scvi_denoised"
    int_f = CNMF_ROOT / name / f"{name}.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
    if not int_f.exists():
        int_f = CNMF_ROOT / name / name / f"{name}.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"

    print(f"[507.compare] uncorrected: {full_f}")
    print(f"[507.compare] scvi-denoised: {int_f}")
    full = pd.read_csv(full_f, sep="\t", index_col=0)
    intg = pd.read_csv(int_f, sep="\t", index_col=0)
    top_full = {p: set(full.loc[p].nlargest(100).index) for p in full.index}
    top_int = {p: set(intg.loc[p].nlargest(100).index) for p in intg.index}

    rows = []
    for pf, sf in top_full.items():
        best_j, best_p = 0.0, None
        for pi, si in top_int.items():
            j = len(sf & si) / max(1, len(sf | si))
            if j > best_j:
                best_j, best_p = j, pi
        rows.append({
            "uncorrected_program": pf,
            "best_integrated_match": best_p,
            "jaccard_top100": best_j,
        })
    df = pd.DataFrame(rows)
    out_f = RD / "integrated_vs_uncorrected_cnmf.tsv"
    df.to_csv(out_f, sep="\t", index=False)
    print(f"[507.compare] wrote {out_f}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["train", "compare"], required=True)
    ap.add_argument("--n-target", type=int, default=200000)
    ap.add_argument("--n-latent", type=int, default=30)
    ap.add_argument("--max-epochs", type=int, default=400)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--dthresh", type=float, default=0.03)
    args = ap.parse_args()
    if args.stage == "train":
        stage_subsample_train(args)
    else:
        stage_compare(args)
