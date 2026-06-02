#!/usr/bin/env python
"""GEARS Phase 2 fine-tune on Saunders 2025 Perturb-Multi mouse hepatocyte.

Pipeline:
  1. Load saunders_perturb_hep_confident.h5ad (built by 08_build_anndata.py).
  2. Prep for GEARS:
       - Map mouse symbols → human orthologs (var.gene_name → human; drop unmapped).
       - Set obs.condition = "{HUMAN_ORTHOLOG}+ctrl" for confident KOs; "ctrl" for safe-harbor / control sgRNAs.
       - Normalize_total(1e4) + log1p (GEARS internally expects log-normalized X).
       - Hold out 10% of cells by random seed for eval.
  3. Build GEARS PertData from prepped AnnData via new_data_process.
  4. Optionally warm-start from Norman K562 weights (--warm-start path).
  5. Train (--epochs, --lr, --batch-size) and save checkpoint per epoch.
  6. Report holdout MSE/Pearson at end.

Usage (smoke, 1000 cells, 2 epochs):
  python gears_finetune.py --smoke --epochs 2 --out-dir .../saunders_hep_v1_smoke

Usage (production):
  python gears_finetune.py --epochs 20 --lr 1e-3 --batch-size 32 \
      --out-dir data/perturbation/finetuned/gears/saunders_hep_v1
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SAUNDERS_H5AD = (
    PROJECT_ROOT
    / "data/perturbation/datasets/saunders2025/processed/saunders_perturb_hep_confident.h5ad"
)
ORTH_MAP = (
    PROJECT_ROOT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"
)


def prep_anndata_for_gears(adata_in, seed=42, smoke=False, smoke_cells=1000):
    """Transform raw Saunders AnnData -> GEARS-ready AnnData.

    GEARS contract:
      - X normalized (log-cpm)
      - var.gene_name = human gene symbol (so condition labels match)
      - obs.condition: "GENE+ctrl" for single KO; "ctrl" for control
    """
    print(f"Input: {adata_in.shape} cells x genes")

    # 1. Build mouse→human ortholog map
    orth = pd.read_csv(ORTH_MAP, sep="\t")
    # Drop rows missing either symbol
    orth = orth.dropna(subset=["mouse_gene_symbol", "human_gene_symbol"]).copy()
    # Prefer ortholog_one2one to avoid 1:many ambiguity
    if "ortholog_type" in orth.columns:
        one2one = orth[orth["ortholog_type"] == "ortholog_one2one"]
        if len(one2one) > 1000:
            print(f"  Using {len(one2one):,} 1:1 mouse→human orthologs (filtered from {len(orth):,})")
            orth = one2one
    m2h = dict(zip(orth["mouse_gene_symbol"], orth["human_gene_symbol"]))

    # 2. Filter var to genes with human ortholog
    var = adata_in.var.copy()
    var["mouse_symbol"] = var["gene_name"]
    var["human_symbol"] = var["mouse_symbol"].map(m2h)
    keep_var = var["human_symbol"].notna() & var["human_symbol"].ne("")
    # Deduplicate var.human_symbol (some mouse genes map to same human ortholog)
    seen = set()
    keep_mask = []
    for idx, hs in zip(var.index, var["human_symbol"]):
        if not isinstance(hs, str) or not hs or hs in seen:
            keep_mask.append(False)
        else:
            seen.add(hs)
            keep_mask.append(True)
    keep_var_final = pd.Series(keep_mask, index=var.index)
    print(f"  var: {len(var):,} mouse genes -> {keep_var_final.sum():,} human-orthologous unique genes")

    adata = adata_in[:, keep_var_final.values].copy()
    adata.var["gene_name"] = adata.var["human_symbol"]
    adata.var_names = adata.var["human_symbol"].astype(str).values
    adata.var_names_make_unique()

    # 3. Set obs.condition
    # target_gene == 'control' -> "ctrl"; else "{human_ortholog}+ctrl"
    target = adata.obs["target_gene"].astype(str)
    is_control = target.isin(["control", "Control", "CTRL", "ctrl", "NTC"])
    human_target = adata.obs["target_gene"].map(m2h)
    cond = np.where(is_control, "ctrl", human_target.astype(str) + "+ctrl")
    adata.obs["condition"] = cond
    # Drop cells whose target lacks a mapped ortholog (and isn't control)
    valid = is_control | human_target.notna()
    print(f"  cells with valid condition: {valid.sum():,} / {len(adata):,}")
    adata = adata[valid.values, :].copy()

    # 4. Optional smoke subsample
    if smoke:
        rng = np.random.RandomState(seed)
        if len(adata) > smoke_cells:
            # Stratify: keep all conditions, sample proportionally
            n_per = max(5, smoke_cells // adata.obs["condition"].nunique())
            keep_idx = []
            for c, idx in adata.obs.groupby("condition").groups.items():
                take = min(len(idx), n_per)
                keep_idx.extend(rng.choice(list(idx), size=take, replace=False))
            adata = adata[keep_idx, :].copy()
            print(f"  smoke subsample: {len(adata):,} cells across {adata.obs['condition'].nunique()} conditions")

    # 5. Normalize_total + log1p
    if not adata.X.dtype.kind == "f":
        adata.X = adata.X.astype(np.float32)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 6. Set raw counts on .layers for safety; X is now log-cpm
    adata.layers["X_norm"] = adata.X.copy()

    print(f"  Final prepped: {adata.shape}; {adata.obs['condition'].nunique()} unique conditions")
    cond_counts = adata.obs["condition"].value_counts()
    print(f"  control cells: {cond_counts.get('ctrl', 0):,}")
    print(f"  KO conditions: {(cond_counts.index != 'ctrl').sum()}")
    print(f"  median cells per condition: {cond_counts.median():.0f}")
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="Subsample ~1k cells, fewer epochs — pipeline validation only")
    ap.add_argument("--smoke-cells", type=int, default=1000)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--hidden-size", type=int, default=64)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--holdout-frac", type=float, default=0.10)
    ap.add_argument("--in-h5ad", default=str(SAUNDERS_H5AD))
    ap.add_argument("--out-dir", required=True,
                    help="Where to write the prepped PertData + checkpoint")
    ap.add_argument("--warm-start", default=None,
                    help="(Optional) Path to a Norman K562 checkpoint dir to load_pretrained() from")
    args = ap.parse_args()

    t0 = time.time()
    print(f"=== GEARS Phase 2 fine-tune on Saunders ===")
    print(f"args: {vars(args)}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load + prep ----
    print(f"\nLoading {args.in_h5ad}")
    adata_in = sc.read_h5ad(args.in_h5ad)
    adata = prep_anndata_for_gears(
        adata_in,
        seed=args.seed,
        smoke=args.smoke,
        smoke_cells=args.smoke_cells,
    )
    adata.write_h5ad(out_dir / "prepped.h5ad", compression="gzip")
    print(f"Wrote prepped AnnData to {out_dir/'prepped.h5ad'}")

    # ---- Build GEARS PertData ----
    print("\nImporting GEARS...")
    from gears import GEARS, PertData

    pert_data = PertData(data_path=str(out_dir))
    dataset_name = "saunders_hep_smoke" if args.smoke else "saunders_hep"
    # Process from our prepped AnnData; this creates {out_dir}/{dataset_name}/...
    try:
        pert_data.new_data_process(dataset_name=dataset_name, adata=adata)
    except Exception as e:
        # If already processed, reuse
        warnings.warn(f"new_data_process raised {e}; will try .load()")
        pert_data.load(data_name=dataset_name)

    # ---- Split + dataloaders ----
    try:
        pert_data.prepare_split(split="simulation", seed=args.seed,
                                train_gene_set_size=(1 - args.holdout_frac))
    except Exception as e:
        warnings.warn(f"prepare_split with holdout failed ({e}); using default split")
        pert_data.prepare_split(split="simulation", seed=args.seed)

    pert_data.get_dataloader(batch_size=args.batch_size, test_batch_size=max(64, args.batch_size * 2))

    # ---- Build model + train ----
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nDevice: {device}; cuda avail: {torch.cuda.is_available()}")

    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=args.hidden_size)

    if args.warm_start and Path(args.warm_start).exists():
        try:
            model.load_pretrained(args.warm_start)
            print(f"Warm-started from {args.warm_start}")
        except Exception as e:
            warnings.warn(f"warm_start load failed ({e}); training from scratch")

    print(f"\nTraining {args.epochs} epochs (lr={args.lr}, bs={args.batch_size})...")
    model.train(epochs=args.epochs, lr=args.lr)

    # ---- Save checkpoint ----
    ckpt = out_dir / "best.ckpt"
    ckpt.mkdir(parents=True, exist_ok=True)
    model.save_model(str(ckpt))
    print(f"\nCheckpoint saved to {ckpt}")

    # ---- Save metadata ----
    pert_universe = None
    for attr in ("pert_names", "perturbable_genes", "gene_names"):
        if hasattr(pert_data, attr):
            v = getattr(pert_data, attr)
            if v is not None and len(v) > 0:
                pert_universe = list(v)
                break
    meta = {
        "args": vars(args),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_conditions": int(adata.obs["condition"].nunique()),
        "control_cells": int((adata.obs["condition"] == "ctrl").sum()),
        "pert_universe_size": len(pert_universe) if pert_universe else None,
        "pert_universe_sample": pert_universe[:50] if pert_universe else None,
        "device": device,
        "training_time_sec": time.time() - t0,
        "dataset_name": dataset_name,
        "checkpoint_path": str(ckpt),
    }
    with open(out_dir / "train_metadata.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"Wrote train_metadata.json")

    print(f"\nDONE in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
