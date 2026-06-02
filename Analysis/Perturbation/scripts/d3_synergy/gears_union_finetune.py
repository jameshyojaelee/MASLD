#!/usr/bin/env python
"""GEARS Phase 2 fine-tune — 3-way UNION-VOCAB (per team-lead full Option C, 2026-05-21).

Combines THREE Perturb-seq datasets into a single PertData with a UNION
perturbation vocabulary, all backed by GEARS' `essential_all_data_pert_genes.pkl`
GO-graph universe (9,975 perturbable genes):

  1. Replogle 2022 K562 essential (~163K cells, 1,092 single perts) — broad
     single-perturbation essentiality vocab.
  2. Norman 2019 K562 combinatorial (~91K cells, ~280 conditions ≈ 105 unique
     perts × singles + combos) — combinatorial training signal that drives the
     GI_predict() head.
  3. Saunders 2025 mouse hepatocyte (~111K confident cells, 181 1:1 human-orth
     perts) — hepatocyte cell-context embedding refinement.

Result:
  - Perturbable predict() universe: 9,975 genes (GEARS GO graph) — supports
    tier1 + tier2 D3 hits via GO embedding extrapolation.
  - Directly-trained perturbations: Replogle (1092) ∪ Norman (~105) ∪ Saunders
    (181) ≈ ~1,300 unique perts; ~280 of these have combinatorial training
    signal from Norman (drives GI_predict calibration).
  - `gene_names` (expression universe) = intersection of all three var.gene_name
    sets (~3K-5K).

Saunders Mouse → Human ortholog mapping uses 1:1 orthologs only from
`data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz`.

Usage (smoke):
  python gears_union_finetune.py --smoke --epochs 2 \\
      --out-dir data/perturbation/finetuned/gears/saunders_norman_union_smoke

Usage (production):
  python gears_union_finetune.py --epochs 20 --lr 1e-3 --batch-size 32 \\
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
from scipy.sparse import vstack, csr_matrix

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SAUNDERS_H5AD = (
    PROJECT_ROOT
    / "data/perturbation/datasets/saunders2025/processed/saunders_perturb_hep_confident.h5ad"
)
REPLOGLE_H5AD = (
    PROJECT_ROOT
    / "Analysis/Perturbation/results/finetuned_checkpoints/replogle_k562_essential/perturb_processed.h5ad"
)
NORMAN_H5AD = (
    PROJECT_ROOT
    / "Analysis/Perturbation/results/finetuned_checkpoints/norman/perturb_processed.h5ad"
)
ORTH_MAP = (
    PROJECT_ROOT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"
)


def prep_saunders(adata_in, m2h, smoke=False, smoke_cells=1000, seed=42):
    """Prep Saunders to Norman-compatible layout (human symbols, condition col)."""
    print(f"[saunders] input: {adata_in.shape}")
    var = adata_in.var.copy()
    var["mouse_symbol"] = var["gene_name"]
    var["human_symbol"] = var["mouse_symbol"].map(m2h)
    # Drop unmapped + dedup
    seen = set()
    keep = []
    for idx, hs in zip(var.index, var["human_symbol"]):
        if not isinstance(hs, str) or not hs or hs in seen:
            keep.append(False)
        else:
            seen.add(hs)
            keep.append(True)
    keep = pd.Series(keep, index=var.index)
    adata = adata_in[:, keep.values].copy()
    adata.var["gene_name"] = adata.var["human_symbol"]
    adata.var_names = adata.var["human_symbol"].astype(str).values
    adata.var_names_make_unique()

    # obs.condition
    target = adata.obs["target_gene"].astype(str)
    is_ctrl = target.isin(["control", "Control", "CTRL", "ctrl", "NTC"])
    human_target = adata.obs["target_gene"].map(m2h)
    cond = np.where(is_ctrl, "ctrl", human_target.astype(str) + "+ctrl")
    adata.obs["condition"] = cond
    valid = is_ctrl | human_target.notna()
    adata = adata[valid.values, :].copy()
    adata.obs["cell_source"] = "saunders"
    print(f"[saunders] after ortholog filter: {adata.shape}; conditions: {adata.obs['condition'].nunique()}")

    if smoke and len(adata) > smoke_cells:
        rng = np.random.RandomState(seed)
        n_per = max(5, smoke_cells // adata.obs["condition"].nunique())
        keep_ix = []
        for c, idx in adata.obs.groupby("condition").groups.items():
            take = min(len(idx), n_per)
            keep_ix.extend(rng.choice(list(idx), size=take, replace=False))
        adata = adata[keep_ix, :].copy()
        print(f"[saunders] smoke subsample: {len(adata):,}")
    return adata


def prep_k562_pertdata(adata_in, source_label, smoke=False, smoke_cells=1000, seed=42):
    """Pass a K562 PertData AnnData through with a cell_source label.

    Both Replogle K562 essential and Norman K562 combinatorial follow the same
    PertData layout: var.gene_name = human HGNC symbols, obs.condition =
    'GENE+ctrl' / 'GENE_A+GENE_B' / 'ctrl'.
    """
    print(f"[{source_label}] input: {adata_in.shape}")
    adata = adata_in.copy()
    adata.obs["cell_source"] = source_label
    if "gene_name" not in adata.var.columns:
        adata.var["gene_name"] = adata.var.index.astype(str)

    if smoke and len(adata) > smoke_cells:
        rng = np.random.RandomState(seed)
        n_per = max(5, smoke_cells // max(1, adata.obs["condition"].nunique()))
        keep_ix = []
        for c, idx in adata.obs.groupby("condition").groups.items():
            take = min(len(idx), n_per)
            keep_ix.extend(rng.choice(list(idx), size=take, replace=False))
        adata = adata[keep_ix, :].copy()
        print(f"[{source_label}] smoke subsample: {len(adata):,}")
    return adata


def union_anndatas(*adatas):
    """Concatenate N AnnDatas on gene-intersection namespace."""
    gene_sets = [set(a.var["gene_name"]) for a in adatas]
    intersect = sorted(set.intersection(*gene_sets))
    src_summary = "  ".join(f"|{a.obs['cell_source'].iloc[0]}|={len(g):,}" for a, g in zip(adatas, gene_sets))
    print(f"\n[union] gene intersection: {src_summary}  ∩={len(intersect):,}")

    def _subset(a, genes):
        m = a.var["gene_name"].isin(genes)
        a2 = a[:, m.values].copy()
        gene_to_pos = {g: i for i, g in enumerate(a2.var["gene_name"].astype(str))}
        order = [gene_to_pos[g] for g in genes if g in gene_to_pos]
        return a2[:, order].copy()

    subsets = [_subset(a, intersect) for a in adatas]
    ref_genes = list(subsets[0].var["gene_name"].astype(str))
    for s in subsets[1:]:
        assert list(s.var["gene_name"].astype(str)) == ref_genes, "Gene order mismatch after intersect"
    print(f"[union] post-intersect: " + "  ".join(f"{s.obs['cell_source'].iloc[0]}={s.shape}" for s in subsets))

    union = ad.concat(subsets, axis=0, join="outer", merge="first", label="dataset")
    union.var = subsets[0].var.copy()  # gene metadata
    print(f"[union] concatenated: {union.shape}; cell_source counts: {union.obs['cell_source'].value_counts().to_dict()}")
    print(f"[union] union conditions: {union.obs['condition'].nunique()}")
    return union


def maybe_normalize(adata):
    """If looks like raw counts (large integer values), normalize. Else assume already log-cpm."""
    if hasattr(adata.X, "toarray"):
        sample = np.asarray(adata.X[:200].toarray()).ravel()
    else:
        sample = np.asarray(adata.X[:200]).ravel()
    sample = sample[sample > 0]
    if len(sample) == 0:
        return adata
    if sample.max() > 30 or (sample == sample.astype(int)).mean() > 0.95:
        # Likely raw counts
        print(f"  X looks like raw counts (max={sample.max():.1f}, integer_frac={(sample==sample.astype(int)).mean():.2f}); applying normalize_total + log1p")
        if not adata.X.dtype.kind == "f":
            adata.X = adata.X.astype(np.float32)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print(f"  X looks already log-normalized (max={sample.max():.2f}); leaving as-is")
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="Subsample each source to ~1k cells, fewer epochs")
    ap.add_argument("--smoke-cells", type=int, default=1000,
                    help="Per-source cap when --smoke")
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--hidden-size", type=int, default=64)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--holdout-frac", type=float, default=0.10)
    ap.add_argument("--saunders-h5ad", default=str(SAUNDERS_H5AD))
    ap.add_argument("--replogle-h5ad", default=str(REPLOGLE_H5AD),
                    help="Path to Replogle K562 essential processed PertData adata")
    ap.add_argument("--norman-h5ad", default=str(NORMAN_H5AD),
                    help="Path to Norman K562 combinatorial processed PertData adata")
    ap.add_argument("--skip-norman", action="store_true",
                    help="Train on Replogle ∪ Saunders only (omit Norman; combinatorial training signal absent)")
    ap.add_argument("--out-dir", required=True,
                    help="Where to write the prepped PertData + checkpoint")
    args = ap.parse_args()

    t0 = time.time()
    print("=== GEARS Phase 2 union-vocab fine-tune (Replogle K562 ∪ Norman K562 ∪ Saunders mouse hep) ===")
    print(f"args: {vars(args)}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Mouse → human 1:1 ortholog map
    print(f"\nLoading ortholog map from {ORTH_MAP}")
    orth = pd.read_csv(ORTH_MAP, sep="\t").dropna(subset=["mouse_gene_symbol","human_gene_symbol"])
    orth_o2o = orth[orth["ortholog_type"]=="ortholog_one2one"]
    m2h = dict(zip(orth_o2o["mouse_gene_symbol"], orth_o2o["human_gene_symbol"]))
    print(f"  loaded {len(m2h):,} 1:1 mouse→human ortholog pairs")

    # ----- Saunders -----
    print(f"\nLoading Saunders: {args.saunders_h5ad}")
    adata_sa = sc.read_h5ad(args.saunders_h5ad)
    adata_sa = prep_saunders(adata_sa, m2h, smoke=args.smoke, smoke_cells=args.smoke_cells, seed=args.seed)

    # ----- Replogle K562 essential -----
    print(f"\nLoading Replogle K562 essential: {args.replogle_h5ad}")
    adata_rep = sc.read_h5ad(args.replogle_h5ad)
    adata_rep = prep_k562_pertdata(adata_rep, "replogle_k562", smoke=args.smoke, smoke_cells=args.smoke_cells, seed=args.seed)

    sources = [adata_sa, adata_rep]
    source_labels = ["saunders", "replogle_k562"]

    if not args.skip_norman:
        print(f"\nLoading Norman K562 combinatorial: {args.norman_h5ad}")
        adata_nor = sc.read_h5ad(args.norman_h5ad)
        adata_nor = prep_k562_pertdata(adata_nor, "norman_k562", smoke=args.smoke, smoke_cells=args.smoke_cells, seed=args.seed)
        sources.append(adata_nor)
        source_labels.append("norman_k562")

    # ----- Union -----
    print("\nNormalizing each (if raw)...")
    for s, lbl in zip(sources, source_labels):
        print(f"  {lbl}:"); maybe_normalize(s)
    union = union_anndatas(*sources)
    union.write_h5ad(out_dir / "prepped_union.h5ad", compression="gzip")
    print(f"\nWrote union AnnData: {out_dir/'prepped_union.h5ad'}")

    # ----- GEARS PertData -----
    print("\nImporting GEARS...")
    from gears import GEARS, PertData

    pert_data = PertData(data_path=str(out_dir))
    dataset_name = "saunders_hep_smoke" if args.smoke else "saunders_hep"
    try:
        pert_data.new_data_process(dataset_name=dataset_name, adata=union)
    except Exception as e:
        warnings.warn(f"new_data_process: {e}; trying .load()")
        pert_data.load(data_name=dataset_name)

    try:
        pert_data.prepare_split(split="simulation", seed=args.seed,
                                train_gene_set_size=(1 - args.holdout_frac))
    except Exception as e:
        warnings.warn(f"prepare_split with holdout failed ({e}); using default")
        pert_data.prepare_split(split="simulation", seed=args.seed)
    pert_data.get_dataloader(batch_size=args.batch_size, test_batch_size=max(64, args.batch_size * 2))

    # ----- Model + train -----
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nDevice: {device}; cuda avail: {torch.cuda.is_available()}")

    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=args.hidden_size)

    print(f"\nTraining {args.epochs} epochs (lr={args.lr}, bs={args.batch_size})...")
    model.train(epochs=args.epochs, lr=args.lr)

    # ----- Save -----
    ckpt = out_dir / "best.ckpt"
    ckpt.mkdir(parents=True, exist_ok=True)
    model.save_model(str(ckpt))
    print(f"\nCheckpoint saved to {ckpt}")

    pert_universe = None
    for attr in ("pert_names", "perturbable_genes", "gene_names"):
        if hasattr(pert_data, attr):
            v = getattr(pert_data, attr)
            if v is not None and len(v) > 0:
                pert_universe = list(v)
                break

    meta = {
        "args": vars(args),
        "n_cells_union": int(union.n_obs),
        "n_genes_union": int(union.n_vars),
        "n_cells_saunders": int((union.obs["cell_source"] == "saunders").sum()),
        "n_cells_replogle_k562": int((union.obs["cell_source"] == "replogle_k562").sum()),
        "n_cells_norman_k562": int((union.obs["cell_source"] == "norman_k562").sum()),
        "n_conditions_union": int(union.obs["condition"].nunique()),
        "pert_universe_size": len(pert_universe) if pert_universe else None,
        "pert_universe_sample": pert_universe[:30] if pert_universe else None,
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
