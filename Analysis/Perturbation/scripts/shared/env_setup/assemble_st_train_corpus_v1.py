#!/usr/bin/env python
"""Assemble a Perturb-seq training corpus for a new STATE Transition (ST) head.

Goal: expand perturbation vocabulary beyond the ~2,024-gene Replogle-essential
set used by the prebuilt HepG2 ST head, so MASLD hits (D3 tier1 pairs, D5
mouse orthologs) get higher coverage during ST-based perturbation inference.

Strategy v1 (what this script does):
  * Reformat the GEARS-processed Replogle K562 essential h5ad
    (1,093 single-gene perts × K562, ~163k cells)
  * Reformat the GEARS-processed Norman 2019 h5ad
    (282 perts × A549, ~91k cells — keeps singles "GENE+ctrl" AND
     combinatorial "GENE1+GENE2" tokens for downstream synergy inference)
  * **Harmonize HVG features across datasets** — cell-load's collate_fn
    stacks `obsm[embed_key]` across cells from different files, so each
    file's X_hvg MUST share the same dimensionality and (for biological
    consistency) the same gene set in the same order. We pick the top-N
    intersection of gene_name across both datasets, ranked by combined
    mean expression.
  * Emit cell-load-compatible h5ad files, one per cell_type:
      data/perturbation/state_train_v1/datasets/replogle/k562.h5ad
      data/perturbation/state_train_v1/datasets/norman/a549.h5ad
  * obs columns: gene, cell_type, gem_group, control (1 for non-targeting)
  * var columns: gene_name (Ensembl ID in index, gene_name col preserved)
  * Control rows have gene = "non-targeting" (cell-load default).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# ---------------------------------------------------------------------------
PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GEARS_DATA = PROJECT_ROOT / "Analysis/Perturbation/results/finetuned_checkpoints"
OUT_ROOT = PROJECT_ROOT / "data/perturbation/state_train_v1/datasets"
N_HVG = 2000  # common HVG count; both files end up with X_hvg shape (n_cells, N_HVG)
# ---------------------------------------------------------------------------


def reformat_condition(cond: str) -> str:
    """Map GEARS-style condition strings to cell-load `gene` tokens.

      'ctrl'             -> 'non-targeting'
      'GENE+ctrl'        -> 'GENE'         (single-gene perturbation)
      'ctrl+GENE'        -> 'GENE'
      'GENE1+GENE2'      -> 'GENE1+GENE2'  (combinatorial — token preserved)
    """
    if cond == "ctrl":
        return "non-targeting"
    parts = [p for p in str(cond).split("+") if p and p != "ctrl"]
    if not parts:
        return "non-targeting"
    if len(parts) == 1:
        return parts[0]
    # combinatorial — preserve as compound token (canonical sort)
    return "+".join(sorted(parts))


def derive_gem_group(adata: ad.AnnData) -> pd.Series:
    """Pull a batch label from `condition_name` ('CELL_COND_R1+R2') or fall
    back to 'gem_group' if absent."""
    if "gem_group" in adata.obs.columns:
        return adata.obs["gem_group"].astype(str)
    if "condition_name" in adata.obs.columns:
        names = adata.obs["condition_name"].astype(str)
        gem = names.str.rsplit("_", n=1).str[-1]
        gem = gem.where(gem.str.len() > 0, "1")
        return gem.fillna("1").astype(str)
    return pd.Series(["1"] * adata.n_obs, index=adata.obs.index, dtype=str)


def gene_name_series(a: ad.AnnData) -> pd.Series:
    """Best-effort gene_name resolution; falls back to Ensembl IDs."""
    if "gene_name" in a.var.columns:
        return a.var["gene_name"].astype(str)
    return pd.Series(a.var.index.astype(str), index=a.var.index)


def compute_shared_hvg(datasets: dict[str, ad.AnnData], n_hvg: int) -> list[str]:
    """Pick the top-`n_hvg` shared genes by combined-mean expression.

    Both inputs must share a meaningful subset of gene_names. We take the
    intersection (set on gene_name string), compute mean expression in each
    file, sum the ranks, and keep the top n.
    """
    # gene_name -> Ensembl index per file
    name_to_idx = {}
    for tag, a in datasets.items():
        gname = gene_name_series(a)
        name_to_idx[tag] = pd.Series(np.arange(a.n_vars), index=gname.values)

    shared = set(name_to_idx[next(iter(datasets))].index)
    for tag in list(datasets)[1:]:
        shared &= set(name_to_idx[tag].index)
    shared_sorted = sorted(shared)
    print(f"[hvg] shared gene_name count: {len(shared_sorted)}")

    # rank-sum across datasets
    rank_sum = pd.Series(0.0, index=shared_sorted)
    for tag, a in datasets.items():
        idx = name_to_idx[tag].loc[shared_sorted].values
        # mean over rows for those columns (sparse-friendly)
        if sp.issparse(a.X):
            means = np.asarray(a.X[:, idx].mean(axis=0)).ravel()
        else:
            means = a.X[:, idx].mean(axis=0)
        # rank descending — higher mean = lower (better) rank
        ranks = pd.Series(-means, index=shared_sorted).rank(method="average")
        rank_sum = rank_sum + ranks

    top = rank_sum.sort_values().head(n_hvg).index.tolist()
    print(f"[hvg] picked top {len(top)} HVGs by combined rank")
    return top


def reformat_one(
    in_path: Path, dataset_name: str, out_dir: Path, hvg_names: list[str]
) -> None:
    """Read GEARS h5ad, reformat, attach shared X_hvg, write per-cell_type h5ads."""
    print(f"[{dataset_name}] reading {in_path}")
    a = ad.read_h5ad(in_path)
    n0 = a.n_obs
    print(f"[{dataset_name}] loaded shape={a.shape}, unique conditions={a.obs['condition'].nunique()}")

    a.obs["gene"] = a.obs["condition"].astype(str).map(reformat_condition)

    if "cell_type" not in a.obs.columns:
        a.obs["cell_type"] = dataset_name
    a.obs["cell_type"] = a.obs["cell_type"].astype(str)

    a.obs["gem_group"] = derive_gem_group(a)
    a.obs["control"] = (a.obs["gene"] == "non-targeting").astype(int)

    if "gene_name" not in a.var.columns:
        a.var["gene_name"] = a.var.index.astype(str)

    if not sp.issparse(a.X):
        a.X = sp.csr_matrix(a.X.astype(np.float32))
    elif a.X.dtype != np.float32:
        a.X = a.X.astype(np.float32)

    # Shared X_hvg — same gene set & same order across all datasets.
    gname_series = gene_name_series(a)
    name_to_pos = pd.Series(np.arange(a.n_vars), index=gname_series.values)
    missing = [g for g in hvg_names if g not in name_to_pos.index]
    assert not missing, (
        f"[{dataset_name}] {len(missing)} shared HVG names missing in this file; "
        "compute_shared_hvg should have intersected first. First 5: "
        f"{missing[:5]}"
    )
    hvg_idx = name_to_pos.loc[hvg_names].values.astype(int)
    X_hvg = a.X[:, hvg_idx]
    if sp.issparse(X_hvg):
        X_hvg = X_hvg.toarray().astype(np.float32)
    a.obsm["X_hvg"] = X_hvg
    print(f"[{dataset_name}] X_hvg built: {X_hvg.shape}")

    obs_keep = ["gene", "cell_type", "gem_group", "control", "condition"]
    if "condition_name" in a.obs.columns:
        obs_keep.append("condition_name")
    a.obs = a.obs[[c for c in obs_keep if c in a.obs.columns]].copy()

    # Drop GEARS uns blobs (heavy + not needed for ST training) — preserve
    # hvg_names so cell-load can read uns/hvg_names when output_space="gene".
    for key in list(a.uns.keys()):
        del a.uns[key]
    a.uns["hvg_names"] = np.asarray(hvg_names, dtype=str)

    out_dir.mkdir(parents=True, exist_ok=True)
    for ct in a.obs["cell_type"].unique():
        sub = a[a.obs["cell_type"] == ct].to_memory()
        safe_ct = ct.replace("/", "_").replace(" ", "_").lower()
        out_path = out_dir / f"{safe_ct}.h5ad"
        sub.write_h5ad(out_path, compression="gzip")
        print(
            f"[{dataset_name}/{ct}] wrote {out_path} "
            f"shape={sub.shape} unique_perts={sub.obs['gene'].nunique()} "
            f"ctrl_cells={int(sub.obs['control'].sum())}"
        )

    vocab = a.obs["gene"].unique()
    print(f"[{dataset_name}] pert_vocab_size={len(vocab)}; "
          f"combinatorial_count={sum('+' in v for v in vocab)}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["replogle_k562_essential", "norman"],
    )
    ap.add_argument("--n-hvg", type=int, default=N_HVG)
    args = ap.parse_args()

    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    plan = {
        "replogle_k562_essential": (
            GEARS_DATA / "replogle_k562_essential/perturb_processed.h5ad",
            OUT_ROOT / "replogle",
        ),
        "norman": (
            GEARS_DATA / "norman/perturb_processed.h5ad",
            OUT_ROOT / "norman",
        ),
    }

    # 1. Read each dataset (memory-resident — they're 1-10 GB each, fine).
    loaded = {}
    for name in args.datasets:
        if name not in plan:
            print(f"[ERROR] unknown dataset {name}; skipping", file=sys.stderr)
            continue
        in_path, _ = plan[name]
        if not in_path.exists():
            print(f"[ERROR] {name}: missing {in_path}", file=sys.stderr)
            continue
        print(f"[load] {name} <- {in_path}")
        loaded[name] = ad.read_h5ad(in_path)
        print(f"[load] {name} shape={loaded[name].shape}")

    if not loaded:
        print("[ERROR] no datasets loaded", file=sys.stderr)
        return 1

    # 2. Pick the cross-dataset HVG list.
    hvg_names = compute_shared_hvg(loaded, args.n_hvg)

    # 3. Reformat + emit.
    for name, a in loaded.items():
        in_path, out_dir = plan[name]
        # we already loaded a in step 1; re-pass through reformat_one's logic
        # but avoid re-reading by stashing into a temp tag — simplest is to
        # write a tiny shim that pulls from in-memory a directly. Inline.
        print(f"[{name}] reformat...")
        # mimic reformat_one but with a passed-in adata
        a.obs["gene"] = a.obs["condition"].astype(str).map(reformat_condition)
        if "cell_type" not in a.obs.columns:
            a.obs["cell_type"] = name
        a.obs["cell_type"] = a.obs["cell_type"].astype(str)
        a.obs["gem_group"] = derive_gem_group(a)
        a.obs["control"] = (a.obs["gene"] == "non-targeting").astype(int)
        if "gene_name" not in a.var.columns:
            a.var["gene_name"] = a.var.index.astype(str)
        if not sp.issparse(a.X):
            a.X = sp.csr_matrix(a.X.astype(np.float32))
        elif a.X.dtype != np.float32:
            a.X = a.X.astype(np.float32)

        gname_series = gene_name_series(a)
        name_to_pos = pd.Series(np.arange(a.n_vars), index=gname_series.values)
        hvg_idx = name_to_pos.loc[hvg_names].values.astype(int)
        X_hvg = a.X[:, hvg_idx]
        if sp.issparse(X_hvg):
            X_hvg = X_hvg.toarray().astype(np.float32)
        a.obsm["X_hvg"] = X_hvg
        print(f"[{name}] X_hvg built: {X_hvg.shape}")

        obs_keep = ["gene", "cell_type", "gem_group", "control", "condition"]
        if "condition_name" in a.obs.columns:
            obs_keep.append("condition_name")
        a.obs = a.obs[[c for c in obs_keep if c in a.obs.columns]].copy()
        for key in list(a.uns.keys()):
            del a.uns[key]
        a.uns["hvg_names"] = np.asarray(hvg_names, dtype=str)

        out_dir.mkdir(parents=True, exist_ok=True)
        for ct in a.obs["cell_type"].unique():
            sub = a[a.obs["cell_type"] == ct].to_memory()
            safe_ct = ct.replace("/", "_").replace(" ", "_").lower()
            out_path = out_dir / f"{safe_ct}.h5ad"
            sub.write_h5ad(out_path, compression="gzip")
            print(
                f"[{name}/{ct}] wrote {out_path} "
                f"shape={sub.shape} unique_perts={sub.obs['gene'].nunique()} "
                f"ctrl_cells={int(sub.obs['control'].sum())}"
            )
        vocab = a.obs["gene"].unique()
        print(f"[{name}] pert_vocab_size={len(vocab)}; "
              f"combinatorial_count={sum('+' in v for v in vocab)}")

    print()
    print("[summary] outputs under", OUT_ROOT)
    for sub in sorted(OUT_ROOT.glob("**/*.h5ad")):
        size_mb = sub.stat().st_size / 1024 / 1024
        print(f"  {sub.relative_to(OUT_ROOT)}  {size_mb:.1f} MB")
    print(f"[summary] X_hvg dim across all files: {len(hvg_names)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
