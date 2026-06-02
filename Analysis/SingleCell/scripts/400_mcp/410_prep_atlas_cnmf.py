#!/usr/bin/env python
"""
410_prep_atlas_cnmf.py — Build confounder-stripped raw-counts AnnData(s) for cNMF.

Outputs:
  results_gpu_v2/mcp/inputs/atlas_cnmf_global.h5ad                (~300K cells x ~20K genes, 5 major cell types)
  results_gpu_v2/mcp/inputs/atlas_cnmf_{celltype}.h5ad            (per-cell-type subsets with raw counts)
  results_gpu_v2/mcp/inputs/donor_metadata.tsv                    (donor-level phenotype table)
  results_gpu_v2/mcp/inputs/prep_summary.json                     (cell counts per cell type, per stage)

Key design choices (mirror bulk Script 95):
  - Use a.raw.X (raw integer counts) as the cNMF input
  - Apply identical confounder strip: chrY, chrM/MT-*, RPS/RPL/MRPS/MRPL, IG[HKL][VJCD],
    hemoglobin cluster, X-inactivation escapees (XIST/TSIX/KDM6A/DDX3X/EIF1AX/UTX/UBA1/
    RLIM/ZFX/RPS4X/EIF2S3/TXLNG), rRNA biotype. No protein-coding-only filter (scRNA should
    keep lncRNA).
  - Subset to 5 core cell types: Hepatocytes, Endothelial cells, Fibroblasts, Macrophages,
    Cholangiocytes. Keep each cell where donor has >=30 cells of that cell type.
  - Disease stage (0=Healthy, 1=Steatosis, 2=Steatohepatitis, 3=Cirrhosis) is propagated from
    hepatocyte_subtype_metadata.csv (sample-level). Samples absent from that file keep
    disease_stage as NaN and are retained (unsupervised factorization still sees them).
  - No batch correction. Dataset / cohort will be a regression covariate downstream.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

ATLAS_PATH = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
HEP_META_PATH = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv"
GENE_META_PATH = PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
OUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CORE_CELLTYPES = ["Hepatocytes", "Endothelial cells", "Fibroblasts", "Macrophages", "Cholangiocytes"]
MIN_CELLS_PER_DONOR_CT = 30

# --- confounder-strip rules (match RNA-seq/95_nmf_clean_ksweep.R verbatim) ------------------
RIBO_RE = re.compile(r"^(RPS|RPL|MRPS|MRPL)[0-9]")
IG_RE = re.compile(r"^IG[HKL][VJCD]")
HB_LIST = {"HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"}
X_ESCAPEE = {"XIST", "TSIX", "KDM6A", "DDX3X", "EIF1AX", "UTX", "UBA1", "RLIM", "ZFX",
             "RPS4X", "EIF2S3", "TXLNG"}
MT_PREFIX = "MT-"


def load_gene_metadata() -> pd.DataFrame:
    return pd.read_csv(GENE_META_PATH, sep="\t")


def build_strip_mask(var: pd.DataFrame, gene_meta: pd.DataFrame) -> pd.Series:
    """Boolean mask: True = DROP this gene."""
    var = var.copy()
    # Merge gene_biotype + chromosome from gencode using symbol
    gm = gene_meta[["gene_name", "chromosome", "gene_biotype"]].drop_duplicates("gene_name")
    var = var.merge(gm, how="left", left_index=True, right_on="gene_name").set_index("gene_name")
    sym = pd.Index(var.index.astype(str))

    drop = pd.Series(False, index=sym)
    # chrY (via gencode)
    drop |= var["chromosome"].eq("chrY").fillna(False).values
    # chrM / MT- prefix
    drop |= sym.str.startswith(MT_PREFIX)
    drop |= var["chromosome"].eq("chrM").fillna(False).values
    # rRNA biotype
    drop |= var["gene_biotype"].isin(["rRNA", "Mt_rRNA", "Mt_tRNA"]).fillna(False).values
    # Ribosomal
    drop |= sym.to_series().str.match(RIBO_RE).fillna(False).values
    # Immunoglobulin chains
    drop |= sym.to_series().str.match(IG_RE).fillna(False).values
    # Hemoglobin
    drop |= sym.isin(HB_LIST)
    # X-escapees
    drop |= sym.isin(X_ESCAPEE)
    return pd.Series(drop.values, index=sym, name="drop")


def load_donor_stage_map() -> pd.DataFrame:
    hep = pd.read_csv(HEP_META_PATH)
    # NOTE: use `len(s.dropna()) > 0` not `s.dropna().any()`: the latter is False for all-zero
    # stage entries (stage=0 = Healthy is falsy in Python), which incorrectly maps Healthy to NaN.
    def _mode_or_nan(s):
        d = s.dropna()
        if len(d) == 0:
            return np.nan
        return d.mode().iloc[0]

    per_sample = hep.groupby("sample").agg(
        disease_stage_coarse=("disease_stage_coarse", _mode_or_nan),
        disease_stage_numeric=("disease_stage_numeric", _mode_or_nan),
        dataset=("dataset", "first"),
        condition=("condition", "first"),
        condition_binary=("condition_binary", "first"),
    ).reset_index()
    return per_sample


def ensure_csr(x):
    if sp.issparse(x):
        return x.tocsr() if not sp.isspmatrix_csr(x) else x
    return sp.csr_matrix(x)


def main(args: argparse.Namespace) -> None:
    print(f"[410] Loading atlas (backed): {ATLAS_PATH}")
    adata = ad.read_h5ad(ATLAS_PATH, backed="r")
    print(f"[410] Full atlas: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # 1. Cell-type subset
    ct_mask = adata.obs["cell_type"].isin(CORE_CELLTYPES).values
    print(f"[410] Core cell types keep: {ct_mask.sum():,} cells")

    # 2. Donor x cell_type >= 30 cells filter
    obs = adata.obs.loc[ct_mask, ["sample", "cell_type", "dataset", "condition"]].copy()
    donor_ct_counts = obs.groupby(["sample", "cell_type"], observed=True).size().rename("n").reset_index()
    keep_pairs = set(map(tuple, donor_ct_counts.loc[donor_ct_counts["n"] >= MIN_CELLS_PER_DONOR_CT, ["sample", "cell_type"]].values))
    obs["keep"] = list(map(lambda s, c: (s, c) in keep_pairs, obs["sample"].astype(str), obs["cell_type"].astype(str)))
    print(f"[410] After donor>=30 per cell-type filter: {obs['keep'].sum():,} cells in {len(keep_pairs):,} (donor, CT) pairs")

    # Combined boolean over all adata cells
    full_keep = np.zeros(adata.shape[0], dtype=bool)
    full_keep[np.where(ct_mask)[0][obs["keep"].values]] = True
    print(f"[410] Final cell-keep count: {full_keep.sum():,}")

    # 3. Gene strip (use .var index which is gene_name)
    gene_meta = load_gene_metadata()
    drop_mask = build_strip_mask(adata.var, gene_meta)
    drop_mask = drop_mask.reindex(adata.var.index).fillna(False)
    keep_genes = (~drop_mask.values)
    print(f"[410] Strip {drop_mask.sum():,} confounder genes; retain {keep_genes.sum():,}")

    # 4. Materialize the subset using raw counts
    print("[410] Materializing subset (raw counts)...")
    sub = adata[full_keep, :].to_memory()
    adata.file.close()
    # Replace X with raw counts
    X_raw = sub.raw.X
    X_raw = ensure_csr(X_raw)
    # Cast to float32 (cNMF expects float; float32 halves memory)
    X_raw = X_raw.astype(np.float32)
    # Apply gene filter
    X_raw = X_raw[:, keep_genes]
    var_kept = sub.var.loc[keep_genes].copy()
    print(f"[410] Subset: {X_raw.shape}")

    # 5. Attach phenotype metadata (stage + donor fields)
    stage_map = load_donor_stage_map()
    obs_full = sub.obs.copy()
    obs_full = obs_full.merge(
        stage_map[["sample", "disease_stage_coarse", "disease_stage_numeric", "condition_binary"]],
        how="left", on="sample", suffixes=("", "_stage"),
    )
    obs_full.index = sub.obs.index

    # 6. Build sparse AnnData with raw counts as X and raw for safe access
    prep = ad.AnnData(
        X=X_raw,
        obs=obs_full,
        var=var_kept,
    )
    # Store raw counts also as layer for cNMF access
    prep.layers["counts"] = prep.X.copy()
    prep.uns["prep_params"] = {
        "min_cells_per_donor_ct": MIN_CELLS_PER_DONOR_CT,
        "core_celltypes": CORE_CELLTYPES,
        "n_dropped_genes": int(drop_mask.sum()),
        "atlas_source": str(ATLAS_PATH),
    }

    # 7. Save global
    global_out = OUT_DIR / "atlas_cnmf_global.h5ad"
    print(f"[410] Writing global: {global_out}")
    prep.write_h5ad(global_out, compression="gzip")

    # 8. Save per-cell-type
    ct_cell_counts = {}
    for ct in CORE_CELLTYPES:
        ct_safe = ct.lower().replace(" ", "_").replace("+", "_").replace("/", "_")
        sub_ct = prep[prep.obs["cell_type"] == ct].copy()
        out = OUT_DIR / f"atlas_cnmf_{ct_safe}.h5ad"
        print(f"[410] Writing {ct}: {sub_ct.shape} -> {out}")
        sub_ct.write_h5ad(out, compression="gzip")
        ct_cell_counts[ct] = int(sub_ct.shape[0])

    # 9. Donor metadata summary
    donor_meta = prep.obs.groupby("sample", observed=True).agg(
        dataset=("dataset", "first"),
        condition=("condition", "first"),
        condition_binary=("condition_binary", "first"),
        disease_stage_coarse=("disease_stage_coarse", "first"),
        disease_stage_numeric=("disease_stage_numeric", "first"),
        n_cells=("cell_type", "size"),
    ).reset_index()
    # Cell-type fractions per donor
    ct_wide = prep.obs.groupby(["sample", "cell_type"], observed=True).size().unstack(fill_value=0)
    ct_wide = ct_wide.divide(ct_wide.sum(axis=1), axis=0).add_prefix("frac_")
    donor_meta = donor_meta.merge(ct_wide.reset_index(), on="sample", how="left")
    donor_meta.to_csv(OUT_DIR / "donor_metadata.tsv", sep="\t", index=False)
    print(f"[410] donor_metadata: {donor_meta.shape[0]} donors")

    # 10. Summary JSON
    summary = {
        "n_cells_total_kept": int(prep.shape[0]),
        "n_genes_total_kept": int(prep.shape[1]),
        "n_genes_dropped": int(drop_mask.sum()),
        "cells_per_celltype": ct_cell_counts,
        "stage_distribution": prep.obs["disease_stage_numeric"].value_counts(dropna=False).to_dict(),
        "datasets": prep.obs["dataset"].value_counts().to_dict(),
        "n_donors": int(prep.obs["sample"].nunique()),
        "n_donors_with_stage": int(prep.obs.loc[prep.obs["disease_stage_numeric"].notna(), "sample"].nunique()),
    }
    with open(OUT_DIR / "prep_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[410] summary: {summary}")
    print("[410] DONE.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    _args = ap.parse_args()
    main(_args)
