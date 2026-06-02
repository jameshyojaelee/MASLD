#!/usr/bin/env python3
"""
Build a GEARS-ready AnnData (.h5ad) from:
  - 32 cellranger h5 GEX matrices (data/perturbation/datasets/saunders2025/raw/h5/)
  - Per-cell sgRNA assignment table (data/.../processed/per_cell_sgrna.tsv.gz)

Output AnnData layout:
  X: gene_expression (sparse, raw counts)
  obs:
    cell_barcode_24    : 24-nt 10x barcode  (16 CB + 8 probe BC)
    sublib             : e.g. "S5_B" — identifies GEM well + probe BC
    sample_bc          : BC001/002/003/004
    sgrna_id           : assigned sgRNA (only for tier=confident)
    target_gene        : the target gene (mouse symbol). 'control' for NCs.
    CR                 : CR1/CR2/CR3 scaffold variant
    tier               : confident/ambiguous_low_dominance/low_total_umis/low_top_umi/none
    n_sgrna_umis       : total sgRNA UMIs
    n_sgrna_reads      : total sgRNA reads
    dominant_frac      : top sgRNA fraction (UMI-weighted)
    human_ortholog     : mouse → human ortholog (where available)
  var:
    gene_id  : Ensembl mouse ID
    gene_name: mouse gene symbol (if available)

Usage:
    python3 08_build_anndata.py \\
        --per-cell data/perturbation/datasets/saunders2025/processed/per_cell_sgrna.tsv.gz \\
        --h5-glob 'data/perturbation/datasets/saunders2025/raw/h5/GSM8478327_*.h5' \\
        --out      data/perturbation/datasets/saunders2025/processed/saunders_perturb_hep.h5ad
"""
import argparse, csv, glob, gzip, os
from collections import defaultdict

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import vstack, csr_matrix


def load_per_cell(path):
    """Load per-cell sgRNA assignments into a dataframe keyed by 24-nt CB."""
    df = pd.read_csv(path, sep="\t", compression="gzip")
    df.set_index("cell_barcode_24", inplace=True)
    return df


def load_h5(h5_path, sublib_name):
    """Load one cellranger h5 into AnnData with sublib annotation."""
    with h5py.File(h5_path, "r") as h:
        m = h["matrix"]
        # Sparse CSC -> CSR (cells × genes)
        data = m["data"][:]
        indices = m["indices"][:]
        indptr = m["indptr"][:]
        shape = m["shape"][:]
        # cellranger writes (features, cells); transpose for AnnData (cells, features)
        from scipy.sparse import csc_matrix
        X = csc_matrix((data, indices, indptr), shape=tuple(shape))
        X = X.T.tocsr()
        barcodes = [b.decode().rstrip("-1").rstrip("-") for b in m["barcodes"][:]]
        features = m["features"]
        gene_ids = [g.decode() for g in features["id"][:]]
        gene_names = [g.decode() for g in features["name"][:]]
        feature_types = [g.decode() for g in features["feature_type"][:]]
    obs = pd.DataFrame({"sublib": sublib_name}, index=barcodes)
    var = pd.DataFrame({"gene_id": gene_ids, "gene_name": gene_names,
                        "feature_type": feature_types}, index=gene_ids)
    return ad.AnnData(X=X, obs=obs, var=var)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-cell", required=True)
    ap.add_argument("--h5-glob", required=True)
    ap.add_argument("--ortholog-map",
                    default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--confident-only", action="store_true",
                    help="Keep only cells with tier=='confident' assignment")
    args = ap.parse_args()

    print(f"Loading per-cell assignments from {args.per_cell}")
    per_cell = load_per_cell(args.per_cell)
    print(f"  {len(per_cell):,} assignments loaded")

    h5_files = sorted(glob.glob(args.h5_glob))
    print(f"Loading {len(h5_files)} h5 files")
    adatas = []
    for fp in h5_files:
        sublib_name = os.path.basename(fp).split("_sample")[0].replace("GSM8478327_", "")
        a = load_h5(fp, sublib_name)
        adatas.append(a)
        print(f"  {sublib_name}: {a.shape}")

    print("Concatenating...")
    adata = ad.concat(adatas, axis=0, join="outer", merge="first", label="sublib")
    print(f"  Combined: {adata.shape}")

    # Annotate with per-cell assignments
    obs = adata.obs.copy()
    obs["cell_barcode_24"] = obs.index
    # Reset index to ensure uniqueness (some CBs may repeat across sublibs).
    # Cast to str: after ad.concat, sublib may be Categorical -> str+str fails.
    obs["cell_id"] = obs["sublib"].astype(str) + "_" + obs["cell_barcode_24"].astype(str)
    obs.set_index("cell_id", inplace=True)

    # Merge with assignments
    obs = obs.merge(per_cell, left_on="cell_barcode_24", right_index=True, how="left",
                    suffixes=("", "_assign"))

    # Optionally filter to confident-only
    if args.confident_only:
        mask = obs["tier"].eq("confident")
        print(f"Filtering to confident-tier cells: {mask.sum():,}/{len(obs):,}")
        keep_idx = obs.index[mask]
        adata = adata[obs.index.isin(keep_idx)].copy()
        obs = obs.loc[keep_idx]

    # Set obs back to adata (handle index alignment)
    adata.obs = obs.loc[adata.obs.index] if adata.obs.index.is_unique else obs.iloc[:len(adata)]

    # Ortholog map (target_gene mouse → human)
    if os.path.exists(args.ortholog_map):
        print(f"Loading ortholog map from {args.ortholog_map}")
        orth = pd.read_csv(args.ortholog_map, sep="\t", compression="gzip" if args.ortholog_map.endswith(".gz") else None)
        # Try common column names
        if "mouse_gene_symbol" in orth.columns and "human_gene_symbol" in orth.columns:
            m2h = dict(zip(orth["mouse_gene_symbol"], orth["human_gene_symbol"]))
        elif "mouse_gene" in orth.columns and "human_gene" in orth.columns:
            m2h = dict(zip(orth["mouse_gene"], orth["human_gene"]))
        elif "MGI symbol" in orth.columns and "HGNC symbol" in orth.columns:
            m2h = dict(zip(orth["MGI symbol"], orth["HGNC symbol"]))
        else:
            print(f"  WARN: unknown ortholog map columns: {orth.columns.tolist()}")
            m2h = {}
        adata.obs["human_ortholog"] = adata.obs["target_gene"].map(m2h)
        n_orth = adata.obs["human_ortholog"].notna().sum()
        print(f"  {n_orth:,} of {len(adata):,} cells have human ortholog mapping")

    print(f"\nFinal AnnData: {adata.shape}")
    print(f"  obs columns: {adata.obs.columns.tolist()}")
    if "tier" in adata.obs.columns:
        print(f"  tier breakdown: {adata.obs['tier'].value_counts().to_dict()}")

    # Coerce object/mixed columns in obs/var to string (h5ad write_vlen_string_array
    # cannot serialize NaN-mixed object columns)
    for c in adata.obs.columns:
        if adata.obs[c].dtype == "object":
            adata.obs[c] = adata.obs[c].astype(str).fillna("").replace("nan", "")
    for c in adata.var.columns:
        if adata.var[c].dtype == "object":
            adata.var[c] = adata.var[c].astype(str).fillna("").replace("nan", "")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    print(f"\nWriting {args.out}")
    adata.write_h5ad(args.out, compression="gzip")
    print("DONE")


if __name__ == "__main__":
    main()
