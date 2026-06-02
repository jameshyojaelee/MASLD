#!/usr/bin/env python3
"""
71_cellex_specificity.py
---------------------------------------------------------------------------
Phase 4B Step 1: CELLEX expression specificity for CELLECT enrichment

Computes per-gene, per-cell-type expression specificity scores (ESmu)
using the CELLEX algorithm.

Input:
  - Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad

Output:
  - RNA-seq/results/causal_inference/cellect/cellex_esmu.csv
  - RNA-seq/results/causal_inference/cellect/cellex_esmu_top10pct.csv

Usage:
  python RNA-seq/71_cellex_specificity.py
---------------------------------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import cellex
import warnings
import time

warnings.filterwarnings("ignore")

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

ATLAS_PATH = os.path.join(
    BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
)

OUT_DIR = os.path.join(BASE, "RNA-seq/results/causal_inference/cellect")
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    print("=== CELLEX Expression Specificity ===")
    print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Atlas: {ATLAS_PATH}")
    print()

    # Load atlas
    print("Loading scRNA atlas...")
    adata = sc.read_h5ad(ATLAS_PATH)
    print(f"  Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

    # CELLEX needs raw counts in a DataFrame
    # Check if .raw exists or if X is raw counts
    if adata.raw is not None:
        print("  Using .raw layer for counts")
        count_matrix = adata.raw.X
        gene_names = list(adata.raw.var_names)
    else:
        print("  Using .X (assuming raw or normalized counts)")
        count_matrix = adata.X
        gene_names = list(adata.var_names)

    cell_types = adata.obs["cell_type"].values

    # Convert to DataFrame for CELLEX
    print("  Converting to DataFrame...")
    # CELLEX expects genes as rows, cells as columns
    # But for large datasets, use the ESObject directly
    import scipy.sparse as sp

    if sp.issparse(count_matrix):
        # For very large datasets, subsample to make CELLEX tractable
        n_cells = count_matrix.shape[0]
        if n_cells > 200000:
            print(f"  Subsampling from {n_cells:,} to 200,000 cells (stratified)...")
            np.random.seed(42)
            idx = []
            for ct in np.unique(cell_types):
                ct_idx = np.where(cell_types == ct)[0]
                n_sample = min(len(ct_idx), max(200, int(200000 * len(ct_idx) / n_cells)))
                idx.extend(np.random.choice(ct_idx, n_sample, replace=False))
            idx = sorted(idx)
            count_matrix = count_matrix[idx]
            cell_types = cell_types[idx]
            print(f"  After subsampling: {len(idx):,} cells")

    # CELLEX ESObject
    print("Computing expression specificity (ESmu)...")
    print("  This may take 30-60 minutes for large datasets...")
    t0 = time.time()

    eso = cellex.ESObject(
        data=count_matrix.T if sp.issparse(count_matrix) else count_matrix.T,
        annotation=pd.Series(cell_types, name="cell_type"),
        verbose=True,
    )
    # Set gene names
    eso.results["gene"] = gene_names

    # Compute ESmu (combined specificity metric)
    eso.compute(verbose=True)

    elapsed = time.time() - t0
    print(f"  CELLEX computation done in {elapsed/60:.1f} minutes")

    # Extract ESmu matrix
    esmu = eso.results["esmu"]
    print(f"  ESmu shape: {esmu.shape}")
    print(f"  Cell types: {list(esmu.columns)}")

    # Save full ESmu
    esmu_file = os.path.join(OUT_DIR, "cellex_esmu.csv")
    esmu.to_csv(esmu_file)
    print(f"  Saved: {esmu_file}")

    # Save top 10% specific genes per cell type (for CELLECT/S-LDSC)
    top_pct = 0.10
    top_genes = {}
    for ct in esmu.columns:
        threshold = esmu[ct].quantile(1 - top_pct)
        ct_genes = esmu[esmu[ct] >= threshold].index.tolist()
        top_genes[ct] = ct_genes
        print(f"  {ct}: {len(ct_genes)} genes in top {top_pct*100:.0f}%")

    # Save as long-format for CELLECT
    rows = []
    for ct, genes in top_genes.items():
        for g in genes:
            rows.append({"cell_type": ct, "gene": g, "esmu": esmu.loc[g, ct]})
    top_df = pd.DataFrame(rows)
    top_file = os.path.join(OUT_DIR, "cellex_esmu_top10pct.csv")
    top_df.to_csv(top_file, index=False)
    print(f"  Saved: {top_file}")

    # Summary stats
    print("\n=== ESmu Summary ===")
    for ct in esmu.columns:
        n_specific = (esmu[ct] > 0.1).sum()
        max_gene = esmu[ct].idxmax()
        max_val = esmu[ct].max()
        print(f"  {ct}: {n_specific} genes with ESmu>0.1, top={max_gene} ({max_val:.3f})")

    print(f"\n=== Done: {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
