#!/usr/bin/env python
"""
figS_scdrs_null_gen.py
======================
Generate hepatocyte-expression-matched null gene sets for the 4 GWAS-coloc
phenotype panels (PolyFun-anchored). Tests whether the "Hepatocytes harbor
GWAS-locus risk" finding in figS_scdrs_gwas Panel A survives a bias-corrected
null: bulk-liver eQTL ascertainment favors hepatocyte-expressed genes by
construction, so we need a null distribution drawn from genes with the SAME
hepatocyte-expression decile profile as each real panel.

Method:
  1. Load the preprocessed scDRS atlas.
  2. Compute per-gene mean expression in Hepatocytes only.
  3. Bin all atlas genes into 10 hep-expression deciles.
  4. For each of the 4 phenotype panels, tally the per-decile distribution
     of the real panel's genes.
  5. For each panel, generate 50 null gene sets matched to that per-decile
     distribution (sampled WITHOUT replacement; weights permuted from real
     panel to keep the weight distribution identical).
  6. Write a single .gs file with 200 null traits, suitable for the existing
     scDRS scoring pipeline (Analysis/SingleCell/scripts/400_mcp/
     run_scdrs_stage2_score_gwas.py works on any .gs file via env var).

Outputs:
  Analysis/SingleCell/results_gpu_v2/disease_signatures/
    GWAS_anchored_hep_matched_nulls.gs
    gene_hep_decile_metadata.csv  (gene -> hep_mean, hep_decile mapping)
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scdrs

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SIG_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
PP_ATLAS_PATH = SIG_DIR / "scdrs_preprocessed_atlas.h5ad"
GS_REAL_PATH  = SIG_DIR / "GWAS_anchored_scdrs.gs"
GS_NULL_PATH  = SIG_DIR / "GWAS_anchored_hep_matched_nulls.gs"
META_PATH     = SIG_DIR / "gene_hep_decile_metadata.csv"

TARGET_PHENOS = [
    "gwas_coloc_nafld_specific_polyfun",
    "gwas_coloc_liver_enzymes_polyfun",
    "gwas_coloc_pdff_polyfun",
    "gwas_coloc_cirrhosis_hcc_polyfun",
]
N_NULLS_PER_PHENO = 50
N_DECILES = 10
RANDOM_SEED = 42


def log(msg: str) -> None:
    print(f"[null_gen] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def main() -> None:
    rng = np.random.default_rng(RANDOM_SEED)

    log(f"loading preprocessed atlas: {PP_ATLAS_PATH}")
    t0 = time.time()
    adata = anndata.read_h5ad(PP_ATLAS_PATH)
    log(f"  atlas: {adata.n_obs:,} x {adata.n_vars:,};  loaded in {time.time()-t0:.1f}s")

    hep_mask = (adata.obs["cell_type"] == "Hepatocytes").values
    n_hep = int(hep_mask.sum())
    log(f"  hepatocyte cells: {n_hep:,}")
    if n_hep == 0:
        raise RuntimeError("No hepatocyte cells found in atlas")

    log("computing per-gene mean expression in hepatocytes")
    t0 = time.time()
    X_hep = adata.X[hep_mask, :]
    # adata.X may be sparse; .mean returns a matrix for sparse
    hep_mean = np.asarray(X_hep.mean(axis=0)).flatten()
    log(f"  hep_mean shape={hep_mean.shape};  took {time.time()-t0:.1f}s")

    gene_df = pd.DataFrame({
        "gene": adata.var_names.astype(str),
        "hep_mean": hep_mean,
    })
    log(f"  hep_mean stats: min={gene_df.hep_mean.min():.4g}  "
        f"median={gene_df.hep_mean.median():.4g}  "
        f"max={gene_df.hep_mean.max():.4g}  "
        f"zero_genes={(gene_df.hep_mean == 0).sum():,}")

    # Bin into N_DECILES. Drop duplicates at the edges (many zeros possible).
    gene_df["hep_decile"] = pd.qcut(
        gene_df["hep_mean"], N_DECILES,
        labels=False, duplicates="drop",
    )
    log(f"  decile counts:")
    for dec, n in gene_df["hep_decile"].value_counts().sort_index().items():
        log(f"    decile {dec}: {n:,} genes")

    gene_df.to_csv(META_PATH, index=False)
    log(f"  wrote {META_PATH}")

    log("loading real .gs file")
    gs_real = scdrs.util.load_gs(str(GS_REAL_PATH))
    log(f"  traits in .gs: {list(gs_real.keys())}")
    missing_targets = [t for t in TARGET_PHENOS if t not in gs_real]
    if missing_targets:
        raise RuntimeError(f"target phenotypes missing from .gs: {missing_targets}")

    log("generating matched nulls per phenotype")
    null_gs = {}
    summary_rows = []
    for pheno in TARGET_PHENOS:
        panel_genes, panel_weights = gs_real[pheno]
        panel_df = pd.DataFrame({"gene": panel_genes, "weight": panel_weights})
        panel_df = panel_df.merge(gene_df, on="gene", how="inner")
        panel_set = set(panel_df["gene"])
        decile_counts = panel_df["hep_decile"].value_counts().to_dict()
        log(f"  {pheno}: panel size = {len(panel_genes)};  "
            f"in-atlas = {len(panel_df)};  "
            f"decile profile = {sorted(decile_counts.items())}")

        n_nulls_built = 0
        for null_i in range(N_NULLS_PER_PHENO):
            sampled_genes = []
            for decile, count in decile_counts.items():
                pool = gene_df[
                    (gene_df["hep_decile"] == decile) &
                    (~gene_df["gene"].isin(panel_set))
                ]["gene"].values
                if len(pool) < count:
                    # very rare; use entire pool (no replacement still)
                    sampled = pool
                else:
                    sampled = rng.choice(pool, count, replace=False)
                sampled_genes.extend(sampled.tolist())

            # Permute the real panel's weights and truncate to sampled length
            shuffled_w = rng.permutation(panel_weights).tolist()
            sampled_weights = shuffled_w[:len(sampled_genes)]

            null_name = f"{pheno}__null{null_i:03d}"
            null_gs[null_name] = (sampled_genes, sampled_weights)
            n_nulls_built += 1
            summary_rows.append({
                "real_trait": pheno,
                "null_index": null_i,
                "n_genes": len(sampled_genes),
            })
        log(f"    built {n_nulls_built} null gene sets for {pheno}")

    log(f"writing {len(null_gs)} null traits to {GS_NULL_PATH}")
    with open(GS_NULL_PATH, "w") as f:
        f.write("TRAIT\tGENESET\n")
        for trait, (genes, weights) in null_gs.items():
            gs_str = ",".join([f"{g}:{w}" for g, w in zip(genes, weights)])
            f.write(f"{trait}\t{gs_str}\n")
    log(f"  wrote {GS_NULL_PATH}  ({sum(1 for _ in open(GS_NULL_PATH)) - 1} traits)")

    pd.DataFrame(summary_rows).to_csv(SIG_DIR / "GWAS_anchored_hep_matched_nulls_summary.csv", index=False)
    log(f"  wrote summary CSV")
    log("done.")


if __name__ == "__main__":
    main()
