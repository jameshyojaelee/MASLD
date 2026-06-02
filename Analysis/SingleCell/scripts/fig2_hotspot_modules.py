#!/usr/bin/env python3
"""
fig3_hotspot_modules.py — Hotspot gene module analysis on hepatocytes.

Finds co-varying gene programs in hepatocyte UMAP space.
Shows the 83% hepatocyte-intrinsic signal is organized into coherent modules.

Outputs to fig2_data/:
  hotspot_gene_modules.csv           — gene, module_id, z_score, local_corr_sum
  hotspot_module_disease_correlation.csv — module_id, n_genes, masld_vs_ctrl_effect, pval, padj, top_hub_genes
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import h5py
import hotspot
from scipy import stats
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC_DIR = os.path.join(BASE, "Analysis/SingleCell")
RESULTS = os.path.join(SC_DIR, "results_gpu_v2")
OUT_DIR = os.path.join(RESULTS, "fig2_data")
os.makedirs(OUT_DIR, exist_ok=True)

H5AD = os.path.join(RESULTS, "integrated_atlas.h5ad")
PB_DE = os.path.join(RESULTS, "pseudobulk_de", "Hepatocytes_de.csv")


def fix_ordered_attrs(h5ad_path):
    """Add missing 'ordered' attribute to categoricals in /obs and /var."""
    with h5py.File(h5ad_path, "a") as f:
        for group_name in ["obs", "var"]:
            if group_name not in f:
                continue
            grp = f[group_name]
            for key in grp.keys():
                item = grp[key]
                if isinstance(item, h5py.Group) and "categories" in item:
                    if "ordered" not in item.attrs:
                        item.attrs["ordered"] = False
                        log.info("Fixed missing 'ordered' attr: /%s/%s", group_name, key)


def build_ensembl_to_symbol_map(adata):
    """Build Ensembl ID → gene symbol mapping from adata.var if available."""
    for col in ["gene_name", "gene_symbols", "gene_short_name", "feature_name", "Symbol"]:
        if col in adata.var.columns:
            m = dict(zip(adata.var.index.tolist(), adata.var[col].tolist()))
            # Also build reverse: symbol → var_names index
            rev = {v: k for k, v in m.items()}
            log.info("Found symbol column '%s' in adata.var", col)
            return m, rev
    return {}, {}


def main():
    log.info("Fixing h5ad ordered attrs...")
    fix_ordered_attrs(H5AD)

    log.info("Loading h5ad...")
    adata = sc.read_h5ad(H5AD)
    log.info("Full atlas: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Subset to hepatocytes
    cell_type_col = None
    for col in ["cell_type", "celltype", "cell_type_coarse", "cell_type_leiden"]:
        if col in adata.obs.columns:
            cell_type_col = col
            break
    if cell_type_col is None:
        raise ValueError("No cell_type column found in adata.obs")

    hep_mask = adata.obs[cell_type_col].str.contains("Hepatocyte", case=False, na=False)
    adata_hep = adata[hep_mask].copy()
    log.info("Hepatocytes: %d cells", adata_hep.n_obs)

    if adata_hep.n_obs == 0:
        raise ValueError(f"No hepatocytes found in column '{cell_type_col}'")

    # Subsample to 100K cells stratified by condition for memory efficiency
    np.random.seed(42)
    target_n = 100_000
    if adata_hep.n_obs > target_n:
        condition_col = None
        for col in ["condition_harmonized", "condition", "disease_status", "group"]:
            if col in adata_hep.obs.columns:
                condition_col = col
                break
        if condition_col:
            # astype(str) avoids Categorical.fillna() rejection of new categories
            conditions = adata_hep.obs[condition_col].astype(str).fillna("Unknown")
            cond_vals = conditions.unique()
            keep_idx = []
            for c in cond_vals:
                c_idx = np.where(conditions == c)[0]
                n_c = int(target_n * len(c_idx) / adata_hep.n_obs)
                keep_idx.extend(np.random.choice(c_idx, min(n_c, len(c_idx)), replace=False).tolist())
            keep_idx = sorted(set(keep_idx))
        else:
            keep_idx = np.random.choice(adata_hep.n_obs, target_n, replace=False).tolist()
        adata_hep = adata_hep[keep_idx].copy()
    log.info("After subsampling: %d cells", adata_hep.n_obs)

    # Use raw counts for Hotspot (DANB model)
    # Subset raw to current obs/var to avoid shape mismatch (raw stores all original genes)
    if adata_hep.raw is not None:
        log.info("Using raw counts from adata.raw")
        adata_hep.X = adata_hep.raw[adata_hep.obs_names, adata_hep.var_names].X
    log.info("Count matrix shape: %d cells x %d genes", adata_hep.n_obs, adata_hep.n_vars)

    # Build Ensembl → symbol map
    _, sym_to_ens = build_ensembl_to_symbol_map(adata_hep)

    # Determine genes to test: prefer DEGs from pseudobulk, fallback to HVGs
    genes_to_test = None
    if os.path.exists(PB_DE):
        de = pd.read_csv(PB_DE)
        # gene column has Ensembl IDs in this dataset
        de_ens = de[de["padj"] < 0.1]["gene"].str.replace(r"\\..*", "", regex=True).tolist()
        log.info("DEGs (padj<0.1): %d Ensembl IDs", len(de_ens))

        # Check if var_names are Ensembl IDs or symbols
        var_names_list = adata_hep.var_names.tolist()
        n_ensg = sum(1 for g in var_names_list[:200] if str(g).startswith("ENSG"))

        if n_ensg > 50:
            # var_names are Ensembl IDs — match directly
            overlap = [g for g in de_ens if g in adata_hep.var_names]
            log.info("Direct Ensembl match: %d / %d DEGs in h5ad", len(overlap), len(de_ens))
        else:
            # var_names are gene symbols — map via gene_name column
            ensg_to_sym = {v: k for k, v in sym_to_ens.items()}  # ensembl → symbol
            de_syms = [ensg_to_sym.get(e, None) for e in de_ens]
            de_syms = [s for s in de_syms if s is not None]
            overlap = [g for g in de_syms if g in adata_hep.var_names]
            log.info("Symbol-mapped DEG match: %d / %d DEGs in h5ad", len(overlap), len(de_ens))

        if len(overlap) >= 100:
            genes_to_test = overlap
        else:
            log.warning("Insufficient DEG overlap (%d); trying padj<0.3", len(overlap))
            de_ens_broad = de[de["padj"] < 0.3]["gene"].str.replace(r"\\..*", "", regex=True).tolist()
            if n_ensg > 50:
                overlap_broad = [g for g in de_ens_broad if g in adata_hep.var_names]
            else:
                de_syms_broad = [ensg_to_sym.get(e, None) for e in de_ens_broad]
                de_syms_broad = [s for s in de_syms_broad if s is not None]
                overlap_broad = [g for g in de_syms_broad if g in adata_hep.var_names]
            if len(overlap_broad) >= 100:
                genes_to_test = overlap_broad
                log.info("Using padj<0.3 overlap: %d genes", len(overlap_broad))

    if genes_to_test is None or len(genes_to_test) < 100:
        # Fallback: use highly variable genes
        log.info("Falling back to HVGs for Hotspot input")
        if "highly_variable" in adata_hep.var.columns:
            genes_to_test = adata_hep.var_names[adata_hep.var["highly_variable"]].tolist()
            log.info("Using %d HVGs from adata.var['highly_variable']", len(genes_to_test))
        else:
            sc.pp.highly_variable_genes(adata_hep, n_top_genes=2000, flavor="seurat_v3", span=1.0)
            genes_to_test = adata_hep.var_names[adata_hep.var["highly_variable"]].tolist()
            log.info("Computed HVGs: %d genes", len(genes_to_test))
        genes_to_test = [g for g in genes_to_test if g in adata_hep.var_names]

    log.info("Genes for Hotspot: %d", len(genes_to_test))
    adata_hep_filt = adata_hep[:, genes_to_test].copy()

    # Ensure PCA embedding exists for KNN
    pca_key = None
    for key in ["X_pca_harmony", "X_pca", "X_scVI"]:
        if key in adata_hep_filt.obsm:
            pca_key = key
            break
    if pca_key is None:
        log.info("Computing PCA for KNN graph...")
        sc.pp.normalize_total(adata_hep_filt, target_sum=1e4)
        sc.pp.log1p(adata_hep_filt)
        sc.pp.pca(adata_hep_filt, n_comps=30)
        pca_key = "X_pca"
    log.info("Using embedding: %s", pca_key)

    # Init Hotspot
    log.info("Initializing Hotspot (DANB model)...")
    hs = hotspot.Hotspot(
        adata_hep_filt,
        layer_key=None,
        model="danb",
        latent_obsm_key=pca_key,
        umi_counts_obs_key=None
    )

    log.info("Creating KNN graph (n_neighbors=30)...")
    hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

    # Checkpoint paths — skip expensive steps if already computed
    ckpt_autocorr = os.path.join(OUT_DIR, "_hs_autocorr_ckpt.csv")
    ckpt_localcorr = os.path.join(OUT_DIR, "_hs_localcorr_ckpt.parquet")

    if os.path.exists(ckpt_autocorr) and os.path.exists(ckpt_localcorr):
        log.info("Loading autocorrelation checkpoint: %s", ckpt_autocorr)
        hs.results = pd.read_csv(ckpt_autocorr, index_col=0)
        log.info("Loading local correlation checkpoint: %s", ckpt_localcorr)
        hs.local_correlation_z = pd.read_parquet(ckpt_localcorr)
    else:
        log.info("Computing autocorrelations (jobs=8)...")
        hs.compute_autocorrelations(jobs=8)
        hs.results.to_csv(ckpt_autocorr)
        log.info("Autocorrelation checkpoint saved")

        # Filter to significant autocorrelated genes
        sig_genes = hs.results.index[hs.results["FDR"] < 0.05].tolist()
        log.info("Significant autocorrelated genes (FDR<0.05): %d", len(sig_genes))
        if len(sig_genes) < 20:
            sig_genes = hs.results.nsmallest(min(200, len(hs.results)), "FDR").index.tolist()
            log.warning("Too few sig genes; using top %d by FDR", len(sig_genes))

        log.info("Computing local correlations (jobs=8)...")
        hs.compute_local_correlations(genes=sig_genes, jobs=8)
        hs.local_correlation_z.to_parquet(ckpt_localcorr)
        log.info("Local correlation checkpoint saved")

    log.info("Creating gene modules...")
    modules = hs.create_modules(
        min_gene_threshold=10,
        core_only=False,
        fdr_threshold=0.05
    )

    # Module assignments
    module_df = modules.reset_index()
    module_df.columns = ["gene", "module_id"]
    n_total = len(module_df)
    module_df = module_df[module_df["module_id"] != -1]
    log.info("Genes in modules (excl. unassigned -1): %d / %d", len(module_df), n_total)

    # Add autocorrelation stats — hs.results columns: C, Z, Pval, FDR (gene in index)
    auto_stats = hs.results.rename(
        columns={"Z": "z_score", "Pval": "pval", "FDR": "fdr"}
    )[["z_score", "pval", "fdr"]].reset_index().rename(
        columns={hs.results.index.name or "index": "gene"}
    )
    module_df = module_df.merge(auto_stats[["gene", "z_score"]], on="gene", how="left")

    # Compute local_corr_sum (connectivity within module)
    if hasattr(hs, "local_correlation_z") and hs.local_correlation_z is not None:
        lc = hs.local_correlation_z
        corr_sums = {}
        for mod_id in module_df["module_id"].unique():
            same_mod = module_df.loc[module_df["module_id"] == mod_id, "gene"].tolist()
            mod_in_lc = [g for g in same_mod if g in lc.columns]
            for gene in mod_in_lc:
                if gene in lc.index:
                    other = [g for g in mod_in_lc if g != gene]
                    corr_sums[gene] = lc.loc[gene, other].sum() if other else 0.0
        module_df["local_corr_sum"] = module_df["gene"].map(corr_sums)

    out_genes = os.path.join(OUT_DIR, "hotspot_gene_modules.csv")
    module_df.to_csv(out_genes, index=False)
    log.info("Saved: %s (%d genes in modules)", out_genes, len(module_df))

    # Compute module scores per cell and disease correlation
    log.info("Computing module-disease correlation...")
    condition_col = None
    for col in ["condition_harmonized", "condition", "disease_status", "group"]:
        if col in adata_hep_filt.obs.columns:
            condition_col = col
            break

    module_results = []
    for mod_id in sorted(module_df["module_id"].unique()):
        mod_genes = module_df.loc[module_df["module_id"] == mod_id, "gene"].tolist()
        if len(mod_genes) < 5:
            continue

        in_var = [g for g in mod_genes if g in adata_hep_filt.var_names]
        if not in_var:
            continue

        X_mod = adata_hep_filt[:, in_var].X
        if hasattr(X_mod, "toarray"):
            X_mod = X_mod.toarray()
        scores = np.asarray(X_mod.mean(axis=1)).ravel()

        if condition_col:
            obs_cond = adata_hep_filt.obs[condition_col].astype(str)
            masld_mask = obs_cond.str.contains("MASLD|NASH|NAFLD|Disease", case=False, na=False)
            ctrl_mask  = obs_cond.str.contains("Healthy|Normal|Control", case=False, na=False)
            masld_scores = scores[masld_mask.values]
            ctrl_scores  = scores[ctrl_mask.values]
            if len(masld_scores) > 10 and len(ctrl_scores) > 10:
                _, pval = stats.mannwhitneyu(masld_scores, ctrl_scores, alternative="two-sided")
                effect = float(np.median(masld_scores) - np.median(ctrl_scores))
            else:
                pval, effect = 1.0, 0.0
        else:
            pval, effect = 1.0, 0.0

        top_hubs = (
            module_df.loc[module_df["module_id"] == mod_id]
            .nlargest(5, "z_score")["gene"]
            .tolist()
        )

        module_results.append({
            "module_id": mod_id,
            "n_genes": len(mod_genes),
            "mean_autocorr_z": module_df.loc[module_df["module_id"] == mod_id, "z_score"].mean(),
            "masld_vs_ctrl_effect": effect,
            "pval": pval,
            "top_hub_genes": ";".join(top_hubs)
        })

    mod_res = pd.DataFrame(module_results)
    if len(mod_res) > 0 and "pval" in mod_res.columns:
        _, padj, _, _ = multipletests(mod_res["pval"].fillna(1.0), method="fdr_bh")
        mod_res["padj"] = padj

    out_mods = os.path.join(OUT_DIR, "hotspot_module_disease_correlation.csv")
    mod_res.to_csv(out_mods, index=False)
    log.info("Saved: %s (%d modules)", out_mods, len(mod_res))
    log.info("=== Hotspot complete ===")


if __name__ == "__main__":
    main()
