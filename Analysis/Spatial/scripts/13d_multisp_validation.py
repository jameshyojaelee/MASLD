#!/usr/bin/env python3
"""
13d_multisp_validation.py — Validate MultiSP spatial multi-omics integration.

Assesses imputation quality (imputed gene-activity vs RNA in hepatocyte spots),
domain stability across parameter sweeps, biological validation against
zonation and fibrosis markers, and modality complementarity (genes where
RNA and chromatin spatial patterns disagree).

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
from scipy import sparse
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, C2L_PREFIX, load_config, load_deconvolved_adata,
    save_csv, print_header, print_step,
    strip_c2l_prefix,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_multisp_adata():
    """Load spatial AnnData with MultiSP domains from 13b."""
    path = RESULTS_DIR / "multisp" / "spatial_with_multisp_domains.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Run 13b first: {path}")
    return sc.read_h5ad(path)


def load_imputed_chromatin():
    """Load imputed chromatin AnnData from 13a."""
    path = RESULTS_DIR / "multisp" / "imputed_chromatin.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Run 13a first: {path}")
    return sc.read_h5ad(path)


def validate_imputation_quality(adata_spatial, adata_chromatin, config):
    """Correlate imputed gene-activity vs actual RNA in hepatocyte-dominant spots.

    Restricts to spots with hepatocyte fraction > hep_abundance_threshold.
    """
    hep_thresh = config["validation"]["hep_abundance_threshold"]

    # Find hepatocyte proportion column
    c2l_cols = [c for c in adata_spatial.obs.columns if c.startswith(C2L_PREFIX)]
    hep_col = None
    for c in c2l_cols:
        name = strip_c2l_prefix(c)
        if "hepatocyte" in name.lower() or "hep" in name.lower():
            hep_col = c
            break

    if hep_col is None:
        print_step("WARNING: Hepatocyte proportion column not found, using all spots")
        hep_mask = np.ones(adata_spatial.n_obs, dtype=bool)
    else:
        # Normalize proportions for the threshold comparison
        props = adata_spatial.obs[c2l_cols].values
        prop_sums = props.sum(axis=1, keepdims=True)
        prop_sums[prop_sums == 0] = 1.0
        props_norm = props / prop_sums
        hep_idx = c2l_cols.index(hep_col)
        hep_mask = props_norm[:, hep_idx] > hep_thresh
        print_step(f"Hepatocyte-dominant spots (>{hep_thresh}): {hep_mask.sum()}/{len(hep_mask)}")

    # Align spots
    common = adata_spatial.obs_names[hep_mask].intersection(adata_chromatin.obs_names)
    if len(common) < 50:
        print_step(f"WARNING: Only {len(common)} hepatocyte spots, relaxing threshold")
        hep_mask = np.ones(adata_spatial.n_obs, dtype=bool)
        common = adata_spatial.obs_names.intersection(adata_chromatin.obs_names)

    rna_sub = adata_spatial[common]
    chrom_sub = adata_chromatin[common]

    # Get matrices
    rna_mat = rna_sub.X
    if sparse.issparse(rna_mat):
        rna_mat = rna_mat.toarray()

    chrom_mat = chrom_sub.X
    if sparse.issparse(chrom_mat):
        chrom_mat = chrom_mat.toarray()

    rna_genes = rna_sub.var_names.tolist()
    chrom_genes = chrom_sub.var_names.tolist()

    # Find overlapping genes
    common_genes = list(set(rna_genes) & set(chrom_genes))
    if len(common_genes) == 0:
        print_step("WARNING: No overlapping genes between RNA and chromatin")
        return pd.DataFrame()

    print_step(f"Testing {len(common_genes)} genes in {len(common)} spots")

    # Compute per-gene correlation
    records = []
    for gene in common_genes:
        rna_idx = rna_genes.index(gene)
        chrom_idx = chrom_genes.index(gene)

        rna_vals = rna_mat[:, rna_idx]
        chrom_vals = chrom_mat[:, chrom_idx]

        if np.std(rna_vals) < 1e-10 or np.std(chrom_vals) < 1e-10:
            continue

        rho, pval = spearmanr(rna_vals, chrom_vals)
        records.append({
            "gene": gene, "spearman_rho": rho, "pval": pval,
            "rna_mean": np.mean(rna_vals), "chrom_mean": np.mean(chrom_vals),
            "rna_nonzero_frac": np.mean(rna_vals > 0),
        })

    result = pd.DataFrame(records)
    if len(result) > 0:
        # BH correction
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(result["pval"].values, method="fdr_bh")
        result["padj_bh"] = padj

        pos_rho = (result["spearman_rho"] > 0).mean()
        mean_rho = result["spearman_rho"].mean()
        median_rho = result["spearman_rho"].median()
        sig = (result["padj_bh"] < 0.05).sum()

        print(f"    Mean rho: {mean_rho:.3f}, median: {median_rho:.3f}")
        print(f"    Positive rho fraction: {pos_rho:.1%}")
        print(f"    Significant genes (padj<0.05): {sig}/{len(result)}")

    return result


def validate_domain_stability(adata):
    """Assess domain stability if multiple parameter sweep results exist."""
    sweep_path = RESULTS_DIR / "multisp" / "multisp_domain_sweep.csv"
    if not sweep_path.exists():
        print_step("No domain sweep results found, skipping stability check")
        return pd.DataFrame()

    sweep = pd.read_csv(sweep_path)
    if len(sweep) < 2:
        print_step("Only 1 parameter configuration tested, skipping stability")
        return pd.DataFrame()

    # If we have the multisp_domain as the final assignment,
    # check if there are alternative domain columns from different runs
    domain_cols = [c for c in adata.obs.columns if "domain" in c.lower()]
    if len(domain_cols) < 2:
        # Report sweep results only
        print_step(f"Domain sweep: {len(sweep)} configurations tested")
        print(f"    Silhouette range: {sweep['silhouette'].min():.3f} - {sweep['silhouette'].max():.3f}")
        return sweep

    # Compute pairwise ARI between domain assignments
    records = []
    for i, col_a in enumerate(domain_cols):
        for col_b in domain_cols[i+1:]:
            labels_a = adata.obs[col_a].astype(str).values
            labels_b = adata.obs[col_b].astype(str).values
            ari = adjusted_rand_score(labels_a, labels_b)
            records.append({
                "domain_set_a": col_a, "domain_set_b": col_b, "ARI": ari,
            })

    stability = pd.DataFrame(records)
    if len(stability) > 0:
        print_step(f"Domain stability (pairwise ARI):")
        for _, row in stability.iterrows():
            print(f"    {row['domain_set_a']} vs {row['domain_set_b']}: ARI={row['ARI']:.3f}")

    return stability


def validate_biology(adata, config):
    """Check domain markers vs known liver zonation and fibrosis markers."""
    zon_config = config["zonation"]
    pp_markers = zon_config["periportal_markers"]
    pc_markers = zon_config["pericentral_markers"]
    fibrosis_markers = ["COL1A1", "COL3A1", "ACTA2", "TGFB1", "PDGFRB"]

    spatial_genes = set(adata.var_names)

    records = []

    # Test each marker category
    for category, markers in [
        ("periportal", pp_markers),
        ("pericentral", pc_markers),
        ("fibrosis", fibrosis_markers),
    ]:
        for marker in markers:
            if marker not in spatial_genes:
                records.append({
                    "marker": marker, "category": category,
                    "domain_enriched": np.nan, "kruskal_pval": np.nan,
                    "max_mean": np.nan, "min_mean": np.nan,
                    "status": "not_in_data",
                })
                continue

            gene_idx = adata.var_names.tolist().index(marker)
            if sparse.issparse(adata.X):
                vals = np.asarray(adata.X[:, gene_idx].todense()).flatten()
            else:
                vals = adata.X[:, gene_idx]

            if np.std(vals) < 1e-10:
                records.append({
                    "marker": marker, "category": category,
                    "domain_enriched": np.nan, "kruskal_pval": np.nan,
                    "max_mean": np.nan, "min_mean": np.nan,
                    "status": "constant",
                })
                continue

            # Kruskal-Wallis across domains
            from scipy.stats import kruskal
            domains = sorted(adata.obs["multisp_domain"].unique())
            groups = [vals[adata.obs["multisp_domain"] == d] for d in domains]
            groups = [g for g in groups if len(g) >= 3]

            if len(groups) < 2:
                continue

            stat, pval = kruskal(*groups)

            # Find the domain with highest mean expression
            means = {}
            for d in domains:
                mask = adata.obs["multisp_domain"] == d
                means[d] = np.mean(vals[mask.values if hasattr(mask, 'values') else mask])
            best_domain = max(means, key=means.get)

            records.append({
                "marker": marker, "category": category,
                "domain_enriched": best_domain,
                "kruskal_pval": pval,
                "max_mean": max(means.values()),
                "min_mean": min(means.values()),
                "fold_change": max(means.values()) / (min(means.values()) + 1e-10),
                "status": "tested",
            })

    result = pd.DataFrame(records)
    tested = result[result["status"] == "tested"]
    if len(tested) > 0:
        sig = tested[tested["kruskal_pval"] < 0.05]
        print_step(f"Biological validation: {len(sig)}/{len(tested)} markers significantly "
                   f"differ across domains (p<0.05)")

        # Check if PP and PC markers segregate to different domains
        pp_domains = set(tested[tested["category"] == "periportal"]["domain_enriched"].dropna())
        pc_domains = set(tested[tested["category"] == "pericentral"]["domain_enriched"].dropna())
        if pp_domains and pc_domains:
            overlap = pp_domains & pc_domains
            print(f"    PP-enriched domains: {pp_domains}")
            print(f"    PC-enriched domains: {pc_domains}")
            print(f"    Separation: {'YES' if not overlap else 'partial' if len(overlap) < min(len(pp_domains), len(pc_domains)) else 'NO'}")

    return result


def assess_modality_complementarity(adata_spatial, adata_chromatin, config):
    """Identify genes where RNA and chromatin spatial patterns disagree.

    Computes Moran's I for each gene in both modalities and flags cases
    where the spatial patterns have opposite directions.
    """
    # Align spots
    common = adata_spatial.obs_names.intersection(adata_chromatin.obs_names)
    rna_sub = adata_spatial[common].copy()
    chrom_sub = adata_chromatin[common].copy()

    # Get overlapping genes
    common_genes = list(set(rna_sub.var_names) & set(chrom_sub.var_names))
    if len(common_genes) == 0:
        print_step("WARNING: No overlapping genes for complementarity analysis")
        return pd.DataFrame()

    # Subsample genes if too many (Moran's I is expensive)
    max_genes = 2000
    if len(common_genes) > max_genes:
        # Prioritize: HVGs + known markers
        hvg_set = set()
        if "highly_variable" in rna_sub.var.columns:
            hvg_set = set(rna_sub.var_names[rna_sub.var["highly_variable"]])
        zon_config = config["zonation"]
        marker_set = set(zon_config["periportal_markers"] + zon_config["pericentral_markers"])
        priority = list((hvg_set | marker_set) & set(common_genes))
        remaining = list(set(common_genes) - set(priority))
        np.random.seed(42)
        n_remaining = max_genes - len(priority)
        if n_remaining > 0 and remaining:
            sampled = list(np.random.choice(remaining, min(n_remaining, len(remaining)), replace=False))
        else:
            sampled = []
        common_genes = priority + sampled
        print_step(f"Subsampled to {len(common_genes)} genes ({len(priority)} priority)")

    # Build spatial neighbors for Moran's I
    sq.gr.spatial_neighbors(rna_sub, coord_type="generic", n_neighs=6)

    # Compute Moran's I for RNA
    rna_gene_idx = [rna_sub.var_names.tolist().index(g) for g in common_genes]
    rna_test = rna_sub[:, common_genes].copy()

    try:
        sq.gr.spatial_autocorr(rna_test, mode="moran", n_perms=100, n_jobs=8)
        rna_morans = rna_test.uns["moranI"].copy()
    except Exception as e:
        print_step(f"WARNING: Moran's I on RNA failed ({e})")
        return pd.DataFrame()

    # Compute Moran's I for chromatin
    chrom_gene_idx = [chrom_sub.var_names.tolist().index(g) for g in common_genes]
    chrom_test = chrom_sub[:, common_genes].copy()

    # Copy spatial connectivity from RNA (same spots)
    chrom_test.obsp["spatial_connectivities"] = rna_sub.obsp["spatial_connectivities"]
    chrom_test.obsp["spatial_distances"] = rna_sub.obsp["spatial_distances"]
    chrom_test.uns["spatial_neighbors"] = rna_sub.uns["spatial_neighbors"]

    try:
        sq.gr.spatial_autocorr(chrom_test, mode="moran", n_perms=100, n_jobs=8)
        chrom_morans = chrom_test.uns["moranI"].copy()
    except Exception as e:
        print_step(f"WARNING: Moran's I on chromatin failed ({e})")
        return pd.DataFrame()

    # Compare Moran's I between modalities
    # Align by gene name
    common_idx = rna_morans.index.intersection(chrom_morans.index)
    if len(common_idx) == 0:
        return pd.DataFrame()

    records = []
    for gene in common_idx:
        rna_I = rna_morans.loc[gene, "I"]
        chrom_I = chrom_morans.loc[gene, "I"]
        rna_pval = rna_morans.loc[gene, "pval_norm"] if "pval_norm" in rna_morans.columns else np.nan
        chrom_pval = chrom_morans.loc[gene, "pval_norm"] if "pval_norm" in chrom_morans.columns else np.nan

        # Determine agreement
        if np.isnan(rna_I) or np.isnan(chrom_I):
            agreement = "NA"
        elif rna_I > 0 and chrom_I > 0:
            agreement = "concordant_spatial"
        elif rna_I < 0 and chrom_I < 0:
            agreement = "concordant_random"
        else:
            agreement = "discordant"

        records.append({
            "gene": gene,
            "rna_morans_I": rna_I, "chrom_morans_I": chrom_I,
            "rna_pval": rna_pval, "chrom_pval": chrom_pval,
            "delta_morans_I": rna_I - chrom_I,
            "agreement": agreement,
        })

    result = pd.DataFrame(records)

    # Summary
    if len(result) > 0:
        agreement_counts = result["agreement"].value_counts()
        print_step("Modality spatial pattern agreement:")
        for cat, n in agreement_counts.items():
            print(f"    {cat}: {n} ({n/len(result):.1%})")

        # Most discordant genes
        discordant = result[result["agreement"] == "discordant"].copy()
        if len(discordant) > 0:
            discordant["abs_delta"] = discordant["delta_morans_I"].abs()
            top_disc = discordant.nlargest(10, "abs_delta")
            print(f"\n    Top 10 discordant genes:")
            for _, row in top_disc.iterrows():
                print(f"      {row['gene']}: RNA_I={row['rna_morans_I']:.3f}, "
                      f"chrom_I={row['chrom_morans_I']:.3f}")

    return result


def main():
    print_header("13d: MultiSP Validation")

    config = load_config()
    output_dir = RESULTS_DIR / "multisp"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load data ────────────────────────────────────────────────────
    print_step("Step 1: Loading data")
    adata_spatial = load_deconvolved_adata()
    adata_chromatin = load_imputed_chromatin()
    adata_multisp = load_multisp_adata()
    print(f"    Spatial: {adata_spatial.n_obs} spots x {adata_spatial.n_vars} genes")
    print(f"    Chromatin: {adata_chromatin.n_obs} spots x {adata_chromatin.n_vars} features")
    print(f"    MultiSP: {adata_multisp.n_obs} spots, {adata_multisp.obs['multisp_domain'].nunique()} domains")

    # ── Step 2: Imputation quality ───────────────────────────────────────────
    print_step("Step 2: Validating imputation quality (hepatocyte-dominant spots)")
    imputation_df = validate_imputation_quality(adata_spatial, adata_chromatin, config)
    if len(imputation_df) > 0:
        save_csv(imputation_df, "validation_imputation.csv", subdir="multisp")

    # ── Step 3: Domain stability ─────────────────────────────────────────────
    print_step("Step 3: Assessing domain stability")
    stability_df = validate_domain_stability(adata_multisp)
    if len(stability_df) > 0:
        save_csv(stability_df, "validation_domains.csv", subdir="multisp")

    # ── Step 4: Biological validation ────────────────────────────────────────
    print_step("Step 4: Biological validation (zonation + fibrosis markers)")
    biology_df = validate_biology(adata_multisp, config)
    if len(biology_df) > 0:
        save_csv(biology_df, "validation_biology.csv", subdir="multisp")

    # ── Step 5: Modality complementarity ─────────────────────────────────────
    print_step("Step 5: Assessing modality complementarity (Moran's I)")
    complementarity_df = assess_modality_complementarity(adata_spatial, adata_chromatin, config)
    if len(complementarity_df) > 0:
        save_csv(complementarity_df, "validation_complementarity.csv", subdir="multisp")

    # ── Final summary ────────────────────────────────────────────────────────
    print(f"\n  Validation summary:")
    if len(imputation_df) > 0:
        tested = imputation_df[imputation_df.get("padj_bh", imputation_df.get("pval", pd.Series(dtype=float))) < 0.05] if "padj_bh" in imputation_df.columns else pd.DataFrame()
        pos = (imputation_df["spearman_rho"] > 0).sum() if "spearman_rho" in imputation_df.columns else 0
        total = len(imputation_df[imputation_df.get("spearman_rho", pd.Series(dtype=float)).notna()]) if "spearman_rho" in imputation_df.columns else 0
        print(f"    Imputation: {pos}/{total} genes with positive RNA-chromatin correlation")
    if len(biology_df) > 0:
        tested_bio = biology_df[biology_df["status"] == "tested"]
        sig_bio = tested_bio[tested_bio["kruskal_pval"] < 0.05] if len(tested_bio) > 0 else pd.DataFrame()
        print(f"    Biology: {len(sig_bio)}/{len(tested_bio)} markers domain-specific (p<0.05)")
    if len(complementarity_df) > 0:
        disc = (complementarity_df["agreement"] == "discordant").sum()
        print(f"    Complementarity: {disc}/{len(complementarity_df)} genes discordant between modalities")

    print_header("13d: Complete")


if __name__ == "__main__":
    main()
