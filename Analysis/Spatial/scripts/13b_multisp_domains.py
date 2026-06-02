#!/usr/bin/env python3
"""
13b_multisp_domains.py — Bimodal spatial domain identification (RNA + chromatin).

Constructs a bimodal embedding from RNA expression (3K HVGs) and imputed
gene-activity (from 13a), then identifies spatial domains using either
MultiSP or a manual spatially-smoothed Leiden fallback. Compares against
the existing 05d RNA-only domains.

SLURM: --partition=gpu --gres=gpu:1 --cpus=8 --mem=64G --time=8:00:00
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
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import normalize

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    save_checkpoint, save_csv, print_header, print_step,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_imputed_chromatin():
    """Load imputed chromatin AnnData from 13a."""
    path = RESULTS_DIR / "multisp" / "imputed_chromatin.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Run 13a first: {path}")
    print_step(f"Loading imputed chromatin: {path}")
    return sc.read_h5ad(path)


def build_bimodal_adata(adata_spatial, adata_chromatin, n_hvgs=3000):
    """Construct bimodal AnnData: RNA as .X, chromatin as .obsm['chromatin'].

    Ensures spots are aligned between the two modalities.
    """
    # Align spots
    common_spots = adata_spatial.obs_names.intersection(adata_chromatin.obs_names)
    if len(common_spots) < adata_spatial.n_obs:
        print(f"    WARNING: {adata_spatial.n_obs - len(common_spots)} spots lost in alignment")
    adata_rna = adata_spatial[common_spots].copy()
    adata_chrom = adata_chromatin[common_spots].copy()

    # Prepare RNA: ensure log-normalized + HVGs
    if "counts" in adata_rna.layers:
        adata_rna.X = adata_rna.layers["counts"].copy()
        sc.pp.highly_variable_genes(adata_rna, n_top_genes=n_hvgs, flavor="seurat_v3")
        sc.pp.normalize_total(adata_rna, target_sum=1e4)
        sc.pp.log1p(adata_rna)
    elif "highly_variable" not in adata_rna.var.columns:
        sc.pp.highly_variable_genes(adata_rna, n_top_genes=n_hvgs)

    # Store chromatin in obsm
    chrom_mat = adata_chrom.X
    if sparse.issparse(chrom_mat):
        chrom_mat = chrom_mat.toarray()
    adata_rna.obsm["chromatin"] = chrom_mat.astype(np.float32)

    print_step(f"Bimodal AnnData: {adata_rna.n_obs} spots, "
               f"RNA={adata_rna.n_vars} genes, chromatin={chrom_mat.shape[1]} features")
    return adata_rna


def try_multisp(adata_bimodal, config):
    """Attempt to run MultiSP for bimodal domain detection.

    Returns (domains, success_flag).
    """
    multisp_config = config["multisp"]
    try:
        import multisp
        print_step("MultiSP package found, running bimodal integration")

        # MultiSP API: spatial+feature fusion with modality-specific VAE
        # The exact API depends on the installed version
        model = multisp.MultiSP(
            adata=adata_bimodal,
            modality_keys={"rna": None, "chromatin": "chromatin"},
            n_latent=multisp_config["n_latent"],
        )
        model.train(max_epochs=multisp_config["epochs"])
        adata_bimodal.obsm["X_multisp"] = model.get_latent_representation()
        print_step(f"MultiSP latent space: {adata_bimodal.obsm['X_multisp'].shape}")
        return True

    except ImportError:
        print_step("MultiSP not available, using manual bimodal fallback")
        return False
    except Exception as e:
        print_step(f"MultiSP failed ({e}), using manual bimodal fallback")
        return False


def manual_bimodal_embedding(adata_bimodal, config):
    """Manual bimodal domain detection: PCA on RNA + PCA on chromatin,
    concatenated with modality weighting, then spatially-smoothed Leiden.
    """
    n_pcs = 30

    # PCA on RNA (HVGs)
    if "highly_variable" in adata_bimodal.var.columns:
        sc.pp.pca(adata_bimodal, n_comps=n_pcs, use_highly_variable=True)
    else:
        sc.pp.pca(adata_bimodal, n_comps=n_pcs)
    rna_pcs = adata_bimodal.obsm["X_pca"].copy()

    # PCA on chromatin
    from sklearn.decomposition import PCA
    chrom_mat = adata_bimodal.obsm["chromatin"]
    # Handle NaN/Inf
    chrom_mat = np.nan_to_num(chrom_mat, nan=0.0, posinf=0.0, neginf=0.0)
    n_chrom_pcs = min(n_pcs, chrom_mat.shape[1] - 1, chrom_mat.shape[0] - 1)
    if n_chrom_pcs < 2:
        print_step("WARNING: Too few chromatin features for PCA, using RNA-only")
        return

    pca_chrom = PCA(n_components=n_chrom_pcs, random_state=42)
    chrom_pcs = pca_chrom.fit_transform(chrom_mat)
    print_step(f"Chromatin PCA: {chrom_pcs.shape[1]} PCs "
               f"({pca_chrom.explained_variance_ratio_.sum():.1%} variance)")

    # Concatenate with equal weighting (scale each modality to unit variance)
    rna_scaled = rna_pcs / (rna_pcs.std(axis=0, keepdims=True) + 1e-8)
    chrom_scaled = chrom_pcs / (chrom_pcs.std(axis=0, keepdims=True) + 1e-8)
    combined = np.hstack([rna_scaled, chrom_scaled])

    adata_bimodal.obsm["X_multisp"] = combined
    print_step(f"Combined embedding: {combined.shape[1]} dimensions "
               f"(RNA={rna_pcs.shape[1]}, chromatin={chrom_pcs.shape[1]})")


def run_spatially_smoothed_leiden(adata, config, embedding_key="X_multisp"):
    """Run spatially-smoothed Leiden clustering (same approach as 05d)."""
    dom_config = config["domains"]
    alpha = dom_config["spatial_weight"]
    n_expr = dom_config["n_expression_neighbors"]
    n_spat = dom_config["n_spatial_neighbors"]

    # Expression-based neighbors on bimodal embedding
    sc.pp.neighbors(adata, n_neighbors=n_expr, use_rep=embedding_key,
                    key_added="bimodal_neighbors")

    # Spatial neighbors
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=n_spat,
                            key_added="spatial_neighbors")

    # Combine connectivity matrices
    expr_conn = adata.obsp["bimodal_neighbors_connectivities"]
    spat_conn = adata.obsp["spatial_neighbors_connectivities"]

    expr_norm = normalize(expr_conn, norm="l1", axis=1)
    spat_norm = normalize(spat_conn, norm="l1", axis=1)
    combined = (1 - alpha) * expr_norm + alpha * spat_norm

    adata.obsp["connectivities"] = combined
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "params": {"method": "bimodal_spatial", "alpha": alpha},
    }

    return adata


def sweep_domains(adata, config, key="multisp_domain"):
    """Sweep n_domains range, select optimal by silhouette score."""
    multisp_config = config["multisp"]
    n_range = multisp_config["n_domains_range"]
    dom_config = config["domains"]
    resolutions = dom_config["resolution_range"]

    # Extend resolution range for finer sweeping
    all_resolutions = sorted(set(resolutions + [r + 0.1 for r in resolutions]
                                 + [r - 0.1 for r in resolutions if r - 0.1 > 0]))

    best_sil = -1
    best_res = 0.5
    best_n = 0
    embedding = adata.obsm.get("X_multisp", adata.obsm.get("X_pca"))

    results = []
    for resolution in all_resolutions:
        sc.tl.leiden(adata, resolution=resolution, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        n_clusters = adata.obs[key].nunique()

        if n_clusters < 2 or n_clusters not in n_range:
            continue

        labels = adata.obs[key].astype(int).values
        sil = silhouette_score(embedding, labels, sample_size=min(5000, len(labels)))
        results.append({
            "resolution": resolution, "n_domains": n_clusters,
            "silhouette": sil,
        })

        if sil > best_sil:
            best_sil = sil
            best_res = resolution
            best_n = n_clusters

    # Apply best resolution
    if best_n > 0:
        sc.tl.leiden(adata, resolution=best_res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        print_step(f"Optimal: {best_n} domains (resolution={best_res:.2f}, silhouette={best_sil:.3f})")
    else:
        # Fallback: use target from config
        target = dom_config["target_n_domains"]
        best_diff = 999
        for resolution in all_resolutions:
            sc.tl.leiden(adata, resolution=resolution, key_added=key,
                         flavor="igraph", n_iterations=2, directed=False)
            n_clusters = adata.obs[key].nunique()
            diff = abs(n_clusters - target)
            if diff < best_diff:
                best_diff = diff
                best_res = resolution
                best_n = n_clusters
            if diff <= 1:
                break
        sc.tl.leiden(adata, resolution=best_res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        print_step(f"Fallback: {best_n} domains (resolution={best_res:.2f}, target={target})")

    return adata, pd.DataFrame(results)


def compare_with_05d(adata, key="multisp_domain"):
    """Compare MultiSP domains with existing 05d RNA-only domains."""
    # Load 05d domain assignments
    dom_path = RESULTS_DIR / "domains" / "spatial_domains.csv"
    if not dom_path.exists():
        print_step("WARNING: 05d domains not found, skipping comparison")
        return pd.DataFrame()

    dom_05d = pd.read_csv(dom_path, index_col=0)
    if "spatial_domain" not in dom_05d.columns:
        print_step("WARNING: spatial_domain column not found in 05d results")
        return pd.DataFrame()

    # Align spots
    common = adata.obs_names.intersection(dom_05d.index)
    if len(common) == 0:
        print_step("WARNING: No overlapping spots between 05d and MultiSP")
        return pd.DataFrame()

    labels_05d = dom_05d.loc[common, "spatial_domain"].astype(str).values
    labels_multisp = adata[common].obs[key].astype(str).values

    ari = adjusted_rand_score(labels_05d, labels_multisp)
    nmi = normalized_mutual_info_score(labels_05d, labels_multisp)

    print(f"    ARI (05d vs MultiSP): {ari:.3f}")
    print(f"    NMI (05d vs MultiSP): {nmi:.3f}")

    comparison = pd.DataFrame([{
        "comparison": "05d_rna_only_vs_multisp_bimodal",
        "ARI": ari, "NMI": nmi,
        "n_spots": len(common),
        "n_domains_05d": len(set(labels_05d)),
        "n_domains_multisp": len(set(labels_multisp)),
    }])
    return comparison


def run_rna_only_ablation(adata, config, key="rna_only_domain"):
    """Ablation: run RNA-only domains with the same pipeline for fair comparison."""
    print_step("Running RNA-only ablation (same Leiden, no chromatin)")

    # PCA on RNA only
    if "highly_variable" in adata.var.columns:
        sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    else:
        sc.pp.pca(adata, n_comps=30)

    # Clear any prior neighbors entry to avoid KeyError from non-standard params
    if "neighbors" in adata.uns:
        del adata.uns["neighbors"]

    # Use RNA PCA for neighbors
    sc.pp.neighbors(adata, n_neighbors=config["domains"]["n_expression_neighbors"],
                    n_pcs=30, key_added="rna_only_neighbors")

    # Spatial neighbors (reuse if already computed)
    if "spatial_neighbors_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic",
                                n_neighs=config["domains"]["n_spatial_neighbors"],
                                key_added="spatial_neighbors")

    alpha = config["domains"]["spatial_weight"]
    expr_conn = adata.obsp["rna_only_neighbors_connectivities"]
    spat_conn = adata.obsp["spatial_neighbors_connectivities"]

    expr_norm = normalize(expr_conn, norm="l1", axis=1)
    spat_norm = normalize(spat_conn, norm="l1", axis=1)
    combined = (1 - alpha) * expr_norm + alpha * spat_norm

    adata.obsp["connectivities"] = combined
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "params": {"method": "rna_only_spatial", "alpha": alpha},
    }

    # Cluster with same target
    target = config["domains"]["target_n_domains"]
    best_res = 0.5
    best_diff = 999
    for resolution in config["domains"]["resolution_range"]:
        sc.tl.leiden(adata, resolution=resolution, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        n_clusters = adata.obs[key].nunique()
        diff = abs(n_clusters - target)
        if diff < best_diff:
            best_diff = diff
            best_res = resolution
        if diff <= 1:
            break

    sc.tl.leiden(adata, resolution=best_res, key_added=key,
                 flavor="igraph", n_iterations=2, directed=False)
    n_domains = adata.obs[key].nunique()
    print(f"    RNA-only ablation: {n_domains} domains (resolution={best_res:.2f})")
    return adata


def main():
    print_header("13b: MultiSP Bimodal Domain Identification")

    config = load_config()
    output_dir = RESULTS_DIR / "multisp"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load data ────────────────────────────────────────────────────
    print_step("Step 1: Loading spatial + imputed chromatin data")
    adata_spatial = load_deconvolved_adata()
    adata_chromatin = load_imputed_chromatin()

    # ── Step 2: Build bimodal AnnData ────────────────────────────────────────
    print_step("Step 2: Building bimodal AnnData (RNA + chromatin)")
    adata = build_bimodal_adata(adata_spatial, adata_chromatin,
                                n_hvgs=config["normalization"]["n_top_genes"])

    # ── Step 3: Compute bimodal embedding ────────────────────────────────────
    print_step("Step 3: Computing bimodal embedding")
    multisp_success = try_multisp(adata, config)
    if not multisp_success:
        manual_bimodal_embedding(adata, config)

    # ── Step 4: Spatially-smoothed Leiden ─────────────────────────────────────
    print_step("Step 4: Running spatially-smoothed Leiden on bimodal embedding")
    adata = run_spatially_smoothed_leiden(adata, config)
    adata, sweep_df = sweep_domains(adata, config)

    if len(sweep_df) > 0:
        save_csv(sweep_df, "multisp_domain_sweep.csv", subdir="multisp")

    # ── Step 5: Compare with 05d domains ─────────────────────────────────────
    print_step("Step 5: Comparing with 05d RNA-only domains")
    comparison_df = compare_with_05d(adata)
    if len(comparison_df) > 0:
        save_csv(comparison_df, "multisp_domain_comparison.csv", subdir="multisp")

    # ── Step 6: RNA-only ablation ────────────────────────────────────────────
    print_step("Step 6: RNA-only ablation test")
    adata = run_rna_only_ablation(adata, config)

    # Compare bimodal vs RNA-only ablation
    labels_bimodal = adata.obs["multisp_domain"].astype(str).values
    labels_rna = adata.obs["rna_only_domain"].astype(str).values
    ari_ablation = adjusted_rand_score(labels_rna, labels_bimodal)
    nmi_ablation = normalized_mutual_info_score(labels_rna, labels_bimodal)
    print(f"    Bimodal vs RNA-only ablation: ARI={ari_ablation:.3f}, NMI={nmi_ablation:.3f}")

    # ── Step 7: Modality variation analysis ──────────────────────────────────
    print_step("Step 7: Modality variation per domain")
    modality_records = []
    for domain in sorted(adata.obs["multisp_domain"].unique()):
        mask = adata.obs["multisp_domain"] == domain
        n_spots = mask.sum()

        # RNA variance in this domain
        rna_sub = adata.X[mask.values] if hasattr(mask, 'values') else adata.X[mask]
        if sparse.issparse(rna_sub):
            rna_sub = rna_sub.toarray()
        rna_var = np.var(rna_sub, axis=0).mean()

        # Chromatin variance in this domain
        chrom_sub = adata.obsm["chromatin"][mask.values if hasattr(mask, 'values') else mask]
        chrom_var = np.var(chrom_sub, axis=0).mean()

        modality_records.append({
            "domain": domain, "n_spots": int(n_spots),
            "mean_rna_variance": rna_var,
            "mean_chromatin_variance": chrom_var,
            "chromatin_rna_var_ratio": chrom_var / (rna_var + 1e-10),
        })

    modality_df = pd.DataFrame(modality_records)
    save_csv(modality_df, "multisp_modality_variation.csv", subdir="multisp")
    print(f"\n  Modality variance by domain:")
    for _, row in modality_df.iterrows():
        print(f"    Domain {row['domain']}: {int(row['n_spots'])} spots, "
              f"RNA_var={row['mean_rna_variance']:.4f}, "
              f"chrom_var={row['mean_chromatin_variance']:.4f}")

    # ── Step 8: Save domain assignments ──────────────────────────────────────
    print_step("Step 8: Saving outputs")
    domain_df = adata.obs[["multisp_domain"]].copy()
    if "rna_only_domain" in adata.obs.columns:
        domain_df["rna_only_domain"] = adata.obs["rna_only_domain"]
    if "sample_id" in adata.obs.columns:
        domain_df["sample_id"] = adata.obs["sample_id"]
    if "condition" in adata.obs.columns:
        domain_df["condition"] = adata.obs["condition"]
    save_csv(domain_df, "multisp_domains.csv", subdir="multisp")

    # Save checkpoint with all embeddings
    save_checkpoint(adata, "spatial_with_multisp_domains.h5ad", subdir="multisp")

    # Summary
    n_domains = adata.obs["multisp_domain"].nunique()
    method = "MultiSP" if multisp_success else "manual_bimodal"
    print(f"\n  Summary:")
    print(f"    Method: {method}")
    print(f"    Domains identified: {n_domains}")
    print(f"    Spots: {adata.n_obs}")
    if len(comparison_df) > 0:
        print(f"    vs 05d: ARI={comparison_df['ARI'].iloc[0]:.3f}, "
              f"NMI={comparison_df['NMI'].iloc[0]:.3f}")
    print(f"    vs RNA-only ablation: ARI={ari_ablation:.3f}, NMI={nmi_ablation:.3f}")

    print_header("13b: Complete")


if __name__ == "__main__":
    main()
