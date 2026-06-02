#!/usr/bin/env python3
"""
13a_atac_pseudobulk_projection.py — Project scATAC pseudobulk onto Visium spots.

Computes pseudobulk gene-activity and peak accessibility per cell type from the
scATAC multiome dataset, then imputes chromatin state for each Visium spot using
cell2location deconvolution proportions as weights.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, C2L_PREFIX,
    load_config, load_deconvolved_adata,
    save_checkpoint, save_csv, print_header, print_step,
    harmonize_cell_types, strip_c2l_prefix,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_scatac_adata(config):
    """Load scATAC processed h5ad from snapatac2 pipeline."""
    atac_path = PROJECT_ROOT / "Analysis/ATAC/Human_Multiome/results/snapatac2/snapatac2_processed.h5ad"
    print_step(f"Loading scATAC data: {atac_path}")
    adata_atac = sc.read_h5ad(atac_path)
    print(f"    Shape: {adata_atac.shape}")
    print(f"    Cell types: {adata_atac.obs['cell_type'].nunique() if 'cell_type' in adata_atac.obs.columns else 'N/A'}")
    return adata_atac


def compute_gene_activity_pseudobulk(adata_atac, config):
    """Compute per-cell-type pseudobulk gene-activity profiles (CPM-normalized).

    If gene_activity is already stored in layers or obsm, use it directly.
    Otherwise, attempt to compute via snapatac2 make_gene_matrix.
    """
    multisp_config = config["multisp"]
    norm_method = multisp_config["gene_activity_normalization"]

    # Determine cell type column
    ct_col = None
    for candidate in ["cell_type", "celltype", "cluster", "leiden"]:
        if candidate in adata_atac.obs.columns:
            ct_col = candidate
            break
    if ct_col is None:
        raise ValueError("Cannot find cell type column in scATAC adata.obs")
    print_step(f"Using cell type column: '{ct_col}' ({adata_atac.obs[ct_col].nunique()} types)")

    cell_types = adata_atac.obs[ct_col].unique().tolist()

    # Check for existing gene activity matrix
    gene_activity_layer = None
    if "gene_activity" in adata_atac.layers:
        gene_activity_layer = "gene_activity"
    elif "GeneActivity" in adata_atac.layers:
        gene_activity_layer = "GeneActivity"
    elif "gene_score" in adata_atac.layers:
        gene_activity_layer = "gene_score"

    if gene_activity_layer is not None:
        print_step(f"Using existing gene activity layer: '{gene_activity_layer}'")
        mat = adata_atac.layers[gene_activity_layer]
    elif "gene_activity" in adata_atac.obsm:
        print_step("Using gene activity from .obsm['gene_activity']")
        # obsm stores as a separate AnnData or matrix; handle both
        ga = adata_atac.obsm["gene_activity"]
        if isinstance(ga, ad.AnnData):
            mat = ga.X
            # Use the gene names from the sub-AnnData
            adata_atac_ga = ga
        else:
            mat = ga
        # Cannot easily build pseudobulk without gene names from obsm matrix
        # Fall back to .X
        if not isinstance(ga, ad.AnnData):
            print_step("WARNING: obsm gene_activity is a raw matrix, falling back to .X")
            mat = adata_atac.X
    else:
        # Try snapatac2 gene matrix computation
        print_step("No pre-computed gene activity found; attempting snapatac2.pp.make_gene_matrix()")
        try:
            import snapatac2 as snap
            gene_mat = snap.pp.make_gene_matrix(adata_atac, gene_anno=None)
            mat = gene_mat.X
            print(f"    Gene matrix shape: {gene_mat.shape}")
        except Exception as e:
            print(f"    WARNING: snapatac2 gene matrix failed ({e}), using .X as fallback")
            mat = adata_atac.X

    # Build pseudobulk per cell type
    if sparse.issparse(mat):
        mat_dense = None  # will convert per-cell-type to save memory
    else:
        mat_dense = mat

    pseudobulk_records = {}
    for ct in cell_types:
        mask = adata_atac.obs[ct_col] == ct
        n_cells = mask.sum()
        if n_cells < 10:
            print(f"    Skipping {ct}: only {n_cells} cells")
            continue

        mask_arr = mask.values if hasattr(mask, 'values') else np.asarray(mask)
        if sparse.issparse(mat):
            ct_sum = np.asarray(mat[mask_arr].sum(axis=0)).flatten()
        else:
            ct_sum = mat_dense[mask_arr].sum(axis=0).flatten()

        # CPM normalization
        if norm_method == "cpm":
            total = ct_sum.sum()
            if total > 0:
                ct_sum = ct_sum / total * 1e6

        pseudobulk_records[ct] = ct_sum
        print(f"    {ct}: {n_cells} cells, {(ct_sum > 0).sum()} expressed features")

    # Build DataFrame: rows = features, columns = cell types
    feature_names = adata_atac.var_names.tolist()
    pseudobulk_df = pd.DataFrame(pseudobulk_records, index=feature_names)

    print_step(f"Pseudobulk gene activity: {pseudobulk_df.shape[0]} features x {pseudobulk_df.shape[1]} cell types")
    return pseudobulk_df, ct_col


def compute_peak_pseudobulk(adata_atac, ct_col, config):
    """Compute per-cell-type pseudobulk peak accessibility (TF-IDF normalized).

    Selects top variable peaks + DA peaks if available.
    """
    multisp_config = config["multisp"]
    n_top = multisp_config["n_top_peaks"]

    # Peak accessibility is in .X for scATAC data
    mat = adata_atac.X
    cell_types = adata_atac.obs[ct_col].unique().tolist()

    # Select top variable peaks
    if "highly_variable" in adata_atac.var.columns:
        peak_mask = adata_atac.var["highly_variable"].values
        print_step(f"Using {peak_mask.sum()} pre-selected variable peaks")
    else:
        # Compute variance and select top N
        if sparse.issparse(mat):
            var_per_peak = np.asarray(mat.power(2).mean(axis=0) - np.power(mat.mean(axis=0), 2)).flatten()
        else:
            var_per_peak = mat.var(axis=0)
        top_idx = np.argsort(var_per_peak)[-n_top:]
        peak_mask = np.zeros(adata_atac.n_vars, dtype=bool)
        peak_mask[top_idx] = True
        print_step(f"Selected top {peak_mask.sum()} variable peaks")

    # Include DA peaks if available
    da_path = PROJECT_ROOT / "Analysis/ATAC/Human_Multiome/results/snapatac2/scatac_da_results.csv"
    if da_path.exists():
        da_df = pd.read_csv(da_path)
        if "peak" in da_df.columns:
            da_peaks = set(da_df["peak"].dropna().unique())
            peak_names = adata_atac.var_names.tolist()
            for i, p in enumerate(peak_names):
                if p in da_peaks:
                    peak_mask[i] = True
            print_step(f"Added DA peaks: total {peak_mask.sum()} peaks selected")

    mat_sub = mat[:, peak_mask]
    peak_names_sub = adata_atac.var_names[peak_mask].tolist()

    # TF-IDF normalization
    if sparse.issparse(mat_sub):
        mat_sub = mat_sub.toarray()

    # Term frequency: normalize each cell
    tf = mat_sub / (mat_sub.sum(axis=1, keepdims=True) + 1e-8)
    # Inverse document frequency: log(N / (1 + df))
    n_cells = mat_sub.shape[0]
    df = (mat_sub > 0).sum(axis=0)
    idf = np.log1p(n_cells / (1 + df))
    tfidf = tf * idf

    # Build pseudobulk per cell type
    pseudobulk_peaks = {}
    for ct in cell_types:
        mask = adata_atac.obs[ct_col] == ct
        n_cells_ct = mask.sum()
        if n_cells_ct < 10:
            continue
        ct_mean = tfidf[mask.values if hasattr(mask, 'values') else mask].mean(axis=0)
        pseudobulk_peaks[ct] = ct_mean

    peak_df = pd.DataFrame(pseudobulk_peaks, index=peak_names_sub)
    print_step(f"Pseudobulk peak accessibility: {peak_df.shape[0]} peaks x {peak_df.shape[1]} cell types")
    return peak_df


def project_onto_spots(pseudobulk_df, adata_spatial, config):
    """Project pseudobulk profiles onto Visium spots using c2l proportions.

    For each spot, imputed_profile = sum(c2l_proportion_ct * pseudobulk_ct).
    """
    # Get c2l proportions — check both prefix variants
    c2l_cols = [c for c in adata_spatial.obs.columns if c.startswith(C2L_PREFIX)]
    if not c2l_cols:
        # Try with "c2l_" prepended (some versions add this)
        full_prefix = "c2l_" + C2L_PREFIX
        c2l_cols = [c for c in adata_spatial.obs.columns if c.startswith(full_prefix)]
    if not c2l_cols:
        # Last resort: try obsm
        if "q05_cell_abundance_w_sf" in adata_spatial.obsm:
            print("  Using obsm['q05_cell_abundance_w_sf'] for proportions")
            obsm_df = pd.DataFrame(
                adata_spatial.obsm["q05_cell_abundance_w_sf"],
                index=adata_spatial.obs_names,
            )
            # Store as obs columns with C2L_PREFIX
            for col in obsm_df.columns:
                adata_spatial.obs[C2L_PREFIX + str(col)] = obsm_df[col].values
            c2l_cols = [c for c in adata_spatial.obs.columns if c.startswith(C2L_PREFIX)]
    if not c2l_cols:
        raise ValueError("No cell2location proportion columns found in spatial adata")

    # Map c2l cell type names to scATAC cell type names
    # strip_c2l_prefix handles both "C2L_PREFIX" and "c2l_C2L_PREFIX"
    c2l_names = []
    for c in c2l_cols:
        name = c
        if name.startswith("c2l_" + C2L_PREFIX):
            name = name[len("c2l_" + C2L_PREFIX):]
        elif name.startswith(C2L_PREFIX):
            name = name[len(C2L_PREFIX):]
        c2l_names.append(name)
    atac_names = pseudobulk_df.columns.tolist()

    # Harmonize: map scATAC names to c2l names
    atac_harmonized = harmonize_cell_types(atac_names, source="atac")
    atac_to_c2l = dict(zip(atac_names, atac_harmonized))

    # Find matching cell types
    matched = {}
    for atac_ct, c2l_ct in atac_to_c2l.items():
        if c2l_ct in c2l_names:
            c2l_col = c2l_cols[c2l_names.index(c2l_ct)]
            matched[atac_ct] = c2l_col
    print_step(f"Matched cell types: {len(matched)}/{len(atac_names)}")
    for atac_ct, c2l_col in matched.items():
        print(f"    {atac_ct} -> {strip_c2l_prefix(c2l_col)}")

    if len(matched) == 0:
        raise ValueError("No cell types could be matched between scATAC and c2l")

    # Normalize c2l proportions for matched cell types only
    prop_cols = list(matched.values())
    props = adata_spatial.obs[prop_cols].values.copy()
    prop_sums = props.sum(axis=1, keepdims=True)
    prop_sums[prop_sums == 0] = 1.0
    props = props / prop_sums

    # Get pseudobulk matrix for matched cell types (features x matched_cts)
    matched_atac_cts = list(matched.keys())

    # CRITICAL: restrict to features that overlap with spatial RNA gene names
    # to avoid 6M peak features → 296GB allocation
    spatial_genes = set(adata_spatial.var_names)
    overlap_features = [f for f in pseudobulk_df.index if f in spatial_genes]
    if len(overlap_features) == 0:
        print("  WARNING: No feature overlap between pseudobulk and spatial RNA")
        print(f"    Pseudobulk features (first 5): {pseudobulk_df.index[:5].tolist()}")
        print(f"    Spatial genes (first 5): {list(spatial_genes)[:5]}")
        # Fall back to top variable features capped at 5000
        feature_var = pseudobulk_df[matched_atac_cts].var(axis=1)
        overlap_features = feature_var.nlargest(5000).index.tolist()
        print(f"    Using top 5000 variable features instead")

    pseudobulk_df = pseudobulk_df.loc[overlap_features]
    print_step(f"Restricted to {len(overlap_features)} overlapping features "
               f"(from {len(pseudobulk_df.index)} total)")

    pb_mat = pseudobulk_df[matched_atac_cts].values  # (n_features, n_matched_cts)

    # Weighted sum: imputed = props @ pb_mat.T -> (n_spots, n_features)
    imputed = props @ pb_mat.T

    print_step(f"Imputed chromatin matrix: {imputed.shape[0]} spots x {imputed.shape[1]} features")
    return imputed, pseudobulk_df.index.tolist()


def validate_imputation(imputed_mat, feature_names, adata_spatial):
    """Validate imputation by correlating imputed gene-activity with RNA expression.

    Tests liver marker genes: ALB, CYP3A4, GLUL.
    """
    from scipy.stats import spearmanr

    markers = ["ALB", "CYP3A4", "GLUL", "CYP2E1", "CYP1A2", "ASS1", "HAL"]
    spatial_genes = adata_spatial.var_names.tolist()

    # Get RNA expression (use log-normalized if available)
    if "counts" in adata_spatial.layers:
        rna_mat = adata_spatial.X
    else:
        rna_mat = adata_spatial.X

    if sparse.issparse(rna_mat):
        rna_dense = rna_mat.toarray()
    else:
        rna_dense = np.array(rna_mat)

    validation_records = []
    for marker in markers:
        # Check if marker exists in both modalities
        if marker not in feature_names:
            validation_records.append({
                "gene": marker, "rho": np.nan, "pval": np.nan,
                "status": "missing_in_chromatin",
            })
            continue
        if marker not in spatial_genes:
            validation_records.append({
                "gene": marker, "rho": np.nan, "pval": np.nan,
                "status": "missing_in_spatial",
            })
            continue

        feat_idx = feature_names.index(marker)
        gene_idx = spatial_genes.index(marker)

        imputed_vals = imputed_mat[:, feat_idx]
        rna_vals = rna_dense[:, gene_idx]

        # Skip if constant
        if np.std(imputed_vals) < 1e-10 or np.std(rna_vals) < 1e-10:
            validation_records.append({
                "gene": marker, "rho": np.nan, "pval": np.nan,
                "status": "constant_expression",
            })
            continue

        rho, pval = spearmanr(imputed_vals, rna_vals)
        validation_records.append({
            "gene": marker, "rho": rho, "pval": pval, "status": "tested",
        })
        print(f"    {marker}: rho={rho:.3f}, p={pval:.2e}")

    return pd.DataFrame(validation_records)


def main():
    print_header("13a: scATAC Pseudobulk Projection onto Visium Spots")

    config = load_config()
    output_dir = RESULTS_DIR / "multisp"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load scATAC data ─────────────────────────────────────────────
    print_step("Step 1: Loading scATAC data")
    adata_atac = load_scatac_adata(config)

    # ── Step 2: Compute pseudobulk gene activity ─────────────────────────────
    print_step("Step 2: Computing pseudobulk gene activity (CPM-normalized)")
    pseudobulk_ga, ct_col = compute_gene_activity_pseudobulk(adata_atac, config)

    # ── Step 3: Compute pseudobulk peak accessibility ────────────────────────
    print_step("Step 3: Computing pseudobulk peak accessibility (TF-IDF)")
    try:
        pseudobulk_peaks = compute_peak_pseudobulk(adata_atac, ct_col, config)
        save_csv(pseudobulk_peaks, "pseudobulk_peak_accessibility.csv", subdir="multisp")
    except Exception as e:
        print(f"    WARNING: Peak pseudobulk failed ({e}), skipping peak modality")
        pseudobulk_peaks = None

    # ── Step 4: Load spatial data with c2l proportions ───────────────────────
    print_step("Step 4: Loading deconvolved spatial data")
    adata_spatial = load_deconvolved_adata()
    print(f"    Spatial data: {adata_spatial.n_obs} spots x {adata_spatial.n_vars} genes")

    # ── Step 5: Project pseudobulk onto spots ────────────────────────────────
    print_step("Step 5: Projecting pseudobulk gene activity onto spots")
    imputed_mat, feature_names = project_onto_spots(pseudobulk_ga, adata_spatial, config)

    # Build imputed chromatin AnnData
    adata_chromatin = ad.AnnData(
        X=imputed_mat.astype(np.float32),
        obs=adata_spatial.obs.copy(),
        var=pd.DataFrame(index=feature_names),
    )
    # Copy spatial coordinates if available
    if "spatial" in adata_spatial.obsm:
        adata_chromatin.obsm["spatial"] = adata_spatial.obsm["spatial"].copy()

    # ── Step 6: Validate imputation ──────────────────────────────────────────
    print_step("Step 6: Validating imputation against RNA expression")
    validation_df = validate_imputation(imputed_mat, feature_names, adata_spatial)
    save_csv(validation_df, "projection_validation.csv", subdir="multisp")

    tested = validation_df[validation_df["status"] == "tested"]
    if len(tested) > 0:
        mean_rho = tested["rho"].mean()
        n_pos = (tested["rho"] > 0).sum()
        print(f"    Mean Spearman rho: {mean_rho:.3f} ({n_pos}/{len(tested)} positive)")

    # ── Step 7: Save outputs ─────────────────────────────────────────────────
    print_step("Step 7: Saving outputs")
    save_checkpoint(adata_chromatin, "imputed_chromatin.h5ad", subdir="multisp")
    save_csv(pseudobulk_ga, "pseudobulk_gene_activity.csv", subdir="multisp")

    # Summary
    print(f"\n  Summary:")
    print(f"    scATAC cells: {adata_atac.n_obs}")
    print(f"    Cell types projected: {pseudobulk_ga.shape[1]}")
    print(f"    Features imputed: {pseudobulk_ga.shape[0]}")
    print(f"    Visium spots: {adata_spatial.n_obs}")
    if pseudobulk_peaks is not None:
        print(f"    Peaks profiled: {pseudobulk_peaks.shape[0]}")

    print_header("13a: Complete")


if __name__ == "__main__":
    main()
