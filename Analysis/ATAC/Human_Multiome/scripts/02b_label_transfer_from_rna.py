#!/usr/bin/env python3
"""
02b_label_transfer_from_rna.py — Transfer cell type labels from donor-matched
snRNA-seq atlas to scATAC-seq data via KNN in PCA space.

For each of 18 donors (GSE244832):
  1. Compute ATAC gene activity matrix from per-donor backed h5ads (SnapATAC2)
  2. Extract matching RNA cells (CellTypist-annotated) from integrated snRNA atlas
  3. Project ATAC gene activity into RNA PCA space (centering + loadings)
  4. Transfer cell_type labels via k=10 nearest-neighbor majority voting

Saves:
  - snapatac2_label_transferred.h5ad with cell_type_transferred column
  - gene_activity_matrix.h5ad (cached for re-use)
  - Confusion matrix, Jaccard, concordance statistics
  - Per-cell annotation comparison table

Pipeline position:
  02_snapatac2_processing.py → **02b_label_transfer_from_rna.py** → 02c_rerun_peaks

Environment:
  micromamba activate snapatac2

Usage:
  cd Analysis/ATAC/Human_Multiome
  python scripts/02b_label_transfer_from_rna.py \
      --atac_dir results/snapatac2 \
      --rna_atlas ../../Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad \
      --donor_pairing ../../data/GSE244832/metadata/donor_pairing.csv \
      --output_dir results/label_transfer
"""

from __future__ import annotations

import argparse
import gc
import json
import logging
import os
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import confusion_matrix
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore", category=FutureWarning)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Matplotlib backend (must be before import)
# ---------------------------------------------------------------------------
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(
        description="Transfer cell type labels from donor-matched snRNA-seq to scATAC-seq"
    )
    p.add_argument(
        "--atac_dir",
        default="results/snapatac2",
        help="SnapATAC2 results directory with per_donor/ and snapatac2_processed.h5ad",
    )
    p.add_argument(
        "--rna_atlas",
        default="../../Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad",
        help="Path to integrated snRNA-seq atlas h5ad (CellTypist-annotated)",
    )
    p.add_argument(
        "--donor_pairing",
        default="../../data/GSE244832/metadata/donor_pairing.csv",
        help="Donor pairing CSV (donor_id, rna_gsm, etc.)",
    )
    p.add_argument(
        "--output_dir",
        default="results/label_transfer",
        help="Output directory for label transfer results",
    )
    p.add_argument(
        "--k_neighbors",
        type=int,
        default=10,
        help="Number of neighbors for KNN label transfer",
    )
    p.add_argument(
        "--confidence_percentile",
        type=float,
        default=95.0,
        help="Percentile threshold for low-confidence flagging (on KNN distance)",
    )
    p.add_argument(
        "--n_hvgs",
        type=int,
        default=3000,
        help="Number of HVGs for label transfer",
    )
    p.add_argument(
        "--n_pcs", type=int, default=30, help="Number of PCs for PCA/ingest"
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Step 1: Compute gene activity matrix from ATAC fragments
# ---------------------------------------------------------------------------
def compute_gene_activity(atac_dir: str) -> ad.AnnData:
    """Compute gene activity matrix from per-donor backed h5ads via SnapATAC2.

    Uses snap.pp.make_gene_matrix() which sums fragment overlaps with gene bodies.
    Requires AnnDataSet with fragment references (from per_donor/*.h5ad).
    """
    import snapatac2 as snap

    per_donor_dir = Path(atac_dir) / "per_donor"
    h5ad_files = sorted(per_donor_dir.glob("*.h5ad"))

    if not h5ad_files:
        raise FileNotFoundError(f"No h5ad files found in {per_donor_dir}")

    log.info("Found %d per-donor h5ad files", len(h5ad_files))

    # Open backed h5ads for AnnDataSet
    # SnapATAC2 AnnDataSet expects list of (name_str, AnnData) tuples
    processed = []
    for f in h5ad_files:
        adata_d = snap.read(str(f), backed="r")
        processed.append((f.stem, adata_d))
        log.info("  %s: %d cells", f.stem, adata_d.n_obs)

    # Create temporary AnnDataSet
    temp_h5ads = str(Path(atac_dir) / "temp_label_transfer.h5ads")
    dataset = snap.AnnDataSet(
        adatas=processed,
        filename=temp_h5ads,
        add_key="donor_id",
    )
    log.info(
        "AnnDataSet created: %d cells x %d features", dataset.n_obs, dataset.n_vars
    )

    # Compute gene activity (uses hg38 built-in annotation, same as original pipeline)
    log.info("Computing gene activity matrix from fragment files...")
    gene_mat = snap.pp.make_gene_matrix(dataset, gene_anno=snap.genome.hg38)
    log.info("Gene activity matrix: %d cells x %d genes", gene_mat.n_obs, gene_mat.n_vars)

    # Convert to in-memory if backed
    if hasattr(gene_mat, "to_memory"):
        gene_mat = gene_mat.to_memory()

    # Clean up
    dataset.close()
    try:
        Path(temp_h5ads).unlink(missing_ok=True)
    except Exception:
        pass
    for _, a in processed:
        try:
            if hasattr(a, "file") and a.file:
                a.file.close()
        except Exception:
            pass

    return gene_mat


# ---------------------------------------------------------------------------
# Step 2: Load RNA reference (atlas subset for GSE244832)
# ---------------------------------------------------------------------------
def load_rna_reference(
    atlas_path: str, donor_pairing_path: str
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Load RNA reference for label transfer.

    First checks for a pre-extracted GSE244832 reference file (created by
    02a_extract_rna_reference.py in rapids_singlecell env to avoid anndata
    version incompatibilities). Falls back to loading the full atlas if the
    pre-extracted file is not found.
    """
    pairing = pd.read_csv(donor_pairing_path)

    # Check for pre-extracted reference (avoids anndata version mismatch)
    pre_extracted = Path(atlas_path).parent / "rna_reference_GSE244832.h5ad"
    # Also check in label_transfer output dir
    alt_pre_extracted = Path("results/label_transfer/rna_reference_GSE244832.h5ad")

    for ref_path in [pre_extracted, alt_pre_extracted]:
        if ref_path.exists():
            log.info("Loading pre-extracted reference: %s", ref_path)
            adata_ref = sc.read_h5ad(str(ref_path))
            log.info(
                "Reference loaded: %d cells, %d cell types",
                adata_ref.n_obs,
                adata_ref.obs["cell_type"].nunique(),
            )
            log.info(
                "Cell type distribution:\n%s",
                adata_ref.obs["cell_type"].value_counts().to_string(),
            )
            return adata_ref, pairing

    # No pre-extracted file — try loading full atlas
    log.info("No pre-extracted reference found, loading full atlas: %s", atlas_path)
    log.info("(If this fails with encoding errors, run 02a_extract_rna_reference.py first)")

    try:
        adata_rna = sc.read_h5ad(atlas_path, backed="r")
    except Exception as e:
        log.error("Failed to load atlas in backed mode: %s", e)
        log.info("Trying non-backed mode...")
        adata_rna = sc.read_h5ad(atlas_path)

    log.info("Atlas: %d cells x %d genes", adata_rna.n_obs, adata_rna.n_vars)
    log.info("Atlas obs columns: %s", list(adata_rna.obs.columns))

    # --- Identify GSE244832 cells ---
    gse_mask = None

    if "dataset" in adata_rna.obs.columns:
        for pattern in ["GSE244832", "gse244832"]:
            candidate = adata_rna.obs["dataset"].str.contains(pattern, case=False, na=False)
            if candidate.sum() > 0:
                gse_mask = candidate
                log.info("Matched %d cells via dataset column ('%s')", gse_mask.sum(), pattern)
                break

    if gse_mask is None or gse_mask.sum() == 0:
        rna_gsms = pairing["rna_gsm"].dropna().unique().tolist()
        log.info("Trying GSM filter with %d accessions...", len(rna_gsms))
        gse_mask = adata_rna.obs["sample"].isin(rna_gsms)
        log.info("Matched %d cells via GSM accessions", gse_mask.sum())

    if gse_mask is None or gse_mask.sum() == 0:
        raise ValueError(
            "No cells found for GSE244832 in atlas. "
            "Check dataset/sample columns and donor_pairing.csv"
        )

    log.info("Loading %d reference cells into memory...", gse_mask.sum())
    adata_ref = adata_rna[gse_mask].to_memory()

    if hasattr(adata_rna, "file") and adata_rna.file:
        adata_rna.file.close()
    del adata_rna
    gc.collect()

    log.info(
        "Reference loaded: %d cells, %d cell types",
        adata_ref.n_obs,
        adata_ref.obs["cell_type"].nunique(),
    )
    log.info("Cell type distribution:\n%s", adata_ref.obs["cell_type"].value_counts().to_string())

    return adata_ref, pairing


# ---------------------------------------------------------------------------
# Step 3: Build donor -> RNA cells mapping
# ---------------------------------------------------------------------------
def build_donor_map(
    adata_ref: ad.AnnData, pairing: pd.DataFrame
) -> dict[str, list[str]]:
    """Map donor_id -> list of RNA cell barcodes in the reference atlas.

    The atlas uses SRR accessions as sample names (not GSM). The donor_pairing
    CSV has `rna_srrs` (semicolon-separated SRR accessions) which is the correct
    matching column. Falls back to `rna_gsm` if SRR matching fails.
    """
    donor_to_cells = {}
    ref_samples = set(adata_ref.obs["sample"].unique())
    log.info("Reference has %d unique sample IDs", len(ref_samples))
    log.info("  Sample format examples: %s", list(ref_samples)[:3])

    for _, row in pairing.iterrows():
        donor_id = row["donor_id"]
        mask = pd.Series(False, index=adata_ref.obs_names)

        # Strategy 1: Match by SRR accessions (rna_srrs column, semicolon-separated)
        rna_srrs = row.get("rna_srrs", "")
        if pd.notna(rna_srrs) and rna_srrs:
            srrs = [s.strip() for s in str(rna_srrs).split(";") if s.strip()]
            mask = adata_ref.obs["sample"].isin(srrs)
            if mask.sum() > 0:
                donor_to_cells[donor_id] = adata_ref.obs_names[mask].tolist()
                log.info("  %s -> %d RNA cells (%d SRR samples)", donor_id, mask.sum(), len(srrs))
                continue

        # Strategy 2: Match by GSM accession
        rna_gsm = row.get("rna_gsm", "")
        if pd.notna(rna_gsm) and rna_gsm:
            gsms = [g.strip() for g in str(rna_gsm).split(";")]
            mask = adata_ref.obs["sample"].isin(gsms)
            if mask.sum() > 0:
                donor_to_cells[donor_id] = adata_ref.obs_names[mask].tolist()
                log.info("  %s -> %d RNA cells (GSM: %s)", donor_id, mask.sum(), rna_gsm)
                continue

        # Strategy 3: Match by rna_sample name (JB_*) — exact match only
        rna_sample = row.get("rna_sample", "")
        if pd.notna(rna_sample) and rna_sample:
            mask = adata_ref.obs["sample"] == str(rna_sample)
            if mask.sum() > 0:
                donor_to_cells[donor_id] = adata_ref.obs_names[mask].tolist()
                log.info("  %s -> %d RNA cells (sample: %s)", donor_id, mask.sum(), rna_sample)
                continue

        log.warning("  %s: No RNA cells found (srrs=%s, gsm=%s)", donor_id, rna_srrs[:30], rna_gsm)

    return donor_to_cells


# ---------------------------------------------------------------------------
# Step 4: Per-donor label transfer
# ---------------------------------------------------------------------------
def transfer_labels_per_donor(
    adata_atac_ga: ad.AnnData,
    adata_ref: ad.AnnData,
    donor_to_rna_cells: dict[str, list[str]],
    n_hvgs: int = 3000,
    n_pcs: int = 30,
    k_neighbors: int = 10,
) -> tuple[pd.Series, pd.Series]:
    """Per-donor label transfer using scanpy.tl.ingest.

    For each donor:
      1. Subset RNA reference and ATAC gene activity to this donor's cells
      2. Find shared genes, normalize, select HVGs
      3. PCA on RNA reference, compute neighbors
      4. scanpy.tl.ingest to project ATAC into RNA PCA space + transfer labels
      5. Compute KNN distances for confidence scoring

    LIMITATION (documented): transfer is per-donor KNN in the RNA reference's
    PCA space (query centered with the reference mean and projected through the
    reference loadings). There is NO joint scaling / batch integration (e.g.
    Harmony, scVI, Seurat anchors) of the ATAC gene-activity modality against
    the RNA modality, so cross-modality distribution shift is not corrected.
    The relabel fraction vs the prior gene-activity annotation is therefore
    reported downstream as a transparency metric.
    """
    transferred_labels = pd.Series(index=adata_atac_ga.obs_names, dtype="object")
    transfer_distances = pd.Series(index=adata_atac_ga.obs_names, dtype="float64")

    # Extract donor_id from ATAC cell names (format: "D01:barcode" or "D01_barcode")
    atac_obs = adata_atac_ga.obs_names.to_series()
    if atac_obs.str.contains(":").any():
        atac_donor_ids = atac_obs.str.split(":").str[0]
    elif atac_obs.str.contains("_").any():
        # Try D01_barcode format
        atac_donor_ids = atac_obs.str.extract(r"(D\d+)")[0]
    else:
        # Fallback: check obs columns
        for col in ["donor_id", "donor", "sample"]:
            if col in adata_atac_ga.obs.columns:
                atac_donor_ids = adata_atac_ga.obs[col].astype(str)
                break
        else:
            raise ValueError(
                "Cannot determine donor_id from cell names or obs columns. "
                f"Sample obs_names: {atac_obs.head(3).tolist()}"
            )

    log.info("ATAC donors found: %s", sorted(atac_donor_ids.unique()))

    for donor_id, rna_cells in donor_to_rna_cells.items():
        log.info("\n--- Donor %s ---", donor_id)

        # Get ATAC cells for this donor
        atac_mask = atac_donor_ids == donor_id
        n_atac = atac_mask.sum()
        if n_atac == 0:
            log.warning("  No ATAC cells for %s, skipping", donor_id)
            continue

        log.info("  ATAC: %d cells, RNA: %d cells", n_atac, len(rna_cells))

        # Subset reference and query
        adata_ref_d = adata_ref[rna_cells].copy()
        adata_query_d = adata_atac_ga[atac_mask].copy()

        # Find shared genes
        shared_genes = sorted(
            adata_ref_d.var_names.intersection(adata_query_d.var_names)
        )
        log.info("  Shared genes: %d", len(shared_genes))

        if len(shared_genes) < 500:
            log.warning("  Too few shared genes (%d), skipping", len(shared_genes))
            continue

        # Subset to shared genes
        adata_ref_d = adata_ref_d[:, shared_genes].copy()
        adata_query_d = adata_query_d[:, shared_genes].copy()

        # Ensure CSR format for sparse
        if sp.issparse(adata_ref_d.X):
            adata_ref_d.X = adata_ref_d.X.tocsr()
        if sp.issparse(adata_query_d.X):
            adata_query_d.X = adata_query_d.X.tocsr()

        # --- Normalize ---
        # RNA: check if already log-normalized
        # Strategy: check uns for normalization log, then fall back to data range
        already_normalized = False
        if "log1p" in adata_ref_d.uns:
            already_normalized = True
            log.info("  RNA reference: log1p found in .uns (already normalized)")
        elif hasattr(adata_ref_d, "raw") and adata_ref_d.raw is not None:
            already_normalized = True
            log.info("  RNA reference: .raw exists (X is likely normalized)")
        else:
            # Heuristic: check if values look like integer counts
            ref_max = (
                adata_ref_d.X.max()
                if not sp.issparse(adata_ref_d.X)
                else adata_ref_d.X.data.max() if adata_ref_d.X.nnz > 0 else 0
            )
            # Check for integer-like values (raw counts) vs continuous (normalized)
            if sp.issparse(adata_ref_d.X):
                sample_vals = adata_ref_d.X.data[:min(10000, len(adata_ref_d.X.data))]
            else:
                sample_vals = adata_ref_d.X.ravel()[:10000]
            frac_integer = np.mean(np.abs(sample_vals - np.round(sample_vals)) < 1e-6)
            if frac_integer > 0.9 and ref_max > 50:
                log.info(
                    "  RNA reference: raw counts detected (max=%.0f, %.0f%% integer)",
                    ref_max, 100 * frac_integer,
                )
            else:
                already_normalized = True
                log.info(
                    "  RNA reference: appears normalized (max=%.1f, %.0f%% integer)",
                    ref_max, 100 * frac_integer,
                )

        if not already_normalized:
            sc.pp.normalize_total(adata_ref_d, target_sum=1e4)
            sc.pp.log1p(adata_ref_d)

        # ATAC gene activity: always normalize (raw fragment counts)
        sc.pp.normalize_total(adata_query_d, target_sum=1e4)
        sc.pp.log1p(adata_query_d)

        # --- HVG selection on reference ---
        n_hvg = min(n_hvgs, len(shared_genes))
        try:
            sc.pp.highly_variable_genes(
                adata_ref_d, n_top_genes=n_hvg, flavor="seurat"
            )
        except Exception:
            # Fallback: use variance-based selection
            sc.pp.highly_variable_genes(
                adata_ref_d, n_top_genes=n_hvg, flavor="cell_ranger"
            )
        hvgs = adata_ref_d.var_names[adata_ref_d.var["highly_variable"]].tolist()
        log.info("  HVGs selected: %d", len(hvgs))

        # Subset both to HVGs
        adata_ref_d = adata_ref_d[:, hvgs].copy()
        adata_query_d = adata_query_d[:, hvgs].copy()

        # --- PCA on reference, project query ---
        n_pc = min(n_pcs, min(adata_ref_d.n_obs, len(hvgs)) - 1)
        sc.tl.pca(adata_ref_d, n_comps=n_pc, svd_solver="arpack")

        # Project query into reference PCA space using the same centering
        # that sc.tl.pca used internally (X - mean) @ loadings
        pca_loadings = adata_ref_d.varm["PCs"]  # (n_genes, n_pcs)
        # Use PCA's own centering vector if available; otherwise X.mean(axis=0)
        # is equivalent since no transformation occurs between pca() and here
        if "pca" in adata_ref_d.uns and "mean" in adata_ref_d.uns["pca"]:
            pca_mean = np.array(adata_ref_d.uns["pca"]["mean"])
        else:
            pca_mean = np.array(adata_ref_d.X.mean(axis=0)).flatten() if sp.issparse(adata_ref_d.X) else adata_ref_d.X.mean(axis=0)

        query_X = adata_query_d.X
        if sp.issparse(query_X):
            query_centered = query_X.toarray() - pca_mean
        else:
            query_centered = np.asarray(query_X) - pca_mean
        query_pca = query_centered @ pca_loadings
        ref_pca = adata_ref_d.obsm["X_pca"]

        # --- KNN label transfer (direct, no scanpy.tl.ingest) ---
        try:
            nn = NearestNeighbors(n_neighbors=k_neighbors, metric="euclidean")
            nn.fit(ref_pca)
            distances, indices = nn.kneighbors(query_pca)

            # Distance-weighted KNN voting (avoids alphabetical tie-breaking bias)
            ref_labels = adata_ref_d.obs["cell_type"].values
            donor_labels = []
            for i in range(len(query_pca)):
                neighbor_labels = ref_labels[indices[i]]
                neighbor_dists = distances[i]
                # Inverse-distance weights (add epsilon to avoid division by zero)
                weights = 1.0 / (neighbor_dists + 1e-8)
                # Weighted vote: sum weights per label, pick highest
                label_weights = {}
                for lbl, w in zip(neighbor_labels, weights):
                    label_weights[lbl] = label_weights.get(lbl, 0.0) + w
                donor_labels.append(max(label_weights, key=label_weights.get))
            donor_labels = pd.Series(donor_labels, index=adata_query_d.obs_names)

            transferred_labels.loc[adata_atac_ga.obs_names[atac_mask]] = (
                donor_labels.values
            )

            # Confidence: mean KNN distance (lower = more confident)
            mean_dist = distances.mean(axis=1)
            transfer_distances.loc[adata_atac_ga.obs_names[atac_mask]] = mean_dist

            log.info("  Transferred labels:\n%s", donor_labels.value_counts().to_string())
            log.info(
                "  Mean KNN distance: %.3f +/- %.3f", mean_dist.mean(), mean_dist.std()
            )

        except Exception as e:
            log.error("  Label transfer failed for %s: %s", donor_id, e)
            import traceback
            traceback.print_exc()
            continue

        # Free memory
        del adata_ref_d, adata_query_d
        gc.collect()

    return transferred_labels, transfer_distances


# ---------------------------------------------------------------------------
# Step 5: Confidence thresholding
# ---------------------------------------------------------------------------
def apply_confidence_threshold(
    labels: pd.Series, distances: pd.Series, percentile: float = 95.0
) -> tuple[float, pd.Series]:
    """Flag low-confidence cells using median + 3*MAD outlier detection.

    Uses robust statistics (median absolute deviation) rather than a fixed
    percentile, so the number of flagged cells reflects genuine outliers
    rather than always being a fixed fraction.
    The percentile parameter is retained as a ceiling: cells beyond both
    the MAD threshold and the percentile are flagged.
    """
    valid = distances.notna() & (distances > 0)
    if valid.sum() == 0:
        return np.inf, pd.Series(False, index=labels.index)

    d = distances[valid]
    median_d = np.median(d)
    mad = np.median(np.abs(d - median_d))
    mad_threshold = median_d + 3 * mad * 1.4826  # 1.4826 scales MAD to SD

    # Also compute percentile as a ceiling
    pct_threshold = np.percentile(d, percentile)

    # Use the more conservative (lower) threshold
    threshold = min(mad_threshold, pct_threshold)
    low_conf = distances > threshold

    log.info(
        "Confidence threshold: %.3f (MAD=%.3f at median+3MAD=%.3f, p%.0f=%.3f)",
        threshold, mad, mad_threshold, percentile, pct_threshold,
    )
    log.info("Low-confidence cells: %d (%.1f%%)", low_conf.sum(), 100 * low_conf.mean())

    return threshold, low_conf


# ---------------------------------------------------------------------------
# Step 6: Compare old vs new annotations
# ---------------------------------------------------------------------------
def _harmonize_label(label: str) -> str:
    """Map cell type labels to shared major-type vocabulary."""
    _MAP = {
        "Hepatocyte": "Hepatocyte", "Hepatocytes": "Hepatocyte",
        "Cholangiocyte": "Cholangiocyte", "Cholangiocytes": "Cholangiocyte",
        "Stellate_Cell": "Mesenchymal", "Fibroblasts": "Mesenchymal",
        "Endothelial": "Endothelial", "LSEC": "Endothelial",
        "Endothelial cells": "Endothelial", "Endothelial_cells": "Endothelial",
        "Kupffer_Cell": "Myeloid", "Kupffer cells": "Myeloid",
        "Macrophage": "Myeloid", "Macrophages": "Myeloid",
        "Mono+mono derived cells": "Myeloid", "Mono+mono_derived_cells": "Myeloid",
        "NK_T_Cell": "Lymphoid", "B_Cell": "Lymphoid", "B cells": "Lymphoid",
        "T cells": "Lymphoid", "T_cells": "Lymphoid",
        "Circulating NK/NKT": "Lymphoid", "Circulating_NK_NKT": "Lymphoid",
        "Resident NK": "Lymphoid", "Resident_NK": "Lymphoid",
        "Plasma_Cell": "Plasma_cell", "Plasma cells": "Plasma_cell",
        "Plasma_cells": "Plasma_cell",
        "cDC1s": "Myeloid", "cDC2s": "Myeloid", "pDCs": "Myeloid",
        "Neutrophils": "Myeloid", "Basophils": "Myeloid",
    }
    return _MAP.get(label, label)


def compare_annotations(
    adata: ad.AnnData,
    old_col: str = "cell_type_gene_activity",
    new_col: str = "cell_type_transferred",
    output_dir: Path | None = None,
) -> dict:
    """Compare gene-activity vs label-transfer annotations.

    Reports both raw-label and harmonized major-type concordance.
    Raw labels differ by naming convention (e.g., Hepatocyte vs Hepatocytes),
    so raw concordance is expected to be ~0. Harmonized concordance is the
    scientifically meaningful metric.
    """
    old_raw = adata.obs[old_col].astype(str)
    new_raw = adata.obs[new_col].astype(str)
    old_harm = old_raw.map(_harmonize_label)
    new_harm = new_raw.map(_harmonize_label)

    # Harmonized concordance (primary metric)
    concordance_harm = (old_harm == new_harm).mean()
    n_changed_harm = (old_harm != new_harm).sum()
    log.info("Harmonized concordance: %.1f%%", 100 * concordance_harm)

    # Raw concordance (expected ~0 due to naming differences)
    concordance_raw = (old_raw == new_raw).mean()

    # --- Raw cross-mapping confusion matrix (old rows → new columns) ---
    old_with_cells = sorted([l for l in old_raw.unique() if (old_raw == l).sum() > 0])
    new_with_cells = sorted([l for l in new_raw.unique() if (new_raw == l).sum() > 0])
    cm_cross = pd.DataFrame(0, index=old_with_cells, columns=new_with_cells)
    for ol in old_with_cells:
        for nl in new_with_cells:
            cm_cross.loc[ol, nl] = ((old_raw == ol) & (new_raw == nl)).sum()

    # --- Harmonized confusion matrix ---
    major_types = sorted(set(old_harm.unique()) | set(new_harm.unique()))
    cm_harm = confusion_matrix(old_harm, new_harm, labels=major_types)
    cm_harm_df = pd.DataFrame(cm_harm, index=major_types, columns=major_types)

    # Per-major-type Jaccard
    jaccard_harm = {}
    for mt in major_types:
        old_set = set(adata.obs_names[old_harm == mt])
        new_set = set(adata.obs_names[new_harm == mt])
        union = len(old_set | new_set)
        jaccard_harm[mt] = len(old_set & new_set) / union if union > 0 else 0.0

    results = {
        "concordance_harmonized": float(concordance_harm),
        "concordance_raw": float(concordance_raw),
        "n_cells": int(len(adata)),
        "n_changed_harmonized": int(n_changed_harm),
        "frac_changed_harmonized": float(n_changed_harm / len(adata)),
        "jaccard_harmonized": {k: float(v) for k, v in jaccard_harm.items()},
        "old_cell_types": sorted(old_raw.unique().tolist()),
        "new_cell_types": sorted(new_raw.unique().tolist()),
    }

    if output_dir is not None:
        output_dir = Path(output_dir)

        # Save cross-mapping and harmonized confusion matrices
        cm_cross.to_csv(output_dir / "confusion_matrix_crossmap.csv")
        cm_harm_df.to_csv(output_dir / "confusion_matrix_harmonized.csv")

        # --- Harmonized heatmap (the meaningful one) ---
        fig, ax = plt.subplots(figsize=(10, 8))
        cm_norm = cm_harm_df.div(cm_harm_df.sum(axis=1), axis=0).fillna(0)
        sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="Blues", ax=ax,
                    xticklabels=True, yticklabels=True, linewidths=0.5)
        ax.set_xlabel("Label Transfer (harmonized)", fontsize=12)
        ax.set_ylabel("Gene Activity (harmonized)", fontsize=12)
        ax.set_title(f"Cell Type Concordance (harmonized = {concordance_harm:.1%})", fontsize=13)
        plt.tight_layout()
        fig.savefig(output_dir / "confusion_matrix_harmonized.pdf", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # --- Cross-mapping heatmap (raw labels, informative for tracing) ---
        fig, ax = plt.subplots(figsize=(12, 8))
        cm_cross_norm = cm_cross.div(cm_cross.sum(axis=1), axis=0).fillna(0)
        sns.heatmap(cm_cross_norm, annot=True, fmt=".2f", cmap="YlOrRd", ax=ax,
                    xticklabels=True, yticklabels=True, linewidths=0.5)
        ax.set_xlabel("Label Transfer cell type", fontsize=12)
        ax.set_ylabel("Gene Activity cell type", fontsize=12)
        ax.set_title("Cross-Mapping: Gene Activity → Label Transfer", fontsize=13)
        plt.tight_layout()
        fig.savefig(output_dir / "confusion_matrix_crossmap.pdf", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # --- Harmonized Jaccard bar plot ---
        fig, ax = plt.subplots(figsize=(8, 5))
        jac_s = pd.Series(jaccard_harm).sort_values()
        colors = ["#2166ac" if v > 0.3 else "#d6604d" for v in jac_s]
        jac_s.plot.barh(ax=ax, color=colors)
        ax.set_xlabel("Jaccard Index (harmonized)")
        ax.set_title(f"Per-Major-Type Concordance (overall = {concordance_harm:.1%})")
        ax.axvline(0.5, color="black", linestyle="--", alpha=0.3)
        for i, (mt, val) in enumerate(jac_s.items()):
            ax.text(val + 0.01, i, f"{val:.3f}", va="center", fontsize=9)
        plt.tight_layout()
        fig.savefig(output_dir / "jaccard_harmonized.pdf", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # --- Summary JSON ---
        with open(output_dir / "comparison_summary.json", "w") as f:
            json.dump(results, f, indent=2)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # File logging
    fh = logging.FileHandler(output_dir / "label_transfer.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    log.addHandler(fh)

    log.info("=" * 60)
    log.info("Label Transfer from Donor-Matched snRNA-seq Atlas")
    log.info("=" * 60)
    log.info("Parameters: %s", vars(args))

    # -----------------------------------------------------------------------
    # Step 1: Gene activity matrix
    # -----------------------------------------------------------------------
    log.info("\n=== Step 1: Gene activity matrix ===")
    ga_cache = Path(args.atac_dir) / "gene_activity_matrix.h5ad"

    if ga_cache.exists():
        log.info("Loading cached gene activity from %s", ga_cache)
        adata_atac_ga = sc.read_h5ad(str(ga_cache))
    else:
        adata_atac_ga = compute_gene_activity(args.atac_dir)
        log.info("Caching gene activity to %s", ga_cache)
        adata_atac_ga.write_h5ad(str(ga_cache))

    log.info(
        "Gene activity: %d cells x %d genes", adata_atac_ga.n_obs, adata_atac_ga.n_vars
    )

    # -----------------------------------------------------------------------
    # Step 2: Load processed ATAC (for metadata)
    # -----------------------------------------------------------------------
    log.info("\n=== Step 2: Load ATAC metadata ===")
    atac_h5ad = Path(args.atac_dir) / "snapatac2_processed.h5ad"
    if not atac_h5ad.exists():
        # Try fixed version
        atac_h5ad = Path(args.atac_dir) / "snapatac2_processed_fixed.h5ad"
    adata_atac = sc.read_h5ad(str(atac_h5ad))
    log.info("ATAC processed: %d cells", adata_atac.n_obs)
    log.info(
        "Original cell types:\n%s", adata_atac.obs["cell_type"].value_counts().to_string()
    )

    # Align gene activity with processed ATAC data
    # Both derive from the same AnnDataSet (same 18 donors, same cell order)
    # but obs_names may differ because processed h5ad applied obs_names_make_unique()
    # while gene activity matrix retains the original AnnDataSet names.
    if adata_atac_ga.n_obs == adata_atac.n_obs:
        log.info(
            "Cell counts match (%d). Using positional alignment (same AnnDataSet source).",
            adata_atac_ga.n_obs,
        )
        # Make gene activity obs_names unique to match processed h5ad
        if not adata_atac_ga.obs_names.is_unique:
            adata_atac_ga.obs_names_make_unique()

        # Verify positional alignment by comparing donor prefixes for ALL cells
        ga_names = adata_atac_ga.obs_names.to_series()
        proc_names = adata_atac.obs_names.to_series()

        def _extract_donor(names):
            if names.str.contains(":").any():
                return names.str.split(":").str[0]
            return names.str.extract(r"(D\d+)")[0]

        ga_donors = _extract_donor(ga_names)
        proc_donors = _extract_donor(proc_names)

        if ga_donors.isna().all() or proc_donors.isna().all():
            log.warning("Cannot extract donor IDs from obs_names for alignment check")
        else:
            match_rate = (ga_donors.values == proc_donors.values).mean()
            if match_rate < 0.99:
                log.error(
                    "ALIGNMENT FAILED: only %.1f%% of cells have matching donor prefixes. "
                    "Positional alignment is unreliable — falling back to obs_names intersection.",
                    100 * match_rate,
                )
                # Fall back to explicit alignment
                shared = ga_names[ga_names.isin(proc_names)]
                log.info("Falling back to obs_names intersection: %d shared", len(shared))
                adata_atac_ga = adata_atac_ga[shared.index].copy()
                for col in adata_atac.obs.columns:
                    if col in ("cell_type", "leiden", "condition", "donor_id", "donor", "batch"):
                        adata_atac_ga.obs[col] = adata_atac.obs.loc[shared, col].values
                # Skip the positional copy below
                adata_atac_ga.obs["_alignment_verified"] = True
            else:
                log.info(
                    "Positional alignment verified: %.1f%% donor prefix match (%d/%d cells)",
                    100 * match_rate, int(match_rate * len(ga_names)), len(ga_names),
                )

        # Copy metadata positionally (same cell order from AnnDataSet)
        if "_alignment_verified" not in adata_atac_ga.obs.columns:
            for col in adata_atac.obs.columns:
                if col in ("cell_type", "leiden", "condition", "donor_id", "donor", "batch"):
                    adata_atac_ga.obs[col] = adata_atac.obs[col].values
    else:
        # Different cell counts — fall back to obs_names intersection
        log.warning(
            "Cell count mismatch: gene_activity=%d, processed=%d. Using obs_names alignment.",
            adata_atac_ga.n_obs, adata_atac.n_obs,
        )
        if not adata_atac_ga.obs_names.is_unique:
            adata_atac_ga.obs_names_make_unique()
        shared_cells = adata_atac_ga.obs_names.intersection(adata_atac.obs_names)
        log.info("Shared cells: %d", len(shared_cells))
        adata_atac_ga = adata_atac_ga[shared_cells].copy()
        for col in adata_atac.obs.columns:
            if col in ("cell_type", "leiden", "condition", "donor_id", "donor", "batch"):
                adata_atac_ga.obs[col] = adata_atac.obs.loc[shared_cells, col].values

    # -----------------------------------------------------------------------
    # Step 3: Load RNA reference
    # -----------------------------------------------------------------------
    log.info("\n=== Step 3: Load RNA reference ===")
    adata_ref, pairing = load_rna_reference(args.rna_atlas, args.donor_pairing)
    donor_to_rna = build_donor_map(adata_ref, pairing)

    if len(donor_to_rna) == 0:
        log.error("No donor-to-RNA mappings found. Cannot proceed.")
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Step 4: Per-donor label transfer
    # -----------------------------------------------------------------------
    log.info("\n=== Step 4: Per-donor label transfer ===")
    transferred_labels, transfer_distances = transfer_labels_per_donor(
        adata_atac_ga,
        adata_ref,
        donor_to_rna,
        n_hvgs=args.n_hvgs,
        n_pcs=args.n_pcs,
        k_neighbors=args.k_neighbors,
    )

    # Free RNA reference
    del adata_ref
    gc.collect()

    # -----------------------------------------------------------------------
    # Step 5: Confidence threshold
    # -----------------------------------------------------------------------
    log.info("\n=== Step 5: Apply confidence threshold ===")
    threshold, low_conf_mask = apply_confidence_threshold(
        transferred_labels, transfer_distances, args.confidence_percentile
    )

    # Final labels: transferred unless low confidence or unassigned
    final_labels = transferred_labels.copy()
    final_labels[low_conf_mask] = "Low_confidence"
    final_labels[transferred_labels.isna()] = "Unassigned"

    # -----------------------------------------------------------------------
    # Step 6: Save results
    # -----------------------------------------------------------------------
    log.info("\n=== Step 6: Save results ===")

    # Update processed ATAC AnnData with transferred labels
    # Labels are indexed by gene_activity obs_names (which may differ from processed).
    # Since both have same cell count and order (from same AnnDataSet), use positional mapping.
    adata_atac.obs["cell_type_gene_activity"] = adata_atac.obs["cell_type"].copy()

    if adata_atac.n_obs == len(final_labels):
        # Positional mapping (same order guaranteed from same AnnDataSet source)
        adata_atac.obs["cell_type_transferred"] = final_labels.values
        adata_atac.obs["cell_type_transfer_distance"] = transfer_distances.values
        adata_atac.obs["cell_type_transfer_raw"] = transferred_labels.values
    else:
        # Fall back to reindex (may produce NaN for unmatched cells)
        adata_atac.obs["cell_type_transferred"] = final_labels.reindex(
            adata_atac.obs_names
        ).values
        adata_atac.obs["cell_type_transfer_distance"] = transfer_distances.reindex(
            adata_atac.obs_names
        ).values
        adata_atac.obs["cell_type_transfer_raw"] = transferred_labels.reindex(
            adata_atac.obs_names
        ).values

    # Save label-transferred AnnData
    out_h5ad = output_dir / "snapatac2_label_transferred.h5ad"
    adata_atac.write_h5ad(str(out_h5ad))
    log.info("Saved: %s", out_h5ad)

    # Per-cell annotation table
    anno_cols = [
        "cell_type_gene_activity",
        "cell_type_transferred",
        "cell_type_transfer_raw",
        "cell_type_transfer_distance",
    ]
    for col in ("donor_id", "donor", "condition", "leiden"):
        if col in adata_atac.obs.columns:
            anno_cols.append(col)
    adata_atac.obs[anno_cols].to_csv(output_dir / "cell_annotations_comparison.csv")

    # -----------------------------------------------------------------------
    # Step 7: Compare old vs new
    # -----------------------------------------------------------------------
    log.info("\n=== Step 7: Compare annotations ===")

    # Exclude low-confidence and unassigned for comparison
    valid_mask = ~adata_atac.obs["cell_type_transferred"].isin(
        ["Low_confidence", "Unassigned"]
    )
    adata_valid = adata_atac[valid_mask].copy()

    comparison = compare_annotations(
        adata_valid,
        old_col="cell_type_gene_activity",
        new_col="cell_type_transferred",
        output_dir=output_dir,
    )

    # -----------------------------------------------------------------------
    # Relabel-fraction transparency (review: KNN transfer can relabel >50% of
    # cells vs their prior gene-activity annotation, with no joint
    # scaling/integration). Report it, do NOT alter the transfer algorithm.
    # Prior label = cell_type_gene_activity; new label = cell_type_transferred.
    # Computed on confidently-transferred cells (Low_confidence/Unassigned
    # already excluded in adata_valid), using harmonized major types so naming
    # conventions (Hepatocyte vs Hepatocytes) do not inflate the fraction.
    # -----------------------------------------------------------------------
    prior_harm = adata_valid.obs["cell_type_gene_activity"].astype(str).map(_harmonize_label)
    new_harm = adata_valid.obs["cell_type_transferred"].astype(str).map(_harmonize_label)
    n_valid = int(len(adata_valid))
    relabel_fraction = float((prior_harm != new_harm).mean()) if n_valid > 0 else float("nan")
    # Low-confidence fraction is over ALL cells (low-confidence is the
    # KNN-distance-flagged minority that never receives a confident label).
    low_conf_fraction = float(
        (adata_atac.obs["cell_type_transferred"] == "Low_confidence").mean()
    )

    comparison["relabel_fraction"] = relabel_fraction
    comparison["n_relabeled_vs_prior"] = (
        int((prior_harm != new_harm).sum()) if n_valid > 0 else 0
    )
    comparison["n_confident_cells"] = n_valid
    comparison["low_confidence_fraction"] = low_conf_fraction
    # Re-emit the summary JSON now that the relabel keys are attached.
    with open(output_dir / "comparison_summary.json", "w") as f:
        json.dump(comparison, f, indent=2)

    log.info(
        "Relabel fraction (harmonized, vs prior gene-activity label): %.3f "
        "(%d / %d confidently-transferred cells)",
        relabel_fraction, comparison["n_relabeled_vs_prior"], n_valid,
    )
    if np.isfinite(relabel_fraction) and relabel_fraction > 0.5:
        log.warning(
            "RELABEL FRACTION > 0.5: %.1f%% of confidently-transferred cells "
            "received a DIFFERENT major cell type than their prior "
            "gene-activity annotation. KNN transfer is per-donor in RNA PCA "
            "space with NO joint scaling/integration of the ATAC modality "
            "(documented limitation); interpret transferred labels with care.",
            100 * relabel_fraction,
        )

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    n_transferred = (~transferred_labels.isna()).sum()
    n_low_conf = low_conf_mask.sum()
    n_unassigned = transferred_labels.isna().sum()

    log.info("\n" + "=" * 60)
    log.info("LABEL TRANSFER SUMMARY")
    log.info("=" * 60)
    log.info("Total cells:             %d", len(adata_atac))
    log.info("Successfully transferred: %d (%.1f%%)", n_transferred, 100 * n_transferred / len(adata_atac))
    log.info("Low confidence (flagged): %d (%.1f%%)", n_low_conf, 100 * n_low_conf / len(adata_atac))
    log.info("Unassigned:              %d (%.1f%%)", n_unassigned, 100 * n_unassigned / len(adata_atac))
    log.info("Harmonized concordance: %.1f%%", 100 * comparison["concordance_harmonized"])
    log.info("Cells changed (harmonized): %d (%.1f%%)", comparison["n_changed_harmonized"], 100 * comparison["frac_changed_harmonized"])
    log.info(
        "\nNew cell type distribution:\n%s",
        adata_atac.obs["cell_type_transferred"].value_counts().to_string(),
    )
    log.info("\nConfidence threshold (distance): %.3f", threshold)
    log.info("Outputs saved to %s", output_dir)

    # Save summary as text
    with open(output_dir / "label_transfer_summary.txt", "w") as f:
        f.write("Label Transfer Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Total cells:              {len(adata_atac)}\n")
        f.write(f"Successfully transferred: {n_transferred}\n")
        f.write(f"Low confidence:           {n_low_conf}\n")
        f.write(f"Unassigned:               {n_unassigned}\n")
        f.write(f"Harmonized concordance:   {comparison['concordance_harmonized']:.3f}\n")
        f.write(f"Raw concordance:          {comparison['concordance_raw']:.3f} (naming mismatch expected)\n")
        f.write(f"Cells changed (harmonized): {comparison['n_changed_harmonized']}\n")
        f.write(f"Relabel fraction (vs prior gene-activity): {comparison['relabel_fraction']:.3f}")
        f.write(" [WARNING >0.5]\n" if (np.isfinite(comparison['relabel_fraction']) and comparison['relabel_fraction'] > 0.5) else "\n")
        f.write(f"Low-confidence fraction:  {comparison['low_confidence_fraction']:.3f}\n")
        f.write(f"Confidence threshold:     {threshold:.3f}\n")
        f.write(f"\nPer-major-type Jaccard (harmonized):\n")
        for ct, jac in sorted(comparison["jaccard_harmonized"].items()):
            f.write(f"  {ct:25s}  {jac:.3f}\n")
        f.write(f"\nNew cell type distribution:\n")
        f.write(adata_atac.obs["cell_type_transferred"].value_counts().to_string())
        f.write("\n")

    log.info("\nDone!")


if __name__ == "__main__":
    main()
