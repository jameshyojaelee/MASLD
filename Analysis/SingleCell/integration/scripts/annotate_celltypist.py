#!/usr/bin/env python3
"""CellTypist annotation of ScaleSC-integrated atlas.

ScaleSC saves adata without X (expression matrix) by design.
This script reconstructs X from the original per-sample h5ad files,
runs CellTypist annotation, validates markers, and exports
proportions + pseudobulk counts.

Requires: rapids_singlecell env with celltypist installed.
"""
import argparse
import logging
import sys
import time
from pathlib import Path

import anndata as ad
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, vstack, csr_matrix

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTEGRATION_DIR = BASE / "Analysis/SingleCell/integration"
RESULTS_DIR = BASE / "Analysis/SingleCell/results_gpu_v2"

CANONICAL_MARKERS = {
    "Hepatocytes": ["ALB", "APOB", "APOC3", "CYP3A4"],
    "Cholangiocytes": ["KRT19", "EPCAM", "KRT7"],
    "Endothelial cells": ["PECAM1", "VWF", "CDH5"],
    "Macrophages": ["CD68", "MARCO", "CD163"],
    "Fibroblasts": ["ACTA2", "COL1A1", "COL3A1"],
    "T cells": ["CD3D", "CD3E", "NKG7"],
}


def reconstruct_expression(scalesc_h5ad: str, input_dir: str) -> ad.AnnData:
    """Reconstruct full adata (with X) from ScaleSC output + per-sample h5ads.

    ScaleSC saves adata without X. We load each per-sample h5ad,
    filter to cells that survived ScaleSC QC, and build the
    expression matrix.
    """
    t0 = time.time()
    log.info(f"Loading ScaleSC output: {scalesc_h5ad}")
    adata_sc = ad.read_h5ad(scalesc_h5ad)
    log.info(f"  {adata_sc.n_obs:,} cells x {adata_sc.n_vars:,} genes")
    log.info(f"  obs columns: {list(adata_sc.obs.columns)}")
    log.info(f"  obsm keys: {list(adata_sc.obsm.keys())}")

    # Get ScaleSC gene list
    sc_genes = set(adata_sc.var_names)
    sc_gene_list = list(adata_sc.var_names)
    log.info(f"  ScaleSC gene set: {len(sc_genes)} genes")

    # Build per-sample barcode sets from ScaleSC output
    # Barcodes are NOT unique across samples, so we need (sample, barcode) pairs
    samples = adata_sc.obs["sample"].unique()
    log.info(f"  {len(samples)} samples in ScaleSC output")

    sample_barcodes = {}
    for sample in samples:
        mask = adata_sc.obs["sample"] == sample
        barcodes = set(adata_sc.obs.index[mask])
        sample_barcodes[sample] = barcodes

    # Reconstruct expression matrix from per-sample h5ads
    input_path = Path(input_dir)
    chunks = []  # list of (sample, adata_subset) in ScaleSC cell order
    cells_found = 0
    cells_missing = 0

    for i, sample in enumerate(samples):
        h5ad_path = input_path / f"{sample}.h5ad"
        if not h5ad_path.exists():
            log.warning(f"  {sample}: h5ad not found at {h5ad_path}")
            cells_missing += len(sample_barcodes[sample])
            continue

        # Load per-sample h5ad
        sdata = ad.read_h5ad(h5ad_path)
        sdata.var_names_make_unique()

        # Filter to barcodes that survived ScaleSC QC
        kept = sample_barcodes[sample]
        mask = sdata.obs.index.isin(kept)
        sdata = sdata[mask]

        if sdata.n_obs == 0:
            log.warning(f"  {sample}: 0 cells matched (expected {len(kept)})")
            cells_missing += len(kept)
            continue

        # Subset to ScaleSC gene set (intersect)
        common_genes = [g for g in sc_gene_list if g in sdata.var_names]
        if len(common_genes) < len(sc_gene_list):
            # Some genes missing in this sample — fill with zeros later
            pass

        sdata = sdata[:, [g for g in sdata.var_names if g in sc_genes]]
        cells_found += sdata.n_obs
        chunks.append((sample, sdata))

        if (i + 1) % 50 == 0 or (i + 1) == len(samples):
            log.info(f"  Loaded {i+1}/{len(samples)} samples ({cells_found:,} cells)")

    log.info(f"  Reconstruction: {cells_found:,} cells found, {cells_missing:,} missing")
    log.info(f"  Loading took {time.time() - t0:.1f}s")

    # Concatenate all chunks
    t1 = time.time()
    log.info("Concatenating expression data...")
    all_adata = [chunk[1] for chunk in chunks]
    adata_concat = ad.concat(all_adata, merge="same")
    log.info(f"  Concatenated: {adata_concat.n_obs:,} cells x {adata_concat.n_vars:,} genes ({time.time() - t1:.1f}s)")

    # Now we need to align the cell order to match ScaleSC output
    # Build a unique cell ID: barcode (obs index) — but barcodes can be duplicated
    # After concat, anndata appends "-0", "-1" etc. to make obs_names unique
    # We need to match by (sample, original_barcode)

    # Build sample column from concat
    concat_samples = adata_concat.obs["sample"].values
    concat_barcodes = adata_concat.obs_names.values

    # ScaleSC obs also has non-unique barcodes. Build composite keys.
    # ScaleSC obs order is the canonical order we need to match
    sc_keys = list(zip(adata_sc.obs["sample"].values, adata_sc.obs_names.values))

    # For the concatenated data, the obs_names may have been suffixed.
    # Let's rebuild from the sample metadata which is still correct
    # We need to match based on sample + original barcode

    # Actually: anndata.concat may rename obs_names to make them unique.
    # Let's use a different approach: iterate through ScaleSC cells and
    # pull from concat by (sample, barcode).

    # Build a lookup: (sample, barcode) -> row index in concat
    log.info("Aligning cell order to ScaleSC output...")
    t2 = time.time()

    # For each sample, get the barcodes in concat order
    concat_lookup = {}
    for idx in range(adata_concat.n_obs):
        sample = concat_samples[idx]
        # The barcode in concat may have been renamed, but the original
        # barcode is from the h5ad. Let's store by position within sample.
        key = (sample, concat_barcodes[idx])
        concat_lookup[key] = idx

    # If barcodes were made unique by concat (appended -0, -1, etc.), we need
    # to handle that. Let's check:
    if adata_concat.obs_names.is_unique:
        # Barcodes were made unique somehow
        # Try matching with the raw barcode (strip suffix)
        pass

    # Alternative approach: since we know the order within each sample
    # is preserved from the h5ad, and ScaleSC preserves order within batches,
    # let's just match by sample and position.
    # Actually, let's do it properly with the sample column:
    # Group concat cells by sample, then for each ScaleSC cell, find it in the right sample group.

    # Simpler: build per-sample barcode→row index maps
    sample_to_rows = {}
    for idx in range(adata_concat.n_obs):
        s = concat_samples[idx]
        if s not in sample_to_rows:
            sample_to_rows[s] = {}
        # Strip any suffix added by concat (e.g., "-0" appended)
        bc = concat_barcodes[idx]
        # anndata.concat with make_names_unique: appends -0, -1, etc.
        # Original barcode already ends with -1 (10x suffix), so concat
        # may append another -N. Let's store the raw bc.
        sample_to_rows[s][bc] = idx

    # Now build reindex array: for each ScaleSC cell, find its row in concat
    sc_obs_names = adata_sc.obs_names.values
    sc_samples = adata_sc.obs["sample"].values

    reindex = np.full(len(sc_obs_names), -1, dtype=int)
    matched = 0
    unmatched_samples = set()

    for i in range(len(sc_obs_names)):
        s = sc_samples[i]
        bc = sc_obs_names[i]
        if s in sample_to_rows and bc in sample_to_rows[s]:
            reindex[i] = sample_to_rows[s][bc]
            matched += 1
        else:
            unmatched_samples.add(s)

    log.info(f"  Matched {matched:,}/{len(sc_obs_names):,} cells ({time.time() - t2:.1f}s)")

    if matched < len(sc_obs_names) * 0.9:
        # Try stripping concat suffixes
        log.warning(f"  Low match rate. Trying to strip concat suffixes...")

        # Rebuild concat lookup without suffix
        sample_to_rows2 = {}
        for idx in range(adata_concat.n_obs):
            s = concat_samples[idx]
            if s not in sample_to_rows2:
                sample_to_rows2[s] = {}
            bc = concat_barcodes[idx]
            # Strip trailing -N that concat added (but keep original -1)
            # Original 10x barcode: ACGT...-1
            # After concat: ACGT...-1-0 or ACGT...-1-1
            parts = bc.rsplit("-", 1)
            if len(parts) == 2 and parts[1].isdigit() and len(parts[1]) <= 3:
                bc_stripped = parts[0]
            else:
                bc_stripped = bc
            sample_to_rows2[s][bc_stripped] = idx

        reindex2 = np.full(len(sc_obs_names), -1, dtype=int)
        matched2 = 0
        for i in range(len(sc_obs_names)):
            s = sc_samples[i]
            bc = sc_obs_names[i]
            if s in sample_to_rows2 and bc in sample_to_rows2[s]:
                reindex2[i] = sample_to_rows2[s][bc]
                matched2 += 1

        if matched2 > matched:
            log.info(f"  After stripping: matched {matched2:,}/{len(sc_obs_names):,}")
            reindex = reindex2
            matched = matched2

    if matched < len(sc_obs_names):
        log.warning(f"  {len(sc_obs_names) - matched:,} cells unmatched — will have zero expression")

    # Reindex the expression matrix
    valid_mask = reindex >= 0
    X_raw = adata_concat.X

    # Build aligned sparse matrix
    from scipy.sparse import lil_matrix
    n_cells = len(sc_obs_names)
    n_genes = adata_concat.n_vars
    log.info(f"Building aligned expression matrix ({n_cells:,} x {n_genes:,})...")
    t3 = time.time()

    # Use fancy indexing on the sparse matrix
    valid_indices = np.where(valid_mask)[0]
    valid_concat_indices = reindex[valid_indices]

    # Create output sparse matrix
    if issparse(X_raw):
        X_aligned = X_raw[valid_concat_indices]
    else:
        X_aligned = csr_matrix(X_raw[valid_concat_indices])

    # If some cells unmatched, insert zero rows
    if matched < n_cells:
        from scipy.sparse import lil_matrix as lil
        X_full = lil((n_cells, n_genes), dtype=X_aligned.dtype)
        X_full[valid_indices] = X_aligned
        X_aligned = csr_matrix(X_full)

    log.info(f"  Aligned matrix: {X_aligned.shape} ({time.time() - t3:.1f}s)")

    # Build the final adata with expression + ScaleSC metadata
    adata_final = ad.AnnData(
        X=X_aligned,
        obs=adata_sc.obs.copy(),
        var=adata_concat.var.copy() if adata_concat.n_vars == len(sc_gene_list) else adata_sc.var.copy(),
        obsm=adata_sc.obsm.copy(),
    )

    # Align gene names if needed
    if adata_final.n_vars != len(sc_gene_list):
        log.warning(f"  Gene count mismatch: concat has {adata_final.n_vars}, ScaleSC has {len(sc_gene_list)}")

    log.info(f"Reconstructed adata: {adata_final.n_obs:,} cells x {adata_final.n_vars:,} genes")
    log.info(f"Total reconstruction time: {time.time() - t0:.1f}s")

    return adata_final


def run_celltypist(adata: ad.AnnData) -> ad.AnnData:
    """Run CellTypist annotation on log-normalized data."""
    # Store raw counts before normalization
    adata.raw = adata.copy()

    # Make obs_names unique (required for CellTypist — barcodes are duplicated across samples)
    log.info("Making obs_names unique for CellTypist...")
    original_obs_names = adata.obs_names.copy()
    adata.obs_names = adata.obs["sample"].astype(str) + "_" + adata.obs_names.astype(str)
    adata.obs_names_make_unique()
    log.info(f"  obs_names now unique: {adata.obs_names.is_unique}")

    # Normalize + log1p for CellTypist
    t0 = time.time()
    log.info("Normalizing for CellTypist (normalize_total + log1p)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    log.info(f"  Normalization: {time.time() - t0:.1f}s")

    # Check value range
    if issparse(adata.X):
        max_val = float(adata.X.data.max()) if adata.X.nnz > 0 else 0.0
    else:
        max_val = float(np.max(adata.X))
    log.info(f"  X max value after log1p: {max_val:.2f}")

    # Load CellTypist model
    t1 = time.time()
    log.info("Loading CellTypist model: Healthy_Human_Liver.pkl")
    model = celltypist.models.Model.load(model="Healthy_Human_Liver.pkl")
    log.info(f"  Model loaded in {time.time() - t1:.1f}s")
    log.info(f"  Cell types: {list(model.cell_types)}")

    # Run CellTypist
    t2 = time.time()
    log.info("Running CellTypist prediction (majority_voting=True)...")
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=True,
    )
    log.info(f"  CellTypist done in {time.time() - t2:.1f}s")

    # Extract predictions directly (avoid to_adata() reindex issues with duplicate obs_names)
    pred_labels = predictions.predicted_labels
    log.info(f"  Prediction columns: {list(pred_labels.columns)}")
    log.info(f"  Prediction shape: {pred_labels.shape}")
    adata.obs["cell_type"] = pred_labels["majority_voting"].values
    adata.obs["cell_type_raw"] = pred_labels["predicted_labels"].values
    adata.obs["cell_type_conf"] = predictions.probability_matrix.max(axis=1).values

    # Summary
    ct_counts = adata.obs["cell_type"].value_counts()
    log.info("Cell type distribution:")
    for ct, n in ct_counts.items():
        pct = 100 * n / adata.n_obs
        log.info(f"  {ct}: {n:,} ({pct:.1f}%)")

    return adata


def validate_markers(adata: ad.AnnData):
    """Log mean expression of canonical markers per cell type."""
    log.info("\n=== Marker Validation ===")
    cell_types = sorted(adata.obs["cell_type"].unique())

    for expected_ct, markers in CANONICAL_MARKERS.items():
        available = [m for m in markers if m in adata.var_names]
        if not available:
            log.warning(f"  {expected_ct}: no markers found")
            continue

        log.info(f"\n  {expected_ct} markers: {available}")
        for ct in cell_types:
            mask = adata.obs["cell_type"] == ct
            if mask.sum() == 0:
                continue
            subset = adata[mask, available]
            if issparse(subset.X):
                mean_expr = np.array(subset.X.mean(axis=0)).flatten()
            else:
                mean_expr = np.mean(subset.X, axis=0)
            vals = ", ".join(f"{m}={v:.2f}" for m, v in zip(available, mean_expr))
            marker_tag = " <-- EXPECTED" if ct == expected_ct else ""
            log.info(f"    {ct}: {vals}{marker_tag}")


def export_proportions(adata: ad.AnnData, output_dir: Path):
    """Export cell-type proportions per sample."""
    output_dir.mkdir(parents=True, exist_ok=True)

    ct_counts = pd.crosstab(adata.obs["sample"], adata.obs["cell_type"])
    ct_props = ct_counts.div(ct_counts.sum(axis=1), axis=0)

    sample_dataset = adata.obs.groupby("sample")["dataset"].first()
    ct_props.insert(0, "dataset", sample_dataset)
    ct_props.index.name = "sample"

    out_path = output_dir / "cell_type_proportions.csv"
    ct_props.to_csv(out_path)
    log.info(f"Proportions: {out_path} ({ct_props.shape[0]} samples x {ct_props.shape[1]-1} types)")

    numeric_cols = ct_props.select_dtypes(include=[np.number])
    mean_props = numeric_cols.mean().sort_values(ascending=False)
    log.info("Mean proportions:")
    for ct, prop in mean_props.items():
        log.info(f"  {ct}: {prop*100:.1f}%")


def export_pseudobulk(adata: ad.AnnData, output_dir: Path, min_cells: int = 10):
    """Export pseudobulk count matrices per cell type from raw counts."""
    pb_dir = output_dir / "pseudobulk"
    pb_dir.mkdir(parents=True, exist_ok=True)

    # Use raw counts for pseudobulk
    if adata.raw is not None:
        log.info("Using adata.raw for pseudobulk counts")
        raw_X = adata.raw.X
        raw_var_names = adata.raw.var_names
    else:
        log.warning("No adata.raw; using expm1(X) as approximate raw counts")
        if issparse(adata.X):
            raw_X = adata.X.expm1()
        else:
            raw_X = np.expm1(adata.X)
        raw_var_names = adata.var_names

    cell_types = sorted(adata.obs["cell_type"].unique())
    for ct in cell_types:
        ct_mask = (adata.obs["cell_type"] == ct).values
        ct_obs = adata.obs[ct_mask]

        sample_counts = ct_obs["sample"].value_counts()
        valid_samples = sample_counts[sample_counts >= min_cells].index

        if len(valid_samples) < 3:
            log.info(f"  {ct}: skipping ({len(valid_samples)} samples with >= {min_cells} cells)")
            continue

        # Sum counts per sample
        pb_dict = {}
        for sample in valid_samples:
            sample_mask = ct_mask & (adata.obs["sample"] == sample).values
            subset_X = raw_X[sample_mask]
            if issparse(subset_X):
                pb_dict[sample] = np.array(subset_X.sum(axis=0)).flatten()
            else:
                pb_dict[sample] = np.sum(subset_X, axis=0)

        pb_df = pd.DataFrame(pb_dict, index=raw_var_names)
        pb_df = pb_df.round(0).astype(int)

        ct_safe = ct.replace(" ", "_").replace("+", "+").replace("/", "_")
        out_path = pb_dir / f"{ct_safe}_pseudobulk.csv"
        pb_df.to_csv(out_path)
        log.info(f"  {ct}: {pb_df.shape[0]:,} genes x {pb_df.shape[1]} samples -> {out_path.name}")


def main():
    parser = argparse.ArgumentParser(description="CellTypist annotation + pseudobulk export")
    parser.add_argument("--scalesc-h5ad", type=str,
                        default=str(INTEGRATION_DIR / "output/human/scalesc_human_integrated.h5ad"))
    parser.add_argument("--input-dir", type=str,
                        default=str(INTEGRATION_DIR / "input_h5ad/human"))
    parser.add_argument("--output-dir", type=str,
                        default=str(RESULTS_DIR))
    parser.add_argument("--skip-pseudobulk", action="store_true")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Reconstruct expression matrix
    adata = reconstruct_expression(args.scalesc_h5ad, args.input_dir)

    # Step 2: CellTypist annotation (normalizes internally)
    adata = run_celltypist(adata)

    # Step 3: Marker validation
    validate_markers(adata)

    # Step 4: Export proportions
    export_proportions(adata, output_dir)

    # Step 5: Export pseudobulk
    if not args.skip_pseudobulk:
        export_pseudobulk(adata, output_dir)

    # Step 6: Save annotated h5ad
    out_h5ad = INTEGRATION_DIR / "output/human/scalesc_human_annotated_celltypist.h5ad"
    log.info(f"Saving annotated h5ad: {out_h5ad}")
    adata.write_h5ad(out_h5ad)

    # Step 7: Create symlink for figure scripts
    symlink_path = output_dir / "integrated_atlas.h5ad"
    if symlink_path.exists() or symlink_path.is_symlink():
        symlink_path.unlink()
    symlink_path.symlink_to(out_h5ad)
    log.info(f"Symlink: {symlink_path} -> {out_h5ad}")

    log.info("Done!")


if __name__ == "__main__":
    main()
