#!/usr/bin/env python3
"""
scVI integration of GSE136103 + GSE192740 human scRNA-seq data.

Inputs:
  - GSE136103 cellranger outputs: data/GSE136103/cellranger/SRR*/outs/filtered_feature_bc_matrix.h5
  - GSE192740 cellranger outputs: Liver_Atlas/cellranger/human/SRR*/outs/filtered_feature_bc_matrix.h5

Outputs:
  - Analysis/SingleCell/results/integrated_atlas.h5ad
  - Analysis/SingleCell/results/cell_type_proportions.csv
  - Analysis/SingleCell/results/pseudobulk/  (per-cell-type matrices)
  - Analysis/SingleCell/results/scvi_model/
"""

import argparse
import csv
import glob
import logging
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import celltypist
from celltypist import models as ct_models

# GPU acceleration support
sys.path.insert(0, str(Path(__file__).parent))
from gpu_utils import init_gpu, get_processor, to_gpu, from_gpu

# --- Seed pinning (added 2026-04-22 per T0.8) -------------------------------
# NOTE: placed after gpu_utils import (which handles PyTorch CUDA init order)
# but before any scvi model construction. init_gpu() is called later in main().
import random
os.environ['PYTHONHASHSEED'] = '42'
random.seed(42)
np.random.seed(42)
import torch
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
try:
    scvi.settings.seed = 42
except Exception:
    pass
# ---------------------------------------------------------------------------

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")


# --- Preparation method metadata ---
# GSE192740 SRR → preparation_method mapping
# Built from SraRunTable.csv + GSE192740_sampleInfo_scRNAseq.tsv
# Digestion method + sample name determine preparation type:
#   Nuclei → "nuclei" (snRNA-seq, captures hepatocytes at ~60%)
#   Ex Vivo + CD45+ → "CD45_positive" (immune-enriched, 0% hepatocytes)
#   Ex Vivo + CD45- → "CD45_negative" (non-immune, may include hepatocytes)
#   Ex Vivo + DCs/monos+macs/Hepatocytes/CD14+CD16+ → "sorted_other"
_GSE192740_PREP_MAP = {
    "SRR17375011": "nuclei",          # ABU8 — H33 nuclei
    "SRR17375012": "CD45_positive",   # CISE06 — H02 CD45+
    "SRR17375013": "CD45_positive",   # CISE06 — H02 CD45+
    "SRR17375014": "CD45_positive",   # CISE07 — H02 CD45+
    "SRR17375015": "CD45_positive",   # CISE07 — H02 CD45+
    "SRR17375016": "CD45_positive",   # CISE08 — H02 CD45+
    "SRR17375017": "CD45_positive",   # CISE08 — H02 CD45+
    "SRR17375018": "CD45_positive",   # CISE09 — H02 CD45+
    "SRR17375019": "CD45_positive",   # CISE09 — H02 CD45+
    "SRR17375020": "CD45_negative",   # CS31 — H02 CD45-
    "SRR17375021": "sorted_other",    # CS32 — H02 DCs
    "SRR17375022": "sorted_other",    # CS33 — H02 DCs
    "SRR17375023": "sorted_other",    # CS34 — H02 monos+macs
    "SRR17375024": "CD45_positive",   # CS37 — H04 CD45+
    "SRR17375025": "CD45_positive",   # CS38 — H04 CD45+
    "SRR17375026": "CD45_positive",   # CS41 — H06 CD45+
    "SRR17375027": "CD45_positive",   # CS42 — H06 CD45+
    "SRR17375028": "CD45_positive",   # CS43 — H06 CD45+
    "SRR17375029": "CD45_positive",   # CS43 — H06 CD45+
    "SRR17375030": "CD45_positive",   # CS44 — H06 CD45+
    "SRR17375031": "CD45_positive",   # CS44 — H06 CD45+
    "SRR17375032": "CD45_positive",   # CS46 — H07 CD45+
    "SRR17375033": "CD45_positive",   # CS71 — H10 CD45+
    "SRR17375034": "CD45_positive",   # CS73 — H11 CD45+
    "SRR17375035": "CD45_positive",   # CS73 — H11 CD45+
    "SRR17375036": "CD45_positive",   # CS81 — H13 CD45+
    "SRR17375037": "CD45_positive",   # CS81 — H13 CD45+
    "SRR17375038": "CD45_positive",   # CS83 — H14 CD45+
    "SRR17375039": "CD45_positive",   # CS85 — H16 CD45+
    "SRR17375040": "CD45_positive",   # CS87 — H16 CD45+
    "SRR17375041": "CD45_positive",   # CS87 — H16 CD45+
    "SRR17375042": "CD45_positive",   # CS101 — H18 CD45+
    "SRR17375043": "CD45_positive",   # CS101 — H18 CD45+
    "SRR17375044": "CD45_positive",   # CS108 — H21 CD45+
    "SRR17375045": "sorted_other",    # CS109 — H21 CD45+ CD14+CD16+
    "SRR17375046": "CD45_positive",   # CS110 — H22 CD45+
    "SRR17375047": "CD45_negative",   # CS111 — H22 CD45-
    "SRR17375048": "CD45_positive",   # CS112 — H23 CD45+
    "SRR17375049": "CD45_positive",   # CS126 — H25 CD45+
    "SRR17375050": "CD45_negative",   # CS127 — H25 CD45-
    "SRR17375051": "nuclei",          # CS161 — H30 nuclei
    "SRR17375052": "nuclei",          # CS162 — H30 nuclei
    "SRR17375053": "nuclei",          # CS164 — H32 nuclei
    "SRR17375054": "nuclei",          # CS166 — H37 nuclei
    "SRR17375055": "nuclei",          # CS167 — H37 nuclei
    "SRR17375056": "nuclei",          # CS169 — H37 nuclei
    "SRR17375057": "nuclei",          # CS170 — H38 nuclei
    "SRR17375058": "nuclei",          # CS171 — H38 nuclei (Ex Vivo in metadata but nuclei isolation)
}


def assign_preparation_method(adata: ad.AnnData) -> ad.AnnData:
    """Add preparation_method column to adata.obs based on dataset and sample ID.

    Categories:
        - enzymatic_dissociation: GSE136103 (Ramachandran) — lyses hepatocytes, <1% hepatocytes
        - nuclei: snRNA-seq samples — captures hepatocytes at ~60%
        - CD45_positive: CD45+ sorted — immune-enriched, 0% hepatocytes by design
        - CD45_negative: CD45- sorted — non-immune fraction
        - sorted_other: DCs, monos+macs, CD14+CD16+ sorted fractions
    """
    prep = pd.Series("unknown", index=adata.obs.index)

    # GSE136103: all samples are enzymatic dissociation
    gse136_mask = adata.obs["dataset"] == "GSE136103"
    prep[gse136_mask] = "enzymatic_dissociation"

    # GSE192740: map by SRR ID
    gse192_mask = adata.obs["dataset"] == "GSE192740"
    if gse192_mask.any():
        for srr, method in _GSE192740_PREP_MAP.items():
            srr_mask = gse192_mask & (adata.obs["sample"] == srr)
            prep[srr_mask] = method

    adata.obs["preparation_method"] = prep.astype("category")

    # Report
    print("  Preparation method distribution:")
    for method, count in adata.obs["preparation_method"].value_counts().items():
        pct = 100 * count / adata.n_obs
        print(f"    {method}: {count} ({pct:.1f}%)")

    return adata


def load_h5_matrices(pattern: str, dataset_id: str) -> list[ad.AnnData]:
    """Load all filtered H5 matrices matching a glob pattern."""
    paths = sorted(glob.glob(pattern))
    adatas = []
    for p in paths:
        try:
            adata = sc.read_10x_h5(p)
            adata.var_names_make_unique()
            srr = Path(p).parent.parent.name  # .../SRR.../outs/filtered...
            adata.obs["sample"] = srr
            adata.obs["dataset"] = dataset_id
            adatas.append(adata)
            print(f"  Loaded {srr}: {adata.n_obs} cells, {adata.n_vars} genes")
        except Exception as e:
            print(f"  WARNING: Failed to load {p}: {e}")
    return adatas


def run_scrublet(adata: ad.AnnData) -> ad.AnnData:
    """Run doublet detection per sample using scanpy's scrublet wrapper.

    Uses sc.pp.scrublet which avoids Annoy's SIMD-incompatible KNN
    and instead uses scanpy's own neighbor graph implementation.
    """
    samples = adata.obs["sample"].unique()
    doublet_scores = np.zeros(adata.n_obs)
    predicted_doublets = np.zeros(adata.n_obs, dtype=bool)

    for sample in samples:
        mask = adata.obs["sample"] == sample
        subset = adata[mask].copy()
        if subset.n_obs < 50:
            continue
        try:
            sc.pp.scrublet(subset, batch_key=None)
            doublet_scores[mask] = subset.obs["doublet_score"].values
            predicted_doublets[mask] = subset.obs["predicted_doublet"].values
        except Exception as e:
            print(f"  Scrublet failed for {sample}: {e}")

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets
    return adata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
    parser.add_argument("--output-dir", default="Analysis/SingleCell/results")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--max-pct-mito", type=float, default=20.0)
    parser.add_argument("--n-top-genes", type=int, default=4000)
    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--n-hidden", type=int, default=128)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--max-epochs", type=int, default=400)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--import-data", type=str, default=None, help="Path to pre-processed h5ad file to skip to scVI directly")
    parser.add_argument("--export-data", type=str, default=None, help="Path to save pre-processed h5ad file before scVI")
    parser.add_argument("--resume-model", type=str, default=None,
                        help="Path to a saved scVI model directory (from CP1). "
                             "Skips training and jumps to latent space + annotation.")
    parser.add_argument("--use-gpu", action="store_true",
                        help="Use rapids_singlecell for GPU-accelerated preprocessing and clustering")
    args = parser.parse_args()

    # --- GPU setup ---
    use_gpu = False
    if args.use_gpu:
        use_gpu = init_gpu()
        if use_gpu:
            print("GPU mode: ON (rapids_singlecell)")
        else:
            print("GPU mode: FAILED — falling back to CPU (scanpy)")
    pp, tl = get_processor(use_gpu)

    root = Path(args.project_root)
    out_dir = root / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "pseudobulk").mkdir(exist_ok=True)
    (out_dir / "scvi_model").mkdir(exist_ok=True)

    if args.import_data:
        print(f"Importing pre-processed data from {args.import_data}...")
        adata = sc.read_h5ad(args.import_data)
        print(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    else:
        # --- Load data ---
        print("Loading GSE136103 (Ramachandran)...")
        adatas_136 = load_h5_matrices(
            str(root / "data/GSE136103/cellranger/*/outs/filtered_feature_bc_matrix.h5"),
            "GSE136103"
        )
        print(f"  Total: {len(adatas_136)} samples")

        print("Loading GSE192740 human (Liver Cell Atlas)...")
        adatas_192 = load_h5_matrices(
            str(root / "Liver_Atlas/cellranger/human/*/outs/filtered_feature_bc_matrix.h5"),
            "GSE192740"
        )
        print(f"  Total: {len(adatas_192)} samples")

        all_adatas = adatas_136 + adatas_192
        print(f"\nConcatenating {len(all_adatas)} samples...")
        adata = ad.concat(all_adatas, join="outer")
        adata.obs_names_make_unique()
        print(f"  Combined: {adata.n_obs} cells, {adata.n_vars} genes")

        # --- Preparation method metadata ---
        print("\nAssigning preparation method metadata...")
        adata = assign_preparation_method(adata)

        # --- QC ---
        print("\nQC filtering...")
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

        # Transfer to GPU for accelerated QC + filtering
        if use_gpu:
            print("  Transferring to GPU...")
            to_gpu(adata)

        pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        n_before = adata.n_obs
        pp.filter_cells(adata, min_genes=args.min_genes)
        pp.filter_genes(adata, min_cells=args.min_cells)

        # Transfer back to CPU for subsetting (obs column access)
        if use_gpu:
            from_gpu(adata)
        adata = adata[adata.obs.pct_counts_mt < args.max_pct_mito, :].copy()
        print(f"  Filtered: {n_before} -> {adata.n_obs} cells")

        # --- Doublet removal (always CPU — scrublet not in rapids) ---
        print("\nRunning scrublet doublet detection...")
        adata = run_scrublet(adata)
        n_doublets = adata.obs["predicted_doublet"].sum()
        print(f"  Detected {n_doublets} doublets ({100*n_doublets/adata.n_obs:.1f}%)")
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        print(f"  After removal: {adata.n_obs} cells")

        # --- Store raw counts for scVI ---
        adata.layers["counts"] = adata.X.copy()

        # --- HVG selection ---
        print("\nSelecting highly variable genes...")
        if use_gpu:
            to_gpu(adata)

        pp.normalize_total(adata, target_sum=1e4)
        pp.log1p(adata)
        pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes,
                                 flavor="seurat_v3", batch_key="dataset",
                                 layer="counts")

        if use_gpu:
            from_gpu(adata)
        print(f"  Selected {adata.var.highly_variable.sum()} HVGs")

        if args.export_data:
            print(f"\nExporting pre-processed data to {args.export_data}...")
            adata.write_h5ad(args.export_data)
            print("Export complete. Exiting.")
            return

    # --- scVI ---
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")

    if args.resume_model:
        print(f"\n[RESUME] Loading trained scVI model from {args.resume_model}...")
        model = scvi.model.SCVI.load(args.resume_model, adata=adata)
        print("  Model loaded successfully. Skipping training.")
    else:
        print("\nSetting up scVI model...")
        model = scvi.model.SCVI(adata, n_latent=args.n_latent, n_hidden=args.n_hidden,
                                 n_layers=args.n_layers)

        print(f"Training scVI (max_epochs={args.max_epochs})...")
        model.train(max_epochs=args.max_epochs, early_stopping=True,
                    early_stopping_patience=20)

        # ===== CHECKPOINT 1: Save model immediately after training =====
        print("\n[CHECKPOINT 1] Saving trained scVI model...")
        model.save(str(out_dir / "scvi_model"), overwrite=True)
        print("  Model saved to:", out_dir / "scvi_model")

    # --- Latent space + clustering ---
    print("\nComputing latent space...")
    adata.obsm["X_scVI"] = model.get_latent_representation()

    if use_gpu:
        print("  GPU-accelerated neighbors + UMAP + leiden...")
        to_gpu(adata)

    pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
    tl.umap(adata)
    tl.leiden(adata, resolution=args.leiden_resolution, key_added="leiden")

    if use_gpu:
        from_gpu(adata)
    print(f"  Found {adata.obs['leiden'].nunique()} clusters")

    # ===== CHECKPOINT 2: Save h5ad with latent space + UMAP + clusters =====
    print("\n[CHECKPOINT 2] Saving h5ad with latent space and clusters...")
    adata.write_h5ad(str(out_dir / "integrated_atlas.h5ad"))
    print("  Saved:", out_dir / "integrated_atlas.h5ad")

    # --- Cell type annotation (CellTypist) ---
    print("\nAnnotating cell types with CellTypist...")
    # Use the best available liver model
    available_models = ct_models.models_description()
    liver_model_name = None
    for candidate in ["Healthy_Human_Liver.pkl", "Immune_All_Low.pkl"]:
        if candidate in available_models["model"].values:
            liver_model_name = candidate
            break
    if liver_model_name is None:
        liver_model_name = "Immune_All_Low.pkl"

    print(f"  Using model: {liver_model_name}")
    ct_model = ct_models.Model.load(model=liver_model_name)

    # CellTypist requires log-normalized data in .X (which we already have)
    predictions = celltypist.annotate(adata, model=ct_model, majority_voting=True)
    adata.obs["cell_type_celltypist"] = predictions.predicted_labels.get("majority_voting", predictions.predicted_labels.iloc[:, 0])
    
    if "conf_score" in predictions.predicted_labels.columns:
        adata.obs["cell_type_confidence"] = predictions.predicted_labels["conf_score"]
    else:
        adata.obs["cell_type_confidence"] = predictions.probability_matrix.max(axis=1)

    # Also keep the original marker-based scoring as validation
    markers = {
        "Hepatocyte": ["ALB", "APOB", "SERPINA1", "HP", "TF"],
        "Cholangiocyte": ["KRT19", "KRT7", "EPCAM", "SOX9"],
        "Endothelial": ["PECAM1", "VWF", "CDH5", "ERG"],
        "Kupffer_cell": ["CD68", "MARCO", "VSIG4", "CD163"],
        "Macrophage": ["CD14", "LYZ", "S100A8", "S100A9"],
        "HSC": ["ACTA2", "COL1A1", "PDGFRB", "DES"],
        "T_cell": ["CD3D", "CD3E", "CD4", "CD8A"],
        "NK_cell": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
        "B_cell": ["CD79A", "MS4A1", "CD19"],
        "Plasma_cell": ["JCHAIN", "MZB1", "IGHA1"],
    }
    for ct, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            sc.tl.score_genes(adata, present, score_name=f"score_{ct}")

    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if score_cols:
        score_df = adata.obs[score_cols]
        adata.obs["cell_type_markers"] = score_df.idxmax(axis=1).str.replace("score_", "")

    # Use CellTypist as the primary annotation
    adata.obs["cell_type"] = adata.obs["cell_type_celltypist"]

    ct_counts = adata.obs["cell_type"].value_counts()
    print("  Cell type distribution (CellTypist):")
    for ct, n in ct_counts.items():
        print(f"    {ct}: {n} ({100*n/adata.n_obs:.1f}%)")

    # Report CellTypist confidence
    low_conf = (adata.obs["cell_type_confidence"] < 0.5).sum()
    print(f"  Low-confidence cells (<0.5): {low_conf} ({100*low_conf/adata.n_obs:.1f}%)")

    # ===== CHECKPOINT 3: Save h5ad with full cell type annotations =====
    print("\n[CHECKPOINT 3] Saving h5ad with cell type annotations...")
    adata.write_h5ad(str(out_dir / "integrated_atlas.h5ad"))
    print("  Saved:", out_dir / "integrated_atlas.h5ad")

    # --- Cell type proportions ---
    print("\nComputing cell type proportions per sample...")
    group_cols = ["sample", "dataset"]
    if "preparation_method" in adata.obs.columns:
        group_cols.append("preparation_method")
    proportions = (adata.obs.groupby(group_cols)["cell_type"]
                   .value_counts(normalize=True)
                   .unstack(fill_value=0))
    proportions.to_csv(out_dir / "cell_type_proportions.csv")

    # Also output proportions grouped by preparation method for diagnosis
    if "preparation_method" in adata.obs.columns:
        prep_proportions = (adata.obs.groupby("preparation_method")["cell_type"]
                            .value_counts(normalize=True)
                            .unstack(fill_value=0))
        prep_proportions.to_csv(out_dir / "cell_type_proportions_by_prep_method.csv")
        print("  Saved proportions by preparation method")

    # --- Pseudobulk (wrapped in try/except to protect checkpoint) ---
    print("\nGenerating pseudobulk matrices...")
    try:
        adata_raw = adata.copy()
        adata_raw.X = adata_raw.layers["counts"]

        for ct in adata.obs["cell_type"].unique():
            ct_mask = adata.obs["cell_type"] == ct
            ct_adata = adata_raw[ct_mask]
            if ct_adata.n_obs < 10:
                continue
            # Sanitize cell type name for filesystem (CellTypist uses / in names like "NK/NKT")
            ct_safe = ct.replace("/", "_").replace(" ", "_")
            # Sum counts per sample
            pb = pd.DataFrame(index=ct_adata.var_names)
            for sample in ct_adata.obs["sample"].unique():
                s_mask = ct_adata.obs["sample"] == sample
                pb[sample] = np.asarray(ct_adata[s_mask].X.sum(axis=0)).flatten()
            pb.to_csv(out_dir / "pseudobulk" / f"{ct_safe}_pseudobulk.csv")
            print(f"  {ct}: {pb.shape[1]} samples")
    except Exception as e:
        print(f"  WARNING: Pseudobulk generation failed: {e}")
        print("  Continuing — model + h5ad are already saved.")

    # --- Diagnostic plots (wrapped in try/except) ---
    print("\nGenerating diagnostic plots...")
    try:
        fig_dir = out_dir / "figures"
        fig_dir.mkdir(exist_ok=True)

        sc.settings.figdir = str(fig_dir)
        sc.set_figure_params(dpi=150, frameon=False)

        sc.pl.umap(adata, color="dataset", save="_by_dataset.png", show=False)
        sc.pl.umap(adata, color="cell_type", save="_by_celltype.png", show=False)
        sc.pl.umap(adata, color="cell_type_confidence", save="_confidence.png", show=False)

        # Marker dotplot
        marker_genes = ["ALB", "APOB", "KRT19", "EPCAM", "PECAM1", "VWF",
                        "CD68", "MARCO", "CD14", "LYZ", "ACTA2", "COL1A1",
                        "CD3D", "CD3E", "NKG7", "GNLY", "CD79A", "MS4A1",
                        "JCHAIN", "MZB1"]
        marker_genes = [g for g in marker_genes if g in adata.var_names]
        if marker_genes:
            sc.pl.dotplot(adata, marker_genes, groupby="cell_type",
                          save="_markers.png", show=False)
    except Exception as e:
        print(f"  WARNING: Plot generation failed: {e}")
        print("  Continuing — model + h5ad are already saved.")

    # ===== FINAL SAVE (redundant but guarantees everything is on disk) =====
    print("\n[FINAL] Saving final outputs...")
    model.save(str(out_dir / "scvi_model"), overwrite=True)
    adata.write_h5ad(str(out_dir / "integrated_atlas.h5ad"))

    # Summary
    print(f"\n=== INTEGRATION COMPLETE ===")
    print(f"Total cells: {adata.n_obs}")
    print(f"Total genes: {adata.n_vars}")
    print(f"Datasets: {adata.obs['dataset'].value_counts().to_dict()}")
    print(f"Cell types: {adata.obs['cell_type'].nunique()}")
    print(f"Outputs saved to: {out_dir}")


if __name__ == "__main__":
    main()
