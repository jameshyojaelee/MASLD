#!/usr/bin/env python3
"""
33_vu_cell2location_deconvolution.py — Cell2location deconvolution of Vu et al.
spatial data using the pre-trained MASLD scRNA reference signatures.

Tests whether F3a / F3b spatial bimodality reflects within-hepatocyte
sub-state biology vs hepatocyte/stellate compositional artifact.

Inputs:
  - Vu spatial AnnData: Analysis/Spatial/results/preprocessed/merged_spatial_vu_raw.h5ad
  - Reference signatures (cell2location, MASLD scRNA atlas):
    Analysis/Spatial/results/cell2location/reference_model/inf_aver.csv
  - Cell types include: Hepatocytes, Cholangiocytes, Endothelial cells, Fibroblasts
    (≈ stellate), Macrophages, T cells, B cells, Plasma cells, Neutrophils, +others

Outputs:
  - Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad
  - RNA-seq/results/granular_staging/vu_deconvolved_substate.csv
    columns: spot_id, sample_id, individual, x, y, c2l_<celltype>..., li_f3a,
             li_f3b, li_f3a_z, li_f3b_z, mixing_axis, hep_dominant_80

SLURM: --partition=gpu --gres=gpu:1 --qos=interactive --mem=64G --cpus=8 --time=24:00:00
       env: spatial
"""

import pathlib
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SPATIAL_DIR = PROJECT_ROOT / "Analysis/Spatial"
GRAN_OUT = PROJECT_ROOT / "RNA-seq/results/granular_staging"
GRAN_OUT.mkdir(parents=True, exist_ok=True)

VU_RAW_H5AD = SPATIAL_DIR / "results/preprocessed/merged_spatial_vu_raw.h5ad"
VU_NORM_H5AD = SPATIAL_DIR / "results/preprocessed/merged_spatial_vu.h5ad"
INF_AVER_CSV = SPATIAL_DIR / "results/cell2location/reference_model/inf_aver.csv"
OUT_MODEL_DIR = SPATIAL_DIR / "results/cell2location/spatial_model_vu"
OUT_MODEL_DIR.mkdir(parents=True, exist_ok=True)
OUT_H5AD = OUT_MODEL_DIR / "spatial_deconvolved_vu.h5ad"
OUT_DECONV_CSV = GRAN_OUT / "vu_deconvolved_substate.csv"

# F3a / F3b signatures (Li et al., re-used from script 32)
LI_F3A = [
    "HSPA5", "DDIT3", "ATF4", "ATF6", "EIF2AK3", "ERN1", "XBP1", "DNAJB9",
    "SREBF1", "SREBF2", "HMGCR", "HMGCS1", "SQLE", "DHCR7", "DHCR24",
    "MVD", "MVK", "FASN", "ACACA", "SCD", "INSIG1",
    "FABP1", "ACOX1", "CPT1A", "PLIN2", "DGAT1",
]
LI_F3B = [
    "IGFBP7", "BGN", "COL1A2", "COL3A1", "TIMP1",
    "COL1A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "FN1", "VCAN", "DCN", "LUM", "LOX", "LOXL1", "LOXL2",
    "MMP2", "MMP9", "TIMP2", "ACTA2", "TAGLN",
    "CDKN1A", "CDKN2A", "GLB1", "SERPINE1",
]


def init_gpu():
    import torch
    if torch.cuda.is_available():
        print(f"  GPU initialized: {torch.cuda.get_device_name(0)}")
    else:
        print("  WARNING: No GPU — cell2location will be very slow on CPU")


def run_cell2location(adata_vis, inf_aver, n_cells_per_location=15,
                      detection_alpha=20, max_epochs=15000,
                      num_posterior_samples=1000):
    """Run cell2location deconvolution on Vu spots.

    Uses pre-trained reference signatures (inf_aver) from the MASLD scRNA atlas.
    """
    import cell2location

    # Intersect genes
    shared_genes = adata_vis.var_names.intersection(inf_aver.index)
    print(f"  Shared genes (Vu × c2l reference): {len(shared_genes)}")
    if len(shared_genes) < 1000:
        raise RuntimeError(
            f"Too few shared genes ({len(shared_genes)}). "
            "Vu var_names probably need symbol mapping."
        )
    adata_vis = adata_vis[:, shared_genes].copy()
    inf_aver = inf_aver.loc[shared_genes, :]

    # Filter very-low-expression genes
    sc.pp.filter_genes(adata_vis, min_cells=5)
    inf_aver = inf_aver.loc[inf_aver.index.isin(adata_vis.var_names), :]
    print(f"  After min_cells=5 filter: {adata_vis.n_vars} genes")

    # Setup cell2location with sample_id as batch
    cell2location.models.Cell2location.setup_anndata(
        adata_vis, batch_key="sample_id"
    )
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=n_cells_per_location,
        detection_alpha=detection_alpha,
    )
    print(f"  Model parameters: "
          f"{sum(p.numel() for p in mod.module.parameters()):,}")

    # Train
    print(f"  Training (max_epochs={max_epochs})...")
    train_kwargs = dict(max_epochs=max_epochs, batch_size=None, train_size=1.0)
    try:
        mod.train(accelerator="gpu", **train_kwargs)
    except TypeError:
        mod.train(use_gpu=True, **train_kwargs)

    # Export posterior
    print(f"  Exporting posterior ({num_posterior_samples} samples)...")
    sample_kwargs = {"num_samples": num_posterior_samples, "batch_size": 2500}
    try:
        adata_vis = mod.export_posterior(adata_vis, sample_kwargs=sample_kwargs)
    except TypeError:
        sample_kwargs["use_gpu"] = True
        adata_vis = mod.export_posterior(adata_vis, sample_kwargs=sample_kwargs)

    # Save model
    mod.save(str(OUT_MODEL_DIR), overwrite=True)
    return adata_vis


def fallback_signature_scoring(adata_norm):
    """Signature-based pseudo-deconvolution as fallback if cell2location fails.

    NOT a deconvolution — emits per-spot z-scores for each cell type
    using established marker panels. Per-spot abundances are derived as
    softmax over the cell-type score vector.
    """
    print("  FALLBACK: signature-based per-spot cell-type scoring (no c2l)")
    panels = {
        "Hepatocytes": ["ALB", "APOB", "CYP3A4", "CYP2E1", "TTR", "TF",
                        "HNF4A", "FABP1", "G6PC", "PCK1"],
        "Cholangiocytes": ["KRT19", "KRT7", "EPCAM", "SOX9", "CFTR"],
        "Endothelial cells": ["PECAM1", "VWF", "CDH5", "FLT1", "ENG", "STAB2"],
        "Fibroblasts": ["ACTA2", "COL1A1", "COL1A2", "COL3A1", "PDGFRB",
                        "DCN", "LUM", "LRAT"],
        "Macrophages": ["CD68", "CD163", "MARCO", "VSIG4", "CD14"],
        "T cells": ["CD3D", "CD3E", "CD3G", "CD8A", "TRAC", "CD4"],
        "Plasma cells": ["JCHAIN", "MZB1", "IGHA1", "IGHG1"],
        "B cells": ["MS4A1", "CD79A", "CD19"],
        "Neutrophils": ["FCGR3B", "CXCR2", "S100A8", "S100A9"],
    }
    score_cols = []
    for ct, genes in panels.items():
        present = [g for g in genes if g in adata_norm.var_names]
        if not present:
            print(f"    {ct}: 0 markers present → score 0")
            adata_norm.obs[f"sig_{ct}"] = 0.0
            score_cols.append(f"sig_{ct}")
            continue
        sc.tl.score_genes(adata_norm, gene_list=present, score_name=f"sig_{ct}",
                          random_state=42)
        score_cols.append(f"sig_{ct}")

    # Convert scores → softmax abundance estimates per spot
    score_mat = adata_norm.obs[score_cols].values
    exp_mat = np.exp(score_mat - score_mat.max(axis=1, keepdims=True))
    abund = exp_mat / exp_mat.sum(axis=1, keepdims=True)
    for i, ct_col in enumerate(score_cols):
        ct = ct_col.replace("sig_", "")
        adata_norm.obs[f"c2l_{ct}"] = abund[:, i]
    return adata_norm, [c.replace("sig_", "c2l_") for c in score_cols]


def add_substate_scores(adata_norm):
    """Score each spot for Li F3a / F3b signatures + mixing axis."""
    f3a_present = [g for g in LI_F3A if g in adata_norm.var_names]
    f3b_present = [g for g in LI_F3B if g in adata_norm.var_names]
    print(f"  F3a genes present: {len(f3a_present)}/{len(LI_F3A)}")
    print(f"  F3b genes present: {len(f3b_present)}/{len(LI_F3B)}")

    sc.tl.score_genes(adata_norm, gene_list=f3a_present, score_name="li_f3a",
                      random_state=42)
    sc.tl.score_genes(adata_norm, gene_list=f3b_present, score_name="li_f3b",
                      random_state=42)
    adata_norm.obs["li_f3a_z"] = adata_norm.obs.groupby(
        "individual", observed=False
    )["li_f3a"].transform(lambda x: (x - x.mean()) / (x.std() + 1e-6))
    adata_norm.obs["li_f3b_z"] = adata_norm.obs.groupby(
        "individual", observed=False
    )["li_f3b"].transform(lambda x: (x - x.mean()) / (x.std() + 1e-6))
    adata_norm.obs["mixing_axis"] = (
        adata_norm.obs["li_f3b_z"] - adata_norm.obs["li_f3a_z"]
    )
    return adata_norm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--use-fallback", action="store_true",
                        help="Skip cell2location, use signature-based scoring")
    parser.add_argument("--max-epochs", type=int, default=15000,
                        help="Max epochs for c2l training (default 15k)")
    parser.add_argument("--force", action="store_true",
                        help="Rerun even if output CSV exists")
    args = parser.parse_args()

    if OUT_DECONV_CSV.exists() and not args.force:
        print(f"  Output exists: {OUT_DECONV_CSV} — pass --force to rerun")
        sys.exit(0)

    # Load normalized Vu (for signature scoring) and raw counts (for c2l)
    print(f"== Loading Vu spatial data ==")
    print(f"  Normalized: {VU_NORM_H5AD}")
    adata_norm = sc.read_h5ad(VU_NORM_H5AD)
    print(f"  Spots: {adata_norm.n_obs}, Genes: {adata_norm.n_vars}")
    print(f"  obs cols: {list(adata_norm.obs.columns)[:20]}")
    if "individual" not in adata_norm.obs.columns:
        raise RuntimeError("No 'individual' column in Vu obs")

    method = "cell2location"
    abund_cols = []

    if not args.use_fallback:
        print(f"\n== Running cell2location on Vu spots ==")
        try:
            init_gpu()
            print(f"  Loading raw counts: {VU_RAW_H5AD}")
            if VU_RAW_H5AD.exists():
                adata_raw = sc.read_h5ad(VU_RAW_H5AD)
            else:
                # If raw not separately stored, derive raw from .layers['counts']
                print(f"  Raw not found → using normalized as input (suboptimal)")
                adata_raw = adata_norm.copy()
                if "counts" in adata_raw.layers:
                    adata_raw.X = adata_raw.layers["counts"].copy()
            print(f"  Raw shape: {adata_raw.n_obs} × {adata_raw.n_vars}")

            print(f"  Loading reference signatures: {INF_AVER_CSV}")
            inf_aver = pd.read_csv(INF_AVER_CSV, index_col=0)
            # Strip the 'means_per_cluster_mu_fg_' prefix from columns
            inf_aver.columns = [
                c.replace("means_per_cluster_mu_fg_", "")
                for c in inf_aver.columns
            ]
            print(f"  Reference: {inf_aver.shape[0]} genes × "
                  f"{inf_aver.shape[1]} cell types")
            print(f"  Cell types: {list(inf_aver.columns)}")

            adata_dec = run_cell2location(
                adata_raw, inf_aver,
                max_epochs=args.max_epochs,
            )
            adata_dec.write_h5ad(OUT_H5AD)
            print(f"  Saved: {OUT_H5AD}")

            # Pull q05 abundance into obs
            abund = adata_dec.obsm.get("q05_cell_abundance_w_sf")
            if abund is None:
                # Different cell2location versions name this differently
                for k in adata_dec.obsm.keys():
                    if "q05" in k and "abundance" in k:
                        abund = adata_dec.obsm[k]
                        break
            if abund is None:
                raise RuntimeError("Could not find q05 abundance matrix in obsm")
            if isinstance(abund, pd.DataFrame):
                ct_cols = list(abund.columns)
                for ct in ct_cols:
                    cleaned = ct.replace("means_per_cluster_mu_fg_", "")
                    adata_norm.obs[f"c2l_{cleaned}"] = abund[ct].values
                    abund_cols.append(f"c2l_{cleaned}")
            else:
                ct_cols = list(inf_aver.columns)
                for i, ct in enumerate(ct_cols):
                    adata_norm.obs[f"c2l_{ct}"] = abund[:, i]
                    abund_cols.append(f"c2l_{ct}")
        except Exception as e:
            import traceback
            print(f"  cell2location FAILED: {e}")
            traceback.print_exc()
            print(f"  Falling back to signature-based scoring")
            method = "signature_softmax"
            adata_norm, abund_cols = fallback_signature_scoring(adata_norm)
    else:
        method = "signature_softmax"
        adata_norm, abund_cols = fallback_signature_scoring(adata_norm)

    # Compute hepatocyte fraction & dominant flag
    hep_col = next(
        (c for c in abund_cols if "Hepatocyt" in c), None
    )
    if hep_col is None:
        raise RuntimeError(f"No hepatocyte abundance col found in {abund_cols}")

    # Normalize abundances to fractions (sum to 1) for proportion calc
    abund_mat = adata_norm.obs[abund_cols].values.astype(float)
    abund_sum = abund_mat.sum(axis=1, keepdims=True)
    abund_sum[abund_sum == 0] = 1.0
    frac_mat = abund_mat / abund_sum
    hep_idx = abund_cols.index(hep_col)
    adata_norm.obs["hep_fraction"] = frac_mat[:, hep_idx]
    adata_norm.obs["hep_dominant_80"] = adata_norm.obs["hep_fraction"] >= 0.80

    print(f"  Method: {method}")
    print(f"  Hepatocyte-dominant (≥80%): "
          f"{adata_norm.obs['hep_dominant_80'].sum()}/{adata_norm.n_obs} spots "
          f"({100 * adata_norm.obs['hep_dominant_80'].mean():.1f}%)")

    # F3a / F3b scoring
    add_substate_scores(adata_norm)

    # Spatial coords
    if "spatial" in adata_norm.obsm:
        adata_norm.obs["x"] = adata_norm.obsm["spatial"][:, 0]
        adata_norm.obs["y"] = adata_norm.obsm["spatial"][:, 1]

    keep_cols = ["sample_id", "individual", "x", "y", "hep_fraction",
                 "hep_dominant_80", "li_f3a", "li_f3b", "li_f3a_z",
                 "li_f3b_z", "mixing_axis"] + abund_cols
    keep_cols = [c for c in keep_cols if c in adata_norm.obs.columns]

    out_df = adata_norm.obs[keep_cols].copy()
    out_df.insert(0, "spot_id", adata_norm.obs_names)
    out_df["method"] = method
    out_df.to_csv(OUT_DECONV_CSV, index=False)
    print(f"  Saved: {OUT_DECONV_CSV} ({len(out_df)} spots)")

    # Save updated AnnData with all annotations for downstream scripts
    out_norm_h5ad = OUT_MODEL_DIR / "vu_deconvolved_substate.h5ad"
    adata_norm.write_h5ad(out_norm_h5ad)
    print(f"  Saved: {out_norm_h5ad}")
    print(f"  Done.")


if __name__ == "__main__":
    main()
