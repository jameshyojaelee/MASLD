"""352a_scanorama_integration.py — third independent integration method.

Addresses critique #9: v1↔v2 fstage LMM Jaccard top-50 is 0.27 across two
scVI-based builds, suggesting heavy method-dependence. To test method-
dependence properly, we need a third integration latent that is independent
of scVI (and ideally independent of Harmony too).

Scanorama (Hie et al. Nat. Biotech. 2019) is a panoramic-stitching integration
algorithm that uses mutual nearest neighbours + SVD; it shares no internal
machinery with scVI (VAE / neural variational inference) or Harmony
(per-cluster centroid adjustment).

Strategy: build a Scanorama-integrated latent over the 7 per-dataset
QC-harmonized hepatocyte h5ads (the same inputs that fed 01_integrate_atlas_v2).
Project the integrated latent back onto cells in the v2 atlas via the same
sample/cell-barcode keys.

Outputs:
  - results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_scanorama.h5ad
    Same shape as hepatocyte_atlas_v2_annotated.h5ad, with `obsm['X_scanorama']`
    (100-dim) replacing `obsm['X_scVI']` as the latent to be used in
    downstream F-stage prediction (352b).
"""
from __future__ import annotations

import gc
import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scanorama

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PER_DATASET_DIR = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset"
ATLAS_ANNOTATED = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad"
ATLAS_OUT = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_scanorama.h5ad"
HVG_TSV = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/scvi/canonical_hvg_genes_v2.tsv"

N_LATENT = 100  # Scanorama default; matches the size that downstream classifier expects
N_HVG_FALLBACK = 3000


def load_hvg_gene_list() -> list[str] | None:
    """Use the same HVG list that scVI v2 trained on, so Scanorama operates on
    the same gene panel — apples-to-apples method comparison."""
    if not HVG_TSV.exists():
        log.warning(f"HVG TSV not found at {HVG_TSV}; will re-select HVGs")
        return None
    df = pd.read_csv(HVG_TSV, sep="\t")
    col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
    hvg = df[col].astype(str).tolist()
    log.info(f"loaded {len(hvg)} canonical v2 HVGs from {HVG_TSV}")
    return hvg


def load_per_dataset_h5ads(hvg: list[str] | None) -> tuple[list[ad.AnnData], list[str]]:
    """Load each of the 7 per-dataset filtered h5ads. We restrict to the
    hepatocyte subset by intersecting cell barcodes with the annotated atlas."""
    log.info(f"loading annotated atlas to get hep barcode set: {ATLAS_ANNOTATED}")
    atlas = ad.read_h5ad(ATLAS_ANNOTATED, backed="r")
    hep_barcodes = set(atlas.obs_names.astype(str).tolist())
    hep_obs_per_dataset = atlas.obs[["sample", "dataset"]].copy()
    hep_obs_per_dataset.index = atlas.obs_names.astype(str)
    log.info(f"  {len(hep_barcodes):,} hep barcodes; {hep_obs_per_dataset['dataset'].nunique()} datasets")

    paths = sorted(PER_DATASET_DIR.glob("*_v2.h5ad"))
    log.info(f"per-dataset files: {len(paths)}")

    adatas: list[ad.AnnData] = []
    ds_labels: list[str] = []
    for p in paths:
        ds_name = p.stem.replace("_v2", "")
        log.info(f"  reading {p.name}")
        a = ad.read_h5ad(p)
        # restrict to hep barcodes that exist in the annotated atlas
        keep = a.obs_names.astype(str).isin(hep_barcodes)
        if keep.sum() == 0:
            log.warning(f"    {ds_name}: 0 hep cells; skipping")
            continue
        a = a[keep, :].copy()
        # Ensure dataset column matches the dataset name (some files had
        # the integer-coded stomp; restore from filename)
        a.obs["dataset"] = ds_name
        log.info(f"    {ds_name}: {a.n_obs:,} hep cells, {a.n_vars} genes")
        adatas.append(a)
        ds_labels.append(ds_name)

    # Intersect genes
    common = set.intersection(*[set(a.var_names.astype(str)) for a in adatas])
    log.info(f"common genes across datasets: {len(common):,}")
    if hvg is not None:
        keep_genes = sorted(set(hvg) & common)
        log.info(f"intersecting with canonical v2 HVG list -> {len(keep_genes)} genes")
        if len(keep_genes) < 500:
            log.warning("HVG intersection too small; falling back to common-genes only")
            keep_genes = sorted(common)
    else:
        keep_genes = sorted(common)
    for i, a in enumerate(adatas):
        adatas[i] = a[:, keep_genes].copy()
    return adatas, ds_labels


def main() -> int:
    if ATLAS_OUT.exists():
        log.info(f"output already exists at {ATLAS_OUT}; skipping (delete to re-run)")
        return 0

    hvg = load_hvg_gene_list()
    adatas, ds_labels = load_per_dataset_h5ads(hvg)
    log.info(f"running scanorama.integrate on {len(adatas)} datasets")

    # Normalize + log1p before Scanorama (it expects log-normalized expression).
    # Use counts layer if present.
    for i, a in enumerate(adatas):
        if "counts" in a.layers:
            a.X = a.layers["counts"].copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        adatas[i] = a

    # Scanorama 1.7.4: correct_scanpy(adatas, return_dimred=True, dimred=K)
    # returns (corrected_adatas_list, genes). The low-dim latent is stored
    # in each corrected adata's obsm['X_scanorama'].
    result = scanorama.correct_scanpy(adatas, return_dimred=True,
                                      dimred=N_LATENT)
    if isinstance(result, tuple) and len(result) == 2:
        corrected, _ = result
    elif isinstance(result, tuple) and len(result) == 3:
        # Older API: (integrated, corrected, genes)
        _, corrected, _ = result
    else:
        corrected = result
    integrated = corrected
    log.info(f"scanorama integration done; {len(integrated)} batches")

    # Concatenate latents in order. Some Scanorama versions write to
    # obsm['X_scanorama']; others to obsm['Scanorama']. Try both.
    def _get_latent(a):
        for k in ("X_scanorama", "Scanorama", "scanorama"):
            if k in a.obsm:
                return np.asarray(a.obsm[k])
        raise KeyError(f"no scanorama obsm key found; keys = {list(a.obsm)}")
    Z_list = [_get_latent(a) for a in integrated]
    obs_list = [a.obs.copy() for a in integrated]
    Z = np.concatenate(Z_list, axis=0)
    obs = pd.concat(obs_list, axis=0)
    log.info(f"concatenated latent shape: {Z.shape}; obs shape: {obs.shape}")

    # Stamp X_scanorama onto a copy of the annotated atlas, aligned by barcode
    log.info(f"re-loading annotated atlas to stamp X_scanorama")
    atlas = ad.read_h5ad(ATLAS_ANNOTATED)
    # Build a barcode -> latent-row mapping
    barcodes_int = []
    for i, a in enumerate(integrated):
        barcodes_int.extend(a.obs_names.astype(str).tolist())
    bc_to_row = {b: r for r, b in enumerate(barcodes_int)}

    # Reorder Z to match atlas.obs_names
    Z_aligned = np.full((atlas.n_obs, N_LATENT), np.nan, dtype=np.float32)
    n_missing = 0
    for r, bc in enumerate(atlas.obs_names.astype(str)):
        idx = bc_to_row.get(bc)
        if idx is None:
            n_missing += 1
        else:
            Z_aligned[r, :] = Z[idx, :]
    log.info(f"alignment: {atlas.n_obs - n_missing:,} cells with X_scanorama, "
             f"{n_missing} missing")

    atlas.obsm["X_scanorama"] = Z_aligned

    log.info(f"writing -> {ATLAS_OUT}")
    atlas.write_h5ad(ATLAS_OUT, compression="gzip")
    log.info("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
