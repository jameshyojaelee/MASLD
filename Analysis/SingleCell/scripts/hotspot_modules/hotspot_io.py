"""Shared paths and loaders for Hotspot module pipeline.

Single source of truth for atlas paths, donor metadata, reference panels.
"""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SC_ROOT = PROJECT_ROOT / "Analysis/SingleCell"
RESULTS = SC_ROOT / "results_gpu_v2/hotspot_modules"
INPUTS = SC_ROOT / "results_gpu_v2/mcp/inputs"

ATLAS_GLOBAL = SC_ROOT / "integration/output/human/scalesc_human_annotated_celltypist.h5ad"

CELLTYPE_SUBSETS = {
    "global":             ATLAS_GLOBAL,
    "hepatocytes":        INPUTS / "atlas_cnmf_hepatocytes.h5ad",
    "macrophages":        INPUTS / "atlas_cnmf_macrophages.h5ad",
    "fibroblasts":        INPUTS / "atlas_cnmf_fibroblasts.h5ad",
    "endothelial_cells":  INPUTS / "atlas_cnmf_endothelial_cells.h5ad",
    "cholangiocytes":     INPUTS / "atlas_cnmf_cholangiocytes.h5ad",
    "tcells":             INPUTS / "atlas_cnmf_tcells.h5ad",
}
RUN_ORDER = list(CELLTYPE_SUBSETS.keys())

DONOR_META = INPUTS / "donor_metadata.tsv"
DONOR_META_EXTENDED = (
    SC_ROOT / "results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)

REF_DIR = PROJECT_ROOT / "data/reference_signatures"
CNMF_GENE_SPECTRA_K16 = (
    SC_ROOT / "results_gpu_v2/mcp/cnmf_runs/global/global.gene_spectra_score.k_16.dt_0_03.txt"
)
CNMF_USAGES_K16 = (
    SC_ROOT / "results_gpu_v2/mcp/cnmf_runs/global/global.usages.k_16.dt_0_03.consensus.txt"
)
BULK_NMF_LABELS = PROJECT_ROOT / "RNA-seq/results/subtypes/program_labels_protonly.csv"
BULK_NMF_RDS = PROJECT_ROOT / "RNA-seq/results/subtypes/nmf_results_cache_clean.rds"

GENE_BIOTYPES = PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
MULTI_EVIDENCE_ATLAS = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"


def load_atlas(cell_type: str, smoke: bool = False) -> ad.AnnData:
    """Load AnnData for a cell-type subset. `smoke=True` returns 5,000-cell sample."""
    if cell_type not in CELLTYPE_SUBSETS:
        raise ValueError(
            f"Unknown cell_type {cell_type!r}. Valid: {RUN_ORDER}"
        )
    path = CELLTYPE_SUBSETS[cell_type]
    adata = ad.read_h5ad(path)
    if smoke:
        n = min(5000, adata.n_obs)
        rng = np.random.default_rng(42)
        idx = rng.choice(adata.n_obs, size=n, replace=False)
        adata = adata[idx].copy()
    return adata


def load_donor_metadata() -> pd.DataFrame:
    """Merge donor_metadata + extended on `sample`."""
    base = pd.read_csv(DONOR_META, sep="\t")
    ext = pd.read_csv(DONOR_META_EXTENDED, sep="\t")
    keep = [
        c for c in ext.columns
        if c not in base.columns or c == "sample"
    ]
    return base.merge(ext[keep], on="sample", how="left")


def load_gene_biotypes() -> pd.DataFrame:
    return pd.read_csv(GENE_BIOTYPES, sep="\t")


def hotspot_outdir(cell_type: str) -> Path:
    p = RESULTS / cell_type
    p.mkdir(parents=True, exist_ok=True)
    return p
