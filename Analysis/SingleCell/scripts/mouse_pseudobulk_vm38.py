#!/usr/bin/env python3
"""
Mouse scRNA-seq Pseudobulk Generation (GENCODE vM38 re-quant).

vM38 variant of mouse_pseudobulk_de.py. Identical marker-based cell-type
annotation and raw-count aggregation, but reads the NEW vM38 Cell Ranger
``filtered_feature_bc_matrix`` directories (one per sample) instead of the
integrated per-sample h5ads. The whole point of the vM38 re-quant is to recover
the ~68 predicted lncRNA loci that vM37 (refdata-gex-GRCm39-2024-A) was missing,
so this script preserves the vM38 **gene_id** (Ensembl ENSMUSG) alongside the
symbol all the way into the pseudobulk CSVs. The downstream substrate
(compute_mouse_hep_specificity_vm38.R) then joins the Cas13 library on the native
vM38 gene_id with no GTF symbol->id remapping.

Differences vs mouse_pseudobulk_de.py:
  - Input: vM38 Cell Ranger MTX (barcodes/features/matrix .tsv/.mtx .gz) per
    sample, NOT integration/input_h5ad/mouse/*.h5ad.
  - Cell Ranger features.tsv col1=gene_id (ENSMUSG), col2=symbol, col3=type.
    We read with var_names="gene_symbols" (so marker scoring is symbol-keyed,
    identical to the original) and keep gene_ids in adata.var['gene_ids'].
  - Output pseudobulk CSVs are keyed on gene_id (index) AND carry a gene_symbol
    column, so the vM38 Ensembl id is recoverable downstream.
  - Output dir: results_gpu_v2/mouse_sc/pseudobulk_vm38/ (parallel to the vM37
    pseudobulk/, which is preserved unchanged).

Marker dict, score_cell_types(), assign_cell_types() and the per-cell-type raw
count aggregation are byte-for-byte the same logic as the original so the only
variable changed is the reference build.

Run (rapids_singlecell env, on a COMPUTE NODE -- NOT login):
  micromamba run -n rapids_singlecell python mouse_pseudobulk_vm38.py \
      [--cellranger-root <dir>] [--manifest <tsv>]

Defaults:
  --cellranger-root  results_gpu_v2/mouse_sc/vm38_requant/cellranger
                     (per-sample outs at <root>/<SRR>/outs/filtered_feature_bc_matrix)
  --manifest         results_gpu_v2/mouse_sc/vm38_requant/sample_manifest.tsv
                     (produced by the re-quant agent; needs srr,dataset,condition cols)
  env overrides: CAS13_VM38_CELLRANGER_ROOT, CAS13_VM38_MANIFEST
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
SC_DIR = BASE / "Analysis/SingleCell"
VM38_DIR = SC_DIR / "results_gpu_v2/mouse_sc/vm38_requant"
DEFAULT_CR_ROOT = VM38_DIR / "cellranger"
DEFAULT_MANIFEST = VM38_DIR / "sample_manifest.tsv"
OUT_DIR = SC_DIR / "results_gpu_v2/mouse_sc/pseudobulk_vm38"

# ---------------------------------------------------------------------------
# Mouse liver cell type marker gene sets  (IDENTICAL to mouse_pseudobulk_de.py)
# ---------------------------------------------------------------------------
MOUSE_MARKERS = {
    "Hepatocytes": [
        "Alb", "Apoe", "Cyp2e1", "Cyp3a11", "Tat", "Arg1", "Ass1",
        "Serpina1a", "Ttr", "Apob", "Fabp1", "Apoa1", "Cyp7a1", "Hmgcr",
    ],
    "Cholangiocytes": [
        "Krt7", "Krt19", "Epcam", "Spp1", "Anxa4", "Sox9", "Cftr",
    ],
    "LSECs": [
        "Pecam1", "Cd34", "Vwf", "Lyve1", "Stab2", "Oit3", "Fcgr2b",
    ],
    "Kupffer_cells": [
        "Csf1r", "Adgre1", "Clec4f", "Timd4", "Vsig4", "Cd68", "Cx3cr1",
    ],
    "Stellate_cells": [
        "Acta2", "Col1a1", "Dcn", "Pdgfrb", "Lrat", "Des", "Reln",
    ],
    "T_cells": [
        "Cd3d", "Cd3e", "Trdc", "Cd8a", "Cd4", "Nkg7",
    ],
    "NK_cells": [
        "Nkg7", "Gzma", "Gzmb", "Ncr1", "Klrb1c",
    ],
    "B_cells": [
        "Cd79a", "Ms4a1", "Cd19", "Ighm",
    ],
    "Monocytes": [
        "Ly6c2", "Ccr2", "Csf3r", "S100a8", "S100a9",
    ],
    "pDCs": [
        "Siglech", "Bst2", "Irf7",
    ],
}


def score_cell_types(counts_csr: sp.csr_matrix, gene_names: list) -> pd.DataFrame:
    """
    Score each cell for each marker gene set using mean normalized expression.
    Returns a DataFrame (cells x cell_types) with scores.

    Method (IDENTICAL to mouse_pseudobulk_de.py):
      1. Normalize counts to 10K (median), log1p
      2. For each cell type, compute mean expression of available markers
    """
    counts_csr = counts_csr.astype(np.float32)
    row_sums = np.asarray(counts_csr.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0  # avoid division by zero
    scale = 10000.0 / row_sums

    from scipy.sparse import diags
    norm_mat = diags(scale) @ counts_csr
    norm_mat.data = np.log1p(norm_mat.data)

    gene_index = {g: i for i, g in enumerate(gene_names)}
    n_cells = counts_csr.shape[0]

    scores = {}
    for ct, markers in MOUSE_MARKERS.items():
        avail = [g for g in markers if g in gene_index]
        if not avail:
            scores[ct] = np.zeros(n_cells, dtype=np.float32)
            continue
        idx = [gene_index[g] for g in avail]
        sub = norm_mat[:, idx].toarray()
        scores[ct] = sub.mean(axis=1).astype(np.float32)

    return pd.DataFrame(scores)


def assign_cell_types(score_df: pd.DataFrame) -> np.ndarray:
    """Assign each cell to the highest-scoring cell type."""
    return score_df.idxmax(axis=1).values


def _read_cellranger_mtx(mtx_dir: Path):
    """
    Read a vM38 Cell Ranger filtered_feature_bc_matrix into (X_csr, gene_symbols,
    gene_ids). var_names="gene_symbols" makes marker scoring symbol-keyed (identical
    to the h5ad path), while gene_ids preserves the native vM38 ENSMUSG ids.
    """
    import scanpy as sc

    adata = sc.read_10x_mtx(
        str(mtx_dir), var_names="gene_symbols", gex_only=True, make_unique=False
    )
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    elif not isinstance(X, sp.csr_matrix):
        X = X.tocsr()

    gene_symbols = list(adata.var_names)
    if "gene_ids" in adata.var:
        gene_ids = list(adata.var["gene_ids"].values)
    else:
        # Fallback: no gene_ids column -> emit symbols as ids (logged loudly).
        log.warning("  features table has no gene_ids column; using symbols as gene_id")
        gene_ids = gene_symbols
    return X, gene_symbols, gene_ids


def process_sample(mtx_dir: Path, sample_id: str) -> dict | None:
    """
    Load one sample's vM38 Cell Ranger matrix, annotate cells, return dict of
    {cell_type: 1D summed-count array}, plus gene_symbols and gene_ids.
    """
    try:
        X, gene_symbols, gene_ids = _read_cellranger_mtx(mtx_dir)
    except Exception as e:
        log.warning(f"  Could not read {mtx_dir}: {e}")
        return None

    n_cells, n_genes = X.shape
    if n_cells == 0 or n_genes == 0:
        log.warning(f"  Empty matrix in {mtx_dir}")
        return None

    log.info(f"  {sample_id}: {n_cells:,} cells x {n_genes:,} genes")

    # Score and assign cell types (symbol-keyed, identical to original).
    score_df = score_cell_types(X, gene_symbols)
    assignments = assign_cell_types(score_df)

    result = {"gene_symbols": gene_symbols, "gene_ids": gene_ids}
    for ct in MOUSE_MARKERS:
        mask = assignments == ct
        if mask.sum() > 0:
            agg = np.asarray(X[mask].sum(axis=0)).ravel()
        else:
            agg = np.zeros(n_genes, dtype=np.float64)
        result[ct] = agg

    ct_counts = {ct: int((assignments == ct).sum()) for ct in MOUSE_MARKERS}
    log.info(
        "  Cell type distribution: "
        + ", ".join(
            f"{k}={v}"
            for k, v in sorted(ct_counts.items(), key=lambda x: -x[1])[:5]
        )
    )
    return result


def _resolve_mtx_dir(cr_root: Path, srr: str) -> Path:
    """
    Locate a sample's filtered_feature_bc_matrix. Canonical Cell Ranger layout is
    <root>/<SRR>/outs/filtered_feature_bc_matrix. Tolerate a couple of common
    variants so a quant agent's exact dir choice doesn't break us.
    """
    candidates = [
        cr_root / srr / "outs" / "filtered_feature_bc_matrix",
        cr_root / srr / "filtered_feature_bc_matrix",
        cr_root / srr,
    ]
    for c in candidates:
        if (c / "matrix.mtx.gz").exists() or (c / "matrix.mtx").exists():
            return c
    return candidates[0]  # canonical; caller logs the miss


def main():
    ap = argparse.ArgumentParser(description="vM38 mouse scRNA pseudobulk generation")
    ap.add_argument(
        "--cellranger-root",
        default=os.environ.get("CAS13_VM38_CELLRANGER_ROOT", str(DEFAULT_CR_ROOT)),
        help="Root holding per-sample Cell Ranger outs (<root>/<SRR>/outs/...)",
    )
    ap.add_argument(
        "--manifest",
        default=os.environ.get("CAS13_VM38_MANIFEST", str(DEFAULT_MANIFEST)),
        help="TSV with srr,dataset,condition columns (vm38 sample manifest)",
    )
    ap.add_argument("--out-dir", default=str(OUT_DIR), help="Output pseudobulk dir")
    args = ap.parse_args()

    cr_root = Path(args.cellranger_root)
    manifest_path = Path(args.manifest)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("=== Mouse scRNA-seq Pseudobulk Generation (vM38) ===")
    log.info(f"Cell Ranger root: {cr_root}")
    log.info(f"Manifest:         {manifest_path}")
    log.info(f"Output dir:       {out_dir}")

    if not manifest_path.exists():
        log.error(f"Manifest not found: {manifest_path}")
        sys.exit(1)

    # Manifest may be TSV or CSV; sniff the separator.
    sep = "\t" if manifest_path.suffix.lower() in (".tsv", ".txt") else ","
    manifest = pd.read_csv(manifest_path, sep=sep)
    # Tolerate species column if present; keep mouse only.
    if "species" in manifest.columns:
        manifest = manifest[manifest["species"] == "mouse"].copy()
    # Required columns
    for col in ("srr", "dataset", "condition"):
        if col not in manifest.columns:
            log.error(f"Manifest missing required column '{col}'. Has: {list(manifest.columns)}")
            sys.exit(1)
    log.info(f"Mouse samples in manifest: {len(manifest)}")
    log.info(f"Conditions: {manifest['condition'].value_counts().to_dict()}")

    # Process each sample.
    all_results = {}  # srr -> {meta + cell_type arrays}
    gene_symbols_ref = None
    gene_ids_ref = None

    for _, row in manifest.iterrows():
        srr = str(row["srr"])
        condition = row["condition"]
        dataset = row["dataset"]
        mtx_dir = _resolve_mtx_dir(cr_root, srr)

        if not ((mtx_dir / "matrix.mtx.gz").exists() or (mtx_dir / "matrix.mtx").exists()):
            log.warning(f"  Missing Cell Ranger matrix for {srr} at {mtx_dir} -- skipping")
            continue

        log.info(f"Processing {srr} (dataset={dataset}, condition={condition})")
        result = process_sample(mtx_dir, srr)
        if result is None:
            continue

        if gene_symbols_ref is None:
            gene_symbols_ref = result["gene_symbols"]
            gene_ids_ref = result["gene_ids"]

        # Enforce identical feature ordering across samples (vM38 ref is fixed, so
        # this should always hold; guard against an accidental ref mismatch).
        if result["gene_ids"] != gene_ids_ref:
            log.warning(f"  Gene IDs differ for {srr} (feature set mismatch) -- skipping")
            continue

        all_results[srr] = {
            "condition": condition,
            "dataset": dataset,
            **{ct: result[ct] for ct in MOUSE_MARKERS},
        }

    if not all_results:
        log.error("No samples processed successfully.")
        sys.exit(1)

    log.info(f"\nSuccessfully processed: {len(all_results)} samples")

    sample_ids = list(all_results.keys())
    conditions = [all_results[s]["condition"] for s in sample_ids]
    datasets = [all_results[s]["dataset"] for s in sample_ids]

    # Save per-cell-type pseudobulk matrices (genes x samples), keyed on gene_id
    # with a gene_symbol column so the native vM38 Ensembl id is recoverable.
    for ct in MOUSE_MARKERS:
        count_arrays = [all_results[s][ct] for s in sample_ids]
        count_mat = np.column_stack(count_arrays)  # genes x samples

        df = pd.DataFrame(count_mat, columns=sample_ids)
        df.insert(0, "gene_symbol", gene_symbols_ref)
        df.insert(0, "gene_id", gene_ids_ref)

        out_path = out_dir / f"{ct}_pseudobulk.csv"
        df.to_csv(out_path, index=False)
        log.info(f"  Saved: {out_path} ({df.shape[0]:,} genes x {len(sample_ids)} samples)")

    # Sample metadata.
    meta_df = pd.DataFrame(
        {"sample": sample_ids, "condition": conditions, "dataset": datasets}
    )
    meta_path = out_dir.parent / "sample_metadata_vm38.csv"
    meta_df.to_csv(meta_path, index=False)
    log.info(f"  Saved: {meta_path}")

    log.info("=== Pseudobulk generation (vM38) complete ===")
    log.info(f"Output directory: {out_dir}")


if __name__ == "__main__":
    main()
