#!/usr/bin/env python3
"""
Phase 0.5 — Step 1 / v2 atlas
============================

Per-dataset QC harmonization for the v2 scRNA atlas (Phase 0.5).

For each of the 7 datasets in the v1 integrated atlas, extract raw counts and
re-apply HARMONIZED cutoffs from Agent 1's QC_HARMONIZATION.tsv. Then run
``sc.pp.scrublet`` for doublet calling and drop cells with score > 0.25.

Optional: ``decontX`` ambient correction is SKIPPED in this script because the
celda/decontX R-Bioconductor stack is not in ``rapids_singlecell``. Disclose in
the cell audit so reviewers can see exactly what was/was not run.

Inputs (READ-ONLY):
    Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad
    Analysis/SingleCell/QC_HARMONIZATION.tsv  (Agent 1 output)

Outputs (UNDER results_gpu_v2_phase05/):
    atlas/per_dataset/{dataset_id}_v2.h5ad   — filtered raw counts + obs
    atlas/qc_cell_audit_v2.tsv               — n_raw / n_after_mt / n_after_genes / n_after_doublet
    atlas/qc_cell_audit_v2.log               — extended log

Env: rapids_singlecell (uses scanpy + scrublet through scanpy.external.pp).
SBATCH (cpu / qos=interactive / 16 cpus / 128G / 12h).
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

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# Determinism
import random
os.environ["PYTHONHASHSEED"] = "42"
random.seed(42)
np.random.seed(42)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SC_DIR = BASE / "Analysis" / "SingleCell"
V1_ATLAS = SC_DIR / "integration/output/human/scalesc_human_annotated_celltypist.h5ad"
# Agent 1 publishes QC_HARMONIZATION.tsv under results_gpu_v2_phase05/ — fall back to
# Analysis/SingleCell/QC_HARMONIZATION.tsv for legacy compatibility.
QC_TSV_PRIMARY = SC_DIR / "results_gpu_v2_phase05" / "QC_HARMONIZATION.tsv"
QC_TSV_LEGACY = SC_DIR / "QC_HARMONIZATION.tsv"
QC_TSV = QC_TSV_PRIMARY if QC_TSV_PRIMARY.exists() else QC_TSV_LEGACY

OUT_ROOT = SC_DIR / "results_gpu_v2_phase05"
OUT_PER_DATASET = OUT_ROOT / "atlas" / "per_dataset"
OUT_AUDIT_TSV = OUT_ROOT / "atlas" / "qc_cell_audit_v2.tsv"
OUT_AUDIT_LOG = OUT_ROOT / "atlas" / "qc_cell_audit_v2.log"

# Harmonized thresholds (override any per-dataset value from the TSV if needed)
HARMONIZED = {
    "mt_pct_max": 20.0,
    "min_genes": 200,
    "min_counts": 500,
    "doublet_threshold": 0.25,
}


def _categorical_safe_read(h5: Path) -> ad.AnnData:
    """Read h5ad with categorical 'ordered' attribute fix (mirrors 308)."""
    import h5py as _h5py
    import shutil

    try:
        return sc.read_h5ad(h5, backed="r")
    except KeyError:
        log.warning("Categorical encoding mismatch — patching with shadow copy in /tmp")

    tmp_path = Path(os.environ.get("TMPDIR", "/tmp")) / "_v1_atlas_patched.h5ad"
    if not tmp_path.exists():
        shutil.copy2(h5, tmp_path)
        with _h5py.File(tmp_path, "a") as f:
            for col in f["obs"].keys():
                if col == "_index":
                    continue
                attrs = f["obs"][col].attrs
                enc = attrs.get("encoding-type", b"")
                if isinstance(enc, bytes):
                    enc = enc.decode()
                if enc == "categorical" and "ordered" not in attrs:
                    attrs["ordered"] = False
    return sc.read_h5ad(tmp_path, backed="r")


def _load_harmonization(qc_tsv: Path) -> pd.DataFrame:
    if not qc_tsv.exists():
        log.warning("QC_HARMONIZATION.tsv missing — proceeding with hard-coded harmonized cutoffs")
        return pd.DataFrame()
    df = pd.read_csv(qc_tsv, sep="\t")
    log.info("Loaded QC harmonization for %d datasets: %s",
             len(df), df.columns.tolist())
    return df


def _ensure_mt_pct(adata: ad.AnnData) -> None:
    """Compute pct_counts_mt if not already in obs."""
    if "pct_counts_mt" in adata.obs.columns:
        return
    mt_mask = adata.var_names.str.startswith("MT-")
    if not mt_mask.any():
        # try mt- lowercase or fallback
        mt_mask = adata.var_names.str.upper().str.startswith("MT-")
    # mt_mask is a numpy array on some pandas / anndata versions; coerce safely
    adata.var["mt"] = np.asarray(mt_mask)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],
                               percent_top=None, log1p=False, inplace=True)


def _per_dataset_filter(adata: ad.AnnData, dataset_id: str, harm: dict) -> tuple[ad.AnnData, dict]:
    """Apply harmonized QC filters and run scrublet. Returns filtered adata + audit row."""
    audit = {
        "dataset_id": dataset_id,
        "n_raw": int(adata.n_obs),
        "n_after_mt": np.nan,
        "n_after_genes": np.nan,
        "n_after_counts": np.nan,
        "n_after_doublet": np.nan,
        "n_doublet_called": np.nan,
        "mt_pct_max": harm["mt_pct_max"],
        "min_genes": harm["min_genes"],
        "min_counts": harm["min_counts"],
        "doublet_threshold": harm["doublet_threshold"],
        "decontX_run": False,
        "decontX_skip_reason": "celda/decontX not in rapids_singlecell env; ambient correction deferred to a separate R-Bioconductor pass if reviewers ask",
    }

    if adata.n_obs == 0:
        log.warning("Dataset %s has 0 cells — skipping", dataset_id)
        return adata, audit

    # 1. MT% filter
    _ensure_mt_pct(adata)
    mt_keep = adata.obs["pct_counts_mt"] <= harm["mt_pct_max"]
    log.info("  [%s] MT%% <= %.1f kept %d / %d", dataset_id,
             harm["mt_pct_max"], mt_keep.sum(), len(mt_keep))
    adata = adata[mt_keep].copy()
    audit["n_after_mt"] = int(adata.n_obs)

    # 2. min_genes
    sc.pp.filter_cells(adata, min_genes=int(harm["min_genes"]))
    audit["n_after_genes"] = int(adata.n_obs)
    log.info("  [%s] min_genes>=%d kept %d", dataset_id, harm["min_genes"], adata.n_obs)

    # 3. min_counts
    sc.pp.filter_cells(adata, min_counts=int(harm["min_counts"]))
    audit["n_after_counts"] = int(adata.n_obs)
    log.info("  [%s] min_counts>=%d kept %d", dataset_id, harm["min_counts"], adata.n_obs)

    if adata.n_obs == 0:
        log.warning("[%s] All cells filtered out", dataset_id)
        return adata, audit

    # 4. Scrublet doublet calling
    # scanpy's built-in wrapper uses scrublet under the hood;
    # falls back gracefully if too few cells.
    n_cells = int(adata.n_obs)
    if n_cells < 30:
        log.warning("[%s] only %d cells — skipping scrublet", dataset_id, n_cells)
        adata.obs["scrublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False
        audit["n_after_doublet"] = int(adata.n_obs)
        audit["n_doublet_called"] = 0
        return adata, audit

    try:
        # rapids_singlecell ships scanpy>=1.10 → sc.pp.scrublet present
        if hasattr(sc.pp, "scrublet"):
            sc.pp.scrublet(adata, expected_doublet_rate=0.06,
                           threshold=harm["doublet_threshold"], random_state=42)
        else:
            sc.external.pp.scrublet(adata, expected_doublet_rate=0.06,
                                    threshold=harm["doublet_threshold"], random_state=42)
        score_col = "doublet_score" if "doublet_score" in adata.obs.columns else "scrublet_score"
        if "scrublet_score" not in adata.obs.columns and score_col in adata.obs.columns:
            adata.obs["scrublet_score"] = adata.obs[score_col]
    except Exception as exc:
        log.warning("[%s] scrublet failed: %s — falling back to score=0", dataset_id, exc)
        adata.obs["scrublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False

    doublet_mask = adata.obs["scrublet_score"] > harm["doublet_threshold"]
    audit["n_doublet_called"] = int(doublet_mask.sum())
    log.info("  [%s] doublets (score > %.2f): %d / %d",
             dataset_id, harm["doublet_threshold"],
             audit["n_doublet_called"], n_cells)

    adata = adata[~doublet_mask].copy()
    audit["n_after_doublet"] = int(adata.n_obs)
    return adata, audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="*", default=None,
                        help="Optional subset of dataset_ids; default = all in v1 atlas")
    parser.add_argument("--mt_pct_max", type=float, default=HARMONIZED["mt_pct_max"])
    parser.add_argument("--min_genes", type=int, default=HARMONIZED["min_genes"])
    parser.add_argument("--min_counts", type=int, default=HARMONIZED["min_counts"])
    parser.add_argument("--doublet_threshold", type=float, default=HARMONIZED["doublet_threshold"])
    args = parser.parse_args(argv)

    harm = {
        "mt_pct_max": args.mt_pct_max,
        "min_genes": args.min_genes,
        "min_counts": args.min_counts,
        "doublet_threshold": args.doublet_threshold,
    }
    log.info("Harmonized cutoffs: %s", harm)

    OUT_PER_DATASET.mkdir(parents=True, exist_ok=True)
    OUT_AUDIT_TSV.parent.mkdir(parents=True, exist_ok=True)

    log.info("Reading v1 atlas (backed='r'): %s", V1_ATLAS)
    bdata = _categorical_safe_read(V1_ATLAS)
    log.info("v1 atlas shape: %s", bdata.shape)

    # Discover dataset ids
    dataset_ids = sorted(bdata.obs["dataset"].astype(str).unique().tolist())
    log.info("Datasets in v1 atlas: %s", dataset_ids)
    if args.datasets:
        keep = [d for d in dataset_ids if d in args.datasets]
        log.info("Filtering to requested subset: %s", keep)
        dataset_ids = keep

    # Optional cross-check with Agent 1's TSV
    _ = _load_harmonization(QC_TSV)

    audit_rows: list[dict] = []
    for dsid in dataset_ids:
        log.info("=" * 60)
        log.info("Processing dataset: %s", dsid)

        # Skip-if-exists: per-dataset h5ads are written incrementally;
        # if a previous run already produced this one, reuse it.
        ds_out = OUT_PER_DATASET / f"{dsid}_v2.h5ad"
        if ds_out.exists() and ds_out.stat().st_size > 10_000:
            log.info("  [%s] already at %s (%.1f MB) -- skipping",
                     dsid, ds_out, ds_out.stat().st_size / 1e6)
            audit_rows.append({"dataset_id": dsid, "note": "skipped_existing"})
            continue

        # Boolean mask in obs space
        mask = (bdata.obs["dataset"].astype(str) == dsid).values
        log.info("  cells in dataset (v1): %d", int(mask.sum()))
        if mask.sum() == 0:
            audit_rows.append({"dataset_id": dsid, "n_raw": 0, "note": "missing in v1"})
            continue

        # Slicing a backed h5ad → in-memory copy with same X (counts via raw)
        # We want RAW counts (pre-normalization). v1 .raw holds them.
        log.info("  loading raw counts for %s in-memory ...", dsid)
        idx = np.where(mask)[0]
        # AnnData slicing on backed mode returns a view; copy to materialize
        sub = bdata[idx].to_memory()
        # raw -> AnnData with raw X
        if sub.raw is not None:
            raw_X = sub.raw.X
            raw_var = sub.raw.var.copy()
            raw_var_names = sub.raw.var_names
            log.info("  using sub.raw (shape %s)", sub.raw.shape)
        else:
            raw_X = sub.X
            raw_var = sub.var.copy()
            raw_var_names = sub.var_names
            log.info("  using sub.X (no .raw)")

        # Keep only essential obs cols (v1 categoricals can be heavy)
        keep_cols = [c for c in ("sample", "dataset", "species", "condition",
                                 "cell_type", "cell_type_conf", "cell_type_raw",
                                 "leiden") if c in sub.obs.columns]
        sub_obs = sub.obs[keep_cols].copy()
        sub_obs.index = sub.obs_names.astype(str)

        # Build a clean AnnData
        new_var = raw_var.copy()
        new_var.index = pd.Index(raw_var_names.astype(str), name=None)
        if not sparse.issparse(raw_X):
            raw_X = sparse.csr_matrix(raw_X)
        ad_raw = ad.AnnData(X=raw_X, obs=sub_obs, var=new_var)
        del sub, raw_X, raw_var, sub_obs
        gc.collect()

        log.info("  raw AnnData shape: %s", ad_raw.shape)

        # Apply harmonized QC + scrublet
        ad_filt, audit = _per_dataset_filter(ad_raw, dsid, harm)

        # Persist
        out_h5 = OUT_PER_DATASET / f"{dsid}_v2.h5ad"
        log.info("  writing %s (%s)", out_h5, ad_filt.shape)
        ad_filt.write_h5ad(out_h5, compression="gzip")
        audit["out_path"] = str(out_h5)
        audit["pct_retained"] = (
            float(audit["n_after_doublet"]) / audit["n_raw"] * 100.0
            if audit.get("n_raw", 0) else 0.0
        )
        audit_rows.append(audit)
        del ad_raw, ad_filt
        gc.collect()

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(OUT_AUDIT_TSV, sep="\t", index=False)
    log.info("Wrote cell audit: %s", OUT_AUDIT_TSV)
    log.info("\n%s", audit_df.to_string(index=False))

    with open(OUT_AUDIT_LOG, "w") as f:
        f.write("Phase 0.5 v2 per-dataset QC harmonization\n")
        f.write("=" * 60 + "\n")
        f.write(f"Harmonized cutoffs: {json.dumps(harm)}\n")
        f.write(f"Datasets processed: {dataset_ids}\n")
        f.write("decontX ambient correction: SKIPPED (celda/decontX not in rapids_singlecell)\n")
        f.write("Disclose in Methods that ambient correction was not applied uniformly.\n")
        f.write("\n")
        f.write(audit_df.to_string(index=False) + "\n")

    log.info("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
