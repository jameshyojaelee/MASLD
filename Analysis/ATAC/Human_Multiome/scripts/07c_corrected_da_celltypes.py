#!/usr/bin/env python3
"""07c_corrected_da_celltypes.py — Donor-level pseudobulk DA for the three
non-hepatocyte compartments shown in Fig S4g (Stellate, Macrophage/Kupffer,
Cholangiocyte).

WHY (2026-07-09)
----------------
Fig S4g panel a previously drew per-CELL SnapATAC2 `diff_test` volcanoes
(scatac_da_results.csv) — pseudoreplication: cells treated as independent
replicates when the experimental unit is the DONOR (18 donors: 13 MASLD [MASL+
MASH] vs 5 NORMAL). 07b already fixed this for hepatocytes; this script applies
the SAME donor-pseudobulk pipeline to the other three panel compartments so the
whole figure is computed at donor resolution.

It is a thin driver: it imports and reuses the validated helpers from
07b_corrected_da_hepatocytes.py (aggregation, tile->peak mapping, the pseudobulk
count export, and the OLS fallback). The only per-compartment inputs are (a) the
obs `cell_type` label(s) to subset and (b) the cell-type-specific peak BED. Each
compartment writes:
  * {label}_pseudobulk_counts.tsv.gz + {label}_pseudobulk_coldata.tsv
      -> edgeR-QLF input for 17c_edger_celltypes.R (the CANONICAL test)
  * scatac_da_corrected_{label}.csv  (the OLS-on-log-CPM arm, for parity with
      scatac_da_corrected_hep.csv; edgeR is the headline)

Hepatocytes are NOT re-run here — scatac_da_corrected_hep{,_edger}.csv already
exist and are canonical; keeping them untouched preserves their bit-for-bit
provenance.

Env: snapatac2 (needs ~/.local scipy for backed-anndata sparse indexing; do NOT
set PYTHONNOUSERSITE). Run via run_celltype_da.sbatch, not on the login node.
"""

import argparse
import importlib.util
import logging
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# utils_pseudobulk (imported by 07b) lives in this dir.
sys.path.insert(0, SCRIPT_DIR)

# 07b's module name starts with a digit -> load it by path.
_spec = importlib.util.spec_from_file_location(
    "da07b", os.path.join(SCRIPT_DIR, "07b_corrected_da_hepatocytes.py")
)
da07b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(da07b)

PEAK_DIR = os.path.join(
    SCRIPT_DIR, "..", "results", "label_transfer", "cell_type_peak_sets_v2"
)

# Compartments to add. `obs` = the obs cell_type value(s) to pool; `display` =
# the label written into the output `cell_type` column (matches Fig S4g ct_order
# exactly so the figure joins without remapping). Kupffer cells are liver-
# resident macrophages -> pooled with Macrophage into the panel's
# "Macrophage/Kupffer" compartment, using the Macrophages peak set.
COMPARTMENTS = [
    dict(label="stellate",      display="Stellate",
         obs=["Stellate_Cell"],                 bed="Fibroblasts_peaks.bed"),
    dict(label="macrophage",    display="Macrophage/Kupffer",
         obs=["Macrophage", "Kupffer_Cell"],    bed="Macrophages_peaks.bed"),
    dict(label="cholangiocyte", display="Cholangiocyte",
         obs=["Cholangiocyte"],                 bed="Cholangiocytes_peaks.bed"),
]


def extract_ct_matrix(adata, ct_mask):
    """Chunked backed read of the peak/tile matrix for the masked cells.

    Returns (X csr [n_ct_cell x n_tile], ct_indices_sorted). Mirrors 07b's
    chunked extraction so backed sparse indexing stays cheap.
    """
    ct_indices = np.where(ct_mask)[0]
    ct_indices_sorted = np.sort(ct_indices)
    chunk_size = 5000
    n_chunks = (len(ct_indices_sorted) + chunk_size - 1) // chunk_size
    X_chunks = []
    for ci, start in enumerate(range(0, len(ct_indices_sorted), chunk_size)):
        chunk_idx = ct_indices_sorted[start:start + chunk_size]
        X_chunk = adata.X[chunk_idx, :]
        if not sparse.issparse(X_chunk):
            X_chunk = sparse.csr_matrix(X_chunk)
        X_chunks.append(X_chunk)
        if ci % 5 == 0:
            log.info("    read chunk %d / %d", ci + 1, n_chunks)
    X = sparse.vstack(X_chunks, format="csr")
    return X, ct_indices_sorted


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True,
                    help="label-transferred h5ad (var_names = 500bp tiles)")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--cell-type-col", default="cell_type")
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--min-donor-detect", type=int, default=3)
    ap.add_argument("--min-total-reads", type=int, default=10)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    t_start = time.time()

    # --- donor -> condition(0/1) map (reuse 07b's exact MASLD/NORMAL sets) ---
    dmeta = pd.read_csv(da07b.DONOR_META_PATH, sep="\t")
    cond_upper = dmeta["condition"].astype(str).str.upper()
    donor_is_masld = {}
    for did, cu in zip(dmeta["donor_id"].astype(str), cond_upper):
        if cu in {"MASLD", "MASH", "MASL", "NASH", "NAFLD", "NAFL"}:
            donor_is_masld[did] = 1.0
        elif cu in {"NORMAL", "HEALTHY", "CONTROL"}:
            donor_is_masld[did] = 0.0
    log.info("Condition map: %d donors (MASLD=%d, NORMAL=%d)",
             len(donor_is_masld),
             sum(v == 1.0 for v in donor_is_masld.values()),
             sum(v == 0.0 for v in donor_is_masld.values()))

    # --- load h5ad backed once; the tile grid is shared across compartments ---
    log.info("Loading %s (backed)...", args.input)
    adata = ad.read_h5ad(args.input, backed="r")
    obs = adata.obs
    for col in (args.cell_type_col, args.donor_col):
        if col not in obs.columns:
            log.error("column %r not in obs (%s)", col, list(obs.columns))
            sys.exit(1)
    tile_names = list(adata.var_names)
    chrom_to_range = da07b.parse_chrom_sizes_from_tiles(tile_names)
    log.info("h5ad: %d cells x %d tiles", adata.n_obs, adata.n_vars)

    summary = []
    for comp in COMPARTMENTS:
        label, display, obs_vals, bed = (
            comp["label"], comp["display"], comp["obs"], comp["bed"]
        )
        log.info("")
        log.info("=== compartment %s (obs %s) ===", display, obs_vals)

        ct_mask = obs[args.cell_type_col].isin(obs_vals).values
        n_ct = int(ct_mask.sum())
        if n_ct == 0:
            log.warning("  no cells for %s; skipping", display)
            continue
        donor_array = obs[args.donor_col].astype(str).values[ct_mask]
        log.info("  %d cells; per-donor counts: %s", n_ct,
                 pd.Series(donor_array).value_counts().sort_index().to_dict())

        X_tiles, ct_idx_sorted = extract_ct_matrix(adata, ct_mask)
        # reorder donor_array to the sorted-index order used for extraction
        sort_perm = np.argsort(np.where(ct_mask)[0], kind="stable")
        donor_array_sorted = donor_array[sort_perm]

        # tiles -> called peaks for THIS compartment's peak set
        peak_bed = os.path.normpath(os.path.join(PEAK_DIR, bed))
        peaks = da07b.load_peaks(peak_bed)
        peak_to_tile_S = da07b.build_peak_to_tile_mapping(
            peaks, chrom_to_range, n_tiles=len(tile_names)
        )
        peak_names = [f"{c}:{s}-{e}" for (c, s, e) in peaks]
        log.info("  %d peaks from %s; map=%s",
                 len(peaks), bed, peak_to_tile_S.shape)

        da = da07b.run_pseudobulk_da(
            X_tiles, donor_array_sorted, peak_names, donor_is_masld,
            peak_to_tile_S,
            min_donor_detect=args.min_donor_detect,
            min_total_reads=args.min_total_reads,
            counts_export_dir=args.output_dir,
            counts_basename=f"{label}_pseudobulk",
        )
        da["cell_type"] = display
        out_csv = os.path.join(args.output_dir, f"scatac_da_corrected_{label}.csv")
        da.to_csv(out_csv, index=False)

        n_sig = int((da["adjusted p-value"] < 0.05).sum())
        log.info("  OLS arm: %d peaks tested, %d sig (padj<0.05) -> %s",
                 len(da), n_sig, out_csv)
        summary.append((display, len(da), n_sig))

    adata.file.close()
    log.info("")
    log.info("=== OLS-arm summary (edgeR via 17c is the headline) ===")
    for disp, n, s in summary:
        log.info("  %-20s %6d peaks, %d sig(padj<0.05)", disp, n, s)
    log.info("Total time: %.1f min", (time.time() - t_start) / 60)


if __name__ == "__main__":
    main()
