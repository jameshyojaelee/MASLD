#!/usr/bin/env python3
"""07b_corrected_da_hepatocytes.py — Corrected DA test for hepatocytes only.

PSEUDOREPLICATION FIX (2026-05-30, review findings F048/F150)
-------------------------------------------------------------
The previous version fit a per-PEAK logistic regression with ONE ROW PER CELL
(~52k hepatocytes from ~18 donors) and no donor random effect. That treats
individual cells as independent replicates when the true experimental unit is
the DONOR, so the effective n was wildly inflated (~41% of peaks "significant").

This rewrite collapses cells -> donors BEFORE testing, mirroring the
already-correct pattern in 14_stage_stratified_da.py:
  * pseudobulk per donor: SUM peak counts across that donor's hepatocytes
    (via utils_pseudobulk.aggregate_cells_to_donors with agg="sum").
  * per-donor library size; convert to log2(CPM+1).
  * per-peak OLS  log2(CPM+1) ~ condition_MASLD  across donors (n ~= 18),
    using the same vectorized OLS helper as 14 (fit_lm_vectorized).
The effective n is now the number of DONORS, not cells.

NOTE: this script does NOT re-call peaks. It reuses the peak matrix already
present in the label-transferred h5ad (var_names ARE peak coordinates), and
simply aggregates the hepatocyte rows to donors.

Per-cell QC covariates from the original model (log10 n_fragment, TSS
enrichment) are dropped: they are cell-level nuisance variables that do not
apply once cells are pooled into a donor pseudobulk. Only donor-level
predictors are admissible here; condition is the sole predictor.

Output format matches the original for compatibility with
08_annotate_peaks_for_l8.py — columns:
  'feature name', 'log2(fold_change)', 'p-value', 'adjusted p-value',
  'cell_type'  (plus a new 'n_donors' column).
"""

import argparse
import logging
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests

# Shared donor-level pseudobulk helper (same directory).
from utils_pseudobulk import aggregate_cells_to_donors

warnings.filterwarnings("ignore", category=RuntimeWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# Donor metadata: the 'donor_id' column (D01..D18) matches the h5ad obs
# 'donor_id'; condition is NORMAL / MASL / MASH. disease = {MASL, MASH};
# healthy = {NORMAL}. (donor_id_atac = MM_### is a different, non-matching key.)
DONOR_META_PATH = (
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv"
)

# The label-transferred h5ad var_names are 500bp genome TILES (~6M), NOT the
# called hepatocyte peaks. We map tiles -> peaks (mirroring 14_stage_stratified_da
# .py) before testing; otherwise DA runs on millions of tiles. (Fix 2026-05-30.)
HEP_PEAK_BED = (
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "Analysis/ATAC/Human_Multiome/results/label_transfer/"
    "cell_type_peak_sets_v2/Hepatocytes_peaks.bed"
)
TILE_SIZE = 500


def parse_chrom_sizes_from_tiles(var_names):
    """Map each chromosome to its [first, last] tile index (tiles are sorted)."""
    chrom_to_range = {}
    cur_chrom, cur_start = None, 0
    for i, nm in enumerate(var_names):
        chrom = nm.split(":", 1)[0]
        if chrom != cur_chrom:
            if cur_chrom is not None:
                chrom_to_range[cur_chrom] = (cur_start, i - 1)
            cur_chrom, cur_start = chrom, i
    chrom_to_range[cur_chrom] = (cur_start, len(var_names) - 1)
    return chrom_to_range


def load_peaks(peak_bed_path):
    peaks = []
    with open(peak_bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                peaks.append((parts[0], int(parts[1]), int(parts[2])))
    return peaks


def build_peak_to_tile_mapping(peaks, chrom_to_range, n_tiles):
    """Sparse (n_peak x n_tile) indicator mapping each peak to its overlapping tiles."""
    rows, cols = [], []
    for p_idx, (chrom, start, end) in enumerate(peaks):
        if chrom not in chrom_to_range:
            continue
        c0, c1 = chrom_to_range[chrom]
        gs = max(c0 + start // TILE_SIZE, c0)
        ge = min(c0 + (end - 1) // TILE_SIZE, c1)
        if ge < gs:
            continue
        for t in range(gs, ge + 1):
            rows.append(p_idx)
            cols.append(t)
    data = np.ones(len(rows), dtype=np.float32)
    return sparse.csr_matrix((data, (rows, cols)),
                             shape=(len(peaks), n_tiles), dtype=np.float32)


def fit_lm_vectorized(Y_log_cpm, x):
    """Vectorized OLS of log2(CPM+1) on x across donors.

    Copied from 14_stage_stratified_da.py (the canonical donor-level DA in this
    pipeline) to keep the two scripts consistent.

    Y_log_cpm: (n_donor, n_peak) log2(CPM+1)
    x: (n_donor,) covariate of interest (0/1 condition indicator here)

    Fits per-peak  Y[:,j] = beta0 + beta1 * x + eps  and returns
    (beta1, se, pval) arrays of shape (n_peak,). For tiny pseudobulk n this OLS
    on log-CPM is the standard limma-voom-style fast approximation (pydeseq2 is
    not available in the snapatac2 env).
    """
    n = Y_log_cpm.shape[0]
    n_peak = Y_log_cpm.shape[1]
    x = x.astype(np.float64)
    x_mean = x.mean()
    x_dev = x - x_mean
    SSx = (x_dev ** 2).sum()
    if SSx <= 0:
        return (np.full(n_peak, np.nan),
                np.full(n_peak, np.nan),
                np.full(n_peak, np.nan))

    # beta1_j = sum_i (x_i - xbar)(y_ij - ybar_j) / SSx
    Y_mean = Y_log_cpm.mean(axis=0)
    Y_dev = Y_log_cpm - Y_mean  # (n, n_peak)
    cov_xy = (x_dev[:, None] * Y_dev).sum(axis=0)
    beta1 = cov_xy / SSx
    beta0 = Y_mean - beta1 * x_mean

    # Residuals + SSE
    yhat = beta0[None, :] + np.outer(x, beta1)
    resid = Y_log_cpm - yhat
    SSE = (resid ** 2).sum(axis=0)
    df_resid = n - 2
    if df_resid <= 0:
        return (beta1,
                np.full(n_peak, np.nan),
                np.full(n_peak, np.nan))
    sigma2 = SSE / df_resid
    se = np.sqrt(sigma2 / SSx)

    # t-stat and two-sided p-value via t-distribution
    from scipy import stats as scipy_stats
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, beta1 / se, np.nan)
    pval = np.full(n_peak, np.nan)
    valid = np.isfinite(tstat)
    pval[valid] = 2.0 * scipy_stats.t.sf(np.abs(tstat[valid]), df=df_resid)
    return beta1, se, pval


def run_pseudobulk_da(X_hep, donor_array, peak_names, donor_is_masld,
                      peak_to_tile_S, min_donor_detect=3, min_total_reads=10,
                      counts_export_dir=None, counts_basename="hep_pseudobulk"):
    """Donor-level pseudobulk DA for hepatocyte peaks.

    Parameters
    ----------
    X_hep : (n_cell, n_peak) sparse/dense hepatocyte peak COUNT matrix.
    donor_array : (n_cell,) donor id per cell (matches donor metadata donor_id).
    peak_names : list of peak coordinate strings (h5ad var_names).
    donor_is_masld : dict donor_id -> 1.0 (MASLD) / 0.0 (NORMAL).
    min_donor_detect, min_total_reads : peak filter (mirrors 14's
        ">= min_total_reads total AND detected in >= min_donor_detect donors").
    counts_export_dir : if set, the FILTERED donor x peak pseudobulk COUNT
        matrix + per-donor coldata (condition 0/1 = MASLD, libsize) are written
        there as hep_pseudobulk_counts.tsv.gz + hep_pseudobulk_coldata.tsv for
        an external edgeR-QLF re-test (17_pseudobulk_de_retest.R). This is the
        Squair et al. 2021 pseudobulk + edgeR empirical-Bayes input; the OLS on
        log2(CPM+1) below is left UNCHANGED so scatac_da_corrected_hep.csv is
        bit-for-bit identical.

    Returns a DataFrame with the columns the original script wrote, plus
    'n_donors'. The effective n is the number of DONORS.
    """
    # 1. Aggregate cells -> per-donor TILE sums, then map tiles -> PEAK counts via
    #    peak_to_tile_S (n_peak x n_tile), mirroring 14_stage_stratified_da.py.
    donors, M_tiles, n_cells = aggregate_cells_to_donors(
        X_hep, donor_array, agg="sum", min_cells=1
    )
    M_counts = np.asarray(peak_to_tile_S @ M_tiles.T).T  # (n_donor, n_peak)
    log.info(
        "  Pseudobulk: %d donors x %d peaks (from %d tiles; cells/donor "
        "median=%.0f, range=[%d, %d])",
        len(donors), M_counts.shape[1], M_tiles.shape[1], float(np.median(n_cells)),
        int(n_cells.min()), int(n_cells.max()),
    )

    # 2. Restrict to donors with a known condition label.
    cond = np.array([donor_is_masld.get(d, np.nan) for d in donors], dtype=float)
    keep_donor = np.isfinite(cond)
    dropped = [d for d, k in zip(donors, keep_donor) if not k]
    if dropped:
        log.warning("  Dropping %d donor(s) lacking condition: %s",
                    len(dropped), dropped)
    donors = [d for d, k in zip(donors, keep_donor) if k]
    M_counts = M_counts[keep_donor, :]
    cond = cond[keep_donor]
    n_donors = len(donors)

    n_masld = int((cond == 1.0).sum())
    n_normal = int((cond == 0.0).sum())
    log.info("  Donors used: %d (MASLD=%d, NORMAL=%d)",
             n_donors, n_masld, n_normal)
    for d, c, nc in zip(donors, cond, n_cells[keep_donor]):
        log.info("    %-10s %-7s %d cells",
                 str(d), "MASLD" if c == 1.0 else "NORMAL", int(nc))
    if n_masld < 2 or n_normal < 2:
        raise RuntimeError(
            f"Need >=2 donors per group; got MASLD={n_masld}, NORMAL={n_normal}"
        )

    # 3. Library size per donor (total peak reads) + log2(CPM+1).
    libsize = M_counts.sum(axis=1)  # (n_donor,)
    log.info("  Donor libsizes: median=%.0f, range=[%.0f, %.0f]",
             float(np.median(libsize)), float(libsize.min()),
             float(libsize.max()))

    # 4. Peak filter: total reads AND detected in >= min_donor_detect donors
    #    (mirrors 14_stage_stratified_da.run_da).
    peak_total = M_counts.sum(axis=0)
    n_donor_detect = (M_counts > 0).sum(axis=0)
    keep_peak = (peak_total >= min_total_reads) & (n_donor_detect >= min_donor_detect)
    n_keep = int(keep_peak.sum())
    log.info(
        "  Peaks passing filter (>=%d reads & detected in >=%d donors): "
        "%d / %d (%.1f%%)",
        min_total_reads, min_donor_detect, n_keep, len(keep_peak),
        100.0 * n_keep / max(len(keep_peak), 1),
    )
    if n_keep == 0:
        raise RuntimeError("No peaks pass the donor-level filter")

    peak_idx = np.where(keep_peak)[0]
    counts_kept = M_counts[:, peak_idx]
    peak_names_kept = [peak_names[i] for i in peak_idx]

    # 4b. Export the FILTERED donor x peak pseudobulk COUNT matrix + coldata for
    #     the external edgeR-QLF empirical-Bayes re-test (Squair et al. 2021).
    #     This is purely additive: it does not alter counts_kept / Y / the OLS.
    #     Counts are rounded to integers (tile->peak mapping can give fractional
    #     sums when a peak spans tile boundaries; edgeR/DGEList wants integers).
    if counts_export_dir is not None:
        os.makedirs(counts_export_dir, exist_ok=True)
        counts_int = np.rint(counts_kept).astype(np.int64)  # (n_donor, n_peak)
        counts_df = pd.DataFrame(
            counts_int, index=list(donors), columns=peak_names_kept
        )
        counts_df.index.name = "donor_id"
        counts_path = os.path.join(counts_export_dir, f"{counts_basename}_counts.tsv.gz")
        counts_df.to_csv(counts_path, sep="\t", compression="gzip")

        coldata_df = pd.DataFrame({
            "donor_id": list(donors),
            "condition": cond.astype(int),  # 1 = MASLD, 0 = NORMAL
            "libsize": libsize.astype(np.int64),
        })
        coldata_path = os.path.join(counts_export_dir, f"{counts_basename}_coldata.tsv")
        coldata_df.to_csv(coldata_path, sep="\t", index=False)
        log.info(
            "  Exported pseudobulk COUNTS for edgeR: %s (%d donors x %d peaks) + %s",
            counts_path, counts_int.shape[0], counts_int.shape[1], coldata_path,
        )

    # 5. log2(CPM+1) using donor libsize.
    cpm = counts_kept / libsize[:, None] * 1e6
    Y = np.log2(cpm + 1.0)

    # 6. Per-peak OLS  log2(CPM+1) ~ condition_MASLD  (n = n_donors).
    log.info("  Fitting vectorized OLS on log2(CPM+1) over %d donors "
             "(n_peak=%d)...", n_donors, n_keep)
    t0 = time.time()
    coefs, ses, pvals = fit_lm_vectorized(Y, cond)
    log.info("  Done fitting in %.2f sec", time.time() - t0)

    # 7. BH FDR among valid p-values.
    padj = np.full(n_keep, np.nan)
    valid = np.isfinite(pvals)
    if valid.sum() > 0:
        _, padj_valid, _, _ = multipletests(pvals[valid], method="fdr_bh")
        padj[valid] = padj_valid

    # coefs are already in log2(CPM+1) units => directly log2 fold-change.
    df = pd.DataFrame({
        "feature name": peak_names_kept,
        "log2(fold_change)": coefs,
        "p-value": pvals,
        "adjusted p-value": padj,
        "se": ses,
        "n_donors": n_donors,
    })
    # Drop peaks with undefined p-value (no variation), matching old behaviour.
    df = df.dropna(subset=["p-value"]).reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Corrected hepatocyte DA via donor pseudobulk (no pseudoreplication)"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Label-transferred h5ad (snapatac2_label_transferred.h5ad); "
             "var_names are peak coordinates",
    )
    parser.add_argument(
        "--per-donor-dir",
        required=False,
        default=None,
        help="(Unused since 2026-05-30) Previously held per-cell QC h5ads; "
             "per-cell covariates are dropped under donor pseudobulk. Kept for "
             "CLI backward-compatibility.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--min-detection",
        type=float,
        default=0.05,
        help="(Deprecated) Per-cell detection-rate prefilter from the logistic "
             "model; ignored under donor pseudobulk. Use --min-donor-detect.",
    )
    parser.add_argument(
        "--min-donor-detect",
        type=int,
        default=3,
        help="Keep peaks detected (>=1 read) in at least this many donors "
             "(mirrors 14_stage_stratified_da; default: 3)",
    )
    parser.add_argument(
        "--min-total-reads",
        type=int,
        default=10,
        help="Keep peaks with at least this many total pseudobulk reads "
             "(mirrors 14_stage_stratified_da; default: 10)",
    )
    parser.add_argument(
        "--cell-type-col",
        default="cell_type",
        help="Column for cell type (default: cell_type)",
    )
    parser.add_argument(
        "--donor-col",
        default="donor_id",
        help="obs column with donor id matching donor metadata (default: donor_id)",
    )
    args = parser.parse_args()

    if args.per_donor_dir is not None:
        log.info("Note: --per-donor-dir is ignored (per-cell QC covariates "
                 "dropped under donor pseudobulk).")
    if abs(args.min_detection - 0.05) > 1e-9:
        log.info("Note: --min-detection is deprecated/ignored under donor "
                 "pseudobulk; using --min-donor-detect=%d.", args.min_donor_detect)

    os.makedirs(args.output_dir, exist_ok=True)
    t_start = time.time()

    # --- Step 1: Load donor metadata -> condition map ---
    log.info("Loading donor metadata: %s", DONOR_META_PATH)
    if not os.path.exists(DONOR_META_PATH):
        log.error("Donor metadata not found: %s", DONOR_META_PATH)
        sys.exit(1)
    dmeta = pd.read_csv(DONOR_META_PATH, sep="\t")
    # JOIN KEY: the h5ad obs 'donor_id' holds D## ids, which match the metadata
    # 'donor_id' column (D01..D18) -- NOT 'donor_id_atac' (= MM_454...). Keying on
    # donor_id_atac would match zero donors and drop the whole cohort. (Fix
    # 2026-05-30; the earlier draft mistakenly used donor_id_atac.)
    if "donor_id" not in dmeta.columns or "condition" not in dmeta.columns:
        log.error("Donor metadata missing required columns "
                  "(donor_id, condition); have %s", list(dmeta.columns))
        sys.exit(1)
    cond_upper = dmeta["condition"].astype(str).str.upper()
    masld_set = {"MASLD", "MASH", "MASL", "NASH", "NAFLD", "NAFL"}
    normal_set = {"NORMAL", "HEALTHY", "CONTROL"}
    donor_is_masld = {}
    for did, cu in zip(dmeta["donor_id"].astype(str), cond_upper):
        if cu in masld_set:
            donor_is_masld[did] = 1.0
        elif cu in normal_set:
            donor_is_masld[did] = 0.0
    log.info("  Condition map: %d donors (MASLD=%d, NORMAL=%d)",
             len(donor_is_masld),
             sum(v == 1.0 for v in donor_is_masld.values()),
             sum(v == 0.0 for v in donor_is_masld.values()))

    # --- Step 2: Load h5ad (backed) + subset to hepatocytes ---
    log.info("Loading %s (backed)...", args.input)
    adata = ad.read_h5ad(args.input, backed="r")
    log.info("Shape: %d x %d", adata.n_obs, adata.n_vars)

    obs = adata.obs
    if args.cell_type_col not in obs.columns:
        log.error("cell-type column %r not in obs (%s)",
                  args.cell_type_col, list(obs.columns))
        sys.exit(1)
    if args.donor_col not in obs.columns:
        log.error("donor column %r not in obs (%s)",
                  args.donor_col, list(obs.columns))
        sys.exit(1)

    all_names = list(adata.obs_names)
    hep_mask = (obs[args.cell_type_col] == "Hepatocyte").values
    hep_indices = np.where(hep_mask)[0]
    log.info("Hepatocytes: %d cells", len(hep_indices))
    if len(hep_indices) == 0:
        log.error("No hepatocytes found in column %r", args.cell_type_col)
        sys.exit(1)

    donor_array = obs[args.donor_col].astype(str).values[hep_indices]
    log.info("  Per-donor hepatocyte counts: %s",
             pd.Series(donor_array).value_counts().sort_index().to_dict())

    # Export a per-cell barcode -> donor_id -> condition map for the SCENIC
    # regulon-activity edgeR/limma re-test (17_pseudobulk_de_retest.R, arm b).
    # The saved regulon_activity_scores.csv is keyed by bare cell barcode with
    # NO donor column, and recovering donor labels otherwise requires reading
    # the 3.5 GB label-transferred h5ad obs in Python (not feasible cleanly in
    # R). We are already in that obs here, so emit the small map (barcodes match
    # the regulon CSV index exactly; verified). Covers all cells so the R join
    # is robust regardless of which cell types the regulon matrix contains.
    cellmap_path = os.path.join(args.output_dir, "cell_donor_condition_map.tsv.gz")
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        cell_map = pd.DataFrame({
            "barcode": list(adata.obs_names),
            "donor_id": obs[args.donor_col].astype(str).values,
            "condition": obs["condition"].astype(str).values
            if "condition" in obs.columns else "NA",
            "cell_type": obs[args.cell_type_col].astype(str).values,
        })
        cell_map.to_csv(cellmap_path, sep="\t", index=False, compression="gzip")
        log.info("  Exported cell->donor->condition map: %s (%d cells)",
                 cellmap_path, len(cell_map))
    except Exception as e:  # never let the map export break the DA run
        log.warning("  Could not write cell->donor map (%s): %s",
                    cellmap_path, e)

    # --- Step 3: Extract hepatocyte peak matrix (chunked, backed read) ---
    log.info("Extracting hepatocyte peak matrix (chunked)...")
    chunk_size = 5000
    hep_indices_sorted = np.sort(hep_indices)
    X_chunks = []
    n_chunks = (len(hep_indices_sorted) + chunk_size - 1) // chunk_size
    for ci, start in enumerate(range(0, len(hep_indices_sorted), chunk_size)):
        chunk_idx = hep_indices_sorted[start:start + chunk_size]
        X_chunk = adata.X[chunk_idx, :]
        if not sparse.issparse(X_chunk):
            X_chunk = sparse.csr_matrix(X_chunk)
        X_chunks.append(X_chunk)
        if ci % 5 == 0:
            log.info("  Read chunk %d / %d", ci + 1, n_chunks)
    X_hep = sparse.vstack(X_chunks, format="csr")
    del X_chunks

    # donor_array is in hep_indices order; reorder it to hep_indices_sorted.
    sort_perm = np.argsort(hep_indices, kind="stable")
    donor_array_sorted = donor_array[sort_perm]

    tile_names = list(adata.var_names)
    adata.file.close()
    log.info("  Hepatocyte tile matrix: %d cells x %d tiles",
             X_hep.shape[0], X_hep.shape[1])

    # Map 500bp tiles -> called hepatocyte peaks (var_names are tiles, not peaks).
    peaks = load_peaks(HEP_PEAK_BED)
    chrom_to_range = parse_chrom_sizes_from_tiles(tile_names)
    peak_to_tile_S = build_peak_to_tile_mapping(
        peaks, chrom_to_range, n_tiles=len(tile_names)
    )
    peak_names = [f"{c}:{s}-{e}" for (c, s, e) in peaks]
    log.info("  Loaded %d hepatocyte peaks; peak-to-tile map = %s",
             len(peaks), peak_to_tile_S.shape)

    # Pseudobulk COUNTS for the external edgeR re-test go next to the corrected
    # DA output (--output-dir is results/snapatac2), so the canonical paths are
    # results/snapatac2/hep_pseudobulk_{counts.tsv.gz,coldata.tsv}, alongside
    # scatac_da_corrected_hep.csv. 17_pseudobulk_de_retest.R reads them there.
    counts_export_dir = args.output_dir

    # --- Step 4: Donor-level pseudobulk DA (peak-level) ---
    da_results = run_pseudobulk_da(
        X_hep, donor_array_sorted, peak_names, donor_is_masld, peak_to_tile_S,
        min_donor_detect=args.min_donor_detect,
        min_total_reads=args.min_total_reads,
        counts_export_dir=counts_export_dir,
    )
    da_results["cell_type"] = "Hepatocyte"

    # --- Step 5: Save (preserve original output path + columns) ---
    out_path = os.path.join(args.output_dir, "scatac_da_corrected_hep.csv")
    da_results.to_csv(out_path, index=False)
    log.info("Saved %d DA results -> %s", len(da_results), out_path)

    # Summary statistics
    n_donors_used = int(da_results["n_donors"].iloc[0]) if len(da_results) else 0
    n_sig = int((da_results["adjusted p-value"] < 0.05).sum())
    n_sig_lfc = int(
        ((da_results["adjusted p-value"] < 0.05) &
         (da_results["log2(fold_change)"].abs() > 0.25)).sum()
    )
    up = int(
        ((da_results["adjusted p-value"] < 0.05) &
         (da_results["log2(fold_change)"] > 0.25)).sum()
    )
    down = int(
        ((da_results["adjusted p-value"] < 0.05) &
         (da_results["log2(fold_change)"] < -0.25)).sum()
    )
    log.info("  Effective n = %d donors", n_donors_used)
    log.info("  Significant (padj<0.05): %d", n_sig)
    log.info("  Significant + |logFC|>0.25: %d (Up=%d, Down=%d)", n_sig_lfc, up, down)

    lfc = da_results["log2(fold_change)"]
    if len(lfc):
        log.info(
            "  logFC: median=%.3f, mean=%.3f, [%.3f, %.3f]",
            lfc.median(), lfc.mean(), lfc.min(), lfc.max(),
        )

    # --- Step 6: Compare with original per-cell logistic results (if present) ---
    orig_path = os.path.join(
        os.path.dirname(args.output_dir), "scatac_da_results.csv"
    )
    if os.path.exists(orig_path):
        log.info("Comparing with original (per-cell) DA results...")
        orig = pd.read_csv(orig_path)
        if "cell_type" in orig.columns:
            orig = orig[orig["cell_type"] == "Hepatocyte"].copy()
        orig = orig.rename(columns={
            "log2(fold_change)": "logFC_orig",
            "adjusted p-value": "padj_orig",
        })
        if {"feature name", "logFC_orig"}.issubset(orig.columns):
            merged = da_results.merge(
                orig[["feature name", "logFC_orig", "padj_orig"]],
                on="feature name",
                how="inner",
            )
            if len(merged) > 10:
                from scipy.stats import spearmanr
                rho, p = spearmanr(
                    merged["log2(fold_change)"], merged["logFC_orig"]
                )
                direction = (
                    np.sign(merged["log2(fold_change)"])
                    == np.sign(merged["logFC_orig"])
                ).mean()
                log.info(
                    "  Overlap: %d peaks, Spearman rho=%.3f (p=%.2e), "
                    "direction concordance=%.1f%%",
                    len(merged), rho, p, 100 * direction,
                )

    log.info("Total time: %.1f min", (time.time() - t_start) / 60)


if __name__ == "__main__":
    main()
