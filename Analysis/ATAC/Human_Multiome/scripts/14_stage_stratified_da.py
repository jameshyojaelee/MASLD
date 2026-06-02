#!/usr/bin/env python3
"""14_stage_stratified_da.py -- Stage-stratified DA per cell type (ATAC B3a).

Computes differential accessibility across the F-stage trajectory using
pseudobulk count aggregation per (donor x cell_type), with a 3-bin ordinal
bucket (F0 / F2+F3 / F4) chosen because F1=0 and F2=1 in this 18-donor
cohort, so the fine-grained 5-stage axis is not feasible.

Models per peak (OLS / limma-voom-style linear model on log2(CPM+1)
pseudobulk, per-donor; NOT a negative binomial GLM -- pydeseq2 is not in the
snapatac2 env, so we substitute the standard small-n log-CPM linear model):
  * Primary ordinal: log2(CPM+1) ~ F_bin (numeric 0,1,2), CPM normalised by
    per-donor library size.
  * Binary transitions: F0 vs F2+F3, F2+F3 vs F4, F0 vs F4 -- same linear
    model with a binary indicator.

Pseudobulk strategy:
  - Load full label-transferred h5ad (in-memory CSR, ~3GB).
  - For each cell type (CT), subset cells and load that CT's peak BED
    (`cell_type_peak_sets_v2/{ATAC_label}_peaks.bed`).
  - Tile-to-peak aggregation: 500bp fixed-step tiles intersect peaks by
    `[start//500, end//500]` index range; we sum tile counts per peak per
    donor. Output: peak-count matrix of shape (n_donor, n_peak).
  - Filter peaks: keep those with >= 10 reads total and >= 1 read in >= 3
    donors.
  - Fit the OLS / limma-voom-style log2(CPM+1) linear model per peak;
    FDR-correct with BH.
  - Annotate to nearest gene via GENCODE v49 GTF (logic mirroring
    08_annotate_peaks_for_l8.py).

Outputs (results/stage_da/):
  - da_stage_ordinal_{CT}.csv
  - da_F0_vs_F2F3_{CT}.csv
  - da_F2F3_vs_F4_{CT}.csv
  - da_F0_vs_F4_{CT}.csv

CAVEAT (measurement error in the F-stage axis):
  The per-donor F-stage used here (`F_stage_augmented`) is largely scVI-INFERRED
  (column `F_stage_source` is `documented` vs `scvi_predicted`; inference LOOCV
  QWK ~0.74-0.76). The DA models treat F-stage as a FIXED, error-free predictor.
  Regression-style measurement error in the predictor attenuates the estimated
  stage effect toward zero, so these stage contrasts are a LOWER BOUND on the
  true accessibility-vs-stage association. The per-run count of documented vs
  inferred donors is logged and written into every output CSV
  (`n_donors_fstage_documented` / `n_donors_fstage_inferred`) for transparency.
  The model is intentionally NOT changed.

Environment: snapatac2 env (PYTHONNOUSERSITE=1).
"""

import os
import sys
import gzip
import time
import logging
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
import anndata as ad
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
PROJECT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATAC_DIR = f"{PROJECT}/Analysis/ATAC/Human_Multiome"
H5AD_PATH = f"{ATAC_DIR}/results/label_transfer/snapatac2_label_transferred.h5ad"
PEAK_DIR = f"{ATAC_DIR}/results/label_transfer/cell_type_peak_sets_v2"
DONOR_META_PATH = f"{ATAC_DIR}/metadata/donor_metadata_curated.tsv"
GTF_PATH = (
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/"
    "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
)
OUT_DIR = f"{ATAC_DIR}/results/stage_da"
os.makedirs(OUT_DIR, exist_ok=True)

TILE_SIZE = 500

# Map of cell_type (h5ad) -> ATAC peak BED basename + group name
CT_GROUPS = {
    "Hep": {
        "h5ad_labels": ["Hepatocyte"],
        "peak_bed": f"{PEAK_DIR}/Hepatocytes_peaks.bed",
    },
    "Mac": {
        "h5ad_labels": ["Macrophage", "Kupffer_Cell"],
        "peak_bed": f"{PEAK_DIR}/Macrophages_peaks.bed",
    },
    "Fib": {
        "h5ad_labels": ["Stellate_Cell"],
        "peak_bed": f"{PEAK_DIR}/Fibroblasts_peaks.bed",
    },
    "Endo": {
        "h5ad_labels": ["Endothelial", "LSEC"],
        "peak_bed": f"{PEAK_DIR}/Endothelial_cells_peaks.bed",
    },
    "Chol": {
        "h5ad_labels": ["Cholangiocyte"],
        "peak_bed": f"{PEAK_DIR}/Cholangiocytes_peaks.bed",
    },
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_donor_metadata():
    """Load and bucket F-stage.

    Also carries `F_stage_source` (documented vs scvi_predicted) so the
    fraction of donors whose stage is INFERRED rather than measured can be
    reported alongside results (the model treats F-stage as fixed; inferred
    stage induces measurement-error attenuation -- see header caveat).
    """
    meta = pd.read_csv(DONOR_META_PATH, sep="\t")
    meta["F_stage_augmented"] = pd.to_numeric(
        meta["F_stage_augmented"], errors="coerce"
    )
    # F_stage_source: 'documented' = measured Kleiner/SAF stage; anything else
    # (e.g. 'scvi_predicted') is treated as inferred.
    if "F_stage_source" in meta.columns:
        meta["F_stage_source"] = meta["F_stage_source"].astype(str)
    else:
        log.warning(
            "F_stage_source column absent from %s; marking all donors 'unknown' "
            "(cannot distinguish documented vs inferred F-stage).",
            DONOR_META_PATH,
        )
        meta["F_stage_source"] = "unknown"

    def bucket(f):
        if pd.isna(f):
            return np.nan
        if f == 0:
            return 0
        if f in (2, 3):
            return 1
        if f == 4:
            return 2
        # F1 = 0 donors in this cohort; treat as missing for safety
        return np.nan

    meta["F_bin"] = meta["F_stage_augmented"].apply(bucket)
    meta = meta.dropna(subset=["F_bin"]).copy()
    meta["F_bin"] = meta["F_bin"].astype(int)
    log.info(
        "Donor F_bin distribution: %s",
        meta["F_bin"].value_counts().sort_index().to_dict(),
    )
    log.info(
        "Donor F_stage_source distribution: %s",
        meta["F_stage_source"].value_counts(dropna=False).to_dict(),
    )
    return meta[
        ["donor_id", "F_stage_augmented", "F_bin", "F_stage_source"]
    ].reset_index(drop=True)


def _count_fstage_source(source_values):
    """Split a series/array of F_stage_source labels into (n_documented,
    n_inferred). Anything not exactly 'documented' (e.g. 'scvi_predicted',
    'unknown') counts as inferred, the conservative choice for the caveat."""
    s = pd.Series(source_values).astype(str).str.strip().str.lower()
    n_documented = int((s == "documented").sum())
    n_inferred = int(len(s) - n_documented)
    return n_documented, n_inferred


def parse_chrom_sizes_from_tiles(var_names):
    """Get max tile coordinate per chrom from the tile index names."""
    log.info("Parsing tile index from var_names (%d tiles)...", len(var_names))
    # Tiles look like 'chr1:0-500'. We need chrom and start.
    # Build mapping: chrom -> (first_idx, last_idx, max_end).
    chrom_to_range = {}
    cur_chrom = None
    cur_start = 0
    for i, nm in enumerate(var_names):
        chrom, rest = nm.split(":", 1)
        if chrom != cur_chrom:
            if cur_chrom is not None:
                # Previous chrom ends at i-1
                chrom_to_range[cur_chrom] = (cur_start, i - 1)
            cur_chrom = chrom
            cur_start = i
    # Close the last chrom
    chrom_to_range[cur_chrom] = (cur_start, len(var_names) - 1)
    return chrom_to_range


def load_peaks(peak_bed_path):
    """Load peaks BED file. Returns list of (chrom, start, end)."""
    peaks = []
    with open(peak_bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            peaks.append((parts[0], int(parts[1]), int(parts[2])))
    log.info("  Loaded %d peaks from %s", len(peaks), peak_bed_path)
    return peaks


def build_peak_to_tile_mapping(peaks, chrom_to_range, n_tiles):
    """Build sparse matrix S (n_peak x n_tile) where S[p,t] = 1 if tile t
    falls within peak p. Then peak_counts = cell_tile_counts @ S.T."""
    log.info("  Building peak <-> tile mapping...")
    rows, cols = [], []
    n_peak = len(peaks)
    n_skipped = 0
    for p_idx, (chrom, start, end) in enumerate(peaks):
        if chrom not in chrom_to_range:
            n_skipped += 1
            continue
        chrom_start_idx, chrom_end_idx = chrom_to_range[chrom]
        # Tile t covers [t*500, (t+1)*500); tile index relative to chrom 0 is
        # t = tile_idx - chrom_start_idx, with tile_start = t*500.
        # We want tiles overlapping [start, end].
        tile_start = start // TILE_SIZE
        tile_end = (end - 1) // TILE_SIZE  # inclusive last tile
        # Global tile indices:
        gs = chrom_start_idx + tile_start
        ge = chrom_start_idx + tile_end
        # Clip to chrom range
        gs = max(gs, chrom_start_idx)
        ge = min(ge, chrom_end_idx)
        if ge < gs:
            continue
        for t in range(gs, ge + 1):
            rows.append(p_idx)
            cols.append(t)
    if n_skipped:
        log.warning("  %d peaks on chroms not in tile index (skipped)",
                    n_skipped)
    data = np.ones(len(rows), dtype=np.float32)
    S = sp.csr_matrix(
        (data, (rows, cols)), shape=(n_peak, n_tiles), dtype=np.float32
    )
    log.info("  Peak->tile mapping: %d entries, avg %.1f tiles/peak",
             S.nnz, S.nnz / n_peak)
    return S


def pseudobulk_by_donor(X_cells, donor_array, peak_to_tile_S):
    """Build per-donor pseudobulk peak count matrix.

    X_cells: csr_matrix (n_cell x n_tile)
    donor_array: array of length n_cell with donor IDs
    peak_to_tile_S: csr_matrix (n_peak x n_tile) of 0/1

    Returns: DataFrame index=donor, columns=peak_idx, values=counts
             and library size per donor (total tile counts).
    """
    donors = sorted(np.unique(donor_array))
    log.info("  Pseudobulking %d cells across %d donors...",
             X_cells.shape[0], len(donors))
    rows = []
    libsizes = []
    for d in donors:
        mask = donor_array == d
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        # Sum tile counts across cells for this donor (1 x n_tile)
        tile_sum = np.asarray(X_cells[mask].sum(axis=0)).ravel().astype(
            np.float32
        )
        libsize = float(tile_sum.sum())
        # Map to peaks: peak_counts = S @ tile_sum (n_peak,)
        peak_counts = peak_to_tile_S @ tile_sum  # (n_peak,)
        rows.append((d, n_cells, libsize, peak_counts))
        libsizes.append(libsize)

    donors_used = [r[0] for r in rows]
    n_cells_per = [r[1] for r in rows]
    libs = [r[2] for r in rows]
    counts = np.vstack([r[3] for r in rows])  # (n_donor, n_peak)

    df_counts = pd.DataFrame(
        counts,
        index=donors_used,
    )
    df_counts.index.name = "donor_id"
    df_meta = pd.DataFrame({
        "n_cells": n_cells_per,
        "libsize": libs,
    }, index=donors_used)
    log.info("  Pseudobulk libsizes: median=%.0f, range=[%.0f, %.0f]",
             np.median(libs), min(libs), max(libs))
    return df_counts, df_meta


def fit_lm_vectorized(Y_log_cpm, x):
    """Vectorized linear regression of log2(CPM+1) on x.

    Y_log_cpm: (n_donor, n_peak) log2(CPM+1)
    x: (n_donor,) covariate of interest (numeric for ordinal, 0/1 for binary)

    Fits per-peak: Y[:,j] = beta0 + beta1 * x + eps
    Returns (beta1, se, pval) arrays of shape (n_peak,).

    For n=5-13 small samples and pseudobulk, OLS on log-CPM is a standard
    fast approximation (limma-voom style for tiny n; we substitute for
    DESeq2/pyDESeq2 since pydeseq2 is not in the snapatac2 env).
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


def run_da(counts_df, meta_df, donor_meta_df, test_type, peaks_coords):
    """Run DA test. test_type is one of: 'ordinal', 'F0_vs_F2F3',
    'F2F3_vs_F4', 'F0_vs_F4'."""
    # Align donors
    donor_meta_idx = donor_meta_df.set_index("donor_id")
    common = [d for d in counts_df.index if d in donor_meta_idx.index]
    if not common:
        log.warning("  No donor overlap for %s", test_type)
        return None
    counts_df = counts_df.loc[common]
    meta_df = meta_df.loc[common]
    F_bin = donor_meta_idx.loc[common, "F_bin"].values.astype(int)

    # Use raw F_stage_augmented for fine-grained contrasts (drops the single
    # F2 donor automatically since every "fine" contrast specifies F0/F3/F4)
    F_fine = donor_meta_idx.loc[common, "F_stage_augmented"].values.astype(float)

    # Build group selector
    if test_type == "ordinal":
        keep = np.ones(len(F_bin), dtype=bool)
        x = F_bin.astype(float)
    elif test_type == "F0_vs_F2F3":
        keep = (F_bin == 0) | (F_bin == 1)
        x = (F_bin == 1).astype(float)
    elif test_type == "F2F3_vs_F4":
        keep = (F_bin == 1) | (F_bin == 2)
        x = (F_bin == 2).astype(float)
    elif test_type == "F0_vs_F4":
        keep = (F_bin == 0) | (F_bin == 2)
        x = (F_bin == 2).astype(float)
    # ── fine-grained contrasts (drop single F2 donor) ──────────────────────
    elif test_type == "F0_vs_F3":
        keep = (F_fine == 0) | (F_fine == 3)
        x = (F_fine == 3).astype(float)
    elif test_type == "F3_vs_F4":
        keep = (F_fine == 3) | (F_fine == 4)
        x = (F_fine == 4).astype(float)
    elif test_type == "ordinal_F0_F3_F4":
        keep = (F_fine == 0) | (F_fine == 3) | (F_fine == 4)
        x = F_fine.copy()
    else:
        raise ValueError(f"Unknown test_type {test_type}")

    counts_df = counts_df.iloc[keep]
    libsize = meta_df.iloc[keep]["libsize"].values
    n_donors = len(counts_df)
    x = x[keep]

    # F-stage provenance for the donors actually entering THIS contrast.
    # `common` is the donor order aligned above; `keep` selects the subset
    # used here. F-stage is largely scVI-inferred (measurement-error caveat in
    # the header) -- carry the documented/inferred split into log + output.
    donors_used = [d for d, k in zip(common, keep) if k]
    if "F_stage_source" in donor_meta_idx.columns:
        n_fstage_documented, n_fstage_inferred = _count_fstage_source(
            donor_meta_idx.loc[donors_used, "F_stage_source"].values
        )
    else:
        n_fstage_documented, n_fstage_inferred = 0, n_donors

    log.info(
        "  [%s] %d donors used (F-stage: %d documented / %d inferred); "
        "group distribution: %s",
        test_type, n_donors, n_fstage_documented, n_fstage_inferred,
        np.unique(x, return_counts=True),
    )
    if n_donors < 4:
        log.warning("  [%s] only %d donors, skipping", test_type, n_donors)
        return None

    # Filter peaks: total reads >= 10 AND >= 1 read in >= 3 donors
    counts_arr = counts_df.values.astype(np.float32)  # (n_donor, n_peak)
    peak_total = counts_arr.sum(axis=0)
    n_donor_detect = (counts_arr > 0).sum(axis=0)
    keep_peak = (peak_total >= 10) & (n_donor_detect >= 3)
    n_keep = int(keep_peak.sum())
    log.info("  [%s] peaks passing filter: %d / %d",
             test_type, n_keep, len(keep_peak))
    if n_keep == 0:
        return None

    peak_indices = np.where(keep_peak)[0]
    counts_kept = counts_arr[:, peak_indices]  # (n_donor, n_keep)

    # Convert to log2(CPM+1) using donor libsize
    cpm = counts_kept / libsize[:, None] * 1e6
    Y = np.log2(cpm + 1.0)

    log.info("  [%s] Fitting vectorized OLS on log2(CPM+1) (n_peak=%d)...",
             test_type, n_keep)
    t0 = time.time()
    coefs, ses, pvals = fit_lm_vectorized(Y, x)
    log.info("  [%s] Done fitting in %.1f sec",
             test_type, time.time() - t0)

    # BH FDR within this (cell_type x test) family. A second, across-family
    # global BH (`padj_global`) is added in main() by pooling p-values across
    # all ~35 contrasts/cell-types of this run (the contrasts are overlapping/
    # nested: ordinal + F0vsF2F3 + F2F3vsF4 + F0vsF4 + fine contrasts, each x
    # several cell types). Both columns are reported: `padj` = per-family
    # (per cell_type x test), `padj_global` = pooled across the whole run.
    valid = np.isfinite(pvals)
    padj = np.full(n_keep, np.nan)
    if valid.sum() > 0:
        _, padj_valid, _, _ = multipletests(pvals[valid], method="fdr_bh")
        padj[valid] = padj_valid

    # Build result frame -- coefs are in log2(CPM+1) units => directly log2FC
    rec = []
    for j, p_idx in enumerate(peak_indices):
        chrom, start, end = peaks_coords[p_idx]
        rec.append({
            "peak_id": f"{chrom}:{start}-{end}",
            "chrom": chrom,
            "start": start,
            "end": end,
            "log2FC": coefs[j],
            "se": ses[j],
            "pvalue": pvals[j],
            "padj": padj[j],
            "test": test_type,
            "n_donors": n_donors,
            "n_donors_fstage_documented": n_fstage_documented,
            "n_donors_fstage_inferred": n_fstage_inferred,
        })
    out_df = pd.DataFrame(rec).sort_values("pvalue", na_position="last")
    n_sig = int((out_df["padj"] < 0.05).sum())
    log.info("  [%s] Significant peaks (padj<0.05): %d", test_type, n_sig)
    return out_df


def parse_tss_from_gtf(gtf_path):
    """Parse gene TSS from GENCODE GTF. Returns DataFrame and chrom->sorted
    (tss_array, gene_array) lookup."""
    log.info("Parsing TSS from GTF: %s", gtf_path)
    records = []
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            if not chrom.startswith("chr"):
                continue
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            attrs = fields[8]
            gene_name = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                    break
            if gene_name is None:
                continue
            tss = start if strand == "+" else end
            records.append((chrom, tss, gene_name))
    df = pd.DataFrame(records, columns=["chrom", "tss", "gene_name"])
    df = df.drop_duplicates(subset="gene_name", keep="first")
    log.info("  Parsed %d TSS", len(df))

    chrom_tss = {}
    for chrom, sub in df.groupby("chrom"):
        sub_sorted = sub.sort_values("tss")
        chrom_tss[chrom] = (
            sub_sorted["tss"].values.astype(np.int64),
            sub_sorted["gene_name"].values,
        )
    return chrom_tss


def annotate_peaks(df, chrom_tss):
    """Annotate each peak with nearest gene."""
    gene_syms = []
    dists = []
    for _, row in df.iterrows():
        chrom = row["chrom"]
        mid = (row["start"] + row["end"]) // 2
        if chrom not in chrom_tss:
            gene_syms.append(None)
            dists.append(np.nan)
            continue
        tss_arr, gene_arr = chrom_tss[chrom]
        idx = np.searchsorted(tss_arr, mid)
        candidates = []
        if idx > 0:
            candidates.append(idx - 1)
        if idx < len(tss_arr):
            candidates.append(idx)
        best_gene = None
        best_dist = np.inf
        for c in candidates:
            dist = abs(int(tss_arr[c]) - mid)
            if dist < best_dist:
                best_dist = dist
                best_gene = gene_arr[c]
        gene_syms.append(best_gene)
        dists.append(best_dist)
    df = df.copy()
    df["gene_symbol"] = gene_syms
    df["distance_to_tss"] = dists
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 70)
    log.info("14_stage_stratified_da.py")
    log.info("=" * 70)

    # 1. Load donor metadata
    donor_meta = load_donor_metadata()

    # 2. Load h5ad
    log.info("Loading label-transferred h5ad (in-memory)...")
    t0 = time.time()
    adata = ad.read_h5ad(H5AD_PATH)
    log.info("  Loaded %s in %.1fs (X type=%s, nnz=%d)",
             adata.shape, time.time() - t0,
             type(adata.X).__name__,
             adata.X.nnz if sp.issparse(adata.X) else "?")

    # 3. Chrom range mapping (once)
    chrom_to_range = parse_chrom_sizes_from_tiles(list(adata.var_names))
    n_tiles = adata.n_vars

    # 4. Annotation TSS (once)
    chrom_tss = parse_tss_from_gtf(GTF_PATH)

    # 5. Per-CT loop
    # Buffer every (cell_type x test) result so a second, across-family global
    # BH (`padj_global`) can be computed by pooling p-values across the whole
    # run before anything is written to disk (review F058). Each entry is
    # (out_path, test_type, ct_name, result_df).
    buffered = []
    for ct_name, info in CT_GROUPS.items():
        log.info("=" * 70)
        log.info("Cell type: %s", ct_name)
        log.info("=" * 70)

        cell_mask = adata.obs["cell_type"].isin(info["h5ad_labels"]).values
        n_cells = int(cell_mask.sum())
        log.info("  %d cells matching labels %s", n_cells, info["h5ad_labels"])
        if n_cells < 100:
            log.warning("  Too few cells, skipping %s", ct_name)
            continue

        # Subset
        ct_adata = adata[cell_mask]
        donor_array = ct_adata.obs["donor_id"].values
        log.info("  Per-donor cell counts: %s",
                 pd.Series(donor_array).value_counts().sort_index().to_dict())

        # Peaks
        peaks = load_peaks(info["peak_bed"])
        if len(peaks) == 0:
            log.warning("  No peaks for %s, skipping", ct_name)
            continue
        S = build_peak_to_tile_mapping(peaks, chrom_to_range, n_tiles)

        # Pseudobulk
        log.info("  Materializing CT X submatrix...")
        # csr submatrix; use .X which is already in memory
        X_ct = ct_adata.X
        if not sp.isspmatrix_csr(X_ct):
            X_ct = X_ct.tocsr()
        counts_df, meta_df = pseudobulk_by_donor(X_ct, donor_array, S)

        # Run DA tests
        for test_type in ["ordinal", "F0_vs_F2F3", "F2F3_vs_F4", "F0_vs_F4",
                          "F0_vs_F3", "F3_vs_F4", "ordinal_F0_F3_F4"]:
            log.info("--- Test: %s ---", test_type)
            res = run_da(
                counts_df, meta_df, donor_meta, test_type, peaks
            )
            if res is None:
                continue
            # Annotate
            res = annotate_peaks(res, chrom_tss)
            res["cell_type"] = ct_name
            out_name = (
                f"da_stage_ordinal_{ct_name}.csv"
                if test_type == "ordinal"
                else f"da_{test_type}_{ct_name}.csv"
            )
            out_path = f"{OUT_DIR}/{out_name}"
            buffered.append((out_path, test_type, ct_name, res))

        # free memory before next CT
        del X_ct, S, counts_df, meta_df

    # 6. Across-family global BH-FDR. The per-family `padj` (per cell_type x
    # test) already lives on each result frame; `padj_global` pools p-values
    # across all buffered contrasts/cell-types of this run, since the ~35
    # contrasts are overlapping/nested and a single global correction guards
    # the run-wide FDR. Both columns are written side by side.
    if buffered:
        all_p = np.concatenate([r[3]["pvalue"].values for r in buffered])
        valid = np.isfinite(all_p)
        padj_global = np.full(all_p.shape, np.nan)
        if valid.sum() > 0:
            _, padj_g_valid, _, _ = multipletests(
                all_p[valid], method="fdr_bh"
            )
            padj_global[valid] = padj_g_valid
        log.info(
            "Global BH across %d buffered contrasts (%d pooled peaks, "
            "%d valid p-values): %d significant at padj_global<0.05",
            len(buffered), all_p.size, int(valid.sum()),
            int(np.nansum(padj_global < 0.05)),
        )
        # Scatter padj_global back to each result frame in pooled order
        offset = 0
        for out_path, test_type, ct_name, res in buffered:
            n = len(res)
            res = res.copy()
            res["padj_global"] = padj_global[offset:offset + n]
            offset += n
            res.to_csv(out_path, index=False)
            log.info(
                "  Saved %d rows -> %s (sig padj<0.05: %d, "
                "padj_global<0.05: %d)",
                n, out_path,
                int((res["padj"] < 0.05).sum()),
                int((res["padj_global"] < 0.05).sum()),
            )

    log.info("=" * 70)
    log.info("14_stage_stratified_da.py DONE")
    log.info("=" * 70)


if __name__ == "__main__":
    main()
