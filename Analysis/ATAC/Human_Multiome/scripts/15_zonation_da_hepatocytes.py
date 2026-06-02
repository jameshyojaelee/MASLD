#!/usr/bin/env python3
"""15_zonation_da_hepatocytes.py -- Zonation-resolved hepatocyte DA (ATAC B3b).

Chromatin-level confirmation of C1's snRNA periportal HNF4A finding.

Approach (mirrors Analysis/SingleCell/scripts/312e_hnf4a_zonation_resolved.py
"Approach A"):
  1. Subset the label-transferred h5ad to hepatocytes only.
  2. Compute per-cell PP and PC zonation scores using gene-activity matrix
     (one summed over promoter + body tiles). Markers are Halpern 2017 +
     Aizarani 2019 with HNF4A dropped from PP (it is one of the TFs we
     test for chromatin rewiring).
  3. zone_score = PP_score - PC_score; tercile split into Periportal / Mid
     / Pericentral.
  4. For each zone, pseudobulk per (donor x zone) and run the same
     OLS / limma-voom-style log2(CPM+1) linear-model stage tests
     (F0 / F2+F3 / F4) as Script 14 (per-donor; NOT a negative binomial GLM).
  5. HNF4A-specific test: which peaks change in Periportal-only across the
     F0_vs_F4 contrast? Compare to RORA/THRB/CEBPB.

Gene activity for zonation scoring:
  We approximate scanpy's score_genes by summing tile counts over a 2kb
  TSS window for each marker gene, then z-scoring per cell.

Outputs (results/stage_da/):
  - da_hep_zonation_stage.csv         (long: peak x zone x test)
  - hep_zone_summary.csv              (per-zone donor counts)
  - hep_zone_periportal_change_summary.csv (peaks changing in PP only,
    annotated to TFs of interest)

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
# Paths
# ---------------------------------------------------------------------------
PROJECT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATAC_DIR = f"{PROJECT}/Analysis/ATAC/Human_Multiome"
H5AD_PATH = f"{ATAC_DIR}/results/label_transfer/snapatac2_label_transferred.h5ad"
HEP_PEAK_BED = f"{ATAC_DIR}/results/label_transfer/cell_type_peak_sets_v2/Hepatocytes_peaks.bed"
DONOR_META_PATH = f"{ATAC_DIR}/metadata/donor_metadata_curated.tsv"
GTF_PATH = (
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/"
    "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
)
OUT_DIR = f"{ATAC_DIR}/results/stage_da"
os.makedirs(OUT_DIR, exist_ok=True)

TILE_SIZE = 500
PROMOTER_WINDOW = 2000  # +/- bp around TSS for gene-activity score

# ---------------------------------------------------------------------------
# Markers (Halpern 2017 + Aizarani 2019; HNF4A excluded from PP to avoid
# circularity)
# ---------------------------------------------------------------------------
PP_MARKERS = [
    "CPS1", "OTC", "ASS1", "ASL", "ARG1",
    "PCK1", "G6PC", "FBP1",
    "ALB", "TTR", "C3", "HAL",
    "SULT1A1", "IGFBP1",
    "SLCO1B1", "SLC10A1",
]

PC_MARKERS = [
    "CYP2E1", "CYP1A2", "CYP3A4",
    "GLUL", "OAT", "AXIN2", "LGR5",
    "CYP7A1", "CYP8B1",
    "ACLY", "SCD",
    "AKR1D1", "ALDH1A1",
]

# TFs of interest (from C1 + GWAS-ATAC plan)
TF_PANEL = ["HNF4A", "RORA", "THRB", "CEBPB", "NR1H4", "FOXA1",
            "FOXA2", "RXRA", "AR", "ESR1", "PPARA", "KLF15"]

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
# Helpers (shared with Script 14)
# ---------------------------------------------------------------------------
def load_donor_metadata():
    meta = pd.read_csv(DONOR_META_PATH, sep="\t")
    meta["F_stage_augmented"] = pd.to_numeric(
        meta["F_stage_augmented"], errors="coerce"
    )

    def bucket(f):
        if pd.isna(f):
            return np.nan
        if f == 0:
            return 0
        if f in (2, 3):
            return 1
        if f == 4:
            return 2
        return np.nan

    meta["F_bin"] = meta["F_stage_augmented"].apply(bucket)
    meta = meta.dropna(subset=["F_bin"]).copy()
    meta["F_bin"] = meta["F_bin"].astype(int)
    log.info("Donor F_bin distribution: %s",
             meta["F_bin"].value_counts().sort_index().to_dict())
    return meta[["donor_id", "F_stage_augmented", "F_bin"]].reset_index(drop=True)


def parse_chrom_sizes_from_tiles(var_names):
    chrom_to_range = {}
    cur_chrom = None
    cur_start = 0
    for i, nm in enumerate(var_names):
        chrom, _ = nm.split(":", 1)
        if chrom != cur_chrom:
            if cur_chrom is not None:
                chrom_to_range[cur_chrom] = (cur_start, i - 1)
            cur_chrom = chrom
            cur_start = i
    chrom_to_range[cur_chrom] = (cur_start, len(var_names) - 1)
    return chrom_to_range


def load_peaks(peak_bed_path):
    peaks = []
    with open(peak_bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            peaks.append((parts[0], int(parts[1]), int(parts[2])))
    log.info("Loaded %d peaks from %s", len(peaks), peak_bed_path)
    return peaks


def build_peak_to_tile_mapping(peaks, chrom_to_range, n_tiles):
    rows, cols = [], []
    n_peak = len(peaks)
    n_skipped = 0
    for p_idx, (chrom, start, end) in enumerate(peaks):
        if chrom not in chrom_to_range:
            n_skipped += 1
            continue
        chrom_start_idx, chrom_end_idx = chrom_to_range[chrom]
        tile_start = start // TILE_SIZE
        tile_end = (end - 1) // TILE_SIZE
        gs = chrom_start_idx + tile_start
        ge = chrom_start_idx + tile_end
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
    return sp.csr_matrix(
        (data, (rows, cols)), shape=(n_peak, n_tiles), dtype=np.float32
    )


def fit_lm_vectorized(Y_log_cpm, x):
    """Vectorized OLS regression of log2(CPM+1) on x.

    Returns (beta1, se, pval) arrays per peak.
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

    Y_mean = Y_log_cpm.mean(axis=0)
    Y_dev = Y_log_cpm - Y_mean
    beta1 = (x_dev[:, None] * Y_dev).sum(axis=0) / SSx
    beta0 = Y_mean - beta1 * x_mean
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
    from scipy import stats as scipy_stats
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, beta1 / se, np.nan)
    pval = np.full(n_peak, np.nan)
    valid = np.isfinite(tstat)
    pval[valid] = 2.0 * scipy_stats.t.sf(np.abs(tstat[valid]), df=df_resid)
    return beta1, se, pval


def parse_tss_from_gtf(gtf_path):
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
            records.append((chrom, tss, strand, gene_name))
    df = pd.DataFrame(records, columns=["chrom", "tss", "strand", "gene_name"])
    df = df.drop_duplicates(subset="gene_name", keep="first")
    log.info("  Parsed %d TSS", len(df))
    return df


def build_chrom_tss_lookup(tss_df):
    chrom_tss = {}
    for chrom, sub in tss_df.groupby("chrom"):
        sub_sorted = sub.sort_values("tss")
        chrom_tss[chrom] = (
            sub_sorted["tss"].values.astype(np.int64),
            sub_sorted["gene_name"].values,
        )
    return chrom_tss


def annotate_peaks(df, chrom_tss):
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
# Zonation scoring via tile-based gene activity
# ---------------------------------------------------------------------------
def gene_to_tile_indices(gene_name, tss_df, chrom_to_range,
                        window=PROMOTER_WINDOW, n_tiles=None):
    """Return the list of tile indices covering [tss-window, tss+window]."""
    row = tss_df[tss_df["gene_name"] == gene_name]
    if len(row) == 0:
        return None
    chrom = row.iloc[0]["chrom"]
    tss = int(row.iloc[0]["tss"])
    if chrom not in chrom_to_range:
        return None
    cs, ce = chrom_to_range[chrom]
    start = max(0, tss - window)
    end = tss + window
    t0 = start // TILE_SIZE
    t1 = (end - 1) // TILE_SIZE
    gs = cs + t0
    ge = cs + t1
    gs = max(gs, cs)
    ge = min(ge, ce)
    if ge < gs:
        return None
    return list(range(gs, ge + 1))


def score_gene_set(X_cells, marker_genes, tss_df, chrom_to_range, label):
    """Sum tile counts over marker promoter windows; z-score; mean across
    valid markers."""
    found, missing = [], []
    cell_n = X_cells.shape[0]
    n_tiles = X_cells.shape[1]
    score = np.zeros(cell_n, dtype=np.float32)
    raw_total = np.zeros(cell_n, dtype=np.float32)  # per-cell raw marker coverage
    valid = 0
    for g in marker_genes:
        idxs = gene_to_tile_indices(g, tss_df, chrom_to_range, n_tiles=n_tiles)
        if idxs is None or len(idxs) == 0:
            missing.append(g)
            continue
        # Sum tile counts within window per cell -> (n_cell,)
        gene_score = np.asarray(
            X_cells[:, idxs].sum(axis=1)
        ).ravel().astype(np.float32)
        # z-score across cells
        mu = gene_score.mean()
        sd = gene_score.std()
        if sd <= 0:
            missing.append(g)
            continue
        score += (gene_score - mu) / sd
        raw_total += gene_score          # accumulate raw coverage for QC mask
        valid += 1
        found.append(g)
    if valid == 0:
        log.warning("  %s: 0 valid markers!", label)
        return None, None
    score /= valid
    log.info("  %s score: %d/%d markers (%s used; %d missing)",
             label, valid, len(marker_genes), ",".join(found[:6]) + ("..." if len(found)>6 else ""), len(missing))
    # Return both the z-scored signature and the per-cell raw coverage so the
    # caller can drop zero-coverage cells (review 2026-05-30, F064): cells with
    # zero counts at every marker tile z-score to a single constant and would
    # otherwise flood one tercile bin.
    return score, raw_total


# ---------------------------------------------------------------------------
# Pseudobulk and DA test
# ---------------------------------------------------------------------------
def pseudobulk_by_donor(X_cells, donor_array, peak_to_tile_S):
    donors = sorted(np.unique(donor_array))
    rows = []
    libsizes = []
    for d in donors:
        mask = donor_array == d
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        tile_sum = np.asarray(
            X_cells[mask].sum(axis=0)
        ).ravel().astype(np.float32)
        libsize = float(tile_sum.sum())
        peak_counts = peak_to_tile_S @ tile_sum
        rows.append((d, n_cells, libsize, peak_counts))
        libsizes.append(libsize)
    donors_used = [r[0] for r in rows]
    n_cells_per = [r[1] for r in rows]
    libs = [r[2] for r in rows]
    counts = np.vstack([r[3] for r in rows])
    df_counts = pd.DataFrame(counts, index=donors_used)
    df_counts.index.name = "donor_id"
    df_meta = pd.DataFrame({"n_cells": n_cells_per, "libsize": libs},
                           index=donors_used)
    return df_counts, df_meta


def run_da(counts_df, meta_df, donor_meta_df, test_type, peaks_coords, zone_tag):
    donor_meta_idx = donor_meta_df.set_index("donor_id")
    common = [d for d in counts_df.index if d in donor_meta_idx.index]
    if not common:
        return None
    counts_df = counts_df.loc[common]
    meta_df = meta_df.loc[common]
    F_bin = donor_meta_idx.loc[common, "F_bin"].values.astype(int)

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
    else:
        raise ValueError(test_type)

    counts_df = counts_df.iloc[keep]
    libsize = meta_df.iloc[keep]["libsize"].values
    n_donors = len(counts_df)
    x = x[keep]
    log.info("  [%s/%s] %d donors, distribution %s",
             zone_tag, test_type, n_donors,
             np.unique(x, return_counts=True))
    if n_donors < 4:
        return None

    counts_arr = counts_df.values.astype(np.float32)
    peak_total = counts_arr.sum(axis=0)
    n_donor_detect = (counts_arr > 0).sum(axis=0)
    keep_peak = (peak_total >= 10) & (n_donor_detect >= 3)
    n_keep = int(keep_peak.sum())
    log.info("  [%s/%s] peaks passing filter: %d / %d",
             zone_tag, test_type, n_keep, len(keep_peak))
    if n_keep == 0:
        return None

    peak_indices = np.where(keep_peak)[0]
    counts_kept = counts_arr[:, peak_indices]
    cpm = counts_kept / libsize[:, None] * 1e6
    Y = np.log2(cpm + 1.0)
    log.info("  [%s/%s] Fitting vectorized OLS log2(CPM+1) on %d peaks...",
             zone_tag, test_type, n_keep)
    t0 = time.time()
    coefs, ses, pvals = fit_lm_vectorized(Y, x)
    log.info("  [%s/%s] OLS done in %.1f sec",
             zone_tag, test_type, time.time() - t0)
    # BH FDR within this (zone x test) family. A second, across-family global
    # BH (`padj_global`) is added in main() by pooling p-values across all
    # zone x test contrasts of this run (3 zones x 4 overlapping/nested tests).
    # Both are reported: `padj` = per-family (per zone x test),
    # `padj_global` = pooled across the whole run.
    valid = np.isfinite(pvals)
    padj = np.full(n_keep, np.nan)
    if valid.sum() > 0:
        _, padj_v, _, _ = multipletests(pvals[valid], method="fdr_bh")
        padj[valid] = padj_v

    rec = []
    for j, p_idx in enumerate(peak_indices):
        chrom, start, end = peaks_coords[p_idx]
        rec.append({
            "peak_id": f"{chrom}:{start}-{end}",
            "chrom": chrom, "start": start, "end": end,
            "log2FC": coefs[j],
            "se": ses[j],
            "pvalue": pvals[j],
            "padj": padj[j],
            "zone": zone_tag,
            "test": test_type,
        })
    out_df = pd.DataFrame(rec).sort_values("pvalue", na_position="last")
    n_sig = int((out_df["padj"] < 0.05).sum())
    log.info("  [%s/%s] sig padj<0.05: %d", zone_tag, test_type, n_sig)
    return out_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 70)
    log.info("15_zonation_da_hepatocytes.py")
    log.info("=" * 70)

    donor_meta = load_donor_metadata()

    log.info("Loading h5ad (in-memory)...")
    t0 = time.time()
    adata = ad.read_h5ad(H5AD_PATH)
    log.info("  Loaded %s in %.1fs", adata.shape, time.time() - t0)

    hep_mask = adata.obs["cell_type"].values == "Hepatocyte"
    n_hep = int(hep_mask.sum())
    log.info("Hepatocytes: %d", n_hep)
    hep = adata[hep_mask]

    # ensure CSR
    X_hep = hep.X if sp.isspmatrix_csr(hep.X) else hep.X.tocsr()

    chrom_to_range = parse_chrom_sizes_from_tiles(list(adata.var_names))
    tss_df = parse_tss_from_gtf(GTF_PATH)
    chrom_tss = build_chrom_tss_lookup(tss_df)

    # -----------------------------------------------------------------------
    # Step 1: per-cell PP / PC score
    # -----------------------------------------------------------------------
    log.info("Scoring PP / PC zonation...")
    pp_score, pp_reads = score_gene_set(X_hep, PP_MARKERS, tss_df, chrom_to_range,
                                        "PP_score")
    pc_score, pc_reads = score_gene_set(X_hep, PC_MARKERS, tss_df, chrom_to_range,
                                        "PC_score")
    if pp_score is None or pc_score is None:
        log.error("Zonation scoring failed; aborting.")
        sys.exit(1)
    zone_score = pp_score - pc_score

    # -----------------------------------------------------------------------
    # Step 2: tercile -> Periportal / Mid / Pericentral
    # Exclude zero-coverage cells first (review 2026-05-30, F064). A cell with
    # no counts at any marker promoter tile z-scores to a single constant
    # zone_score; left in, ~40% of hepatocytes collapsed onto that constant and
    # were all dumped into Pericentral, contaminating every PP-vs-PC contrast.
    # Terciles are computed on covered cells only; uncovered cells are labelled
    # "Unscored" and are excluded from all downstream per-zone DA.
    # -----------------------------------------------------------------------
    covered = (pp_reads + pc_reads) > 0
    n_cov = int(covered.sum())
    log.info("  Zonation coverage: %d/%d hepatocytes have marker-tile signal (%.1f%%); "
             "%d zero-coverage cells excluded from terciles",
             n_cov, covered.size, 100.0 * n_cov / covered.size, covered.size - n_cov)
    if n_cov < 30:
        log.error("Too few covered cells (%d) for zonation terciles; aborting.", n_cov)
        sys.exit(1)
    log.info("  zone_score (covered): mean=%.3f, std=%.3f, q33=%.3f, q67=%.3f",
             zone_score[covered].mean(), zone_score[covered].std(),
             np.quantile(zone_score[covered], 1/3), np.quantile(zone_score[covered], 2/3))

    q_lo, q_hi = np.quantile(zone_score[covered], [1/3, 2/3])
    zone_label = np.full(zone_score.shape, "Unscored", dtype=object)
    zone_label[covered & (zone_score <= q_lo)] = "Pericentral"
    zone_label[covered & (zone_score >= q_hi)] = "Periportal"
    zone_label[covered & (zone_score > q_lo) & (zone_score < q_hi)] = "Mid"
    zone_counts = dict(zip(*np.unique(zone_label, return_counts=True)))
    log.info("  Zone counts: %s", zone_counts)
    # Sanity guard: the three real bins should be roughly balanced terciles.
    real_bins = np.array([zone_counts.get("Periportal", 0),
                          zone_counts.get("Mid", 0),
                          zone_counts.get("Pericentral", 0)], dtype=float)
    if real_bins.min() < 0.5 * (n_cov / 3.0):
        log.warning("  Tercile bins are unbalanced (%s) despite excluding "
                    "zero-coverage cells; check zonation score distribution.",
                    real_bins.astype(int).tolist())

    donor_array = hep.obs["donor_id"].values
    # Save per-cell zone CSV for diagnostics
    pd.DataFrame({
        "barcode": hep.obs_names,
        "donor_id": donor_array,
        "PP_score": pp_score, "PC_score": pc_score,
        "zone_score": zone_score, "zone": zone_label,
    }).to_csv(f"{OUT_DIR}/hep_zone_assignments.csv", index=False)

    # Per-zone donor summary
    summary_rows = []
    for z in ["Periportal", "Mid", "Pericentral"]:
        z_mask = zone_label == z
        n_cells = int(z_mask.sum())
        per_donor = pd.Series(donor_array[z_mask]).value_counts().to_dict()
        summary_rows.append({
            "zone": z, "n_cells": n_cells,
            "n_donors_with_cells": len(per_donor),
            "per_donor_counts": str(sorted(per_donor.items())),
        })
    pd.DataFrame(summary_rows).to_csv(
        f"{OUT_DIR}/hep_zone_summary.csv", index=False
    )

    # -----------------------------------------------------------------------
    # Step 3: peaks + tile->peak mapping (same for all zones)
    # -----------------------------------------------------------------------
    peaks = load_peaks(HEP_PEAK_BED)
    n_tiles = adata.n_vars
    S = build_peak_to_tile_mapping(peaks, chrom_to_range, n_tiles)

    # -----------------------------------------------------------------------
    # Step 4: per-zone DA test
    # -----------------------------------------------------------------------
    all_results = []
    for zone in ["Periportal", "Mid", "Pericentral"]:
        log.info("=" * 60)
        log.info("Zone: %s", zone)
        log.info("=" * 60)
        z_mask = zone_label == zone
        n_cells = int(z_mask.sum())
        if n_cells < 500:
            log.warning("  Zone %s has only %d cells; skipping", zone, n_cells)
            continue
        donor_z = donor_array[z_mask]
        X_z = X_hep[z_mask]
        log.info("  Pseudobulking %d cells across %d donors",
                 n_cells, len(np.unique(donor_z)))
        counts_df, meta_df = pseudobulk_by_donor(X_z, donor_z, S)

        for test_type in ["ordinal", "F0_vs_F2F3", "F2F3_vs_F4", "F0_vs_F4"]:
            res = run_da(counts_df, meta_df, donor_meta, test_type, peaks, zone)
            if res is None:
                continue
            all_results.append(res)

    if not all_results:
        log.error("No DA results produced; aborting.")
        sys.exit(1)

    full = pd.concat(all_results, ignore_index=True)
    full = annotate_peaks(full, chrom_tss)

    # Across-family global BH-FDR. `padj` is per-family (per zone x test);
    # `padj_global` pools p-values across all zone x test contrasts of this run
    # (the contrasts are overlapping/nested), guarding the run-wide FDR. Both
    # columns are reported side by side. The Periportal-only summary below
    # intentionally still keys off the per-family `padj` (zone-specific call).
    glob_valid = np.isfinite(full["pvalue"].values)
    full["padj_global"] = np.nan
    if glob_valid.sum() > 0:
        _, padj_g, _, _ = multipletests(
            full["pvalue"].values[glob_valid], method="fdr_bh"
        )
        full.loc[glob_valid, "padj_global"] = padj_g
    log.info(
        "Global BH across all zone x test contrasts (%d pooled peaks, "
        "%d valid p-values): %d significant at padj_global<0.05",
        len(full), int(glob_valid.sum()),
        int((full["padj_global"] < 0.05).sum()),
    )

    full.to_csv(f"{OUT_DIR}/da_hep_zonation_stage.csv", index=False)
    log.info("Wrote %d rows -> %s",
             len(full), f"{OUT_DIR}/da_hep_zonation_stage.csv")

    # -----------------------------------------------------------------------
    # Step 5: Periportal-only change summary
    # -----------------------------------------------------------------------
    log.info("Computing periportal-only change summary (F0 vs F4)...")
    f0f4 = full[full["test"] == "F0_vs_F4"].copy()
    if len(f0f4) == 0:
        log.warning("No F0_vs_F4 rows; cannot build PP-only summary.")
        return

    # Pivot: rows = peak_id, cols = zone, values = padj/log2FC
    pivot_p = f0f4.pivot_table(index="peak_id", columns="zone",
                                values="padj", aggfunc="first")
    pivot_lfc = f0f4.pivot_table(index="peak_id", columns="zone",
                                 values="log2FC", aggfunc="first")
    pivot_gene = f0f4.groupby("peak_id")["gene_symbol"].first()
    pivot_dist = f0f4.groupby("peak_id")["distance_to_tss"].first()

    # PP-only: PP sig (padj<0.1), Mid and PC not sig (padj>=0.1 or NA)
    pp_sig = pivot_p.get("Periportal").fillna(1.0) < 0.1
    mid_ns = pivot_p.get("Mid", pd.Series(1, index=pivot_p.index)).fillna(1.0) >= 0.1
    pc_ns = pivot_p.get("Pericentral", pd.Series(1, index=pivot_p.index)).fillna(1.0) >= 0.1
    pp_only = pivot_p.index[pp_sig & mid_ns & pc_ns]
    log.info("  PP-only significant peaks (padj<0.1): %d", len(pp_only))

    summary = pd.DataFrame({
        "peak_id": pp_only,
        "gene_symbol": [pivot_gene.get(p, None) for p in pp_only],
        "distance_to_tss": [pivot_dist.get(p, None) for p in pp_only],
        "padj_PP": [pivot_p.loc[p].get("Periportal", np.nan) for p in pp_only],
        "padj_Mid": [pivot_p.loc[p].get("Mid", np.nan) for p in pp_only],
        "padj_PC": [pivot_p.loc[p].get("Pericentral", np.nan) for p in pp_only],
        "log2FC_PP": [pivot_lfc.loc[p].get("Periportal", np.nan) for p in pp_only],
        "log2FC_Mid": [pivot_lfc.loc[p].get("Mid", np.nan) for p in pp_only],
        "log2FC_PC": [pivot_lfc.loc[p].get("Pericentral", np.nan) for p in pp_only],
    })
    summary = summary.sort_values("padj_PP")
    summary.to_csv(f"{OUT_DIR}/hep_zone_periportal_change_summary.csv",
                   index=False)
    log.info("  Wrote PP-only summary -> hep_zone_periportal_change_summary.csv")

    # TF-of-interest spotlight
    log.info("=" * 60)
    log.info("TF spotlight (PP-only peaks near TF genes):")
    for tf in TF_PANEL:
        hits = summary[summary["gene_symbol"] == tf]
        log.info(
            "  %-7s: %d PP-only peaks (best padj=%.3g, max |log2FC|=%.2f)" if len(hits) else "  %-7s: %d PP-only peaks",
            tf, len(hits),
            *( [hits["padj_PP"].min(), hits["log2FC_PP"].abs().max()] if len(hits) else [] )
        )

    # Also dump TFs-of-interest annotation across all peaks/zones for ref
    tf_focus = full[full["gene_symbol"].isin(TF_PANEL)].copy()
    tf_focus.to_csv(f"{OUT_DIR}/hep_zone_tf_focus.csv", index=False)
    log.info("Wrote TF-focus annotation -> hep_zone_tf_focus.csv (%d rows)",
             len(tf_focus))

    log.info("DONE")


if __name__ == "__main__":
    main()
