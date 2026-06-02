#!/usr/bin/env python3
"""
141_atac_masl_vs_mash.py — ATAC-seq Differential Accessibility for Disease Progression
---------------------------------------------------------------------------------------
Runs differential accessibility (DA) for two new contrasts per cell type:
  1. MASL vs MASH  (steatosis -> steatohepatitis progression)
  2. MASL vs NORMAL (steatosis onset)

Method: Logistic regression with covariate adjustment (log10 n_fragment + tsse),
following the approach from 07b_corrected_da_hepatocytes.py but extended to all
cell types with sufficient sample size.

For cell types with too few cells for logistic regression (< 50 per group),
falls back to SnapATAC2 snap.tl.diff_test (Wilcoxon rank-sum), matching the
approach from 07_run_da_test.py.

After DA testing, peaks are mapped to nearest gene TSS (GENCODE v49) and
overlapped with C2 progression DEGs (NAFL-vs-NASH dream, Script 13) to quantify
epigenomic-transcriptomic concordance.

Input:
  - snapatac2_label_transferred.h5ad (label-transferred scATAC from ATAC pipeline)
  - Per-donor h5ad files (for QC covariates: n_fragment, tsse)
  - GENCODE v49 GTF (for peak-to-gene mapping)
  - nafl_vs_nash_dream.csv (C2 progression DEGs from Script 13)

Output (to results/progression/):
  - atac_masl_vs_mash_da.csv       DA peaks for MASL vs MASH per cell type
  - atac_masl_vs_normal_da.csv     DA peaks for MASL vs NORMAL per cell type
  - atac_progression_overlap.csv   Overlap of DA peaks with C2 progression DEGs

SLURM: cpu, 8 CPUs, 128G RAM, 48h
Env:   micromamba activate snapatac2
"""

import argparse
import gzip
import logging
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Defaults ──────────────────────────────────────────────────────────────
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATAC_DIR = os.path.join(BASE, "Analysis/ATAC/Human_Multiome")
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")

DEFAULT_INPUT = os.path.join(ATAC_DIR, "results/label_transfer/snapatac2_label_transferred.h5ad")
DEFAULT_PER_DONOR = os.path.join(ATAC_DIR, "results/snapatac2/per_donor")
DEFAULT_OUTPUT = os.path.join(INTEG, "results/progression")
DEFAULT_GTF = "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
DEFAULT_C2_DEGS = os.path.join(INTEG, "results/disease_signatures/nafl_vs_nash_dream.csv")

# Minimum cells per group for logistic regression; below this, use Wilcoxon
MIN_CELLS_LOGISTIC = 50
# Minimum cells per group for any DA test
MIN_CELLS_ANY = 5
# Minimum peak detection rate in either group
MIN_DETECTION = 0.05
# Promoter window for peak-to-gene mapping (bp)
PROMOTER_WINDOW = 5000


# ── QC covariate recovery ────────────────────────────────────────────────
def recover_qc_covariates(per_donor_dir, cell_names):
    """Recover n_fragment and tsse from per-donor h5ad files."""
    log.info("Recovering QC covariates from per-donor h5ads...")
    records = []
    cell_name_set = set(cell_names)

    donor_files = sorted(f for f in os.listdir(per_donor_dir) if f.endswith(".h5ad"))
    for fname in donor_files:
        path = os.path.join(per_donor_dir, fname)
        adata_d = ad.read_h5ad(path, backed="r")
        obs = adata_d.obs[["n_fragment", "tsse"]].copy()
        obs.index = adata_d.obs_names
        keep = obs.index.isin(cell_name_set)
        if keep.any():
            records.append(obs.loc[keep])
        adata_d.file.close()

    if not records:
        log.warning("  No QC covariates recovered from per-donor h5ads")
        return pd.DataFrame(columns=["n_fragment", "tsse"])

    covariates = pd.concat(records)
    covariates = covariates.loc[covariates.index.isin(cell_name_set)]
    log.info("  Recovered covariates for %d / %d cells", len(covariates), len(cell_names))
    return covariates


# ── Logistic regression DA ───────────────────────────────────────────────
def run_logistic_da(X_sub, design_df, peak_names, condition_col="condition_binary"):
    """Per-peak logistic regression: accessible ~ condition + log10_nfrag + tsse.

    Returns DataFrame with columns matching the 07/07b output format.
    """
    n_peaks = X_sub.shape[1]
    log.info("  Running logistic regression on %d peaks...", n_peaks)

    exog = design_df[["intercept", condition_col, "log10_nfrag", "tsse"]].values

    results = []
    t0 = time.time()
    n_failed = 0

    for i in range(n_peaks):
        if i > 0 and i % 10000 == 0:
            elapsed = time.time() - t0
            rate = i / elapsed
            eta = (n_peaks - i) / rate / 60
            log.info("    %d / %d peaks (%.0f/s, ETA %.1f min)", i, n_peaks, rate, eta)

        col = X_sub[:, i]
        if sparse.issparse(col):
            y = np.asarray(col.todense()).ravel()
        else:
            y = np.asarray(col).ravel()

        y = (y > 0).astype(np.float64)

        try:
            model = GLM(y, exog, family=Binomial())
            fit = model.fit(disp=False, maxiter=25)

            coef = fit.params[1]
            se = fit.bse[1]
            pval = fit.pvalues[1]
            log2fc = coef / np.log(2)

            results.append({
                "feature name": peak_names[i],
                "log2(fold_change)": log2fc,
                "p-value": pval,
                "coef_logodds": coef,
                "se": se,
            })
        except Exception:
            n_failed += 1
            results.append({
                "feature name": peak_names[i],
                "log2(fold_change)": np.nan,
                "p-value": np.nan,
                "coef_logodds": np.nan,
                "se": np.nan,
            })

    elapsed = time.time() - t0
    log.info("    Done: %d peaks in %.1f min (%.0f/s), %d failed",
             n_peaks, elapsed / 60, n_peaks / max(elapsed, 1), n_failed)

    df = pd.DataFrame(results)
    df = df.dropna(subset=["p-value"])

    if len(df) > 0:
        _, padj, _, _ = multipletests(df["p-value"].values, method="fdr_bh")
        df["adjusted p-value"] = padj
    else:
        df["adjusted p-value"] = []

    return df


# ── Wilcoxon fallback DA (via SnapATAC2) ─────────────────────────────────
def run_wilcoxon_da(adata, group1_cells, group2_cells):
    """Run SnapATAC2 diff_test (Wilcoxon rank-sum) for cell types with few cells.

    Follows the approach from 07_run_da_test.py.
    """
    try:
        import snapatac2 as snap
    except ImportError:
        log.error("snapatac2 not available for Wilcoxon fallback")
        return pd.DataFrame()

    try:
        result = snap.tl.diff_test(
            adata,
            cell_group1=group1_cells,
            cell_group2=group2_cells,
        )
        if hasattr(result, "to_pandas"):
            df = result.to_pandas()
        else:
            df = result
        return df
    except Exception as e:
        log.error("  Wilcoxon DA failed: %s", str(e))
        return pd.DataFrame()


# ── GTF parsing for peak-to-gene mapping ─────────────────────────────────
def parse_tss_from_gtf(gtf_path):
    """Parse gene TSS positions from GENCODE GTF.

    Returns DataFrame: gene_name, chrom, tss, strand, gene_type
    """
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
            gene_type = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]
                elif attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]

            if gene_name is None:
                continue

            tss = start if strand == "+" else end
            records.append({
                "gene_name": gene_name,
                "chrom": chrom,
                "tss": tss,
                "strand": strand,
                "gene_type": gene_type or "",
            })

    df = pd.DataFrame(records)
    df = df.drop_duplicates(subset="gene_name", keep="first")
    log.info("  Parsed %d gene TSS positions", len(df))
    return df


def annotate_peaks_with_genes(da_df, tss_df, promoter_window=PROMOTER_WINDOW):
    """Map DA peaks to nearest gene TSS.

    For each peak, annotate with nearest gene symbol and distance.
    """
    if da_df.empty:
        return da_df

    peak_col = "feature name"
    coords = da_df[peak_col].str.extract(r"(chr\w+):(\d+)-(\d+)")
    coords.columns = ["chrom", "start", "end"]
    coords["start"] = pd.to_numeric(coords["start"], errors="coerce")
    coords["end"] = pd.to_numeric(coords["end"], errors="coerce")
    coords["mid"] = ((coords["start"] + coords["end"]) // 2).astype("Int64")

    da_df = da_df.copy()
    da_df["chrom"] = coords["chrom"].values
    da_df["peak_mid"] = coords["mid"].values

    # Build TSS lookup by chromosome
    chrom_tss = {}
    for chrom in tss_df["chrom"].unique():
        ct = tss_df[tss_df["chrom"] == chrom].sort_values("tss")
        chrom_tss[chrom] = ct[["gene_name", "tss"]].values

    gene_symbols = []
    distances = []

    for _, row in da_df.iterrows():
        chrom = row["chrom"]
        mid = row["peak_mid"]

        if chrom not in chrom_tss or pd.isna(mid):
            gene_symbols.append(None)
            distances.append(np.nan)
            continue

        genes_on_chrom = chrom_tss[chrom]
        tss_positions = genes_on_chrom[:, 1].astype(int)

        idx = np.searchsorted(tss_positions, int(mid))
        candidates = []
        if idx > 0:
            candidates.append(idx - 1)
        if idx < len(tss_positions):
            candidates.append(idx)

        best_gene = None
        best_dist = np.inf
        for c in candidates:
            dist = abs(int(tss_positions[c]) - int(mid))
            if dist < best_dist:
                best_dist = dist
                best_gene = genes_on_chrom[c, 0]

        gene_symbols.append(best_gene)
        distances.append(best_dist)

    da_df["gene_symbol"] = gene_symbols
    da_df["distance_to_tss"] = distances
    da_df["within_promoter"] = da_df["distance_to_tss"] <= promoter_window

    n_mapped = da_df["gene_symbol"].notna().sum()
    n_promoter = da_df["within_promoter"].sum()
    log.info("  Peak annotation: %d / %d mapped, %d within %dkb promoter",
             n_mapped, len(da_df), n_promoter, promoter_window // 1000)

    return da_df


# ── Per-contrast DA for all cell types ────────────────────────────────────
def run_contrast_da(adata, obs_df, covariates, all_names, peak_names,
                    group1_labels, group2_labels,
                    contrast_name, cell_type_col="cell_type"):
    """Run DA for a single contrast across all cell types.

    group1 is the numerator (positive logFC = more accessible in group1).
    """
    log.info("=" * 70)
    log.info("Contrast: %s", contrast_name)
    log.info("  Group1 (numerator): %s", group1_labels)
    log.info("  Group2 (denominator): %s", group2_labels)

    cell_types = obs_df[cell_type_col].unique().tolist()
    all_results = []

    for ct in sorted(cell_types):
        ct_mask = obs_df[cell_type_col] == ct
        g1_mask = ct_mask & obs_df["condition"].isin(group1_labels)
        g2_mask = ct_mask & obs_df["condition"].isin(group2_labels)

        n_g1 = int(g1_mask.sum())
        n_g2 = int(g2_mask.sum())

        if n_g1 < MIN_CELLS_ANY or n_g2 < MIN_CELLS_ANY:
            log.info("  %s: SKIP (group1=%d, group2=%d < %d)",
                     ct, n_g1, n_g2, MIN_CELLS_ANY)
            continue

        log.info("  %s: group1=%d, group2=%d", ct, n_g1, n_g2)

        use_logistic = (n_g1 >= MIN_CELLS_LOGISTIC and n_g2 >= MIN_CELLS_LOGISTIC)

        if use_logistic:
            log.info("    Method: logistic regression (covariate-adjusted)")
            # Subset cells for this cell type + contrast
            contrast_mask = g1_mask | g2_mask
            ct_obs = obs_df[contrast_mask].copy()
            ct_names = list(ct_obs.index)

            # Merge QC covariates
            ct_obs = ct_obs.join(covariates, how="left")
            missing = ct_obs["n_fragment"].isna().sum()
            if missing > 0:
                log.info("    %d cells missing covariates, dropping", missing)
                ct_obs = ct_obs.dropna(subset=["n_fragment", "tsse"])
                ct_names = list(ct_obs.index)

            if len(ct_names) < MIN_CELLS_ANY * 2:
                log.info("    Too few cells after covariate join, skipping")
                continue

            ct_obs["condition_binary"] = ct_obs["condition"].isin(group1_labels).astype(float)
            ct_obs["log10_nfrag"] = np.log10(ct_obs["n_fragment"])
            ct_obs["intercept"] = 1.0

            # Sync ct_obs with peak matrix: only keep cells in both obs and matrix
            name_to_idx = {name: i for i, name in enumerate(all_names)}
            # Filter to cells present in the peak matrix
            valid_cell_mask = ct_obs.index.isin(name_to_idx)
            ct_obs = ct_obs[valid_cell_mask].copy()
            # Deduplicate index in case of any duplicate cell names
            if ct_obs.index.duplicated().any():
                log.warning("    Removing %d duplicate cell names", ct_obs.index.duplicated().sum())
                ct_obs = ct_obs[~ct_obs.index.duplicated(keep='first')]
            n_in_matrix = len(ct_obs)
            log.info("    Cells in both obs and peak matrix: %d (filtered %d)",
                     n_in_matrix, int((~valid_cell_mask).sum()))

            if n_in_matrix < 20:
                log.warning("    Too few cells after sync (%d) — skipping", n_in_matrix)
                continue

            # Build sorted index mapping (sorted by global position for consistent row order)
            cell_to_global = {n: name_to_idx[n] for n in ct_obs.index}
            sorted_cells = sorted(cell_to_global.keys(), key=lambda x: cell_to_global[x])
            ct_indices = [cell_to_global[n] for n in sorted_cells]

            chunk_size = 5000
            X_chunks = []
            for start in range(0, len(ct_indices), chunk_size):
                chunk_idx = ct_indices[start:start + chunk_size]
                X_chunk = adata.X[chunk_idx, :]
                if not sparse.issparse(X_chunk):
                    X_chunk = sparse.csr_matrix(X_chunk)
                X_chunks.append(X_chunk)

            X_ct = sparse.vstack(X_chunks, format="csr")

            # Reorder ct_obs to match X_ct row order using the same sorted order
            ct_obs = ct_obs.loc[sorted_cells]

            if X_ct.shape[0] != len(ct_obs):
                log.error("Shape mismatch after sync: X_ct=%d, ct_obs=%d. Skipping cell type.",
                          X_ct.shape[0], len(ct_obs))
                continue

            # Filter peaks by detection rate
            X_bin = (X_ct > 0).astype(np.float32)
            g1_arr = ct_obs["condition_binary"].values == 1.0
            g2_arr = ~g1_arr

            det_g1 = np.asarray(X_bin[g1_arr].mean(axis=0)).ravel()
            det_g2 = np.asarray(X_bin[g2_arr].mean(axis=0)).ravel()
            det_max = np.maximum(det_g1, det_g2)
            keep_peaks = det_max >= MIN_DETECTION

            n_keep = int(keep_peaks.sum())
            log.info("    Peaks passing detection filter: %d / %d", n_keep, len(keep_peaks))

            X_filtered = X_ct[:, keep_peaks]
            peaks_filtered = [peak_names[i] for i in range(len(peak_names)) if keep_peaks[i]]

            del X_ct, X_bin, X_chunks

            # Design matrix
            design_df = ct_obs[["intercept", "condition_binary", "log10_nfrag", "tsse"]].copy()
            design_df = design_df.reset_index(drop=True)

            da_result = run_logistic_da(X_filtered, design_df, peaks_filtered,
                                        condition_col="condition_binary")
            da_result["cell_type"] = ct
            da_result["method"] = "logistic_regression"
            all_results.append(da_result)

            del X_filtered

        else:
            log.info("    Method: Wilcoxon rank-sum (SnapATAC2 diff_test)")
            g1_cells = obs_df.index[g1_mask].tolist()
            g2_cells = obs_df.index[g2_mask].tolist()

            da_result = run_wilcoxon_da(adata, g1_cells, g2_cells)
            if not da_result.empty:
                da_result["cell_type"] = ct
                da_result["method"] = "wilcoxon"
                all_results.append(da_result)

                # Count significant
                for col in ["adjusted p-value", "adjusted_p_value", "padj", "fdr"]:
                    if col in da_result.columns:
                        n_sig = int((da_result[col] < 0.05).sum())
                        log.info("    %d DA peaks (padj<0.05)", n_sig)
                        break
            else:
                log.warning("    Wilcoxon DA produced no results")

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        log.info("  Total DA results for %s: %d rows across %d cell types",
                 contrast_name, len(combined), combined["cell_type"].nunique())
        return combined
    else:
        log.warning("  No DA results for contrast %s", contrast_name)
        return pd.DataFrame()


# ── Overlap with C2 progression DEGs ──────────────────────────────────────
def compute_progression_overlap(da_annotated, c2_degs_path, padj_thresh=0.05,
                                lfc_thresh=0.25):
    """Overlap DA peaks with C2 (NAFL-vs-NASH) progression DEGs.

    For each cell type, counts:
      - DA peaks mapped to progression DEGs
      - Concordant direction (DA up + DEG up, or DA down + DEG down)
    """
    log.info("=" * 70)
    log.info("Computing overlap with C2 progression DEGs")

    if not os.path.exists(c2_degs_path):
        log.warning("  C2 DEG file not found: %s", c2_degs_path)
        return pd.DataFrame()

    c2 = pd.read_csv(c2_degs_path)
    log.info("  C2 DEGs: %d genes total", len(c2))

    # Significant C2 DEGs
    c2_sig = c2[c2["adj.P.Val"] < padj_thresh].copy()
    log.info("  C2 significant (padj<%g): %d genes", padj_thresh, len(c2_sig))

    if "symbol" in c2_sig.columns:
        c2_genes = set(c2_sig["symbol"].dropna())
    elif "gene_symbol" in c2_sig.columns:
        c2_genes = set(c2_sig["gene_symbol"].dropna())
    else:
        log.warning("  Cannot find gene symbol column in C2 DEGs")
        return pd.DataFrame()

    log.info("  C2 unique gene symbols: %d", len(c2_genes))

    # Build C2 direction map
    sym_col = "symbol" if "symbol" in c2_sig.columns else "gene_symbol"
    c2_dir = dict(zip(c2_sig[sym_col], np.sign(c2_sig["logFC"])))

    # Filter DA to significant peaks with gene annotation
    padj_col = None
    for col in ["adjusted p-value", "adjusted_p_value", "padj", "fdr"]:
        if col in da_annotated.columns:
            padj_col = col
            break

    if padj_col is None:
        log.warning("  No p-value column in DA results")
        return pd.DataFrame()

    lfc_col = None
    for col in ["log2(fold_change)", "logFC", "log2FC"]:
        if col in da_annotated.columns:
            lfc_col = col
            break

    if lfc_col is None:
        log.warning("  No logFC column in DA results")
        return pd.DataFrame()

    da_sig = da_annotated[
        (da_annotated[padj_col] < padj_thresh) &
        (da_annotated[lfc_col].abs() > lfc_thresh) &
        da_annotated["gene_symbol"].notna()
    ].copy()

    log.info("  Significant DA peaks with gene annotation: %d", len(da_sig))

    overlap_records = []

    cell_types = da_sig["cell_type"].unique() if "cell_type" in da_sig.columns else ["all"]

    for ct in sorted(cell_types):
        if ct == "all":
            ct_da = da_sig
        else:
            ct_da = da_sig[da_sig["cell_type"] == ct]

        da_genes = set(ct_da["gene_symbol"].dropna())
        overlap_genes = da_genes & c2_genes

        # Direction concordance
        n_concordant = 0
        n_discordant = 0
        concordant_genes = []
        discordant_genes = []

        for gene in overlap_genes:
            da_rows = ct_da[ct_da["gene_symbol"] == gene]
            da_direction = np.sign(da_rows[lfc_col].mean())
            c2_direction = c2_dir.get(gene, 0)

            if da_direction != 0 and c2_direction != 0:
                if da_direction == c2_direction:
                    n_concordant += 1
                    concordant_genes.append(gene)
                else:
                    n_discordant += 1
                    discordant_genes.append(gene)

        overlap_records.append({
            "cell_type": ct,
            "n_da_sig_peaks": len(ct_da),
            "n_da_unique_genes": len(da_genes),
            "n_c2_sig_genes": len(c2_genes),
            "n_overlap_genes": len(overlap_genes),
            "n_concordant": n_concordant,
            "n_discordant": n_discordant,
            "concordance_rate": n_concordant / max(n_concordant + n_discordant, 1),
            "overlap_pct_of_da": len(overlap_genes) / max(len(da_genes), 1),
            "concordant_genes": ";".join(sorted(concordant_genes)[:50]),
            "discordant_genes": ";".join(sorted(discordant_genes)[:50]),
        })

        log.info("  %s: %d DA genes, %d overlap with C2, %d concordant / %d discordant (%.1f%%)",
                 ct, len(da_genes), len(overlap_genes),
                 n_concordant, n_discordant,
                 100 * n_concordant / max(n_concordant + n_discordant, 1))

    overlap_df = pd.DataFrame(overlap_records)
    return overlap_df


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="ATAC-seq DA for MASL-vs-MASH and MASL-vs-NORMAL progression",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", default=DEFAULT_INPUT,
                        help="Label-transferred h5ad")
    parser.add_argument("--per-donor-dir", default=DEFAULT_PER_DONOR,
                        help="Per-donor h5ad directory with QC metrics")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT,
                        help="Output directory")
    parser.add_argument("--gtf", default=DEFAULT_GTF,
                        help="GENCODE GTF for TSS positions")
    parser.add_argument("--c2-degs", default=DEFAULT_C2_DEGS,
                        help="C2 NAFL-vs-NASH dream results CSV")
    parser.add_argument("--cell-type-col", default="cell_type",
                        help="Column name for cell type annotations")
    parser.add_argument("--min-detection", type=float, default=MIN_DETECTION,
                        help="Min detection rate in either group")
    parser.add_argument("--promoter-window", type=int, default=PROMOTER_WINDOW,
                        help="Promoter window in bp for peak-gene mapping")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    t_start = time.time()

    # ── Step 1: Load scATAC data ──────────────────────────────────────────
    log.info("=" * 70)
    log.info("141: ATAC-seq MASL-vs-MASH Differential Accessibility")
    log.info("=" * 70)
    log.info("Loading %s (backed)...", args.input)
    adata = ad.read_h5ad(args.input, backed="r")
    log.info("Shape: %d x %d", adata.n_obs, adata.n_vars)
    log.info("obs columns: %s", list(adata.obs.columns))

    obs = adata.obs.copy()
    all_names = list(adata.obs_names)
    peak_names = list(adata.var_names)

    # Fix non-unique obs_names if needed
    if not adata.obs_names.is_unique:
        n_dup = adata.n_obs - adata.obs_names.nunique()
        log.info("Fixing %d duplicate obs_names (in-memory)", n_dup)
        adata.obs_names_make_unique()
        obs.index = adata.obs_names
        all_names = list(adata.obs_names)

    # ── Step 2: Inspect conditions and cell types ─────────────────────────
    conditions = obs["condition"].unique().tolist()
    cell_types = obs[args.cell_type_col].unique().tolist()
    log.info("Conditions: %s", conditions)
    log.info("Cell types: %s", cell_types)

    for cond in conditions:
        n = int((obs["condition"] == cond).sum())
        log.info("  %s: %d cells", cond, n)

    for ct in sorted(cell_types):
        n = int((obs[args.cell_type_col] == ct).sum())
        log.info("  %s: %d cells", ct, n)

    # Validate expected conditions
    cond_upper = {str(c).upper(): c for c in conditions}
    masl_label = cond_upper.get("MASL")
    mash_label = cond_upper.get("MASH")
    normal_label = cond_upper.get("NORMAL")

    if not masl_label:
        log.error("MASL condition not found. Available: %s", conditions)
        sys.exit(1)
    if not mash_label:
        log.error("MASH condition not found. Available: %s", conditions)
        sys.exit(1)
    if not normal_label:
        log.error("NORMAL condition not found. Available: %s", conditions)
        sys.exit(1)

    log.info("Resolved labels: MASL=%s, MASH=%s, NORMAL=%s",
             masl_label, mash_label, normal_label)

    # ── Step 3: Recover QC covariates for all cells ───────────────────────
    covariates = recover_qc_covariates(args.per_donor_dir, set(all_names))

    # ── Step 4: Run DA — Contrast 1: MASL vs MASH ────────────────────────
    # Positive logFC = more accessible in MASH (progression direction)
    da_masl_mash = run_contrast_da(
        adata=adata,
        obs_df=obs,
        covariates=covariates,
        all_names=all_names,
        peak_names=peak_names,
        group1_labels=[mash_label],
        group2_labels=[masl_label],
        contrast_name="MASH_vs_MASL (progression)",
        cell_type_col=args.cell_type_col,
    )

    # ── Step 5: Run DA — Contrast 2: MASL vs NORMAL ──────────────────────
    # Positive logFC = more accessible in MASL (onset direction)
    da_masl_normal = run_contrast_da(
        adata=adata,
        obs_df=obs,
        covariates=covariates,
        all_names=all_names,
        peak_names=peak_names,
        group1_labels=[masl_label],
        group2_labels=[normal_label],
        contrast_name="MASL_vs_NORMAL (onset)",
        cell_type_col=args.cell_type_col,
    )

    # Close backed file
    adata.file.close()

    # ── Step 6: Save raw DA results ───────────────────────────────────────
    log.info("=" * 70)
    log.info("Saving DA results")

    out_masl_mash = os.path.join(args.output_dir, "atac_masl_vs_mash_da.csv")
    out_masl_normal = os.path.join(args.output_dir, "atac_masl_vs_normal_da.csv")

    if not da_masl_mash.empty:
        da_masl_mash.to_csv(out_masl_mash, index=False)
        log.info("  Saved MASL-vs-MASH DA: %s (%d rows)", out_masl_mash, len(da_masl_mash))

        # Summary
        for padj_col_name in ["adjusted p-value", "adjusted_p_value", "padj", "fdr"]:
            if padj_col_name in da_masl_mash.columns:
                for ct in da_masl_mash["cell_type"].unique():
                    ct_da = da_masl_mash[da_masl_mash["cell_type"] == ct]
                    n_sig = int((ct_da[padj_col_name] < 0.05).sum())
                    log.info("    %s: %d sig DA peaks (padj<0.05)", ct, n_sig)
                break
    else:
        log.warning("  No MASL-vs-MASH DA results to save")

    if not da_masl_normal.empty:
        da_masl_normal.to_csv(out_masl_normal, index=False)
        log.info("  Saved MASL-vs-NORMAL DA: %s (%d rows)", out_masl_normal, len(da_masl_normal))

        for padj_col_name in ["adjusted p-value", "adjusted_p_value", "padj", "fdr"]:
            if padj_col_name in da_masl_normal.columns:
                for ct in da_masl_normal["cell_type"].unique():
                    ct_da = da_masl_normal[da_masl_normal["cell_type"] == ct]
                    n_sig = int((ct_da[padj_col_name] < 0.05).sum())
                    log.info("    %s: %d sig DA peaks (padj<0.05)", ct, n_sig)
                break
    else:
        log.warning("  No MASL-vs-NORMAL DA results to save")

    # ── Step 7: Peak-to-gene annotation ───────────────────────────────────
    log.info("=" * 70)
    log.info("Annotating DA peaks with nearest gene TSS")
    tss_df = parse_tss_from_gtf(args.gtf)

    if not da_masl_mash.empty:
        da_masl_mash_annot = annotate_peaks_with_genes(
            da_masl_mash, tss_df, promoter_window=args.promoter_window)
        da_masl_mash_annot.to_csv(out_masl_mash, index=False)
        log.info("  Updated MASL-vs-MASH with gene annotations")
    else:
        da_masl_mash_annot = pd.DataFrame()

    if not da_masl_normal.empty:
        da_masl_normal_annot = annotate_peaks_with_genes(
            da_masl_normal, tss_df, promoter_window=args.promoter_window)
        da_masl_normal_annot.to_csv(out_masl_normal, index=False)
        log.info("  Updated MASL-vs-NORMAL with gene annotations")
    else:
        da_masl_normal_annot = pd.DataFrame()

    # ── Step 8: Overlap with C2 progression DEGs ──────────────────────────
    overlap_records = []

    if not da_masl_mash_annot.empty:
        overlap_mash = compute_progression_overlap(
            da_masl_mash_annot, args.c2_degs,
            padj_thresh=0.05, lfc_thresh=0.25)
        if not overlap_mash.empty:
            overlap_mash["contrast"] = "MASH_vs_MASL"
            overlap_records.append(overlap_mash)

    if not da_masl_normal_annot.empty:
        overlap_normal = compute_progression_overlap(
            da_masl_normal_annot, args.c2_degs,
            padj_thresh=0.05, lfc_thresh=0.25)
        if not overlap_normal.empty:
            overlap_normal["contrast"] = "MASL_vs_NORMAL"
            overlap_records.append(overlap_normal)

    if overlap_records:
        overlap_df = pd.concat(overlap_records, ignore_index=True)
        out_overlap = os.path.join(args.output_dir, "atac_progression_overlap.csv")
        overlap_df.to_csv(out_overlap, index=False)
        log.info("  Saved progression overlap: %s (%d rows)", out_overlap, len(overlap_df))
    else:
        log.warning("  No overlap results to save")

    # ── Summary ───────────────────────────────────────────────────────────
    elapsed = (time.time() - t_start) / 60
    log.info("=" * 70)
    log.info("141 COMPLETE in %.1f min", elapsed)
    log.info("Output files:")
    log.info("  %s", out_masl_mash)
    log.info("  %s", out_masl_normal)
    if overlap_records:
        log.info("  %s", os.path.join(args.output_dir, "atac_progression_overlap.csv"))


if __name__ == "__main__":
    main()
