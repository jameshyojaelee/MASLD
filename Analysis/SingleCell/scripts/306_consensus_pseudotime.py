#!/usr/bin/env python3
"""
306: Consensus Pseudotime — Align & Merge 3 Methods.

For each of 5 core cell types, aligns pseudotimes from Palantir (301),
DPT (301), and Monocle 3 (304) to a common [0,1] scale using isotonic
regression against disease condition severity, computes pairwise method
concordance, and produces a correlation-weighted consensus pseudotime.
Identifies robust trajectory genes (significant in >=2/3 methods) and
writes the consensus back into each cell-type h5ad.

Inputs (from results_gpu_v2/pseudotime/):
    - palantir_pseudotime_{celltype}.csv   (Script 301)
    - dpt_pseudotime_{celltype}.csv        (Script 301)
    - monocle3_pseudotime_{celltype}.csv   (Script 304)
    - pseudotime_corr_{celltype}.csv       (Script 301: gene-pseudotime Spearman)
    - tradeseq_startend_{celltype}.csv     (Script 304: tradeSeq start-vs-end)
    - {celltype}_subset.h5ad               (Script 300)

Outputs (to results_gpu_v2/pseudotime/):
    - consensus_pseudotime_all.csv         Per-cell consensus for all cell types
    - robust_trajectory_genes.csv          Genes significant in >=2/3 methods
    - method_concordance.csv               Pairwise Spearman between methods

Usage:
    sbatch run_pseudotime_pipeline.sh
"""

import os
import sys
import warnings
import logging
import gc

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SC_DIR = os.path.join(BASE, "Analysis/SingleCell")
RESULTS = os.path.join(SC_DIR, "results_gpu_v2")
PT_DIR = os.path.join(RESULTS, "pseudotime")

CELL_TYPES = [
    "Hepatocytes",
    "Macrophages",
    "Fibroblasts",
    "Endothelial_cells",
    "Cholangiocytes",
]

# ---------------------------------------------------------------------------
# Imports (deferred so startup logging appears promptly)
# ---------------------------------------------------------------------------
from sklearn.isotonic import IsotonicRegression
import scanpy as sc

# ---------------------------------------------------------------------------
# Disease condition ordering
# ---------------------------------------------------------------------------
CONDITION_ORDER = {
    "Healthy": 0,
    "NAFLD": 1,
    "MASLD": 1.5,
    "NASH": 2,
    "Cirrhotic": 3,
}


def condition_score(series):
    """Map condition labels to monotone numeric score."""
    return series.map(CONDITION_ORDER).astype(float)


# ---------------------------------------------------------------------------
# Load pseudotime from a single method CSV
# ---------------------------------------------------------------------------
def _load_pseudotime(path, pt_col):
    """Read a pseudotime CSV and return a Series indexed by cell barcode.

    Handles the two index conventions used across scripts:
      - Palantir / DPT: cell barcode is the DataFrame index (first unnamed col)
      - Monocle 3: cell barcode is in a column named 'cell'
    """
    if not os.path.exists(path):
        return None

    df = pd.read_csv(path)

    # Determine index
    if "cell" in df.columns:
        df = df.set_index("cell")
    elif "Unnamed: 0" in df.columns:
        df = df.set_index("Unnamed: 0")
    elif df.columns[0] not in (pt_col, "condition"):
        # First column is likely the barcode
        df = df.set_index(df.columns[0])

    if pt_col not in df.columns:
        # Try case-insensitive match
        candidates = [c for c in df.columns if "pseudotime" in c.lower()]
        if candidates:
            pt_col = candidates[0]
        else:
            log.warning("Column '%s' not found in %s (cols: %s)",
                        pt_col, path, list(df.columns))
            return None

    series = df[pt_col].astype(float)
    series.index = series.index.astype(str)

    # Carry condition forward if present
    if "condition" in df.columns:
        series.attrs["condition"] = df["condition"]

    return series


def load_all_pseudotimes(ct):
    """Load pseudotime from Palantir, DPT, and Monocle 3 for a cell type.

    Returns
    -------
    dict  method_name -> pd.Series (indexed by cell barcode)
    """
    methods = {
        "palantir": (
            os.path.join(PT_DIR, f"palantir_pseudotime_{ct}.csv"),
            "palantir_pseudotime",
        ),
        "dpt": (
            os.path.join(PT_DIR, f"dpt_pseudotime_{ct}.csv"),
            "dpt_pseudotime",
        ),
        "monocle3": (
            os.path.join(PT_DIR, f"monocle3_pseudotime_{ct}.csv"),
            "monocle3_pseudotime",
        ),
    }

    loaded = {}
    for name, (path, col) in methods.items():
        series = _load_pseudotime(path, col)
        if series is not None:
            n_valid = np.isfinite(series).sum()
            log.info("  %s: %d cells (%d valid)", name, len(series), n_valid)
            loaded[name] = series
        else:
            log.warning("  %s: not available", name)

    return loaded


# ---------------------------------------------------------------------------
# Align pseudotime to [0,1] using isotonic regression on condition score
# ---------------------------------------------------------------------------
def align_pseudotime(pt_series, condition_series):
    """Align a raw pseudotime to [0,1] via isotonic regression.

    The monotone target is the disease condition score.  Isotonic regression
    finds the best monotone non-decreasing fit of pseudotime to the condition
    score, then rescales the fitted values to [0,1].

    Parameters
    ----------
    pt_series : pd.Series
        Raw pseudotime values, indexed by cell barcode.
    condition_series : pd.Series
        Condition labels for the same cells.

    Returns
    -------
    pd.Series  Aligned pseudotime in [0,1], same index.
    """
    common = pt_series.index.intersection(condition_series.index)
    pt = pt_series.loc[common].values.astype(float)
    cond = condition_score(condition_series.loc[common]).values

    valid = np.isfinite(pt) & np.isfinite(cond)
    if valid.sum() < 30:
        log.warning("Too few valid cells (%d) for isotonic alignment — "
                    "falling back to min-max rescale", valid.sum())
        pt_out = pt_series.copy()
        lo, hi = np.nanmin(pt_out), np.nanmax(pt_out)
        if hi > lo:
            pt_out = (pt_out - lo) / (hi - lo)
        return pt_out

    # Fit isotonic regression: pseudotime -> condition score
    ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
    ir.fit(pt[valid], cond[valid])

    # Transform the full pseudotime vector through the fitted function
    all_pt = pt_series.values.astype(float)
    all_valid = np.isfinite(all_pt)
    aligned = np.full(len(all_pt), np.nan)
    aligned[all_valid] = ir.predict(all_pt[all_valid])

    # Rescale to [0, 1]
    lo, hi = np.nanmin(aligned), np.nanmax(aligned)
    if hi > lo:
        aligned = (aligned - lo) / (hi - lo)

    return pd.Series(aligned, index=pt_series.index, name=pt_series.name)


# ---------------------------------------------------------------------------
# Pairwise Spearman concordance
# ---------------------------------------------------------------------------
def pairwise_concordance(aligned_dict):
    """Compute pairwise Spearman correlations between aligned pseudotimes.

    Parameters
    ----------
    aligned_dict : dict  method_name -> pd.Series (aligned pseudotime)

    Returns
    -------
    list[dict]  rows for the concordance table
    """
    methods = sorted(aligned_dict.keys())
    rows = []
    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            m1, m2 = methods[i], methods[j]
            s1, s2 = aligned_dict[m1], aligned_dict[m2]
            common = s1.index.intersection(s2.index)
            v1 = s1.loc[common].values
            v2 = s2.loc[common].values
            valid = np.isfinite(v1) & np.isfinite(v2)
            n = valid.sum()
            if n > 30:
                rho, pval = stats.spearmanr(v1[valid], v2[valid])
            else:
                rho, pval = np.nan, np.nan
            rows.append({
                "method_1": m1,
                "method_2": m2,
                "spearman_rho": rho,
                "pval": pval,
                "n_cells": n,
            })
    return rows


# ---------------------------------------------------------------------------
# Consensus pseudotime: correlation-weighted average
# ---------------------------------------------------------------------------
def consensus_pseudotime(aligned_dict, concordance_rows):
    """Compute correlation-weighted average pseudotime.

    Each method's weight is its mean absolute Spearman correlation with the
    other methods.  This down-weights an outlier method while still including
    its signal.

    Parameters
    ----------
    aligned_dict : dict  method_name -> pd.Series (aligned [0,1])
    concordance_rows : list[dict]  from pairwise_concordance()

    Returns
    -------
    pd.Series  consensus pseudotime, union of all cell barcodes
    dict       method -> weight used
    """
    methods = sorted(aligned_dict.keys())

    # Compute per-method weight = mean |rho| with other methods
    rho_map = {}
    for row in concordance_rows:
        m1, m2, rho = row["method_1"], row["method_2"], row["spearman_rho"]
        if np.isfinite(rho):
            rho_map.setdefault(m1, []).append(abs(rho))
            rho_map.setdefault(m2, []).append(abs(rho))

    weights = {}
    for m in methods:
        if m in rho_map and len(rho_map[m]) > 0:
            weights[m] = np.mean(rho_map[m])
        else:
            weights[m] = 1.0  # fallback if only one method

    # Normalise weights
    w_sum = sum(weights.values())
    if w_sum > 0:
        weights = {m: w / w_sum for m, w in weights.items()}

    log.info("Consensus weights: %s",
             {m: f"{w:.3f}" for m, w in weights.items()})

    # Build aligned DataFrame
    all_idx = sorted(
        set().union(*(s.index for s in aligned_dict.values()))
    )
    df = pd.DataFrame(index=all_idx)
    for m, s in aligned_dict.items():
        df[m] = s.reindex(all_idx)

    # Weighted average (ignoring NaN per cell)
    weighted_sum = np.zeros(len(df))
    weight_total = np.zeros(len(df))
    for m in methods:
        vals = df[m].values.astype(float)
        valid = np.isfinite(vals)
        weighted_sum[valid] += weights[m] * vals[valid]
        weight_total[valid] += weights[m]

    consensus = np.where(weight_total > 0, weighted_sum / weight_total, np.nan)
    return pd.Series(consensus, index=df.index, name="consensus_pseudotime"), weights


# ---------------------------------------------------------------------------
# Robust trajectory genes (significant in >=2/3 methods)
# ---------------------------------------------------------------------------
def find_robust_genes(ct, padj_threshold=0.05):
    """Identify genes significant in >=2 of 3 pseudotime methods.

    Sources:
      1. pseudotime_corr_{ct}.csv — Palantir/DPT Spearman correlations (padj)
      2. tradeseq_startend_{ct}.csv — Monocle 3 tradeSeq startVsEnd (pvalue)

    Palantir and DPT share the same correlation file (Script 301 uses whichever
    was available), so we treat it as "Palantir/DPT" evidence.  For a gene to
    be called robust it must be significant in at least 2 of:
      (a) Palantir/DPT gene-pseudotime correlation (padj < threshold)
      (b) tradeSeq startVsEndTest (padj < threshold OR pvalue < threshold if
          no padj column)
      (c) Monocle 3 pseudotime correlation — re-computed here if monocle3
          pseudotime is available, otherwise not counted

    Returns
    -------
    pd.DataFrame  with columns: gene, cell_type, n_methods_sig, sig_methods,
                  corr_padj, tradeseq_pval, monocle3_corr_padj
    """
    gene_evidence = {}  # gene -> set of method labels

    # --- Source 1: Palantir/DPT gene-pseudotime correlation (Script 301) ---
    corr_path = os.path.join(PT_DIR, f"pseudotime_corr_{ct}.csv")
    corr_df = None
    if os.path.exists(corr_path):
        corr_df = pd.read_csv(corr_path)
        if "padj" in corr_df.columns:
            sig = corr_df.loc[corr_df["padj"] < padj_threshold, "gene"]
            for g in sig:
                gene_evidence.setdefault(g, set()).add("palantir_dpt_corr")
            log.info("  Palantir/DPT corr: %d significant genes", len(sig))
        else:
            log.warning("  Palantir/DPT corr: no padj column in %s", corr_path)
    else:
        log.warning("  Palantir/DPT corr file not found: %s", corr_path)

    # --- Source 2: tradeSeq startVsEndTest (Script 304) ---
    tradeseq_path = os.path.join(PT_DIR, f"tradeseq_startend_{ct}.csv")
    tradeseq_df = None
    if os.path.exists(tradeseq_path):
        tradeseq_df = pd.read_csv(tradeseq_path)
        # tradeSeq output has 'pvalue' and sometimes 'padj'/'FDR'
        pval_col = None
        for cand in ["padj", "FDR", "pvalue"]:
            if cand in tradeseq_df.columns:
                pval_col = cand
                break
        if pval_col is not None:
            sig = tradeseq_df.loc[tradeseq_df[pval_col] < padj_threshold, "gene"]
            for g in sig:
                gene_evidence.setdefault(g, set()).add("tradeseq_startend")
            log.info("  tradeSeq startEnd (%s): %d significant genes",
                     pval_col, len(sig))
        else:
            log.warning("  tradeSeq: no p-value column found (cols: %s)",
                        list(tradeseq_df.columns))
    else:
        log.warning("  tradeSeq file not found: %s", tradeseq_path)

    # --- Source 3: Monocle 3 gene-pseudotime correlation ---
    # Re-compute Spearman correlation of each gene against Monocle 3 pseudotime
    monocle3_path = os.path.join(PT_DIR, f"monocle3_pseudotime_{ct}.csv")
    h5ad_path = os.path.join(PT_DIR, f"{ct}_subset.h5ad")
    monocle3_corr = None

    if os.path.exists(monocle3_path) and os.path.exists(h5ad_path):
        try:
            m3_pt = _load_pseudotime(monocle3_path, "monocle3_pseudotime")
            adata = sc.read_h5ad(h5ad_path)

            common = m3_pt.index.intersection(adata.obs_names)
            pt_vals = m3_pt.loc[common].values.astype(float)
            valid_mask = np.isfinite(pt_vals)

            if valid_mask.sum() > 100:
                from scipy import sparse as sp_sparse
                X = adata[common].X
                if sp_sparse.issparse(X):
                    X = X.toarray()
                X_valid = X[valid_mask]
                pt_valid = pt_vals[valid_mask]

                genes = adata.var_names
                m3_rows = []
                for i in range(X_valid.shape[1]):
                    expr = X_valid[:, i]
                    if np.std(expr) == 0:
                        continue
                    rho, pval = stats.spearmanr(pt_valid, expr)
                    m3_rows.append({"gene": genes[i], "pval": pval})

                if m3_rows:
                    m3_df = pd.DataFrame(m3_rows)
                    from statsmodels.stats.multitest import multipletests
                    _, m3_df["padj"], _, _ = multipletests(
                        m3_df["pval"], method="fdr_bh"
                    )
                    monocle3_corr = m3_df.set_index("gene")["padj"]
                    sig = m3_df.loc[m3_df["padj"] < padj_threshold, "gene"]
                    for g in sig:
                        gene_evidence.setdefault(g, set()).add("monocle3_corr")
                    log.info("  Monocle 3 corr: %d significant genes", len(sig))

            del adata
            gc.collect()

        except Exception as e:
            log.warning("  Monocle 3 gene correlation failed: %s", e)
    else:
        log.info("  Monocle 3 corr: skipped (missing files)")

    # --- Assemble robust gene table ---
    if not gene_evidence:
        log.warning("  No gene evidence collected for %s", ct)
        return pd.DataFrame()

    rows = []
    for gene, methods_set in gene_evidence.items():
        n_sig = len(methods_set)
        # Retrieve per-source statistics
        corr_padj = np.nan
        if corr_df is not None and "padj" in corr_df.columns:
            match = corr_df.loc[corr_df["gene"] == gene, "padj"]
            if len(match) > 0:
                corr_padj = match.values[0]

        ts_pval = np.nan
        if tradeseq_df is not None:
            pval_col_ts = None
            for cand in ["padj", "FDR", "pvalue"]:
                if cand in tradeseq_df.columns:
                    pval_col_ts = cand
                    break
            if pval_col_ts is not None:
                match = tradeseq_df.loc[tradeseq_df["gene"] == gene, pval_col_ts]
                if len(match) > 0:
                    ts_pval = match.values[0]

        m3_padj = np.nan
        if monocle3_corr is not None and gene in monocle3_corr.index:
            m3_padj = monocle3_corr.loc[gene]

        rows.append({
            "gene": gene,
            "cell_type": ct,
            "n_methods_sig": n_sig,
            "sig_methods": ";".join(sorted(methods_set)),
            "corr_padj": corr_padj,
            "tradeseq_pval": ts_pval,
            "monocle3_corr_padj": m3_padj,
        })

    df = pd.DataFrame(rows)
    robust = df[df["n_methods_sig"] >= 2].sort_values(
        "n_methods_sig", ascending=False
    )
    log.info("  Robust trajectory genes (>=2/3): %d / %d total",
             len(robust), len(df))

    return robust


# ---------------------------------------------------------------------------
# Write consensus pseudotime into cell-type h5ad
# ---------------------------------------------------------------------------
def update_h5ad(ct, consensus_series):
    """Add consensus_pseudotime as an obs column to the cell-type h5ad."""
    h5ad_path = os.path.join(PT_DIR, f"{ct}_subset.h5ad")
    if not os.path.exists(h5ad_path):
        log.warning("h5ad not found: %s — skipping update", h5ad_path)
        return

    try:
        adata = sc.read_h5ad(h5ad_path)
        adata.obs["consensus_pseudotime"] = consensus_series.reindex(
            adata.obs_names
        ).values
        n_assigned = np.isfinite(adata.obs["consensus_pseudotime"]).sum()
        log.info("  Updated h5ad with consensus_pseudotime: %d/%d cells",
                 n_assigned, len(adata))
        adata.write_h5ad(h5ad_path)
        del adata
        gc.collect()
    except Exception as e:
        log.error("  Failed to update h5ad for %s: %s", ct, e)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("306: Consensus Pseudotime")
    log.info("=" * 60)

    all_consensus = []
    all_concordance = []
    all_robust_genes = []

    for ct in CELL_TYPES:
        log.info("\n" + "=" * 50)
        log.info("Processing: %s", ct)
        log.info("=" * 50)

        # --- Load pseudotimes from all 3 methods ---
        raw_pt = load_all_pseudotimes(ct)
        if len(raw_pt) == 0:
            log.warning("No pseudotime data for %s — skipping", ct)
            continue

        # --- Obtain condition labels (from any available CSV) ---
        condition_series = None
        for name, series in raw_pt.items():
            if hasattr(series, "attrs") and "condition" in series.attrs:
                condition_series = series.attrs["condition"]
                break

        if condition_series is None:
            # Fallback: load from h5ad
            h5ad_path = os.path.join(PT_DIR, f"{ct}_subset.h5ad")
            if os.path.exists(h5ad_path):
                try:
                    adata_tmp = sc.read_h5ad(h5ad_path, backed="r")
                    condition_series = adata_tmp.obs["condition"]
                    condition_series.index = condition_series.index.astype(str)
                    del adata_tmp
                except Exception as e:
                    log.warning("Could not load condition from h5ad: %s", e)

        if condition_series is None:
            log.warning("No condition labels found for %s — "
                        "using min-max rescale instead of isotonic", ct)

        # --- Align each method to [0,1] via isotonic regression ---
        aligned = {}
        for name, series in raw_pt.items():
            if condition_series is not None:
                aligned[name] = align_pseudotime(series, condition_series)
            else:
                # Simple min-max rescale fallback
                vals = series.values.astype(float)
                lo, hi = np.nanmin(vals), np.nanmax(vals)
                if hi > lo:
                    rescaled = (vals - lo) / (hi - lo)
                else:
                    rescaled = vals
                aligned[name] = pd.Series(
                    rescaled, index=series.index, name=series.name
                )
            log.info("  Aligned %s: mean=%.3f, std=%.3f",
                     name,
                     np.nanmean(aligned[name]),
                     np.nanstd(aligned[name]))

        # --- Pairwise concordance ---
        conc_rows = pairwise_concordance(aligned)
        for row in conc_rows:
            row["cell_type"] = ct
            log.info("  %s vs %s: rho=%.3f (n=%d)",
                     row["method_1"], row["method_2"],
                     row["spearman_rho"], row["n_cells"])
        all_concordance.extend(conc_rows)

        # --- Consensus pseudotime ---
        cons, weights = consensus_pseudotime(aligned, conc_rows)
        n_valid = np.isfinite(cons).sum()
        log.info("Consensus pseudotime: %d valid cells, mean=%.3f, std=%.3f",
                 n_valid, np.nanmean(cons), np.nanstd(cons))

        # Build output DataFrame for this cell type
        ct_df = pd.DataFrame(index=cons.index)
        ct_df["cell_type"] = ct
        ct_df["consensus_pseudotime"] = cons.values
        for m, s in aligned.items():
            ct_df[f"aligned_{m}"] = s.reindex(cons.index).values
        ct_df["n_methods"] = ct_df[
            [c for c in ct_df.columns if c.startswith("aligned_")]
        ].notna().sum(axis=1)
        all_consensus.append(ct_df)

        # --- Update h5ad ---
        update_h5ad(ct, cons)

        # --- Robust trajectory genes ---
        try:
            robust_df = find_robust_genes(ct)
            if len(robust_df) > 0:
                all_robust_genes.append(robust_df)
        except Exception as e:
            log.error("Robust gene identification failed for %s: %s", ct, e)

    # ===================================================================
    # Save combined outputs
    # ===================================================================
    log.info("\n" + "=" * 50)
    log.info("Saving combined outputs")
    log.info("=" * 50)

    # 1. Consensus pseudotime (all cell types)
    if all_consensus:
        full_cons = pd.concat(all_consensus, axis=0)
        out_path = os.path.join(PT_DIR, "consensus_pseudotime_all.csv")
        full_cons.to_csv(out_path)
        log.info("Saved consensus pseudotime: %s -> %s", full_cons.shape, out_path)

        # Per-cell-type summary
        summary = full_cons.groupby("cell_type").agg(
            n_cells=("consensus_pseudotime", "count"),
            n_valid=("consensus_pseudotime", lambda x: np.isfinite(x).sum()),
            mean_pt=("consensus_pseudotime", "mean"),
            n_3methods=("n_methods", lambda x: (x == 3).sum()),
            n_2methods=("n_methods", lambda x: (x == 2).sum()),
            n_1method=("n_methods", lambda x: (x == 1).sum()),
        )
        log.info("Per-cell-type summary:\n%s", summary.to_string())
    else:
        log.error("No consensus pseudotime produced!")

    # 2. Method concordance
    if all_concordance:
        conc_df = pd.DataFrame(all_concordance)
        out_path = os.path.join(PT_DIR, "method_concordance.csv")
        conc_df.to_csv(out_path, index=False)
        log.info("Saved method concordance: %s -> %s", conc_df.shape, out_path)

        # Summary
        mean_rho = conc_df.groupby("cell_type")["spearman_rho"].mean()
        log.info("Mean concordance per cell type:\n%s", mean_rho.to_string())
    else:
        log.warning("No concordance data produced")

    # 3. Robust trajectory genes
    if all_robust_genes:
        robust_all = pd.concat(all_robust_genes, axis=0, ignore_index=True)
        out_path = os.path.join(PT_DIR, "robust_trajectory_genes.csv")
        robust_all.to_csv(out_path, index=False)
        log.info("Saved robust trajectory genes: %s -> %s",
                 robust_all.shape, out_path)

        # Summary per cell type
        ct_summary = robust_all.groupby("cell_type").agg(
            n_robust=("gene", "count"),
            n_all_three=("n_methods_sig", lambda x: (x == 3).sum()),
        )
        log.info("Robust genes per cell type:\n%s", ct_summary.to_string())
    else:
        log.warning("No robust trajectory genes identified")

    log.info("\n=== 306: Consensus Pseudotime COMPLETE ===")


if __name__ == "__main__":
    main()
