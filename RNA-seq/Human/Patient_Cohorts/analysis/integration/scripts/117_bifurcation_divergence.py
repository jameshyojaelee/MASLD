#!/usr/bin/env python3
"""
117_bifurcation_divergence.py
Bifurcation analysis and subtype divergence.

WHY do some patients progress to cirrhosis while others don't?
Uses NMF k=2 subtypes (S1=909, S2=345 female-enriched) to identify:
  1. The pseudotime window where S1 and S2 trajectories diverge
  2. Genes whose expression trajectory differs between subtypes (divergence genes)
  3. Fate probabilities: probability of reaching F4 vs resolving
  4. Clinical correlates of subtypes (sex, age, BMI, NAS components)
  5. SCENIC+ regulon activity in the bifurcation window

Input:
  - results/staging_classifier/nas_embeddings_all_samples.csv
  - results/staging_classifier/modeling_metadata.csv
  - results/staging_classifier/prepared_data.h5
  - results/subtypes/nmf_assignments.csv (NMF k=2 subtypes)

Output (all to results/progression/):
  - bifurcation_analysis.csv         (bifurcation window, metrics)
  - divergence_genes.csv             (full-cohort S1-vs-S2 ranking; DESCRIPTIVE
                                      ONLY — selected using every sample, so it
                                      MUST NOT seed the predictor feature panel)
  - divergence_genes_per_fold.csv    (LEAKAGE-FREE: S1-vs-S2 divergence ranking
                                      recomputed once per LOCO outer fold using
                                      TRAINING cohorts only; long format keyed by
                                      `held_out_fold`. Canonical source for the
                                      `div_*` prognosis feature panel.)
  - fate_probabilities.csv           (per-sample probability of progression)
  - bifurcation_clinical.csv         (clinical correlates of subtypes)
  - bifurcation_summary.csv          (overview statistics)

Leakage fix (2026-06-20, mega-review blocker A — S2-fate selection leak):
  The `div_*` panel feeding the LOCO-CV prognosis predictors (Scripts 150 ->
  160 -> 161/183/183v2) was previously selected from the FULL-cohort S1-vs-S2
  contrast, then frozen inside LOCO-CV — selection saw held-out test folds. The
  per-fold ranking (`divergence_genes_per_fold.csv`) recomputes the S1-vs-S2
  divergence statistic separately for each LOCO fold using only the training
  cohorts, so the held-out cohort never contributes to its own fold's panel.

SLURM: cpu partition, 8 CPUs, 64G RAM, 48h
Env:   micromamba activate rapids_singlecell
"""

import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["RAPIDS_NO_INITIALIZE"] = "1"
os.environ["CUDF_PANDAS"] = "0"

import sys
import time
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import (
    spearmanr, mannwhitneyu, chi2_contingency,
    kruskal, pearsonr, ttest_ind, fisher_exact
)
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score
import h5py

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEED = 42
np.random.seed(SEED)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
STAG = os.path.join(RDIR, "staging_classifier")
OUTDIR = os.path.join(RDIR, "progression")
os.makedirs(OUTDIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(OUTDIR, "117_bifurcation.log")),
    ],
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compute_knn_transition_matrix(embeddings, k=30):
    """Build transition matrix from kNN graph."""
    N = embeddings.shape[0]
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    nn.fit(embeddings)
    distances, indices = nn.kneighbors(embeddings)

    sigma = np.maximum(distances[:, -1], 1e-10)
    W = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        for j_idx in range(k):
            j = indices[i, j_idx]
            d = distances[i, j_idx]
            w = np.exp(-d**2 / (2 * sigma[i] * sigma[j]))
            W[i, j] = w
            W[j, i] = max(W[j, i], w)

    row_sums = np.maximum(W.sum(axis=1, keepdims=True), 1e-10)
    T_mat = W / row_sums
    return T_mat


def compute_fate_probabilities(T_mat, terminal_mask, n_iter=100):
    """
    Compute fate (absorption) probabilities.
    terminal_mask: boolean mask for terminal/absorbing states.
    Returns: probability of being absorbed by terminal states.
    """
    N = T_mat.shape[0]
    n_terminal = terminal_mask.sum()
    n_transient = N - n_terminal

    if n_terminal == 0:
        return np.ones(N) * 0.5

    terminal_idx = np.where(terminal_mask)[0]
    transient_idx = np.where(~terminal_mask)[0]

    # Iterative approach: simulate absorption
    # P(i reaches terminal) = sum_j T[i,j] * P(j reaches terminal) for transient j
    # + sum_j T[i,j] for terminal j
    fate_prob = np.zeros(N)
    fate_prob[terminal_mask] = 1.0

    for _ in range(n_iter):
        new_fate = np.zeros(N)
        new_fate[terminal_mask] = 1.0
        for i in transient_idx:
            new_fate[i] = T_mat[i, :] @ fate_prob
        if np.max(np.abs(new_fate - fate_prob)) < 1e-6:
            break
        fate_prob = new_fate

    return fate_prob


def sliding_window_divergence(pt, subtype, values, window_size=0.15, step=0.05):
    """
    Compute divergence between subtypes in sliding windows along pseudotime.
    Returns: window centers, divergence scores (absolute mean difference).
    """
    centers = np.arange(0, 1 - window_size + step, step)
    divergences = []

    for c in centers:
        mask = (pt >= c) & (pt < c + window_size)
        if mask.sum() < 10:
            divergences.append(np.nan)
            continue

        s1_vals = values[mask & (subtype == "S1")]
        s2_vals = values[mask & (subtype == "S2")]

        if len(s1_vals) < 3 or len(s2_vals) < 3:
            divergences.append(np.nan)
            continue

        # Divergence = |mean(S1) - mean(S2)| / pooled_sd
        pooled_sd = np.sqrt((np.var(s1_vals) + np.var(s2_vals)) / 2)
        if pooled_sd < 1e-10:
            divergences.append(0)
        else:
            divergences.append(np.abs(np.mean(s1_vals) - np.mean(s2_vals)) / pooled_sd)

    return centers + window_size / 2, np.array(divergences)


def detect_divergence_genes(expr, gene_names, s1_mask, s2_mask, min_per_group=5):
    """
    S1-vs-S2 divergence gene detection (Mann-Whitney + Cohen's d) on a *given*
    sample subset defined by s1_mask / s2_mask.

    This is the per-gene selection statistic used to build the `div_*` prognosis
    feature panel. It is factored out so the SAME logic can be applied either to
    the full cohort (descriptive bifurcation report) or to a single LOCO
    training split (leakage-free per-fold feature selection).

    Returns a DataFrame sorted by p-value with columns:
      gene, pval, cohens_d, mean_S1, mean_S2, mean_diff, padj
    Empty DataFrame if either group has < min_per_group samples.
    """
    if s1_mask.sum() < min_per_group or s2_mask.sum() < min_per_group:
        return pd.DataFrame()

    results = []
    for g_idx, gene in enumerate(gene_names):
        s1_expr = expr[s1_mask, g_idx]
        s2_expr = expr[s2_mask, g_idx]

        # Mann-Whitney U test
        try:
            stat, pval = mannwhitneyu(s1_expr, s2_expr, alternative="two-sided")
        except ValueError:
            continue

        # Effect size (Cohen's d)
        pooled_sd = np.sqrt((np.var(s1_expr) + np.var(s2_expr)) / 2)
        if pooled_sd > 1e-10:
            cohens_d = (np.mean(s1_expr) - np.mean(s2_expr)) / pooled_sd
        else:
            cohens_d = 0

        results.append({
            "gene": gene,
            "pval": pval,
            "cohens_d": cohens_d,
            "mean_S1": np.mean(s1_expr),
            "mean_S2": np.mean(s2_expr),
            "mean_diff": np.mean(s1_expr) - np.mean(s2_expr),
        })

    if not results:
        return pd.DataFrame()

    out = pd.DataFrame(results)

    # BH correction
    from scipy.stats import false_discovery_control
    try:
        out["padj"] = false_discovery_control(out["pval"].values, method="bh")
    except Exception:
        # Manual BH
        n = len(out)
        ranks = out["pval"].rank()
        out["padj"] = (out["pval"] * n / ranks).clip(upper=1.0)

    return out.sort_values("pval").reset_index(drop=True)


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    t0 = time.time()
    log.info("=== 117: Bifurcation & Subtype Divergence ===")

    # -------------------------------------------------------------------
    # 1. Load data
    # -------------------------------------------------------------------
    log.info("Loading data...")

    meta = pd.read_csv(os.path.join(STAG, "modeling_metadata.csv"))
    log.info(f"  Metadata: {len(meta)} samples")

    # NAS-VAE embeddings
    emb_path = os.path.join(STAG, "nas_embeddings_all_samples.csv")
    if not os.path.exists(emb_path):
        emb_path = os.path.join(STAG, "embeddings_all_samples.csv")
    emb_df = pd.read_csv(emb_path, index_col=0)
    log.info(f"  Embeddings: {emb_df.shape}")

    # NMF subtype assignments
    nmf_path = os.path.join(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")
    if not os.path.exists(nmf_path):
        log.error(f"NMF assignments not found at {nmf_path}")
        # Try alternative path
        nmf_path = os.path.join(RDIR, "subtypes/nmf_assignments.csv")

    if os.path.exists(nmf_path):
        nmf = pd.read_csv(nmf_path)
        log.info(f"  NMF assignments: {len(nmf)} samples")
        # Determine subtype column name
        subtype_col = None
        for col in ["nmf_subtype", "subtype", "cluster", "nmf_cluster"]:
            if col in nmf.columns:
                subtype_col = col
                break
        if subtype_col is None:
            # Look for k=2 assignment column
            for col in nmf.columns:
                if "k2" in col.lower() or "subtype" in col.lower():
                    subtype_col = col
                    break
        if subtype_col is None:
            subtype_col = nmf.columns[-1]  # Last column as fallback
        log.info(f"  Using subtype column: {subtype_col}")
    else:
        log.warning("NMF assignments not found. Will derive subtypes from embeddings.")
        nmf = None

    # Pseudotime (from Script 94 or 116)
    pt_path = os.path.join(OUTDIR, "consensus_pseudotime.csv")
    if not os.path.exists(pt_path):
        pt_path = os.path.join(STAG, "pseudotime_calibration.csv")
    if os.path.exists(pt_path):
        pt_df = pd.read_csv(pt_path)
        log.info(f"  Pseudotime: {len(pt_df)} samples")
    else:
        log.warning("No pseudotime file found. Will compute from embeddings.")
        pt_df = None

    # Expression data
    h5_path = os.path.join(STAG, "prepared_data.h5")
    has_expr = os.path.exists(h5_path)
    # h5 sample order matches modeling_metadata.csv — use meta sample_id for alignment
    h5_sample_order = list(meta["sample_id"])  # before any subsetting
    if has_expr:
        with h5py.File(h5_path, "r") as h5:
            if "rank_expression" in h5:
                expr_full = h5["rank_expression"][:].T  # (genes x samples) -> (samples x genes)
                gene_names = [g.decode() if isinstance(g, bytes) else g
                              for g in h5["gene_names"][:]]
                log.info(f"  Expression: {expr_full.shape} (samples x genes)")
            elif "raw_logcpm" in h5:
                expr_full = h5["raw_logcpm"][:].T
                gene_names = [g.decode() if isinstance(g, bytes) else g
                              for g in h5["gene_names"][:]]
                log.info(f"  Expression (logcpm): {expr_full.shape} (samples x genes)")
            elif "zscore_expression" in h5:
                expr_full = h5["zscore_expression"][:].T
                gene_names = [g.decode() if isinstance(g, bytes) else g
                              for g in h5["gene_names"][:]]
                log.info(f"  Expression (zscore): {expr_full.shape} (samples x genes)")
            else:
                has_expr = False

    # -------------------------------------------------------------------
    # 2. Align samples and prepare data
    # -------------------------------------------------------------------
    log.info("\nAligning samples...")
    common = set(meta["sample_id"]) & set(emb_df.index)
    if nmf is not None:
        nmf_id_col = "sample_id" if "sample_id" in nmf.columns else nmf.columns[0]
        common &= set(nmf[nmf_id_col])
    common = sorted(common)
    log.info(f"  Common samples: {len(common)}")

    meta = meta.set_index("sample_id").loc[common].reset_index()
    embeddings = emb_df.loc[common].values

    # Align expression to common samples
    if has_expr:
        h5_sid_to_idx = {s: i for i, s in enumerate(h5_sample_order)}
        common_expr_idx = [h5_sid_to_idx[s] for s in common if s in h5_sid_to_idx]
        expr = expr_full[common_expr_idx, :]
        log.info(f"  Expression aligned: {expr.shape}")

    if nmf is not None:
        nmf_aligned = nmf.set_index(nmf_id_col).loc[common]
        subtypes = nmf_aligned[subtype_col].values.astype(str)
        # Standardize naming
        unique_st = sorted(set(subtypes))
        log.info(f"  Subtypes: {unique_st}")
        # Map to S1/S2 if needed
        if len(unique_st) == 2:
            st_map = {unique_st[0]: "S1", unique_st[1]: "S2"}
            subtypes = np.array([st_map[s] for s in subtypes])
    else:
        # Simple k-means clustering as fallback
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=2, random_state=SEED, n_init=10)
        labels = km.fit_predict(embeddings)
        # Assign larger cluster as S1
        if (labels == 0).sum() >= (labels == 1).sum():
            subtypes = np.where(labels == 0, "S1", "S2")
        else:
            subtypes = np.where(labels == 1, "S1", "S2")

    meta["subtype"] = subtypes
    n_s1 = (subtypes == "S1").sum()
    n_s2 = (subtypes == "S2").sum()
    log.info(f"  S1: {n_s1}, S2: {n_s2}")

    # Get pseudotime
    if pt_df is not None:
        pt_id_col = "sample_id" if "sample_id" in pt_df.columns else pt_df.columns[0]
        pt_aligned = pt_df.set_index(pt_id_col).loc[common]
        if "pseudotime_consensus" in pt_aligned.columns:
            pseudotime = pt_aligned["pseudotime_consensus"].values
        elif "pseudotime" in pt_aligned.columns:
            pseudotime = pt_aligned["pseudotime"].values
        else:
            pseudotime = pt_aligned.iloc[:, 0].values
    else:
        # Compute simple pseudotime from diffusion map
        from sklearn.neighbors import NearestNeighbors
        root_mask = (meta["group_binary"] == "Control").values
        nn = NearestNeighbors(n_neighbors=30, metric="euclidean", n_jobs=-1)
        nn.fit(embeddings)
        distances, indices = nn.kneighbors(embeddings)
        sigma = np.maximum(distances[:, -1], 1e-10)
        W = np.zeros((len(common), len(common)))
        for i in range(len(common)):
            for j_idx in range(30):
                j = indices[i, j_idx]
                d = distances[i, j_idx]
                w = np.exp(-d**2 / (2 * sigma[i] * sigma[j]))
                W[i, j] = w
                W[j, i] = max(W[j, i], w)
        row_sums = np.maximum(W.sum(axis=1, keepdims=True), 1e-10)
        T_mat_raw = W / row_sums
        eigenvalues, eigenvectors = eigsh(csr_matrix(W), k=11, which="LM")
        dc = eigenvectors[:, 1:11] * eigenvalues[1:11]
        root_centroid = dc[root_mask].mean(axis=0)
        pseudotime = np.sqrt(((dc - root_centroid)**2).sum(axis=1))
        pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())

    meta["pseudotime"] = pseudotime

    # -------------------------------------------------------------------
    # 3. Find bifurcation window
    # -------------------------------------------------------------------
    log.info("\n=== Bifurcation Window Detection ===")

    # Compute S1/S2 divergence along pseudotime using embedding distance
    # At each pseudotime window, measure how different S1 and S2 are
    # in embedding space (effect size = Cohen's d)
    pt_centers, divergences = sliding_window_divergence(
        pseudotime, subtypes,
        np.linalg.norm(embeddings - embeddings.mean(axis=0), axis=1),
        window_size=0.15, step=0.03
    )

    # Also compute divergence per embedding dimension
    multi_dim_div = np.zeros_like(pt_centers)
    for dim in range(min(10, embeddings.shape[1])):
        _, dim_div = sliding_window_divergence(
            pseudotime, subtypes, embeddings[:, dim],
            window_size=0.15, step=0.03
        )
        multi_dim_div += np.nan_to_num(dim_div)
    multi_dim_div /= min(10, embeddings.shape[1])

    # Find peak divergence
    valid = ~np.isnan(multi_dim_div)
    if valid.any():
        peak_idx = np.nanargmax(multi_dim_div)
        bifurcation_pt = pt_centers[peak_idx]
        peak_div = multi_dim_div[peak_idx]
        log.info(f"  Peak divergence at pseudotime={bifurcation_pt:.3f} (div={peak_div:.3f})")

        # Map to approximate clinical stage
        fib = meta["fibrosis_stage"].values.astype(float)
        fib_mask = ~np.isnan(fib)
        if fib_mask.any():
            near_bif = np.abs(pseudotime - bifurcation_pt) < 0.1
            if near_bif.any() and fib_mask[near_bif].any():
                approx_fib = np.nanmedian(fib[near_bif & fib_mask])
                log.info(f"  Approximate fibrosis stage at bifurcation: F{approx_fib:.1f}")
    else:
        bifurcation_pt = 0.5
        peak_div = 0
        log.warning("  No valid divergence computed")

    # Bifurcation window: pseudotime range where divergence > 50% of peak
    threshold = peak_div * 0.5
    bif_window = pt_centers[multi_dim_div > threshold] if valid.any() else np.array([0.4, 0.6])
    bif_start = bif_window.min() if len(bif_window) > 0 else 0.4
    bif_end = bif_window.max() if len(bif_window) > 0 else 0.6
    log.info(f"  Bifurcation window: [{bif_start:.3f}, {bif_end:.3f}]")

    # -------------------------------------------------------------------
    # 4. Divergence gene detection
    # -------------------------------------------------------------------
    log.info("\n=== Divergence Gene Detection ===")

    in_bif_window = (pseudotime >= bif_start) & (pseudotime <= bif_end)

    if has_expr:
        # ---- Full-cohort divergence (DESCRIPTIVE ONLY) -------------------
        # NOTE (leakage fix 2026-06-20): this full-cohort S1-vs-S2 ranking is
        # retained ONLY for the descriptive bifurcation narrative. It selects
        # genes using EVERY sample, including those that later land in held-out
        # LOCO test folds, so it MUST NOT be used to build the `div_*`
        # prognosis feature panel. The leakage-free selection that the
        # predictor consumes is `divergence_genes_per_fold.csv`, computed
        # below using training cohorts only. (Roadmap A: 117 S2-fate leak.)
        s1_bif = in_bif_window & (subtypes == "S1")
        s2_bif = in_bif_window & (subtypes == "S2")

        log.info(f"  Samples in bifurcation window: S1={s1_bif.sum()}, S2={s2_bif.sum()}")

        div_df = detect_divergence_genes(expr, gene_names, s1_bif, s2_bif, min_per_group=5)

        if len(div_df) > 0:
            n_sig = (div_df["padj"] < 0.05).sum()
            n_strong = ((div_df["padj"] < 0.05) & (div_df["cohens_d"].abs() > 0.5)).sum()
            log.info(f"  [DESCRIPTIVE, full-cohort] Divergence genes (padj<0.05): {n_sig}")
            log.info(f"  [DESCRIPTIVE, full-cohort] Strong divergence (padj<0.05, |d|>0.5): {n_strong}")
            log.info(f"\n  Top 20 divergence genes (full-cohort, descriptive):")
            for _, row in div_df.head(20).iterrows():
                log.info(f"    {row['gene']}: d={row['cohens_d']:.2f}, padj={row['padj']:.2e}")
        else:
            log.warning("  Not enough samples in bifurcation window for divergence analysis")

    # -------------------------------------------------------------------
    # 4b. Fold-internal divergence selection (LEAKAGE-FREE)
    # -------------------------------------------------------------------
    # The `div_*` prognosis feature panel (Scripts 150 -> 160 -> 161/183/183v2)
    # must be selected WITHOUT seeing the held-out LOCO test cohort, otherwise
    # the S2-fate AUROC is inflated by selection leakage (mega-review blocker A:
    # 117_bifurcation_divergence.py:404-438). Here we re-run the divergence
    # detection once per LOCO outer fold, using ONLY the training cohorts
    # (every fold except the held-out one), and emit a long-format ranking that
    # downstream feature engineering (150) consumes per fold. The held-out
    # cohort's samples never contribute to its own fold's selection statistic.
    #
    # Fold definition mirrors Script 161: `loco_fold_fibrosis` in
    # modeling_metadata.csv, where each non-"excluded"/non-NA value names one
    # leave-one-cohort-out fold.
    div_fold_df = pd.DataFrame()
    if has_expr:
        log.info("\n=== Fold-internal Divergence Selection (leakage-free) ===")

        # Align loco_fold_fibrosis to the `common` sample order. `meta` is
        # already reset to `common` (set_index/loc/reset_index above), so its
        # row order matches `expr`, `subtypes`, `pseudotime`.
        if "loco_fold_fibrosis" in meta.columns:
            fold_labels = meta["loco_fold_fibrosis"].astype("object").values
        else:
            log.warning(
                "  'loco_fold_fibrosis' not in modeling_metadata; cannot do "
                "fold-internal selection. div_*_per_fold.csv will be empty."
            )
            fold_labels = np.array([None] * len(meta), dtype=object)

        # Valid outer folds = named cohorts (drop 'excluded' / NaN)
        fold_names = sorted(
            {
                f for f in fold_labels
                if isinstance(f, str) and f not in ("excluded", "nan", "")
            }
        )
        log.info(f"  LOCO outer folds: {fold_names}")

        per_fold_rows = []
        for held_out in fold_names:
            # Training split = all samples whose fold is a valid cohort AND is
            # NOT the held-out cohort. 'excluded'/NaN samples (no F-stage) are
            # left out of selection entirely, matching 161's valid-fold filter.
            is_valid_fold = np.array(
                [isinstance(f, str) and f in fold_names for f in fold_labels]
            )
            train_mask_fold = is_valid_fold & (fold_labels != held_out)

            s1_train = in_bif_window & train_mask_fold & (subtypes == "S1")
            s2_train = in_bif_window & train_mask_fold & (subtypes == "S2")

            fold_div = detect_divergence_genes(
                expr, gene_names, s1_train, s2_train, min_per_group=5
            )
            if len(fold_div) == 0:
                log.warning(
                    f"  Fold held-out={held_out}: insufficient S1/S2 in window "
                    f"(S1={s1_train.sum()}, S2={s2_train.sum()}); skipped."
                )
                continue

            fold_div = fold_div.copy()
            fold_div["held_out_fold"] = held_out
            fold_div["n_train_S1"] = int(s1_train.sum())
            fold_div["n_train_S2"] = int(s2_train.sum())
            fold_div["abs_cohens_d"] = fold_div["cohens_d"].abs()
            # Per-fold rank by |Cohen's d| (the statistic 150 uses to pick top-N)
            fold_div["rank"] = (
                fold_div["abs_cohens_d"].rank(ascending=False, method="first").astype(int)
            )
            n_sig = (fold_div["padj"] < 0.05).sum()
            log.info(
                f"  Fold held-out={held_out}: train S1={s1_train.sum()}, "
                f"S2={s2_train.sum()}; {len(fold_div)} genes, "
                f"{n_sig} sig (padj<0.05)"
            )
            per_fold_rows.append(fold_div)

        if per_fold_rows:
            div_fold_df = pd.concat(per_fold_rows, ignore_index=True)
            log.info(
                f"  Assembled per-fold divergence table: "
                f"{div_fold_df['held_out_fold'].nunique()} folds x "
                f"{div_fold_df.groupby('held_out_fold').size().mean():.0f} genes (mean)"
            )
        else:
            log.warning("  No fold produced a valid divergence ranking.")

    # -------------------------------------------------------------------
    # 5. Fate probabilities
    # -------------------------------------------------------------------
    log.info("\n=== Fate Probabilities ===")

    # Build transition matrix from embeddings
    T_mat = compute_knn_transition_matrix(embeddings, k=30)

    # Terminal states: F4 (progressed) and F0 (resolved/healthy)
    fib = meta["fibrosis_stage"].values
    terminal_f4 = (fib == 4)
    terminal_f0 = (fib == 0) | (meta["group_binary"].values == "Control")

    if terminal_f4.sum() > 0 and terminal_f0.sum() > 0:
        # Probability of reaching F4
        fate_f4 = compute_fate_probabilities(T_mat, terminal_f4, n_iter=100)
        # Probability of reaching F0/healthy
        fate_f0 = compute_fate_probabilities(T_mat, terminal_f0, n_iter=100)

        # Normalize
        total = fate_f4 + fate_f0
        total = np.maximum(total, 1e-10)
        fate_f4_norm = fate_f4 / total
        fate_f0_norm = fate_f0 / total

        meta["fate_prob_F4"] = fate_f4_norm
        meta["fate_prob_resolve"] = fate_f0_norm

        # Compare fate by subtype
        s1_fate = fate_f4_norm[subtypes == "S1"]
        s2_fate = fate_f4_norm[subtypes == "S2"]
        stat, pval = mannwhitneyu(s1_fate, s2_fate, alternative="two-sided")
        log.info(f"  S1 mean P(F4): {s1_fate.mean():.3f}")
        log.info(f"  S2 mean P(F4): {s2_fate.mean():.3f}")
        log.info(f"  S1 vs S2 P(F4): p={pval:.2e}")
    else:
        log.warning("  Not enough F4 or F0 samples for fate probability computation")
        meta["fate_prob_F4"] = np.nan
        meta["fate_prob_resolve"] = np.nan

    # -------------------------------------------------------------------
    # 6. Clinical correlates of subtypes
    # -------------------------------------------------------------------
    log.info("\n=== Clinical Correlates ===")

    clinical_results = []

    # Sex distribution
    if "sex" in meta.columns:
        sex_tab = pd.crosstab(meta["subtype"], meta["sex"])
        if sex_tab.shape == (2, 2):
            _, sex_p = fisher_exact(sex_tab.values)
        else:
            _, sex_p, _, _ = chi2_contingency(sex_tab)
        log.info(f"  Sex distribution:")
        log.info(f"    S1: {sex_tab.loc['S1'].to_dict() if 'S1' in sex_tab.index else 'N/A'}")
        log.info(f"    S2: {sex_tab.loc['S2'].to_dict() if 'S2' in sex_tab.index else 'N/A'}")
        log.info(f"    p={sex_p:.2e}")
        clinical_results.append({"variable": "sex", "test": "fisher/chi2", "pval": sex_p})

    # Continuous variables
    for var in ["age", "fibrosis_stage", "nas_score", "steatosis_grade",
                "lobular_inflammation_grade", "ballooning_grade"]:
        if var in meta.columns:
            s1_vals = pd.to_numeric(meta.loc[subtypes == "S1", var], errors="coerce").dropna()
            s2_vals = pd.to_numeric(meta.loc[subtypes == "S2", var], errors="coerce").dropna()
            if len(s1_vals) > 5 and len(s2_vals) > 5:
                stat, pval = mannwhitneyu(s1_vals, s2_vals, alternative="two-sided")
                log.info(f"  {var}: S1 median={s1_vals.median():.1f}, S2 median={s2_vals.median():.1f}, p={pval:.2e}")
                clinical_results.append({
                    "variable": var, "test": "Mann-Whitney",
                    "S1_median": s1_vals.median(), "S2_median": s2_vals.median(),
                    "pval": pval
                })

    # Disease distribution
    if "group_binary" in meta.columns:
        dis_tab = pd.crosstab(meta["subtype"], meta["group_binary"])
        log.info(f"  Disease distribution:\n{dis_tab}")
        _, dis_p, _, _ = chi2_contingency(dis_tab)
        clinical_results.append({"variable": "disease", "test": "chi2", "pval": dis_p})

    # Dataset distribution (to check if subtypes are batch-driven)
    ds_tab = pd.crosstab(meta["subtype"], meta["dataset"])
    _, ds_p, _, _ = chi2_contingency(ds_tab)
    log.info(f"  Dataset distribution p={ds_p:.2e} (p>>0.05 = NOT batch-driven)")
    clinical_results.append({"variable": "dataset", "test": "chi2", "pval": ds_p})

    clinical_df = pd.DataFrame(clinical_results)

    # -------------------------------------------------------------------
    # 7. Save results
    # -------------------------------------------------------------------
    log.info("\n=== Saving Results ===")

    # Bifurcation analysis
    bif_df = pd.DataFrame({
        "pseudotime_center": pt_centers,
        "divergence_score": multi_dim_div,
    })
    bif_df["in_bifurcation_window"] = (bif_df["pseudotime_center"] >= bif_start) & \
                                       (bif_df["pseudotime_center"] <= bif_end)
    bif_df.to_csv(os.path.join(OUTDIR, "bifurcation_analysis.csv"), index=False)
    log.info(f"  Saved bifurcation_analysis.csv")

    # Divergence genes (full-cohort, DESCRIPTIVE ONLY — see warning at the
    # divergence-detection block; not a leakage-free predictor selection source)
    if has_expr and len(div_df) > 0:
        div_df.to_csv(os.path.join(OUTDIR, "divergence_genes.csv"), index=False)
        log.info(f"  Saved divergence_genes.csv ({len(div_df)} genes) [descriptive]")

    # Fold-internal divergence ranking (LEAKAGE-FREE) — the canonical source for
    # building the `div_*` prognosis feature panel. Long format: one block of
    # rows per held-out LOCO fold, each block selected on training cohorts only.
    if has_expr and len(div_fold_df) > 0:
        div_fold_df.to_csv(
            os.path.join(OUTDIR, "divergence_genes_per_fold.csv"), index=False
        )
        log.info(
            f"  Saved divergence_genes_per_fold.csv "
            f"({div_fold_df['held_out_fold'].nunique()} folds, "
            f"{len(div_fold_df)} rows) [leakage-free]"
        )

    # Fate probabilities + per-sample data
    fate_df = meta[["sample_id", "subtype", "pseudotime", "fibrosis_stage",
                     "nas_score", "dataset", "group_binary",
                     "fate_prob_F4", "fate_prob_resolve"]].copy()
    fate_df.to_csv(os.path.join(OUTDIR, "fate_probabilities.csv"), index=False)
    log.info(f"  Saved fate_probabilities.csv")

    # Clinical correlates
    clinical_df.to_csv(os.path.join(OUTDIR, "bifurcation_clinical.csv"), index=False)
    log.info(f"  Saved bifurcation_clinical.csv")

    # Summary
    summary = {
        "n_samples": len(common),
        "n_S1": n_s1,
        "n_S2": n_s2,
        "bifurcation_pseudotime": bifurcation_pt,
        "bifurcation_window_start": bif_start,
        "bifurcation_window_end": bif_end,
        "peak_divergence_score": peak_div,
        "n_divergence_genes_005": (div_df["padj"] < 0.05).sum() if has_expr and len(div_df) > 0 else 0,
        "n_strong_divergence": ((div_df["padj"] < 0.05) & (div_df["cohens_d"].abs() > 0.5)).sum() if has_expr and len(div_df) > 0 else 0,
        "n_loco_folds_selected": int(div_fold_df["held_out_fold"].nunique()) if has_expr and len(div_fold_df) > 0 else 0,
    }
    pd.DataFrame([summary]).to_csv(os.path.join(OUTDIR, "bifurcation_summary.csv"), index=False)
    log.info(f"  Saved bifurcation_summary.csv")

    elapsed = (time.time() - t0) / 60
    log.info(f"\n=== 117: COMPLETE ({elapsed:.1f} min) ===")


if __name__ == "__main__":
    main()
