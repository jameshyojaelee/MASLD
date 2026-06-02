#!/usr/bin/env python3
"""
301: Palantir + DPT + Disease-Stage DE along Pseudotime.

For each of 5 core cell types, computes:
  1. Palantir pseudotime, branch probabilities, differentiation potential
  2. Scanpy DPT for comparison
  3. Disease-stage DE: sliding-window Wilcoxon tests for "switch genes"
  4. Gene dynamics (Spearman correlations along pseudotime)

Inputs:
    - Cell-type subsets from 300: pseudotime/{celltype}_subset.h5ad

Outputs (to results_gpu_v2/pseudotime/):
    - palantir_pseudotime_{celltype}.csv
    - dpt_pseudotime_{celltype}.csv
    - gene_dynamics_{celltype}.csv
    - pseudotime_corr_{celltype}.csv
    - stage_switch_genes_{celltype}.csv

Usage:
    sbatch run_pseudotime_pipeline.sh
"""

import os
import sys
import warnings
import logging
import gc
from pathlib import Path

# ---------------------------------------------------------------------------
# Seed pinning (pre-import) — must be set before downstream libs import RNG.
# Pins Python hash, random, numpy; scanpy seed pinned post-import below.
# Palantir internal RNG (diffusion-map eigensolver, waypoints) is less
# controllable but we pin what we can for reproducibility.
# ---------------------------------------------------------------------------
import random
os.environ["PYTHONHASHSEED"] = "42"
random.seed(42)
import numpy as np
np.random.seed(42)

import pandas as pd
from scipy import stats, sparse

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
# Imports
# ---------------------------------------------------------------------------
import scanpy as sc
import anndata as ad

# Seed-pin scanpy (best-effort; attribute exists in modern scanpy)
try:
    sc.settings.seed = 42
except Exception:
    pass

# Try importing palantir
try:
    import palantir
    HAS_PALANTIR = True
    log.info("Palantir version: %s", palantir.__version__)
except ImportError:
    HAS_PALANTIR = False
    log.warning("Palantir not installed — will use DPT only. "
                "Install via: pip install palantir")


# ---------------------------------------------------------------------------
# Disease condition ordering (for pseudotime direction validation)
# ---------------------------------------------------------------------------
CONDITION_ORDER = {"Healthy": 0, "NAFLD": 1, "MASLD": 1.5, "NASH": 2, "Cirrhotic": 3}

# Cell-type-specific Kupffer root anchors (canonical resident-KC markers).
# Literature: Remmerie 2020 (PMID 32888418), Tran 2020 (PMID 32937148),
# Jaitin 2019 (PMID 31257031) — in MASLD resident Kupffer cells DEPLETE
# and monocyte-derived LAMs EXPAND, so Kupffer -> LAM goes pseudotime 0 -> 1.
# For the Macrophage subset we must anchor root at max(TIMD4+MARCO+VCAM1).
KUPFFER_ROOT_MARKERS = ("TIMD4", "MARCO", "VCAM1")


def get_condition_score(obs_df):
    """Map condition labels to numeric scores."""
    return obs_df["condition"].map(CONDITION_ORDER).fillna(1.5)


def _gene_column(adata, gene):
    """Return 1-D dense vector of expression for `gene` (raises if missing)."""
    if gene not in adata.var_names:
        raise KeyError(f"Required marker gene {gene!r} missing from adata.var_names")
    idx = adata.var_names.get_loc(gene)
    col = adata.X[:, idx]
    if sparse.issparse(col):
        col = np.asarray(col.todense()).flatten()
    else:
        col = np.asarray(col).flatten()
    return col


def pick_kupffer_root_cell(adata, markers=KUPFFER_ROOT_MARKERS):
    """Return obs_name of the macrophage with highest summed Kupffer-marker
    expression. Raises if any marker is missing from var_names (no silent
    fallback — we want the pipeline to fail loudly if markers are absent).
    """
    missing = [g for g in markers if g not in adata.var_names]
    if missing:
        raise KeyError(
            f"Kupffer root markers missing from macrophage subset: {missing}. "
            f"Cannot anchor pseudotime root — refusing to fall back silently."
        )
    cols = np.stack([_gene_column(adata, g) for g in markers], axis=1)
    # Sum over markers (all log-normalized on same scale via scanpy.normalize)
    score = cols.sum(axis=1)
    root_idx = int(np.argmax(score))
    log.info(
        "Kupffer root anchor: argmax(%s) = cell %s (score=%.3f)",
        "+".join(markers), adata.obs_names[root_idx], float(score[root_idx]),
    )
    return adata.obs_names[root_idx]


# ---------------------------------------------------------------------------
# Palantir trajectory
# ---------------------------------------------------------------------------
def run_palantir(adata, ct_name):
    """Run Palantir pseudotime on a cell-type subset."""
    log.info("Running Palantir for %s (%d cells)...", ct_name, len(adata))

    # Palantir uses PCA embeddings stored in obsm
    pca_key = "X_pca"
    if pca_key not in adata.obsm:
        log.warning("No PCA found — computing...")
        sc.tl.pca(adata, n_comps=50)

    # Diffusion maps
    log.info("Computing diffusion maps...")
    dm_res = palantir.utils.run_diffusion_maps(
        pd.DataFrame(adata.obsm[pca_key], index=adata.obs_names)
    )

    # Multiscale space
    log.info("Determining multiscale space...")
    ms_data = palantir.utils.determine_multiscale_space(dm_res)

    # Root cell selection.
    # For Macrophages: canonical biology requires resident Kupffer cells as
    # the ORIGIN of the Kupffer -> LAM trajectory (literature: Remmerie 2020,
    # Tran 2020, Jaitin 2019). Anchor root at max(TIMD4+MARCO+VCAM1) so
    # pseudotime 0 = Kupffer, 1 = LAM.
    # For all other cell types: keep prior behavior (closest to Healthy centroid
    # in PCA space), since disease-stage progression is the appropriate axis.
    if ct_name == "Macrophages":
        root_cell = pick_kupffer_root_cell(adata)
        root_idx = int(adata.obs_names.get_loc(root_cell))
        log.info(
            "Root cell (Macrophage / Kupffer-anchored): %s (condition=%s)",
            root_cell, adata.obs.loc[root_cell, "condition"],
        )
    else:
        healthy_mask = adata.obs["condition"] == "Healthy"
        if healthy_mask.sum() > 0:
            healthy_pca = adata.obsm[pca_key][healthy_mask.values]
            centroid = healthy_pca.mean(axis=0)
            dists = np.linalg.norm(adata.obsm[pca_key] - centroid, axis=1)
            root_idx = np.argmin(dists)
            root_cell = adata.obs_names[root_idx]
        else:
            # Fallback: use cell with lowest condition score
            cond_scores = get_condition_score(adata.obs)
            root_idx = cond_scores.idxmin()
            root_cell = root_idx
        log.info("Root cell: %s (condition=%s)", root_cell,
                 adata.obs.loc[root_cell, "condition"])

    # Find terminal states: centroids of most diseased clusters
    terminal_states = []
    for cond in ["Cirrhotic", "NASH"]:
        cond_mask = adata.obs["condition"] == cond
        if cond_mask.sum() > 10:
            cond_pca = adata.obsm[pca_key][cond_mask.values]
            centroid = cond_pca.mean(axis=0)
            dists = np.linalg.norm(adata.obsm[pca_key] - centroid, axis=1)
            term_idx = np.argmin(dists)
            terminal_states.append(adata.obs_names[term_idx])

    if not terminal_states:
        # Use cell with highest condition score
        cond_scores = get_condition_score(adata.obs)
        terminal_states = [cond_scores.idxmax()]

    log.info("Terminal states: %s", terminal_states)

    # Run Palantir
    log.info("Running Palantir...")
    pr_res = palantir.core.run_palantir(
        ms_data,
        root_cell,
        terminal_states=terminal_states,
        num_waypoints=min(1200, len(adata) // 5),
    )

    log.info("Palantir pseudotime: median=%.3f, range=[%.3f, %.3f]",
             pr_res.pseudotime.median(),
             pr_res.pseudotime.min(),
             pr_res.pseudotime.max())

    # Validate: pseudotime should correlate with disease severity
    cond_scores = get_condition_score(adata.obs)
    valid = ~cond_scores.isna() & ~pr_res.pseudotime.isna()
    if valid.sum() > 30:
        rho, pval = stats.spearmanr(
            pr_res.pseudotime[valid], cond_scores[valid]
        )
        log.info("Pseudotime-condition correlation: rho=%.3f, p=%.2e", rho, pval)
        if rho < -0.1:
            log.warning("Negative correlation — pseudotime may be reversed!")

    # Macrophage-specific sanity check: Kupffer markers MUST anti-correlate
    # with pseudotime (resident KC = pseudotime 0). If they come out positive
    # the trajectory is inverted and all downstream labels are wrong.
    if ct_name == "Macrophages":
        pt_vals = pr_res.pseudotime.values
        score = np.zeros(len(pt_vals), dtype=float)
        for g in KUPFFER_ROOT_MARKERS:
            score = score + _gene_column(adata, g)
        mask = np.isfinite(pt_vals) & np.isfinite(score)
        if mask.sum() > 30:
            rho_k, pval_k = stats.spearmanr(pt_vals[mask], score[mask])
            log.info(
                "Macrophage sanity: Kupffer-score vs pseudotime rho=%.3f, p=%.2e",
                rho_k, pval_k,
            )
            if rho_k > 0:
                raise RuntimeError(
                    "Macrophage pseudotime direction check FAILED: "
                    f"Kupffer markers (TIMD4+MARCO+VCAM1) correlate POSITIVELY "
                    f"with pseudotime (rho={rho_k:.3f}) — trajectory is inverted. "
                    "Root anchor should force Kupffer -> LAM direction."
                )

    return pr_res, dm_res, ms_data


# ---------------------------------------------------------------------------
# Scanpy DPT
# ---------------------------------------------------------------------------
def run_dpt(adata, ct_name):
    """Run scanpy diffusion pseudotime."""
    log.info("Running scanpy DPT for %s...", ct_name)

    # Compute diffusion map
    sc.tl.diffmap(adata, n_comps=15)

    # Set root: cell closest to Healthy centroid in diffmap space
    healthy_mask = adata.obs["condition"] == "Healthy"
    if healthy_mask.sum() > 0:
        dm_embed = adata.obsm["X_diffmap"]
        healthy_dm = dm_embed[healthy_mask.values]
        centroid = healthy_dm.mean(axis=0)
        dists = np.linalg.norm(dm_embed - centroid, axis=1)
        root_idx = np.argmin(dists)
    else:
        root_idx = 0

    adata.uns["iroot"] = root_idx
    sc.tl.dpt(adata)

    dpt_vals = adata.obs["dpt_pseudotime"].values
    valid = np.isfinite(dpt_vals)
    log.info("DPT: median=%.3f, %d valid cells (%.1f%%)",
             np.nanmedian(dpt_vals), valid.sum(), 100 * valid.sum() / len(dpt_vals))

    return dpt_vals


# ---------------------------------------------------------------------------
# Disease-stage DE along pseudotime (switch genes)
# ---------------------------------------------------------------------------
def find_switch_genes(adata, pseudotime, n_bins=20, flank=2, min_cells_per_bin=50):
    """Find genes that switch on/off at specific pseudotime windows.

    For each bin, tests genes for specific activation vs. flanking bins
    using Wilcoxon rank-sum test.
    """
    log.info("Finding switch genes (%d bins, flank=%d)...", n_bins, flank)

    valid = np.isfinite(pseudotime)
    if valid.sum() < n_bins * min_cells_per_bin:
        log.warning("Too few valid cells (%d) for %d bins", valid.sum(), n_bins)
        n_bins = max(5, valid.sum() // min_cells_per_bin)

    # Equal-frequency binning
    pt_valid = pseudotime[valid]
    bin_edges = np.quantile(pt_valid, np.linspace(0, 1, n_bins + 1))
    bin_labels = np.digitize(pt_valid, bin_edges[1:-1])  # 0 to n_bins-1

    # Map bins to dominant condition
    obs_valid = adata.obs.iloc[np.where(valid)[0]].copy()
    obs_valid["pt_bin"] = bin_labels
    cond_scores = get_condition_score(obs_valid)
    bin_condition = obs_valid.groupby("pt_bin").apply(
        lambda g: g["condition"].mode().iloc[0] if len(g) > 0 else "Unknown"
    )

    # Get expression matrix — keep sparse, only densify per-gene slices
    X = adata.X[valid]
    genes = adata.var_names

    # Use top 5000 variable genes for speed
    hvg_mask = adata.var.get("highly_variable", pd.Series(True, index=genes))
    if hvg_mask.sum() > 5000:
        # Compute variance on sparse matrix efficiently
        if sparse.issparse(X):
            mean = np.array(X.mean(axis=0)).flatten()
            sq_mean = np.array(X.multiply(X).mean(axis=0)).flatten()
            var_per_gene = sq_mean - mean**2
        else:
            var_per_gene = np.var(X, axis=0)
        top_idx = np.argsort(-var_per_gene)[:5000]
        gene_indices = top_idx
    elif hvg_mask.sum() > 0:
        gene_indices = np.where(hvg_mask.values)[0]
    else:
        gene_indices = np.arange(min(5000, X.shape[1]))

    results = []
    for b in range(n_bins):
        in_bin = bin_labels == b
        n_in = in_bin.sum()
        if n_in < min_cells_per_bin:
            continue

        # Flanking bins
        flank_bins = [bb for bb in range(max(0, b - flank), min(n_bins, b + flank + 1))
                      if bb != b]
        in_flank = np.isin(bin_labels, flank_bins)
        n_flank = in_flank.sum()
        if n_flank < min_cells_per_bin:
            continue

        for gi in gene_indices:
            # Extract per-gene slices (works with both sparse and dense)
            col = X[:, gi]
            if sparse.issparse(col):
                col = np.array(col.todense()).flatten()
            elif col.ndim > 1:
                col = col.flatten()
            expr_bin = col[in_bin]
            expr_flank = col[in_flank]

            # Skip genes with no variance
            if np.std(expr_bin) == 0 and np.std(expr_flank) == 0:
                continue

            # Wilcoxon rank-sum (upregulated in bin)
            try:
                stat_up, pval_up = stats.mannwhitneyu(
                    expr_bin, expr_flank, alternative="greater"
                )
                stat_dn, pval_dn = stats.mannwhitneyu(
                    expr_bin, expr_flank, alternative="less"
                )
            except ValueError:
                continue

            # Use whichever direction is more significant
            if pval_up < pval_dn:
                pval, direction = pval_up, "UP"
            else:
                pval, direction = pval_dn, "DOWN"

            if pval < 0.01:  # Pre-filter for speed
                lfc = np.mean(expr_bin) - np.mean(expr_flank)
                results.append({
                    "gene": genes[gi],
                    "bin": b,
                    "bin_condition": bin_condition.get(b, "Unknown"),
                    "direction": direction,
                    "lfc": lfc,
                    "pval": pval,
                    "n_bin": n_in,
                    "n_flank": n_flank,
                    "mean_bin": np.mean(expr_bin),
                    "mean_flank": np.mean(expr_flank),
                })

    if not results:
        log.warning("No switch genes found!")
        return pd.DataFrame()

    df = pd.DataFrame(results)

    # BH FDR correction
    from statsmodels.stats.multitest import multipletests
    _, df["padj"], _, _ = multipletests(df["pval"], method="fdr_bh")

    # Keep FDR < 0.05
    sig = df[df["padj"] < 0.05].sort_values("padj")
    log.info("Switch genes: %d significant (FDR<0.05) out of %d tested",
             len(sig), len(df))

    # Identify transition-specific genes (unique to specific bin transitions)
    if len(sig) > 0:
        n_per_bin = sig.groupby("bin").size()
        log.info("Genes per bin:\n%s", n_per_bin.to_string())

    return sig


# ---------------------------------------------------------------------------
# Gene correlations with pseudotime
# ---------------------------------------------------------------------------
def compute_gene_correlations(adata, pseudotime, n_top=1000):
    """Vectorized Spearman correlation of each gene with pseudotime."""
    log.info("Computing gene-pseudotime correlations (vectorized)...")

    valid = np.isfinite(pseudotime)
    pt = pseudotime[valid]
    n_valid = valid.sum()

    X = adata.X[valid]
    genes = adata.var_names

    # Vectorized Spearman: rank pseudotime once, then correlate with ranked expression
    pt_ranks = stats.rankdata(pt)

    # Process in chunks to avoid OOM with dense conversion
    chunk_size = 1000
    n_genes = X.shape[1]
    corrs = []

    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        X_chunk = X[:, start:end]
        if sparse.issparse(X_chunk):
            X_chunk = X_chunk.toarray()

        for j in range(X_chunk.shape[1]):
            gi = start + j
            expr = X_chunk[:, j]
            if np.std(expr) == 0:
                continue
            expr_ranks = stats.rankdata(expr)
            # Pearson on ranks = Spearman
            rho = np.corrcoef(pt_ranks, expr_ranks)[0, 1]
            # Approximate p-value for Spearman
            t_stat = rho * np.sqrt((n_valid - 2) / (1 - rho**2 + 1e-12))
            pval = 2 * stats.t.sf(np.abs(t_stat), df=n_valid - 2)
            corrs.append({"gene": genes[gi], "spearman_rho": rho, "pval": pval})

    df = pd.DataFrame(corrs)

    if len(df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, df["padj"], _, _ = multipletests(df["pval"], method="fdr_bh")
        df = df.sort_values("spearman_rho", key=abs, ascending=False)

    n_sig = (df["padj"] < 0.05).sum() if "padj" in df.columns else 0
    log.info("Gene correlations: %d genes, %d significant (FDR<0.05)", len(df), n_sig)

    # Top positively and negatively correlated
    n_pos = (df["spearman_rho"] > 0).sum()
    n_neg = (df["spearman_rho"] < 0).sum()
    log.info("Positive: %d, Negative: %d", n_pos, n_neg)

    if len(df) > 0:
        top_pos = df[df["spearman_rho"] > 0].head(5)["gene"].tolist()
        top_neg = df[df["spearman_rho"] < 0].head(5)["gene"].tolist()
        log.info("Top positive: %s", top_pos)
        log.info("Top negative: %s", top_neg)

    return df


# ---------------------------------------------------------------------------
# Gene dynamics (binned expression along pseudotime)
# ---------------------------------------------------------------------------
def compute_gene_dynamics(adata, pseudotime, n_bins=100, n_top_genes=2000):
    """Compute binned mean expression along pseudotime for top variable genes."""
    log.info("Computing gene dynamics (%d bins, %d genes)...", n_bins, n_top_genes)

    valid = np.isfinite(pseudotime)
    pt = pseudotime[valid]
    X = adata.X[valid]

    # Select top variable genes (sparse-compatible variance)
    if sparse.issparse(X):
        mean = np.array(X.mean(axis=0)).flatten()
        sq_mean = np.array(X.multiply(X).mean(axis=0)).flatten()
        gene_var = sq_mean - mean**2
    else:
        gene_var = np.var(X, axis=0)
    top_idx = np.argsort(-gene_var)[:n_top_genes]
    top_genes = adata.var_names[top_idx]

    # Only densify the selected genes
    X_top = X[:, top_idx]
    if sparse.issparse(X_top):
        X_top = X_top.toarray()

    # Bin cells by pseudotime
    bin_edges = np.linspace(pt.min(), pt.max(), n_bins + 1)
    bin_labels = np.digitize(pt, bin_edges[1:-1])

    # Compute mean expression per bin
    dynamics = np.zeros((n_top_genes, n_bins))
    for b in range(n_bins):
        mask = bin_labels == b
        if mask.sum() > 0:
            dynamics[:, b] = X_top[mask].mean(axis=0)

    df = pd.DataFrame(
        dynamics,
        index=top_genes,
        columns=[f"bin_{b}" for b in range(n_bins)],
    )
    df.index.name = "gene"
    log.info("Gene dynamics computed: %s", df.shape)

    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("301: Palantir + DPT + Disease-Stage DE")
    log.info("=" * 60)

    for ct in CELL_TYPES:
        log.info("\n" + "=" * 50)
        log.info("Processing: %s", ct)
        log.info("=" * 50)

        h5ad_path = os.path.join(PT_DIR, f"{ct}_subset.h5ad")
        if not os.path.exists(h5ad_path):
            log.warning("Subset not found: %s — skipping", h5ad_path)
            continue

        adata = sc.read_h5ad(h5ad_path)
        log.info("Loaded %s: %s", ct, adata.shape)

        # --- Palantir ---
        palantir_pt = None
        if HAS_PALANTIR:
            try:
                pr_res, dm_res, ms_data = run_palantir(adata, ct)
                palantir_pt = pr_res.pseudotime.values

                # Build output DataFrame
                pt_df = pd.DataFrame(index=adata.obs_names)
                pt_df["palantir_pseudotime"] = pr_res.pseudotime
                pt_df["palantir_entropy"] = pr_res.entropy
                for col in pr_res.branch_probs.columns:
                    pt_df[f"branch_prob_{col}"] = pr_res.branch_probs[col]
                pt_df["condition"] = adata.obs["condition"].values
                pt_df["leiden_substate"] = adata.obs["leiden_substate"].values

                out_path = os.path.join(PT_DIR, f"palantir_pseudotime_{ct}.csv")
                pt_df.to_csv(out_path)
                log.info("Saved Palantir results to %s", out_path)

            except Exception as e:
                log.error("Palantir failed for %s: %s", ct, e)
                import traceback
                traceback.print_exc()

        # --- Scanpy DPT ---
        try:
            dpt_vals = run_dpt(adata, ct)

            dpt_df = pd.DataFrame(index=adata.obs_names)
            dpt_df["dpt_pseudotime"] = dpt_vals
            dpt_df["condition"] = adata.obs["condition"].values
            dpt_df["leiden_substate"] = adata.obs["leiden_substate"].values

            out_path = os.path.join(PT_DIR, f"dpt_pseudotime_{ct}.csv")
            dpt_df.to_csv(out_path)
            log.info("Saved DPT results to %s", out_path)

        except Exception as e:
            log.error("DPT failed for %s: %s", ct, e)
            dpt_vals = None

        # Choose best pseudotime for downstream analyses
        if palantir_pt is not None:
            pseudotime = palantir_pt
            pt_source = "palantir"
        elif dpt_vals is not None:
            pseudotime = dpt_vals
            pt_source = "dpt"
        else:
            log.error("No pseudotime available for %s — skipping downstream", ct)
            del adata
            gc.collect()
            continue

        log.info("Using %s pseudotime for downstream analyses", pt_source)

        # --- Disease-stage switch genes ---
        try:
            switch_df = find_switch_genes(adata, pseudotime)
            if len(switch_df) > 0:
                out_path = os.path.join(PT_DIR, f"stage_switch_genes_{ct}.csv")
                switch_df.to_csv(out_path, index=False)
                log.info("Saved %d switch genes to %s", len(switch_df), out_path)
        except Exception as e:
            log.error("Switch gene analysis failed for %s: %s", ct, e)

        # --- Gene-pseudotime correlations ---
        try:
            corr_df = compute_gene_correlations(adata, pseudotime)
            out_path = os.path.join(PT_DIR, f"pseudotime_corr_{ct}.csv")
            corr_df.to_csv(out_path, index=False)
            log.info("Saved gene correlations to %s", out_path)
        except Exception as e:
            log.error("Gene correlation failed for %s: %s", ct, e)

        # --- Gene dynamics ---
        try:
            dynamics_df = compute_gene_dynamics(adata, pseudotime)
            out_path = os.path.join(PT_DIR, f"gene_dynamics_{ct}.csv")
            dynamics_df.to_csv(out_path)
            log.info("Saved gene dynamics to %s", out_path)
        except Exception as e:
            log.error("Gene dynamics failed for %s: %s", ct, e)

        # --- Cross-method comparison (if both available) ---
        if palantir_pt is not None and dpt_vals is not None:
            valid = np.isfinite(palantir_pt) & np.isfinite(dpt_vals)
            if valid.sum() > 30:
                rho, pval = stats.spearmanr(palantir_pt[valid], dpt_vals[valid])
                log.info("Palantir vs DPT correlation: rho=%.3f, p=%.2e", rho, pval)

        del adata
        gc.collect()

    log.info("\n=== 301: Palantir + DPT COMPLETE ===")


if __name__ == "__main__":
    main()
