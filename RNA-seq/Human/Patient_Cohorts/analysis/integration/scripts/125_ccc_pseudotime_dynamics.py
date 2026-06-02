#!/usr/bin/env python3
"""
125_ccc_pseudotime_dynamics.py
CCC dynamics along disease pseudotime.

Uses L-R interaction scores from Script 124 to identify:
  1. L-R pairs with monotonic increase/decrease along pseudotime (Spearman)
  2. Pseudotime-binned mean L-R scores (20 bins) per communication axis
  3. Clusters of L-R pairs with shared temporal dynamics (k-means)
  4. Pathway annotation per cluster (TGFb, WNT, NOTCH, TNF, etc.)

Input:
  - results/progression/ccc_lr_scores.csv       (from Script 124)
  - results/progression/ccc_transition_de.csv    (from Script 124, for pathway)
  - results/progression/consensus_pseudotime.csv (from Script 116)
  - results/staging_classifier/modeling_metadata.csv

Output (all to results/progression/):
  - ccc_pseudotime_dynamics.csv  (per LR pair: Spearman, trend, cluster)
  - ccc_lr_clusters.csv          (cluster centers + pathway enrichment)
  - ccc_binned_scores.csv        (pseudotime-binned mean scores)

SLURM: cpu partition, 8 CPUs, 64G RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg125_ccc_pt \
         --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \
         --output=logs/125_ccc_pt_%j.out \
         --error=logs/125_ccc_pt_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 125_ccc_pseudotime_dynamics.py'"
"""

# Prevent GPU initialization on CPU nodes
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
from scipy.stats import spearmanr, pearsonr, mannwhitneyu, kruskal
from scipy.ndimage import uniform_filter1d
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEED = 42
np.random.seed(SEED)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PROG_DIR = os.path.join(INTEG, "results/progression")
STAG_DIR = os.path.join(INTEG, "results/staging_classifier")

N_BINS = 20
MIN_SAMPLES_PER_BIN = 10
MIN_NONZERO_FRAC = 0.20  # Skip LR pairs with >80% zeros
SPEARMAN_PVAL_THRESH = 0.05
MIN_ABS_RHO = 0.15  # Minimum |rho| for "monotonic" classification


def main():
    t_start = time.time()
    log.info("=== 125: CCC Pseudotime Dynamics ===")

    # ── 1. Load data ──────────────────────────────────────────────────────
    log.info("Loading L-R scores from Script 124...")
    scores_path = os.path.join(PROG_DIR, "ccc_lr_scores.csv")
    if not os.path.exists(scores_path):
        log.error("ccc_lr_scores.csv not found. Run Script 124 first.")
        sys.exit(1)
    scores = pd.read_csv(scores_path)
    log.info("  L-R scores: %d rows, %d unique LR pairs, %d axes",
             len(scores), scores["lr_pair"].nunique(), scores["axis"].nunique())

    log.info("Loading consensus pseudotime...")
    pt_path = os.path.join(PROG_DIR, "consensus_pseudotime.csv")
    if not os.path.exists(pt_path):
        log.error("consensus_pseudotime.csv not found. Run Script 116 first.")
        sys.exit(1)
    pseudotime = pd.read_csv(pt_path)
    log.info("  Pseudotime: %d samples", len(pseudotime))

    # Use consensus pseudotime if available, else fallback to NAS-VAE
    if "pseudotime_consensus" in pseudotime.columns:
        pt_col = "pseudotime_consensus"
    elif "pseudotime_nasvae" in pseudotime.columns:
        pt_col = "pseudotime_nasvae"
    else:
        log.error("No pseudotime column found in consensus_pseudotime.csv")
        sys.exit(1)
    log.info("  Using pseudotime column: %s", pt_col)

    log.info("Loading metadata...")
    meta = pd.read_csv(os.path.join(STAG_DIR, "modeling_metadata.csv"))
    log.info("  Metadata: %d samples", len(meta))

    # Load pathway annotations from transition DE (if available)
    de_path = os.path.join(PROG_DIR, "ccc_transition_de.csv")
    pathway_map = {}
    if os.path.exists(de_path):
        de_df = pd.read_csv(de_path)
        if "pathway" in de_df.columns:
            for _, row in de_df[["ligand", "receptor", "pathway"]].drop_duplicates().iterrows():
                lr_key = f"{row['ligand']}::{row['receptor']}"
                if pd.notna(row["pathway"]):
                    pathway_map[lr_key] = row["pathway"]
            log.info("  Loaded pathway annotations for %d L-R pairs", len(pathway_map))

    # ── 2. Merge pseudotime with scores ───────────────────────────────────
    log.info("Merging pseudotime with L-R scores...")
    pt_lookup = pseudotime.set_index("sample_id")[pt_col].to_dict()
    fib_lookup = dict(zip(meta["sample_id"], meta["fibrosis_stage"]))

    scores["pseudotime"] = scores["sample_id"].map(pt_lookup)
    scores["fibrosis_stage"] = scores["sample_id"].map(fib_lookup)

    # Drop samples without pseudotime
    n_before = len(scores)
    scores = scores.dropna(subset=["pseudotime"])
    log.info("  Samples with pseudotime: %d -> %d rows (dropped %d)",
             n_before, len(scores), n_before - len(scores))

    # ── 3. Create pseudotime bins ─────────────────────────────────────────
    log.info("Creating %d pseudotime bins...", N_BINS)
    pt_values = scores["pseudotime"].unique()
    pt_min, pt_max = pt_values.min(), pt_values.max()
    bin_edges = np.linspace(pt_min, pt_max, N_BINS + 1)
    scores["pt_bin"] = pd.cut(scores["pseudotime"], bins=bin_edges,
                               labels=range(N_BINS), include_lowest=True)
    scores["pt_bin"] = scores["pt_bin"].astype(float)

    # Report bin sizes
    sample_bins = scores.groupby("pt_bin")["sample_id"].nunique()
    log.info("  Samples per bin: min=%d, max=%d, median=%d",
             sample_bins.min(), sample_bins.max(), int(sample_bins.median()))

    # ── 4. Compute binned mean scores ─────────────────────────────────────
    log.info("Computing pseudotime-binned mean L-R scores...")
    binned = scores.groupby(["axis", "lr_pair", "pt_bin"]).agg(
        mean_score=("score", "mean"),
        median_score=("score", "median"),
        n_samples=("sample_id", "nunique"),
        nonzero_frac=("score", lambda x: (x > 0).mean())
    ).reset_index()

    log.info("  Binned scores: %d rows", len(binned))

    # Save binned scores
    binned_path = os.path.join(PROG_DIR, "ccc_binned_scores.csv")
    binned.to_csv(binned_path, index=False)
    log.info("  Saved ccc_binned_scores.csv")

    # ── 5. Spearman correlation with pseudotime ───────────────────────────
    log.info("Computing Spearman correlation of each L-R pair with pseudotime...")

    # Work per axis × LR pair
    axis_lr_groups = scores.groupby(["axis", "lr_pair"])
    dynamics_results = []

    for (axis, lr_pair), group in axis_lr_groups:
        # Filter pairs where most samples are zero
        nonzero_frac = (group["score"] > 0).mean()
        if nonzero_frac < MIN_NONZERO_FRAC:
            continue

        pt = group["pseudotime"].values
        sc = group["score"].values

        # Spearman correlation
        rho, pval = spearmanr(pt, sc)

        # Additional metrics
        mean_score = sc.mean()
        std_score = sc.std()

        # Trend classification
        if pval < SPEARMAN_PVAL_THRESH and abs(rho) >= MIN_ABS_RHO:
            trend = "increasing" if rho > 0 else "decreasing"
        else:
            trend = "stable"

        # Compute per-fibrosis-stage means (for trajectory shape)
        stage_means = {}
        for fib in sorted(group["fibrosis_stage"].dropna().unique()):
            stage_mask = group["fibrosis_stage"] == fib
            if stage_mask.sum() >= 5:
                stage_means[f"mean_F{int(fib)}"] = group.loc[stage_mask, "score"].mean()

        lig, rec = lr_pair.split("::")
        pathway = pathway_map.get(lr_pair, np.nan)

        result = {
            "axis": axis,
            "lr_pair": lr_pair,
            "ligand": lig,
            "receptor": rec,
            "spearman_rho": rho,
            "spearman_pval": pval,
            "trend": trend,
            "mean_score": mean_score,
            "std_score": std_score,
            "nonzero_frac": nonzero_frac,
            "n_samples": len(group),
            "pathway": pathway
        }
        result.update(stage_means)
        dynamics_results.append(result)

    dynamics_df = pd.DataFrame(dynamics_results)
    log.info("  Dynamics computed for %d axis-LR combinations", len(dynamics_df))

    # BH correction per axis
    from statsmodels.stats.multitest import multipletests
    dynamics_df["spearman_padj"] = np.nan
    for axis in dynamics_df["axis"].unique():
        mask = dynamics_df["axis"] == axis
        pvals = dynamics_df.loc[mask, "spearman_pval"].values
        if len(pvals) > 0:
            _, padj, _, _ = multipletests(pvals, method="fdr_bh")
            dynamics_df.loc[mask, "spearman_padj"] = padj

    # Re-classify trend using adjusted p-values
    dynamics_df["trend_adj"] = "stable"
    sig_mask = (dynamics_df["spearman_padj"] < SPEARMAN_PVAL_THRESH) & \
               (dynamics_df["spearman_rho"].abs() >= MIN_ABS_RHO)
    dynamics_df.loc[sig_mask & (dynamics_df["spearman_rho"] > 0), "trend_adj"] = "increasing"
    dynamics_df.loc[sig_mask & (dynamics_df["spearman_rho"] < 0), "trend_adj"] = "decreasing"

    n_inc = (dynamics_df["trend_adj"] == "increasing").sum()
    n_dec = (dynamics_df["trend_adj"] == "decreasing").sum()
    n_stable = (dynamics_df["trend_adj"] == "stable").sum()
    log.info("  Trend distribution (BH-adjusted): %d increasing, %d decreasing, %d stable",
             n_inc, n_dec, n_stable)

    # Per-axis summary
    log.info("  Per-axis trend summary:")
    for axis in sorted(dynamics_df["axis"].unique()):
        sub = dynamics_df[dynamics_df["axis"] == axis]
        log.info("    %s: %d inc, %d dec, %d stable (of %d)",
                 axis,
                 (sub["trend_adj"] == "increasing").sum(),
                 (sub["trend_adj"] == "decreasing").sum(),
                 (sub["trend_adj"] == "stable").sum(),
                 len(sub))

    # ── 6. Cluster L-R pairs by pseudotime trajectory shape ──────────────
    log.info("Clustering L-R pairs by pseudotime trajectory...")

    # Build trajectory matrix: for each axis-LR, get the N_BINS-length smoothed score vector
    # Only cluster non-stable pairs (those with significant monotonic trend)
    sig_pairs = dynamics_df[dynamics_df["trend_adj"] != "stable"].copy()
    log.info("  Clustering %d significant (non-stable) axis-LR pairs", len(sig_pairs))

    if len(sig_pairs) < 5:
        log.warning("Too few significant L-R pairs for clustering. Saving dynamics only.")
        dynamics_df.to_csv(os.path.join(PROG_DIR, "ccc_pseudotime_dynamics.csv"), index=False)
        log.info("  Saved ccc_pseudotime_dynamics.csv")

        # Empty cluster file
        empty_clusters = pd.DataFrame(columns=["cluster", "n_pairs", "mean_rho",
                                                "dominant_trend", "top_pathways"])
        empty_clusters.to_csv(os.path.join(PROG_DIR, "ccc_lr_clusters.csv"), index=False)
        log.info("  Saved empty ccc_lr_clusters.csv")

        elapsed = time.time() - t_start
        log.info("=== 125: COMPLETE (%.1f min) ===", elapsed / 60)
        return

    # Build trajectory vectors from binned scores
    trajectory_vectors = []
    pair_labels = []

    for _, row in sig_pairs.iterrows():
        axis = row["axis"]
        lr_pair = row["lr_pair"]
        sub_binned = binned[(binned["axis"] == axis) & (binned["lr_pair"] == lr_pair)]

        if len(sub_binned) < N_BINS * 0.5:
            # Skip pairs with too many missing bins
            continue

        # Fill missing bins with 0 and create trajectory vector
        traj = np.zeros(N_BINS)
        for _, brow in sub_binned.iterrows():
            b = int(brow["pt_bin"])
            if 0 <= b < N_BINS:
                traj[b] = brow["mean_score"]

        # Smooth with window=3 to reduce noise
        traj_smooth = uniform_filter1d(traj, size=3)

        # Z-normalize per pair (we want to cluster by SHAPE, not magnitude)
        if traj_smooth.std() > 1e-10:
            traj_norm = (traj_smooth - traj_smooth.mean()) / traj_smooth.std()
        else:
            traj_norm = traj_smooth - traj_smooth.mean()

        trajectory_vectors.append(traj_norm)
        pair_labels.append((axis, lr_pair))

    if len(trajectory_vectors) < 5:
        log.warning("Too few trajectory vectors for clustering (%d). Saving dynamics only.",
                    len(trajectory_vectors))
        dynamics_df["cluster"] = np.nan
        dynamics_df.to_csv(os.path.join(PROG_DIR, "ccc_pseudotime_dynamics.csv"), index=False)
        empty_clusters = pd.DataFrame(columns=["cluster", "n_pairs", "mean_rho",
                                                "dominant_trend", "top_pathways"])
        empty_clusters.to_csv(os.path.join(PROG_DIR, "ccc_lr_clusters.csv"), index=False)
        elapsed = time.time() - t_start
        log.info("=== 125: COMPLETE (%.1f min) ===", elapsed / 60)
        return

    traj_matrix = np.array(trajectory_vectors)
    log.info("  Trajectory matrix: %d pairs x %d bins", traj_matrix.shape[0], traj_matrix.shape[1])

    # Determine optimal k via silhouette score (k=2..8)
    k_range = range(2, min(9, len(trajectory_vectors)))
    best_k = 2
    best_sil = -1

    for k in k_range:
        km = KMeans(n_clusters=k, random_state=SEED, n_init=10, max_iter=300)
        labels = km.fit_predict(traj_matrix)
        # Need at least 2 clusters with >1 member
        if len(set(labels)) < 2:
            continue
        sil = silhouette_score(traj_matrix, labels)
        log.info("    k=%d: silhouette=%.3f", k, sil)
        if sil > best_sil:
            best_sil = sil
            best_k = k

    log.info("  Optimal k=%d (silhouette=%.3f)", best_k, best_sil)

    # Final clustering
    km_final = KMeans(n_clusters=best_k, random_state=SEED, n_init=20, max_iter=500)
    cluster_labels = km_final.fit_predict(traj_matrix)

    # Map cluster labels back to dynamics_df
    cluster_map = {}
    for i, (axis, lr_pair) in enumerate(pair_labels):
        cluster_map[(axis, lr_pair)] = cluster_labels[i]

    dynamics_df["cluster"] = dynamics_df.apply(
        lambda r: cluster_map.get((r["axis"], r["lr_pair"]), np.nan), axis=1
    )

    # ── 7. Annotate clusters ─────────────────────────────────────────────
    log.info("Annotating clusters...")

    cluster_summaries = []
    for c in range(best_k):
        c_mask = dynamics_df["cluster"] == c
        c_df = dynamics_df[c_mask]

        n_pairs = len(c_df)
        if n_pairs == 0:
            continue

        mean_rho = c_df["spearman_rho"].mean()

        # Dominant trend
        trend_counts = c_df["trend_adj"].value_counts()
        dominant_trend = trend_counts.index[0] if len(trend_counts) > 0 else "mixed"

        # Top pathways
        pw_list = c_df["pathway"].dropna().tolist()
        all_pws = []
        for pw in pw_list:
            all_pws.extend(str(pw).split(";"))
        if all_pws:
            from collections import Counter
            pw_counts = Counter(all_pws)
            top_pws = ";".join([p for p, _ in pw_counts.most_common(5)])
        else:
            top_pws = "NA"

        # Top axes
        axis_counts = c_df["axis"].value_counts()
        top_axes = ";".join(axis_counts.head(3).index.tolist())

        # Cluster center trajectory (mean of z-normalized trajectories)
        c_indices = [i for i, (ax, lr) in enumerate(pair_labels)
                     if cluster_labels[i] == c]
        if c_indices:
            center_traj = traj_matrix[c_indices].mean(axis=0)
            # Classify shape: early-peak, late-peak, monotonic-up, monotonic-down, U-shape
            peak_bin = np.argmax(center_traj)
            trough_bin = np.argmin(center_traj)
            if center_traj[-1] - center_traj[0] > 0.5:
                shape = "monotonic_up"
            elif center_traj[0] - center_traj[-1] > 0.5:
                shape = "monotonic_down"
            elif peak_bin < N_BINS * 0.4:
                shape = "early_peak"
            elif peak_bin > N_BINS * 0.6:
                shape = "late_peak"
            else:
                shape = "mid_peak"
        else:
            shape = "unknown"

        # Representative L-R pairs (top 5 by |rho|)
        top_pairs = c_df.nlargest(5, "spearman_rho", keep="first")["lr_pair"].tolist()

        cluster_summaries.append({
            "cluster": int(c),
            "n_pairs": n_pairs,
            "mean_rho": round(mean_rho, 3),
            "dominant_trend": dominant_trend,
            "trajectory_shape": shape,
            "top_pathways": top_pws,
            "top_axes": top_axes,
            "representative_pairs": ";".join(top_pairs[:5])
        })

        log.info("  Cluster %d: %d pairs, mean_rho=%.3f, trend=%s, shape=%s",
                 c, n_pairs, mean_rho, dominant_trend, shape)
        log.info("    Top pathways: %s", top_pws[:80])
        log.info("    Top pairs: %s", ", ".join(top_pairs[:3]))

    clusters_df = pd.DataFrame(cluster_summaries)

    # ── 8. Identify key pathway dynamics ──────────────────────────────────
    log.info("Identifying pathway-level dynamics...")

    # Map L-R pairs to canonical signaling pathways using keyword matching
    CANONICAL_PATHWAYS = {
        "TGFbeta": ["TGFB", "BMP", "ACVR", "SMAD", "TGFBR", "GDF"],
        "WNT": ["WNT", "FZD", "LRP5", "LRP6", "DKK", "SFRP", "RSPO"],
        "NOTCH": ["NOTCH", "DLL", "JAG", "HES", "HEY"],
        "TNF": ["TNF", "TNFRSF", "TNFSF", "TRAIL", "FASL"],
        "PDGF_FGF": ["PDGF", "FGF", "FGFR", "PDGFR"],
        "VEGF": ["VEGF", "FLT", "KDR", "NRP"],
        "Chemokine": ["CCL", "CXCL", "CCR", "CXCR"],
        "Interleukin": ["IL1", "IL2", "IL4", "IL6", "IL10", "IL13", "IL17", "IL18", "IL33"],
        "ECM_Integrin": ["COL", "FN1", "LAMA", "LAMB", "ITGA", "ITGB", "SDC", "THBS"],
        "Complement": ["C3", "C1Q", "CFB", "CFH", "CR1", "CD46"],
        "EGF": ["EGF", "EGFR", "ERBB", "NRG", "AREG", "HBEGF"],
        "HGF_MET": ["HGF", "MET"],
        "Ephrin": ["EPH", "EFN"],
        "Semaphorin": ["SEMA", "PLXN", "NRP"],
        "ADAM_MMP": ["ADAM", "MMP"],
    }

    # Annotate each pair with canonical pathway
    def classify_canonical_pathway(ligand, receptor):
        genes = ligand.split("_") + receptor.split("_")
        matched = []
        for pw_name, keywords in CANONICAL_PATHWAYS.items():
            for gene in genes:
                if any(gene.upper().startswith(kw.upper()) for kw in keywords):
                    matched.append(pw_name)
                    break
        return ";".join(matched) if matched else "Other"

    dynamics_df["canonical_pathway"] = dynamics_df.apply(
        lambda r: classify_canonical_pathway(r["ligand"], r["receptor"]), axis=1
    )

    # Pathway-level statistics
    log.info("  Canonical pathway summary:")
    pw_summary = []
    for pw_name in CANONICAL_PATHWAYS:
        pw_mask = dynamics_df["canonical_pathway"].str.contains(pw_name, na=False)
        if pw_mask.sum() == 0:
            continue
        pw_df = dynamics_df[pw_mask]
        n_total = len(pw_df)
        n_inc = (pw_df["trend_adj"] == "increasing").sum()
        n_dec = (pw_df["trend_adj"] == "decreasing").sum()
        mean_rho = pw_df["spearman_rho"].mean()

        pw_summary.append({
            "canonical_pathway": pw_name,
            "n_pairs": n_total,
            "n_increasing": n_inc,
            "n_decreasing": n_dec,
            "n_stable": n_total - n_inc - n_dec,
            "mean_spearman_rho": round(mean_rho, 3),
            "dominant_direction": "increasing" if n_inc > n_dec else ("decreasing" if n_dec > n_inc else "mixed")
        })
        if n_inc + n_dec > 0:
            log.info("    %s: %d pairs (%d inc, %d dec, mean_rho=%.3f)",
                     pw_name, n_total, n_inc, n_dec, mean_rho)

    pw_summary_df = pd.DataFrame(pw_summary)

    # ── 9. Save results ──────────────────────────────────────────────────
    log.info("Saving results...")

    # Sort dynamics by significance
    dynamics_df = dynamics_df.sort_values(["spearman_padj", "spearman_rho"],
                                           ascending=[True, False])
    dynamics_df.to_csv(os.path.join(PROG_DIR, "ccc_pseudotime_dynamics.csv"), index=False)
    log.info("  Saved ccc_pseudotime_dynamics.csv (%d rows)", len(dynamics_df))

    clusters_df.to_csv(os.path.join(PROG_DIR, "ccc_lr_clusters.csv"), index=False)
    log.info("  Saved ccc_lr_clusters.csv (%d clusters)", len(clusters_df))

    if len(pw_summary_df) > 0:
        pw_summary_df.to_csv(os.path.join(PROG_DIR, "ccc_pathway_dynamics.csv"), index=False)
        log.info("  Saved ccc_pathway_dynamics.csv (%d pathways)", len(pw_summary_df))

    # ── 10. Biological highlights ─────────────────────────────────────────
    log.info("=== Biological Highlights ===")

    # Top increasing L-R pairs (disease-gained communication)
    inc = dynamics_df[dynamics_df["trend_adj"] == "increasing"].nsmallest(10, "spearman_padj")
    if len(inc) > 0:
        log.info("  Top disease-gained L-R interactions:")
        for _, r in inc.iterrows():
            log.info("    %s [%s] rho=%.3f padj=%.2e (%s)",
                     r["lr_pair"], r["axis"], r["spearman_rho"],
                     r["spearman_padj"], r["canonical_pathway"])

    # Top decreasing L-R pairs (disease-lost communication)
    dec = dynamics_df[dynamics_df["trend_adj"] == "decreasing"].nsmallest(10, "spearman_padj")
    if len(dec) > 0:
        log.info("  Top disease-lost L-R interactions:")
        for _, r in dec.iterrows():
            log.info("    %s [%s] rho=%.3f padj=%.2e (%s)",
                     r["lr_pair"], r["axis"], r["spearman_rho"],
                     r["spearman_padj"], r["canonical_pathway"])

    # Pathway-level narrative
    if len(pw_summary_df) > 0:
        log.info("  Pathway-level dynamics:")
        for _, r in pw_summary_df.iterrows():
            if r["n_increasing"] + r["n_decreasing"] > 0:
                log.info("    %s: %s (%d inc / %d dec)",
                         r["canonical_pathway"], r["dominant_direction"],
                         r["n_increasing"], r["n_decreasing"])

    elapsed = time.time() - t_start
    log.info("=== 125: COMPLETE (%.1f min) ===", elapsed / 60)


if __name__ == "__main__":
    main()
