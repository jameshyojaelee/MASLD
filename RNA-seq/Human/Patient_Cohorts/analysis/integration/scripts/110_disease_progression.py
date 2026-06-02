#!/usr/bin/env python3
"""
110_disease_progression.py
Disease trajectory from cross-sectional data using optimal transport
and pseudotime (simplified TIGON-inspired approach).

Pipeline:
  1. Load NAS-VAE embeddings (64-dim, 1,444 samples) + metadata
  2. Define "timepoints" as NAS score groups or fibrosis stages
  3. Compute optimal transport plans between adjacent timepoints
  4. Extract "transport genes" (genes changing most along transport path)
  5. Compute velocity field (progression direction per sample)
  6. Identify bifurcation points (connects to NMF k=2 subtypes)
  7. Compute UMAP with velocity arrows for visualization

Output (all to results/staging_classifier/):
  - progression_transport_plans.csv    (OT assignments between timepoints)
  - progression_velocity_field.csv     (per-sample velocity vectors)
  - progression_bifurcation.csv        (bifurcation point analysis)
  - progression_umap_data.csv          (UMAP coords + velocities for R plotting)

SLURM: cpu partition, 8 CPUs, 32G RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg110_progression \
         --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
         --output=logs/110_progression_%j.out \
         --error=logs/110_progression_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 110_disease_progression.py'"
"""

# Must be set before ANY other imports to prevent CUDA crash on cpu nodes
import os
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
from scipy.optimize import linear_sum_assignment
from scipy.stats import mannwhitneyu, spearmanr
import h5py

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
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
OUTDIR = os.path.join(RDIR, "staging_classifier")
LOGDIR = os.path.join(OUTDIR, "logs")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

# Input files
EMBED_PATH = os.path.join(OUTDIR, "nas_embeddings_all_samples.csv")
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")
H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
UNIFIED_META = os.path.join(INT, "metadata/unified_metadata.csv")

# NMF subtype assignments (from Phase F)
NMF_PATH = os.path.join(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")

# Outputs
OUT_TRANSPORT = os.path.join(OUTDIR, "progression_transport_plans.csv")
OUT_VELOCITY = os.path.join(OUTDIR, "progression_velocity_field.csv")
OUT_BIFURCATION = os.path.join(OUTDIR, "progression_bifurcation.csv")
OUT_UMAP = os.path.join(OUTDIR, "progression_umap_data.csv")

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(LOGDIR, "110_disease_progression.log")),
    ],
)
log = logging.getLogger(__name__)

print("=" * 70)
print("110: Disease Progression via Optimal Transport + Pseudotime")
print("=" * 70)


# ===================================================================
# 1. Load data
# ===================================================================

def load_embeddings():
    """Load NAS-VAE 64-dim embeddings for all samples."""
    if not os.path.exists(EMBED_PATH):
        log.error("NAS embeddings not found: %s", EMBED_PATH)
        sys.exit(1)

    df = pd.read_csv(EMBED_PATH)
    z_cols = [c for c in df.columns if c.startswith("z")]
    log.info("Loaded embeddings: %d samples x %d dims", len(df), len(z_cols))
    return df, z_cols


def load_metadata():
    """Load modeling metadata with NAS scores and fibrosis stages."""
    meta = pd.read_csv(META_PATH)
    log.info("Loaded metadata: %d samples", len(meta))

    # Also try to load unified metadata for additional annotations
    if os.path.exists(UNIFIED_META):
        umeta = pd.read_csv(UNIFIED_META)
        # Merge extra cols not in modeling metadata
        extra_cols = [c for c in umeta.columns if c not in meta.columns and c != "sample_id"]
        if extra_cols:
            meta = meta.merge(umeta[["sample_id"] + extra_cols], on="sample_id", how="left")
            log.info("  Merged %d extra columns from unified metadata", len(extra_cols))

    return meta


def load_expression():
    """Load gene expression matrix from prepared HDF5 for transport gene analysis."""
    if not os.path.exists(H5_PATH):
        log.warning("Expression HDF5 not found: %s", H5_PATH)
        return None, None, None

    try:
        with h5py.File(H5_PATH, "r") as f:
            keys = list(f.keys())
            log.info("H5 keys: %s", keys)

            # Try common key names
            expr = None
            samples = None
            genes = None

            for k in ["expression", "X", "data", "counts"]:
                if k in f:
                    expr = f[k][:]
                    break
            for k in ["sample_ids", "samples", "obs_names", "row_names"]:
                if k in f:
                    samples = [s.decode() if isinstance(s, bytes) else s for s in f[k][:]]
                    break
            for k in ["gene_names", "genes", "var_names", "col_names"]:
                if k in f:
                    genes = [g.decode() if isinstance(g, bytes) else g for g in f[k][:]]
                    break

            if expr is not None:
                log.info("Loaded expression: %s", expr.shape)
            return expr, samples, genes
    except Exception as e:
        log.warning("Could not read HDF5: %s", e)
        return None, None, None


def load_nmf_subtypes():
    """Load NMF subtype assignments if available."""
    if not os.path.exists(NMF_PATH):
        log.warning("NMF subtype assignments not found: %s", NMF_PATH)
        return None
    try:
        df = pd.read_csv(NMF_PATH)
        log.info("Loaded NMF subtypes: %d samples", len(df))
        return df
    except Exception as e:
        log.warning("Could not load NMF subtypes: %s", e)
        return None


# ===================================================================
# 2. Define timepoints
# ===================================================================

def define_timepoints(meta, mode="nas"):
    """
    Define disease timepoints from cross-sectional metadata.

    mode = "nas":  NAS score groups  (0-1, 2-3, 4-5, 6-8)
    mode = "fib":  Fibrosis stages   (F0, F1, F2, F3, F4)
    mode = "disease": Control -> NAFL -> NASH
    """
    if mode == "nas":
        # NAS 4-group: 0=NAS 0-1, 1=NAS 2-3, 2=NAS 4-5, 3=NAS 6-8
        col = "nas_group4"
        labels = {0: "NAS_0-1", 1: "NAS_2-3", 2: "NAS_4-5", 3: "NAS_6-8"}
    elif mode == "fib":
        col = "fib_stage"
        labels = {0: "F0", 1: "F1", 2: "F2", 3: "F3", 4: "F4"}
    elif mode == "disease":
        # Map diagnosis_harmonized to ordinal
        dmap = {"Control": 0, "NAFL": 1, "NASH": 2, "Borderline": 2}
        meta = meta.copy()
        meta["_disease_ord"] = meta["diagnosis_harmonized"].map(dmap)
        col = "_disease_ord"
        labels = {0: "Control", 1: "NAFL", 2: "NASH"}
    else:
        raise ValueError(f"Unknown mode: {mode}")

    valid = meta[meta[col] >= 0].copy() if col != "_disease_ord" else meta.dropna(subset=[col]).copy()
    valid["timepoint"] = valid[col].astype(int)

    timepoints = sorted(valid["timepoint"].unique())
    tp_labels = {tp: labels.get(tp, f"T{tp}") for tp in timepoints}

    log.info("Mode=%s: %d samples in %d timepoints: %s",
             mode, len(valid), len(timepoints),
             {tp_labels[tp]: (valid["timepoint"] == tp).sum() for tp in timepoints})

    return valid, timepoints, tp_labels


# ===================================================================
# 3. Optimal transport between adjacent timepoints
# ===================================================================

def compute_ot_plans(embed_df, meta_tp, z_cols, timepoints, tp_labels, mode):
    """
    Compute optimal transport plan between each pair of adjacent timepoints
    using Hungarian algorithm on Euclidean distances in embedding space.

    For unequal group sizes, subsample the larger group to match
    (repeated for stability).
    """
    all_plans = []
    N_REPS = 10  # subsampling repetitions for stability

    for i in range(len(timepoints) - 1):
        tp_from = timepoints[i]
        tp_to = timepoints[i + 1]

        # Get samples at each timepoint
        mask_from = meta_tp["timepoint"] == tp_from
        mask_to = meta_tp["timepoint"] == tp_to

        ids_from = meta_tp.loc[mask_from, "sample_id"].values
        ids_to = meta_tp.loc[mask_to, "sample_id"].values

        # Get embeddings
        emb_from = embed_df[embed_df["sample_id"].isin(ids_from)].set_index("sample_id")
        emb_to = embed_df[embed_df["sample_id"].isin(ids_to)].set_index("sample_id")

        X_from = emb_from[z_cols].values
        X_to = emb_to[z_cols].values

        n_from = len(X_from)
        n_to = len(X_to)

        if n_from == 0 or n_to == 0:
            log.warning("Empty timepoint: %s->%s (%d->%d)",
                        tp_labels[tp_from], tp_labels[tp_to], n_from, n_to)
            continue

        n_match = min(n_from, n_to)
        log.info("OT %s -> %s: %d -> %d (matching %d)",
                 tp_labels[tp_from], tp_labels[tp_to], n_from, n_to, n_match)

        # Aggregate transport plans over subsampling repetitions
        pair_counts = {}

        for rep in range(N_REPS):
            rng = np.random.RandomState(SEED + rep)

            idx_from = rng.choice(n_from, n_match, replace=False) if n_from > n_match else np.arange(n_from)
            idx_to = rng.choice(n_to, n_match, replace=False) if n_to > n_match else np.arange(n_to)

            # Cost matrix: Euclidean distances
            cost = cdist(X_from[idx_from], X_to[idx_to], metric="euclidean")

            # Hungarian algorithm
            row_ind, col_ind = linear_sum_assignment(cost)

            for r, c in zip(row_ind, col_ind):
                sid_from = emb_from.index[idx_from[r]]
                sid_to = emb_to.index[idx_to[c]]
                key = (sid_from, sid_to)
                pair_counts[key] = pair_counts.get(key, 0) + 1

        # Record pairs with their matching frequency
        total_cost = 0.0
        n_pairs = 0
        for (sid_from, sid_to), count in pair_counts.items():
            d = np.linalg.norm(
                emb_from.loc[sid_from, z_cols].values -
                emb_to.loc[sid_to, z_cols].values
            )
            all_plans.append({
                "mode": mode,
                "tp_from": tp_labels[tp_from],
                "tp_to": tp_labels[tp_to],
                "sample_from": sid_from,
                "sample_to": sid_to,
                "distance": float(d),
                "match_frequency": count,
                "n_reps": N_REPS,
            })
            total_cost += d
            n_pairs += 1

        if n_pairs > 0:
            log.info("  Mean transport distance: %.4f (unique pairs: %d)",
                     total_cost / n_pairs, n_pairs)

    return pd.DataFrame(all_plans)


# ===================================================================
# 4. Extract transport genes
# ===================================================================

def extract_transport_genes(plans_df, expr, samples, genes, meta_tp, mode, top_n=200):
    """
    Identify genes whose expression changes most along the transport path.

    For each adjacent timepoint pair, compare expression between transported
    samples using Mann-Whitney U test.
    """
    if expr is None or samples is None or genes is None:
        log.warning("No expression data available for transport gene extraction")
        return pd.DataFrame()

    sample_to_idx = {s: i for i, s in enumerate(samples)}
    results = []

    # Get stable transport pairs (matched in >50% of subsampling reps)
    stable = plans_df[plans_df["match_frequency"] >= plans_df["n_reps"] * 0.5]

    transitions = stable.groupby(["tp_from", "tp_to"])

    for (tp_from, tp_to), group in transitions:
        from_ids = group["sample_from"].values
        to_ids = group["sample_to"].values

        # Get expression indices
        from_idx = [sample_to_idx[s] for s in from_ids if s in sample_to_idx]
        to_idx = [sample_to_idx[s] for s in to_ids if s in sample_to_idx]

        if len(from_idx) < 5 or len(to_idx) < 5:
            log.warning("Too few samples for %s -> %s: %d -> %d",
                        tp_from, tp_to, len(from_idx), len(to_idx))
            continue

        expr_from = expr[from_idx, :]
        expr_to = expr[to_idx, :]

        # Gene-level differential expression along transport
        for j, gene in enumerate(genes):
            try:
                stat, pval = mannwhitneyu(
                    expr_to[:, j], expr_from[:, j], alternative="two-sided"
                )
                mean_diff = float(np.mean(expr_to[:, j]) - np.mean(expr_from[:, j]))
            except ValueError:
                continue

            results.append({
                "mode": mode,
                "tp_from": tp_from,
                "tp_to": tp_to,
                "gene": gene,
                "mean_diff": mean_diff,
                "abs_diff": abs(mean_diff),
                "stat": float(stat),
                "pval": float(pval),
                "n_from": len(from_idx),
                "n_to": len(to_idx),
            })

    if not results:
        return pd.DataFrame()

    rdf = pd.DataFrame(results)

    # BH correction per transition
    from statsmodels.stats.multitest import multipletests

    corrected = []
    for (tp_from, tp_to), grp in rdf.groupby(["tp_from", "tp_to"]):
        pvals = grp["pval"].values
        if len(pvals) > 0:
            _, padj, _, _ = multipletests(pvals, method="fdr_bh")
            grp = grp.copy()
            grp["padj"] = padj
            corrected.append(grp)

    if corrected:
        rdf = pd.concat(corrected, ignore_index=True)
        rdf = rdf.sort_values("abs_diff", ascending=False)
        # Report top transport genes per transition
        for (tp_from, tp_to), grp in rdf.groupby(["tp_from", "tp_to"]):
            sig = grp[grp["padj"] < 0.05]
            top = grp.head(10)
            log.info("  %s -> %s: %d sig genes (padj<0.05), top: %s",
                     tp_from, tp_to, len(sig),
                     top["gene"].tolist())

    return rdf


# ===================================================================
# 5. Compute velocity field
# ===================================================================

def compute_velocity_field(embed_df, meta_tp, z_cols, timepoints, tp_labels, mode,
                           k_neighbors=15):
    """
    Estimate velocity (direction of progression) for each sample based on
    nearest neighbors at the next timepoint.

    For each sample at timepoint t:
      1. Find k nearest neighbors at timepoint t+1
      2. Velocity = mean displacement vector to those neighbors
      3. Normalize to unit length

    Samples at the last timepoint get zero velocity.
    """
    results = []

    for tp_idx, tp in enumerate(timepoints):
        mask = meta_tp["timepoint"] == tp
        ids = meta_tp.loc[mask, "sample_id"].values

        emb = embed_df[embed_df["sample_id"].isin(ids)].set_index("sample_id")
        X = emb[z_cols].values

        if tp_idx < len(timepoints) - 1:
            # Get next timepoint embeddings
            tp_next = timepoints[tp_idx + 1]
            mask_next = meta_tp["timepoint"] == tp_next
            ids_next = meta_tp.loc[mask_next, "sample_id"].values
            emb_next = embed_df[embed_df["sample_id"].isin(ids_next)].set_index("sample_id")
            X_next = emb_next[z_cols].values

            if len(X_next) == 0:
                velocities = np.zeros_like(X)
            else:
                # Compute distances to next timepoint
                dists = cdist(X, X_next, metric="euclidean")
                k = min(k_neighbors, len(X_next))

                velocities = np.zeros_like(X)
                for i in range(len(X)):
                    nn_idx = np.argpartition(dists[i], k)[:k]
                    displacement = X_next[nn_idx] - X[i]
                    velocities[i] = displacement.mean(axis=0)
        else:
            velocities = np.zeros_like(X)

        # Store per-sample velocity
        for i, sid in enumerate(emb.index):
            v = velocities[i]
            v_mag = float(np.linalg.norm(v))
            v_norm = v / v_mag if v_mag > 0 else v

            row = {
                "sample_id": sid,
                "mode": mode,
                "timepoint": int(tp),
                "tp_label": tp_labels[tp],
                "velocity_magnitude": v_mag,
            }
            # Store normalized velocity components (first 3 for vis)
            for d in range(min(len(z_cols), 64)):
                row[f"v_{d}"] = float(v_norm[d])
            results.append(row)

    return pd.DataFrame(results)


# ===================================================================
# 6. UMAP for visualization
# ===================================================================

def compute_umap_with_velocity(embed_df, velocity_df, z_cols, meta_tp, mode):
    """
    Compute UMAP embedding and project velocities for visualization.
    """
    from sklearn.decomposition import PCA
    try:
        from umap import UMAP
    except ImportError:
        log.warning("UMAP not available; using PCA projection instead")
        UMAP = None

    # Get all samples in this mode
    valid_ids = meta_tp["sample_id"].values
    emb = embed_df[embed_df["sample_id"].isin(valid_ids)].copy()
    emb = emb.set_index("sample_id")
    X = emb[z_cols].values

    if UMAP is not None:
        log.info("Computing UMAP for %d samples...", len(X))
        reducer = UMAP(n_components=2, n_neighbors=30, min_dist=0.3,
                       metric="euclidean", random_state=SEED)
        umap_coords = reducer.fit_transform(X)
    else:
        log.info("Computing PCA for %d samples...", len(X))
        pca = PCA(n_components=2, random_state=SEED)
        umap_coords = pca.fit_transform(X)

    # Project velocities into UMAP space using finite differences
    # For each sample, project velocity by computing UMAP displacement
    v_cols = [f"v_{d}" for d in range(len(z_cols))]
    vel_lookup = velocity_df.set_index("sample_id")

    results = []
    for i, sid in enumerate(emb.index):
        row = {
            "sample_id": sid,
            "mode": mode,
            "umap_1": float(umap_coords[i, 0]),
            "umap_2": float(umap_coords[i, 1]),
        }

        # Add metadata
        meta_row = meta_tp[meta_tp["sample_id"] == sid]
        if len(meta_row) > 0:
            row["timepoint"] = int(meta_row.iloc[0]["timepoint"])

        # Project velocity into UMAP space via Jacobian approximation
        if sid in vel_lookup.index:
            vel_row = vel_lookup.loc[sid]
            v_mag = vel_row.get("velocity_magnitude", 0.0)
            row["velocity_magnitude"] = float(v_mag)

            if v_mag > 0 and UMAP is not None:
                # Approximate: displace in embedding space, project, compute diff
                v_full = np.array([vel_row.get(f"v_{d}", 0.0)
                                   for d in range(len(z_cols))])
                epsilon = 0.1
                x_displaced = X[i] + epsilon * v_full
                umap_displaced = reducer.transform(x_displaced.reshape(1, -1))[0]
                arrow = (umap_displaced - umap_coords[i]) / epsilon
                row["arrow_u"] = float(arrow[0])
                row["arrow_v"] = float(arrow[1])
            else:
                row["arrow_u"] = 0.0
                row["arrow_v"] = 0.0
        else:
            row["velocity_magnitude"] = 0.0
            row["arrow_u"] = 0.0
            row["arrow_v"] = 0.0

        results.append(row)

    return pd.DataFrame(results)


# ===================================================================
# 7. Bifurcation analysis
# ===================================================================

def analyze_bifurcation(velocity_df, meta_tp, embed_df, z_cols, mode, nmf_df=None):
    """
    Identify bifurcation points where the velocity field diverges.

    Strategy:
      1. For each timepoint, compute velocity field divergence
      2. High divergence = potential bifurcation
      3. If NMF subtypes available, test whether subtypes correspond
         to distinct trajectory branches
    """
    results = []

    timepoints = sorted(meta_tp["timepoint"].unique())

    for tp in timepoints:
        tp_vel = velocity_df[
            (velocity_df["mode"] == mode) &
            (velocity_df["timepoint"] == tp)
        ].copy()

        if len(tp_vel) < 10:
            continue

        # Get velocity vectors
        v_cols = [f"v_{d}" for d in range(len(z_cols))]
        available_vcols = [c for c in v_cols if c in tp_vel.columns]
        V = tp_vel[available_vcols].values

        # Compute velocity divergence as variance of velocity directions
        # (high variance = divergent trajectories = potential bifurcation)
        norms = np.linalg.norm(V, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-10)
        V_unit = V / norms

        # Angular dispersion: mean pairwise cosine similarity
        # (lower = more divergent)
        n = len(V_unit)
        if n > 1:
            # Subsample for efficiency
            max_pairs = 5000
            if n * (n - 1) // 2 > max_pairs:
                idx = np.random.choice(n, min(n, 200), replace=False)
                V_sub = V_unit[idx]
            else:
                V_sub = V_unit

            cos_sim = V_sub @ V_sub.T
            # Upper triangle only
            triu_idx = np.triu_indices(len(V_sub), k=1)
            cos_vals = cos_sim[triu_idx]
            mean_cos = float(np.mean(cos_vals))
            std_cos = float(np.std(cos_vals))
        else:
            mean_cos = 1.0
            std_cos = 0.0

        divergence = 1.0 - mean_cos  # higher = more divergent

        row = {
            "mode": mode,
            "timepoint": int(tp),
            "tp_label": tp_vel["tp_label"].iloc[0] if "tp_label" in tp_vel.columns else f"T{tp}",
            "n_samples": len(tp_vel),
            "mean_velocity_mag": float(tp_vel["velocity_magnitude"].mean()),
            "mean_cosine_similarity": mean_cos,
            "std_cosine_similarity": std_cos,
            "divergence_score": divergence,
        }

        # NMF subtype enrichment at this timepoint
        if nmf_df is not None:
            tp_ids = tp_vel["sample_id"].values
            nmf_tp = nmf_df[nmf_df["sample_id"].isin(tp_ids)]
            if len(nmf_tp) > 0 and "subtype" in nmf_tp.columns:
                subtype_counts = nmf_tp["subtype"].value_counts()
                row["n_subtype1"] = int(subtype_counts.get(1, subtype_counts.get("S1", 0)))
                row["n_subtype2"] = int(subtype_counts.get(2, subtype_counts.get("S2", 0)))

                # Test velocity divergence between subtypes
                s1_ids = nmf_tp[nmf_tp["subtype"].isin([1, "S1"])]["sample_id"].values
                s2_ids = nmf_tp[nmf_tp["subtype"].isin([2, "S2"])]["sample_id"].values

                s1_vel = tp_vel[tp_vel["sample_id"].isin(s1_ids)]
                s2_vel = tp_vel[tp_vel["sample_id"].isin(s2_ids)]

                if len(s1_vel) >= 3 and len(s2_vel) >= 3:
                    V1 = s1_vel[available_vcols].values
                    V2 = s2_vel[available_vcols].values
                    mean_v1 = V1.mean(axis=0)
                    mean_v2 = V2.mean(axis=0)
                    cos_between = float(
                        np.dot(mean_v1, mean_v2) /
                        (np.linalg.norm(mean_v1) * np.linalg.norm(mean_v2) + 1e-10)
                    )
                    row["subtype_velocity_cosine"] = cos_between
                    row["subtype_divergent"] = cos_between < 0.5

                    # Mann-Whitney on velocity magnitudes
                    stat, pval = mannwhitneyu(
                        s1_vel["velocity_magnitude"].values,
                        s2_vel["velocity_magnitude"].values,
                        alternative="two-sided"
                    )
                    row["subtype_velocity_pval"] = float(pval)

        results.append(row)

    rdf = pd.DataFrame(results)

    # Flag the most divergent timepoint as primary bifurcation
    if len(rdf) > 0:
        rdf["is_bifurcation"] = rdf["divergence_score"] == rdf["divergence_score"].max()
        max_div = rdf.loc[rdf["divergence_score"].idxmax()]
        log.info("  Primary bifurcation at %s (divergence=%.4f, n=%d)",
                 max_div.get("tp_label", "?"),
                 max_div["divergence_score"],
                 max_div["n_samples"])

    return rdf


# ===================================================================
# 8. Main
# ===================================================================

def main():
    t0 = time.time()

    # Load data
    log.info("Loading data...")
    embed_df, z_cols = load_embeddings()
    meta = load_metadata()
    expr, expr_samples, expr_genes = load_expression()
    nmf_df = load_nmf_subtypes()

    # Merge embeddings with metadata
    embed_meta = embed_df.merge(meta, on="sample_id", how="inner")
    log.info("Merged: %d samples with both embeddings and metadata", len(embed_meta))

    # Run for multiple trajectory definitions
    all_transport = []
    all_velocity = []
    all_bifurcation = []
    all_umap = []

    for mode in ["nas", "fib", "disease"]:
        log.info("\n" + "=" * 60)
        log.info("Mode: %s", mode.upper())
        log.info("=" * 60)

        try:
            meta_tp, timepoints, tp_labels = define_timepoints(embed_meta, mode=mode)
        except Exception as e:
            log.warning("Could not define timepoints for mode=%s: %s", mode, e)
            continue

        if len(timepoints) < 2:
            log.warning("Need >= 2 timepoints for mode=%s, got %d", mode, len(timepoints))
            continue

        # 3. Optimal transport
        log.info("\nComputing optimal transport plans...")
        transport = compute_ot_plans(embed_df, meta_tp, z_cols, timepoints, tp_labels, mode)
        if len(transport) > 0:
            all_transport.append(transport)
            log.info("Transport plans: %d pairs", len(transport))

        # 4. Transport genes (only if expression available)
        if expr is not None and len(transport) > 0:
            log.info("\nExtracting transport genes...")
            transport_genes = extract_transport_genes(
                transport, expr, expr_samples, expr_genes, meta_tp, mode
            )
            if len(transport_genes) > 0:
                # Save transport genes separately per mode
                tg_path = os.path.join(
                    OUTDIR, f"progression_transport_genes_{mode}.csv"
                )
                transport_genes.to_csv(tg_path, index=False)
                log.info("Saved transport genes: %s (%d genes)", tg_path, len(transport_genes))

        # 5. Velocity field
        log.info("\nComputing velocity field...")
        velocity = compute_velocity_field(
            embed_df, meta_tp, z_cols, timepoints, tp_labels, mode
        )
        if len(velocity) > 0:
            all_velocity.append(velocity)
            log.info("Velocity field: %d samples", len(velocity))

        # 6. UMAP + velocity arrows
        log.info("\nComputing UMAP with velocity arrows...")
        umap_data = compute_umap_with_velocity(
            embed_df, velocity, z_cols, meta_tp, mode
        )
        if len(umap_data) > 0:
            all_umap.append(umap_data)

        # 7. Bifurcation analysis
        log.info("\nAnalyzing bifurcation points...")
        bifurcation = analyze_bifurcation(
            velocity, meta_tp, embed_df, z_cols, mode, nmf_df
        )
        if len(bifurcation) > 0:
            all_bifurcation.append(bifurcation)

    # --- Save combined results ---
    log.info("\n" + "=" * 60)
    log.info("Saving results...")
    log.info("=" * 60)

    if all_transport:
        transport_combined = pd.concat(all_transport, ignore_index=True)
        transport_combined.to_csv(OUT_TRANSPORT, index=False)
        log.info("Saved: %s (%d rows)", OUT_TRANSPORT, len(transport_combined))

    if all_velocity:
        velocity_combined = pd.concat(all_velocity, ignore_index=True)
        velocity_combined.to_csv(OUT_VELOCITY, index=False)
        log.info("Saved: %s (%d rows)", OUT_VELOCITY, len(velocity_combined))

    if all_bifurcation:
        bif_combined = pd.concat(all_bifurcation, ignore_index=True)
        bif_combined.to_csv(OUT_BIFURCATION, index=False)
        log.info("Saved: %s (%d rows)", OUT_BIFURCATION, len(bif_combined))

        # Print bifurcation summary
        log.info("\n--- Bifurcation Summary ---")
        for _, row in bif_combined.iterrows():
            bif_flag = "*" if row.get("is_bifurcation", False) else " "
            log.info("  %s [%s] %s: divergence=%.4f, v_mag=%.4f (n=%d)",
                     bif_flag, row["mode"], row["tp_label"],
                     row["divergence_score"], row["mean_velocity_mag"],
                     row["n_samples"])

    if all_umap:
        umap_combined = pd.concat(all_umap, ignore_index=True)
        umap_combined.to_csv(OUT_UMAP, index=False)
        log.info("Saved: %s (%d rows)", OUT_UMAP, len(umap_combined))

    elapsed = time.time() - t0
    log.info("\n" + "=" * 70)
    log.info("110_disease_progression.py completed in %.1f seconds", elapsed)
    log.info("=" * 70)


if __name__ == "__main__":
    main()
