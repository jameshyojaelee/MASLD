#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_274_sens
#SBATCH --output=logs/net_274_sens_%j.out
#SBATCH --error=logs/net_274_sens_%j.err
# ===========================================================================
# Script 274: Sensitivity & Stability Analyses for the Multi-Evidence Network
# ===========================================================================
# Purpose
# -------
# Address reviewer concerns that the pipeline's pathway-Jaccard cutoff (0.15),
# tail-transform quantile (0.90), Leiden seed, and layer portfolio were chosen
# ad-hoc. This script performs five targeted sensitivity / stability analyses:
#
#   A. Pathway Jaccard cutoff sweep (0.05, 0.10, 0.15, 0.20, 0.30)
#   B. Tail-transform quantile sweep (0.80, 0.85, 0.90, 0.95, 0.98)
#   C. Leiden seed stability (10 seeds, pairwise ARI)
#   D. Edge bootstrap stability (100 iterations, co-membership consensus)
#   E. Layer-drop-one analysis (ppi, coexpr, regulon, lr, genetic, pathway, cerna)
#
# For A/B the composite is re-computed on three layers only (pathway + coexpr +
# ppi) to isolate the pathway transform's effect. For E the composite is
# recomputed on six of the seven active layers.
#
# Metrics
# -------
#   - composite K>=2 count      (multiplicity >= 2 across active layers)
#   - composite P>0.9 count     (high-confidence edges)
#   - clinical-drug fold enrichment (enrichment of clinical-drug-target pairs
#     among high-confidence edges vs. base rate in the node set)
#
# Inputs
# ------
#   - edges_{layer}.csv  (symbol-keyed after 273_normalize_ids.py)
#   - network_nodes.csv
#   - positive_control.csv  (subset where Clinical_Drug is non-empty)
#
# Outputs
# -------
#   figures/misc/network_sensitivity/
#     data/A_jaccard_sweep.csv
#     data/B_quantile_sweep.csv
#     data/C_leiden_seed_ari.csv
#     data/D_bootstrap_consensus.csv
#     data/E_layer_drop_one.csv
#     A_jaccard_sweep.pdf
#     B_quantile_sweep.pdf
#     C_leiden_ari_heatmap.pdf
#     D_consensus_histogram.pdf
#     E_layer_drop_one.pdf
#
# Environment: spatial (pandas, numpy, python-igraph, leidenalg, scikit-learn,
#              matplotlib, tqdm)
# ===========================================================================

from __future__ import annotations

import os
import sys
import time
import warnings
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# -----------------------------------------------------------------------------
# Paths & config
# -----------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
NET_DIR = BASE / "RNA-seq" / "results" / "network"
NODES_CSV = NET_DIR / "network_nodes.csv"
POS_CTRL_CSV = BASE / "results" / "library" / "positive_control.csv"
ID_NORM_SUMMARY = NET_DIR / "id_normalization_summary.csv"

OUT_DIR = BASE / "figures" / "misc" / "network_sensitivity"
DATA_DIR = OUT_DIR / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

ACTIVE_LAYERS = ["ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna"]

# Defaults inherited from 262c_confidence_scores.py
DEFAULT_PATHWAY_JACCARD_MIN = 0.15
DEFAULT_TAIL_QUANTILE = 0.90
REGULON_SCORE_CAP = 0.9

# Sensitivity sweeps
JACCARD_CUTOFFS = [0.05, 0.10, 0.15, 0.20, 0.30]
TAIL_QUANTILES = [0.80, 0.85, 0.90, 0.95, 0.98]

# Composite thresholds for sensitivity metrics
COMPOSITE_P_THRESHOLD = 0.9
COMPOSITE_K_THRESHOLD = 2

# Leiden / bootstrap stability
LEIDEN_RESOLUTION = 1.0
LEIDEN_SEEDS = list(range(10))          # 10 seeds
N_BOOTSTRAP = 100
BOOTSTRAP_FRAC = 0.80
BOOTSTRAP_CONSENSUS_THRESHOLD = 0.70
LEIDEN_EDGE_P_THRESHOLD = 0.5           # Build graph on edges with P>0.5

N_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))

# For A/B we isolate pathway effect by using pathway + the two largest
# non-pathway layers (coexpr, ppi)
BIG_THREE = ["ppi", "coexpr", "pathway"]

t0 = time.time()


def log(msg: str) -> None:
    elapsed = time.time() - t0
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# -----------------------------------------------------------------------------
# Preconditions
# -----------------------------------------------------------------------------
def check_preconditions() -> None:
    if not ID_NORM_SUMMARY.exists():
        print("ERROR: Run 273_normalize_ids.py first")
        sys.exit(1)
    for layer in ACTIVE_LAYERS:
        p = NET_DIR / f"edges_{layer}.csv"
        if not p.exists():
            print(f"ERROR: Missing edge file {p}")
            sys.exit(1)
    if not NODES_CSV.exists():
        print(f"ERROR: Missing {NODES_CSV}")
        sys.exit(1)
    if not POS_CTRL_CSV.exists():
        print(f"ERROR: Missing {POS_CTRL_CSV}")
        sys.exit(1)


# -----------------------------------------------------------------------------
# Core confidence transforms (mirrors 262c with exposed knobs for A/B)
# -----------------------------------------------------------------------------
def tail_transform(scores: np.ndarray, cutoff_q: float) -> np.ndarray:
    scores = np.asarray(scores, dtype=np.float64)
    if scores.size == 0:
        return scores
    q = float(np.quantile(scores, cutoff_q))
    mx = float(np.max(scores))
    if mx - q <= 1e-12:
        return (scores >= q).astype(np.float64)
    return np.clip((scores - q) / (mx - q), 0.0, 1.0)


def confidence_for_layer(
    df: pd.DataFrame,
    name: str,
    pathway_jaccard_min: float = DEFAULT_PATHWAY_JACCARD_MIN,
    tail_quantile: float = DEFAULT_TAIL_QUANTILE,
) -> pd.DataFrame:
    """Compute per-layer confidence score using the same transforms as 262c."""
    df = df.copy()
    df["raw_score"] = (
        pd.to_numeric(df["raw_score"], errors="coerce")
        .fillna(0.0)
        .clip(0.0, 1.0)
    )

    if name in ("ppi", "coexpr", "genetic", "cerna", "lr"):
        df["confidence_score"] = df["raw_score"].values
    elif name == "regulon":
        df["confidence_score"] = df["raw_score"].clip(upper=REGULON_SCORE_CAP).values
    elif name == "pathway":
        df = df[df["raw_score"] >= pathway_jaccard_min].reset_index(drop=True)
        if len(df) > 0:
            df["confidence_score"] = tail_transform(df["raw_score"].values, tail_quantile)
        else:
            df["confidence_score"] = df["raw_score"].values
    else:
        df["confidence_score"] = df["raw_score"].values

    keep = ["gene_a", "gene_b", "confidence_score"]
    return df[keep]


def load_layer_raw(name: str) -> pd.DataFrame:
    """Load raw edge file, canonicalize (gene_a, gene_b) so a<b, drop self-loops."""
    path = NET_DIR / f"edges_{name}.csv"
    df = pd.read_csv(path, usecols=["gene_a", "gene_b", "raw_score"])
    df = df.dropna(subset=["gene_a", "gene_b"])
    df = df[df["gene_a"] != df["gene_b"]]

    # Canonical ordering so (A,B) and (B,A) collide on merge
    a = np.where(df["gene_a"].values < df["gene_b"].values,
                 df["gene_a"].values, df["gene_b"].values)
    b = np.where(df["gene_a"].values < df["gene_b"].values,
                 df["gene_b"].values, df["gene_a"].values)
    df = df.assign(gene_a=a, gene_b=b)

    # Dedup within-layer (keep max raw_score if duplicates arise from canonicalization)
    df = (
        df.groupby(["gene_a", "gene_b"], as_index=False, sort=False)["raw_score"]
          .max()
    )
    return df


def load_confidence_layers(
    layers: list[str],
    pathway_jaccard_min: float = DEFAULT_PATHWAY_JACCARD_MIN,
    tail_quantile: float = DEFAULT_TAIL_QUANTILE,
) -> dict[str, pd.DataFrame]:
    out = {}
    for name in layers:
        raw = load_layer_raw(name)
        conf = confidence_for_layer(
            raw, name,
            pathway_jaccard_min=pathway_jaccard_min,
            tail_quantile=tail_quantile,
        )
        conf = conf.rename(columns={"confidence_score": f"p_{name}"})
        out[name] = conf
    return out


def merge_layers(layer_data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    merged: pd.DataFrame | None = None
    for name, df in layer_data.items():
        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on=["gene_a", "gene_b"], how="outer")
    for name in layer_data:
        col = f"p_{name}"
        if col in merged.columns:
            merged[col] = merged[col].fillna(0.0)
        else:
            merged[col] = 0.0
    return merged


def compute_composite(merged: pd.DataFrame, layers: list[str]) -> pd.DataFrame:
    cols = [f"p_{n}" for n in layers if f"p_{n}" in merged.columns]
    M = merged[cols].values.astype(np.float64)
    log_c = np.log1p(-np.clip(M, 0.0, 1.0 - 1e-15))
    p_composite = 1.0 - np.exp(log_c.sum(axis=1))
    p_composite = np.clip(p_composite, 0.0, 1.0)
    k_mult = (M > 0.5).sum(axis=1).astype(np.int32)
    merged = merged.assign(p_composite=p_composite, k_multiplicity=k_mult)
    return merged


# -----------------------------------------------------------------------------
# Clinical drug target pair set
# -----------------------------------------------------------------------------
def load_clinical_drug_genes() -> set[str]:
    pc = pd.read_csv(POS_CTRL_CSV)
    if "Clinical_Drug" not in pc.columns or "Gene symbol" not in pc.columns:
        log("  WARNING: positive_control.csv missing Clinical_Drug/Gene symbol cols")
        return set()
    clinical = pc[pc["Clinical_Drug"].astype(str).str.strip().replace("nan", "") != ""]
    genes = set(clinical["Gene symbol"].dropna().astype(str).str.strip().tolist())
    genes.discard("")
    return genes


def clinical_fold_enrichment(
    merged: pd.DataFrame,
    clinical_genes: set[str],
    node_universe: set[str],
    composite_mask: np.ndarray,
) -> tuple[float, int, int]:
    """Fold enrichment of clinical-drug-target pairs among high-confidence edges.

    Background rate = expected fraction of random gene-gene pairs (from
    node_universe) that involve a clinical-drug gene.
    """
    if len(clinical_genes) == 0 or len(node_universe) == 0:
        return float("nan"), 0, 0

    n_nodes = len(node_universe)
    n_clin = len(clinical_genes & node_universe)
    if n_nodes < 2 or n_clin == 0:
        return float("nan"), 0, n_clin

    # Background: probability a random unordered pair has >=1 clinical gene
    n_pairs_total = n_nodes * (n_nodes - 1) / 2
    n_pairs_no_clin = (n_nodes - n_clin) * (n_nodes - n_clin - 1) / 2
    bg_rate = 1.0 - (n_pairs_no_clin / n_pairs_total) if n_pairs_total > 0 else 0.0
    if bg_rate <= 0:
        return float("nan"), 0, n_clin

    sub = merged.loc[composite_mask, ["gene_a", "gene_b"]]
    if len(sub) == 0:
        return 0.0, 0, n_clin

    ga = sub["gene_a"].values
    gb = sub["gene_b"].values
    hit_mask = np.array(
        [(a in clinical_genes) or (b in clinical_genes) for a, b in zip(ga, gb)],
        dtype=bool,
    )
    n_hits = int(hit_mask.sum())
    obs_rate = n_hits / len(sub)
    fe = obs_rate / bg_rate
    return float(fe), n_hits, n_clin


# -----------------------------------------------------------------------------
# Analysis A: Pathway Jaccard cutoff sweep (isolated on ppi+coexpr+pathway)
# -----------------------------------------------------------------------------
def analysis_A_jaccard_sweep(clinical_genes: set[str], node_universe: set[str]) -> pd.DataFrame:
    log("\n=== Analysis A: Pathway Jaccard cutoff sweep ===")
    rows = []
    for cutoff in JACCARD_CUTOFFS:
        log(f"  Jaccard cutoff = {cutoff}")
        layer_data = load_confidence_layers(
            BIG_THREE,
            pathway_jaccard_min=cutoff,
            tail_quantile=DEFAULT_TAIL_QUANTILE,
        )
        merged = merge_layers(layer_data)
        merged = compute_composite(merged, BIG_THREE)

        n_k_ge2 = int((merged["k_multiplicity"] >= COMPOSITE_K_THRESHOLD).sum())
        mask_p = (merged["p_composite"] > COMPOSITE_P_THRESHOLD).values
        n_p_gt09 = int(mask_p.sum())
        fe, n_hits, n_clin = clinical_fold_enrichment(
            merged, clinical_genes, node_universe, mask_p
        )
        rows.append({
            "jaccard_cutoff": cutoff,
            "n_edges_total": len(merged),
            "n_composite_K_ge2": n_k_ge2,
            "n_composite_P_gt0p9": n_p_gt09,
            "clinical_fold_enrichment": fe,
            "n_clinical_hits": n_hits,
            "n_clinical_genes_in_nodes": n_clin,
        })
        log(
            f"    edges={len(merged):,}  K>=2={n_k_ge2:,}  "
            f"P>0.9={n_p_gt09:,}  fold-enrich(clinical)={fe:.3f}"
        )
    return pd.DataFrame(rows)


def plot_A(df: pd.DataFrame, out_pdf: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    ax.plot(df["jaccard_cutoff"], df["n_composite_P_gt0p9"], marker="o", label="P>0.9")
    ax.plot(df["jaccard_cutoff"], df["n_composite_K_ge2"], marker="s", label="K>=2")
    ax.set_xlabel("Pathway Jaccard cutoff")
    ax.set_ylabel("Composite edges (log scale)")
    ax.set_yscale("log")
    ax.set_title("Composite edge count vs Jaccard cutoff")
    ax.legend(frameon=False)
    ax.axvline(DEFAULT_PATHWAY_JACCARD_MIN, ls="--", c="grey", alpha=0.6)

    ax = axes[1]
    ax.plot(df["jaccard_cutoff"], df["clinical_fold_enrichment"], marker="o", color="C3")
    ax.set_xlabel("Pathway Jaccard cutoff")
    ax.set_ylabel("Clinical drug target fold enrichment")
    ax.set_title("Fold enrichment (clinical) vs Jaccard cutoff")
    ax.axvline(DEFAULT_PATHWAY_JACCARD_MIN, ls="--", c="grey", alpha=0.6)
    ax.axhline(1.0, ls=":", c="black", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Analysis B: Tail quantile sweep
# -----------------------------------------------------------------------------
def analysis_B_quantile_sweep(clinical_genes: set[str], node_universe: set[str]) -> pd.DataFrame:
    log("\n=== Analysis B: Tail-transform quantile sweep ===")
    rows = []
    for q in TAIL_QUANTILES:
        log(f"  Tail quantile = {q}")
        layer_data = load_confidence_layers(
            BIG_THREE,
            pathway_jaccard_min=DEFAULT_PATHWAY_JACCARD_MIN,
            tail_quantile=q,
        )
        merged = merge_layers(layer_data)
        merged = compute_composite(merged, BIG_THREE)

        n_k_ge2 = int((merged["k_multiplicity"] >= COMPOSITE_K_THRESHOLD).sum())
        mask_p = (merged["p_composite"] > COMPOSITE_P_THRESHOLD).values
        n_p_gt09 = int(mask_p.sum())
        fe, n_hits, n_clin = clinical_fold_enrichment(
            merged, clinical_genes, node_universe, mask_p
        )
        rows.append({
            "tail_quantile": q,
            "n_edges_total": len(merged),
            "n_composite_K_ge2": n_k_ge2,
            "n_composite_P_gt0p9": n_p_gt09,
            "clinical_fold_enrichment": fe,
            "n_clinical_hits": n_hits,
            "n_clinical_genes_in_nodes": n_clin,
        })
        log(
            f"    edges={len(merged):,}  K>=2={n_k_ge2:,}  "
            f"P>0.9={n_p_gt09:,}  fold-enrich(clinical)={fe:.3f}"
        )
    return pd.DataFrame(rows)


def plot_B(df: pd.DataFrame, out_pdf: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    ax.plot(df["tail_quantile"], df["n_composite_P_gt0p9"], marker="o", label="P>0.9")
    ax.plot(df["tail_quantile"], df["n_composite_K_ge2"], marker="s", label="K>=2")
    ax.set_xlabel("Tail quantile")
    ax.set_ylabel("Composite edges (log scale)")
    ax.set_yscale("log")
    ax.set_title("Composite edge count vs tail quantile")
    ax.legend(frameon=False)
    ax.axvline(DEFAULT_TAIL_QUANTILE, ls="--", c="grey", alpha=0.6)

    ax = axes[1]
    ax.plot(df["tail_quantile"], df["clinical_fold_enrichment"], marker="o", color="C3")
    ax.set_xlabel("Tail quantile")
    ax.set_ylabel("Clinical drug target fold enrichment")
    ax.set_title("Fold enrichment (clinical) vs tail quantile")
    ax.axvline(DEFAULT_TAIL_QUANTILE, ls="--", c="grey", alpha=0.6)
    ax.axhline(1.0, ls=":", c="black", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Build composite graph at default parameters for C/D/E
# -----------------------------------------------------------------------------
def build_default_composite(
    layers: list[str] = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Compute composite on the given layer set (default: all active)."""
    if layers is None:
        layers = ACTIVE_LAYERS
    layer_data = load_confidence_layers(
        layers,
        pathway_jaccard_min=DEFAULT_PATHWAY_JACCARD_MIN,
        tail_quantile=DEFAULT_TAIL_QUANTILE,
    )
    merged = merge_layers(layer_data)
    merged = compute_composite(merged, layers)
    return merged, layers


def build_igraph_from_edges(
    edges: pd.DataFrame,
    p_threshold: float = LEIDEN_EDGE_P_THRESHOLD,
) -> "ig.Graph":
    """Build a weighted igraph from composite edges above p_threshold."""
    import igraph as ig

    sub = edges.loc[edges["p_composite"] > p_threshold, ["gene_a", "gene_b", "p_composite"]]
    if len(sub) == 0:
        raise RuntimeError(f"No edges above P>{p_threshold}")

    vertices = pd.unique(np.concatenate([sub["gene_a"].values, sub["gene_b"].values]))
    v2i = {v: i for i, v in enumerate(vertices)}
    src = sub["gene_a"].map(v2i).values
    dst = sub["gene_b"].map(v2i).values
    weights = sub["p_composite"].values.astype(np.float64)

    g = ig.Graph(n=len(vertices), edges=list(zip(src, dst)), directed=False)
    g.vs["name"] = list(vertices)
    g.es["weight"] = weights.tolist()
    return g


# -----------------------------------------------------------------------------
# Analysis C: Leiden seed stability
# -----------------------------------------------------------------------------
def _run_leiden_worker(args):
    """Pickleable worker: run Leiden on a given graph file with a given seed."""
    import igraph as ig
    import leidenalg

    graphml_path, seed, resolution = args
    g = ig.Graph.Read_GraphML(graphml_path)
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=int(seed),
        n_iterations=-1,
    )
    return seed, list(partition.membership)


def analysis_C_leiden_seed(g) -> pd.DataFrame:
    log("\n=== Analysis C: Leiden seed stability ===")

    # Persist the graph to disk so workers can reread it (avoids pickling Graph)
    tmp_graphml = DATA_DIR / "_tmp_leiden_graph.graphml"
    g.write_graphml(str(tmp_graphml))

    tasks = [(str(tmp_graphml), s, LEIDEN_RESOLUTION) for s in LEIDEN_SEEDS]

    memberships: dict[int, list[int]] = {}
    n_workers = min(N_WORKERS, len(tasks))
    log(f"  Running Leiden with {len(tasks)} seeds using {n_workers} workers...")
    with Pool(processes=n_workers) as pool:
        for seed, memb in tqdm(
            pool.imap_unordered(_run_leiden_worker, tasks),
            total=len(tasks),
            desc="leiden-seeds",
        ):
            memberships[seed] = memb

    seeds_sorted = sorted(memberships.keys())
    n = len(seeds_sorted)
    ari_mat = np.eye(n, dtype=np.float64)
    for i, j in combinations(range(n), 2):
        s_i, s_j = seeds_sorted[i], seeds_sorted[j]
        ari = adjusted_rand_score(memberships[s_i], memberships[s_j])
        ari_mat[i, j] = ari_mat[j, i] = ari

    iu = np.triu_indices(n, k=1)
    off_diag = ari_mat[iu]
    log(
        f"  ARI off-diagonal: mean={off_diag.mean():.4f}  "
        f"min={off_diag.min():.4f}  std={off_diag.std():.4f}"
    )

    df = pd.DataFrame(ari_mat, index=seeds_sorted, columns=seeds_sorted)
    df.index.name = "seed"

    # Cleanup tmp graph
    try:
        tmp_graphml.unlink()
    except Exception:
        pass

    return df


def plot_C(ari_df: pd.DataFrame, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(ari_df.values, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(ari_df.columns)))
    ax.set_yticks(range(len(ari_df.index)))
    ax.set_xticklabels(ari_df.columns)
    ax.set_yticklabels(ari_df.index)
    ax.set_xlabel("Seed")
    ax.set_ylabel("Seed")
    iu = np.triu_indices(len(ari_df), k=1)
    off_diag = ari_df.values[iu]
    ax.set_title(
        f"Leiden seed stability (ARI)\n"
        f"mean={off_diag.mean():.3f}, min={off_diag.min():.3f}, std={off_diag.std():.3f}"
    )
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Adjusted Rand Index")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Analysis D: Edge bootstrap stability
# -----------------------------------------------------------------------------
def _bootstrap_worker(args):
    """Worker for D: resample edges, run Leiden, return (iter_id, gene_to_community)."""
    import igraph as ig
    import leidenalg

    (
        iter_id,
        vertices_pkl,
        src_arr,
        dst_arr,
        weights_arr,
        n_edges_keep,
        seed,
        resolution,
    ) = args

    rng = np.random.default_rng(seed)
    # Sample with replacement (80% of edge count) — note: duplicates aggregate as weight sum
    idx = rng.choice(len(src_arr), size=n_edges_keep, replace=True)
    s = src_arr[idx]
    d = dst_arr[idx]
    w = weights_arr[idx]

    # Aggregate duplicate (s,d) into summed weight
    key = np.minimum(s, d).astype(np.int64) * (len(vertices_pkl) + 1) + np.maximum(s, d).astype(np.int64)
    order = np.argsort(key)
    key_s = key[order]
    s_s = np.minimum(s, d)[order]
    d_s = np.maximum(s, d)[order]
    w_s = w[order]

    # Reduce duplicates
    uniq_mask = np.concatenate(([True], key_s[1:] != key_s[:-1]))
    starts = np.nonzero(uniq_mask)[0]
    ends = np.concatenate([starts[1:], [len(key_s)]])
    u_src = s_s[starts]
    u_dst = d_s[starts]
    u_w = np.add.reduceat(w_s, starts)
    # guard (should not trigger but safe)
    _ = ends

    # Drop self-loops if any
    keep = u_src != u_dst
    u_src, u_dst, u_w = u_src[keep], u_dst[keep], u_w[keep]

    g = ig.Graph(n=len(vertices_pkl), edges=list(zip(u_src.tolist(), u_dst.tolist())), directed=False)
    g.es["weight"] = u_w.tolist()

    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=int(seed),
        n_iterations=-1,
    )
    return iter_id, np.asarray(partition.membership, dtype=np.int32)


def analysis_D_bootstrap(g) -> pd.DataFrame:
    log("\n=== Analysis D: Edge bootstrap stability ===")

    import igraph as ig  # noqa: F401  (import check)

    vertices = list(g.vs["name"])
    n_v = len(vertices)
    edges = np.array(g.get_edgelist(), dtype=np.int64)
    src = edges[:, 0]
    dst = edges[:, 1]
    weights = np.asarray(g.es["weight"], dtype=np.float64)
    n_edges = len(src)
    n_keep = int(round(n_edges * BOOTSTRAP_FRAC))
    log(f"  Graph: {n_v:,} vertices, {n_edges:,} edges; bootstrap keep = {n_keep:,}")

    tasks = [
        (i, vertices, src, dst, weights, n_keep, 10_000 + i, LEIDEN_RESOLUTION)
        for i in range(N_BOOTSTRAP)
    ]

    # Co-membership accumulator across bootstraps, summarized per gene by its
    # mode community co-membership fraction (see "consensus_fraction" below).
    # To keep memory bounded we compute, for each gene, the distribution of
    # community IDs it receives across bootstraps and report its most-common
    # relative frequency.

    community_ids = np.zeros((n_v, N_BOOTSTRAP), dtype=np.int32)

    n_workers = min(N_WORKERS, N_BOOTSTRAP)
    log(f"  Running {N_BOOTSTRAP} bootstrap iterations using {n_workers} workers...")
    with Pool(processes=n_workers) as pool:
        for iter_id, memb in tqdm(
            pool.imap_unordered(_bootstrap_worker, tasks),
            total=len(tasks),
            desc="bootstrap",
        ):
            community_ids[:, iter_id] = memb

    # For each gene, compute mode-frequency (dominant community across bootstraps).
    # This is equivalent to the max over communities of "fraction of bootstraps
    # in which gene g belongs to community c" — i.e. co-membership fraction
    # with its own modal community.
    dominant_freq = np.zeros(n_v, dtype=np.float64)
    for i in range(n_v):
        vals, counts = np.unique(community_ids[i], return_counts=True)
        dominant_freq[i] = counts.max() / N_BOOTSTRAP

    n_consensus = int((dominant_freq > BOOTSTRAP_CONSENSUS_THRESHOLD).sum())
    frac_consensus = n_consensus / n_v if n_v > 0 else 0.0
    log(
        f"  Genes with dominant co-membership > {BOOTSTRAP_CONSENSUS_THRESHOLD}: "
        f"{n_consensus:,} / {n_v:,}  ({frac_consensus*100:.1f}%)"
    )

    df = pd.DataFrame({
        "gene": vertices,
        "dominant_community_fraction": dominant_freq,
        "n_bootstraps": N_BOOTSTRAP,
    })
    df["is_consensus"] = df["dominant_community_fraction"] > BOOTSTRAP_CONSENSUS_THRESHOLD
    return df


def plot_D(df: pd.DataFrame, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(df["dominant_community_fraction"].values, bins=40, color="C0", edgecolor="black")
    ax.axvline(BOOTSTRAP_CONSENSUS_THRESHOLD, ls="--", c="red",
               label=f"Consensus cutoff = {BOOTSTRAP_CONSENSUS_THRESHOLD}")
    n_cons = int(df["is_consensus"].sum())
    n_tot = len(df)
    ax.set_xlabel("Dominant-community co-membership fraction")
    ax.set_ylabel("Gene count")
    ax.set_title(
        f"Bootstrap community consensus (N={N_BOOTSTRAP}, {int(BOOTSTRAP_FRAC*100)}%-sample)\n"
        f"Consensus genes: {n_cons:,} / {n_tot:,} ({100*n_cons/max(n_tot,1):.1f}%)"
    )
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Analysis E: Layer-drop-one
# -----------------------------------------------------------------------------
def analysis_E_layer_drop_one(
    clinical_genes: set[str],
    node_universe: set[str],
) -> pd.DataFrame:
    log("\n=== Analysis E: Layer-drop-one ===")
    rows = []

    # Baseline (all layers)
    merged_all, _ = build_default_composite(ACTIVE_LAYERS)
    mask_all = (merged_all["p_composite"] > COMPOSITE_P_THRESHOLD).values
    fe_all, n_hits_all, n_clin = clinical_fold_enrichment(
        merged_all, clinical_genes, node_universe, mask_all
    )
    rows.append({
        "dropped_layer": "NONE_baseline",
        "n_layers_active": len(ACTIVE_LAYERS),
        "n_edges_total": len(merged_all),
        "n_composite_P_gt0p9": int(mask_all.sum()),
        "n_composite_K_ge2": int((merged_all["k_multiplicity"] >= COMPOSITE_K_THRESHOLD).sum()),
        "clinical_fold_enrichment": fe_all,
        "n_clinical_hits": n_hits_all,
        "n_clinical_genes_in_nodes": n_clin,
    })
    log(
        f"  Baseline (all {len(ACTIVE_LAYERS)} layers): "
        f"edges={len(merged_all):,}  P>0.9={int(mask_all.sum()):,}  "
        f"fold-enrich={fe_all:.3f}"
    )

    for drop in ACTIVE_LAYERS:
        layers_kept = [l for l in ACTIVE_LAYERS if l != drop]
        log(f"  Dropping layer '{drop}' -> using {layers_kept}")
        merged, _ = build_default_composite(layers_kept)
        mask = (merged["p_composite"] > COMPOSITE_P_THRESHOLD).values
        fe, n_hits, _ = clinical_fold_enrichment(
            merged, clinical_genes, node_universe, mask
        )
        rows.append({
            "dropped_layer": drop,
            "n_layers_active": len(layers_kept),
            "n_edges_total": len(merged),
            "n_composite_P_gt0p9": int(mask.sum()),
            "n_composite_K_ge2": int((merged["k_multiplicity"] >= COMPOSITE_K_THRESHOLD).sum()),
            "clinical_fold_enrichment": fe,
            "n_clinical_hits": n_hits,
            "n_clinical_genes_in_nodes": n_clin,
        })
        log(
            f"    edges={len(merged):,}  P>0.9={int(mask.sum()):,}  "
            f"fold-enrich={fe:.3f}  (delta vs baseline = {fe - fe_all:+.3f})"
        )

    df = pd.DataFrame(rows)
    df["delta_fold_enrichment_vs_baseline"] = df["clinical_fold_enrichment"] - fe_all
    return df


def plot_E(df: pd.DataFrame, out_pdf: Path) -> None:
    base = df.loc[df["dropped_layer"] == "NONE_baseline"].iloc[0]
    drop_df = df.loc[df["dropped_layer"] != "NONE_baseline"].copy()
    drop_df = drop_df.sort_values("delta_fold_enrichment_vs_baseline")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ax.barh(drop_df["dropped_layer"], drop_df["clinical_fold_enrichment"], color="C0")
    ax.axvline(base["clinical_fold_enrichment"], ls="--", c="red",
               label=f"Baseline FE = {base['clinical_fold_enrichment']:.2f}")
    ax.set_xlabel("Clinical drug target fold enrichment (P>0.9)")
    ax.set_ylabel("Dropped layer")
    ax.set_title("Layer-drop-one: fold enrichment")
    ax.legend(frameon=False)

    ax = axes[1]
    colors = ["C3" if d < 0 else "C2"
              for d in drop_df["delta_fold_enrichment_vs_baseline"].values]
    ax.barh(drop_df["dropped_layer"],
            drop_df["delta_fold_enrichment_vs_baseline"], color=colors)
    ax.axvline(0, ls=":", c="black", alpha=0.5)
    ax.set_xlabel("Delta fold enrichment vs baseline\n(negative = layer was helpful)")
    ax.set_ylabel("Dropped layer")
    ax.set_title("Layer importance (drop-one)")

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    log("=== 274: Sensitivity & Stability Analyses ===")
    log(f"Network dir: {NET_DIR}")
    log(f"Output dir:  {OUT_DIR}")
    log(f"Workers:     {N_WORKERS}")

    check_preconditions()

    # --- Load node universe + clinical gene set (shared) ------------------
    nodes = pd.read_csv(NODES_CSV, usecols=["human_symbol"]).dropna()
    node_universe = set(nodes["human_symbol"].astype(str).tolist())
    node_universe.discard("")
    log(f"Node universe: {len(node_universe):,} symbols")

    clinical_genes = load_clinical_drug_genes()
    log(f"Clinical drug target genes: {len(clinical_genes)}")

    # ===================================================================
    # Analysis A: Jaccard sweep
    # ===================================================================
    df_A = analysis_A_jaccard_sweep(clinical_genes, node_universe)
    df_A.to_csv(DATA_DIR / "A_jaccard_sweep.csv", index=False)
    plot_A(df_A, OUT_DIR / "A_jaccard_sweep.pdf")
    log(f"  Saved: A_jaccard_sweep.{{csv,pdf}}")

    # ===================================================================
    # Analysis B: Quantile sweep
    # ===================================================================
    df_B = analysis_B_quantile_sweep(clinical_genes, node_universe)
    df_B.to_csv(DATA_DIR / "B_quantile_sweep.csv", index=False)
    plot_B(df_B, OUT_DIR / "B_quantile_sweep.pdf")
    log(f"  Saved: B_quantile_sweep.{{csv,pdf}}")

    # ===================================================================
    # Build default composite (all active layers) for C/D
    # ===================================================================
    log("\n--- Building default composite graph (all 7 active layers) ---")
    merged_default, _ = build_default_composite(ACTIVE_LAYERS)
    log(
        f"  Merged composite: {len(merged_default):,} edges; "
        f"P>{LEIDEN_EDGE_P_THRESHOLD} subset = "
        f"{int((merged_default['p_composite'] > LEIDEN_EDGE_P_THRESHOLD).sum()):,}"
    )
    g_default = build_igraph_from_edges(merged_default, LEIDEN_EDGE_P_THRESHOLD)
    log(f"  Graph: {g_default.vcount():,} vertices, {g_default.ecount():,} edges")

    # ===================================================================
    # Analysis C: Leiden seed stability
    # ===================================================================
    ari_df = analysis_C_leiden_seed(g_default)
    ari_df.to_csv(DATA_DIR / "C_leiden_seed_ari.csv")
    plot_C(ari_df, OUT_DIR / "C_leiden_ari_heatmap.pdf")
    log(f"  Saved: C_leiden_seed_ari.csv + C_leiden_ari_heatmap.pdf")

    # ===================================================================
    # Analysis D: Bootstrap consensus
    # ===================================================================
    df_D = analysis_D_bootstrap(g_default)
    df_D.to_csv(DATA_DIR / "D_bootstrap_consensus.csv", index=False)
    plot_D(df_D, OUT_DIR / "D_consensus_histogram.pdf")
    log(f"  Saved: D_bootstrap_consensus.{{csv,pdf}}")

    # ===================================================================
    # Analysis E: Layer-drop-one
    # ===================================================================
    df_E = analysis_E_layer_drop_one(clinical_genes, node_universe)
    df_E.to_csv(DATA_DIR / "E_layer_drop_one.csv", index=False)
    plot_E(df_E, OUT_DIR / "E_layer_drop_one.pdf")
    log(f"  Saved: E_layer_drop_one.{{csv,pdf}}")

    log(f"\n=== 274 Complete === Total time: {time.time()-t0:.1f}s "
        f"({(time.time()-t0)/60:.1f} min)")


if __name__ == "__main__":
    main()
