#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --job-name=net_297_mod_null
#SBATCH --output=logs/net_297_mod_null_%j.out
#SBATCH --error=logs/net_297_mod_null_%j.err
# ===========================================================================
# Script 297: Maslov-Sneppen degree-preserving modularity null for the
#             F2-weighted (binary) and continuous-stage Leiden partitions.
# ===========================================================================
# Purpose:
#   Phase 1b Team-2 C2 / B10.1 audit (RP10 deliberation 2026-05-12)
#   flagged that the headline Q = 0.69 for the F2-weighted Leiden in
#   Script 292 has NO degree-preserving null distribution on disk to
#   support the "lift above chance" claim. Script 260 has
#   `permute_maslov_sneppen()` defined but it is never invoked by the
#   dispatcher (line 267 routes "ppi" to score-shuffle, not edge-rewire).
#
#   This script generates the missing null:
#     1. Load the EXACT graph used by 292 (and optionally 293).
#     2. Run Maslov-Sneppen rewiring N_REWIRES times (default 200).
#        Each rewired graph preserves the degree sequence of every node
#        but randomizes the topology.
#     3. Re-run Leiden on each rewired graph with identical parameters.
#     4. Record Q for each rewired graph.
#     5. Compute empirical p-value: P(Q_null >= Q_observed).
#
#   Output also reports the 95% null Q-distribution band so manuscript
#   can cite "observed Q = 0.69 vs null 95% upper bound = X.XX".
#
# Inputs:
#   RNA-seq/results/network/edges_ppi.csv
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (f2_inflection_logFC)
#
# Outputs (RNA-seq/results/network/q_null/):
#   q_null_distribution.csv   (one row per rewire: Q_rewired)
#   q_null_summary.json       (observed Q, mean, 95% bounds, p-value)
#
# Environment: spatial (igraph, leidenalg, pandas, numpy)
# Run via: sbatch RNA-seq/297_modularity_null.py
# ===========================================================================

import os
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import igraph as ig
import leidenalg

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
NET_DIR = BASE / "RNA-seq" / "results" / "network"
ATLAS = BASE / "RNA-seq" / "results" / "multi_evidence" / "multi_evidence_atlas.csv"
OUT_DIR = NET_DIR / "q_null"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STRING_THRESH = 0.7
LEIDEN_RES = 1.0
SEED_BASE = 42
N_REWIRES = int(os.environ.get("N_REWIRES", "200"))
REWIRE_PASSES = int(os.environ.get("REWIRE_PASSES", "10"))  # Maslov-Sneppen swap count = REWIRE_PASSES * |E|

t0 = time.time()


def log(msg: str) -> None:
    print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)


def f2_weights_from_atlas() -> tuple[dict[str, float], dict[str, float]]:
    """Reproduce Script 292's W_F0F1 and W_F3F4 from atlas."""
    df = pd.read_csv(ATLAS, low_memory=False,
                     usecols=["human_symbol", "f2_inflection_logFC"])
    df = df.dropna(subset=["f2_inflection_logFC"])
    df["w_f01"] = (-df["f2_inflection_logFC"]).clip(lower=0.0)
    df["w_f34"] = (df["f2_inflection_logFC"]).clip(lower=0.0)
    w_f01 = df.set_index("human_symbol")["w_f01"].to_dict()
    w_f34 = df.set_index("human_symbol")["w_f34"].to_dict()
    log(f"  F0-F1 weighting: {len(w_f01):,} genes")
    log(f"  F3-F4 weighting: {len(w_f34):,} genes")
    return w_f01, w_f34


def build_graph(edges_df: pd.DataFrame, w_gene: dict[str, float]) -> ig.Graph:
    """Same construction as Script 292/293."""
    genes = pd.unique(edges_df[["gene_a", "gene_b"]].values.ravel("K"))
    g = ig.Graph()
    g.add_vertices(list(genes))
    name_to_idx = {n: i for i, n in enumerate(g.vs["name"])}

    wa = np.array([w_gene.get(a, 0.0) for a in edges_df["gene_a"]])
    wb = np.array([w_gene.get(b, 0.0) for b in edges_df["gene_b"]])
    mult = np.exp(0.5 * (wa + wb))
    weights = edges_df["raw_score"].to_numpy() * mult

    src = [name_to_idx[a] for a in edges_df["gene_a"]]
    dst = [name_to_idx[b] for b in edges_df["gene_b"]]
    g.add_edges(list(zip(src, dst)))
    g.es["weight"] = weights.tolist()
    g.simplify(combine_edges={"weight": "max"})
    return g


def run_leiden(g: ig.Graph, seed: int) -> float:
    part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=LEIDEN_RES,
        seed=seed,
        n_iterations=-1,
    )
    return g.modularity(part.membership, weights="weight")


def maslov_sneppen_rewire(g: ig.Graph, n_passes: int, rng: np.random.Generator) -> ig.Graph:
    """Perform Maslov-Sneppen rewiring on a copy of g.

    igraph's `rewire(mode='simple')` is a C-speed implementation of
    Maslov-Sneppen degree-preserving rewiring. It SWAPS endpoints across
    edge pairs in place, preserving degree sequence exactly.

    Edge-weight handling: `Graph.rewire()` invalidates the per-edge
    `weight` attribute by replacing edges. We follow the standard
    weighted-modularity-null practice (Maslov-Sneppen Option A): after
    rewiring, attach the ORIGINAL weight distribution to the rewired
    edges in randomized order. This tests "is observed Q higher than
    what you would get from a random topology with the same weight
    distribution?" — the canonical null for community-detection lift.
    """
    g_re = g.copy()
    original_weights = np.asarray(g_re.es["weight"], dtype=float)
    n_iter = n_passes * g_re.ecount()
    g_re.rewire(n=n_iter, mode="simple")
    # Reattach weights in randomized order (Option A: preserve weight
    # distribution, randomize weight->edge assignment).
    perm = rng.permutation(len(original_weights))
    g_re.es["weight"] = original_weights[perm].tolist()
    return g_re


def main() -> None:
    log("=== 297: Maslov-Sneppen Q-null for F2-weighted Leiden ===")
    log(f"Out: {OUT_DIR}")
    log(f"N_REWIRES={N_REWIRES}, REWIRE_PASSES={REWIRE_PASSES}")

    edges_ppi = pd.read_csv(NET_DIR / "edges_ppi.csv")
    edges_high = edges_ppi[edges_ppi["raw_score"] >= STRING_THRESH].copy()
    log(f"PPI edges >= {STRING_THRESH}: {len(edges_high):,}")

    # ---- two F2 weighting variants (match 292) ----
    w_f01, w_f34 = f2_weights_from_atlas()

    for label, w_gene in [("F01_binary", w_f01), ("F34_binary", w_f34)]:
        log(f"---- {label} ----")
        g = build_graph(edges_high, w_gene)
        log(f"  Observed graph: {g.vcount():,} nodes, {g.ecount():,} edges")

        # Observed Q
        Q_obs = run_leiden(g, SEED_BASE)
        log(f"  Observed Q ({label}) = {Q_obs:.4f}")

        # Null distribution
        rng = np.random.default_rng(SEED_BASE)
        Q_null = []
        for i in range(N_REWIRES):
            t_pass = time.time()
            g_re = maslov_sneppen_rewire(g, REWIRE_PASSES, rng)
            Q = run_leiden(g_re, SEED_BASE + 1 + i)
            Q_null.append(Q)
            if (i + 1) % 10 == 0:
                log(f"  rewire {i+1}/{N_REWIRES}: Q={Q:.4f} "
                    f"(dt={time.time()-t_pass:.1f}s)")

        Q_null = np.asarray(Q_null)
        n_extreme = int((Q_null >= Q_obs).sum())
        empirical_p = (n_extreme + 1) / (N_REWIRES + 1)
        summary = {
            "label": label,
            "Q_observed": float(Q_obs),
            "Q_null_mean": float(Q_null.mean()),
            "Q_null_sd": float(Q_null.std(ddof=1)),
            "Q_null_q025": float(np.quantile(Q_null, 0.025)),
            "Q_null_q975": float(np.quantile(Q_null, 0.975)),
            "Q_null_max": float(Q_null.max()),
            "N_REWIRES": N_REWIRES,
            "REWIRE_PASSES": REWIRE_PASSES,
            "empirical_p_one_sided": empirical_p,
        }
        log(f"  Q_null mean={Q_null.mean():.4f}, sd={Q_null.std(ddof=1):.4f}, "
            f"95% upper={np.quantile(Q_null, 0.975):.4f}, max={Q_null.max():.4f}, "
            f"p={empirical_p:.4f}")

        # Save
        out_csv = OUT_DIR / f"q_null_distribution_{label}.csv"
        pd.DataFrame({"rewire": np.arange(N_REWIRES) + 1,
                      "Q_rewired": Q_null}).to_csv(out_csv, index=False)
        log(f"  Wrote {out_csv}")
        out_json = OUT_DIR / f"q_null_summary_{label}.json"
        out_json.write_text(json.dumps(summary, indent=2))
        log(f"  Wrote {out_json}")

    log("=== 297 complete ===")


if __name__ == "__main__":
    main()
