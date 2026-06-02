#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_266_leiden
#SBATCH --output=logs/net_266_leiden_%j.out
#SBATCH --error=logs/net_266_leiden_%j.err
# ===========================================================================
# Script 266: Hierarchical Leiden community detection at 3 resolutions
# ===========================================================================
# Purpose: Run Leiden algorithm on the composite graph at macro (0.5), meso
#          (1.0), and micro (2.0) resolutions. Enforce hierarchical nesting
#          by majority vote so each micro cluster is fully contained within
#          one meso module, and each meso within one macro module.
#
# Input:
#   - RNA-seq/results/network/graphs/composite_graph.graphml
#
# Output:
#   - RNA-seq/results/network/communities/community_assignments.csv
#     (gene, macro_id, meso_id, micro_id)
#   - RNA-seq/results/network/communities/community_sizes.csv
#     (level, id, n_genes, modularity)
#
# Environment: spatial (igraph, leidenalg, pandas, numpy)
# ===========================================================================

import os
import sys
import time
import warnings
from pathlib import Path
from collections import Counter

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
GRAPH_DIR = NET_DIR / "graphs"
COMM_DIR = NET_DIR / "communities"
COMM_DIR.mkdir(parents=True, exist_ok=True)

# Leiden resolution parameters
RESOLUTIONS = {
    "macro": 0.5,
    "meso": 1.0,
    "micro": 2.0,
}

SEED = 42

t0_global = time.time()


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ======================================================================
# 1. Load composite graph
# ======================================================================
def load_composite_graph():
    """Load the composite graph from GraphML."""
    path = GRAPH_DIR / "composite_graph.graphml"
    if not path.exists():
        log(f"ERROR: Composite graph not found: {path}")
        log("  Run Script 265 first.")
        sys.exit(1)
    log(f"Loading composite graph: {path}")
    g = ig.Graph.Read_GraphML(str(path))
    log(f"  Vertices: {g.vcount():,}, Edges: {g.ecount():,}")

    # Verify weights exist
    if "weight" not in g.es.attributes():
        log("ERROR: Graph has no 'weight' edge attribute")
        sys.exit(1)

    # Ensure weights are positive (Leiden requires this)
    weights = np.array(g.es["weight"], dtype=np.float64)
    n_neg = np.sum(weights < 0)
    if n_neg > 0:
        log(f"  WARNING: {n_neg} negative weights found, clipping to 0")
        weights = np.clip(weights, 0, None)
        g.es["weight"] = weights.tolist()

    return g


# ======================================================================
# 2. Run Leiden at a single resolution
# ======================================================================
def run_leiden(g, resolution, level_name):
    """Run Leiden community detection using leidenalg.

    Uses RBConfigurationVertexPartition (modularity with resolution parameter)
    which is the standard choice for weighted networks.
    """
    log(f"\n--- Leiden at resolution {resolution} ({level_name}) ---")
    t_start = time.time()

    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=SEED,
        n_iterations=-1,  # iterate until convergence
    )

    n_communities = len(partition)
    quality = partition.quality()

    # Also compute standard modularity for comparison
    modularity = g.modularity(partition.membership, weights="weight")

    membership = partition.membership
    sizes = [len(part) for part in partition]

    elapsed = time.time() - t_start
    log(f"  Communities: {n_communities}")
    log(f"  Quality (Leiden): {quality:.4f}")
    log(f"  Modularity (standard): {modularity:.4f}")
    log(f"  Size range: {min(sizes)} - {max(sizes)} (median: {np.median(sizes):.0f})")
    log(f"  Singletons: {sum(1 for s in sizes if s == 1)}")
    log(f"  Time: {elapsed:.1f}s")

    return membership, modularity, quality


# ======================================================================
# 3. Force hierarchical nesting by majority vote
# ======================================================================
def enforce_nesting(fine_membership, coarse_membership, fine_name, coarse_name):
    """Ensure each fine community maps to exactly one coarse community.

    For each fine community, assign it to the coarse community that holds
    the majority of its members. If a fine community is split across
    multiple coarse communities, the minority members are reassigned to
    the majority coarse community for that fine cluster.

    Parameters
    ----------
    fine_membership : list
        Per-vertex fine-level community IDs.
    coarse_membership : list
        Per-vertex coarse-level community IDs.

    Returns
    -------
    list : corrected coarse membership (same length as input).
    int : number of vertices reassigned.
    """
    log(f"  Enforcing {fine_name} -> {coarse_name} nesting...")

    n = len(fine_membership)
    corrected_coarse = list(coarse_membership)
    n_reassigned = 0

    # Group vertices by fine community
    fine_communities = {}
    for i in range(n):
        fc = fine_membership[i]
        if fc not in fine_communities:
            fine_communities[fc] = []
        fine_communities[fc].append(i)

    # For each fine community, find the majority coarse community
    for fc, members in fine_communities.items():
        coarse_ids = [coarse_membership[m] for m in members]
        majority_coarse = Counter(coarse_ids).most_common(1)[0][0]

        # Reassign any members not in the majority coarse community
        for m in members:
            if coarse_membership[m] != majority_coarse:
                corrected_coarse[m] = majority_coarse
                n_reassigned += 1

    log(f"    Vertices reassigned: {n_reassigned} / {n} "
        f"({100 * n_reassigned / n:.2f}%)")

    return corrected_coarse, n_reassigned


# ======================================================================
# 4. Relabel communities for clean sequential IDs
# ======================================================================
def relabel_sequential(membership):
    """Relabel community IDs to be sequential 0, 1, 2, ... by decreasing size."""
    counter = Counter(membership)
    # Sort by size descending, then by original ID for ties
    sorted_ids = sorted(counter.keys(), key=lambda x: (-counter[x], x))
    old_to_new = {old: new for new, old in enumerate(sorted_ids)}
    return [old_to_new[m] for m in membership]


# ======================================================================
# Main
# ======================================================================
def main():
    log("=== 266: Hierarchical Leiden Community Detection ===")
    log(f"Output directory: {COMM_DIR}")

    g = load_composite_graph()
    gene_names = g.vs["name"]
    n_genes = len(gene_names)

    # Run Leiden at all three resolutions
    memberships = {}
    modularities = {}
    qualities = {}

    for level, res in RESOLUTIONS.items():
        mem, mod, qual = run_leiden(g, res, level)
        memberships[level] = mem
        modularities[level] = mod
        qualities[level] = qual

    # Enforce hierarchical nesting (micro -> meso -> macro)
    # Step 1: Ensure each micro cluster maps to one meso module
    log("\n--- Enforcing hierarchical nesting ---")
    meso_corrected, n_meso_fix = enforce_nesting(
        memberships["micro"], memberships["meso"], "micro", "meso"
    )
    memberships["meso"] = meso_corrected

    # Step 2: Ensure each meso module maps to one macro module
    macro_corrected, n_macro_fix = enforce_nesting(
        memberships["meso"], memberships["macro"], "meso", "macro"
    )
    memberships["macro"] = macro_corrected

    # Step 3: After macro correction, re-verify micro->meso nesting
    # (macro correction may have shifted meso assignments)
    meso_recheck, n_recheck = enforce_nesting(
        memberships["micro"], memberships["meso"], "micro", "meso (recheck)"
    )
    memberships["meso"] = meso_recheck

    # Relabel all levels to sequential IDs (largest community = 0)
    log("\n--- Relabeling communities ---")
    for level in ["macro", "meso", "micro"]:
        memberships[level] = relabel_sequential(memberships[level])
        n_comm = len(set(memberships[level]))
        log(f"  {level}: {n_comm} communities")

    # Recompute modularity after nesting corrections
    log("\n--- Post-correction modularity ---")
    for level in ["macro", "meso", "micro"]:
        mod = g.modularity(memberships[level], weights="weight")
        modularities[level] = mod
        log(f"  {level} (res={RESOLUTIONS[level]}): Q = {mod:.4f}")

    # --- Save community assignments ---
    log("\n--- Saving outputs ---")
    assign_df = pd.DataFrame({
        "gene": gene_names,
        "macro_id": memberships["macro"],
        "meso_id": memberships["meso"],
        "micro_id": memberships["micro"],
    })
    assign_out = COMM_DIR / "community_assignments.csv"
    assign_df.to_csv(assign_out, index=False)
    log(f"  Assignments: {assign_out} ({len(assign_df):,} genes)")

    # --- Save community sizes ---
    size_rows = []
    for level in ["macro", "meso", "micro"]:
        counter = Counter(memberships[level])
        for cid, count in sorted(counter.items()):
            size_rows.append({
                "level": level,
                "id": cid,
                "n_genes": count,
                "modularity": modularities[level],
                "resolution": RESOLUTIONS[level],
            })
    size_df = pd.DataFrame(size_rows)
    size_out = COMM_DIR / "community_sizes.csv"
    size_df.to_csv(size_out, index=False)
    log(f"  Sizes: {size_out} ({len(size_df)} rows)")

    # --- Nesting verification ---
    log("\n--- Nesting verification ---")
    for fine, coarse in [("micro", "meso"), ("meso", "macro")]:
        fine_to_coarse = {}
        violations = 0
        for i in range(n_genes):
            fc = memberships[fine][i]
            cc = memberships[coarse][i]
            if fc in fine_to_coarse:
                if fine_to_coarse[fc] != cc:
                    violations += 1
            else:
                fine_to_coarse[fc] = cc
        log(f"  {fine} -> {coarse}: {len(fine_to_coarse)} mappings, {violations} violations")

    # Summary
    elapsed = time.time() - t0_global
    log(f"\n=== 266 Complete ===")
    log(f"  Genes: {n_genes:,}")
    log(f"  Macro: {len(set(memberships['macro']))} communities (Q={modularities['macro']:.4f})")
    log(f"  Meso: {len(set(memberships['meso']))} communities (Q={modularities['meso']:.4f})")
    log(f"  Micro: {len(set(memberships['micro']))} communities (Q={modularities['micro']:.4f})")
    log(f"  Nesting fixes: meso={n_meso_fix + n_recheck}, macro={n_macro_fix}")
    log(f"  Total time: {elapsed:.1f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
