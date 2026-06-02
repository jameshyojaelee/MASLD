#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_265_assemble
#SBATCH --output=logs/net_265_assemble_%j.out
#SBATCH --error=logs/net_265_assemble_%j.err
# ===========================================================================
# Script 265: Assemble the multiplex graph from posterior-weighted edges
# ===========================================================================
# Purpose: Build igraph Graph objects from posterior edges (one per layer +
#          one composite), attach node attributes from the atlas, compute
#          per-layer graph statistics, and generate a ForceAtlas2-style layout
#          on the composite graph for the Community Map portal view.
#
# Input:
#   - RNA-seq/results/network/posterior_edges/edges_{name}_posterior.csv
#   - RNA-seq/results/network/composite_edges.csv
#   - RNA-seq/results/network/network_nodes.csv
#
# Output:
#   - RNA-seq/results/network/graphs/composite_graph.graphml
#   - RNA-seq/results/network/graphs/layer_{name}.graphml  (10 files)
#   - RNA-seq/results/network/global_layout.csv  (gene, x, y from FR layout)
#   - RNA-seq/results/network/graph_stats.csv
#
# Environment: spatial (igraph, pandas, pyarrow, numpy)
# ===========================================================================

import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.csv as pq
import igraph as ig

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
NET_DIR = BASE / "RNA-seq" / "results" / "network"
POST_DIR = NET_DIR / "posterior_edges"
GRAPH_DIR = NET_DIR / "graphs"
GRAPH_DIR.mkdir(parents=True, exist_ok=True)

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic",
    "pathway", "spatial", "cosmos", "cerna", "xspecies"
]

# Minimum posterior threshold to include an edge in graphs
EDGE_POSTERIOR_MIN = 0.5

t0_global = time.time()


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ======================================================================
# 1. Load node set and all edge data
# ======================================================================
def load_inputs():
    """Load node attributes, posterior edge parquets, and composite edges."""
    log("=== Loading inputs ===")

    # Node set
    node_path = NET_DIR / "network_nodes.csv"
    if not node_path.exists():
        log(f"ERROR: Node file not found: {node_path}")
        sys.exit(1)
    nodes = pd.read_csv(node_path)
    log(f"  Nodes: {len(nodes):,} genes, columns: {list(nodes.columns)}")

    # Posterior edges per layer
    layer_edges = {}
    for name in LAYER_NAMES:
        path = POST_DIR / f"edges_{name}_posterior.csv"
        if path.exists():
            df = pd.read_parquet(path)
            n_raw = len(df)
            df = df[df["posterior"] >= EDGE_POSTERIOR_MIN].copy()
            layer_edges[name] = df
            log(f"  Layer {name}: {n_raw:,} raw -> {len(df):,} edges (P >= {EDGE_POSTERIOR_MIN})")
        else:
            log(f"  Layer {name}: parquet not found, skipping")

    if not layer_edges:
        log("ERROR: No posterior edge files found. Run Script 262 first.")
        sys.exit(1)

    # Composite edges
    comp_path = NET_DIR / "composite_edges.csv"
    if not comp_path.exists():
        log(f"ERROR: Composite edges not found: {comp_path}")
        sys.exit(1)
    composite = pd.read_csv(comp_path)
    log(f"  Composite edges: {len(composite):,}")

    return nodes, layer_edges, composite


# ======================================================================
# 2. Build igraph Graph from edge DataFrame
# ======================================================================
def build_graph(vertex_names, edges_df, weight_col, graph_name=""):
    """Build an igraph Graph with a fixed vertex set and weighted edges.

    Parameters
    ----------
    vertex_names : list
        Canonical ordered gene list (all vertices present in every graph).
    edges_df : pd.DataFrame
        Must contain gene_a, gene_b, and the column named by weight_col.
    weight_col : str
        Column name for edge weights (e.g. 'posterior' or 'p_composite').
    graph_name : str
        Name attribute for the graph.

    Returns
    -------
    ig.Graph
    """
    name_to_idx = {g: i for i, g in enumerate(vertex_names)}

    # Filter to edges where both endpoints are in the vertex set
    mask_a = edges_df["gene_a"].isin(name_to_idx)
    mask_b = edges_df["gene_b"].isin(name_to_idx)
    valid = edges_df[mask_a & mask_b].copy()

    edge_tuples = list(zip(
        valid["gene_a"].map(name_to_idx),
        valid["gene_b"].map(name_to_idx)
    ))
    weights = valid[weight_col].values.tolist()

    g = ig.Graph(n=len(vertex_names), edges=edge_tuples, directed=False)
    g.vs["name"] = vertex_names
    g.es["weight"] = weights
    g["name"] = graph_name

    # Remove multi-edges, keeping max weight
    g = g.simplify(combine_edges={"weight": "max"})

    return g


# ======================================================================
# 3. Attach node attributes from network_nodes.csv
# ======================================================================
def attach_node_attributes(g, nodes_df):
    """Copy columns from nodes DataFrame to graph vertex attributes."""
    name_to_idx = {v["name"]: v.index for v in g.vs}
    # Determine which columns to attach (skip 'gene' which is the key)
    attr_cols = [c for c in nodes_df.columns if c != "gene"]

    # Initialize with NA / defaults
    for col in attr_cols:
        g.vs[col] = [None] * g.vcount()

    # Map node rows to vertex indices
    for _, row in nodes_df.iterrows():
        gene = row.get("gene") or row.get("human_symbol")
        if gene in name_to_idx:
            idx = name_to_idx[gene]
            for col in attr_cols:
                val = row[col]
                # igraph GraphML cannot store NaN — convert to empty string
                if pd.isna(val):
                    g.vs[idx][col] = ""
                else:
                    g.vs[idx][col] = str(val) if not isinstance(val, (int, float, bool)) else val

    return g


# ======================================================================
# 4. Compute per-layer graph statistics
# ======================================================================
def compute_graph_stats(g, layer_name):
    """Compute summary statistics for a graph."""
    n_v = g.vcount()
    n_e = g.ecount()
    max_edges = n_v * (n_v - 1) / 2 if n_v > 1 else 1
    density = n_e / max_edges if max_edges > 0 else 0.0

    degrees = g.degree()
    mean_deg = np.mean(degrees) if degrees else 0.0
    median_deg = np.median(degrees) if degrees else 0.0
    max_deg = max(degrees) if degrees else 0
    n_isolates = sum(1 for d in degrees if d == 0)

    weights = g.es["weight"] if g.ecount() > 0 else []
    mean_weight = np.mean(weights) if weights else 0.0
    median_weight = np.median(weights) if weights else 0.0

    return {
        "layer": layer_name,
        "n_vertices": n_v,
        "n_edges": n_e,
        "density": density,
        "mean_degree": mean_deg,
        "median_degree": median_deg,
        "max_degree": max_deg,
        "n_isolates": n_isolates,
        "mean_weight": mean_weight,
        "median_weight": median_weight,
    }


# ======================================================================
# 5. Compute layout on composite graph
# ======================================================================
def compute_layout(g):
    """Compute Fruchterman-Reingold layout using composite edge weights.

    ForceAtlas2 is not available in this igraph build, so we use
    Fruchterman-Reingold with weights as an excellent alternative for
    force-directed community-preserving layouts.
    """
    log("  Computing Fruchterman-Reingold layout (weighted)...")

    # FR layout: weights act as spring strengths (higher weight = closer)
    layout = g.layout_fruchterman_reingold(
        weights="weight",
        niter=1000,
        seed=None,  # deterministic via grid init
    )
    return layout


# ======================================================================
# Main
# ======================================================================
def main():
    log("=== 265: Assemble Multiplex Graph ===")
    log(f"Network directory: {NET_DIR}")
    log(f"Graph output: {GRAPH_DIR}")

    nodes, layer_edges, composite = load_inputs()

    # Canonical vertex list — use gene column from nodes
    gene_col = "gene" if "gene" in nodes.columns else "human_symbol"
    vertex_names = sorted(nodes[gene_col].dropna().unique().tolist())
    log(f"\nCanonical vertex set: {len(vertex_names):,} genes")

    # --- Build composite graph ---
    log("\n--- Building composite graph ---")
    weight_col = "p_composite"
    if weight_col not in composite.columns:
        # Fallback: try 'posterior' or 'weight'
        for alt in ["posterior", "weight", "composite_posterior"]:
            if alt in composite.columns:
                weight_col = alt
                break
        else:
            log(f"ERROR: No weight column found in composite edges. Columns: {list(composite.columns)}")
            sys.exit(1)
    log(f"  Using weight column: {weight_col}")

    composite_graph = build_graph(vertex_names, composite, weight_col, graph_name="composite")
    attach_node_attributes(composite_graph, nodes)
    log(f"  Composite graph: {composite_graph.vcount()} vertices, {composite_graph.ecount()} edges")

    # Also store k_multiplicity as edge attribute if available
    if "k_multiplicity" in composite.columns:
        name_to_idx = {g: i for i, g in enumerate(vertex_names)}
        # Build a lookup for multiplicity
        mult_lookup = {}
        for _, row in composite.iterrows():
            a, b = row["gene_a"], row["gene_b"]
            if a in name_to_idx and b in name_to_idx:
                key = (min(name_to_idx[a], name_to_idx[b]), max(name_to_idx[a], name_to_idx[b]))
                mult_lookup[key] = int(row["k_multiplicity"])
        # Assign to edges
        k_vals = []
        for e in composite_graph.es:
            key = (min(e.source, e.target), max(e.source, e.target))
            k_vals.append(mult_lookup.get(key, 1))
        composite_graph.es["k_multiplicity"] = k_vals

    # Save composite graph
    comp_out = GRAPH_DIR / "composite_graph.graphml"
    composite_graph.write_graphml(str(comp_out))
    log(f"  Saved: {comp_out}")

    # --- Build per-layer graphs ---
    log("\n--- Building per-layer graphs ---")
    stats_rows = []

    for name in LAYER_NAMES:
        if name not in layer_edges:
            log(f"  Layer {name}: no edges, skipping graph build")
            continue

        df = layer_edges[name]
        g = build_graph(vertex_names, df, "posterior", graph_name=name)
        attach_node_attributes(g, nodes)

        layer_out = GRAPH_DIR / f"layer_{name}.graphml"
        g.write_graphml(str(layer_out))
        log(f"  Layer {name}: {g.ecount():,} edges -> {layer_out}")

        stats = compute_graph_stats(g, name)
        stats_rows.append(stats)

    # Composite stats
    comp_stats = compute_graph_stats(composite_graph, "composite")
    stats_rows.append(comp_stats)

    # Save graph stats
    stats_df = pd.DataFrame(stats_rows)
    stats_out = NET_DIR / "graph_stats.csv"
    stats_df.to_parquet(stats_out, index=False)
    log(f"\n  Graph stats saved: {stats_out}")
    log("\n" + stats_df.to_string(index=False))

    # --- Compute global layout ---
    log("\n--- Computing global layout on composite graph ---")
    layout = compute_layout(composite_graph)
    coords = np.array(layout.coords)
    layout_df = pd.DataFrame({
        "gene": vertex_names,
        "x": coords[:, 0],
        "y": coords[:, 1],
    })
    layout_out = NET_DIR / "global_layout.csv"
    layout_df.to_parquet(layout_out, index=False)
    log(f"  Layout saved: {layout_out} ({len(layout_df):,} genes)")

    # Summary
    elapsed = time.time() - t0_global
    log(f"\n=== 265 Complete ===")
    log(f"  Composite: {composite_graph.vcount()} vertices, {composite_graph.ecount()} edges")
    log(f"  Layers built: {len([s for s in stats_rows if s['layer'] != 'composite'])}")
    log(f"  Layout: {len(layout_df)} gene positions")
    log(f"  Total time: {elapsed:.1f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
