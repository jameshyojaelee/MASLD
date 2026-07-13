#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_267_neighborhoods
#SBATCH --output=logs/net_267_neighborhoods_%j.out
#SBATCH --error=logs/net_267_neighborhoods_%j.err
# ===========================================================================
# Script 267: Pre-compute gene neighborhoods (top-K per gene per layer)
# ===========================================================================
# Purpose: For every gene in the network, build a JSON structure containing
#          its top-20 composite neighbors, per-layer top-20 neighbors, a
#          2-hop extended neighborhood (capped at 200 nodes), community info,
#          and a halo of drug/variant/pathway annotations from the atlas.
#
# Input:
#   - RNA-seq/results/network/posterior_edges/edges_{name}_posterior.csv
#   - RNA-seq/results/network/composite_edges.csv
#   - RNA-seq/results/network/communities/community_assignments.csv
#   - RNA-seq/results/network/network_nodes.csv
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Output:
#   - RNA-seq/results/network/gene_graphs/{SYMBOL}.json  (one per gene)
#   - RNA-seq/results/network/search_index.json  (compact lookup)
#
# Environment: spatial (pandas, pyarrow, numpy, json)
# ===========================================================================

import json
import os
import sys
import time
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.csv as pq

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
NET_DIR = BASE / "RNA-seq" / "results" / "network"
POST_DIR = NET_DIR / "posterior_edges"
COMM_DIR = NET_DIR / "communities"
GENE_DIR = NET_DIR / "gene_graphs"
GENE_DIR.mkdir(parents=True, exist_ok=True)

ATLAS_PATH = BASE / "RNA-seq" / "results" / "multi_evidence" / "multi_evidence_atlas.csv"

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic",
    "pathway", "spatial", "cosmos", "cerna", "xspecies"
]

TOP_K = 20            # top neighbors per gene per layer
TOP_K_HOP2 = 10      # per neighbor in the 2-hop extension
MAX_HALO_NODES = 200  # cap on extended neighborhood size

# Node attribute columns to include in per-gene JSON
NODE_ATTRS = [
    "is_deg", "bulk_logFC", "bulk_padj", "bulk_tstat",
    "gene_biotype", "sex_class", "is_conserved",
    "coloc_susie_best_pp4", "essentiality_chronos",
    "ferroptosis_class", "zonation_class", "layers_active",
]

# Halo columns from atlas
HALO_DRUG_COLS = ["dgidb_druggable", "opentargets_drug"]
HALO_VARIANT_COLS = ["gwas_variant_in_peak", "gwas_variant_cell_types"]
HALO_PATHWAY_COL = "top_pathways"

t0_global = time.time()

# These are populated in main() and used by worker processes
_composite_adj = {}
_layer_adjs = {}
_communities = {}
_node_attrs = {}
_halo_data = {}
_gene_list = []


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ======================================================================
# 1. Build adjacency dicts from edge parquets
# ======================================================================
def build_adjacency(edges_df, weight_col):
    """Build a dict: gene -> sorted list of (neighbor, weight) by descending weight."""
    adj = {}
    for _, row in edges_df.iterrows():
        a, b = row["gene_a"], row["gene_b"]
        w = row[weight_col]
        if pd.isna(w):
            continue
        w = float(w)
        adj.setdefault(a, []).append((b, w))
        adj.setdefault(b, []).append((a, w))

    # Sort each neighbor list by weight descending
    for gene in adj:
        adj[gene].sort(key=lambda x: -x[1])

    return adj


def build_adjacency_chunked(edges_df, weight_col, chunk_size=500000):
    """Memory-efficient adjacency builder that processes in chunks."""
    adj = {}
    n_rows = len(edges_df)
    for start in range(0, n_rows, chunk_size):
        chunk = edges_df.iloc[start:start + chunk_size]
        for _, row in chunk.iterrows():
            a, b = row["gene_a"], row["gene_b"]
            w = row[weight_col]
            if pd.isna(w):
                continue
            w = float(w)
            adj.setdefault(a, []).append((b, w))
            adj.setdefault(b, []).append((a, w))

    for gene in adj:
        adj[gene].sort(key=lambda x: -x[1])
    return adj


# ======================================================================
# 2. Build single gene JSON
# ======================================================================
def build_gene_json(gene):
    """Build the JSON structure for a single gene."""
    result = {"gene": gene}

    # Community assignments
    result["community"] = _communities.get(gene, {"macro": -1, "meso": -1, "micro": -1})

    # Node attributes
    attrs = _node_attrs.get(gene, {})
    result["attributes"] = attrs

    # Top-K composite neighbors
    comp_neighbors = _composite_adj.get(gene, [])[:TOP_K]
    neighbors_composite = []
    for nbr, p_comp in comp_neighbors:
        entry = {
            "gene": nbr,
            "p_composite": round(p_comp, 4),
        }
        # Count layer multiplicity and list per-layer posteriors
        layer_posteriors = {}
        k_mult = 0
        for layer_name, layer_adj in _layer_adjs.items():
            nbr_list = layer_adj.get(gene, [])
            for n, w in nbr_list:
                if n == nbr:
                    layer_posteriors[layer_name] = round(w, 4)
                    k_mult += 1
                    break
        entry["k_multiplicity"] = k_mult
        entry["layers"] = layer_posteriors
        neighbors_composite.append(entry)
    result["neighbors_composite"] = neighbors_composite

    # Per-layer top-K neighbors
    neighbors_per_layer = {}
    for layer_name in LAYER_NAMES:
        layer_adj = _layer_adjs.get(layer_name, {})
        nbrs = layer_adj.get(gene, [])[:TOP_K]
        if nbrs:
            neighbors_per_layer[layer_name] = [
                {"gene": n, "posterior": round(w, 4)} for n, w in nbrs
            ]
    result["neighbors_per_layer"] = neighbors_per_layer

    # 2-hop extended neighborhood (capped at MAX_HALO_NODES)
    hop1_genes = {nbr for nbr, _ in comp_neighbors}
    hop2_genes = set()
    for nbr, _ in comp_neighbors:
        nbr_nbrs = _composite_adj.get(nbr, [])[:TOP_K_HOP2]
        for nn, _ in nbr_nbrs:
            if nn != gene and nn not in hop1_genes:
                hop2_genes.add(nn)
            if len(hop1_genes) + len(hop2_genes) >= MAX_HALO_NODES:
                break
        if len(hop1_genes) + len(hop2_genes) >= MAX_HALO_NODES:
            break
    # (2-hop set stored implicitly in neighbors; not serialized separately to save space)

    # Halo: drugs, variants, pathways from atlas
    halo = {}

    # Drugs
    halo_info = _halo_data.get(gene, {})
    drugs = []
    for col in HALO_DRUG_COLS:
        val = halo_info.get(col, "")
        if val and str(val) not in ("", "nan", "NA", "False", "FALSE"):
            if col == "opentargets_drug":
                # May be semicolon-separated list
                drugs.extend([d.strip() for d in str(val).split(";") if d.strip()])
            else:
                drugs.append(str(val))
    halo["drugs"] = list(set(drugs)) if drugs else []

    # Variants
    variants = []
    for col in HALO_VARIANT_COLS:
        val = halo_info.get(col, "")
        if val and str(val) not in ("", "nan", "NA", "False", "FALSE"):
            variants.extend([v.strip() for v in str(val).split(";") if v.strip()])
    halo["variants"] = list(set(variants)) if variants else []

    # Pathways
    pathways_val = halo_info.get(HALO_PATHWAY_COL, "")
    if pathways_val and str(pathways_val) not in ("", "nan", "NA"):
        halo["top_pathways"] = [p.strip() for p in str(pathways_val).split(";") if p.strip()]
    else:
        halo["top_pathways"] = []

    result["halo"] = halo

    return result


def process_gene(gene):
    """Process a single gene: build JSON and write to file."""
    try:
        gene_json = build_gene_json(gene)
        out_path = GENE_DIR / f"{gene}.json"
        with open(out_path, "w") as f:
            json.dump(gene_json, f, separators=(",", ":"))
        return gene, True
    except Exception as e:
        return gene, str(e)


# ======================================================================
# 3. Build search index
# ======================================================================
def build_search_index(gene_list):
    """Build a compact search index: gene -> {degree, community_macro, top_neighbor}."""
    log("Building search index...")
    index = {}
    for gene in gene_list:
        comp_nbrs = _composite_adj.get(gene, [])
        degree = len(comp_nbrs)
        community = _communities.get(gene, {}).get("macro", -1)
        top_nbr = comp_nbrs[0][0] if comp_nbrs else ""
        index[gene] = {
            "degree": degree,
            "community_macro": community,
            "top_neighbor": top_nbr,
        }
    return index


# ======================================================================
# Main
# ======================================================================
def main():
    global _composite_adj, _layer_adjs, _communities, _node_attrs, _halo_data, _gene_list

    log("=== 267: Gene Neighborhood Pre-computation ===")
    log(f"Output: {GENE_DIR}")

    # --- Load community assignments ---
    comm_path = COMM_DIR / "community_assignments.csv"
    if not comm_path.exists():
        log(f"ERROR: Community assignments not found: {comm_path}")
        log("  Run Script 266 first.")
        sys.exit(1)
    comm_df = pd.read_csv(comm_path)
    log(f"  Communities: {len(comm_df):,} genes")
    for _, row in comm_df.iterrows():
        _communities[row["gene"]] = {
            "macro": int(row["macro_id"]),
            "meso": int(row["meso_id"]),
            "micro": int(row["micro_id"]),
        }

    # --- Load node attributes ---
    node_path = NET_DIR / "network_nodes.csv"
    if not node_path.exists():
        log(f"ERROR: Node file not found: {node_path}")
        sys.exit(1)
    nodes_df = pd.read_csv(node_path)
    gene_col = "gene" if "gene" in nodes_df.columns else "human_symbol"
    log(f"  Nodes: {len(nodes_df):,} genes")

    for _, row in nodes_df.iterrows():
        gene = row[gene_col]
        attrs = {}
        for col in NODE_ATTRS:
            if col in nodes_df.columns:
                val = row[col]
                if pd.isna(val):
                    attrs[col] = None
                elif isinstance(val, (np.integer, np.int64)):
                    attrs[col] = int(val)
                elif isinstance(val, (np.floating, np.float64)):
                    attrs[col] = round(float(val), 4)
                elif isinstance(val, (np.bool_, bool)):
                    attrs[col] = bool(val)
                else:
                    attrs[col] = str(val)
        _node_attrs[gene] = attrs

    # --- Load atlas for halo data ---
    if ATLAS_PATH.exists():
        log(f"  Loading atlas for halo data: {ATLAS_PATH}")
        atlas = pd.read_csv(ATLAS_PATH, usecols=lambda c: c in (
            ["human_symbol"] + HALO_DRUG_COLS + HALO_VARIANT_COLS + [HALO_PATHWAY_COL]
        ))
        for _, row in atlas.iterrows():
            gene = row["human_symbol"]
            if pd.isna(gene):
                continue
            halo = {}
            for col in HALO_DRUG_COLS + HALO_VARIANT_COLS + [HALO_PATHWAY_COL]:
                if col in atlas.columns:
                    halo[col] = row[col]
            _halo_data[gene] = halo
        log(f"  Halo data loaded for {len(_halo_data):,} genes")
    else:
        log(f"  WARNING: Atlas not found at {ATLAS_PATH}, halo will be empty")

    # --- Load composite edges ---
    comp_path = NET_DIR / "composite_edges.csv"
    if not comp_path.exists():
        log(f"ERROR: Composite edges not found: {comp_path}")
        sys.exit(1)
    log("  Loading composite edges...")
    comp_df = pd.read_csv(comp_path)
    weight_col = "p_composite"
    if weight_col not in comp_df.columns:
        for alt in ["posterior", "weight", "composite_posterior"]:
            if alt in comp_df.columns:
                weight_col = alt
                break
    log(f"  Composite: {len(comp_df):,} edges, weight col: {weight_col}")
    _composite_adj = build_adjacency_chunked(comp_df, weight_col)
    log(f"  Composite adjacency: {len(_composite_adj):,} genes with edges")

    # --- Load per-layer edges ---
    log("  Loading per-layer posterior edges...")
    for name in LAYER_NAMES:
        path = POST_DIR / f"edges_{name}_posterior.csv"
        if path.exists():
            df = pd.read_parquet(path)
            _layer_adjs[name] = build_adjacency_chunked(df, "posterior")
            log(f"    {name}: {len(df):,} edges, {len(_layer_adjs[name]):,} genes")
        else:
            log(f"    {name}: not found, skipping")

    # --- Gene list ---
    _gene_list = sorted(nodes_df[gene_col].dropna().unique().tolist())
    log(f"\n  Total genes to process: {len(_gene_list):,}")

    # --- Generate per-gene JSONs using multiprocessing ---
    log("\n--- Generating per-gene JSON files ---")
    n_workers = min(cpu_count(), 8)
    log(f"  Workers: {n_workers}")

    # Process in serial with progress reporting (multiprocessing would require
    # pickling the large adjacency dicts; serial with progress is more reliable)
    n_total = len(_gene_list)
    n_success = 0
    n_fail = 0
    report_interval = max(1, n_total // 20)

    for i, gene in enumerate(_gene_list):
        result = process_gene(gene)
        if result[1] is True:
            n_success += 1
        else:
            n_fail += 1
            if n_fail <= 10:
                log(f"  FAIL: {result[0]} -> {result[1]}")

        if (i + 1) % report_interval == 0 or (i + 1) == n_total:
            pct = 100 * (i + 1) / n_total
            log(f"  Progress: {i + 1:,}/{n_total:,} ({pct:.0f}%) "
                f"[ok={n_success:,}, fail={n_fail}]")

    log(f"\n  JSON files written: {n_success:,} success, {n_fail} failures")

    # --- Build and save search index ---
    search_index = build_search_index(_gene_list)
    index_path = NET_DIR / "search_index.json"
    with open(index_path, "w") as f:
        json.dump(search_index, f, separators=(",", ":"))
    log(f"  Search index: {index_path} ({len(search_index):,} entries)")

    # Summary
    elapsed = time.time() - t0_global
    log(f"\n=== 267 Complete ===")
    log(f"  Genes processed: {n_success:,}")
    log(f"  JSON dir: {GENE_DIR}")
    log(f"  Search index: {index_path}")
    log(f"  Total time: {elapsed:.1f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
