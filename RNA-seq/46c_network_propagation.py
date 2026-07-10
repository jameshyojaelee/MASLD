#!/usr/bin/env python3
"""
46c_network_propagation.py — Network Propagation with Source Convergence

*** DEPRECATED 2026-07-04 (round-2 audit B5a). DO NOT RE-RUN INTO THE PAPER. ***
The RWR convergence result this script produced ("123 genes at 7/7, OR=4.85")
was RETRACTED 2026-05-22 (E2) as a dead pre-C2 7-source build; the live Fig5
panel was replaced by the C2 convergence distribution (see
scripts/figures/fig5_translation.R panel (e)). This script has 0 active
consumers. Its outputs (network_propagation_scores.csv / network_modules.csv /
network_rank_gainers.csv) were archived under
data/archive/dead_rwr_output_2026-07-04/. Kept here for provenance only.

Random walk with restart (RWR) from 6 source-specific seed vectors on STRING PPI.
Tests whether independent evidence sources converge on the same network modules.

Usage: python 46c_network_propagation.py
Compute: sbatch (cpu partition, ~16GB RAM for sparse matrix ops on ~14K nodes)
"""

import os
import gzip
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import norm as sparse_norm
from collections import defaultdict

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME = os.path.join(BASE, "RNA-seq/results/multi_evidence")
STRING_DIR = os.path.join(BASE, "data/string_ppi")

np.random.seed(42)

# ── 1. Load atlas ───────────────────────────────────────────────────────────
print("Loading atlas...")
atlas = pd.read_csv(os.path.join(ME, "multi_evidence_atlas.csv"))
assert {"bulk_padj", "bulk_logFC"} <= set(atlas.columns), \
    "C2: atlas missing bulk_* — rebuild 27a"
atlas_genes = set(atlas["human_symbol"].dropna().unique())
N_atlas = len(atlas)
print(f"  Atlas: {N_atlas} genes")

# ── 2. Build PPI graph from STRING v12 ──────────────────────────────────────
print("Loading STRING PPI...")
# Load protein info for ENSP → gene symbol mapping
info = pd.read_csv(os.path.join(STRING_DIR, "9606.protein.info.v12.0.txt.gz"),
                   sep="\t", compression="gzip")
ensp_to_gene = dict(zip(info["#string_protein_id"], info["preferred_name"]))

# Load interactions, filter combined_score >= 400
edges = []
n_read = 0
with gzip.open(os.path.join(STRING_DIR, "9606.protein.links.v12.0.txt.gz"), "rt") as f:
    header = f.readline()  # skip header
    for line in f:
        parts = line.strip().split()
        p1, p2, score = parts[0], parts[1], int(parts[2])
        n_read += 1
        if score < 400:
            continue
        g1 = ensp_to_gene.get(p1)
        g2 = ensp_to_gene.get(p2)
        if g1 and g2 and g1 in atlas_genes and g2 in atlas_genes and g1 != g2:
            edges.append((g1, g2, score / 1000.0))  # normalize score to [0,1]

print(f"  Read {n_read:,} interactions, kept {len(edges):,} edges "
      f"(score≥400, both genes in atlas)")

# Build node index
gene_set = set()
for g1, g2, _ in edges:
    gene_set.add(g1)
    gene_set.add(g2)
gene_list = sorted(gene_set)
gene_to_idx = {g: i for i, g in enumerate(gene_list)}
N_nodes = len(gene_list)
print(f"  PPI graph: {N_nodes} nodes, {len(edges)} edges")

# Build sparse adjacency matrix (weighted by combined_score)
row_idx = []
col_idx = []
weights = []
for g1, g2, w in edges:
    i, j = gene_to_idx[g1], gene_to_idx[g2]
    row_idx.extend([i, j])
    col_idx.extend([j, i])
    weights.extend([w, w])

A = sparse.csr_matrix((weights, (row_idx, col_idx)), shape=(N_nodes, N_nodes))

# Degree-normalized transition matrix: W = D^{-1} A
degrees = np.array(A.sum(axis=1)).flatten()
degrees[degrees == 0] = 1  # avoid division by zero
D_inv = sparse.diags(1.0 / degrees)
W = D_inv @ A

print(f"  Transition matrix: {W.shape}, density={W.nnz / W.shape[0]**2:.4f}")

# ── 3. Map atlas features to PPI nodes ──────────────────────────────────────
# Create atlas lookup by gene symbol
atlas_dict = atlas.set_index("human_symbol").to_dict("index")

def get_val(gene, col, default=0.0):
    """Get atlas column value for a gene, defaulting to 0."""
    if gene not in atlas_dict:
        return default
    v = atlas_dict[gene].get(col)
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return default
    return float(v)

def safe_neglog10(x):
    if x <= 0 or np.isnan(x):
        return 0.0
    return -np.log10(max(x, 1e-300))

# Define 6 source seed vectors (continuous weights)
print("Computing source seed vectors...")
seed_vectors = {}

for si, (name, compute_fn) in enumerate({
    "S1": lambda g: safe_neglog10(get_val(g, "bulk_padj", 1.0)) * abs(get_val(g, "bulk_logFC", 0.0)),
    "S2": lambda g: max(
        get_val(g, "coloc_pp4"), get_val(g, "broadaway_coloc_pp4"),
        get_val(g, "best_liver_enzyme_pp4"), get_val(g, "ukbb_alt_coloc_pp4"),
        get_val(g, "sceqtl_coloc_best_pp4"),
        0.5 if (get_val(g, "ieqtl_disease_interaction") == 1) else 0,
        0.5 if (get_val(g, "sceqtl_twas_best_fdr", 1.0) < 0.05) else 0),
    "S3": lambda g: max(0, -get_val(g, "essentiality_chronos", 0.0)),
    "S4": lambda g: (
        abs(get_val(g, "scenic_regulon_activity_diff")) +
        safe_neglog10(get_val(g, "mouse_da_padj", 1.0))),
    "S5": lambda g: get_val(g, "spatial_morans_i") * (1 if get_val(g, "spatial_is_svg") else 0),
    "S6": lambda g: (
        safe_neglog10(get_val(g, "sc_hepatocyte_padj", 1.0)) +
        get_val(g, "liana_n_diff_interactions")),
}.items(), start=1):
    vec = np.array([compute_fn(g) for g in gene_list])
    vec = np.clip(vec, 0, 20)  # clip extreme values
    # Normalize to sum to 1
    vec_sum = vec.sum()
    if vec_sum > 0:
        vec = vec / vec_sum
    seed_vectors[name] = vec
    n_nonzero = np.sum(vec > 0)
    print(f"  {name}: {n_nonzero} non-zero seeds ({n_nonzero/N_nodes*100:.1f}%)")

# ── 4. Random walk with restart ────────────────────────────────────────────
ALPHA = 0.5   # restart probability
TOL = 1e-6
MAX_ITER = 100

def rwr(W, seed, alpha=ALPHA, tol=TOL, max_iter=MAX_ITER):
    """Random walk with restart. Returns converged probability vector."""
    p = seed.copy()
    for it in range(max_iter):
        p_new = alpha * W.T @ p + (1 - alpha) * seed
        # Normalize
        p_sum = p_new.sum()
        if p_sum > 0:
            p_new = p_new / p_sum
        diff = np.abs(p_new - p).sum()
        p = p_new
        if diff < tol:
            break
    return p

print("\nRunning RWR for each source...")
propagated = {}
for name, seed in seed_vectors.items():
    p = rwr(W, seed)
    propagated[name] = p
    top5_idx = np.argsort(-p)[:5]
    top5_genes = [gene_list[i] for i in top5_idx]
    print(f"  {name}: converged, top 5 = {top5_genes}")

# ── 5. Source convergence score ─────────────────────────────────────────────
TOP_PCT = 0.10  # top 10%
k_threshold = int(N_nodes * TOP_PCT)

# For each gene, count how many propagated vectors rank it in top 10%
convergence = np.zeros(N_nodes, dtype=int)
rank_pctiles = np.zeros((N_nodes, 6))

for si, name in enumerate(["S1", "S2", "S3", "S4", "S5", "S6"]):
    p = propagated[name]
    ranks = np.argsort(np.argsort(-p)) + 1  # 1-based rank
    pctile = ranks / N_nodes
    rank_pctiles[:, si] = pctile
    convergence += (ranks <= k_threshold).astype(int)

# Geometric mean of rank percentiles (smoother metric)
# Use rank percentile (lower = better), take geometric mean
geo_mean_pctile = np.exp(np.mean(np.log(np.clip(rank_pctiles, 1e-10, 1)), axis=1))

print(f"\nConvergence score distribution:")
for c in range(7):
    n = np.sum(convergence == c)
    print(f"  {c} sources: {n} genes ({n/N_nodes*100:.1f}%)")

# ── 6. Map back to full atlas + rank gain analysis ─────────────────────────
# Build results for PPI genes
ppi_results = pd.DataFrame({
    "human_symbol": gene_list,
    "network_convergence": convergence,
    "geo_mean_rank_pctile": geo_mean_pctile,
})

# Add per-source propagated ranks
for si, name in enumerate(["S1", "S2", "S3", "S4", "S5", "S6"]):
    ppi_results[f"propagated_{name}_pctile"] = rank_pctiles[:, si]

# Merge with atlas layers_active
atlas_slim = atlas[["human_symbol", "layers_active"]].copy()
ppi_results = ppi_results.merge(atlas_slim, on="human_symbol", how="left")

# Rank by network_convergence (then geo_mean as tiebreaker)
ppi_results["network_rank"] = ppi_results[["network_convergence", "geo_mean_rank_pctile"]].apply(
    lambda x: (-x.iloc[0], x.iloc[1]), axis=1
).rank(method="first").astype(int)

# Rank by layers_active (original naive counting)
ppi_results["naive_rank"] = ppi_results["layers_active"].rank(
    method="first", ascending=False).astype(int)

# Rank gain = naive_rank - network_rank (positive = network propagation improves rank)
ppi_results["rank_gain"] = ppi_results["naive_rank"] - ppi_results["network_rank"]

# Merge back to full atlas (genes not in PPI get NaN)
full_results = atlas[["human_symbol", "ensembl_id", "layers_active"]].merge(
    ppi_results[["human_symbol", "network_convergence", "geo_mean_rank_pctile",
                 "network_rank", "rank_gain"] +
                [f"propagated_{name}_pctile" for name in
                 ["S1", "S2", "S3", "S4", "S5", "S6"]]],
    on="human_symbol", how="left"
)

full_results.to_csv(os.path.join(ME, "network_propagation_scores.csv"), index=False)
print(f"\nNetwork propagation scores saved: {len(full_results)} genes")

# Top rank gainers
rank_gainers = ppi_results.nlargest(50, "rank_gain")[
    ["human_symbol", "network_convergence", "layers_active",
     "network_rank", "naive_rank", "rank_gain"]]
rank_gainers.to_csv(os.path.join(ME, "network_rank_gainers.csv"), index=False)
print(f"Top 50 rank gainers saved")

# ── 7. Module detection (Leiden) on top convergence subgraph ────────────────
TOP_K_MODULES = 200
top_genes_idx = np.argsort(-convergence)[:TOP_K_MODULES]
top_genes = [gene_list[i] for i in top_genes_idx]
top_gene_set = set(top_genes)

try:
    import leidenalg
    import igraph as ig

    # Build subgraph
    sub_edges = []
    sub_weights = []
    top_idx_map = {gene_list[i]: new_i for new_i, i in enumerate(top_genes_idx)}

    for g1, g2, w in edges:
        if g1 in top_idx_map and g2 in top_idx_map:
            sub_edges.append((top_idx_map[g1], top_idx_map[g2]))
            sub_weights.append(w)

    if len(sub_edges) > 0:
        G = ig.Graph(n=len(top_genes), edges=sub_edges, directed=False)
        G.vs["name"] = top_genes
        G.es["weight"] = sub_weights

        partition = leidenalg.find_partition(
            G, leidenalg.ModularityVertexPartition, weights="weight",
            seed=42, n_iterations=10
        )

        module_assignments = pd.DataFrame({
            "human_symbol": top_genes,
            "module": [partition.membership[i] for i in range(len(top_genes))],
            "network_convergence": convergence[top_genes_idx],
        })
        module_assignments = module_assignments.sort_values(["module", "network_convergence"],
                                                             ascending=[True, False])

        # Module sizes
        module_sizes = module_assignments.groupby("module").size().reset_index(name="n_genes")
        print(f"\nLeiden modules (top {TOP_K_MODULES} convergence genes):")
        for _, row in module_sizes.iterrows():
            genes_in_mod = module_assignments[module_assignments["module"] == row["module"]]["human_symbol"].tolist()
            print(f"  Module {int(row['module'])}: {int(row['n_genes'])} genes — "
                  f"{', '.join(genes_in_mod[:5])}{'...' if len(genes_in_mod)>5 else ''}")

        module_assignments.to_csv(os.path.join(ME, "network_modules.csv"), index=False)
        print("Module assignments saved")
    else:
        print("No edges in top-200 subgraph; skipping Leiden")

except ImportError:
    print("\nleidenalg/igraph not available; skipping module detection")
    print("Install with: pip install leidenalg python-igraph")

# ── 8. Conserved enrichment of network convergence genes ───────────────
conserved = set(atlas[atlas["is_conserved"] == True]["human_symbol"])
ppi_gene_set = set(gene_list)

# Among PPI genes, are network-convergent genes enriched for Conserved?
top_conv = set(gene_list[i] for i in np.argsort(-convergence)[:500])
n_conv_cc = len(top_conv & conserved & ppi_gene_set)
n_conv_total = len(top_conv & ppi_gene_set)
n_cc_ppi = len(conserved & ppi_gene_set)

# Fisher's exact test
from scipy.stats import fisher_exact
a = n_conv_cc
b = n_conv_total - n_conv_cc
c = n_cc_ppi - n_conv_cc
d = N_nodes - n_conv_total - c
odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
print(f"\nConserved enrichment in top-500 convergent genes:")
print(f"  {a}/{n_conv_total} convergent are CC ({a/n_conv_total*100:.1f}%)")
print(f"  {n_cc_ppi}/{N_nodes} PPI genes are CC ({n_cc_ppi/N_nodes*100:.1f}%)")
print(f"  OR={odds_ratio:.2f}, Fisher p={pval:.2e}")

# ── Summary ─────────────────────────────────────────────────────────────────
print(f"\n=== RESULTS SUMMARY ===")
print(f"PPI graph: {N_nodes} genes, {len(edges)} edges")
print(f"Genes with convergence >= 5: {np.sum(convergence >= 5)}")
print(f"Genes with convergence >= 6: {np.sum(convergence >= 6)}")
print(f"Genes with convergence == 6: {np.sum(convergence == 6)}")

# Top 20 by convergence
print(f"\nTop 20 by network convergence:")
top20 = ppi_results.nlargest(20, "network_convergence")[
    ["human_symbol", "network_convergence", "layers_active", "rank_gain"]]
print(top20.to_string(index=False))

# Rank gain statistics
gains = ppi_results[ppi_results["rank_gain"] > 500]
print(f"\nGenes gaining >500 ranks from network propagation: {len(gains)}")
if len(gains) > 0:
    print(gains[["human_symbol", "layers_active", "network_convergence",
                  "naive_rank", "network_rank", "rank_gain"]].head(10).to_string(index=False))

# DGIdb druggable enrichment of rank gainers
try:
    dgidb = pd.read_csv(os.path.join(BASE, "RNA-seq/results/drug_repurposing/dgidb_drug_gene_interactions.csv"))
    dgidb_genes = set(dgidb[dgidb["approved"] == True]["gene"].unique())
    n_gainer_drug = len(set(gains["human_symbol"]) & dgidb_genes)
    print(f"  {n_gainer_drug}/{len(gains)} rank gainers are DGIdb druggable")
except Exception:
    pass

print("\nDone.")
