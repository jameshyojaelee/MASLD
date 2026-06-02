"""Audit: hep__24 vs each bulk NMF k=6 program — clean Jaccard, cosine, hyper p."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
W_PATH = ROOT / "data/reference_signatures/bulk_nmf_k6_W.tsv"
MOD_PATH = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv"
OUT_DIR = ROOT / "Analysis/SingleCell/scripts/hotspot_modules/audits"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TOPN = 50

# ---- bulk NMF top-50 per program ----
W = pd.read_csv(W_PATH, sep="\t").set_index("gene")
prog_cols = [c for c in W.columns if c.startswith("P")]
print(f"W: {W.shape} genes x {len(prog_cols)} programs")

# Top-50 table
top_rows = []
top_sets: dict[str, set[str]] = {}
for p in prog_cols:
    s = W[p].sort_values(ascending=False).head(TOPN)
    top_sets[p] = set(s.index)
    for rank, (g, w) in enumerate(s.items(), start=1):
        top_rows.append({"program": p, "rank": rank, "gene": g, "weight": w})
top_df = pd.DataFrame(top_rows)
top_df.to_csv(OUT_DIR / "bulk_nmf_k6_top50.tsv", sep="\t", index=False)
print(f"Wrote {OUT_DIR / 'bulk_nmf_k6_top50.tsv'}")

# ---- HKDC1 placement ----
print("\n=== HKDC1 placement across 6 NMF programs ===")
if "HKDC1" in W.index:
    hkdc1_w = W.loc["HKDC1"]
    print(hkdc1_w.to_string())
    # rank within each program
    for p in prog_cols:
        rank = (W[p] > W.loc["HKDC1", p]).sum() + 1  # 1-based rank
        total = W[p].shape[0]
        in_top50 = "HKDC1" in top_sets[p]
        print(f"  {p}: weight={W.loc['HKDC1', p]:.2f}, rank={rank}/{total}, in_top50={in_top50}")
    # which program has the highest weight for HKDC1?
    best_p = hkdc1_w.idxmax()
    print(f"  HKDC1 max-weight program: {best_p} (weight={hkdc1_w.max():.2f})")
else:
    print("HKDC1 NOT FOUND in bulk NMF W matrix")

# ---- hep mod 24 top-50 ----
mg = pd.read_csv(MOD_PATH, sep="\t")
hep24 = mg.query("module == 24").sort_values("weight", ascending=False)
print(f"\nhep mod 24 total genes: {len(hep24)}")
hep24_top50 = hep24.head(TOPN)
hep24_top50_set = set(hep24_top50["gene"])
print(f"hep mod 24 top-50: {list(hep24_top50_set)[:15]}...")

# Weights as a Series indexed by gene (for cosine)
hep24_w = hep24.set_index("gene")["weight"]

# ---- pairwise comparison ----
# Background: protein-coding gene universe ~ overlap of detected genes
# For a fair hypergeometric, use intersection of both gene-universe spaces
bulk_universe = set(W.index)
hot_universe = set(mg["gene"])  # all hotspot-detected genes
universe = bulk_universe | hot_universe
N_universe = len(universe)
print(f"\nUnion universe size: {N_universe}")

# Also report ~17K protein-coding background per request
N_PC = 17000

rows = []
for p in prog_cols:
    top_bulk = top_sets[p]
    # Symmetric Jaccard @ 50
    inter = hep24_top50_set & top_bulk
    union = hep24_top50_set | top_bulk
    jac = len(inter) / max(len(union), 1)

    # Cosine over union, zero-padded for genes missing in either vector
    union_genes = sorted(set(hep24_w.index) | set(W.index))
    u = hep24_w.reindex(union_genes).fillna(0.0).values
    v = W[p].reindex(union_genes).fillna(0.0).values
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    cos_union = float(np.dot(u, v) / (nu * nv)) if nu * nv > 0 else 0.0

    # Cosine over intersection only (the script's behavior)
    common = list(set(hep24_w.index) & set(W.index))
    if common:
        u_i = hep24_w.reindex(common).values
        v_i = W[p].reindex(common).values
        nu_i, nv_i = np.linalg.norm(u_i), np.linalg.norm(v_i)
        cos_inter = float(np.dot(u_i, v_i) / (nu_i * nv_i)) if nu_i * nv_i > 0 else 0.0
    else:
        cos_inter = 0.0

    # Hypergeometric: overlap of top-50 sets in PC background
    # k=inter, K=50 (bulk top), n=50 (hot top), N=17000
    k = len(inter)
    p_hyper = hypergeom.sf(k - 1, N_PC, TOPN, TOPN)

    rows.append({
        "program": p,
        "jaccard_top50": jac,
        "n_overlap_top50": k,
        "overlap_genes": ",".join(sorted(inter)) if inter else "",
        "cosine_intersection": cos_inter,    # what 504 script computes
        "cosine_zero_padded_union": cos_union,
        "hyper_p_PC17K": p_hyper,
    })

cmp_df = pd.DataFrame(rows)
cmp_df.to_csv(OUT_DIR / "hep24_vs_bulk_nmf_k6_comparison.tsv", sep="\t", index=False)
print("\n=== hep__24 vs each NMF program ===")
print(cmp_df.to_string(index=False))

# ---- Diagnose cosine values: what's driving cosine to be ~0.9? ----
print("\n=== Diagnosis: what's the dominant axis in the W column vectors? ===")
# How sparse / how concentrated is each W column?
for p in prog_cols:
    vec = W[p].values
    # Effective rank: how many genes carry ~most of the L2 mass
    sorted_sq = np.sort(vec ** 2)[::-1]
    cum = np.cumsum(sorted_sq) / sorted_sq.sum()
    n90 = int((cum < 0.90).sum() + 1)  # genes needed for 90% L2 mass
    top1_share = sorted_sq[0] / sorted_sq.sum()
    print(f"  {p}: n_genes_for_90%_L2={n90}, top-1 gene share={top1_share:.3f}, top gene={W[p].idxmax()}")

# Vector overlap on shared gene namespace (intersection) — driver of cosine ~0.9 for all
common = list(set(hep24_w.index) & set(W.index))
print(f"\nIntersection size hep24 weights vs W gene-universe: {len(common)}")
print(f"hep24 weight range: {hep24_w.min():.2f} -> {hep24_w.max():.2f}, n={len(hep24_w)}")

# Compute cosine of W column pairs to see if all NMF programs are themselves nearly collinear
print("\n=== Cosine between bulk NMF programs themselves (intersection / full) ===")
for i, p1 in enumerate(prog_cols):
    for p2 in prog_cols[i+1:]:
        v1 = W[p1].values
        v2 = W[p2].values
        c = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
        print(f"  {p1}-{p2}: {c:.3f}")
