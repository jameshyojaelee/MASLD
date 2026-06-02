#!/usr/bin/env python
"""
342d_consensus_signature.py - Build consensus polyploid signature from Richter + Yin + Katsuda.

Strategy:
- Each source contributes up/down (Richter, Katsuda) or candidate (Yin, direction-agnostic) gene sets.
- Cross-source overlap analysis: how many genes appear in ≥2 sources?
- Consensus signature: gene in ≥2 sources AND consistent direction (where direction available).
- Save per-source + consensus signatures in human gene symbols.
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT  = ROOT / "data/ploidy_signatures"

# --- Load per-source signatures ---
ric_up = pd.read_csv(OUT / "richter2021_polyploid_up.tsv", sep="\t")
ric_dn = pd.read_csv(OUT / "richter2021_polyploid_down.tsv", sep="\t")
kat_up = pd.read_csv(OUT / "katsuda2019_polyploid_up.tsv", sep="\t")
kat_dn = pd.read_csv(OUT / "katsuda2019_polyploid_down.tsv", sep="\t")
yin_all = pd.read_csv(OUT / "yin2024_polyploid_candidates.tsv", sep="\t")
yin_cons = pd.read_csv(OUT / "yin2024_conserved_polyploid.tsv", sep="\t")

ric_up_set = set(ric_up["human_symbol"].dropna().str.upper())
ric_dn_set = set(ric_dn["human_symbol"].dropna().str.upper())
kat_up_set = set(kat_up["human_symbol"].dropna().str.upper())
kat_dn_set = set(kat_dn["human_symbol"].dropna().str.upper())
yin_set    = set(yin_all["human_symbol"].dropna().str.upper())
yin_cons_set = set(yin_cons["human_symbol"].dropna().str.upper())

print(f"Source signatures (human gene symbols):")
print(f"  Richter UP:   {len(ric_up_set)}")
print(f"  Richter DOWN: {len(ric_dn_set)}")
print(f"  Katsuda UP:   {len(kat_up_set)}")
print(f"  Katsuda DOWN: {len(kat_dn_set)}")
print(f"  Yin all:      {len(yin_set)} (direction agnostic)")
print(f"  Yin conserved: {len(yin_cons_set)} (≥2 strain backgrounds)")

# --- Cross-source overlap ---
print(f"\n=== Pairwise overlaps ===")
print(f"  Richter-UP ∩ Katsuda-UP: {len(ric_up_set & kat_up_set)}  -> {sorted(ric_up_set & kat_up_set)}")
print(f"  Richter-DN ∩ Katsuda-DN: {len(ric_dn_set & kat_dn_set)}  -> {sorted(ric_dn_set & kat_dn_set)}")
print(f"  Richter-UP ∩ Katsuda-DN: {len(ric_up_set & kat_dn_set)}  -> {sorted(ric_up_set & kat_dn_set)} (DIRECTION CONFLICT)")
print(f"  Richter-DN ∩ Katsuda-UP: {len(ric_dn_set & kat_up_set)}  -> {sorted(ric_dn_set & kat_up_set)} (DIRECTION CONFLICT)")
print(f"  Richter-UP ∩ Yin: {len(ric_up_set & yin_set)}")
print(f"  Richter-DN ∩ Yin: {len(ric_dn_set & yin_set)}")
print(f"  Katsuda-UP ∩ Yin: {len(kat_up_set & yin_set)}")
print(f"  Katsuda-DN ∩ Yin: {len(kat_dn_set & yin_set)}")

# --- Consensus signature ---
# UP: gene in (Richter-UP ∩ Katsuda-UP) OR appears in both Richter-UP and Yin OR Katsuda-UP and Yin
# Direction-consistent: 2 of 3 sources agree on UP
consensus_up = (ric_up_set & kat_up_set) | (ric_up_set & yin_set) | (kat_up_set & yin_set)
consensus_dn = (ric_dn_set & kat_dn_set) | (ric_dn_set & yin_set) | (kat_dn_set & yin_set)
# Remove direction-conflicts: drop genes that appear in BOTH consensus_up and consensus_dn
conflicts = consensus_up & consensus_dn
consensus_up = consensus_up - conflicts
consensus_dn = consensus_dn - conflicts
print(f"\n=== Consensus signature ===")
print(f"  Consensus UP   (≥2 sources, no conflict): {len(consensus_up)} genes")
print(f"  Consensus DOWN (≥2 sources, no conflict): {len(consensus_dn)} genes")
if conflicts:
    print(f"  Direction conflicts (dropped): {sorted(conflicts)}")
print(f"\n  UP genes:")
print(f"    {sorted(consensus_up)}")
print(f"\n  DOWN genes:")
print(f"    {sorted(consensus_dn)}")

# Build a richer table per consensus gene: which sources, what effect size, etc.
def build_table(genes, ric_df_up, ric_df_dn, kat_df_up, kat_df_dn, yin_df, direction):
    rows = []
    for g in genes:
        row = {"human_symbol": g, "direction": direction}
        ric = ric_df_up if direction == "up" else ric_df_dn
        kat = kat_df_up if direction == "up" else kat_df_dn
        row["in_richter"] = g in set(ric["human_symbol"].dropna().str.upper())
        row["in_katsuda"] = g in set(kat["human_symbol"].dropna().str.upper())
        row["in_yin"]     = g in set(yin_df["human_symbol"].dropna().str.upper())
        row["n_sources"]  = sum([row["in_richter"], row["in_katsuda"], row["in_yin"]])
        ric_row = ric[ric["human_symbol"].str.upper() == g].head(1)
        kat_row = kat[kat["human_symbol"].str.upper() == g].head(1)
        row["richter_LFC"] = float(ric_row["log2FoldChange"].iloc[0]) if len(ric_row) else np.nan
        row["richter_padj"] = float(ric_row["padj"].iloc[0]) if len(ric_row) else np.nan
        row["katsuda_LFC"] = float(kat_row["logFC"].iloc[0]) if len(kat_row) else np.nan
        row["katsuda_padj"] = float(kat_row["adj.P.Val"].iloc[0]) if len(kat_row) else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values("n_sources", ascending=False)

tab_up = build_table(consensus_up, ric_up, ric_dn, kat_up, kat_dn, yin_all, "up")
tab_dn = build_table(consensus_dn, ric_up, ric_dn, kat_up, kat_dn, yin_all, "down")
tab_up.to_csv(OUT / "consensus_polyploid_up.tsv", sep="\t", index=False)
tab_dn.to_csv(OUT / "consensus_polyploid_down.tsv", sep="\t", index=False)
print(f"\nWrote consensus_polyploid_up.tsv ({len(tab_up)} genes)")
print(f"Wrote consensus_polyploid_down.tsv ({len(tab_dn)} genes)")

# --- Permutation null: would we see this overlap by chance? ---
# Hypergeometric: pool of ~20000 human genes; Richter-UP=100, Katsuda-UP=100
# How many expected overlaps by chance? n = K*M/N
GENE_POOL = 20000
expected_up = (len(ric_up_set) * len(kat_up_set)) / GENE_POOL
expected_dn = (len(ric_dn_set) * len(kat_dn_set)) / GENE_POOL
from scipy.stats import hypergeom
obs_up = len(ric_up_set & kat_up_set)
obs_dn = len(ric_dn_set & kat_dn_set)
p_up = hypergeom.sf(obs_up - 1, GENE_POOL, len(ric_up_set), len(kat_up_set))
p_dn = hypergeom.sf(obs_dn - 1, GENE_POOL, len(ric_dn_set), len(kat_dn_set))
print(f"\n=== Hypergeometric test: Richter ∩ Katsuda overlap vs null (gene pool=20k) ===")
print(f"  UP:   observed={obs_up}  expected_by_chance={expected_up:.2f}  hyper-p={p_up:.3g}")
print(f"  DOWN: observed={obs_dn}  expected_by_chance={expected_dn:.2f}  hyper-p={p_dn:.3g}")

# --- Save overlap summary ---
summary = pd.DataFrame([
    {"pair": "Richter-UP ∩ Katsuda-UP", "n_overlap": len(ric_up_set & kat_up_set), "hyper_p": p_up},
    {"pair": "Richter-DN ∩ Katsuda-DN", "n_overlap": len(ric_dn_set & kat_dn_set), "hyper_p": p_dn},
    {"pair": "Richter-UP ∩ Yin (any)", "n_overlap": len(ric_up_set & yin_set)},
    {"pair": "Richter-DN ∩ Yin (any)", "n_overlap": len(ric_dn_set & yin_set)},
    {"pair": "Katsuda-UP ∩ Yin (any)", "n_overlap": len(kat_up_set & yin_set)},
    {"pair": "Katsuda-DN ∩ Yin (any)", "n_overlap": len(kat_dn_set & yin_set)},
    {"pair": "Consensus UP (≥2 sources)", "n_overlap": len(consensus_up)},
    {"pair": "Consensus DOWN (≥2 sources)", "n_overlap": len(consensus_dn)},
])
summary.to_csv(OUT / "signature_overlap_matrix.csv", index=False)
print(f"\nWrote signature_overlap_matrix.csv")
print(summary.to_string(index=False))
