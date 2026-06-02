"""Smoke test for the d5 ortholog wire-up.

(a) Reads top-20 d5_mouse_orthologs.csv rows and verifies mouse_gene populated.
(b) Calls ortholog_mapper.{human_to_mouse, mouse_to_human, map_hit_table}.
(c) Calls concordance_scorer.concordance() with dummy top-K dicts to ensure
    no crash and a valid DataFrame back.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "shared"))

import ortholog_mapper as om  # noqa: E402
import concordance_scorer as cs  # noqa: E402

HITS = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "Analysis/Perturbation/data/hits/d5_mouse_orthologs.csv"
)

print("=== (a) top-20 hits ===")
df = pd.read_csv(HITS).head(20)
print(df[["human_gene", "mouse_gene", "mouse_ensembl", "mouse_gene_source"]])
nonempty = df["mouse_gene"].notna().sum()
print(f"mouse_gene populated: {nonempty}/20")
assert nonempty >= 18, "expected >=18 of top-20 to have a mouse symbol"

print("\n=== (b) ortholog_mapper round trip ===")
for sym in ["TP53", "ALB", "HNF4A", "HKDC1", "THRB"]:
    print(f"  {sym} -> {om.human_to_mouse(sym)}")
for sym in ["Trp53", "Alb", "Hnf4a"]:
    print(f"  {sym} -> {om.mouse_to_human(sym)}")

print("\n=== (b2) map_hit_table ===")
toy = pd.DataFrame({"human_gene": ["TP53", "ALB", "HNF4A", "MADEUPGENE"]})
mapped = om.map_hit_table(toy, column="human_gene")
print(mapped[[c for c in ("human_gene", "mouse_gene", "mouse_gene_symbol", "ortholog_type") if c in mapped.columns]])

print("\n=== (c) concordance_scorer.concordance() dummy run ===")
ortho = df[["human_gene", "mouse_gene"]].dropna()
human_topk = {hg: [f"DOWN_{i}" for i in range(50)] for hg in ortho["human_gene"]}
mouse_topk = {mg: [f"DOWN_{i}" for i in range(50)] for mg in ortho["mouse_gene"]}
out = cs.concordance(human_topk, mouse_topk, ortho, k=50)
print(out.head())
assert (out["jaccard_top_k"] == 1.0).all(), "self-Jaccard should be 1.0"
print(f"rows: {len(out)} — all jaccard==1.0 OK")

print("\nALL SMOKE TESTS PASSED")
