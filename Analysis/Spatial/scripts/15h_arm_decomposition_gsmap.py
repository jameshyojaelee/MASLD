#!/usr/bin/env python
"""Fig4g data-prep: decompose the prioritized MASLD universe into its evidence arms
and test each separately for gsMap spatial-risk enrichment (circularity control).

The prioritized set (prioritized_universe_FINAL.txt = 9,882) = transcriptomic arm
(universe_transcriptomic.txt, 8,088 DE-contrast genes) UNION genetic arm
(universe_genetic.txt, 3,038 COLOC+TWAS+burden+ClinVar genes). gsMap is itself a
genetics-based method, so a naive "genetic genes flagged by a genetics method"
tautology is possible. This partitions into mutually-exclusive sets and runs the
IDENTICAL Fisher enrichment (15f/15g methodology, sample OR) per liver-enzyme trait
per cohort, to show the enrichment is CONVERGENCE-specific (both arms) rather than
driven by either arm alone.

Output: Analysis/Spatial/results/gsmap/arm_decomposition_spatial_risk.csv
Env: spatial (pyarrow). Report-only companion to scripts/figures/fig4g_gsmap_risk_in_tissue.R
"""
import math
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.feather as pf
from scipy.stats import fisher_exact

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GS = BASE / "Analysis/Spatial/results/gsmap"
UDIR = BASE / "Analysis/Spatial/results/universe_validation"
TRAITS = ["ukbb_alt", "ukbb_ast", "ukbb_ggt"]   # the replicating liver-enzyme traits

def rd(f):
    return set(l.strip() for l in open(UDIR / f) if l.strip())

rna_sig = rd("universe_transcriptomic.txt")   # 8,088 transcriptomic arm
coloc   = rd("universe_genetic.txt")          # 3,038 genetic arm
FINAL   = rd("prioritized_universe_FINAL.txt")  # 9,882 union

# groups: label -> (gene set, in_figure). Mutually-exclusive partition + union reference.
GROUPS = [
    ("Convergent",          rna_sig & coloc,  True),   # both arms  -> enriched (driver)
    ("Prioritized union",   FINAL,            True),   # T ∪ G = current 4g -> reference
    ("Transcriptomic-only", rna_sig - coloc,  True),   # T ∖ G      -> null (non-circular test)
    ("Genetic-only",        coloc - rna_sig,  True),   # G ∖ T      -> depleted
    ("Transcriptomic-all",  rna_sig,          False),  # context
    ("Genetic-all",         coloc,            False),  # context
]

def fisher_test(set_a, set_b, universe):
    a = len(set_a & set_b); b = len(set_a - set_b)
    c = len(set_b - set_a); d = len(universe - set_a - set_b)
    if min(a + b, c + d, a + c, b + d) == 0:
        return np.nan, np.nan, (a, b, c, d)
    odds, pval = fisher_exact([[a, b], [c, d]])
    return odds, pval, (a, b, c, d)

def wald_ci(a, b, c, d, z=1.959964):
    if min(a, b, c, d) == 0:
        return np.nan, np.nan
    lo = math.log((a * d) / (b * c)); se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return math.exp(lo - z * se), math.exp(lo + z * se)

# per-dataset universe = union of marker-gene symbols across sample feathers
universe = {}
for ds_dir in sorted(GS.glob("*")):
    if not ds_dir.is_dir() or ds_dir.name == "figures":
        continue
    genes = set()
    for s_dir in sorted(ds_dir.glob("*")):
        fe = s_dir / "latent_to_gene" / f"{s_dir.name}_gene_marker_score.feather"
        if fe.exists():
            genes.update(pf.read_table(fe, columns=["HUMAN_GENE_SYM"]).column("HUMAN_GENE_SYM").to_pylist())
    if genes:
        universe[ds_dir.name] = genes

grl = pd.read_csv(GS / "gene_risk_localization.csv")
risk = {(ds, tr): set(sub["gene"]) for (ds, tr), sub in grl.groupby(["dataset", "trait"])}

rows = []
for label, S, in_fig in GROUPS:
    for trait in TRAITS:
        for ds in ["gse192741", "vu"]:
            OR, p, (a, b, c, d) = fisher_test(S, risk[(ds, trait)], universe[ds])
            lo, hi = wald_ci(a, b, c, d)
            rows.append(dict(group=label, in_figure=in_fig, trait=trait, cohort=ds,
                             n_set=len(S), a=a, fisher_or=OR, ci_low=lo, ci_high=hi,
                             fisher_pval=p, neglog10p=-math.log10(max(p, 1e-300))))

res = pd.DataFrame(rows)
out = GS / "arm_decomposition_spatial_risk.csv"
res.to_csv(out, index=False)
print(f"[15h] saved {out}  ({len(res)} rows)")
print(res[res.in_figure][["group", "trait", "cohort", "n_set", "fisher_or",
                          "ci_low", "ci_high", "fisher_pval"]].to_string(index=False))
