#!/usr/bin/env python
"""
Analysis E: characterize the GENUINE semaglutide-refractory residue (the real white-space).

After the meta-level correction, "refractory" = convergence targets whose disease
dysregulation is NOT reversed and that move WITH disease (not merely underpowered).
RevC warning: do NOT run ORA on the small residue against a whole-genome background
(that resurrects spurious 'residue pathways'). Correct background = the disease-significant
convergence pool the residue is drawn from. Expectation: heterogeneous, no coherent pathway.

Output: RNA-seq/results/glp1ra/glp1ra_residue_targets.csv + summary.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
G = os.path.join(ROOT, "RNA-seq/results/glp1ra")
PART = os.path.join(G, "glp1ra_target_reversal_partition.csv")
DRUG = os.path.join(ROOT, "data/external/drug_targets/drug_target_classification.tsv")
GMT = os.path.join(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt")

t = pd.read_csv(PART)
drug = pd.read_csv(DRUG, sep="\t").rename(columns={"symbol": "human_symbol"})
t = t.merge(drug[["human_symbol", "drug_dev_status", "max_phase_masld", "pharos_tdl"]],
            on="human_symbol", how="left")

# disease-significant convergence pool = correct ORA background
pool = t[(t["conv_tier1"] == True) & (t["mouse_disease_sig"] == True)].copy()
# genuine residue = refractory (moving with disease), highest-confidence = also mouse-concordant
residue = pool[pool["mouse_reversal"] == "refractory"]
residue_conc = residue[residue["mouse_human_concordant"] == True]

# ORA against the disease-sig convergence pool (proper background)
def load_gmt(path, universe):
    s = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                s[p[0]] = set(p[2:]) & universe
    return s

universe = set(pool["human_symbol"])
hall = load_gmt(GMT, universe)
N = len(universe)
q = set(residue["human_symbol"])
rows = []
for name, gs in hall.items():
    if len(gs) < 3:
        continue
    k = len(q & gs)
    if k == 0:
        continue
    rows.append(dict(pathway=name, k=k, K=len(gs), q=len(q),
                     pval=hypergeom.sf(k - 1, N, len(gs), len(q))))
ora = pd.DataFrame(rows).sort_values("pval") if rows else pd.DataFrame()
if len(ora):
    ora["padj"] = multipletests(ora["pval"], method="fdr_bh")[1]

# druggability tiering of the residue
def tier(df):
    dd = df["drug_dev_status"].fillna("unclassified")
    drugged = dd.isin(["masld_approved", "masld_clinical"])
    tractable = df["pharos_tdl"].isin(["Tclin", "Tchem"]) | df["druggability_tier"].notna()
    never = (~dd.isin(["masld_approved", "masld_clinical", "masld_preclinical",
                       "drugged_other_indication"]))
    return len(df), int(drugged.sum()), int(tractable.sum()), int(never.sum())

residue.sort_values("convergence_score", ascending=False)[
    ["human_symbol", "convergence_score", "mdis", "mtrt", "fraction_corrected",
     "mouse_human_concordant", "druggability_tier", "drug_dev_status", "pharos_tdl"]].to_csv(
    os.path.join(G, "glp1ra_residue_targets.csv"), index=False)

n, drugged, tract, never = tier(residue)
L = ["Analysis E: genuine semaglutide-refractory residue (real white-space)", "=" * 60,
     f"disease-sig convergence-Tier1 pool (ORA background): {len(pool)}",
     f"genuine refractory (moving WITH disease): {len(residue)}  "
     f"(mouse-disease-concordant high-conf: {len(residue_conc)})",
     f"  druggability: already MASLD approved/clinical={drugged}; tractable(Tclin/Tchem/druggable)={tract}; "
     f"never-drugged={never}",
     "",
     "Proper-background ORA (residue vs disease-sig convergence pool):",
     (ora[ora["padj"] < 0.1][["pathway", "k", "K", "pval", "padj"]].to_string(index=False)
      if len(ora) and (ora["padj"] < 0.1).any()
      else "  NO Hallmark pathway enriched at padj<0.1 => residue is HETEROGENEOUS (no coherent program)."),
     "",
     "Residue genes (mouse-concordant, highest-confidence), by convergence score:",
     residue_conc.sort_values("convergence_score", ascending=False)[
         ["human_symbol", "mdis", "mtrt", "drug_dev_status", "pharos_tdl"]].head(30).round(2).to_string(index=False)]
summary = "\n".join(L)
with open(os.path.join(G, "glp1ra_residue_summary.txt"), "w") as fh:
    fh.write(summary + "\n")
print(summary)
