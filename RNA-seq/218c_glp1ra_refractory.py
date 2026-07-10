#!/usr/bin/env python
"""
Analysis C: characterize the semaglutide-REFRACTORY convergence white-space.

Question the manuscript needs answered: the ~41% of top convergence targets that
semaglutide does NOT reverse -- what biology is it, is it druggable, and is it
already being chased by next-generation (non-GLP1-only) drugs?

Inputs:
  glp1ra_target_reversal_partition.csv      (Analysis B output)
  data/external/drug_targets/drug_target_classification.tsv  (drug-dev status)
  Analysis/.../genesets/hallmark.gmt         (MSigDB Hallmark for ORA)

Outputs (RNA-seq/results/glp1ra/):
  glp1ra_refractory_targets.csv       refractory Tier-1 targets + drug annotations
  glp1ra_refractory_pathways.csv      Hallmark over-representation (refractory vs reversed)
  glp1ra_refractory_summary.txt
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
t = t.merge(drug[["human_symbol", "drug_dev_status", "max_phase_masld",
                  "n_masld_trials", "pharos_tdl"]], on="human_symbol", how="left")

universe = set(t["human_symbol"].dropna())

# ---- target sets (Tier-1 convergence) ----
c = t[t["conv_tier1"] == True].copy()
refr = c[c["reversal_class"] == "refractory"]
revd = c[c["reversal_class"] == "reversed"]
refr_syms, revd_syms = set(refr["human_symbol"]), set(revd["human_symbol"])

# ---- druggability / mechanism of the refractory white-space ----
def drug_profile(df, name):
    n = len(df)
    st = df["drug_dev_status"].fillna("unclassified").value_counts()
    # "novel white-space" = refractory AND not an approved/clinical MASLD drug target
    drugged_masld = df["drug_dev_status"].isin(["masld_approved", "masld_clinical"])
    tractable = df["pharos_tdl"].isin(["Tclin", "Tchem"]) | df["druggability_tier"].notna()
    return (f"\n[{name}] n={n}\n"
            f"  drug_dev_status: " + ", ".join(f"{k}={v}" for k, v in st.items()) + "\n"
            f"  already MASLD approved/clinical: {int(drugged_masld.sum())}\n"
            f"  pharmacologically tractable (Tclin/Tchem or druggable): {int(tractable.sum())} "
            f"({100*tractable.mean():.0f}%)\n"
            f"  => NOVEL undrugged-for-MASLD white-space: {int((~drugged_masld).sum())}")

# ---- Hallmark over-representation ----
def load_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                sets[p[0]] = set(p[2:]) & universe
    return sets

hall = load_gmt(GMT)
N = len(universe)

def ora(query, label):
    rows = []
    q = query & universe
    for name, gs in hall.items():
        K = len(gs)
        if K < 5:
            continue
        k = len(q & gs)
        if k == 0:
            continue
        p = hypergeom.sf(k - 1, N, K, len(q))
        rows.append(dict(pathway=name, set_source=label, n_overlap=k, set_size=K,
                         query_size=len(q), pval=p,
                         genes=";".join(sorted(q & gs))))
    r = pd.DataFrame(rows)
    if len(r):
        r["padj"] = multipletests(r["pval"], method="fdr_bh")[1]
        r = r.sort_values("pval")
    return r

ora_refr = ora(refr_syms, "refractory")
ora_revd = ora(revd_syms, "reversed")
ora_all = pd.concat([ora_refr, ora_revd], ignore_index=True)
ora_all.to_csv(os.path.join(G, "glp1ra_refractory_pathways.csv"), index=False)

# ---- refractory target table (ranked) ----
refr_out = refr.sort_values("convergence_score", ascending=False)[
    ["human_symbol", "convergence_score", "n_modalities_active", "concordance_state",
     "disease_dir_liver", "druggability_tier", "drug_dev_status", "max_phase_masld",
     "pharos_tdl", "mouse_reversal", "serum_reversal"]]
refr_out.to_csv(os.path.join(G, "glp1ra_refractory_targets.csv"), index=False)

# ---- summary ----
L = ["Analysis C: semaglutide-refractory convergence white-space", "=" * 58]
L.append(drug_profile(refr, "REFRACTORY Tier-1 convergence"))
L.append(drug_profile(revd, "REVERSED Tier-1 convergence (contrast)"))
L.append("\n=== Hallmark pathways enriched in REFRACTORY set (padj<0.1) ===")
sig_refr = ora_refr[ora_refr["padj"] < 0.1].head(15)
L.append(sig_refr[["pathway", "n_overlap", "set_size", "pval", "padj"]].to_string(index=False)
         if len(sig_refr) else "  (none at padj<0.1)")
L.append("\n=== Hallmark pathways enriched in REVERSED set (padj<0.1, contrast) ===")
sig_revd = ora_revd[ora_revd["padj"] < 0.1].head(15)
L.append(sig_revd[["pathway", "n_overlap", "set_size", "pval", "padj"]].to_string(index=False)
         if len(sig_revd) else "  (none at padj<0.1)")
L.append("\n=== Top 25 refractory Tier-1 targets (the residual list) ===")
L.append(refr_out.head(25).to_string(index=False))

summary = "\n".join(L)
with open(os.path.join(G, "glp1ra_refractory_summary.txt"), "w") as fh:
    fh.write(summary + "\n")
print(summary)
print("\n[done] wrote refractory targets + pathways + summary")
