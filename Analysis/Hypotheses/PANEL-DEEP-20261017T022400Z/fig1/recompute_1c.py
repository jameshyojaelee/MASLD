"""Independent recomputation of Figure 1C values from COLOC r2 per-study outputs.

Inputs (read-only):
  GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv   (per study x gene)
  GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv (gene state)
  GWAS/finemapping/config/gwas_trait_tier.tsv, gwas_registry.tsv
Class rule (legend Fig 1C): a gene is 'positive' (436) and supports a trait class if ANY
eligible SuSiE pair with a main tier-1/2 study of that class has PP.H4.susie > 0.5.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = ROOT / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1"
REL = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
CFG = ROOT / "GWAS/finemapping/config"

cols = ["gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "PP.H3.susie", "n_cs_pairs", "method",
        "n_cs_pairs_eligible", "eligibility_source", "ld_reliability"]
a = pd.read_csv(REL / "susie_coloc_all_gwas.csv", usecols=cols, low_memory=False)
tier = pd.read_csv(CFG / "gwas_trait_tier.tsv", sep="\t")
reg = pd.read_csv(CFG / "gwas_registry.tsv", sep="\t")
gl = pd.read_csv(REL / "gene_level_coloc_tier12.csv")

res = {}
res["n_rows"] = len(a)
res["n_studies"] = a["gwas_name"].nunique()
main = tier[(tier["placement"] == "main") & (tier["tier"].isin([1, 2]))]
res["n_main_t12_studies"] = len(main)
res["main_by_class"] = main["tier_label"].value_counts().to_dict()
res["studies_in_rows_missing_from_tier"] = sorted(set(a["gwas_name"]) - set(tier["study_name"]))
res["main_studies_without_rows"] = sorted(set(main["study_name"]) - set(a["gwas_name"]))
a = a.merge(tier[["study_name", "trait", "tier", "tier_label", "placement"]], left_on="gwas_name", right_on="study_name", how="left")
m = a[a["study_name"].isin(main["study_name"])].copy()
res["method_counts_main"] = m["method"].value_counts(dropna=False).to_dict()
res["eligibility_source_counts_main"] = m["eligibility_source"].value_counts(dropna=False).to_dict()
pos = m[m["PP.H4.susie"] > 0.5]
res["pos_pairs_main"] = len(pos)
res["pos_pairs_by_method"] = pos["method"].value_counts(dropna=False).to_dict()
res["pos_pairs_by_eligibility"] = pos["eligibility_source"].value_counts(dropna=False).to_dict()
res["pos_pairs_n_cs_eligible_zero_or_na"] = int((pos["n_cs_pairs_eligible"].fillna(0) == 0).sum())

positive = set(gl.loc[gl["susie_state_t12"] == "positive", "ensembl"])
res["gene_level_positive"] = len(positive)

def classify(p):
    per = p.groupby("ensembl")["tier_label"].agg(lambda s: set(s))
    def lab(s):
        e = "liver_enzyme" in s; d = "direct_MASLD" in s
        return "both" if (e and d) else ("enzyme_only" if e else ("direct_only" if d else "other"))
    return per.map(lab)

cls_raw = classify(pos)
res["raw_any_pp4_gt_0.5_genes"] = int(cls_raw.size)
res["raw_vs_gene_level_extra"] = sorted(set(cls_raw.index) - positive)[:20]
res["raw_vs_gene_level_missing"] = sorted(positive - set(cls_raw.index))[:20]
res["class_counts_raw"] = cls_raw.value_counts().to_dict()
# restrict to the gene-level positive set (any eligible pair)
cls = cls_raw[cls_raw.index.isin(positive)]
res["class_counts_positive436"] = cls.value_counts().to_dict()

# Alternative 1: strongest-driver class (tier of the max PP.H4.susie study)
best = pos[pos["ensembl"].isin(positive)].sort_values("PP.H4.susie", ascending=False).drop_duplicates("ensembl")
res["class_counts_strongest_driver"] = best["tier_label"].value_counts().to_dict()
# Alternative 2: PP.H4.susie > 0.8
pos8 = m[(m["PP.H4.susie"] > 0.8)]
res["class_counts_pp4_gt_0.8"] = classify(pos8).value_counts().to_dict()
# Alternative 3: tested denominator per class: genes with any SuSiE result in a class's studies
tested = m[m["PP.H4.susie"].notna()]
tcls = classify(tested)
res["genes_with_susie_result_by_class_combo"] = tcls.value_counts().to_dict()
res["genes_with_any_direct_susie_result"] = int(tested.loc[tested["tier_label"] == "direct_MASLD", "ensembl"].nunique())
res["genes_with_any_enzyme_susie_result"] = int(tested.loc[tested["tier_label"] == "liver_enzyme", "ensembl"].nunique())
# fraction positive among tested, per class
for lab in ("direct_MASLD", "liver_enzyme"):
    t = tested[tested["tier_label"] == lab]
    g = t.groupby("ensembl")["PP.H4.susie"].max()
    res[f"{lab}_tested_genes"] = int(g.size); res[f"{lab}_positive_genes"] = int((g > 0.5).sum())
# Alternative 4: SuSiE-or-ABF (ABF > 0.5)
posu = m[(m["PP.H4.susie"] > 0.5) | (m["PP.H4.abf"] > 0.5)]
res["class_counts_susie_or_abf"] = classify(posu).value_counts().to_dict()

# per-study positive gene counts and sample size
ps = pos.groupby("gwas_name")["ensembl"].nunique()
tr = tested.groupby("gwas_name")["ensembl"].nunique()
st = main.merge(reg[["study_name", "N_tot", "N_cases", "ancestry"]], on="study_name", how="left")
st["n_genes_susie_tested"] = st["study_name"].map(tr).fillna(0).astype(int)
st["n_genes_pp4_gt_0.5"] = st["study_name"].map(ps).fillna(0).astype(int)
st.to_csv(OUT / "recompute_1c_per_study.tsv", sep="\t", index=False)
res["per_study_sum_by_class"] = st.groupby("tier_label")[["n_genes_susie_tested", "n_genes_pp4_gt_0.5"]].sum().to_dict()

# positive controls: per-study rows
pcsym = ["PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "MBOAT7", "APOE", "MTARC1", "GPAM", "TRIB1", "SERPINA1"]
pc = m[m["gene"].isin(pcsym)][["gene", "gwas_name", "trait", "tier_label", "method", "PP.H4.abf", "PP.H4.susie", "PP.H3.susie", "n_cs_pairs", "n_cs_pairs_eligible", "eligibility_source"]]
pc.sort_values(["gene", "gwas_name"]).to_csv(OUT / "recompute_1c_positive_controls.tsv", sep="\t", index=False)
summ = {}
for s in pcsym:
    q = pc[pc["gene"] == s]
    summ[s] = dict(n_main_rows=len(q), n_susie_finite=int(q["PP.H4.susie"].notna().sum()),
                   max_susie=None if q["PP.H4.susie"].isna().all() else float(q["PP.H4.susie"].max()),
                   max_abf=None if q["PP.H4.abf"].isna().all() else float(q["PP.H4.abf"].max()),
                   max_abf_study=None if q["PP.H4.abf"].isna().all() else q.loc[q["PP.H4.abf"].idxmax(), "gwas_name"],
                   methods=q["method"].value_counts(dropna=False).to_dict())
    ens = m.loc[m["gene"] == s, "ensembl"].unique()
    summ[s]["ensembl"] = list(ens)
    summ[s]["class_in_1c"] = [cls.get(e, "not_in_436") for e in ens]
    summ[s]["gene_level_state"] = gl.loc[gl["ensembl"].isin(ens), "susie_state_t12"].tolist()
res["positive_controls"] = summ
cls.rename("class").to_frame().to_csv(OUT / "recompute_1c_gene_classes.tsv", sep="\t")
json.dump(res, open(OUT / "recompute_1c.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
