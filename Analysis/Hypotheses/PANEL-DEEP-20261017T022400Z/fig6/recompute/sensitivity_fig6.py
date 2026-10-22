#!/usr/bin/env python3
"""Display-choice sensitivity for Figure 6B/6C/6D: recompute the matrix and class counts
under reasonable alternative rules. Inputs are the same upstream files as recompute_fig6.py."""
import os, re, json
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute")
MASTER = os.path.join(ROOT, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")
TIERS = os.path.join(ROOT, "GWAS/finemapping/config/gwas_trait_tier.tsv")
H87 = os.path.join(ROOT, "Analysis/Hypotheses/H87-20261015T130324Z/run-20261015T131921Z/h87_coding_member_gene_map.tsv")
UNTESTABLE = "susie_untestable_insufficient_shared_posterior"

c = pd.read_csv(os.path.join(OUT, "catalog_calls_recomputed.tsv.gz"), sep="\t")
tiers = pd.read_csv(TIERS, sep="\t")
tier = dict(zip(tiers.study_name, tiers.tier))
m = pd.read_csv(MASTER, usecols=["gwas_name", "ensembl", "PP.H4.abf", "PP.H4.susie", "method"], low_memory=False)
m = m[m.ensembl.astype(str).str.match(r"^ENSG\d{11}$")]
m["tier"] = m.gwas_name.map(tier)
for col in ("PP.H4.abf", "PP.H4.susie"):
    m[col] = pd.to_numeric(m[col], errors="coerce")


def gen_state(sub):
    g = sub.groupby("ensembl").agg(s=("PP.H4.susie", "max"), a=("PP.H4.abf", "max"),
                                   u=("method", lambda x: int((x == UNTESTABLE).sum())))
    return g


t12 = gen_state(m[m.tier.isin([1, 2])])
t1 = gen_state(m[m.tier == 1])
n_pos_studies = m[m.tier.isin([1, 2]) & (m["PP.H4.susie"] > 0.5)].groupby("ensembl").gwas_name.nunique()

h = pd.read_csv(H87, sep="\t", low_memory=False)
cod = h[(h.canonical_PA == True) & (h.normalized_susie_pip >= 0.5) & h.tier.isin([1, 2])].ensembl_id.dropna().unique()
cod_all_tiers = h[(h.canonical_PA == True) & (h.normalized_susie_pip >= 0.5)].ensembl_id.dropna().unique()


def classify(gc, rc):
    cls = np.where((gc == "untestable") | (rc == "untestable"), "untested",
          np.where((gc == "supported") & (rc == "supported"), "concordant",
          np.where(gc == "supported", "genetically_anchored",
          np.where(rc == "supported", "established_state_associated", "unresolved"))))
    return cls


def gen_call_from(g, rule):
    key = c.gen_key
    s = g.reindex(key).s.values
    a = g.reindex(key).a.values
    u = g.reindex(key).u.values
    present = np.isin(key, g.index)
    if rule == "susie>0.5":
        pos = s > 0.5
    elif rule == "susie>0.8":
        pos = s > 0.8
    elif rule == "susie_or_abf>0.5":
        pos = (s > 0.5) | (a > 0.5)
    elif rule == "susie>0.5_or_abf>0.8_when_no_susie":
        pos = (s > 0.5) | (np.isnan(s) & (a > 0.8))
    pos = np.nan_to_num(pos, nan=0).astype(bool)
    unt = (~present) | ((u > 0) & ~pos)
    return np.where(pos, "supported", np.where(unt, "untestable", "indeterminate"))


base_rna = c.rna_call.values
q = c.rna_q.values; fc = c.rna_log2fc.values
scen = {}
scen["0_adopted_rule"] = (c.gen_call.values, base_rna)
scen["A_rna_q<0.05_no_fc_floor"] = (c.gen_call.values, np.where(c.rna_call == "untestable", "untestable", np.where(q < 0.05, "supported", "indeterminate")))
scen["B_rna_q<0.05_abs_log2fc>1"] = (c.gen_call.values, np.where(c.rna_call == "untestable", "untestable", np.where((q < 0.05) & (np.abs(fc) > 1), "supported", "indeterminate")))
scen["C_gen_susie>0.8"] = (gen_call_from(t12, "susie>0.8"), base_rna)
scen["D_gen_susie_or_abf>0.5"] = (gen_call_from(t12, "susie_or_abf>0.5"), base_rna)
scen["E_gen_abf>0.8_where_no_susie"] = (gen_call_from(t12, "susie>0.5_or_abf>0.8_when_no_susie"), base_rna)
g_cod = c.gen_call.values.copy()
g_cod = np.where(np.isin(c.ensembl_id, cod) & (g_cod != "supported"), "supported", g_cod)
scen["F_gen_plus_coding_finemap_pip>=0.5_tier12"] = (g_cod, base_rna)
g7 = np.where(c.susie_state.eq("untestable"), "indeterminate", c.gen_call.values)
scen["G_seven_susie_untestable_as_indeterminate"] = (g7, base_rna)
scen["H_gen_tier1_direct_masld_only"] = (gen_call_from(t1, "susie>0.5"), base_rna)

rows, wrows = [], []
for k, (gc, rc) in scen.items():
    cl = classify(np.asarray(gc), np.asarray(rc))
    vc = pd.Series(cl).value_counts()
    mt = pd.crosstab(pd.Series(gc, name="gen"), pd.Series(rc, name="rna"))
    r = {"scenario": k}
    for g_ in ("supported", "indeterminate", "untestable"):
        for r_ in ("supported", "indeterminate", "untestable"):
            r[f"gen_{g_}|rna_{r_}"] = int(mt.loc[g_, r_]) if (g_ in mt.index and r_ in mt.columns) else 0
    for cc in ("untested", "unresolved", "established_state_associated", "genetically_anchored", "concordant"):
        r["class_" + cc] = int(vc.get(cc, 0))
    rows.append(r)
    for gname in ("GLP1R", "HKDC1", "MTARC1", "THRB", "PNPLA3", "TM6SF2"):
        i = np.where(c.symbol.values == gname)[0][0]
        wrows.append({"scenario": k, "gene": gname, "gen_call": gc[i], "rna_call": rc[i], "class": cl[i]})
S = pd.DataFrame(rows)
S.to_csv(os.path.join(OUT, "sensitivity_matrix_and_classes.tsv"), sep="\t", index=False)
W = pd.DataFrame(wrows)
W.to_csv(os.path.join(OUT, "sensitivity_worked_entries.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(S.set_index("scenario").T.to_string())
print(W.pivot(index="gene", columns="scenario", values="class").to_string())
print("coding pip>=0.5 genes (tier1/2):", len(set(cod) & set(c.ensembl_id)), "all tiers:", len(set(cod_all_tiers) & set(c.ensembl_id)))
pos = c[c.gen_call == "supported"].gen_key.map(n_pos_studies)
print("genetics-supported genes by number of Tier-1/2 studies with SuSiE PP.H4>0.5:")
print(pos.value_counts().sort_index().to_string())
for gname in ("HKDC1", "THRB"):
    k = c.loc[c.symbol == gname, "gen_key"].iloc[0]
    print(gname, "studies SuSiE>0.5:", sorted(m[(m.ensembl == k) & m.tier.isin([1, 2]) & (m["PP.H4.susie"] > 0.5)].gwas_name))
json.dump({"n_pos_studies_distribution": pos.value_counts().sort_index().to_dict()}, open(os.path.join(OUT, "sensitivity_extra.json"), "w"), default=str)
