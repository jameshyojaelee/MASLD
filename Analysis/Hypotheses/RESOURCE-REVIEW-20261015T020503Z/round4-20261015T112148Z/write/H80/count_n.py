# H80 prespec counts. Set sizes and sampling frames only; no simulation, no recall.
import glob, os
import numpy as np, pandas as pd

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
R2 = f"{FM}/results/susie_coloc_r2_20260924T053244Z"
OUT = os.path.dirname(os.path.abspath(__file__))

tier = pd.read_csv(f"{FM}/config/gwas_trait_tier.tsv", sep="\t")
t12 = set(tier.loc[tier.tier.isin([1, 2]), "study_name"])
reg = pd.read_csv(f"{FM}/config/gwas_registry.tsv", sep="\t")
print("tier-1/2 studies", len(t12), "in registry", len(t12 & set(reg.study_name)))

cols = ["gwas_name", "ensembl", "chr", "PP.H3.susie", "PP.H4.susie", "method",
        "ancestry", "ld_panel", "n_snps", "n_cs_pairs", "n_cs_pairs_eligible", "eligibility_source"]
m = pd.read_csv(f"{R2}/susie_coloc_all_gwas.csv", usecols=cols)
m = m[m.gwas_name.isin(t12)].copy()
print("tier-1/2 pairs", len(m), "unique (gwas,ensembl)", len(m[["gwas_name", "ensembl"]].drop_duplicates()))
s = m[m["PP.H4.susie"].notna()]
print("SuSiE result", len(s), "PP.H4>0.5", (s["PP.H4.susie"] > 0.5).sum(),
      "PP.H4>0.9", (s["PP.H4.susie"] > 0.9).sum(), "PP.H3>0.8", (s["PP.H3.susie"] > 0.8).sum())
print("method counts, tier-1/2:\n", m.method.value_counts(dropna=False).to_string())
print("SuSiE pairs by registry ld_panel label:\n", s.ld_panel.value_counts(dropna=False).to_string())
print("SuSiE pairs by ancestry:\n", s.ancestry.value_counts(dropna=False).to_string())
s = s.merge(reg[["study_name", "ancestry", "trait_type", "N_tot", "N_cases"]].rename(
    columns={"study_name": "gwas_name", "ancestry": "reg_anc"}), on="gwas_name", how="left")
print("SuSiE pairs by registry ancestry:\n", s.reg_anc.value_counts(dropna=False).to_string())
# Windows: the eQTL side is one fit per gene (EUR polyfun); the GWAS side LD is per ancestry.
win = s[["ensembl", "reg_anc"]].drop_duplicates()
print("unique gene x GWAS-ancestry windows among SuSiE pairs", len(win), "unique genes", s.ensembl.nunique())
print(win.reg_anc.value_counts().to_string())
pos9 = s[s["PP.H4.susie"] > 0.9]
print("PP.H4>0.9 pairs", len(pos9), "genes", pos9.ensembl.nunique(),
      "gene x ancestry", len(pos9[["ensembl", "reg_anc"]].drop_duplicates()))
print(pos9.reg_anc.value_counts().to_string())
tested = s[s["PP.H4.susie"] <= 0.5]
print("PP.H4<=0.5 pairs", len(tested), "gene x ancestry windows",
      len(tested[["ensembl", "reg_anc"]].drop_duplicates()))
print(tested[["ensembl", "reg_anc"]].drop_duplicates().reg_anc.value_counts().to_string())
print("binary SuSiE pairs", (s.trait_type == "binary").sum(), "quant", (s.trait_type != "binary").sum())

# Gene states and eQTL fits (production LD-pinned polyfun fits vs the older eqtl_susie set)
g = pd.read_csv(f"{R2}/tier12/gene_level_coloc_tier12.csv")
print("gene states:\n", g.susie_state_t12.value_counts().to_string())
for d in ["eqtl_susie_polyfun", "eqtl_susie"]:
    fs = sorted(glob.glob(f"{FM}/results/{d}/chr*_summary.tsv"))
    e = pd.concat([pd.read_csv(f, sep="\t") for f in fs])
    e = e.drop_duplicates("ensg", keep="last")
    e["fit"] = np.where(e.status.ne("ok"), "no_fit_" + e.status.astype(str),
                        np.where(~e.converged.astype(str).eq("True"), "not_converged",
                                 np.where(e.n_cs > 0, "cs_ge1", "cs0")))
    x = g.merge(e[["ensg", "fit", "n_snps"]], left_on="ensembl", right_on="ensg", how="left")
    x["fit"] = x.fit.fillna("no_summary_row")
    print(f"== {d}: summary rows {len(e)}")
    print(pd.crosstab(x.susie_state_t12, x.fit, margins=True).to_string())
    ok = x[x.n_snps.notna()]
    print("window SNPs among genes with a fit: median", ok.n_snps.median(),
          "IQR", ok.n_snps.quantile(0.25), ok.n_snps.quantile(0.75), "max", ok.n_snps.max())
    if d == "eqtl_susie_polyfun":
        x.to_csv(f"{OUT}/gene_state_by_polyfun_fit.tsv", sep="\t", index=False)
        w = s.merge(e[["ensg", "n_snps"]].rename(columns={"n_snps": "eqtl_p"}),
                    left_on="ensembl", right_on="ensg", how="left")
        print("SuSiE-pair window eQTL p: median", w.eqtl_p.median(), "p90", w.eqtl_p.quantile(0.9),
              "max", w.eqtl_p.max(), "missing", w.eqtl_p.isna().sum())
s.to_csv(f"{OUT}/susie_pairs_tier12.tsv.gz", sep="\t", index=False)
