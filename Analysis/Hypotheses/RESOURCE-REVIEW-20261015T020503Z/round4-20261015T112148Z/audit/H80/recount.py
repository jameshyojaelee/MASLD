# Audit H80: independent recount of every n stated in PRESPEC_H80.md. Reads only; no simulation, no lead z.
import glob, os, resource, time
import numpy as np, pandas as pd
t0 = time.time()
FM = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
R2 = f"{FM}/results/susie_coloc_r2_20260924T053244Z"
tier = pd.read_csv(f"{FM}/config/gwas_trait_tier.tsv", sep="\t")
reg = pd.read_csv(f"{FM}/config/gwas_registry.tsv", sep="\t")
t12 = set(tier.loc[tier.tier.isin([1, 2]), "study_name"])
print("A tier-1/2 studies:", len(t12), "| in registry:", len(t12 & set(reg.study_name)))
m = pd.read_csv(f"{R2}/susie_coloc_all_gwas.csv", low_memory=False,
                usecols=["gwas_name", "ensembl", "chr", "PP.H3.susie", "PP.H4.susie", "method", "ancestry", "ld_panel", "eligibility_source"])
m = m[m.gwas_name.isin(t12)]
print("B tier-1/2 rows:", len(m), "| unique (gwas,ensembl):", len(m[["gwas_name", "ensembl"]].drop_duplicates()),
      "| rows with empty ensembl:", m.ensembl.isna().sum())
print("C method:", m.method.value_counts().to_dict())
s = m[np.isfinite(m["PP.H4.susie"])].merge(reg[["study_name", "ancestry", "trait_type"]].rename(columns={"study_name": "gwas_name", "ancestry": "reg_anc"}), on="gwas_name")
print("D finite SuSiE:", len(s), "| PP.H4>0.5:", (s["PP.H4.susie"] > 0.5).sum(), "| >0.9:", (s["PP.H4.susie"] > 0.9).sum(),
      "| PP.H3>0.8:", (s["PP.H3.susie"] > 0.8).sum(), "| method of finite rows:", s.method.value_counts().to_dict())
print("E ancestry (master):", s.ancestry.value_counts().to_dict(), "| registry:", s.reg_anc.value_counts().to_dict(),
      "| master==registry:", (s.ancestry == s.reg_anc).all())
print("F ld_panel label:", s.ld_panel.value_counts(dropna=False).to_dict())
w = s[["ensembl", "reg_anc"]].drop_duplicates()
print("G S-eval windows:", len(w), "genes:", s.ensembl.nunique(), w.reg_anc.value_counts().to_dict())
p9 = s[s["PP.H4.susie"] > 0.9]; w9 = p9[["ensembl", "reg_anc"]].drop_duplicates()
print("H S-pos pairs:", len(p9), "genes:", p9.ensembl.nunique(), "windows:", len(w9),
      "| pairs by anc:", p9.reg_anc.value_counts().to_dict(), "| windows by anc:", w9.reg_anc.value_counts().to_dict())
le = s[s["PP.H4.susie"] <= 0.5]; wle = le[["ensembl", "reg_anc"]].drop_duplicates()
print("I PP.H4<=0.5 pairs:", len(le), "windows:", len(wle), wle.reg_anc.value_counts().to_dict())
print("J trait type of SuSiE pairs:", s.trait_type.value_counts().to_dict())
print("K eligibility_source x (PP.H4>0.5):\n", pd.crosstab(s.eligibility_source.fillna("NA"), s["PP.H4.susie"] > 0.5).to_string())
print("K2 carried negatives by ancestry EUR/non-EUR:",
      (le.eligibility_source.eq("not_rerun_negative") & le.reg_anc.eq("EUR")).sum(),
      (le.eligibility_source.eq("not_rerun_negative") & le.reg_anc.ne("EUR")).sum())
rr = s[s.eligibility_source.eq("rerun_shared_posterior_check")]
print("K3 rerun rows EUR:", (rr.reg_anc == "EUR").sum(), "| EUR rerun PP.H4>0.5:", ((rr.reg_anc == "EUR") & (rr["PP.H4.susie"] > 0.5)).sum(),
      "| EUR rerun PP.H4<=0.5:", ((rr.reg_anc == "EUR") & (rr["PP.H4.susie"] <= 0.5)).sum())
g = pd.read_csv(f"{R2}/tier12/gene_level_coloc_tier12.csv")
print("L gene states:", g.susie_state_t12.value_counts().to_dict(), "total", len(g))
fs = sorted(glob.glob(f"{FM}/results/eqtl_susie_polyfun/chr*_summary.tsv"))
e = pd.concat([pd.read_csv(f, sep="\t") for f in fs])
print("M summary rows:", len(e), "| unique ensg:", e.ensg.nunique(), "| ensg with >1 row:", (e.ensg.value_counts() > 1).sum())
rds = {os.path.basename(p).replace("_susie.rds", "") for p in glob.glob(f"{FM}/results/eqtl_susie_polyfun/chr*/ENSG*_susie.rds")}
print("M2 rds files:", len(rds), "| status ok genes (last row):", (e.drop_duplicates('ensg', keep='last').status == 'ok').sum())
for keep in ["last", "first"]:
    el = e.drop_duplicates("ensg", keep=keep).copy()
    el["fit"] = np.select([el.status.ne("ok"), el.converged.astype(str).ne("True"), el.n_cs > 0],
                          ["no_fit_" + el.status.astype(str), "not_converged", "cs_ge1"], "cs0")
    x = g.merge(el[["ensg", "fit", "n_snps", "status"]], left_on="ensembl", right_on="ensg", how="left")
    x["fit"] = x.fit.fillna("no_summary_row")
    print(f"N crosstab keep={keep}:\n", pd.crosstab(x.susie_state_t12, x.fit, margins=True).to_string())
    if keep == "last":
        ok = x[x.status.eq("ok")]
        print("O window SNPs, status ok genes in gene table: n", len(ok), "median", ok.n_snps.median(),
              "IQR", ok.n_snps.quantile(.25), ok.n_snps.quantile(.75), "max", ok.n_snps.max())
        anyfit = x[x.n_snps.notna()]
        print("O2 same over all rows with n_snps (writer's rule): n", len(anyfit), "median", anyfit.n_snps.median(), "max", anyfit.n_snps.max())
        sp = s.merge(el[["ensg", "n_snps"]], left_on="ensembl", right_on="ensg", how="left")
        print("O3 SuSiE-pair eQTL p: median", sp.n_snps.median(), "p90", sp.n_snps.quantile(.9), "max", sp.n_snps.max(), "missing", sp.n_snps.isna().sum())
        print("Q eQTL CS genes in gene table (comparison set, keep=last):", int((x.fit == "cs_ge1").sum()))
        nores = set(x.loc[x.susie_state_t12.eq("no_susie_result") & x.fit.eq("cs0"), "ensembl"])
        notconv = x.loc[x.fit.eq("not_converged") & x.susie_state_t12.isin(["positive", "tested_le_0.5"]), ["ensembl", "susie_state_t12"]]
        notconv.to_csv(os.path.join(os.path.dirname(__file__), "notconverged_summary_genes.tsv"), sep="\t", index=False)
        cs0nores = x.loc[x.susie_state_t12.eq("no_susie_result") & x.fit.eq("cs0"), ["ensembl", "n_snps"]]
        print("O4 cs0 no-result genes n_snps < 50:", (cs0nores.n_snps < 50).sum(), "| min", cs0nores.n_snps.min(),
              "| median", cs0nores.n_snps.median(), "| max", cs0nores.n_snps.max())
mm = m[m.ensembl.isin(nores)].drop(columns="ancestry").merge(reg[["study_name", "ancestry"]], left_on="gwas_name", right_on="study_name")
print("P S-nores genes with a tier-1/2 pair:", mm.ensembl.nunique(), "| pairs:", len(mm), mm.method.value_counts().to_dict(),
      "| EUR windows:", mm.loc[mm.ancestry.eq("EUR"), "ensembl"].nunique())
print("rss_GB", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2), "wall_s", round(time.time() - t0, 1))
