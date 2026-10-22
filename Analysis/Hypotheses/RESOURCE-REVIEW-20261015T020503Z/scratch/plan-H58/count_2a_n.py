"""H58 plan: n for the Fig 2A support count. Counts strata, genes and rows only.
Does not intersect COLOC leads with fine-mapped windows."""
import pandas as pd
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
tier = pd.read_csv(ROOT + "GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
MAIN = set(tier.loc[tier.placement == "main", "study_name"])
print("MAIN strata", len(MAIN))
cs = pd.read_csv(ROOT + "GWAS/finemapping/results/credible_sets.csv", usecols=["study", "susie_cs", "locus"], dtype=str)
print("credible_sets.csv rows", len(cs))
cs = cs[cs.study.isin(MAIN)]
print("MAIN strata with any credible-set row", cs.study.nunique(), "missing:", sorted(MAIN - set(cs.study)))
co = pd.read_csv(ROOT + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
                 usecols=["gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp"])
co = co[co.gwas_name.isin(MAIN) & (co.gene.fillna("").str.strip() != "")]
gmax = co.groupby("gene")["PP.H4.susie"].max()
print("COLOC r2 MAIN strata", co.gwas_name.nunique(), "genes tested", co.gene.nunique(),
      "multi-signal genes (max PP.H4.susie>0.5)", int((gmax > 0.5).sum()),
      "positive gene-study rows", int((co["PP.H4.susie"] > 0.5).sum()))
print("top_snp example", co.top_snp.dropna().head(3).tolist())
