import pandas as pd
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
t=pd.read_csv(B+"/GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t"); MAIN=set(t.loc[t.placement=="main","study_name"])
c=pd.read_csv(B+"/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","gene","PP.H4.susie","top_snp"])
c=c[c.gwas_name.isin(MAIN)&c.gene.notna()&(c.gene!="")&c.top_snp.notna()&(c.top_snp!="")]
c["s4"]=pd.to_numeric(c["PP.H4.susie"],errors="coerce"); c=c[c.s4.notna()]
c["m"]=c.groupby("gene").s4.transform("max"); a=c[c.s4==c.m]
k=a.groupby("gene").apply(lambda d: d[["gwas_name","top_snp"]].drop_duplicates().shape[0])
tg=k[k>1].index.tolist(); print(a[a.gene.isin(tg)].to_string())
s=pd.read_csv("%s/out/fig2b_recheck_shown.tsv"%"/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F2abce",sep="\t")
print("tied genes shown:", sorted(set(tg)&set(s.gene)))
