import pandas as pd
W='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/'
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
tiers=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
t12=set(tiers[tiers.tier.isin([1,2])].study_name)
a=pd.read_csv(C+'susie_coloc_all_gwas.csv',usecols=['gwas_name','ensembl','PP.H3.susie','PP.H4.susie'])
a=a[a.gwas_name.isin(t12)&a.ensembl.notna()&a['PP.H4.susie'].notna()]
g=a.groupby('ensembl').agg(h3=('PP.H3.susie','max'),h4=('PP.H4.susie','max'))
print('genes with finite tier1/2 SuSiE',len(g))
p=g[(g.h3>0.8)&(g.h4<=0.5)]; print('PP.H3>0.8 in >=1 study and best PP.H4<=0.5:',len(p))
m=pd.read_csv(W+'catalog_gene_calls_with_coloc.tsv',sep='\t')
mi=m[m.call_state=='indeterminate']
print('of those, Catalog genetics indeterminate (direct ID):',mi.ensembl_id.isin(p.index).sum())
print('genetics-indeterminate split: PP.H3>0.8 pattern',mi.ensembl_id.isin(p.index).sum(),
      ' finite SuSiE other',(mi.coloc_best_susie_pp4.notna()&~mi.ensembl_id.isin(p.index)).sum(),
      ' no SuSiE',mi.coloc_best_susie_pp4.isna().sum())
