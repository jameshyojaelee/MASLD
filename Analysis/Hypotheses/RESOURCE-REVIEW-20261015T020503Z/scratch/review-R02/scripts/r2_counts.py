import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
R2=B+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
tier=pd.read_csv(B+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
main=tier[tier.placement=='main']
print('main studies',len(main))
cols=['gwas_name','gene','ensembl','chr','PP.H4.abf','PP.H4.susie','method','top_snp','top_snp_PP','ancestry','ld_panel','n_palindromic_kept','n_palindromic_af_tested','n_palindromic_af_dropped','n_ambiguous','n_flip','n_match','n_snps']
d=pd.read_csv(R2,usecols=cols,low_memory=False)
print('r2 rows',len(d),'studies',d.gwas_name.nunique())
m=d[d.gwas_name.isin(main.study_name)].copy()
print('main rows',len(m),'main studies present',m.gwas_name.nunique())
for c in ['PP.H4.abf','PP.H4.susie']: m[c]=pd.to_numeric(m[c],errors='coerce')
g=m.groupby('ensembl').agg(s=('PP.H4.susie','max'),a=('PP.H4.abf','max'))
multi=g.index[g.s>0.5]; single=g.index[(g.a>0.5)&~(g.s>0.5)]
print('multi',len(multi),'single-only',len(single),'union',len(multi)+len(single))
# by symbol
gs=m.groupby('gene').agg(s=('PP.H4.susie','max'),a=('PP.H4.abf','max'))
print('by symbol multi',(gs.s>0.5).sum(),'single-only',((gs.a>0.5)&~(gs.s>0.5)).sum())
# bands
s=g.s[g.s>0.5]; print('multi bands >0.9',(s>0.9).sum(),'0.8-0.9',((s>0.8)&(s<=0.9)).sum(),'0.5-0.8',(s<=0.8).sum())
a=g.a[single]; print('single bands >0.9',(a>0.9).sum(),'0.8-0.9',((a>0.8)&(a<=0.9)).sum(),'0.5-0.8',(a<=0.8).sum())
# ld panel by ancestry
print(m.groupby(['ancestry','ld_panel']).gwas_name.nunique())
print(d.groupby(['ancestry','ld_panel']).gwas_name.nunique())
m.to_pickle('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R02/r2_main.pkl')
