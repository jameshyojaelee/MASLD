import pandas as pd, numpy as np
R='results/susie_coloc_r2_20260924T053244Z/'
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-05/'
t=pd.read_csv('config/gwas_trait_tier.tsv',sep='\t'); t12=t[t.tier.isin([1,2])]
m=pd.read_csv(R+'susie_coloc_all_gwas.csv',usecols=['gwas_name','ensembl','PP.H4.susie','PP.H4.abf'])
m=m[m.gwas_name.isin(t12.study_name)].merge(t12[['study_name','trait']],left_on='gwas_name',right_on='study_name')
c=pd.read_csv(S+'classes_rederived.tsv',sep='\t')
mm=m[m.ensembl.isin(c.ensembl)].merge(c[['ensembl','cls','locus']],on='ensembl')
ggt=mm[mm.trait=='GGT']
a=c[c.cls=='ALTAST_only'].copy()
gt=ggt[ggt.ensembl.isin(a.ensembl)]
susie_tested=gt[gt['PP.H4.susie'].notna()].groupby('ensembl')['PP.H4.susie'].max()
abf=gt.groupby('ensembl')['PP.H4.abf'].max()
a['ggt_susie_tested']=a.ensembl.isin(susie_tested.index)
a['ggt_susie_max']=a.ensembl.map(susie_tested)
a['ggt_abf_max']=a.ensembl.map(abf)
print('ALTAST_only genes:',len(a),'; with any GGT SuSiE PP.H4 value:',a.ggt_susie_tested.sum(),'; GGT ABF PP.H4 > 0.5:',(a.ggt_abf_max>0.5).sum(),'; GGT SuSiE max in (0.2,0.5]:',((a.ggt_susie_max>0.2)&(a.ggt_susie_max<=0.5)).sum())
pos=mm[(mm['PP.H4.susie']>0.5)&mm.trait.isin(['ALT','AST'])&mm.ensembl.isin(a.ensembl)]
src=pos.groupby('ensembl').gwas_name.agg(lambda s: 'MVP_only' if all(x.startswith('MVP') for x in s) else 'has_nonMVP')
print('ALTAST_only genes by source of positive ALT/AST coloc:', src.value_counts().to_dict())
g=c[c.cls=='GGT_only']; gg=ggt[ggt.ensembl.isin(g.ensembl)&(ggt['PP.H4.susie']>0.5)]
print('GGT_only genes by GGT study:', gg.groupby('ensembl').gwas_name.agg(lambda s:'|'.join(sorted(set(s)))).value_counts().to_dict())
alt=mm[mm.trait.isin(['ALT','AST'])&mm.ensembl.isin(g.ensembl)]
print('GGT_only genes: max ALT/AST SuSiE PP.H4 distribution:'); print(alt.groupby('ensembl')['PP.H4.susie'].max().describe())
print('GGT_only genes with ALT/AST ABF PP.H4 > 0.5:', (alt.groupby('ensembl')['PP.H4.abf'].max()>0.5).sum())
