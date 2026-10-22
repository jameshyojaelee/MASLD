import pandas as pd, numpy as np
R='results/susie_coloc_r2_20260924T053244Z/'
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-05/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/'
t=pd.read_csv('config/gwas_trait_tier.tsv',sep='\t'); t12=t[t.tier.isin([1,2])]
matched=[s for s in t12.study_name if s.split('_')[0] in ('UKBB','BBJ','PanUKBB')]
print('matched-ascertainment studies:',matched)
m=pd.read_csv(R+'susie_coloc_all_gwas.csv',usecols=['gwas_name','ensembl','PP.H4.susie'])
c=pd.read_csv(S+'classes_rederived.tsv',sep='\t')
m=m[m.gwas_name.isin(t12.study_name)&m.ensembl.isin(c.ensembl)].merge(t12[['study_name','trait']],left_on='gwas_name',right_on='study_name')
p=m[m['PP.H4.susie']>0.5]
direct=set(p.loc[p.trait.isin(['NAFLD','PDFF','NASH']),'ensembl'])
pm=p[p.gwas_name.isin(matched)]
def k(tr):
    tr=set(tr); g='GGT' in tr; a=bool(tr&{'ALT','AST'})
    return 'GGT_plus_ALTAST' if g and a else ('GGT_only' if g else 'ALTAST_only')
cm=pm.groupby('ensembl').trait.agg(k).rename('cls_matched')
c['cls_matched']=c.ensembl.map(cm).fillna('not_positive_in_matched_studies')
c.loc[c.ensembl.isin(direct),'cls_matched']='direct_any'
print(pd.crosstab(c.cls,c.cls_matched))
f3=pd.read_csv(B+'stage_extension_all_gene_results.tsv',sep='\t',usecols=['gene_id_base','axis','contrast'])
ids=set(f3.gene_id_base)
x=c[c.ensembl.isin(ids)&~c.locus.isin(['L090','L091','L092'])&c.cls_matched.isin(['GGT_only','ALTAST_only','GGT_plus_ALTAST','direct_any'])]
lc=x.groupby('locus').cls_matched.agg(lambda s:'|'.join(sorted(set(s))))
print('pure loci under matched ascertainment (non-MHC, in stage table):'); print(lc[~lc.str.contains(r'\|')].value_counts())
n1=(lc=='GGT_only').sum(); n2=(lc=='ALTAST_only').sum()
sd=0.2332
print('MDE log2FC at SD %.4f: %.4f'%(sd, 2.80158*sd*np.sqrt(1/n1+1/n2)))
