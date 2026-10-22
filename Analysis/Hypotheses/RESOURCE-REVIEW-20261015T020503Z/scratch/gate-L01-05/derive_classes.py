import pandas as pd
R='results/susie_coloc_r2_20260924T053244Z/'
t=pd.read_csv('config/gwas_trait_tier.tsv',sep='\t')
t12=t[t.tier.isin([1,2])]
m=pd.read_csv(R+'susie_coloc_all_gwas.csv',usecols=['gwas_name','gene','ensembl','chr','PP.H4.susie'])
m=m[m.gwas_name.isin(t12.study_name)]
g12=pd.read_csv(R+'gene_level_coloc_tier12.csv')
pos=set(g12.loc[g12.susie_state_t12=='positive','ensembl'])
print('positive genes', len(pos))
mp=m[(m['PP.H4.susie']>0.5)&m.ensembl.isin(pos)].merge(t12[['study_name','trait','tier']],left_on='gwas_name',right_on='study_name')
def cls(s):
    tr=set(s)
    if tr & {'NAFLD','PDFF','NASH'}: return 'direct_any'
    g='GGT' in tr; a=bool(tr & {'ALT','AST'})
    return 'GGT_plus_ALTAST' if (g and a) else ('GGT_only' if g else 'ALTAST_only')
c=mp.groupby('ensembl').trait.agg(cls).rename('cls_rederived')
print(c.value_counts()); print('genes classed', c.size)
L=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L01/gene_trait_class_locus.tsv',sep='\t')
x=L.merge(c,left_on='ensembl',right_index=True,how='left')
print(pd.crosstab(x.cls,x.cls_rederived.fillna('NA')))
# GGT studies count
print(t12.groupby('trait').size())
# H10 loci consistency
h=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv',sep='\t')
rows=[(lid,e) for lid,es in zip(h.locus_id,h.ensembl) for e in es.split(',')]
hl=pd.DataFrame(rows,columns=['locus_h10','ensembl'])
y=L.merge(hl,on='ensembl',how='left')
print('locus id match with H10:', (y.locus==y.locus_h10).sum(), 'of', len(y))
print('MHC loci:', h.loc[h.mhc==True,'locus_id'].tolist())
x.to_csv('%s/classes_rederived.tsv'%'/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-05',sep='\t',index=False)
