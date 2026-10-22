import pandas as pd, numpy as np
from collections import defaultdict
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
tier=pd.read_csv(root+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
MAIN=set(tier.loc[tier.placement=='main','study_name'])
cs=pd.read_csv(root+'GWAS/finemapping/results/credible_sets.csv',usecols=['chromosome','position','study','locus'],dtype={'locus':str,'chromosome':str})
cs=cs[cs.study.isin(MAIN)]
cs['position']=cs.position.astype(float).astype(int)
win=cs.groupby(['study','locus']).agg(chrom=('chromosome','first'),mn=('position','min'),mx=('position','max')).reset_index()
print('windows',len(win))
co=pd.read_csv(root+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',usecols=['gwas_name','gene','PP.H4.susie','PP.H4.abf','top_snp','method'])
co=co[co.gwas_name.isin(MAIN)]
co=co[co.gene.fillna('').str.strip()!='']
pos=co[co['PP.H4.susie']>0.5].copy()
gmax=co.groupby('gene')['PP.H4.susie'].max()
multi=set(gmax[gmax>0.5].index)
print('multi-signal genes',len(multi), 'positive gene-study rows',len(pos))
pos['chrom']=pos.top_snp.str.split(':').str[0]
pos['bp']=pos.top_snp.str.split(':').str[1].astype(float)
print('top_snp missing among positives',pos.top_snp.isna().sum())
wd=defaultdict(list)
for r in win.itertuples(): wd[(r.study,str(r.chrom))].append((r.mn,r.mx))
wa=defaultdict(list)
for r in win.itertuples(): wa[str(r.chrom)].append((r.mn,r.mx))
def inside(lst,bp,pad): return any(a-pad<=bp<=b+pad for a,b in lst)
for pad in (0,250000):
    pos[f'same_{pad}']=[inside(wd[(s,c)],bp,pad) for s,c,bp in zip(pos.gwas_name,pos.chrom,pos.bp)]
    pos[f'any_{pad}']=[inside(wa[c],bp,pad) for c,bp in zip(pos.chrom,pos.bp)]
    g_same=pos.groupby('gene')[f'same_{pad}'].any(); g_any=pos.groupby('gene')[f'any_{pad}'].any()
    print(f'pad={pad}: genes with a positive lead inside same-study window: {g_same.sum()}; inside any window: {g_any.sum()} of {len(g_same)}')
fm_studies=set(win.study)
mvp_eur={'MVP_ALT_EUR','MVP_AST_EUR','MVP_NAFLD_EUR'}
ukbb={'UKBB_ALT','UKBB_AST','UKBB_GGT'}
gs=pos.groupby('gene').gwas_name.apply(set)
only_mvp=sum(1 for s in gs if s<=mvp_eur)
only_mvp_ukbb=sum(1 for s in gs if s<=mvp_eur|ukbb)
print('supported only by MVP EUR (ALT/AST/NAFLD):',only_mvp)
print('supported only by MVP EUR + UKBB ALT/AST/GGT:',only_mvp_ukbb)
print('genes with any positive study lacking fine-mapping (MVP EUR):',sum(1 for s in gs if s & mvp_eur))
print('positive rows by study:'); print(pos.gwas_name.value_counts().to_string())
# driving study per gene (max susie)
idx=co[co.gene.isin(multi)].sort_values(['gene','PP.H4.susie','gwas_name'],ascending=[True,False,True]).groupby('gene').head(1)
print('driving study of 436 genes:'); print(idx.gwas_name.value_counts().to_string())
pos.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-03-coverage/positive_rows_window_overlap.tsv',sep='\t',index=False)
