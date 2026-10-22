# Light read: cluster-level coding fraction vs presence of a shared liver eQTL (COLOC r2), tier1+2 main GWAS
import pandas as pd, numpy as np
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
cs=pd.read_csv(R+'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/trait_directness/corrected_locus_coordinates.tsv',sep='\t')
print('credible sets',len(cs), cs.trait_scope.value_counts().to_dict(), 'gate passed', cs.pip_sum_gate_passed.sum())
cs=cs[cs.pip_sum_gate_passed==True].copy()
cs=cs[cs.trait_scope.isin(['tier1_direct_masld_pdff','tier2_liver_enzyme'])]
# single-linkage clusters of sentinel positions within 1 Mb, per chromosome (pooled scopes)
cs=cs.sort_values(['chrom','pos'])
cid=[];prev=None;k=0
for c,p in zip(cs.chrom,cs.pos):
    if prev is None or c!=prev[0] or p-prev[1]>1_000_000: k+=1
    cid.append(k); prev=(c,p)
cs['cluster']=cid
cs['coding']=cs.protein_altering_pip_mass+cs.canonical_splice_pip_mass
cl=cs.groupby('cluster').agg(chrom=('chrom','first'),lo=('member_min_position','min'),hi=('member_max_position','max'),
   coding_mean=('coding','mean'),n_cs=('credible_set_uid','size'),
   scopes=('trait_scope',lambda s:'|'.join(sorted(set(s)))),studies=('study',lambda s:set(s)))
print('clusters',len(cl), cl.scopes.value_counts().to_dict())
cl['coding_dominant']=cl.coding_mean>=0.5
cl['coding_any']=cl.coding_mean>=0.2
print('coding_mean>=0.5:',int(cl.coding_dominant.sum()),' >=0.2:',int(cl.coding_any.sum()))
# COLOC
cols=['gwas_name','gene','ensembl','chr','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf','PP.H4.susie','top_snp']
d=pd.read_csv(R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',usecols=cols)
d=d[d['PP.H4.abf'].fillna(0).gt(0.5)|d['PP.H4.susie'].fillna(0).gt(0.5)].copy()
d['tpos']=d.top_snp.str.split(':').str[1].astype(float)
d['tchr']=d.top_snp.str.split(':').str[0]
res=[]
for i,r in cl.iterrows():
    m=d[(d.tchr==str(r.chrom))&(d.tpos>=r.lo-250_000)&(d.tpos<=r.hi+250_000)&(d.gwas_name.isin(r.studies))]
    res.append((i,m.ensembl.nunique(), int(m['PP.H4.susie'].fillna(0).gt(0.5).any())))
x=pd.DataFrame(res,columns=['cluster','n_h4_genes','susie_h4']).set_index('cluster')
cl=cl.join(x)
cl['shared_eqtl']=cl.n_h4_genes>0
print(pd.crosstab(cl.coding_dominant,cl.shared_eqtl,margins=True))
print(pd.crosstab(cl.coding_any,cl.shared_eqtl,margins=True))
cl.drop(columns='studies').to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L03/c1_cluster_table.tsv',sep='\t')
# where are PNPLA3 (22:44.32M), MTARC1 (1:220.97M), TM6SF2 (19:19.38M), HSD17B13 (4:88.2M), GCKR (2:27.73M)
for nm,c,p in [('PNPLA3',22,44324727),('MTARC1',1,220970028),('TM6SF2',19,19379549),('HSD17B13',4,88231392),('GCKR',2,27730940),('APOE',19,45411941)]:
    q=cl[(cl.chrom==c)&(cl.lo-1_000_000<=p)&(cl.hi+1_000_000>=p)]
    print(nm, q[['lo','hi','coding_mean','n_cs','scopes','n_h4_genes','susie_h4']].round(3).to_dict('records'))
