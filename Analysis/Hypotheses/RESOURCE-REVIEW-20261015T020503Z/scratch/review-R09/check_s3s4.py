import pandas as pd
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
p3a=pd.read_csv(root+'GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z/tables/mpra_dav_calls_S2_S5.tsv',sep='\t')
print(p3a.groupby(['table','context']).size())
src={}
for cl,cond in [('HepG2','control'),('HepG2','PAOA'),('LX2','control'),('LX2','TGFb')]:
    d=pd.read_csv(root+f'GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction/{cl}.{cond}.official_comparison.tsv.gz',sep='\t')
    src[f'{cl}_{cond}']=d
    print(cl,cond,len(d),'official_fdr<0.05 & active:',int(((d.official_fdr<0.05)&(d.active_variant==True)).sum()),'official_fdr<0.05:',int((d.official_fdr<0.05).sum()))
print()
for (t,c),g in p3a.groupby(['table','context']):
    for k,d in src.items():
        m=g.merge(d,left_on='coordination',right_on='element_id')
        if len(m)==0: continue
        r=m['log2(Fold change)'].corr(m['official_log2FC'])
        exact=(abs(m['log2(Fold change)']-m['official_log2FC'])<1e-9).mean()
        print(t,c,'->',k,'n_join',len(m),'of',len(g),'pearson',round(r,4),'exact_match_frac',round(exact,4))
