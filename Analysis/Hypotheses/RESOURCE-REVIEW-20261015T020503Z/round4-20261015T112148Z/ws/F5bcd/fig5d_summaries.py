"""R11-02 and R06-08 summaries on the 7080 output. Gene is the unit for tests; pairs are summarized for reference.
States: Both sources / One source as in the panel; 'partial' split at a 0.05 posterior-mass floor
(max of the single-source masses); sensitivity at 0.01. Gene-level state = best state over the gene's pairs
within a trait class (sensitivity: worst state). Background test: per gene x lineage x trait class,
mean over pairs of (observed any-source mass - locus-matched expected); Wilcoxon signed-rank over genes;
BH over the 10 lineage x trait-class cells."""
import sys, numpy as np, pandas as pd
from scipy import stats
OUT=sys.argv[1]
r=pd.read_csv('fig5d_locus_background_rows.tsv',sep='\t')
r['maxsrc']=r[['gse244832_accessible_mass','gse281367_accessible_mass']].max(axis=1)
def st(row,floor):
    s=row.evidence_state
    if s in('replicated_accessible','source_dependent','indeterminate'): return s
    return 'partial_ge_floor' if row.maxsrc>=floor else 'trace_lt_floor'
ORDER=['replicated_accessible','source_dependent','partial_ge_floor','trace_lt_floor','indeterminate']
for fl,tag in [(0.05,'f05'),(0.01,'f01')]:
    r['state_'+tag]=r.apply(st,axis=1,floor=fl)
r['rank_f05']=r.state_f05.map({s:i for i,s in enumerate(ORDER)})
r['rank_f01']=r.state_f01.map({s:i for i,s in enumerate(ORDER)})
r.to_csv(OUT+'/fig5d_pair_rows_with_background.tsv',sep='\t',index=False)
# pair-level counts
pc=r.groupby(['trait_class','lineage','state_f05']).size().rename('n').reset_index()
pc['denominator']=pc.trait_class.map(r.drop_duplicates(['trait_class','gwas_name','ensembl']).trait_class.value_counts())
pc['unit']='gene-study pair'
# gene level
rows=[]
for how in ['best','worst']:
    for tag in ['f05','f01']:
        f=('min' if how=='best' else 'max')
        g=r.groupby(['trait_class','lineage','ensembl'])['rank_'+tag].agg(f).reset_index()
        g['state']=g['rank_'+tag].map(dict(enumerate(ORDER)))
        c=g.groupby(['trait_class','lineage','state']).size().rename('n').reset_index()
        c['denominator']=c.trait_class.map(g.drop_duplicates(['trait_class','ensembl']).trait_class.value_counts())
        c['unit']='gene'; c['gene_collapse']=how; c['partial_floor']=tag
        rows.append(c)
gc=pd.concat(rows); gc['proportion']=gc.n/gc.denominator
gc.to_csv(OUT+'/fig5d_gene_level_state_counts.tsv',sep='\t',index=False)
pc.to_csv(OUT+'/fig5d_pair_level_state_counts_floor05.tsv',sep='\t',index=False)
print(gc[(gc.gene_collapse=='best')&(gc.partial_floor=='f05')].pivot_table(index=['trait_class','lineage'],columns='state',values='n',fill_value=0).to_string())
print(pc.pivot_table(index=['trait_class','lineage'],columns='state_f05',values='n',fill_value=0).to_string())
# background
r['diff']=r.obs_any-r.exp_any
gb=r.groupby(['trait_class','lineage','ensembl']).agg(obs=('obs_any','mean'),exp=('exp_any','mean'),diff=('diff','mean'),n_pairs=('gwas_name','size')).reset_index()
gb.to_csv(OUT+'/fig5d_gene_level_observed_vs_locus_expected.tsv',sep='\t',index=False)
rng=np.random.default_rng(20260923)
res=[]
for (tc,lin),x in gb.groupby(['trait_class','lineage']):
    w=stats.wilcoxon(x['diff'],zero_method='wilcox')
    # bootstrap over genes for mean obs, exp, ratio of means
    idx=rng.integers(0,len(x),(2000,len(x)))
    o=x.obs.values[idx].mean(1); e=x.exp.values[idx].mean(1)
    res.append(dict(trait_class=tc,lineage=lin,n_genes=len(x),mean_observed=x.obs.mean(),mean_expected=x.exp.mean(),
        ratio=x.obs.mean()/x.exp.mean(),ratio_lo=np.percentile(o/e,2.5),ratio_hi=np.percentile(o/e,97.5),
        obs_lo=np.percentile(o,2.5),obs_hi=np.percentile(o,97.5),exp_lo=np.percentile(e,2.5),exp_hi=np.percentile(e,97.5),
        n_genes_obs_gt_exp=int((x['diff']>0).sum()),median_diff=x['diff'].median(),wilcoxon_p=w.pvalue))
res=pd.DataFrame(res); res['bh_q_10cells']=stats.false_discovery_control(res.wilcoxon_p)
res['family']='10 lineage x trait-class cells, BH'; res['bootstrap']='2000 gene resamples, seed 20260923'
res.to_csv(OUT+'/fig5d_locus_background_summary.tsv',sep='\t',index=False)
pd.set_option('display.width',250); print(res.drop(columns=['family','bootstrap']).to_string())
# pair-level genome-wide fractions for reference
print('pairs any-lineage max mass <0.05:', (r.groupby(['gwas_name','ensembl']).maxsrc.max()<0.05).sum())
