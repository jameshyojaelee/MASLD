# Does the Fig 7B MPRA outcome table's HepG2_PAOA context match the GEO HepG2 PAOA source calls,
# or the swapped LX2 control set from the p3a deposit?
import pandas as pd, numpy as np
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
ro=pd.read_csv(root+'Analysis/MASLD_Model_Benchmark/executions/gse281364-validated-reconstruction-21066470/validated/replicate_outcomes.tsv.gz',sep='\t')
ro=ro[ro.assay_state=='observed']
ro['act']=np.log2((ro.RNA+1)/(ro.DNA+1))
w=ro.pivot_table(index=['element_id','context_id','experimental_replicate'],columns='allele',values='act').dropna()
w['d']=w['alt']-w['ref']
eff=w['d'].groupby(level=[0,1]).mean().unstack()
print('contexts in outcome table:',list(eff.columns))
for cl,cond in [('HepG2','control'),('HepG2','PAOA'),('LX2','control'),('LX2','TGFb')]:
    s=pd.read_csv(root+f'GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction/{cl}.{cond}.official_comparison.tsv.gz',sep='\t').set_index('element_id')
    for c in eff.columns:
        j=eff[[c]].join(s['official_log2FC'],how='inner').dropna()
        print(f'outcome {c:14s} vs source {cl}.{cond:8s} n={len(j):5d} spearman={j[c].corr(j.official_log2FC,method="spearman"):.4f}')
