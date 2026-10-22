import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R02/'
m=pd.read_pickle(S+'r2_main.pkl')
tier=pd.read_csv(B+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
tr=dict(zip(tier.study_name,tier.trait)); tl=dict(zip(tier.study_name,tier.tier_label))
m=m[m.ensembl.notna()]
best=m[m['PP.H4.susie'].notna()].sort_values(['ensembl','PP.H4.susie'],ascending=[True,False]).groupby('ensembl').head(1)
best=best[best['PP.H4.susie']>0.5]
best['trait']=best.gwas_name.map(tr); best['cls']=best.gwas_name.map(tl)
print(len(best)); print(best.cls.value_counts()); print(best.trait.value_counts())
# tie check: genes whose max SuSiE PP.H4 ties between classes
# Fig1C style: any-class support
g=m[m['PP.H4.susie']>0.5].copy(); g['cls']=g.gwas_name.map(tl)
cs=g.groupby('ensembl').cls.apply(set)
print('enzyme-only',(cs.apply(lambda s:s=={'liver_enzyme'})).sum(),'direct-only',(cs.apply(lambda s:s=={'direct_MASLD'})).sum(),'both',(cs.apply(len)==2).sum())
# context unresolved by ancestry (noncoding members), gated sets
mem=pd.read_csv(B+'GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv',sep='\t',low_memory=False)
mem=mem[mem.pip_sum_gate_passed.astype(str).str.lower()=='true']
mem['unres']=(mem.promoter_context_status!='resolved')|(mem.lineage_accessibility_status!='resolved')
mem['cls']=np.where(mem.trait_scope.str.startswith('tier1'),'direct','enzyme')
nc=mem[mem.consequence_category.isin(['other_noncoding','synonymous_or_utr'])]
t=nc.groupby(['cls','ancestry']).apply(lambda d: pd.Series({'n':len(d),'unres_pip_frac':(d.normalized_susie_pip*d.unres).sum()/d.normalized_susie_pip.sum()}))
print(t)
print(nc.groupby('ancestry').apply(lambda d: d.lineage_accessibility_status.value_counts(normalize=True).round(3).to_dict()))
