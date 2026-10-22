import pandas as pd, warnings
warnings.filterwarnings('ignore')
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/'
m=pd.read_csv(D+'credible_set_members_annotated.tsv',sep='\t',dtype=str,keep_default_na=False)
g=m[m.pip_sum_gate_passed=='true'].copy()
g['w']=g.normalized_susie_pip.astype(float)
g['unres']=g.promoter_context_status=='unresolved_liftover'
# row fraction and pip-weighted fraction by ancestry
t=g.groupby('ancestry').apply(lambda d: pd.Series({'rows':len(d),'unres_rows':d.unres.sum(),'row_frac':d.unres.mean(),'pip_frac':(d.w*d.unres).sum()/d.w.sum()}))
print(t.round(4))
# noncoding only
nc=g[g.consequence_category.isin(['other_noncoding'])]
print(g.consequence_category.value_counts())
t2=nc.groupby('ancestry').apply(lambda d: pd.Series({'rows':len(d),'row_frac':d.unres.mean(),'pip_frac':(d.w*d.unres).sum()/d.w.sum()}))
print('noncoding only'); print(t2.round(4))
# per study
ps=g.groupby(['study','ancestry','trait_scope']).apply(lambda d: pd.Series({'rows':len(d),'row_frac':round(d.unres.mean(),3)})).reset_index()
pd.set_option('display.width',200); pd.set_option('display.max_rows',100)
print(ps.sort_values(['ancestry','study']).to_string())
mvp=ps[ps.study.str.contains('MVP')]
print('MVP row_frac range',mvp.row_frac.min(),mvp.row_frac.max(), 'n', len(mvp))
