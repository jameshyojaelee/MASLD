import pandas as pd, warnings
warnings.filterwarnings('ignore')
pd.set_option('display.width',220)
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/'
m=pd.read_csv(D+'credible_set_members_annotated.tsv',sep='\t',dtype=str,keep_default_na=False)
g=m[m.pip_sum_gate_passed=='true'].copy()
g['unres']=g.promoter_context_status=='unresolved_liftover'
g['site']=g.chromosome+':'+g.position
print('TM6SF2 19:19379549 rows (gated):')
print(g[g.site=='19:19379549'][['study','ancestry','variant_id','susie_pip','promoter_context_status','abc_enhancer_overlap','consequence_category']].to_string())
print('all rows (incl ungated):')
print(m[(m.chromosome=='19')&(m.position=='19379549')][['study','ancestry','variant_id','pip_sum_gate_passed','promoter_context_status']].to_string())
# sites resolved in a EUR study
eur_res=set(g[(g.ancestry=='EUR')&(~g.unres)].site)
ne=g[(g.ancestry!='EUR')&(g.site.isin(eur_res))]
print('non-EUR rows at EUR-resolved sites:',len(ne),' unresolved:',ne.unres.sum())
# by variant_id instead of site
eur_res_v=set(g[(g.ancestry=='EUR')&(~g.unres)].variant_id)
ne2=g[(g.ancestry!='EUR')&(g.variant_id.isin(eur_res_v))]
print('non-EUR rows at EUR-resolved variant_ids (same allele order):',len(ne2),' unresolved:',ne2.unres.sum())
# sites with mixed labels across rows
mix=g.groupby('site').unres.agg(['min','max'])
print('sites total',len(mix),' with mixed labels',((mix['min']==False)&(mix['max']==True)).sum())
# ABC by ancestry
print(pd.crosstab(g.ancestry,g.abc_enhancer_overlap,margins=True))
# trait scope by ancestry rows and credible sets
cs=g.drop_duplicates('credible_set_uid')
print('gated credible-set instances by scope x ancestry'); print(pd.crosstab(cs.trait_scope,cs.ancestry,margins=True))
for sc in cs.trait_scope.unique():
    d=cs[cs.trait_scope==sc]; print(sc,'non-EUR instance frac',round((d.ancestry!='EUR').mean(),4), 'n',len(d))
gr=g
for sc in gr.trait_scope.unique():
    d=gr[gr.trait_scope==sc]; print(sc,'non-EUR member-row frac',round((d.ancestry!='EUR').mean(),4))
