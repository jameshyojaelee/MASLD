import pandas as pd
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
mem=pd.read_csv(B+'GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv',sep='\t',low_memory=False,
  usecols=['study','ancestry','chromosome','position','pip_sum_gate_passed','lineage_accessibility_status','promoter_context_status','consequence_category'])
mem=mem[mem.pip_sum_gate_passed.astype(str).str.lower()=='true']
mem['site']=mem.chromosome.astype(str)+':'+mem.position.astype(str)
print('promoter_context_status by ancestry'); print(mem.groupby('ancestry').promoter_context_status.value_counts(normalize=True).round(3).to_string())
g=mem.groupby('site').lineage_accessibility_status.nunique()
multi=mem[mem.site.isin(g[g>1].index)]
print('sites with >1 lineage status across studies:',(g>1).sum(),'of',len(g),'sites present in >1 row:',(mem.groupby('site').size()>1).sum())
eu=set(mem[(mem.ancestry=='EUR')&(mem.lineage_accessibility_status=='resolved')].site)
ne=mem[(mem.ancestry!='EUR')&(mem.site.isin(eu))]
print('non-EUR rows at sites resolved in a EUR study:',len(ne),'of which unresolved_liftover:',(ne.lineage_accessibility_status=='unresolved_liftover').sum())
