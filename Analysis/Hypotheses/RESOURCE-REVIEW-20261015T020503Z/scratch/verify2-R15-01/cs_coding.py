import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
f='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv'
m=pd.read_csv(f,sep='\t',low_memory=False,dtype={'locus':str})
print('rows',len(m),'studies',m.study.nunique())
g=m[m.pip_sum_gate_passed.astype(str).str.lower()=='true']
print('gate-passed rows',len(g),'studies',g.study.nunique())
cols=['study','trait_scope','variant_id','susie_pip','normalized_susie_pip','consequence_category','pip_sum_gate_passed']
for pos in [(1,220970028),(19,19379549),(22,44324727)]:
    s=m[(m.chromosome.astype(str)==str(pos[0]))&(m.position==pos[1])]
    print('\n== site',pos,'n rows',len(s)); print(s[cols].sort_values('normalized_susie_pip',ascending=False).to_string())
for gate in ['all','gate']:
    d=m if gate=='all' else g
    c=d[d.consequence_category.isin(['protein_altering','canonical_splice'])&(d.normalized_susie_pip>=0.5)]
    print(f'\n[{gate}] coding/splice rows PIP>=0.5:',len(c),'distinct sites:',c.variant_id.nunique(),'(chr:pos):',c[['chromosome','position']].drop_duplicates().shape[0],'studies:',c.study.nunique())
    print(c.groupby('trait_scope').agg(sites=('variant_id','nunique'),studies=('study','nunique')))
    c2=d[d.consequence_category.isin(['protein_altering','canonical_splice'])&(d.susie_pip>=0.5)]
    print(f'[{gate}] same with raw susie_pip>=0.5: sites',c2.variant_id.nunique(),'studies',c2.study.nunique())
