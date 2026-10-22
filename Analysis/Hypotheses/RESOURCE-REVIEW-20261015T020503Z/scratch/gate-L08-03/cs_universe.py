import pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/'
c=pd.read_csv(D+'credible_sets.csv',low_memory=False)
cs=c[c.susie_cs>0]
g=cs.groupby('study').apply(lambda d: d.groupby(['locus','susie_cs']).ngroups).sort_values(ascending=False)
print('CS per study (top 12):'); print(g.head(12).to_string())
a=pd.read_csv(D+'candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_pip_architecture.tsv',sep='\t')
print('architecture CS',len(a),'gated',(a.pip_sum_gate_passed==True).sum(),'studies',a.study.nunique())
print('studies in credible_sets.csv not in architecture:',sorted(set(cs.study)-set(a.study)))
print('CS of those studies:',int(g[[s for s in g.index if s not in set(a.study)]].sum()))
m=pd.read_csv(D+'candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv',sep='\t',low_memory=False)
print('member rows',len(m),'unique variants',m.variant_id.nunique())
