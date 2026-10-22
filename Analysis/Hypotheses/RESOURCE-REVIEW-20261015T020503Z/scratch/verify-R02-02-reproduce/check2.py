import pandas as pd, warnings
warnings.filterwarnings('ignore')
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/'
A='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/gwas_atac/'
m=pd.read_csv(D+'credible_set_members_annotated.tsv',sep='\t',dtype=str,keep_default_na=False)
g=m[m.pip_sum_gate_passed=='true'].copy()
g['unres']=g.promoter_context_status=='unresolved_liftover'
v=pd.read_csv(A+'variant_overlap_summary.csv',dtype=str,keep_default_na=False)
print('liftover table rows',len(v),'empty hg38',(v.chr_hg38=='').sum())
v['vid']=v.chromosome+':'+v.position+':'+v.allele1+':'+v.allele2
v['vid_swap']=v.chromosome+':'+v.position+':'+v.allele2+':'+v.allele1
v['site']=v.chromosome+':'+v.position
studies=set()
for s in v.study_list: studies.update(s.split(';') if ';' in s else s.split(','))
print('studies in liftover table:',sorted(studies))
print('max_pip summary', v.max_pip.astype(float).describe())
# unresolved rows: is the site present with swapped alleles? site present at all?
u=g[g.unres]
g['site']=g.chromosome+':'+g.position
u=g[g.unres]
print('unresolved rows',len(u))
print('  vid matches swapped allele order:',u.variant_id.isin(set(v.vid_swap)).sum())
print('  site present in liftover table:',u.site.isin(set(v.site)).sum())
print('  site absent:',(~u.site.isin(set(v.site))).sum())
# study in liftover table study list?
u2=u.assign(study_in=u.study.isin(studies))
print(pd.crosstab(u2.ancestry,u2.study_in))
r=g[~g.unres]
print('resolved rows study_in', r.study.isin(studies).value_counts().to_dict())
# PIP of unresolved vs resolved
print('susie_pip resolved median', r.susie_pip.astype(float).median(), 'unresolved median', u.susie_pip.astype(float).median())
print('unresolved with susie_pip>=0.01', (u.susie_pip.astype(float)>=0.01).sum(), ' >=0.1', (u.susie_pip.astype(float)>=0.1).sum())
print('resolved with susie_pip<0.01', (r.susie_pip.astype(float)<0.01).sum())
