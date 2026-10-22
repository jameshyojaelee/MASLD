# H55 planning: count n by group. Reads metadata, slide list, pilot exposure, note grades, gene lists.
# Does not relate exposure to expression.
import pandas as pd, numpy as np, gzip
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
R=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/'
S=R+'scratch/plan-H55/'
a=pd.read_csv(B+'data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',sep='\t',low_memory=False)
lv=a[(a.SMTSD=='Liver')&(a.SMAFRZE=='RNASEQ')].copy()
lv['donor']=lv.SAMPID.str.split('-').str[:2].str.join('-')
print('liver RNASEQ samples',len(lv),'donors',lv.donor.nunique())
s=pd.read_csv(B+'Analysis/MASLD_Model_Benchmark/executions/he-molecular-c1-20260924T004124Z/liver_slides_matched.tsv',sep='\t')
print('slide table cols',list(s.columns)); print('slide rows',len(s),'slide donors',s.PatientID.nunique())
both=set(lv.donor)&set(s.PatientID); print('both',len(both))
with gzip.open(S+'inputs/gene_reads_2017-06-05_v8_liver.gct.gz','rt') as f:
    f.readline(); f.readline(); hdr=f.readline().rstrip('\n').split('\t')
cs=[c for c in hdr if c.startswith('GTEX-')]
cd={'-'.join(c.split('-')[:2]):c for c in cs}
print('count gct sample cols',len(cs),'donors',len(cd),'count&both',len(set(cd)&both), 'count samples in RNASEQ freeze',len(set(cs)&set(lv.SAMPID)))
with gzip.open(B+'data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz','rt') as f:
    f.readline(); f.readline(); th=f.readline().rstrip('\n').split('\t')
ts=[c for c in th if c.startswith('GTEX-')]
print('tpm gct sample cols',len(ts),'same sample ids as count gct',set(ts)==set(cs))
sub=pd.read_csv(B+'data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt',sep='\t')
print('subject cols',list(sub.columns))
m=lv[lv.donor.isin(both)].merge(sub,left_on='donor',right_on='SUBJID',how='left')
print('n both',len(m),'nonmissing',{c:int(m[c].notna().sum()) for c in ['SMRIN','SMTSISCH','AGE','SEX','DTHHRDY','SMCENTER','SMEXNCRT','SMATSSCR']})
m=m[m.DTHHRDY.notna()].copy(); m['H']=m.DTHHRDY.astype(int)
print('model n',len(m))
print('center',m.SMCENTER.value_counts().to_dict())
print(pd.crosstab(m.H,m.SMCENTER,margins=True))
print('age',m.AGE.value_counts().sort_index().to_dict(),'sex',m.SEX.value_counts().to_dict())
g=pd.read_csv(R+'scratch/gate-L02-01/note_grades.tsv',sep='\t')
m=m.merge(g[['donor','grade','pct','fib_or_nash']],on='donor',how='left')
print('note grade (model n)',m.grade.value_counts().to_dict(),'stated pct',int(m.pct.notna().sum()),'fib_or_nash',int(m.fib_or_nash.sum()))
print('grade x Hardy'); print(pd.crosstab(m.grade,m.H))
f=pd.read_csv(R+'jobs/3111/morph_fat.tsv',sep='\t')
print('pilot status',f.status.value_counts().to_dict())
m=m.merge(f[['donor','status']],on='donor',how='left')
print('pilot status x grade'); print(pd.crosstab(m.grade,m.status))
print('pilot ok by Hardy x center'); print(pd.crosstab(m.H[m.status=='ok'],m.SMCENTER[m.status=='ok']))
print('ischemic time min quantiles',m.SMTSISCH.quantile([0,.25,.5,.75,1]).tolist())
print('SMATSSCR',m.SMATSSCR.value_counts(dropna=False).to_dict())
m[['donor','SAMPID','H','SMCENTER','AGE','SEX','SMRIN','SMTSISCH','SMEXNCRT','SMATSSCR','grade','pct','fib_or_nash','status']].to_csv(S+'design_meta_counts_only.tsv',sep='\t',index=False)
