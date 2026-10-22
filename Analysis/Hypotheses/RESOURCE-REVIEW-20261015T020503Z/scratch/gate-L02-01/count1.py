import pandas as pd, gzip, re
a=pd.read_csv('data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',sep='\t',low_memory=False)
lv=a[a.SMTSD=='Liver'].copy()
print('liver samples all', len(lv), 'SMAFRZE counts', lv.SMAFRZE.value_counts(dropna=False).to_dict())
lv['donor']=lv.SAMPID.str.split('-').str[:2].str.join('-')
rna=lv[lv.SMAFRZE=='RNASEQ']
print('RNASEQ liver samples',len(rna),'donors',rna.donor.nunique())
s=pd.read_csv('Analysis/MASLD_Model_Benchmark/executions/he-molecular-c1-20260924T004124Z/liver_slides_matched.tsv',sep='\t')
print('slides rows',len(s),'cols',list(s.columns))
pid=[c for c in s.columns if 'atient' in c][0]
print('slide donors',s[pid].nunique())
both=set(rna.donor)&set(s[pid])
print('both',len(both))
# GCT header
with gzip.open('data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz','rt') as f:
    f.readline(); dims=f.readline(); hdr=f.readline().rstrip('\n').split('\t')
print('gct dims',dims.strip(),'hdr first',hdr[:3])
gsamp=hdr[2:]; gd=set('-'.join(x.split('-')[:2]) for x in gsamp)
print('gct samples',len(gsamp),'donors',len(gd),'gct&both',len(gd&both))
# notes on the RNA liver sample
r=rna[rna.donor.isin(both)].copy()
r['n']=r.SMPTHNTS.fillna('').str.lower()
print('notes missing/empty', (r.n=='').sum())
fat=r.n.str.contains('steato|fat')
neg=r.n.str.contains('no significant steatosis|no steatosis|without steatosis|no major')
print('fat mention',fat.sum(),'negated',(fat&neg).sum(),'non-negated',(fat&~neg).sum())
r[['SAMPID','donor','SMPTHNTS','SMRIN','SMTSISCH','SMATSSCR']].to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-01/liver_rna_notes.tsv',sep='\t',index=False)
# covariates
sub=pd.read_csv('data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt',sep='\t')
m=r.merge(sub,left_on='donor',right_on='SUBJID',how='left')
print('n',len(m),'nonmissing SMRIN',m.SMRIN.notna().sum(),'SMTSISCH',m.SMTSISCH.notna().sum(),'AGE',m.AGE.notna().sum(),'SEX',m.SEX.notna().sum(),'DTHHRDY',m.DTHHRDY.notna().sum())
print('DTHHRDY',m.DTHHRDY.value_counts(dropna=False).to_dict())
print('AGE',m.AGE.value_counts().sort_index().to_dict(),'SEX',m.SEX.value_counts().to_dict())
print('SMTSISCH quantiles',m.SMTSISCH.quantile([0,.25,.5,.75,1]).round(0).tolist())
print('SMRIN quantiles',m.SMRIN.quantile([0,.25,.5,.75,1]).round(1).tolist())
m.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-01/liver_rna_meta.tsv',sep='\t',index=False)
