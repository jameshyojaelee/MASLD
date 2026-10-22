import pandas as pd, numpy as np, statsmodels.formula.api as smf
from scipy.stats import norm
g=pd.read_csv('calls-20261009T181453Z/person_genotypes.tsv',sep='\t')
p=pd.read_csv('scores-20261009T181306Z/person_phenotypes.tsv',sep='\t')
cohorts=['GSE130970','GSE135251','GSE162694','GSE174478','GSE193066','GSE213621','GSE240729']
d=p.merge(g,on=['person_id','cohort'])
d=d[d.cohort.isin(cohorts)&d.stage.notna()&d.score_ecm.notna()]
# alt dosage: H3 tested allele counts; map RR/RA/AA to 0/1/2 alt
dose={'RR':0,'RA':1,'AA':2}
for v in ['rs738409','rs58542926','rs2642438','rs1260326']:
    d[v]=d[v+'_call'].map(dose)
print('persons staged with ECM score:',len(d)); print(d.groupby('cohort').size().to_string())
print('PNPLA3 called:',d.rs738409.notna().sum(), ' all4 called:',d[['rs738409','rs58542926','rs2642438','rs1260326']].notna().all(1).sum())
# within-cohort standardize ECM (H3 used cohort reference SD; here cohort SD, approximate)
d['ecm_z']=d.groupby('cohort').score_ecm.transform(lambda x:(x-x.mean())/x.std())
d['stage_c']=d.stage.astype(float)
# path a: genotype->stage ; path b: stage->ECM
a=smf.ols('stage_c~rs738409+C(cohort)+C(sex)',d).fit(cov_type='HC3')
b=smf.ols('ecm_z~stage_c+C(cohort)+C(sex)',d).fit(cov_type='HC3')
t=smf.ols('ecm_z~rs738409+C(cohort)+C(sex)',d).fit(cov_type='HC3')
print('a PNPLA3->stage %.3f se %.3f n %d'%(a.params['rs738409'],a.bse['rs738409'],a.nobs))
print('b stage->ECM SD/stage %.3f se %.3f'%(b.params['stage_c'],b.bse['stage_c']))
ab=a.params['rs738409']*b.params['stage_c']
print('expected total effect a*b = %.3f SD/allele'%ab)
print('observed total PNPLA3->ECM (no stage adj) %.3f se %.3f p %.3f n %d'%(t.params['rs738409'],t.bse['rs738409'],t.pvalues['rs738409'],t.nobs))
se=t.bse['rs738409']
for alpha,side in [(0.05,2),(0.05,1)]:
    z=norm.ppf(1-alpha/side)+norm.ppf(0.8); print('MDE alpha %.2f %d-sided: %.3f SD'%(alpha,side,z*se))
z=norm.ppf(0.975)+norm.ppf(0.8); print('n needed for a*b at 2-sided 0.05: %d'%int(len(d)*(z*se/ab)**2))
# 4-variant score with in-data stage weights (optimistic upper bound on instrument strength)
e=d.dropna(subset=['rs738409','rs58542926','rs2642438','rs1260326']).copy()
w=smf.ols('stage_c~rs738409+rs58542926+rs2642438+rs1260326+C(cohort)+C(sex)',e).fit(cov_type='HC3')
print(w.params.filter(like='rs').round(3).to_dict(), w.bse.filter(like='rs').round(3).to_dict())
e['prs']=sum(w.params[v]*e[v] for v in ['rs738409','rs58542926','rs2642438','rs1260326'])
e['prs_z']=(e.prs-e.prs.mean())/e.prs.std()
a2=smf.ols('stage_c~prs_z+C(cohort)+C(sex)',e).fit(cov_type='HC3')
b2=smf.ols('ecm_z~stage_c+C(cohort)+C(sex)',e).fit(cov_type='HC3')
t2=smf.ols('ecm_z~prs_z+C(cohort)+C(sex)',e).fit(cov_type='HC3')
ab2=a2.params['prs_z']*b2.params['stage_c']
print('PRS(in-data weights, overfit upper bound) n %d: a %.3f, b %.3f, expected total %.3f SD per PRS SD; SE %.3f; MDE 2-sided %.3f; n needed %d'%(len(e),a2.params['prs_z'],b2.params['stage_c'],ab2,t2.bse['prs_z'],z*t2.bse['prs_z'],int(len(e)*(z*t2.bse['prs_z']/ab2)**2)))
print('observed PRS->ECM %.3f p %.3f'%(t2.params['prs_z'],t2.pvalues['prs_z']))
