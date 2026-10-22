"""R11-07: within-MASLD (n=46) partial Spearman of GNMT/MAT1A/CYP2C19 protein with fibrosis and NAS.
Same estimator as 02_proteomics_projection.R partial_spearman (rank, residualize on batch, age_z, bmi_z, sex),
restricted to MASLD participants. Family: 3 proteins x 2 features = 6 tests, BH. Input is the display matrix
panel4c_adjusted_abundance.tsv (batch/age/BMI/sex already removed with group retained)."""
import numpy as np, pandas as pd
from scipy import stats
B='Analysis/Multimodal_Program_Projection/results/proteomics/'
m=pd.read_csv(B+'panel4c_metadata.tsv',sep='\t'); a=pd.read_csv(B+'panel4c_adjusted_abundance.tsv',sep='\t').set_index('gene')
print(m.groupby('group')[['Fibrosis','NAS']].describe().T.to_string())
d=m[m.group=='MASLD'].copy()
d['age_z']=(d.age-d.age.mean())/d.age.std(); d['bmi_z']=(d.bmi-d.bmi.mean())/d.bmi.std()
X=pd.get_dummies(d[['acquisition_batch','sex']].astype(str),drop_first=True).astype(float)
X=np.column_stack([np.ones(len(d)),X.values,d.age_z,d.bmi_z]); rank=np.linalg.matrix_rank(X)
res=lambda y: y-X@np.linalg.lstsq(X,y,rcond=None)[0]
rows=[]
for g in ['GNMT','MAT1A','CYP2C19']:
    p=a.loc[g,d.sample_id].astype(float).values
    for f in ['Fibrosis','NAS']:
        s=d[f].astype(float).values
        r=np.corrcoef(res(stats.rankdata(p)),res(stats.rankdata(s)))[0,1]
        df=len(d)-rank-1; t=r*np.sqrt(df/(1-r*r)); pv=2*stats.t.sf(abs(t),df)
        rows.append(dict(gene=g,feature=f,n=len(d),rho=r,df=df,pvalue=pv))
o=pd.DataFrame(rows); o['padj_BH_6']=stats.false_discovery_control(o.pvalue)
o['family']='3 proteins x 2 features within 46 MASLD participants'
print(o.to_string()); o.to_csv(__import__('sys').argv[1],sep='\t',index=False)
