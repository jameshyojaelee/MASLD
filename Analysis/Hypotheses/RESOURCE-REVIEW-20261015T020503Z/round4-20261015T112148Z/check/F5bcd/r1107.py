import pandas as pd, numpy as np
from scipy import stats
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/results/proteomics/'
m=pd.read_csv(B+'panel4c_metadata.tsv',sep='\t'); a=pd.read_csv(B+'panel4c_adjusted_abundance.tsv',sep='\t').set_index('gene')
d=m[m.group=='MASLD']
Z=pd.get_dummies(d[['acquisition_batch','sex']].astype(str),drop_first=True).astype(float)
Z['age']=d.age.values; Z['bmi']=d.bmi.values; Z.insert(0,'c',1.0); Z=Z.values
R=lambda y: y-Z@np.linalg.pinv(Z)@y
for g in ['GNMT','MAT1A','CYP2C19']:
    p=a.loc[g,d.sample_id].astype(float).values
    for f in ['Fibrosis','NAS']:
        r=stats.pearsonr(R(stats.rankdata(p)),R(stats.rankdata(d[f].values)))[0]
        print(g,f,len(d),round(r,4), 'plain spearman',round(stats.spearmanr(p,d[f])[0],4))
