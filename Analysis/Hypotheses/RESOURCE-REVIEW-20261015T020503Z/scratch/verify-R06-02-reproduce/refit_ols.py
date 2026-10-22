import pandas as pd, numpy as np
from scipy import stats
S='Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv'
M='Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv'
sc=pd.read_csv(S,sep='\t'); md=pd.read_csv(M,sep='\t')
print('meta exclude counts:\n', md.groupby(['exclude_stage_analysis','disease_stage_coarse']).size())
stage_map={'Healthy':0,'Steatosis':1,'Steatohepatitis':2}
md=md[(md.exclude_stage_analysis==False)&(md.disease_stage_coarse.isin(stage_map))].copy()
md['stage']=md.disease_stage_coarse.map(stage_map)
print('eligible donors', md['sample'].nunique(), md.groupby('dataset').size().to_dict())
d=sc[['cell_type','module','donor','score']].merge(md[['sample','dataset','stage']],left_on='donor',right_on='sample')
d=d[np.isfinite(d.score)]
rows=[]
for (ct,mod),g in d.groupby(['cell_type','module']):
    X=pd.get_dummies(g['dataset'],drop_first=True,dtype=float)
    X.insert(0,'stage',g['stage'].astype(float).values); X.insert(0,'const',1.0)
    X=X.values; y=g.score.values; n,p=X.shape
    if np.linalg.matrix_rank(X)<p: rows.append((ct,mod,n,np.nan,np.nan,np.nan)); continue
    XtXi=np.linalg.inv(X.T@X); b=XtXi@X.T@y; e=y-X@b; df=n-p
    s2=e@e/df; se=np.sqrt(s2*XtXi[1,1]); pv=2*stats.t.sf(abs(b[1]/se),df)
    h=np.einsum('ij,jk,ik->i',X,XtXi,X); w=e**2/(1-h)**2
    V=XtXi@(X.T*w)@X@XtXi; se3=np.sqrt(V[1,1]); p3=2*stats.t.sf(abs(b[1]/se3),df)
    rows.append((ct,mod,n,b[1],pv,p3))
r=pd.DataFrame(rows,columns=['cell_type','module','n','beta','p','p_hc3'])
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p)
    q=p[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; out=np.empty(n); out[o]=np.minimum(q,1); return out
r['q']=bh(r.p); r['q_hc3']=bh(r.p_hc3)
print('programs fitted', len(r))
print(r[(r.cell_type=='hepatocytes')&(r.module.isin([8,20]))].to_string())
print('q<0.05:', r[r.q<0.05][['cell_type','module','beta','q']].to_string())
print('HC3 q<0.05 count', (r.q_hc3<0.05).sum())
pe=pd.read_csv('figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv',sep='\t')
c=r.merge(pe[['cell_type','module','beta','pvalue','hc3_pvalue','n_donors']],on=['cell_type','module'])
print('max |beta diff|',(c.beta_x-c.beta_y).abs().max(),'max |p diff|',(c.p-c.pvalue).abs().max(),'max |hc3 p diff|',(c.p_hc3-c.hc3_pvalue).abs().max(),'n mismatch',(c.n!=c.n_donors).sum())
