# Within-dataset label permutation calibration of the Fig 4C primary stage test (117 programs, 62 donors).
import numpy as np, pandas as pd, sys
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=pd.read_csv(root+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv',sep='\t',usecols=['cell_type','module','program_uid','donor','score','dataset'])
M=pd.read_csv(root+'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/qualified_donor_metadata.tsv',sep='\t')
M=M[(M.exclude_stage_analysis.astype(str).str.upper()=='FALSE')&(M.stage_ordinal.isin([0,1,2]))]
print('eligible donors',len(M), M.groupby('dataset').size().to_dict())
S=S.merge(M[['donor','stage_ordinal']],on='donor')
P=S.pivot_table(index='donor',columns='program_uid',values='score')
meta=M.set_index('donor').loc[P.index]
ds=pd.get_dummies(meta.dataset).astype(float).values
stage=meta.stage_ordinal.astype(float).values
prog=P.columns.values; Y=P.values
ok=~np.isnan(Y)
def fit(stage):
    out=np.full(Y.shape[1],np.nan); bet=np.full(Y.shape[1],np.nan)
    for j in range(Y.shape[1]):
        m=ok[:,j]; X=np.column_stack([ds[m],stage[m]]); keep=X.std(0)>0; keep[-1]=True
        # drop empty dataset columns
        X=X[:,np.r_[ (X[:,:-1].sum(0)>0), True]]
        y=Y[m,j]; XtX=np.linalg.pinv(X.T@X); b=XtX@X.T@y; r=y-X@b; df=m.sum()-np.linalg.matrix_rank(X)
        s2=r@r/df; se=np.sqrt(s2*XtX[-1,-1]); out[j]=b[-1]/se; bet[j]=b[-1]
    return out,bet
from scipy import stats
t0,b0=fit(stage)
uid={'Ductular':'hotspot_hepatocytes_48f39dd4d817a10e','ECM':'hotspot_hepatocytes_f05c535ae5bbc0b9'}
for k,u in uid.items():
    j=list(prog).index(u); print(k,'beta',b0[j],'t',t0[j])
rng=np.random.default_rng(20261015)
B=int(sys.argv[1]) if len(sys.argv)>1 else 5000
groups=[np.where(meta.dataset.values==d)[0] for d in meta.dataset.unique()]
T=np.empty((B,Y.shape[1]))
for i in range(B):
    s=stage.copy()
    for g in groups: s[g]=rng.permutation(stage[g])
    T[i],_=fit(s)
res=[]
for k,u in uid.items():
    j=list(prog).index(u)
    pe=(1+(np.abs(T[:,j])>=abs(t0[j])).sum())/(B+1)
    pmax=(1+(np.nanmax(np.abs(T),1)>=abs(t0[j])).sum())/(B+1)
    res.append((k,b0[j],t0[j],pe,pmax))
    print(k,'per-program permutation p',pe,'maxT family-wise p',pmax)
# calibration: fraction of null nominal p<0.05 using t-dist parametric p with df
df=len(meta)-ds.shape[1]-1
pn=2*stats.t.sf(np.abs(T),df)
print('null fraction parametric p<0.05 across programs and perms', np.nanmean(pn<0.05), 'p<0.001', np.nanmean(pn<0.001))
# empirical BH over 117 with empirical per-program p
pemp=(1+(np.abs(T)>=np.abs(t0)).sum(0))/(B+1)
from statsmodels.stats.multitest import multipletests
q=multipletests(pemp,method='fdr_bh')[1]
for k,u in uid.items():
    j=list(prog).index(u); print(k,'empirical p',pemp[j],'BH117 q from empirical p',q[j])
print('n programs',len(prog),'min empirical q',q.min())
