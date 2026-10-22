# Leave-one-dataset-out refit of all 117 programs (OLS score ~ dataset + stage, classical SE), BH over estimable programs.
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=pd.read_csv(root+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv',sep='\t',usecols=['program_uid','donor','score','dataset'])
M=pd.read_csv(root+'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/qualified_donor_metadata.tsv',sep='\t')
M=M[(M.exclude_stage_analysis.astype(str).str.upper()=='FALSE')&(M.stage_ordinal.isin([0,1,2]))]
S=S.merge(M[['donor','stage_ordinal']],on='donor')
P=S.pivot_table(index='donor',columns='program_uid',values='score')
meta=M.set_index('donor').loc[P.index]
def run(mask):
    ps={};bs={}
    for u in P.columns:
        y=P[u].values; m=mask&~np.isnan(y)
        dsv=meta.dataset.values[m]; st=meta.stage_ordinal.values[m].astype(float)
        X=np.column_stack([pd.get_dummies(dsv).astype(float).values, st])
        b,res,rk,_=np.linalg.lstsq(X,y[m],rcond=None); r=y[m]-X@b; df=m.sum()-rk
        cov=np.linalg.pinv(X.T@X)*(r@r/df); t=b[-1]/np.sqrt(cov[-1,-1]); ps[u]=2*stats.t.sf(abs(t),df); bs[u]=b[-1]
    p=pd.Series(ps); q=pd.Series(multipletests(p.values,method='fdr_bh')[1],index=p.index)
    return pd.Series(bs),p,q
D='hotspot_hepatocytes_48f39dd4d817a10e'; E='hotspot_hepatocytes_f05c535ae5bbc0b9'
allmask=np.ones(len(meta),bool)
for lab,mask in [('all 62',allmask)]+[('drop '+d, meta.dataset.values!=d) for d in ['GSE189600','GSE244832','GSE202379']]:
    b,p,q=run(mask)
    print(f"{lab:20s} n={mask.sum()} Ductular beta={b[D]:.4f} p={p[D]:.3g} q={q[D]:.4f} | ECM beta={b[E]:.4f} q={q[E]:.4f} | programs q<0.05: {(q<0.05).sum()}")
