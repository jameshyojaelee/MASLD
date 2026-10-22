"""Re-derive Ductular injury (hepatocytes module 20) primary selection checks.
Model: score ~ stage_ordinal + factor(dataset), OLS and HC3, BH over 117 programs (n=117 fixed)."""
import sys, numpy as np, pandas as pd
from scipy import stats
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
SC=ROOT+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/work/donor_program_scores_primary.tsv"
MD=ROOT+"Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
REF=ROOT+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv"
STAGES={"Healthy":0,"Steatosis":1,"Steatohepatitis":2}
sc=pd.read_csv(SC,sep="\t")[["cell_type","module","donor","score"]]
md=pd.read_csv(MD,sep="\t",keep_default_na=False)
md["excl"]=md.exclude_stage_analysis.astype(str).str.upper()=="TRUE"
el=md[(~md.excl)&md.disease_stage_coarse.isin(STAGES)].copy()
el["stage"]=el.disease_stage_coarse.map(STAGES).astype(float)
el=el.rename(columns={"sample":"donor"})[["donor","dataset","stage","disease_stage_coarse"]]
d=sc.merge(el,on="donor",how="inner")
progs=sc[["cell_type","module"]].drop_duplicates().sort_values(["cell_type","module"]).reset_index(drop=True)
assert len(progs)==117, len(progs)

def fit(sub, stage=None):
    y=sub.score.values; s=sub.stage.values if stage is None else stage
    ds=sub.dataset.values; lev=sorted(set(ds))
    X=np.column_stack([np.ones(len(y)),s]+[(ds==l).astype(float) for l in lev[1:]])
    n,p=X.shape
    if np.linalg.matrix_rank(X)<p or n-p<=0: return (np.nan,)*5
    XtXi=np.linalg.inv(X.T@X); b=XtXi@X.T@y; r=y-X@b; rdf=n-p
    se=np.sqrt(r@r/rdf*XtXi[1,1])
    h=np.einsum("ij,jk,ik->i",X,XtXi,X)
    if np.any(1-h<=1.5e-8): hse=np.nan
    else:
        w=(r/(1-h))**2; V=XtXi@(X.T*w)@X@XtXi; hse=np.sqrt(V[1,1])
    t=b[1]/se; th=b[1]/hse
    return b[1], 2*stats.t.sf(abs(t),rdf), 2*stats.t.sf(abs(th),rdf), t, th

def bh(p, m=117):
    p=np.asarray(p,float); q=np.full(len(p),np.nan); ok=np.isfinite(p)
    pv=p[ok]; o=np.argsort(pv); ranks=np.arange(1,len(pv)+1)
    qq=np.minimum.accumulate((pv[o]*m/ranks)[::-1])[::-1]; qq=np.minimum(qq,1)
    out=np.empty(len(pv)); out[o]=qq; q[ok]=out; return q

def family(dd):
    rows=[]
    for ct,mod in progs.itertuples(index=False):
        sub=dd[(dd.cell_type==ct)&(dd.module==mod)]
        b,p,ph,t,th=fit(sub); rows.append((ct,mod,b,p,ph,len(sub)))
    f=pd.DataFrame(rows,columns=["cell_type","module","beta","p","p_hc3","n"])
    f["q"]=bh(f.p); f["q_hc3"]=bh(f.p_hc3); return f
