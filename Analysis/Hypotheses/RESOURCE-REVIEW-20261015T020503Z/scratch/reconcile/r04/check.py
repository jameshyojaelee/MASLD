import pandas as pd, numpy as np
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H=R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
S=R+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/"
d=pd.read_csv(H,sep="\t",keep_default_na=False)
m=pd.read_csv(S+"qualified_donor_metadata.tsv",sep="\t",keep_default_na=False)
m=m[["donor","disease_stage_coarse","exclude_stage_analysis","stage_ordinal"]].rename(columns={"disease_stage_coarse":"stage_q","exclude_stage_analysis":"excl_q","stage_ordinal":"ord_q"})
d=d.merge(m,on="donor",how="left")
d=d[(d.excl_q.astype(str)=="False")&(d.stage_q.isin(["Healthy","Steatosis","Steatohepatitis"]))].copy()
d["stage"]=d.stage_q.map({"Healthy":0,"Steatosis":1,"Steatohepatitis":2})
pe=pd.read_csv(S+"primary_effects.tsv",sep="\t")
def fit(df):
    f=smf.ols("score ~ stage + C(dataset)",data=df).fit()
    h=f.get_robustcov_results(cov_type="HC3")
    return f.params["stage"],f.bse["stage"],f.pvalues["stage"],len(df),f.df_resid
rows=[]
for (ct,mod),g in d.groupby(["cell_type","module"]):
    if g.stage.nunique()<2: continue
    try:
        b,se,p,n,df=fit(g)
    except Exception as e:
        continue
    rows.append((ct,mod,b,se,p,n))
res=pd.DataFrame(rows,columns=["cell_type","module","beta","se","p","n"])
res["q"]=multipletests(res.p,method="fdr_bh")[1]
cmp=res.merge(pe[["cell_type","module","beta","pvalue","qvalue","n_donors"]],on=["cell_type","module"])
print("programs fitted",len(res),"compared",len(cmp),"max|dbeta|",(cmp.beta_x-cmp.beta_y).abs().max(),"max|dq|",(cmp.q-cmp.qvalue).abs().max())
print("primary q<0.05:",(res.q<0.05).sum())
for mod,name in [(20,"Ductular injury"),(8,"Stromal ECM")]:
    g=d[(d.cell_type=="hepatocytes")&(d.module==mod)]
    print("\n==",name,"n",len(g),g.groupby(["dataset","stage"]).size().to_dict())
    b,se,p,n,df=fit(g); print("pooled FE OLS beta %.6f se %.6f p %.3g n %d"%(b,se,p,n))
    f=smf.ols("score ~ stage + C(dataset)",data=g).fit()
    print("pooled residual SD %.4f"%np.sqrt(f.scale))
    g=g.assign(res=f.resid)
    print("pooled-model residual SD by dataset:", g.groupby("dataset").res.apply(lambda x: np.sqrt((x**2).sum()/max(len(x)-1,1))).round(3).to_dict())
    est={}
    resid={}
    for ds,gg in g.groupby("dataset"):
        if gg.stage.nunique()<2:
            print(ds,"n",len(gg),"one stage; sd",round(gg.score.std(),4)); continue
        ff=smf.ols("score ~ stage",data=gg).fit()
        sd=np.sqrt(ff.scale) if ff.df_resid>0 else float('nan')
        print("%s n %d slope %.4f se %.4f p %.4g resid SD %.4f df %d"%(ds,len(gg),ff.params["stage"],ff.bse["stage"] if ff.df_resid>0 else float('nan'),ff.pvalues["stage"] if ff.df_resid>0 else float('nan'),sd,ff.df_resid))
        if ff.df_resid>0:
            est[ds]=(ff.params["stage"],ff.bse["stage"]); resid[ds]=ff.resid.values
    ks=list(est); b_=np.array([est[k][0] for k in ks]); s_=np.array([est[k][1] for k in ks]); w=1/s_**2
    ivw=(w*b_).sum()/w.sum(); seivw=np.sqrt(1/w.sum()); z=ivw/seivw
    Q=(w*(b_-ivw)**2).sum(); pQ=stats.chi2.sf(Q,len(ks)-1)
    print("IVW (%s) %.4f se %.4f z-p %.3g; Cochran Q %.3f df %d p %.4g"%(",".join(ks),ivw,seivw,2*stats.norm.sf(abs(z)),Q,len(ks)-1,pQ))
    # Brown-Forsythe on within-dataset residuals of per-dataset fits
    print("BF (median-centered Levene) on per-dataset-fit residuals:",stats.levene(*[resid[k] for k in ks],center="median"))
    # BF on pooled-model residuals, datasets with n>=3
    grp=[gg.res.values for ds,gg in g.groupby("dataset") if len(gg)>=3]
    print("BF on pooled-model residuals, datasets n>=3:",stats.levene(*grp,center="median"))
    grp2=[g[g.dataset==k].res.values for k in ks]
    print("BF on pooled-model residuals, GSE202379 vs GSE244832:",stats.levene(*grp2,center="median"))
# drop GSE244832
d2=d[d.dataset!="GSE244832"]
rows=[]
for (ct,mod),g in d2.groupby(["cell_type","module"]):
    if g.stage.nunique()<2: continue
    try: b,se,p,n,df=fit(g)
    except Exception: continue
    if not np.isfinite(p): continue
    rows.append((ct,mod,b,se,p,n))
r2=pd.DataFrame(rows,columns=["cell_type","module","beta","se","p","n"]); r2["q"]=multipletests(r2.p,method="fdr_bh")[1]
print("\nDrop GSE244832: programs",len(r2),"q<0.05",(r2.q<0.05).sum(),"min q",r2.q.min())
print(r2[(r2.cell_type=="hepatocytes")&(r2.module.isin([8,20]))])
