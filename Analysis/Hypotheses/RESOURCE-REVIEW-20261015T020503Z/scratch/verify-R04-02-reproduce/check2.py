import pandas as pd, numpy as np
import statsmodels.formula.api as smf
from scipy import stats
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H=R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
S=R+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/"
d=pd.read_csv(H,sep="\t",keep_default_na=False)
m=pd.read_csv(S+"qualified_donor_metadata.tsv",sep="\t",keep_default_na=False)[["donor","disease_stage_coarse","exclude_stage_analysis"]]
d=d.drop(columns=["disease_stage_coarse","exclude_stage_analysis"]).merge(m,on="donor")
d=d[(d.exclude_stage_analysis.astype(str)=="False")&d.disease_stage_coarse.isin(["Healthy","Steatosis","Steatohepatitis"])].copy()
d["stage"]=d.disease_stage_coarse.map({"Healthy":0,"Steatosis":1,"Steatohepatitis":2})
pe=pd.read_csv(S+"primary_effects.tsv",sep="\t")
for mod,name in [(20,"Ductular"),(8,"ECM")]:
    g=d[(d.cell_type=="hepatocytes")&(d.module==mod)].copy()
    f=smf.ols("score ~ stage + C(dataset)",data=g).fit()
    hc3=f.get_robustcov_results(cov_type="HC3")
    i=list(f.params.index).index("stage")
    # stage centered within dataset
    g["xc"]=g.stage-g.groupby("dataset").stage.transform("mean")
    ssx=(g.xc**2).sum(); g["res"]=f.resid
    share=g.groupby("dataset").apply(lambda t:(t.xc*t.score).sum())/ (g.xc*g.score).sum()
    # dataset-specific residual variance (from pooled model, df-adjusted by n-1 within dataset)
    v=0
    for ds,t in g.groupby("dataset"):
        s2=(t.res**2).sum()/max(len(t)-1,1)
        v+=(t.xc**2).sum()*s2
    se_ds=np.sqrt(v)/ssx
    # WLS with dataset-specific variance (feasible GLS)
    s2map=g.groupby("dataset").res.apply(lambda x:(x**2).sum()/max(len(x)-1,1))
    g["w"]=1/g.dataset.map(s2map).clip(lower=1e-3)
    fw=smf.wls("score ~ stage + C(dataset)",data=g,weights=g.w).fit()
    # DL random effects for two datasets with df>0
    est=[]
    for ds in ["GSE202379","GSE244832"]:
        ff=smf.ols("score ~ stage",data=g[g.dataset==ds]).fit(); est.append((ff.params["stage"],ff.bse["stage"]))
    b=np.array([e[0] for e in est]); s=np.array([e[1] for e in est]); w=1/s**2
    fe=(w*b).sum()/w.sum(); Q=(w*(b-fe)**2).sum(); tau2=max(0,(Q-1)/(w.sum()-(w**2).sum()/w.sum()))
    wr=1/(s**2+tau2); re=(wr*b).sum()/wr.sum(); sere=np.sqrt(1/wr.sum())
    row=pe[(pe.cell_type=="hepatocytes")&(pe.module==mod)].iloc[0]
    print(name,"OLS b %.4f se %.4f | file hc3_se %.4f hc3_p %.4g | HC3 recompute se %.4f | dataset-var SE %.4f | FGLS b %.4f se %.4f p %.3g | DL RE b %.4f se %.4f p %.3g tau2 %.4f"%(f.params["stage"],f.bse["stage"],row.hc3_se,row.hc3_pvalue,hc3.bse[i],se_ds,fw.params["stage"],fw.bse["stage"],fw.pvalues["stage"],re,sere,2*stats.norm.sf(abs(re/sere)),tau2))
    print("  share of OLS slope numerator by dataset:",share.round(3).to_dict())
