"""Recompute Fig 3E NMF activity (F0-centred, adjusted for cohort + inferred sex) from per-unit fixed k6 NMF scores
on the 387 staged units; try candidate score scales and compare with the stage-extension activity table."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
m=pd.read_csv(O+"meta_774.tsv",sep="\t",keep_default_na=False,na_values=["NA",""])
sc=pd.read_csv("RNA-seq/results/manuscript_release/candidates/fixed-nmf-corrected-scores-20261009T141312Z/fixed_nmf_scores_full_stage_k6.tsv",sep="\t")
d=m.merge(sc,on='sample_id')
s=d[d.dataset.isin(['GSE130970','GSE135251','GSE162694'])&d.fibrosis_stage.notna()].copy()
s['stage']='F'+s.fibrosis_stage.astype(int).astype(str)
print("n staged",len(s), s.stage.value_counts().sort_index().to_dict())
ref=pd.read_csv("figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/prespecified_nmf_program_activity.tsv",sep="\t")
ref=ref[ref.axis=='fibrosis']
res=[]
for P in ['P1','P2','P3','P4','P5']:
    for scale in ['raw','z774','z387','zF0']:
        y=s[P].astype(float)
        if scale=='z774': y=(y-d[P].mean())/d[P].std()
        if scale=='z387': y=(y-y.mean())/y.std()
        if scale=='zF0': f0=y[s.stage=='F0']; y=(y-f0.mean())/f0.std()
        s['y']=y
        fit=smf.ols("y ~ C(dataset) + C(inferred_sex) + C(stage, Treatment('F0'))",data=s).fit()
        for st in ['F1','F2','F3','F4']:
            k=f"C(stage, Treatment('F0'))[T.{st}]"
            res.append(dict(program=P,scale=scale,group=st,est=fit.params[k],se=fit.bse[k]))
r=pd.DataFrame(res)
mm=r.merge(ref[['program','group','centered_activity','SE']],on=['program','group'])
mm['diff']=mm.est-mm.centered_activity
print(mm.groupby('scale')['diff'].apply(lambda x: x.abs().max()))
best=mm[mm.scale==mm.groupby('scale')['diff'].apply(lambda x:x.abs().max()).idxmin()]
print(best.round(4).to_string())
best.to_csv(O+"fig3e_recomputed_nmf.tsv",sep="\t",index=False)
# Stage n and SE-based 95% CIs crossing zero
best=best.assign(lo=best.est-1.96*best.se,hi=best.est+1.96*best.se)
print(best[['program','group','est','lo','hi']].round(3).to_string())
