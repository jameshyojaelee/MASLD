"""Recompute lncRNA inset quantities from the 774-unit fit: unadjusted OR (lncRNA vs protein-coding) of canonical call,
direction shares among calls, and an expression-adjusted logistic OR (AveExpr spline only; detection covariate approximated
by fraction of units with CPM>0 is not available here -> AveExpr only, labelled as approximation)."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
deg=pd.read_csv("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
bt=pd.read_csv("figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv",sep="\t",usecols=['gene_id_versioned','gene_type','display_biotype']).drop_duplicates('gene_id_versioned')
d=deg.merge(bt,left_on='gene',right_on='gene_id_versioned',how='left')
print(d.display_biotype.value_counts(dropna=False).head(8).to_dict())
d=d[d.gene_type.isin(['protein_coding','lncRNA'])].copy()
d['call']=((d.padj<0.05)&(d.logFC.abs()>0.5)).astype(int); d['lnc']=(d.gene_type=='lncRNA').astype(int)
t=pd.crosstab(d.lnc,d.call); print(t)
a,b,c,dd=t.loc[1,1],t.loc[1,0],t.loc[0,1],t.loc[0,0]; OR=(a*dd)/(b*c); se=np.sqrt(1/a+1/b+1/c+1/dd)
print("unadjusted OR",round(OR,3),"CI",round(np.exp(np.log(OR)-1.96*se),3),round(np.exp(np.log(OR)+1.96*se),3))
fit=smf.logit("call ~ lnc + bs(AveExpr, df=5)",data=d).fit(disp=0)
print("AveExpr-adjusted OR",round(np.exp(fit.params['lnc']),3), np.exp(fit.conf_int().loc['lnc']).round(3).tolist())
calls=d[d.call==1]; sh=calls.groupby('gene_type').apply(lambda g: pd.Series(dict(n=len(g),up=(g.logFC>0).sum(),share_up=(g.logFC>0).mean())))
print(sh)
u=pd.crosstab(calls.lnc,calls.logFC>0); a,b,c,dd=u.loc[1,True],u.loc[1,False],u.loc[0,True],u.loc[0,False]
print("unmatched direction OR (lnc up odds / pc up odds)",round((a*dd)/(b*c),3))
sh.to_csv(O+"lncrna_recomputed_direction.tsv",sep="\t")
