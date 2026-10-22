"""R08-09: per-study deltas in the Fig 7F grid: GSE189600 vs the other six studies, by model-by-head block;
equal-study vs donor-weighted mean and sign flips."""
import numpy as np, pandas as pd
G="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/tables/panel_7f_exposure_grid.tsv"
OUT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F7b"
g=pd.read_csv(G,sep="\t"); rows=[]
for m,d in g.groupby("model_id",sort=False):
    o=d[d.held_out_study!="GSE189600"]; x=d[d.held_out_study=="GSE189600"].per_study_delta.iloc[0]
    eq=d.per_study_delta.mean(); dw=np.average(d.per_study_delta,weights=d.donors)
    rows.append(dict(model_id=m,n_studies=len(d),GSE189600_delta=x,other6_min=o.per_study_delta.min(),other6_max=o.per_study_delta.max(),
        other6_n_positive=int((o.per_study_delta>0).sum()),equal_study_mean=eq,donor_weighted_mean=dw,sign_flip=bool(np.sign(round(eq,4))!=np.sign(round(dw,4)))))
t=pd.DataFrame(rows); t.to_csv(f"{OUT}/fig7f_per_study_delta_ranges.tsv",sep="\t",index=False,float_format="%.6f")
pd.set_option("display.width",250); print(t.round(4).to_string())
big=t.GSE189600_delta<-0.4
print("GSE189600 large blocks",t[big].model_id.tolist(),t[big].GSE189600_delta.min(),t[big].GSE189600_delta.max())
print("GSE189600 other blocks range",t[~big].GSE189600_delta.min(),t[~big].GSE189600_delta.max())
print("other six studies overall",t.other6_min.min(),t.other6_max.max(),"positive",t.other6_n_positive.sum(),"of",6*len(t))
print("sign flips (4 dp)",t.sign_flip.sum())
