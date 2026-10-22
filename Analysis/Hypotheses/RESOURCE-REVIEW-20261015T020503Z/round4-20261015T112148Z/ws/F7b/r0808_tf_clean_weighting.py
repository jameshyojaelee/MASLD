"""R08-08: TranscriptFormer minus liver encoder on the five clean studies under three study weightings,
point estimates from study_metrics.tsv and intervals from the original saved 10,000 paired draws (seed 20260824 in producer).
No new resampling."""
import numpy as np, pandas as pd
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EVAL=f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation"
G=f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/tables/panel_7f_exposure_grid.tsv"
OUT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F7b"
sm=pd.read_csv(f"{EVAL}/study_metrics.tsv",sep="\t")
pt={(r.model_id,r.held_out_study):r.donor_class_balanced_macro_f1 for r in sm.itertuples()}
grid=pd.read_csv(G,sep="\t")
don=grid.drop_duplicates("held_out_study").set_index("held_out_study").donors
clean=['GSE136103','GSE174748','GSE189600','GSE202379','GSE244832']
schemes={"equal-study mean, 5 clean studies (published)":(clean,"equal"),
         "donor-weighted mean, 5 clean studies":(clean,"donor"),
         "equal-study mean, 4 clean studies without GSE189600":([s for s in clean if s!="GSE189600"],"equal")}
rows=[]
with np.load(f"{EVAL}/bootstrap_distributions.npz") as dr:
  for head in ["linear","two_layer_mlp"]:
    m=f"transcriptformer_tf_sapiens::{head}"; ref=f"scvi_liver_latent::{head}"
    for name,(ss,wt) in schemes.items():
        w=np.array([don[s] if wt=="donor" else 1.0 for s in ss]); w=w/w.sum()
        est=float(np.sum([wi*(pt[m,s]-pt[ref,s]) for wi,s in zip(w,ss)]))
        v=np.sum([wi*(dr[f"{m}::study_macro_f1::{s}"]-dr[f"{ref}::study_macro_f1::{s}"]) for wi,s in zip(w,ss)],axis=0)
        lo,hi=np.quantile(v,[.025,.975])
        rows.append(dict(model_id=m,reference_id=ref,weighting=name,studies=",".join(ss),donors=int(sum(don[s] for s in ss)),
                         delta=est,lo95=lo,hi95=hi,excludes_zero=bool(lo>0 or hi<0),resamples=len(v)))
    for s in clean:
        rows.append(dict(model_id=m,reference_id=ref,weighting=f"single study {s}",studies=s,donors=int(don[s]),delta=pt[m,s]-pt[ref,s],lo95=np.nan,hi95=np.nan,excludes_zero=np.nan,resamples=np.nan))
t=pd.DataFrame(rows); t.to_csv(f"{OUT}/fig7_transcriptformer_clean_subset_weightings.tsv",sep="\t",index=False,float_format="%.6f")
pd.set_option("display.width",250); print(t.round(6).to_string())
