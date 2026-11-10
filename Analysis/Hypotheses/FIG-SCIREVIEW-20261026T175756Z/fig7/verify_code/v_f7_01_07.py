"""Verifier: independent recomputation of 7F study-mean contrasts (F7-01, F7-07). Reads one small TSV."""
import pandas as pd
E = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/"
sm = pd.read_csv(E + "study_metrics.tsv", sep="\t")
print(sm.groupby(["block_id","head_id"]).size().to_string())
f1 = sm.set_index(["block_id","head_id","held_out_study"])["donor_class_balanced_macro_f1"]
don = sm.drop_duplicates("held_out_study").set_index("held_out_study")["donors"]
print(don.to_string(), don.sum())
clean = {"transcriptformer_tf_sapiens":["GSE136103","GSE174748","GSE189600","GSE202379","GSE244832"],
         "uce_4l":["GSE136103","GSE174748","GSE189600","GSE202379","GSE244832"],
         "uce_33l":["GSE136103","GSE174748","GSE189600","GSE202379","GSE244832"],
         "scimilarity_v1_1":["GSE174748","GSE189600","GSE202379","GSE244832","Liver_Atlas"]}
all7 = list(don.index)
out=[]
for ref in ["hvg_pca_task_native","scvi_liver_latent"]:
  for m in ["geneformer_v2_316m","transcriptformer_tf_sapiens","scimilarity_v1_1","uce_4l","uce_33l","scvi_liver_latent"]:
    if m==ref: continue
    for h in ["linear","two_layer_mlp"]:
      for nm, st in [("all7", all7)] + ([("clean", clean[m])] if m in clean else []):
        d = pd.Series({s: f1[(m,h,s)]-f1[(ref,h,s)] for s in st})
        w = don[st]
        out.append(dict(ref=ref, model=m, head=h, set=nm, n_studies=len(st), unweighted=d.mean(), donor_weighted=(d*w).sum()/w.sum()))
o = pd.DataFrame(out)
print(o.round(5).to_string(index=False))
o.to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig7/verify_out/v_f7_01_07.tsv", sep="\t", index=False)
# per-study GSE189600 values
for m in ["uce_4l","uce_33l","hvg_pca_task_native","scvi_liver_latent","transcriptformer_tf_sapiens"]:
  print(m, [round(f1[(m,h,"GSE189600")],4) for h in ["linear","two_layer_mlp"]])
