"""Fig 7F: drawn pooled metric vs unweighted and donor-weighted study means. Light, reads one TSV."""
import pandas as pd
E = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/"
s = pd.read_csv(E + "study_metrics.tsv", sep="\t"); s = s[s.head_id.isin(["linear", "two_layer_mlp"])]
p = s.pivot_table(index=["block_id", "head_id"], columns="held_out_study", values="donor_class_balanced_macro_f1")
don = s.groupby("held_out_study")["donors"].first()
TF = ["GSE136103", "GSE174748", "GSE189600", "GSE202379", "GSE244832"]; SC = ["GSE174748", "GSE189600", "GSE202379", "GSE244832", "Liver_Atlas"]
rows = []
for ref in ["hvg_pca_task_native", "scvi_liver_latent"]:
    for (b, h), r in p.iterrows():
        if b in (ref, "hvg_pca_task_native"): continue
        d = r - p.loc[(ref, h)]
        sets = [("all7", list(p.columns))]
        if b in ("transcriptformer_tf_sapiens", "uce_4l", "uce_33l"): sets.append(("clean", TF))
        if b == "scimilarity_v1_1": sets.append(("clean", SC))
        for nm, st in sets:
            rows.append(dict(reference=ref, model=b, head=h, studies=nm, unweighted_study_mean=d[st].mean(),
                             donor_weighted_study_mean=(d[st] * don[st]).sum() / don[st].sum()))
out = pd.DataFrame(rows); print(out.round(4).to_string(index=False))
out.to_csv(__file__.replace("code/weighting_check.py", "out/weighting_check.tsv"), sep="\t", index=False)
print(p.loc[[("uce_4l","two_layer_mlp"),("uce_33l","two_layer_mlp"),("hvg_pca_task_native","two_layer_mlp")],"GSE189600"])
