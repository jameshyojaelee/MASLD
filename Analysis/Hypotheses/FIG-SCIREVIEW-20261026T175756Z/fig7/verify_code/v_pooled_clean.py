"""Verifier: pooled donor-class-balanced macro F1 on study subsets (7C/7F), written independently.
Validates against model_metrics.tsv all-7 values, then scores clean subsets and SCimilarity with GSE136103 counted clean.
Reads 4 models x 2 heads x 15 prediction shards (about 4 MB each)."""
import numpy as np, pandas as pd
EX = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/"
CL = ["cholangiocyte", "endothelial", "hepatocyte", "immune", "mesenchymal_stromal"]
P = ["probability::" + c for c in CL]
sel = pd.read_csv(EX + "model-data-050-20260824T211318Z-21077049/fixture/selection.tsv", sep="\t",
                  usecols=["row_id", "donor_id", "dataset", "broad_label"]).set_index("row_id")
mm = pd.read_csv(EX + "encoder-field-score-21314196/evaluation/model_metrics.tsv", sep="\t").set_index("model_id")

def pred(block, head):
    acc = None
    for s in (1103, 1201, 1301):
        parts = [pd.read_csv(f"{EX}encoder-field-common-head-21313780/shards/{block}__{head}__seed{s}__fold{f}/predictions.tsv",
                             sep="\t", usecols=["row_id"] + P, index_col="row_id") for f in range(5)]
        m = pd.concat(parts).loc[sel.index, P].to_numpy()
        acc = m if acc is None else acc + m
    return np.array(CL)[acc.argmax(1)]

def f1(obs, prd, donor):
    d = pd.DataFrame({"o": obs, "d": donor})
    cell_n = d.groupby(["o", "d"])["o"].transform("size").to_numpy()
    don_n = d.groupby("o")["d"].transform("nunique").to_numpy()
    K = d["o"].nunique()
    w = 1.0 / (K * don_n * cell_n)
    out = []
    for c in np.unique(obs):
        tp = w[(obs == c) & (prd == c)].sum(); fp = w[(obs != c) & (prd == c)].sum(); fn = w[(obs == c) & (prd != c)].sum()
        out.append(2 * tp / (2 * tp + fp + fn))
    return float(np.mean(out))

def score(p, studies=None):
    m = np.ones(len(sel), bool) if studies is None else sel["dataset"].isin(studies).to_numpy()
    return f1(sel["broad_label"].to_numpy()[m], p[m], sel["donor_id"].to_numpy()[m])

TF5 = ["GSE136103", "GSE174748", "GSE189600", "GSE202379", "GSE244832"]
SC5 = ["GSE174748", "GSE189600", "GSE202379", "GSE244832", "Liver_Atlas"]
SC6 = SC5 + ["GSE136103"]
rows = []
preds = {}
for b in ["transcriptformer_tf_sapiens", "scvi_liver_latent", "scimilarity_v1_1", "hvg_pca_task_native"]:
    for h in ["linear", "two_layer_mlp"]:
        p = pred(b, h); preds[(b, h)] = p
        rows.append(dict(model=b, head=h, all7=score(p), all7_published=mm.loc[f"{b}::{h}", "donor_class_balanced_macro_f1"],
                         TF5=score(p, TF5), SC5=score(p, SC5), SC6=score(p, SC6)))
t = pd.DataFrame(rows); print(t.round(6).to_string(index=False))
for h in ["linear", "two_layer_mlp"]:
    g = lambda a, b, s: t.query("model==@a and head==@h")[s].item() - t.query("model==@b and head==@h")[s].item()
    print(h, "TF-scVI TF5", round(g("transcriptformer_tf_sapiens", "scvi_liver_latent", "TF5"), 6),
          "| SC-HVG SC5", round(g("scimilarity_v1_1", "hvg_pca_task_native", "SC5"), 6), "SC6", round(g("scimilarity_v1_1", "hvg_pca_task_native", "SC6"), 6),
          "| SC-scVI SC5", round(g("scimilarity_v1_1", "scvi_liver_latent", "SC5"), 6), "SC6", round(g("scimilarity_v1_1", "scvi_liver_latent", "SC6"), 6))
print("donors SC6", sel[sel.dataset.isin(SC6)].donor_id.nunique(), "SC5", sel[sel.dataset.isin(SC5)].donor_id.nunique(), "TF5", sel[sel.dataset.isin(TF5)].donor_id.nunique())
t.to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig7/verify_out/v_pooled_clean.tsv", sep="\t", index=False)
