"""Independent recomputation of broad cell-class mapping scores (Figure 7C bottom block, 7F gains).

Inputs (analysis outputs, not figure tables):
  model-data-050-20260824T211318Z-21077049/fixture/selection.tsv           (row_id, donor_id, dataset, broad_label)
  encoder-field-common-head-21313780/shards/<block>__<head>__seed<s>__fold<f>/predictions.tsv
  model-cpu-train-503-21097161/predictions/predictions/hvg_pca_<clf>.tsv     (three-seed ensemble files)
  encoder-field-score-21314196/evaluation/bootstrap_distributions.npz      (producer's paired donor-bootstrap draws)

Metric: donor-class-balanced macro F1. Each cell weight = 1 / (n_classes * n_donors_with_class * n_cells_in_donor_class);
per-class weighted F1, unweighted mean over classes. Prediction = argmax of the mean probability over three seeds.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
OUT = Path(sys.argv[1])
CLASSES = ["cholangiocyte", "endothelial", "hepatocyte", "immune", "mesenchymal_stromal"]
PCOLS = [f"probability::{c}" for c in CLASSES]

sel = pd.read_csv(EX / "model-data-050-20260824T211318Z-21077049/fixture/selection.tsv", sep="\t",
                  usecols=["row_id", "donor_id", "dataset", "broad_label"]).set_index("row_id")
assert len(sel) == 50000


def weights(donor, label):
    df = pd.DataFrame({"d": donor, "l": label})
    n_dc = df.groupby(["d", "l"])["d"].transform("size")
    n_donors_c = df.groupby("l")["d"].transform("nunique")
    k = df["l"].nunique()
    return (1.0 / (k * n_donors_c * n_dc)).to_numpy()


def wmacro_f1(obs, pred, w):
    scores = []
    for c in sorted(set(obs)):
        tp = w[(obs == c) & (pred == c)].sum()
        fp = w[(obs != c) & (pred == c)].sum()
        fn = w[(obs == c) & (pred != c)].sum()
        den = 2 * tp + fp + fn
        scores.append(0.0 if den == 0 else 2 * tp / den)
    return float(np.mean(scores))


def score(pred_series, subset=None):
    s = sel if subset is None else sel[sel["dataset"].isin(subset)]
    p = pred_series.loc[s.index].to_numpy()
    o = s["broad_label"].to_numpy()
    return wmacro_f1(o, p, weights(s["donor_id"].to_numpy(), o))


def load_common(block, head):
    mats = []
    for seed in (1103, 1201, 1301):
        parts = []
        for fold in range(5):
            f = EX / f"encoder-field-common-head-21313780/shards/{block}__{head}__seed{seed}__fold{fold}/predictions.tsv"
            parts.append(pd.read_csv(f, sep="\t", usecols=["row_id"] + PCOLS).set_index("row_id"))
        m = pd.concat(parts)
        assert m.index.is_unique and len(m) == 50000
        mats.append(m.loc[sel.index, PCOLS].to_numpy())
    ens = np.mean(np.stack(mats), axis=0)
    return pd.Series(np.array(CLASSES)[ens.argmax(1)], index=sel.index)


def load_hvg(clf):
    m = pd.read_csv(EX / f"model-cpu-train-503-21097161/predictions/predictions/hvg_pca_{clf}.tsv", sep="\t",
                    usecols=["row_id"] + PCOLS).set_index("row_id")
    seeds = []
    for s in (20260824, 20260825, 20260826):
        x = pd.read_csv(EX / f"model-cpu-train-503-21097161/predictions/predictions/hvg_pca_{clf}--seed-{s}.tsv",
                        sep="\t", usecols=["row_id"] + PCOLS).set_index("row_id")
        seeds.append(x.loc[sel.index, PCOLS].to_numpy())
    ens = np.mean(np.stack(seeds), axis=0)
    ens_file = m.loc[sel.index, PCOLS].to_numpy()
    return pd.Series(np.array(CLASSES)[ens.argmax(1)], index=sel.index), float(np.abs(ens - ens_file).max())


blocks = ["geneformer_v2_316m", "transcriptformer_tf_sapiens", "scimilarity_v1_1", "uce_4l", "uce_33l",
          "scvi_liver_latent", "hvg_pca_task_native"]
studies = sorted(sel["dataset"].unique())
rows, preds = [], {}
for b in blocks:
    for h in ("linear", "two_layer_mlp"):
        p = load_common(b, h)
        preds[f"{b}::{h}"] = p
        r = {"model": f"{b}::{h}", "all_studies": score(p)}
        for st in studies:
            r[st] = score(p, [st])
        rows.append(r)
for clf in ("knn", "linear_svm", "logistic", "elastic_net", "nearest_centroid"):
    p, maxdiff = load_hvg(clf)
    preds[f"hvg_pca_{clf}"] = p
    r = {"model": f"hvg_pca_{clf}", "all_studies": score(p), "seed_mean_vs_ensemble_file_maxdiff": maxdiff}
    for st in studies:
        r[st] = score(p, [st])
    rows.append(r)
sc = pd.DataFrame(rows)
sc.to_csv(OUT / "cellclass_scores_recomputed.tsv", sep="\t", index=False, float_format="%.6f")
print(sc[["model", "all_studies"]].to_string())

# Gains over the task-native HVG+PCA common head, all studies (7F) and the clean-only alternative
expo = {  # per-study exposure state used by the current 7F grid (Oct 14 correction for TranscriptFormer)
    "geneformer_v2_316m": {"GSE136103": "seen", "GSE185477": "seen", "Liver_Atlas": "seen"},
    "transcriptformer_tf_sapiens": {"GSE185477": "seen", "Liver_Atlas": "seen"},
    "scimilarity_v1_1": {"GSE136103": "seen", "GSE185477": "seen"},
    "uce_4l": {"GSE185477": "seen", "Liver_Atlas": "seen"},
    "uce_33l": {"GSE185477": "seen", "Liver_Atlas": "seen"},
    "scvi_liver_latent": {},
}
unknown = {"geneformer_v2_316m": {"GSE174748", "GSE189600", "GSE202379", "GSE244832"}}
draws = np.load(EX / "encoder-field-score-21314196/evaluation/bootstrap_distributions.npz")
grow = []
for b in expo:
    for h in ("linear", "two_layer_mlp"):
        m, ref = f"{b}::{h}", f"hvg_pca_task_native::{h}"
        g_all = score(preds[m]) - score(preds[ref])
        dv = draws[f"{m}::donor_class_balanced_macro_f1"] - draws[f"{ref}::donor_class_balanced_macro_f1"]
        clean = [s for s in studies if s not in expo[b] and s not in unknown.get(b, set())]
        g_clean = (score(preds[m], clean) - score(preds[ref], clean)) if clean else np.nan
        # unweighted study mean on clean studies (draws available per study)
        if clean:
            dsm = np.mean([draws[f"{m}::study_macro_f1::{s}"] - draws[f"{ref}::study_macro_f1::{s}"] for s in clean], axis=0)
            g_clean_sm = float(np.mean([score(preds[m], [s]) - score(preds[ref], [s]) for s in clean]))
            lo_sm, hi_sm = np.quantile(dsm, [0.025, 0.975])
        else:
            g_clean_sm = lo_sm = hi_sm = np.nan
        grow.append({"model": m, "gain_all_point": g_all, "gain_all_lo": np.quantile(dv, 0.025),
                     "gain_all_hi": np.quantile(dv, 0.975), "n_clean_studies": len(clean),
                     "gain_clean_pooled_point": g_clean, "gain_clean_study_mean": g_clean_sm,
                     "gain_clean_study_mean_lo": lo_sm, "gain_clean_study_mean_hi": hi_sm,
                     "vs_scvi_all": score(preds[m]) - score(preds[f"scvi_liver_latent::{h}"])})
gt = pd.DataFrame(grow)
gt.to_csv(OUT / "cellclass_gains_recomputed.tsv", sep="\t", index=False, float_format="%.6f")
print(gt.to_string())
# donor and cell counts per study
cnt = sel.groupby("dataset").agg(cells=("donor_id", "size"), donors=("donor_id", "nunique"))
cnt.to_csv(OUT / "cellclass_study_counts.tsv", sep="\t")
print(cnt, sel["donor_id"].nunique())
