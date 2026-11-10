"""U14 independent re-derivation of the Figure 7C/7F TranscriptFormer numbers.

Written from METHODS.md "Single-cell encoders", the 7C/7F legends and the evaluator's metric definition
(evaluation_receipt.json: "Every class receives equal mass, every donor containing that class receives equal
within-class mass, and cells divide their donor-class share"; study-balanced = arithmetic mean of study values).

Prediction per cell: argmax of the mean class probability over the three screen seeds (1103, 1201, 1301),
read from the out-of-fold prediction shards.
Bootstrap: donors resampled with replacement within each held-out study (a donor drawn twice counts twice);
a draw is redrawn for a study until every class observed in that study is present (class coverage kept).
The same multiplicities are used for every model and head.

usage: cellclass_7c7f.py OUTDIR B SEED [B SEED ...]
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
SEL = EX / "model-data-050-20260824T211318Z-21077049/fixture/selection.tsv"
SHARDS = EX / "encoder-field-common-head-21313780/shards"
SAVED = EX / "encoder-field-score-21314196/evaluation/bootstrap_distributions.npz"
CLASSES = ["cholangiocyte", "endothelial", "hepatocyte", "immune", "mesenchymal_stromal"]
PCOLS = [f"probability::{c}" for c in CLASSES]
MODELS = ["transcriptformer_tf_sapiens", "scvi_liver_latent", "hvg_pca_task_native"]
HEADS = ["linear", "two_layer_mlp"]
TF_CLEAN = ["GSE136103", "GSE174748", "GSE189600", "GSE202379", "GSE244832"]   # TF seen: GSE185477, Liver_Atlas

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
runs = [(int(sys.argv[i]), int(sys.argv[i + 1])) for i in range(2, len(sys.argv), 2)]

sel = pd.read_csv(SEL, sep="\t", usecols=["row_id", "donor_id", "dataset", "broad_label"]).set_index("row_id")
assert len(sel) == 50000 and sel.index.is_unique
K = len(CLASSES)
obs = sel["broad_label"].map({c: i for i, c in enumerate(CLASSES)}).to_numpy()
assert not np.isnan(obs.astype(float)).any()
studies = sorted(sel["dataset"].unique())
print("studies", studies, "donors", sel["donor_id"].nunique(), "cells", len(sel))
assert sel.groupby("donor_id")["dataset"].nunique().max() == 1, "a donor id spans two studies"

donors = sorted(sel["donor_id"].unique())
didx = sel["donor_id"].map({d: i for i, d in enumerate(donors)}).to_numpy()
donor_study = sel.groupby("donor_id")["dataset"].first().loc[donors].to_numpy()
D = len(donors)


def load(model, head):
    mats = []
    for seed in (1103, 1201, 1301):
        parts = [pd.read_csv(SHARDS / f"{model}__{head}__seed{seed}__fold{f}/predictions.tsv", sep="\t",
                             usecols=["row_id"] + PCOLS).set_index("row_id") for f in range(5)]
        m = pd.concat(parts)
        assert m.index.is_unique and len(m) == 50000
        mats.append(m.loc[sel.index, PCOLS].to_numpy())
    return np.mean(np.stack(mats), axis=0).argmax(1)


# per-donor confusion counts, then row-normalised (each donor-class row sums to 1)
R = {}
pres = np.zeros((D, K))
np.add.at(pres, (didx, obs), 1)
ndc = pres.copy()
pres = (pres > 0).astype(float)
for m in MODELS:
    for h in HEADS:
        p = load(m, h)
        C = np.zeros((D, K, K))
        np.add.at(C, (didx, obs, p), 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            Rn = np.where(ndc[:, :, None] > 0, C / ndc[:, :, None], 0.0)
        R[(m, h)] = Rn.reshape(D, K * K)
        print("loaded", m, h, flush=True)


def f1_from(Mult, Rflat):
    """Mult: (B, D) donor multiplicities. Returns (B,) donor-class-balanced macro F1."""
    A = (Mult @ Rflat).reshape(-1, K, K)
    N = Mult @ pres
    ok = N > 0
    A = np.where(ok[:, :, None], A / np.where(ok, N, 1)[:, :, None], 0.0)
    tp = np.einsum("bkk->bk", A)
    fn = A.sum(2) - tp
    fp = A.sum(1) - tp
    den = 2 * tp + fp + fn
    f1 = np.where(den > 0, 2 * tp / np.where(den > 0, den, 1), 0.0)
    return (f1 * ok).sum(1) / ok.sum(1)


def subset_mult(Mult, study_set):
    keep = np.isin(donor_study, study_set).astype(float)
    return Mult * keep


def metrics(Mult):
    res = {}
    for (m, h), Rf in R.items():
        res[(m, h, "pooled_all7")] = f1_from(Mult, Rf)
        res[(m, h, "pooled_tfclean5")] = f1_from(subset_mult(Mult, TF_CLEAN), Rf)
        per = np.stack([f1_from(subset_mult(Mult, [s]), Rf) for s in TF_CLEAN], 1)
        res[(m, h, "studymean_tfclean5")] = per.mean(1)
        per7 = np.stack([f1_from(subset_mult(Mult, [s]), Rf) for s in studies], 1)
        res[(m, h, "studymean_all7")] = per7.mean(1)
    return res


point = metrics(np.ones((1, D)))
pt = {f"{m}::{h}::{k}": float(v[0]) for (m, h, k), v in point.items()}
print(json.dumps(pt, indent=1))

donors_by_study = {s: np.flatnonzero(donor_study == s) for s in studies}
classes_by_study = {s: set(np.flatnonzero(pres[donors_by_study[s]].sum(0) > 0)) for s in studies}


def draw(B, rng):
    Mult = np.zeros((B, D))
    redraws = {}
    for s in studies:
        ids = donors_by_study[s]
        need = np.array(sorted(classes_by_study[s]))
        pick = rng.integers(0, len(ids), (B, len(ids)))
        nre = 0
        while True:
            cov = np.stack([pres[ids[pick[:, j]]] for j in range(len(ids))], 1).sum(1)[:, need] > 0
            bad = ~cov.all(1)
            if not bad.any():
                break
            nre += int(bad.sum())
            pick[bad] = rng.integers(0, len(ids), (int(bad.sum()), len(ids)))
        redraws[s] = nre
        for j in range(len(ids)):
            np.add.at(Mult, (np.arange(B), ids[pick[:, j]]), 1)
    return Mult, redraws


CONTRASTS = [  # (label, model a, model b, head, metric)
    ("7C TF all7 score", "transcriptformer_tf_sapiens", None, "two_layer_mlp", "pooled_all7"),
    ("7C TF all7 score", "transcriptformer_tf_sapiens", None, "linear", "pooled_all7"),
    ("7C TF clean score", "transcriptformer_tf_sapiens", None, "two_layer_mlp", "pooled_tfclean5"),
    ("7C TF clean score", "transcriptformer_tf_sapiens", None, "linear", "pooled_tfclean5"),
    ("7F TF minus scVI, pooled clean", "transcriptformer_tf_sapiens", "scvi_liver_latent", "linear", "pooled_tfclean5"),
    ("7F TF minus scVI, pooled clean", "transcriptformer_tf_sapiens", "scvi_liver_latent", "two_layer_mlp", "pooled_tfclean5"),
    ("Methods TF minus scVI, study mean clean", "transcriptformer_tf_sapiens", "scvi_liver_latent", "linear", "studymean_tfclean5"),
    ("Methods TF minus scVI, study mean clean", "transcriptformer_tf_sapiens", "scvi_liver_latent", "two_layer_mlp", "studymean_tfclean5"),
    ("7F TF gain over HVG+PCA, pooled clean", "transcriptformer_tf_sapiens", "hvg_pca_task_native", "linear", "pooled_tfclean5"),
    ("7F TF gain over HVG+PCA, pooled clean", "transcriptformer_tf_sapiens", "hvg_pca_task_native", "two_layer_mlp", "pooled_tfclean5"),
    ("7F TF gain over HVG+PCA, pooled all7", "transcriptformer_tf_sapiens", "hvg_pca_task_native", "linear", "pooled_all7"),
    ("7F TF gain over HVG+PCA, pooled all7", "transcriptformer_tf_sapiens", "hvg_pca_task_native", "two_layer_mlp", "pooled_all7"),
]


def contrast_vals(res, a, b, h, k):
    v = res[(a, h, k)]
    return v - res[(b, h, k)] if b else v


rows = []
meta = {"inputs": {}, "runs": [], "unit": "donor within held-out study", "tf_clean_studies": TF_CLEAN}
for p in [SEL]:
    meta["inputs"][str(p)] = hashlib.sha256(open(p, "rb").read()).hexdigest()
for B, seed in runs:
    rng = np.random.default_rng(seed)
    Mult, redraws = draw(B, rng)
    res = metrics(Mult)
    for lab, a, b, h, k in CONTRASTS:
        v = contrast_vals(res, a, b, h, k)
        rows.append({"label": lab, "a": a, "b": b or "", "head": h, "metric": k, "B": B, "seed": seed,
                     "point": float(contrast_vals(point, a, b, h, k)[0]),
                     "lo95": np.percentile(v, 2.5), "hi95": np.percentile(v, 97.5), "frac_le0": float(np.mean(v <= 0))})
    meta["runs"].append({"B": B, "seed": seed, "redraws_for_class_coverage": redraws})
    print("done", B, seed, redraws, flush=True)

# cross-check: re-aggregate the evaluator's own saved 10,000 study-level draws (seed 20260824)
z = np.load(SAVED, allow_pickle=False)
for h in HEADS:
    a = np.mean([z[f"transcriptformer_tf_sapiens::{h}::study_macro_f1::{s}"] for s in TF_CLEAN], 0)
    b = np.mean([z[f"scvi_liver_latent::{h}::study_macro_f1::{s}"] for s in TF_CLEAN], 0)
    v = a - b
    rows.append({"label": "saved-draw reaggregation, TF minus scVI, study mean clean", "a": "transcriptformer_tf_sapiens",
                 "b": "scvi_liver_latent", "head": h, "metric": "studymean_tfclean5", "B": len(v), "seed": 20260824,
                 "point": float(contrast_vals(point, "transcriptformer_tf_sapiens", "scvi_liver_latent", h, "studymean_tfclean5")[0]),
                 "lo95": np.percentile(v, 2.5), "hi95": np.percentile(v, 97.5), "frac_le0": float(np.mean(v <= 0))})
    a7 = z[f"transcriptformer_tf_sapiens::{h}::donor_class_balanced_macro_f1"]
    rows.append({"label": "saved-draw TF all7 pooled score", "a": "transcriptformer_tf_sapiens", "b": "", "head": h,
                 "metric": "pooled_all7", "B": len(a7), "seed": 20260824,
                 "point": float(point[("transcriptformer_tf_sapiens", h, "pooled_all7")][0]),
                 "lo95": np.percentile(a7, 2.5), "hi95": np.percentile(a7, 97.5), "frac_le0": float(np.mean(a7 <= 0))})

t = pd.DataFrame(rows)
t.to_csv(out / "cellclass_7c7f_rederived.tsv", sep="\t", index=False, float_format="%.6f")
meta["point"] = pt
json.dump(meta, open(out / "cellclass_7c7f_meta.json", "w"), indent=1)
print(t.to_string())
