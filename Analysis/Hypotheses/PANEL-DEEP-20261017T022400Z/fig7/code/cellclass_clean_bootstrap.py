"""Clean-study-only gains (Figure 7F alternative): for each encoder x head, pooled donor-class-balanced macro F1
on the studies that encoder never saw, minus HVG+PCA with the same head on the same studies.
Interval: study-stratified donor-cluster bootstrap (donors resampled with replacement within each study;
a donor drawn twice counts as two donors in the weights), B = 2000, seed 9873.
Also TranscriptFormer minus scVI on the five TranscriptFormer-clean studies.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(sys.argv[1])
ns = {}
src = open(Path(__file__).with_name("recompute_cellclass.py")).read().split("blocks = [")[0]
exec(src, ns)
sel, CLASSES = ns["sel"], ns["CLASSES"]
cidx = {c: i for i, c in enumerate(CLASSES)}
obs = sel["broad_label"].map(cidx).to_numpy()
donor = pd.factorize(sel["donor_id"])[0]
study = sel["dataset"].to_numpy()
K = len(CLASSES)

models = ["geneformer_v2_316m", "transcriptformer_tf_sapiens", "scimilarity_v1_1", "uce_4l", "uce_33l", "scvi_liver_latent", "hvg_pca_task_native"]
pred = {}
for b in models:
    for h in ("linear", "two_layer_mlp"):
        pred[(b, h)] = ns["load_common"](b, h).map(cidx).to_numpy()

seen = {
    "geneformer_v2_316m": {"GSE136103", "GSE185477", "Liver_Atlas", "GSE174748", "GSE189600", "GSE202379", "GSE244832"},
    "transcriptformer_tf_sapiens": {"GSE185477", "Liver_Atlas"},
    "scimilarity_v1_1": {"GSE136103", "GSE185477"},
    "uce_4l": {"GSE185477", "Liver_Atlas"},
    "uce_33l": {"GSE185477", "Liver_Atlas"},
    "scvi_liver_latent": set(),
}
studies = sorted(set(study))
donors_by_study = {s: np.unique(donor[study == s]) for s in studies}
cells_by_donor = pd.Series(np.arange(len(donor))).groupby(donor).apply(np.array).to_dict()


def f1_macro(o, p, w):
    tp = np.bincount(o[o == p], weights=w[o == p], minlength=K)
    pc = np.bincount(p, weights=w, minlength=K)
    oc = np.bincount(o, weights=w, minlength=K)
    present = np.bincount(o, minlength=K) > 0
    f = np.where(pc + oc > 0, 2 * tp / np.where(pc + oc > 0, pc + oc, 1), 0.0)
    return float(f[present].mean())


def weights_for(inst, o):
    # inst: donor-instance id per cell; weight = 1/(k * n_instances_with_class * n_cells_in_instance_class)
    key = inst * K + o
    n_ic = np.bincount(key)[key]
    k_present = np.unique(o)
    n_inst_c = np.zeros(K)
    for c in k_present:
        n_inst_c[c] = len(np.unique(inst[o == c]))
    return 1.0 / (len(k_present) * n_inst_c[o] * n_ic)


def subset_draw(rng, keep_studies):
    idx, inst = [], []
    j = 0
    for s in keep_studies:
        for d in rng.choice(donors_by_study[s], size=len(donors_by_study[s]), replace=True):
            c = cells_by_donor[d]
            idx.append(c); inst.append(np.full(len(c), j)); j += 1
    return np.concatenate(idx), np.concatenate(inst)


rng = np.random.default_rng(9873)
B = int(__import__("os").environ.get("BOOT_B", "2000"))
tasks = []
for b, s in seen.items():
    clean = [x for x in studies if x not in s]
    if clean:
        tasks.append((b, tuple(clean)))
tasks.append(("transcriptformer_tf_sapiens", tuple(studies)))  # positive control: all-study interval vs producer draws
tf_clean = tuple(x for x in studies if x not in seen["transcriptformer_tf_sapiens"])
rows = []
for b, clean in tasks:
    mask = np.isin(study, clean)
    base_idx = np.nonzero(mask)[0]
    inst0 = donor[base_idx]
    w0 = weights_for(inst0, obs[base_idx])
    for h in ("linear", "two_layer_mlp"):
        pt = f1_macro(obs[base_idx], pred[(b, h)][base_idx], w0) - f1_macro(obs[base_idx], pred[("hvg_pca_task_native", h)][base_idx], w0)
        dr = np.empty(B)
        r2 = np.random.default_rng(9873)  # same draws for every model x head
        for i in range(B):
            idx, inst = subset_draw(r2, clean)
            w = weights_for(inst, obs[idx])
            dr[i] = f1_macro(obs[idx], pred[(b, h)][idx], w) - f1_macro(obs[idx], pred[("hvg_pca_task_native", h)][idx], w)
        rows.append({"model": b, "head": h, "clean_studies": len(clean), "clean_donors": int(len(np.unique(inst0))),
                     "clean_cells": int(len(base_idx)), "gain_clean_pooled": pt,
                     "lo": np.quantile(dr, 0.025), "hi": np.quantile(dr, 0.975)})
        print(rows[-1], flush=True)
# TranscriptFormer minus scVI on the five TranscriptFormer-clean studies
mask = np.isin(study, tf_clean); base_idx = np.nonzero(mask)[0]; w0 = weights_for(donor[base_idx], obs[base_idx])
for h in ("linear", "two_layer_mlp"):
    pt = f1_macro(obs[base_idx], pred[("transcriptformer_tf_sapiens", h)][base_idx], w0) - f1_macro(obs[base_idx], pred[("scvi_liver_latent", h)][base_idx], w0)
    dr = np.empty(B); r2 = np.random.default_rng(9873)
    for i in range(B):
        idx, inst = subset_draw(r2, tf_clean); w = weights_for(inst, obs[idx])
        dr[i] = f1_macro(obs[idx], pred[("transcriptformer_tf_sapiens", h)][idx], w) - f1_macro(obs[idx], pred[("scvi_liver_latent", h)][idx], w)
    rows.append({"model": "transcriptformer_minus_scvi", "head": h, "clean_studies": 5, "clean_donors": int(len(np.unique(donor[base_idx]))),
                 "clean_cells": int(len(base_idx)), "gain_clean_pooled": pt, "lo": np.quantile(dr, 0.025), "hi": np.quantile(dr, 0.975)})
    print(rows[-1], flush=True)
# positive control: all-study gain point with this f1 implementation must equal the recompute table
allidx = np.arange(len(obs)); wall = weights_for(donor, obs)
for h in ("linear", "two_layer_mlp"):
    print("control all-study TF gain", h, f1_macro(obs, pred[("transcriptformer_tf_sapiens", h)], wall) - f1_macro(obs, pred[("hvg_pca_task_native", h)], wall))
pd.DataFrame(rows).to_csv(OUT / "cellclass_clean_bootstrap.tsv", sep="\t", index=False, float_format="%.6f")
