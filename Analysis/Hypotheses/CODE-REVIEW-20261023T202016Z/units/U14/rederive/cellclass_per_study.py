"""U14: per-study donor-class-balanced macro F1 for TranscriptFormer and scVI liver on the five TF-clean studies.

Explains why the pooled clean-study difference (Figure 7F) and the equal-study mean (METHODS.md:481) differ in sign.
Same metric and predictions as cellclass_7c7f.py (3-seed mean probability, argmax).

usage: cellclass_per_study.py OUTDIR
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
SEL = EX / "model-data-050-20260824T211318Z-21077049/fixture/selection.tsv"
SHARDS = EX / "encoder-field-common-head-21313780/shards"
CLASSES = ["cholangiocyte", "endothelial", "hepatocyte", "immune", "mesenchymal_stromal"]
PCOLS = [f"probability::{c}" for c in CLASSES]
TF_CLEAN = ["GSE136103", "GSE174748", "GSE189600", "GSE202379", "GSE244832"]

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
sel = pd.read_csv(SEL, sep="\t", usecols=["row_id", "donor_id", "dataset", "broad_label"]).set_index("row_id")


def load(model, head):
    mats = []
    for seed in (1103, 1201, 1301):
        parts = [pd.read_csv(SHARDS / f"{model}__{head}__seed{seed}__fold{f}/predictions.tsv", sep="\t",
                             usecols=["row_id"] + PCOLS).set_index("row_id") for f in range(5)]
        mats.append(pd.concat(parts).loc[sel.index, PCOLS].to_numpy())
    return np.array(CLASSES)[np.mean(np.stack(mats), axis=0).argmax(1)]


def dcb_macro_f1(df):
    """Every class equal mass; every donor with that class equal share of it; cells split their donor-class share."""
    n = df.groupby(["broad_label", "donor_id"]).size()
    nd = df.groupby("broad_label")["donor_id"].nunique()
    w = 1.0 / (nd.loc[df["broad_label"]].to_numpy() * n.loc[list(zip(df["broad_label"], df["donor_id"]))].to_numpy())
    f1 = []
    for c in sorted(df["broad_label"].unique()):
        tp = w[(df["broad_label"] == c) & (df["pred"] == c)].sum()
        fn = w[(df["broad_label"] == c) & (df["pred"] != c)].sum()
        fp = w[(df["broad_label"] != c) & (df["pred"] == c)].sum()
        f1.append(2 * tp / (2 * tp + fp + fn))
    return float(np.mean(f1))


rows = []
for head in ["linear", "two_layer_mlp"]:
    pred = {m: load(m, head) for m in ["transcriptformer_tf_sapiens", "scvi_liver_latent"]}
    for s in TF_CLEAN + ["pooled_clean5"]:
        keep = sel["dataset"].isin(TF_CLEAN if s == "pooled_clean5" else [s]).to_numpy()
        r = {"head": head, "study": s, "donors": sel.loc[keep, "donor_id"].nunique(), "cells": int(keep.sum())}
        for m, p in pred.items():
            r[m] = dcb_macro_f1(sel.loc[keep].assign(pred=p[keep]))
        r["tf_minus_scvi"] = r["transcriptformer_tf_sapiens"] - r["scvi_liver_latent"]
        rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(out / "tf_minus_scvi_per_study.tsv", sep="\t", index=False, float_format="%.6f")
print(t.to_string())
for head in ["linear", "two_layer_mlp"]:
    x = t[(t["head"] == head) & (t["study"] != "pooled_clean5")]
    print(head, "equal-study mean of differences", round(x["tf_minus_scvi"].mean(), 6))
