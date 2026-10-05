#!/usr/bin/env python3
"""Independent recomputation of the C6 observed statistics (no new test, no null).
Refits the frozen image model as c6 does, saves the TCGA image predictions, and
recomputes each partial Spearman with OLS residuals of scipy ranks and pearsonr.
Usage: python c6_verify.py <gtex_c2_exec_dir> <tcga_exec_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr, rankdata

import c4_ridge as c4
import c5_sealed_test as c5

c2, c6 = Path(sys.argv[1]), Path(sys.argv[2])
out = c6 / "c6_verify"; out.mkdir()
dev = pd.read_csv(c2 / "targets/targets_development.tsv", sep="\t").merge(pd.read_csv(c2 / "embed/slide_features.tsv.gz", sep="\t"), on="SUBJID")
sea = pd.read_csv(c2 / "targets/sealed/targets_sealed.tsv", sep="\t").merge(pd.read_csv(c2 / "embed/sealed/slide_features.tsv.gz", sep="\t"), on="SUBJID")
tc = pd.read_csv(c6 / "targets/targets_tcga.tsv", sep="\t").merge(pd.read_csv(c6 / "embed/slide_features.tsv.gz", sep="\t"), on="SUBJID")
emb = [c for c in dev.columns if c.startswith("emb_")]
neg = ["neg_hypoxia", "neg_immediate_early"]
Xb_dev = c4.baseline(pd.concat([dev, sea], ignore_index=True)).to_numpy(float)[: len(dev)]
cov = sm.add_constant(np.column_stack([rankdata(tc["fat_fraction"]), rankdata(tc["age"]), (tc["gender"] == "male").astype(float)]))
rows, preds = [], {"SUBJID": tc["SUBJID"]}
for version in ("primary", "residualized"):
    saved = pd.read_csv(c6 / f"c6_external_test/external_{version}.tsv", sep="\t").set_index("target")["partial_spearman"]
    for t in saved.index:
        yd, yt = dev[t].to_numpy(float), tc[t].to_numpy(float)
        if version == "residualized":
            yd = sm.OLS(yd, sm.add_constant(dev[neg].to_numpy(float))).fit().resid
            yt = sm.OLS(yt, sm.add_constant(tc[neg].to_numpy(float))).fit().resid
        _, (m1, pca) = c5.fit_predict(Xb_dev, dev[emb].to_numpy(float), yd[:, None], Xb_dev[:1], dev[emb].to_numpy(float)[:1], 256)
        p = m1.predict(pca.transform(tc[emb].to_numpy(float))).ravel()
        preds[f"{version}:{t}"] = p
        r = pearsonr(sm.OLS(rankdata(p), cov).fit().resid, sm.OLS(rankdata(yt), cov).fit().resid)[0]
        rows.append({"version": version, "target": t, "recomputed": r, "saved": saved[t], "abs_diff": abs(r - saved[t])})
tab = pd.DataFrame(rows); tab.to_csv(out / "recomputed_partial_spearman.tsv", sep="\t", index=False)
pd.DataFrame(preds).to_csv(out / "tcga_image_predictions.tsv", sep="\t", index=False)
print(tab.groupby("version")["recomputed"].median().to_dict(), "max abs diff", tab["abs_diff"].max())
