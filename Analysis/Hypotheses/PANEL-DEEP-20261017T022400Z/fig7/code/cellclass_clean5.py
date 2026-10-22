"""TranscriptFormer and other encoders on the five TranscriptFormer-clean studies, pooled donor-class-balanced
macro F1, reusing the recompute functions (imports recompute_cellclass logic without re-running the gain table)."""
import sys, runpy
import numpy as np, pandas as pd
sys.argv = [sys.argv[0], sys.argv[1]]
ns = {}
src = open(sys.argv[0].replace("cellclass_clean5.py", "recompute_cellclass.py")).read()
src = src.split("blocks = [")[0]  # functions and data only
exec(src, ns)
clean = ["GSE136103", "GSE174748", "GSE189600", "GSE202379", "GSE244832"]
rows = []
for h in ("linear", "two_layer_mlp"):
    p = {b: ns["load_common"](b, h) for b in ["transcriptformer_tf_sapiens", "scvi_liver_latent", "hvg_pca_task_native", "uce_4l", "scimilarity_v1_1", "uce_33l"]}
    for b, v in p.items():
        rows.append({"head": h, "model": b, "clean5_pooled": ns["score"](v, clean), "all7": ns["score"](v)})
t = pd.DataFrame(rows)
t.to_csv(sys.argv[1] + "/cellclass_clean5_scores.tsv", sep="\t", index=False, float_format="%.6f")
print(t.to_string())
sel = ns["sel"]; print("clean5 donors", sel[sel.dataset.isin(clean)].donor_id.nunique(), "cells", sel.dataset.isin(clean).sum())
