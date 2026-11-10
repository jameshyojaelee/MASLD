#!/usr/bin/env python
"""CGB8 post hoc sensitivities (DEVIATIONS.md, 2026-10-23), outside the prespecified decision.
S1: primary G vs D model with the Broadaway eQTL count (F5) added as a covariate.
S2: G redefined by COLOC ABF PP4 >= 0.5 (r2 gene-level table) in the same H57 universe."""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cgb8_gene_type import FEATS, SEED, logit, matched, prepare  # noqa: E402

PRIMARY = ["F1_loeuf", "F2_abc_enh", "F3_length", "F4_density"]
ap = argparse.ArgumentParser()
ap.add_argument("--genes", required=True)
ap.add_argument("--eqtl", required=True)
ap.add_argument("--coloc", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=False)
g = pd.read_csv(a.genes, sep="\t")
res = {}

# S1: eQTL count as covariate
u = prepare(g.copy(), a.eqtl)
s1 = logit(u, "G", "D", PRIMARY + ["F5_eqtl_poscontrol"])
s1 = s1[s1.feature.isin(PRIMARY)].copy()
s1["holm_p"] = multipletests(s1.p, method="holm")[1]
s1.to_csv(f"{a.out}/S1_G_vs_D_eqtl_adjusted.tsv", sep="\t", index=False)
res["S1_eqtl_adjusted"] = s1.to_dict(orient="records")

# S2: ABF colocalization set
co = pd.read_csv(a.coloc)
co["g"] = co.ensembl.astype(str).str.split(".").str[0]
abf = co.groupby("g").coloc_best_abf_pp4.max()
abf_genes = set(abf.index[pd.to_numeric(abf, errors="coerce") >= 0.5])
g2 = g.copy()
g2["coloc_gen"] = g2.g.isin(abf_genes)
u2 = prepare(g2, a.eqtl)
res["S2_abf_set_sizes"] = u2.set.value_counts().to_dict()
res["S2_abf_genes_total"] = len(abf_genes)
pc = logit(u2, "G", "D", ["F5_eqtl_poscontrol"])
res["S2_positive_control_F5"] = pc.iloc[0].to_dict()
s2 = logit(u2, "G", "D", PRIMARY)
s2["holm_p"] = multipletests(s2.p, method="holm")[1]
s2.to_csv(f"{a.out}/S2_G_vs_D_abf.tsv", sep="\t", index=False)
res["S2_abf"] = s2.to_dict(orient="records")
rng = np.random.default_rng(SEED)
mt = pd.concat([matched(u2.dropna(subset=list(FEATS) + ["AveExpr"]).reset_index(drop=True), s, rng) for s in ["G", "D"]])
mt.to_csv(f"{a.out}/S2_matched_check.tsv", sep="\t", index=False)
res["S2_matched_loeuf"] = mt[mt.feature == "F1_loeuf"].to_dict(orient="records")
json.dump(res, open(f"{a.out}/results.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
