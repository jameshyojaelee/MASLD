#!/usr/bin/env python
"""BF2 prespecified secondaries (PRESPEC_BF2.md, 'Secondary (descriptive)'), on the primary run's midpoints."""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

RUN, OUT = sys.argv[1], sys.argv[2]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
W = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/gene_system_weights.tsv.gz"
R9D = (f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
       "catalog_r9d-20261016T105059Z/build_out/catalog_r9d_gene_table.parquet")
ZONE = f"{ROOT}/Analysis/Hypotheses/BF1-20261023T111742Z/out-20261023T112740Z/s2_system_zone_map.tsv"

d = pd.read_csv(f"{RUN}/midpoints_point.tsv", sep="\t")
d = d[d.eligible]
m = d.pivot(index="system", columns="set", values="m")
dl = d.pivot(index="system", columns="set", values="delta")
s = pd.DataFrame({"m_avg": m[["GSE162694", "GSE213621"]].mean(1),
                  "sign_162694": np.sign(dl.GSE162694), "sign_213621": np.sign(dl.GSE213621)})
s["direction"] = np.where((s.sign_162694 > 0) & (s.sign_213621 > 0), "rising",
                          np.where((s.sign_162694 < 0) & (s.sign_213621 < 0), "falling", "disagree"))
w = pd.read_csv(W, sep="\t")
hs = pd.read_parquet(R9D, columns=["ensembl_id", "lineage_hep_specificity_log2cpm"]).drop_duplicates("ensembl_id")
w = w.merge(hs, left_on="gene_id_base", right_on="ensembl_id", how="left")
s["hep_spec_mean"] = w.groupby("community_id").lineage_hep_specificity_log2cpm.mean().reindex(s.index)
s["hep_tercile"] = pd.qcut(s.hep_spec_mean, 3, labels=["low", "mid", "high"])
z = pd.read_csv(ZONE, sep="\t").pivot(index="system", columns="layer", values="zone_bias")
s["zone_bias_protein"], s["zone_bias_mrna"] = z.reindex(s.index).protein, z.reindex(s.index).mrna
res = {"n_systems": len(s), "direction_counts": s.direction.value_counts().to_dict()}
hi, lo = s.m_avg[s.hep_tercile == "high"], s.m_avg[s.hep_tercile == "low"]
res["hep_high_vs_low"] = dict(n_high=len(hi), n_low=len(lo), median_high=hi.median(), median_low=lo.median(),
                              wilcoxon_p=stats.mannwhitneyu(hi, lo).pvalue)
ri, fa = s.m_avg[s.direction == "rising"], s.m_avg[s.direction == "falling"]
res["rising_vs_falling"] = dict(n_rising=len(ri), n_falling=len(fa), median_rising=ri.median(),
                                median_falling=fa.median(),
                                wilcoxon_p=stats.mannwhitneyu(ri, fa).pvalue if len(fa) and len(ri) else None)
for L in ["protein", "mrna"]:
    r, p = stats.spearmanr(s[f"zone_bias_{L}"], s.m_avg, nan_policy="omit")
    res[f"zone_bias_{L}_vs_m_avg"] = dict(spearman=r, p=p)
for k in ["GSE162694_within_stage", "GSE213621_within_stage", "GSE130970", "GSE135251"]:
    r, p = stats.spearmanr(m[k], s.m_avg)
    res[f"{k}_vs_m_avg"] = dict(spearman=r, p=p)
s.join(m).to_csv(f"{OUT}/secondary_systems.tsv", sep="\t")
json.dump(res, open(f"{OUT}/secondary.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
