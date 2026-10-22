"""Recompute the Figure 7D severity-transfer block: released instrument Spearman against biopsy stage,
95% participant (or identity) bootstrap interval, and the tie-structure ceiling of each cohort's stage coding.

Inputs (analysis outputs):
  release/masld-severity-v1.2/expected_ood.json                              (which arm is the released instrument)
  executions/release-instrument-20260901T211017Z/sealed/jobC_gse268273_predictions.tsv
  executions/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv
  executions/external-validation-expansion-20260901T190000Z/GSE213621_scores.tsv, GSE276114_scores.tsv
  Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv (fibrotic_stage)
  figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv (genetic identity)
  executions/etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv (stage, etiology)
Tie ceiling: Spearman between the midranked stage and any continuous ordering consistent with stage
(the within-tie order does not change the value). Bootstrap: percentile, B = 10,000, seed 9874, unit resampled; also stage-stratified (stage counts fixed).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
EX = ROOT / "Analysis/MASLD_Model_Benchmark/executions"
OUT = Path(sys.argv[1])
rng = np.random.default_rng(9874)
B = 10000


def sp(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean(); rb -= rb.mean()
    return float((ra * rb).sum() / np.sqrt((ra ** 2).sum() * (rb ** 2).sum()))


def ceiling(y):
    y = np.asarray(y, float)
    order = np.argsort(y, kind="stable")
    r = np.empty(len(y)); r[order] = np.arange(1, len(y) + 1)
    return sp(y, r)


def boot(x, y):
    n = len(x); out = np.empty(B)
    for i in range(B):
        s = rng.integers(0, n, n)
        out[i] = sp(x[s], y[s])
    return np.quantile(out, [0.025, 0.975])


def boot_strat(x, y):
    """Units resampled within stage strata (stage counts held fixed)."""
    strata = [np.nonzero(y == u)[0] for u in np.unique(y)]
    out = np.empty(B)
    for i in range(B):
        s = np.concatenate([st[rng.integers(0, len(st), len(st))] for st in strata])
        out[i] = sp(x[s], y[s])
    return np.quantile(out, [0.025, 0.975])


rows = []


def add(cohort, unit, x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    rho = sp(x, y); c = ceiling(y); lo, hi = boot(x, y); slo, shi = boot_strat(x, y)
    rows.append({"cohort": cohort, "unit": unit, "n": len(x), "spearman": rho, "ci_lo": lo, "ci_hi": hi,
                 "ci_strat_lo": slo, "ci_strat_hi": shi,
                 "tie_ceiling": c, "ratio": rho / c, "stage_levels": int(len(np.unique(y)))})


exp = json.load(open(ROOT / "Analysis/MASLD_Model_Benchmark/release/masld-severity-v1.2/expected_ood.json"))
arm = exp["arm"]
ph = pd.read_csv(EX / "gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv", sep="\t")
pr = pd.read_csv(EX / "release-instrument-20260901T211017Z/sealed/jobC_gse268273_predictions.tsv", sep="\t")
g268 = ph[["row_id", "patient_code", "fibrosis_stage"]].merge(pr[["row_id", arm]], on="row_id")
assert g268.patient_code.is_unique
add("GSE268273", "participant record", g268[arm], g268.fibrosis_stage)

s213 = pd.read_csv(EX / "external-validation-expansion-20260901T190000Z/GSE213621_scores.tsv", sep="\t")
m213 = pd.read_csv(ROOT / "Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv", sep="\t")[["sample_id", "fibrotic_stage"]]
ro = pd.read_csv(ROOT / "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv",
                 sep="\t", low_memory=False)
ro = ro[ro.dataset == "GSE213621"][["sample_id", "analysis_unit_id", "primary_selected"]]
d213 = s213.merge(m213, on="sample_id").merge(ro, on="sample_id", how="left")
code3 = {"F0F1": 0, "F2": 1, "F3F4": 2}
d3 = d213[d213.fibrotic_stage.isin(code3)].copy(); d3["y"] = d3.fibrotic_stage.map(code3)
add("GSE213621", "library record (published set)", d3.latent_severity, d3.y)
nst = d3.groupby("analysis_unit_id")["y"].nunique()
disc = nst[nst > 1].index
keep = d3[~d3.analysis_unit_id.isin(disc)]
prim = keep[keep.primary_selected == True]
assert prim.analysis_unit_id.is_unique
add("GSE213621", "genetic identity, primary library, discordant removed", prim.latent_severity, prim.y)
mean_id = keep.groupby("analysis_unit_id").agg(x=("latent_severity", "mean"), y=("y", "first"))
add("GSE213621", "genetic identity, mean of libraries, discordant removed", mean_id.x, mean_id.y)
print("GSE213621 3-level records", len(d3), "identities", d3.analysis_unit_id.nunique(), "discordant identities", len(disc))

s276 = pd.read_csv(EX / "external-validation-expansion-20260901T190000Z/GSE276114_scores.tsv", sep="\t")
dz = pd.read_csv(EX / "etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv", sep="\t")
d276 = s276.merge(dz[["title", "disease", "stage"]], left_on="sample_id", right_on="title")
assert len(d276) == 177
for dis in ["MASLD", "CVH"]:
    sub = d276[d276.disease == dis]
    add("GSE276114", f"participant, {dis}", sub.latent_severity, sub.stage)
add("GSE276114", "participant, all etiologies (published set)", d276.latent_severity, d276.stage)

t = pd.DataFrame(rows)
t.to_csv(OUT / "severity_transfer_recomputed.tsv", sep="\t", index=False, float_format="%.6f")
print(t.to_string())
print("expected_ood record", exp["ood_spearman"], exp["tie_ceiling"])
