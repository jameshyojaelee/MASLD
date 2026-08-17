#!/usr/bin/env python
"""Falsifiability (minimum detectable effect) layer for the MASLD Gene Catalog.

Attaches an 80%-power MDE to every indeterminate / untestable evidence call in
two validated candidate families:

  A. ATAC context v3 differential accessibility  (jointly testable peaks, Wald)
  B. ATAC context v3 prespecified program projection (117 programs x 2 cohorts, Wald)
  C. Spatial two-program confirmatory family (matched-gene permutation null, Moran)
  D. Spatial 117-program x 7-assay coverage family (observability only -> MDE undefined)

CANDIDATE-ONLY. Reads are read-only; all writes go to this candidate directory.
No promotion, no modification of upstream candidates.
"""

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

SEED = 20260813
np.random.seed(SEED)

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CAND = ROOT / "Analysis/Multimodal_Program_Projection/candidates"
ATAC = CAND / "atac-context-v3-candidate-2026-08-11-r1"
SPAT = CAND / "spatial-resource-candidate-2026-08-11"
SPIMP = CAND / "spatial-impact-figures-candidate-2026-08-11"

OUT = Path(os.environ["MDE_OUT"])
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

Z80 = norm.ppf(0.80)  # 0.8416212335729143

INPUTS = [
    ATAC / "da/da_peak_results.tsv.gz",
    ATAC / "da/da_lineage_summary.tsv",
    ATAC / "da/primary_peak_filter_qc.tsv",
    ATAC / "programs/program_atac_results.tsv",
    SPAT / "effects/spatial_program_effects.tsv",
    SPAT / "coverage/spatial_program_coverage.tsv",
    SPIMP / "data/figS_spatial_matched_null_calibration.tsv",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_implied_alpha(pvals, q=0.05):
    """Effective per-test alpha actually applied by BH within one family.

    If BH rejects at least one hypothesis, the operative per-test alpha is the
    largest p-value that was rejected. If BH rejects nothing, no p-value in the
    family cleared the ladder, and the least stringent rung any test could have
    cleared is rank 1, i.e. q/m. Returns (alpha_eff, n_rejected, m).
    """
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p)]
    m = p.size
    if m == 0:
        return np.nan, 0, 0
    order = np.sort(p)
    ranks = np.arange(1, m + 1)
    passed = order <= (ranks / m) * q
    if passed.any():
        kmax = ranks[passed].max()
        return float(order[kmax - 1]), int(kmax), int(m)
    return float(q / m), 0, int(m)


def mde_wald(se, alpha, power_z=Z80, sided=2):
    z = norm.ppf(1 - alpha / sided)
    return (z + power_z) * np.asarray(se, dtype=float)


# ---------------------------------------------------------------------------
# A. ATAC differential accessibility
# ---------------------------------------------------------------------------
da = pd.read_csv(ATAC / "da/da_peak_results.tsv.gz", sep="\t")
COH = ["gse244832", "gse281367"]

# BH families as fitted upstream: one per lineage x cohort.
alpha_rows = []
alpha_map = {}
for lin, g in da.groupby("lineage"):
    for c in COH:
        a, nrej, m = bh_implied_alpha(g[f"pvalue_{c}"], q=0.05)
        alpha_map[(lin, c)] = a
        alpha_rows.append(
            dict(modality="atac", family=f"da_peak::{lin}::{c}", lineage=lin,
                 cohort=c, n_tests=m, n_bh_rejected=nrej,
                 alpha_eff_bh=a, alpha_nominal=0.05)
        )

for c in COH:
    da[f"mde_bh_{c}"] = np.nan
    da[f"alpha_eff_{c}"] = np.nan
for (lin, c), a in alpha_map.items():
    sel = da["lineage"] == lin
    da.loc[sel, f"alpha_eff_{c}"] = a
    da.loc[sel, f"mde_bh_{c}"] = mde_wald(da.loc[sel, f"standard_error_{c}"], a)
for c in COH:
    da[f"mde_nominal_{c}"] = mde_wald(da[f"standard_error_{c}"], 0.05)

# A `supported` peak must clear the bar in BOTH cohorts, so the operative MDE
# for the supported state is the weaker (larger) of the two per-cohort MDEs.
da["mde_joint_bh"] = da[[f"mde_bh_{c}" for c in COH]].max(axis=1)
da["mde_joint_nominal"] = da[[f"mde_nominal_{c}" for c in COH]].max(axis=1)

# --- prespecified meaningful effect, derived from the supported calls --------
sup = da[da["evidence_state"] == "supported"]
sup_abs = np.concatenate([sup[f"logFC_{c}"].abs().values for c in COH])
delta_atac_peak = float(np.percentile(sup_abs, 10))
delta_atac_peak_minarm = float(
    np.percentile(np.minimum(sup["logFC_gse244832"].abs(), sup["logFC_gse281367"].abs()), 10)
)

def classify(row, delta, mdecol):
    st = row["evidence_state"]
    if st != "indeterminate":
        return st, f"not_eligible_state_{st}"
    m = row[mdecol]
    if not np.isfinite(m):
        return "indeterminate", "mde_not_computable_missing_se"
    if m < delta:
        return "tested_negative", "mde_below_prespecified_meaningful_effect"
    return "indeterminate", "mde_above_prespecified_meaningful_effect"

res = da.apply(lambda r: classify(r, delta_atac_peak, "mde_joint_bh"), axis=1)
da["proposed_state_bh"] = [x[0] for x in res]
da["proposed_reason_bh"] = [x[1] for x in res]
res = da.apply(lambda r: classify(r, delta_atac_peak, "mde_joint_nominal"), axis=1)
da["proposed_state_nominal"] = [x[0] for x in res]

da_out = da[[
    "lineage", "peak_coordinate", "evidence_state", "promoter_genes",
    "logFC_gse244832", "standard_error_gse244832", "pvalue_gse244832", "qvalue_gse244832",
    "logFC_gse281367", "standard_error_gse281367", "pvalue_gse281367", "qvalue_gse281367",
    "alpha_eff_gse244832", "alpha_eff_gse281367",
    "mde_bh_gse244832", "mde_bh_gse281367", "mde_joint_bh",
    "mde_nominal_gse244832", "mde_nominal_gse281367", "mde_joint_nominal",
    "proposed_state_bh", "proposed_reason_bh", "proposed_state_nominal",
]].copy()
da_out.insert(0, "unit_id", da["lineage"] + "::" + da["peak_coordinate"])
da_out.insert(0, "family", "atac_da_peak::" + da["lineage"])
da_out.insert(0, "modality", "atac_accessibility")
da_out["prespecified_meaningful_effect"] = delta_atac_peak
da_out["meaningful_effect_rule"] = "p10_of_abs_logFC_across_both_cohort_arms_of_the_4_supported_peaks"
da_out.to_csv(TAB / "mde_atac_da_peaks.tsv.gz", sep="\t", index=False)

# ---------------------------------------------------------------------------
# B. ATAC program projection (117 programs x 2 cohorts)
# ---------------------------------------------------------------------------
pr = pd.read_csv(ATAC / "programs/program_atac_results.tsv", sep="\t")
pr_alpha = []
pr["alpha_eff_bh"] = np.nan
for (coh, lin), g in pr.groupby(["cohort", "lineage"]):
    a, nrej, m = bh_implied_alpha(g["pvalue"], q=0.05)
    pr.loc[g.index, "alpha_eff_bh"] = a
    pr_alpha.append(dict(modality="atac", family=f"program::{coh}::{lin}",
                         lineage=lin, cohort=coh, n_tests=m, n_bh_rejected=nrej,
                         alpha_eff_bh=a, alpha_nominal=0.05))
alpha_rows += pr_alpha

pr["mde_bh"] = mde_wald(pr["standard_error"], 1.0) * np.nan  # placeholder
pr["mde_bh"] = [
    mde_wald(se, a) if np.isfinite(se) and np.isfinite(a) else np.nan
    for se, a in zip(pr["standard_error"], pr["alpha_eff_bh"])
]
pr["mde_nominal"] = mde_wald(pr["standard_error"], 0.05)

n_supported_pr = int((pr["cross_cohort_state"] == "supported").sum())
if n_supported_pr > 0:
    delta_pr = float(np.percentile(
        pr.loc[pr["cross_cohort_state"] == "supported", "effect"].abs(), 10))
    pr_rule = "p10_of_abs_effect_among_supported_program_rows"
else:
    delta_pr = np.nan
    pr_rule = ("undefined_zero_supported_rows_in_family__no_within_modality_anchor__"
               "no_reclassification_permitted")


def classify_pr(row):
    st = row["cross_cohort_state"]
    if st == "untestable":
        if not np.isfinite(row["standard_error"]):
            return "untestable", "no_se_stored__mde_not_computable"
        return "untestable", "mde_computable_for_the_tested_cohort_arm_only__state_unchanged"
    if st != "indeterminate":
        return st, f"not_eligible_state_{st}"
    if not np.isfinite(delta_pr):
        return "indeterminate", ("mde_reported__no_prespecified_effect_derivable__"
                                 "zero_supported_calls_in_this_modality_family")
    if np.isfinite(row["mde_bh"]) and row["mde_bh"] < delta_pr:
        return "tested_negative", "mde_below_prespecified_meaningful_effect"
    return "indeterminate", "mde_above_prespecified_meaningful_effect"


r = pr.apply(classify_pr, axis=1)
pr["proposed_state_bh"] = [x[0] for x in r]
pr["proposed_reason_bh"] = [x[1] for x in r]
pr_out = pr[[
    "cohort", "lineage", "program_uid", "program_name", "n_normal", "n_mash",
    "n_measured_genes", "program_score_testable", "contrast_testable",
    "lineage_observability_state", "testability_reason",
    "effect", "standard_error", "pvalue", "qvalue",
    "within_source_state", "cross_cohort_state",
    "alpha_eff_bh", "mde_bh", "mde_nominal",
    "proposed_state_bh", "proposed_reason_bh",
]].copy()
pr_out.insert(0, "unit_id", pr["cohort"] + "::" + pr["lineage"] + "::" + pr["program_uid"])
pr_out.insert(0, "family", "atac_program::" + pr["cohort"] + "::" + pr["lineage"])
pr_out.insert(0, "modality", "atac_program_projection")
pr_out.rename(columns={"cross_cohort_state": "current_state"}, inplace=True)
pr_out["prespecified_meaningful_effect"] = delta_pr
pr_out["meaningful_effect_rule"] = pr_rule
pr_out.to_csv(TAB / "mde_atac_program_projection.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# C. Spatial two-program confirmatory family (permutation null, Moran scale)
# ---------------------------------------------------------------------------
sp = pd.read_csv(SPAT / "effects/spatial_program_effects.tsv", sep="\t")
cal = pd.read_csv(SPIMP / "data/figS_spatial_matched_null_calibration.tsv", sep="\t")

# permutation resolution floor: smallest attainable one-sided p with n_null draws
sp["perm_p_floor"] = 1.0 / (sp["n_null_draws"] + 1.0)

sp_alpha = []
sp["alpha_eff_bh"] = np.nan
for (ds, fam), g in sp.groupby(["dataset_id", "multiplicity_family"]):
    a, nrej, m = bh_implied_alpha(g["pvalue"], q=0.05)
    floor = g["perm_p_floor"].min()
    if np.isfinite(a) and np.isfinite(floor):
        a = max(a, floor)
    sp.loc[g.index, "alpha_eff_bh"] = a
    sp_alpha.append(dict(modality="spatial", family=f"spatial::{ds}::{fam}",
                         lineage="", cohort=ds, n_tests=m, n_bh_rejected=nrej,
                         alpha_eff_bh=a, alpha_nominal=0.05))
alpha_rows += sp_alpha

# One-sided upper-tail test against the matched-gene null.
def mde_perm(sd, alpha):
    if not np.isfinite(sd) or not np.isfinite(alpha):
        return np.nan
    return (norm.ppf(1 - alpha) + Z80) * sd

sp["mde_bh"] = [mde_perm(sd, a) for sd, a in zip(sp["matched_null_sd"], sp["alpha_eff_bh"])]
sp["mde_nominal"] = [mde_perm(sd, 0.05) for sd in sp["matched_null_sd"]]

sup_sp = sp[sp["evidence_state"].isin(["supported", "source_dependent"]) &
            sp["matched_null_sd"].notna() &
            (sp["within_source_result"] == "supported")]
delta_sp = float(np.percentile(sup_sp["estimate"].abs(), 10)) if len(sup_sp) else np.nan

DONOR_UNIT_OK = {"donor"}


def classify_sp(row):
    st = row["evidence_state"]
    if st != "indeterminate":
        return st, f"not_eligible_state_{st}"
    if not np.isfinite(row["matched_null_sd"]):
        return "indeterminate", "no_matched_null_sd_stored__mde_not_computable"
    if row["biological_unit"] not in DONOR_UNIT_OK:
        return "indeterminate", (
            f"unit_not_donor__{row['biological_unit']}__"
            "mde_reported_but_reclassification_withheld")
    if np.isfinite(delta_sp) and row["mde_bh"] < delta_sp:
        return "tested_negative", "mde_below_prespecified_meaningful_effect"
    return "indeterminate", "mde_above_prespecified_meaningful_effect"


r = sp.apply(classify_sp, axis=1)
sp["proposed_state_bh"] = [x[0] for x in r]
sp["proposed_reason_bh"] = [x[1] for x in r]
sp_out = sp[[
    "dataset_id", "assay_id", "program_uid", "program_label", "biological_unit",
    "n_biological", "technical_unit", "n_technical", "source_dependence",
    "estimate", "matched_null_sd", "pvalue", "qvalue", "n_null_draws",
    "multiplicity_family", "within_source_result", "evidence_state",
    "perm_p_floor", "alpha_eff_bh", "mde_bh", "mde_nominal",
    "proposed_state_bh", "proposed_reason_bh",
]].copy()
sp_out.insert(0, "unit_id", sp["dataset_id"] + "::" + sp["assay_id"] + "::" + sp["program_uid"])
sp_out.insert(0, "family", "spatial_two_program::" + sp["dataset_id"])
sp_out.insert(0, "modality", "spatial_residual_moran")
sp_out.rename(columns={"evidence_state": "current_state"}, inplace=True)
sp_out["prespecified_meaningful_effect"] = delta_sp
sp_out["meaningful_effect_rule"] = (
    "p10_of_abs_centered_residual_moran_i_among_within_source_supported_rows")
sp_out.to_csv(TAB / "mde_spatial_two_program.tsv", sep="\t", index=False)

# D. spatial coverage family: observability, no effect estimate -> MDE undefined
cov = pd.read_csv(SPAT / "coverage/spatial_program_coverage.tsv", sep="\t")
cov_out = cov[["dataset_id", "assay_id", "program_uid", "program_label", "cell_type",
               "n_genes_measured", "retained_l1_weight", "coverage_status",
               "testability_reason", "biological_unit", "source_dependence"]].copy()
cov_out.insert(0, "modality", "spatial_coverage")
cov_out.insert(1, "family", "spatial_program_by_assay_coverage::" + cov["dataset_id"] + "::" + cov["assay_id"])
cov_out.insert(2, "unit_id", cov["dataset_id"] + "::" + cov["assay_id"] + "::" + cov["program_uid"])
cov_out["standard_error_or_null_sd"] = np.nan
cov_out["alpha_eff"] = np.nan
cov_out["mde"] = np.nan
cov_out["proposed_state"] = cov_out["coverage_status"]
cov_out["proposed_reason"] = (
    "coverage_family_reports_observability_not_an_effect_estimate__"
    "no_test_statistic_stored__mde_undefined")
cov_out.to_csv(TAB / "mde_spatial_coverage_family.tsv.gz", sep="\t", index=False)

pd.DataFrame(alpha_rows).to_csv(TAB / "alpha_eff_by_family.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------
def qs(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict(n=0, median=np.nan, q25=np.nan, q75=np.nan, p90=np.nan, min=np.nan, max=np.nan)
    return dict(n=int(x.size), median=float(np.median(x)),
                q25=float(np.percentile(x, 25)), q75=float(np.percentile(x, 75)),
                p90=float(np.percentile(x, 90)), min=float(x.min()), max=float(x.max()))

summ = []
for lin, g in da.groupby("lineage"):
    ind = g[g["evidence_state"] == "indeterminate"]
    for label, col in [("joint_bh", "mde_joint_bh"), ("joint_nominal", "mde_joint_nominal"),
                       ("gse244832_bh", "mde_bh_gse244832"), ("gse281367_bh", "mde_bh_gse281367")]:
        s = qs(ind[col])
        s.update(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
                 arm=label, n_indeterminate=int(len(ind)),
                 n_proposed_tested_negative=int((ind["proposed_state_bh"] == "tested_negative").sum())
                 if label == "joint_bh" else
                 int((ind["proposed_state_nominal"] == "tested_negative").sum())
                 if label == "joint_nominal" else np.nan,
                 prespecified_meaningful_effect=delta_atac_peak)
        summ.append(s)
for (coh, lin), g in pr.groupby(["cohort", "lineage"]):
    ind = g[g["cross_cohort_state"] == "indeterminate"]
    s = qs(ind["mde_bh"])
    s.update(modality="atac_program_projection", family=f"atac_program::{coh}::{lin}",
             arm="bh", n_indeterminate=int(len(ind)),
             n_proposed_tested_negative=int((ind["proposed_state_bh"] == "tested_negative").sum()),
             prespecified_meaningful_effect=delta_pr)
    summ.append(s)
for ds, g in sp.groupby("dataset_id"):
    ind = g[g["evidence_state"] == "indeterminate"]
    s = qs(ind["mde_bh"])
    s.update(modality="spatial_residual_moran", family=f"spatial_two_program::{ds}",
             arm="bh", n_indeterminate=int(len(ind)),
             n_proposed_tested_negative=int((ind["proposed_state_bh"] == "tested_negative").sum()),
             prespecified_meaningful_effect=delta_sp)
    summ.append(s)
pd.DataFrame(summ).to_csv(TAB / "mde_summary_by_family.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Donors-needed curve (Wald families only)
#   SE scales as sqrt(1/n1 + 1/n2) at fixed dispersion. Rescale each peak's
#   observed SE to a balanced n-per-arm design and recount adequate power.
# ---------------------------------------------------------------------------
NPAIR = {("hepatocyte", "gse244832"): (5, 9), ("hepatocyte", "gse281367"): (6, 6),
         ("stellate", "gse244832"): (5, 9), ("stellate", "gse281367"): (6, 6),
         ("macrophage", "gse244832"): (4, 8), ("macrophage", "gse281367"): (6, 6)}
curve = []
NGRID = [4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200]
for lin, g in da.groupby("lineage"):
    ind = g[g["evidence_state"] == "indeterminate"]
    if len(ind) == 0:
        continue
    for n in NGRID:
        scaled = []
        for c in COH:
            n1, n2 = NPAIR[(lin, c)]
            fac = np.sqrt((2.0 / n) / (1.0 / n1 + 1.0 / n2))
            scaled.append(mde_wald(ind[f"standard_error_{c}"] * fac, alpha_map[(lin, c)]))
        joint = np.nanmax(np.vstack(scaled), axis=0)
        curve.append(dict(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
                          n_donors_per_arm_per_cohort=n,
                          n_indeterminate=int(len(ind)),
                          median_mde=float(np.nanmedian(joint)),
                          n_adequately_powered=int(np.nansum(joint < delta_atac_peak)),
                          frac_adequately_powered=float(np.nanmean(joint < delta_atac_peak)),
                          prespecified_meaningful_effect=delta_atac_peak,
                          note="alpha_eff held at the observed BH-implied value; "
                               "dispersion held fixed; balanced n per arm"))
pd.DataFrame(curve).to_csv(TAB / "donors_needed_curve_atac_da.tsv", sep="\t", index=False)

# per-peak donors needed to reach delta
dn = []
for lin, g in da.groupby("lineage"):
    ind = g[g["evidence_state"] == "indeterminate"]
    if len(ind) == 0:
        continue
    need = []
    for c in COH:
        n1, n2 = NPAIR[(lin, c)]
        cur = mde_wald(ind[f"standard_error_{c}"], alpha_map[(lin, c)])
        ratio = (cur / delta_atac_peak) ** 2
        need.append(2.0 * ratio / (1.0 / n1 + 1.0 / n2))
    need = np.nanmax(np.vstack(need), axis=0)
    s = qs(need)
    s.update(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
             metric="donors_per_arm_per_cohort_needed_for_80pct_power_at_delta")
    dn.append(s)
pd.DataFrame(dn).to_csv(TAB / "donors_needed_per_peak_summary.tsv", sep="\t", index=False)

# power floor: smallest effect that could ever have been called, per family
floor = []
for lin, g in da.groupby("lineage"):
    row = dict(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
               n_tests=int(len(g)))
    for c in COH:
        n1, n2 = NPAIR[(lin, c)]
        row[f"n_normal_{c}"] = n1
        row[f"n_mash_{c}"] = n2
        row[f"alpha_eff_{c}"] = alpha_map[(lin, c)]
        row[f"min_se_{c}"] = float(np.nanmin(g[f"standard_error_{c}"]))
        row[f"power_floor_mde_{c}"] = float(np.nanmin(mde_wald(g[f"standard_error_{c}"], alpha_map[(lin, c)])))
    row["power_floor_joint"] = float(np.nanmin(g["mde_joint_bh"]))
    row["prespecified_meaningful_effect"] = delta_atac_peak
    floor.append(row)
for (coh, lin), g in pr.groupby(["cohort", "lineage"]):
    if g["standard_error"].notna().sum() == 0:
        continue
    floor.append(dict(modality="atac_program_projection", family=f"atac_program::{coh}::{lin}",
                      n_tests=int(g["pvalue"].notna().sum()),
                      n_normal_gse244832=int(g["n_normal"].iloc[0]),
                      n_mash_gse244832=int(g["n_mash"].iloc[0]),
                      alpha_eff_gse244832=float(g["alpha_eff_bh"].iloc[0]),
                      min_se_gse244832=float(np.nanmin(g["standard_error"])),
                      power_floor_mde_gse244832=float(np.nanmin(g["mde_bh"])),
                      power_floor_joint=float(np.nanmin(g["mde_bh"])),
                      prespecified_meaningful_effect=delta_pr))
for ds, g in sp.groupby("dataset_id"):
    if g["matched_null_sd"].notna().sum() == 0:
        continue
    floor.append(dict(modality="spatial_residual_moran", family=f"spatial_two_program::{ds}",
                      n_tests=int(g["pvalue"].notna().sum()),
                      n_normal_gse244832=np.nan, n_mash_gse244832=np.nan,
                      alpha_eff_gse244832=float(g["alpha_eff_bh"].iloc[0]),
                      min_se_gse244832=float(np.nanmin(g["matched_null_sd"])),
                      power_floor_mde_gse244832=float(np.nanmin(g["mde_bh"])),
                      power_floor_joint=float(np.nanmin(g["mde_bh"])),
                      prespecified_meaningful_effect=delta_sp))
pd.DataFrame(floor).to_csv(TAB / "power_floor_by_family.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Reclassification headline counts
# ---------------------------------------------------------------------------
head = dict(
    seed=SEED,
    atac_peaks_total=int(len(da)),
    atac_peaks_indeterminate=int((da["evidence_state"] == "indeterminate").sum()),
    atac_peaks_source_dependent=int((da["evidence_state"] == "source_dependent").sum()),
    atac_peaks_supported=int((da["evidence_state"] == "supported").sum()),
    atac_peaks_to_tested_negative_bh=int(((da["evidence_state"] == "indeterminate") &
                                          (da["proposed_state_bh"] == "tested_negative")).sum()),
    atac_peaks_to_tested_negative_nominal=int(((da["evidence_state"] == "indeterminate") &
                                               (da["proposed_state_nominal"] == "tested_negative")).sum()),
    atac_peak_delta=delta_atac_peak,
    atac_peak_delta_min_arm_variant=delta_atac_peak_minarm,
    atac_program_rows_total=int(len(pr)),
    atac_program_indeterminate=int((pr["cross_cohort_state"] == "indeterminate").sum()),
    atac_program_untestable=int((pr["cross_cohort_state"] == "untestable").sum()),
    atac_program_with_se=int(pr["standard_error"].notna().sum()),
    atac_program_supported=n_supported_pr,
    atac_program_delta=delta_pr,
    atac_program_to_tested_negative=int((pr["proposed_state_bh"] == "tested_negative").sum()),
    spatial_rows_total=int(len(sp)),
    spatial_indeterminate=int((sp["evidence_state"] == "indeterminate").sum()),
    spatial_with_matched_null_sd=int(sp["matched_null_sd"].notna().sum()),
    spatial_delta=delta_sp,
    spatial_to_tested_negative=int((sp["proposed_state_bh"] == "tested_negative").sum()),
    spatial_coverage_rows=int(len(cov)),
    spatial_coverage_mde_defined=0,
)
with open(OUT / "tables/headline_counts.json", "w") as fh:
    json.dump(head, fh, indent=2, default=float)

man = [dict(role="input", path=str(p), sha256=sha256(p), bytes=p.stat().st_size)
       for p in INPUTS]
for p in sorted(TAB.glob("*")):
    man.append(dict(role="output", path=str(p), sha256=sha256(p), bytes=p.stat().st_size))
pd.DataFrame(man).to_csv(OUT / "checksum_manifest.tsv", sep="\t", index=False)

print(json.dumps(head, indent=2, default=float))
print("\nalpha_eff by family:")
print(pd.DataFrame(alpha_rows).to_string(index=False))
print("\nsummary:")
print(pd.DataFrame(summ).to_string(index=False))
print("\ndonors-needed curve:")
print(pd.DataFrame(curve).to_string(index=False))
print("\ndonors needed per peak:")
print(pd.DataFrame(dn).to_string(index=False))
print("\npower floor:")
print(pd.DataFrame(floor).to_string(index=False))
print("\nspatial rows:")
print(sp_out[["unit_id", "current_state", "matched_null_sd", "alpha_eff_bh",
              "mde_bh", "prespecified_meaningful_effect", "proposed_state_bh",
              "proposed_reason_bh"]].to_string(index=False))
