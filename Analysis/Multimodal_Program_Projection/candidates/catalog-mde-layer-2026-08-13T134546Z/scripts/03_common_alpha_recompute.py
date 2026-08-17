#!/usr/bin/env python
"""Recompute the MDE layer at a COMMON alpha and invert primary/sensitivity.

Defect being fixed
------------------
The first pass used the BH-implied alpha, which is a function of how many
rejections each family happened to make. GSE244832 hepatocyte made 1 rejection
(alpha_eff 4.59e-7) and GSE244832 stellate made 820 (alpha_eff 4.43e-3), a
~10,000-fold difference, while their minimum standard errors are nearly equal
(0.0474 vs 0.0461). The BH arm therefore measures discovery count, not power,
and is NOT comparable across families.

New primary
-----------
alpha arms (all applied identically to every family):
  * bonferroni_family : 0.05 / m, m = family size            [PRIMARY]
  * nominal_005       : 0.05                                  [common-alpha comparator]
  * bh_implied        : as-fitted BH threshold                [SENSITIVITY, non-comparable]

meaningful-effect anchors:
  * 0.25, 0.50, 1.00 log2  (external, not selected on this study) [PRIMARY]
  * delta* = 0.796978 log2 (p10 of the 4 supported peaks)         [SENSITIVITY, circular]

Spatial evidence calls are NOT revisited here; spatial MDEs are recomputed as
numbers only and the proposed states from the first pass are left untouched.

CANDIDATE-ONLY. Writes only into the recompute subdirectory.
"""
import hashlib
import json
import os
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
BULK = (ROOT / "results/remediation/bg001/"
        "bg001-fragment-v211-gencode49-20260807T195243Z/frozen_sets/canonical_deg_results.csv")

OUT = Path(os.environ["MDE_OUT"]) / "recompute-common-alpha"
OUT.mkdir(parents=True, exist_ok=True)

Z80 = norm.ppf(0.80)
COH = ["gse244832", "gse281367"]
ANCHORS_PRIMARY = {"log2fc_0.25": 0.25, "log2fc_0.50": 0.50, "log2fc_1.00": 1.00}
DELTA_STAR = 0.796977808032176  # from pass 1; sensitivity only
ANCHORS_ALL = {**ANCHORS_PRIMARY, "delta_star_0.797_SENSITIVITY": DELTA_STAR}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def bh_implied_alpha(pvals, q=0.05):
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p)]
    m = p.size
    if m == 0:
        return np.nan, 0, 0
    o = np.sort(p)
    r = np.arange(1, m + 1)
    ok = o <= (r / m) * q
    if ok.any():
        return float(o[r[ok].max() - 1]), int(r[ok].max()), int(m)
    return float(q / m), 0, int(m)


def mde(se, alpha, sided=2):
    return (norm.ppf(1 - alpha / sided) + Z80) * np.asarray(se, dtype=float)


# ---------------------------------------------------------------- ATAC peaks
da = pd.read_csv(ATAC / "da/da_peak_results.tsv.gz", sep="\t")
alpha_tbl = []
ARMS = {}
for lin, g in da.groupby("lineage"):
    for c in COH:
        a_bh, nrej, m = bh_implied_alpha(g[f"pvalue_{c}"])
        arms = {"bonferroni_family": 0.05 / m, "nominal_005": 0.05, "bh_implied": a_bh}
        ARMS[(lin, c)] = arms
        for name, a in arms.items():
            alpha_tbl.append(dict(modality="atac_accessibility",
                                  family=f"atac_da_peak::{lin}", cohort=c,
                                  alpha_arm=name, alpha=a, family_size_m=m,
                                  n_bh_rejected=nrej,
                                  comparable_across_families=name != "bh_implied"))

for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
    for c in COH:
        da[f"mde_{arm}_{c}"] = np.nan
    for (lin, c), arms in ARMS.items():
        s = da["lineage"] == lin
        da.loc[s, f"mde_{arm}_{c}"] = mde(da.loc[s, f"standard_error_{c}"], arms[arm])
    # supported requires clearing BOTH cohorts -> weaker arm governs
    da[f"mde_{arm}_joint"] = da[[f"mde_{arm}_{c}" for c in COH]].max(axis=1)

for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
    for aname, aval in ANCHORS_ALL.items():
        col = f"state_{arm}__{aname}"
        st = np.where(da["evidence_state"] != "indeterminate", da["evidence_state"],
                      np.where(da[f"mde_{arm}_joint"] < aval,
                               "tested_negative", "indeterminate"))
        da[col] = st

keep = (["lineage", "peak_coordinate", "evidence_state", "promoter_genes"] +
        [f"logFC_{c}" for c in COH] + [f"standard_error_{c}" for c in COH] +
        [f"pvalue_{c}" for c in COH] +
        [c for c in da.columns if c.startswith("mde_") or c.startswith("state_")])
out = da[keep].copy()
out.insert(0, "unit_id", da["lineage"] + "::" + da["peak_coordinate"])
out.insert(0, "family", "atac_da_peak::" + da["lineage"])
out.insert(0, "modality", "atac_accessibility")
out["primary_alpha_arm"] = "bonferroni_family"
out["primary_anchors"] = "0.25|0.50|1.00_log2_external"
out["sensitivity_note"] = ("bh_implied alpha and delta*=0.797 are sensitivity arms; "
                           "bh_implied is not comparable across families")
out.to_csv(OUT / "mde_common_alpha_atac_da_peaks.tsv.gz", sep="\t", index=False)

# ------------------------------------------------------- ATAC program family
pr = pd.read_csv(ATAC / "programs/program_atac_results.tsv", sep="\t")
pr["alpha_bonferroni_family"] = np.nan
pr["alpha_bh_implied"] = np.nan
for (coh, lin), g in pr.groupby(["cohort", "lineage"]):
    a_bh, nrej, m = bh_implied_alpha(g["pvalue"])
    if m == 0:
        continue
    pr.loc[g.index, "alpha_bonferroni_family"] = 0.05 / m
    pr.loc[g.index, "alpha_bh_implied"] = a_bh
    for name, a in [("bonferroni_family", 0.05 / m), ("nominal_005", 0.05),
                    ("bh_implied", a_bh)]:
        alpha_tbl.append(dict(modality="atac_program_projection",
                              family=f"atac_program::{coh}::{lin}", cohort=coh,
                              alpha_arm=name, alpha=a, family_size_m=m,
                              n_bh_rejected=nrej,
                              comparable_across_families=name != "bh_implied"))
pr["mde_bonferroni_family"] = [mde(s, a) if np.isfinite(s) and np.isfinite(a) else np.nan
                               for s, a in zip(pr["standard_error"], pr["alpha_bonferroni_family"])]
pr["mde_nominal_005"] = mde(pr["standard_error"], 0.05)
pr["mde_bh_implied"] = [mde(s, a) if np.isfinite(s) and np.isfinite(a) else np.nan
                        for s, a in zip(pr["standard_error"], pr["alpha_bh_implied"])]
for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
    for aname, aval in ANCHORS_ALL.items():
        pr[f"state_{arm}__{aname}"] = np.where(
            pr["cross_cohort_state"] != "indeterminate", pr["cross_cohort_state"],
            np.where(pr[f"mde_{arm}"] < aval, "tested_negative", "indeterminate"))
prk = (["cohort", "lineage", "program_uid", "program_name", "n_normal", "n_mash",
        "effect", "standard_error", "pvalue", "qvalue", "within_source_state",
        "cross_cohort_state", "alpha_bonferroni_family", "alpha_bh_implied"] +
       [c for c in pr.columns if c.startswith("mde_") or c.startswith("state_")])
po = pr[prk].copy()
po.insert(0, "unit_id", pr["cohort"] + "::" + pr["lineage"] + "::" + pr["program_uid"])
po.insert(0, "family", "atac_program::" + pr["cohort"] + "::" + pr["lineage"])
po.insert(0, "modality", "atac_program_projection")
po["primary_alpha_arm"] = "bonferroni_family"
po.to_csv(OUT / "mde_common_alpha_atac_program.tsv", sep="\t", index=False)

# ------------------------------------------------------------------- spatial
# Numbers only. Evidence calls and proposed states from pass 1 are NOT revisited.
sp = pd.read_csv(SPAT / "effects/spatial_program_effects.tsv", sep="\t")
sp["perm_p_floor"] = 1.0 / (sp["n_null_draws"] + 1.0)
srows = []
for (ds, fam), g in sp.groupby(["dataset_id", "multiplicity_family"]):
    a_bh, nrej, m = bh_implied_alpha(g["pvalue"])
    if m == 0:
        continue
    floor = float(g["perm_p_floor"].min()) if np.isfinite(g["perm_p_floor"]).any() else np.nan
    for _, r in g.iterrows():
        sd = r["matched_null_sd"]
        row = dict(modality="spatial_residual_moran",
                   family=f"spatial_two_program::{ds}",
                   unit_id=f"{ds}::{r['assay_id']}::{r['program_uid']}",
                   program_label=r["program_label"],
                   current_state=r["evidence_state"],
                   biological_unit=r["biological_unit"],
                   matched_null_sd=sd, family_size_m=m,
                   permutation_p_floor=floor)
        for name, a in [("bonferroni_family", 0.05 / m), ("nominal_005", 0.05),
                        ("bh_implied", max(a_bh, floor) if np.isfinite(floor) else a_bh)]:
            aa = a
            if np.isfinite(floor) and aa < floor:
                aa = floor  # permutation resolution caps how small alpha can be
            row[f"alpha_{name}"] = aa
            row[f"mde_{name}"] = ((norm.ppf(1 - aa) + Z80) * sd) if np.isfinite(sd) else np.nan
        srows.append(row)
sp_out = pd.DataFrame(srows)
sp_out["states_unchanged_note"] = ("spatial evidence calls and pass-1 proposed states are "
                                   "not revisited in this recompute; numbers only")
sp_out.to_csv(OUT / "mde_common_alpha_spatial.tsv", sep="\t", index=False)

pd.DataFrame(alpha_tbl).to_csv(OUT / "alpha_arms_by_family.tsv", sep="\t", index=False)

# ------------------------------------------------- reclassification matrix
mat = []
ind = da[da["evidence_state"] == "indeterminate"]
for lin, g in ind.groupby("lineage"):
    for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
        for aname, aval in ANCHORS_ALL.items():
            n = int((g[f"mde_{arm}_joint"] < aval).sum())
            mat.append(dict(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
                            alpha_arm=arm, anchor=aname, anchor_value=aval,
                            n_indeterminate=int(len(g)), n_tested_negative=n,
                            frac=n / len(g),
                            median_mde=float(np.nanmedian(g[f"mde_{arm}_joint"])),
                            p90_mde=float(np.nanpercentile(g[f"mde_{arm}_joint"], 90)),
                            arm_is_primary=(arm == "bonferroni_family"
                                            and aname in ANCHORS_PRIMARY)))
indp = pr[pr["cross_cohort_state"] == "indeterminate"]
for (coh, lin), g in indp.groupby(["cohort", "lineage"]):
    for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
        for aname, aval in ANCHORS_ALL.items():
            n = int((g[f"mde_{arm}"] < aval).sum())
            mat.append(dict(modality="atac_program_projection",
                            family=f"atac_program::{coh}::{lin}",
                            alpha_arm=arm, anchor=aname, anchor_value=aval,
                            n_indeterminate=int(len(g)), n_tested_negative=n,
                            frac=n / len(g),
                            median_mde=float(np.nanmedian(g[f"mde_{arm}"])),
                            p90_mde=float(np.nanpercentile(g[f"mde_{arm}"], 90)),
                            arm_is_primary=(arm == "bonferroni_family"
                                            and aname in ANCHORS_PRIMARY)))
mat = pd.DataFrame(mat)
mat.to_csv(OUT / "reclassification_matrix.tsv", sep="\t", index=False)

# ------------------------------------------------------------- power floor
# delta-free and (for the common arms) alpha-comparable.
fl = []
for lin, g in da.groupby("lineage"):
    for c in COH:
        for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
            fl.append(dict(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
                           cohort=c, alpha_arm=arm,
                           alpha=ARMS[(lin, c)][arm],
                           min_se=float(np.nanmin(g[f"standard_error_{c}"])),
                           median_se=float(np.nanmedian(g[f"standard_error_{c}"])),
                           power_floor_mde=float(np.nanmin(g[f"mde_{arm}_{c}"])),
                           median_mde=float(np.nanmedian(g[f"mde_{arm}_{c}"]))))
    for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
        fl.append(dict(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
                       cohort="joint", alpha_arm=arm, alpha=np.nan,
                       min_se=np.nan, median_se=np.nan,
                       power_floor_mde=float(np.nanmin(g[f"mde_{arm}_joint"])),
                       median_mde=float(np.nanmedian(g[f"mde_{arm}_joint"]))))
for (coh, lin), g in pr.groupby(["cohort", "lineage"]):
    if g["standard_error"].notna().sum() == 0:
        continue
    for arm in ["bonferroni_family", "nominal_005", "bh_implied"]:
        fl.append(dict(modality="atac_program_projection",
                       family=f"atac_program::{coh}::{lin}", cohort=coh, alpha_arm=arm,
                       alpha=float(g[f"alpha_{arm}"].iloc[0]) if arm != "nominal_005" else 0.05,
                       min_se=float(np.nanmin(g["standard_error"])),
                       median_se=float(np.nanmedian(g["standard_error"])),
                       power_floor_mde=float(np.nanmin(g[f"mde_{arm}"])),
                       median_mde=float(np.nanmedian(g[f"mde_{arm}"]))))
pd.DataFrame(fl).to_csv(OUT / "power_floor_common_alpha.tsv", sep="\t", index=False)

# --------------------------------------------------------- donors-needed
NPAIR = {("hepatocyte", "gse244832"): (5, 9), ("hepatocyte", "gse281367"): (6, 6),
         ("stellate", "gse244832"): (5, 9), ("stellate", "gse281367"): (6, 6)}
NGRID = [4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500]
curve = []
for lin, g in ind.groupby("lineage"):
    for arm in ["bonferroni_family", "nominal_005"]:
        for n in NGRID:
            sc = []
            for c in COH:
                n1, n2 = NPAIR[(lin, c)]
                fac = np.sqrt((2.0 / n) / (1.0 / n1 + 1.0 / n2))
                sc.append(mde(g[f"standard_error_{c}"] * fac, ARMS[(lin, c)][arm]))
            joint = np.nanmax(np.vstack(sc), axis=0)
            rec = dict(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
                       alpha_arm=arm, n_donors_per_arm_per_cohort=n,
                       n_indeterminate=int(len(g)),
                       median_mde=float(np.nanmedian(joint)))
            for aname, aval in ANCHORS_ALL.items():
                rec[f"frac_powered__{aname}"] = float(np.nanmean(joint < aval))
                rec[f"n_powered__{aname}"] = int(np.nansum(joint < aval))
            rec["arm_is_primary"] = (arm == "bonferroni_family")
            curve.append(rec)
pd.DataFrame(curve).to_csv(OUT / "donors_needed_curve_primary.tsv", sep="\t", index=False)

# donors needed per peak, per anchor, primary alpha arm
dn = []
for lin, g in ind.groupby("lineage"):
    for aname, aval in ANCHORS_ALL.items():
        need = []
        for c in COH:
            n1, n2 = NPAIR[(lin, c)]
            cur = mde(g[f"standard_error_{c}"], ARMS[(lin, c)]["bonferroni_family"])
            need.append(2.0 * (cur / aval) ** 2 / (1.0 / n1 + 1.0 / n2))
        need = np.nanmax(np.vstack(need), axis=0)
        dn.append(dict(modality="atac_accessibility", family=f"atac_da_peak::{lin}",
                       alpha_arm="bonferroni_family", anchor=aname, anchor_value=aval,
                       median_donors=float(np.nanmedian(need)),
                       q25=float(np.nanpercentile(need, 25)),
                       q75=float(np.nanpercentile(need, 75)),
                       p90=float(np.nanpercentile(need, 90)),
                       max_donors=float(np.nanmax(need)),
                       arm_is_primary=aname in ANCHORS_PRIMARY))
pd.DataFrame(dn).to_csv(OUT / "donors_needed_per_peak_primary.tsv", sep="\t", index=False)

# ------------------------------------------------- bulk reference calibration
bulk = pd.read_csv(BULK)
m_bulk = int(bulk["P.Value"].notna().sum())
a_bulk = {"bonferroni_family": 0.05 / m_bulk, "nominal_005": 0.05,
          "bh_implied": bh_implied_alpha(bulk["P.Value"])[0]}
brows = []
for name, a in a_bulk.items():
    mm = mde(bulk["SE"], a)
    brows.append(dict(modality="bulk_rnaseq_reference", family="bg001_frozen_canonical_deg",
                      alpha_arm=name, alpha=a, family_size_m=m_bulk,
                      min_se=float(np.nanmin(bulk["SE"])),
                      median_se=float(np.nanmedian(bulk["SE"])),
                      power_floor_mde=float(np.nanmin(mm)),
                      median_mde=float(np.nanmedian(mm)),
                      p90_mde=float(np.nanpercentile(mm, 90)),
                      frac_below_0_25=float(np.nanmean(mm < 0.25)),
                      frac_below_0_50=float(np.nanmean(mm < 0.50)),
                      frac_below_1_00=float(np.nanmean(mm < 1.00))))
bt = pd.DataFrame(brows)
pos = bulk[bulk["treat_fdr"] < 0.05]
bt["treat_positive_n"] = int(len(pos))
bt["treat_positive_abs_logFC_p10"] = float(np.percentile(pos["logFC"].abs(), 10)) if len(pos) else np.nan
bt["treat_positive_abs_logFC_median"] = float(np.median(pos["logFC"].abs())) if len(pos) else np.nan
bt["source_path"] = str(BULK)
bt["universe_note"] = (f"{len(bulk)} genes; this is the BG-001 fragment remediation frozen set, "
                       "NOT the promoted canonical release and NOT the 23,370-gene five-cohort "
                       "family; reference calibration only")
bt["modality_transfer_warning"] = ("bulk RNA log2FC and snATAC accessibility log2FC are "
                                   "different measurement scales; any anchor transported "
                                   "between them is an assumption, not an equivalence")
bt.to_csv(OUT / "bulk_reference_calibration.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ headline
prim = mat[(mat["alpha_arm"] == "bonferroni_family") &
           (mat["modality"] == "atac_accessibility")]
head = dict(
    seed=SEED, primary_alpha_arm="bonferroni_family",
    primary_anchors=list(ANCHORS_PRIMARY.values()),
    sensitivity_arms=["nominal_005", "bh_implied", "delta_star_0.797"],
    atac_peaks_indeterminate=int(len(ind)),
    primary_tested_negative_at_0_25=int(prim[prim["anchor"] == "log2fc_0.25"]["n_tested_negative"].sum()),
    primary_tested_negative_at_0_50=int(prim[prim["anchor"] == "log2fc_0.50"]["n_tested_negative"].sum()),
    primary_tested_negative_at_1_00=int(prim[prim["anchor"] == "log2fc_1.00"]["n_tested_negative"].sum()),
    sensitivity_tested_negative_at_delta_star=int(
        prim[prim["anchor"] == "delta_star_0.797_SENSITIVITY"]["n_tested_negative"].sum()),
    pass1_bh_implied_delta_star_count=6590,
    program_indeterminate=int(len(indp)),
    program_primary_tested_negative_at_1_00=int(
        mat[(mat["modality"] == "atac_program_projection") &
            (mat["alpha_arm"] == "bonferroni_family") &
            (mat["anchor"] == "log2fc_1.00")]["n_tested_negative"].sum()),
    bulk_reference_path=str(BULK), bulk_genes=int(len(bulk)),
)
json.dump(head, open(OUT / "headline_common_alpha.json", "w"), indent=2, default=float)

inputs = [ATAC / "da/da_peak_results.tsv.gz", ATAC / "programs/program_atac_results.tsv",
          SPAT / "effects/spatial_program_effects.tsv", BULK]
man = ([dict(role="input", path=str(p), sha256=sha256(p), bytes=p.stat().st_size) for p in inputs] +
       [dict(role="output", path=str(p), sha256=sha256(p), bytes=p.stat().st_size)
        for p in sorted(OUT.glob("*")) if p.is_file()])
pd.DataFrame(man).to_csv(OUT / "checksum_manifest.tsv", sep="\t", index=False)

print(json.dumps(head, indent=2, default=float))
print("\n== alpha arms ==")
print(pd.DataFrame(alpha_tbl).to_string(index=False))
print("\n== reclassification matrix (ATAC peaks) ==")
print(mat[mat["modality"] == "atac_accessibility"].to_string(index=False))
print("\n== reclassification matrix (program) ==")
print(mat[mat["modality"] == "atac_program_projection"].to_string(index=False))
print("\n== power floor ==")
print(pd.DataFrame(fl).to_string(index=False))
print("\n== donors needed per peak ==")
print(pd.DataFrame(dn).to_string(index=False))
print("\n== donors curve (bonferroni arm) ==")
cc = pd.DataFrame(curve)
print(cc[cc["alpha_arm"] == "bonferroni_family"][
    ["family", "n_donors_per_arm_per_cohort", "median_mde",
     "frac_powered__log2fc_0.25", "frac_powered__log2fc_0.50",
     "frac_powered__log2fc_1.00"]].to_string(index=False))
print("\n== bulk reference ==")
print(bt.to_string(index=False))
print("\n== spatial (numbers only) ==")
print(sp_out.to_string(index=False))
