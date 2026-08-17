#!/usr/bin/env python
"""Unified per-row MDE table + sensitivity of the reclassification to delta*.

The prespecified meaningful effect delta* is derived from the `supported` calls
in the same modality. Those calls were themselves selected by passing the same
underpowered test, so delta* is upward-biased by selection. This script shows
how the tested_negative count moves as delta* is varied.
"""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

np.random.seed(20260813)
OUT = Path(os.environ["MDE_OUT"])
TAB = OUT / "tables"

da = pd.read_csv(TAB / "mde_atac_da_peaks.tsv.gz", sep="\t")
pr = pd.read_csv(TAB / "mde_atac_program_projection.tsv", sep="\t")
sp = pd.read_csv(TAB / "mde_spatial_two_program.tsv", sep="\t")
cov = pd.read_csv(TAB / "mde_spatial_coverage_family.tsv.gz", sep="\t")

COLS = ["modality", "family", "unit_id", "current_state", "uncertainty_scale",
        "se_or_null_sd", "alpha_eff", "mde", "prespecified_meaningful_effect",
        "meaningful_effect_rule", "proposed_state", "reason"]


def blk(df, cur, scale, sd, alpha, mde, prop, reason):
    o = pd.DataFrame({
        "modality": df["modality"], "family": df["family"], "unit_id": df["unit_id"],
        "current_state": df[cur], "uncertainty_scale": scale,
        "se_or_null_sd": df[sd] if sd else np.nan,
        "alpha_eff": df[alpha] if alpha else np.nan,
        "mde": df[mde] if mde else np.nan,
        "prespecified_meaningful_effect": df["prespecified_meaningful_effect"]
        if "prespecified_meaningful_effect" in df else np.nan,
        "meaningful_effect_rule": df["meaningful_effect_rule"]
        if "meaningful_effect_rule" in df else "",
        "proposed_state": df[prop], "reason": df[reason]})
    return o[COLS]


da["se_joint_worst"] = da[["standard_error_gse244832", "standard_error_gse281367"]].max(axis=1)
da["alpha_eff_worst"] = da[["alpha_eff_gse244832", "alpha_eff_gse281367"]].min(axis=1)
uni = pd.concat([
    blk(da, "evidence_state", "log2_accessibility_wald_se_two_cohort_joint",
        "se_joint_worst", "alpha_eff_worst", "mde_joint_bh",
        "proposed_state_bh", "proposed_reason_bh"),
    blk(pr, "current_state", "program_score_wald_se",
        "standard_error", "alpha_eff_bh", "mde_bh",
        "proposed_state_bh", "proposed_reason_bh"),
    blk(sp, "current_state", "centered_residual_moran_i_matched_gene_null_sd",
        "matched_null_sd", "alpha_eff_bh", "mde_bh",
        "proposed_state_bh", "proposed_reason_bh"),
    blk(cov, "coverage_status", "none_observability_only",
        "standard_error_or_null_sd", "alpha_eff", "mde",
        "proposed_state", "proposed_reason"),
], ignore_index=True)
uni.to_csv(TAB / "mde_all_rows_unified.tsv.gz", sep="\t", index=False)

# --- delta* sensitivity -----------------------------------------------------
sens = []
ind = da[da["evidence_state"] == "indeterminate"]
d_obs = float(da["prespecified_meaningful_effect"].iloc[0])
for mult, lab in [(0.25, "quarter"), (0.5, "half"), (0.75, "three_quarter"),
                  (1.0, "as_derived"), (1.5, "one_and_half"), (2.0, "double")]:
    d = d_obs * mult
    for lin, g in ind.groupby("family"):
        sens.append(dict(modality="atac_accessibility", family=lin, scenario=lab,
                         delta=d, n_indeterminate=int(len(g)),
                         n_tested_negative=int((g["mde_joint_bh"] < d).sum()),
                         frac=float((g["mde_joint_bh"] < d).mean())))
# absolute anchors on the log2 scale
for d, lab in [(0.25, "log2fc_0.25_limma_treat_style"), (0.5, "log2fc_0.50"),
               (1.0, "log2fc_1.00_twofold")]:
    for lin, g in ind.groupby("family"):
        sens.append(dict(modality="atac_accessibility", family=lin, scenario=lab,
                         delta=d, n_indeterminate=int(len(g)),
                         n_tested_negative=int((g["mde_joint_bh"] < d).sum()),
                         frac=float((g["mde_joint_bh"] < d).mean())))
sens = pd.DataFrame(sens)
sens.to_csv(TAB / "delta_sensitivity_atac_da.tsv", sep="\t", index=False)

head = json.load(open(TAB / "headline_counts.json"))
head["unified_rows"] = int(len(uni))
head["unified_tested_negative"] = int((uni["proposed_state"] == "tested_negative").sum())
json.dump(head, open(TAB / "headline_counts.json", "w"), indent=2, default=float)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


man = pd.read_csv(OUT / "checksum_manifest.tsv", sep="\t")
man = man[man["role"] == "input"]
rows = man.to_dict("records") + [
    dict(role="output", path=str(p), sha256=sha256(p), bytes=p.stat().st_size)
    for p in sorted(TAB.glob("*"))]
pd.DataFrame(rows).to_csv(OUT / "checksum_manifest.tsv", sep="\t", index=False)

print(sens.to_string(index=False))
print("\nunified rows:", len(uni), " tested_negative:", head["unified_tested_negative"])
print(uni.groupby(["modality", "current_state", "proposed_state"]).size().to_string())
