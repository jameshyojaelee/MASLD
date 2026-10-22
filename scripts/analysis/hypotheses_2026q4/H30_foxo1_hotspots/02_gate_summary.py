#!/usr/bin/env python3
"""H30 detectability gate summary: join the blind gate table (depth and background only) with the
774-unit roster and stage, count persons that pass, and apply the prespecified stop rule.

Inputs : GATE_TSV from reads_extraction/h30_power_gate.py (no alt counts), ROSTER_STAGE_TSV
Outputs: OUT_DIR/gate_person_site.tsv, gate_by_stage.tsv, gate_by_cohort.tsv, gate_decision.json
Usage  : 02_gate_summary.py GATE_TSV ROSTER_STAGE_TSV OUT_DIR
Thresholds are those in PRESPEC_H30.md: person passes if min MDV over R21L and S22W <= 0.10;
stop unless >= 30 pass overall, >= 10 in F0 and >= 10 in F3-F4.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

gate_path, roster_path, out_dir = sys.argv[1:4]
os.makedirs(out_dir, exist_ok=False)
THRESH = [0.20, 0.10, 0.05]
PRIMARY = 0.10
MIN_ALL, MIN_LOW, MIN_HIGH = 30, 10, 10

gate = pd.read_csv(gate_path, sep="\t")
assert not any(c.startswith("alt") for c in gate.columns), "gate table must not carry alt counts"
roster = pd.read_csv(roster_path, sep="\t")
assert len(roster) == 774 and roster.sample_id.is_unique

g = gate.merge(roster[["sample_id", "dataset", "fibrosis_stage"]], left_on="library", right_on="sample_id", how="inner")
g = g.drop(columns="sample_id")
g.to_csv(os.path.join(out_dir, "gate_person_site.tsv"), sep="\t", index=False)
missing = sorted(set(roster.sample_id) - set(gate.library))

wide = g.pivot(index="library", columns="hotspot", values=["depth_distinct_fragments", "min_detectable_fraction"])
wide.columns = [f"{a}_{b}" for a, b in wide.columns]
wide = wide.reset_index().merge(roster, left_on="library", right_on="sample_id", how="right")
mdv_cols = [c for c in wide.columns if c.startswith("min_detectable_fraction_")]
dep_cols = [c for c in wide.columns if c.startswith("depth_distinct_fragments_")]
for c in mdv_cols:
    wide[c] = pd.to_numeric(wide[c], errors="coerce")
wide["min_mdv"] = wide[mdv_cols].min(axis=1, skipna=True)
wide["max_depth"] = wide[dep_cols].max(axis=1, skipna=True)


def stage_group(s):
    if pd.isna(s):
        return "missing"
    return "F0" if s == 0 else ("F3-F4" if s >= 3 else "F1-F2")


wide["stage_group"] = wide.fibrosis_stage.map(stage_group)
wide["stage_label"] = wide.fibrosis_stage.map(lambda s: "missing" if pd.isna(s) else f"F{int(s)}")


def table(by):
    rows = []
    for key, d in wide.groupby(by, dropna=False):
        r = {by: key, "persons": len(d),
             "with_library_output": int(d.library.notna().sum()),
             "depth_gt0_any_site": int((d.max_depth > 0).sum()),
             "median_max_depth": float(np.nanmedian(d.max_depth)) if d.max_depth.notna().any() else np.nan,
             "p90_max_depth": float(np.nanpercentile(d.max_depth.dropna(), 90)) if d.max_depth.notna().any() else np.nan,
             "max_max_depth": float(np.nanmax(d.max_depth)) if d.max_depth.notna().any() else np.nan,
             "median_min_mdv": float(np.nanmedian(d.min_mdv)) if d.min_mdv.notna().any() else np.nan}
        for t in THRESH:
            r[f"pass_mdv_le_{t:g}"] = int((d.min_mdv <= t).sum())
        rows.append(r)
    return pd.DataFrame(rows)


by_stage = table("stage_label")
by_group = table("stage_group")
by_cohort = table("dataset")
by_stage.to_csv(os.path.join(out_dir, "gate_by_stage.tsv"), sep="\t", index=False)
by_group.to_csv(os.path.join(out_dir, "gate_by_stage_group.tsv"), sep="\t", index=False)
by_cohort.to_csv(os.path.join(out_dir, "gate_by_cohort.tsv"), sep="\t", index=False)
wide.to_csv(os.path.join(out_dir, "gate_person_wide.tsv"), sep="\t", index=False)

n_all = int((wide.min_mdv <= PRIMARY).sum())
n_low = int(((wide.min_mdv <= PRIMARY) & (wide.stage_group == "F0")).sum())
n_high = int(((wide.min_mdv <= PRIMARY) & (wide.stage_group == "F3-F4")).sum())
n_f4 = int(((wide.min_mdv <= PRIMARY) & (wide.stage_label == "F4")).sum())
decision = {
    "primary_mdv_threshold": PRIMARY, "persons_in_roster": int(len(roster)),
    "roster_libraries_without_gate_row": missing,
    "pass_overall": n_all, "pass_F0": n_low, "pass_F3_F4": n_high, "pass_F4_only": n_f4,
    "rule": {"overall_min": MIN_ALL, "F0_min": MIN_LOW, "F3_F4_min": MIN_HIGH},
    "gate_passes": bool(n_all >= MIN_ALL and n_low >= MIN_LOW and n_high >= MIN_HIGH),
    "persons_with_depth_gt0_any_site": int((wide.max_depth > 0).sum()),
    "pass_at_other_thresholds": {str(t): int((wide.min_mdv <= t).sum()) for t in THRESH},
}
with open(os.path.join(out_dir, "gate_decision.json"), "w") as fh:
    json.dump(decision, fh, indent=1)
print(json.dumps(decision, indent=1))
print(by_stage.to_string(index=False))
print(by_group.to_string(index=False))
print(by_cohort.to_string(index=False))
