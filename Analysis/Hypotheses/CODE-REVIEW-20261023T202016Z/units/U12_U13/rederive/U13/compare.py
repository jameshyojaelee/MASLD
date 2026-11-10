#!/usr/bin/env python3
"""Compare U13 re-derivation runs with the recorded producer tables and across seeds/orders."""
from pathlib import Path

import numpy as np
import pandas as pd

U = Path(__file__).resolve().parent
OUT = U / "out"
P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
R2 = P / ("Analysis/Multimodal_Program_Projection/candidates/"
          "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/full_composition_program_results.tsv")
J1731 = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1731/out/program_level.tsv"
AN4 = P / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN4/20261021T124414Z"
NAME = {"hotspot_hepatocytes_f05c535ae5bbc0b9": "Stromal ECM", "hotspot_hepatocytes_48f39dd4d817a10e": "Ductular injury"}

runs = sorted(d for d in OUT.iterdir() if d.is_dir() and d.name != "smoke_gse")
focal = pd.concat([pd.read_csv(d / "focal_results.tsv", sep="\t").assign(run=d.name)
                   for d in runs if (d / "focal_results.tsv").exists()], ignore_index=True)
focal["program"] = focal["program_id"].map(NAME)
cols = ["run", "dataset", "seed", "order", "scheme", "design", "program", "n_measured", "retained_l1", "observed",
        "null_mean", "excess", "k", "p", "q_bh_two_program", "n_units_excess_pos", "signflip_p_one_sided", "sha_match"]
focal[cols].to_csv(U / "focal_all_runs.tsv", sep="\t", index=False)

# recorded values
r2 = pd.read_csv(R2, sep="\t")
rec = []
for r in r2.itertuples():
    rec.append({"dataset": r.dataset, "scheme": "hepatocytes", "design": "all16", "program": NAME[r.program_id],
                "rec_observed": r.residual_moran_i, "rec_null_mean": r.residual_null_mean, "rec_p": r.residual_pvalue,
                "rec_q": r.residual_padj})
j = pd.read_csv(J1731, sep="\t")
for r in j[j.design == "lineage_only"].itertuples():
    rec.append({"dataset": r.dataset, "scheme": "hepatocytes", "design": "hep_only", "program": NAME[r.program_id],
                "rec_observed": r.residual_moran_i, "rec_null_mean": r.null_mean, "rec_p": r.p,
                "rec_q": r.q_bh_two_program_within_dataset})
fl = pd.read_csv(AN4 / "focal_lineage/focal_lineage_results.tsv", sep="\t")
for r in fl[fl.scheme == "content_lineage"].itertuples():
    rec.append({"dataset": r.dataset, "scheme": "content", "design": "all16", "program": NAME[r.program_id],
                "rec_observed": r.residual_moran_i, "rec_null_mean": r.residual_null_mean, "rec_p": r.residual_pvalue,
                "rec_q": r.residual_padj_two_program})
rec = pd.DataFrame(rec)
cmp = focal.merge(rec, on=["dataset", "scheme", "design", "program"], how="left")
cmp["d_observed"] = cmp["observed"] - cmp["rec_observed"]
cmp["d_null_mean"] = cmp["null_mean"] - cmp["rec_null_mean"]
cmp[cols + ["rec_observed", "rec_null_mean", "rec_p", "rec_q", "d_observed", "d_null_mean"]].to_csv(
    U / "focal_vs_recorded.tsv", sep="\t", index=False)
with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(cmp[["run", "scheme", "design", "program", "observed", "rec_observed", "excess", "p", "rec_p",
               "q_bh_two_program", "rec_q", "sha_match", "n_units_excess_pos", "signflip_p_one_sided"]]
          .round(5).to_string(index=False))

# family
fam = []
for d in runs:
    f = d / "family_results.tsv"
    if f.exists():
        fam.append(pd.read_csv(f, sep="\t").assign(run=d.name))
if fam:
    fam = pd.concat(fam, ignore_index=True)
    fam = fam[fam.testable.astype(bool)]
    an4 = pd.read_csv(AN4 / "combined/family_all_programs.tsv", sep="\t")
    an4 = an4[an4.testable.astype(bool)][["dataset", "program_id", "residual_moran_i", "residual_null_mean",
                                          "residual_pvalue", "residual_padj_family", "rank_excess"]]
    m = fam.merge(an4, on=["dataset", "program_id"], how="left")
    m["call"] = m["q_bh_family"] < 0.05
    m["an4_call"] = m["residual_padj_family"] < 0.05
    rows = []
    for (run, ds), g in m.groupby(["run", "dataset"]):
        rows.append({"run": run, "dataset": ds, "n_testable": len(g), "n_q_lt_0.05": int(g.call.sum()),
                     "an4_n_q_lt_0.05": int(g.an4_call.sum()),
                     "calls_gained_vs_an4": int((g.call & ~g.an4_call).sum()),
                     "calls_lost_vs_an4": int((~g.call & g.an4_call).sum()),
                     "n_sha_match": int(g.sha_match.sum()),
                     "max_abs_d_observed": float(np.max(np.abs(g.observed - g.residual_moran_i))),
                     "n_p_equal_an4": int(np.isclose(g.p, g.residual_pvalue, atol=1e-12).sum()),
                     "stromal_rank": int(g.loc[g.program_id == "hotspot_hepatocytes_f05c535ae5bbc0b9", "rank_excess_x"].iloc[0]),
                     "ductular_rank": int(g.loc[g.program_id == "hotspot_hepatocytes_48f39dd4d817a10e", "rank_excess_x"].iloc[0]),
                     "stromal_q": float(g.loc[g.program_id == "hotspot_hepatocytes_f05c535ae5bbc0b9", "q_bh_family"].iloc[0]),
                     "ductular_q": float(g.loc[g.program_id == "hotspot_hepatocytes_48f39dd4d817a10e", "q_bh_family"].iloc[0])})
    summ = pd.DataFrame(rows)
    summ.to_csv(U / "family_runs_summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))
    flips = m[m.call != m.an4_call][["run", "dataset", "program_id", "program_name", "p", "residual_pvalue",
                                     "q_bh_family", "residual_padj_family"]]
    flips.to_csv(U / "family_call_flips_vs_an4.tsv", sep="\t", index=False)
    print(flips.to_string(index=False))
