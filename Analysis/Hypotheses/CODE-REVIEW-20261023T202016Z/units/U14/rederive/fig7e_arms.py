"""U14 re-derivation of the Figure 7E values from saved per-donor outputs.

Scale (legend): (estimate - null) / (threshold - null), 0 = null, 1 = threshold.
Severity arms: estimate = Delta = mean home advantage in LOW donors + mean home advantage in HIGH donors,
home advantage h = skill of the own-stratum map minus skill of the other-stratum map (PRECOMMIT.md);
null = mean of the label-permutation null; threshold - null = MDE = 2.802 * sd(null) (AMENDMENT_02).
Interval: stratified bootstrap (donors resampled within their stratum), 10,000 draws, percentile 95%.
Per-donor values are saved only for the four arms in results_GSE267145.json / results_GSE296875.json;
for the three post hoc sensitivity arms only the producer's interval ends are saved, so only the scaling is checked.

usage: fig7e_arms.py OUTDIR SEED [SEED ...]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BENCH = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
SEV = BENCH / "crossmodal-severity-stratification-20260831T190048Z/results"
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
seeds = [int(s) for s in sys.argv[2:]]
B = 10000

summary = {(a["cohort"], a["arm"]): a for a in json.load(open(SEV / "ALL_ARMS_SUMMARY.json"))}
rows = []
for f, cohort in [("results_GSE267145.json", "GSE267145"), ("results_GSE296875.json", "GSE296875")]:
    res = json.load(open(SEV / f))
    for arm, st in res["strata"].items():
        a2 = st["arm2"]
        pdn = a2["per_donor"]
        isl = np.asarray(pdn["is_low"], bool)
        sL = np.asarray(pdn["skill_low_trained"], float)
        sH = np.asarray(pdn["skill_high_trained"], float)
        h = np.where(isl, sL - sH, sH - sL)
        assert np.allclose(h, np.asarray(pdn["home_advantage"], float))
        hl, hh = h[isl], h[~isl]
        delta = hl.mean() + hh.mean()
        null = np.asarray(a2["label_permutation_null_draws"], float)
        null_mean, null_sd = null.mean(), null.std(ddof=1)
        s = summary[(cohort, arm)]
        mde = 2.802 * null_sd
        row = {"cohort": cohort, "arm": arm, "precommitted": s["precommitted"], "n_low": int(isl.sum()), "n_high": int((~isl).sum()),
               "delta": delta, "delta_saved": s["delta"], "null_mean": null_mean, "null_mean_saved": s["null_mean"],
               "mde_2.802sd": mde, "mde_saved": s["mde"], "null_sd_ddof0": null.std(ddof=0),
               "perm_p_two_sided_recomputed": float(np.mean(np.abs(null - null_mean) >= abs(delta - null_mean))),
               "p_saved": s["p"], "verdict_saved": s["verdict"],
               "scaled_estimate": (delta - null_mean) / s["mde"]}
        saved_ci = a2["secondary_stratified_bootstrap_precision_only"]
        row["saved_boot_lo"], row["saved_boot_hi"] = saved_ci["boot_ci95_low"], saved_ci["boot_ci95_high"]
        row["scaled_saved_lo"] = (saved_ci["boot_ci95_low"] - s["null_mean"]) / s["mde"]
        row["scaled_saved_hi"] = (saved_ci["boot_ci95_high"] - s["null_mean"]) / s["mde"]
        for seed in seeds:
            rng = np.random.default_rng(seed)
            bl = hl[rng.integers(0, len(hl), (B, len(hl)))].mean(1)
            bh = hh[rng.integers(0, len(hh), (B, len(hh)))].mean(1)
            d = bl + bh
            lo, hi = np.percentile(d, [2.5, 97.5])
            row[f"boot_lo_seed{seed}"], row[f"boot_hi_seed{seed}"] = lo, hi
            row[f"scaled_lo_seed{seed}"] = (lo - s["null_mean"]) / s["mde"]
            row[f"scaled_hi_seed{seed}"] = (hi - s["null_mean"]) / s["mde"]
        rows.append(row)

# post hoc sensitivity arms: scaling of the saved values only
for f, keyarm in [("sensitivity_no_nor_GSE267145.json", "posthoc_NOR_dropped_activity_ge4")]:
    d = json.load(open(SEV / f))
    s = summary[("GSE267145", keyarm)]
    ci = d["stratified_bootstrap_precision_only"]
    null = np.asarray(d["label_permutation_null_draws"], float)
    rows.append({"cohort": "GSE267145", "arm": keyarm, "precommitted": s["precommitted"], "n_low": s["n_low"], "n_high": s["n_high"],
                 "delta": d["interaction_did"]["delta"], "delta_saved": s["delta"], "null_mean": null.mean(),
                 "null_mean_saved": s["null_mean"], "mde_2.802sd": 2.802 * null.std(ddof=1), "mde_saved": s["mde"],
                 "perm_p_two_sided_recomputed": float(np.mean(np.abs(null - null.mean()) >= abs(s["delta"] - null.mean()))),
                 "p_saved": s["p"], "verdict_saved": s["verdict"],
                 "scaled_estimate": (s["delta"] - s["null_mean"]) / s["mde"],
                 "saved_boot_lo": ci["boot_ci95_low"], "saved_boot_hi": ci["boot_ci95_high"],
                 "scaled_saved_lo": (ci["boot_ci95_low"] - s["null_mean"]) / s["mde"],
                 "scaled_saved_hi": (ci["boot_ci95_high"] - s["null_mean"]) / s["mde"]})
d = json.load(open(SEV / "sensitivity_sex_GSE267145.json"))
for arm_json, keyarm in zip(d["arms"], ["posthoc_A_female_only_activity_ge3", "posthoc_B_female_only_NOR_dropped_activity_ge4"]):
    s = summary[("GSE267145", keyarm)]
    ci = arm_json["stratified_bootstrap_precision_only"]
    null = np.asarray(arm_json["label_permutation_null_draws"], float)
    rows.append({"cohort": "GSE267145", "arm": keyarm, "precommitted": s["precommitted"], "n_low": s["n_low"], "n_high": s["n_high"],
                 "delta": arm_json["interaction_did"]["delta"], "delta_saved": s["delta"], "null_mean": null.mean(),
                 "null_mean_saved": s["null_mean"], "mde_2.802sd": 2.802 * null.std(ddof=1), "mde_saved": s["mde"],
                 "perm_p_two_sided_recomputed": float(np.mean(np.abs(null - null.mean()) >= abs(s["delta"] - null.mean()))),
                 "p_saved": s["p"], "verdict_saved": s["verdict"],
                 "scaled_estimate": (s["delta"] - s["null_mean"]) / s["mde"],
                 "saved_boot_lo": ci["boot_ci95_low"], "saved_boot_hi": ci["boot_ci95_high"],
                 "scaled_saved_lo": (ci["boot_ci95_low"] - s["null_mean"]) / s["mde"],
                 "scaled_saved_hi": (ci["boot_ci95_high"] - s["null_mean"]) / s["mde"]})

t = pd.DataFrame(rows)
t.to_csv(out / "fig7e_arms_rederived.tsv", sep="\t", index=False, float_format="%.6f")
with pd.option_context("display.width", 250, "display.max_columns", 40):
    print(t.T.to_string())
