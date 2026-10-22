#!/usr/bin/env python3
"""Independent recomputation of V4 reported numbers from the saved output tables.

Reads per_array_results.tsv and the per-array null matrices, recomputes the
donor summaries with separate code (pandas groupby, t-based interval, a second
sign-flip seed), and prints them next to summary_donor_tests.tsv.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

out = Path(sys.argv[1])
arr = pd.read_csv(out / "per_array_results.tsv", sep="\t", dtype={"individual": str})
summ = pd.read_csv(out / "summary_donor_tests.tsv", sep="\t")
rows = []
for (ds, mode), part in arr.groupby(["dataset", "mode"]):
    null = np.load(out / f"null_per_array__{ds}__{mode}.npy")
    part = part.sort_values("array_id").reset_index(drop=True)
    # recompute centered d per array from the raw null matrix
    d_arr = part["residual_I"].to_numpy() - null.mean(axis=0)
    assert np.allclose(d_arr, part["centered_d"].to_numpy(), atol=1e-12), (ds, mode)
    part["d"] = d_arr
    sets = {"all": part}
    if ds == "HRA007511_HMSMA":
        prim = part[part["condition"] == "primary"]
        top4 = prim.sort_values("d", ascending=False).head(4)["individual"]
        sets = {"primary_32": prim, "all_35": part, "primary_drop_top4": prim[~prim["individual"].isin(top4)]}
    for label, p in sets.items():
        donor = p.groupby("individual")["d"].mean().to_numpy()
        n = len(donor)
        rng = np.random.default_rng(12345)
        flips = rng.choice([-1.0, 1.0], size=(200_000, n)) @ donor / n
        t_lo, t_hi = stats.t.interval(0.95, n - 1, loc=donor.mean(), scale=stats.sem(donor)) if n > 1 else (np.nan, np.nan)
        cols = [part.index[part["array_id"] == a][0] for a in p["array_id"]]
        cohort_null = null[:, cols].mean(axis=1) if ds != "GSE192741" else None
        rows.append({"dataset": ds, "mode": mode, "analysis_set": label, "n_donors": n,
                     "mean_d": donor.mean(), "n_pos": int((donor > 0).sum()),
                     "t95_low": t_lo, "t95_high": t_hi,
                     "signflip_p_seed2": (1 + (flips >= donor.mean()).sum()) / (200_001),
                     "wilcoxon_p_one_sided": stats.wilcoxon(donor, alternative="greater").pvalue if n > 5 else np.nan,
                     "cohort_centered": (p.groupby("individual")["residual_I"].mean().mean() - cohort_null.mean())
                     if cohort_null is not None else np.nan})
re = pd.DataFrame(rows)
pd.set_option("display.width", 250, "display.max_columns", 30)
print(re.to_string(index=False))
print()
print(summ[["dataset", "mode", "analysis_set", "n_donors", "mean_centered_d", "n_donors_d_positive", "boot95_low",
            "boot95_high", "signflip_p_one_sided", "holm_p", "cohort_centered", "cohort_null_sd",
            "cohort_matched_null_p", "n_donors_above_own_null_q95"]].to_string(index=False))
re.to_csv(out / "verification_recompute.tsv", sep="\t", index=False)
