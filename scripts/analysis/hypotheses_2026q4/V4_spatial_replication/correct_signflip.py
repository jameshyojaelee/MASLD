#!/usr/bin/env python3
"""Recompute the donor sign-flip P exactly and re-apply Holm.

Fixes a floating-point comparison in v4_spatial_replication.py (job 22093979):
the all-plus flip mean (BLAS dot / n) could round 1 ulp below np.mean(d), so
the identity pattern was not counted and the Monte Carlo P fell below its
exact floor 2^-n. Exact P: full enumeration for n <= 20; for larger n, if
every donor d > 0, only the identity pattern reaches the observed sum, so
P = 2^-n exactly. Writes a new file; the original summary is not modified.
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

out = Path(sys.argv[1])
donor = pd.read_csv(out / "per_donor_results.tsv", sep="\t", dtype={"individual": str})
summ = pd.read_csv(out / "summary_donor_tests.tsv", sep="\t")
arr = pd.read_csv(out / "per_array_results.tsv", sep="\t", dtype={"individual": str})


def exact_p(d):
    d = np.asarray(d, float); n = len(d); obs = d.sum()
    if n <= 20:
        signs = np.array(list(itertools.product([-1.0, 1.0], repeat=n)))
        return float(np.mean(signs @ d >= obs - 1e-12)), "enumeration"
    if np.all(d > 0):
        return float(2.0 ** -n), "all_positive_identity_only"
    raise RuntimeError("n > 20 with a non-positive donor: exact P needs another method")


rows = []
for _, r in summ.iterrows():
    sub = donor[(donor.dataset == r.dataset) & (donor["mode"] == r["mode"]) & (donor.analysis_set == r.analysis_set)]
    if r.analysis_set == "primary_drop_top4":
        prim = arr[(arr.dataset == r.dataset) & (arr.condition == "primary")]
        keep = set(prim["individual"]) - set(str(r.dropped).split(","))
        d = prim[prim["individual"].isin(keep)].groupby("individual")["centered_d"].mean().to_numpy()
    else:
        d = sub["donor_centered_d"].to_numpy()
    p, method = exact_p(d)
    rows.append({"dataset": r.dataset, "mode": r["mode"], "analysis_set": r.analysis_set, "n_donors": len(d),
                 "mean_centered_d": d.mean(), "signflip_p_mc_original": r.signflip_p_one_sided,
                 "signflip_p_exact": p, "exact_method": method})
res = pd.DataFrame(rows)
fam = ((res.dataset == "HRA007511_HMSMA") & (res.analysis_set == "primary_32")) | \
      ((res.dataset == "Vu_et_al_2025") & (res["mode"] == "c2l"))
p = res.loc[fam, "signflip_p_exact"].to_numpy(); m = len(p); order = np.argsort(p)
adj = np.empty(m); run = 0.0
for k, i in enumerate(order):
    run = max(run, min(1.0, (m - k) * p[i])); adj[i] = run
res.loc[fam, "holm_p_exact"] = adj
res.to_csv(out / "summary_signflip_exact.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(res.to_string(index=False))
