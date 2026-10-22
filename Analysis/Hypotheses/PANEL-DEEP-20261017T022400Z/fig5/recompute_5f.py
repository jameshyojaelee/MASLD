"""Independent recomputation of Figure 5F calibration values from the spatial analysis outputs.

Upstream: Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1731/out/{program_level.tsv, unit_level.tsv, null_donor_draws.npz}
  excess = residual Moran's I - mean of 9,999 matched-gene null draws (collapsed unit)
  P = (k+1)/(n+1), k = draws >= observed; BH over the two programs within each source.
  Matched-null 95% range of the excess = observed - quantiles of the null draws.
  Per-unit excess (GSE192741 donors, Vu arrays) from per-unit null draws; exact one-sided sign-flip P.
Cross-check against F5ef/repair/fig5f_spatial_program_calibration_corrected_r2.tsv afterwards.
"""
import itertools
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
J = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1731/out"
OUT = P / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute"
pl = pd.read_csv(J / "program_level.tsv", sep="\t")
ul = pd.read_csv(J / "unit_level.tsv", sep="\t")
nz = np.load(J / "null_donor_draws.npz", allow_pickle=True)
rows = []
for (ds, prog, des), g in pl.groupby(["dataset", "program_name", "design"]):
    r = g.iloc[0]
    draws = nz[f"{ds}|{prog}|{des}|collapsed"]
    k = int((draws >= r.residual_moran_i).sum())
    pv = (k + 1) / (len(draws) + 1)
    ex = r.residual_moran_i - draws.mean()
    lo, hi = np.quantile(draws - draws.mean(), [0.025, 0.975])
    rows.append({"dataset": ds, "program": prog, "design": des, "moran_i": r.residual_moran_i, "null_mean": draws.mean(), "null_sd": draws.std(ddof=1),
                 "excess": ex, "k": k, "n_draws": len(draws), "p": pv, "null95_lo": lo, "null95_hi": hi,
                 "excess_z": ex / draws.std(ddof=1)})
d = pd.DataFrame(rows)
d["q"] = d.groupby(["dataset", "design"]).p.transform(lambda s: stats.false_discovery_control(s, method="bh"))
print(d.round(5).to_string())
# compare with recorded program_level values
cmp = d.merge(pl, left_on=["dataset", "program", "design"], right_on=["dataset", "program_name", "design"])
print("max |excess diff|", (cmp.excess_x - cmp.excess_y).abs().max(), "k equal", (cmp.k == cmp.k_null_ge_obs).all(), "max |q diff|", (cmp.q - cmp.q_bh_two_program_within_dataset).abs().max())
d.to_csv(OUT / "5f_recomputed_program_level.tsv", sep="\t", index=False)

# per-unit (all16)
u = ul[ul.design == "all16"].copy()
print("unit types", u.groupby(["dataset", "unit_type"]).size().to_dict())
urows = []
for (ds, prog, ut), g in u.groupby(["dataset", "program_name", "unit_type"]):
    for _, r in g.iterrows():
        key = f"{ds}|{prog}|all16|{r.unit_id if ut != 'donor' else r.individual}"
        if key in nz.files:
            dr = nz[key]
            ex = r.residual_moran_i - dr.mean()
        else:
            ex = r.excess
        urows.append({"dataset": ds, "program": prog, "unit_type": ut, "unit_id": r.unit_id, "individual": r.individual, "excess": ex, "excess_recorded": r.excess, "moran_i": r.residual_moran_i})
U = pd.DataFrame(urows)
print(U.round(4).to_string())
U.to_csv(OUT / "5f_recomputed_unit_level.tsv", sep="\t", index=False)

def signflip_p(x):
    x = np.asarray(x); obs = x.mean(); n = len(x)
    cnt = 0; tot = 0
    for s in itertools.product([1, -1], repeat=n):
        tot += 1
        if (np.array(s) * np.abs(x)).mean() >= obs - 1e-15:
            cnt += 1
    return cnt / tot
for (ds, prog, ut), g in U.groupby(["dataset", "program", "unit_type"]):
    if ut == "donor_or_array":
        print(ds, prog, ut, "n", len(g), "positive", int((g.excess > 0).sum()), "mean", round(g.excess.mean(), 4), "exact sign-flip P", round(signflip_p(g.excess.values), 4))
