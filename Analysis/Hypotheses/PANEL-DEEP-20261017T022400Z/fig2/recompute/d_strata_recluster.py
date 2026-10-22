"""Fig 2D ancestry/subclass strata with 1-Mb clusters re-derived inside each stratum.

Input: out_D/D_gated_instances_with_clusters.tsv.gz written by d_recompute.py
(restored sentinel positions, class masses). Own bootstrap, seed 20261017.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
a = pd.read_csv(f"{HERE}/out_D/D_gated_instances_with_clusters.tsv.gz", sep="\t", dtype={"chr": str})

def cluster(df, w=1_000_000):
    lab = pd.Series(index=df.index, dtype=object)
    for c, g in df.groupby("chr"):
        g = g.sort_values("pos")
        cid, last = 0, None
        for i, p in zip(g.index, g.pos):
            if last is not None and p - last > w:
                cid += 1
            lab[i] = f"{c}_{cid}"
            last = p
    return lab

def frac(d):
    return 100 * d.groupby("cl").pa_spl.mean().mean() if len(d) else np.nan

rng = np.random.default_rng(20261017)
def boot(df, B=10000):
    cm = df.groupby(["cl", "cls"]).pa_spl.mean().unstack()
    v = cm[["direct", "enzyme"]].values
    idx = rng.integers(0, len(cm), size=(B, len(cm)))
    out = np.array([100 * (np.nanmean(v[i, 0]) - np.nanmean(v[i, 1])) for i in idx])
    return np.nanpercentile(out, [2.5, 97.5]), np.nanpercentile(out, [100 / 120, 100 - 100 / 120])

rows = []
for name, mask in [("EUR", a.ancestry == "EUR"), ("AMR", a.ancestry == "AMR"), ("nonEUR", a.ancestry != "EUR"),
                   ("EUR_and_AMR_union", a.ancestry.isin(["EUR", "AMR"]))]:
    s = a[mask].copy()
    s["cl"] = cluster(s)
    d, e = s[s.cls == "direct"], s[s.cls == "enzyme"]
    ci95, ci98 = boot(s)
    rows.append({"stratum": name, "clusters_in_stratum": s.cl.nunique(), "direct_pct": frac(d), "direct_sets": len(d),
                 "direct_clusters": d.cl.nunique(), "enzyme_pct": frac(e), "enzyme_sets": len(e),
                 "enzyme_clusters": e.cl.nunique(), "diff_pp": frac(d) - frac(e),
                 "own_lo95": ci95[0], "own_hi95": ci95[1], "own_lo98.333": ci98[0], "own_hi98.333": ci98[1]})
r = pd.DataFrame(rows).set_index("stratum")
fw = {"direct_pct": 16 / 26 * r.loc["EUR", "direct_pct"] + 10 / 26 * r.loc["AMR", "direct_pct"],
      "enzyme_pct": 16 / 26 * r.loc["EUR", "enzyme_pct"] + 10 / 26 * r.loc["AMR", "enzyme_pct"]}
fw["diff_pp"] = fw["direct_pct"] - fw["enzyme_pct"]
r.loc["EUR+AMR fixed 16/26,10/26 (within-stratum clusters)", list(fw)] = list(fw.values())
# EUR+AMR fixed-weight interval: separate-stratum cluster draws (own)
print(r.round(3).to_string())
r.to_csv(f"{HERE}/out_D/D_strata_reclustered.tsv", sep="\t")
# which EUR direct clusters differ between pooled and within-EUR clustering
s = a[a.ancestry == "EUR"].copy(); s["cl"] = cluster(s)
x = s[s.cls == "direct"].groupby("cluster").cl.nunique()
print("pooled direct clusters that split inside EUR:", x[x > 1].to_dict())
