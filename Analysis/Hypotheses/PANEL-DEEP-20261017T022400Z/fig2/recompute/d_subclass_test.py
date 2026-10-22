"""Direct test for the Fig 2D-subclass claim: PDFF minus NAFLD+NASH protein-altering + splice PIP share
(pooled 1-Mb clusters, as in the published subclass rows); cluster bootstrap, seed 20261017, 20,000 draws."""
import numpy as np, pandas as pd
a = pd.read_csv("out_D/D_gated_instances_with_clusters.tsv.gz", sep="\t", dtype={"chr": str})
d = a[a.cls == "direct"].copy(); d["sub"] = np.where(d.trait == "PDFF", "PDFF", "NAFLD+NASH")
frac = lambda x: 100 * x.groupby("cluster").pa_spl.mean().mean()
p, n = d[d["sub"] == "PDFF"], d[d["sub"] == "NAFLD+NASH"]
rng = np.random.default_rng(20261017)
cp = p.groupby("cluster").pa_spl.mean().values; cn = n.groupby("cluster").pa_spl.mean().values
bs = np.array([100 * (cp[rng.integers(0, len(cp), len(cp))].mean() - cn[rng.integers(0, len(cn), len(cn))].mean()) for _ in range(20000)])
ne = n[n.ancestry == "EUR"]
out = pd.DataFrame([{"PDFF_pct": frac(p), "PDFF_sets": len(p), "PDFF_clusters": p.cluster.nunique(), "NAFLDNASH_pct": frac(n), "NAFLDNASH_sets": len(n),
                     "NAFLDNASH_clusters": n.cluster.nunique(), "diff_pp": frac(p) - frac(n), "lo95": np.percentile(bs, 2.5), "hi95": np.percentile(bs, 97.5),
                     "EUR_NAFLDNASH_pct": frac(ne), "EUR_NAFLDNASH_sets": len(ne), "EUR_NAFLDNASH_clusters": ne.cluster.nunique()}])
out.to_csv("out_D/D_subclass_direct_test.tsv", sep="\t", index=False); print(out.round(2).T)
