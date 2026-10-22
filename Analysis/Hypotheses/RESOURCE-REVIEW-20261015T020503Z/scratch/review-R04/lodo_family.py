"""Family-wide (BH117) leave-one-dataset-out and per-dataset slope heterogeneity (review R04)."""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.argv = sys.argv[:1] + ["/dev/null"]
exec(open(__file__.replace("lodo_family.py", "rederive_primary.py")).read().split("rows = {}")[0])

dd = load(META)
out = []
for drop in [None] + sorted(dd.dataset.unique()):
    sub = dd if drop is None else dd[dd.dataset != drop]
    res = []
    for (ct, mod), d in sub.groupby(["cell_type", "module"]):
        if d.groupby("dataset").stage.nunique().max() < 2:
            continue
        f = fit(d)
        f.pop("resid")
        res.append(dict(cell_type=ct, module=mod, **f))
    res = pd.DataFrame(res)
    res["q"] = bh(res.p, 117)
    res["hc3_q"] = bh(res.hc3_p, 117)
    res["dropped"] = drop or "none"
    out.append(res)
    hit = res[res.q < 0.05]
    foc = res[(res.cell_type == "hepatocytes") & res.module.isin([8, 20])][["module", "beta", "q", "hc3_q", "n"]]
    print(f"drop={drop}: programs={len(res)} q<0.05={len(hit)} {list(zip(hit.cell_type, hit.module))} "
          f"| focal: " + "; ".join(f"m{r.module} b={r.beta:.3f} q={r.q:.4f} hc3q={r.hc3_q:.3f} n={r.n}" for r in foc.itertuples()))
pd.concat(out).to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R04/lodo_family_bh117.tsv", sep="\t", index=False)

# Inverse-variance meta of within-dataset slopes (datasets with >=2 stage levels and >=4 donors)
print("\nwithin-dataset slopes, IVW fixed-effect meta and Cochran Q (hepatocyte focal programs)")
for mod in [20, 8]:
    d = dd[(dd.cell_type == "hepatocytes") & (dd.module == mod)]
    sl = []
    for ds, g in d.groupby("dataset"):
        if g.stage.nunique() > 1 and len(g) >= 4:
            r = stats.linregress(g.stage, g.score)
            sl.append((ds, len(g), r.slope, r.stderr))
    s = pd.DataFrame(sl, columns=["dataset", "n", "slope", "se"])
    w = 1 / s.se ** 2
    b = (w * s.slope).sum() / w.sum()
    se = 1 / np.sqrt(w.sum())
    Q = (w * (s.slope - b) ** 2).sum()
    print(s.round(4).to_string(index=False))
    print(f"  module {mod}: IVW beta={b:.4f} se={se:.4f} z p={2*stats.norm.sf(abs(b/se)):.4g}; "
          f"Cochran Q={Q:.3f} df={len(s)-1} p={stats.chi2.sf(Q, len(s)-1):.4g}")
