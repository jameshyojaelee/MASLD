"""Leave-one-donor-out BH117 for the Fig 4C primary model (review R04)."""
import sys
import numpy as np
import pandas as pd
sys.argv = sys.argv[:1] + ["/dev/null"]
exec(open("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R04/rederive_primary.py").read().split("rows = {}")[0])
dd = load(META)
rows = []
for donor in sorted(dd.donor.unique()):
    sub = dd[dd.donor != donor]
    res = []
    for (ct, mod), d in sub.groupby(["cell_type", "module"]):
        f = fit(d); f.pop("resid"); res.append(dict(cell_type=ct, module=mod, **f))
    res = pd.DataFrame(res); res["q"] = bh(res.p, 117)
    r20 = res[(res.cell_type == "hepatocytes") & (res.module == 20)].iloc[0]
    rows.append(dict(dropped=donor, dataset=dd[dd.donor == donor].dataset.iloc[0],
                     stage=int(dd[dd.donor == donor].stage.iloc[0]), beta20=r20.beta, q20=r20.q,
                     n_calls=int((res.q < .05).sum())))
R = pd.DataFrame(rows)
R.to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R04/leave_one_donor_out.tsv", sep="\t", index=False)
print("single-donor deletions:", len(R), " Ductular q>=0.05 in", int((R.q20 >= .05).sum()))
print(R.sort_values("q20", ascending=False).head(8).round(4).to_string(index=False))
