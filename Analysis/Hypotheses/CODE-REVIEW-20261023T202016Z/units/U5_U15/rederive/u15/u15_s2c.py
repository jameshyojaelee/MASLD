#!/usr/bin/env python
"""U15 S2C: merged loci (S2B rule, recomputed here) with a result from all four fine-mappers.
Cross-ancestry loci (SuSiEx, meSuSiE; MVP ALT/AST/NAFLD) are matched to a merged locus when the
locus position lies within 500 kb of the merged window. Best PIP per method per merged locus."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = R / "GWAS/finemapping/results"
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
t = pd.read_csv(R / "GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
MAIN = set(t[t.placement == "main"].study_name)
cr = pd.read_csv(FM / "credible_sets.csv", usecols=["chromosome", "position", "locus", "study", "susie_pip", "carma_pip"])
cr = cr[cr.study.isin(MAIN)]
w = cr.groupby(["study", "locus"]).agg(chr=("chromosome", "first"), lo=("position", "min"), hi=("position", "max")).reset_index()
w["a"], w["b"] = w.lo - 250_000, w.hi + 250_000
mid, rows = {}, []
for ch, g in w.sort_values(["chr", "a"]).groupby("chr"):
    cur = None
    for r in g.itertuples(index=False):
        if cur is None or r.a > cur["b"]:
            if cur is not None:
                rows.append(cur)
            cur = dict(chr=int(ch), a=r.a, b=r.b, lo=r.lo, hi=r.hi, keys=[(r.study, r.locus)])
        else:
            cur["b"] = max(cur["b"], r.b); cur["hi"] = max(cur["hi"], r.hi); cur["lo"] = min(cur["lo"], r.lo)
            cur["keys"].append((r.study, r.locus))
    rows.append(cur)
ml = pd.DataFrame(rows)
for i, r in ml.iterrows():
    for k in r["keys"]:
        mid[k] = i
cr["m"] = [mid[(s, l)] for s, l in zip(cr.study, cr.locus)]
best = cr.groupby("m").agg(susie=("susie_pip", "max"), carma=("carma_pip", "max"))
ml = ml.join(best)
print("merged loci", len(ml), "with CARMA", int(ml.carma.notna().sum()))

sx = pd.read_csv(FM / "susiex_mvp/susiex_cs_mvp.csv")
sx = sx[sx.trait_pair.isin(["ALT", "AST", "NAFLD"])]
sxl = sx.groupby(["locus_id", "trait_pair"]).OVRL_PIP.max().reset_index()
m = sxl.locus_id.str.extract(r"chr(\d+)_(\d+)$").astype(int)
sxl["chr"], sxl["pos"] = m[0], m[1]
me = pd.read_csv(FM / "mesusie_mvp/mesusie_locus_summary_mvp.csv")
me = me[me.trait_pair.isin(["ALT", "AST", "NAFLD"]) & (me.converged == True)].copy()  # noqa: E712
me["pos"] = ((me.window_start + me.window_end) // 2).astype(int)
me["chr"] = me.chr.astype(int)


def match(tab, col, tol=500_000):
    out = []
    for r in ml.itertuples():
        s = tab[(tab.chr == r.chr) & (tab.pos >= r.lo - tol) & (tab.pos <= r.hi + tol)]
        out.append(s[col].max() if len(s) else np.nan)
    return out


ml["susiex"] = match(sxl, "OVRL_PIP")
ml["mesusie"] = match(me, "max_pip")
common = ml.dropna(subset=["susie", "carma", "susiex", "mesusie"])
print("common loci (all four methods)", len(common))


def band(x):
    return pd.cut(x, [-np.inf, 0.5, 0.9, np.inf], right=False, labels=["<0.5", "0.5-0.9", ">=0.9"]).value_counts().to_dict()


for c in ["susie", "carma", "susiex", "mesusie"]:
    b = band(common[c])
    print(c, b, "share >=0.9 %.0f%%" % (100 * b[">=0.9"] / len(common)))
common.drop(columns="keys").to_csv(OUT / "S2C_common_loci_rederived.tsv", sep="\t", index=False)
