"""3x3 SuSiE state cross-tab (enzyme class x direct class) per gene over the 35 main tier-1/2 studies.
State per class: positive (any PP.H4.susie > 0.5), tested_le_0.5 (method == 'susie', all <= 0.5), none (no SuSiE result)."""
import pandas as pd, numpy as np, json
from pathlib import Path
from statsmodels.stats.proportion import proportion_confint
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
REL = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
OUT = ROOT / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1"
tier = pd.read_csv(ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
main = tier[(tier.placement == "main") & tier.tier.isin([1, 2])]
a = pd.read_csv(REL / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "ensembl", "PP.H4.susie", "method"], low_memory=False)
a = a[a.gwas_name.isin(main.study_name)].merge(main[["study_name", "tier_label"]], left_on="gwas_name", right_on="study_name")
na = a[a.ensembl.isna()]
NA_INFO = dict(rows_missing_ensembl=len(na), pos_rows_missing_ensembl=int((na["PP.H4.susie"] > 0.5).sum()))
a = a[a.ensembl.notna()]
t = a[a.method == "susie"]
st = {}
for lab in ("liver_enzyme", "direct_MASLD"):
    g = t[t.tier_label == lab].groupby("ensembl")["PP.H4.susie"].max()
    st[lab] = g
genes = sorted(set(a.ensembl))
def state(lab, e):
    g = st[lab]
    if e not in g.index: return "none"
    return "positive" if g[e] > 0.5 else "tested_le_0.5"
df = pd.DataFrame({"ensembl": genes})
df["enzyme"] = [state("liver_enzyme", e) for e in genes]
df["direct"] = [state("direct_MASLD", e) for e in genes]
ct = pd.crosstab(df.enzyme, df.direct, margins=True)
print(ct)
ct.to_csv(OUT / "recompute_1c_crosstab.tsv", sep="\t")
res = {"crosstab": ct.to_dict(), **NA_INFO}
both_tested = df[(df.enzyme != "none") & (df.direct != "none")]
res["n_tested_both_classes"] = len(both_tested)
res["both_tested_enzyme_pos"] = int((both_tested.enzyme == "positive").sum())
res["both_tested_direct_pos"] = int((both_tested.direct == "positive").sum())
for lab, col in (("liver_enzyme", "enzyme"), ("direct_MASLD", "direct")):
    n = int((df[col] != "none").sum()); k = int((df[col] == "positive").sum())
    lo, hi = proportion_confint(k, n, method="wilson")
    res[f"{lab}_pos_of_tested"] = [k, n, k / n, lo, hi]
    nb = len(both_tested); kb = int((both_tested[col] == "positive").sum())
    lo, hi = proportion_confint(kb, nb, method="wilson")
    res[f"{lab}_pos_of_tested_in_both"] = [kb, nb, kb / nb, lo, hi]
json.dump(res, open(OUT / "recompute_1c_crosstab.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in res.items() if k != "crosstab"}, indent=1))
