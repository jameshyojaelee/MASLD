"""Recompute Figure 5F focal p, q, excess from the stored r2 null draws; sign-flip p from stored per-donor excess;
ranks and content-lineage q from AN4 tables; percent of excess removed."""
import itertools, sys
import numpy as np, pandas as pd
R2, W, J = sys.argv[1], sys.argv[2], sys.argv[3]
res = pd.read_csv(f"{R2}/full_composition_program_results.tsv", sep="\t")
nd = pd.read_csv(f"{R2}/full_composition_null_draws.tsv.gz", sep="\t")
print("null draw columns", nd.columns.tolist())
val = [c for c in nd.columns if "residual" in c and "moran" in c][0]
rows = []
for r in res.itertuples():
    d = nd[(nd.dataset == r.dataset) & (nd.program_id == r.program_id)][val].to_numpy(float)
    k = int((d >= r.residual_moran_i).sum())
    rows.append({"dataset": r.dataset, "program": r.program_name, "B": len(d), "k": k, "p": (k + 1) / (len(d) + 1),
                 "rec_p": r.residual_pvalue, "excess": r.residual_moran_i - d.mean(), "rec_q": r.residual_padj})
t = pd.DataFrame(rows)
for ds, g in t.groupby("dataset"):
    p = g.p.to_numpy(); o = np.argsort(p); a = p[o] * 2 / np.arange(1, 3); a = np.minimum.accumulate(a[::-1])[::-1]
    q = np.empty(2); q[o] = np.minimum(a, 1); t.loc[g.index, "q"] = q
print(t.round(6).to_string(index=False))
def sf(e):
    e = np.asarray(e); s = np.array(list(itertools.product([1, -1], repeat=len(e))))
    return np.mean(s @ np.abs(e) >= e.sum() - 1e-12)
u = pd.read_csv(f"{J}/unit_level.tsv", sep="\t")
u = u[(u.design == "all16") & (u.unit_type == "donor_or_array")]
for (ds, pr), g in u.groupby(["dataset", "program_name"]):
    print(ds, pr, "units", len(g), "excess>0", int((g.excess > 0).sum()), "raw residual I>0", int((g.residual_moran_i > 0).sum()),
          "signflip p", sf(g.excess.to_numpy()))
j = pd.read_csv(f"{J}/program_level.tsv", sep="\t")
for (ds, pr), g in j.groupby(["dataset", "program_name"]):
    lin = g[g.design == "lineage_only"].excess.iloc[0]; a16 = g[g.design == "all16"].excess.iloc[0]
    print(ds, pr, "excess hep-only %.4f all16 %.4f removed %.1f%%" % (lin, a16, 100 * (1 - a16 / lin)))
f = pd.read_csv(f"{W}/combined/family_all_programs.tsv", sep="\t"); f = f[f.testable]
for ds, g in f.groupby("dataset"):
    g = g.assign(r=g.excess.rank(ascending=False, method="min"))
    print(ds, "Stromal rank", g.loc[g.program_id == "hotspot_hepatocytes_f05c535ae5bbc0b9", "r"].iloc[0],
          "Ductular rank", g.loc[g.program_id == "hotspot_hepatocytes_48f39dd4d817a10e", "r"].iloc[0],
          "excess>0", int((g.excess > 0).sum()))
fl = pd.read_csv(f"{W}/focal_lineage/focal_lineage_results.tsv", sep="\t")
fl = fl[fl.scheme == "content_lineage"]
for ds, g in fl.groupby("dataset"):
    p = g.residual_pvalue.to_numpy(); o = np.argsort(p); a = p[o] * 2 / np.arange(1, 3); a = np.minimum.accumulate(a[::-1])[::-1]
    q = np.empty(2); q[o] = np.minimum(a, 1)
    print(ds, list(zip(g.program_name, g.excess.round(4), p, q.round(4), g.residual_padj_two_program)))
