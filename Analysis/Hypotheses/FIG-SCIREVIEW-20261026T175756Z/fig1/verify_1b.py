"""Skeptic re-derivation of Figure 1B (verifier, 2026-10-26). Written without reusing rederive_1b.py.

Inputs: adopted 774-unit bulk fit (deg_results.csv) and COLOC r2 tier-1/2 gene table.
Binning follows PRESPEC_H57 line 82: sort by Ensembl ID, rank AveExpr method=first, 10 equal-count bins.
MH OR and RBG CI computed by hand (not statsmodels).
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, nchypergeom_fisher

R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
gl = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
res = {}

deg["ens"] = deg["gene"].str.split(".").str[0]
deg["DEG"] = (deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.5)
res["bulk_rows"] = len(deg)
res["bulk_unique_ens"] = deg["ens"].nunique()
res["canonical_DEGs"] = int(deg["DEG"].sum())
res["gl_rows"] = len(gl)
res["gl_unique_ens"] = gl["ensembl"].nunique()
res["gl_state_counts"] = gl["susie_state_t12"].value_counts().to_dict()
res["multi_signal_evaluated_t12"] = int(gl["susie_state_t12"].isin(["positive", "tested_le_0.5"]).sum())

st = gl.set_index("ensembl")["susie_state_t12"]
deg["state"] = deg["ens"].map(st).fillna("no_record")
res["bulk_with_record"] = int((deg["state"] != "no_record").sum())
res["DEG_by_state"] = deg.loc[deg["DEG"], "state"].value_counts().to_dict()

# positives outside the bulk table
pos = set(gl.loc[gl.susie_state_t12 == "positive", "ensembl"])
res["positives_total"] = len(pos)
res["positives_in_bulk"] = len(pos & set(deg["ens"]))
ev = set(gl.loc[gl.susie_state_t12.isin(["positive", "tested_le_0.5"]), "ensembl"])
res["evaluated_not_in_bulk"] = len(ev - set(deg["ens"]))

U = deg[deg["state"].isin(["positive", "tested_le_0.5"])].copy()
U["gen"] = U["state"] == "positive"
N, G, D = len(U), int(U.gen.sum()), int(U.DEG.sum())
O = int((U.gen & U.DEG).sum())
res.update(N=N, G=G, D=D, O=O, E_unstrat=round(G * D / N, 4))


def strat(df, key):
    E = 0.0
    pmf = np.array([1.0])
    num = den = 0.0
    sP = sPQ = sQ = 0.0
    sR = sS = 0.0
    for _, s in df.groupby(key):
        n = len(s); g = int(s.gen.sum()); d = int(s.DEG.sum()); a = int((s.gen & s.DEG).sum())
        b = g - a; c = d - a; dd = n - g - d + a
        E += g * d / n
        x = np.arange(0, min(g, d) + 1)
        pmf = np.convolve(pmf, hypergeom.pmf(x, n, g, d))
        Rk = a * dd / n; Sk = b * c / n
        Pk = (a + dd) / n; Qk = (b + c) / n
        sR += Rk; sS += Sk
        sP += Pk * Rk; sPQ += Pk * Sk + Qk * Rk; sQ += Qk * Sk
    orr = sR / sS
    var = sP / (2 * sR ** 2) + sPQ / (2 * sR * sS) + sQ / (2 * sS ** 2)
    lo, hi = np.exp(np.log(orr) - 1.959964 * np.sqrt(var)), np.exp(np.log(orr) + 1.959964 * np.sqrt(var))
    p_min = float(pmf[pmf <= pmf[O] * (1 + 1e-7)].sum())
    cdf = np.cumsum(pmf)
    return dict(E=round(E, 4), OR=round(orr, 4), lo=round(lo, 4), hi=round(hi, 4), p_minlik=round(p_min, 4),
                q025=int(np.searchsorted(cdf, 0.025)), q975=int(np.searchsorted(cdf, 0.975)),
                mass_13_30=round(float(pmf[13:31].sum()), 4), pmf=pmf)


U = U.sort_values("ens").reset_index(drop=True)
U["r_ave"] = U["AveExpr"].rank(method="first")
U["bin10"] = pd.qcut(U["r_ave"], 10, labels=False)
s10 = strat(U, "bin10")
pmf10 = s10.pop("pmf")
res["strat_AveExpr10"] = s10
res["bin_sizes"] = U.groupby("bin10").agg(n=("gen", "size"), g=("gen", "sum"), d=("DEG", "sum")).astype(int).to_dict("list")

U["r_se"] = U["SE"].rank(method="first")
U["binSE"] = pd.qcut(U["r_se"], 10, labels=False)
sSE = strat(U, "binSE"); sSE.pop("pmf")
res["strat_SE10"] = sSE
U["a5"] = pd.qcut(U["r_ave"], 5, labels=False)
U["s5"] = pd.qcut(U["r_se"], 5, labels=False)
U["a5s5"] = U["a5"] * 5 + U["s5"]
sJ = strat(U, "a5s5"); sJ.pop("pmf")
res["strat_Ave5xSE5"] = sJ

# Power: unstratified Fisher noncentral hypergeometric. X = overlap; population N, G 'successes', D drawn.
x = np.arange(max(0, G + D - N), min(G, D) + 1)
p0 = hypergeom.pmf(x, N, G, D)
# equal-tail rejection region at 0.025 each tail
cdf = np.cumsum(p0); sf = np.cumsum(p0[::-1])[::-1]
lo_c = x[cdf <= 0.025].max(); hi_c = x[sf <= 0.025].min()
rej_eq = (x <= lo_c) | (x >= hi_c)
# min-likelihood rejection region: p(x) = sum p0[p0 <= p0[x]] <= 0.05
pminl = np.array([p0[p0 <= v * (1 + 1e-7)].sum() for v in p0])
rej_ml = pminl <= 0.05
res["power_eq_crit"] = [int(lo_c), int(hi_c), round(float(p0[rej_eq].sum()), 4)]
ml_lo = int(x[rej_ml & (x < 25)].max()); ml_hi = int(x[rej_ml & (x > 25)].min())
res["power_minlik_crit"] = [ml_lo, ml_hi, round(float(p0[rej_ml].sum()), 4)]


def power(orv, rej):
    pk = nchypergeom_fisher.pmf(x, N, G, D, orv)
    return float(pk[rej].sum())


grid = np.round(np.arange(1.60, 1.85, 0.01), 2)
res["power_up_eq"] = {str(o): round(power(o, rej_eq), 3) for o in grid}
res["power_up_ml"] = {str(o): round(power(o, rej_ml), 3) for o in grid}
grid_dn = np.round(np.arange(0.40, 0.52, 0.01), 2)
res["power_dn_eq"] = {str(o): round(power(o, rej_eq), 3) for o in grid_dn}
res["power_dn_ml"] = {str(o): round(power(o, rej_ml), 3) for o in grid_dn}

json.dump(res, open(__file__.replace(".py", ".json"), "w"), indent=1, default=int)
print(json.dumps({k: v for k, v in res.items() if not k.startswith("power_up") and not k.startswith("power_dn")}, indent=1, default=int))
print("eq up", res["power_up_eq"])
print("ml up", res["power_up_ml"])
print("eq dn", res["power_dn_eq"])
print("ml dn", res["power_dn_ml"])
