"""Blind re-check of F1 (Figure 1B overlap) from the producing release files.

Independent implementation; does not import H57 or repair code.
Seed 20260923.
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact, nchypergeom_fisher
from scipy.optimize import brentq

SEED = 20260923
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
COLOC = ROOT + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
BULK = ROOT + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
OUT = sys.argv[1]

out = {}
c = pd.read_csv(COLOC)
b = pd.read_csv(BULK)
c["g"] = c["ensembl"].astype(str).str.split(".").str[0]
b["g"] = b["gene"].astype(str).str.split(".").str[0]
out["coloc_rows"], out["coloc_unique"] = len(c), int(c["g"].nunique())
out["bulk_rows"], out["bulk_unique"] = len(b), int(b["g"].nunique())
b["deg"] = (b["padj"] < 0.05) & (b["logFC"].abs() > 0.5)
out["canonical_deg_all"] = int(b["deg"].sum())

u = b.merge(c, on="g", how="inner")
u["gen"] = u["coloc_best_susie_pp4"] >= 0.5
u["gen_gt"] = u["coloc_best_susie_pp4"] > 0.5
u["abf_gen"] = u["coloc_best_abf_pp4"] >= 0.5
out["U_full_N"] = len(u)
out["U_full_state_by_deg"] = pd.crosstab(u["susie_state_t12"], u["deg"]).to_dict()
out["gen_ge_equals_gt"] = bool((u["gen"] == u["gen_gt"]).all())
out["gen_equals_positive_state"] = bool((u["gen"] == (u["susie_state_t12"] == "positive")).all())


def counts(df, gcol="gen"):
    N = len(df); G = int(df[gcol].sum()); D = int(df["deg"].sum())
    O = int((df[gcol] & df["deg"]).sum())
    E = G * D / N
    t = [[O, G - O], [D - O, N - G - D + O]]
    OR, p = fisher_exact(t)
    return dict(N=N, G=G, D=D, O=O, E=E, neither=N - G - D + O, fisher_OR=OR, fisher_p=p,
                nonoverlap_obs=(G - O) / G, nonoverlap_exp_unstrat=1 - D / N,
                marker_gen=E / G, marker_deg=E / D)


def add_bins(df, k):
    df = df.sort_values("g").reset_index(drop=True)
    r = df["AveExpr"].rank(method="first")
    df["bin"] = pd.qcut(r, k, labels=False)
    return df


def strat(df, k, gcol="gen"):
    df = add_bins(df, k)
    tabs = []
    pmf = np.array([1.0])
    E = 0.0
    num = den = 0.0
    # RBG pieces
    sPR = sPS_QR = sQS = 0.0
    R_ = S_ = 0.0
    per_bin = []
    for kk, s in df.groupby("bin"):
        n = len(s); G = int(s[gcol].sum()); D = int(s["deg"].sum())
        a = int((s[gcol] & s["deg"]).sum())
        bb, cc = G - a, D - a
        d = n - G - D + a
        per_bin.append(dict(bin=int(kk), n=n, G=G, D=D, O=a,
                            aveexpr_min=float(s.AveExpr.min()), aveexpr_max=float(s.AveExpr.max())))
        E += G * D / n
        lo, hi = max(0, G + D - n), min(G, D)
        x = np.arange(0, hi + 1)
        pm = hypergeom.pmf(x, n, G, D)
        pmf = np.convolve(pmf, pm)
        Ri = a * d / n; Si = bb * cc / n
        Pi = (a + d) / n; Qi = (bb + cc) / n
        R_ += Ri; S_ += Si
        sPR += Pi * Ri; sPS_QR += Pi * Si + Qi * Ri; sQS += Qi * Si
    O = int((df[gcol] & df["deg"]).sum())
    mh = R_ / S_
    var = sPR / (2 * R_ ** 2) + sPS_QR / (2 * R_ * S_) + sQS / (2 * S_ ** 2)
    lo_, hi_ = np.exp(np.log(mh) - 1.959963984540054 * np.sqrt(var)), np.exp(np.log(mh) + 1.959963984540054 * np.sqrt(var))
    pmf = pmf / pmf.sum()
    p_minlike = float(pmf[pmf <= pmf[O] * (1 + 1e-7)].sum())
    p_double = float(min(1.0, 2 * min(pmf[:O + 1].sum(), pmf[O:].sum())))
    null_mean = float((np.arange(len(pmf)) * pmf).sum())
    G = int(df[gcol].sum()); D = int(df["deg"].sum())
    return dict(k=k, O=O, E_strat=E, null_mean=null_mean, O_over_E=O / E, MH_OR=mh, MH_lo=lo_, MH_hi=hi_,
                p_exact_minlike=p_minlike, p_exact_doubled=p_double,
                marker_gen=E / G, marker_deg=E / D), per_bin, df


ev = u[u["susie_state_t12"].isin(["positive", "tested_le_0.5"])].copy()
ev7 = u[u["susie_state_t12"].isin(["positive", "tested_le_0.5", "untestable"])].copy()
out["U_full"] = counts(u)
out["U_eval"] = counts(ev)
out["U_eval7"] = counts(ev7)
out["no_result_DEGs"] = int((u["deg"] & (u["susie_state_t12"] == "no_susie_result")).sum())
out["no_result_genes"] = int((u["susie_state_t12"] == "no_susie_result").sum())
out["untestable_genes"] = int((u["susie_state_t12"] == "untestable").sum())
out["deg_tested_le"] = int((u["deg"] & (u["susie_state_t12"] == "tested_le_0.5")).sum())

T1, perbin, evb = strat(ev, 10)
out["T1_U_eval_10"] = T1
out["T1_per_bin"] = perbin
for k in (5, 20):
    out[f"U_eval_{k}"] = strat(ev, k)[0]
out["U_eval7_10"] = strat(ev7, 10)[0]
out["U_full_10"] = strat(u, 10)[0]

# D8: ABF map in U_full
out["D8_unstrat"] = counts(u, "abf_gen")
out["D8_10"] = strat(u, 10, "abf_gen")[0]

# Power: unstratified exact equal-tailed (alpha 0.025 each tail) under central hypergeometric,
# power under Fisher noncentral hypergeometric.
N, G, D = out["U_eval"]["N"], out["U_eval"]["G"], out["U_eval"]["D"]
x = np.arange(0, min(G, D) + 1)
cdf = hypergeom.cdf(x, N, G, D)
sf = hypergeom.sf(x - 1, N, G, D)
cL = int(x[cdf <= 0.025].max())
cU = int(x[sf <= 0.025].min())
size = float(hypergeom.cdf(cL, N, G, D) + hypergeom.sf(cU - 1, N, G, D))


def power(OR):
    return float(nchypergeom_fisher.cdf(cL, N, G, D, OR) + nchypergeom_fisher.sf(cU - 1, N, G, D, OR))


out["power"] = dict(cL=cL, cU=cU, size=size,
                    OR80_up=brentq(lambda o: power(o) - 0.8, 1.0, 5.0),
                    OR80_dn=brentq(lambda o: power(o) - 0.8, 0.05, 1.0),
                    power_at_1295=power(1.2947), power_at_0835=power(0.835))

# NC-1: 426 non-colocalized U_eval genes matched 1:1 on AveExpr, greedy nearest neighbour,
# without replacement, for several visiting orders of the genetic genes.
rng = np.random.default_rng(SEED)
gen_df = evb[evb["gen"]]
pool = evb[~evb["gen"]].reset_index(drop=True)
bin_D_over_N = evb.groupby("bin").apply(lambda s: s["deg"].sum() / len(s))
nc = []
orders = [("ascending", np.argsort(gen_df["AveExpr"].values, kind="stable"))]
orders += [(f"random_{i}", rng.permutation(len(gen_df))) for i in range(200)]
pa = pool["AveExpr"].values
for name, ordr in orders:
    used = np.zeros(len(pool), bool)
    picks = []
    gaps = []
    for v in gen_df["AveExpr"].values[ordr]:
        dist = np.abs(pa - v)
        dist[used] = np.inf
        j = int(np.argmin(dist))
        used[j] = True
        picks.append(j)
        gaps.append(dist[j])
    m = pool.iloc[picks]
    O = int(m["deg"].sum())
    E = float(sum(bin_D_over_N[b_] for b_ in m["bin"]))
    nc.append(dict(order=name, O=O, E_strat=E, OE=O / E, max_gap=float(max(gaps))))
ncd = pd.DataFrame(nc)
ncd.to_csv(OUT + "/nc_matched_recheck.tsv", sep="\t", index=False)
out["NC1"] = dict(n_orders=len(ncd), O_values=sorted(ncd["O"].unique().tolist()),
                  E_strat_min=float(ncd.E_strat.min()), E_strat_max=float(ncd.E_strat.max()),
                  OE_mean=float(ncd.OE.mean()), max_gap_median=float(ncd.max_gap.median()))

# Figure 3G universe and display fractions
out["fig1b_display"] = dict(
    gen_both_frac=26 / out["U_eval"]["G"], gen_only_frac=(out["U_eval"]["G"] - 26) / out["U_eval"]["G"],
    deg_both_frac=26 / out["U_eval"]["D"], deg_only_frac=(out["U_eval"]["D"] - 26) / out["U_eval"]["D"])

pd.DataFrame(perbin).to_csv(OUT + "/bins_recheck.tsv", sep="\t", index=False)
json.dump(out, open(OUT + "/recheck.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
