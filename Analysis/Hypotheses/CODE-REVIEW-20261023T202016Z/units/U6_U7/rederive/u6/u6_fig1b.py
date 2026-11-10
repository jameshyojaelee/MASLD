#!/usr/bin/env python3
"""Code review U6: independent Figure 1B overlap from METHODS.md ("Joint testing universe and evidence classes").
Universe: stable-ID join of the 24,115-gene bulk fit and the COLOC r2 Tier-1/2 gene table (14,704), restricted to
susie_state_t12 in {positive, tested_le_0.5}. Genetic = SuSiE PP.H4 > 0.5; DEG = BH q < 0.05 and |log2FC| > 0.5.
Strata: 10 AveExpr quantile bins (ties broken by Ensembl ID). Null: exact convolution of per-bin hypergeometrics;
two-sided p = sum of P(x) with |x - E| >= |O - E|. Effect: Mantel-Haenszel OR with Robins-Breslow-Greenland CI.
Usage: u6_fig1b.py <out_dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
b = pd.read_csv(ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
c = pd.read_csv(ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
b["id"] = b.gene.str.split(".").str[0]; c["id"] = c.ensembl.astype(str).str.split(".").str[0]
b["deg"] = (b.padj < 0.05) & (b.logFC.abs() > 0.5)
j = b.merge(c, on="id", how="inner", validate="one_to_one")
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)
say("bulk", len(b), "COLOC", len(c), "join", len(j))
say("join by SuSiE state", j.susie_state_t12.value_counts().to_dict())
u = j[j.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
u["gen"] = u.coloc_best_susie_pp4 > 0.5
say("U_eval", len(u), "genetic", int(u.gen.sum()), "DEG", int(u.deg.sum()), "both", int((u.gen & u.deg).sum()))
say("canonical DEGs: in U_eval", int(u.deg.sum()), "; join with no SuSiE result (ABF only)", int(j[j.susie_state_t12 == "no_susie_result"].deg.sum()),
    "; join SuSiE-untestable", int(j[j.susie_state_t12 == "untestable"].deg.sum()), "; bulk not in COLOC", int(b[~b.id.isin(c.id)].deg.sum()))
E0 = u.gen.sum() * u.deg.sum() / len(u)
say("unstratified expected", round(E0, 4))

def binning(df, how):
    d = df.sort_values(["AveExpr", "id"]).reset_index(drop=True)
    n = len(d); r = np.arange(n)
    if how == "floor":   d["bin"] = (r * 10) // n
    elif how == "qcut":  d["bin"] = pd.qcut(r, 10, labels=False)
    elif how == "ceil":  d["bin"] = np.ceil((r + 1) * 10 / n).astype(int) - 1
    return d
def analyse(d, tag):
    tabs = d.groupby("bin").apply(lambda g: pd.Series({"N": len(g), "K": int(g.gen.sum()), "n": int(g.deg.sum()), "a": int((g.gen & g.deg).sum())}))
    E = float((tabs.K * tabs.n / tabs.N).sum()); O = int(tabs.a.sum())
    pmf = np.array([1.0])
    for _, t in tabs.iterrows():
        lo = max(0, t.n - (t.N - t.K)); hi = min(t.K, t.n)
        p = np.zeros(int(hi) + 1); x = np.arange(lo, hi + 1)
        p[int(lo):] = stats.hypergeom.pmf(x, t.N, t.K, t.n)
        pmf = np.convolve(pmf, p)
    xs = np.arange(len(pmf))
    p2 = float(pmf[np.abs(xs - E) >= abs(O - E) - 1e-9].sum())
    cdf = np.cumsum(pmf)
    lo95 = int(xs[np.searchsorted(cdf, 0.025)]); hi95 = int(xs[np.searchsorted(cdf, 0.975)])
    # MH OR and RBG CI; per stratum a=gen&deg, b=gen&!deg, c=!gen&deg, d=!gen&!deg
    A = tabs.a; B = tabs.K - tabs.a; C = tabs.n - tabs.a; D = tabs.N - tabs.K - tabs.n + tabs.a; N = tabs.N
    R = A * D / N; S = B * C / N; P = (A + D) / N; Q = (B + C) / N
    OR = R.sum() / S.sum()
    v = (P * R).sum() / (2 * R.sum() ** 2) + ((P * S) + (Q * R)).sum() / (2 * R.sum() * S.sum()) + (Q * S).sum() / (2 * S.sum() ** 2)
    ci = np.exp(np.log(OR) + np.array([-1, 1]) * 1.959963984540054 * np.sqrt(v))
    say(f"[{tag}] bins sizes {tabs.N.astype(int).tolist()} | O {O} E {E:.4f} O/E {O/E:.5f} | exact two-sided p {p2:.4f} | null central 95% {lo95}-{hi95} | MH OR {OR:.4f} ({ci[0]:.4f}-{ci[1]:.4f})")
    return tabs
for how in ("floor", "qcut", "ceil"):
    t = analyse(binning(u, how), how)
t.to_csv(OUT / "fig1b_bins.tsv", sep="\t")
(OUT / "fig1b_log.txt").write_text("\n".join(out) + "\n")
