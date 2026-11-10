#!/usr/bin/env python3
"""Code review U6: Figure 1B expression-matched control (gray diamond), from METHODS.md:
426 non-colocalized U_eval genes matched 1:1 to the 426 colocalized genes on AveExpr by greedy nearest-neighbour
matching (no replacement; no caliper is stated), over 1,000 visiting orders. Count canonical DEGs among controls and the
expectation within the same 10 AveExpr bins. The legend's "exact p = 0.052" has no Methods definition, so three
candidate tests are computed: (1) exact stratified distance p for the control count, (2) exact McNemar on pairs,
(3) exact binomial 26 vs control count. Seeds: 20261023 (visiting orders). Usage: u6_fig1b_nc.py <out_dir>"""
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
u = b.merge(c, on="id", validate="one_to_one")
u = u[u.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
u["gen"] = u.coloc_best_susie_pp4 > 0.5
u = u.sort_values(["AveExpr", "id"]).reset_index(drop=True)
u["bin"] = pd.qcut(np.arange(len(u)), 10, labels=False)
rate = u.groupby("bin").deg.mean()
cases = u[u.gen]; pool = u[~u.gen]
ca = cases.AveExpr.to_numpy(); pa = pool.AveExpr.to_numpy(); pdeg = pool.deg.to_numpy(); pbin = pool.bin.to_numpy()

def match(order):
    used = np.zeros(len(pa), bool); pick = np.empty(len(order), int)
    for k, i in enumerate(order):
        d = np.abs(pa - ca[i]); d[used] = np.inf
        j = int(np.argmin(d)); used[j] = True; pick[k] = j
    out = np.empty(len(ca), int); out[order] = pick
    return out

rng = np.random.default_rng(20261023)
res = []
for r in range(1000):
    m = match(rng.permutation(len(ca)))
    res.append((int(pdeg[m].sum()), float(rate[pbin[m]].sum()), float(np.abs(pa[m] - ca).max()), m))
cnt = np.array([x[0] for x in res]); E = np.array([x[1] for x in res])
lines = [f"cases {len(ca)} pool {len(pa)} pool DEGs {int(pdeg.sum())}",
         f"1000 visiting orders (seed 20261023): control DEG count distribution {pd.Series(cnt).value_counts().to_dict()}",
         f"expected within bins: min {E.min():.3f} max {E.max():.3f}; max |dAveExpr| over runs {max(x[2] for x in res):.4f}"]
m = res[0][3]
# candidate p values for the first visiting order
cnt0 = int(pdeg[m].sum()); mb = pd.Series(pbin[m]).value_counts().sort_index()
tab = u.groupby("bin").agg(N=("deg", "size"), D=("deg", "sum"))
pmf = np.array([1.0])
for bb, row in tab.iterrows():
    k = int(mb.get(bb, 0)); x = np.arange(0, min(k, row.D) + 1)
    pmf = np.convolve(pmf, stats.hypergeom.pmf(x, row.N, row.D, k))
xs = np.arange(len(pmf)); E0 = float((pmf * xs).sum())
p_strat = float(pmf[np.abs(xs - E0) >= abs(cnt0 - E0) - 1e-9].sum())
cdeg = cases.deg.to_numpy()
b10 = int((cdeg & ~pdeg[m]).sum()); b01 = int((~cdeg & pdeg[m]).sum())
p_mcn = stats.binomtest(b10, b10 + b01, 0.5).pvalue
p_bin = stats.binomtest(cnt0, cnt0 + int(cdeg.sum()), 0.5).pvalue
lines += [f"run 1: control DEGs {cnt0}, E within bins (hypergeometric mean with control bin sizes) {E0:.3f}",
          f"(1) exact stratified distance p for control count: {p_strat:.4f}",
          f"(2) exact McNemar: case-only {b10} control-only {b01} both {int((cdeg & pdeg[m]).sum())} p {p_mcn:.4f}",
          f"(3) exact binomial {int(cdeg.sum())} vs {cnt0}: p {p_bin:.4f}"]
print("\n".join(lines)); (OUT / "fig1b_nc_log.txt").write_text("\n".join(lines) + "\n")
