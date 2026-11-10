"""Critic check (1B): expected overlap and MH OR when expression strata are made finer than deciles,
and when strata use the colocalized genes' own expression distribution. Exact stratified two-sided p
uses the same distance-from-expectation rule as H57."""
import numpy as np, pandas as pd
from scipy.stats import hypergeom
R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
gl = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
deg["ens"] = deg.gene.str.split(".").str[0]; deg["DEG"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
U = deg.merge(gl, left_on="ens", right_on="ensembl")
U = U[U.susie_state_t12.isin(["positive", "tested_le_0.5"])].sort_values("ens").reset_index(drop=True)
U["gen"] = U.susie_state_t12 == "positive"
def test(b):
    O = int((U.gen & U.DEG).sum()); E = 0; pmf = np.array([1.0]); num = den = 0
    for k, g in U.groupby(b):
        n = len(g); G = int(g.gen.sum()); D = int(g.DEG.sum()); a = int((g.gen & g.DEG).sum())
        E += G * D / n
        pmf = np.convolve(pmf, hypergeom.pmf(np.arange(min(G, D) + 1), n, G, D))
        bb = G - a; c = D - a; d = n - G - D + a
        num += a * d / n; den += bb * c / n
    x = np.arange(len(pmf)); p = pmf[np.abs(x - E) >= abs(O - E) - 1e-9].sum()
    return O, round(E, 2), round(num / den, 3), round(p, 4)
for q in [10, 20, 50, 100]:
    print(f"{q} equal-count AveExpr bins:", test(pd.qcut(U.AveExpr.rank(method="first"), q, labels=False)))
# DEG rate by expression among non-colocalized genes vs colocalized genes within decile
U["dec"] = pd.qcut(U.AveExpr.rank(method="first"), 10, labels=False)
t = U.groupby(["dec", "gen"]).agg(n=("DEG", "size"), deg=("DEG", "sum"), aveexpr=("AveExpr", "median")).unstack("gen")
print(t.to_string())
