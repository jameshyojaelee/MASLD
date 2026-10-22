"""Gate L01-04 MDE. Locus-stratified one-sided CMH test (score test of conditional logistic
regression with a binary predictor) using the real per-locus counts of expressed COLOC and
non-COLOC window genes from job 3051 locus_testability.tsv. Seed 20261015."""
import numpy as np, pandas as pd
from scipy import stats
J = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3051"
t = pd.read_csv(f"{J}/locus_testability.tsv", sep="\t")
t = t[t.peak_passes_5pct & (t.coloc_genes_expressed_5pct > 0) & (t.noncoloc_expressed_5pct > 0)]
rng = np.random.default_rng(20261015)
B = 4000

def power(sub, p0, OR, alpha=0.05):
    c = sub.coloc_genes_expressed_5pct.to_numpy(); o = sub.noncoloc_expressed_5pct.to_numpy(); n = c + o
    p1 = OR * p0 / (1 - p0 + OR * p0)
    a = rng.binomial(c[None, :], p1, size=(B, len(c)))   # linked COLOC genes per locus
    b = rng.binomial(o[None, :], p0, size=(B, len(c)))   # linked other genes per locus
    m = a + b
    E = c * m / n
    V = c * o * m * (n - m) / (n ** 2 * np.maximum(n - 1, 1))
    num = (a - E).sum(1); den = np.sqrt(V.sum(1))
    z = np.where(den > 0, num / np.where(den > 0, den, 1), 0)
    return float((z > stats.norm.isf(alpha)).mean())

for label, thr in [("in-peak SNP.PP.H4>=0.2", 0.2), ("in-peak SNP.PP.H4>=0.5", 0.5)]:
    sub = t[t.max_pp >= thr]
    print(f"{label}: testable loci {len(sub)}, COLOC genes {sub.coloc_genes_expressed_5pct.sum()}, other genes {sub.noncoloc_expressed_5pct.sum()}")
    for p0 in (0.05, 0.10, 0.20):
        grid = [1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15]
        pw = [power(sub, p0, OR) for OR in grid]
        mde = next((OR for OR, p in zip(grid, pw) if p >= 0.8), None)
        p1 = (mde * p0 / (1 - p0 + mde * p0)) if mde else None
        print(f"   p0={p0}: power by OR {dict(zip(grid, np.round(pw, 2)))} -> MDE OR {mde}, COLOC link rate {None if p1 is None else round(p1, 3)}")
# Candidate's binomial sketch, recomputed
for n in (51, 34, 19):
    for p0 in (0.10,):
        k = int(stats.binom.isf(0.05, n, p0)) + 1
        while stats.binom.sf(k - 1, n, p0) > 0.05: k += 1
        p1 = next(p for p in np.arange(0.10, 0.80, 0.005) if stats.binom.sf(k - 1, n, p) >= 0.8)
        print(f"binomial sketch n={n}, p0={p0}: critical count {k}, size {stats.binom.sf(k-1,n,p0):.4f}, MDE COLOC link rate {p1:.3f}")
