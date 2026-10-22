"""MDE from job 3091 null draws (gate-L03-03). Observed agreement never computed."""
import numpy as np, pandas as pd
from scipy.stats import norm
J = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3091/out/'
d = pd.read_csv(J + 'null_draws.tsv', sep='\t'); s = pd.read_csv(J + 'null_summary.tsv', sep='\t')
za, zb = norm.ppf(1 - 0.05 / 4), norm.ppf(0.8)   # Holm first step (one-sided 0.0125), 80% power
for _, r in s.iterrows():
    x = d[d.lineage == r.lineage].agree.values
    upfrac = r.calls_up / r.calls_expressed
    emp = np.quantile(x, 1 - 0.05 / 4)
    print(f"{r.lineage:18s} B={len(x):4d} up_frac={upfrac:.3f} null mean={x.mean():.3f} sd={x.std(ddof=1):.3f} "
          f"q95={np.quantile(x,0.95):.3f} q98.75={emp:.3f} | MDE normal={x.mean()+(za+zb)*x.std(ddof=1):.3f} "
          f"MDE empirical(q98.75+0.84sd)={emp+zb*x.std(ddof=1):.3f} | draws>=0.65: {np.mean(x>=0.65):.3f} draws<=0.35: {np.mean(x<=0.35):.3f}")
    h, e = np.histogram(x, bins=np.arange(0.1, 0.85, 0.05)); print('   hist', dict(zip(np.round(e[:-1], 2), h)))
# per-gene partial-r MDE (Fisher z), covariates dataset(3 dummies)+sex+stage
for L, n, m in [('Hepatocytes', 70, 23208), ('Endothelial_cells', 65, 13393), ('Fibroblasts', 49, 9968), ('Cholangiocytes', 46, 10970)]:
    k = 4  # nuisance columns besides intercept and stage
    se = 1 / np.sqrt(n - k - 3)
    for a, lab in [(0.05, 'nominal 0.05'), (0.05 / m, 'Bonferroni 0.05/m (first BH rejection)')]:
        z = (norm.ppf(1 - a / 2) + zb) * se
        print(f"{L:18s} n={n} per-gene |partial r| MDE at {lab}: {np.tanh(z):.3f}")
