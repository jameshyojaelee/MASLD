"""Size of RE-DL truncated-HK TOST at the margin boundary (true mu = 0.5), per gene, using each gene's
per-cohort variances and DL tau2 as truth. Also reports median se_nom and the per-gene smallest margin with 80% power at mu=0."""
import numpy as np, pandas as pd
from scipy import stats
exec(open("gate.py").read().split("d = pd.read_csv(COH")[0])
d = pd.read_csv(COH, sep="\t"); d["SE"] = d.logFC / d.t
Y = d.pivot(index="gene_id_versioned", columns="dataset", values="logFC")[COHORTS].values
V = d.pivot(index="gene_id_versioned", columns="dataset", values="SE")[COHORTS].values ** 2
mu, se, se_nom, tau2, *_ = re_meta(Y, V)
print("median se_nom", round(float(np.median(se_nom)), 4), "median SE_HK", round(float(np.median(se)), 4))
tcrit = stats.t.ppf(1 - 0.01772, 4)
mm = (tcrit + stats.norm.ppf(0.9)) * se_nom
print("smallest margin with ~80% power at mu=0 (ignoring HK noise), quantiles 5/25/50/75/95:", np.quantile(mm, [.05,.25,.5,.75,.95]).round(3).tolist())
rng = np.random.default_rng(20260923)
sub = rng.choice(len(Y), 3000, replace=False); B = 400
for a, tc in [("nominal 0.05", stats.t.ppf(0.95, 4)), ("BH-effective 0.0177", tcrit)]:
    rates = []
    for i in sub:
        yy = 0.5 + rng.normal(0, np.sqrt(tau2[i]), (B, 5)) + rng.normal(0, np.sqrt(V[i]), (B, 5))
        m_, s_, *_ = re_meta(yy, np.tile(V[i], (B, 1)))
        rates.append(np.mean((m_ - tc * s_ > -0.5) & (m_ + tc * s_ < 0.5)))
    rates = np.array(rates)
    print(f"true mu=0.5, alpha {a}: mean false-equivalence rate {rates.mean():.4f}; 95th pct over genes {np.quantile(rates,0.95):.4f}; max {rates.max():.4f}")
