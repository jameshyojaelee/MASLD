"""Design-only power for the per-gene half-retained (T_keep) and half-removed (T_lose) tests.
Uses stored unadjusted logFC and SE of the 1,854 calls and the design SE inflation per cell; normal approximation;
assumes every gene has the same true ratio r in every cell. No adjusted outcome is used."""
import numpy as np, pandas as pd
from scipy.stats import norm
P = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
deg = pd.read_csv(P + "deg_results.csv"); c = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)]
b = c.logFC.abs().to_numpy(); se = c.SE.to_numpy()
# (inflation, unadjusted SE ratio vs 774) per cell: A16, C7, D1, DRmain x rebuilt, earlier (n_and_design.out, dr_main_check.out)
cells = [(1.157, 1), (1.087, 1), (1.148, 1), (1.116, 1.049), (1.177, 1), (1.088, 1), (1.140, 1), (1.135, 1.049)]
def bh(p):
    n = len(p); o = np.argsort(p); q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)
for r in (1.0, 0.9, 0.8, 0.7, 0.3, 0.2):
    pk, pl = [], []
    for k, s in cells:
        sek = se * s * np.sqrt(k**2 - 0.75)       # SE of (b_adj - 0.5 b_unadj)
        mu = (r - 0.5) * b                        # its expected value, oriented to the call's sign
        pk.append(norm.sf(mu / sek)); pl.append(norm.cdf(mu / sek))
    pk = np.max(pk, 0); pl = np.max(pl, 0)       # IUT: largest p over the 8 cells (expected-z version)
    print(f"true ratio {r}: expected-z T_keep passes {int((bh(pk) < 0.05).sum())} of 1854; T_lose passes {int((bh(pl) < 0.05).sum())}")
