"""Audit check (design only, no expression data):
1) D in the 537-unit DR-main cell is uncorrelated with the stored 774-unit estimate b_774 (so selecting calls on b_774
   does not change the distribution of D_537), and D in a 774-unit cell is uncorrelated with b_537.
2) Power of the conditional per-gene test t = s(D + 0.5 b_u)/sqrt(s2(u_a^2 - u_u^2)) over 8 cells, BH over 1,854,
   for a common ratio r of adjusted to observed unadjusted effect (normal approximation, same as plan power_keep.py)."""
import numpy as np, pandas as pd
from scipy.linalg import helmert
from scipy.stats import norm
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P = ROOT + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
dz = pd.read_csv(P + "model_design.tsv", sep="\t").set_index("sample_id")
w = pd.read_csv(P + "lvqw_sample_weights.tsv", sep="\t").set_index("sample_id").sample_weight.loc[dz.index].to_numpy()
man = pd.read_csv(ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/sample_manifest.tsv", sep="\t").set_index("sample_id").loc[dz.index]
new = pd.read_csv(ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz", sep="\t").set_index("sample_id")
ct = [c for c in new.columns if c != "dataset"]; F = new.loc[dz.index, ct]
DR5 = {"hep_epi": ["Hepatocytes", "Cholangiocytes"], "stromal": ["Fibroblasts", "Endothelial cells"],
       "mac_mono": ["Macrophages", "Mono+mono derived cells"], "oth_my": ["cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"],
       "lymph": ["T cells", "Circulating NK/NKT", "Resident NK", "B cells", "Plasma cells"]}
def ilr(x, floor):
    x = np.maximum(np.asarray(x, float), floor); return np.log(x) @ (-helmert(x.shape[1]).T)
main = man.dataset.isin(["GSE130970", "GSE162694", "GSE213621"]).to_numpy()
Xu = dz.to_numpy(float)
def hrow(X, wt):
    sw = np.sqrt(wt); Xw = X * sw[:, None]; XtXi = np.linalg.inv(Xw.T @ Xw); return ((XtXi @ Xw.T) * sw[None, :])[-1], XtXi[-1, -1]
# 774 cells
Za = ilr(F.to_numpy(), 1e-6)
h_u774, uu774 = hrow(Xu, w); h_a774, ua774 = hrow(np.hstack([Xu[:, :-1], Za, Xu[:, -1:]]), w)
# DR-main cell (537 units; drop empty GSE135251 dummy and GSE130970 dummy as the baseline)
a = np.column_stack([F[v].sum(1) for v in DR5.values()]); a = a / a.sum(1, keepdims=True); Zd = ilr(a, 1e-6)
cols = [c for c in dz.columns if c not in ("datasetGSE135251", "datasetGSE130970", "group_binaryDisease")] + ["group_binaryDisease"]
X537 = dz.loc[main, cols].to_numpy(float)
h_u537s, uu537 = hrow(X537, w[main]); h_a537s, ua537 = hrow(np.hstack([X537[:, :-1], Zd[main], X537[:, -1:]]), w[main])
h_u537 = np.zeros(len(w)); h_u537[main] = h_u537s; h_a537 = np.zeros(len(w)); h_a537[main] = h_a537s
Winv = 1 / w
c = lambda x, y: (x * y * Winv).sum()
print("design columns DR-main:", cols)
print(f"corr(D_537, b_774) = {c(h_a537 - h_u537, h_u774) / np.sqrt(c(h_a537 - h_u537, h_a537 - h_u537) * c(h_u774, h_u774)):.2e}")
print(f"corr(D_537, b_537) = {c(h_a537 - h_u537, h_u537) / np.sqrt(c(h_a537 - h_u537, h_a537 - h_u537) * c(h_u537, h_u537)):.2e}")
print(f"corr(D_774(A16 rebuilt), b_537) = {c(h_a774 - h_u774, h_u537) / np.sqrt(c(h_a774 - h_u774, h_a774 - h_u774) * c(h_u537, h_u537)):.2e}")
print(f"corr(D_774(A16 rebuilt), b_774) = {c(h_a774 - h_u774, h_u774) / np.sqrt(c(h_a774 - h_u774, h_a774 - h_u774) * c(h_u774, h_u774)):.2e}")
# power of the conditional test (cells: inflation k, unadjusted SE ratio s vs 774), from n_and_design.out / dr_main_check.out
deg = pd.read_csv(P + "deg_results.csv"); cl = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)]
b = cl.logFC.abs().to_numpy(); se = cl.SE.to_numpy()
cells = [(1.157, 1), (1.087, 1), (1.148, 1), (1.116, 1.049), (1.177, 1), (1.088, 1), (1.140, 1), (1.135, 1.049)]
print("\nmedian SE of D per cell (median call SE x s x sqrt(k^2-1)):", [round(float(np.median(se)) * s * np.sqrt(k**2 - 1), 4) for k, s in cells])
def bh(p):
    n = len(p); o = np.argsort(p); q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)
for r in (1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2):
    pk, pl = [], []
    for k, s in cells:
        sed = se * s * np.sqrt(k**2 - 1); mu = (r - 0.5) * b
        pk.append(norm.sf(mu / sed)); pl.append(norm.cdf(mu / sed))
    pk = np.max(pk, 0); pl = np.max(pl, 0)
    print(f"ratio {r}: expected-z T_keep passes {int((bh(pk) < 0.05).sum())} of 1854; T_lose passes {int((bh(pl) < 0.05).sum())}")
