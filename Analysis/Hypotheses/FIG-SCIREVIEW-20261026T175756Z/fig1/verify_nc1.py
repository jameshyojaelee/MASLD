"""Verifier: rebuild the expression-matched non-colocalized control (greedy nearest AveExpr, no replacement) in 3 visiting orders;
list the 10 colocalized genes outside the bulk table."""
import numpy as np, pandas as pd
from scipy.stats import hypergeom
R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
gl = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
deg["ens"] = deg.gene.str.split(".").str[0]; deg["DEG"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
U = deg.merge(gl, left_on="ens", right_on="ensembl")
U = U[U.susie_state_t12.isin(["positive", "tested_le_0.5"])].sort_values("ens").reset_index(drop=True)
U["gen"] = U.susie_state_t12 == "positive"
U["b"] = pd.qcut(U.AveExpr.rank(method="first"), 10, labels=False)
nk = U.groupby("b").size(); Dk = U.groupby("b").DEG.sum()
cand = U[~U.gen]; ca = cand.AveExpr.to_numpy(); cidx = cand.index.to_numpy()
for seed in (1, 2, 3):
    rng = np.random.default_rng(seed); used = np.zeros(len(ca), bool); pick = []
    for i in rng.permutation(U.index[U.gen]):
        d = np.abs(ca - U.at[i, "AveExpr"]); d[used] = np.inf
        m = d.min(); ties = np.where(d == m)[0]; j = ties[np.argmin(ca[ties])]
        used[j] = True; pick.append(cidx[j])
    L = U.loc[pick]; O = int(L.DEG.sum())
    gk = L.groupby("b").size().reindex(range(10), fill_value=0)
    E = float((gk * Dk / nk).sum())
    pmf = np.array([1.0])
    for k in range(10):
        pmf = np.convolve(pmf, hypergeom.pmf(np.arange(min(gk[k], Dk[k]) + 1), nk[k], gk[k], Dk[k]))
    x = np.arange(len(pmf)); p = pmf[np.abs(x - E) >= abs(O - E) - 1e-9].sum()
    print(f"order seed {seed}: O {O}, E_strat {E:.2f}, O/E {O/E:.4f}, exact p {p:.4f}")
pos = gl[gl.susie_state_t12 == "positive"]
out = pos[~pos.ensembl.isin(deg.ens)]
a = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc.csv", usecols=["gene", "ensembl"])
a["ensembl"] = a.ensembl.astype(str).str.split(".").str[0]
print("positives not in bulk:", len(out), sorted(a[a.ensembl.isin(out.ensembl)].gene.tolist()))
ev = gl[gl.susie_state_t12.isin(["positive", "tested_le_0.5"])]
print("evaluated t12:", len(ev), "not in bulk:", int((~ev.ensembl.isin(deg.ens)).sum()))
