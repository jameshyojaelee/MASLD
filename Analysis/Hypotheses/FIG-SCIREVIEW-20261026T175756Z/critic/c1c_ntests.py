"""Critic check (1B design): the colocalized call is 'PP.H4 > 0.5 in ANY tier-1/2 study', so genes
SuSiE-tested in more studies have more chances to be called. Is the number of SuSiE-tested studies
associated with DEG status, and does stratifying on it (with AveExpr deciles) move the 1B result?"""
import numpy as np, pandas as pd
from scipy.stats import hypergeom, mannwhitneyu
R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
gl = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
deg["ens"] = deg.gene.str.split(".").str[0]; deg["DEG"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
U = deg.merge(gl, left_on="ens", right_on="ensembl")
U = U[U.susie_state_t12.isin(["positive", "tested_le_0.5"])].sort_values("ens").reset_index(drop=True)
U["gen"] = U.susie_state_t12 == "positive"
a = pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv", usecols=["gwas_name", "ensembl", "PP.H4.susie"])
a["ensembl"] = a.ensembl.astype(str).str.split(".").str[0]
t12 = set(pd.read_csv(R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv", usecols=["susie_driving_gwas"]).susie_driving_gwas.dropna())
# tier-1/2 study list from the driving and ABF driving columns is incomplete; use all studies appearing for tier12 genes with tier info if present
nt = a[a["PP.H4.susie"].notna()].groupby("ensembl").gwas_name.nunique()
U["ntest"] = U.ens.map(nt).fillna(0).astype(int)
print("studies with finite SuSiE per gene (all studies in file): median", U.ntest.median(), "range", U.ntest.min(), U.ntest.max())
print("ntest by DEG: median DEG", U[U.DEG].ntest.median(), "non-DEG", U[~U.DEG].ntest.median(), "MWU p", mannwhitneyu(U[U.DEG].ntest, U[~U.DEG].ntest).pvalue)
print("positive rate by ntest quartile:", U.groupby(pd.qcut(U.ntest.rank(method='first'), 4, labels=False)).agg(n=("gen", "size"), pos=("gen", "mean"), deg=("DEG", "mean"), med_nt=("ntest", "median")).round(3).to_dict("index"))
def test(b):
    O = int((U.gen & U.DEG).sum()); E = 0; pmf = np.array([1.0]); num = den = 0
    for k, g in U.groupby(b):
        n = len(g); G = int(g.gen.sum()); D = int(g.DEG.sum()); x = int((g.gen & g.DEG).sum())
        E += G * D / n; pmf = np.convolve(pmf, hypergeom.pmf(np.arange(min(G, D) + 1), n, G, D))
        num += x * (n - G - D + x) / n; den += (G - x) * (D - x) / n
    xs = np.arange(len(pmf)); p = pmf[np.abs(xs - E) >= abs(O - E) - 1e-9].sum()
    return O, round(E, 2), round(num / den, 3), round(p, 4)
U["dec"] = pd.qcut(U.AveExpr.rank(method="first"), 10, labels=False)
U["ntq"] = pd.qcut(U.ntest.rank(method="first"), 4, labels=False)
print("deciles only:", test("dec"))
print("ntest quartiles only:", test("ntq"))
print("deciles x ntest quartiles:", test(["dec", "ntq"]))
