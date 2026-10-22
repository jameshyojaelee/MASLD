#!/usr/bin/env python3
"""Re-derive Fig 1B overlap under an exchangeable null and under expression-stratified nulls.
Inputs are the two files named in the interface input_manifest.json (sha256 checked by caller).
Universe and gene sets rebuilt exactly as scripts/manuscript/recompute_genetics_state_interface.py build()."""
import numpy as np, pandas as pd
from scipy import stats
P = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
C = P + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B = P + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
SEED = 20261015; NPERM = 10000

coloc = pd.read_csv(C)
coloc["g"] = coloc["ensembl"].astype(str).str.split(".").str[0]
coloc = coloc.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")
deg = pd.read_csv(B)
deg["g"] = deg["gene"].astype(str).str.split(".").str[0]
deg = deg.drop_duplicates("g")
jt = deg.merge(coloc[["g", "coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_state_t12"]], on="g", how="inner")
jt["gen"] = (pd.to_numeric(jt["coloc_best_susie_pp4"], errors="coerce") >= 0.5).fillna(False).to_numpy()
jt["de"] = ((jt["padj"] < 0.05) & (jt["logFC"].abs() > 0.5)).to_numpy()
N = len(jt); G = jt.gen.sum(); D = jt.de.sum(); O = (jt.gen & jt.de).sum()
tab = [[O, G - O], [D - O, N - G - D + O]]
orr, p = stats.fisher_exact(tab)
print(f"universe {N}  genetic {G}  DEG {D}  overlap {O}  expected {G*D/N:.4f}  OR {orr:.4f}  fisher p {p:.4f}")
print(f"non-overlap obs {100*(G-O)/G:.2f}%  exp {100*(1-D/N):.2f}%")
print(f"NA AveExpr {jt.AveExpr.isna().sum()}  NA SE {jt.SE.isna().sum()}")
for lab, m in [("DEG", jt.de), ("nonDEG", ~jt.de), ("genetic", jt.gen), ("nongenetic", ~jt.gen)]:
    print(f"mean AveExpr {lab}: {jt.loc[m,'AveExpr'].mean():.3f} (median {jt.loc[m,'AveExpr'].median():.3f}, n {m.sum()})")

rng = np.random.default_rng(SEED)
def strat(col, k):
    bins = pd.qcut(jt[col].rank(method="first"), k, labels=False)
    g = jt.groupby(bins).agg(n=("gen", "size"), G=("gen", "sum"), D=("de", "sum"), O=("gen", lambda s: (s & jt.loc[s.index, "de"]).sum()))
    exp = (g.G * g.D / g.n).sum()
    # hypergeometric variance per bin
    var = (g.G * g.D * (g.n - g.G) * (g.n - g.D) / (g.n**2 * (g.n - 1))).sum()
    # permutation: shuffle genetic labels within bins, count overlap
    gen = jt.gen.to_numpy(); de = jt.de.to_numpy(); b = bins.to_numpy()
    idx = [np.where(b == i)[0] for i in range(k)]
    sims = np.empty(NPERM)
    for r in range(NPERM):
        tot = 0
        for ii in idx:
            tot += (rng.permutation(gen[ii]) & de[ii]).sum()
        sims[r] = tot
    mu = sims.mean()
    p2 = (np.sum(np.abs(sims - mu) >= abs(O - mu)) + 1) / (NPERM + 1)
    pup = (np.sum(sims >= O) + 1) / (NPERM + 1); plo = (np.sum(sims <= O) + 1) / (NPERM + 1)
    # Mantel-Haenszel common OR across strata
    tabs = [np.array([[r.O, r.G - r.O], [r.D - r.O, r.n - r.G - r.D + r.O]]) for r in g.itertuples()]
    st = stats.contingency.odds_ratio  # not used
    num = sum(t[0,0]*t[1,1]/t.sum() for t in tabs); den = sum(t[0,1]*t[1,0]/t.sum() for t in tabs)
    print(f"{col} {k:>2} bins: analytic exp {exp:.2f} (sd {np.sqrt(var):.2f}); perm mean {mu:.2f}; "
          f"two-sided p {p2:.3f}; P(>=obs) {pup:.3f} P(<=obs) {plo:.3f}; MH OR {num/den:.3f}")

for k in (5, 10, 20, 50):
    strat("AveExpr", k)
for k in (10,):
    strat("SE", k)

# Logistic regression: DEG ~ genetic + natural spline of AveExpr (4 df)
import statsmodels.formula.api as smf
jt["gen_i"] = jt.gen.astype(int); jt["de_i"] = jt.de.astype(int)
m0 = smf.logit("de_i ~ gen_i", jt).fit(disp=0)
m1 = smf.logit("de_i ~ gen_i + cr(AveExpr, df=5)", jt).fit(disp=0)
for nm, m in [("unadjusted", m0), ("AveExpr spline df5", m1)]:
    b = m.params["gen_i"]; ci = m.conf_int().loc["gen_i"]
    print(f"logit {nm}: OR {np.exp(b):.3f} (95% CI {np.exp(ci[0]):.3f}-{np.exp(ci[1]):.3f}), Wald p {m.pvalues['gen_i']:.3f}")
jt[["g", "symbol", "AveExpr", "SE", "gen", "de"]].to_csv(
    "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R07-03-reproduce/interface_14704_genes.tsv", sep="\t", index=False)
