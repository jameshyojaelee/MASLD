"""F1 reproduction: Fig 1B universes, counts, unstratified and AveExpr-stratified expectation,
exact Fisher 95% CI on the 14,704 universe (R12-06), and the per-gene source table for the
corrected Fig 1B. Independent of the H57 script except for the stated universe rule.
Inputs: r2 gene-level COLOC tier-1/2 table and the adopted 774-unit refit DEG table."""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WS = os.path.dirname(os.path.abspath(__file__))
FIG = ROOT + "/figures/candidates/review-r3-fixes-20261015T112148Z/F1"
C = ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B = ROOT + "/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"

coloc = pd.read_csv(C)
coloc["g"] = coloc["ensembl"].astype(str).str.split(".").str[0]
coloc = coloc.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")
bulk = pd.read_csv(B)
bulk["g"] = bulk["gene"].astype(str).str.split(".").str[0]
full = bulk.drop_duplicates("g").merge(coloc[["g", "coloc_best_susie_pp4", "susie_state_t12"]], on="g", how="inner")
full["susie"] = pd.to_numeric(full.coloc_best_susie_pp4, errors="coerce")
full["gen"] = (full.susie >= 0.5).fillna(False).astype(bool)
full["deg"] = ((full.padj < 0.05) & (full.logFC.abs() > 0.5)).astype(bool)
full = full.sort_values("g").reset_index(drop=True)


def counts(df, name):
    N = len(df); G = int(df.gen.sum()); D = int(df.deg.sum()); O = int((df.gen & df.deg).sum())
    t = [[O, G - O], [D - O, N - G - D + O]]
    orr, p = stats.fisher_exact(t)
    ci = stats.contingency.odds_ratio(t, kind="conditional").confidence_interval(0.95)
    E = G * D / N
    print(f"{name}: N {N} G {G} D {D} O {O} E {E:.4f} sample OR {orr:.4f} cond-MLE CI {ci.low:.3f}-{ci.high:.3f} "
          f"Fisher p {p:.4f}; markers gen {E/G:.4%} deg {E/D:.4%}; observed {O/G:.4%} {O/D:.4%}; "
          f"non-overlap obs {(G-O)/G:.4%} exp {(G-E)/G:.4%}")
    return N, G, D, O


print(pd.crosstab(full.susie_state_t12.fillna("NA"), full.deg))
counts(full, "U_full")
ueval = full[full.susie_state_t12.isin(["positive", "tested_le_0.5"])].reset_index(drop=True)
counts(ueval, "U_eval")
ueval7 = full[full.susie.notna()].reset_index(drop=True)
counts(ueval7, "U_eval7")

# 10 AveExpr quantile bins, ties broken by Ensembl order (H57 "Bins")
bins = pd.qcut(ueval.AveExpr.rank(method="first"), 10, labels=False).to_numpy()
E = 0.0
for k in range(10):
    m = bins == k
    E += ueval.gen[m].sum() * ueval.deg[m].sum() / m.sum()
print(f"U_eval E_strat (10 AveExpr bins) {E:.4f}; markers gen {E/426:.4%} deg {E/220:.4%}; O/E {26/E:.4f}")

never = int(((full.susie_state_t12 == "no_susie_result") & full.deg).sum())
neither = int((~ueval.gen & ~ueval.deg).sum())
print(f"never-evaluated DEGs {never}; U_eval neither {neither}; U_eval7 neither {int((~ueval7.gen & ~ueval7.deg).sum())}")

out = ueval[["g", "AveExpr", "logFC", "padj", "susie", "susie_state_t12", "gen", "deg"]].copy()
out["aveexpr_bin10"] = bins
out.to_csv(FIG + "/fig1b_U_eval_gene_table.tsv", sep="\t", index=False)
