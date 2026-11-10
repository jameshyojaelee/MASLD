"""DIG1 job 55 (exploratory, post hoc to job 53). Are colocalized genes over-represented among the
hepatocyte state genes? Set: COLOC r2 primary SuSiE genes (coloc_best_susie_pp4 > 0.5), mapped to symbols
through the bulk release. State lists from job 53 (depth+injury matched): CCL20 and SOX9 concordant in both
cohorts (mean excess > 0.5 or < -0.5); Flex CXCL10 at q < 0.1 and |excess| > 0.5. Background: genes tested
in that comparison. Fisher's exact test, two-sided, descriptive (no correction family claimed).
Usage: python -I 55_coloc_state_overlap.py <new_out_dir>
"""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)

c = pd.read_csv(COLOC)
b = pd.read_csv(BULK, usecols=["gene", "symbol"])
sym = dict(zip(b.gene.str.replace(r"\..*$", "", regex=True), b.symbol))
c["symbol"] = c.ensembl.map(sym)
p = c[c.coloc_best_susie_pp4 > 0.5]
prim = set(p.symbol.dropna())
info = p.dropna(subset=["symbol"]).set_index("symbol")[["driving_trait", "coloc_best_susie_pp4"]]

D = RUN / "out_chemokine_hepatocyte_state"
X = pd.read_csv(D / "cross_cohort_excess.tsv.gz", sep="\t")
G = pd.read_csv(D / "gene_lfc_matched.tsv.gz", sep="\t")
lists = []
for a in ["CCL20", "SOX9"]:
    o = X[(X.anchor == a) & (X["mode"] == "depth+injury")].dropna(subset=["mean_excess_flex", "mean_excess_snrna"])
    bg = set(o.gene)
    lists.append((f"{a}_up_both", set(o[(o.mean_excess_flex > 0.5) & (o.mean_excess_snrna > 0.5)].gene), bg))
    lists.append((f"{a}_down_both", set(o[(o.mean_excess_flex < -0.5) & (o.mean_excess_snrna < -0.5)].gene), bg))
g = G[(G.cohort == "flex") & (G.anchor == "CXCL10") & (G["mode"] == "depth+injury") & (G.donors == "all")]
bg = set(g.gene)
lists.append(("CXCL10_flex_up_q0.1", set(g[(g.q_bh < 0.1) & (g.mean_excess > 0.5)].gene), bg))
lists.append(("CXCL10_flex_down_q0.1", set(g[(g.q_bh < 0.1) & (g.mean_excess < -0.5)].gene), bg))

rows, hits = [], []
for name, S, bg in lists:
    S = S & bg
    k, K = len(S & prim), len(bg & prim)
    OR, pv = stats.fisher_exact([[k, len(S) - k], [K - k, len(bg) - len(S) - K + k]])
    rows.append(dict(state_list=name, n_genes=len(S), n_background=len(bg), coloc_in_background=K, observed=k,
                     expected=len(S) * K / len(bg), odds_ratio=OR, p_fisher=pv))
    hits += [dict(state_list=name, gene=h, **info.loc[h].to_dict()) for h in sorted(S & prim)]
pd.DataFrame(rows).to_csv(OUT / "coloc_overlap.tsv", sep="\t", index=False)
pd.DataFrame(hits).to_csv(OUT / "coloc_overlap_genes.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
print(f"primary {len(p)} mapped {len(prim)}")
print(pd.DataFrame(rows).round(3).to_string()); print(pd.DataFrame(hits).to_string())
