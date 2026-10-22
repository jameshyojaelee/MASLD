"""Independent recomputation of load-bearing numbers in title/abstract/intro/fig1/fig2.
Read-only on all inputs. Prints key=value lines. Run with PYTHONNOUSERSITE=1."""
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable

R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
C = R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
BD = R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/"
B = BD + "primary_source_controls/deg_results.csv"
CAT = R + "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/"


def out(k, v):
    print(f"{k}\t{v}")


# ---------------- bulk
b = pd.read_csv(B)
deg = (b.padj < 0.05) & (b.logFC.abs() > 0.5)
out("bulk_genes_tested", len(b))
out("bulk_canonical", int(deg.sum()))
out("bulk_canonical_up", int((deg & (b.logFC > 0)).sum()))
out("bulk_canonical_down", int((deg & (b.logFC < 0)).sum()))
cg = pd.read_csv(BD + "primary_source_controls/cohort_group_counts.tsv", sep="\t")
out("bulk_cohorts", len(cg))
out("bulk_controls", int(cg.new_control.sum()))
out("bulk_disease", int(cg.new_disease.sum()))
roster = pd.read_csv(BD + "primary_unit_roster.tsv", sep="\t")
out("bulk_roster_units", roster.analysis_unit_id.nunique())

# ---------------- COLOC gene level
c = pd.read_csv(C)
s = pd.to_numeric(c.coloc_best_susie_pp4, errors="coerce")
a = pd.to_numeric(c.coloc_best_abf_pp4, errors="coerce")
susie_pos = s > 0.5
abf_pos = a > 0.5
out("coloc_rows", len(c))
out("susie_pos", int(susie_pos.sum()))
out("susie_state_positive", int((c.susie_state_t12 == "positive").sum()))
out("abf_pos", int(abf_pos.sum()))
out("union", int((susie_pos | abf_pos).sum()))
out("abf_only", int((abf_pos & ~susie_pos).sum()))
unt = c.susie_state_t12 == "untestable"
out("susie_untestable", int(unt.sum()))
out("untestable_in_abf_only", int((unt & abf_pos & ~susie_pos).sum()))
pos = c[susie_pos]
out("susie_pos_driver_trait", pos.driving_trait.value_counts().to_dict())
enz = pos.driving_trait.isin(["ALT", "AST", "GGT"]).sum()
out("susie_pos_enzyme_driver", int(enz))
out("susie_pos_direct_driver", int(len(pos) - enz))
out("susie_pos_driver_tier", pos.driving_tier.value_counts().to_dict())

# ---------------- interface universes (H57 construction, re-implemented)
b["g"] = b.gene.astype(str).str.split(".").str[0]
c["g"] = c.ensembl.astype(str).str.split(".").str[0]
cc = c.assign(_s=s, _a=a).sort_values(["_s", "_a"], ascending=False).drop_duplicates("g")
jt = b.drop_duplicates("g").merge(cc[["g", "_s", "_a", "susie_state_t12"]], on="g", how="inner")
jt["gen"] = jt._s > 0.5
jt["deg"] = (jt.padj < 0.05) & (jt.logFC.abs() > 0.5)
out("iface_N", len(jt))
out("iface_gen", int(jt.gen.sum()))
out("iface_deg", int(jt.deg.sum()))
out("iface_both", int((jt.gen & jt.deg).sum()))
out("iface_state_x_deg", pd.crosstab(jt.susie_state_t12.fillna("NA"), jt.deg).to_dict())
ue = jt[jt.susie_state_t12.isin(["positive", "tested_le_0.5"])].sort_values("g").reset_index(drop=True)
out("ueval_N", len(ue))
out("ueval_gen", int(ue.gen.sum()))
out("ueval_deg", int(ue.deg.sum()))
out("ueval_both", int((ue.gen & ue.deg).sum()))
nores = jt[jt.susie_state_t12 == "no_susie_result"]
out("no_multisignal_result_N", len(nores))
out("no_multisignal_result_deg", int(nores.deg.sum()))
out("iface_untestable", int((jt.susie_state_t12 == "untestable").sum()))
stateonly = jt[jt.deg & ~jt.gen]
out("state_only", len(stateonly))
out("state_only_no_result", int((stateonly.susie_state_t12 == "no_susie_result").sum()))

# stratified test, 10 AveExpr quantile bins with ties broken by Ensembl order
bins = pd.qcut(ue.AveExpr.rank(method="first"), 10, labels=False).to_numpy()
gen = ue.gen.to_numpy(); dg = ue.deg.to_numpy()
tabs, E, pmf = [], 0.0, np.array([1.0])
for k in np.unique(bins):
    m = bins == k
    n, G, D, O = m.sum(), gen[m].sum(), dg[m].sum(), (gen & dg)[m].sum()
    E += G * D / n
    tabs.append(np.array([[O, G - O], [D - O, n - G - D + O]], dtype=float))
    x = np.arange(0, min(G, D) + 1)
    pmf = np.convolve(pmf, hypergeom.pmf(x, n, D, G))
O = int((gen & dg).sum())
st = StratifiedTable(tabs)
lo, hi = st.oddsratio_pooled_confint()
xx = np.arange(len(pmf))
p = pmf[np.abs(xx - E) >= abs(O - E) - 1e-7].sum()
out("ueval_E_strat10", round(E, 4))
out("ueval_MH_OR", round(st.oddsratio_pooled, 4))
out("ueval_MH_CI", f"{lo:.4f}-{hi:.4f}")
out("ueval_exact_p", round(p, 4))
out("ueval_E_unstrat", round(ue.gen.sum() * ue.deg.sum() / len(ue), 3))

# ---------------- Catalog r8 / r9d
gi = pd.read_parquet(CAT + "passport_gene_index.parquet")
out("catalog_genes", len(gi))
pm = pd.read_parquet(CAT + "passport_program_membership.parquet")
out("catalog_membership_rows", len(pm))
out("catalog_membership_cols", list(pm.columns))
gcol = [x for x in pm.columns if "gene" in x.lower()]
for gc in gcol:
    out(f"catalog_membership_unique_{gc}", pm[gc].nunique())
pcol = [x for x in pm.columns if "program" in x.lower()]
for pc in pcol:
    out(f"catalog_membership_unique_{pc}", pm[pc].nunique())
