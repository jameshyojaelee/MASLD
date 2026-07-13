#!/usr/bin/env python
"""
C2/C3/C5 — Does receptor localization predict semaglutide's reach?

THESIS (pre-registered): semaglutide reverses disease genes it can reach — non-parenchymal
programs + the weight-loss-reversible hepatocyte metabolic program — and leaves
hepatocyte-intrinsic, non-metabolic effectors refractory. Decision rule fixed before results:
  CLEAN   : hep-origin strongly predicts LOWER reversal net of confounds      -> unification headline
  NUANCED : reach = cell-type-of-origin AND weight-loss-metabolic pathway     -> "reach = cell-type + pathway"
  NULL    : no cell-type-of-origin signal                                     -> drop unification, use validation framing

Cell-type-of-origin metric = is_hep_dominant (argmax cell type = Hepatocytes; robust to ambient;
the continuous hep_fraction is ambient-deflated so used only as a secondary check).

C3 selectivity control: are convergence targets reversed more/same/less than matched background?
C5 shortlist: convergence Tier-1 + refractory + hepatocyte-origin + druggable.

Outputs (RNA-seq/results/glp1ra/):
  glp1ra_reach_test_summary.txt, glp1ra_reach_per_gene.csv, glp1ra_residual_target_shortlist.csv
Env: spatial (pandas/scipy/statsmodels).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
import statsmodels.formula.api as smf

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
G = os.path.join(ROOT, "RNA-seq/results/glp1ra")
AXIS = os.path.join(G, "celltype_origin_axis.csv")
PART = os.path.join(G, "glp1ra_target_reversal_partition.csv")
DRUG = os.path.join(ROOT, "data/external/drug_targets/drug_target_classification.tsv")
GMT = os.path.join(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt")

# weight-loss-reversible hepatocyte metabolic sets (the serum-GSEA axis semaglutide moves)
METAB_SETS = {"HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_ADIPOGENESIS",
              "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_XENOBIOTIC_METABOLISM",
              "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_GLYCOLYSIS",
              "HALLMARK_OXIDATIVE_PHOSPHORYLATION"}


def load_metabolic_genes():
    genes = set()
    with open(GMT) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p and p[0] in METAB_SETS:
                genes |= set(p[2:])
    return genes


# ---------------------------------------------------------------- assemble
def to_bool(s):
    return s.map(lambda x: str(x).strip().lower() in ("true", "1", "1.0"))


ax = pd.read_csv(AXIS)
ax["is_hep_dominant"] = to_bool(ax["is_hep_dominant"])
ax["origin_confident"] = to_bool(ax["origin_confident"])
# hepatocyte fold-enrichment (ambient-robust): hep vs the highest non-hep cell type
cpm_cols = [c for c in ax.columns if c.startswith("cpm__")]
hep_col = "cpm__Hepatocytes"
nonhep = [c for c in cpm_cols if c != hep_col]
ax["max_nonhep_cpm"] = ax[nonhep].max(axis=1)
ax["hep_fold_log2"] = np.log2((ax[hep_col] + 1) / (ax["max_nonhep_cpm"] + 1))
ax = ax.rename(columns={"gene": "human_symbol"})

part = pd.read_csv(PART)
for c in ["mouse_disease_sig", "conv_tier1", "is_deg_tier1"]:
    if c in part.columns:
        part[c] = to_bool(part[c])
metab = load_metabolic_genes()

d = part.merge(ax[["human_symbol", "hep_fraction", "hep_fold_log2", "tau",
                   "dominant_celltype", "is_hep_dominant", "origin_confident", "total_cpm"]],
               on="human_symbol", how="left")
d["is_hep_dominant"] = d["is_hep_dominant"].fillna(False).astype(bool)
d["origin_confident"] = d["origin_confident"].fillna(False).astype(bool)
d["in_metabolic"] = d["human_symbol"].isin(metab)
d["abs_logFC"] = d["logFC"].abs()
d["log_expr"] = np.log10(d["AveExpr"].clip(lower=0) + 1) if "AveExpr" in d else np.log10(d["total_cpm"].fillna(0) + 1)

# analysis universe = genes tested for reversal (mouse disease-significant) with a confident origin
u = d[(d["mouse_disease_sig"] == True) & (d["origin_confident"] == True) &
      d["fraction_corrected"].notna()].copy()
u["reversed_bin"] = (u["mouse_reversal"] == "reversed").astype(int)
d.to_csv(os.path.join(G, "glp1ra_reach_per_gene.csv"), index=False)

L = ["C2 — Does receptor localization predict semaglutide's reach?", "=" * 58,
     f"universe = mouse-disease-sig genes with confident cell-type origin: n={len(u)}",
     f"  hepatocyte-dominant: {int(u.is_hep_dominant.sum())} | non-parenchymal-dominant: {int((~u.is_hep_dominant).sum())}",
     f"  in weight-loss-metabolic Hallmark sets: {int(u.in_metabolic.sum())}"]

# --- C2a continuous: hepatocyte specificity vs reversal ---
for xcol, lab in [("hep_fold_log2", "hep fold-enrichment (log2)"), ("hep_fraction", "hep_fraction (ambient-caveat)")]:
    rho, p = spearmanr(u[xcol], u["fraction_corrected"], nan_policy="omit")
    L.append(f"\n[continuous] Spearman({lab}, fraction_corrected) = {rho:+.3f} (p={p:.2e})   "
             f"(negative = hepatocyte-specific genes reverse LESS)")

# --- C2b categorical: reversal by cell-type-of-origin ---
def grp(mask, name):
    s = u[mask]
    vc = s["mouse_reversal"].value_counts()
    return (f"  {name:22s} n={len(s):5d}  median frac_corrected={s.fraction_corrected.median():+.2f}  "
            f"reversed={vc.get('reversed',0)} ({100*vc.get('reversed',0)/max(len(s),1):.0f}%)  "
            f"refractory={vc.get('refractory',0)} ({100*vc.get('refractory',0)/max(len(s),1):.0f}%)")

L.append("\n[categorical] reversal by cell-type-of-origin:")
L.append(grp(u.is_hep_dominant, "hepatocyte-dominant"))
L.append(grp(~u.is_hep_dominant, "non-parenchymal"))
mw = mannwhitneyu(u.loc[u.is_hep_dominant, "fraction_corrected"],
                  u.loc[~u.is_hep_dominant, "fraction_corrected"], alternative="less")
L.append(f"  Mann-Whitney frac_corrected (hep < non-paren): p={mw.pvalue:.2e}")

# --- C2c 2x2 decomposition: cell-type x metabolic-pathway ---
L.append("\n[decomposition] median fraction_corrected by origin x weight-loss-metabolic pathway:")
dec = u.groupby([u.is_hep_dominant, u.in_metabolic])["fraction_corrected"].agg(["size", "median"])
L.append(dec.round(3).to_string())

# --- C2d confound-controlled logistic: reversed ~ hep-origin + effect size + metabolic ---
fit_df = u[u["mouse_reversal"].isin(["reversed", "refractory"])].copy()
fit_df["hep"] = fit_df["is_hep_dominant"].astype(int)
fit_df["metab"] = fit_df["in_metabolic"].astype(int)
try:
    m = smf.logit("reversed_bin ~ hep + abs_logFC + metab + log_expr", data=fit_df).fit(disp=0)
    b, se, pv = m.params["hep"], m.bse["hep"], m.pvalues["hep"]
    L.append(f"\n[logistic] P(reversed) ~ hep_origin + |logFC| + metabolic + log_expr  (n={len(fit_df)})")
    L.append(f"  hep_origin coef = {b:+.3f} (SE {se:.3f}, p={pv:.2e})  OR={np.exp(b):.2f}  "
             f"(<0 = hepatocyte-origin genes reverse LESS, net of confounds)")
    L.append(f"  metabolic coef = {m.params['metab']:+.3f} (p={m.pvalues['metab']:.2e}); "
             f"|logFC| coef = {m.params['abs_logFC']:+.3f}")
except Exception as e:
    L.append(f"\n[logistic] failed: {e}")

# --- pre-registered verdict ---
hep_med = u.loc[u.is_hep_dominant, "fraction_corrected"].median()
np_med = u.loc[~u.is_hep_dominant, "fraction_corrected"].median()
L.append("\n### VERDICT (pre-registered rule) — inspect coef sign + p + decomposition above.")
L.append(f"  hep-dominant median reversal {hep_med:.2f} vs non-parenchymal {np_med:.2f}.")

# ============================ C3 selectivity control ============================
L.append("\n\nC3 — Selectivity: are convergence targets reversed differently from matched background?")
sel = d[(d["mouse_disease_sig"] == True) & d["fraction_corrected"].notna()].copy()
sel["is_conv"] = (sel["conv_tier1"] == True).astype(int)
sel["abs_logFC"] = sel["logFC"].abs()
sel_bin = sel[sel["mouse_reversal"].isin(["reversed", "refractory"])].copy()
sel_bin["reversed_bin"] = (sel_bin["mouse_reversal"] == "reversed").astype(int)
try:
    ms = smf.logit("reversed_bin ~ is_conv + abs_logFC + log_expr", data=sel_bin).fit(disp=0)
    L.append(f"  logistic P(reversed) ~ convergence + |logFC| + log_expr (n={len(sel_bin)}): "
             f"convergence coef = {ms.params['is_conv']:+.3f} (p={ms.pvalues['is_conv']:.2e}, "
             f"OR={np.exp(ms.params['is_conv']):.2f})")
    L.append("  (coef<0 = convergence targets reverse LESS than matched disease genes = they occupy the "
             "residual space; coef~0 = generic; coef>0 = preferentially reversible/validated)")
except Exception as e:
    L.append(f"  logistic failed: {e}")
L.append(f"  raw: convergence Tier-1 median frac_corrected="
         f"{sel.loc[sel.is_conv==1,'fraction_corrected'].median():.2f} (n={int((sel.is_conv==1).sum())}) "
         f"vs non-conv disease genes {sel.loc[sel.is_conv==0,'fraction_corrected'].median():.2f} "
         f"(n={int((sel.is_conv==0).sum())})")

# ============================ C5 shortlist ============================
drug = pd.read_csv(DRUG, sep="\t").rename(columns={"symbol": "human_symbol"})
short = d.merge(drug[["human_symbol", "drug_dev_status", "pharos_tdl"]], on="human_symbol", how="left")
short = short[(short["conv_tier1"] == True) & (short["mouse_reversal"] == "refractory") &
              (short["is_hep_dominant"] == True)]
short["tractable"] = short["pharos_tdl"].isin(["Tclin", "Tchem"]) | short["druggability_tier"].notna()
short = short.sort_values(["tractable", "convergence_score"], ascending=[False, False])
short[["human_symbol", "convergence_score", "dominant_celltype", "hep_fold_log2", "fraction_corrected",
       "in_metabolic", "druggability_tier", "drug_dev_status", "pharos_tdl", "tractable"]].to_csv(
    os.path.join(G, "glp1ra_residual_target_shortlist.csv"), index=False)
L.append(f"\nC5 — residual target shortlist (conv-Tier1 + refractory + hepatocyte-origin): "
         f"{len(short)} genes; {int(short.tractable.sum())} tractable. Top: "
         f"{', '.join(short.head(12)['human_symbol'])}")

summary = "\n".join(L)
with open(os.path.join(G, "glp1ra_reach_test_summary.txt"), "w") as fh:
    fh.write(summary + "\n")
print(summary)
