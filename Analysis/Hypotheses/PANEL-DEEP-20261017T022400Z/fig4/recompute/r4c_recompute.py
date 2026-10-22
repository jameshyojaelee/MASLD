"""Fig 4C independent recomputation from upstream H59 outputs and donor scores.

Inputs (upstream analysis outputs, not the panel source table):
  H59 results_primary_family.tsv  (per-program GLS beta, SE, CI, permutation p/q)
  H59 null_critical_abs_t.tsv     (per-program permutation critical |t| at 0.05/117)
  donor_program_scores_primary.tsv (frozen donor-level Hotspot scores)
  donor_metadata_stage_corrected_donor.tsv (2026-10-14 source-diagnosis stages)
Outputs: r4c_rows.tsv (24 rows, beta, CIs, gray interval, heatmap medians under 3 centring variants)
"""
import pandas as pd, numpy as np, os
root = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out = f"{root}/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
H = f"{root}/Analysis/Hypotheses/H59-20261015T112432Z/out"
fam = pd.read_csv(f"{H}/outcome/results_primary_family.tsv", sep="\t")
crit = pd.read_csv(f"{H}/null_critical_abs_t.tsv", sep="\t")
fam = fam.merge(crit, on="program", how="left")
assert len(fam) == 117

# ---- row choice: per non-T lineage, 3 positive + 3 negative with smallest perm q
rows = []
for ct in ["hepatocytes", "fibroblasts", "macrophages", "cholangiocytes"]:
    x = fam[(fam.cell_type == ct) & np.isfinite(fam.gls_beta)].copy()
    x["abs_t"] = x.gls_t.abs()
    for side, sub in [("positive", x[x.gls_beta > 0]), ("negative", x[x.gls_beta < 0])]:
        sub = sub.sort_values(["perm_q", "perm_p", "abs_t"], ascending=[True, True, False]).head(3)
        sub = sub.assign(side=side)
        rows.append(sub)
rows = pd.concat(rows)
rows["gray_lo"] = rows.gls_beta - rows.null_abs_t_crit * rows.gls_se
rows["gray_hi"] = rows.gls_beta + rows.null_abs_t_crit * rows.gls_se
rows["gray_excludes0"] = (rows.gray_lo > 0) | (rows.gray_hi < 0)
rows["ci_excludes0"] = (rows.gls_ci_low > 0) | (rows.gls_ci_high < 0)
fam["gray_lo"] = fam.gls_beta - fam.null_abs_t_crit * fam.gls_se
fam["gray_hi"] = fam.gls_beta + fam.null_abs_t_crit * fam.gls_se
print("all-117: gray excludes 0:", int(((fam.gray_lo > 0) | (fam.gray_hi < 0)).sum()),
      "; nominal CI excludes 0:", int(((fam.gls_ci_low > 0) | (fam.gls_ci_high < 0)).sum()),
      "; min perm_q", fam.perm_q.min(), "; n perm_q==min", int((fam.perm_q == fam.perm_q.min()).sum()),
      "; selected", int(fam.selected_H59.sum()))

# ---- heatmap: median dataset-centred donor score by stage in GSE202379 + GSE244832
S = pd.read_csv(f"{root}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv", sep="\t")
M = pd.read_csv(f"{root}/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv", sep="\t")
S = S.merge(M[["sample", "disease_stage_coarse", "exclude_stage_analysis"]].rename(
    columns={"sample": "donor", "disease_stage_coarse": "stage_corr", "exclude_stage_analysis": "excl_corr"}), on="donor", how="left")
S = S[S.dataset.isin(["GSE202379", "GSE244832"]) & S.stage_corr.isin(["Healthy", "Steatosis", "Steatohepatitis"]) & (~S.excl_corr.astype(bool))]
g = S.groupby(["cell_type", "module", "dataset"]).score
S["c_mean"] = S.score - g.transform("mean")
S["c_median"] = S.score - g.transform("median")
S["z"] = (S.score - g.transform("mean")) / g.transform("std")
heat = S.groupby(["cell_type", "module", "stage_corr"]).agg(
    med_c_mean=("c_mean", "median"), med_c_median=("c_median", "median"), med_z=("z", "median"),
    mean_c_mean=("c_mean", "mean"), n=("donor", "nunique")).reset_index()
hw = heat.pivot_table(index=["cell_type", "module"], columns="stage_corr",
                      values=["med_c_mean", "med_c_median", "med_z", "mean_c_mean", "n"])
hw.columns = [f"{a}_{b}" for a, b in hw.columns]
hw = hw.reset_index()
rows = rows.merge(hw, on=["cell_type", "module"], how="left")
cols = ["cell_type", "module", "module_name", "side", "n_fit", "gls_beta", "gls_se", "gls_ci_low", "gls_ci_high",
        "null_abs_t_crit", "gray_lo", "gray_hi", "gray_excludes0", "ci_excludes0", "perm_p", "perm_q", "selected_H59",
        "stability_median"] + [c for c in hw.columns if c not in ("cell_type", "module")]
rows[cols].to_csv(f"{out}/r4c_rows.tsv", sep="\t", index=False)
fam.to_csv(f"{out}/r4c_all117_with_gray.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(rows[["cell_type", "module", "module_name", "side", "n_fit", "gls_beta", "gls_ci_low", "gls_ci_high", "gray_lo", "gray_hi",
            "perm_q", "ci_excludes0"]].to_string())
print(rows[["cell_type", "module", "med_c_mean_Healthy", "med_c_mean_Steatosis", "med_c_mean_Steatohepatitis",
            "n_Healthy", "n_Steatosis", "n_Steatohepatitis"]].to_string())
# compare to H59's own E9 row choice (an upstream H59 output)
e9 = pd.read_csv(f"{H}/outcome/E9_fig4c_row_choice.tsv", sep="\t")
mm = rows.merge(e9, on=["cell_type", "module"], how="outer", indicator=True)
print("row choice vs H59 E9:", mm._merge.value_counts().to_dict())
