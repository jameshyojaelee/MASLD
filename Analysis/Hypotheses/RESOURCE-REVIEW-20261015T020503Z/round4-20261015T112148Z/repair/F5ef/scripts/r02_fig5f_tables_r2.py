#!/usr/bin/env python3
"""F5ef repair: relabel the corrected Fig 5F source tables (no statistic is refitted).

Inputs (read only): the first-round corrected tables in figures/candidates/.../F5ef/.
Outputs (new files): figures/candidates/.../F5ef/repair/*_r2.tsv.

Changes:
  * name the sensitivity covariate: the lineage-only design is the cell2location Hepatocytes q05
    abundance plus log1p library size and log1p detected genes for both programs
    (jobs/1731/rerun_spatial.py line 8 and 150), so lineage_only_* becomes hepatocyte_only_*;
  * the matched-gene null's lineage correlation is correlation with the Hepatocytes abundance;
  * the hepatocyte-only P at k=0 is the 1e-4 floor, written as <=1e-4;
  * the overlay section P (k=0) is removed because no section-level P is reported (RESULTS);
  * unit sign-flip values and the between-program contrast are labelled post hoc descriptive:
    no prespecified test compared the two programs.
Every numeric value is checked against the adopted release before writing.
"""
from pathlib import Path

import pandas as pd

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CAND = P / "figures/candidates/review-r3-fixes-20261015T112148Z/F5ef"
OUT = CAND / "repair"
REL = P / "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2"

t = pd.read_csv(CAND / "fig5f_spatial_program_calibration_corrected.tsv", sep="\t")
b = pd.read_csv(CAND / "fig5f_between_program_unit_contrast.tsv", sep="\t")
rp = pd.read_csv(REL / "full_composition_program_results.tsv", sep="\t")

# guard: the plotted values equal the adopted release
m = t.merge(rp[["program_id", "dataset", "residual_moran_i", "residual_null_mean", "residual_null_sd", "residual_padj"]],
            left_on=["program_uid", "source_dataset_id"], right_on=["program_id", "dataset"], suffixes=("", "_rel"))
assert len(m) == 4
for c in ("residual_moran_i", "residual_null_mean", "residual_null_sd"):
    assert (m[c] - m[c + "_rel"]).abs().max() < 1e-6, c
assert (m.bh_q - m.residual_padj).abs().max() < 1e-9
assert ((m.moran_excess - (m.residual_moran_i_rel - m.residual_null_mean_rel)).abs() < 1e-6).all()

t = t.rename(columns={"lineage_only_residual_moran_i": "hepatocyte_only_residual_moran_i",
                      "lineage_only_excess": "hepatocyte_only_excess",
                      "lineage_only_k_null_ge_obs": "hepatocyte_only_k_null_ge_obs"})
t["hepatocyte_only_adjustment"] = ("sensitivity: cell2location Hepatocytes q05 abundance (discovery lineage of both "
                                   "programs), log1p library size and log1p detected genes")
t["hepatocyte_only_p"] = [("<=1e-4 (k=0; resolution floor of 9,999 draws)" if k == 0 else f"{(k + 1) / 10000:.4f}")
                          for k in t.hepatocyte_only_k_null_ge_obs]
t["null_resamples"] = ("genes matched on expression, detection, biotype and correlation with the cell2location "
                       "Hepatocytes abundance; donors and arrays held fixed")
t["unit_signflip_status"] = "post hoc, descriptive; not prespecified; floor 1/2^n_units"
t["overlay_section_rank"] = ["largest excess of 5 GSE192741 sections" if isinstance(s, str) else ""
                             for s in t.overlay_section]
t = t.drop(columns=["overlay_section_p"])
t["between_program_test"] = "none prespecified"
t.to_csv(OUT / "fig5f_spatial_program_calibration_corrected_r2.tsv", sep="\t", index=False)

b["note"] = ["post hoc descriptive; no prespecified test compared the two programs; donors",
             "post hoc descriptive; no prespecified test compared the two programs; arrays are not donors"]
b.to_csv(OUT / "fig5f_between_program_unit_contrast_r2.tsv", sep="\t", index=False)
print(t[["display_name", "source_dataset_id", "evidence_state", "moran_excess", "residual_null_sd", "bh_q",
         "hepatocyte_only_excess", "hepatocyte_only_p"]].to_string(index=False))
print(b.to_string(index=False))
