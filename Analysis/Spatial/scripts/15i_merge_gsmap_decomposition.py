#!/usr/bin/env python
"""Merge the 3-trait (liver-enzyme) and 4-trait (disease-diagnosis) gsMap
convergence-decomposition Fisher-enrichment results into one canonical 7-trait file
for Fig 4g. Supersedes both source files — this becomes the sole source of truth.
No new statistical compute: both inputs already contain final Fisher OR/CI/p values.

Output: Analysis/Spatial/results/gsmap/arm_decomposition_spatial_risk.csv (overwritten,
        7 traits, 56 rows)
Removes: Analysis/Spatial/results/gsmap/arm_decomposition_disease_traits_INVESTIGATION.csv
Env: rnaseq or spatial (pandas only, no other deps)
"""
import pandas as pd
from pathlib import Path

GS = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/results/gsmap")

ENZYME_TRAITS  = {"ukbb_alt", "ukbb_ast", "ukbb_ggt"}
DISEASE_TRAITS = {"finngen_nafld", "finngen_nash", "ghodsian_nafld", "pdff"}
EXPECTED_COLS  = ["group", "in_figure", "trait", "cohort", "n_set", "a",
                  "fisher_or", "ci_low", "ci_high", "fisher_pval", "neglog10p"]

enzymes = pd.read_csv(GS / "arm_decomposition_spatial_risk.csv")
disease = pd.read_csv(GS / "arm_decomposition_disease_traits_INVESTIGATION.csv")

assert list(enzymes.columns) == EXPECTED_COLS, f"enzyme schema drift: {list(enzymes.columns)}"
assert list(disease.columns) == EXPECTED_COLS, f"disease schema drift: {list(disease.columns)}"
assert set(enzymes["trait"].unique()) == ENZYME_TRAITS, f"unexpected enzyme traits: {set(enzymes['trait'].unique())}"
assert set(disease["trait"].unique()) == DISEASE_TRAITS, f"unexpected disease traits: {set(disease['trait'].unique())}"

merged = pd.concat([enzymes, disease], ignore_index=True)

assert len(merged) == len(enzymes) + len(disease), "row count mismatch after concat"
assert set(merged["trait"].unique()) == ENZYME_TRAITS | DISEASE_TRAITS, "trait set incomplete after merge"
in_fig = merged[merged["in_figure"]]
assert len(in_fig) == 7 * 4 * 2, f"expected 56 in_figure rows, got {len(in_fig)}"

merged.to_csv(GS / "arm_decomposition_spatial_risk.csv", index=False)
print(f"[15i] wrote {len(merged)} rows ({merged['trait'].nunique()} traits) to "
      f"{GS / 'arm_decomposition_spatial_risk.csv'}")

investigation_file = GS / "arm_decomposition_disease_traits_INVESTIGATION.csv"
investigation_file.unlink()
print(f"[15i] removed superseded {investigation_file}")
