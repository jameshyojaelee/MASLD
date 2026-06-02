#!/usr/bin/env Rscript
# sex_v3/06h_v5_to_v3_shim.R
# ---------------------------------------------------------------------------
# v5 → v3 schema shim. Rewrites sex_deg_classification_v3.csv so its
# assigned_class / assigned_class_gated / sex_class columns reflect v5
# consensus (class_v5 + evidence_tier).
# Preserves prior v3 columns under *_v3legacy and v4 columns under *_v4legacy.
# Used by downstream atlas (27a) + figures (figS_sex_dimorphism, figS_coloc_sex_direct).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })
SEXV3 <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3"
V3 <- file.path(SEXV3, "sex_deg_classification_v3.csv")
V5 <- file.path(SEXV3, "sex_deg_classification_v5.csv")

v3 <- fread(V3); v5 <- fread(V5)

# Preserve current state under v4legacy/v3legacy aliases
if (!"class_v4_legacy" %in% names(v3) && "class_v4" %in% names(v3)) {
  v3[, class_v4_legacy := class_v4]
}
v3[, assigned_class_v4legacy := assigned_class]

# Join v5 consensus
v3 <- merge(v3, v5[, .(gene, class_v5, evidence_tier,
                       stability_modal_v5 = stability_modal,
                       stability_pass_08_v5 = stability_pass_08,
                       padj_int_5k_v5 = padj_int_5k,
                       beta_int_v5 = beta_int)],
            by = "gene", all.x = TRUE)

# Map class_v5 → v3 schema vocabulary (5-class: Concordant, F_only, M_only,
# Divergent, Uncertain).
map_v5 <- function(x) fcase(
  x %in% c("Female_biased", "Female_biased_M_underpowered"), "F_only",
  x %in% c("Male_biased", "Male_biased_F_underpowered"),     "M_only",
  x %in% c("Divergent", "Divergent_one_sided"),              "divergent",
  x %in% c("Concordant", "Concordant_single_arm",
           "Sex_modifier"),                                   "concordant",
  default = "uncertain")
v3[, assigned_class := map_v5(class_v5)]
v3[, assigned_class_gated := assigned_class]
v3[, sex_class_gated := fcase(
  assigned_class_gated == "F_only",     "Female_biased",
  assigned_class_gated == "M_only",     "Male_biased",
  assigned_class_gated == "divergent",  "Divergent",
  assigned_class_gated == "concordant", "Concordant",
  assigned_class_gated == "uncertain",  "Uncertain",
  default = NA_character_)]
v3[, sex_class := sex_class_gated]

cat("=== v5 → v3 vocabulary rebroadcast ===\n")
print(table(v3$assigned_class_gated, useNA = "ifany"))
cat("\nevidence_tier:\n")
print(table(v3$evidence_tier, useNA = "ifany"))

tmp <- paste0(V3, ".tmp"); fwrite(v3, tmp); file.rename(tmp, V3)
cat("\nWrote refreshed v3 csv (v5 consensus active):", V3, "\n")
cat("  cols:", ncol(v3), "  rows:", nrow(v3), "\n")
