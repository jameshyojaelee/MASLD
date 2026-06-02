#!/usr/bin/env Rscript
# sex_v3/06i_v6_to_v3_shim.R
# ---------------------------------------------------------------------------
# v6 → v3 schema shim. Rewrites sex_deg_classification_v3.csv so its
# assigned_class / assigned_class_gated / sex_class columns reflect the v6
# cross-pillar consensus (class_v6_consensus + class_v6_combined). Preserves
# the v5-derived state (set by 06h) under *_v5legacy columns so downstream
# audit can still recover the fixed-effect candidate list.
#
# Why this script exists (meta-review 2026-05-16, P1-14):
#   `09_atlas_refresh.R` + `27a_assemble_evidence_atlas.R` + figure data
#   loader `load_figure_data.R::load_sex_classification` all read
#   `sex_deg_classification_v3.csv` and consume `assigned_class*`/`sex_class*`
#   columns. Without this shim, the v6 cross-pillar consensus never reaches
#   the atlas or any downstream figure script. 06h emitted the v5 vocabulary;
#   06i overlays v6.
#
# Inputs:
#   sex_deg_classification_v3.csv  (current v3 csv — after 06h v5 shim)
#   sex_deg_classification_v6.csv  (cross-pillar consensus from 18_consensus_v6)
# Output:
#   sex_deg_classification_v3.csv  (rewritten in place, atomic; v6 active)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })

SEXV3 <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3"

# Contrast routing (sex_v3_utils.R::contrast_paths) — overrides SEXV3 / IDIR
# when CONTRAST_NAME != "disease_vs_ctrl"; default preserves legacy layout.
if (!exists("contrast_paths")) source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3/sex_v3_utils.R"))
.cpaths <- contrast_paths()
SEXV3 <- .cpaths$sexv3
IDIR  <- .cpaths$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
V3 <- file.path(SEXV3, "sex_deg_classification_v3.csv")
V6 <- file.path(SEXV3, "sex_deg_classification_v6.csv")

stopifnot(file.exists(V3), file.exists(V6))
v3 <- fread(V3); v6 <- fread(V6)

# Snapshot the v5-derived state (set by 06h) under *_v5legacy aliases so the
# atlas can still recover the fixed-effect candidate list for audit. Idempotent.
if (!"assigned_class_v5legacy" %in% names(v3)) {
  v3[, assigned_class_v5legacy   := assigned_class]
  v3[, sex_class_v5legacy        := sex_class]
}

# Join v6 cross-pillar consensus columns.
v6_cols <- intersect(
  c("gene", "class_v6_consensus", "class_v6_combined",
    "q_emp_per_gene", "padj_int_5k_random", "stability_v6",
    "loco_rho", "I2_pct", "n_pillars_used", "consensus_quality"),
  names(v6))
v3 <- merge(v3, v6[, ..v6_cols], by = "gene", all.x = TRUE)

# Map v6 class_v6_combined → v3 5-class vocabulary
# class_v6_combined examples: "Suggestive_F_biased", "Suggestive_M_biased",
# "Suggestive_Divergent", "Suggestive_NA", "Strong_F_biased",
# "Power_limited", "Uncertain".
map_v6 <- function(x) fcase(
  is.na(x) | x == "Uncertain" | x == "Power_limited", "uncertain",
  grepl("_F_biased$",  x), "F_only",
  grepl("_M_biased$",  x), "M_only",
  grepl("_Divergent$", x), "divergent",
  grepl("_NA$",        x), "uncertain",  # dir-NA = no direction call
  grepl("_Concordant$", x), "concordant",
  default = "uncertain"
)
v3[, assigned_class := map_v6(class_v6_combined)]
v3[, assigned_class_gated := assigned_class]
v3[, sex_class_gated := fcase(
  assigned_class_gated == "F_only",     "Female_biased",
  assigned_class_gated == "M_only",     "Male_biased",
  assigned_class_gated == "divergent",  "Divergent",
  assigned_class_gated == "concordant", "Concordant",
  assigned_class_gated == "uncertain",  "Uncertain",
  default = NA_character_)]
v3[, sex_class := sex_class_gated]

# Evidence-tier column — promote v6 consensus quality so figures and atlas
# can filter on it directly.
v3[, evidence_tier_v6 := class_v6_consensus]

cat("=== v6 → v3 vocabulary rebroadcast ===\n")
cat("assigned_class_gated:\n")
print(table(v3$assigned_class_gated, useNA = "ifany"))
cat("\nclass_v6_consensus:\n")
print(table(v3$class_v6_consensus, useNA = "ifany"))
cat("\nclass_v6_combined:\n")
print(table(v3$class_v6_combined, useNA = "ifany"))
cat("\nconsensus_quality:\n")
print(table(v3$consensus_quality, useNA = "ifany"))

# Atomic write
tmp <- paste0(V3, ".tmp.", Sys.getpid())
fwrite(v3, tmp)
file.rename(tmp, V3)
cat("\nWrote refreshed v3 csv (v6 consensus active):", V3, "\n")
cat("  cols:", ncol(v3), "  rows:", nrow(v3), "\n")
