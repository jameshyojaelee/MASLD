#!/usr/bin/env Rscript
# 217d_hormone_panel_annotation_clean.R -- expanded panel + dual-role flag
#
# Agent: B5b (C5 critique fix)
#
# CHANGE vs 217d: Add `dual_role_flag` column to identify panel members that
# are ALSO targets of panel TFs (avoid double-counting in downstream analyses).
# Specifically CYP3A4 (panel + AR/ESR1/HNF4A target), SULT2A1 (panel +
# AR/HNF4A target), and CYP19A1 (panel + ESR/AR target).
#
# Inputs:
#   - sex_deg_classification.csv
#   - gene_level_coloc.csv
#   - multi_evidence_atlas.csv
#   - RNA-seq/results/stratified_causal/hormone_tf_targetsets_long_clean.csv
#
# Output:
#   - RNA-seq/results/stratified_causal/sex_hormone_panel_expanded_clean.csv

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ORIGINAL_A2 <- c(
  "AR", "ESR1", "ESR2", "PGR", "GPER1",
  "CYP19A1", "SRD5A2", "CYP17A1", "STAR",
  "HSD17B1", "HSD17B2", "HSD17B3", "HSD17B6", "HSD17B7", "HSD17B12",
  "HSD3B1",
  "SHBG", "GHR",
  "KDM6A", "DDX3X", "XIST",
  "DDX3Y", "EIF1AY", "RPS4Y1",
  "BCL6", "STAT5A", "STAT5B", "HNF4A", "CUX2"
)
C2_ADDITIONS <- c("CYP3A4", "CYP1A2", "PRLR", "SULT2A1")
PANEL <- unique(c(ORIGINAL_A2, C2_ADDITIONS))

# Dual-role: panel members that are also TF targets of panel TFs.
# These genes are receptors/biosynthesis that are themselves transcriptionally
# regulated by panel TFs -- double-counted if not flagged.
# (C5 callout: CYP3A4 and SULT2A1; we also flag CYP19A1 since it's an ESR/AR
# pathway transcript.)
DUAL_ROLE_CANDIDATES <- c("CYP3A4", "SULT2A1", "CYP19A1")

cat("=== 217d_clean: Expanded sex-hormone panel + dual_role_flag ===\n")
cat("  N panel genes:", length(PANEL), "\n")

# Load TF-target edges (clean version)
tf_long_path <- file.path(outdir, "hormone_tf_targetsets_long_clean.csv")
panel_tf_targets <- character(0)
if (file.exists(tf_long_path)) {
  tf_long <- fread(tf_long_path)
  # Only count edges sourced from non-self_loop (otherwise every TF would self-flag)
  panel_tf_targets <- unique(tf_long[source != "self_loop", target_gene])
  cat("  Distinct TF-targeted symbols (clean, excl self-loop):",
      length(panel_tf_targets), "\n")
}

# Load sex DEG classification (v3 mashr Bayesian preferred; v2 fallback)
sex_v3_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_file <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
cat("  Source:", basename(dirname(sex_file)), "/", basename(sex_file), "\n")
sex_degs <- fread(sex_file)
# v3 column compatibility shim: v3 names interaction logFC as `beta_interaction`
if (!"interaction_logFC" %in% names(sex_degs) && "beta_interaction" %in% names(sex_degs)) {
  sex_degs[, interaction_logFC := beta_interaction]
}
sex_degs[, ensembl_id := sub("\\..*", "", gene)]

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
sex_degs <- merge(sex_degs, atlas, by = "ensembl_id", all.x = TRUE)

# Load COLOC
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
if ("coloc_best_susie_pp4" %in% names(coloc)) {
  coloc[, coloc_best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                                    coloc_best_susie_pp4, coloc_best_pp4)]
}
coloc <- coloc[gene != "" & !is.na(gene)]
coloc_lookup <- unique(coloc[, .(human_symbol = gene, coloc_best_pp4, coloc_best_gwas)])
coloc_lookup <- coloc_lookup[, .(coloc_best_pp4 = max(coloc_best_pp4, na.rm = TRUE),
                                   coloc_best_gwas = coloc_best_gwas[which.max(coloc_best_pp4)]),
                              by = human_symbol]

panel_dt <- data.table(human_symbol = PANEL)
panel_dt[, in_A2_original := human_symbol %in% ORIGINAL_A2]
panel_dt[, in_C2_addition := human_symbol %in% C2_ADDITIONS]
panel_dt[, dual_role_flag := human_symbol %in% DUAL_ROLE_CANDIDATES |
                              human_symbol %in% panel_tf_targets]

sex_sub <- sex_degs[!is.na(human_symbol) & human_symbol %in% PANEL,
                     .(human_symbol, ensembl_id,
                       logFC_M, logFC_F, padj_M, padj_F,
                       interaction_logFC, interaction_padj,
                       sex_class)]
sex_sub <- sex_sub[order(human_symbol, interaction_padj, na.last = TRUE)]
sex_sub <- sex_sub[!duplicated(human_symbol)]
panel_dt <- merge(panel_dt, sex_sub, by = "human_symbol", all.x = TRUE)
panel_dt <- merge(panel_dt, coloc_lookup, by = "human_symbol", all.x = TRUE)

panel_dt[, category := fcase(
  human_symbol %in% c("AR", "ESR1", "ESR2", "PGR", "GPER1", "PRLR"), "receptor",
  human_symbol %in% c("CYP19A1", "SRD5A2", "CYP17A1", "STAR",
                       "HSD17B1", "HSD17B2", "HSD17B3", "HSD17B6", "HSD17B7",
                       "HSD17B12", "HSD3B1"), "biosynthesis",
  human_symbol %in% c("CYP3A4", "CYP1A2", "SULT2A1"), "hepatic_metabolism",
  human_symbol %in% c("SHBG", "GHR"), "transport",
  human_symbol %in% c("KDM6A", "DDX3X", "XIST",
                       "DDX3Y", "EIF1AY", "RPS4Y1"), "sex_chromosome",
  human_symbol %in% c("BCL6", "STAT5A", "STAT5B", "HNF4A", "CUX2"), "TF_GH_axis",
  default = "other"
)]

setcolorder(panel_dt,
            c("human_symbol", "ensembl_id", "category",
              "in_A2_original", "in_C2_addition", "dual_role_flag",
              "logFC_M", "padj_M", "logFC_F", "padj_F",
              "interaction_logFC", "interaction_padj", "sex_class",
              "coloc_best_pp4", "coloc_best_gwas"))
panel_dt <- panel_dt[order(category, human_symbol)]

fwrite(panel_dt, file.path(outdir, "sex_hormone_panel_expanded_clean.csv"))
cat("  Wrote sex_hormone_panel_expanded_clean.csv (", nrow(panel_dt), "rows)\n")
cat("  Dual-role flagged:",
    paste(panel_dt[dual_role_flag == TRUE, human_symbol], collapse = ", "), "\n")
print(panel_dt[, .(human_symbol, category, dual_role_flag, sex_class,
                    interaction_padj, coloc_best_pp4)])

cat("\n=== 217d_clean complete ===\n")
