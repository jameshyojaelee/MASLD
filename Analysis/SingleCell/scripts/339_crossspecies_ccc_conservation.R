#!/usr/bin/env Rscript
# 339_crossspecies_ccc_conservation.R
#
# Analysis L3 (v1) — Cross-species CCC conservation.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   1. Load LIANA human differential LR pairs.
#   2. Load Cross_Species_Concordance atlas (concordance_atlas_unified.csv) with
#      human_symbol × translatability_tier (High/Medium/Low) + n_mouse_sig +
#      mean_h_lfc + dvc_category (Conserved / Mouse_Specific / Human_Specific /
#      Divergent).
#   3. For each LIANA LR pair, annotate ligand + receptor translatability.
#   4. Classify CCC axis:
#        - Conserved (both ligand + receptor Conserved): highest confidence
#          translational LR pair
#        - Divergent (either gene in Divergent category): species-specific biology
#        - Unclassified: gene not in atlas
#   5. Prioritize conserved MASLD-enriched LR axes.
#
# Env: rnaseq
# Outputs: Analysis/SingleCell/results_gpu_v2/crossspecies_ccc/

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
LIANA   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")
CSC     <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUTDIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/crossspecies_ccc")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading LIANA + Cross-species atlas...")
liana <- fread(file.path(LIANA, "liana_differential_interactions.csv"))
atlas <- fread(file.path(CSC, "concordance_atlas_unified.csv"),
               select = c("human_symbol","translatability_tier","primary_category",
                          "dvc_category","n_concordant","n_discordant","n_diets_sig",
                          "mean_h_lfc","dvc_h_significant","dvc_mean_h_lfc",
                          "translatability_score"))
message(sprintf("  LIANA LR pairs: %d; atlas rows: %d", nrow(liana), nrow(atlas)))

message("[2] Annotating LIANA LR pairs with translatability...")
liana[, ligand_first := sapply(strsplit(ligand_complex, "_"), `[`, 1)]
liana[, receptor_first := sapply(strsplit(receptor_complex, "_"), `[`, 1)]

lig_info <- unique(atlas[, .(ligand_first = human_symbol,
                              ligand_tier = translatability_tier,
                              ligand_category = dvc_category,
                              ligand_score = translatability_score,
                              ligand_h_lfc = dvc_mean_h_lfc)])
rec_info <- unique(atlas[, .(receptor_first = human_symbol,
                              receptor_tier = translatability_tier,
                              receptor_category = dvc_category,
                              receptor_score = translatability_score,
                              receptor_h_lfc = dvc_mean_h_lfc)])

liana <- merge(liana, lig_info, by = "ligand_first", all.x = TRUE)
liana <- merge(liana, rec_info, by = "receptor_first", all.x = TRUE)

message("[3] Classify CCC axes by conservation...")
liana[, axis_conservation := fifelse(
  is.na(ligand_category) | is.na(receptor_category), "Unclassified",
  fifelse(ligand_category == "Conserved" & receptor_category == "Conserved",
          "Fully_Conserved",
  fifelse(ligand_category == "Divergent" | receptor_category == "Divergent",
          "Divergent",
  fifelse(grepl("Conserved", ligand_category) | grepl("Conserved", receptor_category),
          "Partially_Conserved", "Other"))))]

message("[4] Enrichment: MASLD-enriched LR pairs in Conserved ligands vs background?")
ligand_universe <- unique(liana$ligand_first)
conserved_ligands <- unique(atlas[primary_category == "Conserved", human_symbol])
masld_enriched_ligands <- unique(liana[score_diff > 0.2, ligand_first])

hyper <- function(set_a, set_b, bg) {
  set_a <- intersect(set_a, bg); set_b <- intersect(set_b, bg)
  k <- length(intersect(set_a, set_b))
  m <- length(set_a); n <- length(bg) - m
  p <- phyper(k - 1, m, n, length(set_b), lower.tail = FALSE)
  or <- (k / length(set_b)) / (m / length(bg))
  data.table(k = k, n_a = m, n_b = length(set_b), n_bg = length(bg),
             odds_ratio = or, pvalue = p)
}
enr <- hyper(conserved_ligands, masld_enriched_ligands, ligand_universe)
enr[, test := "Conserved enriched for MASLD-up LIANA ligands"]

message("[5] Writing outputs...")
fwrite(liana, file.path(OUTDIR, "liana_crossspecies_annotated.csv"))

# Priority: Fully_Conserved + MASLD-enriched
priority <- liana[axis_conservation == "Fully_Conserved" & score_diff > 0.1][order(-score_diff)]
fwrite(priority, file.path(OUTDIR, "conserved_masld_ccc_priority.csv"))

# Species-divergent pairs (potential model artifact indicators)
divergent <- liana[axis_conservation == "Divergent" & abs(score_diff) > 0.1][order(-abs(score_diff))]
fwrite(divergent, file.path(OUTDIR, "divergent_masld_ccc.csv"))

# Summary
cons_count <- liana[, .N, by = axis_conservation][order(-N)]
masld_enriched_by_cons <- liana[abs(score_diff) > 0.1, .N,
                                by = axis_conservation][order(-N)]

summary_lines <- c(
  sprintf("LIANA LR pairs: %d", nrow(liana)),
  sprintf("Cross-species atlas rows: %d", nrow(atlas)),
  "",
  "=== Axis conservation distribution (all LR pairs) ===",
  capture.output(print(cons_count)),
  "",
  "=== MASLD-enriched LR pairs by conservation (|score_diff|>0.1) ===",
  capture.output(print(masld_enriched_by_cons)),
  "",
  "=== Enrichment: Conserved for MASLD ligands ===",
  capture.output(print(enr)),
  "",
  sprintf("Fully-conserved MASLD-enriched CCC axes: %d", nrow(priority)),
  "",
  "Top 30 conserved MASLD-up CCC axes:",
  capture.output(print(priority[1:min(30, .N),
                                .(source, target, ligand_complex, receptor_complex,
                                  score_diff, ligand_category, receptor_category,
                                  ligand_score, receptor_score)], nrows = 30)),
  "",
  sprintf("Species-divergent LR pairs (|score_diff|>0.1): %d", nrow(divergent)),
  "Top 15 divergent axes:",
  capture.output(print(divergent[1:min(15, .N),
                                 .(source, target, ligand_complex, receptor_complex,
                                   score_diff, ligand_category, receptor_category)],
                       nrows = 15)))
writeLines(summary_lines, file.path(OUTDIR, "crossspecies_ccc_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
