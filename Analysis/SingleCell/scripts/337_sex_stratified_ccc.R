#!/usr/bin/env Rscript
# 337_sex_stratified_ccc.R
#
# Analysis L1 (v1) — Sex-stratified CCC + macrophage trajectories.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   1. Load existing sex-dimorphic DEGs (Script 26: female-biased, male-biased,
#      divergent, concordant).
#   2. Load LIANA differential LR pairs.
#   3. For each LR pair, annotate whether ligand / receptor is sex-dimorphic
#      and which direction (female vs male).
#   4. Enrichment test: are sex-biased DEGs enriched in CCC ligand/receptor
#      gene sets?
#   5. Prioritize sex-specific CCC axes (e.g., Female_biased ligand +
#      MASLD-enriched score_diff > 0 = female-specific progressor axis).
#
# Env: rnaseq
# Outputs: Analysis/SingleCell/results_gpu_v2/sex_ccc/

suppressPackageStartupMessages({
  library(data.table)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_RES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LIANA   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")
OUTDIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/sex_ccc")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading sex-dimorphic DEG classifications (Script 26)...")
sex_file <- file.path(INT_RES, "sex_deg_classification.csv")
if (!file.exists(sex_file)) {
  # try alternate
  sex_file_alt <- list.files(INT_RES, pattern = "sex.*\\.csv", full.names = TRUE)
  sex_file <- sex_file_alt[1]
}
sex <- fread(sex_file)
message(sprintf("  Sex DEG file: %s (%d rows)", basename(sex_file), nrow(sex)))
message(sprintf("  Columns: %s", paste(names(sex), collapse = ", ")))

# Sex file uses ENSG in 'gene' column — join to symbols via bulk dream.
# Preferred column for classification: sex_class (has Female_biased / Male_biased / Divergent / Concordant)
if (!"sex_class" %in% names(sex) && "sex_class_stratified" %in% names(sex)) {
  setnames(sex, "sex_class_stratified", "sex_class")
}
sex[, ensg_base := sub("\\.\\d+$", "", gene)]

# Pull symbols from bulk dream
bulk_sym <- fread(file.path(INT_RES, "dream_results_ashr.csv"),
                  select = c("gene","symbol"))
bulk_sym[, ensg_base := sub("\\.\\d+$", "", gene)]
sex <- merge(sex, unique(bulk_sym[, .(ensg_base, symbol)]),
             by = "ensg_base", all.x = TRUE)

message(sprintf("  Sex class distribution: "))
print(sex[, .N, by = sex_class][order(-N)])

message("[2] Loading LIANA differential LR pairs...")
liana <- fread(file.path(LIANA, "liana_differential_interactions.csv"))
ligands <- unique(liana$ligand_complex)
receptors <- unique(liana$receptor_complex)
split_complex <- function(x) unlist(strsplit(x, "_"))
ligand_genes   <- unique(unlist(lapply(ligands, split_complex)))
receptor_genes <- unique(unlist(lapply(receptors, split_complex)))

message("[3] Annotating LIANA LR pairs with sex-class of ligand + receptor...")
sex_lookup <- unique(sex[, .(symbol, sex_class)])
liana[, ligand_first := sapply(strsplit(ligand_complex, "_"), `[`, 1)]
liana[, receptor_first := sapply(strsplit(receptor_complex, "_"), `[`, 1)]
liana <- merge(liana, sex_lookup[, .(symbol, ligand_sex_class = sex_class)],
               by.x = "ligand_first", by.y = "symbol", all.x = TRUE)
liana <- merge(liana, sex_lookup[, .(symbol, receptor_sex_class = sex_class)],
               by.x = "receptor_first", by.y = "symbol", all.x = TRUE)

# Flag CCC axes
liana[, axis_sex_bias := fifelse(
  !is.na(ligand_sex_class) & ligand_sex_class == "Female_biased", "Female_ligand",
  fifelse(!is.na(ligand_sex_class) & ligand_sex_class == "Male_biased", "Male_ligand",
  fifelse(!is.na(receptor_sex_class) & receptor_sex_class == "Female_biased", "Female_receptor",
  fifelse(!is.na(receptor_sex_class) & receptor_sex_class == "Male_biased", "Male_receptor",
  "Not_sex_biased"))))]

message("[4] Enrichment tests...")
hyper_test <- function(set_a, set_b, bg) {
  set_a <- intersect(set_a, bg); set_b <- intersect(set_b, bg)
  k <- length(intersect(set_a, set_b))
  m <- length(set_a); n <- length(bg) - m
  p <- phyper(k - 1, m, n, length(set_b), lower.tail = FALSE)
  or <- (k / length(set_b)) / (m / length(bg))
  data.table(k = k, n_a = m, n_b = length(set_b), n_bg = length(bg),
             odds_ratio = or, pvalue = p)
}

# Build universe as all expressed genes in bulk dream
bulk <- fread(file.path(INT_RES, "dream_results_ashr.csv"), select = "symbol")
universe <- unique(bulk$symbol)

female_genes <- sex[sex_class == "Female_biased", unique(symbol)]
male_genes   <- sex[sex_class == "Male_biased",   unique(symbol)]
div_genes    <- sex[sex_class %in% c("Divergent","divergent"), unique(symbol)]

enr <- rbind(
  hyper_test(female_genes, ligand_genes, universe)[, test := "Female-biased enriched for CCC ligands"],
  hyper_test(male_genes,   ligand_genes, universe)[, test := "Male-biased enriched for CCC ligands"],
  hyper_test(female_genes, receptor_genes, universe)[, test := "Female-biased enriched for CCC receptors"],
  hyper_test(male_genes,   receptor_genes, universe)[, test := "Male-biased enriched for CCC receptors"],
  hyper_test(div_genes,    c(ligand_genes, receptor_genes), universe)[, test := "Divergent enriched for CCC"],
  fill = TRUE
)
fwrite(enr, file.path(OUTDIR, "sex_ccc_enrichment.csv"))

message("[5] Priority sex-specific MASLD-enriched CCC axes...")
priority <- liana[abs(score_diff) > 0.1 & axis_sex_bias != "Not_sex_biased"]
priority <- priority[order(-score_diff)]
fwrite(priority[, .(source, target, ligand_complex, receptor_complex, score_diff,
                     ligand_sex_class, receptor_sex_class, axis_sex_bias)],
       file.path(OUTDIR, "sex_specific_ccc_axes.csv"))

# Counts by axis sex class
axis_counts <- priority[, .N, by = axis_sex_bias][order(-N)]

summary_lines <- c(
  sprintf("Sex-dimorphic DEGs — Female_biased: %d, Male_biased: %d, Divergent: %d",
          length(female_genes), length(male_genes), length(div_genes)),
  sprintf("LIANA LR pairs: %d", nrow(liana)),
  "",
  "=== Enrichment tests ===",
  capture.output(print(enr[, .(test, k, n_a, n_b, odds_ratio, pvalue)], nrows = 10)),
  "",
  sprintf("Priority sex-specific MASLD-enriched CCC axes: %d", nrow(priority)),
  "",
  "=== Counts by axis sex-bias ===",
  capture.output(print(axis_counts)),
  "",
  "=== Top 30 female-ligand MASLD-enriched axes ===",
  capture.output(print(priority[axis_sex_bias == "Female_ligand"][1:30,
                                .(source, target, ligand_complex, receptor_complex,
                                  score_diff, ligand_sex_class, receptor_sex_class)],
                       nrows = 30)),
  "",
  "=== Top 20 male-ligand MASLD-enriched axes ===",
  capture.output(print(priority[axis_sex_bias == "Male_ligand"][1:20,
                                .(source, target, ligand_complex, receptor_complex,
                                  score_diff, ligand_sex_class, receptor_sex_class)],
                       nrows = 20))
)
writeLines(summary_lines, file.path(OUTDIR, "sex_ccc_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
