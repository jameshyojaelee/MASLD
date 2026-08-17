#!/usr/bin/env Rscript
# ==============================================================================
# Figure 2 phenotype provenance and supplementary eQTL boundary
#
# Purpose
#   Render phenotype provenance as Figure 2J and retain the supplementary
#   deposited-data boundary showing which eQTL power universes are constructible.
#
# Historical method constraints: docs/archive/plans/2026-08-07_paper_program/30_GENETICS_CONTEXT_AND_FIG3.md
#   Plan 30 closed as `coverage_limited_terminal`. Seal GEN_TERMINAL_CLOSURE_READY
#   carries context_rescue_authorized=false and negative_claim_authorized=false.
#   AUTHORIZED here: phenotype provenance and coverage of the prespecified
#   eQTL power universes.
#   PROHIBITED and deliberately absent from this script: context-rescued fractions,
#   enrichment claims or tests, matched-rescue analyses, matched nulls, and any
#   genetic or context NEGATIVE claim. No Fisher/OR/p-value is computed or drawn.
#   Missing coverage is drawn as `indeterminate`, which is NOT a negative.
#
#   Identifier units are NOT interchangeable (Plan 30 lines 49-56):
#     6,564 = source-wide unique source-defined Ensembl eGenes
#     6,583 = mapped gene-symbol annotation rows in gene_observability.tsv
#     6,562 = unique source Ensembl identifiers represented by those rows
#   This script plots SYMBOL ROWS and labels them as such. It never reports
#   6,583 as a count of unique source eGenes.
#
# Inputs: signed 15-row Plan 60 handoff only (plan60_terminal_artifacts.tsv).
# Outputs: Figure 2J and FigS2_eqtl_power_universes.pdf
# Run: sbatch scripts/figures/run_fig2_eqtl_observability.sbatch
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

GEN <- file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates",
                 "program-context-v2-candidate-2026-08-07/genetics_context")
candidate_out <- Sys.getenv("FIG2_CANDIDATE_DIR", "")
OUT <- if (nzchar(candidate_out)) candidate_out else file.path(FIG3_DIR, "panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Seal check — refuse to render if the terminal closure is not the expected one
# ---------------------------------------------------------------------------
seal <- fread(file.path(GEN, "GEN_TERMINAL_CLOSURE_READY"))
stopifnot(seal$status == "coverage_limited_terminal_validated")
stopifnot(identical(tolower(as.character(seal$context_rescue_authorized)), "false"))
stopifnot(identical(tolower(as.character(seal$negative_claim_authorized)),  "false"))
manifest <- fread(file.path(GEN, "plan60_terminal_artifacts.tsv"))
stopifnot(nrow(manifest) == 15L, uniqueN(manifest$artifact_id) == 15L)
message("Seal OK: ", seal$release_id, " / ", seal$status)

# ---------------------------------------------------------------------------
# Colours. Control / no-information grey is locked to #9E9E9E.
# ---------------------------------------------------------------------------
COL_POSITIVE      <- "#4C72B0"   # constructible universe
COL_INDETERMINATE <- "#9E9E9E"   # observability indeterminate (NOT a negative)
COL_MATCHED       <- "#4C72B0"   # eQTL side ancestry-matched (EUR GWAS)
COL_XANC          <- "#DD8452"   # cross-ancestry eQTL-limited

# ===========================================================================
# SUPPLEMENTARY BOUNDARY — the three prespecified eQTL power universes and what
# is constructible. This panel was replaced in the main figure by the
# PIP-composition Figure 2D on 2026-08-13.
# ===========================================================================
pw <- fread(file.path(GEN, "power_stratified_interface.tsv"))
stopifnot(nrow(pw) == 3L)
stopifnot(!any(tolower(as.character(pw$source_negative_authorized)) == "true"))

univ_label <- c(
  all_joint_testable = "1. All jointly testable genes",
  bulk_expressed_and_source_eqtl_tested = "2. Bulk-expressed AND source eQTL-tested",
  source_significant_egene_complete_covariates = "3. Source-significant liver cis-eGenes with\ncomplete covariates")
blocker <- c(
  all_joint_testable = "constructible",
  bulk_expressed_and_source_eqtl_tested =
    "per-gene tested / non-eGene status\nnot deposited (18,322 source\nexpression genes unresolved)",
  source_significant_egene_complete_covariates =
    "source expression and local\ntested-variant density\nnot deposited")

d_pw <- data.table(
  universe = univ_label[pw$universe],
  status   = fifelse(pw$status == "pass", "Constructible", "Coverage-limited: not constructible"),
  n_genes  = suppressWarnings(as.integer(pw$n_genes)),
  blocker  = blocker[pw$universe])
d_pw[blocker == "constructible", blocker := "constructible from the deposited source"]
d_pw[, universe := factor(universe, levels = rev(univ_label))]
d_pw[, status := factor(status, levels = c("Constructible",
                                           "Coverage-limited: not constructible"))]
d_pw[, cell := fifelse(is.na(n_genes), "-", format(n_genes, big.mark = ","))]

p_pw <- ggplot(d_pw, aes(x = 1, y = universe, fill = status)) +
  geom_tile(width = 0.96, height = 0.86, colour = "white", linewidth = 0.3) +
  geom_text(aes(x = 1, label = cell), size = GEOM_TEXT_6PT,
            colour = "black", fontface = "plain") +
  geom_text(aes(x = 1.72, label = blocker), size = GEOM_TEXT_6PT,
            colour = "black", fontface = "plain", hjust = 0, lineheight = 1.05) +
  scale_fill_manual(values = c("Constructible" = COL_POSITIVE,
                               "Coverage-limited: not constructible" = COL_INDETERMINATE),
                    name = NULL) +
  scale_x_continuous(limits = c(0.42, 3.75),
                     breaks = c(1, 1.72),
                     labels = c("Genes", "Why the universe cannot be built"),
                     position = "top", expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(colour = "black", hjust = 0),
        axis.text.y = element_text(colour = "black", hjust = 1, lineheight = 1.05,
                                   margin = margin(r = 3)),
        axis.ticks = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(2.4, "mm"),
        text = element_text(colour = "black", face = "plain"))

save_fig(p_pw, file.path(OUT, "FigS2_eqtl_power_universes.pdf"),
         width = 4.8, height = 1.9)
fwrite(cbind(d_pw[, .(universe = as.character(universe), status, n_genes)],
             pw[, .(reason)]),
       file.path(OUT, "FigS2_eqtl_power_universes_source.csv"))

message(paste0(
  "CAPTION FigS2: Only one of the three prespecified expression-QTL power universes ",
  "can be constructed from the deposited source. The bulk-expressed/eQTL-tested and ",
  "source-significant-with-covariates universes are coverage-limited: the liver eQTL ",
  "source deposits significant leads but not its complete tested-gene universe, source ",
  "expression, or local tested-variant density. Consequently no powered eQTL negative, ",
  "matched-control analysis, or enrichment test is reported anywhere in this figure. ",
  "This is a limit of the deposited data, not evidence that the untested genes lack ",
  "regulatory genetic effects."))

# ===========================================================================
# PANEL Fig2J — phenotype provenance of the 35 Tier-1/2 strata, split by
# whether the regulatory (eQTL) side is ancestry-matched.
# ===========================================================================
reg <- fread(file.path(GEN, "phenotype_registry.tsv"))
stopifnot(nrow(reg) == 35L, all(reg$placement == "main"))
stopifnot(uniqueN(reg$eqtl_panel) == 1L)
EQTL_PANEL <- unique(reg$eqtl_panel)

strat_label <- c(
  alt_ast_or_ggt = "ALT / AST / GGT\n(liver enzyme)",
  direct_masld_mash_diagnosis = "Direct MASLD / MASH\ndiagnosis",
  mri_pdff_or_histologic_steatosis = "MRI-PDFF / histologic\nsteatosis")
anc_label <- c(
  ancestry_matched_eur = "eQTL ancestry-matched (European GWAS)",
  cross_ancestry_eqtl_limited = "Cross-ancestry eQTL-limited (non-European GWAS, European eQTL)")

d_reg <- reg[, .N, by = .(phenotype_stratum, regulatory_ancestry_status)]
d_reg[, stratum := factor(strat_label[phenotype_stratum], levels = rev(strat_label))]
d_reg[, anc := factor(anc_label[regulatory_ancestry_status], levels = anc_label)]
d_reg[, tot := sum(N), by = stratum]

stopifnot(sum(d_reg$N) == 35L)
stopifnot(d_reg[phenotype_stratum == "alt_ast_or_ggt", sum(N)] == 20L)
stopifnot(d_reg[phenotype_stratum == "direct_masld_mash_diagnosis", sum(N)] == 12L)
stopifnot(d_reg[phenotype_stratum == "mri_pdff_or_histologic_steatosis", sum(N)] == 3L)

p_reg <- ggplot(d_reg, aes(x = N, y = stratum, fill = anc)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5, reverse = TRUE),
            size = GEOM_TEXT_6PT, colour = "black", fontface = "plain") +
  scale_fill_manual(values = c(COL_MATCHED, COL_XANC), name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)), breaks = seq(0, 20, 5)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "GWAS strata (35 prespecified Tier-1/2)", y = NULL) +
  theme_masld(base_size = 6) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(2.4, "mm"),
        text = element_text(colour = "black", face = "plain"),
        axis.text = element_text(colour = "black", lineheight = 1.05))

save_fig(p_reg, file.path(OUT, "Fig2J_phenotype_provenance.pdf"),
         width = 4.4, height = 2.0)
fwrite(d_reg[, .(phenotype_stratum, regulatory_ancestry_status, n_studies = N)],
       file.path(OUT, "Fig2J_phenotype_provenance_source.csv"))

message(sprintf(
  paste0("CAPTION Fig2J: Phenotype provenance of the %d prespecified Tier-1/2 GWAS strata. ",
         "%d strata are liver-enzyme traits, %d are direct MASLD/MASH diagnoses and %d are ",
         "MRI-PDFF/histologic steatosis. Every stratum is colocalized against one European ",
         "liver eQTL panel (%s); the %d non-European strata are therefore cross-ancestry ",
         "eQTL-limited, not ancestry-matched regulatory validation."),
  nrow(reg),
  d_reg[phenotype_stratum == "alt_ast_or_ggt", sum(N)],
  d_reg[phenotype_stratum == "direct_masld_mash_diagnosis", sum(N)],
  d_reg[phenotype_stratum == "mri_pdff_or_histologic_steatosis", sum(N)],
  EQTL_PANEL,
  d_reg[regulatory_ancestry_status == "cross_ancestry_eqtl_limited", sum(N)]))

message("DONE. Wrote the supplementary eQTL-boundary panel and Figure 2J to ", OUT)
