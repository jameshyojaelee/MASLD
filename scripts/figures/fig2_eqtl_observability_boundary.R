#!/usr/bin/env Rscript
# ==============================================================================
# Figure 2 — expression-QTL observability boundary (NEW, 2026-08-08)
#
# Purpose
#   Define what the expression-QTL layer of the genetic arm DOES and DOES NOT
#   represent. Three panels:
#     fig2I  positive-only source liver cis-eGene observability across the frozen
#            gene sets; the complement is INDETERMINATE, never negative
#     fig2J  the three prespecified eQTL power universes and their constructibility
#     fig2K  phenotype provenance of the 35 Tier-1/2 GWAS strata, split by whether
#            the regulatory (eQTL) side is ancestry-matched
#
# Historical method constraints: docs/archive/plans/2026-08-07_paper_program/30_GENETICS_CONTEXT_AND_FIG3.md
#   Plan 30 closed as `coverage_limited_terminal`. Seal GEN_TERMINAL_CLOSURE_READY
#   carries context_rescue_authorized=false and negative_claim_authorized=false.
#   AUTHORIZED here: phenotype provenance, POSITIVE-ONLY Broadaway/eGene
#   observability, and the 34/447 descriptive interface.
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
# Outputs: figures/main/fig2_genetics/panels/  (FIG3_DIR == fig2_genetics)
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
OUT <- file.path(FIG3_DIR, "panels")   # FIG3_DIR = figures/main/fig2_genetics
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Seal check — refuse to render if the terminal closure is not the expected one
# ---------------------------------------------------------------------------
seal <- fread(file.path(GEN, "GEN_TERMINAL_CLOSURE_READY"))
stopifnot(seal$status == "coverage_limited_terminal_validated")
stopifnot(identical(tolower(as.character(seal$context_rescue_authorized)), "false"))
stopifnot(identical(tolower(as.character(seal$negative_claim_authorized)),  "false"))
SOURCE_WIDE_EGENES <- as.integer(seal$source_wide_unique_source_defined_egenes)  # 6564
MAPPED_SYMBOL_ROWS <- as.integer(seal$mapped_symbol_annotation_rows)             # 6583
REPRESENTED_ENSG   <- as.integer(seal$represented_unique_source_ensgs)           # 6562
stopifnot(SOURCE_WIDE_EGENES == 6564L, MAPPED_SYMBOL_ROWS == 6583L,
          REPRESENTED_ENSG == 6562L)

manifest <- fread(file.path(GEN, "plan60_terminal_artifacts.tsv"))
stopifnot(nrow(manifest) == 15L, uniqueN(manifest$artifact_id) == 15L)
message("Seal OK: ", seal$release_id, " / ", seal$status)

# ---------------------------------------------------------------------------
# Colours. Control / no-information grey is locked to #9E9E9E.
# ---------------------------------------------------------------------------
COL_POSITIVE      <- "#4C72B0"   # source-positive liver cis-eGene (observed)
COL_INDETERMINATE <- "#9E9E9E"   # observability indeterminate (NOT a negative)
COL_MATCHED       <- "#4C72B0"   # eQTL side ancestry-matched (EUR GWAS)
COL_XANC          <- "#DD8452"   # cross-ancestry eQTL-limited

# ===========================================================================
# PANEL fig2I — positive-only source liver cis-eGene observability
# ===========================================================================
obs <- fread(file.path(GEN, "gene_observability.tsv"),
             select = c("gene_symbol", "joint_testable", "primary_genetic",
                        "established_state_associated", "source_defined_egene",
                        "source_egene_absence_interpretation",
                        "adequate_static_eqtl_negative_authorized"))
lgl <- function(x) tolower(as.character(x)) == "true"

# Hard guard: no row anywhere may authorise a source-eQTL negative.
stopifnot(!any(lgl(obs$adequate_static_eqtl_negative_authorized)))
# Hard guard: the only two absence interpretations are positive / indeterminate.
stopifnot(setequal(unique(obs$source_egene_absence_interpretation),
                   c("source_positive",
                     "indeterminate_complete_tested_universe_not_deposited")))

obs[, `:=`(jt = lgl(joint_testable),
           pg = lgl(primary_genetic),
           es = lgl(established_state_associated),
           eg = lgl(source_defined_egene))]

# Re-derive the frozen contract before plotting anything.
stopifnot(sum(obs$jt) == 14931L)
stopifnot(sum(obs$pg) == 473L, sum(obs$pg & obs$jt) == 447L)
stopifnot(sum(obs$es) == 1915L, sum(obs$es & obs$jt) == 1261L)
stopifnot(sum(obs$pg & obs$es) == 34L, sum(obs$pg & obs$es & obs$jt) == 34L)
stopifnot(sum(obs$eg) == MAPPED_SYMBOL_ROWS)
message("Frozen contract re-derived: 14,931 joint / 447 genetic / 1,261 state / 34 both")

sets <- list(
  "All jointly testable"        = obs[jt == TRUE],
  "Genetically anchored"        = obs[jt == TRUE & pg == TRUE],
  "Established-state associated"= obs[jt == TRUE & es == TRUE],
  "Both"                        = obs[jt == TRUE & pg == TRUE & es == TRUE]
)

d_obs <- rbindlist(lapply(names(sets), function(nm) {
  s <- sets[[nm]]
  data.table(set = nm, n_total = nrow(s),
             n_positive = sum(s$eg), n_indeterminate = sum(!s$eg))
}))
d_obs[, label := sprintf("%s\n(n = %s)", set, trimws(format(n_total, big.mark = ",")))]
d_obs[, label := factor(label, levels = rev(label))]

long_obs <- melt(d_obs, id.vars = c("set", "label", "n_total"),
                 measure.vars = c("n_positive", "n_indeterminate"),
                 variable.name = "state", value.name = "n")
long_obs[, state := factor(
  fifelse(state == "n_positive",
          "Source-positive liver cis-eGene",
          "Observability indeterminate (source tested universe not deposited)"),
  levels = c("Source-positive liver cis-eGene",
             "Observability indeterminate (source tested universe not deposited)"))]
long_obs[, frac := n / n_total]

p_obs <- ggplot(long_obs, aes(x = frac, y = label, fill = state)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = trimws(format(n, big.mark = ","))),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = GEOM_TEXT_6PT, colour = "black", fontface = "plain") +
  scale_fill_manual(values = c(COL_POSITIVE, COL_INDETERMINATE), name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x * 100, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Share of gene set", y = NULL) +
  theme_masld(base_size = 6) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(2.4, "mm"),
        text = element_text(colour = "black", face = "plain"),
        axis.text = element_text(colour = "black"))

save_fig(p_obs, file.path(OUT, "fig2I_eqtl_observability.pdf"),
         width = 4.4, height = 2.5)
fwrite(d_obs, file.path(OUT, "fig2I_eqtl_observability_source.csv"))

message(sprintf(
  paste0("CAPTION fig2I: Positive-only expression-QTL observability across the frozen ",
         "gene sets. Coloured segments are genes with a source-defined liver cis-eGene ",
         "call in the Broadaway liver eQTL meta-analysis (N = 1,183, European); grey ",
         "segments are genes for which observability is INDETERMINATE because the source ",
         "did not deposit its complete tested-gene universe. Grey is not a negative and ",
         "no gene in this figure is called tested_negative. Counts are mapped gene-symbol ",
         "annotation rows (%s rows across the whole table; the source-wide count of unique ",
         "source-defined Ensembl eGenes is %s and the number of unique source Ensembl ",
         "identifiers represented by the symbol rows is %s - these units are not ",
         "interchangeable). Genetically anchored genes are eGene-positive almost by ",
         "construction (%d/%d) because SuSiE colocalization requires an eQTL signal at the ",
         "locus; that is a definitional property of the map, not a result. For the ",
         "established-state set, %d of %d genes have unknown eQTL observability, so the ",
         "question 'does this disease-state gene carry a regulatory genetic effect' is ",
         "unanswerable from the deposited source rather than answered in the negative."),
  format(MAPPED_SYMBOL_ROWS, big.mark = ","),
  format(SOURCE_WIDE_EGENES, big.mark = ","),
  format(REPRESENTED_ENSG, big.mark = ","),
  d_obs[set == "Genetically anchored", n_positive],
  d_obs[set == "Genetically anchored", n_total],
  d_obs[set == "Established-state associated", n_indeterminate],
  d_obs[set == "Established-state associated", n_total]))

# ===========================================================================
# PANEL fig2J — the three prespecified eQTL power universes and what is
# constructible. This is the structural boundary: 2 of 3 cannot be built.
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

save_fig(p_pw, file.path(OUT, "fig2J_eqtl_power_universes.pdf"),
         width = 4.8, height = 1.9)
fwrite(cbind(d_pw[, .(universe = as.character(universe), status, n_genes)],
             pw[, .(reason)]),
       file.path(OUT, "fig2J_eqtl_power_universes_source.csv"))

message(paste0(
  "CAPTION fig2J: Only one of the three prespecified expression-QTL power universes ",
  "can be constructed from the deposited source. The bulk-expressed/eQTL-tested and ",
  "source-significant-with-covariates universes are coverage-limited: the liver eQTL ",
  "source deposits significant leads but not its complete tested-gene universe, source ",
  "expression, or local tested-variant density. Consequently no powered eQTL negative, ",
  "matched-control analysis, or enrichment test is reported anywhere in this figure. ",
  "This is a limit of the deposited data, not evidence that the untested genes lack ",
  "regulatory genetic effects."))

# ===========================================================================
# PANEL fig2K — phenotype provenance of the 35 Tier-1/2 strata, split by
# whether the regulatory (eQTL) side is ancestry-matched.
# ===========================================================================
reg <- fread(file.path(GEN, "phenotype_registry.tsv"))
stopifnot(nrow(reg) == 35L, all(reg$placement == "main"))
stopifnot(uniqueN(reg$eqtl_panel) == 1L)          # single EUR eQTL panel for all 35
EQTL_PANEL <- unique(reg$eqtl_panel)

strat_label <- c(
  alt_ast_or_ggt = "ALT / AST / GGT\n(liver enzyme)",
  direct_masld_mash_diagnosis = "Direct MASLD / MASH\ndiagnosis",
  mri_pdff_or_histologic_steatosis = "MRI-PDFF / histologic\nsteatosis")
anc_label <- c(
  ancestry_matched_eur = "eQTL ancestry-matched (European GWAS)",
  cross_ancestry_eqtl_limited = "Cross-ancestry eQTL-limited (non-European GWAS, European eQTL)")

d_reg <- reg[, .N, by = .(phenotype_stratum, regulatory_ancestry_status)]
d_reg[, stratum := factor(strat_label[phenotype_stratum],
                          levels = rev(strat_label))]
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
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)),
                     breaks = seq(0, 20, 5)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "GWAS strata (35 prespecified Tier-1/2)", y = NULL) +
  theme_masld(base_size = 6) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(2.4, "mm"),
        text = element_text(colour = "black", face = "plain"),
        axis.text = element_text(colour = "black", lineheight = 1.05))

save_fig(p_reg, file.path(OUT, "fig2K_phenotype_provenance.pdf"),
         width = 4.4, height = 2.0)
fwrite(d_reg[, .(phenotype_stratum, regulatory_ancestry_status, n_studies = N)],
       file.path(OUT, "fig2K_phenotype_provenance_source.csv"))

message(sprintf(
  paste0("CAPTION fig2K: Phenotype provenance of the %d prespecified Tier-1/2 GWAS strata. ",
         "%d strata are liver-enzyme traits (ALT/AST/GGT), %d are direct MASLD/MASH ",
         "diagnoses and %d are MRI-PDFF/histologic steatosis, so most of the genetic map ",
         "is anchored on enzyme proxies rather than on diagnosed disease. Every stratum is ",
         "colocalized against a single European liver eQTL panel (%s); the %d non-European ",
         "strata are therefore cross-ancestry eQTL-limited and show locus-level replication, ",
         "NOT ancestry-matched regulatory validation."),
  nrow(reg),
  d_reg[phenotype_stratum == "alt_ast_or_ggt", sum(N)],
  d_reg[phenotype_stratum == "direct_masld_mash_diagnosis", sum(N)],
  d_reg[phenotype_stratum == "mri_pdff_or_histologic_steatosis", sum(N)],
  EQTL_PANEL,
  d_reg[regulatory_ancestry_status == "cross_ancestry_eqtl_limited", sum(N)]))

message("DONE. Wrote 3 panels + 3 source CSVs to ", OUT)
