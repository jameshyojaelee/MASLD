#!/usr/bin/env Rscript
# KEY MESSAGE: Inherited regulatory susceptibility and established RNA-seq
# remodeling are complementary, without evidence of biological antagonism.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
grDevices::pdf.options(useDingbats = FALSE)

base <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
candidate_root <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
stage_root <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
genetics_root <- normalizePath(Sys.getenv("GENETICS_RELEASE_ROOT", ""), mustWork = TRUE)
interface_root <- normalizePath(Sys.getenv("GENETICS_STATE_INTERFACE_ROOT", ""), mustWork = TRUE)
if (!nzchar(candidate_root) || !dir.exists(candidate_root)) {
  stop("FIGURE_CANDIDATE_ROOT must be a fresh directory", call. = FALSE)
}
source(file.path(base, "scripts/figures/publication_theme.R"))
source(file.path(base, "scripts/figures/load_figure_data.R"))  # is_canonical_deg()

# FIG3G_EXPECT_ADOPTED_COUNTS=0 (a corrected bulk fit or COLOC release): the
# adopted-release counts below are replaced by agreement with the interface run.
expect_adopted <- Sys.getenv("FIG3G_EXPECT_ADOPTED_COUNTS", "1") != "0"
model_manifest <- fread(file.path(stage_root, "model_input_manifest.tsv"))
stopifnot(nrow(model_manifest) == 1L,
          !expect_adopted || model_manifest$n_genes == 23370L,
          model_manifest$n_samples == 844L, model_manifest$n_cohorts == 5L)
# The corrected-control refit records the ORIGINAL fit as accepted_reference_path,
# so its own deg_results.csv must be named explicitly.
bulk_path <- Sys.getenv("FIG3G_BULK_DEG", "")
if (!nzchar(bulk_path)) bulk_path <- model_manifest$accepted_reference_path
bulk_path <- normalizePath(bulk_path, mustWork = TRUE)
coloc_path <- file.path(genetics_root, "gene_level_coloc_tier12.csv")
stopifnot(file.exists(coloc_path))

bulk <- fread(bulk_path)
stopifnot(model_manifest$n_genes == nrow(bulk))
n_genes_bulk <- nrow(bulk)
bulk[, gene_id := sub("\\..*$", "", gene)]
bulk <- unique(bulk, by = "gene_id")
coloc <- fread(coloc_path)
coloc[, gene_id := sub("\\..*$", "", ensembl)]
setorder(coloc, -coloc_best_susie_pp4, -coloc_best_abf_pp4)
coloc <- unique(coloc, by = "gene_id")
coloc_keep <- coloc[, .(gene_id, susie_pp4 = coloc_best_susie_pp4)]
# A corrected COLOC release marks a tier-1/2 SuSiE test that failed coloc's
# shared-posterior check as untestable (NA PP.H4). Those genes stay in the jointly
# testable universe and are not genetic, like tested genes at PP.H4 < 0.5.
has_state <- "susie_state_t12" %in% names(coloc)
if (has_state) coloc_keep[, susie_state := coloc$susie_state_t12]
interface <- merge(
  bulk[, .(gene_id, symbol, logFC, padj, treat_lfc, treat_fdr, AveExpr)],
  coloc_keep,
  by = "gene_id", all = FALSE
)
interface[, genetic := is.finite(susie_pp4) & susie_pp4 >= 0.5]
# Main-panel disease-state call: the paper's canonical rule (padj < 0.05 and
# |log2FC| > 0.5). TREAT (lfc 0.25, FDR < 0.05) stays the named sensitivity in the
# Figure S3G expression-matched comparison.
interface[, rna_sig := is_canonical_deg(interface)]
interface[, treat_sig := is.finite(treat_fdr) & treat_fdr < 0.05]
interface[, state := fifelse(genetic & rna_sig, "Both",
                      fifelse(genetic, "Genetics sig.",
                        fifelse(rna_sig, "Transcriptomics sig.", "Neither")))]
scatter <- interface[is.finite(susie_pp4) & is.finite(logFC)]
rho_test <- cor.test(scatter$logFC, scatter$susie_pp4,
                     method = "spearman", exact = FALSE)
rho <- unname(rho_test$estimate)
threshold_free <- fread(file.path(interface_root, "threshold_free_by_arm.tsv"))[
  bulk_arm == "fragment_23370" & pp4_column == "susie"
]
overlap <- fread(file.path(interface_root, "overlap_by_arm.tsv"))[
  bulk_arm == "fragment_23370" & genetic_definition == "susie_pp4_0.5"
]
stopifnot(!has_state ||
            overlap$n_susie_untestable_joint == sum(interface$susie_state == "untestable"))
stopifnot(!expect_adopted || (nrow(interface) == 14112L && nrow(scatter) == 3617L &&
            sum(interface$genetic) == 428L &&
            sum(interface$genetic & interface$treat_sig) == 30L),
          nrow(overlap) == 1L,
          overlap$n_joint_testable == nrow(interface),
          overlap$n_genetic_joint == sum(interface$genetic),
          nrow(threshold_free) == 1L,
          all(interface$treat_lfc == 0.25),
          threshold_free$n_pairs_raw == nrow(scatter),
          abs(rho - threshold_free$spearman_signed_raw) < 1e-10)
scatter[, state := factor(state, levels = c(
  "Neither", "Transcriptomics sig.", "Genetics sig.", "Both"
))]
setorder(scatter, state)

state_colors <- c(
  Neither = "#D8D8D8",
  `Transcriptomics sig.` = masld_colors$human_enriched,
  `Genetics sig.` = masld_colors$down,
  Both = cat_palette[9]
)
p_scatter <- ggplot(scatter, aes(logFC, susie_pp4)) +
  geom_hline(yintercept = 0.5, color = "grey72", linewidth = 0.25,
             linetype = "22") +
  geom_vline(xintercept = 0, color = "grey82", linewidth = 0.25) +
  geom_point(data = scatter[state != "Both"], aes(color = state),
             size = 0.65, alpha = 0.58, stroke = 0) +
  geom_point(data = scatter[state == "Both"], aes(color = state),
             size = 0.65, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = state_colors, drop = FALSE, name = NULL) +
  guides(color = guide_legend(ncol = 2, byrow = TRUE,
                              override.aes = list(size = 1.4, alpha = 1))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  labs(x = "RNA-seq Log2FC",
       y = "COLOC PP.H4") +
  theme_masld_compact() +
  theme(legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 6, face = "plain"),
        plot.margin = margin(2, 2, 0, 2))

matched <- fread(file.path(interface_root, "matched_treat_by_arm.tsv"))[
  bulk_arm == "fragment_23370" & genetic_definition == "susie_pp4_0.5"
]
matched_pairs <- fread(file.path(interface_root, "matched_pairs_by_arm.tsv"))[
  bulk_arm == "fragment_23370" & genetic_definition == "susie_pp4_0.5"
]
stopifnot(nrow(matched) == 1L,
          nrow(matched_pairs) == matched$n_pairs,
          matched$treat_lfc == 0.25,
          matched$both_treat + matched$discordant_case_only +
            matched$discordant_control_only + matched$neither_treat == matched$n_pairs,
          matched$n_case_treat == matched$both_treat + matched$discordant_case_only,
          matched$n_control_treat == matched$both_treat + matched$discordant_control_only,
          !expect_adopted || (matched$n_pairs == 428L &&
            matched$n_case_treat == 30L && matched$n_control_treat == 21L &&
            matched$discordant_case_only == 30L &&
            matched$discordant_control_only == 21L &&
            matched$both_treat == 0L && matched$neither_treat == 377L))
pair_matrix <- data.table(
  control_rna = c(1, 2, 1, 2),
  genetics_rna = c(2, 2, 1, 1),
  n = c(matched$both_treat, matched$discordant_case_only,
        matched$discordant_control_only, matched$neither_treat),
  state = factor(c("Both", "Genetics sig.", "Transcriptomics sig.", "Neither"),
                 levels = names(state_colors)),
  label_color = c("black", "white", "white", "black")
)
p_matrix <- ggplot(pair_matrix, aes(control_rna, genetics_rna)) +
  geom_tile(aes(fill = state), width = 0.92, height = 0.92,
            color = "white", linewidth = 0.45) +
  geom_text(aes(label = n, color = label_color), size = 9 / .pt,
            family = "Helvetica") +
  scale_fill_manual(values = state_colors, guide = "none") +
  scale_color_identity() +
  scale_x_continuous(breaks = c(1, 2), labels = c("TREAT sig.", "Not sig."),
                     expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(breaks = c(1, 2), labels = c("Not sig.", "TREAT sig."),
                     expand = expansion(mult = c(0.03, 0.03))) +
  labs(title = sprintf("Expression-matched pairs (n = %d)", matched$n_pairs),
       subtitle = sprintf("Paired OR = %.2f; P = %.2f", matched$mcnemar_or,
                          matched$mcnemar_p),
       x = "Matched control", y = "Genetics sig. gene") +
  theme_masld_compact() +
  theme(legend.position = "none",
        plot.title = element_text(size = 6, face = "plain", hjust = 0.5),
        plot.subtitle = element_text(size = 6, face = "plain", hjust = 0.5,
                                     margin = margin(b = 1)),
        axis.title = element_text(size = 6, face = "plain"),
        axis.text = element_text(size = 6, face = "plain"),
        plot.margin = margin(1, 4, 2, 4))

main_out_dir <- file.path(candidate_root, "figure3", "panels")
supp_out_dir <- file.path(candidate_root, "figureS3", "panels")
src_dir <- file.path(candidate_root, "source_tables")
dir.create(main_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(main_out_dir, "fig3g_genetic_state_complementarity.pdf")
supp_out <- file.path(supp_out_dir, "figs3g_expression_matched_pairs.pdf")
ggsave(out, p_scatter, width = 2.45, height = 2.45, device = cairo_pdf)
ggsave(supp_out, p_matrix, width = 2.45, height = 2.15, device = cairo_pdf)

fwrite(scatter[, .(gene_id, symbol, logFC, padj, canonical_deg = rna_sig,
                   treat_lfc, treat_fdr, treat_supported = treat_sig,
                   susie_pp4, state)],
       file.path(src_dir, "fig3g_raw_treat_interface.tsv"), sep = "\t")
summary <- data.table(
  bulk_arm = "fragment_23370",
  n_jointly_testable = nrow(interface),
  n_finite_susie_pairs = nrow(scatter),
  spearman_rho = rho,
  spearman_p = rho_test$p.value,
  n_treat_supported_joint = sum(interface$treat_sig),
  n_genetically_anchored_joint = sum(interface$genetic),
  n_genetically_anchored_treat = sum(interface$genetic & interface$treat_sig),
  n_expression_matched_pairs = matched$n_pairs,
  genetic_only_treat_pairs = matched$discordant_case_only,
  matched_control_only_treat_pairs = matched$discordant_control_only,
  neither_treat_pairs = matched$neither_treat,
  both_treat_pairs = matched$both_treat,
  paired_or = matched$mcnemar_or,
  paired_p = matched$mcnemar_p,
  treat_lfc = matched$treat_lfc
)
# "fragment_23370" is the bulk-arm identifier, not a gene count.
summary[, `:=`(
  n_genes_bulk = n_genes_bulk,
  spearman_logfc = "raw log2FC (deg_results logFC), not shrunk",
  state_rule = "canonical padj < 0.05 and abs(log2FC) > 0.5; TREAT is the Figure S3G sensitivity",
  n_canonical_deg_joint = sum(interface$rna_sig),
  n_genetically_anchored_canonical = sum(interface$genetic & interface$rna_sig),
  n_scatter_both = sum(scatter$state == "Both"),
  n_scatter_transcriptomics_only = sum(scatter$state == "Transcriptomics sig."),
  n_scatter_genetics_only = sum(scatter$state == "Genetics sig."),
  n_scatter_neither = sum(scatter$state == "Neither")
)]
if (has_state) {
  # Own state: in n_jointly_testable, never genetic; without a PP.H4 they are
  # outside the scatter and its Spearman rho only.
  summary[, n_susie_untestable_joint := sum(interface$susie_state == "untestable")]
  summary[, n_susie_untestable_without_pp4 :=
            sum(interface$susie_state == "untestable" & !is.finite(interface$susie_pp4))]
}
fwrite(summary, file.path(src_dir, "fig3g_matched_treat_summary.tsv"), sep = "\t")
fwrite(matched_pairs,
       file.path(src_dir, "fig3g_expression_matched_pairs.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(candidate_root, "sessionInfo.txt"))
cat("FIG3G_PROMOTED_RELEASE_RENDER_COMPLETE", out, "\n")
