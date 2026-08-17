#!/usr/bin/env Rscript
# KEY MESSAGE: Inherited regulatory susceptibility and established disease-state
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

model_manifest <- fread(file.path(stage_root, "model_input_manifest.tsv"))
stopifnot(nrow(model_manifest) == 1L, model_manifest$n_genes == 23370L,
          model_manifest$n_samples == 844L, model_manifest$n_cohorts == 5L)
bulk_path <- normalizePath(model_manifest$accepted_reference_path, mustWork = TRUE)
coloc_path <- file.path(genetics_root, "gene_level_coloc_tier12.csv")
stopifnot(file.exists(coloc_path))

bulk <- fread(bulk_path)
bulk[, gene_id := sub("\\..*$", "", gene)]
bulk <- unique(bulk, by = "gene_id")
coloc <- fread(coloc_path)
coloc[, gene_id := sub("\\..*$", "", ensembl)]
setorder(coloc, -coloc_best_susie_pp4, -coloc_best_abf_pp4)
coloc <- unique(coloc, by = "gene_id")
interface <- merge(
  bulk[, .(gene_id, symbol, logFC, treat_lfc, treat_fdr, AveExpr)],
  coloc[, .(gene_id, susie_pp4 = coloc_best_susie_pp4)],
  by = "gene_id", all = FALSE
)
interface[, genetic := is.finite(susie_pp4) & susie_pp4 >= 0.5]
interface[, disease_state := is.finite(treat_fdr) & treat_fdr < 0.05]
interface[, state := fifelse(genetic & disease_state, "Both",
                      fifelse(genetic, "Genetics sig.",
                        fifelse(disease_state, "Transcriptomics sig.", "Neither")))]
scatter <- interface[is.finite(susie_pp4) & is.finite(logFC)]
rho_test <- cor.test(scatter$logFC, scatter$susie_pp4,
                     method = "spearman", exact = FALSE)
rho <- unname(rho_test$estimate)
threshold_free <- fread(file.path(interface_root, "threshold_free_by_arm.tsv"))[
  bulk_arm == "fragment_23370" & pp4_column == "susie"
]
stopifnot(nrow(interface) == 14112L, nrow(scatter) == 3617L,
          nrow(threshold_free) == 1L,
          all(interface$treat_lfc == 0.25),
          sum(interface$genetic) == 428L,
          sum(interface$genetic & interface$disease_state) == 30L,
          threshold_free$n_pairs_raw == nrow(scatter),
          abs(rho - threshold_free$spearman_signed_raw) < 1e-10)
scatter[, state := factor(state, levels = c(
  "Neither", "Transcriptomics sig.", "Genetics sig.", "Both"
))]
setorder(scatter, state)

state_colors <- c(
  Neither = "#D8D8D8",
  `Transcriptomics sig.` = masld_colors$mash,
  `Genetics sig.` = masld_colors$down,
  Both = masld_colors$masl
)
p_scatter <- ggplot(scatter, aes(logFC, susie_pp4)) +
  geom_hline(yintercept = 0.5, color = "grey72", linewidth = 0.25,
             linetype = "22") +
  geom_vline(xintercept = 0, color = "grey82", linewidth = 0.25) +
  geom_point(aes(color = state), size = 0.4, alpha = 0.58, stroke = 0) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("ρ = %.03f; P = %.02f", rho, rho_test$p.value),
           hjust = -0.05, vjust = 1.25, size = 6 / .pt, family = "Helvetica") +
  scale_color_manual(values = state_colors, drop = FALSE, name = NULL) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE,
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
stopifnot(nrow(matched) == 1L, matched$n_pairs == 428L,
          nrow(matched_pairs) == 428L,
          matched$treat_lfc == 0.25,
          matched$n_case_treat == 30L, matched$n_control_treat == 21L,
          matched$discordant_case_only == 30L,
          matched$discordant_control_only == 21L,
          matched$both_treat == 0L, matched$neither_treat == 377L)
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
  scale_x_continuous(breaks = c(1, 2), labels = c("RNA-seq sig.", "Not sig."),
                     expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(breaks = c(1, 2), labels = c("Not sig.", "RNA-seq sig."),
                     expand = expansion(mult = c(0.03, 0.03))) +
  labs(title = "Expression-matched pairs (n = 428)",
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

panel <- p_scatter / p_matrix + plot_layout(heights = c(4.1, 1.7))
out_dir <- file.path(candidate_root, "figure3", "panels")
src_dir <- file.path(candidate_root, "source_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "fig3g_genetic_state_complementarity.pdf")
ggsave(out, panel, width = 4.0, height = 2.85, device = cairo_pdf)

fwrite(scatter[, .(gene_id, symbol, logFC, treat_lfc, treat_fdr,
                   susie_pp4, state)],
       file.path(src_dir, "fig3g_raw_treat_interface.tsv"), sep = "\t")
summary <- data.table(
  bulk_arm = "fragment_23370",
  n_jointly_testable = nrow(interface),
  n_finite_susie_pairs = nrow(scatter),
  spearman_rho = rho,
  spearman_p = rho_test$p.value,
  n_treat_supported_joint = sum(interface$disease_state),
  n_genetically_anchored_joint = sum(interface$genetic),
  n_genetically_anchored_treat = sum(interface$genetic & interface$disease_state),
  n_expression_matched_pairs = matched$n_pairs,
  genetic_only_treat_pairs = matched$discordant_case_only,
  matched_control_only_treat_pairs = matched$discordant_control_only,
  neither_treat_pairs = matched$neither_treat,
  both_treat_pairs = matched$both_treat,
  paired_or = matched$mcnemar_or,
  paired_p = matched$mcnemar_p,
  treat_lfc = matched$treat_lfc
)
fwrite(summary, file.path(src_dir, "fig3g_matched_treat_summary.tsv"), sep = "\t")
fwrite(matched_pairs,
       file.path(src_dir, "fig3g_expression_matched_pairs.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(candidate_root, "sessionInfo.txt"))
cat("FIG3G_PROMOTED_RELEASE_RENDER_COMPLETE", out, "\n")
