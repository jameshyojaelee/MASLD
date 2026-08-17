#!/usr/bin/env Rscript
# KEY MESSAGE: Inherited regulatory evidence and established disease-state expression identify different genes without evidence that genetically anchored genes are depleted for disease responsiveness.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

base <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
preview <- identical(tolower(Sys.getenv("FIG3G_DESIGN_PREVIEW", "false")), "true")
candidate_root <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
if (!preview) {
  stop("Figure 3G is fail-closed until corrected COLOC is promoted; set FIG3G_DESIGN_PREVIEW=true only for a non-promotable layout preview",
       call. = FALSE)
}
if (!nzchar(candidate_root) || !dir.exists(candidate_root)) {
  stop("FIGURE_CANDIDATE_ROOT must be a fresh design-preview directory", call. = FALSE)
}

source(file.path(base, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)

# The preview uses the current official release plus its independently generated
# expression-matched response analysis. It is not a substitute for the pending
# corrected-COLOC rerun and cannot be promoted.
evidence_path <- file.path(base,
  "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv")
bulk_path <- file.path(base,
  "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-POOLED-REPRO/deg_results.csv")
match_root <- file.path(base, "RNA-seq/results/response_amplitude/20260813T083158")

evidence <- fread(evidence_path, select = c(
  "symbol", "ensembl_bulk", "joint_testable", "max_susie_pp4_all", "in_susie_primary"
))
bulk <- fread(bulk_path, select = c("gene", "logFC", "padj", "AveExpr"))
bulk[, gene_id := sub("\\..*$", "", gene)]
evidence[, gene_id := sub("\\..*$", "", ensembl_bulk)]
interface <- merge(
  evidence[joint_testable == TRUE],
  bulk[, .(gene_id, logFC, padj, AveExpr)],
  by = "gene_id", all = FALSE
)
interface[, state := fifelse(in_susie_primary == TRUE & padj < 0.05 & abs(logFC) > 0.5,
                             "Both",
                      fifelse(in_susie_primary == TRUE, "Genetically anchored",
                       fifelse(padj < 0.05 & abs(logFC) > 0.5,
                               "State-associated", "Other")))]

scatter <- interface[is.finite(max_susie_pp4_all) & is.finite(logFC)]
rho <- cor.test(scatter$logFC, scatter$max_susie_pp4_all,
                method = "spearman", exact = FALSE)
message(sprintf("[interface audit] n=%d rho=%.12f",
                nrow(scatter), unname(rho$estimate)))
if (nrow(scatter) != 3713L || abs(unname(rho$estimate) - 0.017257910624) > 1e-8) {
  stop("Design-preview threshold-free interface drift", call. = FALSE)
}
scatter[, state := factor(state, levels = c(
  "Other", "State-associated", "Genetically anchored", "Both"
))]

state_colors <- c(
  Other = "#D8D8D8",
  `State-associated` = "#C9265E",
  `Genetically anchored` = "#1565C0",
  Both = "#2A8C7D"
)
p_scatter <- ggplot(scatter, aes(logFC, max_susie_pp4_all)) +
  geom_hline(yintercept = 0.5, color = "grey70", linewidth = 0.25,
             linetype = "22") +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.25) +
  geom_point(aes(color = state), size = 0.38, alpha = 0.55, stroke = 0) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              color = "black", linewidth = 0.45, span = 0.8) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("\u03c1 = %.02f", unname(rho$estimate)),
           hjust = -0.05, vjust = 1.2, size = 6 / .pt,
           family = "Helvetica") +
  scale_color_manual(values = state_colors, drop = FALSE, name = NULL) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 1.4, alpha = 1))) +
  coord_cartesian(xlim = quantile(scatter$logFC, c(0.005, 0.995), na.rm = TRUE)) +
  labs(x = expression("Disease-state "*log[2]*"FC"), y = "COLOC PP.H4") +
  theme_masld_compact() +
  theme(legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 6, face = "plain"),
        plot.margin = margin(2, 2, 2, 2))

ladder <- fread(file.path(match_root, "matched_ladder.tsv"))
canonical_match <- ladder[arm == "primary_aveexpr" & abs(lfc_floor - 0.5) < 1e-12]
if (nrow(canonical_match) != 1L) stop("Matched-response canonical row drift", call. = FALSE)
matched <- data.table(
  group = factor(c("Genetically anchored", "Expression-matched"),
                 levels = c("Expression-matched", "Genetically anchored")),
  n_deg = c(canonical_match$n_treated_de, canonical_match$n_control_de),
  n_total = c(447L, 447L)
)
matched[, fraction := n_deg / n_total]
message(sprintf("[matched audit] treated=%d/%d control=%d/%d",
                matched$n_deg[[1]], matched$n_total[[1]],
                matched$n_deg[[2]], matched$n_total[[2]]))
if (!identical(matched$n_deg, c(23L, 24L)) || any(matched$n_total != 447L)) {
  stop("Official-release matched response counts drift", call. = FALSE)
}

p_match <- ggplot(matched, aes(fraction, group)) +
  geom_segment(aes(x = 0, xend = fraction, yend = group),
               linewidth = 2.8, color = "#D8D8D8", lineend = "butt") +
  geom_point(aes(shape = group), size = 2.0, stroke = 0.45,
             fill = "white", color = "#2A8C7D") +
  geom_text(aes(x = fraction + 0.004, label = sprintf("%d/%d", n_deg, n_total)),
            hjust = 0, size = 6 / .pt, family = "Helvetica") +
  scale_shape_manual(values = c(`Expression-matched` = 21,
                                `Genetically anchored` = 22), guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.075), breaks = c(0, 0.05)) +
  labs(x = "Disease-state associated", y = NULL) +
  theme_masld_compact() +
  theme(axis.text.y = element_text(size = 6, face = "plain"),
        panel.grid = element_blank(), plot.margin = margin(2, 5, 2, 2))

panel <- p_scatter | p_match + plot_layout(widths = c(1.8, 1))
out_dir <- file.path(candidate_root, "figure3", "panels")
src_dir <- file.path(candidate_root, "source_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "fig3g_genetic_state_complementarity_DESIGN_PREVIEW.pdf")
ggsave(out, panel, width = 4.55, height = 2.15, device = cairo_pdf)

fwrite(scatter[, .(gene_id, symbol, logFC, max_susie_pp4_all, state)],
       file.path(src_dir, "fig3g_threshold_free_interface_DESIGN_PREVIEW.tsv"), sep = "\t")
fwrite(matched,
       file.path(src_dir, "fig3g_matched_response_DESIGN_PREVIEW.tsv"), sep = "\t")
writeLines(c(
  "promotion_state=design_preview_not_promotable",
  "genetics_release=2026-07-15-r2_current_official_not_corrected_rerun",
  "bulk_release=resource-f-five-coloc-v6-candidate-2026-08-10/BULK-POOLED-REPRO",
  "release_mixing_authority=layout_preview_only; scientific interpretation prohibited",
  sprintf("spearman_n=%d", nrow(scatter)),
  sprintf("spearman_rho=%.12f", unname(rho$estimate)),
  "final_renderer_gate=promoted_terminal_corrected_coloc"
), file.path(candidate_root, "DESIGN_PREVIEW_ONLY.txt"))
cat("[saved]", out, "\n")
