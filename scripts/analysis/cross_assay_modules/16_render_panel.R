#!/usr/bin/env Rscript
# 16: the candidate Figure 5 panels.
#
# The panel shows two decisions per cell: does the standardized association
# replicate in the discovery direction (tile fill), and does it exceed
# co-expressed sets of the same size (dot). It carries no combined score, no
# rank and no count of supporting modalities. Every explanatory sentence lives
# in the caption file next to the PDFs, not inside the plotting area.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))
suppressPackageStartupMessages({ library(ggplot2) })
grDevices::pdf.options(useDingbats = FALSE)

contract <- cam_contract()
panel_dir <- file.path(cam_out_root(), "panel")
FIG_VERSION <- Sys.getenv("CAM_FIGURE_VERSION", unset = "figure")
out <- cam_dir(FIG_VERSION)
cam_assert(file.exists(file.path(panel_dir, "READY")), "Run 15 first")

FS <- 6
FAM <- "Helvetica"
INK <- "#222222"
CONTROL_GREY <- "#9E9E9E"
STATE_COLS <- c(
  replicates = "#B2182B", discordant = "#2166AC", tested_negative = "#4D4D4D",
  indeterminate = CONTROL_GREY, untestable = "#E8E8E8", not_applicable = "#EFE6D8",
  coherent = "#F4A582", coherence_unresolved = "#BABABA", source_dependent = "#C79A3E"
)
STATE_LABELS <- c(
  replicates = "association replicates", discordant = "association reversed",
  tested_negative = "tested negative (bound < SESOI)", indeterminate = "not resolved",
  untestable = "untestable (coverage)", not_applicable = "no endpoint of this family",
  coherent = "spatial excess over matched genes", coherence_unresolved = "spatial excess unresolved",
  source_dependent = "source dependent"
)
ASSAY_SHORT <- c(
  bulk_gse268273 = "GSE268273", bulk_gse276114 = "GSE276114",
  proteome_liver_pxd051911 = "Liver proteome", snrna_all_cells = "snRNA (diagnosis)",
  h3k27ac_gse267145 = "H3K27ac", atac_gse296875 = "ATAC"
)
ASSAY_LABEL <- c(
  bulk_gse268273 = "Bulk RNA (GSE268273)",
  bulk_gse276114 = "Bulk RNA (GSE276114)",
  proteome_liver_pxd051911 = "Liver proteome",
  snrna_all_cells = "snRNA (diagnosis proxy)",
  h3k27ac_gse267145 = "H3K27ac promoter",
  atac_gse296875 = "ATAC promoter",
  visium_gse192741 = "Visium (donor)",
  visium_vu = "Visium (array)",
  cosmx_govaere = "CosMx panel",
  proteome_plasma_pxd051911 = "Plasma proteome"
)

base_theme <- theme_minimal(base_size = FS, base_family = FAM) +
  theme(
    text = element_text(size = FS, colour = INK),
    axis.text = element_text(size = FS, colour = INK),
    axis.title = element_text(size = FS, colour = INK),
    legend.text = element_text(size = FS, colour = INK),
    legend.title = element_text(size = FS, colour = INK),
    legend.key.size = unit(3, "mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15, colour = "#EEEEEE"),
    plot.title = element_blank(), plot.subtitle = element_blank()
  )

panel <- fread(file.path(panel_dir, "module_by_assay_states.tsv"))
paired <- if (file.exists(file.path(panel_dir, "paired_transfer_panel.tsv")))
  fread(file.path(panel_dir, "paired_transfer_panel.tsv")) else NULL
examples <- fread(file.path(panel_dir, "prespecified_example_modules.tsv"))

panel[, assay_label := factor(ASSAY_LABEL[assay], levels = unname(ASSAY_LABEL))]
panel[, display_state := factor(display_state, levels = names(STATE_COLS))]
panel[, spec_dot := specificity_state == "exceeds_coexpressed_sets"]

# Row order is fixed before any state is read: by discovery effect size.
order_tab <- unique(panel[, .(module_id, discovery_supported, discovery_effect,
                              score_stability_flag, n_genes, endpoint_family_label)])
setorder(order_tab, -discovery_supported, -discovery_effect)
shown <- order_tab[discovery_supported == TRUE & score_stability_flag == "score_stable"]
cam_say("main panel rows: ", nrow(shown), "; supplementary rows: ", nrow(order_tab))

render_state_grid <- function(dt, module_order, file, width, height, stability_strip = FALSE) {
  d <- dt[module_id %in% module_order]
  d[, module_id := factor(module_id, levels = rev(module_order))]
  present <- names(STATE_COLS)[names(STATE_COLS) %in% unique(as.character(d$display_state))]
  d[, display_state := factor(as.character(display_state), levels = present)]
  p <- ggplot(d, aes(x = assay_label, y = module_id)) +
    geom_tile(aes(fill = display_state), colour = "white", linewidth = 0.12) +
    geom_point(data = d[spec_dot == TRUE], aes(x = assay_label, y = module_id),
               size = 0.35, colour = INK, inherit.aes = FALSE) +
    scale_fill_manual(values = STATE_COLS[present], labels = STATE_LABELS[present],
                      drop = FALSE, name = NULL) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    base_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 0, vjust = 0, size = FS),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(), legend.position = "right",
          plot.margin = margin(t = 14, r = 2, b = 2, l = 2, unit = "mm"))
  if (stability_strip) {
    st <- unique(d[, .(module_id, score_stability_flag)])
    st[, assay_label := factor("score unstable", levels = levels(d$assay_label))]
    p <- p + geom_point(data = st[score_stability_flag == "score_unstable"],
                        aes(x = assay_label, y = module_id), shape = 4, size = 0.6,
                        colour = INK, inherit.aes = FALSE)
  }
  ggsave(file, p, width = width, height = height, units = "in", device = cairo_pdf)
  invisible(p)
}

# 5G: discovery-supported, score-stable modules, by assay.
render_state_grid(panel, shown$module_id,
                  file.path(out, "fig5g_cross_assay_module_states.pdf"), width = 5.0, height = 3.4)

# S5: the whole family, with the score-stability marker as an extra column.
panel_s5 <- copy(panel)
panel_s5[, assay_label := factor(as.character(assay_label), levels = c(unname(ASSAY_LABEL), "score unstable"))]
render_state_grid(panel_s5, order_tab$module_id,
                  file.path(out, "figS5_cross_assay_module_states_all.pdf"), width = 5.2, height = 4.6,
                  stability_strip = TRUE)

# 5H: association q against competitive q, per assay.
cmp <- panel[evidence_axis == "association" & applicable == TRUE & testable == TRUE &
               is.finite(association_q) & is.finite(competitive_q)]
if (nrow(cmp)) {
  cmp[, assay_short := factor(ASSAY_SHORT[assay], levels = unname(ASSAY_SHORT))]
  cmp[, neglog_q := -log10(pmax(association_q, 1e-12))]
  cmp[, neglog_comp := -log10(pmax(competitive_q, 1e-4))]
  present2 <- names(STATE_COLS)[names(STATE_COLS) %in% unique(as.character(cmp$association_state))]
  cmp[, association_state := factor(as.character(association_state), levels = present2)]
  p2 <- ggplot(cmp, aes(x = neglog_q, y = neglog_comp)) +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.2, colour = CONTROL_GREY) +
    geom_vline(xintercept = -log10(0.05), linewidth = 0.2, colour = CONTROL_GREY) +
    geom_point(aes(colour = association_state), size = 0.6) +
    scale_colour_manual(values = STATE_COLS[present2], labels = STATE_LABELS[present2],
                        drop = FALSE, name = NULL) +
    facet_wrap(~ assay_short, nrow = 1) +
    labs(x = "-log10 association BH q (two-sided)", y = "-log10 competitive BH q") +
    base_theme + theme(strip.text = element_text(size = FS, colour = INK),
                       panel.spacing = unit(2.5, "mm"))
  ggsave(file.path(out, "fig5h_association_vs_competitive.pdf"), p2,
         width = 6.4, height = 1.9, units = "in", device = cairo_pdf)
}

# 5I: within-participant transfer, sign shown, adjustment set named per facet.
if (!is.null(paired) && nrow(paired)) {
  pd <- paired[testable_both == TRUE & is.finite(partial_spearman)]
  pd[, pairing_label := fifelse(pairing == "gse267145_rna_vs_h3k27ac",
                                "RNA vs H3K27ac, given histology + sex (99 participants)",
                                "RNA vs ATAC, given well (39 donors)")]
  pd[, sig := !is.na(q_value) & q_value < 0.05]
  pd[, coupling_sign := factor(coupling_sign, levels = c("positive", "inverse"))]
  pd[, lo := fifelse(is.finite(boot_ci_low), boot_ci_low, lowo_min)]
  pd[, hi := fifelse(is.finite(boot_ci_high), boot_ci_high, lowo_max)]
  pd[, ord := seq_len(.N), by = pairing]
  pd <- pd[order(pairing, partial_spearman)][, ord := seq_len(.N), by = pairing]
  p3 <- ggplot(pd, aes(x = partial_spearman, y = ord)) +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = CONTROL_GREY) +
    geom_segment(aes(x = lo, xend = hi, yend = ord), linewidth = 0.15, colour = "#CCCCCC") +
    geom_point(aes(colour = coupling_sign, shape = sig), size = 0.7, stroke = 0.3) +
    scale_colour_manual(values = c(positive = "#B2182B", inverse = "#2166AC"), name = "coupling sign") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = "BH q<0.05", `FALSE` = "not significant"), name = NULL) +
    facet_wrap(~ pairing_label, nrow = 1, scales = "free") +
    labs(x = "partial Spearman (RNA score vs chromatin score)", y = NULL) +
    base_theme +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
          strip.text = element_text(size = FS, colour = INK), panel.spacing = unit(4, "mm"))
  ggsave(file.path(out, "fig5i_paired_within_participant_transfer.pdf"), p3,
         width = 5.2, height = 2.4, units = "in", device = cairo_pdf)
}

# 5J: coverage.
cov <- unique(panel[, .(assay_label, module_id, fraction_measured)])
cam_assert(!any(is.na(cov$fraction_measured)), "5J would drop a measured module")
cov[, assay_label := factor(as.character(assay_label), levels = unname(ASSAY_LABEL))]
p4 <- ggplot(cov, aes(y = assay_label, x = fraction_measured)) +
  geom_vline(xintercept = contract$scoring$testability$min_fraction_measured,
             linewidth = 0.2, colour = CONTROL_GREY) +
  geom_jitter(height = 0.18, width = 0, size = 0.3, alpha = 0.45, colour = INK) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "fraction of module members measured", y = NULL) +
  base_theme
ggsave(file.path(out, "fig5j_module_member_coverage.pdf"), p4,
       width = 4.2, height = 2.3, units = "in", device = cairo_pdf)

# 5K: six prespecified modules, oriented standardized slope with 95% CI per assay.
fk <- panel[evidence_axis == "association" & module_id %in% examples$module_id]
fk[, ci_lo := effect_oriented - stats::qt(0.975, df) * effect_se]
fk[, ci_hi := effect_oriented + stats::qt(0.975, df) * effect_se]
fk[, module_label := paste0(module_id, " (", endpoint_family_label,
                            fifelse(score_stability_flag == "score_unstable", "; unstable", ""), ")")]
fk[, module_label := factor(module_label, levels = unique(module_label[order(endpoint_family_label, -discovery_effect)]))]
fk[, assay_label := factor(as.character(assay_label), levels = rev(unname(ASSAY_LABEL)))]
present3 <- names(STATE_COLS)[names(STATE_COLS) %in% unique(as.character(fk$association_state))]
fk[, association_state := factor(as.character(association_state), levels = present3)]
p5 <- ggplot(fk[applicable == TRUE], aes(y = assay_label, x = effect_oriented)) +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = CONTROL_GREY) +
  geom_vline(xintercept = contract$negatives$sesoi_standardized_slope, linewidth = 0.2,
             colour = CONTROL_GREY, linetype = "dashed") +
  geom_segment(aes(x = ci_lo, xend = ci_hi, yend = assay_label), linewidth = 0.25, colour = INK) +
  geom_point(aes(colour = association_state), size = 0.9) +
  scale_colour_manual(values = STATE_COLS[present3], labels = STATE_LABELS[present3], drop = FALSE, name = NULL) +
  facet_wrap(~ module_label, nrow = 2) +
  labs(x = "oriented standardized slope (SD score per SD endpoint), 95% CI", y = NULL) +
  base_theme + theme(strip.text = element_text(size = FS, colour = INK), panel.spacing = unit(3, "mm"))
ggsave(file.path(out, "fig5k_prespecified_module_forests.pdf"), p5,
       width = 6.4, height = 3.2, units = "in", device = cairo_pdf)

# Source tables.
cam_write_tsv(panel[module_id %in% shown$module_id], file.path(out, "fig5g_source.tsv"))
cam_write_tsv(panel, file.path(out, "figS5_source.tsv"))
if (nrow(cmp)) cam_write_tsv(cmp, file.path(out, "fig5h_source.tsv"))
if (!is.null(paired)) cam_write_tsv(paired, file.path(out, "fig5i_source.tsv"))
cam_write_tsv(cov, file.path(out, "fig5j_source.tsv"))
cam_write_tsv(fk, file.path(out, "fig5k_source.tsv"))

# The caption authority travels with the panels, so a reader of the candidate
# directory never has a figure without the sentences that qualify it, and so
# the validator can check the caption text against the withdrawn phrases.
file.copy(file.path(cam_script_dir(), "CAPTIONS.md"), file.path(out, "CAPTIONS.md"),
          overwrite = TRUE)

pdfs <- list.files(out, pattern = "\\.pdf$", full.names = TRUE)
cam_write_tsv(data.table(file = basename(pdfs), sha256 = vapply(pdfs, cam_sha256, character(1))),
              file.path(out, "figure_checksums.tsv"))
cam_session_info(file.path(cam_out_root(), "sessionInfo.txt"))
writeLines("figures rendered", file.path(out, "READY"))
cam_say("16 complete: ", length(pdfs), " PDFs")
