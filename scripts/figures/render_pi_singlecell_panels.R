# KEY MESSAGE: Donor-level single-cell programs connect stage-associated remodeling to communication and bulk-tissue signatures without redefining the program registry.
# Candidate Figure 4C-E renderer. Sourced by singlecell_module_heatmap.R so the
# established program and communication entrypoint remains authoritative.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

base <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
candidate_root <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
stage_root <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
if (!nzchar(candidate_root) || !dir.exists(candidate_root)) stop("Active candidate root required", call. = FALSE)
source(file.path(base, "scripts/figures/publication_theme.R"))

program_root <- Sys.getenv(
  "PROGRAM_RELEASE_ROOT",
  file.path(base,
    "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
)
registry_path <- file.path(program_root, "program_registry_v2.tsv")
scores_path <- file.path(program_root, "donor_program_scores_primary.tsv")
membership_path <- file.path(program_root, "program_membership_v2.tsv")
registry <- fread(registry_path)
scores <- fread(scores_path)
membership <- fread(membership_path)
if (nrow(registry) != 117L || uniqueN(registry$program_uid) != 117L) stop("Program registry cardinality drift", call. = FALSE)
if (!any(registry$primary_selected)) {
  stop("Corrected registry has no primary-selected programs", call. = FALSE)
}

out_dir <- file.path(candidate_root, "figure4", "panels")
src_dir <- file.path(candidate_root, "source_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

ct_order <- c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells")
ct_labels <- c(hepatocytes = "Hepatocyte", fibroblasts = "Fibroblast",
               macrophages = "Macrophage", cholangiocytes = "Cholangiocyte", tcells = "T cell")
ct_colors <- c(hepatocytes = ct_palette[["Hepatocytes"]], fibroblasts = ct_palette[["Fibroblasts"]],
               macrophages = ct_palette[["Macrophages"]], cholangiocytes = ct_palette[["Cholangiocytes"]],
               tcells = ct_palette[["T cells"]])
registry[, cell_type := factor(cell_type, levels = ct_order)]
registry[, y := as.numeric(cell_type) + seq(-0.28, 0.28, length.out = .N), by = cell_type]
registry[, selected_state := fifelse(primary_selected & hc3_supported, "HC3-supported",
                              fifelse(primary_selected, "Selected", "All programs"))]
registry[, effect_finite := is.finite(primary_beta) & is.finite(primary_se)]
finite_limits <- range(c(registry[effect_finite == TRUE, primary_beta - 1.96 * primary_se],
                         registry[effect_finite == TRUE, primary_beta + 1.96 * primary_se]))
finite_span <- diff(finite_limits)
untestable_x <- finite_limits[[1L]] - 0.12 * finite_span
registry[, plot_beta := fifelse(effect_finite, primary_beta, untestable_x)]
effect_breaks <- pretty(finite_limits, n = 4)

left <- ggplot(registry, aes(plot_beta, y)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_segment(data = registry[effect_finite == TRUE], aes(x = primary_beta - 1.96 * primary_se,
                   xend = primary_beta + 1.96 * primary_se, yend = y,
                   color = cell_type), linewidth = 0.25, alpha = 0.55) +
  geom_point(data = registry[effect_finite == TRUE], aes(fill = cell_type), shape = 21, size = 1.15,
             color = "white", stroke = 0.2) +
  geom_point(data = registry[effect_finite == FALSE], shape = 4, size = 1.15,
             color = "#9E9E9E", stroke = 0.35) +
  geom_point(data = registry[primary_selected & effect_finite], aes(shape = selected_state),
             size = 2.15, fill = NA, color = "black", stroke = 0.5) +
  scale_color_manual(values = ct_colors, guide = "none") +
  scale_fill_manual(values = ct_colors, guide = "none") +
  scale_shape_manual(values = c(`HC3-supported` = 21, Selected = 1), name = NULL) +
  scale_x_continuous(breaks = c(untestable_x, effect_breaks),
                     labels = c("NT", format(effect_breaks, trim = TRUE))) +
  scale_y_continuous(breaks = seq_along(ct_order), labels = ct_labels[ct_order]) +
  labs(x = "Stage effect", y = NULL) + theme_masld_compact() +
  theme(axis.text.y = element_text(size = 6, face = "plain"),
        legend.position = "bottom", legend.key.size = unit(7, "pt"),
        plot.margin = margin(2, 2, 2, 2))

selected_ids <- registry[primary_selected == TRUE, program_uid]
profile <- scores[program_uid %in% selected_ids & exclude_stage_analysis == FALSE &
                    disease_stage_coarse %in% c("Healthy", "Steatosis", "Steatohepatitis")]
profile[, centered_score := as.numeric(scale(score)), by = program_uid]
profile <- profile[, .(
  centered_activity = mean(centered_score),
  SE = sd(centered_score) / sqrt(.N),
  n_donors = uniqueN(donor), n_datasets = uniqueN(dataset)
), by = .(program_uid, disease_stage_coarse)]
profile <- merge(profile, registry[, .(program_uid, module_name, hc3_supported)], by = "program_uid")
profile[, stage := factor(
  disease_stage_coarse,
  levels = coarse_stage_levels[1:3],
  labels = unname(coarse_stage_labels[coarse_stage_levels[1:3]])
)]
profile[, label := sub(" \\([^)]*\\)$", "", module_name)]
profile[, label := factor(label, levels = unique(label))]
right <- ggplot(profile, aes(stage, centered_activity, group = program_uid)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.25) +
  geom_ribbon(aes(ymin = centered_activity - SE, ymax = centered_activity + SE,
                  fill = hc3_supported), alpha = 0.15, color = NA) +
  geom_line(aes(color = hc3_supported), linewidth = 0.65) +
  geom_point(aes(fill = hc3_supported), shape = 21, size = 1.1, color = "white", stroke = 0.2) +
  facet_wrap(~label, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c(`TRUE` = cat_palette[[3]], `FALSE` = "#9E9E9E"), guide = "none") +
  scale_fill_manual(values = c(`TRUE` = cat_palette[[3]], `FALSE` = "#9E9E9E"), guide = "none") +
  labs(x = NULL, y = "Centered activity") + theme_masld_compact() +
  theme(axis.text.x = element_text(size = 6, face = "plain", angle = 30, hjust = 1),
        strip.text = element_text(size = 6, face = "plain"), plot.margin = margin(2, 2, 2, 2))

panel_c <- (left | right) + plot_layout(widths = c(1.15, 1))
ggsave(file.path(out_dir, "fig4c_program_landscape.pdf"), panel_c,
       width = 7.1, height = 2.45, device = cairo_pdf)

# Fixed 13-pair roster. Cells that do not satisfy the roster-specific q gate are
# retained with reduced opacity rather than being reselected.
ccc_path <- Sys.getenv(
  "CCC_TRAJECTORY_FILE",
  file.path(base, "figures/main/fig4_singlecell_programs/source_tables/legacy_from_fig3/ccc_trajectories_data_dc.csv")
)
ccc_audit_path <- Sys.getenv(
  "CCC_AUDIT_FILE",
  file.path(base, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/ccc_dc_before_after_report.tsv")
)
ccc <- fread(ccc_path)
ccc_audit <- fread(ccc_audit_path)[, .(
  headline_label, dc_Estimate, dc_pval, dc_q, dc_n, testable_dc, dc_q05
)]
if (nrow(ccc_audit) != 13L || sum(ccc_audit$dc_q05, na.rm = TRUE) != 0L) {
  stop("Donor-collapsed communication audit drift", call. = FALSE)
}
ccc <- merge(ccc, ccc_audit, by = "headline_label", all.x = TRUE)
ccc <- ccc[disease_stage_coarse %in% c("Healthy", "Steatosis", "Steatohepatitis")]
if (uniqueN(ccc$headline_label) != 13L) stop("Communication roster cardinality drift", call. = FALSE)
ccc[, donor_state := fifelse(dc_q05 == TRUE, "BH-supported",
                      fifelse(testable_dc == TRUE, "Tested, not supported", "Untestable"))]
ccc[, stage := factor(
  disease_stage_coarse,
  levels = coarse_stage_levels[1:3],
  labels = unname(coarse_stage_labels[coarse_stage_levels[1:3]])
)]
ccc[, row_z := as.numeric(scale(mean_score)), by = headline_label]
ccc[, lr_label := gsub("  +", " ", gsub("_", "/", headline_label))]
order_lr <- ccc[, .(
  delta = row_z[stage == "Steatohepatitis"] - row_z[stage == "Healthy"]
), by = lr_label][order(-delta), lr_label]
ccc[, lr_label := factor(lr_label, levels = rev(order_lr))]
panel_d <- ggplot(ccc, aes(stage, lr_label, fill = row_z, alpha = donor_state)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$mash,
                       midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish, name = "L-R z",
                       breaks = c(-1.5, 0, 1.5),
                       guide = guide_colorbar(barwidth = unit(18, "mm"),
                                              barheight = unit(2, "mm"))) +
  scale_alpha_manual(values = c(`BH-supported` = 1, `Tested, not supported` = 0.45,
                                Untestable = 0.22), guide = "none") +
  labs(x = NULL, y = NULL) + theme_masld_compact() +
  theme(axis.text.x = element_text(size = 6, face = "plain", angle = 25, hjust = 1),
        axis.text.y = element_text(size = 6, face = "plain"), legend.position = "bottom",
        legend.title = element_text(size = 6), legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0))
ggsave(file.path(out_dir, "fig4d_communication.pdf"), panel_d,
       width = 2.25, height = 2.45, device = cairo_pdf)

# Bulk projection of every primary-selected program in the supplied release. The display keeps
# HC3 state explicit and reports descriptive member-gene support rather than a
# naive program-level standard error that would ignore gene-gene covariance.
# A full Figure 4E is written only when the same promoted genetics gate used by
# Figure 3F exists.
stage <- fread(file.path(stage_root, "stage_extension_all_gene_results.tsv"))[axis == "fibrosis"]
selected <- registry[primary_selected == TRUE, .(
  program_uid, cell_type, module_name, primary_beta, primary_se,
  primary_qvalue, primary_hc3_qvalue, hc3_supported,
  primary_n_donors, primary_n_datasets
)]
weights <- membership[program_uid %in% selected$program_uid &
                        grepl("^gencode_v49_", mapped_symbol_status) & mapped_symbol != ""]
projection_gene <- merge(weights, stage, by.x = "mapped_symbol", by.y = "gene_name")
projection_gene[, program_effect := weighted.mean(logFC, original_l1_weight, na.rm = TRUE),
                by = .(program_uid, contrast)]
projection <- projection_gene[, .(
  effect = first(program_effect),
  n_genes = uniqueN(mapped_symbol),
  n_bh = uniqueN(mapped_symbol[FDR < 0.05]),
  bh_weight_fraction = weighted.mean(FDR < 0.05, original_l1_weight, na.rm = TRUE),
  direction_agreement = weighted.mean(
    sign(logFC) == sign(first(program_effect)), original_l1_weight, na.rm = TRUE
  )
), by = .(program_uid, contrast)]
projection <- merge(projection, selected, by = "program_uid")
projection[, stage := factor(sub("_vs_F0", "", contrast), levels = paste0("F", 1:4))]
projection[, stage_number := as.integer(stage)]
projection[, series := sub(" \\([^)]*\\)$", "", module_name)]
projection[, series_type := "Program"]
projection[, series_label := sprintf("%s (n=%d)", series, n_genes)]
series_order <- projection[, .(primary_beta = first(primary_beta)), by = series_label][order(primary_beta), series_label]
projection[, series_label := factor(series_label, levels = series_order)]

genetics_root <- Sys.getenv("GENETICS_RELEASE_ROOT", "")
genetics_requested <- nzchar(genetics_root)
promotion_file <- file.path(genetics_root, "PROMOTED.json")
pair_file <- file.path(genetics_root, "figure_signal_pairs.tsv")
manifest_file <- file.path(genetics_root, "output_manifest.tsv")
full_genetics <- nzchar(genetics_root) && all(file.exists(c(promotion_file, pair_file, manifest_file)))
anchor <- data.table(stage = factor(levels = paste0("F", 1:4)), series = character(),
                     series_type = character(), effect = numeric())
if (full_genetics) {
  promoted <- paste(readLines(promotion_file, warn = FALSE), collapse = "")
  full_genetics <- grepl('"status"[[:space:]]*:[[:space:]]*"promoted"', promoted) &&
    grepl('"terminal"[[:space:]]*:[[:space:]]*true', promoted) &&
    grepl('"canonical_promotion_authorized"[[:space:]]*:[[:space:]]*true', promoted)
}
if (full_genetics) {
  genetics_manifest <- fread(manifest_file, colClasses = "character")
  pair_manifest <- genetics_manifest[relative_path == "figure_signal_pairs.tsv"]
  pair_sha <- strsplit(system2("sha256sum", normalizePath(pair_file), stdout = TRUE)[[1L]],
                       "[[:space:]]+")[[1L]][[1L]]
  full_genetics <- nrow(pair_manifest) == 1L && pair_manifest$sha256 == pair_sha &&
    as.numeric(pair_manifest$size_bytes) == file.info(pair_file)$size
}
if (genetics_requested && !full_genetics) {
  stop("Genetics gate closed: Figure 4E requires a promoted terminal checksum-linked signal-pair export",
       call. = FALSE)
}
if (full_genetics) {
  pairs <- fread(pair_file)
  candidates <- unique(pairs[primary_colocalized == TRUE & finalized_figure2_gene == TRUE,
    .(gene_id_versioned, gene_symbol, `PP.H4.susie`)])
  effects <- stage[gene_id_versioned %in% candidates$gene_id_versioned]
  eligible <- effects[, .(finite_all = all(is.finite(logFC)), any_bh = any(FDR < 0.05),
                          max_effect = max(abs(logFC))), by = gene_id_versioned]
  candidates <- merge(candidates, eligible[finite_all & any_bh], by = "gene_id_versioned")
  setorder(candidates, -`PP.H4.susie`, -max_effect, gene_symbol)
  candidates <- candidates[!duplicated(gene_id_versioned)][1:min(2L, .N)]
  anchor <- merge(effects[gene_id_versioned %in% candidates$gene_id_versioned], candidates, by = "gene_id_versioned")
  anchor[, `:=`(stage = factor(sub("_vs_F0", "", contrast), levels = paste0("F", 1:4)),
                 series = gene_symbol, series_type = "Genetic anchor", effect = logFC)]
  stop("Expanded Figure 4E program layout is ready; genetic-anchor row integration remains fail-closed until post-COLOC visual adjudication",
       call. = FALSE)
}
program_rows <- unique(projection[, .(
  series_label, primary_beta, primary_se, primary_hc3_qvalue, hc3_supported
)])
program_rows[, `:=`(
  ci_low = primary_beta - 1.96 * primary_se,
  ci_high = primary_beta + 1.96 * primary_se
)]
p_program <- ggplot(program_rows, aes(primary_beta, series_label)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = series_label),
               linewidth = 0.35, color = "#4D4D4D") +
  geom_point(aes(fill = hc3_supported), shape = 21, size = 1.6,
             color = "black", stroke = 0.4) +
  scale_fill_manual(values = c(`TRUE` = "#2A8C7D", `FALSE` = "white"), guide = "none") +
  labs(x = expression("Single-cell "*beta), y = NULL) +
  theme_masld_compact() +
  theme(axis.text.y = element_text(size = 6, face = "plain"),
        plot.margin = margin(2, 2, 2, 2))

agreement <- projection[stage == "F4", .(series_label, direction_agreement)]
agreement[, stage_number := 4.9]
p_bulk <- ggplot(projection, aes(stage_number, series_label)) +
  geom_point(aes(fill = effect, size = bh_weight_fraction), shape = 21,
             color = "#4D4D4D", stroke = 0.25) +
  geom_text(data = agreement,
            aes(stage_number, series_label,
                label = scales::percent(direction_agreement, accuracy = 1)),
            inherit.aes = FALSE, size = 6 / .pt, family = "Helvetica") +
  scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$mash,
                       midpoint = 0, breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1.0"),
                       name = expression("Bulk "*log[2]*"FC"),
                       guide = guide_colorbar(barwidth = unit(17, "mm"),
                                              barheight = unit(2, "mm"))) +
  scale_size_continuous(range = c(0.7, 3.0), limits = c(0, 1),
                        breaks = c(0.25, 0.75),
                        labels = scales::percent_format(accuracy = 1),
                        name = "BH-supported weight") +
  scale_x_continuous(breaks = c(1:4, 4.9), labels = c(paste0("F", 1:4), "agree"),
                     limits = c(0.55, 5.25)) +
  labs(x = "Fibrosis stage vs F0", y = NULL) +
  theme_masld_compact() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 2, 2, 2))

panel_e <- (p_program | p_bulk) + plot_layout(widths = c(1.15, 1.55))
plot_data <- projection
e_name <- if (full_genetics && nrow(anchor)) "fig4e_bulk_projection.pdf" else "fig4e_bulk_projection_pre_genetics.pdf"
ggsave(file.path(out_dir, e_name), panel_e, width = 4.15, height = 2.45,
       device = cairo_pdf)

fwrite(registry[, .(program_uid, cell_type, module, module_name, primary_beta, primary_se,
                    primary_qvalue, primary_hc3_qvalue, primary_selected, hc3_supported,
                    primary_n_donors, primary_n_datasets, registry_state)],
       file.path(src_dir, "fig4c_all_117_programs.tsv"), sep = "\t")
fwrite(profile, file.path(src_dir, "fig4c_selected_program_profiles.tsv"), sep = "\t")
fwrite(ccc, file.path(src_dir, "fig4d_fixed_13_pair_communication.tsv"), sep = "\t")
fwrite(plot_data[, .(
  program_uid, cell_type, module_name, series, stage, contrast,
  primary_beta, primary_se, primary_qvalue, primary_hc3_qvalue,
  hc3_supported, primary_n_donors, primary_n_datasets,
  effect, n_genes, n_bh, bh_weight_fraction, direction_agreement
)], file.path(src_dir, "fig4e_bulk_projection.tsv"), sep = "\t")
cat("[saved] Figure 4C-D and", e_name, "\n")
