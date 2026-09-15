#!/usr/bin/env Rscript

# KEY MESSAGE: Four externally anchored systems resolve distinct questions about
# the single-cell gene-program registry; each panel shows the system's calibrated estimand
# without converting the systems into a common score.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(project_root, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)

systems_root <- file.path(workstream_root, "systems")
output_root <- Sys.getenv(
  "MASLD_SYSTEM_FIGURE_OUT",
  file.path(project_root, "figures/candidates/program-systems-2026-08-12-v3")
)
if (dir.exists(output_root) || file.exists(output_root)) {
  stop("Refusing to overwrite existing output: ", output_root, call. = FALSE)
}
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(output_root, ".tmp-", Sys.getpid())
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

lineage_levels <- c(
  "hepatocytes", "cholangiocytes", "fibroblasts", "macrophages", "tcells"
)
lineage_labels <- c(
  hepatocytes = "Hepatocyte", cholangiocytes = "Cholangiocyte",
  fibroblasts = "Fibroblast", macrophages = "Macrophage", tcells = "T cell"
)
lineage_colors <- c(
  hepatocytes = "#C9265E", cholangiocytes = "#1565C0",
  fibroblasts = "#00695C", macrophages = "#E08214", tcells = "#7B1FA2"
)

write_source <- function(x, name) {
  fwrite(x, file.path(tmp, paste0(name, "_source.tsv")), sep = "\t", na = "NA")
}

# ---------------------------------------------------------------------- T3 --
# Each vocabulary is compared only with random sets matched to itself. The
# violin preserves all 200 draws; a diamond marks the observed hit rate.
t3_hit_path <- file.path(systems_root, "T3_vocabulary/hit_rates.tsv")
t3_draw_path <- file.path(systems_root, "T3_vocabulary/matched_random_counts.tsv")
t3_hit <- fread(t3_hit_path)[axis == "fibrosis"]
t3_draw <- fread(t3_draw_path)[axis == "fibrosis"]

vocabulary_labels <- c(
  hotspot_117 = "Hotspot 117", hallmark_50 = "Hallmark 50", cnmf_16 = "cNMF 16",
  bnmf_6 = "bulk NMF 6", wgcna_4 = "WGCNA 4", published_panels = "Published panels"
)
vocabulary_order <- c(
  "Published panels", "WGCNA 4", "bulk NMF 6",
  "cNMF 16", "Hallmark 50", "Hotspot 117"
)
t3_hit[, label := factor(vocabulary_labels[vocabulary], levels = vocabulary_order)]
t3_hit[, provenance := fifelse(
  substrate_leak %in% TRUE | reportable_as_support %in% FALSE,
  "Substrate-leaked", "Independent"
)]
t3_draw[, hit_rate := n_supported / n_tested]
t3_draw[, label := factor(vocabulary_labels[vocabulary], levels = vocabulary_order)]

p_t3 <- ggplot(t3_draw, aes(x = hit_rate, y = label)) +
  geom_violin(
    fill = "#E0E0E0", color = NA, scale = "width", width = 0.72,
    trim = TRUE, orientation = "y"
  ) +
  geom_boxplot(
    width = 0.10, outlier.shape = NA, fill = "white", color = "#616161",
    linewidth = STROKE_HAIRLINE, orientation = "y"
  ) +
  geom_point(
    data = t3_hit,
    aes(x = hit_rate, y = label, color = provenance, shape = provenance),
    inherit.aes = FALSE, size = 2.0, stroke = 0.55
  ) +
  scale_color_manual(
    values = c(Independent = "#C9265E", `Substrate-leaked` = "#9E9E9E"),
    breaks = c("Independent", "Substrate-leaked"), name = NULL
  ) +
  scale_shape_manual(
    values = c(Independent = 18, `Substrate-leaked` = 4),
    breaks = c("Independent", "Substrate-leaked"), name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 0.82), breaks = seq(0, 0.8, 0.2),
    labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "Fibrosis-supported fraction", y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(
    legend.position = "top", legend.direction = "horizontal",
    axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )

t3_source <- rbindlist(list(
  t3_hit[, .(
    row_type = "observed", vocabulary, label = as.character(label),
    hit_rate, n_testable, n_supported, enrichment_over_random,
    p_vs_matched_random, substrate_leak, reportable_as_support
  )],
  t3_draw[, .(
    row_type = "matched_random", vocabulary, label = as.character(label),
    hit_rate, n_testable = n_tested, n_supported,
    enrichment_over_random = NA_real_, p_vs_matched_random = NA_real_,
    substrate_leak = NA, reportable_as_support = NA
  )]
), use.names = TRUE, fill = TRUE)
write_source(t3_source, "panel_t3_vocabulary")

# ---------------------------------------------------------------------- T4 --
# Sparse symmetric matrix of calibrated edges. Residual correlation uses a
# clipped diverging scale; opacity identifies edges reproducing in all cohorts.
t4_edge_path <- file.path(systems_root, "T4_interaction/coordinated_edges.tsv")
t4_comm_path <- file.path(systems_root, "T4_interaction/program_communities.tsv")
t4_summary_path <- file.path(systems_root, "T4_interaction/interaction_summary.tsv")
t4_edges <- fread(t4_edge_path)
t4_comm <- fread(t4_comm_path)
t4_summary <- fread(t4_summary_path)
registry <- fread(program_registry)[, .(
  feature_id = program_uid, cell_type, module_name
)]
t4_order <- merge(t4_comm, registry, by = "feature_id", all.x = TRUE)
t4_order[, cell_type := factor(cell_type, levels = lineage_levels)]
setorder(t4_order, community, cell_type, module_name, feature_id)
t4_order[, display_index := .I]

idx_a <- t4_order[, .(program_a = feature_id, index_a = display_index)]
idx_b <- t4_order[, .(program_b = feature_id, index_b = display_index)]
t4_edges <- merge(t4_edges, idx_a, by = "program_a", all.x = TRUE)
t4_edges <- merge(t4_edges, idx_b, by = "program_b", all.x = TRUE)
stopifnot(!anyNA(t4_edges$index_a), !anyNA(t4_edges$index_b))
t4_edges[, all_four_cohorts := n_same_sign_partial == 4L]
t4_display <- rbindlist(list(
  t4_edges[, .(x = index_a, y = index_b, r_partial_k3, all_four_cohorts)],
  t4_edges[, .(x = index_b, y = index_a, r_partial_k3, all_four_cohorts)]
))

clip_lim <- max(abs(quantile(
  t4_display$r_partial_k3, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE
)))
t4_bounds <- t4_order[, .(
  xmin = min(display_index) - 0.5, xmax = max(display_index) + 0.5,
  ymin = min(display_index) - 0.5, ymax = max(display_index) + 0.5,
  midpoint = mean(range(display_index))
), by = community]
t4_order[, cell_type_character := as.character(cell_type)]

p_t4 <- ggplot(t4_display, aes(x = x, y = y)) +
  geom_tile(
    aes(fill = r_partial_k3, alpha = all_four_cohorts),
    width = 0.96, height = 0.96, color = NA
  ) +
  geom_rect(
    data = t4_bounds,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = NA, color = "black", linewidth = STROKE_HAIRLINE
  ) +
  geom_point(
    data = t4_order,
    aes(x = display_index, y = 0.3, color = cell_type_character),
    inherit.aes = FALSE, shape = 15, size = 1.35
  ) +
  geom_point(
    data = t4_order,
    aes(x = 0.3, y = display_index, color = cell_type_character),
    inherit.aes = FALSE, shape = 15, size = 1.35
  ) +
  geom_text(
    data = t4_bounds,
    aes(x = midpoint, y = max(t4_order$display_index) + 2.5,
        label = paste0("C", community)),
    inherit.aes = FALSE, size = PUB_GEOM_TEXT, color = "black"
  ) +
  scale_fill_gradient2(
    low = "#1565C0", mid = "white", high = "#C9265E", midpoint = 0,
    limits = c(-clip_lim, clip_lim), oob = squish,
    name = "Residual r"
  ) +
  scale_alpha_manual(
    values = c(`FALSE` = 0.28, `TRUE` = 1),
    labels = c(`FALSE` = "1–3 cohorts", `TRUE` = "All 4 cohorts"),
    name = "Same direction"
  ) +
  scale_color_manual(
    values = lineage_colors, labels = lineage_labels,
    breaks = lineage_levels, name = NULL
  ) +
  coord_fixed(
    xlim = c(-0.4, max(t4_order$display_index) + 0.5),
    ylim = c(-0.4, max(t4_order$display_index) + 4.5),
    expand = FALSE, clip = "off"
  ) +
  labs(x = "Single-cell gene programs", y = "Single-cell gene programs") +
  theme_masld(base_size = 6) + theme_pub() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    plot.margin = margin(3, 3, 3, 3)
  )

write_source(
  t4_edges[, .(
    program_a, program_b, cell_type_a, cell_type_b, module_name_a, module_name_b,
    r_partial_k3, n_same_sign_partial, n_cohorts, q_partial_k3,
    q_partial_empirical, q_partial_group_null, heterogeneity_I2,
    all_four_cohorts
  )],
  "panel_t4_coordination_edges"
)
write_source(
  t4_order[, .(
    display_index, feature_id, community,
    cell_type = cell_type_character, module_name
  )],
  "panel_t4_program_order"
)

# ---------------------------------------------------------------------- T5 --
# Only the prespecified 27 fibrosis-supported programs are shown. The x-axis is
# the corrected composition-coupled fraction. Full bootstrap CIs remain in the
# source table: six extend below zero and would compress the main visual scale.
t5_path <- file.path(systems_root, "T5_composition/program_composition_decomposition.tsv")
t5 <- fread(t5_path)[fibrosis_supported_linear %in% TRUE]
stopifnot(nrow(t5) == 27L)
setorder(t5, coupling_fraction_corrected, unique_fibrosis)
t5[, display_rank := .I]
t5[, cell_type := factor(cell_type, levels = lineage_levels)]
t5[, robust_label := fifelse(robust_display %in% TRUE, module_name, "")]

p_t5 <- ggplot(t5, aes(x = coupling_fraction_corrected, y = display_rank)) +
  geom_vline(
    xintercept = 0.5, color = "#616161", linetype = "dashed",
    linewidth = STROKE_HAIRLINE
  ) +
  geom_point(
    aes(fill = cell_type), shape = 21, size = 1.9,
    stroke = 0.35, color = "black"
  ) +
  geom_text(
    data = t5[robust_display %in% TRUE],
    aes(label = robust_label), hjust = -0.10,
    size = PUB_GEOM_TEXT, color = "black", fontface = "plain"
  ) +
  scale_fill_manual(
    values = lineage_colors, labels = lineage_labels,
    breaks = lineage_levels, name = NULL
  ) +
  scale_x_continuous(
    limits = c(-0.10, 1.16), breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(0.8, 0.8))) +
  labs(x = "Composition-coupled fraction", y = "27 fibrosis-associated programs") +
  theme_masld(base_size = 6) + theme_pub() +
  theme(
    legend.position = "top", legend.direction = "horizontal",
    axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )

write_source(
  t5[, .(
    display_rank, feature_id, cell_type = as.character(cell_type), module_name,
    robust_display, coupling_class, coupling_fraction_corrected,
    coupling_fraction_corrected_ci_lower,
    coupling_fraction_corrected_ci_upper, unique_fibrosis,
    q_unique_fibrosis, beta_retained_fraction
  )],
  "panel_t5_composition_coupling"
)

# ---------------------------------------------------------------------- T6 --
# Four calibrated sex-modifier views. Deterministic within-row jitter avoids
# overplotting. Status distinguishes an informative null from indeterminacy.
t6_path <- file.path(systems_root, "T6_modifiers/modifier_interaction_meta.tsv")
t6 <- fread(t6_path)[view %chin% c(
  "sex_x_fibrosis", "sex_x_nas",
  "sex_recorded_only_x_fibrosis", "sex_recorded_only_x_nas"
)]
view_labels <- c(
  sex_x_fibrosis = "All donors · fibrosis",
  sex_x_nas = "All donors · NAS",
  sex_recorded_only_x_fibrosis = "Recorded sex · fibrosis",
  sex_recorded_only_x_nas = "Recorded sex · NAS"
)
view_order <- rev(unname(view_labels))
t6[, view_label := factor(view_labels[view], levels = view_order)]
t6[, status := fifelse(support == "informative_null", "Informative null", "Indeterminate")]
setorder(t6, view_label, beta_meta, feature_id)
t6[, jitter_y := seq(-0.16, 0.16, length.out = .N), by = view_label]
t6[, display_y := as.numeric(view_label) + jitter_y]

p_t6 <- ggplot(t6, aes(x = beta_meta, y = display_y)) +
  geom_vline(xintercept = 0, color = "black", linewidth = STROKE_HAIRLINE) +
  geom_vline(
    xintercept = c(-0.2, 0.2), color = "#9E9E9E", linetype = "dashed",
    linewidth = STROKE_HAIRLINE
  ) +
  geom_point(
    aes(color = status, shape = status), size = 1.15,
    alpha = 0.82, stroke = 0.35
  ) +
  scale_color_manual(
    values = c(`Informative null` = "#00695C", Indeterminate = "#9E9E9E"),
    breaks = c("Informative null", "Indeterminate"), name = NULL
  ) +
  scale_shape_manual(
    values = c(`Informative null` = 16, Indeterminate = 1),
    breaks = c("Informative null", "Indeterminate"), name = NULL
  ) +
  scale_y_continuous(
    breaks = seq_along(view_order), labels = view_order,
    expand = expansion(add = c(0.35, 0.35))
  ) +
  scale_x_continuous(
    limits = c(-0.38, 0.38), breaks = seq(-0.3, 0.3, 0.1),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Sex × stage interaction (program-score SD)",
    y = NULL
  ) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(
    legend.position = "top", legend.direction = "horizontal",
    axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )

write_source(
  t6[, .(
    feature_id, view, view_label = as.character(view_label), modifier,
    beta_meta, se_meta, ci_lower, ci_upper, p_value_meta, q_value,
    mde_nominal, mde_family, effect_threshold, status,
    n_cohorts, n_same_direction
  )],
  "panel_t6_sex_modification"
)

# ------------------------------------------------------------------- output --
panels <- list(
  fig_system_t3_vocabulary = list(plot = p_t3, width = 3.55, height = 2.35),
  fig_system_t4_coordination = list(plot = p_t4, width = 5.15, height = 4.15),
  fig_system_t5_composition = list(plot = p_t5, width = 4.30, height = 3.15),
  fig_system_t6_modifiers = list(plot = p_t6, width = 4.25, height = 2.35)
)
for (name in names(panels)) {
  ggsave(
    file.path(tmp, paste0(name, ".pdf")), panels[[name]]$plot,
    width = panels[[name]]$width, height = panels[[name]]$height,
    units = "in", device = cairo_pdf
  )
}

captions <- c(
  paste0(
    "T3 — Vocabulary invariance. Observed fibrosis-supported fractions (diamonds) ",
    "are shown against 200 size- and expression-matched random draws for each ",
    "vocabulary. Hotspot 117 does not exceed its own matched null; Hallmark and ",
    "cNMF do. Gray open diamonds mark substrate-leaked vocabularies and are not ",
    "counted as support."
  ),
  paste0(
    "T4 — Residual program coordination. Calibrated pairs among 113 cell-type-specific ",
    "gene programs are ordered into four communities after saturated histology adjustment ",
    "and removal of three global components. Color is residual correlation; opacity ",
    "distinguishes the 1,658 pairs reproducing in all four cohorts from other ",
    "coordinated pairs. C3 is exclusively cholangiocyte and fibroblast programs."
  ),
  paste0(
    "T5 — Composition coupling. Corrected composition-coupled fractions are shown ",
    "for all 27 fibrosis-associated programs; 20 exceed 50%. Full bootstrap intervals ",
    "are retained in the source table rather than compressed into the main panel. ",
    "Because composition and program scores come from the same RNA matrix, the shared ",
    "fraction is an upper bound on attribution, not causal attribution."
  ),
  paste0(
    "T6 — Effect modification. Program-level sex-by-stage interaction estimates are ",
    "shown for fibrosis and NAS in the full and recorded-sex-only arms. No program ",
    "passes the calibrated family; teal marks informative nulls and gray open marks ",
    "remain indeterminate. Dashed lines show the prespecified ±0.20 effect threshold."
  )
)
writeLines(captions, file.path(tmp, "figure_captions.txt"))
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

inputs <- c(
  t3_hit_path, t3_draw_path,
  t4_edge_path, t4_comm_path, t4_summary_path, program_registry,
  t5_path, t6_path
)
input_manifest <- data.table(
  path = sub(paste0("^", project_root, "/"), "", inputs),
  sha256 = vapply(inputs, sha256_file, character(1))
)
fwrite(input_manifest, file.path(tmp, "input_manifest.tsv"), sep = "\t")

readme <- c(
  "# Resolved program systems figure candidate",
  "",
  "Status: candidate only; not adopted.",
  "Source: recorded BULK-PROGRAM-MAP-v9 outputs.",
  "The 117 cell-type-specific gene programs were defined before cross-cohort testing.",
  "Panels: T3 vocabulary, T4 coordination, T5 composition coupling, T6 modifiers.",
  "Open release dependency: axis map 113/27 versus T3 114/26 remains unreconciled.",
  "No holdout data were read and no system was refit by this figure generator."
)
writeLines(readme, file.path(tmp, "README.md"))

artifacts <- setdiff(list.files(tmp), "output_manifest.tsv")
output_manifest <- data.table(
  artifact = artifacts,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
  size_bytes = file.info(file.path(tmp, artifacts))$size
)
fwrite(output_manifest, file.path(tmp, "output_manifest.tsv"), sep = "\t")

if (!file.rename(tmp, output_root)) {
  stop("Could not publish candidate directory: ", output_root, call. = FALSE)
}
message("RESOLVED_SYSTEM_FIGURES_COMPLETE ", output_root)
