#!/usr/bin/env Rscript
# KEY MESSAGE: Complete pathway families align fibrosis-stage remodeling with stage-adjusted molecular-continuum effects while retaining unsupported and unestimable sets in the declared family.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 65_render_pathway_effect_overviews.R COMPARISON_SOURCE OUTPUT_CANDIDATE",
       call. = FALSE)
}
source_file <- normalizePath(args[[1L]], mustWork = TRUE)
output_candidate <- args[[2L]]
if (file.exists(output_candidate)) {
  stop("Refusing to overwrite output candidate: ", output_candidate, call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
contract <- ml_read_contract()

panel_dir <- file.path(output_candidate, "panels")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(panel_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

family_order <- c("hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc")
family_labels <- c(
  hallmark = "Hallmark", kegg = "KEGG", reactome = "Reactome",
  go_bp = "GO biological process", go_mf = "GO molecular function",
  go_cc = "GO cellular component"
)
expected_sizes <- unlist(contract$pathway_collections[family_order], use.names = TRUE)
support_levels <- c("Neither", "Stage only", "Continuum only", "Both")
support_colors <- c(
  "Neither" = "#9E9E9E", "Stage only" = "#F4A674",
  "Continuum only" = "#7B1FA2", "Both" = "#00695C"
)
support_shapes <- c("Neither" = 1, "Stage only" = 15, "Continuum only" = 17, "Both" = 16)

theme_effect <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, color = "black", face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.title = element_text(size = 6, face = "plain"),
      strip.text = element_text(size = 6, face = "plain"),
      strip.background = element_blank(), plot.title = element_blank(),
      plot.subtitle = element_text(size = 6, face = "plain", hjust = 0),
      plot.caption = element_text(size = 6, face = "plain", hjust = 0),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(2.8, "mm"),
      panel.spacing = grid::unit(1.6, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_panel <- function(plot, filename, width, height) {
  path <- file.path(panel_dir, filename)
  ml_assert(!file.exists(path), paste0("Refusing to overwrite: ", path))
  ggsave(path, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  ml_assert(file.info(path)$size > 1000L, paste0("Unexpectedly small PDF: ", path))
  path
}

all_results <- fread(source_file, na.strings = c("", "NA"))
required <- c(
  "set_id", "beta", "q_value", "stage_supported", "continuum_beta",
  "continuum_q", "continuum_supported", "support_class", "collection",
  "family_size", "bh_family_size"
)
ml_assert(all(required %in% names(all_results)), "Pathway comparison schema drift")
ml_assert(setequal(unique(all_results$collection), family_order),
          "Pathway comparison collections drifted")
all_results[, collection := factor(collection, levels = family_order)]
all_results[, support_class := factor(support_class, levels = support_levels)]
all_results[, coordinate_estimable := is.finite(beta) & is.finite(continuum_beta)]

for (collection_id in family_order) {
  observed <- all_results[collection == collection_id]
  ml_assert(nrow(observed) == expected_sizes[[collection_id]],
            paste0("Complete family size drift: ", collection_id))
  ml_assert(uniqueN(observed$set_id) == nrow(observed),
            paste0("Duplicate pathway IDs: ", collection_id))
  ml_assert(all(observed$family_size == expected_sizes[[collection_id]]),
            paste0("Declared family size drift: ", collection_id))
  ml_assert(all(observed$bh_family_size[is.finite(observed$bh_family_size)] ==
                  expected_sizes[[collection_id]]),
            paste0("BH family size drift: ", collection_id))
}

family_summary <- all_results[, .(
  complete_family = .N,
  coordinate_estimable = sum(coordinate_estimable),
  coordinate_unestimable = sum(!coordinate_estimable),
  stage_supported = sum(stage_supported %in% TRUE, na.rm = TRUE),
  continuum_supported = sum(continuum_supported %in% TRUE, na.rm = TRUE),
  both = sum(support_class == "Both", na.rm = TRUE),
  stage_only = sum(support_class == "Stage only", na.rm = TRUE),
  continuum_only = sum(support_class == "Continuum only", na.rm = TRUE),
  neither = sum(support_class == "Neither", na.rm = TRUE),
  spearman_rho = suppressWarnings(cor(
    beta[coordinate_estimable], continuum_beta[coordinate_estimable],
    method = "spearman"
  )),
  same_direction = sum(sign(beta[coordinate_estimable]) ==
                         sign(continuum_beta[coordinate_estimable])),
  same_direction_fraction = mean(sign(beta[coordinate_estimable]) ==
                                   sign(continuum_beta[coordinate_estimable]))
), by = collection]
family_summary[, family_label := family_labels[as.character(collection)]]
family_summary[, facet_label := sprintf(
  "%s | %d/%d; rho=%.2f",
  family_label, coordinate_estimable, complete_family, spearman_rho
)]

all_results <- merge(
  all_results,
  family_summary[, .(collection, family_label, facet_label)],
  by = "collection", all.x = TRUE, sort = FALSE
)
all_results[, facet_plot_label := sub(" | ", "\n", facet_label, fixed = TRUE)]
all_results[, facet_plot_label := factor(
  facet_plot_label,
  levels = sub(
    " | ", "\n",
    family_summary[match(family_order, collection), facet_label], fixed = TRUE
  )
)]

ml_write_tsv_once(all_results,
                  file.path(source_dir, "complete_pathway_stage_continuum_coordinates.tsv.gz"))
ml_write_tsv_once(family_summary[match(family_order, collection)],
                  file.path(source_dir, "pathway_family_summary.tsv"))

axis_x <- "Fibrosis-stage trend beta\n(pathway SD / stage increment)"
axis_y <- "Stage-adjusted fixed-projection beta\n(pathway SD / continuum SD)"
caption_text <- paste0(
  "Complete-family BH support; unestimable sets remain in the family and source table."
)

compact_breaks <- function(limits) {
  values <- pretty(limits, n = 3)
  values <- values[values >= limits[[1L]] - 1e-10 & values <= limits[[2L]] + 1e-10]
  if (length(values) > 3L) {
    values <- values[unique(round(seq(1, length(values), length.out = 3L)))]
  }
  values
}
compact_labels <- function(values) sprintf("%.1f", values)

base_scatter <- function(data, point_size = 0.7, point_alpha = 0.65) {
  ggplot(
    data[coordinate_estimable == TRUE],
    aes(beta, continuum_beta, color = support_class, shape = support_class)
  ) +
    geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
    geom_abline(slope = 1, intercept = 0, color = "#BDBDBD", linewidth = 0.25,
                linetype = "22") +
    geom_point(size = point_size, alpha = point_alpha, stroke = 0.25) +
    scale_color_manual(values = support_colors, drop = FALSE, name = NULL) +
    scale_shape_manual(values = support_shapes, drop = FALSE, name = NULL) +
    scale_x_continuous(breaks = compact_breaks, labels = compact_labels) +
    scale_y_continuous(breaks = compact_breaks, labels = compact_labels) +
    labs(x = axis_x, y = axis_y, caption = caption_text) +
    theme_effect() +
    guides(
      color = guide_legend(
        nrow = 2, byrow = TRUE,
        override.aes = list(
          size = 1.4, alpha = 1, shape = unname(support_shapes)
        )
      ),
      shape = "none"
    ) +
    theme(legend.position = "bottom")
}

overview <- base_scatter(all_results, point_size = 0.42, point_alpha = 0.42) +
  facet_wrap(~facet_plot_label, ncol = 3, scales = "free")
overview_path <- save_panel(
  overview, "s3_all_pathway_families_stage_vs_continuum_overview.pdf", 5.5, 5.4
)

fixed_hallmark_labels <- setNames(
  c(
    "Fatty acid metabolism", "Oxidative phosphorylation", "EMT",
    "TNF/NF-kB signaling", "Inflammatory response", "TGF-beta signaling"
  ),
  contract$display_hallmarks
)

render_family <- function(collection_id, filename, label_fixed_hallmarks = FALSE) {
  data <- all_results[collection == collection_id]
  summary <- family_summary[collection == collection_id]
  summary_line <- sprintf(
    "%d/%d estimable; rho=%.2f\nBH support: both %d | stage %d | continuum %d | neither %d",
    summary$coordinate_estimable, summary$complete_family, summary$spearman_rho,
    summary$both, summary$stage_only, summary$continuum_only, summary$neither
  )
  finite_range <- range(c(data$beta[data$coordinate_estimable],
                          data$continuum_beta[data$coordinate_estimable], 0), na.rm = TRUE)
  padding <- max(diff(finite_range) * 0.06, 0.02)
  limits <- finite_range + c(-padding, padding)
  point_size <- if (nrow(data) > 3000L) 0.34 else if (nrow(data) > 1000L) 0.42 else
    if (nrow(data) > 300L) 0.58 else 1.05
  point_alpha <- if (nrow(data) > 3000L) 0.34 else if (nrow(data) > 1000L) 0.42 else
    if (nrow(data) > 300L) 0.55 else 0.8
  plot <- base_scatter(data, point_size, point_alpha) +
    coord_equal(xlim = limits, ylim = limits) +
    labs(subtitle = summary_line)
  if (label_fixed_hallmarks) {
    labels <- data[set_id %in% names(fixed_hallmark_labels) & coordinate_estimable == TRUE]
    labels[, display_label := fixed_hallmark_labels[set_id]]
    plot <- plot +
      geom_point(data = labels, size = 1.7, shape = 21, fill = "white",
                 stroke = 0.45, show.legend = FALSE) +
      geom_text(
        data = labels, aes(label = display_label), color = "black",
        size = 6 / ggplot2::.pt, nudge_y = 0.035 * diff(limits),
        check_overlap = TRUE, show.legend = FALSE
      )
  }
  save_panel(plot, filename, 4.3, 4.1)
}

individual_paths <- c(
  hallmark = render_family(
    "hallmark", "s3_hallmark_stage_vs_continuum_effects.pdf", TRUE
  ),
  kegg = render_family("kegg", "s3_kegg_stage_vs_continuum_effects.pdf"),
  reactome = render_family("reactome", "s3_reactome_stage_vs_continuum_effects.pdf"),
  go_bp = render_family("go_bp", "s3_go_bp_stage_vs_continuum_effects.pdf"),
  go_mf = render_family("go_mf", "s3_go_mf_stage_vs_continuum_effects.pdf"),
  go_cc = render_family("go_cc", "s3_go_cc_stage_vs_continuum_effects.pdf")
)

figure_manifest <- data.table(
  figure_id = c("all_pathway_families", names(individual_paths)),
  path = c(overview_path, unname(individual_paths)),
  family = c("all", names(individual_paths)),
  role = c("complete_family_overview", rep("complete_single_family_overview", 6L)),
  inference_source = "existing_two_evaluation_cohort_stage_and_continuum_models",
  result_adaptive_feature_selection = FALSE
)
figure_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

input_manifest <- data.table(
  input = c("pathway_stage_continuum_source", "analysis_contract"),
  path = c(source_file, file.path(script_dir, "00_contract.json"))
)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))
writeLines(c(
  "# Complete pathway stage-versus-continuum overview candidate", "",
  "The overview retains the six complete pathway families: 50 Hallmarks, 658 KEGG sets, 1,787 Reactome sets, 7,583 GO biological-process sets, 1,855 GO molecular-function sets, and 1,042 GO cellular-component sets.", "",
  "Scatter coordinates are available only when both the standardized fibrosis-stage trend and the stage-adjusted fixed-projection coefficient are estimable. Every unestimable set remains explicit in the source table and family denominator.", "",
  "Colors reproduce the existing support classes under complete-family BH correction. The six labeled Hallmarks are the prespecified display roster, not result-selected labels. No model, threshold, membership, or pathway collection was changed."
), file.path(output_candidate, "README.md"))

message("PATHWAY_EFFECT_OVERVIEWS_COMPLETE: ", output_candidate)
