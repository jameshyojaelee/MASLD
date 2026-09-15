#!/usr/bin/env Rscript
# KEY MESSAGE: Only one molecular system exceeds size-matched regions of the same membership geometry, and it does so against the permissive reference only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    paste(
      "Usage: 74_render_system_competitive_null.R",
      "SYSTEM_STRESS_SOURCE OUTPUT_CANDIDATE"
    ),
    call. = FALSE
  )
}

stress_source <- normalizePath(args[[1L]], mustWork = TRUE)
output_candidate <- args[[2L]]
if (file.exists(output_candidate)) {
  stop("Refusing to overwrite output candidate: ", output_candidate, call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
contract <- ml_read_contract()

candidate_root <- normalizePath(
  file.path(ml_project_root(), "figures", "candidates"), mustWork = TRUE
)
candidate_parent <- normalizePath(dirname(output_candidate), mustWork = TRUE)
ml_assert(
  startsWith(paste0(candidate_parent, "/"), paste0(candidate_root, "/")),
  "Competitive-null rendering is restricted to figures/candidates"
)

panel_dir <- file.path(output_candidate, "panels")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(panel_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

system_family_size <- 43L
focal_system <- "system_01"
axes_primary <- as.character(contract$co_primary_axes)
null_levels <- c("connected_subgraph", "size_matched_uniform")
null_titles <- c(
  connected_subgraph = "Connected subgraph (conservative)",
  size_matched_uniform = "Size-matched nodes (permissive)"
)
axis_titles <- c(signature_pc1 = "Cohort PC1", fixed_projection = "Fixed projection")

draws_path <- file.path(
  stress_source, "tables", "system_competitive_null_draws.tsv.gz"
)
summary_path <- file.path(stress_source, "tables", "system_competitive_null.tsv")
checksums_path <- file.path(stress_source, "provenance", "output_checksums.tsv")
required_inputs <- c(draws_path, summary_path, checksums_path)
ml_assert(all(file.exists(required_inputs)), "A competitive-null input is missing")

checksums <- fread(checksums_path)
for (path in c(draws_path, summary_path)) {
  target_basename <- basename(path)
  expected <- checksums[basename(get("path")) == target_basename, unique(sha256)]
  ml_assert(length(expected) == 1L,
            paste0("Stress-test checksum is absent or ambiguous: ", path))
  ml_assert(identical(ml_sha256(path), expected),
            paste0("Stress-test source hash drift: ", path))
}

summaries <- fread(summary_path, na.strings = c("", "NA"))
draws <- fread(draws_path, na.strings = c("", "NA"))
ml_assert(nrow(summaries) == system_family_size * 2L * length(axes_primary),
          "The complete competitive-null family drifted")
ml_assert(setequal(unique(summaries$null_model), null_levels),
          "Competitive null models drifted")
ml_assert(sum(summaries$competitive_q_value < 0.05 &
                summaries$null_model == "connected_subgraph") == 0L,
          "A system now clears the conservative null; revisit the panel message")

summaries[, null_model := factor(null_model, levels = null_levels)]
draws[, null_model := factor(null_model, levels = null_levels)]

theme_supplement <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, color = "black", face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.title = element_text(size = 6, face = "plain"),
      strip.text = element_text(size = 6, color = "black", face = "plain"),
      strip.background = element_blank(),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(2.6, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

null_colors <- c(
  connected_subgraph = "#1565C0",
  size_matched_uniform = "#9E9E9E"
)

# ---------------------------------------------------------------------------
# Panel one: the focal system against its own null distributions
# ---------------------------------------------------------------------------

focal_draws <- draws[community_id == focal_system]
focal_summary <- summaries[community_id == focal_system]
ml_assert(nrow(focal_summary) == 2L * length(axes_primary),
          "The focal system is missing a competitive-null row")
focal_draws[, panel_null := null_titles[as.character(null_model)]]
focal_draws[, panel_axis := axis_titles[axis_id]]
focal_summary[, `:=`(
  panel_null = null_titles[as.character(null_model)],
  panel_axis = axis_titles[axis_id]
)]
focal_draws[, panel_null := factor(panel_null, levels = null_titles[null_levels])]
focal_summary[, panel_null := factor(panel_null, levels = null_titles[null_levels])]

focal_plot <- ggplot(focal_draws, aes(null_absolute_meta_z)) +
  geom_histogram(
    aes(fill = null_model), bins = 40, color = NA, show.legend = FALSE
  ) +
  geom_vline(
    data = focal_summary, aes(xintercept = observed_absolute_meta_z),
    color = "black", linewidth = 0.4
  ) +
  geom_text(
    data = focal_summary,
    aes(x = observed_absolute_meta_z, y = Inf,
        label = sprintf("P = %.4f", competitive_p_value)),
    hjust = 1.08, vjust = 1.6, size = 6 / ggplot2::.pt, color = "black"
  ) +
  facet_grid(panel_null ~ panel_axis, scales = "free_y") +
  scale_fill_manual(values = null_colors) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    expand = expansion(mult = c(0.03, 0.12))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = "|Meta z| of size-matched draws (black line: S01 observed)",
    y = "Draws"
  ) +
  theme_supplement()

focal_path <- file.path(panel_dir, "s4v_system01_competitive_null.pdf")
ggsave(focal_path, focal_plot, width = 4.3, height = 2.9, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
ml_assert(file.info(focal_path)$size > 1000L,
          "Competitive-null focal panel PDF is unexpectedly small")

# ---------------------------------------------------------------------------
# Panel two: the family-wide gate quantity for all 43 systems
# ---------------------------------------------------------------------------

# Support requires both co-primary axes, so the gate quantity is the worse of
# the two BH q values within each system and null model.
gate <- summaries[, .(
  gate_q_value = max(competitive_q_value),
  observed_absolute_meta_z = min(observed_absolute_meta_z),
  n_features = unique(n_features)
), by = .(community_id, null_model)]
ml_assert(nrow(gate) == system_family_size * 2L, "Competitive gate family drifted")

order_frame <- gate[null_model == "connected_subgraph"]
setorder(order_frame, -gate_q_value)
gate[, community_id := factor(community_id, levels = order_frame$community_id)]
gate[, short_label := sprintf(
  "S%02d", as.integer(sub("system_", "", as.character(community_id)))
)]
gate[, short_label := factor(short_label, levels = sprintf(
  "S%02d", as.integer(sub("system_", "", order_frame$community_id))
))]

gate_plot <- ggplot(gate, aes(gate_q_value, short_label)) +
  geom_vline(xintercept = 0.05, color = "#D0D0D0", linewidth = 0.3,
             linetype = "22") +
  geom_point(aes(color = null_model), shape = 16, size = 1.1, stroke = 0) +
  scale_color_manual(
    values = null_colors, labels = null_titles, name = NULL,
    breaks = null_levels
  ) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0.05, 0.25, 0.5, 0.75, 1),
    labels = function(x) formatC(x, format = "f", digits = 2),
    expand = expansion(mult = 0.02)
  ) +
  labs(
    x = "Competitive BH q, worse of the two co-primary axes",
    y = NULL
  ) +
  theme_supplement() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 5, color = "black"),
    panel.grid.major.y = element_line(linewidth = 0.15, color = "#F0F0F0")
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 1.4)))

gate_path <- file.path(panel_dir, "s4w_system_competitive_gate.pdf")
ggsave(gate_path, gate_plot, width = 3.2, height = 3.9, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
ml_assert(file.info(gate_path)$size > 1000L,
          "Competitive gate panel PDF is unexpectedly small")

# ---------------------------------------------------------------------------
# Source tables, manifest, provenance
# ---------------------------------------------------------------------------

focal_source <- focal_summary[, .(
  community_id, null_model, axis_id, n_features,
  observed_absolute_meta_z, null_median_absolute_meta_z,
  null_p95_absolute_meta_z, competitive_p_value, competitive_q_value,
  finite_draws, requested_draws
)]
setorder(focal_source, null_model, axis_id)
ml_write_tsv_once(
  focal_source, file.path(source_dir, "s4v_system01_competitive_null.tsv")
)

gate_source <- gate[, .(
  community_id = as.character(community_id), short_label = as.character(short_label),
  null_model = as.character(null_model), n_features, gate_q_value
)]
setorder(gate_source, null_model, community_id)
ml_write_tsv_once(
  gate_source, file.path(source_dir, "s4w_system_competitive_gate.tsv")
)

figure_manifest <- data.table(
  callout = c("S4V", "S4W"),
  path = c(focal_path, gate_path),
  sha256 = c(ml_sha256(focal_path), ml_sha256(gate_path)),
  complete_system_family = c("1/43 focal", "43/43"),
  connected_null_supported = sprintf(
    "%d/43", sum(gate[null_model == "connected_subgraph", gate_q_value < 0.05])
  ),
  uniform_null_supported = sprintf(
    "%d/43", sum(gate[null_model == "size_matched_uniform", gate_q_value < 0.05])
  ),
  draws_per_system_per_null = unique(summaries$requested_draws),
  inferential_unit = "molecular system",
  embedded_explanatory_footnote = FALSE, promoted = FALSE
)
ml_write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

input_manifest <- data.table(
  path = normalizePath(required_inputs, mustWork = TRUE),
  sha256 = vapply(required_inputs, ml_sha256, character(1))
)
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))

audit <- data.table(
  check = c(
    "stress_input_hashes", "complete_competitive_family",
    "conservative_null_supports_nothing", "focal_system_is_marked_system",
    "both_axes_gate", "panels_nonempty"
  ),
  passed = c(
    TRUE,
    nrow(gate) == system_family_size * 2L,
    sum(gate[null_model == "connected_subgraph", gate_q_value < 0.05]) == 0L,
    sum(gate[null_model == "size_matched_uniform", gate_q_value < 0.05]) == 1L &&
      gate[null_model == "size_matched_uniform" & gate_q_value < 0.05,
           as.character(community_id)] == focal_system,
    TRUE,
    file.info(focal_path)$size > 1000L && file.info(gate_path)$size > 1000L
  ),
  detail = c(
    paste(basename(required_inputs), collapse = ";"),
    "All 43 systems appear under both null models",
    "No system exceeds the connected-subgraph reference at BH q<0.05",
    paste0("Only ", focal_system, " clears the permissive reference"),
    "The gate quantity is the worse BH q across the two co-primary axes",
    paste(basename(c(focal_path, gate_path)), collapse = ";")
  )
)
ml_assert(all(audit$passed), "Competitive-null panel audit failed")
ml_write_tsv_once(audit, file.path(provenance_dir, "audit.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

output_paths <- c(
  focal_path, gate_path,
  file.path(source_dir, "s4v_system01_competitive_null.tsv"),
  file.path(source_dir, "s4w_system_competitive_gate.tsv"),
  file.path(output_candidate, "figure_manifest.tsv"),
  file.path(provenance_dir, "audit.tsv")
)
ml_write_tsv_once(
  data.table(
    path = output_paths,
    sha256 = vapply(output_paths, ml_sha256, character(1))
  ),
  file.path(provenance_dir, "output_checksums.tsv")
)

cat("Competitive-null panels complete:", normalizePath(output_candidate), "\n")
