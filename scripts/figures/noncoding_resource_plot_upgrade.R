#!/usr/bin/env Rscript

# KEY MESSAGE: Noncoding biology strengthens the MASLD Resource through robust
# measurement and mechanism-aware interpretation, not through a separate lncRNA
# significance family or an untested mediation claim.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(digits = 17, scipen = 999)
set.seed(20260812)
grDevices::pdf.options(useDingbats = FALSE)

fail <- function(...) stop(..., call. = FALSE)
is_symlink <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}
sha256_file <- function(path) {
  value <- system2("sha256sum", normalizePath(path, mustWork = TRUE), stdout = TRUE)
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
write_tsv_once <- function(value, path) {
  connection <- file(path, open = "wx")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(
    as.data.frame(value), connection, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = "NA", eol = "\n"
  )
  close(connection)
}
parse_arguments <- function(values) {
  expected <- c(
    "bulk-counts", "biotype-association", "biotype-review", "stage-counts",
    "stage-review", "bulk", "cohort", "loco", "bulk-review",
    "program-content", "program-sensitivity", "program-review",
    "program-registry", "catalog-rules", "catalog-review", "output"
  )
  parsed <- list()
  for (value in values) {
    pieces <- strsplit(value, "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !startsWith(pieces[[1L]], "--")) {
      fail("Arguments must use --name=value syntax: ", value)
    }
    name <- sub("^--", "", pieces[[1L]])
    if (!name %in% expected || name %in% names(parsed)) fail("Unknown argument: ", name)
    parsed[[name]] <- pieces[[2L]]
  }
  missing <- setdiff(expected, names(parsed))
  if (length(missing)) fail("Missing arguments: ", paste(missing, collapse = ", "))
  parsed
}
review_passed <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  grepl('"status"[[:space:]]*:[[:space:]]*"pass"', text) ||
    grepl("(^|[[:space:]])PASS($|[[:space:]])", text)
}
plain_theme <- function() {
  theme_masld(base_size = 6) + theme_pub() +
    theme(
      text = element_text(size = 6, face = "plain", color = "black"),
      plot.title = element_text(size = 6, face = "plain", color = "black"),
      strip.text = element_text(size = 6, face = "plain", color = "black"),
      legend.position = "top",
      legend.direction = "horizontal"
    )
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(
  Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  mustWork = TRUE
)
theme_path <- normalizePath(file.path(project_root, "scripts/figures/publication_theme.R"), mustWork = TRUE)
source(theme_path, local = TRUE)

input_keys <- setdiff(names(arguments), "output")
inputs <- vapply(input_keys, function(key) normalizePath(arguments[[key]], mustWork = TRUE), character(1))
for (key in c("biotype-review", "stage-review", "bulk-review", "program-review", "catalog-review")) {
  if (!review_passed(inputs[[key]])) fail("Input review did not pass: ", key)
}
output_parent <- normalizePath(dirname(arguments[["output"]]), mustWork = TRUE)
output <- file.path(output_parent, basename(arguments[["output"]]))
if (file.exists(output) || is_symlink(output)) fail("Refusing existing output: ", output)

bulk_counts <- fread(inputs[["bulk-counts"]])
association <- fread(inputs[["biotype-association"]])
stage_counts <- fread(inputs[["stage-counts"]])
bulk <- fread(inputs[["bulk"]])
cohort <- fread(inputs[["cohort"]])
loco <- fread(inputs[["loco"]])
program_content <- fread(inputs[["program-content"]])
program_sensitivity <- fread(inputs[["program-sensitivity"]])
program_registry <- fread(inputs[["program-registry"]])
catalog_rules <- fread(inputs[["catalog-rules"]])

if (nrow(bulk_counts) != 3L || sum(bulk_counts$n_tested) != 23370L ||
    sum(bulk_counts$n_treat_positive) != 1616L) fail("Bulk biotype census drift")
if (nrow(association) != 1L || association$comparison[[1L]] != "lncRNA_vs_protein_coding") {
  fail("Biotype association contract drift")
}
if (nrow(stage_counts) != 12L || any(stage_counts$n_fdr_positive != stage_counts$n_up + stage_counts$n_down)) {
  fail("Stage biotype census drift")
}
high <- bulk[high_confidence == TRUE]
if (nrow(bulk) != 6249L || nrow(high) != 271L || nrow(cohort) != 6249L * 5L ||
    nrow(loco) != 6249L * 5L) fail("lncRNA robustness family drift")
if (nrow(program_content) != 117L || sum(program_content$requires_leave_all_lncrna_out_sensitivity) != 49L) {
  fail("Program-content family drift")
}
if (nrow(program_sensitivity) != 49L || sum(program_sensitivity$sensitivity_passed) != 48L ||
    sum(program_sensitivity$stage_direction_preserved_y) != 48L) fail("Program-sensitivity verdict drift")
if (nrow(catalog_rules) != 5L) fail("Catalog experiment-rule family drift")

biotype_levels <- c("protein_coding", "lncRNA", "other")
biotype_labels <- c(protein_coding = "Protein-coding", lncRNA = "lncRNA", other = "Other")
direction_colors <- c(Down = masld_colors$down, Up = masld_colors$up)

# Figure 3 pooled disease-state counts.
bulk_long <- rbindlist(list(
  bulk_counts[, .(display_biotype, n_tested, n_significant = n_treat_down,
                  signed_count = -n_treat_down, direction = "Down")],
  bulk_counts[, .(display_biotype, n_tested, n_significant = n_treat_up,
                  signed_count = n_treat_up, direction = "Up")]
))
bulk_long[, biotype_label := sprintf(
  "%s (%s tested)", biotype_labels[display_biotype],
  format(n_tested, big.mark = ",", scientific = FALSE, trim = TRUE)
)]
bulk_long[, biotype_label := factor(
  biotype_label,
  levels = rev(sprintf(
    "%s (%s tested)", biotype_labels[biotype_levels],
    format(
      bulk_counts[match(biotype_levels, display_biotype)]$n_tested,
      big.mark = ",", scientific = FALSE, trim = TRUE
    )
  ))
)]
bulk_long[, count_label := format(n_significant, big.mark = ",", scientific = FALSE)]
bulk_long[, label_hjust := fifelse(signed_count < 0, 1.12, -0.12)]
bulk_plot <- ggplot(bulk_long, aes(x = signed_count, y = biotype_label, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black") +
  geom_col(width = 0.56, color = NA) +
  geom_text(aes(label = count_label, hjust = label_hjust), size = GEOM_TEXT_6PT, show.legend = FALSE) +
  scale_fill_manual(values = direction_colors, breaks = c("Down", "Up")) +
  scale_x_continuous(
    labels = function(values) format(abs(values), big.mark = ",", scientific = FALSE),
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  labs(x = "TREAT-positive genes (lower in disease ←  → higher in disease)", y = NULL, fill = NULL) +
  plain_theme() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# Figure 3 cross-sectional stage counts.
transition_levels <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
transition_labels <- c(F0_to_F1 = "F0→F1", F1_to_F2 = "F1→F2", F2_to_F3 = "F2→F3", F3_to_F4 = "F3→F4")
stage_long <- rbindlist(list(
  stage_counts[, .(transition, display_biotype, n_tested, n_significant = n_down,
                   signed_count = -n_down, direction = "Down")],
  stage_counts[, .(transition, display_biotype, n_tested, n_significant = n_up,
                   signed_count = n_up, direction = "Up")]
))
stage_long[, transition_label := factor(transition_labels[transition], levels = rev(transition_labels))]
stage_denominators <- stage_counts[, .(n_tested = unique(n_tested)), by = display_biotype]
stage_denominators[, facet_label := sprintf(
  "%s (%s tested)", biotype_labels[display_biotype], format(n_tested, big.mark = ",", scientific = FALSE)
)]
stage_long <- merge(stage_long, stage_denominators, by = c("display_biotype", "n_tested"), all.x = TRUE)
stage_long[, facet_label := factor(
  facet_label,
  levels = stage_denominators[match(biotype_levels, display_biotype)]$facet_label
)]
stage_long[, count_label := format(n_significant, big.mark = ",", scientific = FALSE)]
stage_long[, label_hjust := fifelse(signed_count < 0, 1.10, -0.10)]
stage_plot <- ggplot(stage_long, aes(x = signed_count, y = transition_label, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black") +
  geom_col(width = 0.60, color = NA) +
  geom_text(aes(label = count_label, hjust = label_hjust), size = GEOM_TEXT_6PT, show.legend = FALSE) +
  facet_wrap(~facet_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = direction_colors, breaks = c("Down", "Up")) +
  scale_x_continuous(
    labels = function(values) format(abs(values), big.mark = ",", scientific = FALSE),
    expand = expansion(mult = c(0.24, 0.24))
  ) +
  labs(
    x = "Genes at FDR < 0.05 (lower at later stage ←  → higher at later stage; separate scales)",
    y = "Cross-sectional Kleiner contrast", fill = NULL
  ) +
  plain_theme() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing.x = unit(0.28, "cm"))

# Tested-versus-positive biotype composition.
composition <- rbindlist(list(
  bulk_counts[, .(universe = "All tested genes", display_biotype, n = n_tested)],
  bulk_counts[, .(universe = "TREAT-positive genes", display_biotype, n = n_treat_positive)]
))
composition[, total := sum(n), by = universe]
composition[, percent := 100 * n / total]
composition[, biotype_label := biotype_labels[display_biotype]]
composition[, biotype_y := match(display_biotype, biotype_levels)]
composition[, universe := factor(universe, levels = c("All tested genes", "TREAT-positive genes"))]
composition[, count_label := sprintf("%s / %s", format(n, big.mark = ","), format(total, big.mark = ","))]
composition[, label_y := biotype_y + fifelse(universe == "All tested genes", 0.14, -0.14)]
composition_plot <- ggplot(composition, aes(x = percent, y = biotype_y, group = biotype_label)) +
  geom_line(color = "#BDBDBD", linewidth = 0.35) +
  geom_point(aes(shape = universe, fill = universe), size = 2.2, color = "black", stroke = 0.35) +
  geom_text(
    aes(y = label_y, label = sprintf("%.1f%%", percent)),
    size = GEOM_TEXT_6PT, show.legend = FALSE
  ) +
  scale_shape_manual(values = c("All tested genes" = 21, "TREAT-positive genes" = 22)) +
  scale_fill_manual(values = c("All tested genes" = "white", "TREAT-positive genes" = masld_colors$up)) +
  scale_x_continuous(limits = c(0, 78), breaks = seq(0, 75, 25), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(breaks = seq_along(biotype_levels), labels = biotype_labels[biotype_levels]) +
  labs(x = "Share of each complete gene universe", y = NULL, shape = NULL, fill = NULL) +
  plain_theme() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# High-confidence lncRNA heatmap with analysis blocks.
setorder(high, logFC, gene_id_versioned)
high[, gene_order := seq_len(.N)]
high_ids <- high$gene_id_versioned
cohorts <- sort(unique(cohort$cohort))
heatmap <- rbindlist(list(
  high[, .(gene_id_versioned, analysis_group = "Pooled", analysis = "Pooled", effect = logFC)],
  cohort[gene_id_versioned %in% high_ids, .(
    gene_id_versioned, analysis_group = "Individual cohorts", analysis = sub("^GSE", "", cohort), effect = logFC
  )],
  loco[gene_id_versioned %in% high_ids, .(
    gene_id_versioned, analysis_group = "Leave-one-cohort-out",
    analysis = sub("^GSE", "", excluded_cohort), effect = logFC
  )]
))
heatmap <- merge(heatmap, high[, .(gene_id_versioned, gene_order)], by = "gene_id_versioned")
heatmap[, analysis_group := factor(
  analysis_group, levels = c("Pooled", "Individual cohorts", "Leave-one-cohort-out")
)]
heatmap[, analysis := factor(analysis, levels = c("Pooled", sub("^GSE", "", cohorts)))]
heatmap[, effect_display := pmax(-2, pmin(2, effect))]
if (nrow(heatmap) != 271L * 11L) fail("Heatmap source family drift")
heat_plot <- ggplot(heatmap, aes(x = analysis, y = gene_order, fill = effect_display)) +
  geom_raster() +
  facet_grid(. ~ analysis_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = masld_colors$down, mid = "white", high = masld_colors$up,
    midpoint = 0, limits = c(-2, 2), name = "log2 FC\n(clipped)"
  ) +
  labs(x = "GEO accession", y = "271 high-confidence lncRNAs\n(ordered by pooled effect)") +
  plain_theme() +
  theme(
    legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.spacing.x = unit(0.12, "cm"), panel.border = element_rect(color = "#D0D0D0", fill = NA, linewidth = 0.25)
  )

# Genomic classes for the high-confidence display set.
class_counts <- high[, .N, by = lncrna_genomic_class][order(N)]
class_counts[, percent := 100 * N / sum(N)]
class_counts[, class_label := gsub("_", " ", lncrna_genomic_class)]
class_counts[, class_label := factor(class_label, levels = class_label)]
class_counts[, value_label := sprintf("%d (%.1f%%)", N, percent)]
class_plot <- ggplot(class_counts, aes(x = N, y = class_label)) +
  geom_segment(aes(x = 0, xend = N, yend = class_label), color = "#BDBDBD", linewidth = 0.45) +
  geom_point(color = masld_colors$up, size = 2.0) +
  geom_text(aes(label = value_label), hjust = -0.12, size = GEOM_TEXT_6PT) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = "High-confidence lncRNAs, n (% of 271)", y = NULL) +
  plain_theme() +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank())

# Deterministic example forest.
example_id <- "ENSG00000213062.6"
example_class <- high[gene_id_versioned == example_id, lncrna_genomic_class]
if (length(example_class) != 1L) fail("Deterministic example identity drift")
example <- rbindlist(list(
  bulk[gene_id_versioned == example_id, .(
    source_group = "Pooled", analysis = "Pooled", estimate = logFC,
    lower = logFC - 1.96 * SE, upper = logFC + 1.96 * SE
  )],
  cohort[gene_id_versioned == example_id, .(
    source_group = "Individual cohorts", analysis = cohort,
    estimate = logFC, lower = CI_low, upper = CI_high
  )],
  loco[gene_id_versioned == example_id, .(
    source_group = "Leave-one-cohort-out", analysis = paste0("LOO ", excluded_cohort),
    estimate = logFC, lower = CI_low, upper = CI_high
  )]
))
if (nrow(example) != 11L) fail("Deterministic example family drift")
example[, source_group := factor(
  source_group, levels = c("Pooled", "Individual cohorts", "Leave-one-cohort-out")
)]
example[, analysis := factor(analysis, levels = rev(c(
  "Pooled", cohorts, paste0("LOO ", cohorts)
)))]
forest_plot <- ggplot(example, aes(x = estimate, y = analysis)) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper, color = source_group),
    orientation = "y", width = 0.16, linewidth = 0.35
  ) +
  geom_point(aes(color = source_group, shape = source_group), size = 1.35) +
  scale_color_manual(values = c(
    "Pooled" = "black", "Individual cohorts" = masld_colors$up,
    "Leave-one-cohort-out" = masld_colors$down
  )) +
  scale_shape_manual(values = c("Pooled" = 18, "Individual cohorts" = 16, "Leave-one-cohort-out" = 17)) +
  labs(
    x = "Disease versus control log2 fold change", y = NULL, color = NULL, shape = NULL,
    title = sprintf("%s; %s", example_id, gsub("_", " ", example_class))
  ) +
  plain_theme() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# Complete program lncRNA content.
program_content[, cell_label := factor(
  cell_type,
  levels = c("tcells", "macrophages", "fibroblasts", "cholangiocytes", "hepatocytes"),
  labels = c("T cells", "Macrophages", "Fibroblasts", "Cholangiocytes", "Hepatocytes")
)]
program_content[, l1_percent := 100 * lncrna_l1_fraction]
cell_colors <- c(
  "T cells" = "#7B1FA2", "Macrophages" = "#C9265E", "Fibroblasts" = "#E69F00",
  "Cholangiocytes" = "#0072B2", "Hepatocytes" = "#009E73"
)
program_content_plot <- ggplot(program_content, aes(x = l1_percent, y = cell_label, color = cell_label)) +
  geom_vline(xintercept = 5, color = "#9E9E9E", linewidth = 0.35, linetype = "dashed") +
  geom_jitter(height = 0.15, width = 0, size = 1.25, alpha = 0.85) +
  scale_color_manual(values = cell_colors) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.01, 0.05))) +
  labs(x = "Original program weight carried by uniquely mapped lncRNAs", y = NULL, color = NULL) +
  plain_theme() +
  theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  annotate("text", x = 5, y = 5.45, label = "5% sensitivity trigger", hjust = -0.05, size = GEOM_TEXT_6PT)

# Native leave-all-lncRNA-out sensitivity.
registry_small <- program_registry[, .(program_uid, module_name, robust_display)]
program_sensitivity <- merge(program_sensitivity, registry_small, by = "program_uid", all.x = TRUE)
if (anyNA(program_sensitivity$module_name)) fail("Program registry join failed")
program_sensitivity[, display_state := fifelse(
  !sensitivity_passed, "Direction changed",
  fifelse(robust_display, "Robust display program", "Sensitivity passed")
)]
program_sensitivity[, display_state := factor(
  display_state, levels = c("Sensitivity passed", "Robust display program", "Direction changed")
)]
highlight_uids <- c(
  "hotspot_hepatocytes_f05c535ae5bbc0b9",
  "hotspot_hepatocytes_48f39dd4d817a10e",
  "hotspot_cholangiocytes_576c469934529fdd"
)
program_sensitivity[, point_label := fifelse(program_uid %in% highlight_uids, module_name, "")]
program_sensitivity[, label_hjust := fifelse(
  program_uid %in% c(
    "hotspot_hepatocytes_f05c535ae5bbc0b9",
    "hotspot_hepatocytes_48f39dd4d817a10e"
  ), 1.04, -0.04
)]
program_sensitivity[, label_vjust := fifelse(
  program_uid == "hotspot_hepatocytes_f05c535ae5bbc0b9", -0.65,
  fifelse(program_uid == "hotspot_hepatocytes_48f39dd4d817a10e", 1.45, -0.65)
)]
program_sensitivity_plot <- ggplot(
  program_sensitivity,
  aes(x = original_stage_beta, y = without_lncrna_stage_beta_y)
) +
  geom_abline(intercept = 0, slope = 1, color = "#9E9E9E", linewidth = 0.35, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(aes(color = display_state, shape = display_state), size = 1.75, stroke = 0.4) +
  geom_text(
    data = program_sensitivity[point_label != ""],
    aes(label = point_label, hjust = label_hjust, vjust = label_vjust),
    size = GEOM_TEXT_6PT, show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Sensitivity passed" = "#9E9E9E", "Robust display program" = masld_colors$up,
    "Direction changed" = "#D55E00"
  )) +
  scale_shape_manual(values = c(
    "Sensitivity passed" = 21, "Robust display program" = 23, "Direction changed" = 24
  )) +
  coord_equal() +
  labs(
    x = "Original donor-level stage coefficient",
    y = "Stage coefficient after removing lncRNA members",
    color = NULL, shape = NULL
  ) +
  plain_theme() +
  annotate("text", x = -Inf, y = Inf, label = "48/49 preserve direction", hjust = -0.05, vjust = 1.25, size = GEOM_TEXT_6PT)

# Figure 5 molecular-object routing. This is a schema panel, not a biological result.
route_source <- data.table(
  row = 5:1,
  candidate_object = c(
    "Noncoding regulatory DNA", "Disease-associated lncRNA",
    "Colocalized or ambiguous lncRNA locus", "Protein-coding gene", "Untestable nomination"
  ),
  open_question = c(
    "Which target changes?", "Does the RNA product matter?",
    "DNA element, transcription, or RNA?", "Which context reveals function?", "Can the assay observe it?"
  ),
  next_test = c(
    "Allele-aware reporter or base editing", "RNA depletion with neighboring-gene readout",
    "Paired RNA depletion and locus perturbation", "Context-specific perturbation and rescue",
    "Acquire the missing direct measurement"
  ),
  rule_id = c(
    "EXP_REGULATORY_DNA_V2", "EXP_LNCRNA_RNA_PRODUCT_V2",
    "EXP_LNCRNA_LOCUS_DISAMBIGUATION_V2", "EXP_PROTEIN_STATE_CONTEXT_V2",
    "EXP_MEASURE_MISSING_ASSAY_V2"
  )
)
if (!setequal(route_source$rule_id, catalog_rules$experiment_rule_id)) fail("Catalog route mapping drift")
wrap_text <- function(values, width) {
  vapply(values, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}
route_plot_data <- copy(route_source)
route_plot_data[, candidate_object_plot := wrap_text(candidate_object, 25)]
route_plot_data[, open_question_plot := wrap_text(open_question, 22)]
route_plot_data[, next_test_plot := wrap_text(next_test, 37)]
route_plot <- ggplot(route_plot_data) +
  geom_rect(aes(xmin = 0.0, xmax = 2.0, ymin = row - 0.34, ymax = row + 0.34), fill = "#F2F2F2", color = "#BDBDBD", linewidth = 0.25) +
  geom_rect(aes(xmin = 2.65, xmax = 4.55, ymin = row - 0.34, ymax = row + 0.34), fill = "#F7E4EC", color = "#C9265E", linewidth = 0.25) +
  geom_rect(aes(xmin = 5.2, xmax = 8.0, ymin = row - 0.34, ymax = row + 0.34), fill = "#E8F0F8", color = "#1565C0", linewidth = 0.25) +
  geom_segment(aes(x = 2.05, xend = 2.55, y = row, yend = row), arrow = arrow(length = unit(0.08, "in")), linewidth = 0.3, color = "#777777") +
  geom_segment(aes(x = 4.60, xend = 5.10, y = row, yend = row), arrow = arrow(length = unit(0.08, "in")), linewidth = 0.3, color = "#777777") +
  geom_text(aes(x = 1.0, y = row, label = candidate_object_plot), size = GEOM_TEXT_6PT, lineheight = 0.92) +
  geom_text(aes(x = 3.6, y = row, label = open_question_plot), size = GEOM_TEXT_6PT, lineheight = 0.92) +
  geom_text(aes(x = 6.6, y = row, label = next_test_plot), size = GEOM_TEXT_6PT, lineheight = 0.92) +
  annotate("text", x = c(1.0, 3.6, 6.6), y = 5.72, label = c("Candidate object", "Open question", "Smallest discriminating test"), size = GEOM_TEXT_6PT) +
  coord_cartesian(xlim = c(-0.05, 8.05), ylim = c(0.55, 5.85), clip = "off") +
  theme_void(base_size = 6) +
  theme(
    text = element_text(size = 6, family = "Helvetica", face = "plain", color = "black"),
    plot.margin = margin(3, 3, 3, 3)
  )

tmp <- file.path(output_parent, paste0(".", basename(output), ".tmp.", Sys.getenv("SLURM_JOB_ID", Sys.getpid())))
if (file.exists(tmp)) fail("Temporary output exists: ", tmp)
dir.create(tmp, recursive = FALSE)
if (!dir.exists(tmp)) fail("Could not create candidate output")

plots <- list(
  fig3_pooled_deg_biotype_v2 = list(plot = bulk_plot, width = 3.55, height = 1.95, source = bulk_long),
  fig3_stage_deg_biotype_v2 = list(plot = stage_plot, width = 5.45, height = 2.35, source = stage_long),
  figS_lncrna_universe_composition = list(plot = composition_plot, width = 3.65, height = 2.15, source = composition),
  figS_lncrna_robustness_heatmap_v2 = list(plot = heat_plot, width = 4.75, height = 3.25, source = heatmap),
  figS_lncrna_genomic_classes_v2 = list(plot = class_plot, width = 3.75, height = 2.25, source = class_counts),
  figS_lncrna_example_forest_v2 = list(plot = forest_plot, width = 4.25, height = 3.05, source = example),
  figS_program_lncrna_content = list(plot = program_content_plot, width = 4.15, height = 2.45, source = program_content),
  figS_program_without_lncrna_sensitivity = list(plot = program_sensitivity_plot, width = 4.1, height = 3.35, source = program_sensitivity),
  fig5_noncoding_object_routes = list(plot = route_plot, width = 5.35, height = 3.0, source = route_source)
)
for (name in names(plots)) {
  item <- plots[[name]]
  ggsave(
    file.path(tmp, paste0(name, ".pdf")), item$plot,
    width = item$width, height = item$height, device = cairo_pdf, units = "in"
  )
  write_tsv_once(item$source, file.path(tmp, paste0(name, "_source.tsv")))
}

write_tsv_once(
  data.table(
    input_role = input_keys,
    path = inputs,
    size_bytes = file.info(inputs)$size,
    sha256 = vapply(inputs, sha256_file, character(1))
  ),
  file.path(tmp, "input_manifest.tsv")
)
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"), useBytes = TRUE)
artifacts <- setdiff(list.files(tmp), "output_manifest.tsv")
write_tsv_once(
  data.table(
    relative_path = artifacts,
    size_bytes = file.info(file.path(tmp, artifacts))$size,
    sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1))
  ),
  file.path(tmp, "output_manifest.tsv")
)
if (!file.rename(tmp, output)) fail("Atomic candidate publication failed")
cat("NONCODING_PLOT_UPGRADE_COMPLETE\t", output, "\n", sep = "")
