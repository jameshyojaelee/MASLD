#!/usr/bin/env Rscript
# KEY MESSAGE: The two frozen, disease-associated programs are represented by
# paired biological donors in distinct non-parenchymal compartments, while their
# fixed member genes show cross-sectional bulk tissue-state transport.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE)
grDevices::pdf.options(useDingbats = FALSE)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT <- Sys.getenv(
  "FIG4DE_GEOMETRY_OUT",
  unset = file.path(BASE, "figures/candidates/fig4de-geometry-review-2026-08-17-v5")
)
if (dir.exists(OUT)) {
  stop("Candidate output already exists: ", OUT, call. = FALSE)
}

source(file.path(BASE, "scripts/figures/publication_theme.R"))

dir.create(file.path(OUT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "source_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "captions"), recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(value)) stop("sha256sum failed: ", path, call. = FALSE)
  strsplit(value[[1]], "[[:space:]]+")[[1]][[1]]
}

require_that <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

check_rows <- list()
check <- function(name, condition, observed, expected) {
  check_rows[[length(check_rows) + 1L]] <<- data.table(
    check = name,
    status = if (isTRUE(condition)) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected)
  )
  require_that(condition, paste0(name, ": observed ", observed, "; expected ", expected))
}

score_path <- file.path(
  BASE,
  "Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2/results/hero_lineage_scores.tsv.gz"
)
contrast_path <- file.path(
  BASE,
  "Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2/results/hero_lineage_contrasts.tsv"
)
bulk_path <- file.path(
  BASE,
  "figures/main/fig4_singlecell_programs/source_tables/current_candidate/fig4e_bulk_tissue_state_transport.tsv"
)
theme_path <- file.path(BASE, "scripts/figures/publication_theme.R")

inputs <- data.table(
  role = c("donor_score_input", "lineage_contrast_input", "bulk_transport_input", "theme"),
  path = c(score_path, contrast_path, bulk_path, theme_path)
)
require_that(all(file.exists(inputs$path)), "one or more input files are missing")
inputs[, `:=`(sha256 = vapply(path, sha256, character(1)), bytes = file.info(path)$size)]
fwrite(inputs, file.path(OUT, "input_manifest.tsv"), sep = "\t")

hero <- data.table(
  program_uid = c(
    "hotspot_hepatocytes_f05c535ae5bbc0b9",
    "hotspot_hepatocytes_48f39dd4d817a10e"
  ),
  program_name = c("ECM/IGFBP7", "Ductular-injury/BICC1"),
  target_lineage = c("Fibroblasts", "Cholangiocytes"),
  target_label = c("Fibroblasts", "Cholangiocytes"),
  lineage_plot_label = c(
    "atop('ECM'~'/'~italic(IGFBP7), 'Fibroblasts'~'-'~'hepatocytes')",
    "atop('Ductular injury'~'/'~italic(BICC1), 'Cholangiocytes'~'-'~'hepatocytes')"
  ),
  bulk_plot_label = c("'ECM'~'/'~italic(IGFBP7)", "'Ductular injury'~'/'~italic(BICC1)"),
  program_color = c("#FF6F00", "#00BCD4"),
  program_shape = c(21L, 22L)
)

scores <- fread(score_path, na.strings = c("", "NA"))
program_columns <- hero$program_uid
check("score_columns_present", all(program_columns %in% names(scores)),
      paste(program_columns[program_columns %in% names(scores)], collapse = ";"),
      paste(program_columns, collapse = ";"))
scores <- scores[
  annotation_filter == "all_annotated_cells" & minimum_cells == 50 &
    lineage %chin% c("Hepatocytes", hero$target_lineage),
  c("dataset", "donor", "lineage", program_columns), with = FALSE
]

contrasts <- fread(contrast_path, na.strings = c("", "NA"))
contrasts <- contrasts[
  annotation_filter == "all_annotated_cells" & minimum_cells == 50 &
    program_uid %chin% hero$program_uid
]
check("complete_lineage_family", nrow(contrasts) == 10L, nrow(contrasts), 10L)

donor_rows <- list()
dataset_rows <- list()
pooled_rows <- list()
for (i in seq_len(nrow(hero))) {
  h <- hero[i]
  long <- scores[lineage %chin% c("Hepatocytes", h$target_lineage),
                 .(dataset, donor, lineage, score = get(h$program_uid))]
  wide <- dcast(long, dataset + donor ~ lineage, value.var = "score")
  require_that(
    all(c("Hepatocytes", h$target_lineage) %in% names(wide)),
    paste0("missing paired lineages for ", h$program_name)
  )
  wide <- wide[is.finite(Hepatocytes) & is.finite(get(h$target_lineage))]
  wide[, difference := get(h$target_lineage) - Hepatocytes]
  counts <- wide[, .N, by = dataset]
  eligible <- counts[N >= 3, dataset]
  wide <- wide[dataset %chin% eligible]
  focal <- contrasts[
    program_uid == h$program_uid & comparison_lineage == h$target_lineage
  ]
  check(paste0("focal_contrast_", i), nrow(focal) == 1L, nrow(focal), 1L)
  check(paste0("focal_donor_count_", i), uniqueN(wide$donor) == focal$n_paired_donors,
        uniqueN(wide$donor), focal$n_paired_donors)
  check(paste0("focal_dataset_count_", i), uniqueN(wide$dataset) == focal$n_datasets,
        uniqueN(wide$dataset), focal$n_datasets)
  equal_dataset_beta <- mean(wide[, mean(difference), by = dataset]$V1)
  check(paste0("equal_dataset_beta_", i),
        isTRUE(all.equal(equal_dataset_beta, focal$beta, tolerance = 1e-10)),
        format(equal_dataset_beta, digits = 15), format(focal$beta, digits = 15))
  check(paste0("focal_hc3_q_", i), focal$hc3_qvalue < 0.05,
        format(focal$hc3_qvalue, scientific = TRUE), "<0.05")

  shared <- list(
    program_uid = h$program_uid,
    program_name = h$program_name,
    target_lineage = h$target_lineage,
    reference_lineage = "Hepatocytes",
    biological_unit = "paired biological donor",
    multiple_testing_family = "10 prespecified program-by-lineage contrasts, BH across the complete family",
    unresolved_alternative = "within-lineage state composition, ambient carryover, and shared multicellular response remain possible",
    claim_boundary = "same-atlas donor-paired enrichment relative to hepatocytes; not transcript origin, secretion, cell autonomy, or exclusive lineage specificity"
  )
  donor_rows[[i]] <- cbind(
    wide[, .(dataset, donor, hepatocyte_score = Hepatocytes,
             target_score = get(h$target_lineage), difference)],
    as.data.table(shared),
    display_role = "individual_donor_difference",
    evidence_state = "underlying_donor_observation",
    state_reason = "within-dataset standardized target-lineage minus hepatocyte score"
  )
  dataset_rows[[i]] <- cbind(
    wide[, .(difference = mean(difference), n_paired_donors = .N), by = dataset],
    as.data.table(shared),
    display_role = "dataset_mean",
    evidence_state = "descriptive_dataset_mean",
    state_reason = "unweighted mean of paired biological-donor differences within one eligible dataset"
  )
  pooled_rows[[i]] <- cbind(
    data.table(
      dataset = "Pooled HC3",
      difference = focal$beta,
      hc3_ci_low = focal$hc3_ci_low,
      hc3_ci_high = focal$hc3_ci_high,
      n_paired_donors = focal$n_paired_donors,
      n_datasets = focal$n_datasets,
      hc3_qvalue = focal$hc3_qvalue
    ),
    as.data.table(shared),
    display_role = "equal_dataset_weighted_hc3_estimate",
    evidence_state = "hc3_family_supported",
    state_reason = "equal average of dataset-specific donor-paired differences; robust HC3 uncertainty and BH across ten prespecified contrasts"
  )
}

donor_plot <- rbindlist(donor_rows, use.names = TRUE, fill = TRUE)
dataset_plot <- rbindlist(dataset_rows, use.names = TRUE, fill = TRUE)
pooled_plot <- rbindlist(pooled_rows, use.names = TRUE, fill = TRUE)

fwrite(donor_plot, file.path(OUT, "source_tables/fig4d_donor_paired_differences.tsv"), sep = "\t")
fwrite(dataset_plot, file.path(OUT, "source_tables/fig4d_dataset_means.tsv"), sep = "\t")
fwrite(pooled_plot, file.path(OUT, "source_tables/fig4d_equal_dataset_hc3_estimates.tsv"), sep = "\t")

dataset_levels <- c(
  "GSE174748", "GSE185477", "GSE202379", "GSE244832", "Liver_Atlas", "Pooled HC3"
)
for (frame in list(donor_plot, dataset_plot, pooled_plot)) {
  frame[, dataset_plot_label := factor(dataset, levels = rev(dataset_levels))]
}
donor_plot[, dataset_plot_label := factor(dataset, levels = rev(dataset_levels))]
dataset_plot[, dataset_plot_label := factor(dataset, levels = rev(dataset_levels))]
pooled_plot[, dataset_plot_label := factor(dataset, levels = rev(dataset_levels))]

theme_panel <- function() {
  theme_masld(base_size = 6) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 6, face = "plain", margin = margin(b = 2)),
      legend.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

set.seed(20260817L)
donor_plot <- merge(donor_plot, hero[, .(program_uid, lineage_plot_label, program_color, program_shape)],
                    by = "program_uid", all.x = TRUE, sort = FALSE)
dataset_plot <- merge(dataset_plot, hero[, .(program_uid, lineage_plot_label, program_color, program_shape)],
                      by = "program_uid", all.x = TRUE, sort = FALSE)
pooled_plot <- merge(pooled_plot, hero[, .(program_uid, lineage_plot_label, program_color, program_shape)],
                     by = "program_uid", all.x = TRUE, sort = FALSE)

lineage_order <- c(
  "Fibroblasts", "Cholangiocytes", "Endothelial cells", "Macrophages", "T cells"
)
forest <- merge(
  contrasts,
  hero[, .(
    program_uid, target_lineage, bulk_plot_label, program_color, program_shape,
    program_key = program_name
  )],
  by = "program_uid", all.x = TRUE, sort = FALSE
)
forest[, lineage_label := factor(comparison_lineage, levels = rev(lineage_order))]
forest[, target := comparison_lineage == target_lineage]
forest[, `:=`(
  biological_unit = "paired biological donor",
  multiple_testing_family = "10 prespecified program-by-lineage contrasts, BH across the complete family",
  evidence_state = fifelse(hc3_qvalue < 0.05, "hc3_family_supported", "tested_nonsignificant"),
  state_reason = fifelse(
    hc3_qvalue < 0.05,
    "lineage-minus-hepatocyte contrast passes HC3 BH q<0.05 in the complete ten-comparison family",
    "lineage-minus-hepatocyte contrast does not pass HC3 BH q<0.05 in the complete ten-comparison family"
  ),
  unresolved_alternative = "within-lineage state composition, ambient carryover, and shared multicellular response remain possible",
  claim_boundary = "same-atlas descriptive lineage enrichment relative to hepatocytes; not exclusive specificity, transcript origin, secretion, or cell autonomy"
)]
fwrite(forest, file.path(OUT, "source_tables/fig4d_compact_lineage_forest.tsv"), sep = "\t")

p_d <- ggplot(forest, aes(y = lineage_label)) +
  geom_vline(xintercept = 0, linewidth = 0.30, color = "#707070") +
  geom_segment(
    data = forest[target == FALSE],
    aes(x = hc3_ci_low, xend = hc3_ci_high, yend = lineage_label),
    linewidth = 0.35, color = "#8F8F8F"
  ) +
  geom_point(
    data = forest[target == FALSE], aes(x = beta),
    shape = 21, size = 1.65, stroke = 0.25, fill = "white", color = "#707070"
  ) +
  geom_segment(
    data = forest[target == TRUE],
    aes(x = hc3_ci_low, xend = hc3_ci_high, yend = lineage_label, color = program_key),
    linewidth = 0.55
  ) +
  geom_point(
    data = forest[target == TRUE],
    aes(x = beta, fill = program_key, shape = program_key),
    size = 2.30, stroke = 0.30, color = "black"
  ) +
  facet_wrap(
    ~bulk_plot_label, nrow = 1,
    labeller = labeller(bulk_plot_label = label_parsed)
  ) +
  scale_color_manual(values = setNames(hero$program_color, hero$program_name), guide = "none") +
  scale_fill_manual(values = setNames(hero$program_color, hero$program_name), guide = "none") +
  scale_shape_manual(values = setNames(hero$program_shape, hero$program_name), guide = "none") +
  scale_y_discrete(drop = TRUE) +
  scale_x_continuous(
    limits = c(-0.12, 2.68), breaks = c(0, 1, 2),
    expand = c(0, 0)
  ) +
  labs(x = "Program-score difference versus hepatocytes", y = NULL) +
  theme_panel() +
  theme(
    axis.ticks.y = element_blank(),
    panel.spacing.x = unit(0.16, "inches")
  )

bulk <- fread(bulk_path, na.strings = c("", "NA"))
check("bulk_rows", nrow(bulk) == 8L, nrow(bulk), 8L)
check("bulk_heroes", setequal(bulk$program_uid, hero$program_uid),
      paste(sort(unique(bulk$program_uid)), collapse = ";"),
      paste(hero$program_uid, collapse = ";"))
check("bulk_stages", setequal(bulk$stage, c("F1", "F2", "F3", "F4")),
      paste(sort(unique(bulk$stage)), collapse = ";"), "F1;F2;F3;F4")
bulk <- merge(bulk, hero[, .(program_uid, bulk_plot_label, program_color, program_shape)],
              by = "program_uid", all.x = TRUE, sort = FALSE)
bulk[, stage_label := factor(paste0(stage, " vs F0"),
                             levels = c("F1 vs F0", "F2 vs F0", "F3 vs F0", "F4 vs F0"))]
bulk[, program_row := factor(
  program_name,
  levels = c("ECM/IGFBP7", "Ductular-injury/BICC1")
)]
bulk[, display_role := "fixed_member_gene_stage_contrast"]
fwrite(bulk, file.path(OUT, "source_tables/fig4e_unconnected_stage_effects.tsv"), sep = "\t")

# Dot plot, not stem+ball: lollipops are barred by house style. Dropping the
# stem keeps both encodings intact (x = effect, point area = BH weight fraction),
# which a bar chart would cost.
p_e <- ggplot(bulk, aes(x = effect, y = program_row, color = program_name, fill = program_name)) +
  geom_point(
    aes(size = bh_weight_fraction, shape = program_name),
    stroke = 0.30, color = "black"
  ) +
  facet_grid(. ~ stage_label) +
  scale_color_manual(values = setNames(hero$program_color, hero$program_name), guide = "none") +
  scale_fill_manual(values = setNames(hero$program_color, hero$program_name), guide = "none") +
  scale_shape_manual(values = setNames(hero$program_shape, hero$program_name), guide = "none") +
  scale_size_area(
    limits = c(0, 1), max_size = 3.0,
    breaks = c(0.25, 0.50, 0.75), labels = percent_format(accuracy = 1),
    name = "BH-supported\nfixed weight"
  ) +
  scale_x_continuous(
    limits = c(0, 1.20), breaks = c(0, 0.6, 1.2), expand = c(0, 0)
  ) +
  scale_y_discrete(labels = c(
    "ECM/IGFBP7" = expression("ECM / " * italic(IGFBP7)),
    "Ductular-injury/BICC1" = expression("Ductular injury / " * italic(BICC1))
  )) +
  labs(x = expression("Fixed-program weighted " * log[2] * "FC"), y = NULL) +
  theme_panel() +
  theme(
    panel.spacing.x = unit(0.16, "inches"),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.height = unit(0.10, "inches"),
    legend.key.width = unit(0.12, "inches")
  ) +
  guides(size = guide_legend(
    nrow = 1, title.position = "left",
    override.aes = list(shape = 21, fill = "white", color = "#707070")
  ))

save_panel <- function(plot, path, width, height) {
  ggsave(path, plot, width = width, height = height, device = cairo_pdf)
}

save_panel(p_d, file.path(OUT, "panels/fig4d_lineage_enrichment_forest.pdf"), 4.35, 2.20)
save_panel(p_e, file.path(OUT, "panels/fig4e_independent_stage_contrasts.pdf"), 4.60, 2.05)

caption <- c(
  "Geometry-review candidate only. D, Dataset-adjusted donor-paired lineage-minus-hepatocyte program-score contrasts. Points are model estimates and intervals are HC3 95% confidence intervals. The colored point marks the predefined biologically matched lineage; BH correction spans all ten displayed contrasts. This is descriptive lineage enrichment, not exclusive specificity, transcript origin, secretion, or cell autonomy.",
  "E, Each unconnected point is a fixed L1-weighted member-gene effect for one cross-sectional fibrosis-stage-versus-F0 bulk contrast. Point area is the fraction of fixed program weight meeting BH q<0.05 in the complete 23,370-gene family for that contrast. The display is not a longitudinal trajectory and does not provide a program-level significance test."
)
writeLines(caption, file.path(OUT, "captions/fig4de_geometry_review.md"))

fwrite(rbindlist(check_rows), file.path(OUT, "validation.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
