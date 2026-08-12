#!/usr/bin/env Rscript

# KEY MESSAGE: High-confidence disease-associated lncRNAs retain direction
# across cohorts and leave-one-cohort-out fits, while genomic context remains a
# descriptive annotation rather than a functional claim.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(digits = 17, scipen = 999)
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
  expected <- c("bulk", "cohort", "loco", "review", "output")
  parsed <- list()
  for (value in values) {
    pieces <- strsplit(value, "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !startsWith(pieces[[1L]], "--")) {
      fail("Arguments must use --name=value: ", value)
    }
    name <- sub("^--", "", pieces[[1L]])
    if (!name %in% expected || name %in% names(parsed)) fail("Unknown argument: ", name)
    parsed[[name]] <- pieces[[2L]]
  }
  missing <- setdiff(expected, names(parsed))
  if (length(missing)) fail("Missing arguments: ", paste(missing, collapse = ","))
  parsed
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(
  Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  mustWork = TRUE
)
theme_path <- normalizePath(file.path(project_root, "scripts/figures/publication_theme.R"), mustWork = TRUE)
source(theme_path, local = TRUE)

inputs <- vapply(c("bulk", "cohort", "loco", "review"), function(key) {
  normalizePath(arguments[[key]], mustWork = TRUE)
}, character(1))
review_text <- paste(readLines(inputs[["review"]], warn = FALSE), collapse = "\n")
if (!grepl('"status"[[:space:]]*:[[:space:]]*"pass"', review_text)) {
  fail("Bulk lncRNA candidate lacks an independent passing review")
}
output_parent <- normalizePath(dirname(arguments[["output"]]), mustWork = TRUE)
output <- file.path(output_parent, basename(arguments[["output"]]))
if (file.exists(output) || is_symlink(output)) fail("Refusing existing output")

bulk <- fread(inputs[["bulk"]])
cohort <- fread(inputs[["cohort"]])
loco <- fread(inputs[["loco"]])
high <- bulk[high_confidence == TRUE]
if (nrow(bulk) != 6249L || nrow(high) != 271L ||
    nrow(cohort) != 6249L * 5L || nrow(loco) != 6249L * 5L) {
  fail("Validated lncRNA family census drift")
}

setorder(high, logFC, gene_id_versioned)
high[, gene_order := seq_len(.N)]
high_ids <- high$gene_id_versioned
heatmap <- rbindlist(list(
  high[, .(
    gene_id_versioned, analysis = "Pooled", effect = logFC
  )],
  cohort[gene_id_versioned %in% high_ids, .(
    gene_id_versioned, analysis = cohort, effect = logFC
  )],
  loco[gene_id_versioned %in% high_ids, .(
    gene_id_versioned, analysis = paste0("Leave out ", excluded_cohort), effect = logFC
  )]
))
heatmap <- merge(heatmap, high[, .(gene_id_versioned, gene_order)], by = "gene_id_versioned")
cohorts <- sort(unique(cohort$cohort))
heat_levels <- c("Pooled", cohorts, paste0("Leave out ", cohorts))
heatmap[, analysis := factor(analysis, levels = heat_levels)]
heatmap[, effect_display := pmax(-2, pmin(2, effect))]
if (nrow(heatmap) != 271L * 11L || anyNA(heatmap$analysis)) fail("Heatmap family drift")

heat_plot <- ggplot(heatmap, aes(x = analysis, y = gene_order, fill = effect_display)) +
  geom_raster() +
  scale_fill_gradient2(
    low = masld_colors$down, mid = "white", high = masld_colors$up,
    midpoint = 0, limits = c(-2, 2), name = "log2 FC\n(clipped)"
  ) +
  scale_x_discrete(labels = function(x) sub("GSE", "", sub("Leave out ", "LOO ", x))) +
  labs(x = NULL, y = "271 high-confidence lncRNAs\n(ordered by pooled effect)") +
  theme_masld(base_size = 6) + theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), legend.position = "right"
  )

class_counts <- high[, .N, by = lncrna_genomic_class][order(N)]
class_counts[, class_label := gsub("_", " ", lncrna_genomic_class)]
class_counts[, class_label := factor(class_label, levels = class_label)]
class_plot <- ggplot(class_counts, aes(x = N, y = class_label)) +
  geom_col(fill = masld_colors$up, width = 0.65) +
  geom_text(aes(label = N), hjust = -0.15, size = GEOM_TEXT_6PT) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "High-confidence lncRNAs", y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

example_id <- "ENSG00000213062.6"
example <- rbindlist(list(
  bulk[gene_id_versioned == example_id, .(
    analysis = "Pooled", estimate = logFC, lower = logFC - 1.96 * SE,
    upper = logFC + 1.96 * SE
  )],
  cohort[gene_id_versioned == example_id, .(
    analysis = cohort, estimate = logFC, lower = CI_low, upper = CI_high
  )],
  loco[gene_id_versioned == example_id, .(
    analysis = paste0("Leave out ", excluded_cohort),
    estimate = logFC, lower = CI_low, upper = CI_high
  )]
))
if (nrow(example) != 11L) fail("Deterministic example family drift")
example[, analysis := factor(analysis, levels = rev(heat_levels))]
example[, source_type := fifelse(grepl("Leave out", analysis), "Leave-one-cohort-out",
                                 fifelse(analysis == "Pooled", "Pooled", "Individual cohort"))]
forest_plot <- ggplot(example, aes(x = estimate, y = analysis, color = source_type)) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, linewidth = 0.35) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c(
    "Pooled" = "black", "Individual cohort" = masld_colors$up,
    "Leave-one-cohort-out" = masld_colors$down
  )) +
  labs(
    x = "Disease versus control log2 fold change", y = NULL, color = NULL,
    title = example_id
  ) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(legend.position = "top", axis.line.y = element_blank(), axis.ticks.y = element_blank())

tmp <- file.path(output_parent, paste0(".", basename(output), ".tmp.", Sys.getpid()))
if (file.exists(tmp)) fail("Temporary output exists")
dir.create(tmp, recursive = FALSE)
ggsave(file.path(tmp, "figS_lncrna_robustness_heatmap.pdf"), heat_plot,
       width = 4.4, height = 3.2, device = cairo_pdf, units = "in")
ggsave(file.path(tmp, "figS_lncrna_genomic_classes.pdf"), class_plot,
       width = 3.35, height = 2.2, device = cairo_pdf, units = "in")
ggsave(file.path(tmp, "figS_lncrna_example_forest.pdf"), forest_plot,
       width = 4.0, height = 3.0, device = cairo_pdf, units = "in")
write_tsv_once(heatmap, file.path(tmp, "figS_lncrna_robustness_heatmap_source.tsv"))
write_tsv_once(class_counts, file.path(tmp, "figS_lncrna_genomic_classes_source.tsv"))
write_tsv_once(example, file.path(tmp, "figS_lncrna_example_forest_source.tsv"))
write_tsv_once(
  data.table(
    role = c("bulk", "cohort", "leave_one_cohort_out", "independent_review", "theme"),
    path = c(inputs, theme_path),
    size_bytes = file.info(c(inputs, theme_path))$size,
    sha256 = vapply(c(inputs, theme_path), sha256_file, character(1))
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
