#!/usr/bin/env Rscript

# KEY MESSAGE: The same complete gene-level models capture protein-coding and
# lncRNA remodeling without creating a separate lncRNA significance family.

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
  expected <- c("bulk-counts", "bulk-review", "stage-counts", "stage-review", "output")
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
  if (length(missing)) fail("Missing arguments: ", paste(missing, collapse = ","))
  parsed
}
review_passed <- function(path) {
  lines <- readLines(path, warn = FALSE)
  any(grepl('"status"[[:space:]]*:[[:space:]]*"pass"', lines))
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(
  Sys.getenv(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
  ),
  mustWork = TRUE
)
theme_path <- normalizePath(file.path(project_root, "scripts/figures/publication_theme.R"), mustWork = TRUE)
source(theme_path, local = TRUE)

bulk_counts_path <- normalizePath(arguments[["bulk-counts"]], mustWork = TRUE)
bulk_review_path <- normalizePath(arguments[["bulk-review"]], mustWork = TRUE)
stage_counts_path <- normalizePath(arguments[["stage-counts"]], mustWork = TRUE)
stage_review_path <- normalizePath(arguments[["stage-review"]], mustWork = TRUE)
output_parent <- normalizePath(dirname(arguments[["output"]]), mustWork = TRUE)
output <- file.path(output_parent, basename(arguments[["output"]]))
inputs <- c(
  bulk_counts_path, bulk_review_path, stage_counts_path, stage_review_path, theme_path
)
if (file.exists(output) || is_symlink(output)) fail("Refusing existing output: ", output)
if (any(vapply(inputs, is_symlink, logical(1)))) fail("Symlinked figure input is prohibited")
if (!review_passed(bulk_review_path) || !review_passed(stage_review_path)) {
  fail("Figure inputs lack independent passing reviews")
}

bulk <- fread(bulk_counts_path)
required_bulk <- c(
  "display_biotype", "n_tested", "n_treat_positive", "n_treat_up", "n_treat_down"
)
if (!identical(names(bulk), required_bulk) || nrow(bulk) != 3L ||
    sum(bulk$n_tested) != 23370L || sum(bulk$n_treat_positive) != 1616L ||
    any(bulk$n_treat_positive != bulk$n_treat_up + bulk$n_treat_down)) {
  fail("Validated pooled biotype count schema or census drift")
}
stage <- fread(stage_counts_path)
required_stage <- c(
  "transition", "contrast", "display_biotype", "n_tested",
  "n_fdr_positive", "n_up", "n_down"
)
if (!identical(names(stage), required_stage) || nrow(stage) != 12L ||
    any(stage$n_fdr_positive != stage$n_up + stage$n_down) ||
    any(stage[, .N, by = transition]$N != 3L)) {
  fail("Validated stage biotype count schema or census drift")
}

biotype_levels <- c("protein_coding", "lncRNA", "other")
biotype_labels <- c(
  protein_coding = "Protein-coding", lncRNA = "lncRNA", other = "Other"
)
direction_colors <- c(Down = masld_colors$down, Up = masld_colors$up)

bulk_long <- rbindlist(list(
  bulk[, .(
    display_biotype, n_tested, n_significant = n_treat_down,
    signed_count = -n_treat_down, direction = "Down"
  )],
  bulk[, .(
    display_biotype, n_tested, n_significant = n_treat_up,
    signed_count = n_treat_up, direction = "Up"
  )]
))
bulk_long[, display_biotype := factor(
  display_biotype, levels = rev(biotype_levels), labels = rev(biotype_labels)
)]
bulk_long[, count_label := format(n_significant, big.mark = ",", scientific = FALSE)]
bulk_long[, label_hjust := fifelse(signed_count < 0, 1.12, -0.12)]

bulk_plot <- ggplot(bulk_long, aes(x = signed_count, y = display_biotype, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black") +
  geom_col(width = 0.58, color = NA) +
  geom_text(
    aes(label = count_label, hjust = label_hjust),
    size = GEOM_TEXT_6PT, color = "black", show.legend = FALSE
  ) +
  scale_fill_manual(values = direction_colors, breaks = c("Down", "Up")) +
  scale_x_continuous(
    labels = function(values) format(abs(values), big.mark = ",", scientific = FALSE),
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  labs(
    x = "TREAT-positive genes (down ←  → up)",
    y = NULL,
    fill = NULL
  ) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

transition_levels <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
transition_labels <- c(
  F0_to_F1 = "F0→F1", F1_to_F2 = "F1→F2",
  F2_to_F3 = "F2→F3", F3_to_F4 = "F3→F4"
)
stage_long <- rbindlist(list(
  stage[, .(
    transition, contrast, display_biotype, n_tested,
    n_significant = n_down, signed_count = -n_down, direction = "Down"
  )],
  stage[, .(
    transition, contrast, display_biotype, n_tested,
    n_significant = n_up, signed_count = n_up, direction = "Up"
  )]
))
stage_long[, transition := factor(
  transition, levels = rev(transition_levels), labels = rev(transition_labels)
)]
stage_long[, display_biotype := factor(
  display_biotype, levels = biotype_levels, labels = biotype_labels
)]
stage_long[, count_label := format(n_significant, big.mark = ",", scientific = FALSE)]
stage_long[, label_hjust := fifelse(signed_count < 0, 1.10, -0.10)]

stage_plot <- ggplot(stage_long, aes(x = signed_count, y = transition, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black") +
  geom_col(width = 0.62, color = NA) +
  geom_text(
    aes(label = count_label, hjust = label_hjust),
    size = GEOM_TEXT_6PT, color = "black", show.legend = FALSE
  ) +
  facet_wrap(~display_biotype, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = direction_colors, breaks = c("Down", "Up")) +
  scale_x_continuous(
    labels = function(values) format(abs(values), big.mark = ",", scientific = FALSE),
    expand = expansion(mult = c(0.22, 0.22))
  ) +
  labs(
    x = "Genes at FDR < 0.05 (down ←  → up; separate scales)",
    y = "Cross-sectional Kleiner contrast",
    fill = NULL
  ) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.x = unit(0.3, "cm")
  )

tmp <- file.path(
  output_parent,
  paste0(".", basename(output), ".tmp.", Sys.getenv("SLURM_JOB_ID", Sys.getpid()))
)
if (file.exists(tmp) || is_symlink(tmp)) fail("Refusing existing temporary output: ", tmp)
dir.create(tmp, recursive = FALSE)
if (!dir.exists(tmp)) fail("Failed to create temporary figure output")

bulk_pdf <- file.path(tmp, "fig3_pooled_deg_biotype.pdf")
stage_pdf <- file.path(tmp, "fig3_stage_deg_biotype.pdf")
ggsave(
  bulk_pdf, bulk_plot, width = 3.35, height = 1.75,
  device = cairo_pdf, units = "in"
)
ggsave(
  stage_pdf, stage_plot, width = 5.2, height = 2.25,
  device = cairo_pdf, units = "in"
)
write_tsv_once(bulk_long, file.path(tmp, "fig3_pooled_deg_biotype_source.tsv"))
write_tsv_once(stage_long, file.path(tmp, "fig3_stage_deg_biotype_source.tsv"))
write_tsv_once(
  data.table(
    input_role = c(
      "pooled_biotype_counts", "pooled_independent_review",
      "stage_biotype_counts", "stage_independent_review", "publication_theme"
    ),
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
if (!file.rename(tmp, output)) fail("Candidate figure publication failed: ", tmp)
cat("NONCODING_FIG3_PANELS_COMPLETE\t", output, "\n", sep = "")
