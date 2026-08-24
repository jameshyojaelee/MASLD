#!/usr/bin/env Rscript
# KEY MESSAGE: Two fixed single-cell program member sets show increasing
# cross-sectional bulk fibrosis-stage effects without a program-significance
# encoding.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

grDevices::pdf.options(useDingbats = FALSE)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT <- Sys.getenv(
  "FIG4E_UNIFORM_OUT",
  unset = file.path(BASE, "figures/candidates/fig4e-uniform-points-2026-08-17-v1")
)
if (dir.exists(OUT)) stop("Candidate output already exists: ", OUT, call. = FALSE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))
dir.create(file.path(OUT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "source_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "captions"), recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(value)) stop("sha256sum failed: ", path, call. = FALSE)
  strsplit(value[[1]], "[[:space:]]+")[[1]][[1]]
}

source_path <- file.path(
  BASE,
  "figures/main/fig4_singlecell_programs/source_tables/current_candidate/fig4e_bulk_tissue_state_transport.tsv"
)
theme_path <- file.path(BASE, "scripts/figures/publication_theme.R")
inputs <- data.table(
  role = c("bulk_transport_input", "theme"),
  path = c(source_path, theme_path)
)
if (!all(file.exists(inputs$path))) stop("missing input", call. = FALSE)
inputs[, `:=`(sha256 = vapply(path, sha256, character(1)), bytes = file.info(path)$size)]
fwrite(inputs, file.path(OUT, "input_manifest.tsv"), sep = "\t")

d <- fread(source_path, na.strings = c("", "NA"))
stopifnot(
  nrow(d) == 8L,
  setequal(d$stage, c("F1", "F2", "F3", "F4")),
  uniqueN(d$program_uid) == 2L,
  all(is.finite(d$effect)),
  all(c(
    "biological_unit", "multiple_testing_family", "evidence_state",
    "state_reason", "unresolved_alternative", "claim_boundary"
  ) %in% names(d))
)
if (!file.copy(source_path, file.path(OUT, "source_tables/fig4e_bulk_tissue_state_transport.tsv"))) {
  stop("source-table copy failed", call. = FALSE)
}

program_colors <- c(
  "ECM/IGFBP7" = "#FF6F00",
  "Ductular-injury/BICC1" = "#00BCD4"
)
program_shapes <- c(
  "ECM/IGFBP7" = 21L,
  "Ductular-injury/BICC1" = 22L
)

d[, stage_index := match(stage, c("F1", "F2", "F3", "F4"))]
d[, program_key := fifelse(
  program_name == "ECM/IGFBP7", "ECM/IGFBP7", "Ductular-injury/BICC1"
)]
d[, gene_label := fifelse(
  program_key == "ECM/IGFBP7", "italic(IGFBP7)", "italic(BICC1)"
)]
end_labels <- d[stage == "F4"]
end_labels[, label_y := effect + fifelse(program_key == "ECM/IGFBP7", -0.015, 0.015)]

theme_panel <- function() {
  theme_masld(base_size = 6) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(3, 4, 3, 3)
    )
}

p <- ggplot(
  d,
  aes(
    x = stage_index, y = effect, group = program_key,
    color = program_key, shape = program_key
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#BDBDBD") +
  geom_line(linewidth = 0.65) +
  geom_point(
    aes(fill = program_key),
    size = 2.25, stroke = 0.3, color = "black"
  ) +
  geom_text(
    data = end_labels,
    aes(x = 4.08, y = label_y, label = gene_label),
    inherit.aes = FALSE,
    parse = TRUE, hjust = 0, size = 6 / ggplot2::.pt, family = "Helvetica"
  ) +
  scale_color_manual(values = program_colors, guide = "none") +
  scale_fill_manual(values = program_colors, guide = "none") +
  scale_shape_manual(values = program_shapes, guide = "none") +
  scale_x_continuous(
    breaks = 1:4, labels = c("F1", "F2", "F3", "F4"),
    limits = c(0.85, 4.65), expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.2, 0.3), limits = c(-0.03, 1.20),
    expand = c(0, 0)
  ) +
  labs(
    x = "Fibrosis stage versus F0",
    y = expression("Bulk member-gene " * log[2] * "FC")
  ) +
  theme_panel()

panel_path <- file.path(OUT, "panels/fig4e_bulk_tissue_state_transport.pdf")
ggsave(panel_path, p, width = 2.72, height = 2.43, device = cairo_pdf)

writeLines(
  paste0(
    "Figure 4E candidate. Fixed L1-weighted member-gene bulk log2 fold changes ",
    "for F1–F4 versus F0. Uniform points identify the eight descriptive ",
    "program-by-contrast summaries; no point-size or program-significance ",
    "encoding is shown. This is cross-sectional tissue-state transport, not ",
    "longitudinal within-person change or a program-level significance test."
  ),
  file.path(OUT, "captions/fig4e.md")
)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
