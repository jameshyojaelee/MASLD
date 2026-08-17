#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
})

SCRIPT_PATH <- normalizePath(sub(
  "^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
), mustWork = TRUE)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID
)
EXPECTED <- normalizePath(CANDIDATE, mustWork = TRUE)
if (!identical(normalizePath(CANDIDATE, mustWork = TRUE), EXPECTED)) {
  stop("Unsafe candidate root")
}

READY <- file.path(CANDIDATE, "integration/INTEGRATION_READY")
if (!file.exists(READY)) stop("INTEGRATION_READY is required")
ready <- fread(READY)
if (nrow(ready) == 0L || any(ready$status != "READY") ||
    any(ready$release_id != RELEASE_ID)) {
  stop("Invalid INTEGRATION_READY seal")
}
for (i in seq_len(nrow(ready))) {
  artifact <- if (startsWith(ready$artifact[i], "Analysis/")) {
    file.path(BASE, ready$artifact[i])
  } else {
    file.path(CANDIDATE, ready$artifact[i])
  }
  if (!file.exists(artifact)) stop("Missing integration artifact: ", artifact)
  observed <- digest(artifact, algo = "sha256", file = TRUE, serialize = FALSE)
  if (!identical(observed, ready$sha256[i])) stop("Integration hash drift: ", artifact)
}

OUT <- file.path(CANDIDATE, "integration/figures")
if (dir.exists(OUT)) stop("Refusing to overwrite Figure 4F integration directory")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
grDevices::pdf.options(useDingbats = FALSE)
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))

states <- fread(file.path(CANDIDATE, "integration/fig4f_atac_v3_state_contract.tsv"))
if (nrow(states) != 234L || uniqueN(states, by = c("cohort", "program_uid")) != 234L) {
  stop("Figure 4F state contract must contain 234 unique cohort-program rows")
}

coverage <- states[, .(
  numerator = sum(n_promoter_measured_genes),
  denominator = sum(n_program_genes),
  proportion = sum(n_promoter_measured_genes) / sum(n_program_genes),
  n_programs = .N
), by = .(cohort, lineage)]
coverage[, `:=`(
  metric = "Promoter-gene coverage",
  note = "summed measured genes / summed frozen program memberships",
  display_label = sprintf("%.0f%%\n%d programs", 100 * proportion, n_programs)
)]

score <- states[, .(
  numerator = sum(program_score_state == "testable"),
  denominator = .N,
  proportion = mean(program_score_state == "testable"),
  n_programs = .N
), by = .(cohort, lineage)]
score[, `:=`(
  metric = "Program score testable",
  note = "at least 8 genes and at least 20% frozen L1 weight retained",
  display_label = paste0(numerator, "/", denominator)
)]

contrast <- states[, .(
  numerator = sum(contrast_state == "testable"),
  denominator = .N,
  proportion = mean(contrast_state == "testable"),
  n_programs = .N
), by = .(cohort, lineage)]
contrast[, `:=`(
  metric = "MASH contrast testable",
  note = "program score testable and at least 4 eligible donors per condition",
  display_label = paste0(numerator, "/", denominator)
)]

source_table <- rbindlist(list(coverage, score, contrast), use.names = TRUE)
source_table[, release_id := RELEASE_ID]
setcolorder(source_table, c(
  "release_id", "cohort", "lineage", "metric", "numerator", "denominator",
  "proportion", "n_programs", "display_label", "note"
))
setorder(source_table, cohort, lineage, metric)
source_path <- file.path(OUT, "fig4f_atac_v3_contract_source.tsv")
fwrite(source_table, source_path, sep = "\t")

source_table[, lineage := factor(
  lineage,
  levels = rev(c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")),
  labels = rev(c("Hepatocyte", "Stellate", "Macrophage", "Cholangiocyte", "T/NK"))
)]
source_table[, metric := factor(
  metric,
  levels = c("Promoter-gene coverage", "Program score testable", "MASH contrast testable")
)]

theme_v3 <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      plot.title = element_text(size = 6, face = "plain"),
      plot.subtitle = element_text(size = 6, face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, face = "plain", color = "black"),
      strip.text = element_text(size = 6, face = "plain"),
      strip.background = element_blank(),
      legend.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.position = "bottom"
    )
}

plot <- ggplot(source_table, aes(metric, lineage, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = display_label), size = 5.5 / .pt, lineheight = 0.85) +
  facet_wrap(~ cohort, nrow = 1L) +
  scale_fill_gradientn(
    colors = purple_gradient,
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "Coverage / testable"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "ATAC context for all 117 frozen programs",
    subtitle = "GSE281367 T/NK uses a combined NK/T source label and remains source-dependent"
  ) +
  theme_v3() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

pdf_path <- file.path(OUT, "fig4f_atac_v3_contract.pdf")
device <- if (capabilities("cairo")) cairo_pdf else pdf
ggsave(pdf_path, plot, width = 5.6, height = 3.2, device = device)
session_path <- file.path(OUT, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_path)

manifest <- data.table(
  release_id = RELEASE_ID,
  panel = "fig4f_atac_v3_contract",
  pdf = basename(pdf_path),
  source_table = basename(source_path),
  pdf_sha256 = digest(pdf_path, algo = "sha256", file = TRUE, serialize = FALSE),
  source_sha256 = digest(source_path, algo = "sha256", file = TRUE, serialize = FALSE),
  renderer_sha256 = digest(SCRIPT_PATH, algo = "sha256", file = TRUE, serialize = FALSE)
)
fwrite(manifest, file.path(OUT, "figure_manifest.tsv"), sep = "\t")
message("[FIGURE] Wrote explicit candidate Figure 4F ATAC contract: ", pdf_path)
