#!/usr/bin/env Rscript

# KEY MESSAGE: ATAC measurement and disease-contrast testability are distinct
# states, and dynamic accessibility replication is reported across all peaks.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE <- Sys.getenv(
  "ATAC_V3_CANDIDATE_ROOT",
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID)
)
EXPECTED <- normalizePath(
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID),
  mustWork = FALSE
)
if (!identical(normalizePath(CANDIDATE, mustWork = FALSE), EXPECTED)) stop("Unsafe candidate root")
OUT <- file.path(CANDIDATE, "figures")
if (dir.exists(OUT)) stop("Refusing to overwrite candidate figure directory")
dir.create(file.path(OUT, "source_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "panels"), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf.options(useDingbats = FALSE)
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))

theme_v3 <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      plot.title = element_text(size = 6, face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, face = "plain", color = "black"),
      legend.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      strip.text = element_text(size = 6, face = "plain"),
      strip.background = element_blank(),
      legend.position = "bottom"
    )
}
save_panel <- function(plot, path, width, height) {
  device <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(path, plot, width = width, height = height, device = device)
}

programs <- fread(file.path(CANDIDATE, "programs/program_atac_results.tsv"))
programs[, display_state := fcase(
  !program_score_testable, "not measured adequately",
  !contrast_testable, "measured; contrast untestable",
  qvalue < 0.05 & sensitivity_sign_agree, "MASH-associated",
  default = "measured; indeterminate"
)]
program_source <- programs[, .(n_programs = .N), by = .(cohort, lineage, display_state)]
all_grid <- CJ(
  cohort = unique(programs$cohort),
  lineage = c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"),
  display_state = c(
    "not measured adequately", "measured; contrast untestable",
    "measured; indeterminate", "MASH-associated"
  )
)
program_source <- merge(all_grid, program_source, by = c("cohort", "lineage", "display_state"), all.x = TRUE)
program_source[is.na(n_programs), n_programs := 0L]
program_source[, denominator := sum(n_programs), by = .(cohort, lineage)]
setorder(program_source, cohort, lineage, display_state)
program_source_path <- file.path(OUT, "source_tables", "fig4f_atac_state_counts.tsv")
fwrite(program_source, program_source_path, sep = "\t")

program_source[, lineage := factor(
  lineage,
  levels = rev(c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")),
  labels = rev(c("Hepatocyte", "Stellate", "Macrophage", "Cholangiocyte", "T/NK"))
)]
program_source[, display_state := factor(
  display_state,
  levels = c(
    "not measured adequately", "measured; contrast untestable",
    "measured; indeterminate", "MASH-associated"
  )
)]
p_program <- ggplot(program_source, aes(display_state, lineage, fill = n_programs)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = paste0(n_programs, "/", denominator)), size = 6 / .pt) +
  facet_wrap(~ cohort, nrow = 1L) +
  scale_fill_gradientn(colors = purple_gradient, name = "Programs") +
  labs(x = NULL, y = NULL, title = "Frozen-program ATAC observability and contrast state") +
  theme_v3() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
program_pdf <- file.path(OUT, "panels", "fig4f_atac_context.pdf")
save_panel(p_program, program_pdf, 5.4, 2.6)

da <- fread(file.path(CANDIDATE, "da/da_peak_results.tsv.gz"))
panel_records <- list(list(pdf = program_pdf, source = program_source_path, panel = "fig4f_atac_context"))
for (lineage_value in c("hepatocyte", "stellate", "macrophage")) {
  source <- da[lineage == lineage_value]
  source_path <- file.path(OUT, "source_tables", paste0("supp_da_effects_", lineage_value, ".tsv.gz"))
  fwrite(source, source_path, sep = "\t")
  plot <- ggplot(source, aes(logFC_gse244832, logFC_gse281367, color = evidence_state)) +
    geom_hline(yintercept = 0, color = "#BDBDBD", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "#BDBDBD", linewidth = 0.25) +
    ggrastr::rasterise(geom_point(size = 0.35, alpha = 0.45), dpi = 300) +
    scale_color_manual(values = c(
      supported = "#00695C", discordant = "#C9265E",
      source_dependent = "#7B1FA2", indeterminate = "#9E9E9E"
    ), drop = FALSE) +
    coord_equal() +
    labs(
      x = "GSE244832 MASH-normal log2 accessibility",
      y = "GSE281367 MASH-normal log2 accessibility",
      color = NULL,
      title = paste("Cross-cohort peak effects:", lineage_value)
    ) +
    theme_v3()
  pdf_path <- file.path(OUT, "panels", paste0("supp_da_effects_", lineage_value, ".pdf"))
  save_panel(plot, pdf_path, 3.5, 3.4)
  panel_records[[length(panel_records) + 1L]] <- list(
    pdf = pdf_path, source = source_path, panel = paste0("supp_da_effects_", lineage_value)
  )
}

manifest <- rbindlist(lapply(panel_records, function(value) data.table(
  release_id = RELEASE_ID,
  panel = value$panel,
  pdf = substring(value$pdf, nchar(CANDIDATE) + 2L),
  source_table = substring(value$source, nchar(CANDIDATE) + 2L),
  pdf_sha256 = digest(value$pdf, algo = "sha256", file = TRUE, serialize = FALSE),
  source_sha256 = digest(value$source, algo = "sha256", file = TRUE, serialize = FALSE)
)))
fwrite(manifest, file.path(OUT, "figure_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("[FIGURE] Wrote candidate-only non-genetic panels: ", OUT)
