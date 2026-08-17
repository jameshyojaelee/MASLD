#!/usr/bin/env Rscript

# Figure 5A: input firewall for the physical-context analyses.
# This is a provenance schematic, not a quantitative validation funnel.
# KEY MESSAGE: Fixed upstream inputs feed distinct assay-native Figure 5 branches without a global validation score.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
grDevices::pdf.options(useDingbats = FALSE)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FAM <- "Helvetica"
FS <- 6
INK <- "#202124"
MUTED <- "#6B7280"
GENE <- "#6A51A3"
PROGRAM <- "#C9265E"
PROTEIN <- "#00838F"
BOX_FILL <- "#FAFAFA"

branches <- data.table(
  y = c(2.45, 1.45, 0.45),
  source_title = c("Fig. 2 + Fig. 3", "Program registry", "Fixed before refit"),
  source_sub = c("prioritized genes", "prespecified display set", "25 proteins"),
  contract = c("gene-level context", "no reselection", "descriptive set"),
  readout = c(
    "Protein triage\nopen chromatin",
    "ATAC  ·  DIA-MS\nspatial",
    "Protein  ·  histology"
  ),
  panels = c("B · D", "E · F", "C"),
  colour = c(GENE, PROGRAM, PROTEIN)
)

write_sidecar <- branches[, .(
  branch = c("prioritized_gene_context", "prespecified_program_projection", "fixed_protein_display"),
  source_title,
  source_sub,
  contract,
  readout = gsub("\n", "; ", readout, fixed = TRUE),
  panels,
  interpretation = c(
    "context_not_global_validation",
    "assay_native_projection_without_reselection",
    "selection_conditioned_descriptive_display"
  )
)]
CANDIDATE_ROOT <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
OUTPUT_DIR <- if (nzchar(CANDIDATE_ROOT)) file.path(CANDIDATE_ROOT, "figure5", "panels") else file.path(FIG5_CONTEXT_DIR, "panels")
sidecar <- file.path(OUTPUT_DIR, "data", "fig5a_input_firewall.tsv")
dir.create(dirname(sidecar), recursive = TRUE, showWarnings = FALSE)
fwrite(write_sidecar, sidecar, sep = "\t", quote = FALSE)

p <- ggplot() +
  annotate("text", x = 0.62, y = 3.08, label = "FIXED INPUT", hjust = 0.5,
           family = FAM, size = FS/.pt, colour = MUTED) +
  annotate("text", x = 3.70, y = 3.08, label = "ASSAY-NATIVE READOUT", hjust = 0.5,
           family = FAM, size = FS/.pt, colour = MUTED) +
  annotate("segment", x = 2.00, xend = 2.00, y = 0.00, yend = 2.90,
           linewidth = 0.35, linetype = "22", colour = "#BDBDBD") +
  annotate("text", x = 2.00, y = 3.08, label = "PRESPECIFIED", hjust = 0.5,
           family = FAM, size = FS/.pt, colour = MUTED) +
  geom_curve(
    data = branches,
    aes(x = 1.24, y = y, xend = 2.53, yend = y, colour = colour),
    curvature = 0, linewidth = 0.55,
    arrow = grid::arrow(length = grid::unit(0.045, "in"), type = "closed"),
    show.legend = FALSE
  ) +
  geom_rect(
    data = branches,
    aes(xmin = 0.05, xmax = 1.22, ymin = y - 0.34, ymax = y + 0.34),
    fill = BOX_FILL, colour = branches$colour, linewidth = 0.55
  ) +
  geom_text(
    data = branches,
    aes(x = 0.635, y = y + 0.10, label = source_title),
    family = FAM, fontface = "plain", size = FS/.pt, colour = INK
  ) +
  geom_text(
    data = branches,
    aes(x = 0.635, y = y - 0.12, label = source_sub),
    family = FAM, size = FS/.pt, colour = INK
  ) +
  geom_text(
    data = branches,
    aes(x = 0.635, y = y - 0.25, label = contract),
    family = FAM, size = FS/.pt, colour = MUTED
  ) +
  geom_rect(
    data = branches,
    aes(xmin = 2.55, xmax = 4.70, ymin = y - 0.34, ymax = y + 0.34),
    fill = BOX_FILL, colour = "#C7C7C7", linewidth = 0.45
  ) +
  geom_segment(
    data = branches,
    aes(x = 2.55, xend = 2.55, y = y - 0.34, yend = y + 0.34, colour = colour),
    linewidth = 1.35, show.legend = FALSE
  ) +
  geom_text(
    data = branches,
    aes(x = 2.72, y = y, label = readout),
    hjust = 0, family = FAM, size = FS/.pt, colour = INK
  ) +
  geom_label(
    data = branches,
    aes(x = 4.44, y = y, label = panels),
    family = FAM, fontface = "plain", size = FS/.pt,
    fill = "white", colour = INK, linewidth = 0,
    label.padding = grid::unit(0.08, "lines")
  ) +
  coord_cartesian(xlim = c(0, 4.78), ylim = c(-0.02, 3.18), expand = FALSE, clip = "off") +
  scale_colour_identity() +
  theme_void() +
  theme(
    text = element_text(family = FAM, face = "plain"),
    plot.margin = margin(2, 2, 2, 2)
  )

out <- file.path(OUTPUT_DIR, "fig5a_input_firewall.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 3.04, height = 2.12, device = grDevices::cairo_pdf)
message("[fig5a] saved input firewall: ", out)
message(
  "CAPTION (Fig. 5A): Fixed-input analysis firewall. Figure 2/3 prioritized genes feed gene-level ",
  "protein and chromatin context; the prespecified program display set feeds assay-native ATAC, DIA-MS, ",
  "and spatial projections without reselection; and the fixed 25-protein set remains a separate, ",
  "selection-conditioned descriptive display. No global validation endpoint or cross-modal score is used."
)
