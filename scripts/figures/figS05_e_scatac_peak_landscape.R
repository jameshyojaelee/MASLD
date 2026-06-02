#!/usr/bin/env Rscript
# figS05_e_scatac_peak_landscape.R
# Per-cell-type scATAC peak landscape:
#   Panel A: horizontal stacked bar — N peaks (promoter / genic / distal)
#   Panel B: horizontal bar — N significant DA peaks F0_vs_F4 (up / down)
#
# Inputs:
#   Analysis/ATAC/Human_Multiome/results/snapatac2/peak_landscape_per_ct.tsv
#   Analysis/ATAC/Human_Multiome/results/stage_da/da_F0_vs_F4_<CT>.csv
# Output:
#   figures/supplementary/figS05_epigenomic_spatial/figS05_e_scatac_peak_landscape.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

LANDSCAPE_TSV <- file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/snapatac2/peak_landscape_per_ct.tsv")
DA_DIR        <- file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/stage_da")
PADJ          <- 0.05

# Map full peak-BED cell-type names to the 5-letter DA acronyms used in DA files
ct_da_lookup <- c(
  "Hepatocytes"        = "Hep",
  "Macrophages"        = "Mac",
  "Fibroblasts"        = "Fib",
  "Endothelial_cells"  = "Endo",
  "Cholangiocytes"     = "Chol",
  "T_cells"            = NA,
  "Plasma_cells"       = NA,
  "Resident_NK"        = NA,
  "Circulating_NK_NKT" = NA
)

# Display labels (collapse underscores)
ct_display <- c(
  "Hepatocytes"        = "Hepatocytes",
  "Macrophages"        = "Macrophages",
  "Fibroblasts"        = "Fibroblasts",
  "Endothelial_cells"  = "Endothelial",
  "Cholangiocytes"     = "Cholangiocytes",
  "T_cells"            = "T cells",
  "Plasma_cells"       = "Plasma cells",
  "Resident_NK"        = "Resident NK",
  "Circulating_NK_NKT" = "Circulating NK/NKT"
)

# ── Panel A data: peak landscape ──────────────────────────────────────────────
land <- fread(LANDSCAPE_TSV)
land[, ct_label := ct_display[cell_type]]
# Order CTs by total peaks descending (Hep largest typically)
setorder(land, -n_peaks_total)
ct_order <- land$ct_label

land_long <- melt(
  land,
  id.vars      = c("cell_type", "ct_label", "n_peaks_total"),
  measure.vars = c("n_promoter", "n_genic", "n_distal"),
  variable.name = "annotation",
  value.name    = "n"
)
land_long[, annotation := factor(
  annotation,
  levels = c("n_promoter", "n_genic", "n_distal"),
  labels = c("Promoter", "Genic", "Distal")
)]
land_long[, ct_label := factor(ct_label, levels = rev(ct_order))]

# Label position: right edge of bar
land[, ct_label_f := factor(ct_label, levels = rev(ct_order))]

panel_a <- ggplot(land_long, aes(y = ct_label, x = n / 1000, fill = annotation)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    name   = NULL,
    values = c(Promoter = "#9E9E9E",
               Genic    = "#1565C0",
               Distal   = "#C9265E")
  ) +
  scale_x_continuous(
    name   = "Peaks (thousands)",
    labels = comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(tag = "a", y = NULL) +
  theme_masld(base_size = 9) +
  theme(
    axis.title.x    = element_text(size = 9, face = "bold"),
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 7),
    plot.tag        = element_text(size = 11, face = "bold"),
    plot.tag.position = c(0.02, 0.97),
    legend.position = "bottom",
    legend.text     = element_text(size = 7),
    legend.key.size = unit(0.30, "cm"),
    legend.margin   = margin(t = -2),
    panel.grid      = element_blank()
  )

# ── Panel B data: F0_vs_F4 DA peaks ───────────────────────────────────────────
da_rows <- rbindlist(lapply(names(ct_da_lookup), function(ct_full) {
  ct_da <- ct_da_lookup[[ct_full]]
  if (is.na(ct_da)) return(NULL)
  f <- file.path(DA_DIR, sprintf("da_F0_vs_F4_%s.csv", ct_da))
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(
    fread(f, select = c("log2FC", "padj")),
    error = function(e) NULL
  )
  if (is.null(d) || !"padj" %in% names(d)) return(NULL)
  d <- d[!is.na(padj) & padj < PADJ]
  data.table(
    cell_type = ct_full,
    ct_label  = ct_display[[ct_full]],
    n_up      = sum(d$log2FC > 0),
    n_down    = sum(d$log2FC < 0)
  )
}), fill = TRUE)

da_long <- melt(
  da_rows,
  id.vars       = c("cell_type", "ct_label"),
  measure.vars  = c("n_up", "n_down"),
  variable.name = "direction",
  value.name    = "n"
)
da_long[, direction := factor(
  direction,
  levels = c("n_up", "n_down"),
  labels = c("Opening (up)", "Closing (down)")
)]
# Keep same CT ordering as panel A (so visually consistent), but only show DA-available CTs
da_ct_order <- ct_order[ct_order %in% da_rows$ct_label]
da_long[, ct_label := factor(ct_label, levels = rev(da_ct_order))]

# Totals for label
da_tot <- da_rows[, .(total = n_up + n_down), by = ct_label]
da_tot[, ct_label := factor(ct_label, levels = rev(da_ct_order))]

panel_b <- ggplot(da_long, aes(y = ct_label, x = n, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    name   = NULL,
    values = c("Opening (up)"   = "#C9265E",
               "Closing (down)" = "#1565C0")
  ) +
  scale_x_continuous(
    name   = "DA peaks (F0 vs F4)",
    labels = comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(tag = "b", y = NULL) +
  theme_masld(base_size = 9) +
  theme(
    axis.title.x    = element_text(size = 9, face = "bold"),
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 7),
    plot.tag        = element_text(size = 11, face = "bold"),
    plot.tag.position = c(0.02, 0.97),
    legend.position = "bottom",
    legend.text     = element_text(size = 7),
    legend.key.size = unit(0.30, "cm"),
    legend.margin   = margin(t = -2),
    panel.grid      = element_blank()
  )

# ── Compose ───────────────────────────────────────────────────────────────────
combo <- (panel_a | panel_b) + plot_layout(widths = c(1, 1))

out_pdf <- file.path(FIGS05_DIR, "figS05_e_scatac_peak_landscape.pdf")
ggsave(out_pdf, combo, width = 6.0, height = 3.0, device = cairo_pdf)
message("Wrote: ", out_pdf)

# Save underlying data
out_csv_a <- sub("\\.pdf$", "_panelA.csv", out_pdf)
out_csv_b <- sub("\\.pdf$", "_panelB.csv", out_pdf)
fwrite(land, out_csv_a)
fwrite(da_rows, out_csv_b)
message("Wrote: ", out_csv_a)
message("Wrote: ", out_csv_b)
