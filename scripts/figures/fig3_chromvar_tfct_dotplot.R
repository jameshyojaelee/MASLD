# Fig 3 panel — chromVAR TF × cell-type activity dot plot
# Curated MASLD regulatory landscape: 18 biologically-meaningful TFs grouped
# by function (hepatocyte identity / nuclear-receptor drug targets / metabolic
# / stress-inflammation / fibrogenic-EMT) across 4 MASLD-core cell types
# (Hepatocytes, Cholangiocytes, Fibroblasts, Macrophages).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# Donor-level limma table (per cell type): the pseudoreplication-corrected fix.
# Replaces the retired per-CELL Mann-Whitney table (chromvar_tf_activity.csv),
# whose 4,832 "sig" hits were fabricated by treating each cell as an independent
# replicate of an <=18-donor contrast. The donor-level limma yields 110 sig
# (adj.P.Val < 0.05) genome-wide, ALL in the Low_confidence cell type; none of
# the 4 displayed cell types carry a chromVAR-significant motif (see megareview
# A6.1). Columns: cell_type, TF, logFC, P.Value, adj.P.Val, n_donors.
CHROMVAR_CSV <- file.path(ATAC_DIR, "results/chromvar_v2/chromvar_limma_per_ct.csv")
OUT_PDF      <- file.path(FIGS05_DIR, "figS05_scatac_chromvar_dotplot.pdf")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Curated TF list (display name -> chromVAR motif name; case-handled)
# Organised top-to-bottom: hep identity -> NR drug targets -> metabolic ->
# stress/inflammation -> fibrogenic/EMT.
# -----------------------------------------------------------------------------
TF_GROUPS <- list(
  `Hepatocyte-related` =
    c(HNF4A  = "HNF4A",  HNF1A = "HNF1A",
      FOXA1  = "FOXA1",  FOXA2 = "FOXA2"),
  `NR drug targets` =
    c(THRB   = "THRB",   NR1H4 = "Nr1H4",
      PPARA  = "Ppara",  PPARG = "Pparg::Rxra",
      RORA   = "RORA",   RXRA  = "Rxra"),
  `Metabolic` =
    c(MLXIPL = "MLXIPL", `NR1H3 (LXRα)` = "Nr1h3", KLF15 = "KLF15"),
  `Stress / inflammation` =
    c(XBP1   = "XBP1",   ETS1  = "ETS1"),
  `Fibrogenic / EMT` =
    c(SNAI1  = "SNAI1",  ZEB1  = "ZEB1",  TCF4  = "TCF4")
)

BOLD_TFS <- c("THRB", "NR1H4", "RORA", "HNF4A", "PPARA", "PPARG")

CT_KEEP  <- c("Hepatocytes", "Cholangiocytes", "Fibroblasts", "Macrophages")

# Flatten into a lookup table: display_name | motif_name | group
tf_lookup <- rbindlist(lapply(names(TF_GROUPS), function(g) {
  v <- TF_GROUPS[[g]]
  data.table(display = names(v), motif = unname(v), group = g)
}))

# -----------------------------------------------------------------------------
# Load + filter
# -----------------------------------------------------------------------------
dt <- fread(CHROMVAR_CSV)
# Donor-level limma columns -> harmonise to the names used downstream.
setnames(dt, c("TF", "logFC", "adj.P.Val"),
             c("tf_name", "logFC_deviation", "padj"))
dt <- dt[cell_type %in% CT_KEEP]
dt <- dt[order(padj)][!duplicated(dt[, .(tf_name, cell_type)])]

plot_dt <- merge(dt, tf_lookup, by.x = "tf_name", by.y = "motif",
                 all.y = TRUE)
missing <- plot_dt[is.na(cell_type), unique(display)]
if (length(missing) > 0) {
  message("WARNING: TFs missing from chromVAR data: ",
          paste(missing, collapse = ", "))
}
plot_dt <- plot_dt[!is.na(cell_type)]

# Display-name order (top-to-bottom on y-axis = first group at top)
y_order <- unlist(lapply(TF_GROUPS, names), use.names = FALSE)
y_order <- intersect(y_order, plot_dt$display)
plot_dt[, display := factor(display, levels = rev(y_order))]

plot_dt[, cell_type := factor(cell_type, levels = CT_KEEP)]
plot_dt[, neglog10_padj := pmin(-log10(pmax(padj, 1e-300)), 6)]
plot_dt[, lfc_clip      := pmax(pmin(logFC_deviation, 1.5), -1.5)]
plot_dt[, is_sig        := padj < 0.05]

# -----------------------------------------------------------------------------
# Y-axis label aesthetics: bold for user-specified TFs (THRB, NR1H4, ...)
# -----------------------------------------------------------------------------
y_levels <- levels(plot_dt$display)
# strip any parenthetical for matching (e.g., "NR1H3 (LXRα)" -> "NR1H3")
y_stem   <- sub(" \\(.*", "", y_levels)
y_face   <- ifelse(y_stem %in% BOLD_TFS, "italic", "plain")

# -----------------------------------------------------------------------------
# Group separators (horizontal lines between TF functional groups)
# -----------------------------------------------------------------------------
group_of <- tf_lookup$group[match(y_levels, tf_lookup$display)]
# group changes between row i and i+1 in the y-axis order (rev: top->bottom)
breaks_between <- which(diff(as.integer(factor(group_of, levels = unique(group_of)))) != 0)
hline_y <- length(y_levels) - breaks_between + 0.5

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
p <- ggplot(plot_dt,
            aes(x = cell_type, y = display,
                colour = lfc_clip, size = neglog10_padj)) +
  geom_hline(yintercept = hline_y,
             colour = "gray85", linewidth = 0.3) +
  geom_point(data = plot_dt[is_sig == TRUE],
             shape = 21, fill = NA, colour = "black",
             stroke = 0.35, show.legend = FALSE) +
  geom_point(alpha = 0.95) +
  scale_colour_gradient2(
    low      = "#1565C0",
    mid      = "#9E9E9E",
    high     = "#C9265E",
    midpoint = 0,
    limits   = c(-1.5, 1.5),
    oob      = scales::squish,
    name     = "Motif deviation Δ (MASLD−CT)"
  ) +
  scale_size_continuous(
    range  = c(0.5, 3.6),
    limits = c(0, 6),
    breaks = c(2, 4, 6),
    labels = c("2", "4", expression("">=6)),
    name   = expression(-log[10]*" padj")
  ) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme_pub() +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1, vjust = 1,
                                    size = PUB_AXIS_TEXT, colour = "black"),
    axis.text.y      = element_text(size = PUB_AXIS_TEXT, face = y_face,
                                    colour = "black"),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin    = margin(l = 4, r = 0),
    legend.title     = element_text(size = PUB_LEGEND_TIT, face = "plain"),
    legend.text      = element_text(size = PUB_LEGEND),
    plot.margin      = margin(4, 6, 2, 4)
  ) +
  guides(
    colour = guide_colourbar(title.position = "top",
                             barwidth  = unit(0.22, "cm"),
                             barheight = unit(2.0, "cm"),
                             order = 1),
    size   = guide_legend(title.position = "top",
                          override.aes = list(colour = "#444444"),
                          order = 2)
  )

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
# RETIRED 2026-07-14 (chromVAR panel consolidation): the significance dimension this
# dotplot carried (donor-level limma adj.P.Val) is now folded into the single
# figS05_scatac_chromvar_celltype_progression heatmap as asterisks. Output suppressed.
# ggsave(OUT_PDF, p,
#        width = 4.0, height = 3.6,
#        device = cairo_pdf)

# -----------------------------------------------------------------------------
# Verification
# -----------------------------------------------------------------------------
n_tfs <- length(unique(plot_dt$display))
n_cts <- length(unique(plot_dt$cell_type))
message(sprintf("[fig3_chromvar_tfct_dotplot] wrote %s", OUT_PDF))
message(sprintf("  TFs in plot      : %d  (%s)",
                n_tfs, paste(rev(levels(plot_dt$display)), collapse = ", ")))
message(sprintf("  Cell types       : %d  (%s)",
                n_cts, paste(levels(plot_dt$cell_type), collapse = ", ")))
message(sprintf("  Bold TFs         : %s",
                paste(intersect(BOLD_TFS, y_stem), collapse = ", ")))

# Source data CSV (donor-level limma schema)
out_csv <- sub("\\.pdf$", ".csv", OUT_PDF)
fwrite(plot_dt[, .(display, group, cell_type,
                   logFC = logFC_deviation, P.Value, adj_P_Val = padj,
                   n_donors)],
       out_csv)
n_sig_plot <- plot_dt[is_sig == TRUE, .N]
message(sprintf("  Source CSV       : %s", out_csv))
message(sprintf("  Donor-level sig  : %d of %d dots (adj.P.Val < 0.05)",
                n_sig_plot, plot_dt[, .N]))
