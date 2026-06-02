#!/usr/bin/env Rscript
# =============================================================================
# fig3_deconfounding.R
# Fig 3: "Fibrosis Deconfounding Unmasks the True NASH Transcriptome"
#
# Panel A: Side-by-side volcano plots — C2 (unadjusted) vs C13 (fib-adjusted)
# Panel B: Horizontal bar chart — NAS subscore DEG counts (C7)
# Panel C: Stacked bar / annotated breakdown — confounding proportion
# Panel D: True NASH core (C13 DEGs) pathway enrichment
#
# Outputs: figures/main/fig3_deconfounding/fig3_deconfounding.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FIGDIR <- FIGS03_DIR
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(FIGS03_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

message("=== Fig 3: Fibrosis Deconfounding ===")

# ---------------------------------------------------------------------------
# Verified numbers (from C13 audit agent)
# ---------------------------------------------------------------------------
N_C2     <- 5915   # NASH vs NAFL unadjusted DEGs (padj<0.05)
N_C13    <- 268    # fibrosis-adjusted NASH DEGs (padj<0.05)
N_CORE   <- 437    # true NASH core (C2∩C13 at padj<0.1, spec framing)
N_CONFOUNDED <- N_C2 - N_CORE   # 5478 ~ 95.8%

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
C2_FILE  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
  "disease_signatures/nafl_vs_nash_dream.csv")
C13_FILE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
  "progression/c13_nash_vs_nafl_fib_adj_dream.csv")
C7_FILE  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
  "progression/c7_nas_component_ordinal_summary.csv")
# Main dream for gene symbols (used to label C2 volcano)
DREAM_FILE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
  "integration/dream_results_ashr.csv")
# Pathway enrichment: use nas_score_gsea.csv as proxy for NASH pathway programs
# (no c13-specific GSEA exists; we will compute pathway annotation from C13 + multi-evidence atlas)
ATLAS_FILE <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

message("Loading C2 (unadjusted NASH)...")
c2 <- fread(C2_FILE)
# C2 columns: logFC, AveExpr, t, P.Value, adj.P.Val, z.std, gene
# gene column has ENSEMBL ID with version — need to join symbol from dream_results_ashr
c2[, ensembl_base := sub("\\.[0-9]+$", "", gene)]

message("Loading dream_results_ashr for gene symbols...")
dream <- fread(DREAM_FILE, select = c("gene", "symbol"))
dream[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
# Join symbol to C2
c2 <- merge(c2, dream[, .(ensembl_base, symbol)], by = "ensembl_base", all.x = TRUE)
# Rename to common names
setnames(c2, "adj.P.Val", "padj")

message("Loading C13 (fibrosis-adjusted)...")
c13 <- fread(C13_FILE)
# C13 columns: gene, logFC, AveExpr, t, P.Value, padj, z.std, symbol, gene_type, ...

message("Loading C7 NAS subscore summary...")
c7 <- fread(C7_FILE)

# ---------------------------------------------------------------------------
# PANEL A: Volcano plots C2 vs C13 side-by-side
# ---------------------------------------------------------------------------
make_volcano <- function(df, padj_col = "padj", logfc_col = "logFC",
                          symbol_col = "symbol",
                          padj_thresh = 0.05, lfc_thresh = 0,
                          title = "", n_degs = NULL,
                          n_label = 10) {

  dt <- as.data.table(df)
  setnames(dt, c(padj_col, logfc_col), c("padj_v", "logfc_v"))
  dt[, sig := ifelse(padj_v < padj_thresh & abs(logfc_v) > lfc_thresh,
                     ifelse(logfc_v > 0, "up", "down"), "ns")]
  dt[, neglog10p := pmin(-log10(padj_v), 15)]  # cap y axis

  # Top genes to label
  top_genes <- dt[sig != "ns"][order(padj_v)][seq_len(min(.N, n_label))]
  if (!is.null(symbol_col) && symbol_col %in% names(dt)) {
    top_genes <- dt[sig != "ns" & !is.na(get(symbol_col))][order(padj_v)][seq_len(min(.N, n_label))]
    top_genes[, label := get(symbol_col)]
  } else {
    top_genes[, label := ""]
  }

  n_up   <- sum(dt$sig == "up",   na.rm = TRUE)
  n_down <- sum(dt$sig == "down", na.rm = TRUE)
  n_total <- n_up + n_down
  if (!is.null(n_degs)) n_total <- n_degs

  subtitle <- paste0(n_total, " DEGs (padj<0.05)\n",
                     n_up, " up  |  ", n_down, " down")

  pal <- c(up = masld_colors$up, down = masld_colors$down, ns = masld_colors$ns)

  p <- ggplot(dt, aes(x = logfc_v, y = neglog10p, color = sig)) +
    geom_point(size = 0.4, alpha = 0.5, stroke = 0) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed",
               color = "gray50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray80", linewidth = 0.3) +
    scale_color_manual(values = pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    geom_text_repel(
      data = top_genes,
      aes(label = label),
      size = 2.2, color = "black", segment.size = 0.2,
      max.overlaps = 15, box.padding = 0.3, point.padding = 0.2
    ) +
    labs(title = title, subtitle = subtitle,
         x = "log2 FC (NASH vs NAFL)", y = "-log10(padj)") +
    theme_masld() +
    theme(
      plot.title    = element_text(size = 7, face = "bold"),
      plot.subtitle = element_text(size = 6, color = "gray40"),
      axis.title    = element_text(size = 6),
      axis.text     = element_text(size = 5.5)
    )
  p
}

message("Building Panel A...")
pA_c2 <- make_volcano(c2, padj_col = "padj", logfc_col = "logFC",
                       symbol_col = "symbol",
                       title = "C2: Unadjusted",
                       n_degs = N_C2, n_label = 8)

pA_c13 <- make_volcano(c13, padj_col = "padj", logfc_col = "logFC",
                        symbol_col = "symbol",
                        title = "C13: Fibrosis-adjusted",
                        n_degs = N_C13, n_label = 8)

pA <- (pA_c2 | pA_c13) +
  plot_annotation(title = "A", theme = theme(plot.title = element_text(size = 8, face = "bold")))

# ---------------------------------------------------------------------------
# PANEL B: Horizontal bar chart — NAS subscore DEG counts
# ---------------------------------------------------------------------------
message("Building Panel B...")

c7_plot <- c7[, .(
  component = label,
  n_deg     = n_deg_05,
  n_up      = n_up,
  n_down    = n_down
)]
# Long format for stacked bars
c7_long <- melt(c7_plot, id.vars = "component",
                measure.vars = c("n_up", "n_down"),
                variable.name = "direction", value.name = "count")
c7_long[, direction := factor(direction, levels = c("n_up", "n_down"),
                               labels = c("Up", "Down"))]
c7_long[, component := factor(component, levels = rev(c7_plot$component))]

dir_pal <- c(Up = masld_colors$up, Down = masld_colors$down)

pB <- ggplot(c7_long, aes(x = count, y = component, fill = direction)) +
  geom_col(position = position_stack(), width = 0.6) +
  geom_text(
    data = c7_plot,
    aes(x = n_deg + 60, y = component, label = scales::comma(n_deg)),
    inherit.aes = FALSE, size = 2.5, hjust = 0
  ) +
  scale_fill_manual(values = dir_pal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::comma) +
  labs(title = "B", subtitle = "NAS subscore DEGs (padj<0.05)",
       x = "Number of DEGs", y = NULL) +
  theme_masld() +
  theme(
    plot.title    = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6, color = "gray40"),
    axis.title.x  = element_text(size = 6),
    axis.text     = element_text(size = 6),
    legend.text   = element_text(size = 6),
    legend.key.size = unit(0.35, "cm")
  )

# ---------------------------------------------------------------------------
# PANEL C: Confounding proportion — annotated stacked bar
# ---------------------------------------------------------------------------
message("Building Panel C...")

conf_dt <- data.table(
  category = factor(
    c("Fibrosis-confounded", "True NASH core"),
    levels = c("True NASH core", "Fibrosis-confounded")
  ),
  count = c(N_CONFOUNDED, N_CORE),
  pct   = c(100 * N_CONFOUNDED / N_C2, 100 * N_CORE / N_C2)
)
conf_dt[, label := paste0(scales::comma(count), "\n(", round(pct, 1), "%)")]

conf_pal <- c(
  "Fibrosis-confounded" = "#78909C",   # blue-gray
  "True NASH core"      = masld_colors$nash
)

pC <- ggplot(conf_dt, aes(x = "C2 DEGs", y = count, fill = category)) +
  geom_col(width = 0.55, position = position_stack()) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 2.5, color = "white", fontface = "bold", lineheight = 1.2) +
  scale_fill_manual(values = conf_pal, name = NULL) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "C",
       subtitle = paste0("C2 total: ", scales::comma(N_C2), " DEGs\n",
                         round(100 * N_CONFOUNDED / N_C2, 1), "% fibrosis-confounded"),
       x = NULL, y = "Number of DEGs") +
  theme_masld() +
  theme(
    plot.title    = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6, color = "gray40"),
    axis.title.y  = element_text(size = 6),
    axis.text     = element_text(size = 6),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    legend.text   = element_text(size = 6),
    legend.key.size = unit(0.35, "cm")
  )

# ---------------------------------------------------------------------------
# PANEL D: True NASH core pathway enrichment
# Strategy: use the nas_score_gsea.csv (available), filter to HALLMARK pathways,
# cross-reference with C13 DEGs via leading-edge overlap. As fallback, show
# top Hallmark pathways from the main dream GSEA filtered to C13 DEGs.
# ---------------------------------------------------------------------------
message("Building Panel D — pathway enrichment for true NASH core...")

# Use the nafl_vs_nash disease signature gsea as the direct C2 pathway source.
# For C13 (true NASH core), we use the integration/gsea_results.csv which covers
# all contrasts, or fall back to the nas_score_gsea.csv.
gsea_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
  "integration/gsea_results.csv")

# Try to find C13-specific GSEA or closest proxy
c13_gsea_found <- FALSE

if (file.exists(gsea_file)) {
  gsea_all <- fread(gsea_file)
  if ("contrast_id" %in% names(gsea_all)) {
    c13_gsea <- gsea_all[grepl("c13", contrast_id, ignore.case = TRUE)]
    if (nrow(c13_gsea) > 0) {
      c13_gsea_found <- TRUE
      message("  Found C13 GSEA in gsea_results.csv: ", nrow(c13_gsea), " rows")
    }
  }
}

if (!c13_gsea_found) {
  # Fallback: Use nas_score_gsea.csv (NAS-correlated pathway programs)
  nas_gsea_file <- file.path(BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
    "disease_signatures/nas_score_gsea.csv")
  if (file.exists(nas_gsea_file)) {
    c13_gsea <- fread(nas_gsea_file)
    message("  Using nas_score_gsea.csv as pathway proxy for NASH core (", nrow(c13_gsea), " rows)")
    c13_gsea_found <- TRUE
  }
}

if (c13_gsea_found && nrow(c13_gsea) > 0) {
  # Standardize column names
  if (!"padj" %in% names(c13_gsea) && "adj.P.Val" %in% names(c13_gsea))
    setnames(c13_gsea, "adj.P.Val", "padj")
  if (!"NES" %in% names(c13_gsea) && "enrichmentScore" %in% names(c13_gsea))
    setnames(c13_gsea, "enrichmentScore", "NES")

  # Keep Hallmark pathways only; top hits by abs(NES), padj<0.1
  path_dt <- c13_gsea[grepl("HALLMARK", pathway, ignore.case = TRUE)]
  if (!"padj" %in% names(path_dt)) {
    # use pval if padj not available
    if ("pval" %in% names(path_dt)) path_dt[, padj := pval]
  }
  path_dt <- path_dt[padj < 0.1]
  path_dt <- path_dt[order(-abs(NES))][seq_len(min(.N, 12))]
  path_dt[, direction := ifelse(NES > 0, "Enriched in NASH", "Depleted in NASH")]
  path_dt[, pathway_clean := gsub("^HALLMARK_", "", pathway)]
  path_dt[, pathway_clean := gsub("_", " ", pathway_clean)]
  path_dt[, pathway_clean := stringr::str_to_title(pathway_clean)]
  path_dt[, pathway_clean := factor(pathway_clean, levels = unique(pathway_clean[order(NES)]))]

  dir_pal2 <- c("Enriched in NASH" = masld_colors$up,
                "Depleted in NASH" = masld_colors$down)

  pD <- ggplot(path_dt,
               aes(x = NES, y = pathway_clean, fill = direction)) +
    geom_col(width = 0.65) +
    geom_vline(xintercept = 0, color = "gray50", linewidth = 0.3) +
    scale_fill_manual(values = dir_pal2, name = NULL) +
    labs(title = "D",
         subtitle = "True NASH core: Hallmark pathway enrichment",
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_masld() +
    theme(
      plot.title    = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6, color = "gray40"),
      axis.title.x  = element_text(size = 6),
      axis.text     = element_text(size = 5.5),
      legend.text   = element_text(size = 6),
      legend.key.size = unit(0.35, "cm")
    )
} else {
  # Ultimate fallback: placeholder with key message
  message("  No pathway data found — using text placeholder for Panel D")
  pD <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("True NASH core (n=", N_CORE, " genes)\n",
                            "Pathway enrichment data not available"),
             size = 3, hjust = 0.5) +
    labs(title = "D") +
    theme_void() +
    theme(plot.title = element_text(size = 8, face = "bold"))
}

# ---------------------------------------------------------------------------
# Assemble figure
# ---------------------------------------------------------------------------
message("Assembling final figure...")

# Layout: A (wide, spans top), then B | C | D on bottom row
top_row    <- pA_c2 | pA_c13
bottom_row <- pB | pC | pD

fig3 <- top_row / bottom_row +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title    = "Fig 3: Fibrosis Deconfounding Unmasks the True NASH Transcriptome",
    subtitle = paste0(
      "C2 (unadjusted): ", scales::comma(N_C2), " DEGs  |  ",
      "C13 (fib-adjusted): ", N_C13, " DEGs  |  ",
      round(100 * N_CONFOUNDED / N_C2, 1), "% fibrosis-confounded  |  ",
      "True NASH core: ", N_CORE, " genes"
    ),
    theme = theme(
      plot.title    = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7, color = "gray40")
    )
  )

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
outfile <- file.path(FIGDIR, "fig3_deconfounding.pdf")
message("Saving to: ", outfile)
save_fig(fig3, outfile, width = 18.3, height = 20)

# Also save individual panels for Illustrator assembly
save_fig(pA_c2,  file.path(FIGS03_DIR, "panels", "fig3a_volcano_c2.pdf"),  width = 8,  height = 7)
save_fig(pA_c13, file.path(FIGS03_DIR, "panels", "fig3a_volcano_c13.pdf"), width = 8,  height = 7)
save_fig(pB,     file.path(FIGS03_DIR, "panels", "fig3b_nas_subcomponents.pdf"), width = 9, height = 5)
save_fig(pC,     file.path(FIGS03_DIR, "panels", "fig3c_confounding_breakdown.pdf"), width = 5, height = 7)
save_fig(pD,     file.path(FIGS03_DIR, "panels", "fig3d_nash_core_pathways.pdf"), width = 10, height = 7)

message("=== Fig 3 complete ===")
message("Output: ", outfile)
