##############################################################################
# DEG Threshold Landscape — pairwise similarity panels
#
# Two pairwise size-ratio matrices (no canonical anchor):
#   (a2) Across padj thresholds at fixed |LFC| = 0.5
#   (a3) Across |LFC| thresholds at fixed padj = 0.05
#
# Cell (i,j) = min(|A|,|B|) / max(|A|,|B|): what fraction of the larger
# list does the smaller list represent? Avoids the Jaccard triviality for
# nested sets (Jaccard dominated by set size; overlap coefficient = 1
# everywhere since stricter sets are always subsets of looser ones).
# Diagonal shows DEG count.
# Outputs → figures_STAR_s0/supplementary/figS_sensitivity/panels/
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(BASE,
  "figures_STAR_s0/supplementary/figS_sensitivity/panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Load dream results
# ─────────────────────────────────────────────────────────────────────────────
message("Loading dream results...")
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "dream_results_ashr.csv"),
  select = c("gene", "logFC", "padj", "symbol"))
message("  ", nrow(dream), " genes loaded")

padj_vals <- c(1e-3, 0.005, 0.01, 0.025, 0.05, 0.1)
lfc_vals  <- c(0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5)

size_ratio <- function(a, b) {
  min(length(unique(a)), length(unique(b))) /
  max(length(unique(a)), length(unique(b)))
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build lower-triangle size-ratio data.table from a named list of sets
# ─────────────────────────────────────────────────────────────────────────────
make_ratio_dt <- function(sets) {
  n     <- length(sets)
  nms   <- names(sets)
  sizes <- vapply(sets, function(s) length(unique(s)), integer(1))

  rows <- vector("list", n * n)
  k <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      val <- if (i == j) NA_real_ else size_ratio(sets[[i]], sets[[j]])
      rows[[k]] <- data.table(xi = i, xj = j, val = val)
      k <- k + 1L
    }
  }
  dt <- rbindlist(rows)
  dt[xi < xj, val := NA_real_]   # blank upper triangle
  dt[, is_diag := xi == xj]
  dt[, label := fcase(
    is_diag,    format(sizes[xi], big.mark = ","),
    !is.na(val), sprintf("%.2f", val),
    default = ""
  )]
  dt[, x_fac := factor(nms[xi], levels = nms)]
  dt[, y_fac := factor(nms[xj], levels = rev(nms))]
  dt
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: draw the triangle size-ratio heatmap
# ─────────────────────────────────────────────────────────────────────────────
draw_ratio_heatmap <- function(dt, x_lab, y_lab, title, canonical_x = NULL) {
  p <- ggplot(dt, aes(x = x_fac, y = y_fac, fill = val)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label, fontface = ifelse(is_diag, "bold", "plain")),
              size = 2.2, color = "gray10") +
    scale_fill_gradient(low = "#F5F9FB", high = "#C65B2E",
                        limits = c(0, 1), na.value = "white",
                        name = "Size\nratio") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_masld() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.8, "cm"),
          legend.key.width  = unit(0.25, "cm"))

  if (!is.null(canonical_x)) {
    n_lvls <- length(levels(dt$x_fac))
    cx <- which(levels(dt$x_fac) == canonical_x)
    p <- p +
      annotate("rect",
               xmin = cx - 0.5, xmax = cx + 0.5,
               ymin = 0.5, ymax = n_lvls + 0.5,
               color = "#FFB300", fill = NA, linewidth = 0.9)
  }
  p
}

# ═══════════════════════════════════════════════════════════════════════════
# Panel a2: size ratio across padj thresholds at fixed |LFC| = 0.5
# ═══════════════════════════════════════════════════════════════════════════
message("Panel a2: size ratio across padj thresholds (|LFC| = 0.5)...")

padj_labels <- ifelse(padj_vals < 0.001,
                      format(padj_vals, scientific = TRUE),
                      as.character(padj_vals))

sets_padj <- setNames(
  lapply(padj_vals, function(pa) dream[padj < pa & abs(logFC) >= 0.5, symbol]),
  padj_labels
)

dt_padj <- make_ratio_dt(sets_padj)

p_a2 <- draw_ratio_heatmap(
  dt_padj,
  x_lab       = "padj threshold",
  y_lab       = "padj threshold",
  title       = "DEG-set size ratio across padj thresholds (|LFC| \u2265 0.5)",
  canonical_x = "0.05"
)

out_a2 <- file.path(PANEL_DIR, "deg_landscape_a_sizeratio_padj.pdf")
save_fig(p_a2, out_a2, width = fig_half_width, height = fig_half_width)
message("Saved: ", out_a2)

# ═══════════════════════════════════════════════════════════════════════════
# Panel a3: size ratio across |LFC| thresholds at fixed padj = 0.05
# ═══════════════════════════════════════════════════════════════════════════
message("Panel a3: size ratio across |LFC| thresholds (padj = 0.05)...")

lfc_labels <- sprintf("%.1f", lfc_vals)

sets_lfc <- setNames(
  lapply(lfc_vals, function(lf) dream[padj < 0.05 & abs(logFC) >= lf, symbol]),
  lfc_labels
)

dt_lfc <- make_ratio_dt(sets_lfc)

p_a3 <- draw_ratio_heatmap(
  dt_lfc,
  x_lab       = "|log\u2082FC| threshold",
  y_lab       = "|log\u2082FC| threshold",
  title       = "DEG-set size ratio across |LFC| thresholds (padj < 0.05)",
  canonical_x = "0.5"
)

out_a3 <- file.path(PANEL_DIR, "deg_landscape_a_sizeratio_lfc.pdf")
save_fig(p_a3, out_a3, width = fig_half_width, height = fig_half_width)
message("Saved: ", out_a3)

message("\n=== figS_deg_threshold_landscape_pairwise.R complete ===")
