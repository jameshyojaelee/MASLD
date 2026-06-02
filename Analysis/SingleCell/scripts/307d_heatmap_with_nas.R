#!/usr/bin/env Rscript
#' 307d: Bulk-SC concordance heatmap with BOTH fibrosis AND NAS transitions.
#'
#' Extends Panel D from 307c by including NAS transition signatures
#' (NAS01→NAS24, NAS24→NAS5, NAS5→NAS68) alongside fibrosis stages.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PT_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

theme_src <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(theme_src)) source(theme_src)

CT_LABELS <- c(Hepatocytes = "Hepatocytes", Macrophages = "Macrophages",
               Fibroblasts = "Fibroblasts", Endothelial_cells = "Endothelial",
               Cholangiocytes = "Cholangiocytes")

cat("=== 307d: Bulk-SC Heatmap with NAS Transitions ===\n")

# =========================================================================
# Load concordance data
# =========================================================================
conc <- fread(file.path(PT_DIR, "bulk_sc_concordance.csv"))
cat(sprintf("Loaded %d rows from bulk_sc_concordance.csv\n", nrow(conc)))

# Filter to palantir; fall back to dpt if absent
conc_sel <- conc[pseudotime_method == "palantir"]
if (nrow(conc_sel) == 0) {
  cat("No palantir rows; falling back to dpt\n")
  conc_sel <- conc[pseudotime_method == "dpt"]
}
cat(sprintf("Using %d rows (method: %s)\n", nrow(conc_sel),
            conc_sel[1, pseudotime_method]))

# =========================================================================
# Clean signature names
# =========================================================================
conc_sel[, sig_clean := gsub("transition_", "", signature)]

# NAS group labels: NAS01 -> "NAS 0-1", NAS24 -> "NAS 2-4", etc.
conc_sel[, sig_clean := gsub("NAS01", "NAS 0-1", sig_clean)]
conc_sel[, sig_clean := gsub("NAS24", "NAS 2-4", sig_clean)]
conc_sel[, sig_clean := gsub("NAS5",  "NAS 5",   sig_clean)]
conc_sel[, sig_clean := gsub("NAS68", "NAS 6-8", sig_clean)]

# Arrows and direction symbols
conc_sel[, sig_clean := gsub("_to_", " \u2192 ", sig_clean)]
conc_sel[, sig_clean := gsub("_DOWN", " \u2193", sig_clean)]
conc_sel[, sig_clean := gsub("_UP",   " \u2191", sig_clean)]

# =========================================================================
# Define signature order: fibrosis first, then NAS, each UP then DOWN
# =========================================================================
fib_order <- c(
  "F0 \u2192 F1 \u2191", "F0 \u2192 F1 \u2193",
  "F1 \u2192 F2 \u2191", "F1 \u2192 F2 \u2193",
  "F2 \u2192 F3 \u2191", "F2 \u2192 F3 \u2193",
  "F3 \u2192 F4 \u2191", "F3 \u2192 F4 \u2193"
)
nas_order <- c(
  "NAS 0-1 \u2192 NAS 2-4 \u2191", "NAS 0-1 \u2192 NAS 2-4 \u2193",
  "NAS 2-4 \u2192 NAS 5 \u2191",   "NAS 2-4 \u2192 NAS 5 \u2193",
  "NAS 5 \u2192 NAS 6-8 \u2191",   "NAS 5 \u2192 NAS 6-8 \u2193"
)
sig_order <- c(fib_order, nas_order)

# Keep only signatures in our defined order
conc_sel <- conc_sel[sig_clean %in% sig_order]
cat(sprintf("Retained %d rows after signature filtering\n", nrow(conc_sel)))

# Factor: reversed so first item appears at top of heatmap
conc_sel[, sig_clean := factor(sig_clean, levels = rev(sig_order))]

# =========================================================================
# Cell type labels
# =========================================================================
conc_sel[, ct_label := CT_LABELS[cell_type]]
ct_level_order <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                     "Endothelial", "Cholangiocytes")
conc_sel[, ct_label := factor(ct_label, levels = ct_level_order)]

# =========================================================================
# Build heatmap
# =========================================================================
# Horizontal separator between fibrosis and NAS sections.
# In the reversed factor, NAS signatures occupy bottom positions (1-6),
# fibrosis occupies top positions (7-14). The boundary sits at 6.5.
n_nas <- length(nas_order)
hline_y <- n_nas + 0.5

p <- ggplot(conc_sel, aes(x = ct_label, y = sig_clean, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho),
                color = ifelse(abs(spearman_rho) > 0.3, "high", "low")),
            size = 2.0, show.legend = FALSE) +
  scale_color_manual(values = c(high = "white", low = "gray30")) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0, limits = c(-0.65, 0.65),
                       name = "Spearman \u03c1") +
  geom_hline(yintercept = hline_y, linewidth = 0.6, color = "gray40") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.6, "cm"),
        plot.title = element_text(size = 8, face = "bold")) +
  labs(x = NULL, y = NULL,
       title = "Bulk transition signatures vs\nsc pseudotime")

# =========================================================================
# Save
# =========================================================================
out <- file.path(FIG_DIR, "figS_bulk_sc_heatmap_full.pdf")
save_fig(p, out, width = 5, height = 5)
cat("Saved:", out, "\n")

cat("=== 307d COMPLETE ===\n")
