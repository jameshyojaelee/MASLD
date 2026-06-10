#!/usr/bin/env Rscript
# figS_stage_transitions.R
# Stage transition analysis: consecutive pairwise DEG counts and mean |LFC|
# Reproduces Govaere-style transition analysis using dream mega-analysis
#
# 4 panels: NAS mean|LFC|, Fibrosis mean|LFC|, NAS DEG counts, Fibrosis DEG counts

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

SIGS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")
PADJ_THRESH <- 0.05  # Tier-2 progression convention (padj<0.05, no LFC). Was 0.1 (Govaere-matched) before 2026-06-09 harmonization.
LFC_THRESH  <- 0     # No LFC threshold for transition counts (Tier-2 convention)

# --- Load consecutive dream results ---
nas_consec <- fread(file.path(SIGS, "nas_consecutive_dream.csv"))
fib_consec <- fread(file.path(SIGS, "fibrosis_consecutive_dream.csv"))

cat("NAS consecutive:", nrow(nas_consec), "rows,", length(unique(nas_consec$contrast)), "contrasts\n")
cat("Fibrosis consecutive:", nrow(fib_consec), "rows,", length(unique(fib_consec$contrast)), "contrasts\n")

# --- Panel A: NAS Expression Magnitude (Mean |LFC|) ---
nas_mag <- nas_consec[, .(mean_abs_lfc = mean(abs(logFC), na.rm = TRUE)), by = contrast]
# Create transition labels (e.g., "0→1")
nas_mag[, transition := gsub("NAS(\\d+)_vs_NAS(\\d+)", "\\2→\\1", contrast)]
# Order by transition
nas_mag[, order_idx := as.integer(gsub(".*→(\\d+)", "\\1", transition))]
nas_mag <- nas_mag[order(order_idx)]
nas_mag[, transition := factor(transition, levels = transition)]

p_a <- ggplot(nas_mag, aes(x = transition, y = mean_abs_lfc)) +
  geom_col(fill = "#FF8A65", width = 0.7) +
  geom_text(aes(label = round(mean_abs_lfc, 3)), vjust = -0.5, size = 2.2) +
  labs(x = NULL, y = "Mean |log2FC|",
       title = "NAS: Expression Magnitude (Mean |LFC|)") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6, angle = 0))

# --- Panel B: Fibrosis Expression Magnitude (Mean |LFC|) ---
fib_mag <- fib_consec[, .(mean_abs_lfc = mean(abs(logFC), na.rm = TRUE)), by = contrast]
fib_mag[, transition := gsub("F(\\d+)_vs_F(\\d+)", "F\\2→F\\1", contrast)]
fib_mag[, order_idx := as.integer(gsub(".*→F(\\d+)", "\\1", transition))]
fib_mag <- fib_mag[order(order_idx)]
fib_mag[, transition := factor(transition, levels = transition)]

p_b <- ggplot(fib_mag, aes(x = transition, y = mean_abs_lfc)) +
  geom_col(fill = "#F48FB1", width = 0.7) +
  geom_text(aes(label = round(mean_abs_lfc, 3)), vjust = -0.5, size = 2.2) +
  labs(x = NULL, y = "Mean |log2FC|",
       title = "Fibrosis: Expression Magnitude (Mean |LFC|)") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6, angle = 0))

# --- Panel C: NAS Significant Genes (mirrored bar) ---
nas_sig <- nas_consec[padj < PADJ_THRESH]
nas_deg_up <- nas_sig[logFC > 0, .(n = .N, direction = "Upregulated"), by = contrast]
nas_deg_down <- nas_sig[logFC < 0, .(n = -.N, direction = "Downregulated"), by = contrast]
nas_deg <- rbindlist(list(nas_deg_up, nas_deg_down))

# Add transition labels
nas_deg[, transition := gsub("NAS(\\d+)_vs_NAS(\\d+)", "\\2→\\1", contrast)]
nas_deg[, order_idx := as.integer(gsub(".*→(\\d+)", "\\1", transition))]
nas_deg <- nas_deg[order(order_idx)]

# Ensure all transitions present
all_nas_trans <- nas_mag$transition
for (tr in levels(all_nas_trans)) {
  for (dir in c("Upregulated", "Downregulated")) {
    if (!any(nas_deg$transition == tr & nas_deg$direction == dir)) {
      nas_deg <- rbindlist(list(nas_deg, data.table(
        contrast = NA, n = 0, direction = dir, transition = tr, order_idx = NA)))
    }
  }
}
nas_deg[, transition := factor(transition, levels = levels(all_nas_trans))]

# Labels for bar counts
nas_deg[, label := as.character(abs(n))]
nas_deg[n == 0, label := ""]

p_c <- ggplot(nas_deg, aes(x = transition, y = n, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = label,
                vjust = ifelse(n >= 0, -0.3, 1.3)),
            size = 2, color = "gray20") +
  scale_fill_manual(values = c(Upregulated = "#EC407A", Downregulated = "#FFAB91"),
                    name = NULL) +
  labs(x = NULL, y = paste0("Significant Genes (FDR < ", PADJ_THRESH, ")"),
       title = "NAS: Significant Genes") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5))

# --- Panel D: Fibrosis Significant Genes (mirrored bar) ---
fib_sig <- fib_consec[padj < PADJ_THRESH]
fib_deg_up <- fib_sig[logFC > 0, .(n = .N, direction = "Upregulated"), by = contrast]
fib_deg_down <- fib_sig[logFC < 0, .(n = -.N, direction = "Downregulated"), by = contrast]
fib_deg <- rbindlist(list(fib_deg_up, fib_deg_down))

fib_deg[, transition := gsub("F(\\d+)_vs_F(\\d+)", "F\\2→F\\1", contrast)]
fib_deg[, order_idx := as.integer(gsub(".*→F(\\d+)", "\\1", transition))]
fib_deg <- fib_deg[order(order_idx)]

all_fib_trans <- fib_mag$transition
for (tr in levels(all_fib_trans)) {
  for (dir in c("Upregulated", "Downregulated")) {
    if (!any(fib_deg$transition == tr & fib_deg$direction == dir)) {
      fib_deg <- rbindlist(list(fib_deg, data.table(
        contrast = NA, n = 0, direction = dir, transition = tr, order_idx = NA)))
    }
  }
}
fib_deg[, transition := factor(transition, levels = levels(all_fib_trans))]

fib_deg[, label := as.character(abs(n))]
fib_deg[n == 0, label := ""]

p_d <- ggplot(fib_deg, aes(x = transition, y = n, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = label,
                vjust = ifelse(n >= 0, -0.3, 1.3)),
            size = 2, color = "gray20") +
  scale_fill_manual(values = c(Upregulated = "#EC407A", Downregulated = "#FFAB91"),
                    name = NULL) +
  labs(x = NULL, y = paste0("Significant Genes (FDR < ", PADJ_THRESH, ")"),
       title = "Fibrosis: Significant Genes") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5))

# --- Assembly ---
fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "Cross-Cohort Transition Analysis (Integrated Mega-Analysis)",
    theme = theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  )

OUT_PDF <- file.path(FIGS02_DIR, "figS_stage_transitions.pdf")
save_fig_tall(fig, OUT_PDF, width = fig_full_width, height = 7, dpi = 300)
cat("Stage transitions saved to:", OUT_PDF, "\n")
