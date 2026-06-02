#!/usr/bin/env Rscript
# figS_progression_landscape.R
# Cross-cohort disease progression landscape (NAS × Fibrosis)
# Reproduces Hoang/Govaere-style 2D heatmaps using dream mega-analysis signatures
#
# 4 panels + summary annotation:
#   A: Sample distribution
#   B: NAS gene signature intensity
#   C: Fibrosis gene signature intensity
#   D: NAS/Fibrosis signature ratio

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

# --- Load data ---
cells <- fread(file.path(SIGS, "progression_landscape_cells.csv"))
sample_dist <- fread(file.path(SIGS, "progression_sample_distribution.csv"))
samples <- fread(file.path(SIGS, "progression_landscape_samples.csv"))

cat("Landscape cells:", nrow(cells), "\n")
cat("Total samples:", sum(sample_dist$N), "\n")

# Signature correlation
cor_test <- cor.test(samples$nas_sig_score, samples$fib_sig_score)
cat("Signature correlation: r =", round(cor_test$estimate, 3),
    ", p =", format(cor_test$p.value, digits = 3), "\n")

# Create full NAS × Fibrosis grid (include empty cells as NA)
full_grid <- CJ(nas_group = 0:7, fib_stage = 0:4)
cells_full <- merge(full_grid, cells, by = c("nas_group", "fib_stage"), all.x = TRUE)
dist_full <- merge(full_grid, sample_dist, by = c("nas_group", "fib_stage"), all.x = TRUE)
dist_full[is.na(N), N := 0]

# Factor labels
cells_full[, nas_label := factor(nas_group)]
cells_full[, fib_label := factor(paste0("F", fib_stage), levels = paste0("F", 4:0))]
dist_full[, nas_label := factor(nas_group)]
dist_full[, fib_label := factor(paste0("F", fib_stage), levels = paste0("F", 4:0))]

# Shared theme for heatmap panels
hm_theme <- theme_masld() +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 7, face = "bold"))

# --- Panel A: Sample Distribution ---
p_a <- ggplot(dist_full, aes(x = nas_label, y = fib_label, fill = N)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = dist_full[N > 0], aes(label = N),
            size = 2.2, color = "gray20") +
  scale_fill_gradient(low = "#FFF9C4", high = "#7B1FA2",
                      name = "n samples", na.value = "gray95") +
  labs(x = "NAS Score", y = "Fibrosis Stage",
       title = "A. Sample Distribution") +
  hm_theme

# --- Panel B: NAS Gene Signature ---
p_b <- ggplot(cells_full[!is.na(nas_sig_mean)],
              aes(x = nas_label, y = fib_label, fill = nas_sig_mean)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "inferno",
                        name = expression("Mean log"[2]*"(CPM+1)"),
                        na.value = "gray95") +
  labs(x = "NAS Score", y = "Fibrosis Stage",
       title = "B. NAS Gene Signature") +
  hm_theme

# --- Panel C: Fibrosis Gene Signature ---
p_c <- ggplot(cells_full[!is.na(fib_sig_mean)],
              aes(x = nas_label, y = fib_label, fill = fib_sig_mean)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colours = c("#FFF0F5", "#FFB6C1", "#FF69B4", "#C71585", "#800020"),
    name = expression("Mean log"[2]*"(CPM+1)"),
    na.value = "gray95") +
  labs(x = "NAS Score", y = "Fibrosis Stage",
       title = "C. Fibrosis Gene Signature") +
  hm_theme

# --- Panel D: NAS/Fibrosis Signature Ratio ---
ratio_range <- cells_full[!is.na(ratio), range(ratio)]
max_dev <- max(abs(ratio_range - 1))

p_d <- ggplot(cells_full[!is.na(ratio)],
              aes(x = nas_label, y = fib_label, fill = ratio)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828",
                        midpoint = 1, name = "Ratio",
                        limits = c(1 - max_dev, 1 + max_dev)) +
  labs(x = "NAS Score", y = "Fibrosis Stage",
       title = "D. NAS/Fibrosis Signature Ratio") +
  hm_theme

# --- Summary text panel ---
summary_text <- paste0(
  "Dataset Summary\n",
  "Total samples: ", sum(sample_dist$N), "\n",
  "NAS range: 0-", max(cells$nas_group), "\n",
  "Fibrosis range: F0-F4\n",
  "Datasets: 5 (GSE130970, GSE135251,\n",
  "  GSE162694, GSE174478, GSE193066)\n\n",
  "Signature Analysis\n",
  "Top 50 NAS genes\n",
  "Top 50 Fibrosis genes\n",
  "Correlation: r = ", round(cor_test$estimate, 3),
  "\n  p = ", format(cor_test$p.value, digits = 3), "\n\n",
  "Key Observations\n",
  "NAS signature increases with\n  NAS score across all F-stages\n",
  "Fibrosis signature increases\n  primarily with F-stage\n",
  "Signatures are positively\n  correlated"
)

p_summary <- ggplot() +
  annotate("text", x = 0, y = 0, label = summary_text,
           hjust = 0, vjust = 0.5, size = 2.2, family = "mono",
           lineheight = 1.2) +
  theme_void() +
  xlim(-0.1, 3) + ylim(-1, 1)

# --- Assembly ---
fig <- ((p_a | p_b) / (p_c | p_d)) | p_summary
fig <- fig +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "Cross-Cohort Disease Progression Landscape (Integrated Mega-Analysis)",
    theme = theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  )

OUT_PDF <- file.path(FIGS02_DIR, "figS02_progression_landscape.pdf")
save_fig_tall(fig, OUT_PDF, width = 10, height = 7, dpi = 300)
cat("Progression landscape saved to:", OUT_PDF, "\n")
