#!/usr/bin/env Rscript
# M05_mouse_volcano_plots.R
# ---------------------------------------------------------------------------
# Publication-quality volcano plots for per-diet DE and dream mega-analysis.
# Uses Sanjana Lab mouse color scheme (Cyan/Green family).
# Input:  results/per_diet/*_de_results.csv, results/meta_analysis/lvqw_pooled_results.csv
# Output: results/volcanos_per_diet.pdf, results/volcano_dream_mega.pdf
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration"
RDIR <- file.path(BASE, "results")

# --- Publication theme ---
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Mouse color scheme (Cyan/Green family); NS = control gray
volcano_colors <- c("NS" = "#9E9E9E", "Up" = "#00BCD4", "Down" = "#1B5E20")
LFSR_THRESH <- 0.05  # ashr local false sign rate gate (shrunk effects)
LFC_THRESH  <- 0.3   # 2026-06-27: 0.5 -> 0.3 (track new human canonical; config mouse_lfc=0.3)
# Gating now uses ashr-shrunk effect sizes (shrunk_logFC) and lfsr, not raw logFC/padj.
# Per-diet *_de_results.csv carry shrunk_logFC + lfsr (ashr applied upstream);
# pooled lvqw_pooled_results.csv carries the same.

# ===== Per-Diet Volcano Plots =====
cat("Loading per-diet DE results...\n")
diet_files <- list.files(file.path(RDIR, "per_diet"),
                          pattern = "_de_results\\.csv$", full.names = TRUE)
per_diet <- rbindlist(lapply(diet_files, function(f) {
  dt <- fread(f)
  # Extract diet name from filename
  dt[, diet := gsub("_de_results\\.csv$", "", basename(f))]
  dt
}))

cat("Total genes across diets:", nrow(per_diet), "\n")
cat("Diets:", paste(unique(per_diet$diet), collapse = ", "), "\n")

# Classify genes (ashr-shrunk effect + lfsr)
per_diet[, category := fcase(
  lfsr < LFSR_THRESH & shrunk_logFC > LFC_THRESH, "Up",
  lfsr < LFSR_THRESH & shrunk_logFC < -LFC_THRESH, "Down",
  default = "NS"
)]
per_diet[, category := factor(category, levels = c("NS", "Up", "Down"))]

# Cap -log10(lfsr) for visualization
per_diet[, neg_log10_lfsr := pmin(-log10(lfsr), 50)]

# --- Multi-panel volcano: all diets ---
cat("\nGenerating per-diet volcano plots...\n")
cairo_pdf(file.path(RDIR, "volcanos_per_diet.pdf"), width = 14, height = 10)

# Label top 10 genes per diet (by lfsr among significant)
top_labels <- per_diet[category != "NS",
                        .SD[order(lfsr)][1:min(10, .N)],
                        by = diet]

# Order diets for display (auto-discovered from loaded data)
per_diet[, diet := factor(diet, levels = sort(unique(diet)))]

p <- ggplot(per_diet, aes(x = shrunk_logFC, y = neg_log10_lfsr, color = category)) +
  geom_point(data = per_diet[category == "NS"], size = 0.3, alpha = 0.3) +
  geom_point(data = per_diet[category != "NS"], size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(LFSR_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = volcano_colors) +
  facet_wrap(~ diet, scales = "free", ncol = 3) +
  labs(
    title = "Per-Diet Volcano Plots — Mouse MASLD Models (Disease vs Control)",
    x = expression(ashr ~ "shrunk" ~ log[2] ~ "fold change"),
    y = expression(-log[10] ~ "lfsr"),
    color = ""
  ) +
  theme_pub
print(p)
dev.off()
# Caption (no in-plot subtitle, house style)
message("CAPTION volcanos_per_diet.pdf: per-diet ashr-shrunk volcanoes; ",
        "gate lfsr<", LFSR_THRESH, " & |shrunk_logFC|>", LFC_THRESH,
        "; dashed guides at +/-", LFC_THRESH, ".")
cat("Saved: volcanos_per_diet.pdf\n")

# Print per-diet DEG counts
for (d in levels(per_diet$diet)) {
  sub <- per_diet[diet == d]
  n_up <- sum(sub$category == "Up")
  n_down <- sum(sub$category == "Down")
  cat(sprintf("  %s: Up=%d Down=%d\n", d, n_up, n_down))
}

# ===== Pooled (cohort-adjusted) Volcano =====
cat("\nLoading pooled limma-voom-qw results...\n")
dream <- fread(file.path(RDIR, "meta_analysis/lvqw_pooled_results.csv"))

dream[, category := fcase(
  lfsr < LFSR_THRESH & shrunk_logFC > LFC_THRESH, "Up",
  lfsr < LFSR_THRESH & shrunk_logFC < -LFC_THRESH, "Down",
  default = "NS"
)]
dream[, category := factor(category, levels = c("NS", "Up", "Down"))]
dream[, neg_log10_lfsr := pmin(-log10(lfsr), 100)]

# Top 20 labeled genes
top_dream <- dream[category != "NS"][order(lfsr)][1:min(20, sum(dream$category != "NS"))]

n_up <- sum(dream$category == "Up")
n_down <- sum(dream$category == "Down")
# Derive counts from data (not hardcoded)
n_diets <- length(unique(per_diet$diet))
# Read QC-passing sample count from metadata if available, else from de_summary
qc_file <- file.path(BASE, "qc/sample_qc_report.csv")
if (file.exists(qc_file)) {
  qc_meta   <- fread(qc_file)
  n_samples <- sum(qc_meta$pass_qc == TRUE, na.rm = TRUE)
} else {
  meta_file <- file.path(RDIR, "meta_matched.rds")
  if (file.exists(meta_file)) {
    n_samples <- nrow(readRDS(meta_file))
  } else {
    n_samples <- NA_integer_
  }
}

cairo_pdf(file.path(RDIR, "volcano_dream_mega.pdf"), width = 9, height = 7)
p2 <- ggplot(dream, aes(x = shrunk_logFC, y = neg_log10_lfsr, color = category)) +
  geom_point(data = dream[category == "NS"], size = 0.5, alpha = 0.2) +
  geom_point(data = dream[category != "NS"], size = 1, alpha = 0.6) +
  geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(LFSR_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_text_repel(
    data = top_dream,
    aes(label = gene),
    size = 2.5, max.overlaps = 20,
    segment.color = "grey60", segment.size = 0.3
  ) +
  scale_color_manual(values = volcano_colors) +
  labs(
    title = "Pooled Volcano — Mouse MASLD (Disease vs Control)",
    x = expression(ashr ~ "shrunk" ~ log[2] ~ "fold change"),
    y = expression(-log[10] ~ "lfsr"),
    color = ""
  ) +
  theme_pub
print(p2)
dev.off()
# Caption (no in-plot subtitle, house style)
message(sprintf("CAPTION volcano_dream_mega.pdf: pooled limma-voom-qw, %s samples, %d diets; ",
                ifelse(is.na(n_samples), "?", as.character(n_samples)), n_diets),
        sprintf("gate lfsr<%.2f & |shrunk_logFC|>%.1f; Up=%d Down=%d; guides at +/-%.1f.",
                LFSR_THRESH, LFC_THRESH, n_up, n_down, LFC_THRESH))

cat(sprintf("Pooled volcano saved. Up: %d, Down: %d\n", n_up, n_down))
cat("Saved: volcano_dream_mega.pdf\n")
