#!/usr/bin/env Rscript
# fig1d_candidates.R — Cumulative discovery (rarefaction) curve for Fig 1 panel (d)
# Shows why multi-cohort integration matters: naive union inflates with noise,
# but dream recovers weak-but-consistent cross-cohort signals invisible to any single study.
# Thresholds match UpSet panel (e): padj < 0.1 & |LFC| > 0.5
# Output: figures/main/fig1_atlas_overview/fig1d_rarefaction.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# FIG1_DIR is defined in load_figure_data.R (figures/main/fig1_atlas_overview)

PANEL_DIR <- file.path(FIG1_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Data prep — thresholds matching UpSet panel (e)
# ============================================================================
message("Loading data...")
dream <- load_dream_results()
per_study <- load_per_study_de()

COMPARE_STUDIES <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
STUDY_SHORT <- c(GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
                 GSE162694 = "Bril",   GSE213621 = "Chen")

# padj-only for rarefaction (LFC cutoff asymmetrically penalizes
# integrated analysis where averaged LFCs shrink toward zero)
padj_thr <- 0.1

# Build per-study DEG lists (padj-only)
deg_lists <- lapply(COMPARE_STUDIES, function(ds) {
  dt <- per_study[dataset == ds & padj < padj_thr]
  unique(sub("\\..*", "", dt$gene))
})
names(deg_lists) <- STUDY_SHORT[COMPARE_STUDIES]

# Dream DEGs (padj-only)
dream_degs <- unique(sub("\\..*", "", dream[dream_padj < padj_thr, gene]))

message(sprintf("Dream DEGs (padj<%.1f): %s", padj_thr, comma(length(dream_degs))))
for (nm in names(deg_lists)) {
  message(sprintf("  %s: %s DEGs", nm, comma(length(deg_lists[[nm]]))))
}

# ============================================================================
# Rarefaction curve with permutations
# ============================================================================
message("Computing rarefaction curves...")
set.seed(42)
n_perms <- 120
study_names_5 <- names(deg_lists)

perm_curves <- rbindlist(lapply(seq_len(n_perms), function(p) {
  ord <- sample(study_names_5)
  cumul <- integer(5)
  running_union <- character(0)
  for (k in seq_along(ord)) {
    running_union <- union(running_union, deg_lists[[ord[k]]])
    cumul[k] <- length(running_union)
  }
  data.table(perm = p, n_cohorts = 1:5, n_degs = cumul)
}))

# Summary stats per n_cohorts
ribbon_dt <- perm_curves[, .(mean_degs = mean(n_degs),
                              min_degs  = min(n_degs),
                              max_degs  = max(n_degs)),
                          by = n_cohorts]

n_dream <- length(dream_degs)
union_mean_5 <- round(ribbon_dt[n_cohorts == 5, mean_degs])

# Individual study DEG counts for jittered points at x=1
study_counts <- data.table(
  study = names(deg_lists),
  n_degs = sapply(deg_lists, length)
)
study_counts[, x_jitter := 1 + seq(-0.15, 0.15, length.out = nrow(study_counts))]

# ============================================================================
# Plot
# ============================================================================
message("Generating rarefaction plot...")

p_rarefaction <- ggplot() +
  # Permutation ribbon (min/max)
  geom_ribbon(data = ribbon_dt,
              aes(x = n_cohorts, ymin = min_degs, ymax = max_degs),
              fill = "gray70", alpha = 0.15) +
  # Mean curve
  geom_line(data = ribbon_dt,
            aes(x = n_cohorts, y = mean_degs),
            linewidth = 0.9, color = masld_colors$up) +
  geom_point(data = ribbon_dt,
             aes(x = n_cohorts, y = mean_degs),
             size = 1.6, color = masld_colors$up, shape = 16) +
  # Individual study points at x=1 (jittered)
  geom_point(data = study_counts,
             aes(x = x_jitter, y = n_degs),
             size = 1.2, color = masld_colors$up, alpha = 0.7, shape = 16) +
  geom_text(data = study_counts,
            aes(x = x_jitter, y = n_degs, label = study),
            size = 1.6, color = "gray40", vjust = -0.8, hjust = 0.5) +
  # Dream horizontal line
  geom_hline(yintercept = n_dream, linetype = "dashed", linewidth = 0.5,
             color = masld_colors$down) +
  # Dream label
  annotate("text", x = 1.1, y = n_dream,
           label = paste0("Integrated: ", comma(n_dream)),
           size = 2.2, color = masld_colors$down, hjust = 0, vjust = -0.6,
           fontface = "bold") +
  # Union mean label at x=5
  annotate("text", x = 5, y = union_mean_5,
           label = paste0(comma(union_mean_5), "\n(naive union)"),
           size = 2, color = "gray40", hjust = -0.1, vjust = 0.5,
           lineheight = 0.85) +
  # Permutation count annotation
  annotate("text", x = 1, y = max(ribbon_dt$max_degs) * 1.02,
           label = paste0("n = ", n_perms, " permutations"),
           size = 1.7, color = "gray50", hjust = 0) +
  scale_x_continuous(breaks = 1:5, labels = 1:5) +
  scale_y_continuous(labels = comma) +
  labs(x = "Number of cohorts combined",
       y = paste0("Unique DEGs (padj < ", padj_thr, ")"),
       title = "Cumulative discovery curve vs. integrated analysis") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold"))

save_fig(p_rarefaction, file.path(PANEL_DIR, "fig1d_rarefaction.pdf"),
         width = fig_half_width, height = 3)
message("  Saved fig1d_rarefaction.pdf")

# ============================================================================
# Done
# ============================================================================
fp <- file.path(PANEL_DIR, "fig1d_rarefaction.pdf")
if (file.exists(fp)) {
  message(sprintf("Output: %s (%s)", fp,
                  utils:::format.object_size(file.size(fp), "auto")))
} else {
  message("ERROR: Output file not created")
}
