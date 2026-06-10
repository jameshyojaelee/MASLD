#!/usr/bin/env Rscript
# 15c_aggregate.R
# ---------------------------------------------------------------------------
# Aggregate per-iteration outputs from 15c_one_iter.R + regenerate the
# original power_saturation_results.csv and figS_power_saturation.pdf.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

PROJECT  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE     <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
INT      <- file.path(BASE, "analysis/integration")
RDIR     <- file.path(INT, "results/integration")
OUT_DIR  <- file.path(PROJECT, "RNA-seq/results/audit_sensitivity/power_saturation")
ITER_DIR <- file.path(OUT_DIR, "iter")
FIG_DIR  <- file.path(PROJECT, "figures/supplementary/figS_methods_validation/sensitivity")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT, "scripts/figures/publication_theme.R"))

TARGET_NS <- c(200L, 400L, 600L, 800L, 1000L, 1200L)
N_ITER    <- 10L
PADJ_THRESH <- 0.1

# --- Load partial CSVs ---
files <- list.files(ITER_DIR, pattern = "^iter_N\\d+_i\\d+\\.csv$", full.names = TRUE)
cat(sprintf("Found %d iter CSVs in %s\n", length(files), ITER_DIR))
if (!length(files)) stop("No per-iter CSVs to aggregate.")

results <- rbindlist(lapply(files, fread))
cat(sprintf("Aggregated %d rows\n", nrow(results)))

# --- Compute full-model reference (same as monolithic script) ---
full_dream <- fread(file.path(RDIR, "dream_results.csv"))
full_deg_genes <- full_dream[padj < PADJ_THRESH, gene]
full_n_degs <- length(full_deg_genes)

pc <- fread(file.path(PROJECT, "results/library/positive_control.csv"))
pc_genes <- intersect(pc$`Gene symbol`, full_dream$gene)
pc_full_recovery <- sum(pc_genes %in% full_deg_genes) / length(pc_genes)

# Compute N_FULL the same way as 15c_one_iter
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(PROJECT, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]
na_sex <- is.na(matched_sex)
N_FULL <- sum(!na_sex)
cat(sprintf("N_FULL = %d\n", N_FULL))

# Append full-model reference row
results_with_full <- rbind(
  results,
  data.table(
    target_n        = N_FULL,
    actual_n        = N_FULL,
    iteration       = 0L,
    seed            = NA_integer_,
    n_degs          = full_n_degs,
    n_shared_degs   = full_n_degs,
    spearman_rho    = 1.0,
    dir_concordance = 1.0,
    jaccard         = 1.0,
    pc_recovery     = pc_full_recovery,
    n_genes_tested  = nrow(full_dream)
  ),
  fill = TRUE
)

out_csv <- file.path(OUT_DIR, "power_saturation_results.csv")
fwrite(results_with_full, out_csv)
cat("Saved:", out_csv, "\n")

# --- Summary ---
summary_dt <- results[, .(
  mean_degs        = mean(n_degs),
  sd_degs          = sd(n_degs),
  mean_rho         = mean(spearman_rho, na.rm = TRUE),
  sd_rho           = sd(spearman_rho, na.rm = TRUE),
  mean_jaccard     = mean(jaccard),
  sd_jaccard       = sd(jaccard),
  mean_pc          = mean(pc_recovery),
  sd_pc            = sd(pc_recovery),
  mean_dir_conc    = mean(dir_concordance, na.rm = TRUE),
  sd_dir_conc      = sd(dir_concordance, na.rm = TRUE),
  n_iters          = .N
), by = target_n][order(target_n)]
print(summary_dt)

# --- Figure (4 panels, same as monolithic script) ---
plot_dt <- copy(results)

p_a <- ggplot(plot_dt, aes(x = target_n, y = n_degs)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$up, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$up, size = 1.5) +
  stat_summary(fun.data = function(x) data.frame(ymin = mean(x) - 1.96 * sd(x),
                                                  ymax = mean(x) + 1.96 * sd(x)),
               geom = "ribbon", alpha = 0.2, fill = masld_colors$up) +
  geom_hline(yintercept = full_n_degs, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  annotate("text", x = max(TARGET_NS), y = full_n_degs,
           label = sprintf("Full model (N=%d): %s DEGs", N_FULL,
                           formatC(full_n_degs, format = "d", big.mark = ",")),
           hjust = 1, vjust = -0.5, size = 2, color = "gray40") +
  scale_x_continuous(breaks = c(TARGET_NS, N_FULL),
                     labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Sample size (N)", y = "DEGs (padj < 0.1)", title = "DEG discovery") +
  theme_masld()

p_b <- ggplot(plot_dt, aes(x = target_n, y = spearman_rho)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$down, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$down, size = 1.5) +
  stat_summary(fun.data = function(x) data.frame(ymin = mean(x) - 1.96 * sd(x),
                                                  ymax = mean(x) + 1.96 * sd(x)),
               geom = "ribbon", alpha = 0.2, fill = masld_colors$down) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_x_continuous(breaks = c(TARGET_NS, N_FULL),
                     labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  coord_cartesian(ylim = c(NA, 1.0)) +
  labs(x = "Sample size (N)", y = "Spearman rho vs full model",
       title = "Effect size correlation") +
  theme_masld()

p_c <- ggplot(plot_dt, aes(x = target_n, y = jaccard)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$conserved, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$conserved, size = 1.5) +
  stat_summary(fun.data = function(x) data.frame(ymin = mean(x) - 1.96 * sd(x),
                                                  ymax = mean(x) + 1.96 * sd(x)),
               geom = "ribbon", alpha = 0.2, fill = masld_colors$conserved) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_x_continuous(breaks = c(TARGET_NS, N_FULL),
                     labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  coord_cartesian(ylim = c(0, 1.0)) +
  labs(x = "Sample size (N)", y = "Jaccard similarity vs full model",
       title = "DEG set stability") +
  theme_masld()

p_d <- ggplot(plot_dt, aes(x = target_n, y = pc_recovery)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$human_enriched, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$human_enriched, size = 1.5) +
  stat_summary(fun.data = function(x) data.frame(ymin = mean(x) - 1.96 * sd(x),
                                                  ymax = mean(x) + 1.96 * sd(x)),
               geom = "ribbon", alpha = 0.2, fill = masld_colors$human_enriched) +
  geom_hline(yintercept = pc_full_recovery, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  annotate("text", x = max(TARGET_NS), y = pc_full_recovery,
           label = sprintf("Full model: %.0f%%", pc_full_recovery * 100),
           hjust = 1, vjust = -0.5, size = 2, color = "gray40") +
  scale_x_continuous(breaks = c(TARGET_NS, N_FULL),
                     labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Sample size (N)", y = "Positive control recovery",
       title = "Known target recovery") +
  theme_masld()

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 8, face = "bold")))

fig_path <- file.path(FIG_DIR, "figS_power_saturation.pdf")
save_fig(fig, fig_path, width = fig_full_width, height = 5)
cat("Saved figure:", fig_path, "\n")
cat("Done.\n")
