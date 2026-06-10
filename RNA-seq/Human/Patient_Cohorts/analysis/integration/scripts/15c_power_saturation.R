#!/usr/bin/env Rscript
# 15c_power_saturation.R
# ---------------------------------------------------------------------------
# Power saturation analysis for the MASLD dream mega-analysis.
#
# Demonstrates diminishing returns of sample size by running dream on
# stratified random subsamples at N = {200, 400, 600, 800, 1000, 1200},
# 10 iterations each. Compares each subsample to the full-model (N=1,444)
# dream results.
#
# Metrics: DEG count, Spearman rho, direction concordance, Jaccard,
#          positive control recovery.
#
# Output:
#   - RNA-seq/results/audit_sensitivity/power_saturation/power_saturation_results.csv
#   - figures/supplementary/sensitivity/figS_power_saturation.pdf
# ---------------------------------------------------------------------------

t0 <- proc.time()

# ---------------------------------------------------------------------------
# Library loading + reformulas injection (MUST match script 05)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE    <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
INT     <- file.path(BASE, "analysis/integration")
RDIR    <- file.path(INT, "results/integration")

OUT_DIR <- file.path(PROJECT, "RNA-seq/results/audit_sensitivity/power_saturation")
FIG_DIR <- file.path(PROJECT, "figures/supplementary/figS_methods_validation/sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)

# Publication theme
source(file.path(PROJECT, "scripts/figures/publication_theme.R"))

# ---------------------------------------------------------------------------
# Parallel setup
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
TARGET_NS   <- c(200L, 400L, 600L, 800L, 1000L, 1200L)
N_ITER      <- 10L
PADJ_THRESH <- 0.1

# ============================================================================
# 1. Load full-model dream results (reference)
# ============================================================================
cat("\n===== Loading full-model dream results =====\n")

full_dream <- fread(file.path(RDIR, "dream_results.csv"))
cat(sprintf("Full model: %d genes, %d DEGs (padj < %.1f)\n",
            nrow(full_dream),
            sum(full_dream$padj < PADJ_THRESH, na.rm = TRUE),
            PADJ_THRESH))

full_deg_genes <- full_dream[padj < PADJ_THRESH, gene]
full_n_degs    <- length(full_deg_genes)

# Build named LFC vector for full model
full_lfc <- setNames(full_dream$logFC, full_dream$gene)

# ============================================================================
# 2. Load positive controls
# ============================================================================
cat("\n===== Loading positive controls =====\n")

pc <- fread(file.path(PROJECT, "results/library/positive_control.csv"))
pc_genes <- pc$`Gene symbol`
# Only keep positive controls that are in the dream results
pc_genes <- intersect(pc_genes, full_dream$gene)
cat(sprintf("Positive controls in dream results: %d / %d\n",
            length(pc_genes), nrow(pc)))

# Which positive controls are DEGs in the full model?
pc_full_recovery <- sum(pc_genes %in% full_deg_genes) / length(pc_genes)
cat(sprintf("Full-model positive control recovery: %.1f%%\n", pc_full_recovery * 100))

# ============================================================================
# 3. Load merged DGE and prepare metadata
# ============================================================================
cat("\n===== Loading merged DGE =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Exclude datasets with irrecoverable confounds (same as script 05)
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "
")
cat("Samples available:", ncol(dge_mega), "\n")

# Load metadata with inferred sex
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# Build info data.frame
info <- data.frame(
  sample_id    = colnames(dge_mega),
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

# Drop samples with NA sex
na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) {
  cat("Dropping", sum(na_sex), "samples with NA sex\n")
  dge_mega <- dge_mega[, !na_sex]
  info     <- info[!na_sex, , drop = FALSE]
}

N_FULL <- ncol(dge_mega)
cat(sprintf("Total samples for analysis: %d\n", N_FULL))

cat("\nDataset distribution:\n")
print(table(info$dataset, info$group_binary))

# ============================================================================
# 4. Stratified subsampling + dream
# ============================================================================
cat("\n===== Starting power saturation analysis =====\n")
cat(sprintf("Target Ns: %s\n", paste(TARGET_NS, collapse = ", ")))
cat(sprintf("Iterations per N: %d\n", N_ITER))
cat(sprintf("Total dream runs: %d\n", length(TARGET_NS) * N_ITER))

# Dream formula (same as script 05 base model)
form <- ~ group_binary + inferred_sex + (1|dataset)

# Pre-compute dataset-level sampling frames
datasets <- levels(info$dataset)
ds_indices <- lapply(datasets, function(d) which(info$dataset == d))
names(ds_indices) <- datasets
ds_sizes <- sapply(ds_indices, length)

cat("\nDataset sizes:\n")
print(ds_sizes)

# Results collector
results_list <- list()
run_counter <- 0L
total_runs <- length(TARGET_NS) * N_ITER

for (target_n in TARGET_NS) {
  cat(sprintf("\n===== N = %d =====\n", target_n))

  if (target_n > N_FULL) {
    cat("Skipping: target N exceeds available samples\n")
    next
  }

  for (iter in seq_len(N_ITER)) {
    run_counter <- run_counter + 1L
    seed <- 42000L + target_n + iter
    set.seed(seed)

    cat(sprintf("  [%d/%d] N=%d, iter=%d (seed=%d)... ",
                run_counter, total_runs, target_n, iter, seed))

    # --- Stratified sampling: proportional to dataset size ---
    # Allocate samples per dataset proportionally, ensuring at least 2 per
    # dataset (need >=2 for random effect estimation) and at least 1 per
    # group_binary level per dataset
    frac <- target_n / N_FULL
    draw_per_ds <- setNames(pmax(round(ds_sizes * frac), 2L), datasets)

    # Adjust total to match target_n (trim from largest datasets)
    total_allocated <- sum(draw_per_ds)
    if (total_allocated != target_n) {
      diff <- total_allocated - target_n
      # Sort datasets by size descending, trim/add from largest
      ds_order <- names(sort(ds_sizes, decreasing = TRUE))
      i <- 1L
      while (diff != 0) {
        d <- ds_order[((i - 1L) %% length(ds_order)) + 1L]
        if (diff > 0 && draw_per_ds[d] > 2L) {
          draw_per_ds[d] <- draw_per_ds[d] - 1L
          diff <- diff - 1L
        } else if (diff < 0 && draw_per_ds[d] < ds_sizes[d]) {
          draw_per_ds[d] <- draw_per_ds[d] + 1L
          diff <- diff + 1L
        }
        i <- i + 1L
        if (i > length(ds_order) * abs(diff) + 100L) break  # safety valve
      }
    }

    # Within each dataset, stratified sample preserving group_binary ratio
    sampled_idx <- integer(0)
    skip_this <- FALSE

    for (d in datasets) {
      d_idx <- ds_indices[[d]]
      d_info <- info[d_idx, ]
      n_draw <- draw_per_ds[d]

      # Cap at actual dataset size
      n_draw <- min(n_draw, length(d_idx))

      ctrl_idx <- d_idx[d_info$group_binary == "Control"]
      dis_idx  <- d_idx[d_info$group_binary == "Disease"]

      if (length(ctrl_idx) == 0 || length(dis_idx) == 0) {
        # Dataset has only one group — draw all from that group
        sampled_idx <- c(sampled_idx, sample(d_idx, n_draw))
        next
      }

      # Proportional split
      ctrl_frac <- length(ctrl_idx) / length(d_idx)
      n_ctrl <- max(round(n_draw * ctrl_frac), 1L)
      n_dis  <- n_draw - n_ctrl
      n_ctrl <- min(n_ctrl, length(ctrl_idx))
      n_dis  <- min(n_dis, length(dis_idx))

      sampled_idx <- c(sampled_idx,
                       sample(ctrl_idx, n_ctrl),
                       sample(dis_idx, n_dis))
    }

    # Subset DGE
    dge_sub <- dge_mega[, sampled_idx]
    info_sub <- info[sampled_idx, , drop = FALSE]
    rownames(info_sub) <- colnames(dge_sub)

    # Drop factor levels with no observations
    info_sub$dataset <- droplevels(info_sub$dataset)

    actual_n <- ncol(dge_sub)

    # Check we have both levels of group_binary
    if (length(unique(info_sub$group_binary)) < 2) {
      cat("SKIP (missing group level)\n")
      next
    }

    # Check we have >1 dataset level for random effect
    if (nlevels(info_sub$dataset) < 2) {
      cat("SKIP (single dataset)\n")
      next
    }

    # --- Run dream ---
    tryCatch({
      v_sub <- suppressWarnings(
        voomWithDreamWeights(dge_sub, form, info_sub, BPPARAM = param)
      )
      fit_sub <- suppressWarnings(
        dream(v_sub, form, info_sub, BPPARAM = param)
      )

      res_sub <- topTable(fit_sub, coef = "group_binaryDisease",
                          number = Inf, sort.by = "none")
      res_sub$gene <- rownames(res_sub)
      res_dt <- as.data.table(res_sub)

      # --- Compute metrics ---
      sub_deg_genes <- res_dt[adj.P.Val < PADJ_THRESH, gene]
      n_degs <- length(sub_deg_genes)

      # Spearman rho (over all shared genes)
      shared_genes <- intersect(res_dt$gene, full_dream$gene)
      sub_lfc <- setNames(res_dt$logFC, res_dt$gene)
      rho <- cor(full_lfc[shared_genes], sub_lfc[shared_genes],
                 method = "spearman", use = "complete.obs")

      # Direction concordance (among genes DEG in BOTH)
      shared_degs <- intersect(sub_deg_genes, full_deg_genes)
      if (length(shared_degs) > 0) {
        same_dir <- sum(sign(full_lfc[shared_degs]) == sign(sub_lfc[shared_degs]))
        dir_concordance <- same_dir / length(shared_degs)
      } else {
        dir_concordance <- NA_real_
      }

      # Jaccard
      union_degs <- length(union(sub_deg_genes, full_deg_genes))
      jaccard <- if (union_degs > 0) {
        length(shared_degs) / union_degs
      } else {
        0
      }

      # Positive control recovery
      pc_recovered <- sum(pc_genes %in% sub_deg_genes) / length(pc_genes)

      results_list[[length(results_list) + 1]] <- data.table(
        target_n      = target_n,
        actual_n      = actual_n,
        iteration     = iter,
        seed          = seed,
        n_degs        = n_degs,
        n_shared_degs = length(shared_degs),
        spearman_rho  = rho,
        dir_concordance = dir_concordance,
        jaccard       = jaccard,
        pc_recovery   = pc_recovered,
        n_genes_tested = nrow(res_dt)
      )

      cat(sprintf("N=%d, DEGs=%d, rho=%.3f, Jaccard=%.3f, PC=%.1f%%\n",
                  actual_n, n_degs, rho, jaccard, pc_recovered * 100))

    }, error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    })
  }
}

# ============================================================================
# 5. Compile results
# ============================================================================
cat("\n===== Compiling results =====\n")

results <- rbindlist(results_list)

# Add the full model as a reference row
results_with_full <- rbind(
  results,
  data.table(
    target_n      = N_FULL,
    actual_n      = N_FULL,
    iteration     = 0L,
    seed          = NA_integer_,
    n_degs        = full_n_degs,
    n_shared_degs = full_n_degs,
    spearman_rho  = 1.0,
    dir_concordance = 1.0,
    jaccard       = 1.0,
    pc_recovery   = pc_full_recovery,
    n_genes_tested = nrow(full_dream)
  )
)

# Save results
out_csv <- file.path(OUT_DIR, "power_saturation_results.csv")
fwrite(results_with_full, out_csv)
cat("Saved results:", out_csv, "\n")

# Summary statistics per target_n
summary_dt <- results[, .(
  mean_degs        = mean(n_degs),
  sd_degs          = sd(n_degs),
  ci_lo_degs       = mean(n_degs) - 1.96 * sd(n_degs),
  ci_hi_degs       = mean(n_degs) + 1.96 * sd(n_degs),
  mean_rho         = mean(spearman_rho, na.rm = TRUE),
  sd_rho           = sd(spearman_rho, na.rm = TRUE),
  mean_jaccard     = mean(jaccard),
  sd_jaccard       = sd(jaccard),
  mean_pc          = mean(pc_recovery),
  sd_pc            = sd(pc_recovery),
  mean_dir_conc    = mean(dir_concordance, na.rm = TRUE),
  sd_dir_conc      = sd(dir_concordance, na.rm = TRUE)
), by = target_n]

cat("\n===== Summary =====\n")
print(summary_dt)

# ============================================================================
# 6. Generate figure
# ============================================================================
cat("\n===== Generating figure =====\n")

# Prepare plot data
plot_dt <- copy(results)
plot_dt[, target_n_f := factor(target_n)]

# --- Panel (a): DEG count vs N ---
p_a <- ggplot(plot_dt, aes(x = target_n, y = n_degs)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$up, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$up, size = 1.5) +
  stat_summary(fun.data = function(x) {
    data.frame(ymin = mean(x) - 1.96 * sd(x),
               ymax = mean(x) + 1.96 * sd(x))
  }, geom = "ribbon", alpha = 0.2, fill = masld_colors$up) +
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

# --- Panel (b): Spearman rho vs N ---
p_b <- ggplot(plot_dt, aes(x = target_n, y = spearman_rho)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$down, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$down, size = 1.5) +
  stat_summary(fun.data = function(x) {
    data.frame(ymin = mean(x) - 1.96 * sd(x),
               ymax = mean(x) + 1.96 * sd(x))
  }, geom = "ribbon", alpha = 0.2, fill = masld_colors$down) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_x_continuous(breaks = c(TARGET_NS, N_FULL),
                     labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  coord_cartesian(ylim = c(NA, 1.0)) +
  labs(x = "Sample size (N)", y = "Spearman rho vs full model",
       title = "Effect size correlation") +
  theme_masld()

# --- Panel (c): Jaccard vs N ---
p_c <- ggplot(plot_dt, aes(x = target_n, y = jaccard)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$conserved, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$conserved, size = 1.5) +
  stat_summary(fun.data = function(x) {
    data.frame(ymin = mean(x) - 1.96 * sd(x),
               ymax = mean(x) + 1.96 * sd(x))
  }, geom = "ribbon", alpha = 0.2, fill = masld_colors$conserved) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_x_continuous(breaks = c(TARGET_NS, N_FULL),
                     labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  coord_cartesian(ylim = c(0, 1.0)) +
  labs(x = "Sample size (N)", y = "Jaccard similarity vs full model",
       title = "DEG set stability") +
  theme_masld()

# --- Panel (d): Positive control recovery vs N ---
p_d <- ggplot(plot_dt, aes(x = target_n, y = pc_recovery)) +
  stat_summary(fun = mean, geom = "line", color = masld_colors$human_enriched, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", color = masld_colors$human_enriched, size = 1.5) +
  stat_summary(fun.data = function(x) {
    data.frame(ymin = mean(x) - 1.96 * sd(x),
               ymax = mean(x) + 1.96 * sd(x))
  }, geom = "ribbon", alpha = 0.2, fill = masld_colors$human_enriched) +
  geom_hline(yintercept = pc_full_recovery, linetype = "dashed", color = "gray40",
             linewidth = 0.4) +
  annotate("text", x = max(TARGET_NS), y = pc_full_recovery,
           label = sprintf("Full model: %.0f%%", pc_full_recovery * 100),
           hjust = 1, vjust = -0.5, size = 2, color = "gray40") +
  scale_x_continuous(breaks = c(TARGET_NS, N_FULL),
                     labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Sample size (N)", y = "Positive control recovery",
       title = "Known target recovery") +
  theme_masld()

# --- Compose 2x2 ---
fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = 8, face = "bold"))
  )

# Save
fig_path <- file.path(FIG_DIR, "figS_power_saturation.pdf")
save_fig(fig, fig_path, width = fig_full_width, height = 5)
cat("Saved figure:", fig_path, "\n")

# ============================================================================
# 7. Print saturation assessment
# ============================================================================
cat("\n===== Saturation Assessment =====\n")

# Marginal gain from N=1000 to N=1200
if (all(c(1000, 1200) %in% summary_dt$target_n)) {
  gain_1000_1200 <- summary_dt[target_n == 1200, mean_degs] -
                    summary_dt[target_n == 1000, mean_degs]
  gain_800_1000  <- summary_dt[target_n == 1000, mean_degs] -
                    summary_dt[target_n == 800, mean_degs]
  gain_600_800   <- summary_dt[target_n == 800, mean_degs] -
                    summary_dt[target_n == 600, mean_degs]

  cat(sprintf("Marginal DEG gain (600->800):   +%.0f DEGs (+%.1f%%)\n",
              gain_600_800,
              gain_600_800 / summary_dt[target_n == 600, mean_degs] * 100))
  cat(sprintf("Marginal DEG gain (800->1000):  +%.0f DEGs (+%.1f%%)\n",
              gain_800_1000,
              gain_800_1000 / summary_dt[target_n == 800, mean_degs] * 100))
  cat(sprintf("Marginal DEG gain (1000->1200): +%.0f DEGs (+%.1f%%)\n",
              gain_1000_1200,
              gain_1000_1200 / summary_dt[target_n == 1000, mean_degs] * 100))
}

# Rho at N=1200
if (1200 %in% summary_dt$target_n) {
  cat(sprintf("\nAt N=1200: rho=%.3f, Jaccard=%.3f, PC recovery=%.1f%%\n",
              summary_dt[target_n == 1200, mean_rho],
              summary_dt[target_n == 1200, mean_jaccard],
              summary_dt[target_n == 1200, mean_pc] * 100))
}

cat(sprintf("\nFull model (N=%d): %s DEGs, PC recovery=%.1f%%\n",
            N_FULL,
            formatC(full_n_degs, format = "d", big.mark = ","),
            pc_full_recovery * 100))

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\nTotal runtime: %.1f hours (%.0f minutes)\n",
            elapsed / 3600, elapsed / 60))
