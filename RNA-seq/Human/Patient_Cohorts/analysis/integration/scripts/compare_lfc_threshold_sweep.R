#!/usr/bin/env Rscript
# compare_lfc_threshold_sweep.R
# ---------------------------------------------------------------------------
# Build a side-by-side comparison of Pillar A stability metrics across 4 LFC
# cutoffs (0, 0.1, 0.3, 0.5). Emits a CSV summary table and a 4-panel PDF.
# Run AFTER aggregate_pillar_A.R has been run with each LFC_THR.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
AUDIT <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Threshold catalog (extended 2026-05-13 to include 0.75, 1, 2)
sweep <- data.table(
  thr   = c(0,        0.1,    0.3,    0.5,                     0.75,     1,        2),
  tag   = c("nolfc",  "lfc1", "lfc3", "lfc5",                  "lfc075", "lfc100", "lfc200"),
  label = c("|LFC|>0\n(padj<0.05 only)",
            "|LFC|>0.1", "|LFC|>0.3",
            "|LFC|>0.5\n(Tier 1)",
            "|LFC|>0.75", "|LFC|>1.0", "|LFC|>2.0")
)

# --- Read each per-gene file and compute headline metrics ---
collect <- function(tag, thr) {
  f <- file.path(AUDIT, paste0("pillar_A_stability_", tag, ".csv"))
  if (!file.exists(f)) {
    warning(sprintf("Missing %s — skipping", f)); return(NULL)
  }
  pA <- fread(f)
  deg <- pA[is_canonical_DEG == TRUE]
  data.table(
    threshold     = thr,
    tag           = tag,
    n_DEG         = nrow(deg),
    mean_pi_hat   = round(mean(deg$cpss_pi_hat, na.rm = TRUE), 3),
    pi_hat_ge_07  = sum(deg$cpss_pi_hat >= 0.7, na.rm = TRUE),
    pct_pi_ge_07  = round(100 * mean(deg$cpss_pi_hat >= 0.7, na.rm = TRUE), 1),
    mean_boot     = round(mean(deg$bootstrap_freq, na.rm = TRUE), 3),
    boot_ge_09    = sum(deg$bootstrap_freq >= 0.9, na.rm = TRUE),
    pct_boot_ge_09 = round(100 * mean(deg$bootstrap_freq >= 0.9, na.rm = TRUE), 1),
    mean_kfold10  = round(mean(deg$kfold10_recur, na.rm = TRUE), 2),
    mean_kfold50  = round(mean(deg$kfold50_recur, na.rm = TRUE), 2),
    mean_pi_nonDEG = round(mean(pA[is_canonical_DEG == FALSE]$cpss_pi_hat, na.rm = TRUE), 3))
}

cmp <- rbindlist(lapply(seq_len(nrow(sweep)),
                         function(i) collect(sweep$tag[i], sweep$thr[i])))
cmp[, label := sweep$label[match(threshold, sweep$thr)]]
fwrite(cmp, file.path(AUDIT, "pillar_A_lfc_threshold_sweep.csv"))
cat("=== Comparison table ===\n"); print(cmp)

# --- Per-iter ρ data, one source per threshold ---
iter_long <- rbindlist(lapply(sweep$tag, function(tag) {
  f <- file.path(AUDIT, paste0("pillar_A_iter_summary_", tag, ".csv"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f); dt[, threshold := sweep$thr[match(tag, sweep$tag)]]
  dt[, label := sweep$label[match(threshold, sweep$thr)]]
  dt
}))

# ============================================================================
# Figure: 4-panel comparison
# ============================================================================
.short_levels <- c("0", "0.1", "0.3", "0.5", "0.75", "1", "2")
.short_map    <- setNames(.short_levels, as.character(sweep$thr))
cmp[, label_short := factor(.short_map[as.character(threshold)], levels = .short_levels)]

p1 <- ggplot(cmp, aes(x = label_short, y = n_DEG, fill = label_short)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(n_DEG)), vjust = -0.3, size = 3) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "A — DEG count by |logFC| cutoff",
       x = "|logFC| threshold", y = "DEGs (padj<0.05 + |LFC|>x)") +
  theme_pub() + theme(legend.position = "none")

p2 <- ggplot(cmp, aes(x = label_short, y = mean_pi_hat, fill = label_short)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_pi_hat)), vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.05))) +
  labs(title = expression("B — mean CPSS "*hat(pi)*" of DEGs"),
       x = "|logFC| threshold", y = expression(hat(pi))) +
  theme_pub() + theme(legend.position = "none")

p3 <- ggplot(cmp, aes(x = label_short, y = pct_pi_ge_07, fill = label_short)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_pi_ge_07)), vjust = -0.3, size = 3) +
  geom_hline(yintercept = 70, linetype = "dashed", colour = "grey30") +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.05))) +
  labs(title = expression("C — % DEGs with "*hat(pi)*" "*"" >= ""*" 0.7"),
       x = "|logFC| threshold", y = "% DEGs",
       caption = "Shah-Samworth π_thr = 0.7") +
  theme_pub() + theme(legend.position = "none")

# Recovery summary: combine 10-fold and 50-fold
p4 <- ggplot(cmp, aes(x = label_short)) +
  geom_col(aes(y = mean_kfold10 / 10 * 100, fill = "K=10 (10% held out)"),
           position = position_dodge(width = 0.8), width = 0.35) +
  geom_col(aes(y = mean_kfold50 / 50 * 100, fill = "K=50 (2% held out)"),
           position = position_nudge(x = 0.4), width = 0.35) +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("K=10 (10% held out)" = "#4B0082",
                                "K=50 (2% held out)"  = "#E14B9D")) +
  labs(title = "D — K-fold recurrence (% of folds DEG was re-selected)",
       x = "|logFC| threshold", y = "% recurrence per fold", fill = "") +
  theme_pub() + theme(legend.position = "bottom")

# Per-iter ρ violin
if (nrow(iter_long) > 0) {
  iter_long[, label_short := factor(.short_map[as.character(threshold)],
                                     levels = .short_levels)]
  iter_long[, src_lbl := factor(source, levels = c("cpss", "bootstrap", "kfold_K10", "kfold_K50"),
                                 labels = c("CPSS", "Bootstrap", "K-fold-10", "K-fold-50"))]
  p5 <- ggplot(iter_long[!is.na(rho)], aes(x = label_short, y = rho, fill = label_short)) +
    geom_violin(scale = "width", alpha = 0.6, colour = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
    facet_wrap(~ src_lbl, nrow = 1) +
    geom_hline(yintercept = 0.85, linetype = "dashed", colour = "grey30") +
    labs(title = "E — Per-iter Spearman ρ by source × threshold",
         x = "|logFC| threshold", y = "ρ vs full dream") +
    theme_pub() + theme(legend.position = "none")
} else {
  p5 <- ggplot() + theme_void()
}

combined <- (p1 | p2 | p3) / p4 / p5 + plot_layout(heights = c(1, 1, 1.3))
out_pdf <- file.path(OUT_DIR, "figS_robustness_lfc_sweep.pdf")
ggsave(out_pdf, combined, width = 14, height = 11, device = cairo_pdf)
cat(sprintf("Saved: %s\n", out_pdf))
fwrite(iter_long[!is.na(rho), .(threshold, source, rho, jaccard, recovery_pct, n_sel)],
       file.path(AUDIT, "pillar_A_lfc_sweep_iter_long.csv"))
cat(sprintf("Saved: %s\n", file.path(AUDIT, "pillar_A_lfc_sweep_iter_long.csv")))
