#!/usr/bin/env Rscript
# =============================================================================
# Fig 1 candidate panels (f) and (g): "Why Integration Matters"
# Generates 3 standalone PDF plots for visual comparison:
#   1. Integration-rescued genes histogram + known gene annotation strip
#   2. Cohort accumulation / rarefaction curve (union-based)
#   3. Per-gene forest plot strip (effect size concordance)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ─── Shared: author name mapping ─────────────────────────────────────────────
STUDY_NAMES <- c(
  GSE213621  = "Chen",
  GSE135251  = "Govaere",
  GSE130970  = "Hoang",
  GSE162694  = "Bril",
  GSE174478  = "Kawamura",
  GSE193066  = "Hoshida",
  GSE240729  = "Verschuren",
  GSE126848  = "Suppli",
  GSE167523  = "Kozumi"
)
# PRJNA512027 (Gerhard 2018) excluded from cohort presentation: L0/S0
# library-prep batch perfectly confounded with diagnosis.

# ─── Load data ───────────────────────────────────────────────────────────────
cat("Loading data...\n")
dream   <- load_dream_results()
per_study <- load_per_study_de()

dream[, gene_clean := sub("\\..*", "", gene)]
per_study[, gene_clean := sub("\\..*", "", gene)]

# Studies for comparison — PRJNA512027 already excluded from cohort
# presentation upstream; GSE167523 (NAFL-only) also excluded here.
COMPARE_STUDIES <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694",
                     "GSE174478", "GSE193066", "GSE213621", "GSE240729")
ps <- per_study[dataset %in% COMPARE_STUDIES]
N_STUDIES <- length(COMPARE_STUDIES)
cat(sprintf("Using %d cohorts for per-study comparison\n", N_STUDIES))


# =============================================================================
# IDEA 1: Integration-Rescued Genes Histogram
# =============================================================================
cat("\n── Idea 1: Integration-rescued genes histogram ──\n")

dream_degs <- dream[dream_padj < 0.1]
cat(sprintf("Dream DEGs (padj<0.1): %s\n", comma(nrow(dream_degs))))

gene_n_tested <- ps[, .(n_tested = uniqueN(dataset)), by = gene_clean]
genes_ok <- gene_n_tested[n_tested >= 3, gene_clean]
dream_degs_f <- dream_degs[gene_clean %in% genes_ok]
cat(sprintf("After >=3 study filter: %s DEGs\n", comma(nrow(dream_degs_f))))

# Per-study detection at padj < 0.1 (same threshold as dream)
ps_sig <- ps[gene_clean %in% dream_degs_f$gene_clean & padj < 0.1,
             .(n_sig = uniqueN(dataset)), by = gene_clean]
dream_degs_f <- merge(dream_degs_f, ps_sig, by = "gene_clean", all.x = TRUE)
dream_degs_f[is.na(n_sig), n_sig := 0L]

# Summary stats
n_zero <- sum(dream_degs_f$n_sig == 0)
pct_zero <- round(100 * n_zero / nrow(dream_degs_f), 1)
n_low  <- sum(dream_degs_f$n_sig <= 1)
pct_low <- round(100 * n_low / nrow(dream_degs_f), 1)
n_half <- sum(dream_degs_f$n_sig <= N_STUDIES / 2)
pct_half <- round(100 * n_half / nrow(dream_degs_f), 1)
median_sig <- median(dream_degs_f$n_sig)
cat(sprintf("  0-study: %s (%s%%); 0-1: %s (%s%%); <=4: %s (%s%%); median=%d\n",
            comma(n_zero), pct_zero, comma(n_low), pct_low,
            comma(n_half), pct_half, median_sig))

# ── Known MASLD genes (>= 10 required) ──────────────────────────────────────
known_masld <- c(
  # Fibrosis / ECM
  "COL1A1", "COL3A1", "ACTA2", "FAP", "THY1", "TIMP1",
  # Inflammation / immune
  "TREM2", "SPP1", "CD68",
  # Lipid metabolism
  "CIDEC", "PLIN2", "FASN", "SCD",
  # Bile acid / hepatocyte function
  "CYP7A1", "CYP1A2", "SLC27A5", "ACSM3",
  # Other established MASLD genes
  "LPL", "GDF15", "FABP1"
)

# Get bin assignments for known genes (must be dream DEGs)
gene_bins <- dream_degs_f[symbol %in% known_masld][!duplicated(symbol),
                           .(symbol, n_sig)]
gene_bins <- gene_bins[order(n_sig, symbol)]
cat(sprintf("  Known MASLD genes in dream DEGs: %d / %d\n",
            nrow(gene_bins), length(known_masld)))
for (i in seq_len(nrow(gene_bins)))
  cat(sprintf("    %s: %d/%d cohorts\n",
              gene_bins$symbol[i], gene_bins$n_sig[i], N_STUDIES))

# ── Build histogram ──────────────────────────────────────────────────────────
grad_cols <- colorRampPalette(c("#0D47A1", "#7B1FA2", "#C9265E"))(N_STUDIES + 1)
bar_dt <- dream_degs_f[, .N, by = n_sig][order(n_sig)]
max_bar <- max(bar_dt$N)

# Gene annotation strip: position known genes above their bars
gene_bins[, rank := seq_len(.N), by = n_sig]
y_base  <- max_bar * 1.03
y_step  <- max_bar * 0.045
gene_bins[, y_pos := y_base + (rank - 1) * y_step]

p1 <- ggplot(dream_degs_f, aes(x = factor(n_sig, levels = 0:N_STUDIES))) +
  geom_bar(aes(fill = factor(n_sig, levels = 0:N_STUDIES)), show.legend = FALSE) +
  scale_fill_manual(values = setNames(grad_cols, as.character(0:N_STUDIES)),
                    drop = FALSE) +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(drop = FALSE) +
  # Gene annotation: diamonds + repelled labels
  geom_point(data = gene_bins,
             aes(x = factor(n_sig, levels = 0:N_STUDIES), y = y_pos),
             shape = 18, size = 1.2, color = "#880E4F") +
  geom_text_repel(data = gene_bins,
                  aes(x = factor(n_sig, levels = 0:N_STUDIES),
                      y = y_pos, label = symbol),
                  size = 1.8, fontface = "italic", color = "#880E4F",
                  segment.size = 0.15, segment.color = "#880E4F60",
                  box.padding = 0.15, point.padding = 0.1,
                  max.overlaps = 30, direction = "both",
                  min.segment.length = 0,
                  nudge_y = max_bar * 0.02, seed = 42) +
  # Median dashed line
  geom_vline(xintercept = median_sig + 1, linetype = "dotted",
             color = "gray40", linewidth = 0.35) +
  annotate("text", x = median_sig + 1.4, y = max_bar * 0.55,
           label = sprintf("median = %d", median_sig),
           hjust = 0, size = 2, color = "gray40") +
  coord_cartesian(ylim = c(0, max(gene_bins$y_pos) + max_bar * 0.12),
                  clip = "off") +
  labs(x = "# individual cohorts detecting gene (padj < 0.1)",
       y = "Number of integrated DEGs (padj < 0.1)") +
  theme_masld()

save_fig(p1, file.path(FIG1_DIR, "fig1_idea1_rescued_genes.pdf"),
         width = fig_half_width, height = 4.2)
cat("  Saved fig1_idea1_rescued_genes.pdf\n")


# =============================================================================
# IDEA 2: Cohort Accumulation Curve (simplified)
# =============================================================================
cat("\n── Idea 2: Cohort accumulation curve ──\n")

# Per-study DEG sets at padj < 0.1 (matching dream)
study_deg_sets <- lapply(COMPARE_STUDIES, function(s) {
  ps[dataset == s & padj < 0.1, unique(gene_clean)]
})
names(study_deg_sets) <- COMPARE_STUDIES

cat("Per-study DEG counts (padj<0.1):\n")
for (s in COMPARE_STUDIES) {
  cat(sprintf("  %s (%s): %s\n", STUDY_NAMES[s], s,
              comma(length(study_deg_sets[[s]]))))
}

# Dream DEGs used as reference
dream_gene_set <- dream_degs[gene_clean %in% genes_ok, unique(gene_clean)]
dream_total <- length(dream_gene_set)

# Accumulation: for each random ordering, track how many dream DEGs are in
# the growing per-study union
set.seed(42)
N_PERMS <- 100

accum_dt <- rbindlist(lapply(seq_len(N_PERMS), function(perm) {
  ord <- sample(COMPARE_STUDIES)
  running <- character()
  rbindlist(lapply(seq_along(ord), function(k) {
    running <<- union(running, study_deg_sets[[ord[k]]])
    data.table(perm = perm, k = k,
               n_recovered = sum(dream_gene_set %in% running))
  }))
}))

accum_summ <- accum_dt[, .(mean_n = mean(n_recovered),
                            lo = quantile(n_recovered, 0.025),
                            hi = quantile(n_recovered, 0.975)), by = k]

cat(sprintf("  Dream total: %s DEGs\n", comma(dream_total)))
cat(sprintf("  Recovered at k=1: %s (%.0f%%); k=%d: %s (%.0f%%)\n",
            comma(round(accum_summ[k == 1, mean_n])),
            100 * accum_summ[k == 1, mean_n] / dream_total,
            N_STUDIES,
            comma(round(accum_summ[k == N_STUDIES, mean_n])),
            100 * accum_summ[k == N_STUDIES, mean_n] / dream_total))

p2 <- ggplot() +
  # Dream total as horizontal reference
  geom_hline(yintercept = dream_total, linetype = "dashed",
             color = "#880E4F", linewidth = 0.5) +
  annotate("text", x = 0.6, y = dream_total,
           label = sprintf("Integrated analysis\n%s DEGs", comma(dream_total)),
           hjust = 0, vjust = -0.3, size = 2.2, fontface = "bold",
           color = "#880E4F", lineheight = 0.85) +
  # Permutation traces (faint)
  geom_line(data = accum_dt, aes(x = k, y = n_recovered, group = perm),
            alpha = 0.06, linewidth = 0.2, color = "gray50") +
  # 95% ribbon
  geom_ribbon(data = accum_summ, aes(x = k, ymin = lo, ymax = hi),
              alpha = 0.15, fill = masld_colors$up) +
  # Mean line
  geom_line(data = accum_summ, aes(x = k, y = mean_n),
            linewidth = 1, color = masld_colors$up) +
  geom_point(data = accum_summ, aes(x = k, y = mean_n),
             size = 2, color = masld_colors$up, shape = 16) +
  # Endpoint labels
  annotate("text", x = N_STUDIES + 0.15,
           y = accum_summ[k == N_STUDIES, mean_n],
           label = sprintf("%s\n(%.0f%%)",
                           comma(round(accum_summ[k == N_STUDIES, mean_n])),
                           100 * accum_summ[k == N_STUDIES, mean_n] / dream_total),
           hjust = 0, vjust = 0.5, size = 2, color = masld_colors$up,
           lineheight = 0.85) +
  annotate("text", x = 1 - 0.15,
           y = accum_summ[k == 1, mean_n],
           label = sprintf("%s\n(%.0f%%)",
                           comma(round(accum_summ[k == 1, mean_n])),
                           100 * accum_summ[k == 1, mean_n] / dream_total),
           hjust = 1, vjust = 0.5, size = 2, color = "gray50",
           lineheight = 0.85) +
  scale_x_continuous(breaks = 1:N_STUDIES,
                     expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "Number of cohorts combined",
       y = "Integrated DEGs overlapping\nper-study DEG union") +
  theme_masld()

save_fig(p2, file.path(FIG1_DIR, "fig1_idea2_rarefaction.pdf"),
         width = fig_half_width, height = 3.5)
cat("  Saved fig1_idea2_rarefaction.pdf\n")


# =============================================================================
# IDEA 3: Per-Gene Forest Plot Strip (expanded, author names)
# =============================================================================
cat("\n── Idea 3: Per-gene forest plot ──\n")

# Expanded candidate list (~20 genes across MASLD biology)
candidates <- c(
  # Fibrosis / ECM (upregulated)
  "COL1A1", "COL3A1", "COL4A1", "ACTA2", "FAP", "THY1", "TIMP1", "VIM",
  # Inflammation / immune (upregulated)
  "TREM2", "SPP1", "CD68", "CCL2",
  # Lipid metabolism
  "CIDEC", "PLIN2", "FASN", "SCD",
  # Bile acid / hepatocyte function (downregulated)
  "CYP7A1", "CYP1A2", "SLC27A5", "ACSM3", "HMGCS2",
  # Other
  "LPL", "GDF15", "FABP1"
)

# Keep genes present in dream results with a valid t-stat
forest_pool <- dream[symbol %in% candidates & !is.na(t) & abs(t) > 0]
forest_pool <- forest_pool[order(dream_padj)][!duplicated(symbol)]

# If fewer than 20, supplement with top dream DEGs
if (nrow(forest_pool) < 20) {
  top_up <- dream[dream_padj < 1e-6 & dream_logFC > 0.5 &
                    !symbol %in% forest_pool$symbol][order(-abs(t))][1:5]
  top_dn <- dream[dream_padj < 1e-6 & dream_logFC < -0.5 &
                    !symbol %in% forest_pool$symbol][order(-abs(t))][1:5]
  extras <- rbind(top_up, top_dn, fill = TRUE)
  extras <- extras[!is.na(symbol)][!duplicated(symbol)]
  forest_pool <- rbind(forest_pool, extras, fill = TRUE)
}
forest_pool <- forest_pool[seq_len(min(22, nrow(forest_pool)))]
forest_genes <- forest_pool$symbol
cat(sprintf("  Selected %d genes: %s\n", length(forest_genes),
            paste(forest_genes, collapse = ", ")))

# Dream estimates + CI
fd <- forest_pool[, .(symbol, gene, dream_logFC, dream_padj, t,
                       se = abs(dream_logFC / t))]
fd[, ci_lo := dream_logFC - 1.96 * se]
fd[, ci_hi := dream_logFC + 1.96 * se]
fd[, sig_label := fifelse(dream_padj < 0.001, "***",
                   fifelse(dream_padj < 0.01, "**",
                   fifelse(dream_padj < 0.1, "*", "ns")))]

# Add I² from meta-analysis (addresses reviewer feedback)
meta <- tryCatch(load_meta_results(), error = function(e) NULL)
if (!is.null(meta)) {
  meta[, gene_clean := sub("\\..*", "", gene)]
  fd[, gene_clean := sub("\\..*", "", gene)]
  fd <- merge(fd, meta[, .(gene_clean, meta_I2)],
              by = "gene_clean", all.x = TRUE)
  fd[, i2_label := fifelse(!is.na(meta_I2),
                            sprintf("I\u00B2=%.0f%%", meta_I2), "")]
  cat("  I² annotations added for", sum(!is.na(fd$meta_I2)), "/", nrow(fd), "genes\n")
} else {
  fd[, i2_label := ""]
  fd[, meta_I2 := NA_real_]
  cat("  WARNING: meta-analysis results not found; skipping I² annotations\n")
}

# Per-study data (include t-stat for inverse-variance weighting)
fps <- ps[symbol %in% forest_genes, .(symbol, logFC, padj, t, dataset)]
fps <- fps[order(padj)][!duplicated(paste(symbol, dataset))]

# Compute inverse-variance weight: 1/SE² where SE = |logFC/t|
fps[, se := fifelse(abs(t) > 0, abs(logFC / t), NA_real_)]
fps[, inv_var := fifelse(!is.na(se) & se > 0, 1 / se^2, NA_real_)]
# Cap extreme weights at 99th percentile for visual clarity
inv_var_cap <- quantile(fps$inv_var, 0.99, na.rm = TRUE)
fps[, inv_var_capped := pmin(inv_var, inv_var_cap, na.rm = TRUE)]

# Map dataset to author name
fps[, author := STUDY_NAMES[dataset]]

# Sort genes: upregulated at top, downregulated at bottom
gene_order <- fd[order(-dream_logFC), symbol]
fd[, symbol := factor(symbol, levels = gene_order)]
fps[, symbol := factor(symbol, levels = gene_order)]

# Author-level colors (ordered by sample size, largest = darkest)
author_order <- c("Chen", "Govaere", "Hoang", "Bril",
                  "Kawamura", "Hoshida", "Verschuren", "Suppli")
author_cols <- setNames(
  c("#0D47A1", "#1565C0", "#42A5F5", "#7B1FA2",
    "#AD1457", "#E91E63", "#F48FB1", "#C2185B"),
  author_order
)
fps[, author := factor(author, levels = author_order)]

# X limits for significance labels
x_range <- range(c(fps$logFC, fd$ci_lo, fd$ci_hi), na.rm = TRUE)
x_pad <- diff(x_range) * 0.06
x_sig_pos <- x_range[2] + x_pad

p3 <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60",
             linewidth = 0.3) +
  # Shaded up/down background
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#FCE4EC", alpha = 0.25) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#E3F2FD", alpha = 0.25) +
  # Per-study dots sized by inverse-variance weight (1/SE²)
  geom_point(data = fps,
             aes(x = logFC, y = symbol, color = author, size = inv_var_capped),
             alpha = 0.6,
             position = position_jitter(height = 0.2, width = 0, seed = 42), shape = 16) +
  scale_size_continuous(range = c(0.5, 2.5), guide = "none") +
  # Dream CI
  geom_errorbar(data = fd,
                aes(xmin = ci_lo, xmax = ci_hi, y = symbol),
                width = 0.3, linewidth = 0.5, color = "#880E4F",
                orientation = "y") +
  # Dream estimate (diamond)
  geom_point(data = fd, aes(x = dream_logFC, y = symbol),
             shape = 18, size = 2.5, color = "#880E4F") +
  # Significance stars + I² annotation
  geom_text(data = fd, aes(x = x_sig_pos, y = symbol, label = sig_label),
            hjust = 0, size = 1.8, color = "#880E4F") +
  geom_text(data = fd[nchar(i2_label) > 0],
            aes(x = x_sig_pos + x_pad * 1.5, y = symbol, label = i2_label),
            hjust = 0, size = 1.5, color = "gray40") +
  scale_color_manual(values = author_cols, name = "Cohort",
                     drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(x = expression(log[2]~"fold-change (Disease vs Control)"),
       y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(face = "italic", size = 5.5),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 5.5),
        plot.margin = margin(3, 42, 3, 3))

save_fig(p3, file.path(FIG1_DIR, "fig1_idea3_forest_plot.pdf"),
         width = fig_half_width, height = 5.5)
cat("  Saved fig1_idea3_forest_plot.pdf\n")


cat("\n=== All 3 panels generated in figures/ ===\n")
cat("  1. fig1_idea1_rescued_genes.pdf\n")
cat("  2. fig1_idea2_rarefaction.pdf\n")
cat("  3. fig1_idea3_forest_plot.pdf\n")
