#!/usr/bin/env Rscript
# diagnose_coloc_snp_counts.R
# Diagnostic: SNP count distribution across COLOC windows
# Question: what fraction of hits survive at different n_snps thresholds (50/100/200)?
# Run: Rscript diagnose_coloc_snp_counts.R
# Output: printed tables + figures/supplementary/coloc_snp_diagnostic.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

RESULTS_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/susie_coloc"
FIGURE_OUT  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/supplementary/coloc_snp_diagnostic.pdf"

# ── 1. Load all-GWAS COLOC results ────────────────────────────────────────────
cat("Loading susie_coloc_all_gwas.csv ...\n")
all_coloc <- fread(file.path(RESULTS_DIR, "susie_coloc_all_gwas.csv"))
cat(sprintf("  %d gene-GWAS pairs, %d unique genes, %d GWAS\n",
            nrow(all_coloc),
            n_distinct(all_coloc$gene),
            n_distinct(all_coloc$gwas_name)))

# ── 2. SNP count distribution across ALL tested windows ───────────────────────
cat("\n=== SNP count distribution (all tested windows) ===\n")
quants <- quantile(all_coloc$n_snps, probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99))
print(round(quants))

cat(sprintf("\nWindows with n_snps < 10:  %d (%.1f%%)\n",
            sum(all_coloc$n_snps < 10),  mean(all_coloc$n_snps < 10)*100))
cat(sprintf("Windows with n_snps < 50:  %d (%.1f%%)\n",
            sum(all_coloc$n_snps < 50),  mean(all_coloc$n_snps < 50)*100))
cat(sprintf("Windows with n_snps < 100: %d (%.1f%%)\n",
            sum(all_coloc$n_snps < 100), mean(all_coloc$n_snps < 100)*100))
cat(sprintf("Windows with n_snps < 200: %d (%.1f%%)\n",
            sum(all_coloc$n_snps < 200), mean(all_coloc$n_snps < 200)*100))

# ── 3. Sensitivity: hit counts at each PP.H4 threshold × n_snps threshold ─────
cat("\n=== Sensitivity analysis: hits surviving each (PP.H4, n_snps) threshold ===\n")

pp_thresholds   <- c(0.5, 0.8, 0.9)
snp_thresholds  <- c(10, 50, 100, 200, 500)

sensitivity <- expand.grid(pp4 = pp_thresholds, min_snps = snp_thresholds) |>
  rowwise() |>
  mutate(
    n_gene_gwas_pairs = sum(all_coloc$PP.H4.abf >= pp4 & all_coloc$n_snps >= min_snps),
    n_unique_genes    = n_distinct(all_coloc$gene[all_coloc$PP.H4.abf >= pp4 & all_coloc$n_snps >= min_snps])
  ) |>
  ungroup()

# Pivot to readable table
for (pp in pp_thresholds) {
  cat(sprintf("\n  PP.H4 > %.1f\n", pp))
  sub <- sensitivity[sensitivity$pp4 == pp, c("min_snps", "n_unique_genes", "n_gene_gwas_pairs")]
  # Add pct retained relative to min_snps=10
  base_genes <- sub$n_unique_genes[sub$min_snps == 10]
  sub$pct_genes_retained <- round(sub$n_unique_genes / base_genes * 100, 1)
  print(as.data.frame(sub), row.names = FALSE)
}

# ── 4. SNP count distribution WITHIN top hits ─────────────────────────────────
cat("\n=== SNP counts for high-confidence hits (PP.H4 > 0.5) ===\n")
top_hits <- all_coloc[all_coloc$PP.H4.abf >= 0.5, ]
cat(sprintf("  %d gene-GWAS pairs with PP.H4 > 0.5\n", nrow(top_hits)))
cat(sprintf("  %d unique genes\n", n_distinct(top_hits$gene)))
cat("  n_snps quantiles:\n")
print(round(quantile(top_hits$n_snps, probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99))))

cat(sprintf("\n  PP.H4>0.5 hits with n_snps <  50: %d (%.1f%%)\n",
            sum(top_hits$n_snps < 50),  mean(top_hits$n_snps < 50)*100))
cat(sprintf("  PP.H4>0.5 hits with n_snps < 100: %d (%.1f%%)\n",
            sum(top_hits$n_snps < 100), mean(top_hits$n_snps < 100)*100))
cat(sprintf("  PP.H4>0.5 hits with n_snps < 200: %d (%.1f%%)\n",
            sum(top_hits$n_snps < 200), mean(top_hits$n_snps < 200)*100))

# ── 5. Genes lost at n_snps >= 100 that have PP.H4 > 0.5 ─────────────────────
cat("\n=== Genes with PP.H4 > 0.5 that would be LOST at n_snps >= 100 ===\n")
lost_genes <- top_hits[top_hits$n_snps < 100, ] |>
  group_by(gene) |>
  summarise(
    best_pp4   = max(PP.H4.abf),
    max_n_snps = max(n_snps),
    n_gwas     = n(),
    gwas_names = paste(gwas_name, collapse = "; ")
  ) |>
  arrange(desc(best_pp4))

cat(sprintf("  %d unique genes would be lost (have PP.H4>0.5 only in windows with <100 SNPs)\n",
            nrow(lost_genes)))
if (nrow(lost_genes) > 0) {
  cat("  Top 20 by PP.H4:\n")
  print(as.data.frame(head(lost_genes, 20)), row.names = FALSE)
}

# ── 6. Per-GWAS breakdown ─────────────────────────────────────────────────────
cat("\n=== Per-GWAS: median n_snps and hit counts at PP.H4 > 0.5 ===\n")
per_gwas <- all_coloc |>
  group_by(gwas_name) |>
  summarise(
    n_windows    = n(),
    median_snps  = median(n_snps),
    p10_snps     = quantile(n_snps, 0.10),
    hits_pp4_05  = sum(PP.H4.abf >= 0.5),
    hits_pp4_05_snps100 = sum(PP.H4.abf >= 0.5 & n_snps >= 100),
    pct_retained = round(hits_pp4_05_snps100 / pmax(hits_pp4_05, 1) * 100, 1)
  ) |>
  arrange(desc(hits_pp4_05))
print(as.data.frame(per_gwas), row.names = FALSE)

# ── 7. Figures ────────────────────────────────────────────────────────────────
cat("\nGenerating figures ...\n")
dir.create(dirname(FIGURE_OUT), recursive = TRUE, showWarnings = FALSE)

pdf(FIGURE_OUT, width = 12, height = 10)

# Panel A: n_snps distribution for all windows vs top hits
df_plot <- rbind(
  data.frame(n_snps = all_coloc$n_snps, group = "All windows"),
  data.frame(n_snps = top_hits$n_snps,  group = "PP.H4 > 0.5")
)
p1 <- ggplot(df_plot, aes(x = pmin(n_snps, 1000), fill = group)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(50, 100, 200), linetype = "dashed",
             color = c("orange", "red", "darkred")) +
  annotate("text", x = c(55, 105, 205), y = Inf, vjust = 1.5,
           label = c("50", "100", "200"), color = c("orange", "red", "darkred"), size = 3) +
  scale_x_continuous("n_snps in COLOC window (capped at 1000)") +
  scale_fill_manual(values = c("All windows" = "steelblue", "PP.H4 > 0.5" = "firebrick")) +
  labs(title = "SNP count distribution in COLOC windows",
       subtitle = "Dashed lines: candidate minimum thresholds",
       fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
print(p1)

# Panel B: Sensitivity — unique genes retained vs min_snps threshold
sens_long <- sensitivity |>
  mutate(pp4_label = paste0("PP.H4 > ", pp4))
p2 <- ggplot(sens_long, aes(x = factor(min_snps), y = n_unique_genes,
                             group = pp4_label, color = pp4_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_label(aes(label = n_unique_genes), show.legend = FALSE, size = 3, nudge_y = 2) +
  scale_color_manual(values = c("PP.H4 > 0.5" = "steelblue",
                                 "PP.H4 > 0.8" = "darkorange",
                                 "PP.H4 > 0.9" = "firebrick")) +
  labs(x = "Minimum n_snps threshold",
       y = "Unique COLOC genes",
       title = "Sensitivity of COLOC gene set to minimum SNP threshold",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
print(p2)

# Panel C: Per-GWAS median SNP count
p3 <- ggplot(per_gwas, aes(x = reorder(gwas_name, median_snps), y = median_snps)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = c(100, 200), linetype = "dashed",
             color = c("red", "darkred")) +
  coord_flip() +
  labs(x = NULL, y = "Median n_snps per window",
       title = "Per-GWAS median SNP count in COLOC windows",
       subtitle = "Red dashed: 100 and 200 thresholds") +
  theme_bw(base_size = 9)
print(p3)

dev.off()
cat(sprintf("Figures saved to %s\n", FIGURE_OUT))
cat("\nDone.\n")
