#!/usr/bin/env Rscript
# 206b_hoang_method_investigation.R
# Investigate why Hoang (GSE130970) concordance is r=0.71-0.77 while
# Govaere (GSE135251) is r=0.93-0.97.
#
# Strategy: reproduce Hoang's exact ordinal regression (CLM) on our
# re-quantified counts, compute their range_log2FC metric, and compare
# to their published results. This isolates method vs quantification.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(edgeR)
  library(ordinal)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INTB <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
INT  <- file.path(INTB, "results")
PUB  <- file.path(BASE, "data/published_degs")
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 206b: Hoang Method Investigation ===\n")

# ── Load data ───────────────────────────────────────────────────────────────
cat("Loading data...\n")
counts <- readRDS(file.path(INT, "integration/merged_counts_raw.rds"))
sample_meta <- readRDS(file.path(INT, "integration/meta_matched.rds"))
qc <- fread(file.path(INTB, "qc/sample_qc_report.csv"))
sample_meta <- sample_meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Gene annotation
annot <- fread(file.path(INT, "gene_annotation/human_ensg_to_symbol.tsv"))
annot_pc <- annot[gene_type == "protein_coding"]
sym2ens <- setNames(annot_pc$gene_base, annot_pc$symbol)
ens2sym <- setNames(annot_pc$symbol, annot_pc$gene_base)

# GSE130970 samples
m130 <- sample_meta[dataset == "GSE130970" & !is.na(nas_score) & !is.na(fibrosis_stage)]
cat(sprintf("  GSE130970 samples: %d\n", nrow(m130)))

# Load published Hoang results
hoang_nas <- as.data.table(read_excel(file.path(PUB, "GSE130970/MOESM2.xlsx"),
                                       sheet = "NAS ordinal regression"))
hoang_fib <- as.data.table(read_excel(file.path(PUB, "GSE130970/MOESM2.xlsx"),
                                       sheet = "fibrosis ordinal regression"))
hoang_nas[, gene_base := sym2ens[gene_symbol]]
hoang_fib[, gene_base := sym2ens[gene_symbol]]

# ── Step 1: Prepare log2CPM matrix (matching Hoang's normalization) ─────────
cat("\nPreparing log2CPM matrix...\n")
idx <- colnames(counts) %in% m130$sample_id
dge <- DGEList(counts = counts[, idx])
dge <- calcNormFactors(dge, method = "TMM")
keep <- filterByExpr(dge, group = factor(m130$nas_score))
dge <- dge[keep, , keep.lib.sizes = FALSE]
log2cpm <- cpm(dge, log = TRUE, prior.count = 1)  # log2(CPM + prior)
cat(sprintf("  %d genes × %d samples\n", nrow(log2cpm), ncol(log2cpm)))

# Match sample order to metadata
m130_ordered <- m130[match(colnames(log2cpm), m130$sample_id)]

# ── Step 2: Compute range_log2FC (Hoang's metric) on our counts ────────────
cat("\nComputing range_log2FC on our counts...\n")

compute_range_log2fc <- function(log2cpm_mat, scores) {
  # "Difference in mean log2(CPM) between top 2 and bottom 2 ordinal levels"
  ulevels <- sort(unique(scores))
  n_lev <- length(ulevels)
  top2 <- ulevels[(n_lev - 1):n_lev]
  bot2 <- ulevels[1:2]
  top_idx <- which(scores %in% top2)
  bot_idx <- which(scores %in% bot2)
  top_mean <- rowMeans(log2cpm_mat[, top_idx, drop = FALSE])
  bot_mean <- rowMeans(log2cpm_mat[, bot_idx, drop = FALSE])
  top_mean - bot_mean
}

our_range_nas <- compute_range_log2fc(log2cpm, m130_ordered$nas_score)
our_range_fib <- compute_range_log2fc(log2cpm, m130_ordered$fibrosis_stage)

# ── Step 3: Run ordinal regression (CLM) on our counts ──────────────────────
cat("\nRunning ordinal regression (CLM) on our counts...\n")
cat("  This matches Hoang's exact statistical method.\n")

run_clm <- function(log2cpm_mat, meta_df, score_col) {
  scores_ordered <- factor(meta_df[[score_col]], ordered = TRUE)
  age <- meta_df$age
  sex <- meta_df$sex

  n_genes <- nrow(log2cpm_mat)
  results <- data.table(
    gene = rownames(log2cpm_mat),
    coefficient = rep(NA_real_, n_genes),
    std_err = rep(NA_real_, n_genes),
    z_score = rep(NA_real_, n_genes),
    p_value = rep(NA_real_, n_genes)
  )

  converged <- 0
  for (i in seq_len(n_genes)) {
    if (i %% 5000 == 0) cat(sprintf("    %d / %d genes...\n", i, n_genes))
    expr_val <- log2cpm_mat[i, ]
    df <- data.frame(score = scores_ordered, expr = expr_val, age = age, sex = sex)
    tryCatch({
      fit <- clm(score ~ expr + age + sex, data = df)
      s <- summary(fit)
      # Extract coefficient for expr
      coef_row <- s$coefficients["expr", , drop = FALSE]
      results[i, `:=`(
        coefficient = coef_row[1, "Estimate"],
        std_err = coef_row[1, "Std. Error"],
        z_score = coef_row[1, "z value"],
        p_value = coef_row[1, "Pr(>|z|)"]
      )]
      converged <- converged + 1
    }, error = function(e) NULL, warning = function(w) NULL)
  }
  cat(sprintf("    Converged: %d / %d (%.1f%%)\n", converged, n_genes,
              100 * converged / n_genes))
  results[, gene_base := sub("\\.\\d+$", "", gene)]
  results[, adj_P := p.adjust(p_value, method = "BH")]
  results
}

cat("  NAS ordinal CLM...\n")
clm_nas <- run_clm(log2cpm, m130_ordered, "nas_score")

cat("  Fibrosis ordinal CLM...\n")
clm_fib <- run_clm(log2cpm, m130_ordered, "fibrosis_stage")

# Add range_log2FC to CLM results
clm_nas[, range_log2FC := our_range_nas[match(gene, rownames(log2cpm))]]
clm_fib[, range_log2FC := our_range_fib[match(gene, rownames(log2cpm))]]

# ── Step 4: Compare to published results ────────────────────────────────────
cat("\n=== Comparison Results ===\n")

compare <- function(pub, ours, pub_lfc, our_lfc, label) {
  m <- merge(pub[!is.na(gene_base), .(gene_base, pub_lfc = get(pub_lfc))],
             ours[!is.na(gene_base), .(gene_base, our_lfc = get(our_lfc))],
             by = "gene_base")
  m <- m[is.finite(pub_lfc) & is.finite(our_lfc)]
  r <- cor(m$pub_lfc, m$our_lfc)
  rho <- cor(m$pub_lfc, m$our_lfc, method = "spearman")
  dir <- mean(sign(m$pub_lfc) == sign(m$our_lfc)) * 100
  cat(sprintf("  %-50s r=%.3f  rho=%.3f  dir=%.1f%%  n=%d\n",
              label, r, rho, dir, nrow(m)))
  list(r = r, rho = rho, dir = dir, n = nrow(m), data = m, label = label)
}

# Load our limma results for reference
limma_fib <- fread(file.path(INT, "disease_signatures/fibrosis_ordinal_per_study.csv"))
limma_fib[, gene_base := sub("\\.\\d+$", "", gene)]
limma_fib <- limma_fib[dataset == "GSE130970"]

cat("\n--- NAS ordinal ---\n")
n1 <- compare(hoang_nas, clm_nas, "range_log2FC", "range_log2FC",
              "Published range_log2FC vs Our CLM range_log2FC")
n2 <- compare(hoang_nas, clm_nas, "coefficient", "coefficient",
              "Published coefficient vs Our CLM coefficient")
n3 <- compare(hoang_nas, clm_nas, "range_log2FC", "coefficient",
              "Published range_log2FC vs Our CLM coefficient")

cat("\n--- Fibrosis ordinal ---\n")
f1 <- compare(hoang_fib, clm_fib, "range_log2FC", "range_log2FC",
              "Published range_log2FC vs Our CLM range_log2FC")
f2 <- compare(hoang_fib, clm_fib, "coefficient", "coefficient",
              "Published coefficient vs Our CLM coefficient")
f3 <- compare(hoang_fib, limma_fib, "range_log2FC", "logFC",
              "Published range_log2FC vs Our limma logFC (baseline)")

cat("\n--- Variance decomposition ---\n")
cat(sprintf("  Govaere (method matched, limma vs limma):     r ~ 0.95\n"))
cat(sprintf("  Hoang range_log2FC same-metric CLM:           r = %.3f  ← quantification-only gap\n", f1$r))
cat(sprintf("  Hoang coefficient same-method CLM:            r = %.3f  ← quantification-only gap\n", f2$r))
cat(sprintf("  Hoang range_log2FC vs limma logFC (baseline):  r = %.3f  ← method + quantification\n", f3$r))
cat(sprintf("  Gap explained by method: %.3f\n", f1$r - f3$r))

# ── Step 5: Discordance analysis ────────────────────────────────────────────
cat("\n=== Discordance Analysis ===\n")

# Using fibrosis range_log2FC comparison (most comparable)
disc <- f1$data[sign(pub_lfc) != sign(our_lfc)]
conc <- f1$data[sign(pub_lfc) == sign(our_lfc)]
cat(sprintf("  Concordant: %d (%.1f%%), Discordant: %d (%.1f%%)\n",
            nrow(conc), 100*nrow(conc)/nrow(f1$data),
            nrow(disc), 100*nrow(disc)/nrow(f1$data)))

# Expression level analysis
f1$data[, symbol := ens2sym[gene_base]]
f1$data[, abs_pub := abs(pub_lfc)]
f1$data[, abs_our := abs(our_lfc)]
f1$data[, concordant := sign(pub_lfc) == sign(our_lfc)]

# Split by expression quartile
f1_with_expr <- merge(f1$data,
  data.table(gene_base = sub("\\.\\d+$", "", rownames(log2cpm)),
             mean_expr = rowMeans(log2cpm)),
  by = "gene_base")
f1_with_expr[, expr_q := cut(mean_expr, quantile(mean_expr, 0:4/4),
                              include.lowest = TRUE,
                              labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))]
disc_by_q <- f1_with_expr[, .(
  n = .N,
  n_disc = sum(!concordant),
  pct_disc = 100 * mean(!concordant)
), by = expr_q]
cat("\n  Direction discordance by expression quartile:\n")
print(disc_by_q[order(expr_q)])

# Top discordant genes
top_disc <- f1_with_expr[concordant == FALSE][order(-abs(pub_lfc - our_lfc))][1:20]
cat("\n  Top 20 most discordant genes (fibrosis range_log2FC):\n")
print(top_disc[, .(symbol, gene_base, pub_lfc = round(pub_lfc, 3),
                    our_lfc = round(our_lfc, 3), mean_expr = round(mean_expr, 1))])

# ── Step 6: Figure ──────────────────────────────────────────────────────────
cat("\nGenerating investigation figure...\n")

make_scatter <- function(data, xlab, ylab, title) {
  r <- cor(data$pub_lfc, data$our_lfc)
  rho <- cor(data$pub_lfc, data$our_lfc, method = "spearman")
  dir <- mean(sign(data$pub_lfc) == sign(data$our_lfc)) * 100
  anno <- sprintf("r = %.3f\nrho = %.3f\ndir = %.1f%%\nn = %s",
                  r, rho, dir, formatC(nrow(data), big.mark = ","))
  ggplot(data, aes(x = pub_lfc, y = our_lfc)) +
    geom_point(alpha = 0.15, size = 0.3, color = masld_colors$ns) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = masld_colors$up) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted",
                color = "grey40") +
    annotate("text", x = -Inf, y = Inf, label = anno,
             hjust = -0.1, vjust = 1.3, size = 2, color = "grey30") +
    labs(x = xlab, y = ylab, title = title) +
    theme_masld()
}

# Panel a: Fib range_log2FC - same metric, CLM vs CLM
pa <- make_scatter(f1$data,
  "Hoang Fib range_log2FC", "Our CLM range_log2FC",
  "Same metric (range_log2FC)")

# Panel b: Fib coefficient - same method, CLM vs CLM
pb <- make_scatter(f2$data,
  "Hoang Fib coefficient", "Our CLM coefficient",
  "Same method (CLM coefficient)")

# Panel c: NAS range_log2FC - same metric
pc <- make_scatter(n1$data,
  "Hoang NAS range_log2FC", "Our CLM range_log2FC",
  "NAS: same metric (range_log2FC)")

# Panel d: NAS coefficient - same method
pd <- make_scatter(n2$data,
  "Hoang NAS coefficient", "Our CLM coefficient",
  "NAS: same method (CLM coefficient)")

# Panel e: Discordance by expression quartile
pe <- ggplot(disc_by_q, aes(x = expr_q, y = pct_disc)) +
  geom_col(fill = masld_colors$up, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d)", pct_disc, n_disc)),
            vjust = -0.2, size = 2) +
  labs(x = "Expression quartile", y = "Direction discordance (%)",
       title = "Discordance by expression level") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_masld()

# Panel f: Variance decomposition summary
decomp <- data.table(
  comparison = c("Govaere\nlimma vs limma", "Hoang CLM\nsame metric",
                 "Hoang CLM\nsame method", "Hoang limma\nbaseline"),
  r = c(0.955, f1$r, f2$r, f3$r),
  driver = c("Quant only", "Quant only", "Quant only", "Method + Quant")
)
decomp[, comparison := factor(comparison, levels = rev(comparison))]
pf <- ggplot(decomp, aes(x = r, y = comparison, fill = driver)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", r)), hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c("Quant only" = masld_colors$down,
                               "Method + Quant" = masld_colors$up),
                    name = "Gap driver") +
  labs(x = "Pearson r", y = NULL, title = "Variance decomposition") +
  xlim(0, 1.1) +
  theme_masld()

fig <- (pa | pb | pc) / (pd | pe | pf) +
  plot_annotation(tag_levels = "a",
                  title = "Hoang concordance investigation") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, file.path(OUTDIR, "hoang_investigation.pdf"),
         width = fig_full_width, height = 5.5)
cat("  Saved hoang_investigation.pdf\n")

# Save detailed results
fwrite(clm_nas[!is.na(coefficient)], file.path(OUTDIR, "clm_nas_results.csv"))
fwrite(clm_fib[!is.na(coefficient)], file.path(OUTDIR, "clm_fib_results.csv"))

# Summary table
summary_dt <- data.table(
  analysis = c("NAS range_log2FC CLM-vs-CLM", "NAS coefficient CLM-vs-CLM",
               "Fib range_log2FC CLM-vs-CLM", "Fib coefficient CLM-vs-CLM",
               "Fib range_log2FC vs limma (baseline)"),
  r = c(n1$r, n2$r, f1$r, f2$r, f3$r),
  rho = c(n1$rho, n2$rho, f1$rho, f2$rho, f3$rho),
  dir_conc = c(n1$dir, n2$dir, f1$dir, f2$dir, f3$dir),
  n = c(n1$n, n2$n, f1$n, f2$n, f3$n)
)
fwrite(summary_dt, file.path(OUTDIR, "hoang_investigation_summary.csv"))
cat("  Saved summary table\n")

cat("\n=== Done ===\n")
