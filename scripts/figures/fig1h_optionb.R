#!/usr/bin/env Rscript
# fig1h_optionb.R — Alternative Fig 1h framing.
#
# Hypothesis: integration-only DEGs (significant in dream but missed by every
# cohort's canonical DEG threshold) carry real biological signal. If they do,
# the per-cohort logFC signs should *agree* with the integrated direction far
# more often than for a random non-DEG control set.
#
# Comparison:
#   Group A: Integration-only — dream padj < 0.05 AND no cohort meets canonical
#            (padj < 0.05 AND |logFC| > 0.5 AND same direction as dream).
#   Group B: Random non-DEG — random size-matched sample of dream padj >= 0.05.
#
# Metric per gene: count of 5 cohorts where sign(per-cohort logFC) ==
#                  sign(dream logFC). Range 0..5.
#
# Plot: grouped bars (% of genes per concordance level) for the two groups,
#       plus 4+/5 fraction comparison and Fisher's exact p-value.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Outputs (fig1h_optionb.pdf + fig1h_optionb_table.csv) relocated to the
# Figure-3 RNA-seq dir (FIG2_DIR = figures/main/fig3_RNAseq, back-compat constant name).
PANEL_DIR <- file.path(FIG2_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_INT      <- 0.05
PADJ_COHORT   <- 0.05
LFC_CANONICAL <- 0.5
SET_SEED      <- 42

FIVE_COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

# ----------------------------------------------------------------------------
# Load + merge
# ----------------------------------------------------------------------------
message("Loading data...")
# Per-cohort sign concordance is computed from A2's OWN percohort_lvqw files plus
# the canonical integrated direction — this reproduces A2's n_concord EXACTLY
# (verified 100% per-gene identical; 46.7% concordant in >=4/5). The newer
# standalone per-study CSVs (QW-flipped 2026-06-11) drift to ~52%, so we
# deliberately do NOT use load_per_study_de() here.
A2_DIR <- file.path(INT_RESULTS, "integration_only_C2")
can <- load_dream_results()                        # canonical_deg_results.csv (C2)
can[, gene_clean := sub("\\..*", "", gene)]
sgn <- can[, .(gene_clean, symbol, bulk_logFC, bulk_padj,
               int_sign = sign(bulk_logFC))]
# AveExpr (for A2's AveExpr-decile-matched control) read straight from canonical
ave <- fread(file.path(INT_RESULTS, "canonical_deg_results.csv"),
             select = c("gene", "AveExpr"))
ave[, gene_clean := sub("\\..*", "", gene)]
sgn <- merge(sgn, ave[, .(gene_clean, AveExpr)], by = "gene_clean", all.x = TRUE)
for (co in FIVE_COHORTS) {
  pc <- fread(file.path(A2_DIR, sprintf("percohort_lvqw_%s.csv", co)))
  pc[, gene_clean := sub("\\..*", "", gene)]
  pc <- pc[, .(gene_clean, s = sign(logFC))]
  setnames(pc, "s", co)
  sgn <- merge(sgn, pc, by = "gene_clean", all.x = TRUE)
}
sign_mat <- as.matrix(sgn[, ..FIVE_COHORTS])
sgn[, n_sign_concordant := rowSums(sign_mat == int_sign, na.rm = TRUE)]
panel_df <- sgn

# ----------------------------------------------------------------------------
# Define groups
# ----------------------------------------------------------------------------
# Group A: the CANONICAL integration-only set (A2_integration_only_genes_C2.csv,
# 5,926 genes) with A2's own per-cohort concordance — matches the manuscript
# exactly (integration_only_n = 5926, dirconcord_io_pct = 46.7).
io_gids <- sub("\\..*", "", fread(file.path(A2_DIR,
                  "A2_integration_only_genes_C2.csv"))$gid)
group_A <- panel_df[gene_clean %in% io_gids]
cat(sprintf("[sanity] group_A n=%d (A2 canonical=5926); >=4/5 concordant=%.1f%% (A2=46.7%%)\n",
            nrow(group_A), 100 * mean(group_A$n_sign_concordant >= 4)))

# Group B: AveExpr-DECILE-size-matched non-DEGs — replicates A2's matched control
# exactly (A2 matches on AveExpr deciles of the io set -> 16.1% / OR 4.55; a plain
# random non-DEG sample undershoots to ~12.6%).
non_deg_pool <- panel_df[!is.na(bulk_padj) & bulk_padj >= PADJ_INT & !is.na(AveExpr)]
io_ave <- group_A$AveExpr
brk    <- unique(quantile(io_ave, probs = seq(0, 1, 0.1), na.rm = TRUE))
io_bin <- cut(io_ave, breaks = brk, include.lowest = TRUE)
ns_bin <- cut(non_deg_pool$AveExpr, breaks = brk, include.lowest = TRUE)
set.seed(SET_SEED)
matched <- character(0)
for (b in levels(io_bin)) {
  need <- sum(io_bin == b, na.rm = TRUE)
  pool <- non_deg_pool$gene_clean[which(ns_bin == b)]
  if (length(pool) == 0 || need == 0) next
  matched <- c(matched, sample(pool, need, replace = length(pool) < need))
}
group_B <- non_deg_pool[match(matched, gene_clean)]

message(sprintf("Group A (Integration-only): %s genes", comma(nrow(group_A))))
message(sprintf("Group B (Random non-DEG):   %s genes (sampled from %s pool)",
                comma(nrow(group_B)), comma(nrow(non_deg_pool))))

group_A[, group := "Integration-only DEGs"]
group_B[, group := "Size-matched non-DEGs"]
combined <- rbindlist(list(group_A, group_B), use.names = TRUE)

# ----------------------------------------------------------------------------
# Tabulate concordance distributions
# ----------------------------------------------------------------------------
combined[, group := factor(group, levels = c("Integration-only DEGs", "Size-matched non-DEGs"))]
dist_dt <- combined[, .N, by = .(group, n_sign_concordant)]
dist_dt[, total := sum(N), by = group]
dist_dt[, pct := 100 * N / total]
setorder(dist_dt, group, n_sign_concordant)
print(dcast(dist_dt, n_sign_concordant ~ group, value.var = "pct"))

# 4+/5 fraction
high_conc <- combined[, .(n_high = sum(n_sign_concordant >= 4),
                          n_low  = sum(n_sign_concordant <  4)),
                      by = group]
print(high_conc)

# Fisher's exact test
mat <- as.matrix(high_conc[, .(n_high, n_low)])
rownames(mat) <- high_conc$group
fisher_p  <- fisher.test(mat)$p.value
fisher_or <- fisher.test(mat)$estimate
message(sprintf("Fisher's exact (4+/5 vs <4/5): OR=%.2f, p=%.2e",
                fisher_or, fisher_p))

frac_int_high <- high_conc[group == "Integration-only DEGs",     n_high / (n_high + n_low)]
frac_rnd_high <- high_conc[group == "Size-matched non-DEGs",   n_high / (n_high + n_low)]
n_hi_A <- high_conc[group == "Integration-only DEGs", n_high]
n_hi_B <- high_conc[group == "Size-matched non-DEGs", n_high]

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
group_colors <- c(
  "Integration-only DEGs"      = "#C9265E",
  "Size-matched non-DEGs"    = "#9E9E9E"
)

# Ensure all 6 levels (0..5) appear
all_levels <- CJ(group = levels(combined$group),
                 n_sign_concordant = 0:5)
dist_full <- merge(all_levels, dist_dt[, .(group, n_sign_concordant, N)],
                   by = c("group", "n_sign_concordant"), all.x = TRUE)
dist_full[is.na(N), N := 0]
dist_full[, x_factor := factor(n_sign_concordant, levels = 0:5)]

p_str <- if (fisher_p < 1e-300) "p<1e-300" else
  paste0("p=", formatC(fisher_p, format = "e", digits = 1))
label_text <- sprintf(
  "≥4 / 5 cohorts agree: %.1f%% vs %.1f%% (OR=%.1f, %s)",
  100 * frac_int_high, 100 * frac_rnd_high,
  fisher_or, p_str
)

p <- ggplot(dist_full,
            aes(x = x_factor, y = N, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.13))) +
  labs(
    x = "Cohorts with concordant log2FC sign (of 5)",
    y = "Number of genes",
    title = "Integration-rescued DEGs share effect direction across cohorts"
  ) +
  theme_masld(base_size = 10) +
  theme(
    plot.title    = element_text(size = 11, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 8.5, color = "gray30", margin = margin(b = 2)),
    plot.margin   = margin(2, 3, 2, 2),
    axis.title.x  = element_text(margin = margin(t = 1)),
    axis.title.y  = element_text(margin = margin(r = 1)),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA),
    legend.margin = margin(1, 3, 1, 3),
    legend.key.height = unit(0.3, "cm")
  )

save_fig(p, file.path(PANEL_DIR, "integration_discovery_optionb.pdf"),
         width = fig_half_width * 1.55, height = 3.4)

# Per-gene table
fwrite(combined[, .(gene_clean, symbol, bulk_logFC, bulk_padj,
                    n_sign_concordant, group)],
       file.path(PANEL_DIR, "integration_discovery_optionb_table.csv"))

fp <- file.path(PANEL_DIR, "integration_discovery_optionb.pdf")
if (file.exists(fp)) {
  message(sprintf("\nOutput: %s (%s)", fp,
                  utils:::format.object_size(file.size(fp), "auto")))
}
