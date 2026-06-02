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

PANEL_DIR <- file.path(FIG1_DIR, "panels")
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
dream <- load_dream_results()
ps    <- load_per_study_de()[dataset %in% FIVE_COHORTS]

dream[, gene_clean := sub("\\..*", "", gene)]
ps[,    gene_clean := sub("\\..*", "", gene)]

# Wide table: one row per gene, cohort logFCs as columns
ps_wide <- dcast(ps[, .(gene_clean, dataset, logFC, padj)],
                 gene_clean ~ dataset, value.var = c("logFC", "padj"))

panel_df <- merge(dream[, .(gene_clean, symbol, dream_logFC, dream_padj)],
                  ps_wide, by = "gene_clean", all.x = TRUE)

# ----------------------------------------------------------------------------
# Per-gene metrics
# ----------------------------------------------------------------------------
lfc_cols  <- paste0("logFC_", FIVE_COHORTS)
padj_cols <- paste0("padj_",  FIVE_COHORTS)

# n_sign_concordant: cohorts with logFC sign matching dream
panel_df[, n_sign_concordant := rowSums(
  sign(.SD[, ..lfc_cols, with = FALSE]) ==
    sign(panel_df$dream_logFC), na.rm = TRUE
)]

# n_cohorts_canonical: cohorts where padj<0.05 AND |logFC|>0.5 AND same dir
canonical_mat <- mapply(function(lfc_col, padj_col) {
  lfc  <- panel_df[[lfc_col]]
  padj <- panel_df[[padj_col]]
  !is.na(lfc) & !is.na(padj) & padj < PADJ_COHORT &
    abs(lfc) > LFC_CANONICAL &
    sign(lfc) == sign(panel_df$dream_logFC) & sign(lfc) != 0
}, lfc_cols, padj_cols, SIMPLIFY = TRUE)
panel_df[, n_cohorts_canonical := rowSums(canonical_mat, na.rm = TRUE)]

# ----------------------------------------------------------------------------
# Define groups
# ----------------------------------------------------------------------------
# Group A: integration-only (sig integrated, no cohort canonical)
group_A <- panel_df[!is.na(dream_padj) & dream_padj < PADJ_INT &
                    n_cohorts_canonical == 0 &
                    !is.na(dream_logFC) & dream_logFC != 0]

# Group B: random non-DEG, size-matched
non_deg_pool <- panel_df[!is.na(dream_padj) & dream_padj >= PADJ_INT &
                         !is.na(dream_logFC) & dream_logFC != 0]

set.seed(SET_SEED)
group_B <- non_deg_pool[sample(.N, min(nrow(group_A), .N))]

message(sprintf("Group A (Integration-only): %s genes", comma(nrow(group_A))))
message(sprintf("Group B (Random non-DEG):   %s genes (sampled from %s pool)",
                comma(nrow(group_B)), comma(nrow(non_deg_pool))))

group_A[, group := "Integration-only (dream padj<0.05, no cohort canonical)"]
group_B[, group := "Non-DEG control (dream padj≥0.05, size-matched)"]
combined <- rbindlist(list(group_A, group_B), use.names = TRUE)

# ----------------------------------------------------------------------------
# Tabulate concordance distributions
# ----------------------------------------------------------------------------
combined[, group := factor(group, levels = c("Integration-only (dream padj<0.05, no cohort canonical)", "Non-DEG control (dream padj≥0.05, size-matched)"))]
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

frac_int_high <- high_conc[group == "Integration-only (dream padj<0.05, no cohort canonical)",     n_high / (n_high + n_low)]
frac_rnd_high <- high_conc[group == "Non-DEG control (dream padj≥0.05, size-matched)",   n_high / (n_high + n_low)]

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
group_colors <- c(
  "Integration-only (dream padj<0.05, no cohort canonical)"      = "#C9265E",
  "Non-DEG control (dream padj≥0.05, size-matched)"    = "#9E9E9E"
)

# Ensure all 6 levels (0..5) appear
all_levels <- CJ(group = levels(combined$group),
                 n_sign_concordant = 0:5)
dist_full <- merge(all_levels, dist_dt[, .(group, n_sign_concordant, pct)],
                   by = c("group", "n_sign_concordant"), all.x = TRUE)
dist_full[is.na(pct), pct := 0]
dist_full[, x_factor := factor(n_sign_concordant, levels = 0:5)]

p_str <- if (fisher_p < 1e-300) "p<1e-300" else
  paste0("p=", formatC(fisher_p, format = "e", digits = 1))
label_text <- sprintf(
  "≥4 / 5 cohorts agree: %.1f%% vs %.1f%% (OR=%.1f, %s)",
  100 * frac_int_high, 100 * frac_rnd_high,
  fisher_or, p_str
)

p <- ggplot(dist_full,
            aes(x = x_factor, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_continuous(labels = function(v) paste0(v, "%"),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = "Cohorts with concordant log2FC sign (of 5)",
    y = "Genes (%)",
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

save_fig(p, file.path(PANEL_DIR, "fig1h_optionb.pdf"),
         width = fig_half_width * 1.55, height = 3.4)

# Per-gene table
fwrite(combined[, .(gene_clean, symbol, dream_logFC, dream_padj,
                    n_sign_concordant, n_cohorts_canonical, group)],
       file.path(PANEL_DIR, "fig1h_optionb_table.csv"))

fp <- file.path(PANEL_DIR, "fig1h_optionb.pdf")
if (file.exists(fp)) {
  message(sprintf("\nOutput: %s (%s)", fp,
                  utils:::format.object_size(file.size(fp), "auto")))
}
