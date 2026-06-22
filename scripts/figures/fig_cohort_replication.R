#!/usr/bin/env Rscript
# Cross-cohort replication analysis for MASLD DEGs
# Shows how gene counts change as we vary the minimum number of cohorts
# in which a gene must be significantly DE (padj < 0.1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(project_root, "scripts/figures/publication_theme.R"))
source(file.path(project_root, "scripts/figures/load_figure_data.R"))

out_dir   <- FIGS01_DIR
panel_dir <- file.path(out_dir, "panels")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load per-study DE results — ONLY mega-eligible cohorts (Disease vs Control)
# =============================================================================
# Of the 9 cohorts in the cohort presentation, 4 are excluded from this
# Disease-vs-Control replication panel (different contrast, no controls):
# GSE167523, GSE174478, GSE193066, GSE240729.
# PRJNA512027 (Gerhard 2018) is excluded from the cohort presentation
# entirely: L0/S0 library-prep batch is perfectly confounded with diagnosis
# and produces 45,146 DEGs (85% of genes), which inflates replication counts
# artificially. The pipeline still uses it in fibrosis-vs-healthy.
mega_cohorts <- c("GSE126848", "GSE135251", "GSE130970",
                  "GSE213621", "GSE162694")
n_cohorts <- length(mega_cohorts)

per_study_dir <- file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")
files <- list.files(per_study_dir, pattern = "_de_results\\.csv$", full.names = TRUE)
cat("Loading per-study DE files...\n")
all_de <- rbindlist(lapply(files, fread))
all_de <- all_de[dataset %in% mega_cohorts]
cat("  Total rows (6 mega cohorts):", nrow(all_de), "\n")
cat("  Unique genes:", uniqueN(all_de$gene), "\n")
cat("  Cohorts:", paste(sort(unique(all_de$dataset)), collapse = ", "), "\n")

# =============================================================================
# 2. Load dream mega-analysis results
# =============================================================================
dream <- fread(file.path(project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
cat("Dream genes:", nrow(dream), "| DEGs (padj<0.1):", sum(dream$padj < 0.1, na.rm = TRUE), "\n")

# =============================================================================
# 3. Per-gene replication stats
# =============================================================================
gene_rep <- all_de[, .(
  n_tested     = .N,
  n_sig        = sum(adj.P.Val < 0.1, na.rm = TRUE),
  n_up         = sum(adj.P.Val < 0.1 & logFC > 0, na.rm = TRUE),
  n_down       = sum(adj.P.Val < 0.1 & logFC < 0, na.rm = TRUE),
  n_nominal    = sum(P.Value < 0.05, na.rm = TRUE),
  n_nominal_up = sum(P.Value < 0.05 & logFC > 0, na.rm = TRUE),
  n_nominal_dn = sum(P.Value < 0.05 & logFC < 0, na.rm = TRUE),
  mean_logFC   = mean(logFC, na.rm = TRUE)
), by = gene]

# Merge with dream
gene_rep <- merge(gene_rep,
  dream[, .(gene, bulk_logFC = logFC, bulk_padj = padj)],
  by = "gene", all.x = TRUE)
gene_rep[, bulk_sig := !is.na(bulk_padj) & bulk_padj < 0.1]

cat("\n--- Replication summary ---\n")
cat("Genes tested in all", n_cohorts, "cohorts:", sum(gene_rep$n_tested == n_cohorts), "\n")
cat("Genes tested in >=4 cohorts:", sum(gene_rep$n_tested >= 4), "\n")
cat("Dream DEGs significant in 0 per-study cohorts:",
    sum(gene_rep$bulk_sig & gene_rep$n_sig == 0), "\n")

# =============================================================================
# 4. Panel A: Cumulative gene counts at each threshold
# =============================================================================
thresholds <- 1:n_cohorts

# Adjusted p-value threshold
cum_adj <- rbindlist(lapply(thresholds, function(t) {
  data.table(
    min_cohorts = t,
    direction = c("Upregulated", "Downregulated", "Either direction"),
    n_genes = c(
      sum(gene_rep$n_up >= t),
      sum(gene_rep$n_down >= t),
      sum(gene_rep$n_sig >= t)
    ),
    threshold = "Per-study padj < 0.1"
  )
}))

# Nominal p-value threshold
cum_nom <- rbindlist(lapply(thresholds, function(t) {
  data.table(
    min_cohorts = t,
    direction = c("Upregulated", "Downregulated", "Either direction"),
    n_genes = c(
      sum(gene_rep$n_nominal_up >= t),
      sum(gene_rep$n_nominal_dn >= t),
      sum(gene_rep$n_nominal >= t)
    ),
    threshold = "Per-study P < 0.05"
  )
}))

cum_all <- rbind(cum_adj, cum_nom)

# Print table
cat("\n--- Cumulative gene counts (padj < 0.1) ---\n")
dcast_tbl <- dcast(cum_adj, min_cohorts ~ direction, value.var = "n_genes")
print(dcast_tbl)

# Panel A plot — padj < 0.1 only
pA <- ggplot(cum_adj, aes(x = min_cohorts, y = n_genes, color = direction)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  # Only label "Either direction" to avoid overlap at crowded points
  geom_text(data = cum_adj[min_cohorts %in% c(1, 3, 5) & direction == "Either direction"],
            aes(label = format(n_genes, big.mark = ",")),
            vjust = -0.8, hjust = 0.5, size = PUB_GEOM_TEXT, show.legend = FALSE,
            fontface = "bold") +
  scale_x_continuous(breaks = 1:n_cohorts, labels = 1:n_cohorts) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c(
    "Upregulated"      = fig1_colors$up,
    "Downregulated"    = fig1_colors$down,
    "Either direction" = "#424242"
  )) +
  labs(
    x = "Minimum cohorts with significant DE",
    y = "Number of genes",
    title = "Genes reaching per-study significance across cohorts",
    subtitle = "Per-study padj < 0.1 | No logFC threshold",
    color = "Direction"
  ) +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom")

# =============================================================================
# 5. Panel B: Marginal genes gained per step decrease
# =============================================================================
marginal <- cum_adj[direction != "Either direction"][order(direction, min_cohorts)]
marginal[, genes_gained := c(NA, -diff(n_genes)), by = direction]

# From threshold N to N-1: how many genes are gained
pB <- ggplot(marginal[!is.na(genes_gained)],
             aes(x = factor(min_cohorts), y = genes_gained, fill = direction)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = format(genes_gained, big.mark = ",")),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = PUB_GEOM_TEXT, fontface = "bold") +
  scale_fill_manual(values = c(
    "Upregulated"   = fig1_colors$up,
    "Downregulated" = fig1_colors$down
  )) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Cohort threshold relaxed from N \u2192 N\u22121",
    y = "Genes gained",
    title = "Marginal genes gained per threshold decrease",
    subtitle = "padj < 0.1 | Each bar = additional genes from relaxing by one cohort",
    fill = "Direction"
  ) +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom")

# =============================================================================
# 6. Panel C: Dream DEG replication distribution (Fig 1F).
# Uses the canonical project DEG thresholds — same convention as Fig 1h:
#   integrated: padj < 0.05 AND |logFC| > 0.5  (consensus_degs.csv definition)
#   per-study:  padj < 0.05 AND same direction as integrated
# Note: panels A/B/D above retain padj < 0.1 (no LFC, no direction) on purpose;
# they answer different questions and live in the supplementary figure.
# =============================================================================
PADJ_INT_C    <- 0.05
LFC_INT_C     <- 0.5
PADJ_COHORT_C <- 0.05

# Per-gene direction-aware concordance count using per-study DE.
# We must recompute from raw per-study rows because gene_rep$n_up/n_down
# count direction within a cohort, not concordance with the integrated effect.
ps_with_dream <- merge(
  all_de[, .(gene, dataset, ps_logFC = logFC, ps_padj = adj.P.Val)],
  dream[, .(gene, bulk_logFC = logFC, bulk_padj = padj)],
  by = "gene", all.x = FALSE
)
ps_with_dream[, concordant_sig := !is.na(ps_padj) & ps_padj < PADJ_COHORT_C]
fig1f_concord <- ps_with_dream[, .(n_cohorts_concordant = sum(concordant_sig, na.rm = TRUE)),
                               by = gene]

# Integrated DEGs at canonical threshold
dream_genes_C <- merge(
  dream[!is.na(padj) & padj < PADJ_INT_C & abs(logFC) > LFC_INT_C,
        .(gene, bulk_logFC = logFC, bulk_padj = padj)],
  fig1f_concord, by = "gene", all.x = TRUE)
dream_genes_C[is.na(n_cohorts_concordant), n_cohorts_concordant := 0L]
dream_genes_C[, dir := fifelse(bulk_logFC > 0,
                               "Integrated upregulated",
                               "Integrated downregulated")]

cat("\nFig 1F integrated DEGs (padj <", PADJ_INT_C, ", |logFC| >", LFC_INT_C,
    "):", nrow(dream_genes_C), "\n")
cat("  Replication counts (per-study padj <", PADJ_COHORT_C, "):\n")
print(dream_genes_C[, .N, by = n_cohorts_concordant][order(n_cohorts_concordant)])

# Distribution for stacked bar
rep_dist_C <- dream_genes_C[, .(count = .N), by = .(n_cohorts_concordant, dir)]
totals_C   <- dream_genes_C[, .(total = .N), by = n_cohorts_concordant]

pC <- ggplot(rep_dist_C,
             aes(x = factor(n_cohorts_concordant), y = count, fill = dir)) +
  geom_col(position = position_stack(), width = 0.7) +
  geom_text(
    data = totals_C,
    aes(x = factor(n_cohorts_concordant), y = total,
        label = format(total, big.mark = ","), fill = NULL),
    vjust = -0.4, size = 2.0, fontface = "bold", inherit.aes = FALSE
  ) +
  # Liang et al. 2025 Fig 1E palette: soft magenta + warm peach.
  # Lower-chroma than the volcano's #C9265E / #1565C0 pair so the panel
  # reads quieter at fig1f scale.
  scale_fill_manual(values = c(
    "Integrated upregulated"   = "#D7367A",
    "Integrated downregulated" = "#FFB481"
  )) +
  scale_x_discrete(limits = as.character(0:n_cohorts)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = paste0("Cohorts significant (padj < ", PADJ_COHORT_C, ")"),
    y = "Number of integrated DEGs",
    title = "Per-study replication of integrated DEGs",
    subtitle = paste0(
      format(nrow(dream_genes_C), big.mark = ","),
      " integrated DEGs (padj < ", PADJ_INT_C,
      ", |logFC| > ", LFC_INT_C, ")"
    ),
    fill = "Direction"
  ) +
  theme_masld() + theme_pub() +
  # Bump fonts a notch above theme_pub baseline so labels read at fig1f scale.
  theme(legend.position = "bottom",
        plot.title    = element_text(size = PUB_TITLE + 1, face = "bold"),
        plot.subtitle = element_text(size = PUB_SUBTITLE + 1, color = "gray30"),
        axis.title    = element_text(size = PUB_AXIS_TITLE + 1),
        axis.text     = element_text(size = PUB_AXIS_TEXT + 2),
        legend.title  = element_text(size = PUB_LEGEND_TIT + 1, face = "bold"),
        legend.text   = element_text(size = PUB_LEGEND + 2))

# =============================================================================
# 7. Panel D: Dream-rescued genes (dream-sig but per-study non-sig)
# =============================================================================
# Genes significant in dream but in 0 individual cohorts
rescued <- gene_rep[bulk_sig == TRUE & n_sig == 0]
cat("Dream-rescued genes (0 per-study sig):", nrow(rescued), "\n")

# For these genes, show nominal replication
rescued_dist <- rescued[, .(count = .N), by = n_nominal]

pD <- ggplot(rescued_dist, aes(x = factor(n_nominal), y = count)) +
  geom_col(fill = fig1_colors$mixed, width = 0.7) +
  geom_text(aes(label = count), vjust = -0.4, size = PUB_GEOM_TEXT, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Number of cohorts with nominal significance (P < 0.05)",
    y = "Number of genes",
    title = "Integration-rescued genes: nominal replication",
    subtitle = paste0(format(nrow(rescued), big.mark = ","),
                      " integrated DEGs with 0 per-study padj < 0.1"),
    fill = NULL
  ) +
  theme_masld() + theme_pub()

# =============================================================================
# 8. Composite figure — supplementary composite uses A/B/D only.
# Panel C (replication_dist) is the canonical Fig 1F and is rendered separately
# (see Fig 1 panel promotion below); including it here would duplicate that
# panel within the supplement.
# =============================================================================
fig <- (pA | pB) / (pD | plot_spacer()) +
  plot_annotation(
    title = "Cross-cohort replication of MASLD disease-vs-control DEGs",
    subtitle = paste0(n_cohorts, " Disease-vs-Control cohorts | ",
                      format(uniqueN(gene_rep$gene), big.mark = ","),
                      " genes tested | Panel C = Fig 1F"),
    tag_levels = list(c("A", "B", "C")),
    theme = theme(
      plot.title    = element_text(size = PUB_TITLE + 1, face = "bold", family = "Helvetica"),
      plot.subtitle = element_text(size = PUB_SUBTITLE, color = "gray35", family = "Helvetica")
    )
  )

ggsave(file.path(out_dir, "cohort_replication.pdf"), fig,
       width = 9, height = 7.5, device = cairo_pdf)
cat("\nSaved:", file.path(out_dir, "cohort_replication.pdf"), "\n")

# Individual panels → panels/ subdirectory
# Panel C (replication_dist) is NOT written here — it is the canonical fig1f
# and is promoted to FIG1_DIR/panels/fig1f.pdf below. Writing a second copy
# in this supp panels dir caused the same panel to appear in two manuscript
# figures.
ggsave(file.path(panel_dir, "cohort_replication_A_cumulative.pdf"), pA, width = 5, height = 4, device = cairo_pdf)
ggsave(file.path(panel_dir, "cohort_replication_B_marginal.pdf"), pB, width = 5, height = 4, device = cairo_pdf)
ggsave(file.path(panel_dir, "cohort_replication_D_rescued.pdf"), pD, width = 5, height = 4, device = cairo_pdf)
cat("Saved individual panels to:", panel_dir, "\n")

# Promote pC to main Fig 1F (per replication being a headline integration claim).
# Save both the rendered PDF and the ggplot RDS so fig1_atlas_overview_v2.R can
# embed the same panel in the combined Fig 1.
fig1_panel_dir <- file.path(FIG1_DIR, "panels")
dir.create(fig1_panel_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig1_panel_dir, "fig1f.pdf"), pC,
       width = 3.6, height = 3.0, device = cairo_pdf)
saveRDS(pC, file.path(out_dir, "cohort_replication_pC.rds"))
cat("Promoted pC -> ", file.path(fig1_panel_dir, "fig1f.pdf"), "\n", sep = "")

# =============================================================================
# 9. Summary table output
# =============================================================================
summary_table <- cum_adj[, .(min_cohorts, direction, n_genes)]
summary_table <- dcast(summary_table, min_cohorts ~ direction, value.var = "n_genes")
setnames(summary_table, c("min_cohorts", "Down", "Either", "Up"))
summary_table[, `:=`(
  pct_dream_up = round(Up / sum(dream$padj < 0.1 & dream$logFC > 0, na.rm = TRUE) * 100, 1),
  pct_dream_dn = round(Down / sum(dream$padj < 0.1 & dream$logFC < 0, na.rm = TRUE) * 100, 1)
)]
fwrite(summary_table, file.path(out_dir, "replication_summary.csv"))
cat("\nSaved: replication_summary.csv\n")

# Per-gene replication table
fwrite(gene_rep, file.path(out_dir, "per_gene_replication.csv"))
cat("Saved: per_gene_replication.csv\n")

cat("\nDone.\n")
