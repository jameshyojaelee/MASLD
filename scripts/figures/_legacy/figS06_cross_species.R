#!/usr/bin/env Rscript
# ==========================================================================
# figS06_cross_species.R — Supplementary Figure 6: Cross-Species Concordance
# Comprehensive view of 5 mouse diet models + cross-species concordance
#
# Panels:
#   a: Per-diet DEG counts (bidirectional bar)
#   b: LFC correlation heatmap across diet models
#   c: Diet model concordance with human dream (scatter per diet)
#   d: Consensus tier composition (stacked bar showing Tier1/Tier2/NS per diet)
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS06_DIR, "figS06_cross_species.pdf")
dir.create(file.path(FIGS06_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# Load mouse data
# ==========================================================================
mouse_dir <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results")

# Per-diet DE
per_diet_files <- list.files(file.path(mouse_dir, "per_diet"),
                             pattern = "_de_results\\.csv$|^[A-Z]+\\.csv$",
                             full.names = TRUE)
# Filter out summary
per_diet_files <- per_diet_files[!grepl("summary", per_diet_files)]

per_diet <- rbindlist(lapply(per_diet_files, function(f) {
  dt <- fread(f)
  if (!"diet_model" %in% names(dt)) {
    dm <- gsub("(_de_results)?\\.csv$", "", basename(f))
    dt[, diet_model := dm]
  }
  dt
}), fill = TRUE)

if ("adj.P.Val" %in% names(per_diet)) setnames(per_diet, "adj.P.Val", "padj")

# Mouse consensus
consensus_path <- file.path(mouse_dir, "mouse_consensus_degs.csv")
consensus <- if (file.exists(consensus_path)) fread(consensus_path) else NULL

# Mouse meta per diet
meta_path <- file.path(mouse_dir, "meta_analysis/meta_per_diet.csv")
meta <- if (file.exists(meta_path)) fread(meta_path) else NULL

# Human dream for concordance
dream <- load_dream_results()

# ==========================================================================
# Panel a: Per-diet DEG counts
# ==========================================================================
diets <- unique(per_diet$diet_model)
diets <- diets[diets %in% names(diet_colors)]

deg_counts <- per_diet[diet_model %in% diets, .(
  n_up = sum(padj < 0.05 & logFC > 0, na.rm = TRUE),
  n_down = sum(padj < 0.05 & logFC < 0, na.rm = TRUE)
), by = diet_model]

deg_long <- melt(deg_counts, id.vars = "diet_model",
                 variable.name = "direction", value.name = "n")
deg_long[, direction := fifelse(direction == "n_up", "Up", "Down")]
deg_long[direction == "Down", n := -n]

totals <- deg_counts[, .(total = n_up + n_down), by = diet_model][order(total)]
deg_long[, diet_model := factor(diet_model, levels = totals$diet_model)]

p_a <- ggplot(deg_long, aes(x = n, y = diet_model, fill = direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(Up = masld_colors$up, Down = masld_colors$down)) +
  scale_x_continuous(labels = function(x) format(abs(x), big.mark = ",")) +
  theme_masld() +
  theme(legend.position = "bottom") +
  labs(x = "Number of DEGs (padj < 0.05)", y = NULL, fill = NULL,
       title = "Per-diet model DE")

# ==========================================================================
# Panel b: LFC correlation heatmap across diets
# ==========================================================================
# Build gene x diet LFC matrix
diet_lfc <- dcast(per_diet[diet_model %in% diets],
                  gene ~ diet_model, value.var = "logFC", fun.aggregate = mean)
genes <- diet_lfc$gene
diet_mat <- as.matrix(diet_lfc[, -1])
rownames(diet_mat) <- genes

# Correlation matrix
cor_mat <- cor(diet_mat, use = "pairwise.complete.obs", method = "spearman")
cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("Diet1", "Diet2", "rho"))

p_b <- ggplot(cor_dt, aes(x = Diet1, y = Diet2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                       midpoint = 0, limits = c(-1, 1), name = expression(rho)) +
  coord_fixed() +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Cross-diet LFC correlation")

# ==========================================================================
# Panel c: Each diet vs human (C2 canonical) — ortholog-joined LFC concordance
# ==========================================================================
# C2 fix (2026-06-21): canonical_deg_results.csv carries no mouse_symbol and the
# mouse per-diet tables carry no human symbol, so the old shared-symbol merge
# dropped to "mapping not available". Join through the 1:1 ortholog map
# (replicates concordance script 01 / fig4g): human C2 logFC -> mouse_gene_id ->
# mouse per-diet logFC.
ortho_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/ortholog_mapping.tsv")
if (!is.null(dream) && !is.null(per_diet) && "bulk_logFC" %in% names(dream) && file.exists(ortho_path)) {
  ortho <- fread(ortho_path)
  hd <- copy(dream); hd[, gene_base := gsub("\\..*", "", gene)]
  h_mapped <- merge(hd[, .(gene_base, human_logFC = bulk_logFC)],
                    ortho[, .(human_gene_id, mouse_gene_id)],
                    by.x = "gene_base", by.y = "human_gene_id")
  pd <- copy(per_diet[diet_model %in% diets]); pd[, mouse_base := gsub("\\..*", "", gene)]
  diet_human <- merge(pd[, .(mouse_base, logFC, diet_model)],
                      h_mapped[, .(mouse_gene_id, human_logFC)],
                      by.x = "mouse_base", by.y = "mouse_gene_id", allow.cartesian = TRUE)

  if (nrow(diet_human) > 0) {
    cors <- diet_human[, .(
      rho = cor(logFC, human_logFC, use = "complete.obs", method = "spearman"),
      n = .N
    ), by = diet_model]

    p_c <- ggplot(diet_human, aes(x = human_logFC, y = logFC)) +
      rasterize_layer(geom_point(size = 0.1, alpha = 0.15, color = "gray50", shape = 16)) +
      geom_smooth(method = "lm", linewidth = 0.4, color = masld_colors$up, se = FALSE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  linewidth = 0.2, color = "gray40") +
      geom_text(data = cors, aes(label = sprintf("rho=%.2f\nn=%s", rho, format(n, big.mark = ","))),
                x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 2, inherit.aes = FALSE) +
      facet_wrap(~diet_model, nrow = 1) +
      coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
      theme_masld() +
      labs(x = expression("Human integrated log"[2]*"FC (C2)"),
           y = expression("Mouse diet log"[2]*"FC"),
           title = "Mouse-human LFC concordance by diet model")
  } else {
    p_c <- placeholder("No overlapping orthologs between mouse and human")
  }
} else {
  p_c <- placeholder("Ortholog map or human/mouse C2 data missing")
}

# ==========================================================================
# Panel d: Consensus tier composition
# ==========================================================================
# C2 fix (2026-06-21): the consensus column is `tier` (Significant/Not_significant),
# not the legacy `consensus_tier`; accept either so the panel renders a real bar.
tier_col <- if (!is.null(consensus)) intersect(c("consensus_tier", "tier"), names(consensus))[1] else NA_character_
if (!is.null(consensus) && !is.na(tier_col)) {
  tier_summary <- consensus[, .N, by = c(tier_col)]
  setnames(tier_summary, tier_col, "grp")
  tier_summary[, pct := N / sum(N) * 100]

  p_d <- ggplot(tier_summary, aes(x = reorder(grp, -N), y = N, fill = grp)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = paste0(format(N, big.mark = ","), "\n(",
                                 sprintf("%.1f%%", pct), ")")),
              vjust = -0.3, size = 2.2, color = "black") +
    scale_fill_manual(values = c(Significant = masld_colors$up,
                                 Not_significant = "#9E9E9E",
                                 Consensus = tryCatch(masld_colors$conserved,
                                                      error = function(e) "#1B7C6F")),
                      guide = "none", na.value = "#9E9E9E") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_masld() +
    labs(x = NULL, y = "Mouse genes",
         title = "Mouse consensus significance (across diet models)")
} else if (!is.null(meta)) {
  # Fallback: show per-diet significant gene counts
  meta_sig <- meta[, .(
    n_sig = sum(meta_padj < 0.1, na.rm = TRUE),
    n_total = .N
  )]
  p_d <- placeholder(paste0("Meta: ", meta_sig$n_sig, " / ", meta_sig$n_total, " DEGs"))
} else {
  p_d <- placeholder("Mouse consensus data not available")
}

# ==========================================================================
# Assemble
# ==========================================================================
fig <- ((p_a | p_b) / (p_c) / (p_d)) +
  plot_layout(heights = c(1.2, 1, 0.8)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, OUT, width = fig_full_width, height = 9)
message("Mouse overview figure saved to ", OUT)
