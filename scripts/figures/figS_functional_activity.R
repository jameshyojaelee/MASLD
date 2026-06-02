#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Functional Activity and TF-Mediated Mechanistic Paths
#
# Panels:
#   a: TF activity bar plot (top 15 activated + 15 repressed, ULM scores)
#   b: PROGENy pathway activity (14 pathways, colored by significance)
#   c: decoupleR vs SCENIC+ TF activity concordance scatter
#   d: Statistically significant bridging TFs (GWAS enrichment)
#   e: Top GWAS→TF→DE mechanistic paths (Sankey-like flow)
#   f: Drug target TF regulatory context
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

FUNC_DIR <- file.path(ME, "functional_activity")
COSMOS_DIR <- file.path(ME, "cosmos_mechanistic")
OUT <- file.path(FIGS08_DIR, "figS_functional_activity.pdf")
dir.create(file.path(FIGS08_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== Supplementary Figure: Functional Activity ===\n")

# ==========================================================================
# Load data
# ==========================================================================
tf_scores <- fread(file.path(FUNC_DIR, "tf_activity_scores.csv"))
pw_scores <- fread(file.path(FUNC_DIR, "pathway_activity_scores.csv"))
scenic_comp <- fread(file.path(FUNC_DIR, "decoupler_vs_scenic_comparison.csv"))
tf_bridge <- fread(file.path(COSMOS_DIR, "tf_gwas_de_bridge.csv"))
tf_paths <- fread(file.path(COSMOS_DIR, "gwas_tf_de_paths.csv"))
tf_conv <- fread(file.path(FUNC_DIR, "tf_convergence_summary.csv"))

cat("Data loaded:\n")
cat("  TF scores:", nrow(tf_scores), "TFs\n")
cat("  Pathway scores:", nrow(pw_scores), "pathways\n")
cat("  SCENIC+ comparison:", nrow(scenic_comp), "TFs\n")
cat("  TF bridges:", nrow(tf_bridge), "TFs\n")
cat("  Mechanistic paths:", nrow(tf_paths), "\n\n")

# ==========================================================================
# Panel (a): TF activity bar plot — top 15 activated + 15 repressed
# ==========================================================================
cat("Panel a: TF activity bar plot...\n")

# Select top 15 activated and 15 repressed by absolute score
top_act <- tf_scores[direction == "activated"][order(-abs_score)][1:15]
top_rep <- tf_scores[direction == "repressed"][order(-abs_score)][1:15]
top_tfs <- rbind(top_rep, top_act)
top_tfs[, tf := factor(tf, levels = rev(top_tfs$tf))]
top_tfs[, sig_label := fifelse(padj < 0.05, "*", "")]
top_tfs[, sig_label := fifelse(padj < 0.001, "***",
                        fifelse(padj < 0.01, "**", sig_label))]

pa <- ggplot(top_tfs, aes(x = score, y = tf, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sig_label,
                x = score + sign(score) * 0.3),
            size = 2, hjust = ifelse(top_tfs$score > 0, 0, 1)) +
  scale_fill_manual(values = c("activated" = masld_colors$up,
                                "repressed" = masld_colors$down),
                     name = "Direction") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(x = "TF activity score (ULM)", y = NULL,
       title = "Disease-dysregulated TFs (decoupleR)") +
  theme_masld() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 7))

# ==========================================================================
# Panel (b): Pathway activity — all 14 PROGENy pathways
# ==========================================================================
cat("Panel b: Pathway activity...\n")

pw_scores[, pathway := factor(pathway, levels = pathway[order(score)])]
pw_scores[, sig := padj < 0.05]

# Color by significance: significant = red/blue by direction, NS = gray
pw_scores[, fill_group := fifelse(sig == TRUE, direction, "ns")]

pb <- ggplot(pw_scores, aes(x = score, y = pathway)) +
  geom_col(aes(fill = fill_group), width = 0.7) +
  scale_fill_manual(values = c("activated" = masld_colors$up,
                                "repressed" = masld_colors$down,
                                "ns" = "gray75"),
                     labels = c("Activated (sig)", "Repressed (sig)", "Not significant"),
                     name = NULL) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(x = "Pathway activity score (ULM)", y = NULL,
       title = "Signaling pathway activities (PROGENy)") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 7)) +
  annotate("text", x = max(pw_scores$score) * 0.7, y = 3,
           label = paste0(pw_scores[padj < 0.05, .N], "/14 sig\n(padj < 0.05)"),
           size = 2, color = "gray40")

# ==========================================================================
# Panel (c): decoupleR vs SCENIC+ concordance
# ==========================================================================
cat("Panel c: decoupleR vs SCENIC+ comparison...\n")

rho_val <- cor(scenic_comp$decoupler_score, scenic_comp$scenic_activity_diff,
               method = "spearman", use = "complete.obs")
n_conc <- sum(sign(scenic_comp$decoupler_score) == sign(scenic_comp$scenic_activity_diff))
binom_p <- binom.test(n_conc, nrow(scenic_comp), p = 0.5)$p.value

pc <- ggplot(scenic_comp, aes(x = scenic_activity_diff, y = decoupler_score)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray80") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray80") +
  geom_point(aes(color = sign(decoupler_score) == sign(scenic_activity_diff)),
             size = 2, shape = 16) +
  ggrepel::geom_text_repel(aes(label = tf), size = 2, max.overlaps = 15) +
  scale_color_manual(values = c("TRUE" = masld_colors$conserved,
                                 "FALSE" = masld_colors$up),
                      labels = c("Discordant", "Concordant"),
                      name = "Direction") +
  labs(x = "SCENIC+ regulon activity difference",
       y = "decoupleR TF activity (ULM)",
       title = "Bulk-inferred vs chromatin-level TF activity") +
  annotate("text", x = min(scenic_comp$scenic_activity_diff) * 0.8,
           y = max(scenic_comp$decoupler_score) * 0.9,
           label = paste0("rho = ", round(rho_val, 2), "\n",
                          n_conc, "/", nrow(scenic_comp), " concordant\n",
                          "binomial p = ", signif(binom_p, 2)),
           size = 2, hjust = 0) +
  theme_masld() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 7))

# ==========================================================================
# Panel (d): Statistically significant bridging TFs
# ==========================================================================
cat("Panel d: Bridging TFs...\n")

# Show top 20 bridging TFs ranked by GWAS enrichment
bridge_top <- tf_bridge[bridges_gwas_to_de == TRUE][order(gwas_hypergeom_pval)][1:min(20, sum(tf_bridge$bridges_gwas_to_de))]
bridge_top[, tf := factor(tf, levels = rev(tf))]
bridge_top[, sig := gwas_hypergeom_padj < 0.1]

pd <- ggplot(bridge_top, aes(x = gwas_enrichment, y = tf)) +
  geom_col(aes(fill = sig), width = 0.6) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("TRUE" = masld_colors$up,
                                "FALSE" = "gray70"),
                     labels = c("padj >= 0.1", "padj < 0.1"),
                     name = "GWAS enrichment") +
  geom_text(aes(label = paste0(n_gwas_targets, "/", round(expected_gwas)),
                x = gwas_enrichment + 0.15),
            size = 1.8, hjust = 0, color = "black") +
  labs(x = "GWAS target fold enrichment (obs/exp annotated)",
       y = NULL,
       title = "TFs bridging GWAS-causal to DE genes") +
  theme_masld() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 7)) +
  coord_cartesian(xlim = c(0, max(bridge_top$gwas_enrichment) * 1.3))

# ==========================================================================
# Panel (e): Top mechanistic paths (dot plot)
# ==========================================================================
cat("Panel e: Mechanistic paths...\n")

# Show top 15 paths by TF activity × n_de_targets
path_top <- tf_paths[order(-abs(tf_activity) * n_de_targets)][1:min(15, nrow(tf_paths))]
# Shorten path descriptions
path_top[, short_desc := paste0(
  sub(";.*", " +", gwas_gene), " → ", tf, " → ",
  sub(";.*", " +", top_de_targets))]
path_top[, short_desc := factor(short_desc, levels = rev(short_desc))]

pe <- ggplot(path_top, aes(x = n_de_targets, y = short_desc)) +
  geom_point(aes(size = abs(tf_activity),
                  color = tf_activity > 0), shape = 16) +
  scale_color_manual(values = c("TRUE" = masld_colors$up,
                                 "FALSE" = masld_colors$down),
                      labels = c("Repressed TF", "Activated TF"),
                      name = "TF direction") +
  scale_size_continuous(range = c(1.5, 5), name = "|TF activity|") +
  labs(x = "DE targets downstream of TF", y = NULL,
       title = expression("GWAS" %->% "TF" %->% "DE mechanistic paths")) +
  theme_masld() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 5),
        plot.title = element_text(size = 7))

# ==========================================================================
# Panel (f): Drug target TF regulatory context
# ==========================================================================
cat("Panel f: Drug target TF context...\n")

drug_targets <- c("THRB", "NR1H4", "PPARA", "PPARG", "PNPLA3", "TM6SF2", "GLP1R")
drug_conv <- tf_conv[tf %in% drug_targets |
                      (n_druggable_targets > 3 & convergence_layers >= 4)]

# If drug targets not in TF list, show top convergent TFs instead
if (nrow(drug_conv) < 5) {
  drug_conv <- tf_conv[convergence_layers >= 4][order(-n_druggable_targets)][1:min(15, sum(tf_conv$convergence_layers >= 4))]
}

if (nrow(drug_conv) > 0) {
  drug_conv[, tf := factor(tf, levels = rev(tf[order(n_druggable_targets)]))]

  pf <- ggplot(drug_conv, aes(x = n_druggable_targets, y = tf)) +
    geom_segment(aes(xend = 0, yend = tf), linewidth = 0.3, color = "gray80") +
    geom_point(aes(size = n_deg_targets, color = direction), alpha = 0.8, shape = 16) +
    geom_text(aes(label = n_gwas_targets, x = n_druggable_targets + 0.5),
              size = 2, hjust = 0) +
    scale_color_manual(values = c("activated" = masld_colors$up,
                                   "repressed" = masld_colors$down),
                        name = "TF direction") +
    scale_size_continuous(range = c(2, 6), name = "DE targets") +
    labs(x = "Druggable targets in regulon",
         y = NULL,
         title = "High-convergence TFs: druggable + GWAS context") +
    theme_masld() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 7)) +
    annotate("text", x = max(drug_conv$n_druggable_targets) * 0.5,
             y = 1.5, label = "numbers = GWAS targets",
             size = 1.6, color = "gray50", fontface = "italic")
} else {
  pf <- placeholder("No high-convergence TFs")
}

# ==========================================================================
# Assemble figure
# ==========================================================================
cat("Assembling figure...\n")

fig <- (pa | pb) /
       (pc | pd) /
       (pe | pf) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(fig, OUT, width = fig_full_width, height = 9)
cat("Saved:", OUT, "\n")

# Also save individual panels for flexibility
panel_dir <- file.path(FIGS08_DIR, "panels")
save_fig(pa, file.path(panel_dir, "panel_a_tf_activity.pdf"), width = fig_half_width, height = 3.5)
save_fig(pb, file.path(panel_dir, "panel_b_pathway_activity.pdf"), width = fig_half_width, height = 2.5)
save_fig(pc, file.path(panel_dir, "panel_c_scenic_comparison.pdf"), width = fig_half_width, height = 3)
save_fig(pd, file.path(panel_dir, "panel_d_bridging_tfs.pdf"), width = fig_half_width, height = 3.5)
save_fig(pe, file.path(panel_dir, "panel_e_mechanistic_paths.pdf"), width = fig_half_width, height = 3.5)
save_fig(pf, file.path(panel_dir, "panel_f_tf_convergence.pdf"), width = fig_half_width, height = 3)

cat("\nDone. All panels saved to:", panel_dir, "\n")
