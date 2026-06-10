#!/usr/bin/env Rscript
# figS_healthy_control_audit.R
#
# Supplementary figure: healthy-control characterization (sensitivity + resilience).
#
# Eight panels:
#   (a) Cohort flow + 3-way table
#   (b) Subclinical screen: NMF P_Pro_inflammatory_score × NAS for control pool, tertiles
#   (c) progression-metric DE robustness: original vs cleaned controls scatter, per-cohort LOO
#   (d) Effect-size attenuation barplot (sex / COLOC / plasma)
#   (e) Resilient transcriptome volcanoes (vs lean-healthy + vs obese-MASLD)
#   (f) Pathway NES heatmap
#   (g) Olink progression-metric Z boxplot
#   (h) Genetic anchor: 3-gene expression + TWAS concordance
#
# Outputs PDF panels into figures/supplementary/figS_methods_validation/healthy_control_audit/.
#
# Spec: docs/superpowers/specs/2026-04-27-healthy-control-audit-design.md
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

HCDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/healthy_control_audit")
DEDIR <- file.path(HCDIR, "resilience_de")
GENDIR <- file.path(HCDIR, "gwas_overlap")
OUTDIR <- FIGS_HCAUDIT_DIR
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Color palette for the 3 groups (consistent across panels)
GRP_COLORS <- c(lean_healthy = "#377eb8",
                resilient    = "#4daf4a",
                obese_MASLD  = "#e41a1c")
GRP_LABELS <- c(lean_healthy = "Lean healthy",
                resilient    = "Resilient",
                obese_MASLD  = "Obese-MASLD")

theme_panel <- theme_classic(base_size = 9) +
  theme(legend.position = "right",
        plot.title = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))

# ---------------------------------------------------------------------------
# Panel (a): Cohort flow + 3-way table
# ---------------------------------------------------------------------------
cat("[Panel a] Cohort flow + 3-way table\n")
defs <- fread(file.path(HCDIR, "controls_definitions.csv"))
tbl_a <- defs[is_resilient | is_lean_healthy | is_obese_MASLD,
               .N, by = .(dataset, group3 = fcase(
                 is_resilient,    "resilient",
                 is_lean_healthy, "lean_healthy",
                 is_obese_MASLD,  "obese_MASLD"))]
tbl_a[, group3 := factor(group3, levels = c("lean_healthy", "resilient", "obese_MASLD"))]

p_a <- ggplot(tbl_a, aes(x = dataset, y = N, fill = group3)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = GRP_COLORS, labels = GRP_LABELS) +
  labs(title = "(a) 3-way pool composition by cohort",
       x = NULL, y = "Subjects", fill = "Group") +
  theme_panel +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTDIR, "panel_a_cohort_breakdown.pdf"), p_a,
       width = 4.5, height = 3, device = cairo_pdf)

# ---------------------------------------------------------------------------
# Panel (b): Subclinical screen scatter
# ---------------------------------------------------------------------------
cat("[Panel b] Subclinical screen scatter\n")
sc <- fread(file.path(HCDIR, "subclinical_screen.csv"))
# Color by p1_tertile
sc[, p1_tertile := factor(p1_tertile,
                            levels = c("clean", "intermediate", "suspect", "unscored"))]
TERT_COLORS <- c(clean = "#41ab5d", intermediate = "#ffd92f",
                  suspect = "#e08214", unscored = "grey70")
p_b <- ggplot(sc[!is.na(P_Pro_inflammatory_score)],
              aes(x = nas_score, y = P_Pro_inflammatory_score, color = p1_tertile)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.85, size = 1.4) +
  scale_color_manual(values = TERT_COLORS, drop = FALSE) +
  facet_wrap(~ dataset, nrow = 2) +
  labs(title = "(b) Subclinical screen of histology-healthy controls",
       x = "NAS score (clinical)", y = "NMF P_Pro_inflammatory_score",
       color = "Tertile") +
  theme_panel
ggsave(file.path(OUTDIR, "panel_b_subclinical_screen.pdf"), p_b,
       width = 5.5, height = 3, device = cairo_pdf)

# ---------------------------------------------------------------------------
# Panel (c): progression-metric DE robustness scatter
# ---------------------------------------------------------------------------
cat("[Panel c] progression-metric DE robustness scatter\n")
sens_fp <- file.path(HCDIR, "sensitivity/sensitivity_logFC_scatter.csv")
sens_summary_fp <- file.path(HCDIR, "sensitivity/sensitivity_summary.csv")
if (file.exists(sens_fp)) {
  sens <- fread(sens_fp)
  # Top 1000 by |orig logFC| for the focal scatter
  sens_top <- sens[order(-abs(logFC_orig))][1:min(1000, .N)]
  rho <- if (file.exists(sens_summary_fp)) {
    s <- fread(sens_summary_fp)
    s[metric == "spearman_rho_top1000", as.numeric(value)]
  } else NA

  p_c <- ggplot(sens_top, aes(x = logFC_orig, y = logFC_clean)) +
    geom_hex(bins = 50) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_fill_continuous(type = "viridis") +
    annotate("text", x = -Inf, y = Inf,
              label = sprintf("Spearman ρ = %.3f", rho),
              hjust = -0.1, vjust = 1.5, size = 3) +
    labs(title = sprintf("(c) F2+ vs F0 logFC: original vs cleaned controls\nTop 1000 by |orig logFC|"),
         x = "logFC (orig: all F0 controls)",
         y = "logFC (clean: subclinical-suspect dropped)") +
    theme_panel
  ggsave(file.path(OUTDIR, "panel_c_progression_robustness.pdf"), p_c,
         width = 4.5, height = 3.5, device = cairo_pdf)
} else {
  # Placeholder — atlas-only
  cat("  [warn] Script 232 sensitivity output not present; using atlas placeholder\n")
  atlas_f2 <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
                     select = c("ensembl_id", "human_symbol",
                                "f2_inflection_logFC", "f2_inflection_padj"))
  atlas_sig <- atlas_f2[!is.na(f2_inflection_padj) &
                          f2_inflection_padj < 0.05 &
                          abs(f2_inflection_logFC) > 0.3]
  p_c <- ggplot(atlas_sig, aes(x = f2_inflection_logFC)) +
    geom_histogram(fill = "#377eb8", alpha = 0.7, bins = 60) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "(c) Placeholder: progression-metric logFC distribution",
         x = "Progression-inflection logFC", y = "Genes") +
    theme_panel
  ggsave(file.path(OUTDIR, "panel_c_progression_robustness.pdf"), p_c,
         width = 4.5, height = 3, device = cairo_pdf)
}

# ---------------------------------------------------------------------------
# Panel (d): Effect-size attenuation barplot
# ---------------------------------------------------------------------------
cat("[Panel d] Effect-size attenuation barplot\n")
if (file.exists(sens_summary_fp)) {
  s <- fread(sens_summary_fp)
  # Build attenuation data from 232 outputs
  rho_top   <- s[metric == "spearman_rho_top1000", as.numeric(value)]
  rho_all   <- s[metric == "spearman_rho_all", as.numeric(value)]
  med_atten <- s[metric == "median_attenuation_ratio", as.numeric(value)]
  retention <- s[metric == "retention_rate", as.numeric(value)] / 100
  n_sig_orig <- s[metric == "n_sig_orig", as.numeric(value)]
  n_sig_clean <- s[metric == "n_sig_clean", as.numeric(value)]
  attn <- data.table(
    metric  = c("Spearman ρ (top 1000)",
                "Spearman ρ (all)",
                "Median |LFC clean/orig|",
                "Sig-gene retention",
                "FDR<0.05 retained / orig"),
    value   = c(rho_top, rho_all, med_atten, retention, n_sig_clean / max(n_sig_orig, 1)),
    target  = c(0.85, 0.85, 0.9, 0.8, 0.8)
  )
  fwrite(attn, file.path(OUTDIR, "panel_d_attenuation_data.csv"))

  p_d <- ggplot(attn, aes(x = metric, y = value)) +
    geom_col(fill = "#984ea3", width = 0.6) +
    geom_hline(aes(yintercept = target), color = "red", linetype = "dashed") +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.5, size = 2.5) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(title = "(d) F2+vs-F0 DE robustness to control cleaning",
         x = NULL, y = "Value (red = robustness target)") +
    theme_panel +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
  ggsave(file.path(OUTDIR, "panel_d_effect_attenuation.pdf"), p_d,
         width = 4, height = 3, device = cairo_pdf)
} else {
  # Placeholder
  cat("  [warn] sensitivity_summary.csv not present; using placeholder\n")
  attn <- data.table(metric = "pending", value = NA_real_)
  p_d <- ggplot(attn, aes(x = metric, y = value)) +
    geom_col(fill = "grey80") +
    labs(title = "(d) Placeholder: pending Script 232",
         x = NULL, y = NULL) +
    theme_panel
  ggsave(file.path(OUTDIR, "panel_d_effect_attenuation.pdf"), p_d,
         width = 3, height = 2, device = cairo_pdf)
}

# ---------------------------------------------------------------------------
# Panel (e): Resilient transcriptome volcanoes
# ---------------------------------------------------------------------------
cat("[Panel e] Resilient transcriptome volcanoes\n")
de_files <- list(
  GSE126848  = file.path(DEDIR, "resilient_de_pairwise_GSE126848.csv"),
  CrossCohort = file.path(DEDIR, "resilient_de_pairwise_crosscohort.csv")
)
volcano_plots <- list()
for (label in names(de_files)) {
  fp <- de_files[[label]]
  if (!file.exists(fp)) {
    cat(sprintf("  [skip] %s pending\n", fp))
    next
  }
  d <- fread(fp)
  for (cn in unique(d$contrast)) {
    dd <- d[contrast == cn]
    dd[, neglog10p := -log10(pmax(adj.P.Val, 1e-300))]
    dd[, sig := !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 0.5]
    p_v <- ggplot(dd, aes(x = logFC, y = neglog10p, color = sig)) +
      geom_point(alpha = 0.4, size = 0.5) +
      scale_color_manual(values = c(`FALSE` = "grey80", `TRUE` = "#e41a1c"),
                          guide = "none") +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey30") +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey30") +
      labs(title = sprintf("(e) %s: %s", label, cn),
           x = "logFC", y = "-log10 FDR") +
      theme_panel
    volcano_plots[[paste(label, cn, sep = "_")]] <- p_v
  }
}
if (length(volcano_plots) > 0) {
  pe <- wrap_plots(volcano_plots, ncol = 2)
  ggsave(file.path(OUTDIR, "panel_e_volcanos.pdf"), pe,
         width = 7, height = 2.5 * ceiling(length(volcano_plots) / 2),
         device = cairo_pdf)
}

# ---------------------------------------------------------------------------
# Panel (f): Pathway NES heatmap
# ---------------------------------------------------------------------------
cat("[Panel f] Pathway NES heatmap\n")
gsea_fp <- file.path(DEDIR, "resilient_gsea.csv")
if (file.exists(gsea_fp)) {
  gs <- fread(gsea_fp)
  gs[, contrast_short := sub("group", "", contrast)]
  gs[, comparison := paste(analysis, contrast_short, sep = " | ")]

  # Top pathways by total significance (sum -log10 padj across comparisons)
  gs[, abs_NES := abs(NES)]
  setorder(gs, -abs_NES)
  top_paths <- gs[padj < 0.05, .SD[1:1], by = pathway]
  top_paths <- top_paths[order(-abs_NES)][1:min(40, .N)]

  gs_top <- gs[pathway %in% top_paths$pathway]
  gs_top[, pathway_short := sub("HALLMARK_", "", pathway)]
  gs_top[, pathway_short := sub("REACTOME_", "R: ", pathway_short)]
  gs_top[, pathway_short := substr(pathway_short, 1, 50)]

  p_f <- ggplot(gs_top, aes(x = comparison, y = reorder(pathway_short, NES, FUN = mean),
                              fill = NES)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#377eb8", mid = "white", high = "#e41a1c",
                          midpoint = 0, limits = c(-3, 3), oob = squish) +
    geom_text(aes(label = ifelse(padj < 0.05, "*", "")), size = 2.5) +
    labs(title = "(f) Pathway NES across resilient contrasts",
         x = NULL, y = NULL, fill = "NES") +
    theme_panel +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))
  ggsave(file.path(OUTDIR, "panel_f_gsea_heatmap.pdf"), p_f,
         width = 6, height = 7, device = cairo_pdf)
} else {
  cat("  [skip] resilient_gsea.csv pending\n")
}

# ---------------------------------------------------------------------------
# Panel (g): Olink progression-metric Z boxplot
# ---------------------------------------------------------------------------
cat("[Panel g] Olink progression-metric Z boxplot\n")
# Path retained for upstream CSV compatibility; in-script identifiers neutralized.
olink_fp <- file.path(HCDIR, "resilient_olink_f2switch.csv")
if (file.exists(olink_fp)) {
  ok <- fread(olink_fp)
  ok_labeled <- ok[!is.na(disease_broad)]
  ok_labeled[, disease_label := factor(
    disease_broad,
    levels = c("Healthy", "CLD_early", "CLD_F3", "CLD_F4"),
    labels = c("Healthy\n(n=40)", "CLD F0-2\n(n=39)", "CLD F3\n(n=24)", "CLD F4\n(n=114)")
  )]
  # Upstream CSV column (legacy name on disk) is renamed to neutral identifier here.
  setnames(ok_labeled, old = grep("switch_z$", names(ok_labeled), value = TRUE),
           new = "progression_z")

  # KW p-value
  kw_p <- kruskal.test(progression_z ~ disease_broad, data = ok_labeled)$p.value
  kw_lbl <- sprintf("KW p = %.3f", kw_p)

  pal_g <- c("Healthy\n(n=40)"  = "#2c7bb6",
             "CLD F0-2\n(n=39)" = "#74add1",
             "CLD F3\n(n=24)"   = "#f46d43",
             "CLD F4\n(n=114)"  = "#d73027")

  p_g <- ggplot(ok_labeled[!is.na(disease_label)],
                aes(x = disease_label, y = progression_z, fill = disease_label)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
    scale_fill_manual(values = pal_g, guide = "none") +
    annotate("text", x = 2.5, y = max(ok_labeled$progression_z, na.rm = TRUE) * 0.92,
             label = kw_lbl, size = 3) +
    labs(title = "(g) Plasma progression score by liver disease stage",
         subtitle = "Yang et al. 2025 (GSE276114); n=217 subjects",
         x = NULL, y = "Progression Z-score") +
    theme_panel
  ggsave(file.path(OUTDIR, "panel_g_olink_progression.pdf"), p_g,
         width = 5, height = 3.5, device = cairo_pdf)
}

# ---------------------------------------------------------------------------
# Panel (h): Genetic anchor — 3-gene boxplot + TWAS concordance
# ---------------------------------------------------------------------------
cat("[Panel h] Genetic anchor\n")
hand_fp <- file.path(GENDIR, "hand_curated_expression_long.csv")
gen_fp <- file.path(GENDIR, "resilient_genetic_enrichment.csv")
if (file.exists(hand_fp)) {
  hd <- fread(hand_fp)
  hd[, group3 := factor(group3, levels = c("lean_healthy", "resilient", "obese_MASLD"))]
  p_h1 <- ggplot(hd, aes(x = group3, y = log2_cpm, fill = group3)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~ symbol, scales = "free_y") +
    scale_fill_manual(values = GRP_COLORS, labels = GRP_LABELS, guide = "none") +
    labs(title = "(h.i) Curated protective alleles: bulk expression",
         x = NULL, y = "log2 CPM") +
    theme_panel +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(OUTDIR, "panel_h_protective_genes.pdf"), p_h1,
         width = 5, height = 3, device = cairo_pdf)
}
if (file.exists(gen_fp)) {
  ge <- fread(gen_fp)
  cat("Genetic enrichment summary:\n")
  print(ge)
  fwrite(ge, file.path(OUTDIR, "panel_h_genetic_enrichment_table.csv"))
}

# ---------------------------------------------------------------------------
# Composite figure (all panels)
# ---------------------------------------------------------------------------
cat("\nAttempting composite figure...\n")
panels <- list()
if (exists("p_a")) panels$a <- p_a
if (exists("p_b")) panels$b <- p_b
if (exists("p_c")) panels$c <- p_c
if (exists("p_d")) panels$d <- p_d
if (length(volcano_plots) > 0) panels$e <- wrap_plots(volcano_plots, ncol = 2)
if (exists("p_f")) panels$f <- p_f
if (exists("p_g")) panels$g <- p_g
if (exists("p_h1")) panels$h <- p_h1

if (length(panels) >= 4) {
  comp <- wrap_plots(panels, ncol = 2)
  ggsave(file.path(OUTDIR, "figS_healthy_control_audit_composite.pdf"),
         comp, width = 11, height = 14, device = cairo_pdf, limitsize = FALSE)
  cat("  Composite written.\n")
}

cat("\nAll panels written to:", OUTDIR, "\n")
cat("Done.\n")
