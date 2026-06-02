##############################################################################
# Script 212: Progression Figures + Atlas Integration
#
# Loads pre-computed outputs from Scripts 210 + 211:
#   - progression_coloc_enrichment.csv  (9 rows: per-transition Fisher's exact)
#   - onset_vs_progression_genes.csv    (5,121 rows: gene-level classification)
#   - progression_driver_genetics.csv   (top drivers with causal annotation)
#   - drug_target_progression_classification.csv (drug targets with prog class)
#
# Generates 4 figure panels:
#   Panel A: Per-transition COLOC enrichment (grouped bars, OR + 95% CI)
#   Panel B: Scatter — tau (x) vs COLOC PP.H4 (y), quadrant labels
#   Panel C: Onset vs progression gene comparison + GWAS phenotype
#   Panel D: Drug target progression classification (dot plot)
#
# Atlas integration: adds progression_coloc_class, progression_causal_score,
#   stage_specific_coloc_gwas to multi_evidence_atlas.csv
#
# Outputs:
#   figures/supplementary/figS04_coloc/figS_progression_*.pdf
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (updated)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ---------------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------------
STRAT_RES <- file.path(BASE, "RNA-seq/results/stratified_causal")
PANEL_DIR <- file.path(FIG_OUT, "supplementary/figS04_coloc")
dir.create(STRAT_RES, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 1. Load pre-computed Script 210/211 outputs
# ==========================================================================
cat("=== Script 212: Progression Figures + Atlas Integration ===\n")

# --- Required input files ---
enrich_f     <- file.path(STRAT_RES, "progression_coloc_enrichment.csv")
onset_prog_f <- file.path(STRAT_RES, "onset_vs_progression_genes.csv")
driver_f     <- file.path(STRAT_RES, "progression_driver_genetics.csv")
drug_f       <- file.path(STRAT_RES, "drug_target_progression_classification.csv")

required_files <- c(enrich_f, onset_prog_f, driver_f, drug_f)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Script 212 requires pre-computed outputs from Scripts 210/211.\n",
       "Missing files:\n  ",
       paste(missing, collapse = "\n  "),
       "\nPlease run Scripts 210 and 211 first.")
}

enrich      <- fread(enrich_f)
onset_prog  <- fread(onset_prog_f)
drivers     <- fread(driver_f)
drug_class  <- fread(drug_f)

cat(sprintf("  Enrichment: %d transitions\n", nrow(enrich)))
cat(sprintf("  Onset/progression genes: %d genes\n", nrow(onset_prog)))
cat(sprintf("  Progression drivers: %d drivers\n", nrow(drivers)))
cat(sprintf("  Drug target classification: %d entries\n", nrow(drug_class)))

# --- Also load COLOC gene-level and per-GWAS for atlas integration ---
coloc_f <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
coloc_gl <- NULL
if (file.exists(coloc_f)) {
  coloc_gl <- fread(coloc_f)
  cat(sprintf("  COLOC gene-level: %d genes\n", nrow(coloc_gl)))
}

coloc_gwas_f <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
coloc_gwas <- NULL
if (file.exists(coloc_gwas_f)) {
  coloc_gwas <- fread(coloc_gwas_f)
  cat(sprintf("  Per-GWAS COLOC: %d entries\n", nrow(coloc_gwas)))
}

gm <- load_gene_map()


# ==========================================================================
# 2. Panel A: Per-transition COLOC enrichment (grouped bars, OR + 95% CI)
# ==========================================================================
cat("\n--- Panel A: Per-transition COLOC enrichment ---\n")

# Enrichment file columns: transition, n_transition_genes, n_coloc_genes,
#   n_overlap, odds_ratio, ci_lower, ci_upper, pvalue, universe_size, padj

# Define display order — show transition-style labels (F0->F1, etc.)
# Data has both transition (F0_to_F1) and OVR (F0_vs_rest) formats
trans_display <- c(
  "F0_to_F1"  = "F0\u2192F1",
  "F1_to_F2"  = "F1\u2192F2",
  "F2_to_F3"  = "F2\u2192F3",
  "F3_to_F4"  = "F3\u2192F4",
  "F0_vs_rest" = "F0 (OVR)",
  "F1_vs_rest" = "F1 (OVR)",
  "F2_vs_rest" = "F2 (OVR)",
  "F3_vs_rest" = "F3 (OVR)",
  "F4_vs_rest" = "F4 (OVR)"
)

enrich[, display_label := fifelse(
  transition %in% names(trans_display),
  trans_display[transition],
  transition
)]

# Separate transitions from OVR for cleaner plot
is_transition <- grepl("_to_", enrich$transition)
enrich[, analysis_type := fifelse(grepl("_to_", transition), "Transition", "OVR")]

# Order: transitions first, then OVR
enrich[, display_label := factor(display_label,
  levels = c("F0\u2192F1", "F1\u2192F2", "F2\u2192F3", "F3\u2192F4",
             "F0 (OVR)", "F1 (OVR)", "F2 (OVR)", "F3 (OVR)", "F4 (OVR)"))]

# Significance labels
enrich[, sig_label := ""]
enrich[padj < 0.05, sig_label := "*"]
enrich[padj < 0.01, sig_label := "**"]
enrich[padj < 0.001, sig_label := "***"]

# Handle infinite OR (F1->F2 has 0 overlap, OR=0)
enrich[is.infinite(ci_upper) | ci_upper > 50, ci_upper := 50]

# Color by significance
enrich[, sig_fill := fifelse(padj < 0.05, "FDR < 0.05", "NS")]

p_a <- ggplot(enrich[!is.na(display_label)],
              aes(x = display_label, y = odds_ratio, fill = sig_fill)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = pmin(ci_upper, 5)),
                width = 0.2, linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  geom_text(aes(label = sig_label,
                y = pmin(ci_upper, 5) + 0.08),
            size = 2.5, vjust = 0) +
  scale_fill_manual(values = c("FDR < 0.05" = masld_colors$up,
                               "NS" = masld_colors$ns),
                    name = NULL) +
  coord_cartesian(ylim = c(0, max(2, max(enrich$odds_ratio, na.rm = TRUE) * 1.3))) +
  labs(x = NULL,
       y = "COLOC enrichment (OR)",
       title = "COLOC gene enrichment per fibrosis transition") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")

save_fig(p_a, file.path(PANEL_DIR, "figS_progression_enrichment_panelA.pdf"),
         width = fig_half_width, height = 3.2)
cat("  Saved Panel A\n")


# ==========================================================================
# 3. Panel B: Scatter — tau (x) vs COLOC PP.H4 (y), quadrant labels
# ==========================================================================
cat("\n--- Panel B: Tau vs COLOC PP.H4 scatter ---\n")

# onset_prog columns: symbol, progression_class, coloc_best_pp4, tau,
#   progression_causal_score, gwas_phenotype_class, etc.

scatter_dt <- onset_prog[!is.na(coloc_best_pp4) & coloc_best_pp4 > 0]

# Handle tau: may be NA for some genes
if (!"tau" %in% names(scatter_dt) || all(is.na(scatter_dt$tau))) {
  cat("  WARNING: tau column missing or all NA; using n_stages_sig as proxy\n")
  scatter_dt[, tau := n_stages_sig / 5]
}
scatter_dt <- scatter_dt[!is.na(tau)]

# Quadrant thresholds
tau_thresh <- 0.5
pp4_thresh <- 0.5

# Quadrant colors aligned with progression_class from Script 210
class_colors <- c(
  onset       = "#42A5F5",  # Light blue
  mid_stage   = "#7B1FA2",  # Violet
  progression = "#880E4F",  # Dark magenta
  pan_stage   = "#00695C",  # Teal
  not_causal  = "#E0E0E0"
)

# Label top genes per quadrant (high progression_causal_score)
causal_genes <- scatter_dt[progression_class %in% c("onset", "progression",
                                                     "mid_stage", "pan_stage")]
label_dt <- causal_genes[order(-progression_causal_score)][1:min(25, .N)]

p_b <- ggplot(scatter_dt, aes(x = tau, y = coloc_best_pp4)) +
  # Quadrant shading
  annotate("rect", xmin = tau_thresh, xmax = Inf, ymin = pp4_thresh, ymax = Inf,
           fill = "#F3E5F5", alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = tau_thresh, ymin = pp4_thresh, ymax = Inf,
           fill = "#E8F5E9", alpha = 0.3) +
  # Points
  geom_point(aes(color = progression_class), size = 0.5, alpha = 0.5) +
  # Threshold lines
  geom_vline(xintercept = tau_thresh, linetype = "dashed", color = "gray50",
             linewidth = 0.3) +
  geom_hline(yintercept = pp4_thresh, linetype = "dashed", color = "gray50",
             linewidth = 0.3) +
  # Labels for top genes
  geom_text_repel(data = label_dt, aes(label = symbol, color = progression_class),
                  size = 1.6, max.overlaps = 20, segment.size = 0.15,
                  min.segment.length = 0, show.legend = FALSE) +
  # Quadrant annotations
  annotate("text", x = 0.85, y = 0.97, label = "Stage-specific\ncausal",
           size = 1.8, fontface = "bold", color = "gray30") +
  annotate("text", x = 0.15, y = 0.97, label = "Pan-stage\ncausal",
           size = 1.8, fontface = "bold", color = "gray30") +
  annotate("text", x = 0.85, y = 0.03, label = "Stage-specific\nnon-causal",
           size = 1.6, color = "gray50") +
  annotate("text", x = 0.15, y = 0.03, label = "Pan-stage\nnon-causal",
           size = 1.6, color = "gray50") +
  scale_color_manual(values = class_colors, name = "Class") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression("Stage specificity (" * tau * ")"),
       y = "COLOC PP.H4 (best GWAS)",
       title = "Stage specificity vs genetic causality") +
  theme_masld() +
  theme(legend.position = "right",
        legend.key.size = unit(0.2, "cm"))

save_fig(p_b, file.path(PANEL_DIR, "figS_progression_enrichment_panelB.pdf"),
         width = fig_half_width + 0.5, height = 3.8)
cat("  Saved Panel B\n")


# ==========================================================================
# 4. Panel C: Onset vs progression gene comparison + GWAS phenotype
# ==========================================================================
cat("\n--- Panel C: Onset vs progression comparison ---\n")

# Count COLOC+ genes per progression class
coloc_pos <- onset_prog[is_coloc == TRUE | coloc_best_pp4 > 0.5]
class_counts <- coloc_pos[, .N, by = progression_class]
class_counts <- class_counts[progression_class != "not_causal" & !is.na(progression_class)]

class_order <- c("onset", "mid_stage", "pan_stage", "progression")
class_counts[, progression_class := factor(progression_class,
                                           levels = intersect(class_order, progression_class))]

p_c_bars <- ggplot(class_counts[!is.na(progression_class)],
                   aes(x = progression_class, y = N, fill = progression_class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = N), vjust = -0.3, size = 2) +
  scale_fill_manual(values = class_colors, guide = "none") +
  labs(x = "Progression class",
       y = "N genes (COLOC PP.H4 > 0.5)",
       title = "Gene classification") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# GWAS phenotype breakdown if gwas_phenotype_class is available
p_c_gwas <- NULL
if ("gwas_phenotype_class" %in% names(coloc_pos)) {
  gwas_cross <- coloc_pos[!is.na(progression_class) &
                           progression_class != "not_causal" &
                           !is.na(gwas_phenotype_class),
                          .N, by = .(progression_class, gwas_phenotype_class)]

  if (nrow(gwas_cross) > 0) {
    gwas_cross[, progression_class := factor(progression_class,
                                             levels = class_order)]

    gwas_pal <- c(
      "nafld_specific"          = "#F48FB1",
      "cirrhosis_specific"      = "#880E4F",
      "enzyme_specific"         = "#42A5F5",
      "shared_cirrhosis_enzyme" = "#7B1FA2",
      "steatosis_specific"      = "#00695C",
      "multi_class"             = "#F57F17",
      "hcc_specific"            = "#4A148C",
      "no_coloc"                = "#E0E0E0"
    )

    p_c_gwas <- ggplot(gwas_cross[!is.na(progression_class)],
                       aes(x = progression_class, y = N, fill = gwas_phenotype_class)) +
      geom_col(position = "stack", width = 0.6) +
      scale_fill_manual(values = gwas_pal, name = "GWAS phenotype") +
      labs(x = "Progression class",
           y = "N genes",
           title = "GWAS phenotype by class") +
      theme_masld() +
      theme(axis.text.x = element_text(angle = 25, hjust = 1),
            legend.position = "right",
            legend.text = element_text(size = 5))
  }
}

if (!is.null(p_c_gwas)) {
  p_c <- p_c_bars + p_c_gwas + plot_layout(widths = c(1, 1.4))
} else {
  p_c <- p_c_bars
}

save_fig(p_c, file.path(PANEL_DIR, "figS_progression_enrichment_panelC.pdf"),
         width = fig_full_width, height = 3.5)
cat("  Saved Panel C\n")


# ==========================================================================
# 5. Panel D: Drug target progression classification
# ==========================================================================
cat("\n--- Panel D: Drug target progression classification ---\n")

plot_drugs <- copy(drug_class)
plot_drugs <- plot_drugs[!is.na(target_gene) & target_gene != ""]

# Build display label: drug name (target gene)
if ("drug" %in% names(plot_drugs)) {
  plot_drugs[, label := fifelse(
    !is.na(drug) & drug != "" & drug != target_gene,
    paste0(drug, " (", target_gene, ")"),
    target_gene
  )]
} else {
  plot_drugs[, label := target_gene]
}

# Ensure tau and coloc columns exist
if (!"tau" %in% names(plot_drugs)) plot_drugs[, tau := 0]
if (!"coloc_best_pp4" %in% names(plot_drugs)) plot_drugs[, coloc_best_pp4 := 0]
if (!"progression_class" %in% names(plot_drugs)) {
  # Use target_class if available (from Script 211 output)
  if ("target_class" %in% names(plot_drugs)) {
    plot_drugs[, progression_class := target_class]
  } else {
    plot_drugs[, progression_class := "not_causal"]
  }
}

# Sort by COLOC PP.H4 descending
plot_drugs <- plot_drugs[order(-coloc_best_pp4)]
plot_drugs[, label := factor(label, levels = rev(unique(label)))]

# Dot plot: y = drug label, x = tau, size = PP.H4, color = class
# Map target_class to our color scheme
drug_class_colors <- c(
  onset_target              = "#42A5F5",
  mid_stage_target          = "#7B1FA2",
  progression_target        = "#880E4F",
  pan_stage_target          = "#00695C",
  enzyme_associated_target  = "#64B5F6",
  not_genetically_causal    = "#E0E0E0",
  low_confidence_causal     = "#F48FB1",
  high_confidence_causal    = "#C2185B",
  # Fallback aliases for non-target genes
  onset                     = "#42A5F5",
  progression               = "#880E4F",
  pan_stage                 = "#00695C",
  not_causal                = "#E0E0E0"
)

p_d <- ggplot(plot_drugs, aes(x = tau, y = label)) +
  geom_point(aes(size = coloc_best_pp4,
                 color = progression_class)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50",
             linewidth = 0.3) +
  scale_size_continuous(range = c(1.5, 5), name = "COLOC PP.H4",
                        limits = c(0, 1)) +
  scale_color_manual(values = drug_class_colors,
                     name = "Progression class", drop = TRUE) +
  labs(x = expression("Stage specificity (" * tau * ")"),
       y = NULL,
       title = "Drug target progression classification") +
  theme_masld() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 5.5))

save_fig(p_d, file.path(PANEL_DIR, "figS_progression_enrichment_panelD.pdf"),
         width = fig_half_width + 1, height = 3.5)
cat("  Saved Panel D\n")


# ==========================================================================
# 6. Combined 4-panel figure
# ==========================================================================
cat("\n--- Assembling combined figure ---\n")

p_combined <- (p_a | p_b) / p_c / p_d +
  plot_layout(heights = c(1, 0.85, 0.9)) +
  plot_annotation(
    title = "Progression-stratified causal architecture",
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 9, face = "bold"))
  )

save_fig_tall(p_combined,
              file.path(PANEL_DIR, "figS_progression_combined.pdf"),
              width = fig_full_width, height = 11)
cat("  Saved combined figure\n")


# ==========================================================================
# 7. Atlas Integration — add progression columns
# ==========================================================================
cat("\n--- Atlas integration ---\n")

atlas_path <- file.path(ME, "multi_evidence_atlas.csv")
if (!file.exists(atlas_path)) {
  cat("  WARNING: multi_evidence_atlas.csv not found; skipping atlas integration.\n")
} else {
  atlas_full <- fread(atlas_path)
  n_before <- ncol(atlas_full)

  # Remove old columns if re-running
  old_cols <- c("progression_coloc_class", "progression_causal_score",
                "stage_specific_coloc_gwas")
  drop <- intersect(old_cols, names(atlas_full))
  if (length(drop) > 0) atlas_full[, (drop) := NULL]

  # Prepare onset_prog for merge: use symbol to match atlas human_symbol
  prog_merge <- onset_prog[, .(
    symbol,
    progression_coloc_class = progression_class,
    progression_causal_score = progression_causal_score
  )]

  # Add gwas_phenotype_class if available
  if ("gwas_phenotype_class" %in% names(onset_prog)) {
    prog_merge <- onset_prog[, .(
      symbol,
      progression_coloc_class = progression_class,
      progression_causal_score = progression_causal_score,
      stage_specific_coloc_gwas = gwas_phenotype_class
    )]
  }

  prog_merge <- prog_merge[!is.na(symbol) & symbol != ""]
  prog_merge <- prog_merge[!duplicated(symbol)]

  # Merge on human_symbol
  atlas_full <- merge(atlas_full, prog_merge,
                      by.x = "human_symbol", by.y = "symbol",
                      all.x = TRUE)

  # If stage_specific_coloc_gwas was not added from onset_prog, try per-GWAS COLOC
  if (!"stage_specific_coloc_gwas" %in% names(atlas_full) && !is.null(coloc_gwas)) {
    coloc_gwas[, gwas_category := fcase(
      grepl("NAFLD|NASH", gwas_name, ignore.case = TRUE), "NAFLD/NASH",
      grepl("Cirrhosis", gwas_name, ignore.case = TRUE), "Cirrhosis",
      grepl("HCC", gwas_name, ignore.case = TRUE), "HCC",
      grepl("ALT|AST|GGT", gwas_name, ignore.case = TRUE), "Liver_enzyme",
      grepl("PDFF", gwas_name, ignore.case = TRUE), "PDFF",
      default = "Other"
    )]
    best_cat <- coloc_gwas[PP.H4.abf > 0.5,
                           .(stage_specific_coloc_gwas = gwas_category[which.max(PP.H4.abf)]),
                           by = gene]
    atlas_full <- merge(atlas_full, best_cat,
                        by.x = "human_symbol", by.y = "gene", all.x = TRUE)
  }

  n_after <- ncol(atlas_full)
  cat(sprintf("  Atlas: %d -> %d columns (+%d)\n", n_before, n_after, n_after - n_before))

  # Coverage report
  n_class <- sum(!is.na(atlas_full$progression_coloc_class) &
                   atlas_full$progression_coloc_class != "not_causal", na.rm = TRUE)
  n_score <- sum(!is.na(atlas_full$progression_causal_score) &
                   atlas_full$progression_causal_score > 0, na.rm = TRUE)
  n_gwas  <- if ("stage_specific_coloc_gwas" %in% names(atlas_full))
               sum(!is.na(atlas_full$stage_specific_coloc_gwas)) else 0
  cat(sprintf("  progression_coloc_class: %d genes classified\n", n_class))
  cat(sprintf("  progression_causal_score > 0: %d genes\n", n_score))
  cat(sprintf("  stage_specific_coloc_gwas: %d genes annotated\n", n_gwas))

  fwrite(atlas_full, atlas_path)
  cat(sprintf("  Updated atlas: %s\n", atlas_path))
}


# ==========================================================================
# Summary
# ==========================================================================
cat("\n=== Script 212 complete ===\n")
cat(sprintf("  Panel A (enrichment bars): %s\n",
            file.path(PANEL_DIR, "figS_progression_enrichment_panelA.pdf")))
cat(sprintf("  Panel B (tau vs PP.H4): %s\n",
            file.path(PANEL_DIR, "figS_progression_enrichment_panelB.pdf")))
cat(sprintf("  Panel C (onset vs progression): %s\n",
            file.path(PANEL_DIR, "figS_progression_enrichment_panelC.pdf")))
cat(sprintf("  Panel D (drug targets): %s\n",
            file.path(PANEL_DIR, "figS_progression_enrichment_panelD.pdf")))
cat(sprintf("  Combined figure: %s\n",
            file.path(PANEL_DIR, "figS_progression_combined.pdf")))
