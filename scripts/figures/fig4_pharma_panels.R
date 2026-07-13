##############################################################################
# Figure 4 panels (a)-(b): Pharmacotranscriptomics
#
#   (a) MOA-level LINCS reversal summary — bar chart of drug mechanism
#       classes ranked by mean composite reversal score
#   (b) Clinical drug validation tile — 10 MASLD therapeutics × evidence
#
# Exports: p_a, p_b
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

dir.create(file.path(FIG4_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIG4_DIR, "panels", "fig4_pharma_panels.pdf")

# Pre-initialize placeholders
p_a <- placeholder("Panel a: LINCS data not found")
p_b <- placeholder("Panel b: Drug validation not found")

# Known MASLD-relevant MOA classes (for highlighting)
MASLD_MOAS <- c("PPAR receptor agonist", "FXR agonist", "HMGCR inhibitor",
                "Insulin sensitizer", "Retinoid receptor agonist",
                "Retinoid receptor ligand", "Fatty acid synthase inhibitor",
                "PPAR receptor antagonist", "FXR antagonist",
                "Lipase inhibitor")

# ==========================================================================
# Panel (a): MOA-level LINCS reversal summary
#   Aggregate compounds by mechanism of action, show mean composite score
#   per MOA class. Highlight MASLD-relevant classes.
# ==========================================================================
cat("Panel a: MOA-level LINCS summary...\n")

lincs <- load_lincs_compounds()

if (!is.null(lincs) && "composite_score" %in% names(lincs)) {

  lk <- copy(lincs)

  # Consolidate moa columns (moa.x -> moa.y -> MOAss)
  if ("moa" %in% names(lk)) {
    lk[, moa_final := as.character(moa)]
  } else {
    lk[, moa_final := as.character(NA)]
    if ("moa.x" %in% names(lk))
      lk[!is.na(moa.x) & nchar(as.character(moa.x)) > 0,
         moa_final := as.character(moa.x)]
    if ("moa.y" %in% names(lk))
      lk[is.na(moa_final) & !is.na(moa.y) & nchar(as.character(moa.y)) > 0,
         moa_final := as.character(moa.y)]
    if ("MOAss" %in% names(lk))
      lk[is.na(moa_final) & !is.na(MOAss) & nchar(as.character(MOAss)) > 0,
         moa_final := as.character(MOAss)]
  }

  # Aggregate by MOA class
  moa_stats <- lk[!is.na(moa_final), .(
    n_compounds  = .N,
    mean_score   = mean(composite_score, na.rm = TRUE),
    max_score    = max(composite_score, na.rm = TRUE),
    median_score = median(composite_score, na.rm = TRUE)
  ), by = moa_final]

  # Require at least 2 compounds per MOA for robustness
  moa_stats <- moa_stats[n_compounds >= 2]
  setorder(moa_stats, -mean_score)
  moa_top <- moa_stats[1:min(12, .N)]

  # Flag MASLD-relevant
  moa_top[, masld_relevant := moa_final %in% MASLD_MOAS]

  # Clean display names: truncate long names
  moa_top[, display := sub("^(.{35}).*$", "\\1...", moa_final)]
  moa_top[nchar(moa_final) <= 35, display := moa_final]
  moa_top[, display := factor(display, levels = rev(display))]

  # Colors: MASLD-relevant = magenta, other = blue
  bar_fill <- c(`TRUE` = masld_colors$up, `FALSE` = "#42A5F5")

  p_a <- ggplot(moa_top, aes(x = mean_score, y = display, fill = masld_relevant)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0("n=", n_compounds)),
              hjust = -0.15, size = GEOM_TEXT_6PT, color = "black") +
    scale_fill_manual(values = bar_fill,
                      labels = c("Other", "MASLD-relevant"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = "Mean composite reversal score",
         y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 6),
          legend.position = c(0.78, 0.18),
          legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
          legend.key.size = unit(0.25, "cm"))

  cat("  ", nrow(moa_top), "MOA classes plotted (n>=2 compounds each)\n")
  message("[caption] Drug mechanism classes reversing MASLD signature")
}

# ==========================================================================
# Panel (b): Clinical drug validation tile
#   10 MASLD therapeutics × evidence layers
#   Fill: Strong (dark magenta), Moderate (bright magenta), Absent (gray)
# ==========================================================================
cat("Panel b: Clinical drug validation tile...\n")

drug_val <- load_drug_validation_table()

if (!is.null(drug_val) && "drug_name" %in% names(drug_val)) {

  evidence_cols <- setdiff(names(drug_val),
                           c("drug_name", "drug_class", "target", "mechanism"))

  if (length(evidence_cols) > 0) {

    val_long <- melt(drug_val, id.vars = "drug_name",
                     measure.vars = evidence_cols,
                     variable.name = "evidence", value.name = "level")
    val_long[, level := factor(level, levels = c("Strong", "Moderate", "Absent"))]

    # Clean evidence labels
    val_long[, evidence_label := gsub("_", " ", evidence)]

    # Order drugs by number of Strong + Moderate hits (best at top)
    drug_strength <- val_long[level %in% c("Strong", "Moderate"), .N, by = drug_name]
    drug_order <- drug_strength[order(-N)]$drug_name
    # Add any drugs not in the ranking (all Absent)
    drug_order <- c(drug_order, setdiff(unique(val_long$drug_name), drug_order))
    val_long[, drug_name := factor(drug_name, levels = rev(drug_order))]

    p_b <- ggplot(val_long, aes(x = evidence_label, y = drug_name, fill = level)) +
      geom_tile(color = "white", linewidth = 0.4) +
      scale_fill_manual(values = drug_evidence_colors,
                        name = "Evidence", na.value = "gray90") +
      labs(x = NULL, y = NULL) +
      theme_masld() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6),
            legend.position = "right",
            legend.key.size = unit(0.3, "cm"))

    cat("  ", uniqueN(val_long$drug_name), "drugs ×",
        uniqueN(val_long$evidence), "evidence layers\n")
    message("[caption] Known MASLD drug target recovery")
  }
}

# ==========================================================================
# Standalone assembly
# ==========================================================================
cat("Composing pharma panels...\n")

fig <- p_a | p_b
fig <- fig +
  plot_annotation(tag_levels = list(c("a", "b"))) &
  theme(plot.tag = element_text(size = 9, face = "plain"))

save_fig(fig, OUT, width = fig_full_width, height = 3.5)
cat("Saved:", OUT, "\n")

# Individual panels for Illustrator
save_fig(p_a, file.path(FIG4_DIR, "panels", "a_moa_lincs.pdf"),
         width = fig_half_width, height = 3.5)
save_fig(p_b, file.path(FIG4_DIR, "panels", "b_drug_validation.pdf"),
         width = fig_half_width, height = 3.5)

cat("Done. Individual panels in:", file.path(FIG4_DIR, "panels"), "\n")
