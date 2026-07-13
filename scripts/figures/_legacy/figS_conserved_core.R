# Supplementary Figure: Conserved Enrichment Analysis — 6 panels
#   (a) Forest plot of ORs (fixed)
#   (b) -log10(p) significance bars
#   (c) 4×5 concordance heatmap (human sig × mouse diet)
#   (d) n_concordant histogram (Conserved threshold at ≥3)
#   (e) Cross-anchor tier breakdown bar chart
#   (f) Discordance barcode: top 50 Species_Discordant genes

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

OUT        <- file.path(FIGS06_DIR, "figS_conserved.pdf")
dir.create(file.path(FIGS06_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
CONC_DIR   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
matrix_path  <- file.path(CONC_DIR, "gene_concordance_matrix_20x.csv")
pergene_path <- file.path(CONC_DIR, "gene_concordance_per_gene.csv")
allsig_path  <- file.path(CONC_DIR, "gene_concordance_all_signatures.csv")
atlas_path   <- file.path(CONC_DIR, "concordance_atlas_unified.csv")

# Pre-initialize placeholders
p_a <- placeholder("Panel a: enrichment data not found")
p_b <- placeholder("Panel b: enrichment data not found")
p_c <- placeholder("Panel c: concordance matrix not found")
p_d <- placeholder("Panel d: per-gene concordance not found")
p_e <- placeholder("Panel e: atlas not found")
p_f <- placeholder("Panel f: per-gene concordance not found")

# ==========================================================================
# Panels (a) and (b): Enrichment forest plot + significance bars
# ==========================================================================
enr <- load_conserved_enr()

if (!is.null(enr) && nrow(enr) > 0) {
  enr[OR == 0 | is.na(OR), OR := 0.01]
  enr[, neg_log10p := -log10(p_value)]
  enr[, sig := p_value < 0.05]
  enr[, test_label := gsub("_", " ", test)]
  enr[, test_label := factor(test_label, levels = rev(test_label[order(OR)]))]

  # CI: clip p_value to avoid Inf from qnorm
  enr[, log_or := log(OR)]
  enr[, p_clip := pmax(p_value, 1e-300)]
  enr[, se_log_or := abs(log_or / qnorm(p_clip / 2))]
  enr[is.infinite(se_log_or) | is.nan(se_log_or) | se_log_or > 5, se_log_or := 0.5]
  enr[, ci_lo := exp(log_or - 1.96 * se_log_or)]
  enr[, ci_hi := pmin(exp(log_or + 1.96 * se_log_or), 200)]

  p_a <- ggplot(enr, aes(x = OR, y = test_label, color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "gray50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.4) +
    geom_point(size = 2.5, shape = 16) +
    geom_text(aes(label = sprintf("OR=%.1f", OR)), hjust = -0.15, size = 2, color = "black") +
    scale_color_manual(values = c(`TRUE` = masld_colors$conserved, `FALSE` = masld_colors$ns),
                       labels = c("NS", "p < 0.05"), name = NULL) +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.4)),
                  breaks = c(1, 2, 5, 10, 20, 50)) +
    labs(x = "Odds Ratio (log scale)", y = NULL, title = "Conserved enrichment") +
    theme_masld() +
    theme(legend.position = "inside", legend.position.inside = c(0.85, 0.15),
          legend.background = element_blank())

  p_b <- ggplot(enr, aes(x = neg_log10p, y = test_label, fill = sig)) +
    geom_col(width = 0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.3, color = "gray50") +
    geom_text(aes(label = sprintf("p=%.1e", p_value)), hjust = -0.1, size = 2) +
    scale_fill_manual(values = c(`TRUE` = masld_colors$conserved,
                                  `FALSE` = masld_colors$ns), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
    labs(x = expression(-log[10]~p), y = NULL, title = "Significance") +
    theme_masld()
}

# ==========================================================================
# Panel (c): 4×5 concordance rate heatmap (human_signature × diet)
# ==========================================================================
if (file.exists(matrix_path)) {
  mat <- fread(matrix_path)

  # Clean labels
  mat[, sig_label := fcase(
    human_signature == "nafl_vs_nash",    "NAFL→NASH",
    human_signature == "disease_vs_ctrl", "Disease vs Ctrl",
    human_signature == "nafl_specific",   "NAFL vs Ctrl",
    human_signature == "fibrosis",        "Fibrosis",
    default = human_signature
  )]
  mat[, sig_label := factor(sig_label,
    levels = c("NAFL→NASH", "Disease vs Ctrl", "NAFL vs Ctrl", "Fibrosis"))]

  p_c <- ggplot(mat, aes(x = diet, y = sig_label, fill = concordance_pct)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", concordance_pct)), size = 2.2) +
    scale_fill_gradient2(
      low = masld_colors$down, mid = "white", high = masld_colors$conserved,
      midpoint = 50, limits = c(0, 100),
      name = "Concordance\n(%)"
    ) +
    labs(x = "Mouse diet model", y = "Human signature",
         title = "Cross-species concordance (% concordant DEGs)") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ==========================================================================
# Panel (d): Histogram of n_concordant diets per gene (NAFL-vs-NASH anchor)
# ==========================================================================
if (file.exists(pergene_path)) {
  pg <- fread(pergene_path)

  # Count genes per n_concordant value
  hist_dt <- pg[, .N, by = n_concordant]
  setorder(hist_dt, n_concordant)
  hist_dt[, threshold := n_concordant >= 3]

  p_d <- ggplot(pg, aes(x = n_concordant)) +
    geom_histogram(aes(fill = n_concordant >= 3), binwidth = 1,
                   color = "white", linewidth = 0.3) +
    geom_vline(xintercept = 2.5, linetype = "dashed",
               color = masld_colors$up, linewidth = 0.6) +
    annotate("text", x = 3.1, y = Inf, label = "Conserved\nthreshold (≥3)",
             vjust = 1.5, hjust = 0, size = 2, color = masld_colors$up) +
    scale_fill_manual(
      values = c(`FALSE` = masld_colors$not_sig, `TRUE` = masld_colors$conserved),
      labels = c("< 3 diets", "≥ 3 diets (Conserved)"),
      name = NULL
    ) +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(labels = comma) +
    labs(x = "Number of concordant mouse diets",
         y = "Number of genes",
         title = "Concordance breadth per gene (NAFL→NASH anchor)") +
    theme_masld() +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.85),
          legend.background = element_blank(), legend.key = element_blank())
}

# ==========================================================================
# Panel (e): Cross-anchor tier breakdown (from unified atlas)
# ==========================================================================
if (file.exists(atlas_path)) {
  atlas <- fread(atlas_path)

  tier_order <- c("Dual_Conserved", "NASH_Conserved", "MASLD_Conserved", "Other")
  tier_labels <- c(
    Dual_Conserved  = "Dual Conserved\n(NAFL→NASH & Disease-vs-Ctrl)",
    NASH_Conserved  = "NASH Conserved\n(NAFL→NASH only)",
    MASLD_Conserved = "MASLD Conserved\n(Disease-vs-Ctrl only)",
    Other           = "Other"
  )
  tier_colors <- c(
    Dual_Conserved  = masld_colors$conserved,
    NASH_Conserved  = "#42A5F5",
    MASLD_Conserved = "#66BB6A",
    Other           = masld_colors$not_sig
  )

  # Only look at Conserved genes
  core <- atlas[primary_category == "Conserved"]
  core[!cross_anchor_tier %in% tier_order, cross_anchor_tier := "Other"]
  tier_counts <- core[, .N, by = cross_anchor_tier]
  tier_counts[, cross_anchor_tier := factor(cross_anchor_tier, levels = rev(tier_order))]
  tier_counts[, label := paste0(N, "\n(", round(N / sum(N) * 100, 1), "%)")]

  p_e <- ggplot(tier_counts, aes(x = N, y = cross_anchor_tier, fill = cross_anchor_tier)) +
    geom_col(width = 0.6, color = NA) +
    geom_text(aes(label = label), hjust = -0.05, size = 2) +
    scale_fill_manual(values = tier_colors, guide = "none",
                      labels = tier_labels) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(tier_labels[x], width = 25)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35)), labels = comma) +
    labs(x = "Number of Conserved genes", y = NULL,
         title = "Cross-anchor breakdown of Conserved") +
    theme_masld()
}

# ==========================================================================
# Panel (f): Discordance barcode — top 50 Species_Discordant genes
# ==========================================================================
if (file.exists(allsig_path) && file.exists(pergene_path)) {
  allsig <- fread(allsig_path)
  pg     <- fread(pergene_path)

  # Top 30 genes by total discordance, using the primary (NAFL-vs-NASH) anchor
  # Rank by sum of n_discordant across all signatures
  disc_rank <- allsig[, .(total_discordant = sum(n_discordant),
                           total_concordant = sum(n_concordant)), by = human_symbol]
  disc_rank <- disc_rank[total_discordant > 0][order(-total_discordant)]
  top30 <- disc_rank[1:min(30, .N), human_symbol]

  # We need raw LFC/padj per dataset for true Up/Down/NS status.
  # allsig only has summaries. Let's load the raw human and mouse results.
  human_files <- list(
    `Human (NAFL→NASH)`    = file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nafl_vs_nash_dream.csv"),
    `Human (Disease vs Ctrl)` = file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
    `Human (NAFL vs Ctrl)` = file.path(CONC_DIR, "nafl_vs_ctrl_dream.csv"),
    `Human (Fibrosis)`     = file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_dream.csv")
  )

  mouse_diets <- c("MCD", "HFD", "CDAHFD", "FPC", "LIDPAD")
  mouse_dir <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")

  # Load all and extract the top30
  all_rows <- list()

  # Human
  for (ds_name in names(human_files)) {
    if (file.exists(human_files[[ds_name]])) {
      dt <- fread(human_files[[ds_name]])
      if (!"padj" %in% names(dt) && "adj.P.Val" %in% names(dt)) setnames(dt, "adj.P.Val", "padj")
      
      # Clean existing symbol columns to avoid .x / .y collisions
      if ("symbol" %in% names(dt)) dt[, symbol := NULL]
      
      dt <- add_symbols(dt, "gene")
      sub <- dt[symbol %in% top30, .(human_symbol = symbol, dataset = ds_name, lfc = logFC, padj = padj)]
      all_rows[[ds_name]] <- sub
    }
  }

  # Mouse
  for (diet in mouse_diets) {
    f <- file.path(mouse_dir, paste0(diet, "_de_results.csv"))
    if (file.exists(f)) {
      dt <- fread(f)
      if (!"padj" %in% names(dt) && "adj.P.Val" %in% names(dt)) setnames(dt, "adj.P.Val", "padj")
      dt[, ensembl_clean := sub("\\..*", "", gene)]
      # Map to human symbol
      ortho <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/human_mouse_ortholog_comparison.csv"))
      ortho <- unique(ortho[!is.na(human_symbol) & human_symbol != "", .(mouse_gene_id = sub("\\..*", "", mouse_gene_id), human_symbol)])
      dt <- merge(dt, ortho, by.x = "ensembl_clean", by.y = "mouse_gene_id")
      sub <- dt[human_symbol %in% top30, .(human_symbol, dataset = diet, lfc = logFC, padj = padj)]
      all_rows[[diet]] <- sub
    }
  }

  rows <- rbindlist(all_rows, fill = TRUE)

  # Classify regulation state
  rows[, reg_state := fcase(
    !is.na(padj) & padj < 0.05 & lfc > 0, "Upregulated",
    !is.na(padj) & padj < 0.05 & lfc < 0, "Downregulated",
    default = "Not Significant"
  )]

  # Order genes by most discordant overall
  rows[, human_symbol := factor(human_symbol, levels = top30)]  # plot left-to-right

  # Order datasets
  ds_order <- c(
    "Human (NAFL→NASH)", "Human (Disease vs Ctrl)", "Human (NAFL vs Ctrl)", "Human (Fibrosis)",
    "LIDPAD", "FPC", "CDAHFD", "HFD", "MCD"
  )
  rows[, dataset := factor(dataset, levels = rev(ds_order))]

  reg_colors <- c(
    Upregulated    = "#D81B60",   # magenta
    Downregulated  = "#00897B",   # teal
    `Not Significant` = "#E0E0E0" # light gray
  )

  p_f <- ggplot(rows, aes(x = human_symbol, y = dataset, fill = reg_state)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = reg_colors, name = NULL) +
    labs(x = "Symbol", y = NULL,
         title = "Discordance Barcode: Top 30 Species_Discordant Genes") +
    coord_fixed(ratio = 1) +
    theme_masld() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = 5.5, face = "italic"),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      panel.grid = element_blank()
    )
}

# ==========================================================================
# Assemble: 4-row layout
# Row 1: (a) Forest  | (b) -log10(p) bars
# Row 2: (c) Concordance heatmap | (d) n_concordant histogram
# Row 3: (e) Cross-anchor tier breakdown (full width)
# Row 4: (f) Discordance barcode (full width)
# ==========================================================================
figS4 <- (p_a | p_b) /
          (p_c | p_d) /
          p_e /
          p_f +
  plot_layout(heights = c(0.8, 1, 0.7, 0.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(figS4, OUT, height = 12.5)
message("FigS4 saved to ", OUT)
