##############################################################################
# Supplementary Figure 9: Multi-Ancestry COLOC
#
# Three high-priority supplementary analyses for Fig 3:
#   (a) Cross-ancestry PP.H4 scatter (UKBB ALT vs BBJ ALT)
#   (b) COLOC prior sensitivity heatmap (p12 = 1e-6 to 1e-4)
#   (c) UpSet plot for COLOC GWAS gene overlap (4 European GWAS)
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

PANEL_DIR <- FIGS09_DIR
dir.create(file.path(FIGS09_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# Panel (a): Cross-ancestry PP.H4 scatter — UKBB ALT vs BBJ ALT
# ==========================================================================
cat("Panel (a): Cross-ancestry PP.H4 scatter\n")
p_a <- placeholder("(a) Cross-ancestry scatter")

ukbb_alt_f <- file.path(CAUSAL, "broadaway_ukbb", "coloc_results.csv")
bbj_alt_f  <- file.path(CAUSAL, "bbj_alt", "coloc_results.csv")

if (file.exists(ukbb_alt_f) && file.exists(bbj_alt_f)) {
  ukbb <- fread(ukbb_alt_f)
  bbj  <- fread(bbj_alt_f)

  ukbb <- add_symbols(ukbb, "gene")
  bbj  <- add_symbols(bbj, "gene")

  # Best PP.H4 per gene in each
  ukbb_best <- ukbb[, .(ukbb_pp4 = max(PP.H4, na.rm = TRUE)), by = symbol]
  bbj_best  <- bbj[, .(bbj_pp4 = max(PP.H4, na.rm = TRUE)), by = symbol]

  # Inner join on shared genes
  cross <- merge(ukbb_best, bbj_best, by = "symbol")

  # Load dream for DEG status
  me <- load_multi_evidence()
  if (!is.null(me)) {
    deg_genes <- me[is_dream_deg(me) & abs(bulk_logFC) > 0.3, human_symbol]
    cross[, is_deg := symbol %in% deg_genes]
  } else {
    cross[, is_deg := FALSE]
  }

  # Cross-validated: PP.H4 > 0.5 in BOTH
  cross[, cross_validated := ukbb_pp4 > 0.5 & bbj_pp4 > 0.5]
  n_cross <- sum(cross$cross_validated)

  # Count AFR and CSA exploratory genes from cross-ancestry comparison or COLOC files
  n_afr <- 0L
  n_csa <- 0L
  ca_file <- file.path(CAUSAL, "cross_ancestry", "cross_ancestry_comparison.csv")
  if (file.exists(ca_file)) {
    ca_dt <- tryCatch(fread(ca_file), error = function(e) NULL)
    if (!is.null(ca_dt) && "class" %in% names(ca_dt)) {
      n_afr <- sum(ca_dt$class == "AFR_exploratory", na.rm = TRUE)
      n_csa <- sum(ca_dt$class == "CSA_exploratory", na.rm = TRUE)
      message("  AFR exploratory: ", n_afr, ", CSA exploratory: ", n_csa)
    }
  }
  # Fallback: count from COLOC result files directly
  if (n_afr == 0L) {
    afr_f <- file.path(CAUSAL, "panukbb_afr_alt", "coloc_results.csv")
    if (file.exists(afr_f)) {
      afr_dt <- tryCatch(fread(afr_f), error = function(e) NULL)
      if (!is.null(afr_dt) && "PP.H4" %in% names(afr_dt))
        n_afr <- uniqueN(afr_dt[PP.H4 > 0.5, gene])
    }
  }
  if (n_csa == 0L) {
    csa_f <- file.path(CAUSAL, "panukbb_csa_alt", "coloc_results.csv")
    if (file.exists(csa_f)) {
      csa_dt <- tryCatch(fread(csa_f), error = function(e) NULL)
      if (!is.null(csa_dt) && "PP.H4" %in% names(csa_dt))
        n_csa <- uniqueN(csa_dt[PP.H4 > 0.5, gene])
    }
  }

  subtitle_text <- paste0(n_cross, " EUR+EAS cross-ancestry validated; AFR (",
                           n_afr, ") & CSA (", n_csa, ") exploratory")
  message("[caption] Cross-ancestry COLOC replication: ", subtitle_text)

  # Label top cross-validated genes
  known_genes <- c("THRB", "DGAT2", "HSD17B13", "SLC39A8", "SORT1",
                   "CELSR2", "CDK6", "RORA", "MARC1", "EFHD1")
  label_dt <- cross[cross_validated == TRUE |
                     (symbol %in% known_genes & (ukbb_pp4 > 0.5 | bbj_pp4 > 0.5))]
  label_dt <- head(label_dt[order(-ukbb_pp4 - bbj_pp4)], 20)

  p_a <- ggplot(cross, aes(x = ukbb_pp4, y = bbj_pp4)) +
    # Background quadrant shading
    annotate("rect", xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1,
             fill = "#E8F5E9", alpha = 0.5) +
    rasterize_layer(
      geom_point(aes(color = is_deg, shape = cross_validated),
                 size = 0.8, alpha = 0.5)
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "gray60") +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "gray60") +
    geom_text_repel(data = label_dt, aes(label = symbol),
                    size = GEOM_TEXT_6PT, max.overlaps = 15, segment.size = 0.15,
                    min.segment.length = 0, fontface = "italic",
                    color = "black") +
    scale_color_manual(values = c("TRUE" = masld_colors$up, "FALSE" = "gray60"),
                       labels = c("TRUE" = "DEG", "FALSE" = "Non-DEG"),
                       name = "DEG status") +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                       labels = c("TRUE" = "Cross-validated", "FALSE" = "Single ancestry"),
                       name = NULL) +
    annotate("text", x = 0.75, y = 0.05,
             label = paste0(n_cross, " cross-ancestry\nvalidated genes"),
             size = GEOM_TEXT_6PT, color = "#2E7D32", fontface = "plain") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "UKBB ALT PP.H4 (European)",
         y = "BBJ ALT PP.H4 (East Asian)") +
    theme_masld() +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"))
}

# ==========================================================================
# Panel (b): COLOC prior sensitivity — p12 variation
#   Load pre-computed sensitivity results for AFR and EUR sources
#   Show PP.H4 trajectories across p12 values, colored by fragility
# ==========================================================================
cat("Panel (b): COLOC prior sensitivity\n")
p_b <- placeholder("(b) Prior sensitivity")

sens_file <- file.path(CAUSAL, "prior_sensitivity", "prior_sensitivity_results.csv")
if (file.exists(sens_file)) {
  tryCatch({
    sens_all <- fread(sens_file)
    message("  Loaded prior sensitivity: ", nrow(sens_all), " rows, ",
            uniqueN(sens_all$gene), " genes, ",
            uniqueN(sens_all$source), " sources")

    # Gene column already has symbols — no add_symbols() needed
    # Determine fragility per gene per source: robust if PP.H4 > 0.5 at ALL p12
    fragility <- sens_all[, .(
      min_pp4 = min(PP.H4, na.rm = TRUE),
      max_pp4 = max(PP.H4, na.rm = TRUE),
      default_pp4 = PP.H4[which.min(abs(p12 - 5e-6))]
    ), by = .(source, source_N, gene)]
    fragility[, fragility := fifelse(
      min_pp4 > 0.5, "Robust",
      fifelse(max_pp4 > 0.5, "Fragile", "Below")
    )]

    # Build facet labels with fragility counts
    facet_info <- fragility[, .(
      n_total = .N,
      n_fragile = sum(fragility == "Fragile"),
      n_robust = sum(fragility == "Robust"),
      N = source_N[1]
    ), by = source]
    facet_info[, source_label := fifelse(
      grepl("AFR", source),
      paste0("AFR (N=", formatC(N, format = "d", big.mark = ","), "): ",
             n_fragile, "/", n_total, " fragile"),
      paste0("EUR (N=", formatC(N, format = "d", big.mark = ","), "): ",
             n_fragile, "/", n_total, " fragile")
    )]

    # Merge fragility + facet labels back
    sens_all <- merge(sens_all, fragility[, .(source, gene, fragility)],
                      by = c("source", "gene"))
    sens_all <- merge(sens_all, facet_info[, .(source, source_label)],
                      by = "source")

    # Filter genes to show: AFR = all genes; EUR = top 30 by default PP.H4
    afr_genes <- fragility[grepl("AFR", source)]$gene
    eur_top <- fragility[grepl("UKBB", source)]
    setorder(eur_top, -default_pp4)
    eur_genes <- head(eur_top, 30)$gene

    sens_plot <- sens_all[
      (grepl("AFR", source) & gene %in% afr_genes) |
      (grepl("UKBB", source) & gene %in% eur_genes)
    ]

    frag_colors <- c("Robust" = "#00695C", "Fragile" = "#C2185B", "Below" = "#BDBDBD")

    p_b <- ggplot(sens_plot, aes(x = p12, y = PP.H4, group = gene,
                                  color = fragility)) +
      geom_line(linewidth = 0.3, alpha = 0.6) +
      geom_point(size = 0.5, shape = 16, alpha = 0.6) +
      geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.2,
                 color = "gray60") +
      scale_x_log10(labels = scales::scientific) +
      scale_color_manual(values = frag_colors, name = "Fragility") +
      facet_wrap(~source_label, scales = "free_y") +
      labs(x = expression(p[12]~"prior"),
           y = "PP.H4") +
      theme_masld() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.25, "cm"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            strip.text = element_text(size = 6))

    message("  Panel (b) built: ", nrow(sens_plot), " data points")
    message("[caption] COLOC prior sensitivity (p12 variation)")
  }, error = function(e) {
    message("  Panel (b) error: ", conditionMessage(e))
    p_b <<- placeholder("(b) Prior sensitivity — error")
  })
} else {
  message("  Prior sensitivity file not found: ", sens_file)
}

# ==========================================================================
# Panel (c): UpSet plot — COLOC gene overlap across 4 European GWAS
# ==========================================================================
cat("Panel (c): COLOC UpSet plot\n")
p_c <- placeholder("(c) UpSet overlap")

broadaway_dirs_c <- c(
  "UKBB ALT" = "broadaway_ukbb",
  "UKBB AST" = "broadaway_ukbb_ast",
  "UKBB GGT" = "broadaway_ukbb_ggt",
  "PDFF"     = "broadaway_pdff"
)

gwas_gene_sets <- list()
for (gwas_label in names(broadaway_dirs_c)) {
  f <- file.path(CAUSAL, broadaway_dirs_c[[gwas_label]], "coloc_results.csv")
  if (file.exists(f)) {
    tmp <- fread(f)
    if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
      tmp <- add_symbols(tmp, "gene")
      gwas_gene_sets[[gwas_label]] <- unique(tmp[PP.H4 > 0.5, symbol])
    }
  }
}

if (length(gwas_gene_sets) >= 2) {
  # Build membership matrix
  all_genes <- unique(unlist(gwas_gene_sets))
  member_mat <- data.table(gene = all_genes)
  for (g in names(gwas_gene_sets)) {
    member_mat[, (g) := gene %in% gwas_gene_sets[[g]]]
  }

  # Compute intersection sizes
  gwas_names <- names(gwas_gene_sets)
  # Create intersection key
  member_mat[, int_key := do.call(paste, c(lapply(gwas_names, function(g) {
    ifelse(get(g), g, "")
  }), sep = "|"))]
  member_mat[, int_key := gsub("\\|+", " & ", int_key)]
  member_mat[, int_key := gsub("^ & | & $", "", int_key)]
  member_mat[, int_key := trimws(int_key)]

  # Simplified: count per combination
  member_mat[, n_gwas := rowSums(.SD), .SDcols = gwas_names]

  # Create combo column
  member_mat[, combo := apply(.SD, 1, function(x) {
    paste(gwas_names[x], collapse = "\n")
  }), .SDcols = gwas_names]

  combo_counts <- member_mat[, .N, by = combo]
  setorder(combo_counts, -N)

  # Use a horizontal bar chart as a simpler UpSet-style visualization
  # (full UpSet requires UpSetR which may not be installed)
  combo_counts[, combo := factor(combo, levels = rev(combo_counts$combo))]

  # Color by number of GWAS in the intersection
  combo_counts[, n_gwas_combo := sapply(as.character(combo), function(x) {
    length(strsplit(x, "\n")[[1]])
  })]

  p_c <- ggplot(combo_counts, aes(x = N, y = combo, fill = factor(n_gwas_combo))) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = N), hjust = -0.15, size = GEOM_TEXT_6PT, color = "black") +
    scale_fill_manual(values = c("1" = "#90CAF9", "2" = "#42A5F5",
                                  "3" = "#1565C0", "4" = "#0D47A1"),
                      name = "# GWAS") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "Number of genes (PP.H4 > 0.5)",
         y = "GWAS combination") +
    theme_masld() +
    theme(axis.text.y = element_text(size = 6),
          legend.position = "inside",
          legend.position.inside = c(0.8, 0.8),
          legend.key.size = unit(0.25, "cm"))
}
message("[caption] COLOC gene overlap across 4 European GWAS")

# ==========================================================================
# Assemble
# ==========================================================================
fig_supp <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

OUT <- file.path(PANEL_DIR, "figS09_multi_ancestry_coloc.pdf")
save_fig_tall(fig_supp, OUT, height = 9)

# Individual panels
save_fig(p_a, file.path(PANEL_DIR, "panels", "figSa_cross_ancestry_scatter.pdf"), height = 4)
save_fig(p_b, file.path(PANEL_DIR, "panels", "figSb_prior_sensitivity.pdf"), height = 4)
save_fig(p_c, file.path(PANEL_DIR, "panels", "figSc_upset_overlap.pdf"), height = 4)

message("Supplementary causal panels saved to ", file.path(PANEL_DIR, "panels"))
