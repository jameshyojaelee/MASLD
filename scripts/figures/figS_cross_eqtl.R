##############################################################################
# Supplementary Figure: Cross-eQTL Validation of COLOC Results
#
# Two panels (1 row x 2 cols):
#   (a) Broadaway vs GTEx PP.H4 scatter (UKBB ALT)
#       X: Broadaway PP.H4 (N=1,183), Y: GTEx PP.H4 (N=208)
#       Color by validation category, quadrant shading, gene labels
#   (b) Overlap of PP.H4 > 0.5 genes across eQTL sources
#       Broadaway bulk (N=1,183), GTEx bulk (N=208), sc-eQTL (N=312)
#       Euler/Venn for 2 sources, horizontal bar chart for 3
#
# Data:
#   cross_eqtl_comparison.csv  — per-gene broadaway_pp4, gtex_pp4, sceqtl_pp4
#   cross_eqtl_summary.csv     — summary statistics
#   cross_eqtl_validated_genes.csv — cross-validated gene list
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

PANEL_DIR <- file.path(FIGS04_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

CROSS_DIR <- file.path(CAUSAL, "cross_eqtl_comparison")

# Pre-initialize all panels as placeholders
p_a <- placeholder("(a) Broadaway vs GTEx PP.H4 scatter")
p_b <- placeholder("(b) eQTL source overlap")

# ==========================================================================
# Panel (a): Broadaway vs GTEx PP.H4 scatter (UKBB ALT)
# ==========================================================================
cat("Panel (a): Broadaway vs GTEx PP.H4 scatter\n")

comp_file <- file.path(CROSS_DIR, "cross_eqtl_comparison.csv")
if (file.exists(comp_file)) {
  tryCatch({
    comp <- fread(comp_file)
    message("  Loaded cross_eqtl_comparison.csv: ", nrow(comp), " rows, ",
            ncol(comp), " cols")

    # Ensure required columns
    has_broadaway <- "broadaway_pp4" %in% names(comp)
    has_gtex      <- "gtex_pp4" %in% names(comp)

    if (has_broadaway && has_gtex) {
      # Remove rows where both are NA
      comp <- comp[!is.na(broadaway_pp4) | !is.na(gtex_pp4)]
      comp[is.na(broadaway_pp4), broadaway_pp4 := 0]
      comp[is.na(gtex_pp4), gtex_pp4 := 0]

      # Add gene symbols if needed
      gene_col <- intersect(c("gene", "symbol", "human_symbol"), names(comp))[1]
      if (!is.na(gene_col) && gene_col != "symbol") {
        comp <- add_symbols(comp, gene_col)
      }
      # If symbol still missing, create from whatever gene column exists
      if (!"symbol" %in% names(comp)) {
        if ("gene" %in% names(comp)) {
          comp[, symbol := gene]
        } else {
          comp[, symbol := paste0("gene_", .I)]
        }
      }

      # Thresholds: Broadaway PP.H4 > 0.5, GTEx PP.H4 > 0.3 (lower due to smaller N)
      broadaway_thresh <- 0.5
      gtex_thresh      <- 0.3

      # Classify validation categories
      comp[, category := fifelse(
        broadaway_pp4 > broadaway_thresh & gtex_pp4 > gtex_thresh,
        "Cross-validated",
        fifelse(broadaway_pp4 > broadaway_thresh & gtex_pp4 <= gtex_thresh,
                "Broadaway only",
                fifelse(broadaway_pp4 <= broadaway_thresh & gtex_pp4 > gtex_thresh,
                        "GTEx only",
                        "Neither"))
      )]

      n_cross <- sum(comp$category == "Cross-validated")
      n_broadaway_only <- sum(comp$category == "Broadaway only")
      n_gtex_only <- sum(comp$category == "GTEx only")

      # Spearman correlation
      rho <- cor(comp$broadaway_pp4, comp$gtex_pp4, method = "spearman",
                 use = "complete.obs")

      # Category colors
      cat_colors <- c(
        "Cross-validated" = masld_colors$conserved,  # "#00695C" teal
        "Broadaway only"  = masld_colors$up,               # "#C2185B" magenta
        "GTEx only"       = masld_colors$down,             # "#1565C0" blue
        "Neither"         = masld_colors$ns                 # "#BDBDBD" gray
      )

      # Ensure category ordering for legend
      comp[, category := factor(category,
                                levels = c("Cross-validated", "Broadaway only",
                                           "GTEx only", "Neither"))]

      # Select genes to label: top cross-validated + notable genes
      known_genes <- c("THRB", "DGAT2", "HSD17B13", "SLC39A8", "SORT1",
                        "CELSR2", "CDK6", "RORA", "MARC1", "EFHD1",
                        "PNPLA3", "TM6SF2")
      label_dt <- comp[category == "Cross-validated" |
                          (symbol %in% known_genes &
                           (broadaway_pp4 > broadaway_thresh | gtex_pp4 > gtex_thresh))]
      label_dt <- head(label_dt[order(-(broadaway_pp4 + gtex_pp4))], 20)

      p_a <- ggplot(comp, aes(x = broadaway_pp4, y = gtex_pp4)) +
        # Top-right quadrant shading (cross-validated zone)
        annotate("rect",
                 xmin = broadaway_thresh, xmax = 1,
                 ymin = gtex_thresh, ymax = 1,
                 fill = "#E8F5E9", alpha = 0.5) +
        # Diagonal reference line
        geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                    linewidth = 0.2, color = "gray40") +
        # Threshold reference lines
        geom_vline(xintercept = broadaway_thresh, linetype = "dashed",
                   linewidth = 0.2, color = "gray60") +
        geom_hline(yintercept = gtex_thresh, linetype = "dashed",
                   linewidth = 0.2, color = "gray60") +
        # Points (rasterized for PDF performance)
        rasterize_layer(
          geom_point(aes(color = category), size = 0.8, alpha = 0.6)
        ) +
        # Gene labels
        geom_text_repel(data = label_dt, aes(label = symbol),
                        size = 1.8, max.overlaps = 15, segment.size = 0.15,
                        min.segment.length = 0, fontface = "italic",
                        color = "black") +
        # Annotation: rho and cross-validated count
        annotate("text", x = 0.02, y = 0.95,
                 label = paste0("Spearman \u03c1 = ", sprintf("%.3f", rho),
                                "\n", n_cross, " cross-validated"),
                 size = 2.2, hjust = 0, fontface = "bold", color = "gray30") +
        # Threshold annotations
        annotate("text", x = broadaway_thresh + 0.02, y = 0.02,
                 label = "PP.H4 = 0.5", size = 1.8, hjust = 0,
                 color = "gray50", fontface = "italic") +
        annotate("text", x = 0.02, y = gtex_thresh + 0.02,
                 label = "PP.H4 = 0.3", size = 1.8, hjust = 0,
                 color = "gray50", fontface = "italic") +
        scale_color_manual(values = cat_colors, name = NULL, drop = FALSE) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(x = "Broadaway eQTL PP.H4 (N = 1,183)",
             y = "GTEx eQTL PP.H4 (N = 208)",
             title = "Cross-eQTL COLOC validation (UKBB ALT)") +
        theme_masld() +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.25, "cm"))

      message("  Panel (a) built: ", nrow(comp), " genes, ",
              n_cross, " cross-validated, rho = ", sprintf("%.3f", rho))
    } else {
      message("  Missing required columns (broadaway_pp4, gtex_pp4)")
    }
  }, error = function(e) {
    message("  Panel (a) error: ", conditionMessage(e))
    p_a <<- placeholder("(a) Broadaway vs GTEx scatter — error")
  })
} else {
  message("  cross_eqtl_comparison.csv not found: ", comp_file)
}

# ==========================================================================
# Panel (b): Overlap of PP.H4 > 0.5 genes across eQTL sources
# ==========================================================================
cat("Panel (b): eQTL source overlap\n")

if (file.exists(comp_file)) {
  tryCatch({
    comp <- fread(comp_file)

    # Determine which eQTL sources are available
    has_broadaway <- "broadaway_pp4" %in% names(comp)
    has_gtex      <- "gtex_pp4" %in% names(comp)
    has_sceqtl    <- "sceqtl_pp4" %in% names(comp)

    # Add symbols if needed
    gene_col <- intersect(c("gene", "symbol", "human_symbol"), names(comp))[1]
    if (!is.na(gene_col) && gene_col != "symbol") {
      comp <- add_symbols(comp, gene_col)
    }
    if (!"symbol" %in% names(comp)) {
      if ("gene" %in% names(comp)) comp[, symbol := gene]
    }

    # Build gene sets (PP.H4 > 0.5 for each source)
    eqtl_sets <- list()
    if (has_broadaway) {
      eqtl_sets[["Broadaway\n(N=1,183)"]] <-
        comp[broadaway_pp4 > 0.5 & !is.na(broadaway_pp4), symbol]
    }
    if (has_gtex) {
      eqtl_sets[["GTEx\n(N=208)"]] <-
        comp[gtex_pp4 > 0.5 & !is.na(gtex_pp4), symbol]
    }
    if (has_sceqtl) {
      eqtl_sets[["sc-eQTL\n(N=312)"]] <-
        comp[sceqtl_pp4 > 0.5 & !is.na(sceqtl_pp4), symbol]
    }

    n_sources <- length(eqtl_sets)
    message("  Available eQTL sources: ", n_sources)

    if (n_sources >= 2) {
      # Build membership matrix
      all_genes <- unique(unlist(eqtl_sets))
      member_mat <- data.table(gene = all_genes)
      source_names <- names(eqtl_sets)
      for (s in source_names) {
        member_mat[, (s) := gene %in% eqtl_sets[[s]]]
      }

      # Count number of sources per gene
      member_mat[, n_sources := rowSums(.SD), .SDcols = source_names]

      # Create combination label
      member_mat[, combo := apply(.SD, 1, function(x) {
        paste(source_names[x], collapse = " + ")
      }), .SDcols = source_names]

      # Count per combination
      combo_counts <- member_mat[, .N, by = combo]
      setorder(combo_counts, -N)

      # Number of sources per combo (for color)
      combo_counts[, n_src := sapply(as.character(combo), function(x) {
        length(strsplit(x, " \\+ ")[[1]])
      })]

      combo_counts[, combo := factor(combo, levels = rev(combo_counts$combo))]

      # Color by number of sources
      src_fill_colors <- c("1" = "#90CAF9", "2" = "#1565C0", "3" = "#0D47A1")

      p_b <- ggplot(combo_counts, aes(x = N, y = combo,
                                       fill = factor(n_src))) +
        geom_bar(stat = "identity", width = 0.7) +
        geom_text(aes(label = N), hjust = -0.15, size = 2, color = "gray30") +
        scale_fill_manual(values = src_fill_colors, name = "# eQTL\nsources") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
        labs(x = "Number of genes (PP.H4 > 0.5)",
             y = "eQTL source combination",
             title = "COLOC gene overlap across eQTL datasets") +
        theme_masld() +
        theme(axis.text.y = element_text(size = 5),
              legend.position = "inside",
              legend.position.inside = c(0.82, 0.82),
              legend.key.size = unit(0.25, "cm"))

      # Summary annotation
      n_any <- length(all_genes)
      n_multi <- sum(member_mat$n_sources >= 2)
      message("  Panel (b) built: ", n_any, " total genes, ",
              n_multi, " validated in 2+ sources")
    } else if (n_sources == 1) {
      message("  Only 1 eQTL source available — cannot build overlap panel")
      p_b <- placeholder("(b) eQTL overlap — only 1 source available")
    } else {
      message("  No eQTL source columns found in comparison data")
    }
  }, error = function(e) {
    message("  Panel (b) error: ", conditionMessage(e))
    p_b <<- placeholder("(b) eQTL overlap — error")
  })
} else {
  message("  cross_eqtl_comparison.csv not found: ", comp_file)
}

# ==========================================================================
# Assemble composite figure
# ==========================================================================
message("Assembling composite figure...")

fig <- p_a | p_b
fig <- fig +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Save composite
OUT <- file.path(FIGS04_DIR, "figS_cross_eqtl.pdf")
save_fig(fig, OUT, width = fig_full_width, height = 4)
message("Saved composite: ", OUT)

# Save individual panels
save_fig(p_a, file.path(PANEL_DIR, "panel_a_broadaway_vs_gtex_scatter.pdf"),
         width = fig_half_width, height = 4)
save_fig(p_b, file.path(PANEL_DIR, "panel_b_eqtl_source_overlap.pdf"),
         width = fig_half_width, height = 4)

message("=== figS_cross_eqtl.R complete ===")
