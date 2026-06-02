# DEPRECATED — output superseded by renumbered figures
##############################################################################
# Figure 3: Causal Architecture (9 panels)
#   Row 1: (a) COLOC Manhattan (multi-GWAS PP.H4 vs chromosome)
#           (b) Top COLOC genes evidence heatmap (PP.H4 across GWAS)
#   Row 2: (c) sc-eQTL COLOC x cell-type dot plot
#           (d) Causal method coverage bar chart
#   Row 3: (e) ieQTL x DEG directional concordance scatter
#           (f) Multi-ancestry COLOC summary
#   Row 4: (g) TWAS cross-GWAS concordance (Ghodsian vs Chen)
#           (h) [retired — MR forest panel ditched 2026-04-22]
#           (i) sc-TWAS cell-type dot plot
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

OUT <- file.path(FIG3_DIR, "fig3_compact.pdf")
dir.create(file.path(FIG3_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# Pre-initialize all panels with placeholders
p_a <- placeholder("(a) COLOC Manhattan")
p_b <- placeholder("(b) Top COLOC genes heatmap")
p_c <- placeholder("(c) sc-eQTL COLOC x cell type")
p_d <- placeholder("(d) Causal method coverage")
p_e <- placeholder("(e) ieQTL x DEG concordance")
p_f <- placeholder("(f) Multi-ancestry COLOC")
p_g <- placeholder("(g) TWAS Ghodsian vs Chen")
p_h <- placeholder("(h) [retired — MR ditched 2026-04-22]")
p_i <- placeholder("(i) sc-TWAS cell-type dot plot")

# Cell type styling (consistent across panels c and e)
celltype_colors <- c(
  Hepatocyte     = "#0D47A1",
  Endothelial    = "#80CBC4",
  Cholangiocyte  = "#CE93D8",
  Stellate       = "#F8BBD0"
)

celltype_clean_map <- c(
  hepatocyte       = "Hepatocyte",
  endothelial_cell = "Endothelial",
  cholangiocyte    = "Cholangiocyte",
  stellate_cell    = "Stellate"
)

# GWAS-specific colors (consistent across panels)
gwas_colors <- c(
  "ALT" = "#1565C0",
  "AST" = "#D84315",
  "GGT" = "#6A1B9A",
  "PDFF" = "#2E7D32"
)

# ==========================================================================
# Panel (a): COLOC Manhattan — PP.H4 vs chromosome across 4 GWAS
#   Alternating chromosome shading, point size ~ #GWAS, labels prioritized
# ==========================================================================
broadaway_dirs <- c(
  "ALT" = "broadaway_ukbb",
  "AST" = "broadaway_ukbb_ast",
  "GGT" = "broadaway_ukbb_ggt",
  "PDFF" = "broadaway_pdff"
)

coloc_parts_a <- list()
for (gwas_label in names(broadaway_dirs)) {
  f <- file.path(CAUSAL, broadaway_dirs[[gwas_label]], "coloc_results.csv")
  if (file.exists(f)) {
    tmp <- fread(f)
    if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
      tmp[, gwas := gwas_label]
      tmp <- add_symbols(tmp, "gene")
      coloc_parts_a[[gwas_label]] <- tmp
    }
  }
}

if (length(coloc_parts_a) > 0) {
  coloc_all_a <- rbindlist(coloc_parts_a, fill = TRUE)

  # Parse chromosome
  coloc_all_a[, chr_num := suppressWarnings(as.integer(gsub("chr", "", chr)))]
  coloc_all_a <- coloc_all_a[!is.na(chr_num) & chr_num >= 1 & chr_num <= 22]

  # Per-gene #GWAS with PP.H4 > 0.5 (for size encoding)
  gene_n_gwas <- coloc_all_a[PP.H4 > 0.5, .(n_gwas_sig = uniqueN(gwas)), by = symbol]
  coloc_all_a <- merge(coloc_all_a, gene_n_gwas, by = "symbol", all.x = TRUE)
  coloc_all_a[is.na(n_gwas_sig), n_gwas_sig := 0]

  # Keep genes with PP.H4 > 0.3 for plotting (show distribution); label at > 0.5
  coloc_plot_a <- coloc_all_a[PP.H4 > 0.3]

  if (nrow(coloc_plot_a) > 0) {
    # Jitter within chromosome
    coloc_plot_a[, x_pos := chr_num + runif(.N, -0.3, 0.3)]

    # Alternating chromosome background shading rectangles
    chr_shading <- data.table(
      chr = 1:22,
      xmin = (1:22) - 0.45,
      xmax = (1:22) + 0.45,
      ymin = 0.3, ymax = 1.02,
      fill = rep(c("even", "odd"), 11)
    )

    # Genes to label — prioritize 3/3 liver enzyme + known targets
    known_masld_genes <- c("THRB", "DGAT2", "HSD17B13", "SLC39A8", "SORT1",
                           "CELSR2", "CDK6", "RORA", "MARC1", "PNPLA3",
                           "TM6SF2", "MBOAT7", "GCKR", "SAMM50", "EFHD1",
                           "FABP1", "HKDC1", "SPTLC3")
    coloc_best_per_gene <- coloc_plot_a[, .SD[which.max(PP.H4)], by = symbol]

    # Label: only well-known MASLD/liver genes with strong colocalization
    label_dt_a <- coloc_best_per_gene[
      symbol %in% known_masld_genes & PP.H4 > 0.5
    ]

    # Summary counts
    n_total_genes <- uniqueN(coloc_all_a[PP.H4 > 0.5, symbol])
    n_high_genes  <- uniqueN(coloc_all_a[PP.H4 > 0.8, symbol])
    n_replicated_3 <- gene_n_gwas[n_gwas_sig >= 3, .N]

    # Clamp point size: use n_gwas_sig for genes above threshold, else fixed small
    coloc_plot_a[, pt_size := fifelse(n_gwas_sig >= 1, 0.5 + n_gwas_sig * 0.3, 0.5)]

    p_a <- ggplot(coloc_plot_a, aes(x = x_pos, y = PP.H4, color = gwas)) +
      # Alternating chromosome shading
      geom_rect(data = chr_shading[fill == "even"],
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "gray95", alpha = 0.5) +
      geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.2,
                 color = "gray50") +
      geom_hline(yintercept = 0.5, linetype = "dotted", linewidth = 0.2,
                 color = "gray70") +
      rasterize_layer(
        geom_point(aes(size = pt_size), alpha = 0.6, shape = 16)
      ) +
      geom_label_repel(
        data = label_dt_a,
        aes(label = symbol),
        size = 1.7, max.overlaps = 15,
        label.padding = 0.1, segment.size = 0.15,
        min.segment.length = 0, fontface = "italic",
        color = "black", show.legend = FALSE
      ) +
      scale_color_manual(values = gwas_colors, name = "GWAS") +
      scale_size_identity() +
      scale_x_continuous(breaks = 1:22, labels = 1:22,
                         limits = c(0.5, 22.5),
                         expand = expansion(mult = 0.01)) +
      scale_y_continuous(limits = c(0.3, 1.02), expand = c(0, 0)) +
      labs(x = "Chromosome",
           y = "PP.H4 (colocalization probability)",
           title = "Broadaway eQTL \u00d7 UKBB GWAS colocalization",
           subtitle = paste0(n_total_genes, " genes PP.H4 > 0.5; ",
                             n_high_genes, " PP.H4 > 0.8; ",
                             n_replicated_3, " in 3+ GWAS")) +
      theme_masld() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.2, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            plot.subtitle = element_text(size = 5, color = "gray40"))
  }
}

# ==========================================================================
# Panel (b): Top COLOC genes evidence heatmap — PP.H4 across GWAS
#   20 genes, no numeric labels, side annotation bar for # GWAS
# ==========================================================================
if (length(coloc_parts_a) > 0) {
  coloc_all_b <- rbindlist(coloc_parts_a, fill = TRUE)

  # Best PP.H4 per gene x GWAS
  coloc_gene_gwas <- coloc_all_b[, .(PP.H4 = max(PP.H4, na.rm = TRUE)),
                                   by = .(symbol, gwas)]

  # Filter to genes with PP.H4 > 0.5 in at least one GWAS
  gene_max_b <- coloc_gene_gwas[, .(max_pp4 = max(PP.H4),
                                     n_gwas_coloc = sum(PP.H4 > 0.5)),
                                  by = symbol]
  gene_max_b <- gene_max_b[max_pp4 > 0.5]
  gene_max_b <- gene_max_b[!grepl("^AC[0-9]|^AL[0-9]|^RP[0-9]|^LOC[0-9]|^LINC|-DT$", symbol)]
  setorder(gene_max_b, -n_gwas_coloc, -max_pp4)

  # Top 20 genes (reduced from 30 for readability at print size)
  top_b_genes <- head(gene_max_b, 20)$symbol
  gene_order_b <- rev(top_b_genes)

  # Full grid: gene x GWAS
  gwas_levels_b <- c("ALT", "AST", "GGT", "PDFF")
  heat_grid_b <- CJ(symbol = top_b_genes,
                     gwas = gwas_levels_b)
  heat_grid_b <- merge(heat_grid_b,
                        coloc_gene_gwas[, .(symbol, gwas, PP.H4)],
                        by = c("symbol", "gwas"), all.x = TRUE)
  heat_grid_b[is.na(PP.H4), PP.H4 := 0]
  heat_grid_b[, symbol := factor(symbol, levels = gene_order_b)]
  heat_grid_b[, gwas := factor(gwas, levels = gwas_levels_b)]

  # Highlight known MASLD genes in bold
  known_b <- c("HSD17B13", "FABP1", "RORA", "HKDC1", "SPTLC3",
               "EFHD1", "PNPLA3", "TM6SF2", "THRB", "PPARG", "MBOAT7",
               "DGAT2", "SLC39A8", "SORT1", "CELSR2", "CDK6", "MARC1",
               "GCKR", "SAMM50", "CHEK2")
  y_faces_b <- ifelse(gene_order_b %in% known_b, "bold.italic", "italic")
  names(y_faces_b) <- gene_order_b

  # N genes replicated across GWAS
  n_replicated <- gene_max_b[n_gwas_coloc >= 2, .N]

  # Side annotation bar: # GWAS per gene
  n_gwas_side_b <- gene_max_b[symbol %in% top_b_genes, .(symbol, n_gwas_coloc)]
  n_gwas_side_b[, symbol := factor(symbol, levels = gene_order_b)]

  # Main heatmap (no numeric text labels — rely on color gradient)
  p_b_main <- ggplot(heat_grid_b, aes(x = gwas, y = symbol)) +
    geom_tile(aes(fill = PP.H4), color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#9C27B0",
                        limits = c(0, 1), name = "PP.H4") +
    labs(x = NULL, y = NULL,
         title = paste0("COLOC across GWAS (",
                        n_replicated, " genes in 2+ GWAS)")) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5, face = y_faces_b),
          plot.title = element_text(size = 7))

  # Side annotation strip: # GWAS with PP.H4 > 0.5
  p_b_strip <- ggplot(n_gwas_side_b, aes(x = "#GWAS", y = symbol,
                                           fill = n_gwas_coloc)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = n_gwas_coloc), size = 1.8, color = "white",
              fontface = "bold") +
    scale_fill_gradient(low = "#CE93D8", high = "#4A148C",
                        limits = c(1, 4), name = "#GWAS\nsig") +
    labs(x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(size = 4.5, angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(2, 2, 2, 1),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5))

  p_b <- wrap_elements(full =
    p_b_main + p_b_strip +
    plot_layout(widths = c(5, 1))
  )
}

# ==========================================================================
# Panel (c): sc-eQTL COLOC x cell-type — DOT PLOT
#   Size = PP.H4, color = cell type, 15 genes
# ==========================================================================
ct_coloc <- data.table()

sceqtl_coloc <- load_sceqtl_coloc()
if (!is.null(sceqtl_coloc) && nrow(sceqtl_coloc) > 0) {
  pp4_col <- intersect(c("PP.H4_all", "PP.H4"), names(sceqtl_coloc))[1]
  if (!is.na(pp4_col)) {
    ct_coloc <- rbind(ct_coloc,
                      sceqtl_coloc[, .(gene, cell_type, PP.H4 = get(pp4_col),
                                       gwas_source = "Ghodsian")],
                      fill = TRUE)
  }
}

ukbb_coloc_c <- load_sceqtl_coloc_ukbb()
if (!is.null(ukbb_coloc_c) && nrow(ukbb_coloc_c) > 0) {
  pp4_col_uk <- intersect(c("PP.H4_all", "PP.H4"), names(ukbb_coloc_c))[1]
  if (!is.na(pp4_col_uk)) {
    ct_coloc <- rbind(ct_coloc,
                      ukbb_coloc_c[, .(gene, cell_type, PP.H4 = get(pp4_col_uk),
                                       gwas_source = "UKBB ALT")],  # sc-eQTL panel uses separate labels
                      fill = TRUE)
  }
}

if (nrow(ct_coloc) > 0) {
  ct_best_c <- ct_coloc[, .(PP.H4 = max(PP.H4, na.rm = TRUE),
                             best_source = gwas_source[which.max(PP.H4)]),
                          by = .(gene, cell_type)]
  gene_max_c <- ct_best_c[, .(max_pp4 = max(PP.H4),
                               n_cell_types = sum(PP.H4 > 0.3)),
                             by = gene]
  setorder(gene_max_c, -max_pp4)
  # Top 15 genes (reduced from 20 for clarity)
  top_c_genes <- head(gene_max_c[max_pp4 > 0.5], 15)$gene
  gene_order_c <- rev(gene_max_c[gene %in% top_c_genes, gene])

  cell_types_ordered <- c("hepatocyte", "endothelial_cell",
                          "cholangiocyte", "stellate_cell")
  dot_grid_c <- CJ(gene = top_c_genes, cell_type = cell_types_ordered)
  dot_grid_c <- merge(dot_grid_c,
                       ct_best_c[, .(gene, cell_type, PP.H4)],
                       by = c("gene", "cell_type"), all.x = TRUE)
  dot_grid_c[is.na(PP.H4), PP.H4 := 0]
  dot_grid_c[, cell_type_clean := celltype_clean_map[cell_type]]
  dot_grid_c[, cell_type_clean := factor(cell_type_clean,
    levels = c("Hepatocyte", "Endothelial", "Cholangiocyte", "Stellate"))]
  dot_grid_c[, gene := factor(gene, levels = gene_order_c)]

  # Bold known MASLD genes
  known_masld_c <- c("HSD17B13", "PNPLA3", "TM6SF2", "MBOAT7", "GCKR",
                     "CIDEC", "PPARG", "COL1A1", "THRB", "EFHD1",
                     "FABP1", "RORA", "HKDC1", "SPTLC3")
  y_faces_c <- ifelse(gene_order_c %in% known_masld_c, "bold.italic", "italic")
  names(y_faces_c) <- gene_order_c

  p_c <- ggplot(dot_grid_c[PP.H4 > 0], aes(x = cell_type_clean, y = gene)) +
    geom_point(aes(size = PP.H4, color = cell_type_clean), alpha = 0.8, shape = 16) +
    scale_size_continuous(range = c(0.5, 4), limits = c(0, 1),
                          breaks = c(0.3, 0.5, 0.8, 1.0),
                          name = "PP.H4") +
    scale_color_manual(values = celltype_colors, name = "Cell type") +
    labs(x = NULL, y = NULL,
         title = paste0("sc-eQTL COLOC (MASLD eQTL x Ghodsian/UKBB, ",
                        length(top_c_genes), " genes)")) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          axis.text.y = element_text(size = 5, face = y_faces_c),
          legend.position = "right",
          legend.key.size = unit(0.25, "cm"),
          plot.margin = margin(2, 2, 2, 2))
}

# ==========================================================================
# Panel (d): Causal method coverage bar chart
#   Consolidated: ~5 bars (TWAS, COLOC-Broadaway union, sc-TWAS,
#   sc-eQTL COLOC, ieQTL) + S3 summary bar
#   Removed: Ghodsian COLOC (0 hits); MR bar removed 2026-04-22 (MR ditched).
# ==========================================================================
me <- load_multi_evidence()

if (!is.null(me) && nrow(me) > 0) {
  dream_degs <- me[dream_padj < 0.1 & abs(dream_logFC) > 0.5, human_symbol]

  # --- COLOC (Broadaway): union across 4 GWAS ---
  coloc_cols <- c("ukbb_alt_coloc_pp4", "ast_coloc_pp4",
                  "ggt_coloc_pp4", "pdff_coloc_pp4")
  coloc_cols_present <- intersect(coloc_cols, names(me))
  if (length(coloc_cols_present) > 0) {
    me[, coloc_any_pp4 := do.call(pmax, c(.SD, na.rm = TRUE)),
       .SDcols = coloc_cols_present]
    me[is.infinite(coloc_any_pp4), coloc_any_pp4 := NA]
  } else {
    me[, coloc_any_pp4 := NA_real_]
  }

  # --- sc-TWAS: load directly from source (not yet in atlas) ---
  sceqtl_twas_f <- file.path(CAUSAL, "sceqtl_twas", "sceqtl_twas_all_results.csv")
  n_sc_twas_tested <- 0L; n_sc_twas_sig <- 0L; n_sc_twas_deg <- 0L
  if (file.exists(sceqtl_twas_f)) {
    sc_twas_raw <- fread(sceqtl_twas_f)
    # Best FDR per gene (across all cell types, GWAS, conditions, methods)
    sc_twas_best <- sc_twas_raw[, .(best_fdr = min(fdr, na.rm = TRUE)),
                                  by = exposure]
    sc_twas_best <- sc_twas_best[is.finite(best_fdr)]
    # Map Ensembl IDs to symbols for overlap
    gm <- load_gene_map()
    sc_twas_best[, ensembl_clean := sub("\\..*", "", exposure)]
    sc_twas_best <- merge(sc_twas_best, gm, by = "ensembl_clean", all.x = TRUE)
    n_sc_twas_tested <- nrow(sc_twas_best)
    n_sc_twas_sig <- sc_twas_best[best_fdr < 0.05, .N]
    n_sc_twas_deg <- sc_twas_best[best_fdr < 0.05 & symbol %in% dream_degs, .N]
  }

  # --- S3 summary: any causal evidence ---
  # Gene has S3 evidence if any of: TWAS, COLOC, sc-eQTL COLOC, ieQTL, sc-TWAS
  # (MR removed 2026-04-22 — MR ditched from paper)
  me[, has_s3 := FALSE]
  if ("twas_pval" %in% names(me))
    me[!is.na(twas_pval) & twas_pval < 0.05, has_s3 := TRUE]
  me[!is.na(coloc_any_pp4) & coloc_any_pp4 > 0.5, has_s3 := TRUE]
  if ("sceqtl_coloc_best_pp4" %in% names(me))
    me[!is.na(sceqtl_coloc_best_pp4) & sceqtl_coloc_best_pp4 > 0.5, has_s3 := TRUE]
  if ("ieqtl_disease_interaction" %in% names(me))
    me[ieqtl_disease_interaction %in% c(TRUE, "TRUE") &
       !is.na(ieqtl_interaction_pval) & ieqtl_interaction_pval < 1e-5, has_s3 := TRUE]

  method_stats <- rbindlist(list(
    data.table(
      method = "TWAS (GTEx)",
      n_tested  = me[!is.na(twas_pval), .N],
      n_sig     = me[!is.na(twas_pval) & twas_pval < 0.05, .N],
      n_deg_overlap = me[!is.na(twas_pval) & twas_pval < 0.05 &
                          human_symbol %in% dream_degs, .N]
    ),
    # "MR (GTEx)" bar removed 2026-04-22 — MR ditched from paper.
    data.table(
      method = "COLOC (Broadaway, 4 GWAS)",
      n_tested  = me[!is.na(coloc_any_pp4), .N],
      n_sig     = me[!is.na(coloc_any_pp4) & coloc_any_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(coloc_any_pp4) & coloc_any_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "sc-TWAS (4 cell types)",
      n_tested  = n_sc_twas_tested,
      n_sig     = n_sc_twas_sig,
      n_deg_overlap = n_sc_twas_deg
    ),
    data.table(
      method = "COLOC (sc-eQTL)",
      n_tested  = me[!is.na(sceqtl_coloc_best_pp4), .N],
      n_sig     = me[!is.na(sceqtl_coloc_best_pp4) & sceqtl_coloc_best_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(sceqtl_coloc_best_pp4) & sceqtl_coloc_best_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "ieQTL (disease-interaction)",
      n_tested  = me[!is.na(ieqtl_interaction_pval), .N],
      n_sig     = me[ieqtl_disease_interaction %in% c(TRUE, "TRUE") &
                     !is.na(ieqtl_interaction_pval) & ieqtl_interaction_pval < 1e-5, .N],
      n_deg_overlap = me[ieqtl_disease_interaction %in% c(TRUE, "TRUE") &
                         !is.na(ieqtl_interaction_pval) & ieqtl_interaction_pval < 1e-5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "S3: Any causal evidence",
      n_tested  = nrow(me),
      n_sig     = me[has_s3 == TRUE, .N],
      n_deg_overlap = me[has_s3 == TRUE & human_symbol %in% dream_degs, .N]
    )
  ))

  # Stacked bar: sig + overlap breakdown
  method_stats[, n_sig_nondeg := n_sig - n_deg_overlap]
  long_d <- melt(method_stats[, .(method, n_deg_overlap, n_sig_nondeg)],
                 id.vars = "method",
                 measure.vars = c("n_deg_overlap", "n_sig_nondeg"),
                 variable.name = "category", value.name = "n")
  long_d[, category := fifelse(category == "n_deg_overlap",
                                "Sig & DEG overlap", "Sig (not DEG)")]
  long_d[, category := factor(category,
                               levels = c("Sig & DEG overlap", "Sig (not DEG)"))]

  # Order: method type grouping
  # "MR (GTEx)" removed from method_order 2026-04-22 — MR ditched from paper.
  method_order <- c("TWAS (GTEx)",
                    "COLOC (Broadaway, 4 GWAS)", "sc-TWAS (4 cell types)",
                    "COLOC (sc-eQTL)", "ieQTL (disease-interaction)",
                    "S3: Any causal evidence")
  long_d[, method := factor(method, levels = rev(method_order))]
  method_stats[, method_f := factor(method, levels = rev(method_order))]

  # S3 summary bar: distinct fill to separate from individual methods
  s3_pct <- round(100 * method_stats[method == "S3: Any causal evidence", n_sig] /
                    method_stats[method == "S3: Any causal evidence", n_tested], 1)

  p_d <- ggplot(long_d, aes(x = n, y = method, fill = category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(data = method_stats,
              aes(x = n_sig + max(method_stats$n_sig) * 0.04,
                  y = method_f, label = paste0(n_sig, " / ", n_tested)),
              inherit.aes = FALSE, size = 1.7, hjust = 0, color = "gray30") +
    scale_fill_manual(
      values = c("Sig & DEG overlap" = masld_colors$up,
                 "Sig (not DEG)"     = masld_colors$twas),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(x = "Number of genes", y = NULL,
         title = "Causal method coverage (S3)",
         subtitle = paste0("Significant / tested; S3 total: ",
                           s3_pct, "% of atlas genes")) +
    theme_masld() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.75, 0.15),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(0.25, "cm"),
          plot.subtitle = element_text(size = 5))

  # Clean up temporary columns
  me[, c("coloc_any_pp4", "has_s3") := NULL]
}

# ==========================================================================
# Panel (e): ieQTL x DEG directional concordance scatter
#   + quadrant annotations + both Pearson r and Spearman rho
# ==========================================================================
ieqtl <- load_ieqtl()

if (!is.null(ieqtl) && nrow(ieqtl) > 0) {
  ieqtl_plot <- ieqtl[!is.na(interaction_beta) & !is.na(dream_logFC)]

  if (nrow(ieqtl_plot) > 0) {
    # Best (most significant) interaction per gene
    ieqtl_best <- ieqtl_plot[, .SD[which.min(interaction_pval)], by = gene]
    ieqtl_best[, cell_type_clean := celltype_clean_map[cell_type]]
    ieqtl_best[is.na(cell_type_clean), cell_type_clean := cell_type]

    # Quadrant classification
    ieqtl_best[, quadrant := fifelse(
      interaction_beta > 0 & dream_logFC > 0, "Concordant Up",
      fifelse(interaction_beta < 0 & dream_logFC < 0, "Concordant Down",
              fifelse(interaction_beta > 0 & dream_logFC < 0, "Discordant",
                      "Discordant"))
    )]

    # Quadrant counts
    q_counts <- ieqtl_best[, .N, by = quadrant]
    n_ie <- nrow(ieqtl_best)
    n_concordant <- ieqtl_best[grepl("Concordant", quadrant), .N]
    pct_conc <- round(100 * n_concordant / n_ie, 1)

    # Both correlation metrics
    rho_ie <- cor(ieqtl_best$interaction_beta, ieqtl_best$dream_logFC,
                  method = "spearman", use = "complete.obs")
    r_ie <- cor(ieqtl_best$interaction_beta, ieqtl_best$dream_logFC,
                method = "pearson", use = "complete.obs")

    # Quadrant annotation positions (corners of the plot)
    x_range <- range(ieqtl_best$dream_logFC, na.rm = TRUE)
    y_range <- range(ieqtl_best$interaction_beta, na.rm = TRUE)
    x_pad <- diff(x_range) * 0.05
    y_pad <- diff(y_range) * 0.05

    # Compute per-quadrant counts directly from ieqtl_best
    n_disc_ul <- ieqtl_best[interaction_beta > 0 & dream_logFC < 0, .N]
    n_disc_lr <- ieqtl_best[interaction_beta < 0 & dream_logFC > 0, .N]
    n_conc_ur <- ieqtl_best[interaction_beta > 0 & dream_logFC > 0, .N]
    n_conc_ll <- ieqtl_best[interaction_beta < 0 & dream_logFC < 0, .N]

    quad_labels_simple <- data.table(
      x = c(x_range[2] - x_pad, x_range[1] + x_pad,
            x_range[1] + x_pad, x_range[2] - x_pad),
      y = c(y_range[2] - y_pad, y_range[1] + y_pad,
            y_range[2] - y_pad, y_range[1] + y_pad),
      label = c(paste0("Conc. Up\n", n_conc_ur),
                paste0("Conc. Down\n", n_conc_ll),
                paste0("Disc.\n", n_disc_ul),
                paste0("Disc.\n", n_disc_lr)),
      hjust = c(1, 0, 0, 1),
      vjust = c(1, 1, 0, 0)
    )

    # Genes to label
    label_genes_e <- c("THRB", "HSD17B13", "PNPLA3", "FABP1", "RORA",
                       "SCD", "FASN", "COL1A1", "SPP1", "ACSL4")
    top_sig_e <- head(ieqtl_best[order(interaction_pval)], 10)
    label_dt_e <- ieqtl_best[gene %in% label_genes_e | gene %in% top_sig_e$gene]
    label_dt_e <- label_dt_e[, .SD[which.min(interaction_pval)], by = gene]
    if (nrow(label_dt_e) > 15) label_dt_e <- head(label_dt_e[order(interaction_pval)], 15)

    p_e <- ggplot(ieqtl_best, aes(x = dream_logFC, y = interaction_beta,
                                    color = cell_type_clean)) +
      rasterize_layer(
        geom_point(size = 0.5, alpha = 0.5, shape = 16)
      ) +
      geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dashed",
                 color = "gray50") +
      geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed",
                 color = "gray50") +
      # Quadrant annotations
      geom_text(data = quad_labels_simple,
                aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
                inherit.aes = FALSE, size = 1.5, color = "gray50",
                fontface = "italic", lineheight = 0.85) +
      geom_text_repel(data = label_dt_e,
                      aes(label = gene), size = 1.7,
                      max.overlaps = 20, segment.size = 0.2,
                      min.segment.length = 0, box.padding = 0.3,
                      color = "black") +
      scale_color_manual(values = celltype_colors, name = "Cell type") +
      annotate("text", x = Inf, y = -Inf,
               label = paste0("rho = ", round(rho_ie, 3),
                              "; r = ", round(r_ie, 3),
                              "\nn = ", comma(n_ie),
                              "\n", pct_conc, "% concordant"),
               hjust = 1.1, vjust = -0.3, size = 2, fontface = "italic",
               lineheight = 0.85) +
      labs(x = "Integrated logFC (disease vs control)",
           y = "ieQTL interaction beta",
           title = "ieQTL x DEG directional concordance") +
      theme_masld() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.2, "cm"),
            legend.margin = margin(0, 0, 0, 0))
  }
}

# ==========================================================================
# Panel (f): Multi-ancestry COLOC summary
#   Bars for PP.H4 > 0.8; cross-ancestry overlap
#   BBJ = East Asian; PanUKBB = AFR/CSA; Broadaway = primary European
# ==========================================================================
ancestry_sources <- list()

# BBJ liver enzyme COLOC (East Asian x Broadaway eQTLs)
bbj_dirs <- c("BBJ ALT" = "bbj_alt", "BBJ AST" = "bbj_ast", "BBJ GGT" = "bbj_ggt")
bbj_genes_08 <- list()
for (gwas_label in names(bbj_dirs)) {
  f <- file.path(CAUSAL, bbj_dirs[[gwas_label]], "coloc_results.csv")
  if (file.exists(f)) {
    tmp <- fread(f)
    if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
      tmp <- add_symbols(tmp, "gene")
      ancestry_sources[[gwas_label]] <- data.table(
        gwas = gwas_label,
        ancestry = "East Asian (BBJ)",
        n_tested = nrow(tmp),
        n_coloc = sum(tmp$PP.H4 > 0.8, na.rm = TRUE)
      )
      bbj_genes_08[[gwas_label]] <- tmp[PP.H4 > 0.8, unique(symbol)]
    }
  }
}

# RESTORED 2026-04-09: FinnGen COLOC restored (FinnGen_NAFLD, FinnGen_NASH, FinnGen_HCC verified R12)
finngen_dirs <- c("FinnGen NAFLD" = "finngen_nafld",
                  "FinnGen NASH"  = "finngen_nash",
                  "FinnGen HCC"   = "finngen_hcc")
finngen_genes_08 <- list()
for (gwas_label in names(finngen_dirs)) {
  f <- file.path(CAUSAL, finngen_dirs[[gwas_label]], "coloc_results.csv")
  if (file.exists(f)) {
    tmp <- fread(f)
    if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
      tmp <- add_symbols(tmp, "gene")
      ancestry_sources[[gwas_label]] <- data.table(
        gwas = gwas_label,
        ancestry = "European (FinnGen)",
        n_tested = nrow(tmp),
        n_coloc = sum(tmp$PP.H4 > 0.8, na.rm = TRUE)
      )
      finngen_genes_08[[gwas_label]] <- tmp[PP.H4 > 0.8, unique(symbol)]
    }
  }
}

# Pan-UKBB AFR COLOC (African, exploratory — underpowered)
panukbb_afr_dirs <- c("PanUKBB AFR ALT" = "panukbb_afr_alt",
                      "PanUKBB AFR AST" = "panukbb_afr_ast",
                      "PanUKBB AFR GGT" = "panukbb_afr_ggt")
for (gwas_label in names(panukbb_afr_dirs)) {
  f <- file.path(CAUSAL, panukbb_afr_dirs[[gwas_label]], "coloc_results.csv")
  if (file.exists(f)) {
    tmp <- fread(f)
    if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
      tmp <- add_symbols(tmp, "gene")
      ancestry_sources[[gwas_label]] <- data.table(
        gwas = gwas_label,
        ancestry = "African (PanUKBB)",
        n_tested = nrow(tmp),
        n_coloc = sum(tmp$PP.H4 > 0.8, na.rm = TRUE)
      )
    }
  }
}

# Pan-UKBB CSA COLOC (Central/South Asian, exploratory — underpowered)
panukbb_csa_dirs <- c("PanUKBB CSA ALT" = "panukbb_csa_alt",
                      "PanUKBB CSA AST" = "panukbb_csa_ast",
                      "PanUKBB CSA GGT" = "panukbb_csa_ggt")
for (gwas_label in names(panukbb_csa_dirs)) {
  f <- file.path(CAUSAL, panukbb_csa_dirs[[gwas_label]], "coloc_results.csv")
  if (file.exists(f)) {
    tmp <- fread(f)
    if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
      tmp <- add_symbols(tmp, "gene")
      ancestry_sources[[gwas_label]] <- data.table(
        gwas = gwas_label,
        ancestry = "C/S Asian (PanUKBB)",
        n_tested = nrow(tmp),
        n_coloc = sum(tmp$PP.H4 > 0.8, na.rm = TRUE)
      )
    }
  }
}

# Broadaway European COLOC (primary) for reference
broadaway_eu_dirs <- c("UKBB ALT" = "broadaway_ukbb",
                       "UKBB AST" = "broadaway_ukbb_ast",
                       "UKBB GGT" = "broadaway_ukbb_ggt",
                       "PDFF"     = "broadaway_pdff")
eu_genes_08 <- list()
for (gwas_label in names(broadaway_eu_dirs)) {
  f <- file.path(CAUSAL, broadaway_eu_dirs[[gwas_label]], "coloc_results.csv")
  if (file.exists(f)) {
    tmp <- fread(f)
    if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
      tmp <- add_symbols(tmp, "gene")
      ancestry_sources[[gwas_label]] <- data.table(
        gwas = gwas_label,
        ancestry = "European (UKBB)",
        n_tested = nrow(tmp),
        n_coloc = sum(tmp$PP.H4 > 0.8, na.rm = TRUE)
      )
      eu_genes_08[[gwas_label]] <- tmp[PP.H4 > 0.8, unique(symbol)]
    }
  }
}

if (length(ancestry_sources) > 0) {
  anc_dt <- rbindlist(ancestry_sources)

  # Cross-ancestry overlap: genes in BOTH European AND East Asian at PP.H4 > 0.8
  eu_all <- unique(unlist(eu_genes_08))
  bbj_all <- unique(unlist(bbj_genes_08))
  cross_ancestry_overlap <- intersect(eu_all, bbj_all)
  n_cross_ancestry <- length(cross_ancestry_overlap)

  # Strip source prefix from GWAS label to get short trait name
  anc_dt[, trait := gsub("^(UKBB|BBJ|PanUKBB AFR|PanUKBB CSA)\\s+", "", gwas)]

  # Source grouping label (used as facet strip)
  # RESTORED 2026-04-09: FinnGen source mapping restored
  anc_dt[, source := fcase(
    ancestry == "European (UKBB)",     "UKBB",
    ancestry == "European (FinnGen)",  "FinnGen",
    ancestry == "East Asian (BBJ)",    "BBJ",
    ancestry == "African (PanUKBB)",   "PanUKBB AFR",
    ancestry == "C/S Asian (PanUKBB)", "PanUKBB CSA"
  )]

  # Order: UKBB first, then FinnGen, BBJ, PanUKBB AFR, PanUKBB CSA
  source_order <- c("UKBB", "FinnGen", "BBJ", "PanUKBB AFR", "PanUKBB CSA")
  source_order <- source_order[source_order %in% anc_dt$source]
  anc_dt[, source := factor(source, levels = source_order)]

  # Order traits within each source
  gwas_order <- c(names(broadaway_eu_dirs), names(finngen_dirs), names(bbj_dirs),
                  names(panukbb_afr_dirs), names(panukbb_csa_dirs))
  gwas_order <- gwas_order[gwas_order %in% anc_dt$gwas]
  anc_dt[, gwas := factor(gwas, levels = rev(gwas_order))]
  # Trait order within facets (reversed for coord_flip-like horizontal bars)
  trait_order <- anc_dt[order(gwas), unique(trait)]
  anc_dt[, trait := factor(trait, levels = rev(trait_order))]

  # Ancestry colors
  # RESTORED 2026-04-09: FinnGen color restored
  ancestry_fill <- c(
    "European (UKBB)"      = "#1565C0",
    "European (FinnGen)"   = "#2E7D32",
    "East Asian (BBJ)"     = "#E91E63",
    "African (PanUKBB)"    = "#FF8F00",
    "C/S Asian (PanUKBB)"  = "#FFB74D"
  )

  p_f <- ggplot(anc_dt, aes(x = n_coloc, y = trait, fill = ancestry)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = n_coloc),
              hjust = -0.15, size = 1.3, color = "gray30", show.legend = FALSE) +
    facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = ancestry_fill, name = NULL,
                      guide = guide_legend(nrow = 2)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Genes (PP.H4 > 0.8)",
         y = NULL,
         title = "Multi-ancestry COLOC (Broadaway eQTLs)") +
    theme_masld() +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.15, "cm"),
          legend.text = element_text(size = 4.5),
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing.x = unit(0.1, "cm"),
          plot.title = element_text(size = 7),
          plot.margin = margin(2, 2, 2, 2),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          axis.title.x = element_text(size = 6),
          strip.placement = "outside",
          strip.text.y.left = element_text(size = 5, angle = 0, hjust = 1),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.1, "cm"))
}

# ==========================================================================
# Panel (g): TWAS — Ghodsian vs Chen z-score concordance scatter
#   Each dot is a gene tested in both GWAS; color = significance category
# ==========================================================================
twas_combined <- load_twas_combined()

if (!is.null(twas_combined) && nrow(twas_combined) > 0) {
  # Pivot: one row per gene with z-scores from both GWAS
  twas_gh <- twas_combined[gwas == "ghodsian", .(gene, symbol, z_ghodsian = zscore,
                                                   fdr_ghodsian = fdr)]
  twas_ch <- twas_combined[gwas == "chen", .(gene, z_chen = zscore, fdr_chen = fdr)]
  twas_both <- merge(twas_gh, twas_ch, by = "gene")

  if (nrow(twas_both) > 10) {
    # Significance categories
    twas_both[, sig_cat := fcase(
      fdr_ghodsian < 0.05 & fdr_chen < 0.05, "Both FDR<0.05",
      fdr_ghodsian < 0.05 | fdr_chen < 0.05, "One FDR<0.05",
      default = "NS"
    )]
    twas_both[, sig_cat := factor(sig_cat,
      levels = c("Both FDR<0.05", "One FDR<0.05", "NS"))]

    # Correlation
    rho <- cor(twas_both$z_ghodsian, twas_both$z_chen, method = "spearman",
               use = "complete.obs")
    n_both_sig <- twas_both[sig_cat == "Both FDR<0.05", .N]

    # Label top concordant hits
    twas_both[, abs_z_mean := (abs(z_ghodsian) + abs(z_chen)) / 2]
    known_twas <- c("THRB", "HSD17B13", "PNPLA3", "TM6SF2", "MBOAT7", "GCKR",
                     "MARC1", "SORT1", "EFHD1", "SLC39A8", "FABP1")
    label_genes <- twas_both[sig_cat == "Both FDR<0.05"][order(-abs_z_mean)][1:min(15, .N)]
    label_genes <- rbind(label_genes,
                          twas_both[symbol %in% known_twas & sig_cat != "NS"],
                          fill = TRUE)
    label_genes <- label_genes[!duplicated(gene)]

    sig_colors <- c("Both FDR<0.05" = masld_colors$up,
                     "One FDR<0.05"  = masld_colors$twas,
                     "NS"            = "grey80")

    p_g <- ggplot(twas_both, aes(x = z_ghodsian, y = z_chen, color = sig_cat)) +
      geom_point(size = 0.3, alpha = 0.5) +
      geom_point(data = twas_both[sig_cat != "NS"], size = 0.6, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50",
                   linewidth = 0.3) +
      geom_text_repel(data = label_genes,
                       aes(label = symbol), size = 1.8, max.overlaps = 20,
                       segment.size = 0.2, min.segment.length = 0.3,
                       color = "gray20", fontface = "italic") +
      scale_color_manual(values = sig_colors, name = NULL) +
      labs(x = "TWAS z-score (Ghodsian NAFLD)",
           y = "TWAS z-score (Chen NAFLD)",
           title = paste0("TWAS cross-GWAS concordance (rho=",
                           round(rho, 2), ", n=", n_both_sig, " both sig)")) +
      theme_masld() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.15, "cm"),
            legend.text = element_text(size = 5),
            plot.title = element_text(size = 7),
            axis.title = element_text(size = 6))

    cat("Panel g: TWAS Ghodsian vs Chen scatter — done\n")
  }
}

# ==========================================================================
# Panel (h): MR forest plot — REMOVED 2026-04-22 (MR ditched from paper)
# Panel h is no longer assembled into Fig 3; TWAS + COLOC + INTACT cover the
# causal axis (panels a-g, i). Placeholder kept as NULL so p_h references in
# the patchwork layout below resolve without error.
# ==========================================================================
p_h <- patchwork::plot_spacer()
cat("Panel h: MR forest REMOVED (MR ditched 2026-04-22) — spacer inserted\n")

# ==========================================================================
# Panel (i): sc-TWAS cell-type dot plot
#   Gene x cell type, size = -log10(FDR), color = cell type
#   MASLD-condition eQTLs, best GWAS per gene-celltype
# ==========================================================================
sceqtl_twas_f <- file.path(CAUSAL, "sceqtl_twas", "sceqtl_twas_all_results.csv")

if (file.exists(sceqtl_twas_f)) {
  sc_twas_raw <- fread(sceqtl_twas_f)

  # Use MASLD-condition eQTLs (more instruments, biologically relevant)
  sc_twas_masld <- sc_twas_raw[condition == "masld"]

  if (nrow(sc_twas_masld) > 0) {
    # Best FDR per gene x cell type (across GWAS and MR methods)
    sc_twas_best <- sc_twas_masld[, .(best_fdr = min(fdr, na.rm = TRUE),
                                       best_b = b[which.min(fdr)]),
                                    by = .(exposure, cell_type)]
    sc_twas_best <- sc_twas_best[is.finite(best_fdr)]
    sc_twas_best[, symbol := exposure]  # exposure is already gene symbol

    # Per-gene best FDR (across cell types) for ranking
    gene_rank <- sc_twas_best[, .(min_fdr = min(best_fdr),
                                   n_ct_sig = sum(best_fdr < 0.05)),
                                by = symbol]
    setorder(gene_rank, min_fdr)

    # Top 20 genes
    top_sc_genes <- head(gene_rank[min_fdr < 0.05], 20)$symbol
    if (length(top_sc_genes) < 10) {
      top_sc_genes <- head(gene_rank, 20)$symbol
    }

    gene_order_i <- rev(gene_rank[symbol %in% top_sc_genes, symbol])

    ct_ordered <- c("hepatocyte", "endothelial_cell", "cholangiocyte", "stellate_cell")
    dot_grid_i <- CJ(symbol = top_sc_genes, cell_type = ct_ordered)
    dot_grid_i <- merge(dot_grid_i,
                         sc_twas_best[, .(symbol, cell_type, best_fdr, best_b)],
                         by = c("symbol", "cell_type"), all.x = TRUE)
    dot_grid_i[is.na(best_fdr), best_fdr := 1]
    dot_grid_i[, neg_log_fdr := pmin(-log10(best_fdr), 10)]
    dot_grid_i[, cell_type_clean := celltype_clean_map[cell_type]]
    dot_grid_i[, cell_type_clean := factor(cell_type_clean,
      levels = c("Hepatocyte", "Endothelial", "Cholangiocyte", "Stellate"))]
    dot_grid_i[, symbol := factor(symbol, levels = gene_order_i)]

    # Bold known genes
    known_sc <- c("HSD17B13", "PNPLA3", "TM6SF2", "MBOAT7", "GCKR",
                   "CIDEC", "PPARG", "COL1A1", "THRB", "EFHD1",
                   "FABP1", "RORA", "HKDC1", "SPTLC3")
    y_faces_i <- ifelse(gene_order_i %in% known_sc, "bold.italic", "italic")
    names(y_faces_i) <- gene_order_i

    p_i <- ggplot(dot_grid_i[best_fdr < 1], aes(x = cell_type_clean, y = symbol)) +
      geom_point(aes(size = neg_log_fdr, color = cell_type_clean), alpha = 0.8, shape = 16) +
      scale_size_continuous(range = c(0.5, 4),
                             breaks = c(1, 2, 5, 10),
                             name = expression(-log[10](FDR))) +
      scale_color_manual(values = celltype_colors, name = "Cell type") +
      labs(x = NULL, y = NULL,
           title = paste0("sc-TWAS by cell type (",
                           length(top_sc_genes), " genes, MASLD eQTLs)")) +
      theme_masld() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5, face = y_faces_i),
            legend.position = "right",
            legend.key.size = unit(0.25, "cm"),
            plot.title = element_text(size = 7),
            plot.margin = margin(2, 2, 2, 2))

    cat("Panel i: sc-TWAS cell-type dot plot — done\n")
  }
}

# ==========================================================================
# Assemble: expanded layout (9 panels, 3 columns x 3 rows)
# ==========================================================================
row1 <- p_a + p_b + plot_layout(widths = c(1.2, 0.8))
row2 <- p_c + p_d + plot_layout(widths = c(0.9, 1.1))
row3 <- p_e + p_f + plot_layout(widths = c(1, 1))
row4 <- p_g + p_h + p_i + plot_layout(widths = c(1, 0.8, 0.8))

fig3 <- (row1 / row2 / row3 / row4) +
  plot_layout(heights = c(1.2, 1.2, 1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(fig3, OUT, height = 14)

# Save individual panels
panel_dir <- file.path(FIG3_DIR, "panels")
save_fig(p_a, file.path(panel_dir, "fig3a_coloc_manhattan.pdf"), height = 3.5)
save_fig(p_b, file.path(panel_dir, "fig3b_coloc_heatmap.pdf"), height = 3.5)
save_fig(p_c, file.path(panel_dir, "fig3c_sceqtl_dotplot.pdf"), width = fig_half_width, height = 3.5)
save_fig(p_d, file.path(panel_dir, "fig3d_causal_coverage.pdf"), height = 3.5)
save_fig(p_e, file.path(panel_dir, "fig3e_ieqtl_concordance.pdf"), height = 3.5)
save_fig(p_f, file.path(panel_dir, "fig3f_multiancestry.pdf"), width = fig_half_width, height = 2.8)
save_fig(p_g, file.path(panel_dir, "fig3g_twas_concordance.pdf"), height = 3.5)
# Panel h (MR forest) removed 2026-04-22 — MR ditched from paper.
save_fig(p_i, file.path(panel_dir, "fig3i_sctwas_dotplot.pdf"), width = fig_half_width, height = 3.5)

message("Fig 3 saved to ", OUT)
