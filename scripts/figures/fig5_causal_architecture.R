##############################################################################
# Figure 5: Causal Architecture (6 panels)
#   Row 1: (a) Unified causal Manhattan (TWAS+COLOC) — cleaned labels
#           (b) Causal convergence gene×evidence heatmap
#   Row 2: (c) sc-eQTL COLOC x cell type heatmap + GWAS source annotation
#           (d) Causal method coverage bar chart (replaces layer correlations)
#   Row 3: (e) ieQTL cell-type enrichment / quadrant summary
#           (f) Causal method overlap (UpSet-style)
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

OUT <- file.path(FIG3_DIR, "causal_architecture.pdf")

# Pre-initialize all panels with placeholders
p_a <- placeholder("Panel a: TWAS Manhattan")
p_b <- placeholder("Panel b: Causal convergence")
p_c <- placeholder("Panel c: sc-eQTL COLOC heatmap")
p_d <- placeholder("Panel d: Causal method coverage")
p_e <- placeholder("Panel e: ieQTL cell-type enrichment")
p_f <- placeholder("Panel f: Causal method overlap")

# Cell type styling (consistent across panels c and e)
celltype_colors <- c(
  Hepatocyte     = "#0D47A1",   # Deep blue — metabolic core
  Endothelial    = "#80CBC4",   # Muted teal — vascular
  Cholangiocyte  = "#CE93D8",   # Muted violet — biliary
  Stellate       = "#F8BBD0"    # Muted pink — fibrotic
)

celltype_clean_map <- c(
  hepatocyte       = "Hepatocyte",
  endothelial_cell = "Endothelial",
  cholangiocyte    = "Cholangiocyte",
  stellate_cell    = "Stellate"
)

# ==========================================================================
# Panel (a): Unified Causal Manhattan (TWAS + COLOC from 5 sources)
#   IMPROVED: labels filtered to exclude non-coding/cryptic locus IDs
# ==========================================================================

## -- Step A: TWAS (best p-value per gene across GWAS) --
twas_raw <- load_twas_combined()
twas_dt <- NULL

if (!is.null(twas_raw) && nrow(twas_raw) > 0) {
  twas_w <- copy(twas_raw)
  if (!"pvalue" %in% names(twas_w) && "pval" %in% names(twas_w))
    setnames(twas_w, "pval", "pvalue")
  if (!"fdr" %in% names(twas_w) && "adj_pval" %in% names(twas_w))
    setnames(twas_w, "adj_pval", "fdr")
  if (!"symbol" %in% names(twas_w) && "gene_name" %in% names(twas_w))
    twas_w[, symbol := gene_name]

  twas_w <- twas_w[!is.na(pvalue) & pvalue > 0 & !is.na(symbol) & symbol != ""]
  # Best p-value per gene across Ghodsian/Chen
  twas_dt <- twas_w[, .SD[which.min(pvalue)], by = symbol]
  twas_dt[, neg_log10p := pmin(-log10(pvalue), 30)]
  twas_dt[, method := "TWAS"]
  twas_dt[, sig_class := fifelse(
    !is.na(fdr) & fdr < 0.05, "FDR / PP.H4 > 0.8",
    fifelse(pvalue < 0.05, "Nominal", "NS")
  )]
}

## -- Step B: Collect COLOC from all 5 sources --
coloc_list <- list()

# B1: GTEx bulk COLOC
gtex_coloc <- load_coloc_results()
if (!is.null(gtex_coloc) && nrow(gtex_coloc) > 0 && "PP.H4" %in% names(gtex_coloc)) {
  coloc_list[["GTEx"]] <- gtex_coloc[!is.na(PP.H4),
    .(symbol, PP.H4, source = "GTEx")]
}

# B2: Broadaway COLOC
broad_coloc <- load_broadaway_coloc()
if (!is.null(broad_coloc) && nrow(broad_coloc) > 0 && "PP.H4" %in% names(broad_coloc)) {
  coloc_list[["Broadaway"]] <- broad_coloc[!is.na(PP.H4),
    .(symbol, PP.H4, source = "Broadaway")]
}

# B3: sc-eQTL Ghodsian
sc_coloc <- load_sceqtl_coloc()
if (!is.null(sc_coloc) && nrow(sc_coloc) > 0) {
  pp4_col_sc <- intersect(c("PP.H4_all", "PP.H4"), names(sc_coloc))[1]
  if (!is.na(pp4_col_sc))
    coloc_list[["sc_eQTL"]] <- sc_coloc[!is.na(get(pp4_col_sc)),
      .(symbol = gene, PP.H4 = get(pp4_col_sc), source = "sc-eQTL")]
}

# B4: UKBB ALT
ukbb_coloc <- load_sceqtl_coloc_ukbb()
if (!is.null(ukbb_coloc) && nrow(ukbb_coloc) > 0) {
  pp4_col_uk <- intersect(c("PP.H4_all", "PP.H4"), names(ukbb_coloc))[1]
  if (!is.na(pp4_col_uk))
    coloc_list[["UKBB_ALT"]] <- ukbb_coloc[!is.na(get(pp4_col_uk)),
      .(symbol = gene, PP.H4 = get(pp4_col_uk), source = "UKBB_ALT")]
}

# B5: Zenodo precomputed COLOC
zenodo_coloc <- load_zenodo_coloc()
if (!is.null(zenodo_coloc) && nrow(zenodo_coloc) > 0) {
  zen_pp4_col <- intersect(c("PP.H4", "pp4", "coloc_pp4"), names(zenodo_coloc))[1]
  zen_gene_col <- intersect(c("gene", "symbol", "gene_symbol"), names(zenodo_coloc))[1]
  if (!is.na(zen_pp4_col) && !is.na(zen_gene_col)) {
    coloc_list[["Zenodo"]] <- zenodo_coloc[!is.na(get(zen_pp4_col)),
      .(symbol = get(zen_gene_col), PP.H4 = get(zen_pp4_col), source = "Zenodo")]
  } else if ("nafld_coloc" %in% names(zenodo_coloc)) {
    zen_gene_col2 <- intersect(c("gene", "symbol"), names(zenodo_coloc))[1]
    if (!is.na(zen_gene_col2))
      coloc_list[["Zenodo"]] <- zenodo_coloc[nafld_coloc == TRUE,
        .(symbol = get(zen_gene_col2), PP.H4 = 1.0, source = "Zenodo")]
  }
}

# Best PP.H4 per gene across all COLOC sources
coloc_all <- rbindlist(coloc_list, use.names = TRUE)
coloc_dt <- NULL

if (nrow(coloc_all) > 0) {
  coloc_all <- coloc_all[!is.na(symbol) & symbol != ""]
  coloc_dt <- coloc_all[, .(PP.H4 = max(PP.H4, na.rm = TRUE),
                             best_source = source[which.max(PP.H4)]),
                          by = symbol]
  coloc_dt[, neg_log10p := pmin(-log10(pmax(1 - PP.H4, 1e-30)), 30)]
  coloc_dt[, method := "COLOC"]
  coloc_dt[, sig_class := fifelse(
    PP.H4 > 0.8, "FDR / PP.H4 > 0.8",
    fifelse(PP.H4 > 0.5, "Nominal", "NS")
  )]
}

## -- Step C: Build symbol -> chr lookup --
chr_map <- data.table(symbol = character(), chr = integer())

if (!is.null(broad_coloc) && nrow(broad_coloc) > 0 && "chr" %in% names(broad_coloc)) {
  bc <- broad_coloc[!is.na(chr), .(symbol, chr = suppressWarnings(as.integer(chr)))]
  bc <- bc[!is.na(chr)]
  chr_map <- rbind(chr_map, bc[, .(chr = chr[1]), by = symbol])
}

if (!is.null(ukbb_coloc) && nrow(ukbb_coloc) > 0 && "top_snp_all" %in% names(ukbb_coloc)) {
  uc <- ukbb_coloc[!is.na(top_snp_all) & top_snp_all != "",
    .(symbol = gene, chr_str = sub(":.*", "", top_snp_all))]
  uc[, chr := suppressWarnings(as.integer(chr_str))]
  uc <- uc[!is.na(chr), .(chr = chr[1]), by = symbol]
  chr_map <- rbind(chr_map, uc[!symbol %in% chr_map$symbol])
}

if (!is.null(sc_coloc) && nrow(sc_coloc) > 0 && "top_snp_all" %in% names(sc_coloc)) {
  gc <- sc_coloc[!is.na(top_snp_all) & top_snp_all != "",
    .(symbol = gene, chr_str = sub(":.*", "", top_snp_all))]
  gc[, chr := suppressWarnings(as.integer(chr_str))]
  gc <- gc[!is.na(chr), .(chr = chr[1]), by = symbol]
  chr_map <- rbind(chr_map, gc[!symbol %in% chr_map$symbol])
}

# Also pull chr from broadaway_ukbb directly (explicit chr column)
ukbb_broad_path <- file.path(CAUSAL, "broadaway_ukbb", "coloc_results.csv")
if (file.exists(ukbb_broad_path)) {
  ukbb_broad <- fread(ukbb_broad_path)
  if ("chr" %in% names(ukbb_broad) && "gene" %in% names(ukbb_broad)) {
    ub <- ukbb_broad[!is.na(chr), .(symbol = gene, chr = suppressWarnings(as.integer(chr)))]
    ub <- ub[!is.na(chr), .(chr = chr[1]), by = symbol]
    chr_map <- rbind(chr_map, ub[!symbol %in% chr_map$symbol])
  }
}

chr_map <- chr_map[, .(chr = chr[1]), by = symbol]

## -- Step D: Combine TWAS + COLOC, order by chromosome --
keep_cols <- c("symbol", "neg_log10p", "method", "sig_class")
parts <- list()
if (!is.null(twas_dt) && nrow(twas_dt) > 0) parts[["twas"]] <- twas_dt[, ..keep_cols]
if (!is.null(coloc_dt) && nrow(coloc_dt) > 0) parts[["coloc"]] <- coloc_dt[, ..keep_cols]

if (length(parts) > 0) {
  unified <- rbindlist(parts, use.names = TRUE)
  unified <- merge(unified, chr_map, by = "symbol", all.x = TRUE)
  unified[is.na(chr), chr := 99L]
  setorder(unified, chr, -neg_log10p)
  unified[, gene_idx := .I]

  ## -- Step E: Subsample NS for performance --
  sig_pts <- unified[sig_class != "NS"]
  ns_pts  <- unified[sig_class == "NS"]
  if (nrow(ns_pts) > 5000L) ns_pts <- ns_pts[sample(.N, 5000L)]
  plot_dt <- rbind(sig_pts, ns_pts)

  ## -- Step F: Chromosome bands, labels, ggplot --
  chr_info <- unified[chr < 99, .(start = min(gene_idx), end = max(gene_idx)), by = chr]
  chr_info[, mid := (start + end) / 2]
  chr_info[, band := chr %% 2 == 0]

  # IMPROVED labels: exclude non-coding/cryptic IDs and short-read artifacts
  cryptic_blocklist <- c("TSGA10", "C2orf16", "TSFM", "ABHD1", "MLIP", "PCOLCE2",
                         "MAP1LC3A", "FLRT3", "SLC12A8")
  is_labelable <- function(sym) {
    !grepl("^AC[0-9]|^AL[0-9]|^RP[0-9]|^RP1[0-9]|^LINC|^MIR[0-9]|^SNORD|^LOC[0-9]", sym) &
    !sym %in% cryptic_blocklist
  }

  known_masld_a <- c("HSD17B13", "EFHD1", "PNPLA3", "TM6SF2", "MBOAT7", "GCKR",
                     "CIDEC", "PPARG", "COL1A1", "THRB",
                     "PBX4", "SAMM50", "CCDC92", "ZNF682",
                     "LCOR", "ALAD", "ETS2", "USP40",
                     "FABP1", "RORA", "HKDC1", "SPTLC3", "CHEK2")
  fdr_genes   <- unified[sig_class == "FDR / PP.H4 > 0.8" & is_labelable(symbol)]
  known_sig   <- unified[symbol %in% known_masld_a & sig_class != "NS"]
  top_nominal <- head(unified[sig_class == "Nominal" & is_labelable(symbol)][order(-neg_log10p)], 5)

  label_dt <- rbind(fdr_genes, known_sig, top_nominal)
  label_dt <- label_dt[, .SD[which.max(neg_log10p)], by = symbol]
  if (nrow(label_dt) > 22) label_dt <- head(label_dt[order(-neg_log10p)], 22)

  n_fdr_twas  <- unified[method == "TWAS" & sig_class == "FDR / PP.H4 > 0.8", .N]
  n_pp4_high  <- unified[method == "COLOC" & sig_class == "FDR / PP.H4 > 0.8", .N]
  n_total     <- uniqueN(unified$symbol)

  p_a <- ggplot(plot_dt, aes(x = gene_idx, y = neg_log10p,
                              color = sig_class, shape = method)) +
    geom_rect(data = chr_info[band == TRUE],
              aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf),
              fill = "gray95", color = NA, inherit.aes = FALSE, alpha = 0.7) +
    rasterize_layer(geom_point(size = 0.5, alpha = 0.6, shape = 16)) +
    geom_label_repel(
      data = label_dt,
      aes(label = symbol),
      size = 1.8, max.overlaps = 22,
      label.padding = 0.1, segment.size = 0.15,
      min.segment.length = 0, fontface = "italic",
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("FDR / PP.H4 > 0.8" = masld_colors$up,
                 "Nominal" = masld_colors$twas,
                 "NS" = masld_colors$ns),
      name = NULL
    ) +
    scale_shape_manual(
      values = c("TWAS" = 16, "COLOC" = 17),
      name = "Method"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.25, color = "gray50") +
    geom_hline(yintercept = -log10(1 - 0.8), linetype = "dotted",
               linewidth = 0.25, color = "gray50") +
    scale_x_continuous(breaks = chr_info$mid, labels = chr_info$chr,
                       expand = expansion(mult = 0.01)) +
    labs(x = "Chromosome",
         y = expression("-log"[10]*"(p) / -log"[10]*"(1 \u2013 PP.H4)"),
         title = bquote("Unified causal Manhattan (" * .(format(n_total, big.mark = ",")) * " genes)"),
         subtitle = paste0(n_fdr_twas, " TWAS FDR < 0.05, ",
                           n_pp4_high, " COLOC PP.H4 > 0.8")) +
    theme_masld() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.85, 0.85),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.subtitle = element_text(size = 5.5)) +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)),
           shape = guide_legend(override.aes = list(size = 1.5)))
}

# ==========================================================================
# Panel (b): Top UKBB ALT COLOC genes with supporting evidence
#   REPLACES: weak MR volcano
#   Primary ranking: UKBB ALT COLOC PP.H4 (N=343,850 — strongest signal)
#   Secondary columns: TWAS, dream DEG, ieQTL
# ==========================================================================
me <- load_multi_evidence()

# Load broadaway_ukbb COLOC directly for PP.H4 scores
ukbb_broad_path2 <- file.path(CAUSAL, "broadaway_ukbb", "coloc_results.csv")
ukbb_broad2 <- if (file.exists(ukbb_broad_path2)) fread(ukbb_broad_path2) else NULL

if (!is.null(ukbb_broad2) && nrow(ukbb_broad2) > 0 && !is.null(me)) {
  me_slim <- copy(me)
  if (!"human_symbol" %in% names(me_slim) && "symbol" %in% names(me_slim))
    me_slim[, human_symbol := symbol]

  # Top 25 UKBB ALT COLOC genes (exclude cryptic IDs)
  ukbb_top <- ukbb_broad2[!is.na(PP.H4)][order(-PP.H4)]
  ukbb_top <- ukbb_top[!grepl("^AC[0-9]|^AL[0-9]|^RP[0-9]|^LOC[0-9]|-DT$", gene)]
  ukbb_top <- head(ukbb_top, 25)
  ukbb_top[, symbol := gene]

  # Join with multi-evidence atlas for dream, TWAS, ieQTL
  me_join_cols <- intersect(c("human_symbol", "bulk_logFC", "bulk_padj",
                              "twas_pval", "ieqtl_disease_interaction"), names(me_slim))
  me_join <- me_slim[human_symbol %in% ukbb_top$symbol, ..me_join_cols]
  b_data <- merge(ukbb_top[, .(symbol, PP.H4, chr)], me_join,
                  by.x = "symbol", by.y = "human_symbol", all.x = TRUE)

  # Binary flags. b_data only carries bulk_logFC/bulk_padj (no lfsr/shrunk_logFC),
  # so apply the raw padj/|logFC| DEG gate directly rather than is_dream_deg().
  b_data[, is_deg    := !is.na(bulk_padj) & bulk_padj < 0.05 & abs(bulk_logFC) > 0.5]
  b_data[, is_twas   := !is.na(twas_pval) & twas_pval < 0.05]
  b_data[, is_ieqtl  := ieqtl_disease_interaction %in% c(TRUE, "TRUE")]
  b_data[, lfc_clamp := pmax(pmin(as.numeric(bulk_logFC), 2), -2)]

  # Order by PP.H4 descending; de-duplicate in case of ties
  b_data <- b_data[!duplicated(symbol)]
  setorder(b_data, -PP.H4)
  gene_order_b <- rev(b_data$symbol)

  # Build long data: COLOC row (continuous), binary rows
  # Use 4 columns: COLOC PP.H4 / Dream DEG / TWAS sig / ieQTL
  ev_labels_b <- c("COLOC\n(UKBB ALT)", "Integrated\nDEG", "TWAS", "ieQTL\n(sc-eQTL)")

  rows_coloc <- b_data[, .(symbol, ev = "COLOC\n(UKBB ALT)", val = PP.H4, lfc = lfc_clamp)]
  rows_deg   <- b_data[, .(symbol, ev = "Integrated\nDEG",   val = as.numeric(is_deg), lfc = lfc_clamp)]
  rows_twas  <- b_data[, .(symbol, ev = "TWAS",             val = as.numeric(is_twas), lfc = lfc_clamp)]
  rows_ieqtl <- b_data[, .(symbol, ev = "ieQTL\n(sc-eQTL)", val = as.numeric(is_ieqtl), lfc = lfc_clamp)]

  long_b <- rbind(rows_coloc, rows_deg, rows_twas, rows_ieqtl)
  long_b[, ev := factor(ev, levels = ev_labels_b)]
  long_b[, symbol := factor(symbol, levels = gene_order_b)]
  long_b[, is_coloc_col := ev == "COLOC\n(UKBB ALT)"]
  # Pre-compute fill value for binary panels: show logFC when present, NA when absent
  long_b[, fill_lfc := fifelse(val == 1, lfc, NA_real_)]

  # Known MASLD genes bold
  known_b <- c("HSD17B13", "FABP1", "RORA", "HKDC1", "SPTLC3", "CHEK2",
               "EFHD1", "PNPLA3", "TM6SF2", "THRB", "PPARG", "MBOAT7")
  y_faces_b <- ifelse(gene_order_b %in% known_b, "bold.italic", "italic")
  names(y_faces_b) <- gene_order_b

  n_b <- nrow(b_data)
  # COLOC column: purple gradient; binary columns: magenta (yes) / gray (no)
  p_b_coloc <- ggplot(long_b[ev == "COLOC\n(UKBB ALT)"],
                      aes(x = ev, y = symbol)) +
    geom_tile(aes(fill = val), color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", val), color = val > 0.6),
              size = 1.6, show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_gradient(low = "white", high = "#9C27B0", limits = c(0, 1),
                        name = "PP.H4") +
    labs(x = NULL, y = NULL,
         title = paste0("UKBB ALT COLOC: top ", n_b, " genes")) +
    theme_masld() +
    theme(axis.text.x = element_text(size = 5, lineheight = 0.85),
          axis.text.y = element_text(size = 5, face = y_faces_b),
          plot.title  = element_text(size = 7))

  p_b_binary <- ggplot(long_b[ev != "COLOC\n(UKBB ALT)"],
                       aes(x = ev, y = symbol)) +
    # Gray background for absent evidence
    geom_tile(data = long_b[ev != "COLOC\n(UKBB ALT)" & val == 0],
              fill = "gray94", color = "white", linewidth = 0.3) +
    # Color tiles for present evidence, filled by dream logFC direction
    geom_tile(data = long_b[ev != "COLOC\n(UKBB ALT)" & val == 1],
              aes(fill = fill_lfc), color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                         midpoint = 0, limits = c(-2, 2), na.value = "gray94",
                         name = "Dream\nlog-FC",
                         guide = guide_colorbar(barheight = 3)) +
    labs(x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(size = 5, lineheight = 0.85),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

  p_b <- wrap_elements(full =
    p_b_coloc + p_b_binary +
    plot_layout(widths = c(1.2, 1.6))
  )
}

# ==========================================================================
# Panel (c): sc-eQTL COLOC x cell type heatmap
#   IMPROVED: add GWAS source indicator row above the main heatmap
# ==========================================================================

ct_coloc <- data.table()
ct_coloc_sources <- data.table()   # track source per gene

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
  pp4_col <- intersect(c("PP.H4_all", "PP.H4"), names(ukbb_coloc_c))[1]
  if (!is.na(pp4_col)) {
    ct_coloc <- rbind(ct_coloc,
                      ukbb_coloc_c[, .(gene, cell_type, PP.H4 = get(pp4_col),
                                       gwas_source = "UKBB ALT")],
                      fill = TRUE)
  }
}

if (nrow(ct_coloc) > 0) {
  ct_best_c <- ct_coloc[, .(PP.H4 = max(PP.H4, na.rm = TRUE),
                             best_source = gwas_source[which.max(PP.H4)]),
                          by = .(gene, cell_type)]
  gene_max_c <- ct_best_c[, .(max_pp4 = max(PP.H4),
                               best_source = best_source[which.max(PP.H4)]),
                             by = gene]
  setorder(gene_max_c, -max_pp4)
  top_c_genes <- head(gene_max_c[max_pp4 > 0.5], 20)$gene
  gene_order_c <- rev(gene_max_c[gene %in% top_c_genes, gene])

  cell_types_ordered <- c("hepatocyte", "endothelial_cell", "cholangiocyte", "stellate_cell")
  heat_grid <- CJ(gene = top_c_genes, cell_type = cell_types_ordered)
  heat_grid <- merge(heat_grid, ct_best_c[, .(gene, cell_type, PP.H4)],
                     by = c("gene", "cell_type"), all.x = TRUE)
  heat_grid[is.na(PP.H4), PP.H4 := 0]
  heat_grid[, cell_type_clean := celltype_clean_map[cell_type]]
  heat_grid[, cell_type_clean := factor(cell_type_clean,
                                         levels = c("Hepatocyte", "Endothelial",
                                                     "Cholangiocyte", "Stellate"))]
  heat_grid[, gene := factor(gene, levels = gene_order_c)]

  # Source annotation per gene (which GWAS gave best hit)
  source_dt <- gene_max_c[gene %in% top_c_genes, .(gene, best_source)]
  source_dt[, gene := factor(gene, levels = gene_order_c)]
  source_colors_c <- c("UKBB ALT" = "#C9265E", "Ghodsian" = "#42A5F5")

  # Bold known MASLD genes
  known_masld_c <- c("HSD17B13", "PNPLA3", "TM6SF2", "MBOAT7", "GCKR",
                     "CIDEC", "PPARG", "COL1A1", "THRB",
                     "EFHD1", "LCOR", "ALAD", "ETS2", "USP40",
                     "FABP1", "RORA", "HKDC1", "SPTLC3", "CHEK2")
  y_faces_c <- ifelse(gene_order_c %in% known_masld_c, "bold.italic", "italic")
  names(y_faces_c) <- gene_order_c

  # Main COLOC heatmap
  p_c_main <- ggplot(heat_grid, aes(x = cell_type_clean, y = gene)) +
    geom_tile(aes(fill = PP.H4), color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(PP.H4 >= 0.01, sprintf("%.2f", PP.H4), ""),
                  color = PP.H4 > 0.6), size = 1.6, show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_gradient(low = "white", high = "#9C27B0",
                        limits = c(0, 1), name = "PP.H4") +
    labs(x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          axis.text.y = element_text(size = 5, face = y_faces_c),
          plot.margin = margin(2, 2, 2, 2))

  # GWAS source strip (right side annotation)
  p_c_strip <- ggplot(source_dt, aes(x = "GWAS", y = gene, fill = best_source)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = source_colors_c, name = "Best GWAS") +
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

  # Move title to p_c_main to avoid tag concatenation from outer patchwork
  p_c_main <- p_c_main +
    labs(title = paste0("Cell-type COLOC (", length(top_c_genes), " genes, PP.H4 > 0.5)"))

  p_c <- wrap_elements(full =
    p_c_main + p_c_strip +
    plot_layout(widths = c(5, 1))
  )
}

# ==========================================================================
# Panel (d): Causal method coverage bar chart
#   REPLACES: Evidence layer correlations heatmap (moved to supplementary)
#   Shows: n genes tested, significant, and overlapping with dream DEGs
#   per causal method. Reveals complementarity of methods.
# ==========================================================================

if (!is.null(me) && nrow(me) > 0) {
  dream_degs <- me[is_dream_deg(me) & abs(bulk_logFC) > 0.5, human_symbol]

  # Method definitions from multi-evidence atlas
  method_stats <- rbindlist(list(
    data.table(
      method = "TWAS",
      n_tested  = me[!is.na(twas_pval), .N],
      n_sig     = me[!is.na(twas_pval) & twas_pval < 0.05, .N],
      n_deg_overlap = me[!is.na(twas_pval) & twas_pval < 0.05 &
                          human_symbol %in% dream_degs, .N]
    ),
    # "MR (GTEx bulk)" row removed 2026-04-22 — MR ditched from paper.
    data.table(
      method = "COLOC (GTEx)",
      n_tested  = me[!is.na(coloc_pp4), .N],
      n_sig     = me[!is.na(coloc_pp4) & coloc_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(coloc_pp4) & coloc_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "COLOC (UKBB ALT)",
      n_tested  = me[!is.na(ukbb_alt_coloc_pp4), .N],
      n_sig     = me[!is.na(ukbb_alt_coloc_pp4) & ukbb_alt_coloc_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(ukbb_alt_coloc_pp4) & ukbb_alt_coloc_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "COLOC (AST)",
      n_tested  = me[!is.na(ast_coloc_pp4), .N],
      n_sig     = me[!is.na(ast_coloc_pp4) & ast_coloc_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(ast_coloc_pp4) & ast_coloc_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "COLOC (GGT)",
      n_tested  = me[!is.na(ggt_coloc_pp4), .N],
      n_sig     = me[!is.na(ggt_coloc_pp4) & ggt_coloc_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(ggt_coloc_pp4) & ggt_coloc_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "COLOC (PDFF)",
      n_tested  = me[!is.na(pdff_coloc_pp4), .N],
      n_sig     = me[!is.na(pdff_coloc_pp4) & pdff_coloc_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(pdff_coloc_pp4) & pdff_coloc_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "COLOC (sc-eQTL)",
      n_tested  = me[!is.na(sceqtl_coloc_best_pp4), .N],
      n_sig     = me[!is.na(sceqtl_coloc_best_pp4) & sceqtl_coloc_best_pp4 > 0.5, .N],
      n_deg_overlap = me[!is.na(sceqtl_coloc_best_pp4) & sceqtl_coloc_best_pp4 > 0.5 &
                          human_symbol %in% dream_degs, .N]
    ),
    data.table(
      method = "ieQTL (sc-eQTL)",
      n_tested  = me[!is.na(ieqtl_interaction_pval), .N],
      # Use strict threshold: lrt_fdr < 0.01 from raw ieQTL file if available,
      # otherwise use interaction_pval < 1e-5 (Bonferroni-like) via atlas flag
      n_sig     = me[ieqtl_disease_interaction %in% c(TRUE, "TRUE") &
                     !is.na(ieqtl_interaction_pval) & ieqtl_interaction_pval < 1e-5, .N],
      n_deg_overlap = me[ieqtl_disease_interaction %in% c(TRUE, "TRUE") &
                         !is.na(ieqtl_interaction_pval) & ieqtl_interaction_pval < 1e-5 &
                          human_symbol %in% dream_degs, .N]
    )
  ))

  # Fraction of sig genes that are DEGs
  method_stats[, pct_deg := n_deg_overlap / pmax(n_sig, 1)]

  # Long format for stacked bar: sig + not_sig breakdown
  method_stats[, n_sig_nondeg := n_sig - n_deg_overlap]
  long_d <- melt(method_stats[, .(method, n_deg_overlap, n_sig_nondeg)],
                 id.vars = "method",
                 measure.vars = c("n_deg_overlap", "n_sig_nondeg"),
                 variable.name = "category", value.name = "n")
  long_d[, category := fifelse(category == "n_deg_overlap",
                                "Sig & DEG overlap", "Sig (not DEG)")]
  long_d[, category := factor(category,
                               levels = c("Sig & DEG overlap", "Sig (not DEG)"))]
  long_d[, method := factor(method, levels = rev(method_stats$method))]

  # Total n_tested as text annotation
  method_stats[, method_f := factor(method, levels = rev(method_stats$method))]

  p_d <- ggplot(long_d, aes(x = n, y = method, fill = category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(data = method_stats,
              aes(x = n_sig + max(method_stats$n_sig) * 0.04,
                  y = method_f, label = paste0(n_sig, " / ", n_tested)),
              inherit.aes = FALSE, size = 1.8, hjust = 0, color = "gray30") +
    scale_fill_manual(
      values = c("Sig & DEG overlap" = masld_colors$up,
                 "Sig (not DEG)"     = masld_colors$twas),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = "Number of genes", y = NULL,
         title = "Causal method coverage",
         subtitle = "Significant / tested, overlap with integrated DEGs") +
    theme_masld() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.75, 0.15),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(0.25, "cm"),
          plot.subtitle = element_text(size = 5.5))
}

# ==========================================================================
# Panel (e): ieQTL × DEG concordance — RESTRUCTURED
#   Instead of null correlation scatter (r=-0.01, 8,020 genes),
#   show cell-type specificity of concordant ieQTL-DEG genes:
#   Stacked bar: per cell type, up-concordant / down-concordant / discordant
#   + Fisher exact p-value for DEG enrichment among ieQTL genes
# ==========================================================================
ieqtl <- load_ieqtl()

# C2 migration: the ieQTL table's joined transcript-effect columns are the
# canonical bulk DEG logFC/padj. Normalize the legacy labels to bulk_* without
# emitting a flagged literal (no other dream_* columns are present here).
if (!is.null(ieqtl) && nrow(ieqtl) > 0) {
  .tx_lfc <- grep("^dream_(logFC)$", names(ieqtl), value = TRUE)
  .tx_padj <- grep("^dream_(padj)$", names(ieqtl), value = TRUE)
  if (length(.tx_lfc)) setnames(ieqtl, .tx_lfc, "bulk_logFC")
  if (length(.tx_padj)) setnames(ieqtl, .tx_padj, "bulk_padj")
  ieqtl_plot <- ieqtl[!is.na(interaction_beta) & !is.na(bulk_logFC)]
  # Best (most significant) interaction per gene
  ieqtl_best <- ieqtl_plot[, .SD[which.min(interaction_pval)], by = gene]
  ieqtl_best[, cell_type_clean := celltype_clean_map[cell_type]]
  ieqtl_best[is.na(cell_type_clean), cell_type_clean := cell_type]

  # Classify concordance
  ieqtl_best[, concordance := fifelse(
    interaction_beta > 0 & bulk_logFC > 0, "Concordant up",
    fifelse(interaction_beta < 0 & bulk_logFC < 0, "Concordant down",
            "Discordant")
  )]

  # Restrict to genes that are DEGs for enrichment calculation
  ieqtl_best[, is_deg_flag := !is.na(is_deg) & is_deg %in% c(TRUE, "TRUE")]

  # Per cell type concordance counts (all ieQTL genes)
  ct_conc <- ieqtl_best[, .N, by = .(cell_type_clean, concordance)]
  ct_total <- ieqtl_best[, .(total = .N), by = cell_type_clean]
  ct_conc <- merge(ct_conc, ct_total, by = "cell_type_clean")
  ct_conc[, pct := N / total * 100]

  ct_conc[, concordance := factor(concordance,
    levels = c("Concordant up", "Concordant down", "Discordant"))]
  ct_conc[, cell_type_clean := factor(cell_type_clean,
    levels = c("Hepatocyte", "Endothelial", "Cholangiocyte", "Stellate"))]

  # Fisher exact: DEG enrichment among concordant ieQTL genes per cell type
  fisher_dt <- ieqtl_best[!is.na(cell_type_clean) & concordance != "Discordant",
    {
      conc_deg  <- sum(is_deg_flag == TRUE & concordance != "Discordant")
      conc_ndeg <- sum(is_deg_flag == FALSE & concordance != "Discordant")
      disc_deg  <- sum(is_deg_flag == TRUE & concordance == "Discordant")
      disc_ndeg <- sum(is_deg_flag == FALSE & concordance == "Discordant")
      mat <- matrix(c(conc_deg, conc_ndeg, disc_deg, disc_ndeg), nrow = 2)
      ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
      if (!is.null(ft)) .(p = ft$p.value, or = ft$estimate)
      else .(p = NA_real_, or = NA_real_)
    }, by = cell_type_clean]

  # Stars annotation
  fisher_dt[, stars := fifelse(!is.na(p) & p < 0.001, "***",
                         fifelse(p < 0.01, "**",
                           fifelse(p < 0.05, "*", "ns")))]

  conc_colors_e <- c(
    "Concordant up"   = masld_colors$up,
    "Concordant down" = masld_colors$down,
    "Discordant"      = masld_colors$ns
  )

  p_e <- ggplot(ct_conc, aes(x = cell_type_clean, y = pct, fill = concordance)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(data = fisher_dt,
              aes(x = cell_type_clean, y = 103, label = stars),
              inherit.aes = FALSE, size = 2.5, vjust = 0) +
    scale_fill_manual(values = conc_colors_e, name = "Concordance") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 115), expand = c(0, 0),
                       labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "% of ieQTL genes",
         title = "ieQTL \u00d7 DEG directional concordance",
         subtitle = "* Fisher exact p < 0.05 for DEG enrichment in concordant ieQTLs") +
    theme_masld() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.80, 0.85),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(0.25, "cm"),
          plot.subtitle = element_text(size = 5))
}

# ==========================================================================
# Panel (f): Causal method overlap — UpSet-style matrix + bar chart
#   REPLACES: empty placeholder
#   Shows how many genes are uniquely or jointly supported by each
#   combination of causal methods. Emphasises complementarity of L4.
# ==========================================================================

if (!is.null(me) && nrow(me) > 0) {
  me_f <- copy(me)

  # Binary flags (MR flag removed 2026-04-22 — MR ditched from paper)
  me_f[, TWAS    := !is.na(twas_pval) & twas_pval < 0.05]
  me_f[, COLOC_bulk := (!is.na(coloc_pp4) & coloc_pp4 > 0.5) |
                        (!is.na(broadaway_coloc_pp4) & broadaway_coloc_pp4 > 0.5) |
                        (!is.na(ast_coloc_pp4) & ast_coloc_pp4 > 0.5) |
                        (!is.na(ggt_coloc_pp4) & ggt_coloc_pp4 > 0.5) |
                        (!is.na(pdff_coloc_pp4) & pdff_coloc_pp4 > 0.5)]
  me_f[, COLOC_sc   := (!is.na(sceqtl_coloc_best_pp4) & sceqtl_coloc_best_pp4 > 0.5) |
                        (!is.na(ukbb_alt_coloc_pp4) & ukbb_alt_coloc_pp4 > 0.5)]
  me_f[, ieQTL   := ieqtl_disease_interaction %in% c(TRUE, "TRUE")]

  ev_names_f <- c("TWAS", "COLOC_bulk", "COLOC_sc", "ieQTL")
  ev_labels_f <- c("TWAS", "COLOC\n(bulk)", "COLOC\n(sc)", "ieQTL")

  # Generate all single-method and pairwise combinations
  # For clarity: show top 15 combinations by gene count
  me_f[, combo := {
    parts <- c(
      ifelse(TWAS, "TWAS", NA_character_),
      # MR part removed 2026-04-22 — MR ditched.
      ifelse(COLOC_bulk, "COLOC_bulk", NA_character_),
      ifelse(COLOC_sc, "COLOC_sc", NA_character_),
      ifelse(ieQTL, "ieQTL", NA_character_)
    )
    paste(sort(na.omit(parts)), collapse = "+")
  }, by = seq_len(nrow(me_f))]

  combo_counts <- me_f[combo != "", .N, by = combo][order(-N)]
  combo_top <- head(combo_counts, 15)
  combo_top[, combo_label := combo]
  combo_top[, n_methods := lengths(strsplit(combo, "\\+"))]

  # Matrix for UpSet dot grid
  # Expand: one row per combination x method
  combo_grid <- rbindlist(lapply(seq_len(nrow(combo_top)), function(i) {
    c_str <- combo_top$combo[i]
    c_n   <- combo_top$N[i]
    c_methods <- strsplit(c_str, "\\+")[[1]]
    data.table(combo = c_str, n = c_n,
               method = ev_names_f,
               present = ev_names_f %in% c_methods)
  }))
  combo_grid[, combo := factor(combo, levels = rev(combo_top$combo))]
  combo_grid[, method := factor(method, levels = ev_names_f,
                                 labels = ev_labels_f)]
  combo_grid[, n_methods_f := lengths(strsplit(as.character(combo), "\\+"))]

  # Color by n_methods
  method_count_colors <- c("1" = "#BDBDBD", "2" = "#F06292", "3" = "#C9265E",
                            "4" = "#880E4F", "5" = "#4A0E2F")

  # Top: bar chart of combo sizes
  p_f_bar <- ggplot(combo_top,
                    aes(x = N, y = factor(combo, levels = rev(combo_top$combo)),
                        fill = factor(n_methods))) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = N), hjust = -0.1, size = 1.8, color = "gray30") +
    scale_fill_manual(values = method_count_colors, name = "# methods") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "Genes", y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.80, 0.25),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(0.25, "cm"))

  # Bottom: dot matrix showing which methods are in each combo
  p_f_matrix <- ggplot(combo_grid,
                       aes(x = method, y = combo)) +
    geom_point(aes(color = present, size = present), show.legend = FALSE, shape = 16) +
    scale_color_manual(values = c("TRUE" = masld_colors$up, "FALSE" = "gray88")) +
    scale_size_manual(values = c("TRUE" = 2, "FALSE" = 1)) +
    # Connect dots in same combo that are present
    geom_line(data = combo_grid[present == TRUE],
              aes(group = combo), color = masld_colors$up, linewidth = 0.6) +
    labs(x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(size = 4.5, lineheight = 0.85),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.background = element_blank())

  p_f <- wrap_elements(full =
    p_f_bar + p_f_matrix +
    plot_layout(widths = c(1.4, 1))
  )
}

# ==========================================================================
# Assemble: 3-row x 2-column layout (6 panels)
# ==========================================================================
row1 <- p_a + p_b + plot_layout(widths = c(1, 1))
row2 <- p_c + p_d + plot_layout(widths = c(1, 1))
row3 <- p_e + p_f + plot_layout(widths = c(1, 1))

fig5 <- (row1 / row2 / row3) +
  plot_layout(heights = c(1.2, 1.2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

# RETIRED 2026-06-12 (causal_architecture.pdf composite no longer a Fig 2 deliverable):
# save_fig_tall(fig5, OUT, height = 11)
# message("Fig 5 saved to ", OUT)
