#!/usr/bin/env Rscript
# =============================================================================
# figS07_ncrna.R
# Supplementary Figure 7: Non-coding RNA Landscape of the MASLD Transcriptome
#
# 8 panels characterizing ncRNA biology in the MASLD atlas:
# (a) Biotype breakdown in DEGs
# (b) lncRNA vs protein-coding effect size
# (c) Top MASLD ceRNA subnetwork
# (d) Cell-type-specific lncRNA heatmap
# (e) lncRNA synteny conservation
# (f) ATAC accessibility at ncRNA vs PC promoters
# (g) Known MASLD lncRNA validation (forest plot)
# (h) ncRNA multi-evidence upset
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Input paths
ncrna_dir    <- file.path(BASE, "RNA-seq/results/ncrna")
atlas_path   <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
sc_dir       <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/elatus")

# Output
out_dir <- FIGS07_DIR
dir.create(file.path(out_dir, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("=== figS_ncrna: Non-coding RNA Landscape ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# Load data
atlas <- fread(atlas_path)

# ncRNA color palette
ncrna_colors <- c(
  lncRNA   = "#7B1FA2",  # Violet
  miRNA    = "#E91E63",  # Pink
  snoRNA   = "#00695C",  # Teal
  snRNA    = "#42A5F5",  # Light blue
  misc_RNA = "#F48FB1",  # Soft pink
  scaRNA   = "#78909C",  # Blue-gray
  protein_coding = "#BDBDBD"  # Gray
)

# =============================================================================
# Panel (a): ncRNA Biotype Breakdown in DEGs
# =============================================================================
cat("--- Panel (a): Biotype breakdown ---\n")

landscape_path <- file.path(ncrna_dir, "ncrna_landscape_summary.csv")
if (file.exists(landscape_path)) {
  landscape <- fread(landscape_path)
  landscape <- landscape[gene_biotype != "protein_coding"]

  # Prepare data for grouped bar plot
  bar_data <- melt(landscape[, .(gene_biotype, n_tested, n_deg)],
                    id.vars = "gene_biotype",
                    variable.name = "category",
                    value.name = "count")
  bar_data[, category := fifelse(category == "n_tested", "Tested", "DEG (padj<0.05)")]
  bar_data[, gene_biotype := factor(gene_biotype,
    levels = landscape[order(-n_total)]$gene_biotype)]

  pa <- ggplot(bar_data, aes(x = gene_biotype, y = count, fill = category)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Tested" = "#BDBDBD", "DEG (padj<0.05)" = "#7B1FA2")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "Number of genes", fill = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = c(0.7, 0.85))

  # Add fraction labels
  frac_labels <- landscape[, .(gene_biotype, label = sprintf("%.0f%%", 100 * frac_deg),
                                y = n_tested)]
  pa <- pa + geom_text(data = frac_labels, aes(x = gene_biotype, y = y, label = label),
                        inherit.aes = FALSE, size = 2, vjust = -0.3)

  save_fig(pa, file.path(out_dir, "panels", "panel_a.pdf"), width = fig_half_width, height = 3)
  cat("  Saved panel_a.pdf\n")
} else {
  cat("  WARNING: landscape_summary.csv not found\n")
  pa <- placeholder("Panel (a)\nBiotype breakdown")
}

# =============================================================================
# Panel (b): lncRNA vs Protein-Coding Effect Size
# =============================================================================
cat("--- Panel (b): lncRNA vs PC effect size ---\n")

atlas[, is_deg := is_dream_deg(atlas)]
degs <- atlas[is_deg == TRUE & gene_biotype %in% c("protein_coding", "lncRNA")]
degs[, biotype_label := fifelse(gene_biotype == "protein_coding", "Protein-coding", "lncRNA")]

if (nrow(degs) > 10) {
  pb <- ggplot(degs, aes(x = biotype_label, y = abs(dream_logFC), fill = biotype_label)) +
    geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.3) +
    scale_fill_manual(values = c("Protein-coding" = "#BDBDBD", "lncRNA" = "#7B1FA2")) +
    scale_y_continuous(limits = c(0, quantile(abs(degs$dream_logFC), 0.99))) +
    labs(x = NULL, y = "|log2 FC|", fill = NULL) +
    theme_masld() +
    theme(legend.position = "none")

  # Add Wilcoxon p-value
  wt <- wilcox.test(abs(degs[gene_biotype == "lncRNA"]$dream_logFC),
                     abs(degs[gene_biotype == "protein_coding"]$dream_logFC))
  pb <- pb + annotate("text", x = 1.5, y = quantile(abs(degs$dream_logFC), 0.97),
                       label = sprintf("p = %.2e", wt$p.value), size = 2.5)

  save_fig(pb, file.path(out_dir, "panels", "panel_b.pdf"), width = fig_half_width * 0.6, height = 3)
  cat("  Saved panel_b.pdf\n")
} else {
  pb <- placeholder("Panel (b)\nEffect size comparison")
}

# =============================================================================
# Panel (c): Top MASLD ceRNA Subnetwork
# =============================================================================
cat("--- Panel (c): ceRNA network ---\n")

cerna_hubs_path <- file.path(ncrna_dir, "cerna_hubs.csv")
cerna_network_path <- file.path(ncrna_dir, "cerna_network_masld.csv")

if (file.exists(cerna_hubs_path) && file.exists(cerna_network_path)) {
  hubs <- fread(cerna_hubs_path)
  cerna <- fread(cerna_network_path)

  if (nrow(hubs) > 0 && requireNamespace("igraph", quietly = TRUE)) {
    library(igraph)

    # Take top 30 hub lncRNAs + their connections
    top_lnc <- head(hubs[node_type == "lncRNA"][order(-degree)], 30)$node
    sub_cerna <- cerna[lncrna %in% top_lnc]

    if (nrow(sub_cerna) > 0) {
      # Build edges for plotting
      edges <- rbind(
        data.table(from = sub_cerna$lncrna, to = paste0("miR:", sub_cerna$shared_mirnas)),
        data.table(from = paste0("miR:", sub_cerna$shared_mirnas), to = sub_cerna$mrna)
      )
      # Simplify: use direct lncRNA-mRNA edges
      edges_simple <- data.table(from = sub_cerna$lncrna, to = sub_cerna$mrna)
      edges_simple <- unique(edges_simple)

      if (nrow(edges_simple) > 0) {
        g <- graph_from_data_frame(head(edges_simple, 200), directed = FALSE)
        V(g)$type <- ifelse(V(g)$name %in% top_lnc, "lncRNA", "mRNA")
        V(g)$color <- ifelse(V(g)$type == "lncRNA", "#7B1FA2", "#42A5F5")

        # Plot to PDF directly
        pdf(file.path(out_dir, "panels", "panel_c.pdf"), width = fig_half_width, height = fig_half_width)
        par(mar = c(0.5, 0.5, 0.5, 0.5))
        plot(g, vertex.size = ifelse(V(g)$type == "lncRNA", 6, 3),
             vertex.color = V(g)$color,
             vertex.label = ifelse(V(g)$type == "lncRNA", V(g)$name, NA),
             vertex.label.cex = 0.4, vertex.label.color = "black",
             vertex.frame.color = NA,
             edge.color = alpha("#BDBDBD", 0.3), edge.width = 0.3,
             layout = layout_with_fr(g))
        legend("bottomleft", legend = c("lncRNA hub", "mRNA target"),
               col = c("#7B1FA2", "#42A5F5"), pch = 16, cex = 0.6, bty = "n")
        dev.off()
        cat("  Saved panel_c.pdf\n")
      }
    }
  }
} else {
  cat("  WARNING: ceRNA network files not found\n")
  pc <- placeholder("Panel (c)\nceRNA network")
  save_fig(pc, file.path(out_dir, "panels", "panel_c.pdf"), width = fig_half_width, height = fig_half_width)
}

# =============================================================================
# Panel (d): Cell-Type-Specific lncRNA Heatmap
# =============================================================================
cat("--- Panel (d): Cell-type lncRNA heatmap ---\n")

ct_expr_path <- file.path(sc_dir, "lncrna_celltype_expression.csv")
ct_spec_path <- file.path(sc_dir, "lncrna_celltype_specific.csv")

if (file.exists(ct_spec_path) && file.exists(ct_expr_path) &&
    requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  library(ComplexHeatmap)
  library(circlize)

  ct_spec <- fread(ct_spec_path)
  ct_expr <- fread(ct_expr_path)

  # Try to find gene column
  gene_col <- intersect(names(ct_spec), c("gene", "human_symbol", "gene_name"))
  tau_col <- intersect(names(ct_spec), c("tau", "specificity_index"))

  if (length(gene_col) > 0 && length(tau_col) > 0) {
    setnames(ct_spec, gene_col[1], "gene")
    setnames(ct_spec, tau_col[1], "tau")

    # Top 30 most specific lncRNAs
    top_specific <- head(ct_spec[order(-tau)], 30)

    # Build expression matrix from ct_expr
    ct_col <- intersect(names(ct_expr), c("cell_type", "celltype"))
    expr_col <- intersect(names(ct_expr), c("mean_expr", "mean_expression", "detection_rate"))

    if (length(ct_col) > 0 && length(expr_col) > 0) {
      gene_col2 <- intersect(names(ct_expr), c("gene", "human_symbol"))
      if (length(gene_col2) > 0) {
        setnames(ct_expr, gene_col2[1], "gene")
        setnames(ct_expr, ct_col[1], "cell_type")
        setnames(ct_expr, expr_col[1], "value")

        expr_mat_dt <- ct_expr[gene %in% top_specific$gene]
        # Deduplicate before dcast
        expr_mat_dt <- expr_mat_dt[, .(value = mean(value, na.rm = TRUE)), by = .(gene, cell_type)]
        expr_mat <- dcast(expr_mat_dt, gene ~ cell_type, value.var = "value", fill = 0)
        gene_names <- expr_mat$gene
        expr_mat <- as.matrix(expr_mat[, -1])
        rownames(expr_mat) <- gene_names

        # Z-score normalize rows
        expr_mat_z <- t(scale(t(expr_mat)))
        expr_mat_z[is.nan(expr_mat_z)] <- 0

        # Heatmap
        col_fun <- colorRamp2(c(-2, 0, 2), c("#1565C0", "white", "#C2185B"))
        ht <- Heatmap(expr_mat_z,
                       name = "Z-score",
                       col = col_fun,
                       row_names_gp = gpar(fontsize = 5),
                       column_names_gp = gpar(fontsize = 6),
                       column_names_rot = 45,
                       show_row_dend = TRUE,
                       show_column_dend = TRUE,
                       width = unit(fig_half_width * 0.8, "inches"),
                       height = unit(4, "inches"))

        pdf(file.path(out_dir, "panels", "panel_d.pdf"), width = fig_half_width, height = 4.5)
        draw(ht)
        dev.off()
        cat("  Saved panel_d.pdf\n")
      }
    }
  }
} else {
  cat("  WARNING: Cell-type expression files not found or ComplexHeatmap unavailable\n")
  pd <- placeholder("Panel (d)\nCell-type lncRNA heatmap")
  save_fig(pd, file.path(out_dir, "panels", "panel_d.pdf"), width = fig_half_width, height = 4.5)
}

# =============================================================================
# Panel (e): lncRNA Synteny Conservation
# =============================================================================
cat("--- Panel (e): Synteny conservation ---\n")

synteny_path <- file.path(ncrna_dir, "lncrna_synteny_conservation.csv")
if (file.exists(synteny_path)) {
  synteny <- fread(synteny_path)

  if ("synteny_status" %in% names(synteny)) {
    # Stacked bar: synteny status breakdown
    syn_counts <- synteny[, .N, by = synteny_status]
    syn_counts[, frac := N / sum(N)]
    syn_counts[, synteny_status := factor(synteny_status,
      levels = c("Synteny_conserved", "Synteny_partial", "Synteny_absent",
                 "No_flanking_ortholog", "No_flanking_match", "Not_assessed"))]

    syn_colors <- c(
      Synteny_conserved = "#00695C",
      Synteny_partial = "#42A5F5",
      Synteny_absent = "#F48FB1",
      No_flanking_ortholog = "#BDBDBD",
      No_flanking_match = "#E0E0E0",
      Not_assessed = "#F5F5F5"
    )

    pe_bar <- ggplot(syn_counts[!is.na(synteny_status)],
                      aes(x = "lncRNAs", y = N, fill = synteny_status)) +
      geom_col(width = 0.5) +
      scale_fill_manual(values = syn_colors, name = "Conservation") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = NULL, y = "Number of lncRNAs") +
      theme_masld() +
      theme(legend.position = "right", legend.text = element_text(size = 6),
            legend.key.size = unit(0.3, "cm"))

    # PhastCons boxplot: DEG vs non-DEG
    if ("phastcons_tss" %in% names(synteny) && "is_deg" %in% names(synteny)) {
      phast_data <- synteny[!is.na(phastcons_tss) & !is.na(is_deg)]
      phast_data[, deg_status := fifelse(is_deg, "DEG", "Non-DEG")]

      pe_phast <- ggplot(phast_data, aes(x = deg_status, y = phastcons_tss, fill = deg_status)) +
        geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
        scale_fill_manual(values = c("DEG" = "#7B1FA2", "Non-DEG" = "#BDBDBD")) +
        labs(x = NULL, y = "PhastCons (TSS ± 500bp)") +
        theme_masld() +
        theme(legend.position = "none")

      pe <- pe_bar + pe_phast + plot_layout(widths = c(1, 1))
    } else {
      pe <- pe_bar
    }

    save_fig(pe, file.path(out_dir, "panels", "panel_e.pdf"), width = fig_half_width, height = 3)
    cat("  Saved panel_e.pdf\n")
  }
} else {
  cat("  WARNING: synteny file not found\n")
  pe <- placeholder("Panel (e)\nSynteny conservation")
  save_fig(pe, file.path(out_dir, "panels", "panel_e.pdf"), width = fig_half_width, height = 3)
}

# =============================================================================
# Panel (f): ATAC Accessibility ncRNA vs PC
# =============================================================================
cat("--- Panel (f): ATAC accessibility ---\n")

atac_comp_path <- file.path(ncrna_dir, "ncrna_atac_comparison.csv")
if (file.exists(atac_comp_path)) {
  atac_comp <- fread(atac_comp_path)

  # Paired bar chart
  atac_long <- melt(atac_comp, id.vars = "metric",
                     variable.name = "biotype", value.name = "value")
  atac_long[, biotype := fifelse(biotype == "ncRNA", "ncRNA", "Protein-coding")]
  atac_long[, metric := gsub("_frac$", "", metric)]
  atac_long[, metric := gsub("_", " ", metric)]

  pf <- ggplot(atac_long, aes(x = metric, y = value, fill = biotype)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("ncRNA" = "#7B1FA2", "Protein-coding" = "#BDBDBD")) +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "Fraction", fill = NULL) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          legend.position = c(0.7, 0.85))

  save_fig(pf, file.path(out_dir, "panels", "panel_f.pdf"), width = fig_half_width, height = 3)
  cat("  Saved panel_f.pdf\n")
} else {
  cat("  WARNING: atac comparison file not found\n")
  pf <- placeholder("Panel (f)\nATAC accessibility")
  save_fig(pf, file.path(out_dir, "panels", "panel_f.pdf"), width = fig_half_width, height = 3)
}

# =============================================================================
# Panel (g): Known MASLD lncRNA Validation (Forest Plot)
# =============================================================================
cat("--- Panel (g): Known MASLD lncRNA forest plot ---\n")

known_path <- file.path(ncrna_dir, "known_masld_lncrna_validation.csv")
if (file.exists(known_path)) {
  known <- fread(known_path)
  known <- known[in_atlas == TRUE & !is.na(dream_logFC)]

  if (nrow(known) > 0) {
    # Compute approximate 95% CI from t-statistic
    # SE = logFC / t; CI = logFC +/- 1.96*SE
    known[, se := abs(dream_logFC / dream_tstat)]
    known[, ci_lo := dream_logFC - 1.96 * se]
    known[, ci_hi := dream_logFC + 1.96 * se]
    known[, sig_label := fifelse(is_dream_deg(known), "*", "")]
    known[, concordant := (dream_direction == literature_direction)]

    # Deduplicate genes (keep first occurrence = highest |logFC| per direction)
    known <- known[!duplicated(gene)]
    # Order by literature direction then LFC
    setorder(known, literature_direction, -dream_logFC)
    known[, gene := factor(gene, levels = rev(unique(gene)))]

    pg <- ggplot(known, aes(x = dream_logFC, y = gene)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.3) +
      geom_point(aes(color = concordant, shape = is_deg), size = 2) +
      scale_color_manual(values = c("TRUE" = "#00695C", "FALSE" = "#E91E63"),
                          labels = c("Discordant", "Concordant"),
                          name = "Direction") +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                          labels = c("NS", "DEG"),
                          name = "Significance") +
      labs(x = "Integrated log2 FC", y = NULL) +
      theme_masld() +
      theme(legend.position = "right", legend.text = element_text(size = 6))

    save_fig(pg, file.path(out_dir, "panels", "panel_g.pdf"), width = fig_half_width, height = 3)
    cat("  Saved panel_g.pdf\n")
  }
} else {
  cat("  WARNING: known MASLD lncRNA validation file not found\n")
  pg <- placeholder("Panel (g)\nKnown MASLD lncRNA validation")
  save_fig(pg, file.path(out_dir, "panels", "panel_g.pdf"), width = fig_half_width, height = 3)
}

# =============================================================================
# Panel (h): ncRNA Multi-Evidence UpSet
# =============================================================================
cat("--- Panel (h): Multi-evidence UpSet ---\n")

ncrna_biotypes <- c("lncRNA", "miRNA", "snoRNA", "snRNA", "misc_RNA", "scaRNA")
ncrna_atlas <- atlas[gene_biotype %in% ncrna_biotypes & is_deg == TRUE]

if (nrow(ncrna_atlas) > 0 && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  library(ComplexHeatmap)

  # Build binary evidence matrix
  evidence_mat <- data.table(
    gene = ncrna_atlas$human_symbol,
    S1_Human_Bulk = is_dream_deg(ncrna_atlas),
    S2_Mouse = !is.na(ncrna_atlas$mouse_meta_padj) & ncrna_atlas$mouse_meta_padj < 0.1,
    # mr_sig term removed 2026-04-22 — MR ditched from paper.
    S3_Genetic = (!is.na(ncrna_atlas$twas_pval) & ncrna_atlas$twas_pval < 0.05) |
                  (!is.na(ncrna_atlas$coloc_pp4) & ncrna_atlas$coloc_pp4 > 0.5),
    S5_Epigenomic = !is.na(ncrna_atlas$hepatocyte_da_padj) & ncrna_atlas$hepatocyte_da_padj < 0.05,
    S7_SingleCell = !is.na(ncrna_atlas$sc_hepatocyte_padj) & ncrna_atlas$sc_hepatocyte_padj < 0.05
  )

  # Filter to genes with at least one evidence source
  evidence_mat[, n_sources := rowSums(.SD, na.rm = TRUE), .SDcols = patterns("^S")]
  evidence_mat <- evidence_mat[n_sources >= 1]

  if (nrow(evidence_mat) > 5) {
    # Build combination matrix
    mat <- as.matrix(evidence_mat[, .(S1_Human_Bulk, S2_Mouse, S3_Genetic,
                                       S5_Epigenomic, S7_SingleCell)])
    rownames(mat) <- evidence_mat$gene
    mode(mat) <- "integer"

    comb_mat <- make_comb_mat(mat)

    pdf(file.path(out_dir, "panels", "panel_h.pdf"), width = fig_half_width, height = 3)
    ht <- UpSet(comb_mat,
                set_order = colnames(mat),
                comb_order = order(-comb_size(comb_mat)),
                top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE,
                                                       numbers_gp = gpar(fontsize = 5)),
                right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE,
                                                           numbers_gp = gpar(fontsize = 6)),
                row_names_gp = gpar(fontsize = 6))
    draw(ht)
    dev.off()
    cat("  Saved panel_h.pdf\n")
  } else {
    cat("  WARNING: Too few ncRNA DEGs with evidence for UpSet plot\n")
  }
} else {
  cat("  WARNING: No ncRNA DEGs or ComplexHeatmap unavailable\n")
  ph <- placeholder("Panel (h)\nMulti-evidence UpSet")
  save_fig(ph, file.path(out_dir, "panels", "panel_h.pdf"), width = fig_half_width, height = 3)
}

# =============================================================================
# Assemble Composite Figure
# =============================================================================
cat("\n--- Assembling composite figure ---\n")

# The panels are saved individually. Attempt assembly with patchwork
# for panels that are ggplot objects.
panel_files <- list.files(file.path(out_dir, "panels"), pattern = "panel_.*\\.pdf$", full.names = TRUE)
cat(sprintf("  Panel PDFs generated: %d / 8\n", length(panel_files)))

# Note: full assembly may need manual adjustment in Illustrator.
# Provide a basic ggplot-only assembly for panels a, b, e, f, g
tryCatch({
  assembled <- (pa | pb) / (pe) / (pf | pg) +
    plot_annotation(tag_levels = list(c("a", "b", "e", "f", "g")))
  save_fig(assembled, file.path(out_dir, "figS_ncrna_partial.pdf"),
           width = fig_full_width, height = 10)
  cat("  Saved figS_ncrna_partial.pdf (ggplot panels only)\n")
}, error = function(e) {
  cat(sprintf("  Assembly skipped: %s\n", e$message))
})

cat("\n=== figS_ncrna Complete ===\n")
cat("End:", format(Sys.time()), "\n")
cat("Outputs in:", out_dir, "\n")
