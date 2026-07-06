#!/usr/bin/env Rscript
# fig2_coloc_summary.R  (2026-06-19)  — Fig 2 comprehensive COLOC summary
# Two DRAFT formats summarizing ALL significant GWAS-eQTL colocalizations
# (genes, genomic location, ancestry, coding/non-coding, DEG direction):
#   Draft A: gene x ancestry PP.H4 heatmap (ComplexHeatmap) with chr/class/DEG
#            row annotations; genes ordered by genomic position.
#   Draft B: genomic "coloc Manhattan" — genes plotted at chromosomal position,
#            y = best PP.H4.susie, colored by coding/non-coding, sized by
#            #ancestries, top genes labeled.
# Out: figures/main/fig2_genetics/panels/coloc_summary_heatmap.pdf
#      figures/main/fig2_genetics/panels/coloc_summary_manhattan.pdf  (+ source CSV)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
  library(ComplexHeatmap); library(circlize); library(grid)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

# ── assemble per-gene table ──────────────────────────────────────────────────
sc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
# MAIN (Tier-1/2, liver-specific) restriction (2026-07-06): keep only placement=="main"
# strata (NAFLD/NASH/PDFF + ALT/AST/GGT); the Tier-3/4 supp strata (MVP Cirrhosis/
# ChronLiver/Albumin/Platelet) move to a supplementary full-portfolio figure.
MAIN_STUDIES <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[
  placement == "main", study_name]
sc <- sc[gwas_name %in% MAIN_STUDIES]
# ancestry from the GWAS registry (Tier-1/2 portfolio incl. MVP NAFLD/ALT/AST AMR/AFR/EAS/
# EUR strata) — NOT the retired grepl() heuristic, which misrouted every MVP stratum into EUR.
sc[, ancestry := as.character(gwas_ancestry(gwas_name))]
sc[, pos := as.integer(sub("^[0-9]+:", "", top_snp))]
sc <- sc[!is.na(PP.H4.susie)]
# best PP.H4 per gene x ancestry
ba <- sc[, .(pp4 = max(PP.H4.susie)), by = .(gene, ancestry)]
mat_dt <- dcast(ba, gene ~ ancestry, value.var = "pp4")
# overall best + its chr/pos (genomic location)
best <- sc[, .SD[which.max(PP.H4.susie)], by = gene][, .(gene, chr, pos, best = PP.H4.susie)]
g <- merge(mat_dt, best, by = "gene")
for (a in GWAS_ANCESTRY_LEVELS) if (!a %in% names(g)) g[[a]] <- NA_real_
g[, n_anc := rowSums(.SD > 0.9, na.rm = TRUE), .SDcols = GWAS_ANCESTRY_LEVELS]
g[, n_anc5 := rowSums(.SD > 0.5, na.rm = TRUE), .SDcols = GWAS_ANCESTRY_LEVELS]

# class (coding / non-coding)
cva <- fread(file.path(BASE, "RNA-seq/results/coloc_variant_classes/coloc_variant_annotation.csv"))
cls <- cva[!is.na(gene_symbol) & gene_symbol != ""][order(-pp4_best),
        .(coarse_class = coarse_class[1]), by = gene_symbol]
g <- merge(g, cls, by.x = "gene", by.y = "gene_symbol", all.x = TRUE)
g[is.na(coarse_class), coarse_class := "non-coding"]
# DEG direction
deg <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
# dedupe: a symbol can map to >1 Ensembl row — without this the merge duplicates
# coloc genes (9 of them) and inflates the counts (368 -> 378, 210 -> 213).
degu <- deg[order(padj)][!duplicated(symbol), .(symbol, logFC, padj)]
g <- merge(g, degu, by.x = "gene", by.y = "symbol", all.x = TRUE)
g[, deg_dir := fifelse(!is.na(padj) & padj < 0.05 & logFC > 0, "Up",
              fifelse(!is.na(padj) & padj < 0.05 & logFC < 0, "Down", "n.s."))]
g[, chr := as.integer(chr)]
setorder(g, chr, pos)

fwrite(g[best > 0.5, .(gene, chr, pos, EUR = round(EUR,3), AFR = round(AFR,3),
        AMR = round(AMR,3), EAS = round(EAS,3), SAS = round(SAS,3),
        best = round(best,3), n_anc_gt0.9 = n_anc,
        class = coarse_class, deg_dir)],
       file.path(PANEL_DIR, "coloc_summary_source.csv"))

# ════════════════════════════════════════════════════════════════════════════
# Draft A — gene x ancestry heatmap (genes >0.9 in >=1 ancestry)
# ════════════════════════════════════════════════════════════════════════════
hi <- g[best > 0.9]
setorder(hi, chr, pos)
M <- as.matrix(hi[, .(EUR, AFR, AMR, EAS, SAS)]); rownames(M) <- hi$gene
col_pp4 <- colorRamp2(c(0.5, 0.75, 1.0), c("#E8F0EF", "#5BA89C", "#00695C"))
cls_col <- c("coding" = "#C9265E", "non-coding" = "#1565C0")
deg_col <- c("Up" = "#C0392B", "Down" = "#2471A3", "n.s." = "#D5D8DC")
# label only multi-ancestry (>=2 of the 5 ancestries) genes to keep it readable
lab_idx <- which(hi$n_anc >= 2)
left_anno <- rowAnnotation(
  Class = hi$coarse_class, DEG = hi$deg_dir,
  col = list(Class = cls_col, DEG = deg_col),
  annotation_name_gp = gpar(fontsize = 7), simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(Class = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)),
                                 DEG = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6))))
right_mark <- rowAnnotation(mark = anno_mark(at = lab_idx, labels = hi$gene[lab_idx],
                            labels_gp = gpar(fontsize = 5.5, fontface = "italic"), link_width = unit(4, "mm")))
ht <- Heatmap(M, name = "PP.H4", col = col_pp4, na_col = "grey94",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 8, fontface = "bold"), column_names_rot = 0,
  column_names_centered = TRUE,
  row_title = sprintf("%d high-confidence genes (PP.H4 > 0.9; subset of %d at > 0.5), by chromosome",
                      nrow(M), nrow(g[best > 0.5])),
  row_title_gp = gpar(fontsize = 7),
  left_annotation = left_anno, right_annotation = right_mark,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6),
                              at = c(0.5, 0.75, 1.0)),
  width = unit(2.6, "cm"))

pdf(file.path(PANEL_DIR, "coloc_summary_heatmap.pdf"), width = 5.2, height = 8.6)
draw(ht, column_title = "Cross-ancestry colocalization landscape",
     column_title_gp = gpar(fontsize = 9, fontface = "bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
dev.off()
cat(sprintf("[summary A] heatmap: %d genes >0.9 (%d multi-ancestry labeled)\n", nrow(M), length(lab_idx)))

# ════════════════════════════════════════════════════════════════════════════
# Draft B — genomic "coloc Manhattan" (all genes >0.5)
# ════════════════════════════════════════════════════════════════════════════
hg19 <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
          146364022,141213431,135534747,135006516,133851895,115169878,107349540,
          102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
off <- c(0, cumsum(as.numeric(hg19)))[1:22]; names(off) <- 1:22
gm <- g[best > 0.5 & chr %in% 1:22]
gm[, gx := off[as.character(chr)] + pos]
gm[, class := factor(coarse_class, levels = c("non-coding","coding"))]
axis_df <- data.table(chr = 1:22, center = off + hg19/2)
axis_df <- axis_df[chr %in% unique(gm$chr)]
# label: tri/multi-ancestry + top PP.H4 + the coding ones (rare, notable)
gm[, lab := fifelse(n_anc >= 2 | coarse_class == "coding", gene, "")]
band <- data.table(chr = 1:22, xmin = off, xmax = off + hg19)[chr %in% 1:22][chr %% 2 == 0]

pB <- ggplot(gm, aes(gx, best)) +
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = 0.48, ymax = 1.02),
            inherit.aes = FALSE, fill = "grey95") +
  geom_hline(yintercept = c(0.8, 0.9), linetype = "22", linewidth = 0.25, color = "grey70") +
  geom_point(aes(color = class, size = n_anc5), alpha = 0.85) +
  geom_text_repel(aes(label = lab, color = class), size = 2, fontface = "italic",
                  max.overlaps = 20, segment.size = 0.15, min.segment.length = 0.2,
                  box.padding = 0.2, show.legend = FALSE) +
  scale_color_manual(values = c("non-coding" = "#1565C0", "coding" = "#C9265E"),
                     name = "Lead variant") +
  scale_size_continuous(range = c(0.8, 2.6), breaks = c(1,2,3,4,5), name = "# ancestries") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr, expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0.48, 1.02), breaks = c(0.5,0.7,0.9), expand = c(0, 0)) +
  labs(x = "Chromosome", y = expression("Best PP.H"[4]*" (SuSiE-coloc)"),
       title = sprintf("Colocalization landscape: %d SuSiE-coloc effector genes (PP.H4 > 0.5)", nrow(gm))) +
  theme_masld(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.direction = "horizontal", legend.box = "horizontal",
        legend.key.size = unit(0.2, "cm"), legend.title = element_text(size = 6),
        legend.text = element_text(size = 5.5),
        axis.text.x = element_text(size = 6))

save_fig(pB, file.path(PANEL_DIR, "coloc_summary_manhattan.pdf"),
         width = fig_full_width * 0.92, height = 3.0)
cat(sprintf("[summary B] manhattan: %d genes >0.5 plotted\n", nrow(gm)))
