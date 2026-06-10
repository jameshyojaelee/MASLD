#!/usr/bin/env Rscript
# 05_render_fig3_main.R
# Three compact main-fig3 panels — F1, F2, F6 from the sketches:
#   fig3g — COLOC genes by GWAS-variant class (F1)
#   fig3h — Tier-1 DEG → COLOC funnel (F2)
#   fig3i — DEG × COLOC scatter, fig3e style (F6)
#
# Consistent typography across panels:
#   title       = 11 bold
#   subtitle    = 9  grey30
#   axis title  = 10
#   axis text   = 9
#   legend      = 8.5 (title) / 8 (text)
#   data label  = 3.4 bold
#   gene label  = 2.4 italic (text-repel)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
  library(scales)
})

BASE   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT    <- file.path(BASE, "figures/main/fig3_regulatory_architecture/panels")
DATDIR <- file.path(BASE, "RNA-seq/results/coloc_variant_classes")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

fine_levels <- c("coding","spliceSite","fiveUTR","threeUTR",
                 "promoter","intron","intergenic")
# Pretty labels (apostrophes for UTRs)
fine_label <- c(
  "coding"     = "coding",
  "spliceSite" = "splice site",
  "fiveUTR"    = "5' UTR",
  "threeUTR"   = "3' UTR",
  "promoter"   = "promoter",
  "intron"     = "intron",
  "intergenic" = "intergenic")
fine_pal <- c(
  "coding"     = "#C2185B",
  "spliceSite" = "#880E4F",
  "fiveUTR"    = "#9CCC65",
  "threeUTR"   = "#7CB342",
  "promoter"   = "#1565C0",
  "intron"     = "#42A5F5",
  "intergenic" = "#90A4AE")
coarse_pal <- c("coding" = "#C2185B", "non-coding" = "#42A5F5",
                "no COLOC" = "gray80")

base_theme <- theme_masld() +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 9),
    legend.title  = element_text(size = 8.5),
    legend.text   = element_text(size = 8),
    legend.key.size = unit(0.32, "cm"),
    strip.text    = element_text(size = 10, face = "bold"))

save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, name), p, width = w, height = h, device = cairo_pdf)
  cat(" wrote", name, "\n")
}

# ===========================================================================
# fig3g — All COLOC genes split by GWAS-variant class (F1)
# ===========================================================================
f1 <- fread(file.path(DATDIR, "panel_F1_data.csv"))
f1[, fine_class := factor(fine_class, levels = fine_levels)]
f1_totals <- f1[, .(total = sum(n_genes)), by = method]

fig3g <- ggplot(f1, aes(x = method, y = n_genes, fill = fine_class)) +
  geom_col(width = 0.62) +
  geom_text(data = f1_totals, aes(x = method, y = total, label = total),
            inherit.aes = FALSE, vjust = -0.35, size = 3.4,
            fontface = "bold") +
  scale_fill_manual(values = fine_pal, labels = fine_label, drop = FALSE,
                    name = "GWAS variant class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::label_comma()) +
  labs(x = NULL, y = "# genes",
       title = "COLOC genes by where the COLOC top-H4 SNP falls (PP.H4 ≥ 0.5)") +
  base_theme +
  theme(legend.position = "right")

save_panel(fig3g, "coloc_variant_class_split.pdf", w = 5.6, h = 4.2)

# ===========================================================================
# fig3h — Tier-1 DEG → COLOC cascade (F2)
# ===========================================================================
funnel <- fread(file.path(DATDIR, "panel_F2_funnel.csv"))
hp <- fread(file.path(DATDIR, "panel_F2_hypergeom.csv"))
funnel[, step := factor(step, levels = step)]
funnel[, group := fcase(
  grepl("Tier", step),       "DEG total",
  grepl("DEGs",  step),      "DEG ∩ COLOC",
  grepl("Coding", step),     "Coding",
  grepl("Non-coding", step), "Non-coding")]
funnel[, group := factor(group,
   levels = c("DEG total","DEG ∩ COLOC","Coding","Non-coding"))]
fill_pal <- c("DEG total"     = "#37474F",
              "DEG ∩ COLOC"   = "#1565C0",
              "Coding"        = "#C2185B",
              "Non-coding"    = "#42A5F5")

fig3h <- ggplot(funnel, aes(x = count, y = step, fill = group)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = scales::label_comma()(count)),
            hjust = -0.18, size = 3.4, fontface = "bold") +
  scale_x_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.18))) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = fill_pal, guide = "none") +
  labs(x = "# genes", y = NULL,
       title = sprintf(
         "Tier-1 DEGs → COLOC: %d/%d DEG-COLOC have GWAS variant in CDS (vs %d/%d in all COLOC)",
         hp$n_deg_coloc_coding[1], hp$n_deg_coloc[1],
         hp$n_all_coloc_coding[1], hp$n_all_coloc[1])) +
  base_theme +
  theme(axis.text.y = element_text(size = 8.5),
        plot.title = element_text(face = "bold", size = 10.5))

save_panel(fig3h, "deg_coloc_funnel.pdf", w = 7.0, h = 3.5)

# Drop the combined biotype panel from the previous iteration
old_combined <- file.path(OUT, "fig3g_deg_coloc_biotypes.pdf")
if (file.exists(old_combined)) {
  file.remove(old_combined)
  cat(" removed fig3g_deg_coloc_biotypes.pdf (replaced by separate fig3g + fig3h)\n")
}

# ===========================================================================
# fig3i — All 797 COLOC genes plotted by their C2 (limma-voom-qw) log2FC × PP.H4
#          Label every coding-variant gene + the top non-coding hits.
#          Dashed vlines at logFC = ±0.5 (Tier-1 DEG threshold).
# ===========================================================================
deg_all <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("logFC","padj","symbol"))
deg_all <- deg_all[symbol != "" & !is.na(symbol)]
setorder(deg_all, symbol, padj)
deg_all <- deg_all[!duplicated(symbol)]

defC_g <- fread(file.path(DATDIR, "coloc_variant_annotation.csv"),
                select = c("gene_symbol","pp4_best","coarse_class","fine_class"))
setorder(defC_g, gene_symbol, -pp4_best)
defC_g <- defC_g[!duplicated(gene_symbol)]

# Inner join: COLOC genes ⨝ dream results; require padj < 0.05 (classic DEG
# significance), but no |logFC| filter — that lets readers see the LFC range
# of significant COLOC genes including the ones below the Tier-1 0.5 cutoff.
coloc_pts <- merge(defC_g, deg_all, by.x = "gene_symbol", by.y = "symbol",
                   all.x = FALSE)
coloc_pts <- coloc_pts[!is.na(padj) & padj < 0.05]
coloc_pts[, coarse_class := factor(coarse_class,
   levels = c("coding","non-coding"))]
coloc_pts[, is_tier1_DEG := abs(logFC) > 0.5]

n_coding    <- coloc_pts[coarse_class == "coding", .N]
n_noncoding <- coloc_pts[coarse_class == "non-coding", .N]
cat(sprintf("fig3i scatter: %d COLOC genes at padj<0.05 (%d coding, %d non-coding)\n",
            nrow(coloc_pts), n_coding, n_noncoding))

# Always label the 35 coding-variant genes (the user wants to see them all)
coding_labels <- coloc_pts[coarse_class == "coding"]
# Plus top 12 non-coding by |logFC| × PP.H4 to anchor the visual story
noncoding_top <- coloc_pts[coarse_class == "non-coding"]
noncoding_top[, score := abs(logFC) * pp4_best]
noncoding_top <- noncoding_top[order(-score)][seq_len(min(12, .N))]
labelled <- rbind(coding_labels, noncoding_top, fill = TRUE)

fig3i <- ggplot(coloc_pts, aes(x = logFC, y = pp4_best, color = coarse_class)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray55") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             linewidth = 0.25, color = "gray55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
  # Order so coding points draw on top of non-coding
  geom_point(data = coloc_pts[coarse_class == "non-coding"],
             size = 1.6, alpha = 0.55) +
  geom_point(data = coloc_pts[coarse_class == "coding"],
             size = 2.2, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = labelled, aes(label = gene_symbol),
    size = 2.3, max.overlaps = Inf, force = 3.5, force_pull = 0.8,
    box.padding = 0.32, point.padding = 0.15,
    segment.size = 0.2, segment.color = "gray55",
    min.segment.length = 0, fontface = "italic",
    bg.color = "white", bg.r = 0.12,
    show.legend = FALSE, max.iter = 10000, seed = 1) +
  scale_color_manual(values = coarse_pal, drop = FALSE,
                     name = "GWAS variant class",
                     breaks = c("coding","non-coding")) +
  scale_y_continuous(limits = c(0.45, 1.05),
                     breaks = c(0.5, 0.7, 0.9, 1.0),
                     labels = c("0.5","0.7","0.9","1")) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  labs(x = expression("Bulk dream log"[2]*"FC (MASLD vs control)"),
       y = "Best COLOC PP.H4",
       title = "DEG vs. COLOC",
       subtitle = sprintf(
         paste0("%d COLOC genes (padj < 0.05; PP.H4 ≥ 0.5); dashed lines: |logFC| = 0.5 Tier-1 cutoff\n",
                "Color = where the COLOC top-H4 SNP falls: ",
                "coding (CDS) vs non-coding (intron / UTR / promoter / intergenic)"),
         nrow(coloc_pts))) +
  base_theme +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 8.5, color = "gray30",
                                      lineheight = 1.1)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), nrow = 1))

save_panel(fig3i, "deg_coloc_scatter.pdf", w = 6.4, h = 5.4)

# Drop the temporary fig3h_deg_coloc_scatter.pdf from prior iteration
old_h_scatter <- file.path(OUT, "fig3h_deg_coloc_scatter.pdf")
if (file.exists(old_h_scatter)) {
  file.remove(old_h_scatter)
  cat(" removed fig3h_deg_coloc_scatter.pdf (scatter is now fig3i again)\n")
}

cat("\n[05_render_fig3_main] DONE — wrote 3 panels to", OUT, "\n")
