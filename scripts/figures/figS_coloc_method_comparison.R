#!/usr/bin/env Rscript
# figS_coloc_method_comparison.R — Supplementary: ABF COLOC vs SuSiE-COLOC vs INTACT
#
# 4 panels:
# (a) ABF vs SuSiE-COLOC PP.H4 scatter (method concordance)
# (b) ABF vs INTACT score scatter (COLOC vs Bayesian-filtered)
# (c) UpSet plot: 3-way overlap at thresholds
# (d) Venn-style summary with gene counts
#
# Output: figures/supplementary/figS04_coloc/figS_coloc_method_comparison.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

outdir <- FIGS04_DIR
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGS04_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# Load data
# ===========================================================================
cat("Loading data...\n")

atlas <- fread(file.path(BASE,
  "RNA-seq/results/gwas_rna_integration/gwas_rna_scrna_integrated.csv"),
  select = c("gene", "coloc_susie_best_pp4", "intact_score_bulk"))

coloc_gl <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc_gl <- coloc_gl[gene != "" & !is.na(gene)]

# Merge ABF COLOC (best PP.H4 across 24 GWAS)
dt <- merge(atlas, coloc_gl[, .(gene, abf_pp4 = coloc_best_pp4)], by = "gene", all.x = TRUE)
dt <- dt[!is.na(abf_pp4) | !is.na(coloc_susie_best_pp4) | !is.na(intact_score_bulk)]

cat("  Genes with any method:", nrow(dt), "\n")

# Priority genes for labeling
priority <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
              "THRB", "NR1H4", "MARC1")

# ===========================================================================
# Panel A: ABF COLOC vs SuSiE-COLOC
# ===========================================================================
cat("Panel A: ABF vs SuSiE-COLOC...\n")

dt_as <- dt[!is.na(abf_pp4) & !is.na(coloc_susie_best_pp4)]
rho_as <- cor(dt_as$abf_pp4, dt_as$coloc_susie_best_pp4, method = "spearman")

dt_as[, concordance := fcase(
  abf_pp4 > 0.9 & coloc_susie_best_pp4 > 0.9, "Both high",
  abf_pp4 > 0.9 & coloc_susie_best_pp4 <= 0.9, "ABF only",
  abf_pp4 <= 0.9 & coloc_susie_best_pp4 > 0.9, "SuSiE only",
  default = "Both low"
)]
dt_as[, show_label := gene %in% priority |
        (concordance %in% c("ABF only", "SuSiE only") & (abf_pp4 > 0.95 | coloc_susie_best_pp4 > 0.95))]

p_a <- ggplot(dt_as, aes(x = abf_pp4, y = coloc_susie_best_pp4, color = concordance)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 0.9, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0.9, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = c("Both high" = "#D32F2F", "ABF only" = "#FF9800",
                                "SuSiE only" = "#7B1FA2", "Both low" = "#E0E0E0"),
                     name = NULL) +
  geom_text_repel(data = dt_as[show_label == TRUE], aes(label = gene),
                  color = "black", size = 1.8, fontface = "italic",
                  max.overlaps = 15, segment.size = 0.15, seed = 42,
                  box.padding = 0.4, force = 2) +
  annotate("text", x = 0.05, y = 0.95,
           label = sprintf("rho = %.3f\nn = %s", rho_as, format(nrow(dt_as), big.mark = ",")),
           size = 2, hjust = 0, vjust = 1, color = "grey30") +
  labs(x = "ABF COLOC PP.H4 (best across 24 GWAS)",
       y = "SuSiE-COLOC PP.H4 (fine-mapping)",
       title = "ABF vs SuSiE-COLOC concordance") +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 5))

# ===========================================================================
# Panel B: ABF COLOC vs INTACT
# ===========================================================================
cat("Panel B: ABF vs INTACT...\n")

dt_ai <- dt[!is.na(abf_pp4) & !is.na(intact_score_bulk)]
rho_ai <- cor(dt_ai$abf_pp4, dt_ai$intact_score_bulk, method = "spearman")

dt_ai[, concordance := fcase(
  abf_pp4 > 0.9 & intact_score_bulk > 0.5, "Both high",
  abf_pp4 > 0.9 & intact_score_bulk <= 0.5, "COLOC only",
  abf_pp4 <= 0.9 & intact_score_bulk > 0.5, "INTACT only",
  default = "Both low"
)]
dt_ai[, show_label := gene %in% priority |
        (concordance %in% c("COLOC only", "INTACT only") & (abf_pp4 > 0.95 | intact_score_bulk > 0.7))]

p_b <- ggplot(dt_ai, aes(x = abf_pp4, y = intact_score_bulk, color = concordance)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0.9, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = c("Both high" = "#D32F2F", "COLOC only" = "#FF9800",
                                "INTACT only" = "#4CAF50", "Both low" = "#E0E0E0"),
                     name = NULL) +
  geom_text_repel(data = dt_ai[show_label == TRUE], aes(label = gene),
                  color = "black", size = 1.8, fontface = "italic",
                  max.overlaps = 15, segment.size = 0.15, seed = 42,
                  box.padding = 0.4, force = 2) +
  annotate("text", x = 0.05, y = 0.95,
           label = sprintf("rho = %.3f\nn = %s", rho_ai, format(nrow(dt_ai), big.mark = ",")),
           size = 2, hjust = 0, vjust = 1, color = "grey30") +
  labs(x = "ABF COLOC PP.H4 (best across 24 GWAS)",
       y = "Multi-INTACT score (bulk + scTWAS × COLOC)",
       title = "ABF COLOC vs Multi-INTACT concordance") +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 5))

# ===========================================================================
# Panel C: 3-way set sizes as stacked bar / UpSet-style
# ===========================================================================
cat("Panel C: 3-way overlap...\n")

abf_set <- dt[!is.na(abf_pp4) & abf_pp4 > 0.9]$gene
susie_set <- dt[!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.9]$gene
intact_set <- dt[!is.na(intact_score_bulk) & intact_score_bulk > 0.5]$gene

# Compute all 7 intersection categories
dt[, in_abf := gene %in% abf_set]
dt[, in_susie := gene %in% susie_set]
dt[, in_intact := gene %in% intact_set]
dt[, any_method := in_abf | in_susie | in_intact]

intersections <- dt[any_method == TRUE, .(
  category = fcase(
    in_abf & in_susie & in_intact, "All three",
    in_abf & in_susie & !in_intact, "ABF + SuSiE",
    in_abf & !in_susie & in_intact, "ABF + INTACT",
    !in_abf & in_susie & in_intact, "SuSiE + INTACT",
    in_abf & !in_susie & !in_intact, "ABF only",
    !in_abf & in_susie & !in_intact, "SuSiE only",
    !in_abf & !in_susie & in_intact, "INTACT only"
  )
), by = gene]

int_counts <- intersections[, .N, by = category]
int_counts[, category := factor(category,
  levels = c("All three", "ABF + SuSiE", "ABF + INTACT", "SuSiE + INTACT",
             "ABF only", "SuSiE only", "INTACT only"))]
setorder(int_counts, -N)

# Color by overlap level
int_counts[, overlap_level := fcase(
  grepl("All three", category), "3 methods",
  grepl("\\+", category), "2 methods",
  default = "1 method"
)]

p_c <- ggplot(int_counts, aes(x = reorder(category, N), y = N, fill = overlap_level)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.2, size = 2.2) +
  scale_fill_manual(values = c("3 methods" = "#D32F2F", "2 methods" = "#FF9800",
                               "1 method" = "#42A5F5"),
                    name = "Overlap") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  coord_flip() +
  labs(x = NULL, y = "Number of genes",
       title = sprintf("3-way overlap (ABF>0.9, SuSiE>0.9, INTACT>0.5)\n%d unique genes",
                       sum(dt$any_method))) +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.95, 0.3), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.key.size = unit(0.2, "cm"))

# ===========================================================================
# Panel D: Method summary table as plot
# ===========================================================================
cat("Panel D: Summary...\n")

summary_dt <- data.table(
  Method = c("ABF COLOC\n(24 GWAS × Broadaway eQTL)",
             "SuSiE-COLOC\n(fine-mapping, multi-signal)",
             "INTACT\n(scTWAS × COLOC)"),
  Threshold = c("PP.H4 > 0.9", "PP.H4 > 0.9", "Score > 0.5"),
  N_genes = c(length(abf_set), length(susie_set), length(intact_set)),
  Description = c("Single-variant colocalization\nacross 24 GWAS",
                   "Multi-signal fine-mapping\nthen colocalization",
                   "TWAS filtered by\ncolocalization confirmation")
)

summary_dt[, Method := factor(Method, levels = rev(Method))]

p_d <- ggplot(summary_dt, aes(x = Method, y = N_genes, fill = Method)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(N_genes, " genes\n(", Threshold, ")")),
            hjust = -0.05, size = 2, lineheight = 0.9) +
  scale_fill_manual(values = c("#1565C0", "#7B1FA2", "#D32F2F"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  coord_flip() +
  labs(x = NULL, y = "Genes above threshold",
       title = "Method comparison summary") +
  theme_masld(base_size = 7)

# ===========================================================================
# Combine and save
# ===========================================================================
cat("Assembling...\n")

p_combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d"))) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(file.path(outdir, "figS_coloc_method_comparison.pdf"),
       p_combined, width = 170, height = 150, units = "mm", device = cairo_pdf)
cat("Saved: figS_coloc_method_comparison.pdf\n")
cat("Done.\n")
