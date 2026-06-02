#!/usr/bin/env Rscript
# 04_render_figures.R
# F1-F8 panels into figures/sketches/coloc_coding_noncoding/. PDF only.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
  library(scales);     library(ggalluvial)
})

BASE   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT    <- file.path(BASE, "figures/sketches/coloc_coding_noncoding")
DATDIR <- file.path(BASE, "RNA-seq/results/coloc_variant_classes")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# ---- Color palette ---------------------------------------------------------
coarse_pal <- c("coding" = "#C2185B", "non-coding" = "#42A5F5")
fine_pal   <- c("coding"     = "#C2185B",
                "spliceSite" = "#880E4F",
                "fiveUTR"    = "#7E57C2",
                "threeUTR"   = "#5E35B1",
                "promoter"   = "#1565C0",
                "intron"     = "#42A5F5",
                "intergenic" = "#90A4AE")
fine_levels <- c("coding","spliceSite","fiveUTR","threeUTR",
                 "promoter","intron","intergenic")
def_pal <- c("A_lead_gwas"   = "#37474F",
             "B_finemapped"  = "#1565C0",
             "C_coloc_top"   = "#C2185B")
ancestry_pal <- c(EUR = "#1565C0", EAS = "#C2185B",
                  AFR = "#388E3C", SAS = "#6A1B9A")

save_panel <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT, name), p, width = w, height = h, device = cairo_pdf)
  cat(" wrote", name, "\n")
}

# ---- F1: All COLOC genes — counts by method × class ------------------------
f1 <- fread(file.path(DATDIR, "panel_F1_data.csv"))
f1[, fine_class := factor(fine_class, levels = fine_levels)]
p1a <- ggplot(f1, aes(x = method, y = n_genes, fill = fine_class)) +
  geom_col(width = 0.65) +
  geom_text(data = f1[, .(n = sum(n_genes)), by = method],
            aes(x = method, y = n, label = n), inherit.aes = FALSE,
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = fine_pal, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "# genes", fill = "Fine class",
       subtitle = "All COLOC genes (PP.H4 ≥ 0.5) by method") +
  theme_masld() + theme(legend.position = "right",
                        plot.subtitle = element_text(size = 11))
p1b <- ggplot(f1, aes(x = method, y = n_genes, fill = fine_class)) +
  geom_col(position = "fill", width = 0.65) +
  scale_fill_manual(values = fine_pal, drop = FALSE) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion", fill = "Fine class",
       subtitle = "Proportions") +
  theme_masld() + theme(legend.position = "none",
                        plot.subtitle = element_text(size = 11))
save_panel((p1a / p1b) +
             plot_annotation(title = "F1 — COLOC genes split by variant class",
                             theme = theme(plot.title = element_text(face = "bold"))),
           "panel_F1_coloc_genes_class_split.pdf", w = 8, h = 7)

# ---- F2: DEG funnel --------------------------------------------------------
funnel <- fread(file.path(DATDIR, "panel_F2_funnel.csv"))
hp <- fread(file.path(DATDIR, "panel_F2_hypergeom.csv"))
funnel[, step := factor(step, levels = step)]
funnel[, group := fcase(
  grepl("Tier", step), "DEG total",
  grepl("DEGs",  step), "DEG ∩ COLOC",
  grepl("Coding", step),     "Coding subset",
  grepl("Non-coding", step), "Non-coding subset")]
funnel[, group := factor(group, levels = unique(group))]
fill_pal <- c("DEG total"           = "#37474F",
              "DEG ∩ COLOC"         = "#1565C0",
              "Coding subset"       = "#C2185B",
              "Non-coding subset"   = "#42A5F5")
p2 <- ggplot(funnel, aes(x = count, y = step, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::label_comma()(count)),
            hjust = -0.18, size = 4, fontface = "bold") +
  scale_x_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.18))) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = fill_pal) +
  labs(x = "# genes", y = NULL,
       title = "F2 — DEG cascade: Tier-1 → COLOC → coding/non-coding",
       subtitle = sprintf(
         "%d DEGs ∩ COLOC at PP.H4 ≥ 0.5 — hypergeometric p = %.3g",
         hp$n_deg_coloc[1], hp$hypergeom_p[1])) +
  theme_masld() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray30"))
save_panel(p2, "panel_F2_deg_coloc_funnel.pdf", w = 9, h = 4.5)

# ---- F3: Lead causal vs COLOC top concordance ------------------------------
ac <- fread(file.path(DATDIR, "panel_F3_AC.csv"))
bc <- fread(file.path(DATDIR, "panel_F3_BC.csv"))
ac[, A_match_class := factor(A_match_class,
   levels = c("coding","non-coding","no_match"))]
bc[, B_match_class := factor(B_match_class,
   levels = c("coding","non-coding","no_match"))]
p3a <- ggplot(ac, aes(x = C_coarse, y = n, fill = A_match_class)) +
  geom_col(width = 0.65, position = "stack") +
  scale_fill_manual(values = c("coding" = "#C2185B",
                                "non-coding" = "#42A5F5",
                                "no_match" = "#CFD8DC")) +
  labs(x = "Definition C (COLOC top SNP) class", y = "# (gene, study) pairs",
       fill = "Definition A (GWAS lead)\nmatch class",
       subtitle = "A vs C — matched only when same variant_key") +
  theme_masld() +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 10))
p3b <- ggplot(bc, aes(x = C_coarse, y = n, fill = B_match_class)) +
  geom_col(width = 0.65, position = "stack") +
  scale_fill_manual(values = c("coding" = "#C2185B",
                                "non-coding" = "#42A5F5",
                                "no_match" = "#CFD8DC")) +
  labs(x = "Definition C (COLOC top SNP) class", y = "# (gene, study) pairs",
       fill = "Definition B (top finemapped)\nmatch class",
       subtitle = "B vs C — matched only when same variant_key") +
  theme_masld() +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 10))
save_panel((p3a | p3b) +
             plot_annotation(title = "F3 — Lead causal vs COLOC top SNP",
                             theme = theme(plot.title = element_text(face = "bold"))),
           "panel_F3_lead_vs_coloc_top.pdf", w = 11, h = 5)

# ---- F4: Per-ancestry ------------------------------------------------------
f4 <- fread(file.path(DATDIR, "panel_F4_ancestry.csv"))
f4[, ancestry := factor(ancestry, levels = c("EUR","EAS","AFR","SAS"))]
f4[, definition := factor(definition,
   levels = c("A_lead_gwas","B_finemapped","C_coloc_top"),
   labels = c("A: GWAS lead","B: Top finemapped","C: COLOC top"))]
p4 <- ggplot(f4, aes(x = ancestry, y = n, fill = coarse_class)) +
  geom_col(position = "fill", width = 0.65) +
  facet_grid(~ definition) +
  scale_fill_manual(values = coarse_pal) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion", fill = "Class",
       title = "F4 — Coding vs non-coding by ancestry × definition") +
  theme_masld() +
  theme(plot.title = element_text(face = "bold"))
save_panel(p4, "panel_F4_ancestry_class.pdf", w = 11, h = 4.5)

# ---- F5: Per-trait ---------------------------------------------------------
f5 <- fread(file.path(DATDIR, "panel_F5_trait.csv"))
trait_levels <- c("NAFLD/NASH","HCC","Cirrhosis","PDFF","Liver enzymes","Other")
f5[, trait_cat := factor(trait_cat, levels = trait_levels)]
f5[, definition := factor(definition,
   levels = c("A_lead_gwas","B_finemapped","C_coloc_top"),
   labels = c("A: GWAS lead","B: Top finemapped","C: COLOC top"))]
p5 <- ggplot(f5, aes(x = trait_cat, y = n, fill = coarse_class)) +
  geom_col(position = "fill", width = 0.65) +
  facet_grid(~ definition) +
  scale_fill_manual(values = coarse_pal) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion", fill = "Class",
       title = "F5 — Coding vs non-coding by trait × definition") +
  theme_masld() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 9))
save_panel(p5, "panel_F5_trait_class.pdf", w = 12, h = 4.5)

# ---- F6: DEG × COLOC PP4 — fig3e style scatter ----------------------------
# Match the main fig3e aesthetic: gray atlas-DEG background, colored
# coding/non-coding foreground, text-repel italic labels with white halo,
# PP4 ref lines at 0.5 and 0.9, vline at logFC=0.
f6_fg <- fread(file.path(DATDIR, "panel_F6_lfc_pp4.csv"))
# Assemble background: all Tier-1 DEGs with their best PP4 across canonical
# COLOC tests (gene-level rollup). DEGs that don't appear in defC get PP4=0
# (i.e., land on the bottom of the plot).
deg <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv"),
  select = c("logFC","padj","symbol"))
deg <- deg[!is.na(padj) & padj < 0.05 & abs(logFC) > 0.5 &
            symbol != "" & !is.na(symbol)]
setorder(deg, symbol, padj)
deg <- deg[!duplicated(symbol)]
defC_g <- fread(file.path(DATDIR, "coloc_variant_annotation.csv"),
                select = c("gene_symbol","pp4_best","coarse_class","fine_class"))
setorder(defC_g, gene_symbol, -pp4_best)
defC_g <- defC_g[!duplicated(gene_symbol)]
deg <- merge(deg, defC_g, by.x = "symbol", by.y = "gene_symbol", all.x = TRUE)
deg[is.na(pp4_best),     pp4_best     := 0]
deg[is.na(coarse_class), coarse_class := "no COLOC"]

f6_top <- f6_fg[order(-abs(deg_logFC) * pp4_best)][seq_len(min(15, .N))]

p6 <- ggplot(deg, aes(x = logFC, y = pp4_best, color = coarse_class)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
  # Background: non-COLOC DEGs (PP4 < 0.5) — light gray
  geom_point(data = deg[pp4_best < 0.5], inherit.aes = FALSE,
             aes(x = logFC, y = pp4_best),
             color = "gray80", size = 0.4, alpha = 0.4, shape = 16) +
  # Foreground: COLOC DEGs colored by coarse class
  geom_point(data = deg[pp4_best >= 0.5],
             aes(color = coarse_class), size = 2.2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = merge(f6_top[, .(symbol = gene_symbol)],
                 deg, by = "symbol")[, .(symbol, x = logFC, y = pp4_best)],
    aes(x = x, y = y, label = symbol), inherit.aes = FALSE,
    size = 2.4, max.overlaps = Inf,
    force = 3, force_pull = 1,
    box.padding = 0.3, point.padding = 0.15,
    segment.size = 0.2, segment.color = "gray55",
    min.segment.length = 0, fontface = "italic",
    bg.color = "white", bg.r = 0.12,
    show.legend = FALSE,
    max.iter = 8000, seed = 1) +
  scale_color_manual(values = c("coding" = "#C2185B",
                                 "non-coding" = "#42A5F5",
                                 "no COLOC" = "gray70"),
                     name = "Variant class",
                     breaks = c("coding","non-coding")) +
  scale_y_continuous(limits = c(-0.02, 1.05),
                     breaks = c(0, 0.5, 0.9, 1.0),
                     labels = c("0", "0.5", "0.9", "1")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = expression("DEG log"[2]*"FC (Tier 1: padj<0.05, |logFC|>0.5)"),
       y = "Best COLOC PP.H4",
       title = "DEG × COLOC by variant class") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 2.5), nrow = 1))
save_panel(p6, "panel_F6_lfc_vs_pp4.pdf", w = 8, h = 6.5)

# ---- F7: Canonical examples -----------------------------------------------
f7 <- fread(file.path(DATDIR, "panel_F7_canonical.csv"))
f7[, gene_symbol := factor(gene_symbol, levels = rev(unique(gene_symbol)))]
f7[, fine_class := factor(fine_class, levels = fine_levels)]
p7 <- ggplot(f7, aes(x = pp4_best, y = gene_symbol)) +
  geom_segment(aes(xend = 0, yend = gene_symbol), color = "gray80",
               linewidth = 0.4, na.rm = TRUE) +
  geom_point(aes(color = fine_class, size = pmax(abs(deg_logFC), 0.3,
                                                  na.rm = TRUE)),
             alpha = 0.9, na.rm = TRUE) +
  geom_text(data = f7[is.na(pp4_best)], aes(x = 0, label = "n/a"),
            hjust = -0.2, size = 3, color = "gray50", fontface = "italic") +
  scale_color_manual(values = fine_pal, drop = FALSE, na.value = "gray70") +
  scale_size_continuous(range = c(2, 7), name = "|logFC|") +
  scale_x_continuous(limits = c(0, 1.05)) +
  labs(x = "Best COLOC PP.H4", y = NULL, color = "Fine class",
       title = "F7 — Canonical MASLD examples (coding vs non-coding)",
       subtitle = "PNPLA3 / TM6SF2 = expected coding; RORA / HKDC1 = expected intron") +
  theme_masld() +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray30"))
save_panel(p7, "panel_F7_canonical_examples.pdf", w = 8, h = 5)

# ---- F8: CS-member coding fraction ----------------------------------------
f8_path <- file.path(DATDIR, "panel_F8_cs_coding_fraction.csv")
if (file.exists(f8_path)) {
  f8 <- fread(f8_path)
  f8[, has_DEG_gene := !is.na(coloc_genes)]   # all rows have coloc_genes from C
  p8 <- ggplot(f8, aes(x = pip_weighted_coding)) +
    geom_histogram(bins = 25, fill = "#42A5F5", color = "white",
                   linewidth = 0.3) +
    facet_wrap(~ "All COLOC-positive credible sets") +
    labs(x = "PIP-weighted coding fraction within CS", y = "# credible sets",
         title = "F8 — Credible-set coding-fraction (uncertainty)",
         subtitle = sprintf("%d CS-pair credible sets across %d studies",
                            nrow(f8), uniqueN(f8$study))) +
    theme_masld() +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, color = "gray30"))
  save_panel(p8, "panel_F8_cs_coding_fraction.pdf", w = 8, h = 4.5)
} else {
  cat("F8 data file not found; skipping\n")
}

cat("\n[04_render_figures] DONE\n")
