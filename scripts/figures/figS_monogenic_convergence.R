#!/usr/bin/env Rscript
# figS_monogenic_convergence.R
# Common-variant SuSiE-COLOC × rare-variant monogenic (Okur et al. 2026) convergence.
#
# 4-panel supplementary figure:
#   (a) Forest plot: panel-gene best PP.H4 across any GWAS, ranked descending.
#   (b) DEG signal for panel genes (dream LFC vs -log10 padj, panel genes highlighted).
#   (c) Complementarity: scatter of all atlas genes (max-PP4 vs |LFC|), panel genes overlaid.
#   (d) MAFLD-relevant subset: 4 genes with P/LP variants in the n=25 MAFLD/MASH subset
#       (APOB, PPARG, LMNA, NR0B2) — atlas evidence axes.
#
# Output: figures/supplementary/figS04_coloc/figS_monogenic_convergence.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(PROJ, "scripts/figures/load_figure_data.R"))
source(file.path(PROJ, "scripts/figures/publication_theme.R"))

ATLAS_PATH    <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
PANEL_PATH    <- file.path(PROJ, "data/external/okur2026_monogenic/okur2026_panel_genes_partial.tsv")
VARIANTS_PATH <- file.path(PROJ, "data/external/okur2026_monogenic/okur2026_pl_variants_pdf.tsv")
OUT_PDF       <- file.path(FIGS04_DIR, "figS_monogenic_convergence.pdf")

cat("Loading...\n")
atlas    <- fread(ATLAS_PATH)
panel    <- fread(PANEL_PATH, fill = TRUE)
variants <- fread(VARIANTS_PATH, fill = TRUE)

panel_genes <- unique(panel$gene)
pl_genes    <- unique(variants$gene)
mafld_genes <- unique(variants[grepl("MAFLD|MASH", indication, ignore.case=TRUE), gene])

# Atlas slice for panel genes
hits <- atlas[human_symbol %in% panel_genes]
hits[, max_pp4 := pmax(coloc_susie_best_pp4, coloc_abf_best_pp4, na.rm = TRUE)]
hits[is.na(max_pp4) | is.infinite(max_pp4), max_pp4 := 0]
hits[, has_pl  := human_symbol %in% pl_genes]
hits[, has_mafld := human_symbol %in% mafld_genes]
hits[, neg_log_padj := -log10(pmax(bulk_padj, 1e-50))]

cat("  Panel genes in atlas:", nrow(hits), "/", length(panel_genes), "\n")
cat("  P/LP-cohort genes in atlas:", sum(hits$has_pl), "\n")
cat("  MAFLD-flagged genes in atlas:", sum(hits$has_mafld), "\n")

# ===========================================================================
# Panel (a) — Forest plot of best PP.H4 across panel genes
# ===========================================================================
fa_dat <- hits[order(-max_pp4)]
fa_dat[, gene_label := factor(human_symbol, levels = rev(human_symbol))]
# Annotate which atlas-source contributes the max
fa_dat[, source_lab := dplyr::case_when(
  !is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 == max_pp4 & max_pp4 > 0 ~ "SuSiE",
  !is.na(coloc_abf_best_pp4)   & coloc_abf_best_pp4   == max_pp4 & max_pp4 > 0 ~ "ABF",
  TRUE ~ "below_thresh")]

p_a <- ggplot(fa_dat, aes(x = max_pp4, y = gene_label, color = has_pl)) +
  geom_segment(aes(x = 0, xend = max_pp4, yend = gene_label), linewidth = 0.5) +
  geom_point(aes(size = has_mafld)) +
  geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed",
             color = c("grey50", "grey25"), linewidth = 0.3) +
  scale_x_continuous("Best PP.H4 (any GWAS, ABF or SuSiE)",
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual("P/LP variant\nin cohort",
                     values = c("TRUE" = "#D62728", "FALSE" = "grey50"),
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  scale_size_manual("MAFLD/MASH\nP/LP",
                    values = c("TRUE" = 3, "FALSE" = 1.5),
                    labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  labs(y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "right")

# ===========================================================================
# Panel (b) — DEG signal for panel genes (volcano-style with labels)
# ===========================================================================
fb_dat <- atlas[, .(human_symbol, bulk_logFC, bulk_padj)]
fb_dat[, in_panel := human_symbol %in% panel_genes]
fb_dat[, neg_log_padj := -log10(pmax(bulk_padj, 1e-50))]
fb_dat <- fb_dat[!is.na(bulk_logFC) & !is.na(bulk_padj)]

# Label panel genes
fb_label <- fb_dat[in_panel & bulk_padj < 0.05]

p_b <- ggplot() +
  geom_point(data = fb_dat[in_panel == FALSE],
             aes(x = bulk_logFC, y = neg_log_padj),
             color = "grey80", alpha = 0.3, size = 0.4) +
  geom_point(data = fb_dat[in_panel == TRUE],
             aes(x = bulk_logFC, y = neg_log_padj),
             color = "#D62728", size = 1.6) +
  ggrepel::geom_text_repel(data = fb_label,
             aes(x = bulk_logFC, y = neg_log_padj, label = human_symbol),
             size = GEOM_TEXT_6PT, max.overlaps = 30, box.padding = 0.3) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  labs(x = "dream logFC (NASH vs control)",
       y = "-log10 dream padj") +
  theme_masld(base_size = 6) +
  coord_cartesian(xlim = c(-3, 3))

# ===========================================================================
# Panel (c) — Complementarity scatter (max-PP4 vs |LFC|, panel genes overlaid)
# ===========================================================================
fc_dat <- atlas[, .(human_symbol, coloc_susie_best_pp4, coloc_abf_best_pp4,
                    bulk_logFC, bulk_padj)]
fc_dat[, max_pp4 := pmax(coloc_susie_best_pp4, coloc_abf_best_pp4, na.rm = TRUE)]
fc_dat[is.na(max_pp4) | is.infinite(max_pp4), max_pp4 := 0]
fc_dat[, abs_lfc := abs(bulk_logFC)]
fc_dat[, in_panel := human_symbol %in% panel_genes]
fc_dat <- fc_dat[!is.na(abs_lfc)]

p_c <- ggplot() +
  geom_point(data = fc_dat[in_panel == FALSE],
             aes(x = max_pp4, y = abs_lfc),
             color = "grey80", alpha = 0.3, size = 0.4) +
  geom_point(data = fc_dat[in_panel == TRUE & max_pp4 > 0.05],
             aes(x = max_pp4, y = abs_lfc),
             color = "#D62728", size = 1.6) +
  ggrepel::geom_text_repel(data = fc_dat[in_panel == TRUE & (max_pp4 > 0.2 | abs_lfc > 0.5)],
             aes(x = max_pp4, y = abs_lfc, label = human_symbol),
             size = GEOM_TEXT_6PT, max.overlaps = 30, box.padding = 0.3) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "Max PP.H4 (any GWAS)",
       y = "|dream logFC| (NASH vs control)") +
  theme_masld(base_size = 6) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 3.5))

# ===========================================================================
# Panel (d) — MAFLD-subset evidence summary (4 genes from n=25 MAFLD subset)
# ===========================================================================
fd_genes <- mafld_genes
fd_dat <- atlas[human_symbol %in% fd_genes,
                .(human_symbol,
                  best_PP4 = pmax(coloc_susie_best_pp4, coloc_abf_best_pp4, na.rm=TRUE),
                  liver_enz_PP4 = best_liver_enzyme_pp4,
                  finngen_nafld_PP4 = finngen_nafld_coloc_pp4,
                  dream_LFC = bulk_logFC,
                  bulk_padj = bulk_padj)]
fd_long <- melt(fd_dat, id.vars = "human_symbol",
                measure.vars = c("best_PP4", "liver_enz_PP4", "finngen_nafld_PP4"),
                variable.name = "evidence", value.name = "PP4")
fd_long[, evidence := factor(evidence, levels = c("best_PP4", "liver_enz_PP4", "finngen_nafld_PP4"),
                              labels = c("Best PP.H4", "Liver enzyme PP.H4", "FinnGen NAFLD PP.H4"))]

p_d <- ggplot(fd_long, aes(x = evidence, y = human_symbol, fill = PP4)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", PP4)), size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_gradient2(low = "white", mid = "#FFE5B4", high = "#D62728",
                       midpoint = 0.25, limits = c(0, 1), na.value = "grey85",
                       name = "PP.H4") +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# ===========================================================================
# Compose
# ===========================================================================
composite <- (p_a | p_b) / (p_c | p_d)

message("[caption] Figure S | Convergence of common-variant SuSiE-COLOC and rare-variant monogenic evidence (Okur et al. 2026). ",
        sprintf("Partial 90-gene panel (n=%d genes from PDF text + Table 3); full panel pending Zenodo embargo (2026-05-01). ",
                length(panel_genes)),
        "(a) Common-variant COLOC across panel genes; (b) DEG signal for panel genes; ",
        "(c) rare/common evidence complementarity; (d) MAFLD/MASH subset evidence summary.")

dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PDF, composite, width = fig_full_width, height = 5.91, device = cairo_pdf)
cat("Wrote:", OUT_PDF, "\n")
