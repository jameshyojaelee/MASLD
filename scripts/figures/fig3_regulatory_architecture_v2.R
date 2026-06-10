#!/usr/bin/env Rscript
##############################################################################
# fig3_regulatory_architecture_v2.R  (2026-04-15; SuSiE-COLOC refresh 2026-04-21;
#                                     2026-04-29 layout cascade)
# Figure 3 | Multi-ancestry regulatory architecture
#   Panels (binding: docs/manuscript/FIGURE_PLAN_REVISED.md §Figure 3):
#     3a  Twin Manhattan: max PIP (20 validated GWAS, top) + SuSiE-COLOC PP.H4
#           (20 validated GWAS, bottom). The full SuSiE-COLOC portfolio is 23
#           (14 EUR + 3 EAS + 3 AFR + 3 SAS Pan-UKBB); 20 validated shown here
#           (14 EUR + 3 BBJ EAS). AFR/SAS Pan-UKBB are exploratory (see figS09).
#           Canonical counts (2026-04-21 rebuild, gene_level_coloc.csv):
#             SuSiE PP.H4.susie > 0.5: 368 (289 > 0.8, 210 > 0.9)
#             ABF fallback > 0.5: 618 (282 > 0.8, 186 > 0.9)
#     3b  Per-ancestry COLOC eGene counts (23-GWAS portfolio): grouped bars
#           at PP.H4>0.5/0.8/0.9, prefer SuSiE then ABF fallback per gene.
#           Source: GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
#     3c  RORA / GGT chr15:60883281 cross-ancestry locus LD-zoom (external panel
#           rendered by run_fig3c_rora_locus_zoom.sh; figS09_locus_zoom.R writes
#           directly to panels/fig3c.pdf).
#     3d  GWAS-ATAC PIP vs |alleleDiff| scatter, SCENIC+ disease-regulon TFs
#           highlighted (was 3c pre-2026-04-29).
#     3e  High-PIP (>=0.8) variant × TF heatmap (was 3d pre-2026-04-29).
#     3f  Drug-target genetic validation scatter: bulk_logFC vs best SuSiE PP4
#           across 23 GWAS; 8 clinical anchors + 11 novel druggable genes labelled.
#           (Was 3e pre-2026-04-29.)
#
#   Demoted to figS_therapeutics 2026-04-29 (was 3f):
#     drug-target finemapping scatter (bulk_logFC vs SuSiE-X max PIP)
#     -> figures/supplementary/figS_therapeutics/panels/fig3_drug_finemapping_demoted.pdf
#
#   Outputs: figures/main/fig3_regulatory_architecture/panels/fig3{a,b,d,e,f}.pdf
#            figures/main/fig3_regulatory_architecture/fig3_regulatory_architecture.pdf
#            (fig3c.pdf written by external RORA launcher)
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
cat("[fig3] PANEL_DIR:", PANEL_DIR, "\n")

save_panel <- function(p, name, width = fig_half_width, height = 3.2) {
  out <- file.path(PANEL_DIR, name)
  save_fig(p, out, width = width, height = height)
  cat("[fig3] Saved:", out, "\n")
  invisible(p)
}

# ---------------------------------------------------------------------------
# GWAS classification (20 VALIDATED target GWAS: 14 EUR + 3 BBJ EAS).
#   The full SuSiE-COLOC portfolio as of 2026-04-21 is 23 GWAS:
#     14 EUR  (UKBB ALT/AST/GGT, 6 NAFLD cohorts, 3 PDFF, 2 FinnGen R12
#               NAFLD/NASH), fine-mapped against UKBB EUR LD (N = 337K)
#     3  EAS  (BBJ ALT/AST/GGT), 1000 Genomes EAS LD
#     3  AFR  (Pan-UKBB ALT/AST/GGT N ≈ 6.6K), 1000 Genomes AFR LD
#     3  SAS  (Pan-UKBB CSA ALT/AST/GGT N ≈ 8.9K), 1000 Genomes SAS LD
#   The main Fig 3 shows the 20 validated subset; AFR/SAS exploratory go to figS09.
#   EUR_14 contains the 2 FinnGen R12 studies. The 32514122_*_EAS naming is
#   historical; these ARE BBJ data. No Korean GWAS are publicly available.
#   (2026-06-07 portfolio refactor: dropped all cirrhosis/HCC GWAS +
#    PanUKBB sex-stratified EUR duplicates; 34 -> 23 studies.)
# ---------------------------------------------------------------------------
EUR_17 <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2021_34128465_PDFF_EUR",  "2021_34841290_NAFLD_EUR",
  "2021_34957434_PDFF_EUR",  "2022_36402844_PDFF_EUR",
  "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR", "FinnGen_NAFLD", "FinnGen_NASH",
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT"
)
FINNGEN_3       <- c("FinnGen_NAFLD", "FinnGen_NASH")
BBJ_ENZYME_3    <- c("BBJ_ALT", "BBJ_AST", "BBJ_GGT")
BBJ_5           <- BBJ_ENZYME_3

gwas_category <- function(study) {
  fcase(
    study %in% BBJ_ENZYME_3,  "BBJ liver enzyme (EAS)",
    study %in% FINNGEN_3, "FinnGen R12",
    grepl("NAFLD|NASH|Cirrhosis|HCC", study), "EUR disease",
    grepl("PDFF", study), "EUR PDFF",
    grepl("UKBB_(ALT|AST|GGT)$", study), "EUR liver enzyme",
    default = "Other"
  )
}

cat_colors <- c(
  "EUR disease"              = "#C9265E",
  "EUR liver enzyme"         = "#880E4F",
  "EUR PDFF"                 = "#7B1FA2",
  "FinnGen R12"              = "#1565C0",
  "BBJ liver enzyme (EAS)"   = "#F4511E"
)

# ===========================================================================
# Panel 3a: Twin Manhattan — max PIP (top) + max SuSiE-COLOC PP.H4 (bottom) across 20 GWAS
#   x = genome position (chr, bp); y = max PIP at locus (SuSiE/CARMA recommended)
#   color = GWAS category; points = 1 per (study, locus)
# ===========================================================================
cat("[fig3] Panel 3a: Locus-PIP Manhattan ...\n")

fmap <- fread(file.path(BASE,
  "GWAS/finemapping/results/combined_finemapping.csv"))
fmap[, gwas_cat := gwas_category(study)]
# Restrict to 20 target studies: 14 EUR (incl. 2 FinnGen) + 3 BBJ EAS (all Japanese)
target_22 <- c(EUR_17, BBJ_5)
# 23-GWAS COLOC refactor (2026-06-06) dropped cirrhosis/HCC GWAS; the displayed
# fine-mapping panel is the 17 MASLD/enzyme/PDFF studies: 14 EUR (incl. 2 FinnGen) + 3 BBJ EAS.
stopifnot(length(target_22) == 17L)
fmap <- fmap[study %in% target_22]
# Top variant per study × locus
top_loci <- fmap[, .SD[which.max(max_pip)], by = .(study, locus, chromosome)]
top_loci <- top_loci[!is.na(max_pip) & max_pip > 0]

# Compute cumulative genomic position for plotting
chrom_lengths <- data.table(
  chrom = 1:22,
  len = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
          159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
          114364328, 107043718, 101991189, 90338345,  83257441,  80373285,
          58617616,  64444167,  46709983,  50818468)
)
chrom_lengths[, cum_start := cumsum(c(0, len[-length(len)]))]
top_loci <- merge(top_loci,
                  chrom_lengths[, .(chromosome = chrom, cum_start, len)],
                  by = "chromosome", all.x = TRUE)
top_loci[, x_cum := cum_start + position]

# Chromosome tick positions
chrom_ticks <- chrom_lengths[, .(chrom, mid = cum_start + len / 2)]

# Sanity counts for caption
n_total_loci <- nrow(top_loci)
n_studies    <- length(unique(top_loci$study))
cat(sprintf("  [3a] %d locus-study points, %d studies\n", n_total_loci, n_studies))

# Highlight top genes/loci by PIP (anchor gene names from COLOC)
coloc_anno <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc_anno <- coloc_anno[gene != "" & !is.na(gene)]
# Well-known MASLD loci to label
label_coords <- data.table(
  gene    = c("PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "THRB", "MARC1", "MBOAT7", "APOE", "HNF1A"),
  chrom   = c(22, 19, 4, 2, 3, 1, 19, 19, 12),
  pos_hg19 = c(44324727, 19379549, 88231392, 27730940, 24174589, 220974135,
               54677189, 45411941, 121416650)
)
label_coords <- merge(label_coords,
                     chrom_lengths[, .(chrom, cum_start)], by = "chrom")
label_coords[, x_cum := cum_start + pos_hg19]
# Match to top_loci by nearest position within 1Mb; use max_pip for y
label_coords[, y := sapply(seq_len(.N), function(i) {
  sub <- top_loci[chromosome == label_coords$chrom[i] &
                  abs(position - label_coords$pos_hg19[i]) < 1e6]
  if (nrow(sub) == 0) return(NA_real_)
  max(sub$max_pip, na.rm = TRUE)
})]
label_coords <- label_coords[!is.na(y) & y > 0.3]

# Build PP.H4 track: parse top_snp "chr:pos" to genomic coords, keep max PP.H4
# per (study, chromosome, ~locus) so the lower rail mirrors the upper rail.
sc_all <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
sc_all[, study := gwas_name]
sc_all <- sc_all[study %in% target_22 &
                 !is.na(PP.H4.abf) &
                 !is.na(top_snp) & top_snp != ""]
sc_all[, c("chromosome", "position") := tstrsplit(top_snp, ":", fixed = TRUE,
                                                  type.convert = TRUE)]
sc_all <- sc_all[!is.na(chromosome) & chromosome %in% 1:22 & !is.na(position)]
sc_all[, gwas_cat := gwas_category(study)]
# Collapse to one point per (study, chromosome, ~500kb bin): best PP.H4
sc_all[, locus_bin := floor(position / 5e5)]
top_h4 <- sc_all[, .SD[which.max(PP.H4.abf)],
                 by = .(study, chromosome, locus_bin)]
top_h4 <- merge(top_h4,
                chrom_lengths[, .(chromosome = chrom, cum_start)],
                by = "chromosome", all.x = TRUE)
top_h4[, x_cum := cum_start + position]
n_h4_points <- nrow(top_h4)
cat(sprintf("  [3a] %d PP.H4 points (upper rail: %d loci)\n",
            n_h4_points, n_total_loci))

# Labeling y-coord for PP.H4 rail (reuse label_coords genes)
label_h4 <- copy(label_coords)
label_h4[, y_h4 := sapply(seq_len(.N), function(i) {
  sub <- top_h4[chromosome == label_h4$chrom[i] &
                abs(position - label_h4$pos_hg19[i]) < 1e6]
  if (!nrow(sub)) return(NA_real_)
  max(sub$PP.H4.abf, na.rm = TRUE)
})]
label_h4 <- label_h4[!is.na(y_h4) & y_h4 > 0.3]

common_x <- list(
  scale_color_manual(values = cat_colors, name = "GWAS"),
  scale_x_continuous(breaks = chrom_ticks$mid, labels = chrom_ticks$chrom,
                     expand = expansion(mult = c(0.005, 0.005))),
  scale_y_continuous(limits = c(0, 1.02),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)),
  theme_masld(),
  theme(legend.position = "top",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
)

p3a_top <- ggplot(top_loci, aes(x = x_cum, y = max_pip, color = gwas_cat)) +
  rasterize_layer(geom_point(size = 0.7, alpha = 0.75, stroke = 0)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray40") +
  geom_label_repel(
    data = label_coords,
    aes(x = x_cum, y = y, label = gene),
    inherit.aes = FALSE,
    size = 2.2, fontface = "italic",
    label.size = 0.25, label.padding = unit(0.1, "lines"),
    box.padding = 0.35, point.padding = 0.2,
    segment.size = 0.5, segment.color = "black",
    min.segment.length = 0,
    max.overlaps = 20, seed = 42,
    fill = "white"
  ) +
  common_x +
  labs(y = "Max fine-map PIP (20 ancestry-matched GWAS)",
       title = "Cross-ancestry fine-mapping + COLOC (14 EUR + 3 EAS SuSiE-COLOC; 6 Pan-UKBB ABF-only)") +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

p3a_bot <- ggplot(top_h4, aes(x = x_cum, y = PP.H4.abf, color = gwas_cat)) +
  rasterize_layer(geom_point(size = 0.7, alpha = 0.75, stroke = 0)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray40") +
  geom_label_repel(
    data = label_h4,
    aes(x = x_cum, y = y_h4, label = gene),
    inherit.aes = FALSE,
    size = 2.2, fontface = "italic",
    label.size = 0.25, label.padding = unit(0.1, "lines"),
    box.padding = 0.35, point.padding = 0.2,
    segment.size = 0.5, segment.color = "black",
    min.segment.length = 0,
    max.overlaps = 20, seed = 42,
    fill = "white"
  ) +
  common_x +
  labs(x = "Chromosome", y = "Max COLOC PP.H4 (20 ancestry-matched GWAS)") +
  theme(axis.text.x = element_text(size = 5),
        legend.position = "none")

p3a <- p3a_top / p3a_bot +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "top")

save_panel(p3a, "fig3a.pdf", width = fig_full_width, height = 5.2)

# ===========================================================================
# Panel 3b: Per-ancestry COLOC eGene counts at PP.H4 > 0.5 / 0.8 / 0.9
#   23-GWAS portfolio split: 14 EUR + 3 EAS + 3 AFR + 3 SAS.
#   Per gene per ancestry: take best-available PP4 (prefer SuSiE, fall back to
#   ABF when SuSiE did not converge). Bars grouped by ancestry × threshold.
#   Sidecar CSV preserves caption-ready counts.
# ===========================================================================
cat("[fig3] Panel 3b: Per-ancestry COLOC eGene counts ...\n")

sc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))

AFR_3 <- c("PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT")
SAS_3 <- c("PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT")

ancestry_for_gwas <- function(g) {
  fcase(
    g %in% EUR_17, "EUR",
    g %in% BBJ_5,  "EAS",
    g %in% AFR_3,  "AFR",
    g %in% SAS_3,  "SAS",
    default       = NA_character_
  )
}
sc[, ancestry := ancestry_for_gwas(gwas_name)]
sc <- sc[!is.na(ancestry)]

# Best-available PP4: prefer SuSiE, fall back to ABF
sc[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
sc[is.infinite(pp4_best), pp4_best := NA_real_]

# Per-gene-per-ancestry max PP4 across that ancestry's GWAS
per_gene_anc <- sc[!is.na(pp4_best),
                   .(max_pp4 = max(pp4_best, na.rm = TRUE)),
                   by = .(gene, ancestry)]

thresholds <- c(0.5, 0.8, 0.9)
ancestry_levels <- c("EUR", "EAS", "AFR", "SAS")
ancestry_n_gwas <- c(EUR = 14L, EAS = 3L, AFR = 3L, SAS = 3L)

bar_dt <- rbindlist(lapply(thresholds, function(thr) {
  per_gene_anc[, .(n_genes = sum(max_pp4 > thr)),
               by = ancestry][, threshold := thr]
}))
bar_dt[, ancestry := factor(ancestry, levels = ancestry_levels)]
bar_dt[, threshold_label := factor(sprintf("PP4 > %.1f", threshold),
                                    levels = sprintf("PP4 > %.1f", thresholds))]
bar_dt[, n_gwas := ancestry_n_gwas[as.character(ancestry)]]
bar_dt[, ancestry_label := sprintf("%s (n=%d GWAS)", ancestry, n_gwas)]
bar_dt[, ancestry_label := factor(ancestry_label,
                                   levels = sprintf("%s (n=%d GWAS)",
                                                    ancestry_levels,
                                                    ancestry_n_gwas[ancestry_levels]))]

ancestry_colors_3b <- c(
  "EUR" = "#C9265E",  # magenta (matches EUR disease in 3a)
  "EAS" = "#F4511E",  # orange (matches BBJ liver enzyme in 3a)
  "AFR" = "#00695C",  # teal
  "SAS" = "#7B1FA2"   # purple
)
threshold_alphas <- c("PP4 > 0.5" = 0.45,
                       "PP4 > 0.8" = 0.75,
                       "PP4 > 0.9" = 1.00)

p3b <- ggplot(bar_dt, aes(x = ancestry_label, y = n_genes,
                           fill = ancestry, alpha = threshold_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75,
           color = "white", linewidth = 0.2) +
  geom_text(aes(label = n_genes),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 1.9, color = "gray25") +
  scale_fill_manual(values = ancestry_colors_3b, guide = "none") +
  scale_alpha_manual(values = threshold_alphas, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL,
       y = "Colocalised eGenes (best PP.H4 across ancestry's GWAS)",
       title = "Per-ancestry COLOC support",
       subtitle = "Best SuSiE PP.H4 per gene; ABF fallback when SuSiE did not converge") +
  theme_masld() +
  theme(legend.position = "top",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6, color = "gray35"),
        axis.text.x = element_text(size = 6.5))

save_panel(p3b, "fig3b.pdf", width = fig_half_width, height = 3.0)
fwrite(bar_dt[, .(ancestry, n_gwas, threshold, n_genes)],
       file.path(FIG3_DIR, "fig3b_ancestry_coloc_counts.csv"))

# Caption sidecar (kept for backward compatibility with figS09 / NUMBERS.md)
sc_eur <- sc[ancestry == "EUR"]
per_gene_eur_abf <- sc_eur[!is.na(PP.H4.abf),
                           .(max_pp4 = max(PP.H4.abf)),
                           by = gene][max_pp4 > 0.1]
fwrite(data.table(threshold = thresholds,
                  n_genes   = sapply(thresholds, function(t) sum(per_gene_eur_abf$max_pp4 > t)),
                  n_gwas    = 14L,
                  ancestry  = "EUR",
                  method    = "ABF",
                  note      = "EUR-only ABF subset for legacy caption; full per-ancestry counts in fig3b_ancestry_coloc_counts.csv"),
       file.path(FIG3_DIR, "fig3a_pph4_thresholds.csv"))

# ===========================================================================
# Supplementary inset (writes to figS09): Cross-ancestry locus zoom at
#   chr22:24334710 (GGT1).  4-panel stack: UKBB_GGT GWAS, BBJ_GGT GWAS,
#   MESuSiE PIP_EUR, PIP_EAS.  EUR lead 24334710, EAS lead 24463564 (~129 kb
#   apart).  Slot 3c in main Fig 3 is now the RORA / GGT chr15 LD-zoom written
#   externally by run_fig3c_rora_locus_zoom.sh.
# ===========================================================================
cat("[fig3] figS09 inset: Cross-ancestry locus zoom (chr22 GGT1) ...\n")

LZ_LOCUS <- "locus_GGT_chr22_24334710"
FM_BASE  <- file.path(BASE, "GWAS/finemapping")
loci_tbl <- fread(file.path(FM_BASE, "results/susiex/shared_loci.csv"))
lz <- loci_tbl[locus_id == LZ_LOCUS][1]
stopifnot(nrow(lz) == 1L)
LZ_CHR   <- lz$chr
LZ_WSTART <- lz$window_start
LZ_WEND   <- lz$window_end
EUR_LEAD <- lz$eur_lead_pos
EAS_LEAD <- lz$eas_lead_pos
cat(sprintf("  [3c] chr%s:%s-%s | EUR lead=%s | EAS lead=%s\n",
            LZ_CHR, LZ_WSTART, LZ_WEND, EUR_LEAD, EAS_LEAD))

lz_eur_ss <- fread(file.path(FM_BASE, "data/sumstats/UKBB_GGT_reformatted_hg19.tsv"),
                   select = c("chromosome", "position", "pval"))
lz_eur_ss <- lz_eur_ss[chromosome == LZ_CHR & position >= LZ_WSTART & position <= LZ_WEND]
lz_eur_ss[, mlog10p := -log10(pmax(as.numeric(pval), 1e-300))]

lz_eas_ss <- fread(file.path(FM_BASE, "data/sumstats/BBJ_GGT_reformatted_hg19.tsv"),
                   select = c("chromosome", "position", "pval"))
lz_eas_ss <- lz_eas_ss[chromosome == LZ_CHR & position >= LZ_WSTART & position <= LZ_WEND]
lz_eas_ss[, mlog10p := -log10(pmax(as.numeric(pval), 1e-300))]

lz_pips <- fread(file.path(FM_BASE, "results/mesusie/mesusie_variant_summary.csv"))
lz_pips <- lz_pips[locus_id == LZ_LOCUS]

GWS_LINE <- -log10(5e-8)
COL_EUR  <- "#1565C0"   # match cat_colors EUR-FinnGen
COL_EAS  <- "#F4511E"   # match BBJ color

lz_common <- list(
  scale_x_continuous(limits = c(LZ_WSTART, LZ_WEND), expand = c(0, 0)),
  theme_masld(base_size = 6),
  theme(panel.grid.minor = element_blank(),
        plot.title  = element_blank(),
        plot.margin = margin(t = 1, r = 6, b = 1, l = 2))
)
lz_no_x <- theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())

p3c_1 <- ggplot(lz_eur_ss, aes(position, mlog10p)) +
  rasterize_layer(geom_point(color = "grey80", size = 0.5, alpha = 0.5), dpi = 500) +
  geom_hline(yintercept = GWS_LINE, linetype = "longdash", linewidth = 0.3, color = "black") +
  geom_point(data = lz_eur_ss[position == EUR_LEAD], fill = COL_EUR,
             color = "black", size = 2, shape = 21, stroke = 0.35) +
  labs(y = expression(UKBB~GGT~-log[10](italic(p)))) +
  lz_common + lz_no_x

p3c_2 <- ggplot(lz_eas_ss, aes(position, mlog10p)) +
  rasterize_layer(geom_point(color = "grey80", size = 0.5, alpha = 0.5), dpi = 500) +
  geom_hline(yintercept = GWS_LINE, linetype = "longdash", linewidth = 0.3, color = "black") +
  geom_point(data = lz_eas_ss[position == EAS_LEAD], fill = COL_EAS,
             color = "black", size = 2, shape = 21, stroke = 0.35) +
  labs(y = expression(BBJ~GGT~-log[10](italic(p)))) +
  lz_common + lz_no_x

# Identify top-PIP variants (fine-mapped leads ≠ GWAS leads)
top_eur <- lz_pips[which.max(pip_eur)]
top_eas <- lz_pips[which.max(pip_eas)]
cat(sprintf("  [3c] Top PIP EUR: %.3f at %d | Top PIP EAS: %.3f at %d\n",
            top_eur$pip_eur, top_eur$pos, top_eas$pip_eas, top_eas$pos))

# Credible set coloring for PIP panels
lz_pips[, cs_col := fifelse(in_cs, fifelse(pip_eur > pip_eas, "EUR CS", "EAS CS"), "not in CS")]

p3c_3 <- ggplot(lz_pips, aes(pos, pip_eur)) +
  geom_segment(data = lz_pips[in_cs == TRUE],
               aes(x = pos, xend = pos, y = 0, yend = pip_eur),
               color = "grey70", linewidth = 0.3, alpha = 0.5) +
  geom_point(aes(color = cs_col), size = 1.2, alpha = 0.7) +
  geom_point(data = top_eur, fill = COL_EUR,
             color = "black", size = 2.5, shape = 21, stroke = 0.4) +
  geom_text(data = top_eur, aes(label = sprintf("PIP=%.2f", pip_eur)),
            hjust = -0.1, vjust = 0.5, size = 2, color = COL_EUR, fontface = "bold") +
  scale_color_manual(values = c("EUR CS" = COL_EUR, "EAS CS" = COL_EAS,
                                "not in CS" = "grey80"), guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1)) +
  labs(y = "PIP EUR") +
  lz_common + lz_no_x

p3c_4 <- ggplot(lz_pips, aes(pos, pip_eas)) +
  geom_segment(data = lz_pips[in_cs == TRUE],
               aes(x = pos, xend = pos, y = 0, yend = pip_eas),
               color = "grey70", linewidth = 0.3, alpha = 0.5) +
  geom_point(aes(color = cs_col), size = 1.2, alpha = 0.7) +
  geom_point(data = top_eas, fill = COL_EAS,
             color = "black", size = 2.5, shape = 21, stroke = 0.4) +
  geom_text(data = top_eas, aes(label = sprintf("PIP=%.2f", pip_eas)),
            hjust = -0.1, vjust = 0.5, size = 2, color = COL_EAS, fontface = "bold") +
  scale_color_manual(values = c("EUR CS" = COL_EUR, "EAS CS" = COL_EAS,
                                "not in CS" = "grey80"), guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1)) +
  labs(x = sprintf("chr%s position (hg19, bp)", LZ_CHR), y = "PIP EAS") +
  lz_common

p3c <- (p3c_1 / p3c_2 / p3c_3 / p3c_4) +
  plot_layout(heights = c(3, 3, 2, 2)) +
  plot_annotation(
    title    = "GGT1 locus (chr22:24.23–24.56 Mb)",
    subtitle = sprintf("Cross-ancestry fine-mapping · EUR lead %s bp, EAS lead %s bp (\u0394 = %d kb)",
                       format(EUR_LEAD, big.mark = ","),
                       format(EAS_LEAD, big.mark = ","),
                       abs(round((EUR_LEAD - EAS_LEAD) / 1e3))),
    theme    = theme(plot.title = element_text(size = 7, face = "bold"),
                     plot.subtitle = element_text(size = 5.5, color = "grey30"))
  )

# Save to figS09 (multi-ancestry supplement) instead of fig3 panels.
save_fig(p3c,
         file.path(FIGS09_DIR, "figS09_chr22_ggt1_cross_ancestry_locus.pdf"),
         width = fig_half_width, height = 4.0)
cat("[fig3] 3c (GGT1 locus) moved to figS09.\n")

# NOTE: the following panels have been retired from Fig 3:
#   - old 3c (cross-ancestry locus zoom) -> moved to figS09 (line above)
#   - old 3d (cross-ancestry elastic-net classifier) -> REMOVED (negative result)
#   - old 3e (disease-regulon TF bar chart) -> merged into new Panel 3d
#                                              (per-TF peak + motif bars)
# Cleanup stale panel files from prior runs (post-2026-04-29 cascade keeps
# fig3{a,b,c,d,e,f}.pdf; only legacy g/h slots are stale):
for (stale in c("fig3g.pdf", "fig3h.pdf")) {
  p <- file.path(PANEL_DIR, stale)
  if (file.exists(p)) {
    file.remove(p); cat("[fig3] Removed stale", stale, "\n")
  }
}

# NOTE: the old "disease-regulon TF peak + motif bars" panel has been dropped.
# Reason: motif_disruption_scores.csv only contains variant-TF rows where a motif
# match was detected in a peak, so "# in peak" and "# motif-disrupting" counts
# were identical per TF (same metric in two bars). The per-TF information is
# already captured by Panel 3d (scatter) and Panel 3e (heatmap) below.

# ===========================================================================
# Panels 3d + 3e: GWAS-ATAC motif disruption at MASLD disease-regulon TFs
#   3d — scatter of all PIP>=0.2 variant-TF pairs (max_pip vs |alleleDiff|);
#        SCENIC+ disease-regulon hits labeled as TF@gene
#   3e — high-PIP (>=0.8) variant × TF heatmap, collapsed to TFs that are
#        in the disease regulon OR hit by >=2 variants
#
#   Source: GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#           + gwas_atac_variant_annotation.csv
# ===========================================================================
cat("[fig3] Panels 3d + 3e: GWAS-ATAC motif disruption ...\n")


library(ggrepel)

ATAC_DIR <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
motif_raw <- fread(file.path(ATAC_DIR, "motif_disruption_scores.csv"))
ms <- unique(motif_raw[, .(SNP_id, tf_name, alleleDiff, max_pip,
                            motif_in_disease_regulon, activity_padj)])

va_all <- fread(file.path(ATAC_DIR, "gwas_atac_variant_annotation.csv"))
va_all[, has_gene := !is.na(nearest_gene) & nearest_gene != ""]
setorder(va_all, -has_gene, -max_pip)
va_uniq <- va_all[, .SD[1], by = variant_id][,
  .(variant_id, nearest_gene, distance_to_tss)]

merged <- merge(ms, va_uniq, by.x = "SNP_id", by.y = "variant_id", all.x = TRUE)
credible <- merged[max_pip >= 0.2]
credible[, locus      := gsub("^(\\d+):(\\d+):.+$", "chr\\1:\\2", SNP_id)]
credible[, gene_label := ifelse(!is.na(nearest_gene) & nearest_gene != "",
                                 nearest_gene, "intergenic")]
credible[, abs_diff   := abs(alleleDiff)]

REG_COL_3 <- "#C9265E"
OTH_COL_3 <- "#BDBDBD"

# --- Panel 3d (was 3c pre-2026-04-29): PIP vs disruption scatter -----------
# Label ONLY disease-regulon hits; dedup repeated labels.
credible[, label_3c := ""]
credible[motif_in_disease_regulon == TRUE,
         label_3c := paste0(tf_name, " @ ", gene_label)]
credible[duplicated(label_3c), label_3c := ""]

p3c <- ggplot(credible, aes(x = max_pip, y = abs_diff)) +
  geom_vline(xintercept = 0.8, linetype = "dashed",
             linewidth = 0.3, colour = "#888888") +
  geom_vline(xintercept = 0.5, linetype = "dotted",
             linewidth = 0.25, colour = "#AAAAAA") +
  geom_point(data = credible[motif_in_disease_regulon == FALSE],
             colour = OTH_COL_3, alpha = 0.55, shape = 16, size = 1.3) +
  geom_point(data = credible[motif_in_disease_regulon == TRUE],
             colour = REG_COL_3, alpha = 0.95, shape = 16, size = 2.2) +
  geom_text_repel(
    aes(label = label_3c),
    colour = REG_COL_3,
    size = 2.1, segment.size = 0.2, min.segment.length = 0.1,
    box.padding = 0.15, point.padding = 0.1,
    max.overlaps = Inf, seed = 42, force = 0.5
  ) +
  annotate("text", x = 0.8, y = max(credible$abs_diff) * 1.05,
           label = "PIP 0.8", hjust = -0.1, size = 2, colour = "#666666") +
  annotate("text", x = 0.5, y = max(credible$abs_diff) * 1.05,
           label = "PIP 0.5", hjust = -0.1, size = 2, colour = "#888888") +
  scale_x_continuous(limits = c(0.18, 1.05),
                     breaks = c(0.2, 0.5, 0.8, 1.0)) +
  labs(x = "GWAS max PIP", y = "|motif alleleDiff|",
       title = "Credible MASLD variants disrupt TF binding motifs",
       subtitle = sprintf(paste0(
         "%d pairs (PIP >= 0.2); magenta = SCENIC+ disease-regulon TFs ",
         "(%d pairs, %d TFs)"),
         nrow(credible),
         sum(credible$motif_in_disease_regulon),
         length(unique(credible[motif_in_disease_regulon==TRUE]$tf_name)))) +
  theme_masld() +
  theme(plot.margin = margin(4, 4, 4, 4))

save_panel(p3c, "fig3d.pdf", width = 5.0, height = 2.8)

# --- Panel 3e (was 3d pre-2026-04-29): high-PIP variant × TF heatmap -------
hi <- credible[max_pip >= 0.8]
tf_stats <- hi[, .(n_vars = uniqueN(SNP_id),
                    regulon = any(motif_in_disease_regulon)),
                by = tf_name]
keep_tfs <- tf_stats[n_vars >= 2 | regulon == TRUE]$tf_name
hi_f <- hi[tf_name %in% keep_tfs]

d_per_var <- unique(va_all[variant_id %in% unique(hi_f$SNP_id),
                            .(SNP_id = variant_id,
                              distance_to_tss = distance_to_tss)])
d_per_var <- d_per_var[, .(d_tss = min(distance_to_tss, na.rm = TRUE)),
                        by = SNP_id]
max_d_kb_3h <- round(max(d_per_var$d_tss, na.rm = TRUE) / 1000, 1)

hi_f[, var_label := paste0(locus, "  (", gene_label, ")")]
var_order <- hi_f[, .(mp = max(max_pip)), by = var_label][order(-mp)]$var_label
hi_f[, var_label := factor(var_label, levels = rev(var_order))]

tf_order_hi <- tf_stats[tf_name %in% keep_tfs
                        ][order(-regulon, -n_vars)]$tf_name
hi_f[, tf_name := factor(tf_name, levels = tf_order_hi)]

grid_all <- CJ(var_label = levels(hi_f$var_label),
               tf_name   = levels(hi_f$tf_name))
hi_grid <- merge(grid_all,
                 hi_f[, .(var_label, tf_name, alleleDiff,
                          motif_in_disease_regulon)],
                 by = c("var_label", "tf_name"), all.x = TRUE)
hi_grid[, var_label := factor(var_label, levels = levels(hi_f$var_label))]
hi_grid[, tf_name   := factor(tf_name,   levels = levels(hi_f$tf_name))]

regulon_tfs <- tf_stats[regulon == TRUE]$tf_name

p3d <- ggplot(hi_grid, aes(x = tf_name, y = var_label)) +
  # light grid backdrop with thin borders; keeps cells square-ish
  geom_tile(fill = "#F3F3F3", colour = "white", linewidth = 0.2) +
  geom_tile(data = hi_grid[!is.na(alleleDiff)],
            aes(fill = alleleDiff),
            colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#1565C0", mid = "#FFFFFF", high = "#C9265E",
    midpoint = 0, limits = c(-2.5, 2.5),
    breaks = c(-2, 0, 2), name = "alleleDiff"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL,
       y = sprintf("Variant (nearest gene within %.0f kb)", ceiling(max_d_kb_3h)),
       title = "High-confidence MASLD variants (PIP >= 0.8) disrupting TF motifs",
       subtitle = sprintf(paste0(
         "%d variants × %d TFs (disease-regulon + TFs hit by >=2 variants); ",
         "bold = SCENIC+ MASLD disease-regulon TF"),
         length(unique(hi_f$var_label)),
         length(unique(hi_f$tf_name)))) +
  theme_masld() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 0, vjust = 0,
                                face = ifelse(levels(hi_f$tf_name)
                                              %in% regulon_tfs,
                                              "bold", "plain"),
                                colour = ifelse(levels(hi_f$tf_name)
                                                %in% regulon_tfs,
                                                REG_COL_3, "black"),
                                size = 6.5),
    axis.text.y     = element_text(family = "mono", size = 5.8),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    legend.position = "right",
    legend.key.width  = unit(0.2, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.text     = element_text(size = 6),
    legend.title    = element_text(size = 6.5),
    plot.title      = element_text(hjust = 0, size = 8),
    plot.subtitle   = element_text(hjust = 0, size = 6,
                                    margin = margin(b = 2)),
    plot.title.position = "plot",
    plot.margin     = margin(4, 4, 4, 4)
  )

save_panel(p3d, "fig3e_tf_heatmap.pdf", width = 5.6, height = 3.0)

# ===========================================================================
# Panel 3f (+ demoted figS_therapeutics inset): Drug-target genetic validation
#   (moved 2026-04-29 from fig4_validation.R panel 4e/4e_pip; finemapping
#    scatter demoted to figS_therapeutics 2026-04-29)
#
#   Two-dimensional narrative panel (transcription x genetics) over a
#   19-gene universe: 8 clinical anchors + 11 novel druggable genes with
#   SuSiE PP4 >= 0.5 across the 23-GWAS portfolio.
#     - Failed PPARs + GLP1R + SCD cluster at PP4 ~ 0 (no genetic support)
#     - THRB (FDA) = sole genetically supported approved target
#     - 11 novels (HKDC1, GAS6, PKN3, CTSD, F2RL1, ST14, ALDH1B1, PKM,
#       VDR, ADH4, CHI3L1) sit at PP4 >= 0.5 with concordant DE.
#
#   Sources: RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#            RNA-seq/results/multi_evidence/convergence_evidence.csv
# ===========================================================================
cat("[fig3] Panel 3f + demoted finemapping scatter ...\n")

# --- Gene panel + clinical tiers ---
clinical_tiers <- data.table(
  gene = c("THRB", "NR1H4", "PPARA", "PPARG",
           "GLP1R", "SCD", "HSD17B13", "DGAT2"),
  drug = c("Resmetirom", "Obeticholic acid", "Elafibranor", "Lanifibranor",
           "Semaglutide", "Aramchol", "Rapirosiran", "ION224"),
  status = c("FDA approved", "Phase 3 (rejected)", "Phase 3 (failed)",
             "Phase 3", "FDA approved", "Phase 3",
             "Phase 3 planned", "Phase 2b")
)
clinical_tiers[, tier := fcase(
  grepl("FDA", status),     "FDA approved",
  grepl("Phase 3", status), "Phase 3",
  grepl("Phase 2", status), "Phase 2/earlier",
  default = "Other"
)]
# Three-category taxonomy (2026-05-03):
#   Drug target       = active clinical/preclinical drug program in MASLD
#   MASLD-associated  = appears in any published MASLD gene panel (Govaere
#                       2020 STM, Govaere 2020 Cell Metab fibrosis,
#                       SteatoSITE 15, Piras Key 30, Feng Top 27) or curated
#                       fig5 landmark set
#   Novel             = passes atlas filter but in neither prior list
#
# Atlas-driven anchor selection: all DEGs at PP4>=0.5, |logFC|>0.3, FDR<0.05
# are labeled (replaces hand-curated novel list).

drug_target_genes <- c("THRB","RORA","PNPLA3","NR1H4","PPARA",
                       "PPARG","GLP1R","SCD","HSD17B13","DGAT2")

# Published MASLD gene panels — canonical definitions in
# scripts/figures/figS_published_panel_benchmark.R. Inlined here so this
# panel renders without sourcing the benchmark script.
govaere_25 <- c("AKR1B10","DUSP6","GDF15","THBS2","A2M","CDH2",
                "COL1A1","COL3A1","COL4A1","COL4A2","COL6A3","DCN",
                "FBN1","FSTL1","IGFBP7","LUM","MFAP4","MMP2",
                "POSTN","SPARC","SPP1","TAGLN","THY1","TIMP1","VCAN")
govaere_fibrosis <- c("ACTA2","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2",
                      "CCN2","IL8","DCN","FN1","IL1B","IL6","LOX","LOXL2",
                      "MMP2","MMP9","PDGFRB","SERPINE1","TGFB1","TGFB2",
                      "TIMP1","TIMP2","TNF","VIM","CCN4")
steatosite_15 <- c("MT1F","KCNH7","COL25A1","RASD2","CCN2","STC1",
                   "GDNF","PRRX1","FGF7","LCNL1","DPEP1","CHRDL2",
                   "LHX6","POU4F1","CDH16")
piras_key <- c("HSD17B13","PNPLA3","TM6SF2","MARC1","MBOAT7",
               "CIDEB","GPAM","APOH","SERPINA6","GHR",
               "CYP7A1","ABCB4","SLC22A1","SLC27A5","ACSM3",
               "ADH4","ADH1B","CYP2E1","CYP4A11","CES1",
               "AKR1D1","SLC10A1","HGF","CRP","LCAT",
               "CETP","ALDOB","PCK1","GLS2","HAO1")
feng_top <- c("EFHD1","MLIP","TREM2","SPP1","GPNMB","CCL2",
              "CCL20","CXCL1","CXCL6","IL1B","IL32",
              "COL1A1","COL1A2","COL3A1","FN1","LOX","LOXL2",
              "ACTA2","PDGFRB","TGFB1","SERPINE1","MMP9",
              "AKR1B10","GDF15","THY1","THBS2","LUM")
landmark_masld_genes <- c("LPL","FABP4","MMP9","CFLAR")
masld_associated_genes <- unique(c(govaere_25, govaere_fibrosis,
                                   steatosite_15, piras_key, feng_top,
                                   landmark_masld_genes))

atlas_lite <- fread(
  file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "bulk_logFC", "bulk_padj",
             "coloc_susie_best_pp4", "susiex_max_pip"))
setnames(atlas_lite, "human_symbol", "gene")

bayev_lite <- fread(
  file.path(BASE, "RNA-seq/results/multi_evidence/convergence_evidence.csv"),
  select = c("human_symbol", "tier", "concordance_state",
             "convergence_score", "convergence_rank"))
setnames(bayev_lite, "human_symbol", "gene")
setnames(bayev_lite, "tier", "tier_46d")

# Atlas-filtered hits + always-shown landmarks (drug targets + curated MASLD).
# PP4>=0.5 entry, then per-band cap on Novels (see below) so labels distribute
# across the PP4 range without crowding at the top.
atlas_anchor_genes <- atlas_lite[
  !is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 >= 0.5 &
  !is.na(bulk_padj) & bulk_padj < 0.05 &
  !is.na(bulk_logFC) & abs(bulk_logFC) > 0.3, gene]
all_anchors <- unique(c(drug_target_genes, landmark_masld_genes,
                        atlas_anchor_genes))

gene_tbl <- atlas_lite[gene %in% all_anchors]
gene_tbl <- merge(gene_tbl, bayev_lite, by = "gene", all.x = TRUE)
gene_tbl[, coloc_pp4 := coloc_susie_best_pp4]
gene_tbl[is.na(coloc_pp4),   coloc_pp4   := 0]
gene_tbl[is.na(bulk_logFC), bulk_logFC := 0]

# Background: all atlas DEGs (FDR<0.05) shown as light-gray points so the
# anchors are read against the global cloud, not in isolation.
atlas_bg <- atlas_lite[!is.na(bulk_logFC) & !is.na(bulk_padj) &
                       bulk_padj < 0.05]
atlas_bg[, coloc_pp4 := coloc_susie_best_pp4]
atlas_bg[is.na(coloc_pp4), coloc_pp4 := 0]
atlas_bg <- atlas_bg[!gene %in% gene_tbl$gene]

classify_gene <- function(g) {
  fcase(
    g %in% drug_target_genes,      "Drug target",
    g %in% masld_associated_genes, "MASLD-associated",
    default                       = "Novel"
  )
}

gene_tbl[, fig5_group := classify_gene(gene)]

# Distribute Novel labels across the PP4 range so the panel isn't dominated
# by ~30 novels stacked at PP4>=0.9. Cap per-band to match fig3g density.
gene_tbl[, pp4_band := fcase(
  coloc_pp4 >= 0.9, "high",
  coloc_pp4 >= 0.7, "mid",
  coloc_pp4 >= 0.5, "low",
  default          = "below")]
NOVEL_CAP_PER_BAND <- c(high = 6, mid = 5, low = 5)

novel_pool <- gene_tbl[fig5_group == "Novel"]
novel_pool[, abs_lfc := abs(bulk_logFC)]
setorder(novel_pool, -coloc_pp4, -abs_lfc)
novel_pool[, rank_in_band := seq_len(.N), by = pp4_band]
novel_keep <- novel_pool[
  pp4_band %in% names(NOVEL_CAP_PER_BAND) &
  rank_in_band <= NOVEL_CAP_PER_BAND[pp4_band], gene]

gene_tbl <- gene_tbl[fig5_group != "Novel" | gene %in% novel_keep]
gene_tbl[, c("pp4_band") := NULL]

gene_tbl[, fig5_group := factor(fig5_group, levels = c(
  "Drug target", "Novel", "MASLD-associated"))]

tier_colors_3e <- c(
  "Drug target"      = "#880E4F",
  "Novel"            = "#00695C",
  "MASLD-associated" = "#7B1FA2"
)

# All anchors get labeled (filter is the entry criterion, not the label gate).
gene_tbl[, do_label := TRUE]
gene_tbl[, is_sig := factor(ifelse(!is.na(bulk_padj) & bulk_padj < 0.05,
                                    "sig", "ns"), levels = c("sig", "ns"))]

p3e <- ggplot(gene_tbl, aes(x = bulk_logFC, y = coloc_pp4, color = fig5_group)) +
  geom_point(data = atlas_bg, inherit.aes = FALSE,
             aes(x = bulk_logFC, y = coloc_pp4),
             color = "gray80", size = 0.4, alpha = 0.4, shape = 16) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
  geom_point(aes(shape = is_sig), size = 2.2, alpha = 0.9) +
  ggrepel::geom_text_repel(data = gene_tbl[do_label == TRUE],
                   aes(label = gene),
                   size = 2.1, max.overlaps = Inf,
                   force = 3, force_pull = 1,
                   box.padding = 0.3, point.padding = 0.15,
                   segment.size = 0.2, segment.color = "gray55",
                   min.segment.length = 0, fontface = "italic",
                   bg.color = "white", bg.r = 0.12,
                   show.legend = FALSE,
                   max.iter = 8000, seed = 1) +
  scale_color_manual(values = tier_colors_3e, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(limits = c(-0.02, 1.05),
                     breaks = c(0, 0.5, 0.9, 1.0),
                     labels = c("0", "0.5", "0.9", "1")) +
  labs(x = expression("Transcript log"[2]*"FC (MASLD vs control)"),
       y = "Best SuSiE colocalization PP4",
       title = "COLOC vs. RNA-seq DEG") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.25, "cm")) +
  scale_shape_manual(values = c("sig" = 16, "ns" = 1), name = "DEG Significance",
                     labels = c("sig" = "FDR < 0.05", "ns" = "n.s.")) +
  guides(color = guide_legend(override.aes = list(size = 2.5), nrow = 1),
         shape = guide_legend(nrow = 1))

# Finemapping PIP version
novel_genes_pip <- c("ALDH2", "AMN", "ANKRD9", "ATF6B", "BAK1", "BRD2",
                     "BTN2A2", "BTN3A3", "C9orf43", "CLCNKA")
novel_tiers_pip <- data.table(gene = novel_genes_pip, drug = NA_character_,
                              status = "Novel",
                              tier   = "Novel")
gene_meta_pip <- rbind(clinical_tiers, novel_tiers_pip)
gene_meta_pip[, tier := factor(tier,
  levels = c("FDA approved", "Phase 3", "Phase 2/earlier", "Novel"))]

gene_tbl_pip <- merge(gene_meta_pip, atlas_lite, by = "gene", all.x = TRUE)
gene_tbl_pip <- merge(gene_tbl_pip, bayev_lite, by = "gene", all.x = TRUE)
gene_tbl_pip[is.na(susiex_max_pip), susiex_max_pip := 0]
gene_tbl_pip[is.na(bulk_logFC), bulk_logFC := 0]

gene_tbl_pip[, fig5_group := classify_gene(gene)]
gene_tbl_pip[, fig5_group := factor(fig5_group, levels = c(
  "Drug target", "Novel", "MASLD-associated"))]

gene_tbl_pip[, do_label_pip := tier != "Novel" | susiex_max_pip >= 0.5]
gene_tbl_pip[, is_sig := factor(ifelse(!is.na(bulk_padj) & bulk_padj < 0.05,
                                        "sig", "ns"), levels = c("sig", "ns"))]

p3f <- ggplot(gene_tbl_pip, aes(x = bulk_logFC, y = susiex_max_pip,
                                color = fig5_group)) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "gray55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
  geom_point(aes(shape = is_sig), size = 2.2, alpha = 0.9) +
  geom_label_repel(data = gene_tbl_pip[do_label_pip == TRUE],
                   aes(label = gene),
                   size = 2.0, max.overlaps = 30,
                   label.padding = 0.1, segment.size = 0.2,
                   min.segment.length = 0, fontface = "italic",
                   fill = alpha("white", 0.85), show.legend = FALSE) +
  scale_color_manual(values = tier_colors_3e, name = NULL, drop = FALSE) +
  scale_y_continuous(limits = c(-0.02, 1.05),
                     breaks = c(0, 0.5, 0.9, 1.0),
                     labels = c("0", "0.5", "0.9", "1")) +
  annotate("text", x = Inf, y = 0.9, label = "PIP = 0.9 (strong)",
           hjust = 1.05, vjust = -0.3, size = 1.9, color = "gray30") +
  annotate("text", x = Inf, y = 0.5, label = "PIP = 0.5 (canonical)",
           hjust = 1.05, vjust = -0.3, size = 1.9, color = "gray45") +
  labs(x = expression("Transcript log"[2]*"FC (MASLD vs control)"),
       y = "Best SuSiE-X finemapping PIP (across 23 GWAS)",
       title = "Finemapping vs. RNA-seq") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm")) +
  scale_shape_manual(values = c("sig" = 16, "ns" = 1), name = "Significance",
                     labels = c("sig" = "FDR < 0.05", "ns" = "n.s.")) +
  guides(color = guide_legend(override.aes = list(size = 2.5), nrow = 2),
         shape = guide_legend(nrow = 1))

# Persist panel data for reviewers / caption text
fwrite(gene_tbl[, .(gene, fig5_group, bulk_logFC, bulk_padj,
                    coloc_pp4)],
       file.path(BASE, "RNA-seq/results/drug_repurposing/fig3f_scatter_data.csv"))

save_panel(p3e, "fig3e.pdf", width = fig_half_width, height = 4.0)

# Demote drug-target finemapping (PIP) scatter to figS_therapeutics.
THERA_PANELS <- file.path(FIGS_THERA_DIR, "panels")
dir.create(THERA_PANELS, showWarnings = FALSE, recursive = TRUE)
save_fig(p3f, file.path(THERA_PANELS, "fig3_drug_finemapping_demoted.pdf"),
         width = fig_half_width, height = 4.0)
cat("[fig3] Demoted drug-target finemapping scatter -> figS_therapeutics/panels/fig3_drug_finemapping_demoted.pdf\n")

# NOTE: fig3c.pdf (RORA / GGT chr15 locus LD-zoom) is generated externally
# by scripts/figures/run_fig3c_rora_locus_zoom.sh (figS09_locus_zoom.R writes
# directly via its 4th CLI arg). Do NOT delete fig3c.pdf.
for (stale in c("fig3g.pdf", "fig3h.pdf")) {
  p <- file.path(PANEL_DIR, stale)
  if (file.exists(p)) {
    file.remove(p); cat("[fig3] Removed stale", stale, "\n")
  }
}
# Old fig3b.pdf (RORA at the wrong slot) is retired by the cascade — remove.
old_b <- file.path(PANEL_DIR, "fig3b.pdf")
if (file.exists(old_b)) {
  bak <- file.path(PANEL_DIR, ".fig3b_pre_cascade.pdf")
  file.rename(old_b, bak)
  cat("[fig3] Old fig3b.pdf (RORA at wrong slot) renamed ->", bak, "\n")
}

# ===========================================================================
# Composite figure (6 displayed panels)
#   Row 1: 3a twin Manhattan (full width)
#   Row 2: 3b per-ancestry COLOC bars | 3c RORA locus LD-zoom (external PDF)
#   Row 3: 3d GWAS-ATAC PIP vs disruption | 3e high-PIP TF heatmap
#   Row 4: 3f drug-target logFC × COLOC PP4 (half-width, full row)
# Cascade (2026-04-29):
#   - old 3b RETIRED (figS09)             -> new 3b: per-ancestry COLOC bars
#   - old 3c GGT1 locus zoom              -> figS09
#   - external RORA panel                 -> moved fig3b -> fig3c
#   - old 3c GWAS-ATAC scatter            -> new 3d
#   - old 3d TF heatmap                   -> new 3e
#   - old 3e drug COLOC scatter           -> new 3f
#   - old 3f drug PIP scatter             -> demoted to figS_therapeutics
# ===========================================================================
cat("[fig3] Assembling composite (6 panels: a, b, c, d, e, f) ...\n")

p3a_wrapped <- wrap_elements(full = p3a)

# Inline the external RORA PDF as a raster so the composite includes it.
rora_pdf <- file.path(PANEL_DIR, "fig3c.pdf")
# magick::image_read_pdf needs the pdftools/poppler backend; if it is unavailable
# (rnaseq env lacks pdftools) fall through to the placeholder rather than abort,
# since the individual panels are the deliverable and fig3c is rasterised externally.
rora_img <- if (file.exists(rora_pdf) && requireNamespace("magick", quietly = TRUE)) {
  tryCatch(magick::image_read_pdf(rora_pdf, density = 300), error = function(e) NULL)
} else NULL
if (!is.null(rora_img)) {
  rora_grob <- grid::rasterGrob(rora_img, interpolate = TRUE)
  rora_panel <- wrap_elements(full = rora_grob)
} else {
  cat("[fig3] WARNING: fig3c.pdf (RORA) not yet generated; composite will use a placeholder.\n")
  rora_panel <- wrap_elements(full = grid::textGrob(
    "fig3c.pdf (RORA) pending\nrun run_fig3c_rora_locus_zoom.sh",
    gp = grid::gpar(fontsize = 10, col = "grey50")))
}

composite <- p3a_wrapped /
             (p3b | rora_panel) /
             (p3c | p3d) /
             p3e +
  plot_layout(heights = c(1.35, 1.0, 1.0, 1.2)) +
  plot_annotation(
    title = "Figure 3 | Multi-ancestry regulatory architecture of MASLD (23-GWAS portfolio: 14 EUR + 3 EAS + 3 AFR + 3 SAS)",
    tag_levels = list(c("a", "b", "c", "d", "e", "f"))
  ) &
  theme(plot.tag = element_text(size = 9, face = "bold"))

out_composite <- file.path(FIG3_DIR, "fig3_regulatory_architecture.pdf")
save_fig(composite, out_composite,
         width = fig_full_width, height = 13.5)
cat("[fig3] Composite saved:", out_composite, "\n")

cat("[fig3] DONE (6 displayed panels a,b,c,d,e,f; drug-PIP scatter demoted; GGT1 locus -> figS09).\n")
