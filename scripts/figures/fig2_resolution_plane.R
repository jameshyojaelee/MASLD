#!/usr/bin/env Rscript
# fig2_resolution_plane.R  (2026-06-12)
# Fig 2 HERO panel — "the resolution plane".
#
# A single 2-D mechanism plane for every fine-mapped MASLD effector gene:
#   x = within-ancestry causal-variant fine-mapping PIP  (how well-localized is the variant?)
#   y = SuSiE-coloc PP.H4                                 (does the signal act through expression?)
# colored by lead-variant class (coding vs non-coding), so the manuscript thesis
# (para 9) is visible at a glance:
#   * expression-mediated class  -> top-right  (high PIP, high PP.H4)  = RESOLVED by coloc
#   * coding / inactivating class-> high PIP, LOW PP.H4 = colocalization-BLIND (e.g. PNPLA3)
#
# x source = max within-ancestry SuSiE/recommended PIP among the gene's coloc-lead
#            variants (combined_finemapping.csv). NOT the contested cross-ancestry
#            SuSiEx max_pip (the "PNPLA3 0.987" artifact) and NOT the coloc-test PIP.
# y source = max PP.H4.susie across the 23-GWAS portfolio (susie_coloc_all_gwas.csv).
#
# Output: figures/main/fig2_genetics/panels/resolution_plane.pdf
#         figures/main/fig2_genetics/panels/resolution_plane_source.csv
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")          # FIG3_DIR == figures/main/fig2_genetics
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Per-gene colocalization (y) + coloc-lead variants
# ---------------------------------------------------------------------------
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
setnames(coloc, c("PP.H4.susie", "PP.H4.abf"), c("pp4_susie", "pp4_abf"), skip_absent = TRUE)
coloc[, pp4_susie := suppressWarnings(as.numeric(pp4_susie))]
coloc[, pp4_abf   := suppressWarnings(as.numeric(pp4_abf))]
coloc[, top_snp_PP := suppressWarnings(as.numeric(top_snp_PP))]

gene_y <- coloc[, .(
  pp4_susie = max(pp4_susie, na.rm = TRUE),
  pp4_abf   = max(pp4_abf,   na.rm = TRUE)
), by = gene]
gene_y[!is.finite(pp4_susie), pp4_susie := NA_real_]
gene_y[!is.finite(pp4_abf),   pp4_abf   := NA_real_]

# all coloc-lead variants a gene presents across the portfolio (chr:pos, hg19)
gene_leads <- coloc[!is.na(top_snp) & top_snp != "",
                    .(top_snp = unique(top_snp)), by = gene]

# ---------------------------------------------------------------------------
# 2. Within-ancestry causal PIP (x) from combined finemapping
# ---------------------------------------------------------------------------
fm <- fread(file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv"),
            select = c("chromosome", "position", "susie_pip", "recommended_pip"))
fm[, vkey := paste0(chromosome, ":", position)]
fm[, pip := pmax(suppressWarnings(as.numeric(recommended_pip)),
                 suppressWarnings(as.numeric(susie_pip)), na.rm = TRUE)]
fm_v <- fm[is.finite(pip), .(pip = max(pip)), by = vkey]

gene_leads <- merge(gene_leads, fm_v, by.x = "top_snp", by.y = "vkey", all.x = TRUE)
gene_x <- gene_leads[, .(fm_pip = if (all(is.na(pip))) NA_real_ else max(pip, na.rm = TRUE)),
                     by = gene]

# Fallback PIP sources so every coloc gene is placed (genes whose best coloc came
# from a GWAS not SuSiE-fine-mapped lack fm_pip): cross-ancestry SuSiEx max PIP,
# then the coloc-lead variant PIP. x = best available causal-localization confidence.
sx <- fread(file.path(BASE, "GWAS/finemapping/results/susiex/susiex_gene_summary.csv"),
            select = c("GeneSymbol", "susiex_max_pip"))
setnames(sx, c("GeneSymbol", "susiex_max_pip"), c("gene", "susiex_pip"))
sx[, susiex_pip := suppressWarnings(as.numeric(susiex_pip))]
gene_tp <- coloc[is.finite(top_snp_PP), .(coloc_lead_pip = max(top_snp_PP)), by = gene]

gene_x <- merge(gene_x, sx, by = "gene", all = TRUE)
gene_x <- merge(gene_x, gene_tp, by = "gene", all = TRUE)
gene_x[, causal_pip := pmax(fm_pip, susiex_pip, coloc_lead_pip, na.rm = TRUE)]
gene_x[!is.finite(causal_pip), causal_pip := NA_real_]
gene_x[, x_source := fifelse(is.na(fm_pip), fifelse(!is.na(susiex_pip), "susiex", "coloc_lead"),
                             "within_ancestry_susie")]

# ---------------------------------------------------------------------------
# 3. Lead-variant class (coding / non-coding) per gene
# ---------------------------------------------------------------------------
ann <- fread(file.path(BASE,
  "RNA-seq/results/coloc_variant_classes/coloc_variant_annotation.csv"))
ann[, pp4_best := suppressWarnings(as.numeric(pp4_best))]
ann <- ann[!is.na(gene_symbol) & gene_symbol != ""]
gene_cls <- ann[order(-pp4_best)][, .(coarse_class = coarse_class[1]), by = gene_symbol]
setnames(gene_cls, "gene_symbol", "gene")

# Curated label roster (panel_F7_canonical + the MASLD coding-locus literature).
# Label only the genes discussed in the Fig 2 text (keeps the panel uncluttered).
coding_genes  <- c("PNPLA3","HSD17B13","MBOAT7","MTARC1")
straddle_genes<- c("TM6SF2","GCKR")
expr_genes    <- c("RORA","GGT1","EPHA2")
roster <- c(coding_genes, straddle_genes, expr_genes)

# ---------------------------------------------------------------------------
# 4. Assemble plot table
# ---------------------------------------------------------------------------
dt <- merge(gene_y, gene_x, by = "gene", all = TRUE)
dt <- merge(dt, gene_cls, by = "gene", all.x = TRUE)

# class = LEAD-VARIANT class (not outcome): coding vs non-coding lead. Curated
# roster wins; else coarse_class of the gene's top coloc variant; else non-coding
# (96% of coloc-led genes are non-coding). The resolved-vs-blind OUTCOME is shown
# by the quadrant shading, not the colour.
dt[, class := fifelse(gene %in% coding_genes, "Coding",
              fifelse(gene %in% straddle_genes, "Coding + regulatory",
              fifelse(gene %in% expr_genes, "Regulatory (non-coding)",
              fifelse(coarse_class == "coding", "Coding",
                      "Regulatory (non-coding)"))))]
dt[is.na(class), class := "Regulatory (non-coding)"]
dt[, class := factor(class, levels = c("Regulatory (non-coding)","Coding + regulatory","Coding"))]

# y floor at 0 for display; genes with no coloc test but in roster -> pp4 = 0
dt[is.na(pp4_susie), pp4_susie := 0]
dt[, pp4_plot := pmax(pp4_susie, 0)]

# colocalization OUTCOME (the thesis axis): COLOUR now encodes whether the signal
# is resolved by colocalization, NOT the lead-variant class (which moves to SHAPE).
# This makes the para-9 partition visible in colour space (blue resolved cloud
# top, magenta coloc-blind band bottom) instead of in near-invisible shading.
dt[, outcome := fifelse(pp4_plot > 0.55, "Expression-resolved",
                fifelse(pp4_plot < 0.45, "Colocalization-blind", "Straddles threshold"))]
dt[, outcome := factor(outcome, levels = c("Expression-resolved", "Straddles threshold",
                                           "Colocalization-blind"))]

# Every gene with an x (causal_pip) OR in the roster is a candidate point.
have_xy <- dt[!is.na(causal_pip) | gene %in% roster]
have_xy <- have_xy[!is.na(causal_pip)]               # need an x to place a point

# ---- DECLUTTER: do NOT plot all ~18,975 genes (a saturated black blob at the
# origin). Keep the labeled roster ALWAYS, plus only background genes that carry
# real signal. The combined_finemapping universe is already locus-restricted, so
# nearly every gene has a high PIP (median ~0.43, 1,483 sit at exactly 1.0) and a
# loose pip gate keeps ~18k. Two complementary signal gates preserve the
# informative L-shape while removing the dense mid-PIP fuzz:
#   * pp4_susie  > 0.05  -> the colocalization-bearing genes (the horizontal
#                           low-PP.H4 band that rises into the resolved cloud).
#   * causal_pip >= 0.99 -> a representative high-confidence vertical band at
#                           PIP~1 (the coding / coloc-blind structure).
PP4_GATE <- 0.05
PIP_GATE <- 0.99
have_xy[, is_signal := (pp4_susie > PP4_GATE) | (causal_pip >= PIP_GATE) | gene %in% roster]
plot_dt <- have_xy[is_signal == TRUE]

lab_dt <- plot_dt[gene %in% roster]
n_bg    <- sum(!plot_dt$gene %in% roster)
message(sprintf("[resolution_plane] %d genes plotted (%d background + %d roster); gates pp4>%.2f OR pip>=%.2f",
                nrow(plot_dt), n_bg, nrow(lab_dt), PP4_GATE, PIP_GATE))
message(sprintf("  (declutter: dropped %d/%d low-signal genes from the origin fuzz)",
                nrow(have_xy) - nrow(plot_dt), nrow(have_xy)))
missing_roster <- setdiff(roster, lab_dt$gene)
if (length(missing_roster)) message("  roster genes without a causal-PIP hit (dropped): ",
                                     paste(missing_roster, collapse = ", "))

# ---------------------------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------------------------
# COLOUR = colocalization outcome (the thesis); SHAPE = lead-variant class.
out_resolved <- "#1565C0"   # blue    — expression-resolved (PP.H4 > 0.55)
out_straddle <- "#7B1FA2"   # violet  — straddles the 0.5 threshold
out_blind    <- "#C9265E"   # magenta — colocalization-blind (PP.H4 < 0.45)
outcome_cols <- c("Expression-resolved" = out_resolved,
                  "Straddles threshold" = out_straddle,
                  "Colocalization-blind" = out_blind)
class_shapes <- c("Regulatory (non-coding)" = 21, "Coding + regulatory" = 23, "Coding" = 24)

# Two display layers: a muted background cloud (light grey, small, raster) and a
# vivid labeled-roster layer (full class colour, larger, white stroke) so the
# canonical genes + the partition pop. Split the table.
bg_dt  <- plot_dt[!gene %in% roster]

# Directional repel TARGETS into OPEN regions. The roster piles up in two dense
# corners (expression genes at PIP~1 / PP.H4~0.85-1.0; coding genes at PIP~1 /
# PP.H4~0). Both sit at x~1.0, so nudging right pushes labels off-scale. We anchor
# each label to a fixed open coordinate (nudge_* = target - point) so leaders stay
# short and labels never overlap; repel then resolves residual collisions:
#   * expression / non-coding labels -> left of the high-PP.H4 cluster (open band
#     around x~0.66), fanned vertically over the upper plot.
#   * coding / inactivating labels   -> the empty lower-middle (x~0.62), fanned.
#   * coding+expression (straddle)   -> just left of their true position.
lab_dt[, c("tx", "ty") := .(causal_pip, pp4_plot)]   # default target = point

# expression band: spread targets across the open upper-left (x 0.55-0.78,
# y descending 0.98 -> 0.62) ordered by their true PP.H4 so leaders don't cross
expr_lab <- lab_dt[class == "Regulatory (non-coding)"][order(-pp4_plot)]
if (nrow(expr_lab)) {
  ty_seq <- seq(0.985, 0.62, length.out = nrow(expr_lab))
  tx_seq <- seq(0.50,  0.74, length.out = nrow(expr_lab))
  lab_dt[expr_lab, on = "gene", `:=`(tx = tx_seq, ty = ty_seq)]
}
# coding band: spread targets across the empty lower-middle (x 0.50-0.78,
# y rising 0.06 -> 0.32) ordered by true PIP
cod_lab <- lab_dt[class == "Coding"][order(causal_pip)]
if (nrow(cod_lab)) {
  ty_seq <- seq(0.06, 0.34, length.out = nrow(cod_lab))
  tx_seq <- seq(0.50, 0.80, length.out = nrow(cod_lab))
  lab_dt[cod_lab, on = "gene", `:=`(tx = tx_seq, ty = ty_seq)]
}
# straddle genes: nudge left and spread on y around the 0.5 threshold
str_lab <- lab_dt[class == "Coding + regulatory"][order(pp4_plot)]
if (nrow(str_lab)) {
  lab_dt[str_lab, on = "gene",
         `:=`(tx = pmax(0.40, causal_pip - 0.14),
              ty = seq(0.30, 0.62, length.out = nrow(str_lab)))]
}
lab_dt[, nudge_x := tx - causal_pip]
lab_dt[, nudge_y := ty - pp4_plot]

p <- ggplot(plot_dt, aes(x = causal_pip, y = pp4_plot)) +
  # gentle regime hint (colour now carries the regime, so keep shading faint)
  annotate("rect", xmin = 0, xmax = 1.03, ymin = 0.5, ymax = 1.03,
           fill = out_resolved, alpha = 0.04) +
  annotate("rect", xmin = 0, xmax = 1.03, ymin = -0.02, ymax = 0.5,
           fill = out_blind, alpha = 0.04) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  # --- background layer, COLOURED BY OUTCOME (rasterized) — shows the partition
  #     across the whole population, not just the labeled roster ---
  rasterize_layer(geom_point(data = bg_dt, aes(color = outcome), size = 0.5,
                             alpha = 0.32, shape = 16)) +
  # --- vivid labeled-roster layer: fill = outcome, SHAPE = lead-variant class ---
  geom_point(data = lab_dt, aes(fill = outcome, shape = class), size = 2.2,
             color = "white", stroke = 0.4) +
  ggrepel::geom_text_repel(
    data = lab_dt, aes(label = gene, color = outcome),
    size = 2.9, fontface = "italic", max.overlaps = Inf, show.legend = FALSE,
    nudge_x = lab_dt$nudge_x, nudge_y = lab_dt$nudge_y,
    box.padding = 0.45, point.padding = 0.2, min.segment.length = 0,
    force = 1.5, force_pull = 1.0, max.iter = 50000, max.time = 3,
    direction = "both", segment.size = 0.2, segment.color = "grey55",
    segment.alpha = 0.9, seed = 1, bg.color = "white", bg.r = 0.12) +
  scale_fill_manual(values = outcome_cols, name = "Colocalization outcome") +
  scale_color_manual(values = outcome_cols, guide = "none") +
  scale_shape_manual(values = class_shapes, name = "Lead variant") +
  scale_x_continuous(limits = c(0, 1.03), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(-0.02, 1.03), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Causal-variant fine-mapping PIP (within-ancestry SuSiE)",
       y = expression("Colocalization posterior  PP.H"[4]*"  (SuSiE-coloc)"),
       title = "SuSiE PIP vs. PPH4") +
  theme_masld(base_size = 9) +
  theme(legend.position = c(0.015, 0.55), legend.justification = c(0, 0.5),
        plot.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.30, "cm"),
        legend.spacing.y = unit(0.02, "cm"),
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA)) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 3.2, shape = 21)),
         shape = guide_legend(order = 2, override.aes = list(size = 2.6, fill = "grey55")))

save_fig(p, file.path(PANEL_DIR, "resolution_plane.pdf"),
         width = fig_col_width * 1.02, height = 4.05)

# ---------------------------------------------------------------------------
# 6. Freeze source CSV (so the cited numbers are reproducible). Record ALL genes
#    with an x/y (not only the plotted subset); `plotted` flags which passed the
#    declutter gates, `is_label` flags the curated roster.
# ---------------------------------------------------------------------------
plotted_genes <- plot_dt$gene
out <- have_xy[, .(gene, causal_pip = round(causal_pip, 4), x_source,
                   pp4_susie = round(pp4_susie, 4), pp4_abf = round(pp4_abf, 4),
                   class, outcome,
                   is_label = gene %in% roster,
                   plotted  = gene %in% plotted_genes)][order(-is_label, -plotted, -pp4_susie)]
fwrite(out, file.path(PANEL_DIR, "resolution_plane_source.csv"))
cat(sprintf("[resolution_plane] wrote panel + source CSV (%d genes total, %d plotted)\n",
            nrow(out), sum(out$plotted)))
