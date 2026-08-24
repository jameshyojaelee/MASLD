#!/usr/bin/env Rscript
# KEY MESSAGE: Variant localization and multi-signal colocalization resolve
# complementary parts of the inherited liver-trait architecture.
# fig2_resolution_plane.R  (2026-06-12)
# Fig 2 HERO panel — "the resolution plane".
#
# A single 2-D mechanism plane for every fine-mapped MASLD effector gene:
#   x = within-ancestry causal-variant fine-mapping PIP  (how well-localized is the variant?)
#   y = SuSiE-coloc PP.H4
# genes colored by lead-variant class (coding vs non-coding), so the
# manuscript thesis (para 9) is visible at a glance:
#   * expression-mediated class  -> top-right  (high PIP, high PP.H4)  = RESOLVED by coloc
#   * coding / inactivating class-> high PIP, LOW PP.H4 = colocalization-BLIND (e.g. PNPLA3)
#
# x source = max within-ancestry SuSiE/recommended PIP among the gene's coloc-lead
#            variants (combined_finemapping.csv). NOT the contested cross-ancestry
#            SuSiEx max_pip (the "PNPLA3 0.987" artifact) and NOT the coloc-test PIP.
# y source = max PP.H4.susie across the 35 Tier-1/2 (liver-specific) GWAS (susie_coloc_all_gwas.csv).
#
# Output: figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc.pdf
#         figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- Sys.getenv("FIG2_CANDIDATE_DIR", file.path(FIG3_DIR, "panels"))
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
PDF_REL <- "main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc.pdf"
SIZE_FILE <- file.path(BASE, "figures/layout_specs/figure2_panel_sizes.tsv")
EXPECTED_SIZE <- c(width_in = 1.80, height_in = 2.35)

# ---------------------------------------------------------------------------
# 1. Per-gene colocalization (y) + coloc-lead variants
# ---------------------------------------------------------------------------
COLOC_INPUT <- Sys.getenv(
  "FIG2_COLOC_INPUT",
  file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
coloc <- fread(COLOC_INPUT)
# MAIN (Tier-1/2, liver-specific) restriction (2026-07-06): keep only placement=="main"
# strata (NAFLD/NASH/PDFF + ALT/AST/GGT); Tier-3/4 supp strata move to a supplementary
# full-portfolio figure. y (max PP.H4.susie) is thus computed over the MAIN strata only.
MAIN_STUDIES <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[
  placement == "main", study_name]
coloc <- coloc[gwas_name %in% MAIN_STUDIES]
setnames(coloc, c("PP.H4.susie", "PP.H4.abf"), c("pp4_susie", "pp4_abf"), skip_absent = TRUE)
coloc <- coloc[!is.na(gene) & gene != ""]
coloc[, pp4_susie := suppressWarnings(as.numeric(pp4_susie))]
coloc[, pp4_abf   := suppressWarnings(as.numeric(pp4_abf))]
coloc[, top_snp_PP := suppressWarnings(as.numeric(top_snp_PP))]
coloc <- coloc[!is.na(top_snp) & top_snp != ""]
coloc[, c("snp_chr", "snp_pos") := tstrsplit(top_snp, ":", fixed = TRUE, type.convert = TRUE)]

# ---------------------------------------------------------------------------
# 2. Within-ancestry causal PIP (x) from combined finemapping
# ---------------------------------------------------------------------------
fm_canonical <- fread(file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv"),
                      select = c("study", "chromosome", "position", "susie_pip", "recommended_pip"))
fm_canonical[, x_source := "canonical_finemapping"]

fm_targeted_path <- file.path(BASE, "GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv")
if (file.exists(fm_targeted_path)) {
  fm_targeted <- fread(fm_targeted_path,
                       select = c("study", "chromosome", "position", "susie_pip", "recommended_pip"))
  fm_targeted[, x_source := "coloc_targeted_finemapping"]
  fm <- rbindlist(list(fm_canonical, fm_targeted), use.names = TRUE)
} else {
  message("WARNING: combined_finemapping_coloc_targeted.csv not found; using canonical fine-mapping only")
  fm <- fm_canonical
}
fm[, causal_pip := pmax(suppressWarnings(as.numeric(recommended_pip)),
                        suppressWarnings(as.numeric(susie_pip)), na.rm = TRUE)]
fm[!is.finite(causal_pip), causal_pip := NA_real_]
fm <- fm[is.finite(causal_pip)]

# Match COLOC top SNPs to fine-mapping by GWAS + variant coordinate. The
# targeted table recovers COLOC-positive variants missed by genome-wide clumping.
fm[, source_rank := fifelse(x_source == "canonical_finemapping", 1L, 2L)]
fm <- fm[order(study, chromosome, position, -causal_pip, source_rank)]
fm_unique <- unique(fm, by = c("study", "chromosome", "position"))

coloc_match <- merge(
  coloc,
  fm_unique[, .(study, chromosome, position, causal_pip, x_source)],
  by.x = c("gwas_name", "snp_chr", "snp_pos"),
  by.y = c("study", "chromosome", "position"),
  all.x = FALSE
)

message(sprintf("[resolution_plane] COLOC top-SNP fine-mapping coverage: %d / %d rows",
                nrow(coloc_match), nrow(coloc)))

# y = method-specific max PP.H4 across the GWAS portfolio.
max_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}
label_at_max <- function(value, label) {
  if (all(is.na(value))) NA_character_ else label[which.max(fifelse(is.na(value), -Inf, value))]
}
gene_y <- coloc[, .(
  pp4_susie = max_or_na(pp4_susie),
  pp4_abf   = max_or_na(pp4_abf),
  susie_gwas = label_at_max(pp4_susie, gwas_name),
  abf_gwas   = label_at_max(pp4_abf,   gwas_name)
), by = gene]

# x = best matched fine-mapping PIP among a gene's COLOC lead variants. This is
# intentionally gene-level, preserving the mechanism-plane contrast where coding
# loci can be high-PIP but low-PP.H4.
gene_x <- coloc_match[order(gene, -causal_pip, x_source),
                      .SD[1], by = gene][,
                      .(gene, causal_pip, x_source, x_top_snp = top_snp)]

# ---------------------------------------------------------------------------
# 3. Lead-variant class (coding / non-coding) per gene
# ---------------------------------------------------------------------------
VARIANT_CLASS_DIR <- Sys.getenv(
  "FIG2_VARIANT_CLASS_DIR",
  file.path(BASE, "RNA-seq/results/coloc_variant_classes"))
ann <- fread(file.path(VARIANT_CLASS_DIR, "coloc_variant_annotation.csv"))
ann[, pp4_best := suppressWarnings(as.numeric(pp4_best))]
ann <- ann[!is.na(gene_symbol) & gene_symbol != ""]
gene_cls <- ann[order(-pp4_best)][, .(coarse_class = coarse_class[1]), by = gene_symbol]
setnames(gene_cls, "gene_symbol", "gene")

# Curated class overrides remain broader than the label roster so unlabeled
# literature anchors retain their correct lead-variant class.
coding_genes  <- c("PNPLA3","HSD17B13","MBOAT7","MTARC1")
straddle_genes<- c("TM6SF2","GCKR")
expr_genes    <- c("RORA","GGT1","EPHA2","HKDC1")
# Five labels is the maximum readable at the contracted 1.80-in width.
roster <- c("HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA")

# ---------------------------------------------------------------------------
# 4. Assemble plot table
# ---------------------------------------------------------------------------
dt <- merge(gene_y, gene_x, by = "gene", all.x = TRUE)
dt <- merge(dt, gene_cls, by = "gene", all.x = TRUE)

# class = LEAD-VARIANT class (not outcome): coding vs non-coding lead. Curated
# roster wins; else coarse_class of the gene's top coloc variant; else non-coding
# (96% of coloc-led genes are non-coding). The resolved-vs-blind OUTCOME is shown
# by the quadrant shading, not the colour.
dt[, class := fifelse(gene %in% coding_genes, "Coding",
              fifelse(gene %in% straddle_genes, "Both",
              fifelse(gene %in% expr_genes, "Non-coding",
              fifelse(coarse_class == "coding", "Coding",
                      "Non-coding"))))]
dt[is.na(class), class := "Non-coding"]
dt[, class := factor(class, levels = c("Non-coding","Coding","Both"))]

# One SuSiE-coloc point per gene. ABF is intentionally not plotted or used as a
# fallback here; genes with no SuSiE posterior are absent from this panel.
plot_long <- dt[!is.na(pp4_susie) & is.finite(causal_pip),
                .(gene, pp4_gwas = susie_gwas, x_top_snp, causal_pip,
                  x_source, class, pp4 = pmax(pp4_susie, 0))]
plot_long[, `:=`(method = "SuSiE",
                 x_plot = causal_pip,
                 y_plot = pmax(pp4, 0))]

# Keep the figure focused on signal-bearing points plus the curated roster.
PP4_GATE <- 0.05
PIP_GATE <- 0.99
plot_long[, is_signal := (pp4 > PP4_GATE) | (causal_pip >= PIP_GATE) | gene %in% roster]
plot_dt <- plot_long[is_signal == TRUE]
bg_dt   <- plot_dt[!gene %in% roster]
roster_dt <- plot_dt[gene %in% roster]

# Draw the dominant regulatory class first so rare coding points are not hidden
# underneath coincident blue points in the rasterized background layer.
bg_dt <- bg_dt[order(class)]

lab_dt <- plot_dt[gene %in% roster][, {
  idx <- which.max(fifelse(is.na(y_plot), -Inf, y_plot))
  .SD[idx]
}, by = gene]

lab_dt[, c("tx", "ty") := .(x_plot, y_plot)]

message(sprintf("[resolution_plane] %d SuSiE-coloc points plotted (%d background + %d roster anchors); gates pp4>%.2f OR pip>=%.2f | x = fine-mapping PIP at matched COLOC top SNP",
                nrow(plot_dt), nrow(plot_dt[!gene %in% roster]), nrow(lab_dt), PP4_GATE, PIP_GATE))
message(sprintf("  (declutter: dropped %d/%d low-signal points from the origin fuzz)",
                nrow(plot_long) - nrow(plot_dt), nrow(plot_long)))
missing_roster <- setdiff(roster, lab_dt$gene)
if (length(missing_roster)) message("  roster genes without a causal-PIP hit (dropped): ",
                                     paste(missing_roster, collapse = ", "))

# ---------------------------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------------------------
# Every gene carries its lead-variant class color; roster genes are larger so
# their labels remain easy to associate with their points.
class_cols <- c("Non-coding" = "#1565C0",
                "Both" = "#7B1FA2",
                "Coding" = "#C9265E")

# Directional repel TARGETS into OPEN regions. The roster piles up in two dense
# corners (expression genes at PIP~1 / PP.H4~0.85-1.0; coding genes at PIP~1 /
# PP.H4~0). Both sit at x~1.0, so nudging right pushes labels off-scale. We anchor
# each label to a fixed open coordinate (nudge_* = target - point) so leaders stay
# short and labels never overlap; repel then resolves residual collisions:
#   * expression / non-coding labels -> left of the high-PP.H4 cluster (open band
#     around x~0.66), fanned vertically over the upper plot.
#   * coding / inactivating labels   -> the empty lower-middle (x~0.62), fanned.
#   * coding+expression (straddle)   -> just left of their true position.
# expression band: spread targets across the open upper-left (x 0.55-0.78,
# y descending 0.98 -> 0.62) ordered by their true PP.H4 so leaders don't cross
expr_lab <- lab_dt[class == "Non-coding"][order(-y_plot)]
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
str_lab <- lab_dt[class == "Both"][order(y_plot)]
if (nrow(str_lab)) {
  lab_dt[str_lab, on = "gene",
         `:=`(tx = pmax(0.40, x_plot - 0.14),
              ty = seq(0.28, 0.45, length.out = nrow(str_lab)))]
}
lab_dt[, nudge_x := tx - x_plot]
lab_dt[, nudge_y := ty - y_plot]

p <- ggplot(plot_dt, aes(x = x_plot, y = y_plot)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             color = "grey45") +
  # background layer: unlabeled genes use the same class colors as the legend
  rasterize_layer(geom_point(data = bg_dt, aes(color = class), size = 0.45,
                             alpha = 0.8, shape = 16)) +
  # labeled roster: class colors, unified point shape, no outlines
  geom_point(data = roster_dt, aes(color = class), size = 1.6, shape = 16,
             alpha = 0.95, stroke = 0) +
  ggrepel::geom_text_repel(
    data = lab_dt, aes(label = gene),
    size = 6 / ggplot2::.pt, fontface = "plain", max.overlaps = Inf, show.legend = FALSE,
    nudge_x = lab_dt$nudge_x, nudge_y = lab_dt$nudge_y,
    box.padding = 0.25, point.padding = 0.15, min.segment.length = 0,
    force = 1.5, force_pull = 1.0, max.iter = 50000, max.time = 3,
    direction = "both", segment.size = 0.15, segment.color = "grey55",
    segment.alpha = 0.9, seed = 1, bg.color = "white", bg.r = 0.12) +
  scale_color_manual(values = class_cols, name = NULL,
                     breaks = c("Non-coding", "Coding", "Both")) +
  # Both axes keep the full 0-1 probability range (plus the same small
  # point/label clearance at each edge), but the plane is NOT forced square:
  # the panel fills the 1.80 x 2.35-in slot, so it renders portrait (~1:1.32).
  # Nothing here is read off a 45-degree diagonal - the content is the point
  # mass at x=0, x=1 and y=1 plus the PP.H4=0.5 threshold - so the vertical stretch costs no interpretation
  # and recovers ~0.5 in of otherwise stranded height below the x-axis title.
  scale_x_continuous(limits = c(-0.02, 1.03), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(-0.02, 1.03), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Fine-mapping PIP",
       y = "Multi-signal COLOC PP.H4") +
  coord_cartesian(clip = "off") +
  theme_masld(base_size = 6) +
  theme(legend.position = c(0.48, 0.64),
        legend.justification = "center",
        legend.direction = "horizontal",
        plot.title = element_blank(),
        plot.margin = margin(8, 1, 1, 3, "pt"),
        axis.ticks.length = unit(2, "pt"),
        axis.title = element_text(size = 6, face = "plain"),
        axis.title.y = element_text(margin = margin(r = 0.5, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 1, unit = "pt")),
        axis.text = element_text(size = 6, face = "plain", color = "black"),
        axis.text.y = element_text(margin = margin(r = 0.5, unit = "pt")),
        axis.text.x = element_text(margin = margin(t = 1, unit = "pt")),
        legend.key.size = unit(5, "pt"),
        legend.spacing.x = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 6, face = "plain"),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA)) +
  guides(color = guide_legend(order = 1,
                              override.aes = list(size = 1.8, alpha = 0.95),
                              nrow = 1, byrow = TRUE))

source(file.path(BASE, "figures/layout_specs/regenerate_panels.R"))
sizes <- read_sizes(SIZE_FILE)
contract <- sizes[sizes$pdf == PDF_REL, , drop = FALSE]
stopifnot(nrow(contract) == 1L)
if (!isTRUE(all.equal(unname(as.numeric(unlist(contract[1, c("width_in", "height_in")]))),
                            unname(EXPECTED_SIZE), tolerance = 1e-12))) {
  stop("Figure 2B size contract must be 1.80 x 2.35 in before rendering")
}
if (nzchar(Sys.getenv("FIG2_CANDIDATE_DIR"))) {
  contract$pdf <- basename(PDF_REL)
  save_panel(p, basename(PDF_REL), contract, PANEL_DIR)
} else {
  save_panel(p, PDF_REL, sizes, file.path(BASE, "figures"))
}

# ---------------------------------------------------------------------------
# 6. Freeze source CSV (so the cited numbers are reproducible). Record ALL
#    plotted gene-method points (not only the labeled roster).
# ---------------------------------------------------------------------------
out <- copy(plot_long)
out[, `:=`(
  causal_pip = round(causal_pip, 4),
  pp4 = round(pp4, 4),
  x_plot = round(x_plot, 4),
  y_plot = round(y_plot, 4),
  is_label = FALSE,
  plotted = is_signal
)]
out[lab_dt, on = .(gene, method), is_label := TRUE]
out <- out[, .(gene, method, pp4_gwas, x_top_snp, causal_pip, pp4, x_plot, y_plot, x_source,
               class, is_label, plotted)][order(-is_label, gene, method)]
fwrite(out, file.path(PANEL_DIR, "Fig2B_PIP_vs_SuSiE-coloc_source.csv"))
cat(sprintf("[resolution_plane] wrote panel + source CSV (%d gene-method rows, %d plotted)\n",
            nrow(out), sum(out$plotted)))
