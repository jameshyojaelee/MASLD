#!/usr/bin/env Rscript
# composite_crossancestry.R — CANDIDATE composite (cross-ancestry portability)
# LEFT : co-expression matrix of cross-ancestry COLOC genes (pie glyphs), grouped by
#        ancestry-portability class. RIGHT: EUR PP.H4 vs EAS PP.H4 lollipop (SAME units;
#        the gap = ancestry portability / population-specific genetic support).
# Population genetics x transcriptomics. Output: FIG2_DIR/panels/composite_crossancestry.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/composite_helpers.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels"); DATA_DIR <- file.path(PANEL_DIR, "data")

ca <- fread(file.path(BASE, "RNA-seq/results/causal_inference/cross_ancestry/cross_ancestry_comparison.csv"))
# NEUTRAL selection: colocalized (PP.H4>0.5) in EITHER ancestry AND tested in >=2 ancestries. Do NOT
# pre-split into EUR-/EAS-specific classes — that split mechanically forced the EUR-vs-EAS lollipop to
# an artefactual r=-0.465. Now the lollipop shows the honest concordant/discordant mix.
ca <- ca[n_ancestry_tested >= 2 & (best_pp4_eur > 0.5 | best_pp4_eas > 0.5)]
ca[, pmax_pp4 := pmax(best_pp4_eur, best_pp4_eas, na.rm = TRUE)]
ca <- ca[order(-pmax_pp4)][1:min(.N, 24)]

co <- bulk_coexpr(ca$gene, BASE)
ca <- ca[gene %in% co$present]
cormat <- co$cormat[ca$gene, ca$gene]
# group by data-driven CO-EXPRESSION cluster (NOT ancestry class) -> group membership is independent of
# the PP.H4 lollipop (removes circularity); the matrix's blocks are genuine co-expression modules.
k <- 3   # ward.D2 for balanced modules (co-expression among genetically-selected genes is weak)
cl <- cutree(hclust(as.dist(1 - cormat), method = "ward.D2"), k = k)
gl <- paste0("Co-expr module ", 1:k)
group_of <- setNames(gl[cl], names(cl))
ord <- order_by_group_then_clust(cormat, group_of, gl)
cat(sprintf("[cross-ancestry] %d genes, %d co-expr modules; r(EUR PP.H4, EAS PP.H4)=%.2f (was -0.465)\n",
    nrow(ca), k, cor(ca$best_pp4_eur, ca$best_pp4_eas)))

mat <- pie_glyph_matrix(cormat, ord, group_of, group_levels = gl, italic_items = TRUE)
loll <- lollipop_panel(ca, "gene", ord, est1 = "best_pp4_eur", est2 = "best_pp4_eas",
                       est1_lab = "EUR PP.H4", est2_lab = "EAS PP.H4",
                       xlab = "Colocalization posterior (PP.H4)", group_of = group_of, ref0 = FALSE)
assemble_composite(mat, loll, widths = c(2.4, 1),
  out_pdf = file.path(PANEL_DIR, "composite_crossancestry.pdf"), width = 7.4, height = 4.4,
  caption = sprintf(paste0("CAPTION (cross-ancestry portability composite): LEFT = bulk co-expression ",
    "of %d colocalized genes (PP.H4>0.5 in EUR and/or EAS, tested in >=2 ancestries), grouped into ",
    "data-driven co-expression modules (NOT ancestry class, so grouping is independent of the PP.H4 ",
    "lollipop). RIGHT = European vs East-Asian colocalization posterior (PP.H4) per gene; the gap = ",
    "population-specific genetic support (EUR eQTL panel used for both, so hold non-EUR to a higher bar). ",
    "Portability is also shown at locus level in Fig 2 (tri-ancestry colocalization)."),
    nrow(ca)))
fwrite(ca[, .(gene, module = group_of[gene], best_pp4_eur = round(best_pp4_eur,3), best_pp4_eas = round(best_pp4_eas,3))],
       file.path(DATA_DIR, "composite_crossancestry.csv"))
