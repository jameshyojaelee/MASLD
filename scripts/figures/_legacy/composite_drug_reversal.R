#!/usr/bin/env Rscript
# composite_drug_reversal.R — CANDIDATE composite (⑦ drug-target opportunity)
# LEFT : co-expression matrix of druggable disease targets (pie glyphs), grouped by
#        function. RIGHT: disease dysregulation vs druggability (number of DGIdb drugs,
#        log-scaled) per target; the gap = strongly dysregulated biology that is poorly
#        drugged (the pharmacological opportunity). Transcriptomics x pharmacology.
# (LINCS reversal / OpenTargets phase were binary/empty in the summary, so druggability
#  = DGIdb drug count is the continuous "captured" axis.)
# Output: FIG2_DIR/panels/composite_drug_reversal.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/composite_helpers.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels"); DATA_DIR <- file.path(PANEL_DIR, "data")

grp <- list(
  "Metabolic"    = c("FASN","SCD","LPL","FABP4","CA12","AKR1B10","SLC22A12","ACACA"),
  "Inflammatory" = c("CCL20","CXCL10","MMP9","OLR1","SPP1","LGALS3","IL1B","CCR2"),
  "Fibrotic"     = c("COL1A1","COL3A1","COL10A1","COMP","LOXL2","PDGFRB","MMP2","TIMP1"))
gl <- names(grp)
group_of <- setNames(rep(gl, lengths(grp)), unlist(grp, use.names = FALSE))

d <- fread(file.path(BASE, "RNA-seq/results/drug_repurposing/pharmacotranscriptomics_summary.csv"))
d <- d[symbol %in% names(group_of) & !is.na(bulk_logFC) & !is.na(dgidb_n_drugs)]
d <- unique(d, by = "symbol")

co <- bulk_coexpr(d$symbol, BASE)
d <- d[symbol %in% co$present]
group_of <- group_of[d$symbol]; gl <- gl[gl %in% group_of]
cormat <- co$cormat[d$symbol, d$symbol]
d[, dysregulation := abs(bulk_logFC) / max(abs(bulk_logFC))]
d[, druggability := log1p(dgidb_n_drugs) / max(log1p(dgidb_n_drugs))]
ord <- order_by_group_then_clust(cormat, group_of, gl)
cat(sprintf("[drug-opportunity] %d targets across %d functional groups\n", nrow(d), length(gl)))

mat <- pie_glyph_matrix(cormat, ord, group_of, group_levels = gl, italic_items = TRUE)
loll <- lollipop_panel(d, "symbol", ord, est1 = "dysregulation", est2 = "druggability",
                       est1_lab = "disease dysregulation", est2_lab = "DGIdb drug count (log)",
                       xlab = "Strength (0-1)", group_of = group_of, ref0 = FALSE)
assemble_composite(mat, loll, widths = c(2.4, 1),
  out_pdf = file.path(PANEL_DIR, "composite_drug_reversal.pdf"), width = 7.4, height = 4.4,
  caption = sprintf(paste0("CAPTION (drug-target opportunity composite): LEFT = bulk co-expression of ",
    "%d druggable disease targets (pie fill proportional to |Spearman|), grouped by function. RIGHT = ",
    "disease dysregulation (|bulk log2FC|, scaled) vs the number of DGIdb-annotated drugs (log, scaled). ",
    "CAVEAT: DGIdb drug count is an ANNOTATION/effort measure, NOT true druggability — it is biased toward ",
    "well-studied genes (e.g. kinases like PDGFRB), so a low count (e.g. TIMP1, collagens) reflects hard-to-",
    "drug biology as much as unmet opportunity. A continuous LINCS reversal score was unavailable for these ",
    "targets (only top-50 compounds are scored)."), nrow(d)))
fwrite(d[, .(symbol, group = group_of[symbol], bulk_logFC = round(bulk_logFC,3), dgidb_n_drugs)],
       file.path(DATA_DIR, "composite_drug_reversal.csv"))
