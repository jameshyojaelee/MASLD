#!/usr/bin/env Rscript
# composite_transcription_genetics.R — CANDIDATE composite (④ transcription vs genetics)
# LEFT : co-expression matrix of disease genes (pie glyphs), grouped by function.
# RIGHT: transcription evidence (|bulk logFC|, scaled 0-1) vs genetic evidence (COLOC PP.H4)
#        lollipop — the 3E disconnect as a total-vs-captured gap.
# Output: FIG2_DIR/panels/composite_transcription_genetics.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/composite_helpers.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels"); DATA_DIR <- file.path(PANEL_DIR, "data")

# DE-CIRCULARIZED: the high-COLOC genes (THRB/RORA/HKDC1/CYP3A4/GCKR/TM6SF2) are distributed
# into their FUNCTIONAL groups, NOT isolated in a "Colocalized" group (which made group membership
# predict the PP.H4 y-axis at R^2=0.69). Now group ⊥ PP.H4.
grp <- list(
  "Lipid & metabolism"    = c("FASN","SCD","DGAT2","ACACA","HKDC1","GCKR","TM6SF2","THRB"),
  "Inflammation & immune" = c("CRP","LCN2","S100A9","CXCL9","CCL20","ANXA1","RORA"),
  "ECM & fibrosis"        = c("COL1A1","COL3A1","LUM","SPARC","TIMP1","LGALS3"),
  "Detox & redox"         = c("GSTM1","AKR1B10","CYP2E1","SOD2","CYP3A4"))
gl <- names(grp)
group_of <- setNames(rep(gl, lengths(grp)), unlist(grp, use.names = FALSE))

at <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
            select = c("human_symbol", "bulk_logFC"))
at <- at[human_symbol %in% names(group_of) & is.finite(bulk_logFC)]
at <- unique(at, by = "human_symbol")
# canonical PP.H4 (gene_level_coloc.csv), matching Fig 3E — NOT the atlas coloc column
glc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
             select = c("gene", "coloc_best_pp4"))
glc <- glc[gene != "" & !is.na(gene)][, .(coloc_pp4 = max(coloc_best_pp4, na.rm = TRUE)), by = gene]
at <- merge(at, glc, by.x = "human_symbol", by.y = "gene", all.x = TRUE)
at[is.na(coloc_pp4), coloc_pp4 := 0]

co <- bulk_coexpr(at$human_symbol, BASE)
at <- at[human_symbol %in% co$present]
group_of <- group_of[at$human_symbol]
gl <- gl[gl %in% group_of]
cormat <- co$cormat[at$human_symbol, at$human_symbol]
# two evidence tracks on a shared 0-1 scale. Transcription uses a FIXED cap (|logFC|=2 -> 1),
# NOT max-normalization (a moving target that shifts if one outlier gene is added/removed).
at[, transcription := pmin(abs(bulk_logFC) / 2, 1)]
at[, genetics := pmin(coloc_pp4, 1)]
ord <- order_by_group_then_clust(cormat, group_of, gl)
cat(sprintf("[transcription-genetics] %d genes across %d groups\n", nrow(at), length(gl)))

mat <- pie_glyph_matrix(cormat, ord, group_of, group_levels = gl, italic_items = TRUE)
loll <- lollipop_panel(at, "human_symbol", ord, est1 = "transcription", est2 = "genetics",
                       est1_lab = "transcription (|log2FC|/2, capped)", est2_lab = "genetics (PP.H4)",
                       xlab = "Evidence strength (0-1)", group_of = group_of, ref0 = FALSE)
assemble_composite(mat, loll, widths = c(2.4, 1),
  out_pdf = file.path(PANEL_DIR, "composite_transcription_genetics.pdf"), width = 7.4, height = 4.4,
  caption = sprintf(paste0("CAPTION (transcription vs genetics composite): LEFT = bulk co-expression of ",
    "%d disease genes (pie fill proportional to |Spearman|), grouped by FUNCTION (high-COLOC genes ",
    "distributed into functional groups, so group membership does not encode PP.H4). RIGHT = transcription ",
    "evidence (|bulk log2FC|/2, fixed cap) vs genetic evidence (COLOC PP.H4) per gene. Within each ",
    "functional group, strong dysregulation and genetic support rarely coincide (the genetics<->expression ",
    "disconnect; complements the gene-level Fig 3E scatter)."), nrow(at)))
fwrite(at[, .(human_symbol, group = group_of[human_symbol], bulk_logFC = round(bulk_logFC,3),
        coloc_pp4 = round(coloc_pp4,3))], file.path(DATA_DIR, "composite_transcription_genetics.csv"))
