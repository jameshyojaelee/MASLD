#!/usr/bin/env Rscript
# fig2_coding_noncoding.R  (2026-06-12)  — Fig 2C
# Where does the colocalizing causal variant fall? Coding vs non-coding split.
#
# Answers the draft's Fig 2C placeholder "[How many noncoding vs coding variants?]".
# Single vertical stacked bar over the de-duplicated multi- or single-signal COLOC set
# (best PP.H4 > 0.5 per gene), segmented by the VEP CONSEQUENCE of that gene's strongest
# COLOC lead SNP (NOT a TSS-distance cut). Descriptive locus annotation (MAIN, Tier-1/2):
# only a small minority of colocalizing genes have a CODING strongest lead SNP; for the large
# majority the top lead variant falls in non-coding sequence. This describes where the
# strongest lead variant sits -- it is NOT a genome-wide inherited-risk architecture claim
# The noncoding-vs-coding split of total heritability is not estimated here;
# coding architecture is quantified separately from this gene-level lead annotation.
#
# 2026-06-25: switched from SuSiE-only to the de-duped SuSiE-OR-ABF union so the
#             denominator matches Fig2E (ancestry partition, also SuSiE/ABF). Per-gene
#             consequence = the fine_class of the gene's max-PP.H4 colocalization.
# 2026-07-04: portfolio rebuilt with MVP (23->50 GWAS).
# 2026-08-17: membership and lead variants are re-derived directly from the promoted
#             corrected long COLOC release; no July annotation table sets the denominator.
# Data: GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv plus
#       hg19 TxDb lead-variant classification performed in this script.
# Out : figures/main/fig2_genetics/panels/Fig2C_coding_noncoding_split.pdf (+ source CSV)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(VariantAnnotation); library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(GenomicRanges)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/fig2_consequence_palette.R"))
PANEL_DIR <- Sys.getenv("FIG2_CANDIDATE_DIR", file.path(FIG3_DIR, "panels"))
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
PDF_REL <- "main/fig2_genetics/panels/Fig2C_coding_noncoding_split.pdf"
SIZE_FILE <- file.path(BASE, "figures/layout_specs/figure2_panel_sizes.tsv")
EXPECTED_SIZE <- c(width_in = 1.10, height_in = 2.35)
COLOC_INPUT <- Sys.getenv(
  "FIG2_COLOC_INPUT",
  file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
VARIANT_CLASS_DIR <- Sys.getenv("FIG2_VARIANT_CLASS_DIR", "")

# Re-derive the named Tier-1/2 union directly from the promoted long COLOC
# release. The prior coloc_variant_annotation.csv froze the July 1,030-gene
# membership and must not determine the corrected-release denominator.
MAIN_STUDIES <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[
  placement == "main", study_name]
av <- fread(COLOC_INPUT,
  select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp"))
av <- av[gwas_name %in% MAIN_STUDIES & !is.na(gene) & gene != ""]
av[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
av[!is.finite(pp4_best), pp4_best := NA_real_]
setorder(av, gene, -pp4_best, gwas_name)
best <- av[, .SD[1], by = gene][is.finite(pp4_best) & pp4_best > 0.5]
setnames(best, "gene", "gene_symbol")
stopifnot(uniqueN(best$gene_symbol) == nrow(best))

# Classify the promoted lead variants in place with the same hg19 TxDb and
# precedence used by coloc_coding_noncoding/02_annotate_variants.R. This avoids
# silently joining corrected-release variants to a stale July annotation.
parts <- tstrsplit(best$top_snp, ":", fixed = TRUE)
variant_dt <- unique(data.table(
  variant_key = best$top_snp,
  chr = suppressWarnings(as.integer(parts[[1]])),
  pos = suppressWarnings(as.integer(parts[[2]]))))
if (variant_dt[is.na(chr) | is.na(pos), .N])
  stop("Malformed promoted COLOC top_snp coordinates in ", COLOC_INPUT)
class_rank <- c(coding = 1L, spliceSite = 2L, fiveUTR = 3L,
                threeUTR = 4L, promoter = 5L, intron = 6L, intergenic = 7L)
if (nzchar(VARIANT_CLASS_DIR)) {
  class_file <- file.path(VARIANT_CLASS_DIR, "variant_classification.csv")
  if (!file.exists(class_file)) stop("FIG2_VARIANT_CLASS_DIR lacks variant_classification.csv: ", class_file)
  variant_class <- fread(class_file, select = c("variant_key", "fine_class"))
  variant_class <- unique(variant_class, by = "variant_key")
  missing_variant <- setdiff(variant_dt$variant_key, variant_class$variant_key)
  if (length(missing_variant))
    stop(length(missing_variant), " promoted lead variants are absent from ", class_file)
} else {
  gr <- GRanges(seqnames = paste0("chr", variant_dt$chr),
                ranges = IRanges(start = variant_dt$pos, end = variant_dt$pos),
                variant_key = variant_dt$variant_key)
  seqlevelsStyle(gr) <- "UCSC"
  loc <- VariantAnnotation::locateVariants(
    gr, TxDb.Hsapiens.UCSC.hg19.knownGene,
    AllVariants(promoter = PromoterVariants(upstream = 2000, downstream = 200)))
  location_to_class <- function(x) {
    x <- as.character(x)
    fifelse(x %chin% names(class_rank), x, "intergenic")
  }
  variant_class <- data.table(
    variant_key = mcols(gr)$variant_key[loc$QUERYID],
    fine_class = location_to_class(loc$LOCATION))
  variant_class[, rank := class_rank[fine_class]]
  setorder(variant_class, variant_key, rank)
  variant_class <- variant_class[, .SD[1], by = variant_key][, rank := NULL]
  missing_variant <- setdiff(variant_dt$variant_key, variant_class$variant_key)
  if (length(missing_variant)) {
    variant_class <- rbind(
      variant_class,
      data.table(variant_key = missing_variant, fine_class = "intergenic"))
  }
}
best <- merge(best, variant_class, by.x = "top_snp", by.y = "variant_key",
              all.x = TRUE, sort = FALSE)
if (best[is.na(fine_class), .N]) stop("Not every promoted lead variant was classified")

# Merge the single spliceSite gene (1 gene) into coding to simplify the legend and bar
best[fine_class == "spliceSite", fine_class := "coding"]
best[fine_class %in% c("fiveUTR", "threeUTR"), fine_class := "UTR"]   # clump 5'/3' UTR into one group
f1 <- best[, .(n_genes = .N), by = fine_class]
f1[, `:=`(method = "Multi- or single-signal COLOC",
          coarse_class = fifelse(fine_class == "coding", "coding", "non-coding"))]

# consequence order — coding first so that coding lands at the very top of the vertical stack
lev <- c("coding","UTR","promoter","intron","intergenic")
labs <- c(intergenic="intergenic", intron="intron", promoter="promoter",
          UTR="UTR", coding="coding")
f1[, fine_class := factor(fine_class, levels = lev)]
f1[, method := factor(method, levels = "Multi- or single-signal COLOC")]
present_lev <- lev[lev %in% as.character(f1$fine_class)]   # legend shows only classes with data

fig2_consequence_cols <- FIG2_CONSEQUENCE_PALETTE[
  c("intergenic", "intron", "promoter", "UTR", "coding")]

tot  <- f1[, .(n = sum(n_genes)), by = method]
codg <- f1[fine_class %in% c("coding","spliceSite"), .(coding = sum(n_genes)), by = method]
ann  <- merge(tot, codg, by = "method")
ann[, pct := round(100 * coding / n, 1)]

ntot <- ann$n; ncod <- ann$coding; nnon <- ntot - ncod
legend_counts <- setNames(f1$n_genes, as.character(f1$fine_class))
legend_labs <- paste0(unname(labs[present_lev]), "\n", legend_counts[present_lev])
if (ntot > 1050) stop("Figure 2C y-axis contract requires review: union exceeds 1,050 genes")

p <- ggplot(f1, aes(x = method, y = n_genes, fill = fine_class)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.18) +
  # total above the vertical bar (black)
  geom_text(data = tot, aes(x = method, y = n, label = n),
            inherit.aes = FALSE, vjust = -0.55, hjust = 0.5,
            size = GEOM_TEXT_6PT, fontface = "plain") +
  scale_fill_manual(values = fig2_consequence_cols, labels = legend_labs,
                    name = NULL, breaks = present_lev) +
  scale_y_continuous(limits = c(0, 1100), breaks = c(0, 500, 1000),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = c(0.18, 0.18))) +
  labs(x = NULL, y = "Colocalizing genes") +
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(plot.margin = margin(10, 3, 3, 8, "pt"),
        legend.position = "right",
        legend.justification = "center",
        legend.key.width = unit(5, "pt"),
        legend.key.height = unit(5, "pt"),
        legend.text = element_text(size = 6, lineheight = 0.85),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "pt"),
        legend.spacing.y = unit(1, "pt"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

message(sprintf("[caption] Fig2C variant class: non-coding %d (%.0f%%) vs coding %d (%.0f%%); best PP.H4 > 0.5.",
                nnon, 100 - ann$pct, ncod, ann$pct))
source(file.path(BASE, "figures/layout_specs/regenerate_panels.R"))   # save_panel(): exact contract size + cairo_pdf
sizes <- read_sizes(SIZE_FILE)
contract <- sizes[sizes$pdf == PDF_REL, , drop = FALSE]
stopifnot(nrow(contract) == 1L)
if (!isTRUE(all.equal(unname(as.numeric(unlist(contract[1, c("width_in", "height_in")]))),
                            unname(EXPECTED_SIZE), tolerance = 1e-12))) {
  stop("Figure 2C size contract must be 1.10 x 2.35 in before rendering")
}
if (nzchar(Sys.getenv("FIG2_CANDIDATE_DIR"))) {
  contract$pdf <- basename(PDF_REL)
  save_panel(p, basename(PDF_REL), contract, PANEL_DIR)
} else {
  save_panel(p, PDF_REL, sizes, file.path(BASE, "figures"))
}

out <- merge(f1[, .(method, consequence = as.character(fine_class), coarse_class, n_genes)],
             ann[, .(method, method_total = n, coding_led = coding, coding_pct = pct)],
             by = "method")
fwrite(out, file.path(PANEL_DIR, "Fig2C_coding_noncoding_source.csv"))
cat(sprintf("[fig2C] wrote Fig2C_coding_noncoding_split.pdf  (multi- or single-signal COLOC coding-led = %d/%d = %.1f%%)\n",
            ann$coding, ann$n, ann$pct))
message(sprintf("CAPTION (Fig2C): Consequence class of each colocalizing gene's strongest COLOC lead SNP (best PP.H4 > 0.5 under multi- or single-signal COLOC; de-duplicated %d-gene named union). Only %d (%.1f%%) have a coding strongest lead SNP; for the remaining %.0f%% it falls in non-coding sequence (descriptive locus annotation, not a genome-wide noncoding-vs-coding heritability split).",
                ann$n, ann$coding, ann$pct, 100 - ann$pct))
