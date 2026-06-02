#!/usr/bin/env Rscript
# figS_thrb_locus.R
# THRB locus hero panel — stacked-track browser composite.
#
# KEY MESSAGE: THRB (resmetirom target; FDA-approved 2024 for MASH) sits in an
# active, hepatocyte-accessible regulatory locus that colocalizes specifically
# with UKBB GGT (SuSiE-COLOC PP.H4 = 1.000; lead variant chr3:24520283), with
# every other liver GWAS (NAFLD, ALT, AST, PDFF, cirrhosis, HCC) showing
# negligible PP.H4 — i.e. the GGT signal is biochemically specific, not generic
# liver disease.
#
# Tracks (top -> bottom):
#   1. THRB gene model (GENCODE v49)
#   2. ENCODE adult-liver H3K27ac (ENCSR230IMS / ENCFF905FLR, FC-over-control)
#   3. ENCODE adult-liver H3K27me3 (ENCSR266OIG / ENCFF653YIU, FC-over-control)
#   4. Hepatocyte scATAC peak track (Pipeline L8 cell-type peak set)
#   5. SuSiE-COLOC top-SNP lollipop (UKBB_GGT, chr3:24520283, PP.H4 = 1.000)
#   6. Multi-GWAS PP.H4 strip across 16 GWAS at the THRB locus
#
# Output: figures/supplementary/figS_therapeutics/figS_thrb_locus.pdf
# (`FIGS_THERA_DIR` already exists for THRB track + drug validation.)
#
# Run:  micromamba run -n rnaseq Rscript scripts/figures/figS_thrb_locus.R
#
# Notes / gaps documented in final-report message accompanying this script:
#   - No fine-mapped credible-set variant in our table sits in the THRB cis-window
#     for any GWAS we ran SuSiE on. THRB COLOC uses the coloc.abf fallback over the
#     full eGene cis-window (5,381 SNPs); the "top SNP" is the variant with the
#     largest single-SNP PP4 contribution from coloc.abf, NOT a SuSiE-PIP credible
#     set variant. We therefore plot it as a lollipop with its PP_top_snp from the
#     abf fallback, not a PIP.
#   - No motif-disruption hit FALLS AT the THRB locus — the 10 "THRB" rows in
#     motif_disruption_scores.csv refer to variants ELSEWHERE that disrupt the
#     THRB *motif* (i.e. THRB as a TF binding site). We omit this track to avoid
#     misrepresentation.

suppressPackageStartupMessages({
  library(Gviz)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(PROJ, "scripts/figures/publication_theme.R"))
source(file.path(PROJ, "scripts/figures/load_figure_data.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
CHR        <- "chr3"
GENE_START <- 24117153L   # GENCODE v49 THRB gene_start (HAVANA, minus strand)
GENE_END   <- 24495756L   # GENCODE v49 THRB gene_end
TOP_SNP    <- 24520283L   # UKBB_GGT top SNP (abf_fallback PP=1.000)
WIN_START  <- min(GENE_START, TOP_SNP) - 50000L   # ~24,067,153
WIN_END    <- max(GENE_END,   TOP_SNP) + 30000L   # ~24,550,283

cat(sprintf("Window: %s:%d-%d  (%.0f kb)\n", CHR, WIN_START, WIN_END,
            (WIN_END - WIN_START) / 1e3))

ROI <- GRanges(CHR, IRanges(WIN_START, WIN_END))

# ── 1. Gene model from GTF (GENCODE v49) ──────────────────────────────────────
GTF_PATH <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
cat("Loading GENCODE v49 GTF...\n")
gtf <- rtracklayer::import(GTF_PATH,
                           which = ROI,
                           feature.type = c("exon", "CDS", "UTR", "gene", "transcript"))

# Keep only THRB features
thrb <- gtf[!is.na(gtf$gene_name) & gtf$gene_name == "THRB"]

# Use the canonical MANE Select / Ensembl canonical transcript if present;
# fall back to longest CDS-bearing transcript.
canonical_tag <- function(x) {
  tags <- x$tag
  if (is.null(tags)) return(FALSE)
  any(grepl("MANE_Select", as.character(tags))) ||
    any(grepl("Ensembl_canonical", as.character(tags)))
}
transcripts <- thrb[thrb$type == "transcript"]
tx_canonical_id <- NULL
for (i in seq_along(transcripts)) {
  tx <- transcripts[i]
  tags <- as.character(unlist(tx$tag))
  if (any(c("MANE_Select", "Ensembl_canonical") %in% tags)) {
    tx_canonical_id <- tx$transcript_id
    break
  }
}
if (is.null(tx_canonical_id)) {
  # fall back: longest transcript
  tx_lengths <- width(transcripts)
  tx_canonical_id <- transcripts$transcript_id[which.max(tx_lengths)]
}
cat("THRB canonical transcript:", tx_canonical_id, "\n")

# Build GeneRegionTrack data.frame
exons <- thrb[thrb$type == "exon" & thrb$transcript_id == tx_canonical_id]
gene_df <- data.frame(
  chromosome = as.character(seqnames(exons)),
  start      = start(exons),
  end        = end(exons),
  width      = width(exons),
  strand     = as.character(strand(exons)),
  feature    = "protein_coding",
  gene       = "THRB",
  exon       = paste0("exon", seq_along(exons)),
  transcript = exons$transcript_id,
  symbol     = "THRB",
  stringsAsFactors = FALSE
)

gene_track <- GeneRegionTrack(
  range       = gene_df,
  genome      = "hg38",
  chromosome  = CHR,
  name        = "THRB",
  transcriptAnnotation = "symbol",
  collapseTranscripts  = "longest",
  shape       = "smallArrow",
  fill        = "#1A237E",
  col         = "#1A237E",
  background.title = "white",
  col.title   = "black",
  fontcolor.title = "black",
  cex.title   = 0.9
)

# ── 2-3. ENCODE bigWig tracks ─────────────────────────────────────────────────
ENCODE_DIR <- file.path(PROJ, "data/external/encode_liver")
H3K27AC_BW  <- file.path(ENCODE_DIR, "ENCFF905FLR_H3K27ac_ENCSR230IMS.bigWig")
H3K27ME3_BW <- file.path(ENCODE_DIR, "ENCFF653YIU_H3K27me3_ENCSR266OIG.bigWig")
stopifnot(file.exists(H3K27AC_BW), file.exists(H3K27ME3_BW))

cat("Importing ENCODE bigWigs over ROI...\n")
ac_gr  <- rtracklayer::import(H3K27AC_BW,  which = ROI, format = "bigWig")
me3_gr <- rtracklayer::import(H3K27ME3_BW, which = ROI, format = "bigWig")
cat(sprintf("  H3K27ac: %d intervals, max=%.2f\n",  length(ac_gr),  max(score(ac_gr))))
cat(sprintf("  H3K27me3: %d intervals, max=%.2f\n", length(me3_gr), max(score(me3_gr))))

ac_track <- DataTrack(
  range  = ac_gr,
  type   = "histogram",
  name   = "H3K27ac\n(liver, FC)",
  col.histogram = "#C2185B",
  fill.histogram = "#C2185B",
  background.title = "white",
  col.title = "black",
  fontcolor.title = "black",
  col.axis = "black",
  cex.title = 0.8,
  cex.axis  = 0.6,
  ylim = c(0, ceiling(max(score(ac_gr))))
)

me3_track <- DataTrack(
  range  = me3_gr,
  type   = "histogram",
  name   = "H3K27me3\n(liver, FC)",
  col.histogram = "#1565C0",
  fill.histogram = "#1565C0",
  background.title = "white",
  col.title = "black",
  fontcolor.title = "black",
  col.axis = "black",
  cex.title = 0.8,
  cex.axis  = 0.6,
  ylim = c(0, max(2, ceiling(max(score(me3_gr)))))
)

# ── 4. Hepatocyte scATAC peaks ────────────────────────────────────────────────
ATAC_PEAKS_BED <- file.path(PROJ,
  "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2/Hepatocytes_peaks.bed")
atac <- fread(ATAC_PEAKS_BED, header = FALSE,
              col.names = c("chrom", "start", "end"))
atac_roi <- atac[chrom == CHR & start >= WIN_START & end <= WIN_END]
cat(sprintf("Hepatocyte ATAC peaks in window: %d\n", nrow(atac_roi)))

atac_track <- AnnotationTrack(
  start  = atac_roi$start,
  end    = atac_roi$end,
  chromosome = CHR,
  genome     = "hg38",
  name       = "scATAC\n(Hep)",
  fill       = "#00695C",
  col        = "#00695C",
  background.title = "white",
  col.title  = "black",
  fontcolor.title = "black",
  cex.title  = 0.8,
  stacking   = "dense"
)

# ── 5. SuSiE-COLOC variant lollipop ───────────────────────────────────────────
# Read per-GWAS THRB COLOC rows; build a lollipop track with the top SNP
# (from abf_fallback) for the UKBB_GGT-converged signal.
ALL_GWAS_CSV <- file.path(PROJ,
  "GWAS/finemapping/results/susie_coloc_polyfun/susie_coloc_all_gwas_polyfun.csv")
all_gwas <- fread(ALL_GWAS_CSV)
thrb_rows <- all_gwas[gene == "THRB"]
cat("THRB COLOC rows:", nrow(thrb_rows), "\n")
print(thrb_rows[, .(gwas_name, PP.H4.abf, top_snp, top_snp_PP)])

# The single dominant signal: UKBB_GGT, top_snp chr3:24520283, PP4=0.9999
top_gwas_row <- thrb_rows[gwas_name == "UKBB_GGT"]
top_snp_pos  <- as.integer(sub("^3:", "", top_gwas_row$top_snp))
top_snp_pp4  <- top_gwas_row$PP.H4.abf

# Build a track: one strong lollipop at TOP_SNP with height = PP.H4 for GGT
# plus weaker lollipops for other GWAS lead SNPs at this locus (height = PP.H4
# of that GWAS, mostly < 0.05, so visually flat).
lolli_dt <- thrb_rows[, .(
  pos      = as.integer(sub("^3:", "", top_snp)),
  pp4      = PP.H4.abf,
  gwas_lab = gwas_name
)]
# Drop the canonical UKBB_GGT row here; we'll plot it separately for emphasis
lolli_other <- lolli_dt[gwas_lab != "UKBB_GGT"]

# DataTrack for lollipop heights — points + drop-lines
# Use a DataTrack with type=c("p","h") to get vertical "lollipop" sticks
lolli_top_track <- DataTrack(
  data      = top_snp_pp4,
  start     = top_snp_pos,
  end       = top_snp_pos,
  chromosome = CHR,
  genome    = "hg38",
  name      = "GGT\nPP.H4",
  type      = c("p", "h"),
  ylim      = c(0, 1.05),
  col       = "#C2185B",
  pch       = 19,
  cex       = 1.6,
  lwd       = 2,
  background.title = "white",
  col.title  = "black",
  fontcolor.title = "black",
  col.axis   = "black",
  cex.title  = 0.8,
  cex.axis   = 0.6
)

# Other GWAS lead-SNPs — gray lollipops, height = their PP.H4 (mostly ~0.005-0.035)
lolli_other_track <- DataTrack(
  data      = lolli_other$pp4,
  start     = lolli_other$pos,
  end       = lolli_other$pos,
  chromosome = CHR,
  genome    = "hg38",
  name      = "Other GWAS\ntop-SNP PP.H4",
  type      = c("p", "h"),
  ylim      = c(0, 0.05),
  col       = "#9E9E9E",
  pch       = 19,
  cex       = 0.9,
  lwd       = 1,
  background.title = "white",
  col.title  = "black",
  fontcolor.title = "black",
  col.axis   = "black",
  cex.title  = 0.7,
  cex.axis   = 0.55
)

# ── 6. Multi-GWAS PP.H4 strip — annotate THRB locus signal specificity ───────
# We build a HeatmapTrack-style band by encoding each GWAS as a colored
# AnnotationTrack interval at TOP_SNP +/- 5kb, with color = PP.H4.
# Since Gviz lacks a true horizontal heatmap row track, we instead emit a
# DataTrack with type="heatmap" centered at TOP_SNP for each GWAS at a fixed
# y-coordinate; this renders cleanly.
# Approach: collapse to a single GRanges with one row per GWAS at top_snp+/-2kb.

heat_dt <- thrb_rows[, .(gwas_name, PP.H4.abf)][order(-PP.H4.abf)]
heat_dt[, color := vapply(PP.H4.abf, function(x) {
  # Diverging-from-zero: white -> magenta
  pal <- colorRampPalette(c("#FFFFFF", "#FCE4EC", "#F8BBD0", "#F48FB1",
                            "#EC407A", "#C2185B", "#880E4F"))(101)
  pal[round(x * 100) + 1]
}, FUN.VALUE = character(1))]
print(heat_dt)

# Build a vector of values + tracks — easier: use IdeogramTrack-style discrete
# row of colored ticks. We'll render with AnnotationTrack(stacking="dense")
# per-GWAS (16 rows). Simpler: combine into a single 1-row strip per GWAS.
# To avoid 16 stacked tracks (excessive), collapse to 3 logical groups:
#   GGT (specific), other UKBB enzymes (ALT,AST), all NAFLD/PDFF/HCC/cirrhosis
heat_dt[, group := fcase(
  gwas_name == "UKBB_GGT",                     "UKBB GGT",
  gwas_name %in% c("UKBB_ALT", "UKBB_AST"),    "Other UKBB liver enzymes",
  default                                      = "NAFLD / PDFF / cirrhosis / HCC"
)]
heat_dt[, group := factor(group,
  levels = c("UKBB GGT",
             "Other UKBB liver enzymes",
             "NAFLD / PDFF / cirrhosis / HCC"))]

# Group-mean PP.H4
group_summary <- heat_dt[, .(
  mean_pp4 = mean(PP.H4.abf),
  max_pp4  = max(PP.H4.abf),
  n        = .N
), by = group][order(group)]
print(group_summary)

# Build 3 single-row strip tracks (one per GWAS category) — each fills the
# entire window in a colour proportional to that category's max PP.H4. Track
# titles on the left identify the category; the magenta gradient encodes
# PP.H4 magnitude.
pp4_to_color <- function(pp4) {
  pal <- colorRampPalette(c("#FFFFFF", "#FCE4EC", "#F48FB1",
                            "#EC407A", "#C2185B", "#880E4F"))(101)
  pal[round(pp4 * 100) + 1]
}

make_strip <- function(label, pp4) {
  AnnotationTrack(
    start      = WIN_START,
    end        = WIN_END,
    chromosome = CHR,
    genome     = "hg38",
    name       = sprintf("%s (%.3f)", label, pp4),
    fill       = pp4_to_color(pp4),
    col        = "#BDBDBD",
    background.title = "white",
    col.title  = "black",
    fontcolor.title = "black",
    cex.title  = 0.55,
    stacking   = "dense"
  )
}

strip_ggt   <- make_strip("GGT",
                          group_summary[group == "UKBB GGT", max_pp4])
strip_other <- make_strip("ALT/AST",
                          group_summary[group == "Other UKBB liver enzymes", max_pp4])
strip_nafld <- make_strip("NAFLD/HCC/Cirr",
                          group_summary[group == "NAFLD / PDFF / cirrhosis / HCC", max_pp4])

# ── Genome axis + ideogram ────────────────────────────────────────────────────
gaxis <- GenomeAxisTrack(
  add53 = TRUE, add35 = TRUE,
  col = "black", fontcolor = "black", cex = 0.7
)

# Highlight track over the top SNP
ht <- HighlightTrack(
  trackList  = list(gene_track, ac_track, me3_track, atac_track,
                    lolli_top_track, lolli_other_track,
                    strip_ggt, strip_other, strip_nafld),
  start      = TOP_SNP - 1500,
  end        = TOP_SNP + 1500,
  chromosome = CHR,
  col        = "#C2185B",
  fill       = "#FCE4EC",
  inBackground = TRUE,
  alpha = 0.25
)

# ── Render ────────────────────────────────────────────────────────────────────
OUT_PDF <- file.path(FIGS_THERA_DIR, "figS_thrb_locus.pdf")
cat("Writing PDF:", OUT_PDF, "\n")

pdf(OUT_PDF, width = 7.5, height = 10, useDingbats = FALSE)
plotTracks(
  list(gaxis, ht),
  from        = WIN_START,
  to          = WIN_END,
  chromosome  = CHR,
  sizes       = c(0.06, 0.10, 0.17, 0.15, 0.07, 0.14, 0.11, 0.06, 0.06, 0.06),
  cex.main    = 0.85,
  main = sprintf(
    "THRB locus  |  chr3:%.2f-%.2f Mb  |  UKBB-GGT COLOC PP.H4 = %.3f  (top SNP chr3:%s)",
    WIN_START / 1e6, WIN_END / 1e6,
    top_snp_pp4,
    format(top_snp_pos, big.mark = ",")
  ),
  background.panel = "white"
)
dev.off()

cat("\nDone.\n")
cat("PDF:", OUT_PDF, "\n")
cat(sprintf("PDF size: %.1f KB\n", file.info(OUT_PDF)$size / 1024))
