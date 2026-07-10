#!/usr/bin/env Rscript
# spotcheck_guides_plotgardener.R
# plotgardener coverage figure for the Cas13 guide library: renders EVERY targeting
# guide in the final library (data/guides/library.csv, non-targeting controls
# excluded), one compact track per gene, with the gene's selected guides overlaid.
#
# Each gene is drawn as a COLLAPSED (union) gene model — CDS / UTR / lncRNA exons
# colour-coded, introns with strand arrows — with its guides shown as red position
# markers carrying de-overlapped gNN labels, exon-region calls, TIGER / Cas13Design
# scores, and inter-guide spacing (mature-mRNA nt). A plotgardener genome-coordinate
# ruler anchors each track. Design mirrors the matplotlib spot-check figure
# (spotcheck2_gene_tracks.pdf) but at library scale with real annotation structure.
#
# Guide genomic coordinates are read from the coordinate table produced by
# build_all_guide_coords.py (which reuses spotcheck_100_guides.analyse()), so they use
# the same trusted sequence-search coordinates as the matplotlib QC figure.
#
# Inputs:
#   data/guide_qc_all_library.tsv                (gene, guide, genomic_pos, mrna_pos, region, scores)
#   data/guides/library.csv                      (gene_symbol_mouse -> gene_id_mouse map)
#   cellranger GRCm39-vM38 genes.gtf.gz          (TxDb source; same GTF as the QC step)
# Output:
#   figures/guide_qc/spotcheck2_gene_tracks_plotgardener.pdf  (paginated)
#
# Env:
#   SPOTCHECK_TSV    coordinate table (default data/guide_qc_all_library.tsv)
#   GENES_PER_PAGE   default 6
#   TEST_N           if set (integer), render only the first N genes (debugging)

suppressPackageStartupMessages({
  library(plotgardener)
  library(GenomicFeatures)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(rtracklayer)
  library(grid)
})

# Build a TxDb straight from the GTF via the low-level makeTxDb() constructor.
build_txdb_from_gtf <- function(gtf) {
  message("[txdb] importing GTF via rtracklayer (heavy, one-time) ...")
  gr <- rtracklayer::import(gtf, format = "gtf")
  tx <- gr[gr$type == "transcript"]
  transcripts <- data.frame(
    tx_id = seq_along(tx), tx_name = as.character(tx$transcript_id),
    tx_chrom = as.character(seqnames(tx)), tx_strand = as.character(strand(tx)),
    tx_start = start(tx), tx_end = end(tx), stringsAsFactors = FALSE)
  txmap <- setNames(transcripts$tx_id, transcripts$tx_name)
  ex <- gr[gr$type == "exon"]
  spl <- data.frame(
    tx_id = unname(txmap[as.character(ex$transcript_id)]),
    exon_rank = as.integer(ex$exon_number),
    exon_start = start(ex), exon_end = end(ex), stringsAsFactors = FALSE)
  cds <- gr[gr$type == "CDS"]
  if (length(cds)) {
    cds_df <- data.frame(
      tx_id = unname(txmap[as.character(cds$transcript_id)]),
      exon_rank = as.integer(cds$exon_number),
      cds_start = start(cds), cds_end = end(cds), stringsAsFactors = FALSE)
    spl <- merge(spl, cds_df, by = c("tx_id", "exon_rank"), all.x = TRUE)
  }
  spl <- spl[order(spl$tx_id, spl$exon_rank), ]
  genes <- data.frame(tx_id = transcripts$tx_id,
                      gene_id = as.character(tx$gene_id), stringsAsFactors = FALSE)
  message("[txdb] makeTxDb: ", nrow(transcripts), " transcripts, ", nrow(spl), " exons")
  GenomicFeatures::makeTxDb(transcripts = transcripts, splicings = spl, genes = genes)
}

ROOT   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CAS13  <- file.path(ROOT, "Cas13_Library_Design")
TSV    <- Sys.getenv("SPOTCHECK_TSV", file.path(CAS13, "data/guide_qc_all_library.tsv"))
LIBCSV <- file.path(CAS13, "data/guides/library.csv")
GTF    <- "/gpfs/commons/home/jameslee/reference_genome/refdata-cellranger-GRCm39-vM38/genes/genes.gtf.gz"
CACHE  <- file.path(CAS13, "data/guides/cache/txdb_vM38_cellranger.sqlite")
OUTPDF <- file.path(CAS13, "figures/guide_qc/spotcheck2_gene_tracks_plotgardener.pdf")

GUIDE_LEN      <- 23L
GENES_PER_PAGE <- as.integer(Sys.getenv("GENES_PER_PAGE", "6"))
TEST_N         <- suppressWarnings(as.integer(Sys.getenv("TEST_N", "")))

# palette (matches spotcheck2_gene_tracks.pdf)
CDS_C <- "#3B6FB6"; UTR_C <- "#BBBBBB"; LNC_C <- "#8E6FC7"
INTRON_C <- "#9AA0A6"; RED <- "#C0143C"; LEAD_C <- "#B9414E"; GREY <- "#777777"

# ---------------------------------------------------------------------------
# 1. TxDb + collapsed (union) gene models
# ---------------------------------------------------------------------------
if (file.exists(CACHE)) {
  message("[txdb] loading cached ", CACHE); txdb <- loadDb(CACHE)
} else {
  message("[txdb] building from ", GTF, " (one-time, slow)")
  txdb <- build_txdb_from_gtf(GTF); saveDb(txdb, CACHE)
}

message("[models] building union exon / CDS models per gene ...")
exbg  <- exonsBy(txdb, by = "gene")
cdsbg <- cdsBy(txdb, by = "gene")
n_iso <- lengths(transcriptsBy(txdb, by = "gene"))
strip <- function(x) sub("\\.\\d+$", "", x)
names(exbg)  <- strip(names(exbg))
names(cdsbg) <- strip(names(cdsbg))
names(n_iso) <- strip(names(n_iso))

asm <- assembly(Genome = "mm39_vM38", TxDb = txdb, OrgDb = org.Mm.eg.db,
                gene.id.column = "GENEID", display.column = "GENEID")

# ---------------------------------------------------------------------------
# 2. symbol -> gene_id_mouse map + guide coordinate table
# ---------------------------------------------------------------------------
lib <- read.csv(LIBCSV, stringsAsFactors = FALSE)
sym2gid_raw <- tapply(lib$gene_id_mouse, lib$gene_symbol_mouse, function(x) x[1])
sym2gid <- setNames(strip(as.character(sym2gid_raw)), names(sym2gid_raw))

qc <- read.delim(TSV, stringsAsFactors = FALSE)
qc <- qc[!is.na(suppressWarnings(as.numeric(qc$genomic_pos))), ]
qc$genomic_pos <- as.numeric(qc$genomic_pos)
qc$mrna_pos    <- suppressWarnings(as.numeric(qc$mrna_pos))
genes_order <- unique(qc$gene)
if (!is.na(TEST_N)) genes_order <- head(genes_order, TEST_N)
message("[guides] ", nrow(qc), " guides over ", length(unique(qc$gene)),
        " genes; rendering ", length(genes_order), " genes")

# ---------------------------------------------------------------------------
# 3. layout (inches). Text uses plotgardener top-origin y; custom shapes use
#    grid in a device-filling viewport (bottom-origin), converted via yt().
# ---------------------------------------------------------------------------
PAGE_W <- 9.0; MARGIN_L <- 0.95; BLOCK_W <- 7.35
TOP <- 0.45; BLOCK_H <- 1.95
PAGE_H <- TOP + GENES_PER_PAGE * BLOCK_H + 0.30
RX <- MARGIN_L + BLOCK_W                      # right edge of the track area

# within-block offsets from block top
OFF_TITLE <- 0.02; OFF_SCORE <- 0.22; OFF_REGION <- 0.35; OFF_LABEL <- 0.48
OFF_MARK  <- 0.60; OFF_MODEL <- 1.04; OFF_SPACE <- 1.34; OFF_RULER <- 1.62

yt <- function(y) unit(PAGE_H - y, "inches")   # top-origin inch -> grid native
xi <- function(x) unit(x, "inches")

# map a genomic position to a page x (inches) within [cs, ce]
x2page <- function(pos, cs, ce) MARGIN_L + (pmin(pmax(pos, cs), ce) - cs) / (ce - cs) * BLOCK_W

# 1-D label de-overlap: keep order, enforce a minimum gap, clamp to [lo, hi]
spread <- function(x, mingap, lo, hi) {
  if (length(x) <= 1) return(pmin(pmax(x, lo), hi))
  o <- order(x); xs <- x[o]
  for (i in 2:length(xs)) if (xs[i] < xs[i - 1] + mingap) xs[i] <- xs[i - 1] + mingap
  if (xs[length(xs)] > hi) {
    xs[length(xs)] <- hi
    for (i in (length(xs) - 1):1) if (xs[i] > xs[i + 1] - mingap) xs[i] <- xs[i + 1] - mingap
  }
  if (xs[1] < lo) {
    xs[1] <- lo
    for (i in 2:length(xs)) if (xs[i] < xs[i - 1] + mingap) xs[i] <- xs[i - 1] + mingap
  }
  out <- numeric(length(x)); out[o] <- xs; out
}

fmt_score <- function(v) { n <- suppressWarnings(as.numeric(v)); ifelse(is.na(n), "NA", sprintf("%.2f", n)) }

# small strand chevrons along an intron line
draw_strand_arrows <- function(xL, xR, y, strand) {
  if (xR - xL < 0.25) return(invisible())
  xsq <- seq(xL, xR, length.out = 14); xsq <- xsq[-c(1, length(xsq))]
  dx <- if (strand == "+") 0.035 else -0.035
  for (xx in xsq) {
    grid.lines(x = xi(c(xx - dx, xx)), y = yt(c(y - 0.028, y)), gp = gpar(col = "#C4C8CC", lwd = 0.5))
    grid.lines(x = xi(c(xx - dx, xx)), y = yt(c(y + 0.028, y)), gp = gpar(col = "#C4C8CC", lwd = 0.5))
  }
}

# ---------------------------------------------------------------------------
# 4. one gene block
# ---------------------------------------------------------------------------
plot_one_gene <- function(gene, block_y) {
  gid <- sym2gid[[gene]]
  gg  <- qc[qc$gene == gene, , drop = FALSE]
  bt  <- if ("biotype" %in% names(gg)) gg$biotype[1] else ""
  if (is.null(gid) || is.na(gid) || is.null(exbg[[gid]])) {
    plotText(label = sprintf("%s  (%s)  — no gene model", gene, bt),
             x = MARGIN_L, y = block_y + OFF_MODEL, just = c("left", "center"),
             fontsize = 7, fontcolor = GREY, default.units = "inches"); return(invisible())
  }

  uex <- reduce(exbg[[gid]])
  chrom  <- as.character(seqnames(uex))[1]
  strand <- as.character(strand(uex))[1]
  ex_s <- start(uex); ex_e <- end(uex)
  gmin <- min(ex_s); gmax <- max(ex_e)
  pad  <- max(round((gmax - gmin) * 0.03), 120L)
  cs <- gmin - pad; ce <- gmax + pad
  cds_rng <- if (!is.null(cdsbg[[gid]])) c(min(start(cdsbg[[gid]])), max(end(cdsbg[[gid]]))) else NULL
  is_lnc  <- is.null(cds_rng)
  niso    <- if (!is.na(n_iso[gid])) n_iso[[gid]] else NA

  ymodel <- block_y + OFF_MODEL

  # --- title ---
  plotText(label = sprintf("%s", gene), x = MARGIN_L, y = block_y + OFF_TITLE,
           just = c("left", "top"), fontsize = 8.5, fontface = "bold.italic",
           default.units = "inches")
  cc  <- if ("control_class" %in% names(gg)) gg$control_class[1] else "target"
  cctag <- if (!is.na(cc) && cc != "target") sprintf("  ·  [%s]", cc) else ""
  meta <- sprintf("  %s  ·  %s:%s-%s  ·  %s strand  ·  %s nt  ·  %s%s  ·  %d guides%s",
                  bt, chrom, format(gmin, big.mark = ","), format(gmax, big.mark = ","),
                  strand, format(gmax - gmin, big.mark = ","),
                  if (is.na(niso)) "?" else niso, if (!is.na(niso) && niso == 1) " isoform" else " isoforms",
                  nrow(gg), cctag)
  plotText(label = meta, x = MARGIN_L + nchar(gene) * 0.078 + 0.10, y = block_y + OFF_TITLE + 0.006,
           just = c("left", "top"), fontsize = 6, fontcolor = GREY, default.units = "inches")

  # --- collapsed gene model: intron line + strand arrows + exon boxes ---
  xL <- x2page(gmin, cs, ce); xR <- x2page(gmax, cs, ce)
  grid.lines(x = xi(c(xL, xR)), y = yt(c(ymodel, ymodel)), gp = gpar(col = INTRON_C, lwd = 0.9))
  draw_strand_arrows(xL, xR, ymodel, strand)
  for (i in seq_along(ex_s)) {
    xs <- x2page(ex_s[i], cs, ce); xe <- x2page(ex_e[i], cs, ce)
    if (is_lnc) {
      grid.rect(x = xi((xs + xe) / 2), y = yt(ymodel), width = xi(max(xe - xs, 0.004)),
                height = unit(0.115, "inches"), just = "centre",
                gp = gpar(fill = LNC_C, col = NA, alpha = 0.75))
    } else {
      a <- cds_rng[1]; b <- cds_rng[2]
      s <- ex_s[i]; e <- ex_e[i]
      # 5' UTR portion
      if (s < min(e, a)) grid.rect(x = xi((xs + x2page(min(e, a), cs, ce)) / 2), y = yt(ymodel),
                                   width = xi(max(x2page(min(e, a), cs, ce) - xs, 0.003)),
                                   height = unit(0.075, "inches"), just = "centre",
                                   gp = gpar(fill = UTR_C, col = NA))
      # CDS portion
      cS <- max(s, a); cE <- min(e, b)
      if (cS < cE) grid.rect(x = xi((x2page(cS, cs, ce) + x2page(cE, cs, ce)) / 2), y = yt(ymodel),
                             width = xi(max(x2page(cE, cs, ce) - x2page(cS, cs, ce), 0.003)),
                             height = unit(0.135, "inches"), just = "centre",
                             gp = gpar(fill = CDS_C, col = NA, alpha = 0.9))
      # 3' UTR portion
      if (max(s, b) < e) grid.rect(x = xi((x2page(max(s, b), cs, ce) + xe) / 2), y = yt(ymodel),
                                   width = xi(max(xe - x2page(max(s, b), cs, ce), 0.003)),
                                   height = unit(0.075, "inches"), just = "centre",
                                   gp = gpar(fill = UTR_C, col = NA))
    }
  }

  # --- guides: markers on the model + de-overlapped annotation stack ---
  gg <- gg[order(gg$genomic_pos), , drop = FALSE]
  gx_true <- x2page(gg$genomic_pos, cs, ce)
  gx_lab  <- spread(gx_true, mingap = 0.52, lo = MARGIN_L + 0.10, hi = RX - 0.10)
  labs <- sub("^.*_", "", gg$guide)
  reg  <- if ("region_recomputed" %in% names(gg)) gg$region_recomputed else gg$region_pool
  for (i in seq_len(nrow(gg))) {
    xt <- gx_true[i]; xl <- gx_lab[i]
    # marker: short vertical line from just above the model up into the guide row
    grid.lines(x = xi(c(xt, xt)), y = yt(c(ymodel - 0.085, block_y + OFF_MARK)),
               gp = gpar(col = RED, lwd = 1.1))
    grid.rect(x = xi(xt), y = yt(block_y + OFF_MARK), width = unit(0.022, "inches"),
              height = unit(0.075, "inches"), just = "centre", gp = gpar(fill = RED, col = NA))
    # leader from marker top to the (nudged) annotation column
    grid.lines(x = xi(c(xt, xl)), y = yt(c(block_y + OFF_MARK - 0.02, block_y + OFF_LABEL + 0.06)),
               gp = gpar(col = LEAD_C, lwd = 0.4))
    plotText(label = labs[i], x = xl, y = block_y + OFF_LABEL, just = c("center", "bottom"),
             fontsize = 5.6, fontface = "bold", fontcolor = "black", default.units = "inches")
    plotText(label = reg[i], x = xl, y = block_y + OFF_REGION, just = c("center", "bottom"),
             fontsize = 4.6, fontcolor = "#555555", default.units = "inches")
    plotText(label = sprintf("T %s · C %s", fmt_score(gg$tiger[i]), fmt_score(gg$cas13design[i])),
             x = xl, y = block_y + OFF_SCORE, just = c("center", "bottom"),
             fontsize = 4.0, fontcolor = RED, default.units = "inches")
  }

  # --- inter-guide spacing (mature-mRNA nt) between adjacent on-reference guides ---
  onref <- gg[!is.na(gg$mrna_pos) & (is.na(gg$on_reference_tx) | gg$on_reference_tx == "True"), , drop = FALSE]
  if (nrow(onref) >= 2) {
    onref <- onref[order(onref$mrna_pos), , drop = FALSE]
    ox <- x2page(onref$genomic_pos, cs, ce)
    ys <- block_y + OFF_SPACE
    for (i in seq_len(nrow(onref) - 1)) {
      x1 <- ox[i]; x2 <- ox[i + 1]; if (x2 - x1 < 0.02) next
      grid.segments(xi(x1), yt(ys), xi(x2), yt(ys), arrow = arrow(ends = "both", length = unit(0.03, "inches")),
                    gp = gpar(col = "#AAAAAA", lwd = 0.5))
      # stagger consecutive labels above / below the arc line so short adjacent arcs don't collide
      ly <- if (i %% 2 == 1) ys - 0.015 else ys + 0.10
      lj <- if (i %% 2 == 1) "bottom" else "top"
      plotText(label = sprintf("%d nt", onref$mrna_pos[i + 1] - onref$mrna_pos[i]),
               x = (x1 + x2) / 2, y = ly, just = c("center", lj),
               fontsize = 4.0, fontcolor = GREY, default.units = "inches")
    }
  }

  # --- genome coordinate ruler ---
  tryCatch(
    plotGenomeLabel(chrom = chrom, chromstart = cs, chromend = ce, assembly = asm,
                    x = MARGIN_L, y = block_y + OFF_RULER, length = BLOCK_W,
                    just = c("left", "top"), default.units = "inches", fontsize = 6, scale = "bp"),
    error = function(e) message("  [ruler-fail] ", gene, ": ", conditionMessage(e)))
  invisible()
}

# compact one-row legend at the bottom of each page
draw_legend <- function() {
  y <- PAGE_H - 0.20
  items <- list(list(UTR_C, "UTR exon", 0.075), list(CDS_C, "CDS exon", 0.135),
                list(LNC_C, "lncRNA exon", 0.115))
  x <- MARGIN_L
  for (it in items) {
    grid.rect(x = xi(x + 0.06), y = yt(y), width = unit(0.12, "inches"),
              height = unit(it[[3]], "inches"), just = "centre", gp = gpar(fill = it[[1]], col = NA))
    plotText(label = it[[2]], x = x + 0.15, y = y, just = c("left", "center"),
             fontsize = 6, default.units = "inches"); x <- x + 0.15 + 0.017 * nchar(it[[2]]) * 6 + 0.18
  }
  grid.lines(x = xi(c(x + 0.02, x + 0.16)), y = yt(c(y, y)), gp = gpar(col = INTRON_C, lwd = 0.9))
  plotText(label = "intron", x = x + 0.22, y = y, just = c("left", "center"),
           fontsize = 6, default.units = "inches"); x <- x + 0.22 + 0.30
  grid.rect(x = xi(x + 0.05), y = yt(y), width = unit(0.03, "inches"),
            height = unit(0.14, "inches"), just = "centre", gp = gpar(fill = RED, col = NA))
  plotText(label = "guide (23 nt) · T=TIGER C=Cas13Design", x = x + 0.12, y = y,
           just = c("left", "center"), fontsize = 6, fontcolor = GREY, default.units = "inches")
}

# ---------------------------------------------------------------------------
# 5. render (paginated)
# ---------------------------------------------------------------------------
dir.create(dirname(OUTPDF), showWarnings = FALSE, recursive = TRUE)
pdf(OUTPDF, width = PAGE_W, height = PAGE_H)
n_pages <- ceiling(length(genes_order) / GENES_PER_PAGE)
for (p in seq_len(n_pages)) {
  pageCreate(width = PAGE_W, height = PAGE_H, default.units = "inches", showGuides = FALSE)
  # device-filling viewport for custom grid shapes (bottom-origin inches)
  pushViewport(viewport(x = unit(0, "npc"), y = unit(0, "npc"), just = c("left", "bottom"),
                        width = unit(PAGE_W, "inches"), height = unit(PAGE_H, "inches"),
                        xscale = c(0, PAGE_W), yscale = c(0, PAGE_H)))
  idx <- ((p - 1) * GENES_PER_PAGE + 1):min(p * GENES_PER_PAGE, length(genes_order))
  for (j in seq_along(idx)) plot_one_gene(genes_order[idx[j]], TOP + (j - 1) * BLOCK_H)
  draw_legend()
  popViewport()
  if (p %% 10 == 0 || p == n_pages) message("[page] ", p, "/", n_pages)
}
invisible(dev.off())
message("figure -> ", OUTPDF)
