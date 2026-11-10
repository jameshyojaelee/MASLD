#!/usr/bin/env Rscript
# CONF1 2026-11-09: TNFSF10 locus panel replacing Fig 2G/2H (RORA, FABP1). Copy of
# panel-deep-fixes-20261017T035633Z/fig2/code/fig2GH_locus.R; only OUT, UP and CONFIG changed.
# KEY MESSAGE: the strongest enzyme-only locus (ALT up, MASLD null) colocalizes with a liver TNFSF10 eQTL.
# Edited copy of code/originals/figS09_locus_zoom_disease.R (LZ_MATCHED_MAIN=1 layout).
# Kept: track order (GWAS / SuSiE PIP / eQTL / genes / ruler), the EUR LD r2 colouring
# (polyfun_eur reference, #E6E6E5 -> #1B2F5B), the COLOC row rule (PP.H4.susie > 0.5 in the
# r2 aggregate, same GWAS x eQTL gene), the hg19 knownGene transcript model, the log2FC ramp
# (#1565C0 -> white -> #C9265E, saturating at |log2FC| = 0.3; gray if padj >= 0.05 or |log2FC| < 0.05).
# Changes (deep review 2026-10-17, PANEL_SPECS 2G, 2H):
#   G: the PIP track drawn at render time by an inline L = 1 susie_rss fit is removed;
#      MVP_NAFLD_EUR has no stored fine-mapping record, and the slot says so (untested, not negative).
#   G/H: SNP.PP.H4 is printed as text, not plotted on the PIP axis.
#   G/H: GWAS label with build, N and lead p; eQTL source and N; x ticks in Mb.
#   G/H: gene names black italic; gene body fill from the 774-unit bulk fit
#        (bulk-corrected-controls-20261008T162308Z primary_source_controls deg_results.csv);
#        genes outside that DE gene universe drawn as outlines only.
#   G/H: diamonds hollow so the points under them stay visible.
#   H:   lead moved to the GWAS lead 2:88423510 (was 2:88424066, the COLOC top variant);
#        r2 recomputed to that lead; window +/-250 kb, the same width as G.
#   H:   PIP points drawn without stems; the GWS points are no longer redrawn in black.
# LD r2 to the lead reads only the lead's row of the LD matrix (the producer read the full matrix).
# Usage: Rscript fig2GH_locus.R   (renders both)
suppressPackageStartupMessages({
  library(data.table); library(plotgardener); library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db); library(AnnotationDbi); library(grid)
})
set.seed(20261017)
OUT  <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/fig2_tnfsf10"
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
UP   <- file.path(OUT, "upstream")
LD_BASE <- file.path(PROJ, "GWAS/finemapping/data/ld_ref/polyfun_eur")
COLOC_INPUT <- file.path(PROJ, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")
DEG_FILE <- file.path(PROJ, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
REG <- fread(file.path(PROJ, "GWAS/finemapping/config/gwas_registry.tsv"))

CONFIG <- list(
  TNFSF10 = list(panel = "2_TNFSF10", study = "UKBB_ALT", chr = 3L, lead = 172294500L, eqtl = "TNFSF10",
                 win_start = 172294500L - 250000L, win_end = 172294500L + 250000L,
                 gwas_file = file.path(UP, "gwas_window_UKBB_ALT_chr3.tsv"),
                 eqtl_file = file.path(UP, "eqtl_TNFSF10.tsv"), pip_mode = "coloc_targeted",
                 gwas_name = "UKBB ALT GWAS, EUR (hg19)")
)

coloc_all <- fread(COLOC_INPUT, select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp", "top_snp_PP", "method"))
deg <- fread(DEG_FILE, select = c("symbol", "logFC", "padj"))
deg <- deg[!is.na(symbol) & nzchar(symbol)][order(padj)][!duplicated(symbol)]
pip_all <- fread(file.path(UP, "combined_finemapping_coloc_targeted_GH_windows.tsv.gz"),
                 select = c("chromosome", "position", "study", "susie_pip", "susie_cs", "susie_converged"))
stored_studies <- unique(c(pip_all$study, fread(file.path(UP, "combined_finemapping_GH_windows.tsv.gz"), select = "study")$study))

LFC_SAT <- 0.3; LFC_MIN <- 0.05
deg_palette <- colorRampPalette(c("#1565C0", "#90CAF9", "#FFFFFF", "#F48FB1", "#C9265E"))(101)
gene_fill <- function(sym) {
  h <- deg[symbol == sym]
  if (!nrow(h)) return(list(fill = "#FFFFFF", state = "not in bulk DE gene universe", lfc = NA, padj = NA))
  if (is.na(h$padj) || h$padj >= 0.05 || abs(h$logFC) < LFC_MIN)
    return(list(fill = "#9E9E9E", state = "n.s.", lfc = h$logFC, padj = h$padj))
  v <- pmin(pmax(h$logFC / LFC_SAT, -1), 1)
  list(fill = deg_palette[round((v + 1) * 50) + 1], state = if (h$logFC > 0) "up" else "down", lfc = h$logFC, padj = h$padj)
}
ld_pal_vec <- colorRampPalette(c("#E6E6E5", "#1B2F5B"))(100)
LEAD_COL <- "#C2185B"; COLOC_COL <- "#C79100"; EQTL_COL <- "#E07B39"; CS_COL <- "#1565C0"

r2_row_to_lead <- function(chr, lead, ws, we) {
  blocks <- fread(file.path(LD_BASE, "approx_LD_blocks.txt"))
  cc <- chr
  b <- blocks[blocks$chr == cc & blocks$start <= lead & blocks$stop >= lead]
  stopifnot(nrow(b) == 1)
  stem <- file.path(LD_BASE, paste0("chr", chr), paste0(b$start, ".", b$stop), paste0(b$start, ".", b$stop))
  bim <- fread(paste0(stem, ".bim"), col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  k <- which.min(abs(bim$pos - lead))
  row <- scan(pipe(sprintf("sed -n '%dp;%dq' %s", k, k, shQuote(paste0(stem, ".ld")))), quiet = TRUE)
  stopifnot(length(row) == nrow(bim))
  dt <- data.table(position = bim$pos, r2 = pmin(row^2, 1))[position >= ws & position <= we]
  list(r2 = dt[, .(r2 = max(r2)), by = position], lead_ld_pos = bim$pos[k], lead_rsid = bim$rsid[k])
}

render <- function(GENE) {
  cfg <- CONFIG[[GENE]]
  CHR <- cfg$chr; WS <- cfg$win_start; WE <- cfg$win_end; LEAD <- cfg$lead
  gw <- fread(cfg$gwas_file)[chromosome == CHR & position >= WS & position <= WE]
  gw[, p := as.numeric(pval)]
  gw <- gw[is.finite(p) & p > 0][, .(p = min(p)), by = position]
  ld <- r2_row_to_lead(CHR, LEAD, WS, WE)
  gw <- merge(gw, ld$r2, by = "position", all.x = TRUE)
  gw[is.na(r2), r2 := 0]
  gw[, lp := -log10(p)]
  lead_p <- gw[position == LEAD]$p
  stopifnot(length(lead_p) == 1)
  eq <- fread(cfg$eqtl_file)[GeneSymbol == cfg$eqtl & POS >= WS & POS <= WE]
  eq <- eq[, .(p = min(as.numeric(PVAL))), by = .(position = POS)][, lp := -log10(p)]
  eqtl_n <- unique(fread(cfg$eqtl_file, select = "N")$N)
  hit <- coloc_all[gwas_name == cfg$study & gene == cfg$eqtl & is.finite(PP.H4.susie) & PP.H4.susie > 0.5]
  stopifnot(nrow(hit) == 1)
  ctop <- as.integer(sub("^[0-9]+:", "", hit$top_snp))
  gwas_N <- REG[study_name == cfg$study]$N_tot
  gwas_cases <- REG[study_name == cfg$study]$N_cases
  if (cfg$pip_mode == "none") {
    stopifnot(!cfg$study %in% stored_studies)
    pips <- data.table(position = integer(), susie_pip = double(), susie_cs = integer())
  } else {
    pips <- pip_all[study == cfg$study & chromosome == CHR & position >= WS & position <= WE & susie_converged == TRUE &
                      !is.na(susie_pip) & susie_pip > 0, .(position, susie_pip, susie_cs)]
  }
  # ---- page geometry (inches from top)
  MARGIN_L <- 0.62; PLOT_W <- 2.50; MARGIN_R <- 0.14; PAGE_W <- MARGIN_L + PLOT_W + MARGIN_R
  y0 <- 0.05
  y1 <- 0.36; h1 <- 0.74
  y3 <- y1 + h1 + 0.08; h3 <- 0.36
  y4 <- y3 + h3 + 0.08; h4 <- 0.80
  y5 <- y4 + h4 + 0.06; h5 <- 0.55
  y6 <- y5 + h5 + 0.02
  PAGE_H <- y6 + 0.30
  pdf_file <- file.path(OUT, sprintf("Fig%s_%s_%s.pdf", cfg$panel, GENE, cfg$study))
  cairo_pdf(pdf_file, width = PAGE_W, height = PAGE_H, family = "Helvetica")
  pageCreate(width = PAGE_W, height = PAGE_H, default.units = "inches", showGuides = FALSE, xgrid = 0, ygrid = 0)
  G <- function(y_top) PAGE_H - y_top
  ux <- function(x) unit(x, "inches")
  xin <- function(pos) MARGIN_L + (pmax(WS, pmin(WE, pos)) - WS) / (WE - WS) * PLOT_W
  txt <- function(label, x, y, just = c("left", "center"), face = "plain", rot = 0)
    grid.text(label, ux(x), ux(G(y)), just = just, rot = rot,
              gp = gpar(fontsize = 6, fontfamily = "Helvetica", fontface = face, col = "black"))
  track_axis <- function(y, h, ymax, ticks, title) {
    yb <- G(y + h); yt <- G(y)
    to_y <- function(v) yb + v / ymax * (yt - yb)
    grid.segments(ux(MARGIN_L), ux(yb), ux(MARGIN_L), ux(yt), gp = gpar(col = "black", lwd = 0.5))
    grid.segments(ux(MARGIN_L), ux(yb), ux(MARGIN_L + PLOT_W), ux(yb), gp = gpar(col = "black", lwd = 0.5))
    for (v in ticks) {
      grid.segments(ux(MARGIN_L - 0.04), ux(to_y(v)), ux(MARGIN_L), ux(to_y(v)), gp = gpar(col = "black", lwd = 0.5))
      grid.text(format(v), ux(MARGIN_L - 0.06), ux(to_y(v)), just = "right",
                gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
    }
    grid.text(title, ux(MARGIN_L - 0.42), ux((yb + yt) / 2), rot = 90,
              gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
    to_y
  }
  marks <- list()
  addm <- function(...) marks[[length(marks) + 1]] <<- data.table(...)
  # ---- header
  n_txt <- if (gwas_cases > 0) sprintf("N = %s (%s cases)", format(gwas_N, big.mark = ","), format(gwas_cases, big.mark = ",")) else
    sprintf("N = %s", format(gwas_N, big.mark = ","))
  txt(sprintf("%s, %s", cfg$gwas_name, n_txt), MARGIN_L, y0 + 0.06)
  sig_note <- if (lead_p < 5e-8) "genome-wide significant" else "not genome-wide significant"
  txt(sprintf("lead %d:%s, p = %s, %s", CHR, format(LEAD, big.mark = ","), formatC(lead_p, format = "e", digits = 1), sig_note),
      MARGIN_L, y0 + 0.17)
  # ---- GWAS track
  ymax1 <- max(ceiling(max(gw$lp)) + 1, 9)
  ticks1 <- pretty(c(0, max(gw$lp, 7.3)))
  ticks1 <- ticks1[ticks1 <= ymax1]
  to_y1 <- track_axis(y1, h1, ymax1, ticks1, "GWAS −log10 p")
  sig_y <- to_y1(-log10(5e-8))
  grid.segments(ux(MARGIN_L), ux(sig_y), ux(MARGIN_L + PLOT_W), ux(sig_y), gp = gpar(col = "black", lwd = 0.4, lty = 2))
  txt("p = 5e-8", MARGIN_L + PLOT_W, G(sig_y) - 0.035, just = c("right", "bottom"))
  setorder(gw, r2)
  gw[, fill := ld_pal_vec[pmin(100, pmax(1, ceiling(r2 * 100)))]]
  grid.points(ux(xin(gw$position)), ux(to_y1(gw$lp)), pch = 19, size = ux(0.03), gp = gpar(col = gw$fill))
  addm(track = "gwas", position = gw$position, value = gw$lp, p = gw$p, r2 = gw$r2, colour = gw$fill)
  grid.points(ux(xin(LEAD)), ux(to_y1(-log10(lead_p))), pch = 23, size = ux(0.085), gp = gpar(col = LEAD_COL, fill = NA, lwd = 0.8))
  addm(track = "gwas", position = LEAD, value = -log10(lead_p), p = lead_p, mark = "lead GWAS variant", colour = LEAD_COL)
  cgw <- gw[position == ctop]
  if (nrow(cgw)) {
    grid.points(ux(xin(ctop)), ux(to_y1(cgw$lp)), pch = 23, size = ux(0.11), gp = gpar(col = COLOC_COL, fill = NA, lwd = 0.8))
    addm(track = "gwas", position = ctop, value = cgw$lp, p = cgw$p, mark = "COLOC top variant", colour = COLOC_COL)
  }
  # ---- PIP track
  if (cfg$pip_mode == "none") {
    yb <- G(y3 + h3); yt <- G(y3)
    grid.segments(ux(MARGIN_L), ux(yb), ux(MARGIN_L), ux(yt), gp = gpar(col = "black", lwd = 0.5))
    txt(sprintf("no stored fine-mapping record for %s", gsub("_", " ", cfg$study)), MARGIN_L + PLOT_W / 2, y3 + h3 / 2,
        just = c("center", "center"))
    grid.text("SuSiE PIP", ux(MARGIN_L - 0.42), ux((yb + yt) / 2), rot = 90,
              gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
    addm(track = "pip", mark = "no stored fine-mapping record", label = cfg$study)
  } else {
    to_y3 <- track_axis(y3, h3, 1, c(0, 0.5, 1), "SuSiE PIP")
    in_cs <- !is.na(pips$susie_cs) & pips$susie_cs > 0
    grid.points(ux(xin(pips$position[!in_cs])), ux(to_y3(pips$susie_pip[!in_cs])), pch = 21, size = ux(0.028),
                gp = gpar(col = CS_COL, fill = "white", lwd = 0.4))
    grid.points(ux(xin(pips$position[in_cs])), ux(to_y3(pips$susie_pip[in_cs])), pch = 19, size = ux(0.045),
                gp = gpar(col = CS_COL))
    addm(track = "pip", position = pips$position, value = pips$susie_pip, in_cs = in_cs, colour = CS_COL)
    lp_pip <- pips[position == LEAD]$susie_pip
    if (length(lp_pip)) grid.points(ux(xin(LEAD)), ux(to_y3(lp_pip)), pch = 23, size = ux(0.085),
                                    gp = gpar(col = LEAD_COL, fill = NA, lwd = 0.8))
    cp_pip <- pips[position == ctop]$susie_pip
    if (length(cp_pip)) grid.points(ux(xin(ctop)), ux(to_y3(cp_pip)), pch = 23, size = ux(0.11),
                                    gp = gpar(col = COLOC_COL, fill = NA, lwd = 0.8))
    if (length(lp_pip)) addm(track = "pip", position = LEAD, value = lp_pip, mark = "lead GWAS variant", colour = LEAD_COL)
    if (length(cp_pip)) addm(track = "pip", position = ctop, value = cp_pip, mark = "COLOC top variant", colour = COLOC_COL)
  }
  # ---- eQTL track
  ymax4 <- ceiling(max(eq$lp) * 1.55)
  to_y4 <- track_axis(y4, h4, ymax4, pretty(c(0, max(eq$lp))), "eQTL −log10 p")
  grid.points(ux(xin(eq$position)), ux(to_y4(eq$lp)), pch = 19, size = ux(0.03), gp = gpar(col = EQTL_COL))
  addm(track = "eqtl", position = eq$position, value = eq$lp, p = eq$p, colour = EQTL_COL)
  ceq <- eq[position == ctop]
  if (nrow(ceq)) {
    grid.points(ux(xin(ctop)), ux(to_y4(ceq$lp)), pch = 23, size = ux(0.11), gp = gpar(col = COLOC_COL, fill = NA, lwd = 0.8))
    addm(track = "eqtl", position = ctop, value = ceq$lp, p = ceq$p, mark = "COLOC top variant", colour = COLOC_COL)
  }
  txt(sprintf("Broadaway liver eQTL, %s (EUR, N = %s)", cfg$eqtl, format(eqtl_n, big.mark = ",")), MARGIN_L + 0.04, y4 + 0.02,
      just = c("left", "top"))
  txt(sprintf("COLOC PP.H4 = %.3f (SuSiE)", hit$PP.H4.susie), MARGIN_L + 0.04, y4 + 0.11, just = c("left", "top"))
  txt(sprintf("top variant %s, SNP.PP.H4 = %.2f", hit$top_snp, hit$top_snp_PP), MARGIN_L + 0.04, y4 + 0.20,
      just = c("left", "top"))
  addm(track = "text", label = c(sprintf("COLOC PP.H4 = %.2f", hit$PP.H4.susie), sprintf("SNP.PP.H4 = %.2f", hit$top_snp_PP),
                                 sprintf("GWAS N = %d", gwas_N), sprintf("eQTL N = %d", eqtl_n),
                                 sprintf("lead p = %.3g", lead_p)),
       value = c(hit$PP.H4.susie, hit$top_snp_PP, gwas_N, eqtl_n, lead_p))
  # ---- gene track
  locus_gr <- GRanges(paste0("chr", CHR), IRanges(WS, WE))
  gwin <- subsetByOverlaps(suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                                  filter = list(tx_chrom = paste0("chr", CHR)))), locus_gr)
  ex_by_gene <- suppressMessages(exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene"))
  gmod <- list()
  for (i in seq_along(gwin)) {
    gid <- as.character(gwin$gene_id[i])
    sym <- suppressMessages(mapIds(org.Hs.eg.db, gid, "SYMBOL", "ENTREZID"))
    if (is.na(sym) || is.null(ex_by_gene[[gid]])) next
    exr <- reduce(ex_by_gene[[gid]])
    if (max(end(exr)) < WS || min(start(exr)) > WE) next
    gf <- gene_fill(sym)
    gmod[[length(gmod) + 1]] <- list(sym = sym, strand = as.character(strand(exr))[1], gs = min(start(exr)), ge = max(end(exr)),
                                     es = start(exr), ee = end(exr), fill = gf$fill, state = gf$state, lfc = gf$lfc, padj = gf$padj)
  }
  gmod <- gmod[order(sapply(gmod, `[[`, "gs"))]
  row_end <- numeric(0)
  for (j in seq_along(gmod)) {
    gs <- max(WS, gmod[[j]]$gs); placed <- FALSE
    for (r in seq_along(row_end)) if (gs > row_end[r] + 0.02 * (WE - WS)) {
      gmod[[j]]$row <- r; row_end[r] <- min(WE, gmod[[j]]$ge); placed <- TRUE; break }
    if (!placed) { row_end <- c(row_end, min(WE, gmod[[j]]$ge)); gmod[[j]]$row <- length(row_end) }
  }
  nrows <- max(length(row_end), 1)
  stopifnot(nrows <= 2)
  row_y <- function(r) y5 + 0.12 + (r - 1) * 0.24
  for (g in gmod) {
    yc <- row_y(g$row); yg <- G(yc)
    x0 <- xin(g$gs); x1 <- xin(g$ge)
    outline <- g$state == "not in bulk DE gene universe"
    grid.lines(ux(c(x0, x1)), ux(c(yg, yg)), gp = gpar(col = if (outline) "black" else g$fill, lwd = 0.8))
    for (k in seq_along(g$es)) {
      ex0 <- xin(g$es[k]); ex1 <- xin(g$ee[k])
      grid.rect(ux(ex0), ux(yg), width = ux(max(ex1 - ex0, 0.008)), height = ux(0.08), just = c("left", "center"),
                gp = gpar(col = if (outline) "black" else NA, lwd = 0.3, fill = g$fill))
    }
    addm(track = "gene", label = g$sym, start = g$gs, end = g$ge, row = g$row, colour = g$fill, state = g$state,
         lfc = g$lfc, padj = g$padj, strand = g$strand)
  }
  # gene symbols: centred above the body of the same row; if that overlaps a label already
  # placed above this row, put it below the body instead
  placed <- list()
  gmod <- gmod[order(sapply(gmod, function(g) (max(WS, g$gs) + min(WE, g$ge)) / 2))]
  for (g in gmod) {
    w <- convertWidth(grobWidth(textGrob(g$sym, gp = gpar(fontsize = 6, fontface = "italic", fontfamily = "Helvetica"))),
                      "inches", valueOnly = TRUE)
    xm <- (xin(g$gs) + xin(g$ge)) / 2
    xm <- min(max(xm, MARGIN_L + w / 2), MARGIN_L + PLOT_W + MARGIN_R - w / 2 - 0.01)
    yc <- row_y(g$row)
    cand <- c(yc - 0.085, yc + 0.085)
    ok <- vapply(cand, function(yy) all(vapply(placed, function(P) abs(P$y - yy) > 0.05 ||
      (xm - w / 2 > P$x1 + 0.03) || (xm + w / 2 < P$x0 - 0.03), logical(1))), logical(1))
    yy <- if (ok[1]) cand[1] else cand[2]
    placed[[length(placed) + 1]] <- list(y = yy, x0 = xm - w / 2, x1 = xm + w / 2)
    txt(g$sym, xm, yy, just = c("center", "center"), face = "italic")
  }
  grid.text("Genes", ux(MARGIN_L - 0.42), ux(G(y5 + h5 / 2)), rot = 90,
            gp = gpar(fontsize = 6, fontfamily = "Helvetica", col = "black"))
  # ---- Mb ruler
  brks <- pretty(c(WS, WE) / 1e6, n = 5)
  brks <- brks[brks * 1e6 >= WS & brks * 1e6 <= WE]
  grid.segments(ux(MARGIN_L), ux(G(y6)), ux(MARGIN_L + PLOT_W), ux(G(y6)), gp = gpar(col = "black", lwd = 0.5))
  for (b in brks) {
    grid.segments(ux(xin(b * 1e6)), ux(G(y6)), ux(xin(b * 1e6)), ux(G(y6 + 0.04)), gp = gpar(col = "black", lwd = 0.5))
    txt(format(b, nsmall = 1), xin(b * 1e6), y6 + 0.10, just = c("center", "center"))
  }
  txt(sprintf("chr%d position (Mb, hg19)", CHR), MARGIN_L + PLOT_W / 2, y6 + 0.22, just = c("center", "center"))
  addm(track = "ruler", position = as.integer(round(brks * 1e6)), value = brks)
  dev.off()
  src <- rbindlist(marks, fill = TRUE)
  src[, `:=`(panel = cfg$panel, gene_locus = GENE, gwas = cfg$study, chr = CHR, win_start = WS, win_end = WE,
             lead_pos = LEAD, lead_ld_pos = ld$lead_ld_pos, lead_rsid = ld$lead_rsid, coloc_top = ctop,
             x_in = fifelse(is.na(position), NA_real_, MARGIN_L + (position - WS) / (WE - WS) * PLOT_W),
             plot_x0_in = MARGIN_L, plot_w_in = PLOT_W)]
  fwrite(src, file.path(OUT, "source_tables", paste0(cfg$panel, ".tsv")), sep = "\t")
  cat(sprintf("%s: %s GWAS %d, PIP %d, eQTL %d, genes %d; lead p %.3g; LD row %d (%s); COLOC top %d PP.H4 %.4f SNP.PP.H4 %.4f\n",
              cfg$panel, pdf_file, nrow(gw), nrow(pips), nrow(eq), length(gmod), lead_p, ld$lead_ld_pos, ld$lead_rsid,
              ctop, hit$PP.H4.susie, hit$top_snp_PP))
}
for (g in names(CONFIG)) render(g)
writeLines(capture.output(sessionInfo()), file.path(OUT, "logs/sessionInfo_fig2GH.txt"))
