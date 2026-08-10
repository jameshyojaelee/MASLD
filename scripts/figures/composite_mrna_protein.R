#!/usr/bin/env Rscript
# KEY MESSAGE: A fixed, previously selected set of liver protein processes shows
# coherent same-cohort mRNA/protein remodeling, presented descriptively rather
# than as independent validation of the Figure 3 modules.
# composite_mrna_protein.R — FIG 4c descriptive same-cohort refit, full-row proteomics composite. One shared 25-protein axis
# (rows), grouped into 5 previously selected MASLD protein processes (lipid-droplet, ECM/basement-membrane, immune-
# inflammation, detox, amino-acid-catabolism), read left->right through complementary views:
#   PANEL 1 (heatmap)  : protein abundance z-score across PXD051911 liver patients ORDERED BY SEVERITY
#                        (fibrosis then NAS), with NAS + fibrosis colour strips on top.
#   PANEL 2 (pie)      : protein-vs-HISTOLOGY correlation (Spearman of abundance vs each pathologist
#                        score: steatosis/ballooning/inflammation/fibrosis/NAS) — feature specificity
#                        (lipid->steatosis, ECM->fibrosis).
#   PANEL 3 (paired dots): mRNA vs protein log2FC (amplification / buffering).
#   PANEL 4 (NES bars) : pathway-level fgsea NES for all 6 programs, incl. Mito/OXPHOS which is the
#                        2nd-strongest signal but has no honest gene-level representatives (shown NES-only).
# Gene labels + program colour-chip live ONCE on the far left (shared). Strips occupy rows n+1,n+2 on the
# heatmap; the pie & dot panel leave that band blank so ALL gene rows align across the first three panels
# (identical y-limit c(0.5, ytop); the dot panel's internal y-scale is replaced to match).
# Everything from PXD051911 liver DIA-MS (n=58). The exact 25 display rows are
# frozen from the protected original panel, but all protein estimates, pathway
# statistics, patient abundances, and histology associations are rebuilt with
# acquisition-batch-aware source-compatible normalization.
# Acquisition batch remains in every adjusted model and diagnostic sidecar but
# is intentionally not displayed as an annotation strip.
# Output: FIG4_DIR/panels/fig4c_mrna_protein_composite.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(42)
FAM <- "Helvetica"   # PROJECT STANDARD: every text element in this figure is 6pt Helvetica, face "plain"
                     # (italics reserved for gene/module symbols only) — set explicitly since geom_text()/
                     # annotate("text",...) do NOT inherit family from theme_void()/theme(text=...).
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/composite_helpers.R"))
DATA_DIR <- file.path(FIG4_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
PROTEIN_CONTEXT <- file.path(PROGRAM_CONTEXT_DIR, "proteomics")

# ── corrected logFC pairs + nuisance-adjusted liver DIA-MS matrix ─────────────
con <- fread(file.path(PROTEIN_CONTEXT, "panel4c_mrna_protein.tsv"))
ab <- fread(file.path(PROTEIN_CONTEXT, "panel4c_adjusted_abundance.tsv"))
abmat <- as.matrix(ab[, -"gene"])
rownames(abmat) <- ab$gene
storage.mode(abmat) <- "double"

# ── programs + genes: fixed descriptive row contract ─────────────────────────────────────────────────────
# These original 5 x 5 rows were selected in PXD051911 using protein significance
# and mRNA/protein direction concordance. They are retained as a descriptive
# process display and are never reselected during the batch-aware refit.
sel <- fread(file.path(PROTEIN_CONTEXT, "fixed_panel4c_rows.tsv"))
setorder(sel, prog_order, gene_order)
ord <- sel$gene; genes <- ord; group_of <- setNames(sel$program, sel$gene)
gl <- unique(sel$program)                                  # program display order (top -> bottom)
stopifnot(all(ord %in% rownames(abmat)))                   # every selected protein has an abundance row
cand <- con[gene %in% genes]
cand[, both_sig := bulk_sig & protein_sig]
n <- length(ord); ytop <- n + 4.2                          # genes 1..n; NAS/fibrosis strips n+2:n+3
g <- unname(group_of[ord]); bnd <- which(g[-1] != g[-n])
# 2026-07-09: reverted to the original custom hues -- the Okabe-Ito swap (+ exact
# RdBu heatmap poles) read as "hideous" / no longer like a regular heatmap in
# review. Keeping ink/stroke/font fixes; palette reverted pending a better direction.
prog_cols <- c("Lipid droplet" = "#3B6EA5", "ECM / BM" = "#B05A5A", "Immune / inflammation" = "#6E5AA5",
               "Detox" = "#4E9E76", "AA catabolism" = "#C79A3E")
gcol <- prog_cols[gl]

# ── batch-aware metadata + partial histology associations ────────────────────
mh <- fread(file.path(PROTEIN_CONTEXT, "panel4c_metadata.tsv"))
setnames(mh, "sample_id", "raw")
mh <- mh[raw %in% colnames(abmat)]
feat_levels <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
hist <- fread(file.path(PROTEIN_CONTEXT, "panel4c_histology_partial.tsv"))
cmat <- as.matrix(dcast(hist[gene %in% ord], gene ~ feature, value.var = "rho")
                  [match(ord, gene), ..feat_levels])
rownames(cmat) <- ord

# ── PANEL 1: severity heatmap ─────────────────────────────────────────────────
sev <- copy(mh); setorder(sev, Fibrosis, NAS); pts <- sev$raw; np <- length(pts)
z <- t(scale(t(abmat[ord, pts, drop = FALSE]))); z <- pmax(pmin(z, 2), -2)
hm <- CJ(gi = seq_len(n), pj = seq_len(np)); hm[, `:=`(zz = z[cbind(gi, pj)], ypos = n - gi + 1)]
nas_pal <- colorRampPalette(c("#EFEDF5", "#54278F"))(9); fib_pal <- colorRampPalette(c("#FEE6CE", "#8C2D04"))(5)
strips <- rbind(
  data.table(pj = seq_len(np), ypos = n + 3, hex = fib_pal[sev$Fibrosis + 1]),
  data.table(pj = seq_len(np), ypos = n + 2, hex = nas_pal[pmin(sev$NAS, 8) + 1])
)
bands <- rbindlist(lapply(gl, function(p) { idx <- which(g == p)
  data.table(program = factor(p, levels = gl), ymin = n - max(idx) + 0.5, ymax = n - min(idx) + 1.5,
             sidehex = gcol[p], cy = n - mean(idx) + 1) }))
# 2-line wrap for long program names: keeps the label's rendered HEIGHT well inside its own band
# (~5 rows tall) regardless of string length, unlike a single diagonal line whose height scales with
# string length and spills into neighbouring brackets once >~2 words.
wrap_lab <- c("Lipid droplet" = "Lipid\ndroplet", "ECM / BM" = "ECM / BM",
              "Immune / inflammation" = "Immune /\ninflammation",
              "Detox" = "Detox", "AA catabolism" = "AA\ncatabolism")
bands[, label := wrap_lab[as.character(program)]]
# 2026-07-09: pathway-level fgsea NES now rides on the LEFT program brackets (arrow-bracket encoding),
# replacing the old far-right NES bar column (former Panel 4). Read + merge NES onto `bands` here so `hp`
# can draw it.
enr <- fread(file.path(PROTEIN_CONTEXT, "panel4c_process_enrichment.tsv"))[representable == TRUE]
bands <- merge(bands, enr[, .(program = factor(program, levels = gl), NES, padj)], by = "program", sort = FALSE)
setorder(bands, program)                                   # restore gl (top -> bottom) order after merge
# program grouping brackets: "[" opening right. Order across the left gutter =
# label | bracket | arrow | NES/q | (gap) | gene rows -- i.e. the bracket sits between the
# pathway name (to its left) and the NES arrow (to its right), not next to the genes.
XB <- -20.0
brk <- rbind(bands[, .(x = XB,        xend = XB,        y = ymin + 0.2, yend = ymax - 0.2, col = sidehex)],
             bands[, .(x = XB,        xend = XB + 0.45, y = ymax - 0.2, yend = ymax - 0.2, col = sidehex)],
             bands[, .(x = XB,        xend = XB + 0.45, y = ymin + 0.2, yend = ymin + 0.2, col = sidehex)])
# ── arrow-bracket NES glyphs: one filled triangle per program just RIGHT of the bracket spine. Apex UP =
# enriched (+NES) / DOWN = depleted (-NES); SIZE proportional to |NES| FROM ZERO (so the narrow 2.0-2.7
# spread is never exaggerated); exact NES printed in black. SCALE_MAG=FALSE -> fixed-size (sign+number only).
SCALE_MAG <- TRUE
XG <- XB + 2.3                                            # arrow centre x, RIGHT of the bracket; NES/q sit to the arrow's right
HH_MAX <- 0.58; WW_MAX <- 1.2                            # half-height (rows) / half-width (x-units): a row renders ~1.7x
                                                          # taller than an x-unit here, so width > height => proper arrowhead
bands[, sz := if (SCALE_MAG) abs(NES) / max(abs(NES)) else 1]
bands[, s  := sign(NES)]                                   # +1 apex-up (enriched) / -1 apex-down (depleted)
chev <- bands[, .(x = c(XG, XG - WW_MAX * sz, XG + WW_MAX * sz),
                  y = c(cy + s * HH_MAX * sz, cy - s * HH_MAX * sz, cy - s * HH_MAX * sz)),
                by = program]
# hp absorbs the freed nesp column (2.65 -> 3.30 width units). KEEP_TILES=TRUE holds the heatmap block
# pixel-identical by pushing xlim left in proportion (right edge np+0.5 fixed) -> freed space becomes gutter;
# KEEP_TILES=FALSE lets the tiles grow to fill (xlim unchanged) -> smaller gutter, no blank left margin.
KEEP_TILES <- FALSE   # chosen layout = minimal whitespace: heatmap tiles grow ~18% to fill the freed nesp
                      # column so there is no blank left margin (pies + lollipop stay pixel-identical either way).
                      # TRUE would instead hold the tiles pixel-identical and leave the freed width as blank gutter.
HP_W_OLD <- 2.65; HP_W_NEW <- 3.30
xL <- if (KEEP_TILES) (np + 0.5) - (np + 28.5) * (HP_W_NEW / HP_W_OLD) else (XB - 12.5)

# manual horizontal colour-bar drawn (clip=off) right under its OWN panel, so each legend hugs its panel
# instead of being pushed to the figure bottom by the shared x-axis reservation. Tick/title offsets
# tightened (was 1.05) to trim the below-zero footprint now that the panel is taller (more pt/data-unit).
grad_bar <- function(x0, x1, y0, cols, brks, blim, title, bh = 0.5) {
  cr <- grDevices::colorRamp(cols); u <- seq(0, 1, length.out = 64)
  d <- data.table(x = x0 + u * (x1 - x0), y = y0, hex = grDevices::rgb(cr(u), maxColorValue = 255))
  bx <- x0 + (brks - blim[1]) / (blim[2] - blim[1]) * (x1 - x0)
  list(geom_tile(data = d, aes(x = x, y = y), fill = d$hex, width = (x1 - x0) / 63, height = bh),
       annotate("text", x = bx, y = y0 - 0.85, label = as.character(brks), size = 6/.pt, family = FAM),
       annotate("text", x = (x0 + x1) / 2, y = y0 - 1.9, label = title, size = 6/.pt, family = FAM))
}
hp <- ggplot() +
  # heatmap + severity strips rasterised at 600 dpi (crisp print-res image, no low-DPI
  # preflight flag) while all text / lines / pies stay vector. geom_raster is already a
  # single bitmap but embeds at native 58x25 (~23 dpi at panel size); rasterise upsamples
  # the solid-colour cells losslessly (interpolate stays off, so cell edges stay sharp).
  ggrastr::rasterise(geom_raster(data = hm, aes(x = pj, y = ypos, fill = zz)), dpi = 600) +
  ggrastr::rasterise(geom_tile(data = strips, aes(x = pj, y = ypos), fill = strips$hex, height = 0.86), dpi = 600) +
  # 2026-07-09 house re-skin: 0.8 (~1.7pt rendered, flagged in the style handoff)
  # -> STROKE_HAIRLINE (0.25pt); the color-coded bracket itself still identifies
  # the program without needing a heavy line.
  geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = yend),
               colour = brk$col, linewidth = 0.35/.pt, lineend = "round") +   # slightly > STROKE_HAIRLINE (0.25) per request
  geom_polygon(data = chev, aes(x = x, y = y, group = program), fill = house_ink, colour = NA) +   # arrows = black (data glyph); brackets keep program colour
  geom_text(data = bands, aes(x = XG + 1.8, y = cy + 0.55, label = sprintf("%+.1f", NES)),
            hjust = 0, size = 6/.pt, family = FAM, colour = house_ink) +
  geom_text(data = bands, aes(x = XG + 1.8, y = cy - 0.65,
                              label = paste0("q=", gsub("e-0", "e-", sprintf("%.0e", padj)))),
            hjust = 0, size = 6/.pt, family = FAM, colour = house_ink) +
  annotate("text", x = XG + 1.8, y = n + 1.2, label = "adjusted NES", hjust = 0, size = 6/.pt, family = FAM, colour = house_ink) +
  geom_text(data = bands, aes(x = XB - 0.8, y = cy, label = label), hjust = 1, vjust = 0.5,
            lineheight = 0.82, size = 6/.pt, family = FAM, colour = house_ink) +
  # 2026-07-09: reverted to the original moderate red/white/blue -- the exact RdBu
  # extremes (#67001f/#053061) read as too dark/heavy, no longer "a regular heatmap".
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2, 2),
                       na.value = "#D9D9D9", guide = "none") +
  grad_bar(np * 0.30, np * 0.55, -0.6, c("#2166AC", "white", "#B2182B"), c(-2, 0, 2), c(-2, 2), "Protein z-score", bh = 0.5) +
  geom_hline(yintercept = n - bnd + 0.5, linetype = "dashed", linewidth = 0.22, colour = "grey55") +
  annotate("text", x = -0.25, y = n:1, label = ord, hjust = 1, size = 6/.pt, fontface = "italic", family = FAM, colour = house_ink) +
  annotate("text", x = -0.25, y = c(n + 2, n + 3), label = c("NAS", "Fibrosis"), hjust = 1, size = 6/.pt, family = FAM, colour = house_ink) +
  annotate("segment", x = 0.5, xend = np, y = n + 3.65, yend = n + 3.65,
           arrow = arrow(length = unit(0.045, "in"), type = "closed"), linewidth = 0.3, colour = "grey45") +
  annotate("text", x = 0.5, y = n + 3.90, label = "increasing severity", hjust = 0, vjust = 0, size = 6/.pt, colour = "grey35", family = FAM) +
  coord_cartesian(xlim = c(xL, np + 0.5), ylim = c(0.5, ytop), expand = FALSE, clip = "off") +
  theme_void() +
  theme(text = element_text(family = FAM, face = "plain"),
        plot.title = element_text(size = 6, family = FAM, face = "plain", hjust = 0, margin = margin(0, 0, 2, 0)),
        legend.position = "none", plot.margin = margin(1, 0, 1, 2))

# ── PANEL 2: protein vs histology — PAC-MAN pie matrix drawn as MANUAL POLYGONS under coord_cartesian.
# Why manual: coord_fixed (needed for circular pies) letterboxes inside patchwork and breaks row-alignment.
# Under plain cartesian the panel FILLS its cell exactly like the heatmap/lollipop (alignment guaranteed);
# each wedge's x/y radius is scaled by RX/RY so the pies still render circular (RY/RX = cell aspect; tune).
nf <- length(feat_levels)
RY <- 0.46; ASP <- 2.36; RX <- RY / ASP                    # RY≈half a row; ASP = (row height : col width),
                                                            # retuned 1.94 -> 2.36 for the 5.5in canvas (height unchanged at 2.83in),
                                                            # which narrows each column while row height is fixed
mkcirc <- function(cx, cy, np = 44) { th <- seq(0, 2 * pi, length.out = np)
  data.table(x = cx + RX * cos(th), y = cy + RY * sin(th)) }
mkwedge <- function(cx, cy, frac) { if (!is.finite(frac)) frac <- 1
  if (frac < 1e-6) return(NULL)
  t <- seq(0, frac * 2 * pi, length.out = max(4, round(44 * frac) + 2)); a <- pi/2 - t  # top, clockwise
  rbind(data.table(x = cx, y = cy), data.table(x = cx + RX * cos(a), y = cy + RY * sin(a))) }
circ <- list(); wed <- list(); k <- 0L
for (gi in seq_len(n)) for (fj in seq_len(nf)) { k <- k + 1L
  cy <- n - gi + 1; r <- cmat[ord[gi], feat_levels[fj]]
  circ[[k]] <- cbind(mkcirc(fj, cy), gid = k)
  w <- mkwedge(fj, cy, abs(r)); if (!is.null(w)) wed[[k]] <- cbind(w, gid = k, rr = r) }
circ <- rbindlist(circ); wed <- rbindlist(wed)
bbands <- rbindlist(lapply(gl, function(p) { idx <- which(g == p)
  data.table(ymin = n - max(idx) + 0.5, ymax = n - min(idx) + 1.5, tint = paste0(gcol[p], "33")) }))
pp <- ggplot() +
  geom_rect(data = bbands, aes(xmin = 0.5, xmax = nf + 0.5, ymin = ymin, ymax = ymax), fill = bbands$tint) +
  geom_hline(yintercept = n - bnd + 0.5, linetype = "dashed", linewidth = 0.22, colour = "grey55") +
  geom_polygon(data = circ, aes(x = x, y = y, group = gid), fill = "white", colour = "grey80", linewidth = 0.12) +
  geom_polygon(data = wed, aes(x = x, y = y, group = gid, fill = rr), colour = NA) +
  scale_fill_gradient2(low = PIE_NEG, mid = "white", high = PIE_POS, midpoint = 0, limits = c(-1, 1),
                       na.value = "#D9D9D9",
                       guide = "none") +
  grad_bar(0.9, nf + 0.1, -0.6, c(PIE_NEG, "white", PIE_POS), c(-1, 0, 1), c(-1, 1), "Partial Spearman rho", bh = 0.5) +
  annotate("text", x = seq_len(nf), y = n + 1.9, label = feat_levels, angle = 45, hjust = 0, size = 6/.pt, family = FAM) +
  coord_cartesian(xlim = c(0.5, nf + 0.5), ylim = c(0.5, ytop), expand = FALSE, clip = "off") +
  theme_void() +
  theme(text = element_text(family = FAM, face = "plain"),
        plot.title = element_text(size = 6, family = FAM, face = "plain", hjust = 0, margin = margin(0, 0, 2, 0)),
        legend.position = "none", plot.margin = margin(1, 10, 1, 1))

# ── PANEL 3: mRNA vs protein paired dots (y-scale replaced so 24 rows align) ──
# House style bans lollipop geometry, so this is composite_helpers::lollipop_panel()
# with the two zero-anchored stems dropped. Every plotted value is unchanged: the
# two estimates stay at the same x, the same shape/size/stroke, the same
# significance-driven fill, and the same x=0 reference line and group separators.
dot_panel <- function(dt, item_col, ord, est1, est2, est1_lab, est2_lab, xlab,
                      col1 = "#8EC7E2", col2 = "#E8853A", group_of = NULL,
                      ref0 = TRUE, sig1 = NULL, sig2 = NULL) {
  d <- copy(dt)[match(ord, get(item_col))]
  n <- length(ord)
  d[, ypos := n:1]                                   # row 1 at top, aligned to matrix
  d[, e1 := as.numeric(get(est1))]; d[, e2 := as.numeric(get(est2))]
  d[, sig1_plot := if (is.null(sig1)) TRUE else as.logical(get(sig1))]
  d[, sig2_plot := if (is.null(sig2)) TRUE else as.logical(get(sig2))]
  d[is.na(sig1_plot), sig1_plot := FALSE]
  d[is.na(sig2_plot), sig2_plot := FALSE]

  p <- ggplot(d)
  if (ref0) p <- p + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85")
  if (!is.null(group_of)) {                          # group separators mirror the matrix
    g <- unname(group_of[ord]); bnd <- which(g[-1] != g[-n])
    p <- p + geom_hline(yintercept = n - bnd + 0.5, linetype = "dashed",
                        linewidth = 0.25, colour = "grey45")
  }
  p +
    geom_point(aes(x = e1, y = ypos, colour = "e1"), shape = 21,
               fill = ifelse(d$sig1_plot, col1, "white"), stroke = 0.35, size = 1.25) +
    geom_point(aes(x = e2, y = ypos, colour = "e2"), shape = 21,
               fill = ifelse(d$sig2_plot, col2, "white"), stroke = 0.35, size = 1.25) +
    scale_colour_manual(values = c(e1 = col1, e2 = col2),
                        labels = c(e1 = est1_lab, e2 = est2_lab), name = NULL) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    labs(x = xlab) +
    theme_masld_compact() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(), legend.position = "top",
          legend.text = element_text(size = 6), legend.key.size = unit(0.3, "lines"),
          legend.margin = margin(0, 0, 0, 0), plot.margin = margin(18, 2, 2, 2))
}

lxmin <- min(cand$bulk_logFC, cand$protein_logFC, na.rm = TRUE)   # left edge for the inline legend
lxmax <- max(cand$bulk_logFC, cand$protein_logFC, na.rm = TRUE)
sig_key_x <- lxmax - 1.85   # 1.65 -> 1.85: keeps the key box off the panel edge at the 5.5in canvas
lp <- dot_panel(cand, "gene", ord, est1 = "bulk_logFC", est2 = "protein_logFC",
                est1_lab = "mRNA log2FC", est2_lab = "protein log2FC",
                xlab = "log2 FC", group_of = group_of,   # 2026-07-09: plotmath expression() was the DejaVuSans leak source (plain string is font-safe)
                sig1 = "bulk_sig", sig2 = "protein_sig") +
  scale_y_continuous(limits = c(0.5, ytop), expand = c(0, 0)) +
  guides(colour = "none") +                                       # legend drawn inside the panel instead
  annotate("rect", xmin = lxmin - 0.18, xmax = lxmin + 2.18, ymin = n + 1.65, ymax = n + 3.98,
           fill = "white", colour = "grey70", linewidth = 0.2) +
  # 2026-07-09: reverted to original, matching lollipop_panel()'s reverted col1/col2.
  annotate("point", x = lxmin + 0.05, y = c(n + 3.2, n + 2.2), shape = 21,
           colour = c("#8EC7E2", "#E8853A"), fill = c("#8EC7E2", "#E8853A"), stroke = 0.35, size = 1.2) +
  annotate("text", x = lxmin + 0.22, y = n + 3.2, label = "mRNA log2FC", hjust = 0, size = 6/.pt, family = FAM, colour = house_ink) +
  annotate("text", x = lxmin + 0.22, y = n + 2.2, label = "protein log2FC", hjust = 0, size = 6/.pt, family = FAM, colour = house_ink) +
  annotate("rect", xmin = sig_key_x - 0.16, xmax = sig_key_x + 1.82, ymin = 0.95, ymax = 2.95,
           fill = "white", colour = "grey70", linewidth = 0.2) +
  annotate("point", x = sig_key_x + 0.05, y = c(2.35, 1.55), shape = 21,
           colour = "grey40", fill = c("grey40", "white"), stroke = 0.35, size = 1.2) +
  annotate("text", x = sig_key_x + 0.22, y = 2.35, label = "padj < 0.05",
           hjust = 0, size = 6/.pt, family = FAM) +
  annotate("text", x = sig_key_x + 0.22, y = 1.55, label = "padj >= 0.05",
           hjust = 0, size = 6/.pt, family = FAM) +
  ggtitle("mRNA vs protein") +
  theme(text = element_text(family = FAM, face = "plain"),
        plot.title = element_text(size = 6, family = FAM, face = "plain", hjust = 0, margin = margin(0, 0, 2, 0)),
        axis.title.x = element_text(size = 6, family = FAM, face = "plain", margin = margin(t = 1)),
        axis.text.x = element_text(size = 6, family = FAM, margin = margin(t = 7)),
        legend.position = "none", plot.margin = margin(1, 4, 1, 0))

# PANEL 4 (former far-right fgsea NES bar column) REMOVED 2026-07-09 — pathway-level NES now rides on the
# LEFT program brackets as up/down arrow glyphs (see `chev` build above; enr read moved up to feed `bands`).
# Mito/OXPHOS stays documented in the caption. hp absorbs the freed 0.65 width unit (2.65 -> 3.30).
# 2026-08-08: canvas narrowed 6.63 -> 5.5 in (house maximum) at unchanged height; RY/ASP below retuned
# so the pies stay circular at the new cell aspect.

# NOT collected: each colorbar stays under its OWN panel (z-score under heatmap, histology-r under pie),
# so the two similar blue-red scales are spatially tied to their data and separated (mRNA key is inline).
comp <- hp + pp + lp + plot_layout(widths = c(3.30, 0.62, 1.02))
out <- file.path(FIG4_DIR, "panels", "fig4c_mrna_protein_composite.pdf")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
ggsave(out, comp, width = 5.5, height = 2.83, device = grDevices::cairo_pdf)
cat("[fig4c] saved:", out, "\n[fig4c]", n, "proteins x", np, "patients;", length(gl), "gene-row programs (NES on program brackets)\n")
message("CAPTION (Fig 4c full-row proteomics composite): the shared 25-protein axis (five selected MASLD ",
        "process sets) was originally chosen in PXD051911 using protein significance and mRNA/protein direction ",
        "concordance, then fixed before this batch-aware re-estimation in the same liver DIA-MS cohort ",
        "(n=58). Protein-group intensities were log2 transformed and quantile normalized. MASLD-vs-control protein ",
        "effects adjust for acquisition batch, age, BMI, and sex. LEFT = nuisance-adjusted abundance z-score across ",
        "patients ordered by fibrosis then NAS, with NAS/fibrosis strips and gray denoting missing abundance. Acquisition ",
        "batch is retained in the statistical model and diagnostic sidecars but is not displayed as a biological annotation. ",
        "MIDDLE = partial Spearman associations ",
        "with histology after residualizing the same nuisance variables (BH across 25 x 5 tests). RIGHT = canonical ",
        "bulk mRNA log2FC (TREAT significance) versus adjusted protein log2FC. Program brackets show fgsea NES after BH correction across ",
        "the full Hallmark, Reactome, and curated-set family. These process sets are not Hotspot modules, and their ",
        "selection-conditioned row statistics are descriptive rather than independent validation.")

fwrite(cand[, .(gene, bulk_logFC = round(bulk_logFC, 3), bulk_padj = signif(bulk_padj, 4),
                 bulk_sig, protein_logFC = round(protein_logFC, 3),
                 protein_padj = signif(protein_padj, 4), protein_sig, both_sig,
                 selection_conditioned, interpretation,
                 module = group_of[gene])], file.path(DATA_DIR, "composite_mrna_protein_corrected.csv"))
fwrite(as.data.table(cmat, keep.rownames = "gene"), file.path(DATA_DIR, "composite_protein_histology_partial_corrected.csv"))
