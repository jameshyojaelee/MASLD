#!/usr/bin/env Rscript
# composite_helpers.R — shared grammar for the "Abdellaoui-style" composite figures:
#   LEFT  = grouped pie-glyph correlation matrix (pac-man arcs, red +/blue -,
#           thematic blocks with dashed separators)
#   RIGHT = paired-lollipop companion (two estimates of the SAME quantity per row;
#           the gap = "total vs captured"). This is the ONE sanctioned exception to
#           the no-lollipop rule.
# Reused by every composite panel so they share an identical look. Source AFTER
# publication_theme.R. PDF only, 6 pt.

suppressPackageStartupMessages({
  library(ggplot2); library(ggforce); library(data.table); library(patchwork)
})

# 2026-07-09: reverted to the original moderate red/blue -- the exact RdBu poles
# read as too dark/heavy in review ("doesn't look like a regular heatmap").
PIE_POS <- "#B2182B"   # red   = positive correlation
PIE_NEG <- "#2166AC"   # blue  = negative correlation

# ── grouped pie-glyph correlation matrix ─────────────────────────────────────
# cormat: symmetric matrix (rownames = items). ord: item order (top→bottom / left→right).
# group_of: named vector item -> group label. italic_items: TRUE if items are gene/module symbols.
pie_glyph_matrix <- function(cormat, ord, group_of, group_levels = NULL,
                             rad = 0.46, italic_items = TRUE, block_tint = "#F2F5F9",
                             corr_title = "Correlation") {
  n <- length(ord)
  stopifnot(all(ord %in% rownames(cormat)))
  g <- unname(group_of[ord])
  if (is.null(group_levels)) group_levels <- unique(g)

  # lower-triangle cells (row i at top → yi = n-i+1)
  cells <- rbindlist(lapply(seq_len(n), function(i) {
    if (i == 1) return(NULL)
    j <- seq_len(i - 1)
    data.table(xi = j, yi = n - i + 1, corr = cormat[ord[i], ord[j]])
  }))
  cells <- cells[!is.na(corr)]
  cells[, corr := pmax(pmin(corr, 1), -1)]
  cells[, `:=`(r0 = 0, rr = rad)]                    # radius as mapped aesthetics (ggforce needs them)

  # group boundaries (after ordered position k) → dashed separators
  bnd <- which(g[-1] != g[-n])
  # on-diagonal block rectangles (item index range per group)
  blocks <- rbindlist(lapply(group_levels, function(gl) {
    idx <- which(g == gl); if (!length(idx)) return(NULL)
    a <- min(idx); b <- max(idx)
    data.table(group = gl, xmin = a - 0.5, xmax = b + 0.5,
               ymin = n - b + 0.5, ymax = n - a + 0.5)
  }))
  # 2026-07-09: reverted to the original custom palette (Okabe-Ito swap read as
  # part of the "hideous" combination in review).
  pal_solid <- c("#3B6EA5", "#B05A5A", "#4E9E76", "#C79A3E", "#8A6BA8", "#5AA6A6", "#A0743A")
  gcol <- setNames(rep(pal_solid, length.out = length(group_levels)), group_levels)
  blocks[, cx := (xmin + xmax) / 2][, cyd := n - cx + 1]
  blocks[, tint := paste0(gcol[group], "24")]
  blocks[, group := factor(group, levels = group_levels)]

  face <- if (italic_items) "italic" else "plain"
  p <- ggplot() +
    # block tint on the diagonal (per group)
    geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = blocks$tint, colour = NA) +
    # invisible layer to expose a group colour legend (fill scale is used by the pies)
    geom_point(data = blocks, aes(x = cx, y = cyd, colour = group), size = 0.01, alpha = 0) +
    # empty reference circle then filled pac-man arc
    ggforce::geom_circle(data = cells, aes(x0 = xi, y0 = yi, r = rr),
                         colour = "grey85", fill = "white", linewidth = 0.12, inherit.aes = FALSE) +
    ggforce::geom_arc_bar(data = cells[abs(corr) > 1e-6],
                          aes(x0 = xi, y0 = yi, r0 = r0, r = rr, start = 0,
                              end = abs(corr) * 2 * pi, fill = corr),
                          colour = NA) +
    scale_fill_gradient2(low = PIE_NEG, mid = "white", high = PIE_POS, midpoint = 0,
                         limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                         name = corr_title,
                         guide = guide_colourbar(order = 1, barwidth = 0.3, barheight = 2.2)) +
    scale_colour_manual(values = gcol, breaks = group_levels, name = NULL,
                        guide = guide_legend(order = 2, override.aes = list(size = 1.4, alpha = 1))) +
    # dashed group separators
    geom_vline(xintercept = bnd + 0.5, linetype = "dashed", linewidth = 0.25, colour = "grey45") +
    geom_hline(yintercept = n - bnd + 0.5, linetype = "dashed", linewidth = 0.25, colour = "grey45") +
    # row labels (left) + column labels (top, angled)
    # 2026-07-09: family="Helvetica" added explicitly -- annotate()/geom_text() do
    # NOT inherit family from theme(text=...); the missing arg here was silently
    # falling back to DejaVuSans (the font bug flagged for fig4c, root-caused here).
    annotate("text", x = 0.35, y = n:1, label = ord, hjust = 1, size = 6/ggplot2::.pt,
             fontface = face, family = "Helvetica", colour = house_ink) +
    annotate("text", x = 1:n, y = n + 0.7, label = ord, angle = 55, hjust = 0,
             size = 6/ggplot2::.pt, fontface = face, family = "Helvetica", colour = house_ink) +
    coord_fixed(xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5), clip = "off") +
    theme_void() +
    theme(legend.position = "left",
          legend.title = element_text(size = 6), legend.text = element_text(size = 6),
          legend.key.size = unit(0.28, "lines"), legend.spacing = unit(2, "pt"),
          plot.margin = margin(18, 2, 2, 26))
  list(plot = p, ord = ord, n = n, group_of = group_of, blocks = blocks)
}

# ── rectangular gene x feature pie-glyph matrix (e.g. protein-vs-histology corr) ──
# cmat: gene (rows) x feature (cols) numeric matrix in [-1,1]. gene_ord: rows top->bottom
# (must match the lollipop's ord). feature_levels: column order left->right. group_of / group_levels:
# gene -> program row blocks. Same pac-man pie aesthetic as pie_glyph_matrix, but rectangular.
pie_feature_matrix <- function(cmat, gene_ord, feature_levels, group_of, group_levels,
                               rad = 0.46, corr_title = "Correlation", show_row_labels = TRUE,
                               collabel_y = NULL) {
  n <- length(gene_ord); nf <- length(feature_levels)
  g <- unname(group_of[gene_ord])
  cells <- rbindlist(lapply(seq_len(n), function(i) {
    data.table(yi = n - i + 1, xi = seq_len(nf), corr = as.numeric(cmat[gene_ord[i], feature_levels]))
  }))
  cells <- cells[is.finite(corr)]
  cells[, corr := pmax(pmin(corr, 1), -1)][, `:=`(r0 = 0, rr = rad)]

  bnd <- which(g[-1] != g[-n])                                    # program boundaries
  blocks <- rbindlist(lapply(group_levels, function(gg) {
    idx <- which(g == gg); if (!length(idx)) return(NULL)
    data.table(group = gg, ymin = n - max(idx) + 0.5, ymax = n - min(idx) + 0.5)
  }))
  pal_solid <- c("#3B6EA5", "#B05A5A", "#4E9E76", "#C79A3E", "#8A6BA8", "#5AA6A6", "#A0743A")
  gcol <- setNames(rep(pal_solid, length.out = length(group_levels)), group_levels)
  blocks[, tint := paste0(gcol[group], "24")][, group := factor(group, levels = group_levels)]
  blocks[, cy := (ymin + ymax) / 2]

  ggplot() +
    geom_rect(data = blocks, aes(xmin = 0.5, xmax = nf + 0.5, ymin = ymin, ymax = ymax),
              fill = blocks$tint, colour = NA) +
    geom_point(data = blocks, aes(x = 1, y = cy, colour = group), size = 0.01, alpha = 0) +
    ggforce::geom_circle(data = cells, aes(x0 = xi, y0 = yi, r = rr),
                         colour = "grey85", fill = "white", linewidth = 0.12, inherit.aes = FALSE) +
    ggforce::geom_arc_bar(data = cells[abs(corr) > 1e-6],
                          aes(x0 = xi, y0 = yi, r0 = r0, r = rr, start = 0,
                              end = abs(corr) * 2 * pi, fill = corr), colour = NA) +
    scale_fill_gradient2(low = PIE_NEG, mid = "white", high = PIE_POS, midpoint = 0,
                         limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), name = corr_title,
                         guide = guide_colourbar(order = 1, barwidth = 0.3, barheight = 2.2)) +
    scale_colour_manual(values = gcol, breaks = group_levels, name = NULL,
                        guide = guide_legend(order = 2, override.aes = list(size = 1.4, alpha = 1))) +
    geom_hline(yintercept = n - bnd + 0.5, linetype = "dashed", linewidth = 0.25, colour = "grey45") +
    { if (show_row_labels) annotate("text", x = 0.35, y = n:1, label = gene_ord, hjust = 1,
                                    size = 6/ggplot2::.pt, fontface = "italic",
                                    family = "Helvetica", colour = house_ink) } +
    annotate("text", x = seq_len(nf), y = if (is.null(collabel_y)) n + 0.7 else collabel_y,
             label = feature_levels, angle = 45, hjust = 0, size = 6/ggplot2::.pt,
             family = "Helvetica", colour = house_ink) +
    coord_fixed(xlim = c(0.5, nf + 0.5), ylim = c(0.5, n + 0.5), clip = "off") +
    theme_void() +
    theme(legend.position = "left",
          legend.title = element_text(size = 6), legend.text = element_text(size = 6),
          legend.key.size = unit(0.28, "lines"), legend.spacing = unit(2, "pt"),
          plot.margin = margin(18, 2, 2, if (show_row_labels) 26 else 4))
}

# ── paired-lollipop companion ────────────────────────────────────────────────
# dt: per-item table with item column + two estimate columns. ord must match the matrix.
lollipop_panel <- function(dt, item_col, ord, est1, est2, est1_lab, est2_lab, xlab,
                           col1 = "#8EC7E2", col2 = "#E8853A", group_of = NULL,   # 2026-07-09: reverted to original
                           group_levels = NULL, ref0 = TRUE, sig1 = NULL, sig2 = NULL) {
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
  # group separators to mirror the matrix
  if (!is.null(group_of)) {
    g <- unname(group_of[ord]); bnd <- which(g[-1] != g[-n])
    p <- p + geom_hline(yintercept = n - bnd + 0.5, linetype = "dashed",
                        linewidth = 0.25, colour = "grey45")
  }
  p <- p +
    geom_segment(aes(x = 0, xend = e1, y = ypos, yend = ypos), colour = col1, linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = e2, y = ypos, yend = ypos), colour = col2, linewidth = 0.3) +
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
  p
}

# ── assemble matrix | lollipop, row-aligned ──────────────────────────────────
assemble_composite <- function(mat, loll, widths = c(2.4, 1), out_pdf, height = 4.4,
                               width = 7.2, caption = NULL) {
  comp <- mat$plot + loll + patchwork::plot_layout(widths = widths)
  ggplot2::ggsave(out_pdf, comp, width = width, height = height, device = grDevices::cairo_pdf)
  if (!is.null(caption)) message(caption)
  cat("[composite] Saved:", out_pdf, "\n")
  invisible(comp)
}

# ── bulk co-expression matrix (symbol x symbol) from the C2 corrected logCPM ─
# Maps gene symbols -> versioned ENSG via canonical_deg_results, returns Spearman
# co-expression among the symbols that map. base = MASLD_PROJECT_ROOT.
bulk_coexpr <- function(symbols, base) {
  m <- readRDS(file.path(base,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/corrected_logcpm_c2_nogerhard.rds"))
  rownames(m) <- sub("\\..*", "", rownames(m))
  map <- fread(file.path(base,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
    select = c("gene", "symbol"))
  map[, ensg := sub("\\..*", "", gene)]
  map <- map[symbol %in% symbols & ensg %in% rownames(m)][!duplicated(symbol)]
  sm <- m[map$ensg, , drop = FALSE]; rownames(sm) <- map$symbol
  cm <- cor(t(sm), method = "spearman"); cm[is.na(cm)] <- 0
  list(cormat = cm, present = map$symbol)
}

# ── disease-residualized protein co-abundance (partial correlation) ──────────
# abmat: gene x sample numeric matrix (rownames = symbols, colnames = sample IDs).
# genes: subset to correlate. meta: data.table with a `raw` column matching
# colnames(abmat) plus the covariates referenced by `formula` (default
# disease_group + age + bmi + sex). Each protein's log2 abundance is regressed on
# those covariates (na.exclude keeps sample alignment) and the Spearman correlation
# is taken among the RESIDUALS — so the shared disease axis is removed and the
# remaining structure reflects genuine co-regulation, not "everything moves in disease".
partial_coexpr <- function(abmat, genes, meta,
                           formula = y ~ disease_group + age + bmi + sex, min_n = 6) {
  g <- genes[genes %in% rownames(abmat)]
  sub <- abmat[g, , drop = FALSE]
  m <- meta[match(colnames(sub), meta$raw)]            # align covariates to abmat columns
  design_df <- as.data.frame(m)[, setdiff(all.vars(formula), "y"), drop = FALSE]
  logmat <- log2(sub); logmat[!is.finite(logmat)] <- NA
  resid <- t(apply(logmat, 1, function(y) {
    if (sum(is.finite(y)) < min_n) return(rep(NA_real_, length(y)))
    df  <- data.frame(y = as.numeric(y), design_df)
    fit <- tryCatch(lm(formula, data = df, na.action = na.exclude), error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, length(y)))
    as.numeric(residuals(fit))                          # na.exclude pads NA → length preserved
  }))
  rownames(resid) <- g; colnames(resid) <- colnames(sub)
  cm <- cor(t(resid), use = "pairwise.complete.obs", method = "spearman")
  cm[is.na(cm)] <- 0
  cm
}

# ── diagnostic: mean |corr| within vs between pathway blocks (circularity check) ─
block_contrast <- function(cormat, group_of) {
  it <- rownames(cormat); gv <- unname(group_of[it])
  lt <- which(lower.tri(cormat), arr.ind = TRUE)
  same <- gv[lt[, 1]] == gv[lt[, 2]]
  v <- abs(cormat[lower.tri(cormat)])
  c(within = mean(v[same], na.rm = TRUE), between = mean(v[!same], na.rm = TRUE))
}

# ── helper: order items within groups by hierarchical clustering of |corr| ───
order_by_group_then_clust <- function(cormat, group_of, group_levels) {
  unlist(lapply(group_levels, function(gl) {
    it <- names(group_of)[group_of == gl]; it <- it[it %in% rownames(cormat)]
    if (length(it) <= 2) return(it)
    cm <- cormat[it, it]; cm[is.na(cm)] <- 0
    d <- as.dist(1 - cm); it[hclust(d, "average")$order]
  }), use.names = FALSE)
}
