# KEY MESSAGE: what a model could have reached was measured before what it did reach was
# judged; the distance between the two marks on each row is the whole argument.

source(file.path(Sys.getenv("FIG7_SCRIPTS"), "_fig7_common.R"))
suppressPackageStartupMessages(library(patchwork))

b <- tab("panel_7d_bounds.tsv")
b$value <- num(b$value)
stopifnot(!any(is.na(b$value)))

BA <- unique(b$block[startsWith(b$block, "a.")]); stopifnot(length(BA) == 1L)
BB <- unique(b$block[startsWith(b$block, "b.")]); stopifnot(length(BB) == 1L)
BC <- unique(b$block[startsWith(b$block, "c.")]); stopifnot(length(BC) == 1L)
BD <- unique(b$block[startsWith(b$block, "d.")]); stopifnot(length(BD) == 1L)

FA <- "signed eQTL, AUC"
FB <- "severity transfer, Spearman"
FC <- "LD-block join, Pearson r"
FD <- "routing, macro Spearman gain"

get1 <- function(blk, row, quantity) {
  v <- b$value[b$block == blk & b$row == row & b$quantity == quantity]
  stopifnot(length(v) == 1L)
  v
}

# ---------------------------------------------------------------- a. signed eQTL AUC
STRATA <- unique(b$row[b$block == BA])
SLAB   <- c("all held-out rows" = "all rows")
slab   <- function(s) ifelse(s %in% names(SLAB), unname(SLAB[s]), s)
aa <- do.call(rbind, lapply(STRATA, function(s) {
  data.frame(facet = FA, row = slab(s),
             x    = c(get1(BA, s, "achieved AUC (best embedding head)"),
                      get1(BA, s, "sequence-attributable bound"),
                      get1(BA, s, "measured chance")),
             mark = c("achieved", "measured bound", "reference"),
             lo   = c(get1(BA, s, "achieved AUC 95% CI low"), NA, NA),
             hi   = c(get1(BA, s, "achieved AUC 95% CI high"), NA, NA),
             stringsAsFactors = FALSE)
}))

# ------------- b. released severity instrument in its three transfer cohorts (review R08-12)
# Corrected 2026-10-15: the block now draws the cohorts the legend's 0.68-0.79 refers to,
# achieved Spearman (filled, 95% CI) beside each cohort's own tie-structure ceiling (hollow).
BLAB <- c("GSE268273 fibrosis 0-4, n=109 participant records"   = "GSE268273, 109",
          "GSE213621 3-level stage, n=299 participant records"  = "GSE213621, 299",
          "GSE276114 3-level stage, n=177 participant records"  = "GSE276114, 177")
QACH <- "achieved Spearman (released v1.1/v1.2 weights, no refit)"
bb <- do.call(rbind, lapply(names(BLAB), function(s) {
  data.frame(facet = FB, row = unname(BLAB[s]),
             x    = c(get1(BB, s, QACH), get1(BB, s, "tie-structure ceiling")),
             mark = c("achieved", "measured bound"),
             lo   = c(get1(BB, s, "achieved Spearman 95% CI low"), NA),
             hi   = c(get1(BB, s, "achieved Spearman 95% CI high"), NA),
             stringsAsFactors = FALSE)
}))
stopifnot(nrow(bb) == 2 * length(BLAB), length(unique(b$row[b$block == BB])) == length(BLAB))
bb0 <- data.frame(row = names(BLAB),
                  ach = vapply(names(BLAB), function(s) get1(BB, s, QACH), 0),
                  ceil = vapply(names(BLAB), function(s) get1(BB, s, "tie-structure ceiling"), 0))

# --------------------------------------------------------- c. LD-block join structure
CLAB <- c("per-gene mean beta, the row's own block HELD OUT" = "own block held out",
          "per-gene mean beta, the row's own block LEFT IN"  = "own block left in")
cc0 <- b[b$block == BC & b$quantity == "Pearson r", ]
stopifnot(nrow(cc0) == length(CLAB))
cc <- data.frame(facet = FC, row = short(cc0$row, CLAB), x = cc0$value,
                 mark = ifelse(cc0$kind == "honest", "achieved", "measured bound"),
                 lo = NA_real_, hi = NA_real_, stringsAsFactors = FALSE)

# ------------------------------------------------------------------- d. routing oracle
DLAB <- c("the routing rule actually fitted"             = "rule as fitted",
          "ORACLE router with perfect post-hoc knowledge" = "oracle router",
          "pre-declared worthwhile-improvement floor"     = "pre-declared floor")
dd0 <- b[b$block == BD, ]
stopifnot(nrow(dd0) == length(DLAB))
DMARK <- c("achieved" = "achieved", "ceiling" = "measured bound", "threshold" = "reference")
dd <- data.frame(facet = FD, row = short(dd0$row, DLAB), x = dd0$value,
                 mark = short(dd0$kind, DMARK), lo = NA_real_, hi = NA_real_,
                 stringsAsFactors = FALSE)

d <- rbind(aa, bb, cc, dd)
d$facet <- factor(d$facet, levels = c(FA, FB, FC, FD))
d$mark  <- factor(d$mark, levels = c("achieved", "measured bound", "reference"))

SHP  <- c("achieved" = 21, "measured bound" = 21, "reference" = 22)
FILC <- c("achieved" = COL_UP, "measured bound" = "white", "reference" = COL_NULL)
LINC <- c("achieved" = COL_UP, "measured bound" = COL_NULL, "reference" = COL_NULL)

# Each block is measured in its own metric, so each block gets its own x scale.
# facet_grid can free a scale across columns only, so the blocks are separate plots.
block_plot <- function(dd, keys) {
  lv <- unique(dd$row)
  dd$key <- factor(dd$row, levels = rev(lv))
  ggplot(dd, aes(x = x, y = key)) +
    geom_errorbarh(data = dd[!is.na(dd$lo), ], aes(xmin = lo, xmax = hi),
                   height = 0, linewidth = 0.4, colour = COL_UP) +
    geom_point(aes(shape = mark, fill = mark, colour = mark), size = 1.2, stroke = 0.4) +
    scale_shape_manual(values = SHP, drop = FALSE) +
    scale_fill_manual(values = FILC, drop = FALSE, guide = "none") +
    scale_colour_manual(values = LINC, drop = FALSE, guide = "none") +
    scale_x_continuous(n.breaks = 4) +
    facet_grid(rows = vars(facet)) +
    labs(x = NULL, y = NULL) +
    guides(shape = guide_legend(nrow = 1,
             override.aes = list(fill = unname(FILC), colour = unname(LINC)))) +
    theme_fig7() +
    theme(strip.text.y = element_text(size = FS, colour = "black", face = "plain",
                                      angle = 0, hjust = 0, margin = margin(1, 1, 1, 4)),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = if (keys) "bottom" else "none",
          plot.margin = margin(1, 3, 1, 3))
}

LV    <- levels(d$facet)
parts <- lapply(seq_along(LV),
                function(i) block_plot(d[d$facet == LV[i], ], keys = i == length(LV)))
hts   <- vapply(LV, function(f) length(unique(d$row[d$facet == f])) + 1.7, numeric(1))
hts[length(hts)] <- hts[length(hts)] + 1.2
p <- patchwork::wrap_plots(parts, ncol = 1, heights = hts)

save_panel(p, "fig7D_measured_bounds_and_join_structure.pdf", width = 3.6, height = 3.15)

caption("Fig. 7D. Measured bounds beside the values actually achieved. Filled magenta = ",
        "achieved (whisker = 95% CI where the producing file records one); hollow = a bound ",
        "measured from the data itself; gray tick = a reference fixed before the answer. Each ",
        "block carries its own metric and therefore its own x scale.")
caption("  ", FA, " - ", b$note[b$block == BA][1],
        ". The bound drawn is the sequence-attributable bound: at most 1.89% of held-out ",
        "rows can be causal, so no sequence model can exceed it. The producing file also ",
        "records wider sign-reliability ceilings for the same strata (empirical Bayes ",
        paste(sprintf("%.4f", b$value[b$block == BA & b$quantity == "sign-reliability ceiling, empirical Bayes"]), collapse = ", "),
        "; flat prior ",
        paste(sprintf("%.4f", b$value[b$block == BA & b$quantity == "sign-reliability ceiling, flat prior"]), collapse = ", "),
        "; full oracle ",
        paste(sprintf("%.4f", b$value[b$block == BA & b$quantity == "sign-reliability ceiling, full oracle"]), collapse = ", "),
        ", in row order ", paste(slab(STRATA), collapse = ", "), "). Stratum sizes: ",
        paste(unique(paste0(slab(b$row[b$block == BA]), " ", b$note[b$block == BA])), collapse = "; "),
        ". Source ", unique(b$source[b$block == BA]), ".")
caption("  ", FB, " - the released severity instrument (no refit) in the three cohorts it ",
        "never trained on: filled = achieved Spearman with its 95% CI, hollow = the largest ",
        "Spearman an untied predictor could reach given that cohort's own stage coding. ",
        paste(sprintf("%s Spearman %.4f, ceiling %.4f, ratio %.4f", short(bb0$row, BLAB), bb0$ach,
                      bb0$ceil, bb0$ach / bb0$ceil), collapse = "; "),
        ". Source ", unique(b$source[b$block == BB]), ". Ceilings for other cohorts (GSE267145, ",
        "GSE130970, GSE193066, GSE83452) bound other questions and are listed in the supplement.")
caption("  ", FC, " - what an honest split removes. ",
        paste(sprintf("%s %.4f (%s)", short(cc0$row, CLAB), cc0$value, cc0$note), collapse = "; "),
        ". Two further guards in the same producing files are not drawn because they are ",
        "fractions, not correlations: held-out rows whose gene has any training data ",
        sprintf("%.5f", get1(BC, "held-out rows whose gene has any training data", "fraction")),
        " (", b$note[b$block == BC & b$quantity == "fraction"], "); cross-fold r2 violations ",
        "at the 0.01 threshold ",
        sprintf("%.4f", get1(BC, "cross-fold r2 violations at the 0.01 threshold", "fraction of examined pairs")),
        " (", b$note[b$block == BC & b$quantity == "fraction of examined pairs"], "). Sources ",
        paste(unique(b$source[b$block == BC]), collapse = "; "), ".")
caption("  ", FD, " - ",
        paste(sprintf("%s %.5f (%s)", short(dd0$row, DLAB), dd0$value, dd0$note), collapse = "; "),
        ". Source ", unique(dd0$source), ".")
