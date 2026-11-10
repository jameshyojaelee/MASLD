# Critical-fix candidate 2026-10-26 (FIG-SCIREVIEW F7-03, V7-01): treatment response and routing move from
# tested negative to untestable in the source table; the unused tested-negative key is removed. Values unchanged.
# Candidate redraw 2026-10-17 (panel-deep-fixes), edited copy of
# Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7a/scripts/fig7e_verdict_vocabulary_and_refusal_thresholds.R
#
# Changes against the installed panel (REVIEW.md, 7E):
#   * one x scale for every row: (estimate - null) / (threshold - null), so 0 is the null and
#     1 is the threshold the verdict uses; rows with a recorded interval draw it on that scale;
#   * all seven severity arms from ALL_ARMS_SUMMARY.json, each labelled precommitted or post hoc
#     from the producing JSON's own `precommitted` field, with n per arm;
#   * signed eQTL is drawn as untestable, under the panel's own rule: the sequence-attributable
#     bound (open diamond) lies below the threshold, so no sequence model could reach it;
#   * routing draws the fitted rule and the oracle bound; treatment response its interval;
#   * the program-h2 AST/ALT rows (0 of 117 programs, a count, not on this scale) move to the
#     legend text; the zone row uses the MDE threshold rather than the significance bar.
# Verdicts are copied from the producing files, never recomputed here.
# All drawn values come from source_tables/fig7E.tsv.

source(file.path(Sys.getenv("FIG7_CODE"), "_fig7_common.R"))
suppressPackageStartupMessages(library(patchwork))

d <- src_tab("7E")
KEYLAB <- c("supported" = "supported", "untestable" = "untestable", "tested negative" = "tested negative",
            "sequence-attributable bound" = "bound measured from the data",
            "oracle router bound" = "bound measured from the data")
d$key <- short(d$cls, KEYLAB)
KEYS <- data.frame(key = c("supported", "untestable", "bound measured from the data"),
                   glyph = c("circle_filled", "circle_open", "diamond_open"),
                   colour = c(COL_UP, COL_NULL, "#000000"))
stopifnot(all(paste(d$key, d$glyph, d$colour) %in% paste(KEYS$key, KEYS$glyph, KEYS$colour)))
d$strip_lab <- d$block
VL <- list(list(x = 0, lty = "solid", col = COL_NULL), list(x = 1, lty = "dashed", col = "black"))
blocks <- unique(d$block[order(d$block_order)])
parts <- lapply(seq_along(blocks), function(i) {
  b <- blocks[i]
  forest_block(d[d$block == b, ],
               if (i == length(blocks)) "(estimate - null) / (threshold - null): 0 = null, 1 = threshold" else NULL,
               vlines = VL, legend = i == length(blocks), xlim = c(-1.15, 2.7),
               breaks = c(-1, 0, 1, 2), keys = KEYS)
})
hts <- vapply(blocks, function(b) length(unique(d$row_label[d$block == b])) + 1.6, 0)
hts[length(hts)] <- hts[length(hts)] + 2.2
p <- wrap_plots(parts, ncol = 1, heights = hts)
save_panel(p, "fig7E_verdict_vocabulary_and_refusal_thresholds.pdf", width = 4.6, height = 4.4)

for (b in blocks) {
  x <- d[d$block == b, ]
  caption("Fig. 7E ", b, ": ",
          paste(sprintf("%s [%s] %.4f%s (%s)", x$row_label, x$mark, x$x,
                        ifelse(is.na(x$lo), "", sprintf(" [%.4f, %.4f]", x$lo, x$hi)), x$unit), collapse = "; "), ".")
}
