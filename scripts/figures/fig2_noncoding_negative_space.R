#!/usr/bin/env Rscript
# KEY MESSAGE: Most fine-mapped noncoding posterior mass lacks interpretable
# deposited liver context; unresolved loci require element- and target-resolving
# experiments rather than nearest-gene assignment.

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
INPUT <- Sys.getenv(
  "FIG2_NONCODING_CONTEXT_INPUT",
  file.path(BASE, "figures/candidates/noncoding-resource-2026-08-12/noncoding-dna-v5",
            "fig2_noncoding_regulatory_context_source.tsv")
)
OUT_DIR <- Sys.getenv("FIG2_NEGATIVE_SPACE_DIR", file.path(BASE, "figures/candidates"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

d <- fread(INPUT)
required <- c("trait_scope", "trait_scope_label", "value", "n_credible_sets",
              "total_noncoding_pip_mass", "metric", "context_label")
stopifnot(all(required %in% names(d)))
stopifnot(all(d$value >= 0 & d$value <= 1))

scope_order <- c("tier1_direct_masld_pdff", "tier2_liver_enzyme")
scope_short <- c(
  tier1_direct_masld_pdff = "Direct MASLD / liver fat",
  tier2_liver_enzyme = "Liver enzymes"
)
unresolved <- d[metric == "unresolved_context_mass"]
unresolved[, scope_order_index := match(trait_scope, scope_order)]
setorder(unresolved, scope_order_index)
unresolved[, scope_order_index := NULL]
stopifnot(nrow(unresolved) == 2L)

context_metrics <- c("promoter_mass", "accessible_mass", "abc_mass")
context_labels <- c(
  promoter_mass = "Within 2 kb of TSS",
  accessible_mass = "Open in any liver lineage",
  abc_mass = "Liver ABC enhancer"
)
context <- d[metric %in% context_metrics]
context[, metric := factor(metric, levels = context_metrics)]
context[, trait_scope := factor(trait_scope, levels = scope_order)]
setorder(context, metric, trait_scope)
stopifnot(nrow(context) == 6L)

# Re-derive the exact values printed in the panel from the immutable source.
expected <- c(
  tier1_direct_masld_pdff.promoter_mass = 0.0787249975787247,
  tier2_liver_enzyme.promoter_mass = 0.0445912158230019,
  tier1_direct_masld_pdff.accessible_mass = 0.0731733768551714,
  tier2_liver_enzyme.accessible_mass = 0.0415012498200219,
  tier1_direct_masld_pdff.unresolved_context_mass = 0.521129683271921,
  tier2_liver_enzyme.unresolved_context_mass = 0.604308028252816
)
observed <- setNames(
  d[paste(trait_scope, metric, sep = ".") %in% names(expected), value],
  d[paste(trait_scope, metric, sep = ".") %in% names(expected),
    paste(trait_scope, metric, sep = ".")]
)
stopifnot(isTRUE(all.equal(observed[names(expected)], expected, tolerance = 1e-12)))

PAGE_W <- 5.50
PAGE_H <- 2.85
out_pdf <- file.path(OUT_DIR, "Fig2_companion_noncoding_negative_space.pdf")
cairo_pdf(out_pdf, width = PAGE_W, height = PAGE_H, family = "Helvetica",
          onefile = FALSE)
grid.newpage()
pushViewport(viewport(xscale = c(0, PAGE_W), yscale = c(0, PAGE_H),
                      default.units = "inches"))

u <- function(x) unit(x, "native")
txt <- function(label, x, y, ..., size = 6, just = "center", col = "black") {
  grid.text(label, x = u(x), y = u(y), just = just,
            gp = gpar(fontfamily = "Helvetica", fontsize = size, col = col, ...))
}
rect <- function(x, y, w, h, fill, col = NA, lwd = 0.5) {
  grid.rect(x = u(x), y = u(y), width = u(w), height = u(h),
            just = c("left", "center"), gp = gpar(fill = fill, col = col, lwd = lwd))
}
seg <- function(x0, y0, x1, y1, col = "black", lwd = 0.5, lty = 1) {
  grid.segments(u(x0), u(y0), u(x1), u(y1),
                gp = gpar(col = col, lwd = lwd, lty = lty))
}

COL_UNRES <- "#C2185B"
COL_REST <- "#E6E6E5"
COL_DIRECT <- "#1B4F8A"
COL_ENZYME <- "#8FB9D8"
COL_TEXT_GRAY <- "#4D4D4D"
COL_ROUTE <- "#F2F2F2"

# Left: explicitly foreground unresolved mass.
txt("Unresolved deposited context", 0.25, 2.62, just = "left", size = 6.5,
    fontface = "plain")
bar_x <- 0.25
bar_w <- 2.35
bar_h <- 0.28
bar_y <- c(2.12, 1.58)
for (i in seq_len(nrow(unresolved))) {
  value <- unresolved$value[i]
  scope <- unresolved$trait_scope[i]
  label <- scope_short[[scope]]
  ncs <- unresolved$n_credible_sets[i]
  txt(sprintf("%s  (n = %s credible sets)", label, format(ncs, big.mark = ",")),
      bar_x, bar_y[i] + 0.23, just = "left", size = 6)
  rect(bar_x, bar_y[i], bar_w, bar_h, COL_REST)
  rect(bar_x, bar_y[i], bar_w * value, bar_h, COL_UNRES)
  txt(sprintf("%.2f%% unresolved", 100 * value),
      bar_x + bar_w * value - 0.06, bar_y[i], just = "right",
      size = 6.5, col = "white", fontface = "plain")
  txt("remainder", bar_x + bar_w * (value + (1 - value) / 2), bar_y[i],
      size = 6, col = COL_TEXT_GRAY)
}
seg(bar_x + bar_w * 0.5, 1.35, bar_x + bar_w * 0.5, 2.32,
    col = "white", lwd = 0.7, lty = 2)
txt("0", bar_x, 1.31, just = "center", size = 6, col = COL_TEXT_GRAY)
txt("50", bar_x + bar_w * 0.5, 1.31, just = "center", size = 6, col = COL_TEXT_GRAY)
txt("100%", bar_x + bar_w, 1.31, just = "center", size = 6, col = COL_TEXT_GRAY)

# Right: deposited maps capture a minority. Metrics are nonexclusive.
right_x0 <- 3.48
right_x1 <- 5.27
label_x <- 2.86
txt("Positive deposited context", label_x, 2.62, just = "left", size = 6.5)
axis_y <- 1.28
seg(right_x0, axis_y, right_x1, axis_y, col = "black", lwd = 0.55)
for (tick in c(0, 5, 10)) {
  xx <- right_x0 + (tick / 10) * (right_x1 - right_x0)
  seg(xx, axis_y, xx, axis_y - 0.05, lwd = 0.5)
  txt(sprintf("%d%%", tick), xx, axis_y - 0.13, size = 6, col = COL_TEXT_GRAY)
}
metric_y <- c(promoter_mass = 2.19, accessible_mass = 1.82, abc_mass = 1.45)
for (metric_name in context_metrics) {
  yy <- metric_y[[metric_name]]
  txt(context_labels[[metric_name]], label_x, yy, just = "left", size = 6)
  rows <- context[as.character(metric) == metric_name]
  for (j in seq_len(nrow(rows))) {
    value_pct <- 100 * rows$value[j]
    xx <- right_x0 + (value_pct / 10) * (right_x1 - right_x0)
    dy <- if (as.character(rows$trait_scope[j]) == scope_order[1]) 0.065 else -0.065
    if (as.character(rows$trait_scope[j]) == scope_order[1]) {
      grid.points(u(xx), u(yy + dy), pch = 21, size = unit(0.095, "inches"),
                  gp = gpar(fill = COL_DIRECT, col = "#333333", lwd = 0.5))
    } else {
      grid.points(u(xx), u(yy + dy), pch = 22, size = unit(0.090, "inches"),
                  gp = gpar(fill = COL_ENZYME, col = "#333333", lwd = 0.5))
    }
    txt(sprintf("%.2f", value_pct), xx + 0.06, yy + dy,
        just = "left", size = 6, col = COL_TEXT_GRAY)
  }
}
grid.points(u(3.58), u(2.57), pch = 21, size = unit(0.085, "inches"),
            gp = gpar(fill = COL_DIRECT, col = "#333333", lwd = 0.5))
txt("Direct MASLD / liver fat", 3.68, 2.57, just = "left", size = 6)
grid.points(u(4.66), u(2.57), pch = 22, size = unit(0.082, "inches"),
            gp = gpar(fill = COL_ENZYME, col = "#333333", lwd = 0.5))
txt("Liver enzymes", 4.76, 2.57, just = "left", size = 6)

# Bottom: route unresolved loci to the assays that resolve the missing link.
seg(0.25, 1.08, 5.27, 1.08, col = "#BDBDBD", lwd = 0.6)
txt("Unresolved loci", 0.25, 0.84, just = "left", size = 6.5,
    col = COL_UNRES, fontface = "plain")
txt("route to", 1.18, 0.84, just = "left", size = 6, col = COL_TEXT_GRAY)
route_x <- c(1.72, 2.92, 4.12)
route_w <- 1.05
route_labels <- c("Allele-aware\nelement testing",
                  "Cell-state-specific\nchromatin",
                  "3D target\nlinking")
for (i in seq_along(route_x)) {
  rect(route_x[i], 0.82, route_w, 0.42, COL_ROUTE, col = "#BDBDBD", lwd = 0.5)
  txt(route_labels[i], route_x[i] + route_w / 2, 0.82, size = 6)
}
txt("Do not force nearest-gene labels", 0.25, 0.47, just = "left", size = 6,
    col = COL_TEXT_GRAY)
txt("Context categories overlap; percentages must not be summed.",
    5.27, 0.47, just = "right", size = 6, col = COL_TEXT_GRAY)
txt("Regulatory DNA context and colocalized transcript biotype are separate evidence fields.",
    0.25, 0.20, just = "left", size = 6, col = "black")

popViewport()
dev.off()

source_out <- copy(d)
source_out[, displayed_percent := 100 * value]
source_out[, panel_role := fifelse(
  metric == "unresolved_context_mass", "dominant_unresolved_bar",
  fifelse(metric %in% context_metrics, "positive_context_point", "source_only")
)]
fwrite(source_out, file.path(OUT_DIR, "Fig2_companion_noncoding_negative_space_source.tsv"),
       sep = "\t")

caption <- paste(
  "Fine-mapped noncoding posterior mass remains largely outside interpretable deposited liver context.",
  "Context is unresolved for 52.11% of direct MASLD/liver-fat mass and 60.43% of liver-enzyme mass.",
  "Only 7.87% and 4.46%, respectively, lies within 2 kb of a GENCODE TSS, while 7.32% and 4.15% is open in at least one available liver lineage.",
  "Positive regulatory-context categories overlap and are not additive; credible-set instances are descriptive study-by-locus units.",
  "Unresolved loci are routed to allele-aware element testing, cell-state-specific chromatin, and 3D target-linking rather than forced nearest-gene assignment.",
  "Regulatory-DNA context is kept separate from colocalized-transcript biotype."
)
writeLines(caption, file.path(OUT_DIR, "Fig2_companion_noncoding_negative_space_caption.txt"))
cat("wrote:", out_pdf, "\n")
