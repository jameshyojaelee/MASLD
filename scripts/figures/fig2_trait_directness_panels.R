#!/usr/bin/env Rscript
# Panel A KEY MESSAGE: the closer a liver phenotype sits to the molecular lesion,
# the more of its fine-mapped heritability is protein-coding.
# Panel B KEY MESSAGE (revised 2026-08-17 after family multiplicity correction):
# within the noncoding remainder, direct-trait posterior mass is promoter-proximal
# and better resolved overall. It is NOT that the gap is specific to distal
# enhancer-gene maps.
#
# The earlier message, "essentially absent from liver ABC enhancers", is RETIRED.
# Across the family of six direct-vs-enzyme tests on these groups, ABC depletion
# has permutation p = 0.398 (BH q = 0.398) and generic liver accessibility
# p = 0.185 (q = 0.222); neither survives. What survives is promoter proximity
# (3.92x, q = 0.0042), unresolved context (0.53x, q = 6.0e-4) and resolved-no-
# context (2.06x, q = 0.0085). The ABC row is still drawn, because showing the
# non-surviving comparison beside the surviving ones is the honest display, but
# it may not be described as a depletion.
# Source: GWAS/finemapping/results/fig2_multiplicity/20260817T161128Z/
#
# Both panels use 1-Mb genomic-region groups as the unit, not credible-set
# instances. The grouping prevents repeated recovery of one region from
# dominating the summary; it does not assert statistical independence.

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/fig2_consequence_palette.R"))
STATS_DIR <- Sys.getenv("FIG2_TRAIT_DIRECTNESS_STATS")
OUT_DIR <- Sys.getenv("FIG2_TRAIT_DIRECTNESS_OUT")
stopifnot(nzchar(STATS_DIR), nzchar(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

est <- fread(file.path(STATS_DIR, "pip_architecture_estimates.tsv"))
con <- fread(file.path(STATS_DIR, "trait_class_contrasts.tsv"))
ind <- fread(file.path(STATS_DIR, "portfolio_independence.tsv"))
cnc <- fread(file.path(STATS_DIR, "direct_trait_coding_concentration.tsv"))

DIRECT <- "tier1_direct_masld_pdff"
ENZYME <- "tier2_liver_enzyme"
stopifnot(all(c(DIRECT, ENZYME) %in% est$trait_scope))

COL_DIRECT <- "#C9265E"   # Liang deep magenta: direct disease phenotype
COL_ENZYME <- "#1565C0"   # Figure 2 liver-enzyme blue
COL_RULE <- "#BDBDBD"
FS <- 6

grDevices::pdf.options(useDingbats = FALSE)

get_est <- function(scope, cat) {
  row <- est[trait_scope == scope & category == cat]
  stopifnot(nrow(row) == 1L)
  row
}
get_con <- function(cat) {
  row <- con[category == cat]
  stopifnot(nrow(row) == 1L)
  row
}

open_page <- function(path, w, h) {
  cairo_pdf(path, width = w, height = h, family = "Helvetica", onefile = FALSE)
  grid.newpage()
  pushViewport(viewport(xscale = c(0, w), yscale = c(0, h),
                        default.units = "inches"))
}
u <- function(x) unit(x, "native")
txt <- function(label, x, y, size = FS, just = "center", col = "black",
                fontface = "plain", rot = 0) {
  grid.text(label, x = u(x), y = u(y), just = just, rot = rot,
            gp = gpar(fontfamily = "Helvetica", fontsize = size, col = col,
                      fontface = fontface))
}
bar <- function(x, y, w, h, fill, col = NA, lwd = 0.4) {
  grid.rect(x = u(x), y = u(y), width = u(w), height = u(h),
            just = c("left", "center"),
            gp = gpar(fill = fill, col = col, lwd = lwd))
}
seg <- function(x0, y0, x1, y1, col = "black", lwd = 0.5, lty = 1) {
  grid.segments(u(x0), u(y0), u(x1), u(y1), gp = gpar(col = col, lwd = lwd, lty = lty))
}

fmt_p <- function(p) if (p < 1e-4) "P < 1e-4" else sprintf("P = %.3g", p)

row_scope <- c(DIRECT, ENZYME)
row_col <- c(COL_DIRECT, COL_ENZYME)

ci_row <- function(r, y, col, to_x, tick_h = 0.033, pt = 0.080) {
  seg(to_x(r$locus_boot_lo), y, to_x(r$locus_boot_hi), y, col = col, lwd = 0.9)
  seg(to_x(r$locus_boot_lo), y - tick_h, to_x(r$locus_boot_lo), y + tick_h,
      col = col, lwd = 0.9)
  seg(to_x(r$locus_boot_hi), y - tick_h, to_x(r$locus_boot_hi), y + tick_h,
      col = col, lwd = 0.9)
  grid.points(u(to_x(r$instance_weighted_fraction)), u(y), pch = 3,
              size = unit(0.05, "inches"), gp = gpar(col = "black", lwd = 0.5))
  grid.points(u(to_x(r$locus_weighted_fraction)), u(y), pch = 21,
              size = unit(pt, "inches"),
              gp = gpar(fill = col, col = "black", lwd = 0.4))
}

# ---------------------------------------------------------------- Panel A ----
# Complete mutually exclusive consequence composition. This answers what the
# protein-altering share is relative to, without mixing in lead-SNP counts.
PA_W <- 1.85
PA_H <- 2.10
pa_pdf <- file.path(OUT_DIR, "Fig2D_pip_architecture_by_trait_directness.pdf")
open_page(pa_pdf, PA_W, PA_H)

comp_cats <- c("protein_altering_pip_mass", "canonical_splice_pip_mass",
               "synonymous_or_utr_pip_mass", "other_noncoding_pip_mass")
comp_lab <- c("Protein-altering", "Canonical splice", "Synonymous/UTR",
              "Other noncoding")
# Match the Figure 2C sequence palette wherever the classes correspond.
# Canonical splice uses the same dark blue reserved for a distinct sequence
# annotation class in C; the other three mappings are exact semantic matches.
fig2_sequence_palette <- FIG2_CONSEQUENCE_PALETTE[
  c("protein_altering", "canonical_splice", "synonymous_or_utr", "other_noncoding")]
comp_fill <- unname(fig2_sequence_palette)

plot_x0 <- 0.36
plot_y0 <- 0.43
plot_y1 <- 1.93
plot_h <- plot_y1 - plot_y0
bar_w <- 0.20
bar_x <- c(0.59, 0.96)
for (i in seq_along(row_scope)) {
  values <- vapply(comp_cats, function(category) {
    get_est(row_scope[i], category)$locus_weighted_fraction
  }, numeric(1))
  stopifnot(abs(sum(values) - 1) < 1e-8)
  yy <- plot_y0
  for (j in rev(seq_along(values))) {
    hh <- plot_h * values[j]
    grid.rect(x = u(bar_x[i]), y = u(yy), width = u(bar_w), height = u(hh),
              just = c("center", "bottom"),
              gp = gpar(fill = comp_fill[j], col = "white", lwd = 0.45))
    if (values[j] >= 0.05) {
      txt(sprintf("%.1f", 100 * values[j]), bar_x[i], yy + hh / 2, size = 6,
          col = if (j == 1) "white" else "#333333")
    }
    yy <- yy + hh
  }
}

# Keep the two small enzyme segments readable without filling the bars with text.
enzyme_values <- vapply(comp_cats, function(category) {
  get_est(ENZYME, category)$locus_weighted_fraction
}, numeric(1))
enzyme_mid <- vapply(seq_along(enzyme_values), function(j) {
  below <- if (j < length(enzyme_values)) sum(enzyme_values[(j + 1):length(enzyme_values)]) else 0
  plot_y0 + plot_h * (below + enzyme_values[j] / 2)
}, numeric(1))
for (j in c(1, 3)) {
  label_y <- if (j == 1) 1.84 else 1.70
  seg(bar_x[2] + bar_w / 2, enzyme_mid[j], 1.09, label_y,
      col = "#555555", lwd = 0.45)
  txt(sprintf("%.1f", 100 * enzyme_values[j]), 1.11, label_y,
      just = "left", size = 6, col = "#333333")
}

seg(plot_x0, plot_y0, plot_x0, plot_y1, lwd = 0.5)
for (tick in c(0, 50, 100)) {
  yy <- plot_y0 + plot_h * tick / 100
  seg(plot_x0 - 0.04, yy, plot_x0, yy, lwd = 0.5)
  txt(sprintf("%d", tick), plot_x0 - 0.07, yy, just = "right", size = 6)
}
txt("Fine-mapping probability (%)", 0.075, (plot_y0 + plot_y1) / 2,
    size = 6, rot = 90)
txt("Direct MASLD /\nliver fat", bar_x[1], 0.30, size = 6)
txt("Liver\nenzymes", bar_x[2], 0.30, size = 6)

# One right-hand legend column. These four classes are mutually exclusive and
# exhaustive; moving it beside the bars preserves the plot height at 2.10 in.
legend_x <- rep(1.15, 4)
legend_y <- c(1.55, 1.28, 1.01, 0.74)
for (j in seq_along(comp_cats)) {
  bar(legend_x[j], legend_y[j], 5 / 72, 5 / 72, comp_fill[j])
  txt(comp_lab[j], legend_x[j] + 0.095, legend_y[j], just = "left", size = 6,
      col = "#4D4D4D")
}
popViewport()
dev.off()

# ---------------------------------------------------------------- Panel B ----
# One visual claim: deposited liver maps leave most noncoding mass unresolved.
PB_W <- 5.50
PB_H <- 2.70
pb_pdf <- file.path(OUT_DIR, "Fig2_companion_noncoding_annotation_gap.pdf")
open_page(pb_pdf, PB_W, PB_H)

txt("Unresolved context", 0.25, 2.31, just = "left", size = 6.5)

eb_x <- 0.25
eb_w <- 2.35
eb_y <- c(1.94, 1.35)
for (i in seq_along(row_scope)) {
  ru <- get_est(row_scope[i], "unresolved_context_mass")
  value <- ru$instance_weighted_fraction
  txt(if (row_scope[i] == DIRECT) "Direct MASLD / liver fat" else "Liver enzymes",
      eb_x, eb_y[i] + 0.24, just = "left", size = 6)
  bar(eb_x, eb_y[i], eb_w, 0.27, "#E6E6E5")
  bar(eb_x, eb_y[i], eb_w * value, 0.27, row_col[i])
  txt(sprintf("%.2f%% unresolved", 100 * value),
      eb_x + eb_w * value - 0.06, eb_y[i], just = "right",
      size = 6.5, col = "white")
}
for (tick in c(0, 50, 100)) {
  xx <- eb_x + eb_w * tick / 100
  txt(sprintf("%s", if (tick == 100) "100%" else tick), xx, 0.99,
      size = 6, col = "#4D4D4D")
}

bx0 <- 3.63
bx1 <- 5.25
to_bx <- function(v) bx0 + (v / 0.10) * (bx1 - bx0)
ctx_cats <- c("promoter_mass", "accessible_mass", "abc_mass")
ctx_lab <- c("Within 2 kb of TSS", "Open in liver", "Liver ABC enhancer")
cy <- c(2.02, 1.61, 1.20)
txt("Deposited context", 2.88, 2.31, just = "left", size = 6.5)
for (k in seq_along(ctx_cats)) {
  txt(ctx_lab[k], 2.88, cy[k], just = "left", size = 6)
  for (i in seq_along(row_scope)) {
    r <- get_est(row_scope[i], ctx_cats[k])
    yy <- cy[k] + c(0.065, -0.065)[i]
    value <- r$instance_weighted_fraction
    grid.points(u(to_bx(value)), u(yy), pch = if (i == 1) 21 else 22,
                size = unit(0.085, "inches"),
                gp = gpar(fill = row_col[i], col = "#333333", lwd = 0.5))
    txt(sprintf("%.2f%%", 100 * value), to_bx(value) + 0.06, yy,
        just = "left", size = 6, col = "#4D4D4D")
  }
}
bay <- 0.92
seg(bx0, bay, bx1, bay, lwd = 0.5)
for (tick in c(0, 5, 10)) {
  xx <- to_bx(tick / 100)
  seg(xx, bay, xx, bay - 0.05, lwd = 0.5)
  txt(sprintf("%d%%", tick), xx, bay - 0.13, size = 6)
}
grid.points(u(3.02), u(0.61), pch = 21, size = unit(0.075, "inches"),
            gp = gpar(fill = COL_DIRECT, col = "#333333", lwd = 0.4))
txt("Direct", 3.10, 0.61, just = "left", size = 6)
grid.points(u(3.66), u(0.61), pch = 22, size = unit(0.075, "inches"),
            gp = gpar(fill = COL_ENZYME, col = "#333333", lwd = 0.4))
txt("Enzyme", 3.74, 0.61, just = "left", size = 6)

seg(0.25, 0.46, 5.25, 0.46, col = COL_RULE, lwd = 0.5)
txt("Test element  •  map cell state  •  link target in 3D", 0.25, 0.27,
    just = "left", size = 6, col = "#4D4D4D")
txt("Do not force nearest-gene labels", 5.25, 0.27, just = "right",
    size = 6, col = COL_DIRECT)
txt("Regulatory DNA context ≠ colocalized transcript biotype  |  context categories overlap",
    0.25, 0.09, just = "left", size = 6, col = "#4D4D4D")

popViewport()
dev.off()

for (f in c("pip_architecture_estimates.tsv", "trait_class_contrasts.tsv",
            "portfolio_independence.tsv", "direct_trait_coding_concentration.tsv",
            "input_manifest.tsv", "environment.txt")) {
  file.copy(file.path(STATS_DIR, f), file.path(OUT_DIR, f), overwrite = TRUE)
}

cap_a <- paste(
  "Distribution of fine-mapping posterior inclusion probability across four mutually",
  "exclusive variant-consequence classes. Within each trait class, the four segments sum",
  "to 100%: protein-altering, canonical splice, synonymous or UTR, and other noncoding.",
  "Displayed locus-balanced estimates are 22.4%, 0.3%, 4.7%, and 72.6% for direct",
  "MASLD/liver-fat traits and 3.4%, 0.0%, 3.6%, and 93.0% for liver-enzyme traits.",
  "To prevent a region recovered in many studies from being counted repeatedly, credible-set",
  "sentinels within 1 Mb on the same chromosome were grouped as one genomic region and",
  "study-specific instances were averaged within that region before pooling. This grouping",
  "does not assert LD or statistical independence. The 24 direct-trait and 259 enzyme-trait",
  "regions are descriptive genomic units, not biological replicates. Cluster-bootstrap",
  "intervals and the credible-set-weighted sensitivity estimator remain in the source tables.",
  "The panel describes allocation of fine-mapping probability, not heritability, effect size,",
  "variant count, or gene-causality probability."
)
cap_b <- paste(
  "Deposited regulatory context of fine-mapped noncoding posterior mass. Descriptive",
  "credible-set-weighted fractions are shown because the estimand is aggregate posterior",
  "mass: 52.11% of direct MASLD/liver-fat mass and 60.43% of liver-enzyme mass has",
  "unresolved deposited context. Only 7.87% and 4.46%, respectively, lies within 2 kb",
  "of a GENCODE TSS; 7.32% and 4.15% is open in an available liver lineage; and 0.23%",
  "and 1.57% overlaps a liver ABC enhancer. Positive context categories overlap and must",
  "not be summed. Across the family of six direct-versus-enzyme tests on these groups,",
  "corrected by Benjamini-Hochberg, the promoter-proximal difference survives (3.92-fold,",
  "q=0.0042) as does the unresolved-context difference (0.53-fold, q=6.0e-4); the liver",
  "ABC enhancer difference does not (q=0.398), nor does generic liver accessibility",
  "(q=0.222). The ABC and accessibility rows are shown for completeness and must not be",
  "described as depletion or enrichment. Groups are operational 1-Mb clusters, not",
  "independent loci: no LD was consulted in forming them. Contrasts and sensitivity",
  "estimates remain in the source tables. Regulatory-DNA context is separate from",
  "colocalized-transcript biotype; unresolved loci should not be assigned automatically",
  "to the nearest gene."
)
writeLines(cap_a, file.path(OUT_DIR, "Fig2D_pip_architecture_by_trait_directness_caption.txt"))
writeLines(cap_b, file.path(OUT_DIR, "Fig2_companion_noncoding_annotation_gap_caption.txt"))

message("Panel A caption: ", cap_a)
message("Panel B caption: ", cap_b)
cat("wrote:", pa_pdf, "\n")
cat("wrote:", pb_pdf, "\n")
