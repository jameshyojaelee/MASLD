#!/usr/bin/env Rscript
# Candidate panel — evidence-class spatial organisation across 35 HMSMA Visium arrays.
#
# Shows excess spatial autocorrelation (Moran's I) over an expression- and
# detection-matched gene null, by evidence class.
#
# The unit of inference is the MATCHED GENE SET with the 35 arrays held fixed, not
# the array. HMSMA has no array-to-donor key, so arrays cannot be treated as
# independent biological replicates and no across-array t-test or CI is drawn.
# Per-array points are descriptive; the gray band is the matched-gene null.
#
# Source: Analysis/Spatial/results/hmsma_class_spatial_v2/
# Estimand matches the existing spatial arm (GSE192741, Vu):
#   residual_spatial_autocorrelation_against_matched_gene_null
#
# Class colours are taken verbatim from scripts/figures/fig1_fig4_evidence_classes.R
# so this panel reads as part of the same evidence-class family.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

IN  <- Sys.getenv("HMSMA_CLASS_SPATIAL_DIR",
  file.path(BASE, "Analysis/Spatial/results/hmsma_class_spatial_v3"))
OUT <- file.path(FIG5_CONTEXT_DIR, "panels")
dir.create(file.path(OUT, "data"), showWarnings = FALSE, recursive = TRUE)

per <- fread(file.path(IN, "per_array_class_delta.tsv"))
sm  <- fread(file.path(IN, "class_summary.tsv"))
ct  <- fread(file.path(IN, "contrast_summary.tsv"))

# structural checks only — never assert a specific effect size or p-value here,
# or a rerun that legitimately moves the numbers will silently ship a stale caption
stopifnot(nrow(per) > 0, nrow(ct) == 1L)
stopifnot(all(c("genetic_only", "disease_state_only", "convergent") %in% sm$cls))
stopifnot(all(is.finite(sm$p_matched_gene_null)), is.finite(ct$p_matched_gene_null))
stopifnot(uniqueN(per$sample) == max(sm$n_arrays))

lab <- c(genetic_only = "Genetic-only",
         convergent = "Convergent",
         disease_state_only = "Disease-state-only")
col <- c(genetic_only = "#1565C0",
         convergent = "#00695C",
         disease_state_only = "#C9265E")

ord <- c("genetic_only", "convergent", "disease_state_only")
per[, cls := factor(cls, levels = ord)]
sm[,  cls := factor(cls, levels = ord)]
per <- per[!is.na(delta)]

fmt_p <- function(p, n_null) {
  floor_p <- 1 / (n_null + 1)
  ifelse(p <= floor_p, sprintf("< %.0e", floor_p), sprintf("%.2g", p))
}

p <- ggplot(per, aes(x = delta, y = cls)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#9E9E9E") +
  geom_segment(data = sm, inherit.aes = FALSE,
               aes(y = cls, yend = cls, x = null_q050, xend = null_q950),
               linewidth = 2.2, colour = "#9E9E9E", lineend = "butt") +
  geom_jitter(aes(colour = cls), height = 0.16, size = 0.5, alpha = 0.55, stroke = 0) +
  geom_point(data = sm, inherit.aes = FALSE,
             aes(x = null_q500, y = cls), shape = 23, size = 1.0,
             fill = "#616161", colour = "#616161") +
  geom_point(data = sm, inherit.aes = FALSE,
             aes(x = observed_delta, y = cls), size = 1.5, colour = "black") +
  scale_colour_manual(values = col, guide = "none") +
  scale_y_discrete(labels = lab) +
  labs(x = "Excess spatial autocorrelation over matched-gene null (Moran's I)",
       y = NULL) +
  theme_masld(base_size = 6) +
  theme(panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"))

ggsave(file.path(OUT, "fig5_hmsma_class_spatial.pdf"), p,
       width = 3.5, height = 1.9, useDingbats = FALSE)

fwrite(per, file.path(OUT, "data", "fig5_hmsma_class_spatial_source.csv"))
fwrite(sm,  file.path(OUT, "data", "fig5_hmsma_class_spatial_summary.csv"))
fwrite(ct,  file.path(OUT, "data", "fig5_hmsma_class_spatial_contrast.csv"))

d_ds <- sm[cls == "disease_state_only"]
d_go <- sm[cls == "genetic_only"]
d_cv <- sm[cls == "convergent"]

# Verdict language follows the p-value rather than a fixed interpretation, and
# uses the Resource's evidence-state vocabulary: failing to exceed the matched
# null is indeterminate, never negative. The ratio of class deltas is
# deliberately NOT reported: the genetic-only delta sits on zero, so a ratio is
# unstable and reads as an effect size it cannot support.
verdict <- function(row) {
  # plural subject in the caption ("... genes %s")
  if (row$p_matched_gene_null < 0.05) {
    "exceed matched-gene background"
  } else {
    "do not exceed matched-gene background (indeterminate, not negative)"
  }
}

message(sprintf(
paste0("Evidence-class spatial organisation, %d HMSMA Visium arrays. Coloured points ",
"are per-array excess Moran's I over an expression- and detection-matched null drawn ",
"from the 'neither' class; gray bars are the 5th-95th percentiles of that matched-gene ",
"null and gray diamonds its median; black points are the observed cohort value. ",
"Disease-state-only genes %s (%+.4f, p %s, %d/%d arrays positive). Genetic-only genes ",
"%s (%+.4f, p %s, %d/%d arrays positive). Convergent genes %s (%+.4f, p %s), a median ",
"of %.0f genes per array, shown for completeness only. The within-array paired ",
"difference between disease-state-only and genetic-only is %+.4f (p %s, %d/%d arrays ",
"positive), so the two classes are separated by matched-gene calibration even though ",
"both carry positive raw autocorrelation in most arrays. The unit of inference is the ",
"matched gene set with arrays held fixed; HMSMA has no array-to-donor key, so arrays ",
"are not independent donors, no across-array test or confidence interval is drawn, and ",
"these p-values are conditional on these arrays rather than generalizable to new ",
"donors. Arrays are unlabelled for disease state, so this is within-tissue organisation ",
"pooled across the cohort, not a disease contrast; on-tissue status is a UMI-threshold ",
"proxy and coordinates are Visium v1 array positions rather than image-registered."),
max(sm$n_arrays),
verdict(d_ds), d_ds$observed_delta, fmt_p(d_ds$p_matched_gene_null, d_ds$n_null),
d_ds$n_arrays_delta_positive, d_ds$n_arrays,
verdict(d_go), d_go$observed_delta, fmt_p(d_go$p_matched_gene_null, d_go$n_null),
d_go$n_arrays_delta_positive, d_go$n_arrays,
verdict(d_cv), d_cv$observed_delta, fmt_p(d_cv$p_matched_gene_null, d_cv$n_null),
median(per[cls == "convergent", n]),
ct$observed_difference, fmt_p(ct$p_matched_gene_null, ct$n_null),
ct$n_arrays_difference_positive, ct$n_arrays))
