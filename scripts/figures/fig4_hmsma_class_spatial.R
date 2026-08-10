#!/usr/bin/env Rscript
# Candidate panel — evidence-class spatial organisation across 35 HMSMA Visium samples.
#
# Shows per-sample excess spatial autocorrelation (Moran's I) over an
# expression/detection-matched null, by evidence class. The unit of replication is
# the SAMPLE: one point per sample per class, with mean and 95% CI.
#
# Source: Analysis/Spatial/results/hmsma_class_spatial/ (job 19663173)
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

IN  <- file.path(BASE, "Analysis/Spatial/results/hmsma_class_spatial")
OUT <- file.path(FIG4_DIR, "panels")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

per <- fread(file.path(IN, "per_sample_class_delta.tsv"))
sm  <- fread(file.path(IN, "class_summary.tsv"))

stopifnot(nrow(per) > 0, uniqueN(per$sample) == 35L)

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

# assert the headline numbers rather than trusting the file
d_ds <- sm[cls == "disease_state_only"]
d_go <- sm[cls == "genetic_only"]
stopifnot(d_ds$n_samples_positive == 35L, d_ds$n_samples == 35L)
ratio <- d_ds$mean_delta / d_go$mean_delta
stopifnot(ratio > 6.5, ratio < 7.5)

p <- ggplot(per, aes(x = delta, y = cls, colour = cls)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#9E9E9E") +
  geom_jitter(height = 0.16, size = 0.5, alpha = 0.55, stroke = 0) +
  geom_errorbar(data = sm, inherit.aes = FALSE,
                aes(y = cls, xmin = ci_lo, xmax = ci_hi),
                orientation = "y", width = 0, linewidth = 0.5, colour = "black") +
  geom_point(data = sm, inherit.aes = FALSE,
             aes(x = mean_delta, y = cls), size = 1.4, colour = "black") +
  scale_colour_manual(values = col, guide = "none") +
  scale_y_discrete(labels = lab) +
  labs(x = "Excess spatial autocorrelation over matched null (Moran's I)",
       y = NULL) +
  theme_masld(base_size = 6) +
  theme(panel.grid.major.y = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"))

ggsave(file.path(OUT, "fig4_hmsma_class_spatial.pdf"), p,
       width = 3.5, height = 1.9, useDingbats = FALSE)

fwrite(per, file.path(OUT, "data", "fig4_hmsma_class_spatial_source.csv"))

message(sprintf(
paste0("Evidence-class spatial organisation, 35 HMSMA Visium samples. Points are ",
"per-sample excess Moran's I over an expression- and detection-matched null drawn ",
"from the 'neither' class; black points and bars are the mean and 95%% CI across ",
"samples (sample is the unit of replication). Disease-state-only genes show %.1fx ",
"the excess of genetic-only genes (paired within sample, %d/%d samples, ",
"p = 3.1e-09). Genetic-only genes are above their matched null (%d/%d samples, ",
"p = 1.2e-05) but far below disease-state-only: a gradient, not a dissociation. ",
"The convergent class is a median of 18 genes per sample and is shown for ",
"completeness only. Samples are unlabelled for disease state, so this is ",
"within-tissue organisation pooled across the cohort, not a disease contrast; ",
"on-tissue status is a UMI-threshold proxy and coordinates are Visium v1 array ",
"positions rather than image-registered."),
ratio, 35L, 35L, d_go$n_samples_positive, 35L))
