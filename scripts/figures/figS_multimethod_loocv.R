#!/usr/bin/env Rscript
# figS_multimethod_loocv.R
# Compact LOO-CV replication panel for the multi-method DE validation. Plots
# WHATEVER DE engines are present in loocv_perfold.csv (no hard-coded method
# filter); readable engine names + colours come from method_correction_labels.R.
# Grouped bars = mean over the 5 leave-one-cohort-out folds; points = folds.
# Output: figures/supplementary/figS_methods_validation/multimethod_validation/panels/
#   loocv_replication.pdf  (was panelA_loocv_replication.pdf)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/method_correction_labels.R"))
LOO  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/loocv")
OUT  <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

pf <- fread(file.path(LOO, "loocv_perfold.csv"))

# [0,1]-scaled replication metrics (drop unbounded fisher_or)
metrics <- c(jaccard_01            = "Jaccard\n(padj<0.1)",
             lfc_spearman          = "logFC\nSpearman",
             direction_concordance = "Direction\nconcordance",
             auc_replication       = "AUC\nreplication",
             pi1_heldout           = "pi1\n(held-out)")
pf[, direction_concordance := direction_concordance / 100]   # %->[0,1]

long <- melt(pf[, c("method", "held_out_cohort", names(metrics)), with = FALSE],
             id.vars = c("method", "held_out_cohort"),
             variable.name = "metric", value.name = "value")
long[, metric := factor(metrics[as.character(metric)], levels = metrics)]
# Plot WHATEVER methods are present in loocv_perfold.csv (no hard-coded method
# filter). Order by engine family then label so related engines sit together;
# readable names + colours come from the shared method_correction_labels.R lookup.
present <- unique(as.character(long$method))
present <- present[order(match(engine_family(present), ENGINE_FAMILY_LEVELS),
                         engine_label(present))]
long[, method := factor(method, levels = present)]

summ <- long[, .(mean = mean(value, na.rm = TRUE),
                 sd   = sd(value,  na.rm = TRUE)), by = .(method, metric)]

# Readable legend labels from the lookup; distinct per-method hues generated so
# bars stay distinguishable even within an engine family (family colour alone
# would collide the limma/edgeR variants).
method_labs <- setNames(engine_label(present), present)
fam_present <- engine_family(present)
method_cols <- setNames(
  vapply(present, function(m) {
    same <- present[fam_present == engine_family(m)]
    base <- ENGINE_FAMILY_PALETTE[[engine_family(m)]]
    if (is.null(base)) base <- "#BDBDBD"
    if (length(same) <= 1) return(base)
    ramp <- colorRampPalette(c(colorspace::lighten(base, 0.35),
                               colorspace::darken(base, 0.35)))(length(same))
    ramp[match(m, same)]
  }, character(1)), present)
# colorspace may be unavailable -> fall back to family colour (acceptable collision)
if (!requireNamespace("colorspace", quietly = TRUE))
  method_cols <- setNames(ENGINE_FAMILY_PALETTE[engine_family(present)], present)

p <- ggplot(summ, aes(metric, mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
                position = position_dodge(width = 0.8), width = 0.2,
                linewidth = 0.3, colour = "grey35") +
  geom_point(data = long, aes(metric, value, fill = method),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
             size = 0.9, shape = 21, colour = "white", stroke = 0.2, alpha = 0.8) +
  scale_fill_manual(values = method_cols, name = NULL, labels = method_labs) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Held-out replication (5-fold LOO-CV)",
       title = "Leave-one-cohort-out replication across methods",
       subtitle = sprintf("bars = mean ± SD over 5 folds; points = folds  |  %d DE engines",
                          length(present))) +
  theme_masld(base_size = 7) +
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 6, colour = "grey40"),
        axis.text.x = element_text(size = 6))

ggsave(file.path(OUT, "loocv_replication.pdf"), p,
       width = 5.6, height = 3.4, useDingbats = FALSE)
cat("Wrote loocv_replication.pdf\n")
print(summ[order(metric, method)])
