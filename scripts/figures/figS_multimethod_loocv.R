#!/usr/bin/env Rscript
# figS_multimethod_loocv.R
# Compact LOO-CV replication panel for the multi-method DE validation
# (dream vs DESeq2; metafor folds in when enabled). Grouped bars = mean over
# the 5 leave-one-cohort-out folds; points = individual folds. Paired-Wilcoxon
# across folds annotated (all n.s. => methods replicate equivalently).
# Output: figures/supplementary/figS_methods_validation/multimethod_validation/panels/
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
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
long[, method := factor(method, levels = c("dream", "deseq2", "metafor"))]

summ <- long[, .(mean = mean(value, na.rm = TRUE),
                 sd   = sd(value,  na.rm = TRUE)), by = .(method, metric)]

method_cols <- c(dream = "#1b9e77", deseq2 = "#66a61e", metafor = "#1f78b4")

p <- ggplot(summ, aes(metric, mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
                position = position_dodge(width = 0.8), width = 0.2,
                linewidth = 0.3, colour = "grey35") +
  geom_point(data = long, aes(metric, value, fill = method),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
             size = 0.9, shape = 21, colour = "white", stroke = 0.2, alpha = 0.8) +
  scale_fill_manual(values = method_cols, name = NULL,
                    labels = c(dream = "dream", deseq2 = "DESeq2", metafor = "metafor")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Held-out replication (5-fold LOO-CV)",
       title = "Leave-one-cohort-out replication across methods",
       subtitle = "bars = mean ± SD over 5 folds; points = folds  |  dream≈DESeq2 (n.s.); metafor: fewer but higher-confidence DEGs") +
  theme_masld(base_size = 7) +
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 6, colour = "grey40"),
        axis.text.x = element_text(size = 6))

ggsave(file.path(OUT, "panelA_loocv_replication.pdf"), p,
       width = 5.6, height = 3.4, useDingbats = FALSE)
cat("Wrote panelA_loocv_replication.pdf\n")
print(summ[order(metric, method)])
