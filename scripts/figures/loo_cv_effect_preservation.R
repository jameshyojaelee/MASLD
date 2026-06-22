#!/usr/bin/env Rscript
# KEY MESSAGE: leave-one-cohort-out refits barely move the integrated effect sizes —
# per-gene full-model log2FC vs each 4-cohort leave-one-out log2FC hug the y=x
# diagonal (pooled Spearman rho 0.95; among DEGs ~99.9% keep their sign), with only
# the largest held-out cohort (Chen) loosening the cloud.
# ============================================================================
# loo_cv_effect_preservation.R — Figure 3 LOO-CV effect-size stability (scatter)
#   full 5-cohort integrated shrunk_logFC vs each leave-one-out 4-cohort fit.
#   Source: .../results/integration/{canonical_deg_results.csv,
#            loo_cv_C2/lvqw_loo_C2_<cohort>.csv}.
# Output: FIG2_DIR/panels/figs3c_loo_cv_effect_preservation.pdf
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
C2  <- file.path(INT, "loo_cv_C2")

# Full 5-cohort integrated effect sizes (ashr-shrunk)
full <- fread(file.path(INT, "canonical_deg_results.csv"))[, .(gene, full = shrunk_logFC)]

# Each leave-one-cohort-out 4-cohort refit; pool all gene-fold pairs
folds <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
           GSE162694 = "GSE162694",   GSE213621 = "GSE213621")
dt <- rbindlist(lapply(names(folds), function(g) {
  f <- fread(file.path(C2, sprintf("lvqw_loo_C2_%s.csv", g)))[, .(gene, loo = shrunk_logFC)]
  m <- merge(full, f, by = "gene")
  m[, held_out := folds[g]]
  m
}))
dt <- dt[is.finite(full) & is.finite(loo)]

# Pooled statistics for annotation
rho_all  <- cor(dt$full, dt$loo, method = "spearman")
deg      <- dt[abs(full) > 0.5]                       # full-model effect-size DEGs
dir_deg  <- 100 * mean(sign(deg$full) == sign(deg$loo))
n_pairs  <- nrow(dt)
lim      <- as.numeric(quantile(abs(dt$full), 0.999, na.rm = TRUE))
lim      <- max(lim, 3.0)

# DEG-set-overlap means (recovery, Jaccard) from the LOO-CV summary, so the panel
# carries all four robustness numbers the manuscript cites.
sm        <- fread(file.path(C2, "loo_cv_C2_summary.csv"))
rec_mean  <- mean(sm$pct_full_recovered_01)
jac_mean  <- mean(sm$full_vs_train_jaccard_01)

ann <- sprintf("Spearman ρ = %.2f  ·  %.0f%% of DEGs keep sign\nrecovery %.1f%%  ·  Jaccard %.2f  ·  n = %s",
               rho_all, dir_deg, rec_mean, jac_mean, format(n_pairs, big.mark = ","))

p <- ggplot(dt, aes(full, loo)) +
  geom_hline(yintercept = 0, color = "gray85", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "gray85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray45", linewidth = 0.35) +
  geom_hex(bins = 80) +
  scale_fill_gradient(low = "#E8E2EF", high = masld_colors$mash,
                      trans = "log10", name = "gene-fold\npairs") +
  annotate("text", x = -lim * 0.98, y = lim * 0.96, hjust = 0, vjust = 1,
           label = ann, size = 2.5, fontface = "plain", color = "gray20",
           lineheight = 0.95) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = expression("Full 5-cohort log"[2]*"FC"),
       y = expression("Leave-one-out log"[2]*"FC"),
       title = "Integrated effect sizes are preserved under leave-one-cohort-out") +
  theme_masld(base_size = 9) +
  theme(plot.title = element_text(size = 9.5, face = "bold", margin = margin(b = 4)),
        legend.position = "right",
        legend.key.width = unit(0.3, "cm"),
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 6.5))

save_fig(p, file.path(PANEL_DIR, "figs3c_loo_cv_effect_preservation.pdf"),
         width = fig_half_width * 1.5, height = 3.3)

fwrite(data.table(metric = c("pooled_spearman_rho", "deg_direction_pct", "n_gene_fold_pairs"),
                  value  = c(round(rho_all, 4), round(dir_deg, 2), n_pairs)),
       file.path(DATA_DIR, "loo_cv_effect_preservation_stats.csv"))

cat(sprintf("Wrote figs3c_loo_cv_effect_preservation.pdf\n  pooled rho=%.3f  DEG-direction=%.1f%%  n_pairs=%d\n",
            rho_all, dir_deg, n_pairs))
