#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# T1 finalize — dv_results.csv → dv_* atlas columns + headline 2×2 + figure.
#   · dv_direction: fan-out (more variable in disease) / canalize (less).
#   · additive atlas slot `dv_atlas_columns.tsv` (keyed by human_symbol, 27a-merge).
#   · Jaccard(variance-only DV set, Tier-1 mean DEGs) — the orthogonality headline.
#   · heterogeneity-vs-mean-shift scatter (dv_t vs bulk_tstat, coloured by partition).
# Lightweight (reads one csv) — safe off-node.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
P  <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
d  <- fread(file.path(P, "dv_results.csv"))

d[, dv_direction := fcase(dv_sig & dv_t > 0, "fan_out",
                          dv_sig & dv_t < 0, "canalize", default = "ns")]
atlas_cols <- d[, .(human_symbol = symbol, gene_ensembl = gene,
                    dv_t, dv_logfc, dv_p_perm, dv_q_perm, dv_p_pool, dv_q_pool,
                    dv_nb_logratio = nb_logratio, dv_nb_concordant = nb_concordant,
                    dv_loco_signfrac = loco_signfrac, dv_loco_pass = loco_pass,
                    dv_sig, dv_direction, dv_partition = partition)]
atlas_cols <- atlas_cols[!is.na(human_symbol) & human_symbol != ""][order(-dv_sig, dv_q_pool)]
fwrite(atlas_cols, file.path(P, "dv_atlas_columns.tsv"), sep = "\t")

# headline numbers
n_tier1 <- d[is_mean_deg == TRUE, .N]; n_dv <- d[dv_sig == TRUE, .N]
vonly <- d[dv_sig == TRUE & is_mean_deg == FALSE, unique(gene_base)]
tier1 <- d[is_mean_deg == TRUE, unique(gene_base)]
dvset <- d[dv_sig == TRUE, unique(gene_base)]
jac <- length(intersect(dvset, tier1)) / length(union(dvset, tier1))
part <- d[, .N, by = partition][order(-N)]
fanc <- d[dv_sig == TRUE, .N, by = dv_direction]

cat("── T1 headline ──\n")
print(part)
cat(sprintf("\nTier-1 mean DEGs        = %d\nDV-sig genes           = %d\nvariance-ONLY (novel)  = %d\n",
            n_tier1, n_dv, d[partition == "variance_only", .N]))
cat(sprintf("Jaccard(DV-sig, Tier-1) = %.3f  (low ⇒ orthogonal axes)\n", jac))
cat("DV direction among hits:\n"); print(fanc)

summ <- data.table(metric = c("tier1_mean_degs","dv_sig","variance_only","mean_and_variance",
                              "jaccard_dv_tier1","fan_out","canalize","n_universe"),
                   value = c(n_tier1, n_dv, d[partition=="variance_only", .N],
                            d[partition=="mean_and_variance", .N], round(jac,4),
                            fanc[dv_direction=="fan_out", N], fanc[dv_direction=="canalize", N],
                            nrow(d)))
fwrite(summ, file.path(P, "t1_summary.csv"))

# ── headline scatter: heterogeneity (DV) vs mean-shift, coloured by partition ──
pal <- c(variance_only = "#C2185B", mean_and_variance = "#6A1B9A",
         mean_only = "#1565C0", neither = "grey80")
dp <- d[!is.na(bulk_tstat)]
g <- ggplot(dp, aes(bulk_tstat, dv_t, color = partition)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey50") +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = "Mean-shift (bulk DEG t-statistic)", y = "Differential variability (DV t)",
       title = sprintf("Variance carries orthogonal information: %d variance-only genes",
                       d[partition=="variance_only", .N])) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
ggsave(file.path(P, "t1_heterogeneity_vs_meanshift.pdf"), g, width = 5.2, height = 4, device = cairo_pdf)
cat(sprintf("\nWrote dv_atlas_columns.tsv (%d genes), t1_summary.csv, t1_heterogeneity_vs_meanshift.pdf\n",
            nrow(atlas_cols)))
