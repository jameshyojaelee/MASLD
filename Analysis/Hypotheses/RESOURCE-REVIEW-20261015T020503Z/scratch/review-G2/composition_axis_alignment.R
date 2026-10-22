# Light diagnostic (login node, no model fit): does the released F1-vs-F0 effect vector align with the
# within-cohort gene-wise association of expression with the complexity covariates?
suppressPackageStartupMessages({library(edgeR); library(data.table)})
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
d <- readRDS(file.path(REL, "merged_dge.rds"))
s <- fread(file.path(RUN, "jobs/1721/stage_roster_with_complexity.tsv"))
lc <- cpm(d[, s$sample_id], log = TRUE, prior.count = 2)   # TMM factors carried in d$samples
X0 <- model.matrix(~ dataset + inferred_sex, s)
res <- function(M, X) M - t(X %*% qr.solve(X, t(M)))
R <- res(lc, X0)
assoc <- sapply(c("z_loglib", "z_gini", "z_nf"), function(v) { z <- res(matrix(s[[v]], nrow = 1), X0)[1, ]
  as.numeric(R %*% z) / sum(z^2) })     # per-gene slope of residual logCPM on residual covariate
rownames(assoc) <- rownames(d)
ext <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv",
             select = c("gene_id_versioned", "logFC", "FDR", "contrast", "axis"))[axis == "fibrosis"]
out <- rbindlist(lapply(c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0"), function(cn) {
  e <- ext[contrast == cn]; a <- assoc[e$gene_id_versioned, ]
  sig <- e$FDR < 0.05
  data.table(contrast = cn, n_sig = sum(sig),
    rho_all_lib = cor(e$logFC, a[, 1], method = "spearman"), rho_all_gini = cor(e$logFC, a[, 2], method = "spearman"),
    rho_all_nf = cor(e$logFC, a[, 3], method = "spearman"),
    rho_sig_gini = cor(e$logFC[sig], a[sig, 2], method = "spearman"), rho_sig_nf = cor(e$logFC[sig], a[sig, 3], method = "spearman"),
    rho_sig_lib = cor(e$logFC[sig], a[sig, 1], method = "spearman"))
}))
print(out, digits = 3)
fwrite(out, file.path(RUN, "scratch/review-G2/composition_axis_alignment.tsv"), sep = "\t")
