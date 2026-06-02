#!/usr/bin/env Rscript
# Aggregate sex-permutation array results into final empirical p-values.
suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_v6_permutation_null"
rds_files <- list.files(OUT_DIR, pattern = "^perm_task_.*\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "task files\n")
stopifnot(length(rds_files) > 0)

results <- lapply(rds_files, readRDS)
gene_ids <- results[[1]]$gene_ids
t_obs <- results[[1]]$t_obs
total_perm <- sum(sapply(results, `[[`, "n_perm"))
n_extreme_total <- Reduce(`+`, lapply(results, `[[`, "n_extreme"))

p_emp <- (n_extreme_total + 1L) / (total_perm + 1L)

out <- data.table(
  gene = gene_ids,
  t_obs = t_obs,
  n_extreme = n_extreme_total,
  n_perm = total_perm,
  p_emp = p_emp
)
out[, padj_bh := p.adjust(p_emp, method = "BH")]

sig_005 <- sum(out$padj_bh < 0.05, na.rm = TRUE)
sig_010 <- sum(out$padj_bh < 0.10, na.rm = TRUE)
cat("Total permutations:", total_perm, "\n")
cat("Genes:", nrow(out), "\n")
cat("Significant at BH<0.05:", sig_005, "\n")
cat("Significant at BH<0.10:", sig_010, "\n")

fwrite(out, file.path(OUT_DIR, "sex_v6_permutation_null.csv"))
cat("Saved:", file.path(OUT_DIR, "sex_v6_permutation_null.csv"), "\n")
