#!/usr/bin/env Rscript
# compare_deg_tables_vs_degx.R
# Head-to-head concordance between the fresh re-run DEG tables
# (figS_disease_signature_sweep_deg_tables.R; 27,638-gene sweep universe) and the
# pre-existing degx Rexact disease_vs_control battery (~/degx/runs/Rexact/).
# Per overlapping method: gene-universe overlap, logFC Spearman/Pearson on shared
# genes, padj<0.05 DEG counts side by side, Jaccard of DEG sets, direction agreement.

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEGDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/disease_signature_sweep/deg_tables")
DEGX   <- path.expand("~/degx/runs/Rexact")
SIG <- 0.05

# my re-run method  ->  degx Rexact file stem (disease_vs_control)
map <- list(
  limma_voom       = "limma_voom",
  dream            = "dream",
  deseq2           = "deseq2_wald",
  edger_qlf        = "edger_qlf",
  metafor_re       = "metafor_voom_re",
  combatseq_deseq2 = "combatseq_deseq2"
)

rows <- rbindlist(lapply(names(map), function(m) {
  f_mine <- file.path(DEGDIR, sprintf("deg_%s_indep.csv", m))
  f_degx <- file.path(DEGX, sprintf("Rexact_disease_vs_control_%s.csv", map[[m]]))
  if (!file.exists(f_mine)) { message("missing re-run: ", f_mine); return(NULL) }
  if (!file.exists(f_degx)) { message("missing degx:   ", f_degx); return(NULL) }

  a <- fread(f_mine)[, .(gene, logFC, padj)]          # re-run
  b <- fread(f_degx)[, .(gene, logFC, padj)]          # degx
  mab <- merge(a, b, by = "gene", suffixes = c("_me", "_dx"))

  sig_me  <- a[padj < SIG & !is.na(padj), gene]
  sig_dx  <- b[padj < SIG & !is.na(padj), gene]
  inter   <- length(intersect(sig_me, sig_dx))
  uni     <- length(union(sig_me, sig_dx))

  # direction agreement among genes sig in BOTH
  both <- mab[gene %in% intersect(sig_me, sig_dx)]
  dir_agree <- if (nrow(both) > 0) mean(sign(both$logFC_me) == sign(both$logFC_dx)) else NA_real_

  data.table(
    method          = m,
    degx_method     = map[[m]],
    n_genes_rerun   = nrow(a),
    n_genes_degx    = nrow(b),
    n_shared        = nrow(mab),
    logFC_spearman  = round(cor(mab$logFC_me, mab$logFC_dx, method = "spearman", use = "complete.obs"), 4),
    logFC_pearson   = round(cor(mab$logFC_me, mab$logFC_dx, use = "complete.obs"), 4),
    n_sig_rerun     = length(sig_me),
    n_sig_degx      = length(sig_dx),
    n_sig_both      = inter,
    jaccard_sig     = round(inter / uni, 4),
    direction_agree = round(dir_agree, 4)
  )
}))

out <- file.path(DEGDIR, "comparison_indep_vs_degx.csv")
fwrite(rows, out)

cat("\n=========== DEG re-run vs degx Rexact — concordance ===========\n")
print(rows)
cat("\nGene universe: re-run", unique(rows$n_genes_rerun), "vs degx", unique(rows$n_genes_degx), "\n")
cat("Wrote:", out, "\n")
