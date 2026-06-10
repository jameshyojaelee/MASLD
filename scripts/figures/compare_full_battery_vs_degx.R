#!/usr/bin/env Rscript
# compare_full_battery_vs_degx.R
# After the parallel all-19-methods re-run (run_degx_battery_array.sh), gather the
# 19 per-method shard tables into deg_tables/ as deg_<method>_degx.csv (consistent
# naming) and compare each, method by method, against the stored canonical degx run
# (~/degx/runs/Rexact/). Both come from the SAME degx engine on the SAME data, so
# this is a reproducibility check: deterministic methods match exactly; stochastic
# ones (sva) may drift.

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEGD   <- file.path(BASE, "figures/supplementary/figS_methods_validation/disease_signature_sweep/deg_tables")
SHARDS <- file.path(BASE, "figures/supplementary/figS_methods_validation/disease_signature_sweep/work/degx_rerun/shards")
DEGX   <- path.expand("~/degx/runs/Rexact")
SIG <- 0.05

methods <- c("deseq2_wald","deseq2_lrt","deseq2_apeglm","deseq2_ashr",
             "edger_qlf","edger_qlf_robust","edger_exact",
             "limma_voom","limma_trend","limma_voom_qw","dream",
             "metafor_deseq2_re","metafor_deseq2_fe","metafor_deseq2_hk",
             "metafor_voom_re","metafor_voom_fe","metafor_voom_hk",
             "combatseq_deseq2","sva_limma")

# --- gather shard tables into deg_tables/ as deg_<method>_degx.csv -----------------
for (m in methods) {
  src <- file.path(SHARDS, m, sprintf("Rexact_disease_vs_control_%s.csv", m))
  if (file.exists(src)) file.copy(src, file.path(DEGD, sprintf("deg_%s_degx.csv", m)), overwrite = TRUE)
}

# --- compare each method: fresh re-run vs stored degx -----------------------------
rows <- rbindlist(lapply(methods, function(m) {
  f_new <- file.path(DEGD, sprintf("deg_%s_degx.csv", m))
  f_old <- file.path(DEGX, sprintf("Rexact_disease_vs_control_%s.csv", m))
  if (!file.exists(f_new)) return(data.table(method = m, status = "rerun_missing"))
  if (!file.exists(f_old)) return(data.table(method = m, status = "degx_missing"))
  a <- fread(f_new)[, .(gene, logFC, padj)]
  b <- fread(f_old)[, .(gene, logFC, padj)]
  mab <- merge(a, b, by = "gene", suffixes = c("_new", "_old"))
  sa <- a[padj < SIG & !is.na(padj), gene]; sb <- b[padj < SIG & !is.na(padj), gene]
  inter <- length(intersect(sa, sb)); uni <- length(union(sa, sb))
  both <- mab[gene %in% intersect(sa, sb)]
  dir_ag <- if (nrow(both)) mean(sign(both$logFC_new) == sign(both$logFC_old)) else NA_real_
  data.table(
    method = m, status = "ok",
    n_rerun = nrow(a), n_degx = nrow(b), n_shared = nrow(mab),
    logFC_spearman = round(suppressWarnings(cor(mab$logFC_new, mab$logFC_old, method="spearman", use="complete.obs")), 4),
    logFC_pearson  = round(suppressWarnings(cor(mab$logFC_new, mab$logFC_old, use="complete.obs")), 4),
    n_sig_rerun = length(sa), n_sig_degx = length(sb),
    jaccard_sig = if (uni) round(inter/uni, 4) else NA_real_,
    direction_agree = round(dir_ag, 4)
  )
}), fill = TRUE)

out <- file.path(DEGD, "comparison_degxrerun_vs_degx.csv")
fwrite(rows, out)
cat("\n=========== All-19-methods re-run vs stored degx ===========\n")
print(rows)
cat("\nGathered", sum(rows$status == "ok", na.rm = TRUE), "of", length(methods), "methods. Wrote:", out, "\n")
