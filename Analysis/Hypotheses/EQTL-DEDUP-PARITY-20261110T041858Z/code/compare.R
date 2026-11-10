#!/usr/bin/env Rscript
# Old (COLOC r2) vs refit comparison for the 90 affected primary genes.
suppressPackageStartupMessages(library(data.table))
run <- Sys.getenv("RUN_DIR")
fm <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
r2 <- file.path(fm, "results/susie_coloc_r2_20260924T053244Z")
genes <- fread(file.path(run, "inputs/affected_primary_genes.tsv"))
filt <- fread(file.path(run, "inputs/gene_study_filter.tsv"))
studies <- unique(filt$gwas_name)
out <- file.path(run, "results"); dir.create(out, showWarnings = FALSE)

# ---- old gene-study rows (r2 master) ----
old <- fread(file.path(r2, "susie_coloc_all_gwas.csv"),
             select = c("ensembl", "gwas_name", "PP.H4.abf", "PP.H4.susie", "method"))
old <- old[ensembl %in% genes$ensembl & gwas_name %in% studies]
setnames(old, c("PP.H4.abf", "PP.H4.susie", "method"), c("old_abf", "old_susie", "old_method"))

# ---- new gene-study rows ----
new_files <- list.files(file.path(run, "coloc"), pattern = "^susie_coloc_chr[0-9]+\\.csv$",
                        recursive = TRUE, full.names = TRUE)
new <- rbindlist(lapply(new_files, fread), fill = TRUE)
new <- new[ensembl %in% genes$ensembl, .(ensembl, gwas_name, new_abf = PP.H4.abf,
                                         new_susie = PP.H4.susie, new_method = method)]
gs <- merge(filt[, .(gwas_name, ensembl, gene)], old, by = c("gwas_name", "ensembl"), all.x = TRUE)
gs <- merge(gs, new, by = c("gwas_name", "ensembl"), all.x = TRUE)
expected_runs <- unique(filt[, .(gwas_name, chr)])
done_runs <- unique(data.table(
  gwas_name = basename(dirname(new_files)),
  chr = as.integer(sub("^susie_coloc_chr([0-9]+)\\.csv$", "\\1", basename(new_files)))))
missing_runs <- fsetdiff(expected_runs, done_runs)
fwrite(missing_runs, file.path(out, "missing_study_chr_runs.tsv"), sep = "\t")

# ---- eQTL refit status and coverage vs old fit ----
summ <- rbindlist(lapply(list.files(file.path(run, "eqtl_susie"), "_summary.tsv$", full.names = TRUE),
                         fread), fill = TRUE)
loss_files <- list.files(file.path(run, "identity_losses"), "^eqtl_chr", full.names = TRUE)
eloss <- if (length(loss_files)) rbindlist(lapply(loss_files, fread)) else data.table()
cov <- rbindlist(lapply(seq_len(nrow(genes)), function(i) {
  g <- genes[i]; ch <- g$chr
  of <- readRDS(file.path(fm, "results/eqtl_susie_polyfun", paste0("chr", ch), paste0(g$ensembl, "_susie.rds")))
  oid <- colnames(of$lbf_variable); opos <- as.integer(sub("^[0-9]+:", "", oid))
  cs_idx <- unique(unlist(of$sets$cs))
  nf_path <- file.path(run, "eqtl_susie", paste0("chr", ch), paste0(g$ensembl, "_susie.rds"))
  res <- data.table(ensembl = g$ensembl, old_columns = length(oid),
    old_duplicate_columns = sum(duplicated(oid)), old_unique_positions = uniqueN(opos),
    old_converged = isTRUE(of$converged), old_n_cs = length(of$sets$cs))
  if (file.exists(nf_path)) {
    nf <- readRDS(nf_path); nid <- colnames(nf$lbf_variable)
    npos <- as.integer(sub("^[0-9]+:([0-9]+):.*$", "\\1", nid))
    lost <- setdiff(unique(opos), npos)
    res[, `:=`(new_columns = length(nid), new_duplicate_ids = sum(duplicated(nid)),
      new_positions_with_two_events = sum(duplicated(npos)),
      positions_lost = length(lost), positions_gained = length(setdiff(npos, opos)),
      lost_positions = paste(lost, collapse = ";"),
      lost_positions_in_old_cs = sum(opos[cs_idx] %in% lost),
      old_cs_alpha_on_lost = if (length(of$sets$cs)) max(vapply(names(of$sets$cs), function(nm) {
        l <- as.integer(sub("L", "", nm)); sum(of$alpha[l, opos %in% lost]) }, 0)) else 0,
      new_converged = isTRUE(nf$converged), new_n_cs = length(nf$sets$cs))]
  }
  res
}), fill = TRUE)
if (nrow(eloss)) {
  el <- eloss[, .(rejected_rows = .N, rejected_multiple_reference = sum(reason == "multiple_reference_records"),
                  rejected_multiple_summary = sum(reason == "multiple_summary_rows"),
                  min_p_rejected = min(p)), by = .(ensembl = locus)]
  cov <- merge(cov, el, by = "ensembl", all.x = TRUE)
}
# eQTL lead p for the gene (all rows of the source table) for context
lead <- rbindlist(lapply(unique(genes$chr), function(ch) {
  e <- fread(sprintf("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/broadaway_eqtl/chr%d_marginal_summary_results.tsv", ch),
             select = c("ENSG", "PVAL"))
  e[ENSG %in% genes$ensembl, .(eqtl_lead_p = min(PVAL, na.rm = TRUE)), by = .(ensembl = ENSG)]
}))
cov <- merge(cov, lead, by = "ensembl", all.x = TRUE)
cov <- merge(cov, summ[, .(ensembl = ensg, refit_status = status, refit_n_snps = n_snps)],
             by = "ensembl", all.x = TRUE)

# ---- gene states ----
untest <- function(m) !is.na(m) & grepl("^susie_untestable", m)
gene_tab <- gs[, .(
  old_best_susie = suppressWarnings(max(old_susie[old_method == "susie"], na.rm = TRUE)),
  old_driver = gwas_name[old_method == "susie"][which.max(old_susie[old_method == "susie"])][1],
  new_best_susie = suppressWarnings(max(new_susie[new_method %in% "susie"], na.rm = TRUE)),
  new_driver = gwas_name[new_method %in% "susie"][which.max(new_susie[new_method %in% "susie"])][1],
  new_studies_susie = sum(new_method %in% "susie"),
  new_studies_gt05 = sum(new_method %in% "susie" & new_susie > 0.5),
  old_studies_gt05 = sum(old_method %in% "susie" & old_susie > 0.5),
  new_untestable_studies = sum(untest(new_method)),
  new_methods = paste(sort(unique(na.omit(new_method))), collapse = ";"),
  old_best_abf = max(old_abf, na.rm = TRUE), new_best_abf = max(new_abf, na.rm = TRUE),
  new_rows = sum(!is.na(new_method))), by = .(ensembl, gene)]
for (v in c("old_best_susie", "new_best_susie")) set(gene_tab, which(!is.finite(gene_tab[[v]])), v, NA_real_)
gene_tab <- merge(genes[, .(ensembl, chr, in_refit)], gene_tab, by = "ensembl", all.x = TRUE)
gene_tab[, gene := genes$gene[match(ensembl, genes$ensembl)]]
gene_tab <- merge(gene_tab, cov, by = "ensembl", all.x = TRUE)
gene_tab[, eqtl_usable := refit_status %in% "ok" & new_converged %in% TRUE]
gene_tab[, new_state := fifelse(in_refit == FALSE, "untestable_unresolved_source_event",
  fifelse(!is.na(new_best_susie) & new_best_susie > 0.5, "positive",
  fifelse(!eqtl_usable | new_untestable_studies > 0, "untestable", "lost")))]
gene_tab[, new_state_gt08 := fifelse(new_state == "positive" & new_best_susie > 0.8, "positive_gt0.8", new_state)]
gene_tab[, driver_changed := new_state == "positive" & new_driver != old_driver]
setorder(gene_tab, new_state, -old_best_susie)
fwrite(gene_tab, file.path(out, "eqtl_dedup_gene_comparison.tsv"), sep = "\t")
fwrite(gs, file.path(out, "eqtl_dedup_gene_study_comparison.tsv"), sep = "\t")
if (nrow(eloss)) fwrite(eloss, file.path(out, "eqtl_identity_losses_all.tsv"), sep = "\t")
gwas_loss_files <- setdiff(list.files(file.path(run, "identity_losses"), full.names = TRUE), loss_files)
if (length(gwas_loss_files)) {
  gl <- rbindlist(lapply(gwas_loss_files, fread))
  fwrite(gl, file.path(out, "gwas_identity_losses_all.tsv"), sep = "\t")
  fwrite(gl[, .(rows = .N), by = .(assay, reason)], file.path(out, "gwas_identity_losses_summary.tsv"), sep = "\t")
}
print(gene_tab[, .N, by = new_state]); print(gene_tab[, .N, by = new_state_gt08])
print(gs[, .N, by = .(old_method, new_method)][order(-N)])
cat("missing study-chr runs:", nrow(missing_runs), "\n")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo_compare.txt"))
