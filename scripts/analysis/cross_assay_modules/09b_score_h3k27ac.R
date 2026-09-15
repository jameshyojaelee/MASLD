#!/usr/bin/env Rscript
# 09b: score the modules on promoter H3K27ac, and on the matched RNA.
#
# GSE267145 is the only cohort in this Resource where two assays are measured
# on the same participants alongside a severity grade, so it carries the
# panel's only within-person transfer test. This script produces both scores;
# step 11 correlates them.
#
# ENDPOINTS. The NAS component sum is the activity-family endpoint. Recorded
# fibrosis is 71 zeros out of 99 and tops out at stage 3; it is fitted as the
# fibrosis-family endpoint and its cells will mostly print an upper bound.
#
# NORMALISATION. Promoter window counts are divided by the sample's ORIGINAL
# region library size, not by the window sum.
#
# v2: no direction at scoring; endpoints z-scored; df carried.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))

contract <- cam_contract()
out <- cam_dir("assays")
chrom <- file.path(cam_out_root(), "chromatin")
cam_assert(dir.exists(file.path(chrom, "h3k27ac_gene_windows")), "Run 09a first")

membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)
module_ids <- names(modules)

py <- "/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
load_npy <- function(path) {
  tsv <- paste0(path, ".tsv")
  if (!file.exists(tsv)) {
    system2(py, c("-c", shQuote(sprintf(
      "import numpy as np; a=np.load('%s'); np.savetxt('%s', a, delimiter='\\t', fmt='%%.6f')",
      path, tsv))))
  }
  as.matrix(fread(tsv))
}

axis <- fread(file.path(chrom, "h3k27ac_gene_windows", "gene_two_window_axis.tsv"))
G <- load_npy(file.path(chrom, "h3k27ac_gene_windows", "gene_two_window_counts.npy"))
lib <- as.vector(load_npy(file.path(chrom, "h3k27ac_gene_windows",
                                    "library_size_original_region_sum.npy")))
samples <- readLines(file.path(chrom, "h3k27ac_input", "samples.txt"))
cam_assert(nrow(G) == length(samples), "Chromatin matrix rows and sample list disagree")

promoter <- which(axis$kind == "promoter" & axis$window_observed)
G <- G[, promoter, drop = FALSE]
gene_ids <- axis$gene_id[promoter]
gencode <- fread(cam_input("gencode_metadata", contract))
sym <- gencode$gene_name[match(ml_base_gene_id(gene_ids), gencode$ensembl_base)]
keep <- !is.na(sym)
G <- G[, keep, drop = FALSE]; sym <- sym[keep]

cpm <- log2(t(G) / rep(lib, each = ncol(G)) * 1e6 + 1)
rownames(cpm) <- sym
cpm <- rowsum(cpm, group = rownames(cpm)) / as.vector(table(sym)[sort(unique(sym))])
colnames(cpm) <- samples
cam_say("promoter H3K27ac: ", nrow(cpm), " genes x ", ncol(cpm), " samples")

join <- fread(cam_input("gse267145_join", contract))
participant <- sub("_.*$", "", samples)
jm <- join[match(participant, participant_id)]
cam_assert(sum(is.na(jm$participant_id)) == 0L,
           "Some H3K27ac samples do not join to a participant")
jm[, nas_component_sum := suppressWarnings(as.numeric(steatosis) + as.numeric(ballooning) +
                                             as.numeric(lobular_inflammation))]
jm[, fibrosis_numeric := suppressWarnings(as.numeric(fibrosis))]
cam_say("NAS component sum range ", paste(range(jm$nas_component_sum, na.rm = TRUE), collapse = "-"),
        "; fibrosis table ", paste(names(table(jm$fibrosis_numeric)),
                                   table(jm$fibrosis_numeric), collapse = " "))

# Matched RNA from the same participants, for the paired test in step 11.
rna_raw <- fread(cam_input("gse267145_rna", contract))
setnames(rna_raw, 1L, "feature")
rna_mat <- as.matrix(rna_raw[, -1L]); rownames(rna_mat) <- rna_raw$feature
rna_part <- sub("_.*$", "", colnames(rna_mat))
shared <- intersect(rna_part, participant)
rna_mat <- rna_mat[, rna_part %in% shared, drop = FALSE]
rna_sym <- gencode$gene_name[match(ml_base_gene_id(rownames(rna_mat)), gencode$ensembl_base)]
rk <- !is.na(rna_sym)
rna_mat <- rna_mat[rk, , drop = FALSE]; rna_sym <- rna_sym[rk]
rna_cpm <- log2(t(t(rna_mat) / colSums(rna_mat)) * 1e6 + 1)
rownames(rna_cpm) <- rna_sym
rna_cpm <- rowsum(rna_cpm, group = rownames(rna_cpm)) /
  as.vector(table(rna_sym)[sort(unique(rna_sym))])
cam_say("matched RNA: ", nrow(rna_cpm), " genes x ", ncol(rna_cpm), " samples")

score_all <- function(mat) {
  z <- cam_zscore_rows(mat)
  S <- t(vapply(seq_along(module_ids), function(i) {
    idx <- intersect(modules[[i]], rownames(z))
    if (!length(idx)) return(rep(NA_real_, ncol(z)))
    standardize_vector(colMeans(z[idx, , drop = FALSE], na.rm = TRUE))
  }, numeric(ncol(mat))))
  rownames(S) <- module_ids
  list(scores = S, measured = rownames(z))
}
chrom_sc <- score_all(cpm)
rna_sc <- score_all(rna_cpm)
colnames(chrom_sc$scores) <- participant
colnames(rna_sc$scores) <- sub("_.*$", "", colnames(rna_cpm))
cam_write_rds(chrom_sc$scores, file.path(chrom, "h3k27ac_module_scores.rds"))
cam_write_rds(rna_sc$scores, file.path(chrom, "gse267145_rna_module_scores.rds"))
# Coverage on each side, so step 11 can apply the testability rule to both.
writeLines(chrom_sc$measured, file.path(chrom, "h3k27ac_measured_genes.txt"))
writeLines(rna_sc$measured, file.path(chrom, "gse267145_rna_measured_genes.txt"))
# Gene-level z restricted to module members, for the shared-member sensitivity
# in step 11 (both scores recomputed on the genes measured in BOTH assays).
mod_genes <- unique(membership$gene_symbol)
zc <- cam_zscore_rows(cpm); zc <- zc[rownames(zc) %in% mod_genes, , drop = FALSE]
colnames(zc) <- participant
zr <- cam_zscore_rows(rna_cpm); zr <- zr[rownames(zr) %in% mod_genes, , drop = FALSE]
colnames(zr) <- sub("_.*$", "", colnames(rna_cpm))
cam_write_rds(zc, file.path(chrom, "h3k27ac_module_gene_z.rds"))
cam_write_rds(zr, file.path(chrom, "gse267145_rna_module_gene_z.rds"))

endpoints <- list(
  nas_component_sum = list(family = "activity", note = "steatosis+ballooning+lobular inflammation, z-scored"),
  fibrosis_numeric = list(family = "fibrosis", note = "fibrosis 0-3, 71 of 99 at zero, z-scored"))
rows <- list()
for (ep in names(endpoints)) {
  yz <- cam_endpoint_z(jm[[ep]], data.frame(sex = jm$sex))
  for (i in seq_along(module_ids)) {
    tb <- cam_testability(modules[[i]], chrom_sc$measured, contract)
    d <- data.table(score = chrom_sc$scores[i, ], y = yz, sex = jm$sex)
    fit <- if (tb$testable) cam_fit_score(d, "score", "y", "sex") else cam_fit_empty()
    rows[[length(rows) + 1L]] <- data.table(
      assay = "h3k27ac_gse267145", module_id = module_ids[i],
      endpoint = ep, endpoint_family = endpoints[[ep]]$family, endpoint_note = endpoints[[ep]]$note,
      unit = "participant", n_units = fit$n, df = fit$df,
      testable = tb$testable, n_members = tb$n_members,
      n_measured = tb$n_measured, fraction_measured = tb$fraction_measured,
      beta = fit$beta, se = fit$se, p_two_sided = fit$p_two_sided)
  }
}
res <- rbindlist(rows)
cam_assert_no_prohibited_columns(res, contract)
cam_write_tsv(res, file.path(out, "h3k27ac_results.tsv"))
cam_write_json(list(
  participants = ncol(cpm), promoter_genes = nrow(cpm),
  matched_rna_participants = ncol(rna_cpm), matched_rna_genes = nrow(rna_cpm),
  modules_testable = sum(res[endpoint == "nas_component_sum", testable]),
  modules_testable_rna_side = sum(vapply(modules, function(M)
    cam_testability(M, rna_sc$measured, contract)$testable, logical(1))),
  fibrosis_distribution = as.list(table(jm$fibrosis_numeric)),
  promoter_bp = 1000, direction_applied_at_scoring = FALSE
), file.path(out, "h3k27ac_summary.json"))
writeLines("h3k27ac scored", file.path(out, "READY_h3k27ac"))
cam_say("09b complete")
