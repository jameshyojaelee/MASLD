#!/usr/bin/env Rscript
# 14: the competitive reference, in every assay that carries a statistical test.
#
# For each testable module, the 10,000 size-matched connected subgraphs drawn
# once in step 05 are scored and fitted exactly as the module was, in that
# assay, on the endpoint of the module's family. The question is not "is this
# association real" but "would a coherent set of this size have done as well
# here". It is one family per assay, BH-corrected over the 160 modules, and it
# never decides the association state; it is reported as a separate
# specificity state.
#
# TWO SENSITIVITIES TRAVEL WITH IT. The 2,000 draws closest in within-set
# coherence to the module (coherence measured on the discovery substrate, where
# the sets were drawn) give a coherence-matched p; the calibration of the
# unmatched machinery under label permutation was recorded in 05.
#
# THE SCORING MUST BE IDENTICAL. Each loader below re-derives the assay's
# observed coefficient and REFUSES to continue unless it reproduces the value
# the assay script already wrote, to 1e-8.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))
suppressPackageStartupMessages({ library(edgeR); library(igraph) })

contract <- cam_contract()
out <- cam_dir("nulls")
assays_dir <- file.path(cam_out_root(), "assays")

membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
registry <- fread(file.path(cam_out_root(), "modules", "module_registry.tsv"))
labels <- fread(file.path(cam_out_root(), "discovery", "module_discovery_labels.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)
module_ids <- names(modules)
n_family <- length(module_ids)
family_of <- setNames(labels$endpoint_family, labels$module_id)[module_ids]
size_of <- setNames(registry$n_genes, registry$module_id)[module_ids]
sets_by_size <- readRDS(file.path(out, "competitive_sets_by_size.rds"))
matched_idx <- readRDS(file.path(out, "competitive_coherence_matched_index.rds"))
fam_map <- contract$endpoint_families

py <- "/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
load_npy <- function(path) {
  tsv <- paste0(path, ".tsv")
  if (!file.exists(tsv)) system2(py, c("-c", shQuote(sprintf(
    "import numpy as np; a=np.load('%s'); np.savetxt('%s', a, delimiter='\\t', fmt='%%.6f')",
    path, tsv))))
  as.matrix(fread(tsv))
}
collapse_mean <- function(m, group) {
  rowsum(m, group = group) / as.vector(table(group)[sort(unique(group))])
}

# --- loaders: one definition of each assay's scoring substrate ---------------
# Each returns the gene-by-unit matrix, a named list of z-scored endpoints by
# family, the covariate frame, and the block used for standardisation.
loaders <- list()

loaders$bulk_gse268273 <- function() {
  ga <- fread(cam_input("gse268273_gene_axis", contract))
  pa <- fread(cam_input("gse268273_participant_axis", contract))
  ph <- fread(cam_input("gse268273_phenotype", contract))[match(pa$row_id, row_id)]
  mat <- t(as.matrix(fread(file.path(assays_dir, "gse268273_counts_readback.tsv"))))
  rownames(mat) <- ga$gencode_v49_gene_name; colnames(mat) <- pa$row_id
  expr <- edgeR::cpm(rowsum(mat, group = rownames(mat)), log = TRUE, prior.count = 1)
  covar <- data.frame(sex = ph$sex)
  list(expr = expr, covar = covar, block = rep("all", ncol(expr)),
       y = list(fibrosis = cam_endpoint_z(suppressWarnings(as.numeric(ph$fibrosis_stage)), covar)))
}

loaders$bulk_gse276114 <- function() {
  raw <- fread(cam_input("gse276114_counts", contract)); setnames(raw, 1L, "gene")
  des <- fread(cam_input("gse276114_design", contract))
  m <- as.matrix(raw[, -1L]); rownames(m) <- raw$gene
  des <- des[match(colnames(m), matrix_column_name)]
  keep <- des$disease == "MASLD"
  # CPM on the full 177-column matrix, then subset, exactly as 06 does.
  expr_all <- edgeR::cpm(rowsum(m, group = rownames(m)), log = TRUE, prior.count = 1)
  expr <- expr_all[, keep, drop = FALSE]
  y <- unname(c("F0-2" = 0, "F3" = 1, "F4" = 2)[des$disease_group[keep]])
  list(expr = expr, covar = NULL, block = rep("all", sum(keep)),
       y = list(fibrosis = cam_endpoint_z(y, NULL)))
}

loaders$proteome_liver_pxd051911 <- function() {
  prot <- fread(cam_input("proteome_liver", contract))
  meta <- fread(cam_input("proteome_meta", contract))
  universe <- fread(file.path(cam_out_root(), "universe", "universe.tsv"))$gene_symbol
  vc <- names(prot)[4:ncol(prot)]
  mat <- as.matrix(prot[, ..vc])
  gl <- lapply(strsplit(prot$Genes, ";", fixed = TRUE), function(g) trimws(g[nzchar(trimws(g))]))
  gene <- vapply(gl, function(g) {
    if (length(g) == 1L) return(g); inu <- g[g %in% universe]
    if (length(inu) == 1L) inu else NA_character_ }, character(1))
  keep <- rowMeans(!is.na(mat)) >= contract$universe$proteome_min_detection_fraction & !is.na(gene)
  mat <- log2(mat[keep, , drop = FALSE]); gene <- gene[keep]
  present <- matrix(as.numeric(!is.na(mat)), nrow(mat), ncol(mat))
  expr <- rowsum(replace(mat, is.na(mat), 0), group = gene) / rowsum(present, group = gene)
  expr[!is.finite(expr)] <- NA_real_
  pm <- meta[match(colnames(expr), liver_proteomics_filename)]
  covar <- data.frame(age = suppressWarnings(as.numeric(pm$alder)),
                      bmi = suppressWarnings(as.numeric(pm$bmi)), sex = pm$gender)
  list(expr = expr, covar = covar, block = rep("all", ncol(expr)),
       y = list(activity = cam_endpoint_z(suppressWarnings(as.numeric(pm$nafld_activity_score)), covar),
                fibrosis = cam_endpoint_z(suppressWarnings(as.numeric(sub("^F", "", pm$kleiner_fibrosis_grade))), covar)))
}

loaders$h3k27ac_gse267145 <- function() {
  chrom <- file.path(cam_out_root(), "chromatin")
  axis <- fread(file.path(chrom, "h3k27ac_gene_windows", "gene_two_window_axis.tsv"))
  G <- load_npy(file.path(chrom, "h3k27ac_gene_windows", "gene_two_window_counts.npy"))
  lib <- as.vector(load_npy(file.path(chrom, "h3k27ac_gene_windows", "library_size_original_region_sum.npy")))
  samples <- readLines(file.path(chrom, "h3k27ac_input", "samples.txt"))
  pr <- which(axis$kind == "promoter" & axis$window_observed)
  G <- G[, pr, drop = FALSE]
  gencode <- fread(cam_input("gencode_metadata", contract))
  sym <- gencode$gene_name[match(ml_base_gene_id(axis$gene_id[pr]), gencode$ensembl_base)]
  k <- !is.na(sym); G <- G[, k, drop = FALSE]; sym <- sym[k]
  cpm <- log2(t(G) / rep(lib, each = ncol(G)) * 1e6 + 1)
  rownames(cpm) <- sym
  expr <- collapse_mean(cpm, sym); colnames(expr) <- samples
  join <- fread(cam_input("gse267145_join", contract))
  jm <- join[match(sub("_.*$", "", samples), participant_id)]
  covar <- data.frame(sex = jm$sex)
  nas <- suppressWarnings(as.numeric(jm$steatosis) + as.numeric(jm$ballooning) + as.numeric(jm$lobular_inflammation))
  list(expr = expr, covar = covar, block = rep("all", ncol(expr)),
       y = list(activity = cam_endpoint_z(nas, covar),
                fibrosis = cam_endpoint_z(suppressWarnings(as.numeric(jm$fibrosis)), covar)))
}

loaders$snrna_all_cells <- function() {
  source(file.path(CAM_PROJECT_ROOT, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
  ac <- readRDS(file.path(cam_out_root(), "universe", "snrna_all_cells_pseudobulk.rds"))
  s2d <- build_srr_to_donor_map(CAM_PROJECT_ROOT)
  donor <- ifelse(colnames(ac) %in% names(s2d), s2d[colnames(ac)], colnames(ac))
  dm <- t(rowsum(t(ac), group = donor))
  expr <- edgeR::cpm(dm, log = TRUE, prior.count = 1)
  sn <- unique(fread(cam_input("snrna_metadata", contract))[, .(sample, dataset, disease_stage_coarse)])
  sn[, donor := ifelse(sample %in% names(s2d), s2d[sample], sample)]
  dmeta <- sn[, .(dataset = names(sort(table(dataset), decreasing = TRUE))[1],
                  stage = names(sort(table(disease_stage_coarse), decreasing = TRUE))[1]),
              by = donor][match(colnames(expr), donor)]
  covar <- data.frame(dataset = dmeta$dataset)
  y <- unname(c(Healthy = 0, Steatosis = 1, Steatohepatitis = 2, Cirrhosis = 3)[dmeta$stage])
  list(expr = expr, covar = covar, block = rep("all", ncol(expr)),
       y = list(fibrosis = cam_endpoint_z(y, covar)))
}

loaders$atac_gse296875 <- function() {
  ap <- file.path(cam_out_root(), "chromatin", "atac_promoter")
  axis <- readLines(file.path(ap, "gene_axis.txt"))
  A <- load_npy(file.path(ap, "gene_promoter_counts.npy"))
  lib <- as.vector(load_npy(file.path(ap, "library_size_all_peaks.npy")))
  donors <- fread(file.path(cam_input("gse296875_donor_dir", contract), "donor_axis.tsv"))
  ep <- fread(cam_input("gse296875_endpoints", contract))[match(donors$donor_id, donor_id)]
  expr <- log2(t(A) / rep(lib, each = ncol(A)) * 1e6 + 1)
  rownames(expr) <- axis; colnames(expr) <- as.character(donors$donor_id)
  fib <- as.numeric(ep$fibrosis_any == "true" | ep$fibrosis_any == TRUE)
  fib[ep$fibrosis_observed %in% c("false", FALSE)] <- NA_real_
  # The well is a fixed-effect covariate, as in 10b.
  list(expr = expr, covar = data.frame(well = donors$well_id), block = rep("all", ncol(expr)),
       y = list(fibrosis = cam_endpoint_z(fib)))
}

# --- generic scoring and fitting -------------------------------------------
score_sets_in <- function(expr, block, sets) {
  z <- cam_zscore_rows(expr)
  S <- t(vapply(seq_along(sets), function(i) {
    idx <- intersect(sets[[i]], rownames(z))
    if (!length(idx)) return(rep(NA_real_, ncol(z)))
    colMeans(z[idx, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(expr))))
  for (b in unique(block)) {
    j <- which(block == b)
    S[, j] <- t(apply(S[, j, drop = FALSE], 1L, standardize_vector))
  }
  list(scores = S, measured = rownames(z))
}

fit_sets <- function(S, y, covar) {
  ok <- is.finite(y)
  if (!is.null(covar)) ok <- ok & stats::complete.cases(covar)
  X <- NULL
  if (!is.null(covar)) {
    cv <- covar[ok, , drop = FALSE]
    cv <- cv[, vapply(cv, function(v) length(unique(v)) > 1L, logical(1)), drop = FALSE]
    if (ncol(cv) > 0L) X <- stats::model.matrix(~ ., data = cv)[, -1, drop = FALSE]
  }
  Sf <- S[, ok, drop = FALSE]
  bad <- !stats::complete.cases(Sf)
  Sf[bad, ] <- 0
  fit <- cam_fast_assoc(Sf, y[ok], X)
  fit$beta[bad] <- NA_real_; fit$se[bad] <- NA_real_; fit$p_two_sided[bad] <- NA_real_
  fit
}

observed_files <- list(
  bulk_gse268273 = "external_bulk_results.tsv", bulk_gse276114 = "external_bulk_results.tsv",
  proteome_liver_pxd051911 = "proteome_results.tsv", h3k27ac_gse267145 = "h3k27ac_results.tsv",
  snrna_all_cells = "snrna_results.tsv", atac_gse296875 = "atac_results.tsv")

null_rows <- list(); check_rows <- list(); diag_rows <- list()
for (assay in names(loaders)) {
  cam_say("competitive null: ", assay)
  L <- loaders[[assay]]()
  families <- names(L$y)
  obs_sc <- score_sets_in(L$expr, L$block, modules)
  obs_fit <- lapply(families, function(f) fit_sets(obs_sc$scores, L$y[[f]], L$covar))
  names(obs_fit) <- families

  # Module rows for this assay: the family endpoint, testability from the
  # assay table, and the guard that the loader reproduces the reported beta.
  this_assay <- assay
  saved <- fread(file.path(assays_dir, observed_files[[assay]]))
  saved <- saved[saved$assay == this_assay]
  endpoint_of <- vapply(module_ids, function(mid) {
    e <- fam_map[[family_of[[mid]]]][[assay]]; if (is.null(e)) NA_character_ else e }, character(1))
  saved_beta <- saved_testable <- rep(NA, n_family)
  for (i in seq_len(n_family)) {
    if (is.na(endpoint_of[i])) next
    r <- saved[module_id == module_ids[i] & endpoint == endpoint_of[i]]
    if (nrow(r) == 1L) { saved_beta[i] <- r$beta; saved_testable[i] <- r$testable }
  }
  recomputed <- vapply(seq_len(n_family), function(i) {
    if (is.na(endpoint_of[i])) return(NA_real_)
    obs_fit[[family_of[[i]]]]$beta[i] }, numeric(1))
  cmp <- data.table(module_id = module_ids, recomputed = recomputed, saved = as.numeric(saved_beta))
  cmp <- cmp[is.finite(recomputed) & is.finite(saved)]
  worst <- if (nrow(cmp)) max(abs(cmp$recomputed - cmp$saved)) else NA_real_
  check_rows[[length(check_rows) + 1L]] <- data.table(
    assay = assay, n_compared = nrow(cmp), max_abs_difference = worst)
  cam_assert(is.finite(worst) && worst < 1e-8, paste0(
    "Null loader for ", assay, " does not reproduce the reported observed beta (max diff ",
    signif(worst, 3), "); the null would not match the statistic"))
  cam_say("  loader reproduces observed beta to ", signif(worst, 3), " over ", nrow(cmp), " modules")

  for (s in sort(unique(size_of))) {
    sets <- sets_by_size[[as.character(s)]]
    idx <- which(size_of == s & isTRUE_vec(saved_testable) & !is.na(endpoint_of))
    if (!length(sets) || !length(idx)) next
    dsc <- score_sets_in(L$expr, L$block, sets)
    nullz_by_family <- lapply(families, function(f) {
      d <- fit_sets(dsc$scores, L$y[[f]], L$covar); abs(d$beta / d$se) })
    names(nullz_by_family) <- families
    for (i in idx) {
      f <- family_of[[i]]
      of <- obs_fit[[f]]
      obs_z <- abs(of$beta[i] / of$se[i])
      nullz <- nullz_by_family[[f]]
      mi <- matched_idx[[module_ids[i]]]
      null_rows[[length(null_rows) + 1L]] <- data.table(
        assay = assay, module_id = module_ids[i], endpoint = endpoint_of[i], endpoint_family = f,
        observed_abs_z = obs_z,
        competitive_null_mean_abs_z = mean(nullz, na.rm = TRUE),
        competitive_null_q95_abs_z = unname(stats::quantile(nullz, 0.95, na.rm = TRUE)),
        competitive_exceedance_count = sum(nullz >= obs_z, na.rm = TRUE),
        competitive_n_draws = sum(is.finite(nullz)),
        competitive_connected_p = cam_empirical_p_ge(obs_z, nullz),
        coherence_matched_n_draws = sum(is.finite(nullz[mi$idx])),
        coherence_matched_p = cam_empirical_p_ge(obs_z, nullz[mi$idx]),
        coherence_gap = mi$gap)
    }
    diag_rows[[length(diag_rows) + 1L]] <- data.table(
      assay = assay, module_size = s, n_drawn = length(sets), n_modules_tested = length(idx))
  }
}

nulls <- rbindlist(null_rows)
nulls[, competitive_connected_q := ml_complete_bh(competitive_connected_p, n_family), by = assay]
nulls[, coherence_matched_q := ml_complete_bh(coherence_matched_p, n_family), by = assay]
cam_assert_no_prohibited_columns(nulls, contract)
cam_write_tsv(nulls, file.path(out, "competitive_nulls_per_assay.tsv"))
cam_write_tsv(rbindlist(check_rows), file.path(out, "loader_reproduction_check.tsv"))
cam_write_tsv(rbindlist(diag_rows), file.path(out, "competitive_draw_diagnostics.tsv"))
by_assay <- nulls[, .(
  n_tested = .N,
  n_p05_uncorrected_descriptive = sum(competitive_connected_p < 0.05, na.rm = TRUE),
  expected_at_0p05_if_independent_descriptive = 0.05 * .N,
  n_bh_q05 = sum(competitive_connected_q < 0.05, na.rm = TRUE),
  n_coherence_matched_bh_q05 = sum(coherence_matched_q < 0.05, na.rm = TRUE),
  median_coherence_gap = stats::median(coherence_gap, na.rm = TRUE)), by = assay]
cam_write_tsv(by_assay, file.path(out, "competitive_counts_by_assay.tsv"))
cam_write_json(list(
  assays = names(loaders), draws_per_size = length(sets_by_size[[1]]),
  n_sizes = length(sets_by_size),
  max_loader_difference = max(rbindlist(check_rows)$max_abs_difference),
  by_assay = by_assay,
  note = "counts at uncorrected p<0.05 are descriptive; the state uses the BH-corrected q"
), file.path(out, "competitive_null_summary.json"))
writeLines("competitive nulls done", file.path(out, "READY"))
cam_say("14 complete")
