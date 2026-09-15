#!/usr/bin/env Rscript
# 06: the two external bulk cohorts.
#
# These are the only assays that ask the same question in the same modality as
# discovery, so they are the cleanest test of whether a module's association
# is a property of liver biology or of the five cohorts it was built on.
# Neither shares an accession with a discovery cohort.
#
# GSE276114 is a mixed-etiology cohort: 81 participants with MASLD, 82 with
# chronic viral hepatitis, 14 with alcohol-related disease. The primary test
# uses the MASLD participants only. The all-etiology fit with an etiology term
# is kept as a secondary, and both are reported.
#
# v2: scores carry no direction (direction enters once, in 15); the endpoint
# is z-scored so the coefficient is a standardized slope; the fit's residual
# degrees of freedom travel with the row for the equivalence bound.

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/analysis/cross_assay_modules/lib_cross_assay_modules.R"))
suppressPackageStartupMessages(library(edgeR))

contract <- cam_contract()
out <- cam_dir("assays")
cam_assert(file.exists(file.path(cam_out_root(), "discovery", "READY")), "Run 05 first")

membership <- fread(file.path(cam_out_root(), "modules", "module_membership.tsv"))
modules <- split(membership$gene_symbol, membership$module_id)
module_ids <- names(modules)
n_family <- length(module_ids)

# Shared scoring: gene z within the cohort, equal-weight mean over measured
# members, standardised across participants. No direction.
score_matrix <- function(expr, members_list) {
  z <- cam_zscore_rows(expr)
  S <- t(vapply(seq_along(members_list), function(i) {
    idx <- intersect(members_list[[i]], rownames(z))
    if (!length(idx)) return(rep(NA_real_, ncol(z)))
    standardize_vector(colMeans(z[idx, , drop = FALSE]))
  }, numeric(ncol(z))))
  rownames(S) <- names(members_list)
  list(scores = S, measured = rownames(z))
}

results <- list()
row_of <- function(assay, i, endpoint, family, note, tb, fit) data.table(
  assay = assay, module_id = module_ids[i], endpoint = endpoint,
  endpoint_family = family, endpoint_note = note, unit = "participant",
  n_units = fit$n, df = fit$df, testable = tb$testable,
  n_members = tb$n_members, n_measured = tb$n_measured,
  fraction_measured = tb$fraction_measured,
  beta = fit$beta, se = fit$se, p_two_sided = fit$p_two_sided)

# --- GSE268273 --------------------------------------------------------------
cam_say("GSE268273")
npy <- cam_input("gse268273_counts", contract)
gene_axis <- fread(cam_input("gse268273_gene_axis", contract))
part_axis <- fread(cam_input("gse268273_participant_axis", contract))
pheno <- fread(cam_input("gse268273_phenotype", contract))

tmp_tsv <- file.path(out, "gse268273_counts_readback.tsv")
if (!file.exists(tmp_tsv)) {
  py <- "/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
  code <- sprintf(
    "import numpy as np; a=np.load('%s'); np.savetxt('%s', a, delimiter='\\t', fmt='%%.6f')",
    npy, tmp_tsv)
  system2(py, c("-c", shQuote(code)))
}
mat <- as.matrix(fread(tmp_tsv))
cam_assert(nrow(mat) == nrow(part_axis) && ncol(mat) == nrow(gene_axis),
           paste0("GSE268273 matrix is ", nrow(mat), "x", ncol(mat),
                  " but axes say ", nrow(part_axis), "x", nrow(gene_axis)))
expr268 <- t(mat)
rownames(expr268) <- gene_axis$gencode_v49_gene_name
colnames(expr268) <- part_axis$row_id
expr268 <- rowsum(expr268, group = rownames(expr268))
expr268 <- edgeR::cpm(expr268, log = TRUE, prior.count = 1)

pheno268 <- pheno[match(part_axis$row_id, row_id)]
cam_assert(sum(!is.na(pheno268$row_id)) == nrow(part_axis),
           "GSE268273 phenotype does not cover every delivered participant")
fib268 <- suppressWarnings(as.numeric(pheno268$fibrosis_stage))
sex268 <- pheno268$sex
y268 <- cam_endpoint_z(fib268, data.frame(sex = sex268))
cam_say("  participants ", ncol(expr268), ", fibrosis recorded ", sum(is.finite(y268)))

sc <- score_matrix(expr268, modules)
for (i in seq_len(n_family)) {
  tb <- cam_testability(modules[[i]], sc$measured, contract)
  d <- data.table(score = sc$scores[i, ], y = y268, sex = sex268)
  fit <- if (tb$testable) cam_fit_score(d, "score", "y", "sex") else cam_fit_empty()
  results[[length(results) + 1L]] <- row_of("bulk_gse268273", i, "fibrosis_stage_0_4",
    "fibrosis", "fibrosis stage 0-4, z-scored", tb, fit)
}

# --- GSE276114 --------------------------------------------------------------
cam_say("GSE276114")
raw <- fread(cam_input("gse276114_counts", contract))
setnames(raw, 1L, "gene")
design <- fread(cam_input("gse276114_design", contract))
m276 <- as.matrix(raw[, -1L]); rownames(m276) <- raw$gene
design <- design[match(colnames(m276), matrix_column_name)]
cam_assert(sum(is.na(design$matrix_column_name)) == 0L,
           "GSE276114 design does not cover every count column")
m276 <- rowsum(m276, group = rownames(m276))
expr276 <- edgeR::cpm(m276, log = TRUE, prior.count = 1)
stage276 <- c("F0-2" = 0, "F3" = 1, "F4" = 2)[design$disease_group]
cam_say("  participants ", ncol(expr276), "; MASLD ", sum(design$disease == "MASLD"))

for (arm in c("masld_only", "all_etiology")) {
  keep <- if (arm == "masld_only") design$disease == "MASLD" else rep(TRUE, nrow(design))
  sub <- expr276[, keep, drop = FALSE]
  sc2 <- score_matrix(sub, modules)
  covar <- if (arm == "all_etiology") data.frame(etiology = design$disease[keep]) else NULL
  yz <- cam_endpoint_z(unname(stage276[keep]), covar)
  for (i in seq_len(n_family)) {
    tb <- cam_testability(modules[[i]], sc2$measured, contract)
    d <- data.table(score = sc2$scores[i, ], y = yz, etiology = design$disease[keep])
    covars <- if (arm == "all_etiology") "etiology" else character(0)
    fit <- if (tb$testable) cam_fit_score(d, "score", "y", covars) else cam_fit_empty()
    results[[length(results) + 1L]] <- row_of(
      if (arm == "masld_only") "bulk_gse276114" else "bulk_gse276114_all_etiology", i,
      "stage_group_F0_2_F3_F4", "fibrosis",
      "grouped fibrosis F0-2 < F3 < F4, z-scored; a coarser instrument than stage 0-4", tb, fit)
  }
}

res <- rbindlist(results)
cam_assert_no_prohibited_columns(res, contract)
cam_write_tsv(res, file.path(out, "external_bulk_results.tsv"))
cam_write_json(list(
  gse268273_participants = ncol(expr268), gse268273_genes = nrow(expr268),
  gse276114_participants = ncol(expr276),
  gse276114_masld_participants = sum(design$disease == "MASLD"),
  gse276114_etiology_counts = as.list(table(design$disease)),
  modules_testable_gse268273 = sum(res[assay == "bulk_gse268273", testable]),
  modules_testable_gse276114 = sum(res[assay == "bulk_gse276114", testable]),
  direction_applied_at_scoring = FALSE
), file.path(out, "external_bulk_summary.json"))
writeLines("external bulk scored", file.path(out, "READY_external_bulk"))
cam_say("06 complete")
