#!/usr/bin/env Rscript
# H54 plan: count programs testable and coordinated edges testable in the two new
# cohorts under the exact T4 scoring path, with genes restricted to the frozen
# discovery gene universe (rows of holdout_dge.rds). No edge statistic is computed.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(Matrix)})
script_dir <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/program_observability_map"
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "t4_interaction_lib.R"))
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
T4 <- normalizePath(file.path(sealed_input_root, "..", "systems", "T4_interaction"))
communities <- fread(file.path(T4, "program_communities.tsv"))
coord <- fread(file.path(T4, "coordinated_edges.tsv"))
membership <- fread(program_membership)
registry <- fread(program_registry)
features <- communities$feature_id
weights <- t4_membership_weights(membership, features)
universe <- sub("\\..*$", "", rownames(readRDS(holdout_dge)))
cat("frozen gene universe rows:", length(universe), "\n")
lin <- registry$cell_type[match(features, registry$program_uid)]
cat("programs per lineage:\n"); print(table(lin, useNA = "ifany"))
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))

run <- function(label, counts_path, meta) {
  cts <- fread(counts_path, header = TRUE)
  ids <- sub("\\..*$", "", cts[[1]])
  M <- as.matrix(cts[, -1, with = FALSE]); rownames(M) <- ids
  M <- M[, meta$matrix_column, drop = FALSE]
  colnames(M) <- meta$sample_id
  keep <- ids %in% universe
  cat(sprintf("%s: rows %d, in universe %d, counts share in universe median %.4f\n",
              label, length(ids), sum(keep), median(colSums(M[keep, ]) / colSums(M))))
  dge <- calcNormFactors(DGEList(M[keep, , drop = FALSE]), method = "TMM")
  lib <- cpm(dge, log = TRUE, prior.count = 1)
  sm <- collapse_symbols(lib, rownames(lib), annotation)
  panel <- intersect(sort(unique(weights$gene_symbol)), rownames(sm))
  blk <- t4_build_block(sm[panel, meta$sample_id, drop = FALSE], meta, weights, features,
                        program_weight_coverage)
  cv <- blk$coverage
  testable <- cv[cohort_testable == TRUE, feature_id]
  ok_edge <- coord$program_a %in% testable & coord$program_b %in% testable
  dropped <- cv[cohort_testable == FALSE]
  cat(sprintf("%s: n=%d, design cols %d, residual df %d, panel symbols %d of %d, variable %d; programs testable %d of %d; min coverage %.3f; coordinated edges with both programs testable %d of %d\n",
              label, blk$n, blk$n_par, blk$df, length(panel), uniqueN(weights$gene_symbol),
              length(blk$genes), length(testable), length(features), min(cv$cohort_weight_coverage),
              sum(ok_edge), nrow(coord)))
  if (nrow(dropped)) {
    dropped[, community := communities$community[match(feature_id, communities$feature_id)]]
    dropped[, lineage := lin[match(feature_id, features)]]
    print(dropped[, .(feature_id, lineage, community, cohort_weight_coverage)])
  }
  invisible(NULL)
}

a <- fread(file.path(root, "data/GSE281797/metadata/sample_phenotypes_resolved.tsv"))
a91 <- a[sex_concordant == TRUE]
mk <- function(d, sexcol) d[, .(sample_id = as.character(sample_id), matrix_column = as.character(matrix_column),
                                 dataset = "new", fibrosis_stage = as.integer(fibrosis),
                                 nas_score = as.numeric(nas), sex_final = get(sexcol))]
run("GSE281797 all-stage n91", file.path(root, "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz"), mk(a91, "sex_expressed"))
run("GSE281797 F0-F1 n66", file.path(root, "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz"), mk(a91[fibrosis <= 1], "sex_expressed"))
run("GSE281797 all n94 (deposited sex)", file.path(root, "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz"), mk(a, "sex_deposited"))
r <- fread(file.path(root, "data/GSE225740/metadata/sample_phenotypes_resolved.tsv"))
r[, sample_id := paste0("pt", pt_id)]
r87 <- r[!is.na(fibrosis) & !is.na(nas)]
run("GSE225740 all-stage n87", file.path(root, "data/GSE225740/GSE225740_counts_ensembl_stable.tsv.gz"), mk(r87, "sex_expressed"))
run("GSE225740 F0-F1 n53", file.path(root, "data/GSE225740/GSE225740_counts_ensembl_stable.tsv.gz"), mk(r87[fibrosis <= 1], "sex_expressed"))
