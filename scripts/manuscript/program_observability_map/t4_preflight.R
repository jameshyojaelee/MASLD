#!/usr/bin/env Rscript

# T4 preflight: build the real blocks, time one permutation replicate, and print
# the overlap geometry that decides how large the testable family will be.
# Nothing is written and nothing is sealed. Its only job is to size the run and
# to surface a structural problem before a compute job is queued.

suppressPackageStartupMessages({
  library(edgeR); library(data.table); library(Matrix)
})
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "t4_interaction_lib.R"))
set.seed(seed)

dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()

registry <- fread(program_registry); membership <- fread(program_membership)
scored <- score_programs(symbol_matrix, meta, membership, registry, program_weight_coverage)
features <- scored$coverage[testable == TRUE, feature_id]
cat("testable programs:", length(features), "\n")

discovery_meta <- meta[dataset %in% discovery_cohorts & !is.na(fibrosis_stage) & !is.na(nas_score)]
weights <- t4_membership_weights(membership, features)
panel <- intersect(sort(unique(weights$gene_symbol)), rownames(symbol_matrix))
cat("program gene panel:", length(panel), "of", uniqueN(weights$gene_symbol), "\n")
panel_matrix <- symbol_matrix[panel, discovery_meta$sample_id, drop = FALSE]
rm(symbol_matrix); gc()

t0 <- Sys.time()
blocks <- lapply(discovery_cohorts, function(cc)
  t4_build_block(panel_matrix, discovery_meta[dataset == cc], weights, features,
                 program_weight_coverage))
names(blocks) <- discovery_cohorts
cat("blocks built in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
print(rbindlist(lapply(blocks, function(b) data.table(
  cohort = b$cohort, n = b$n, n_par = b$n_par, df = b$df, genes = length(b$genes),
  scored = sum(b$coverage$cohort_testable)))))

pairs <- t4_pair_table(features, weights)
cat("pairs:", nrow(pairs), " sharing >=1 gene:", pairs[n_shared > 0L, .N], "\n")
print(summary(pairs$jaccard))
print(summary(pairs[n_shared > 0L, weight_cosine]))
cat("pairs with weight cosine > 0.2:", pairs[weight_cosine > 0.2, .N], "\n")

t0 <- Sys.time()
ops <- lapply(blocks, t4_disjoint_operators, pairs = pairs, features = features,
              min_retained_weight = program_weight_coverage)
cat("operators built in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
elig <- Reduce(`&`, lapply(discovery_cohorts, function(cc) {
  e <- rep(TRUE, nrow(pairs)); e[ops[[cc]]$pair_index] <- ops[[cc]]$eligible; e }))
cat("membership-entangled pairs at 0.80 retention:", sum(!elig),
    " -> family size", sum(elig), "\n")

t0 <- Sys.time()
obs <- t4_one_pass(blocks, pairs, ops, want_matrix = TRUE)
cat("observed pass:", round(as.numeric(Sys.time() - t0, units = "secs"), 2), "s\n")
cat("observed |r_full| quantiles:\n"); print(quantile(abs(tanh(obs$z_full)), c(.5,.9,.99,1), na.rm=TRUE))
cat("observed globals:\n"); print(obs$globals)

leak <- t4_stage_leakage(blocks)
print(leak[, .(max_raw_r2 = max(raw_stage_r2), mean_r2 = mean(stage_r2),
               chance = chance_stage_r2[[1L]]), by = cohort])

t0 <- Sys.time()
for (b in 1:5) {
  Zl <- lapply(blocks, function(bl) t4_permute_genes(bl$Z))
  t4_one_pass(blocks, pairs, ops, Zlist = Zl)
}
per_rep <- as.numeric(Sys.time() - t0, units = "secs") / 5
cat("per replicate:", round(per_rep, 3), "s\n")
cat("7000 replicates on 16 cores:", round(7000 * per_rep / 16 / 60, 1), "min\n")
cat("PREFLIGHT_OK\n")
