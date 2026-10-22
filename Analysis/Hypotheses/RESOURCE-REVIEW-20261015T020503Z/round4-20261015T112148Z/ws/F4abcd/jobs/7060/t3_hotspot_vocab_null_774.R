#!/usr/bin/env Rscript
# Review R10: re-derive the Hotspot (unweighted) matched-random vocabulary null
# quoted in docs/manuscript/draft/fig4.md:193-195 ("fibrosis 26/114 versus null
# mean 18.445, nominal P=0.134; NAS 4 versus 4.225, P=0.512").
#
# Arm A reproduces the stored August T3 value on the 469-record working roster
# (positive control for the port). Arms B and C apply the identical code to the
# current Figure 4E rosters: B = 774 units, three cohorts (364 units);
# C = 774 units plus the unchanged GSE174478 records (457), the arm installed in
# Figure 4E. Only the hotspot_117 vocabulary is run. Copied library code is
# unchanged. Seeds follow t3_run_vocabulary.R (t3_seed 20260812; hotspot is the
# first vocabulary, so draw seed = 20260812 + 100000, then + feature index).

suppressPackageStartupMessages({
  library(edgeR); library(data.table); library(Matrix); library(matrixStats); library(metafor)
})
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
Sys.setenv(MASLD_PROJECT_ROOT = repo)
lib <- file.path(getwd(), "lib")
source(file.path(lib, "config.R"))
source(file.path(lib, "analysis_lib.R"))
source(file.path(lib, "t3_vocabulary_lib.R"))
out_dir <- file.path(getwd(), "out")
mc <- 4L
n_random <- 200L
draw_seed <- 20260812L + 100000L * 1L
set.seed(20260812L)

ann <- fread(gene_annotation, select = c("gene_id", "gene_name"))
hs <- fread(program_membership)
hs <- hs[!is.na(mapped_symbol) & mapped_symbol != ""]
vocab <- unique(data.table(feature_id = hs$program_uid,
                           gene_symbol = toupper(trimws(hs$mapped_symbol)),
                           direction = "up"), by = c("feature_id", "gene_symbol"))
vocab <- vocab[!is.na(gene_symbol) & gene_symbol != ""]

to_symbols <- function(dge) {
  lc <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  collapse_symbols(lc, rownames(lc), ann)
}
nonconstant <- function(sm, meta, cohorts, univ) {
  for (co in cohorts) {
    j <- which(meta$dataset == co)
    s <- matrixStats::rowSds(sm[univ, j, drop = FALSE])
    univ <- univ[is.finite(s) & s > 0]
  }
  univ
}
run_arm <- function(label, S_obs, designs, zlist, universe, bins, testable) {
  t0 <- Sys.time()
  fit_obs <- t3_fit_vocabulary(S_obs, designs, mc_cores = mc)
  obs <- fit_obs[, .(n_supported = sum(q_value < 0.05, na.rm = TRUE),
                     n_estimable = sum(is.finite(p_value_meta))), by = axis]
  members <- split(vocab[feature_id %in% testable, gene_symbol],
                   vocab[feature_id %in% testable, feature_id])
  members <- lapply(members, function(g) intersect(g, universe))
  draws <- t3_matched_random_sets(members, universe, bins, n_random, seed = draw_seed)
  flat <- unlist(draws, recursive = FALSE, use.names = FALSE)
  names(flat) <- paste(rep(names(draws), each = n_random),
                       rep(seq_len(n_random), times = length(draws)), sep = "||")
  sid <- unlist(lapply(designs, `[[`, "sample_id"), use.names = FALSE)
  S_rand <- t3_score_sets_fast(zlist, flat)[, sid, drop = FALSE]
  fit_rand <- t3_fit_vocabulary(S_rand, designs, mc_cores = mc)
  fit_rand[, replicate := as.integer(sub("^.*\\|\\|", "", feature_id))]
  fit_rand[, q_within := p.adjust(p_value_meta, method = "BH", n = .N), by = replicate]
  counts <- fit_rand[, .(n_supported = sum(q_within < 0.05, na.rm = TRUE),
                         n_tested = uniqueN(feature_id)), by = .(replicate, axis)]
  res <- rbindlist(lapply(c("fibrosis", "nas"), function(ax) {
    o <- obs[axis == ax, n_supported]
    nv <- counts[axis == ax, n_supported]
    data.table(arm = label, axis = ax, n_units = sum(vapply(designs, `[[`, integer(1), "n")),
               n_cohorts = length(designs), n_testable = length(testable),
               bh_family = nrow(fit_obs), observed_supported = o,
               null_mean = mean(nv), null_q025 = quantile(nv, 0.025, names = FALSE),
               null_q975 = quantile(nv, 0.975, names = FALSE), null_max = max(nv),
               n_draws = length(nv), p_vs_matched_random = t3_empirical_p_ge(o, nv),
               resolution_floor = 1 / (length(nv) + 1), universe = length(universe),
               seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }))
  fwrite(res, file.path(out_dir, paste0("summary_", label, ".tsv")), sep = "\t")
  fwrite(fit_obs, file.path(out_dir, paste0("observed_fit_", label, ".tsv")), sep = "\t")
  fwrite(counts, file.path(out_dir, paste0("null_counts_", label, ".tsv")), sep = "\t")
  print(res)
  rm(S_rand, fit_rand); gc()
  res
}

# ---- Arm A: stored 469-record working roster (positive control) ----------
message("[A] loading old substrate")
old_meta <- as.data.table(readRDS(nonholdout_meta))
old_dge <- readRDS(nonholdout_dge)
stopifnot(identical(colnames(old_dge), old_meta$sample_id))
sym_old <- to_symbols(old_dge); rm(old_dge); gc()
discA <- old_meta[dataset %in% discovery_cohorts & !is.na(fibrosis_stage) & !is.na(nas_score)]
stopifnot(nrow(discA) == 469L)
univA <- nonconstant(sym_old, old_meta, discovery_cohorts, rownames(sym_old))
sgA <- score_gene_sets(sym_old, old_meta, vocab, gene_set_min_genes, gene_set_coverage, FALSE)
testA <- sgA$coverage[testable == TRUE, feature_id]
message("[A] testable ", length(testA), " universe ", length(univA))
resA <- run_arm("A_old469_reproduction",
                sgA$scores[testA, discA$sample_id, drop = FALSE],
                t3_designs(discA, discovery_cohorts),
                t3_cohort_z(sym_old, old_meta, discovery_cohorts, univA),
                univA, t3_expression_deciles(sym_old), testA)

# ---- New 774-unit substrate ----------------------------------------------
message("[B] loading 774-unit substrate")
new_root <- file.path(repo, "RNA-seq/results/manuscript_release/candidates",
                      "bulk-corrected-controls-20261008T162308Z/primary_source_controls")
new_meta <- as.data.table(readRDS(file.path(new_root, "meta_matched.rds")))
new_dge <- readRDS(file.path(new_root, "merged_dge.rds"))
stopifnot(nrow(new_meta) == 774L, setequal(colnames(new_dge), new_meta$sample_id))
new_meta <- new_meta[match(colnames(new_dge), sample_id)]
sym_new <- to_symbols(new_dge); rm(new_dge); gc()
three <- c("GSE130970", "GSE135251", "GSE162694")
keep <- c("sample_id", "dataset", "sex_final", "fibrosis_stage", "nas_score")
discB <- new_meta[dataset %in% three & !is.na(fibrosis_stage) & !is.na(nas_score), ..keep]
stopifnot(nrow(discB) == 364L)
bins_new <- t3_expression_deciles(sym_new)
univB <- nonconstant(sym_new, new_meta, three, rownames(sym_new))
sgB <- score_gene_sets(sym_new, new_meta, vocab, gene_set_min_genes, gene_set_coverage, FALSE)
testB <- sgB$coverage[testable == TRUE, feature_id]
message("[B] testable ", length(testB), " universe ", length(univB))
resB <- run_arm("B_774_three_cohorts",
                sgB$scores[testB, discB$sample_id, drop = FALSE],
                t3_designs(discB, three),
                t3_cohort_z(sym_new, new_meta, three, univB),
                univB, bins_new, testB)

# ---- Arm C: 774 units plus unchanged GSE174478 records (Figure 4E arm) ----
message("[C] combining with GSE174478 from the old matrix")
m174 <- old_meta[dataset == "GSE174478"]
sym174 <- sym_old[, m174$sample_id, drop = FALSE]
univC <- intersect(univB, rownames(sym174))
univC <- nonconstant(sym174, m174, "GSE174478", univC)
sg174 <- score_gene_sets(sym174, m174, vocab, gene_set_min_genes, gene_set_coverage, FALSE)
d174 <- m174[!is.na(fibrosis_stage) & !is.na(nas_score), ..keep]
discC <- rbind(discB, d174)
stopifnot(nrow(discC) == 457L)
S_C <- cbind(sgB$scores[testB, discB$sample_id, drop = FALSE],
             sg174$scores[testB, d174$sample_id, drop = FALSE])
fourC <- c(three, "GSE174478")
zC <- c(t3_cohort_z(sym_new, new_meta, three, univC),
        t3_cohort_z(sym174, m174, "GSE174478", univC))
message("[C] testable ", length(testB), " universe ", length(univC),
        " GSE174478-unscored testable programs ",
        sum(!testB %in% sg174$coverage[testable == TRUE, feature_id]))
resC <- run_arm("C_774_plus_GSE174478", S_C, t3_designs(discC, fourC), zC,
                univC, bins_new, testB)

fwrite(rbindlist(list(resA, resB, resC)), file.path(out_dir, "summary_all_arms.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("done")
