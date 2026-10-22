#!/usr/bin/env Rscript
# N4: within-person agreement of the 117 frozen Hotspot program scores against recorded
# fibrosis stage. Implements <out>/PRESPEC.md. Persons are the unit; libraries of one
# person are never treated as independent. Counts are the HISTORICAL fragment counts;
# results must be recomputed after the bulk recount.
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
n1_dir <- if (length(args) >= 2) args[[2]] else "none"
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
source(file.path(repo, "scripts/analysis/hypotheses_2026q4/N4_within_person_agreement/n4_null_eligibility.R"))
SEED <- 20261007L
B_BOOT <- 20000L
B_NULL <- 2000L
MIN_TEST <- 15L
MIN_REPORT <- 5L
writeLines(sprintf('{"bootstrap_seed": %d, "null_seed": %d, "spearman_null_seed": %d, "B_boot": %d, "B_null": %d}',
                   SEED, SEED + 1L, SEED + 2L, B_BOOT, B_NULL), file.path(out, "seeds.json"))

counts_dir <- file.path(repo, "results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T194845Z/frozen_sets/read_counts")
id_dir <- file.path(repo, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified")
roster <- fread(file.path(id_dir, "resource_identity_roster.tsv"))
bulk <- fread(file.path(id_dir, "bulk_roster.tsv"))
membership <- fread(file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"))
registry <- fread(file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv"))
annotation <- fread(file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))
pairs <- fread(file.path(out, "pair_type_table.tsv"))
stopifnot(nrow(registry) == 117L, uniqueN(membership$program_uid) == 117L)

# ---- GSE193066 eligibility from N1 -------------------------------------------------
pairs[, n1_eligible := NA]
ok_runs <- character()
n1_note <- "N1 outcome not supplied; GSE193066 pairs listed, not analysed"
if (n1_dir != "none" && file.exists(file.path(n1_dir, "summary.json"))) {
  n1_summary <- paste(readLines(file.path(n1_dir, "summary.json")), collapse = " ")
  accepted <- grepl('"method_accepted": true', n1_summary, fixed = TRUE)
  calls <- fread(file.path(n1_dir, "gse193066_library_calls.tsv"))
  ok_runs <- calls[n1_status == "not_supported_as_mixed", run]
  for (i in which(pairs$cohort == "GSE193066" & pairs$definition == "source_declared_person")) {
    runs <- strsplit(pairs$libraries[i], ",")[[1]]
    pairs$n1_eligible[i] <- accepted && isTRUE(pairs$genotype_same_individual[i]) && all(runs %in% ok_runs)
  }
  n1_note <- sprintf("N1 method accepted: %s; GSE193066 source pairs eligible: %d of %d", accepted,
                     sum(pairs$n1_eligible, na.rm = TRUE), sum(pairs$cohort == "GSE193066" & pairs$definition == "source_declared_person"))
}
message(n1_note)
pairs[, analysed := (cohort != "GSE193066" & definition == "genotype_component") |
        (cohort == "GSE193066" & definition == "source_declared_person" & n1_eligible %in% TRUE)]
fwrite(pairs, file.path(out, "pair_type_table_with_n1.tsv"), sep = "\t")

# ---- scoring -------------------------------------------------------------------------
base_id <- function(x) sub("\\.[0-9]+.*$", "", x)
annotation[, symbol := toupper(gene_name)]
id2sym <- annotation[symbol != "" & !is.na(symbol), .(symbol = symbol[1]), by = ensembl_base]
memb <- membership[!is.na(mapped_symbol) & mapped_symbol != "",
                   .(weight = sum(as.numeric(original_l1_weight))), by = .(program_uid, symbol = toupper(mapped_symbol))]
total_w <- membership[, .(total = sum(as.numeric(original_l1_weight))), by = program_uid]

cohorts <- unique(pairs[analysed == TRUE, cohort])
scores <- list(); testability <- list(); logcpm_sym <- list(); ref_libs <- list()
for (co in cohorts) {
  ct <- fread(file.path(counts_dir, co, "gene_counts.txt"), skip = 1)
  m <- as.matrix(ct[, 7:ncol(ct)])
  colnames(m) <- sub("^.*?(SRR[0-9]+).*$", "\\1", basename(colnames(m)))
  rownames(m) <- ct[[1]]
  lc <- cpm(DGEList(m), log = TRUE, prior.count = 1)   # library size = all counted genes
  lc <- lc[rowSums(m) > ncol(m), , drop = FALSE]
  sym <- id2sym$symbol[match(base_id(rownames(lc)), id2sym$ensembl_base)]
  ok <- !is.na(sym)
  s <- rowsum(lc[ok, , drop = FALSE], sym[ok], reorder = FALSE)
  s <- s / as.numeric(table(sym[ok])[rownames(s)])
  ref <- if (co == "GSE193066") roster[dataset == co & first_biopsy %in% c(TRUE, "True", "TRUE"), sample_id]
         else bulk[dataset == co & primary_selected %in% c(TRUE, "True", "TRUE"), sample_id]
  ref <- intersect(ref, colnames(s))
  stopifnot(length(ref) >= 30L)
  mu <- rowMeans(s[, ref]); sdv <- apply(s[, ref], 1, sd)
  keep <- is.finite(sdv) & sdv > 0
  z <- (s[keep, ] - mu[keep]) / sdv[keep]
  sc <- matrix(NA_real_, nrow(registry), ncol(z), dimnames = list(registry$program_uid, colnames(z)))
  for (p in registry$program_uid) {
    mm <- memb[program_uid == p & symbol %in% rownames(z)]
    frac <- sum(mm$weight) / total_w[program_uid == p, total]
    test <- nrow(mm) >= 10L && frac >= 0.8
    testability[[length(testability) + 1L]] <- data.table(cohort = co, program_uid = p, n_observed_genes = nrow(mm),
                                                         observed_l1_fraction = frac, testable = test)
    if (!test) next
    v <- as.numeric(crossprod(mm$weight / sum(mm$weight), z[mm$symbol, , drop = FALSE]))
    sc[p, ] <- (v - mean(v[match(ref, colnames(z))])) / sd(v[match(ref, colnames(z))])
  }
  scores[[co]] <- sc; logcpm_sym[[co]] <- s; ref_libs[[co]] <- ref
  message(co, ": libraries ", ncol(sc), ", reference persons ", length(ref), ", programs scored ", sum(!is.na(sc[, 1])))
}
testability <- rbindlist(testability)
fwrite(testability, file.path(out, "program_testability.tsv"), sep = "\t")

# ---- long table of analysed libraries -----------------------------------------------
lib <- pairs[analysed == TRUE, .(sample_id = strsplit(libraries, ",")[[1]]), by = .(cohort, person_id, pair_type)]
lib <- merge(lib, roster[, .(sample_id, fibrosis_stage, nas_score)], by = "sample_id", all.x = TRUE, sort = FALSE)
lib[, in_counts := mapply(function(s, co) s %in% colnames(scores[[co]]), sample_id, cohort)]
fwrite(lib, file.path(out, "analysed_libraries.tsv"), sep = "\t")
lib <- lib[in_counts == TRUE]
lib <- lib[person_id %in% lib[, .N, by = person_id][N >= 2, person_id]]

feature_matrix <- function(d) {
  y <- t(do.call(cbind, lapply(seq_len(nrow(d)), function(i) scores[[d$cohort[i]]][, d$sample_id[i]])))
  colnames(y) <- registry$program_uid
  cbind(y, stage = d$fibrosis_stage, nas = d$nas_score)
}

# one-way random ICC(1,1) with unbalanced n0, vectorised over features (columns of y)
icc_oneway <- function(y, g) {
  g <- as.integer(factor(g))
  n_i <- tabulate(g); k <- length(n_i); N <- length(g)
  s <- rowsum(y, g); q <- rowsum(y^2, g)
  grand <- colSums(s) / N
  ssw <- colSums(q - s^2 / n_i)
  ssb <- colSums(n_i * (sweep(s / n_i, 2, grand))^2)
  msb <- ssb / (k - 1); msw <- ssw / (N - k)
  n0 <- (N - sum(n_i^2) / N) / (k - 1)
  (msb - msw) / (msb + (n0 - 1) * msw)
}

# bootstrap over persons using per-person sufficient statistics
icc_boot <- function(y, g, B, seed) {
  g <- as.integer(factor(g))
  n_i <- tabulate(g); k <- length(n_i)
  s <- rowsum(y, g); q <- rowsum(y^2, g)
  set.seed(seed)
  cnt <- t(replicate(B, tabulate(sample.int(k, k, replace = TRUE), k)))   # B x persons
  N <- as.numeric(cnt %*% n_i)
  grand <- (cnt %*% s) / N
  ssw <- cnt %*% (q - s^2 / n_i)
  ssb <- cnt %*% (s^2 / n_i) - N * grand^2
  n0 <- (N - as.numeric(cnt %*% n_i^2) / N) / (k - 1)
  msb <- ssb / (k - 1); msw <- ssw / (N - k)
  (msb - msw) / (msb + (n0 - 1) * msw)
}

bh <- function(p) p.adjust(p, method = "BH")

weighted_kappa_sym <- function(a, b) {   # linear weights, double entry so library order does not matter
  x <- c(a, b); y <- c(b, a)
  lev <- sort(unique(x)); if (length(lev) < 2) return(NA_real_)
  w <- 1 - abs(outer(lev, lev, "-")) / (max(lev) - min(lev))
  tab <- table(factor(x, lev), factor(y, lev)) / length(x)
  po <- sum(w * tab); pe <- sum(w * outer(rowSums(tab), colSums(tab)))
  (po - pe) / (1 - pe)
}

# ---- random-pair null: different persons, same cohort, same stage labels -------------
null_pool <- list()
multi_libs <- roster[, if (.N > 1) .SD, by = screen_component]$sample_id
for (co in cohorts) {
  pool <- roster[dataset == co & sample_id %in% ref_libs[[co]] & !(sample_id %in% multi_libs) & !is.na(fibrosis_stage),
                 .(sample_id, fibrosis_stage)]
  if (co == "GSE193066") pool <- roster[dataset == co & sample_id %in% ref_libs[[co]] & !is.na(fibrosis_stage) &
                                          !(source_person_id %in% sub("^GSE193066::", "", pairs[cohort == co, person_id])),
                                        .(sample_id, fibrosis_stage, source_person_id, a1_component)]
  if (co == "GSE193066") {
    # Equal N1 non-mixture eligibility is necessary but does not qualify identity.
    observed_runs <- lib[cohort == co, sample_id]
    observed_components <- roster[sample_id %in% observed_runs, a1_component]
    pool <- pool[sample_id %in% ok_runs & !is.na(a1_component) & a1_component != "NA" &
                   !(a1_component %in% observed_components)]
    # Retain only unambiguous one-to-one source/component units; never select one
    # arbitrary library to resolve conflicting person labels.
    pool <- pool[!(source_person_id %in% pool[, .N, by = source_person_id][N > 1, source_person_id]) &
                   !(a1_component %in% pool[, .N, by = a1_component][N > 1, a1_component])]
  }
  null_pool[[co]] <- pool
}
fwrite(rbindlist(lapply(names(null_pool), function(co) copy(null_pool[[co]])[, cohort := co]), fill = TRUE),
       file.path(out, "eligible_null_libraries.tsv"), sep = "\t")
draw_null <- function(d, substitutions) n4_draw_null(d, null_pool)

# ---- analysis per group ----------------------------------------------------------------
groups <- unique(rbind(lib[, .(group = paste0(pair_type, "|pooled"), pair_type, cohort_filter = "all")],
                       lib[, .(group = paste0(pair_type, "|", cohort), pair_type, cohort_filter = cohort)]))
res_all <- list(); stage_all <- list(); summ <- list()
for (gi in seq_len(nrow(groups))) {
  gname <- groups$group[gi]
  d <- lib[pair_type == groups$pair_type[gi] & (groups$cohort_filter[gi] == "all" | cohort == groups$cohort_filter[gi])]
  n_person <- uniqueN(d$person_id)
  status <- if (n_person >= MIN_TEST) "tested" else if (n_person >= MIN_REPORT) "UNDERPOWERED_estimates_only" else "listed_only"
  first2 <- d[, .SD[1:2], by = person_id]
  a <- first2[, fibrosis_stage[1], by = person_id]$V1; b <- first2[, fibrosis_stage[2], by = person_id]$V1
  st <- data.table(group = gname, n_persons = n_person, n_libraries = nrow(d), status = status,
                   cohorts = paste(sort(unique(d$cohort)), collapse = ","),
                   stage_exact_agreement = sum(a == b), stage_discordant_ge1 = sum(a != b),
                   stage_weighted_kappa = weighted_kappa_sym(a, b))
  if (status == "listed_only") { stage_all[[gname]] <- st; next }
  y <- feature_matrix(d)
  has_nas <- all(is.finite(y[, "nas"]))
  if (!has_nas) y <- y[, colnames(y) != "nas", drop = FALSE]
  testable <- colSums(!is.finite(y)) == 0
  yt <- y[, testable, drop = FALSE]
  obs <- icc_oneway(yt, d$person_id)
  boot <- icc_boot(yt, d$person_id, B_BOOT, SEED)
  ci <- apply(boot, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
  dstar <- boot - boot[, "stage"]
  dci <- apply(dstar, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
  p_boot <- (1 + colSums(is.na(dstar) | dstar <= 0)) / (B_BOOT + 1)
  st[, `:=`(stage_icc = obs[["stage"]], stage_icc_lo = ci[1, "stage"], stage_icc_hi = ci[2, "stage"],
            nas_icc = if (has_nas) obs[["nas"]] else NA_real_,
            nas_icc_lo = if (has_nas) ci[1, "nas"] else NA_real_, nas_icc_hi = if (has_nas) ci[2, "nas"] else NA_real_)]
  # null
  capacity <- n4_null_capacity(d, null_pool)
  capacity[, group := gname]
  fwrite(capacity, file.path(out, paste0("null_capacity_", gi, ".tsv")), sep = "\t")
  null_testable <- all(capacity$sufficient)
  st[, null_status := if (null_testable) "exact_stage_eligible_retrospective_diagnostic" else "untestable_insufficient_exact_stage_eligible_units"]
  st[, identity_state := if ("GSE193066" %in% d$cohort) "GSE193066_unqualified_retrospective_only" else "genetic_clusters_not_verified_people"]
  set.seed(SEED + 1L)
  subs <- 0L
  null_icc <- matrix(NA_real_, B_NULL, ncol(yt), dimnames = list(NULL, colnames(yt)))
  for (bi in if (null_testable) seq_len(B_NULL) else integer()) {
    dn <- draw_null(d, subs); subs <- dn$substitutions
    dd <- copy(d)[, sample_id := dn$sample_id]
    yn <- feature_matrix(dd)[, colnames(yt), drop = FALSE]
    null_icc[bi, ] <- icc_oneway(yn, d$person_id)
  }
  p_null <- if (null_testable) (1 + colSums(sweep(null_icc, 2, obs, ">="), na.rm = TRUE)) / (B_NULL + 1) else rep(NA_real_, ncol(yt))
  r <- data.table(group = gname, status = status, n_persons = n_person, feature = colnames(y), testable = testable)
  r[testable == TRUE, `:=`(icc = obs, icc_lo = ci[1, ], icc_hi = ci[2, ],
                           diff_vs_stage = obs - obs[["stage"]], diff_lo = dci[1, ], diff_hi = dci[2, ],
                           p_boot_diff_le0 = p_boot,
                           null_icc_median = apply(null_icc, 2, median, na.rm = TRUE),
                           null_icc_p95 = apply(null_icc, 2, quantile, 0.95, na.rm = TRUE), p_null = p_null)]
  r[, is_program := feature %in% registry$program_uid]
  if (status == "tested") {
    r[is_program == TRUE, `:=`(p_boot_family = fifelse(testable, p_boot_diff_le0, 1), p_null_family = fifelse(testable, p_null, 1))]
    # Bootstrap tail fractions are descriptive diagnostics, not calibrated tests.
    r[is_program == TRUE, `:=`(bh_q_exceeds_stage = NA_real_, bh_q_exceeds_null = if (null_testable) bh(p_null_family) else NA_real_)]
  }
  r <- merge(r, registry[, .(feature = program_uid, cell_type, module, module_name)], by = "feature", all.x = TRUE, sort = FALSE)
  res_all[[gname]] <- r; stage_all[[gname]] <- st
  pr <- r[is_program == TRUE & testable == TRUE][order(-icc)]
  med <- pr[which.min(abs(icc - median(icc)))]
  summ[[gname]] <- data.table(
    group = gname, status = status, null_status = st$null_status, identity_state = st$identity_state,
    n_persons = n_person, programs_testable = nrow(pr),
    stage_icc = obs[["stage"]], stage_icc_lo = ci[1, "stage"], stage_icc_hi = ci[2, "stage"],
    best_program = pr$feature[1], best_name = paste(pr$cell_type[1], pr$module_name[1]),
    best_icc = pr$icc[1], best_lo = pr$icc_lo[1], best_hi = pr$icc_hi[1],
    median_program = med$feature, median_icc = med$icc, median_lo = med$icc_lo, median_hi = med$icc_hi,
    n_point_estimate_above_stage = sum(pr$icc > obs[["stage"]]),
    n_exceed_stage_bh05 = NA_integer_,
    n_exceed_null_bh05 = if (status == "tested" && null_testable) sum(r$bh_q_exceeds_null < 0.05, na.rm = TRUE) else NA_integer_,
    null_median_of_program_medians = median(pr$null_icc_median), null_p95_median = median(pr$null_icc_p95),
    null_stage_substitutions = subs, min_attainable_p_boot = 1 / (B_BOOT + 1), min_attainable_p_null = 1 / (B_NULL + 1))
}
fwrite(rbindlist(stage_all, fill = TRUE), file.path(out, "stage_agreement_by_group.tsv"), sep = "\t")
if (length(res_all)) {
  fwrite(rbindlist(res_all, fill = TRUE), file.path(out, "icc_by_group_and_feature.tsv"), sep = "\t")
  fwrite(rbindlist(summ, fill = TRUE), file.path(out, "summary_by_group.tsv"), sep = "\t")
}

# ---- per-person scores and within-person differences ----------------------------------
yl <- feature_matrix(lib)
fwrite(cbind(lib[, .(sample_id, cohort, person_id, pair_type)], as.data.table(yl)), file.path(out, "program_scores_analysed_libraries.tsv.gz"), sep = "\t")

# ---- descriptive: genome-wide Spearman within real and random pairs --------------------
set.seed(SEED + 2L)
sp <- list()
for (co in unique(lib$cohort)) {
  rk <- apply(logcpm_sym[[co]], 2, rank)
  d <- lib[cohort == co]
  for (p in unique(d$person_id)) {
    r <- d[person_id == p, sample_id]
    sp[[length(sp) + 1L]] <- data.table(cohort = co, kind = "real_pair", pair_type = d[person_id == p, pair_type[1]], id = p,
                                        a = r[1], b = r[2], spearman = cor(rk[, r[1]], rk[, r[2]]),
                                        stage_a = d[sample_id == r[1], fibrosis_stage], stage_b = d[sample_id == r[2], fibrosis_stage])
    for (j in 1:20) {
      if (!all(n4_null_capacity(d[person_id == p][1:2], null_pool)$sufficient)) next
      dn <- draw_null(d[person_id == p][1:2], 0L)$sample_id
      sp[[length(sp) + 1L]] <- data.table(cohort = co, kind = "random_pair_same_cohort_stage", pair_type = d[person_id == p, pair_type[1]],
                                          id = paste0(p, "::null", j), a = dn[1], b = dn[2], spearman = cor(rk[, dn[1]], rk[, dn[2]]),
                                          stage_a = d[person_id == p, fibrosis_stage[1]], stage_b = d[person_id == p, fibrosis_stage[2]])
    }
  }
}
sp <- rbindlist(sp)
fwrite(sp, file.path(out, "pair_genomewide_spearman.tsv"), sep = "\t")
fwrite(sp[, .(n = .N, median = median(spearman), q05 = quantile(spearman, 0.05), q95 = quantile(spearman, 0.95),
              min = min(spearman), max = max(spearman)), by = .(cohort, pair_type, kind)],
       file.path(out, "pair_genomewide_spearman_summary.tsv"), sep = "\t")

writeLines(c(n1_note, "", capture.output(sessionInfo())), file.path(out, "sessionInfo.txt"))
message("N4 done")
