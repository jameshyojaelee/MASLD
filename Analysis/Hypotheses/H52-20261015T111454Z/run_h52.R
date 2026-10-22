#!/usr/bin/env Rscript
# H52: does program-pair covariance found in three discovery cohorts recur, with
# the same signs and the same community structure, in GSE213621?
#
# Implements PRESPEC_H52.md (copied into this run directory; sha256 verified
# below, the script refuses to run on any changed prespecification or input).
# Pipeline = job 3811 feasibility.R functions with the T4 library sourced
# unchanged from the repository (read only). mod_of(), rewire_once() and dens()
# are taken unchanged from the 2026-09-08 script A_coordination_replication.R by
# parsing that hash-verified file and evaluating only those three assignments;
# their edge-weight column `z_holdout` carries the GSE213621 z_partial_k1.
#
# Modes. Final run: no H52_* environment overrides; writes to <run>/out, which
# must be empty. Debug run (H52_DEBUG=1): every GSE213621 expression value is
# permuted per gene across units BEFORE anything is computed on GSE213621, so no
# real GSE213621 edge statistic is ever produced; outputs go to
# <run>/tmp/debug_<H52_TAG>/ and draw counts may be reduced.

suppressPackageStartupMessages({library(data.table); library(edgeR); library(Matrix)})
options(warn = 1, width = 250, datatable.print.nrows = 200)

RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H52-20261015T111454Z"
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SEED <- 20260923L
EXPECTED_PRESPEC <- "5bab40f0efdc2794c0cb19b7d2e8ad511c0871082198d62a06fab6c2b6d3c8c6"

DEBUG <- identical(Sys.getenv("H52_DEBUG"), "1")
env_int <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) as.integer(v) else d }
overrides <- c("H52_NBOOT", "H52_NDRAW", "H52_NNC", "H52_NLIN", "H52_STOP_AFTER", "H52_TAG")
if (!DEBUG && any(nzchar(Sys.getenv(overrides)))) stop("H52_* overrides are allowed only with H52_DEBUG=1", call. = FALSE)
N_DRAW <- env_int("H52_NDRAW", 1000L)   # every F1, PC, F2 and descriptive null
N_BOOT <- env_int("H52_NBOOT", 1000L)   # unit bootstrap resamples
N_NC <- env_int("H52_NNC", 100L)        # permutations per negative control
N_LIN <- env_int("H52_NLIN", 200L)      # lineage-matched draws per resample / permutation
STOP_AFTER <- Sys.getenv("H52_STOP_AFTER", "")

OUT <- if (DEBUG) file.path(RUN, "tmp", paste0("debug_", Sys.getenv("H52_TAG", "x"))) else file.path(RUN, "out")
dir.create(OUT, showWarnings = FALSE, recursive = DEBUG)
if (!DEBUG && length(list.files(OUT, all.files = TRUE, no.. = TRUE)) > 0L)
  stop("out/ is not empty; nothing is overwritten", call. = FALSE)
cat("H52 run; DEBUG =", DEBUG, "; out =", OUT, "; draws", N_DRAW, "boot", N_BOOT, "nc", N_NC, "lin", N_LIN, "\n")
T0 <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1fs] ", as.numeric(difftime(Sys.time(), T0, units = "secs"))), ..., "\n", sep = "")

# ---------------------------------------------------------------- diagnostics
DIAG <- list()
diag_add <- function(id, item, value, expected = "") {
  DIAG[[length(DIAG) + 1L]] <<- data.table(diagnostic = id, item = item,
                                          value = paste(format(value, digits = 6), collapse = " "),
                                          expected = as.character(expected))
}
diag_write <- function() fwrite(rbindlist(DIAG), file.path(OUT, "diagnostics.tsv"), sep = "\t")
f4 <- function(x) sprintf("%.4f", x)

stop_rule <- function(rule, reason, numbers = character()) {
  if (length(DIAG)) diag_write()
  txt <- c("# RESULT H52", "",
           if (DEBUG) "**DEBUG RUN on permuted GSE213621 data. Not a result.**" else NULL, "",
           sprintf("Verdict: **untestable** (stop rule %s).", rule), "",
           sprintf("Reason: %s", reason), "",
           "Numbers reached:", "", paste0("- ", numbers), "",
           "No GSE213621 statistic was computed after this point. The title clause stays where it is now:",
           "historical diagnostics only, per title.md:26-38.",
           "", sprintf("Written %s by run_h52.R.", format(Sys.time(), tz = "UTC", usetz = TRUE)))
  writeLines(txt, file.path(OUT, "RESULT_H52.md"))
  cat("STOP RULE", rule, ":", reason, "\n")
  quit(save = "no", status = 0)
}

# ------------------------------------------------- 0. hashes and sessionInfo
sha <- function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE)[[1]], "[[:space:]]+")[[1]][[1]]
rel <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
cand <- "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10"
T4d <- file.path(cand, "workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction")
hot <- "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
rr <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
hold <- "RNA-seq/results/gse193066_holdout/holdout-20260907T105813Z"
INPUTS <- data.table(path = c(
  file.path(rel, "merged_dge.rds"), file.path(rel, "meta_matched.rds"),
  file.path(T4d, "program_communities.tsv"), file.path(T4d, "coordinated_edges.tsv"),
  file.path(T4d, "program_interaction_edges.tsv.gz"), file.path(T4d, "interaction_summary.tsv"),
  file.path(hot, "program_membership_v2.tsv"), file.path(hot, "program_registry_v2.tsv"),
  file.path(cand, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  "scripts/manuscript/program_observability_map/t4_interaction_lib.R",
  "scripts/manuscript/program_observability_map/analysis_lib.R",
  file.path(rr, "jobs/3811/out/approx_rediscovery_coarse_371.tsv.gz"),
  file.path(rr, "jobs/3811/out/result_B.tsv"), file.path(rr, "jobs/3811/out/result_B2.tsv"),
  file.path(rr, "jobs/3811/out/result_C.tsv"), file.path(rr, "jobs/3811/out/result_D.tsv"),
  file.path(rr, "jobs/3811/out/result_D2.tsv"), file.path(rr, "jobs/3811/feasibility.R"),
  file.path(rr, "scratch/validity-L10-01/out/loco_composition_proxy.tsv"),
  file.path(rr, "scratch/validity-L10-01/out/loco_s1_by_subset.tsv"),
  file.path(rr, "scratch/validity-L10-01/out/boot_GSE135251.tsv"),
  file.path(rr, "scratch/plan-H52/benchmark_s1.R"), file.path(rr, "scratch/plan-H52/benchmark_s1.tsv"),
  file.path(hold, "prespec/PRESPEC_A_coordination.md"), file.path(hold, "A_coordination_replication.R"),
  file.path(hold, "results/A_coordination_results.tsv"),
  file.path(rr, "jobs/3811/out/result_A.tsv"), file.path(rr, "jobs/3811/out/result_B3.tsv"),
  file.path(rr, "jobs/3811/out/result_E.tsv"), file.path(rr, "jobs/3131/gene_permutation_null.tsv"),
  file.path(rr, "scratch/validity-L10-01/out_boot.log"), file.path(rr, "scratch/validity-L10-01/out_comp.log")),
  expected = c(
  "73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b", "429e27376ad421f5f5fae19c75f5e760030e727fab65c5768ef562a9939edfc3",
  "f4b1b6e51b6e0c3edbdcdcff80654b34ea192957e918c8de4cd5571171bc0480", "af6a58989639aa382dbdb7a39d3c155fbb933773a5e1200b87eb3c17684716a3",
  "645fce85f9d9af20656de53594857999485c3de9729627e6e1aa1e294b769776", "dd0f9f37e22d6d89f570eed7f047c64c2a6a0028bca1507fa71405ffc696d778",
  "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b", "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7",
  "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189",
  "2d64504184a47a13d0c3c9cd2ed05b763e5d5b903a744fb3461dcc6b1d00de24",
  "a6439390506398be9283eb335b3866c45e9a510fa87aec07f1bf7fc21accb71b",
  "1e54abddf5a3587e7565e810bbd58486d048379f6e7a68325c8e44fa3f6291c8",
  "41261c09c72406ed388077024a4b9b19d81f9a28b938d01a387560e8dfbe4a58", "1da96dfea76b84b84576bd4c411507bcd3786f5b1fd4e38e0d4f372e676875c5",
  "7d464e9036ce7f890f70e8651cb6879760134f7458138f9e3cd30c35f56be609", "88cc4ee5b2f223cd6ab86434992d3104aae3620ffdd65b125c69adb6ac4ae733",
  "6f4e606c56ecbf296ea0cf20b6d16055e29ccd58fc89cec7c56479f17750e7ec", "35aa1f59dda0390634452b679dbfa2c32f0f06ad1e5496b167180a852ba478f3",
  "f15239a4dc462ef66e0a5e85dae7ed50fb67850f7df7cd168b4e9562beb23d1c",
  "4f1043621c853afa0eabb3cf8a5b245eae2534ba08efd7ba7aadcb4c0ef80503",
  "274d98cb5f8cd17095dd8701c730e23a1ec89a179517b8b571b7eb75cd4c0df0",
  "68bf81a8879288b73e790bf436fbf6ba08926211772149dbf4254c610c153b5e", "a690cc1b968f85ab3e73c653789dc0dbb1355e941cab7d39b06e0560f1ceec3b",
  "fd00249637a1aa900f9c3504684e77b8a3f337e734e47528caba7e1f7b6a53e4", "f8d2f9bc90e94eacd813f9fe544c94acb7fbf101f7def81ca2849d0a7935fbea",
  "89ec85216af0c55c74077caf60b1394bb69a900d76d4d25f63d38cb4bc60bb14",
  "640b42b82436c6cb96114d6f810cf7f26f319ef8ae68e43d51ef3d4d5474e173", "e92238d344dcf8441c60c5317e710664a1506df1488f53220adc80ca3c2da0ba",
  "ced64d0b952fbf1b902eda98b854260a1d1e9b3065f08554c3b7b842e97f7294", "6f0ee1d345ce80950503b26c1c1a9b265bdb237d1822aa36c53893031117cce0",
  "9bd2f8a312bcaa42f45492eb4e0ac0b0ea3fb7132338916cf9ee346893aff020", "aebee69752bf4c7d993c519bf69abb5fe8bd555a91ce38ed59936e6311ec4e7d"))
P <- function(p) file.path(ROOT, p)

prespec_got <- sha(file.path(RUN, "PRESPEC_H52.md"))
prespec_rec <- strsplit(readLines(file.path(RUN, "PRESPEC_H52.sha256"))[[1]], "[[:space:]]+")[[1]][[1]]
if (!identical(prespec_got, EXPECTED_PRESPEC) || !identical(prespec_rec, EXPECTED_PRESPEC))
  stop_rule(1, "prespecification sha256 differs", c(sprintf("got %s, recorded %s, expected %s", prespec_got, prespec_rec, EXPECTED_PRESPEC)))
INPUTS[, got := vapply(path, function(p) if (file.exists(P(p))) sha(P(p)) else "MISSING", character(1))]
INPUTS[, ok := got == expected]
fwrite(INPUTS, file.path(OUT, "input_hashes.tsv"), sep = "\t")
if (!all(INPUTS$ok)) stop_rule(1, "input sha256 differs from the prespecification table",
                               INPUTS[ok == FALSE, sprintf("%s: got %s expected %s", path, got, expected)])
tick("hashes verified: prespec + ", nrow(INPUTS), " inputs")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

# --------------------------------------------------------------- libraries
source(P("scripts/manuscript/program_observability_map/analysis_lib.R"))
source(P("scripts/manuscript/program_observability_map/t4_interaction_lib.R"))
A_SCRIPT <- P(file.path(hold, "A_coordination_replication.R"))
for (e in parse(A_SCRIPT, keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
      as.character(e[[2L]]) %in% c("mod_of", "rewire_once", "dens")) eval(e, envir = globalenv())
}
stopifnot(exists("mod_of"), exists("rewire_once"), exists("dens"))

# ------------------------------------------------------------- frozen inputs
annotation <- fread(P(file.path(cand, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")),
                    select = c("gene_id", "gene_name"))
communities <- fread(P(file.path(T4d, "program_communities.tsv")))
features <- communities$feature_id
membership <- fread(P(file.path(hot, "program_membership_v2.tsv")))
registry <- fread(P(file.path(hot, "program_registry_v2.tsv")), select = c("program_uid", "cell_type", "module_name"))
weights <- t4_membership_weights(membership, features)
pairs0 <- t4_pair_table(features, weights)
key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
pairs0[, ekey := key(program_a, program_b)]
NP <- nrow(pairs0)
frozen_coord <- fread(P(file.path(T4d, "coordinated_edges.tsv")))[, ekey := key(program_a, program_b)]
frozen_all <- fread(P(file.path(T4d, "program_interaction_edges.tsv.gz")),
                    select = c("program_a", "program_b", "r_partial", "membership_entangled", "module_name_a", "module_name_b"))
frozen_all[, ekey := key(program_a, program_b)]
frozen_3811 <- fread(P(file.path(rr, "jobs/3811/out/approx_rediscovery_coarse_371.tsv.gz")))

part <- setNames(communities$community, communities$feature_id)
lin <- registry[match(communities$feature_id, registry$program_uid), cell_type]
names(lin) <- communities$feature_id
stopifnot(!anyNA(lin))
c3 <- communities[community == 3, feature_id]
c3_lin <- table(lin[c3])
modname <- unique(rbind(frozen_all[, .(f = program_a, m = module_name_a)], frozen_all[, .(f = program_b, m = module_name_b)]))
modname <- setNames(modname$m, modname$f)[features]
TECH_RX <- "Translation|ribosome|Oxidative phosphorylation|Mitochondrial|Immediate-early|AP-1|housekeeping"
tech <- features[grepl(TECH_RX, modname)]
stopifnot(length(tech) == 11L, all(part[tech] == 1L))
pairs0[, `:=`(la = lin[program_a], lb = lin[program_b], ca = part[program_a], cb = part[program_b])]
pairs0[, `:=`(cross = la != lb, tech_any = program_a %in% tech | program_b %in% tech)]
cat("frozen: programs", length(features), "communities", paste(table(part), collapse = "/"), "c3 lineage", paste(names(c3_lin), c3_lin), "\n")

# ------------------------------------------------------------- 774 roster
meta <- as.data.table(readRDS(P(file.path(rel, "meta_matched.rds"))))
stopifnot(nrow(meta) == 774L, uniqueN(meta$analysis_unit_id) == 774L, !anyNA(meta$group_binary), !anyNA(meta$dataset))
DISC <- c("GSE130970", "GSE135251", "GSE162694")
md <- meta[dataset %in% DISC & group_binary == "Disease"]
g <- meta[dataset == "GSE213621" & group_binary == "Disease"]
bin3 <- function(f) ifelse(f <= 1, 1L, ifelse(f == 2, 2L, 3L))
stopifnot(!anyNA(md$fibrosis_stage), !anyNA(g$fibrosis_stage), !anyNA(md$sex_final), !anyNA(g$sex_final))
md[, fibrosis_stage := bin3(fibrosis_stage)]
cnt_d <- md[, .N, by = dataset][order(dataset)]
cnt_g <- g[, .N, by = fibrosis_stage][order(fibrosis_stage)]
stopifnot(identical(cnt_d$N, c(65L, 197L, 109L)), nrow(md) == 371L, nrow(g) == 258L,
          identical(cnt_g$N, c(87L, 94L, 77L)), identical(cnt_g$fibrosis_stage, 1:3),
          all(g[, table(condition, fibrosis_stage)][cbind(c("Fibrosis_F0F1", "Fibrosis_F2", "Fibrosis_F3F4"), c("1", "2", "3"))] == c(87, 94, 77)),
          all(g$sex_source == "inferred_kmeans"))
stopifnot(md[, .N, by = .(dataset, fibrosis_stage)][order(dataset, fibrosis_stage)]$N ==
            c(44, 9, 12, 81, 51, 65, 65, 27, 17))
stopifnot(all(c("SRR21622926") %in% g$sample_id), "SRR9036348" %in% md$sample_id)
diag_add("roster", "discovery disease units by cohort (GSE130970/GSE135251/GSE162694)", cnt_d$N, "65 197 109")
diag_add("roster", "GSE213621 disease units by bin (F0-F1/F2/F3-F4)", cnt_g$N, "87 94 77")

# ---------------------------------------------------- 1. log-CPM, library metrics
dge <- readRDS(P(file.path(rel, "merged_dge.rds")))
stopifnot(all(c(md$sample_id, g$sample_id) %in% colnames(dge$counts)))
ids_all <- c(md$sample_id, g$sample_id)
cnt <- dge$counts[, ids_all]
libsize <- colSums(cnt)
prp <- sweep(cnt, 2, libsize, "/")
ent <- apply(prp, 2, function(p) { p <- p[p > 0]; -sum(p * log(p)) })
top100 <- apply(prp, 2, function(p) sum(sort(p, decreasing = TRUE)[1:100]))
rm(prp, cnt); invisible(gc())
libdt <- data.table(sample_id = ids_all, log_lib = log10(libsize[ids_all]), entropy = ent[ids_all], top100 = top100[ids_all])
md <- merge(md, libdt, by = "sample_id", sort = FALSE); g <- merge(g, libdt, by = "sample_id", sort = FALSE)
lc <- edgeR::cpm(dge, log = TRUE, prior.count = 1)     # full 774, as job 3811
rm(dge); invisible(gc())
MK <- list(chol = c("KRT7", "EPCAM", "SOX9", "KRT23", "MMP7", "HNF1B"), fib = c("COL1A1", "COL1A2", "COL3A1", "LUM"),
           endo = c("VWF", "CDH5"), bcell = c("MS4A1", "CD79A"))
stopifnot(!any(unlist(MK) %in% weights$gene_symbol))
sym_lc <- annotation$gene_name[match(sub("\\..*$", "", rownames(lc)), sub("\\..*$", "", annotation$gene_id))]
mk_rows <- lapply(MK, function(s) which(sym_lc %in% s))
lcm <- lc[unlist(mk_rows), , drop = FALSE]
mk_rows <- split(seq_len(nrow(lcm)), rep(names(MK), lengths(mk_rows)))[names(MK)]
sm <- collapse_symbols(lc, rownames(lc), annotation); rm(lc); invisible(gc())
panel <- intersect(sort(unique(weights$gene_symbol)), rownames(sm))
sm <- sm[panel, , drop = FALSE]
cat("panel symbols", length(panel), "of", uniqueN(weights$gene_symbol), "\n")
diag_add("pipeline", "panel symbols present in 774 universe", length(panel), "3296")

if (DEBUG) {   # permute every GSE213621 value per gene across its 258 units; nothing real from GSE213621 is used
  set.seed(777L)
  sm[, g$sample_id] <- t4_permute_genes(sm[, g$sample_id, drop = FALSE])
  lcm[, g$sample_id] <- t4_permute_genes(lcm[, g$sample_id, drop = FALSE])
  cat("DEBUG: GSE213621 expression permuted per gene\n")
}

# ----------------------------------------------------------- pipeline helpers
COARSE <- ~ sex_final + fibrosis_factor
COMP <- ~ sex_final + fibrosis_factor + prx_chol + prx_fib + prx_endo + prx_bcell
LIB <- ~ sex_final + fibrosis_factor + log_lib + entropy + top100
ACC <- ~ sex_final + fibrosis_factor + acc_q
make_block <- function(mat, mc, formula) { T4_RESIDUAL_FORMULA <<- formula; t4_build_block(mat, mc, weights, features, 0.80) }
add_proxies <- function(mc, cols = mc$sample_id) {   # within-cohort z over the supplied units
  for (n in names(MK)) {
    z <- t(scale(t(lcm[mk_rows[[n]], cols, drop = FALSE])))
    set(mc, j = paste0("prx_", n), value = colMeans(z, na.rm = TRUE))
  }
  mc
}
for (cc in DISC) md[dataset == cc, c("prx_chol", "prx_fib", "prx_endo", "prx_bcell") :=
                      as.list(add_proxies(copy(md[dataset == cc]))[, .(prx_chol, prx_fib, prx_endo, prx_bcell)])]

approx_discovery <- function(blocks) {    # job 3811, plus the k = 3 columns for the descriptive sensitivity
  pairs <- copy(pairs0)
  ops <- lapply(blocks, t4_disjoint_operators, pairs = pairs, features = features, min_retained_weight = 0.80)
  elig <- Reduce(`&`, lapply(ops, function(op) { e <- rep(TRUE, nrow(pairs)); e[op$pair_index] <- op$eligible; e }))
  obs <- t4_one_pass(blocks, pairs, ops, want_matrix = TRUE)
  w_c <- vapply(blocks, function(b) b$df - 3L, numeric(1))
  wsum <- function(zc) rowSums(matrix(w_c, nrow(zc), ncol(zc), byrow = TRUE) * is.finite(zc))
  pairs[, `:=`(entangled = !elig,
               z_disj = obs$z_disjoint, z_part = obs$z_partial_k1,
               zs_disj = obs$z_disjoint * sqrt(wsum(obs$zc_disjoint)),
               zs_part = obs$z_partial_k1 * sqrt(wsum(obs$zc_partial[[1L]])),
               z_part3 = obs$z_partial_k3, zs_part3 = obs$z_partial_k3 * sqrt(wsum(obs$zc_partial[[2L]])))]
  fam <- pairs[, !entangled & is.finite(zs_disj) & is.finite(zs_part)]
  pairs[, `:=`(q_disj = NA_real_, q_part = NA_real_)]
  pairs[fam, q_disj := p.adjust(2 * pnorm(-abs(zs_disj)), "BH")]
  pairs[fam, q_part := p.adjust(2 * pnorm(-abs(zs_part)), "BH")]
  pairs[, coordinated := fam & q_disj < 0.05 & q_part < 0.05 & z_part > 0]
  pairs[, in_family := fam]
  fam3 <- pairs[, !entangled & is.finite(zs_disj) & is.finite(zs_part3)]
  pairs[, `:=`(q_disj3 = NA_real_, q_part3 = NA_real_)]
  pairs[fam3, q_disj3 := p.adjust(2 * pnorm(-abs(zs_disj)), "BH")]
  pairs[fam3, q_part3 := p.adjust(2 * pnorm(-abs(zs_part3)), "BH")]
  pairs[, `:=`(in_family3 = fam3, coordinated3 = fam3 & q_disj3 < 0.05 & q_part3 < 0.05 & z_part3 > 0)]
  list(pairs = pairs, globals = obs$globals)
}
rep_pass <- function(block, Z = NULL, ops = NULL) {   # job 3811 replicate_z(), vector aligned to pairs0
  if (is.null(ops)) ops <- t4_disjoint_operators(block, pairs0, features, min_retained_weight = 0.80)
  obs <- t4_one_pass(list(block), pairs0, list(ops), Zlist = if (is.null(Z)) NULL else list(Z))
  list(z = obs$z_partial_k1, z3 = obs$z_partial_k3, globals = obs$globals, ops = ops)
}
edge_sets <- function(dp, zcol = "z_part", fam = "in_family", co = "coordinated", keep = rep(TRUE, NP)) {
  f <- which(dp[[fam]] & keep)
  list(co = f[dp[[co]][f]], nc = f[!dp[[co]][f]], zd = dp[[zcol]])
}
s1_eval <- function(S, z, n_draws = 0L) {
  zd <- S$zd
  co <- S$co[is.finite(z[S$co]) & is.finite(zd[S$co])]
  nc <- S$nc[is.finite(z[S$nc]) & is.finite(zd[S$nc])]
  ag_co <- sign(z[co]) == sign(zd[co]); ag_nc <- sign(z[nc]) == sign(zd[nc])
  ncp <- nc[zd[nc] > 0]
  s1 <- mean(ag_co); pool <- mean(ag_nc); poolpos <- mean(z[ncp] > 0)
  out <- data.table(n_coord = length(co), n_noncoord = length(nc), n_noncoord_pos = length(ncp),
                    draw_overlap = length(co) / length(nc), s1 = s1, pool_concordance = pool, excess = s1 - pool,
                    pool_pos_concordance = poolpos, signmatched_excess = s1 - poolpos,
                    null_mean = NA_real_, null_sd = NA_real_, null_q95 = NA_real_, p = NA_real_,
                    n_draws = n_draws, null_feasible = length(nc) >= length(co))
  if (n_draws > 0L && out$null_feasible) {
    set.seed(SEED)
    nul <- replicate(n_draws, mean(ag_nc[sample.int(length(nc), length(co))]))
    out[, `:=`(null_mean = mean(nul), null_sd = sd(nul), null_q95 = unname(quantile(nul, 0.95)),
               p = (1 + sum(nul >= s1)) / (1 + n_draws))]
  }
  out
}
set_graph <- function(z) {   # the positive graph; column z_holdout carries this cohort's z_partial_k1
  keep <- is.finite(z) & z > 0
  ge <- data.table(program_a = pairs0$program_a[keep], program_b = pairs0$program_b[keep], z_holdout = z[keep])
  assign("g_edges", ge, envir = globalenv())
  ge
}
lin_pick <- function() unlist(lapply(names(c3_lin), function(L) {
  pool <- names(lin)[lin == L]; sample(pool, min(c3_lin[[L]], length(pool))) }))
fin_count <- function(z, set) { k <- pairs0$program_a %in% set & pairs0$program_b %in% set; sum(is.finite(z[k])) }
graph_tests <- function(z, n_draws, label) {   # S2a, S2b, S3 with the 2026-09-08 nulls; seed before each
  ge <- set_graph(z)
  s2 <- mod_of(ge, part)
  set.seed(SEED); nul_a <- replicate(n_draws, mod_of(ge, setNames(sample(part), names(part))))
  set.seed(SEED); nul_b <- replicate(n_draws, mod_of(rewire_once(ge, 10L * nrow(ge)), part))
  set_graph(z)
  s3 <- dens(c3)
  set.seed(SEED)
  nul3 <- replicate(n_draws, { pick <- lin_pick(); c(dens(pick), fin_count(z, pick)) })
  pp <- function(nul, obs) (1 + sum(nul >= obs, na.rm = TRUE)) / (1 + sum(is.finite(nul)))
  list(tab = data.table(set = label,
                        statistic = c("S2a_modularity_vs_label_permutation", "S2b_modularity_vs_rewiring", "S3_community3_density"),
                        observed = c(s2, s2, s3),
                        null_mean = c(mean(nul_a, na.rm = TRUE), mean(nul_b, na.rm = TRUE), mean(nul3[1, ])),
                        null_sd = c(sd(nul_a, na.rm = TRUE), sd(nul_b, na.rm = TRUE), sd(nul3[1, ])),
                        p = c(pp(nul_a, s2), pp(nul_b, s2), pp(nul3[1, ], s3)), n_draws = n_draws,
                        n_positive_edges = nrow(ge)),
       c3_eligible = fin_count(z, c3), null_eligible_mean = mean(nul3[2, ]))
}
s3_lin_excess <- function(z, n) {  # S3 minus the mean of n lineage-matched draws (no seed: caller's stream)
  set_graph(z); s3 <- dens(c3)
  c(s3 = s3, s3_null = mean(replicate(n, dens(lin_pick()))))
}

# ------------------------------------------- 1b. discovery blocks and rediscovery
tick("building discovery blocks")
b_coarse <- lapply(split(md, by = "dataset"), function(mc) make_block(sm, mc, COARSE))[DISC]
for (b in b_coarse) diag_add("D2", sprintf("%s n / model columns / df / programs testable", b$cohort),
                             c(b$n, b$n_par, b$df, sum(b$coverage$cohort_testable)),
                             c(GSE162694 = "109 4 105 113", GSE130970 = "65 4 61 113", GSE135251 = "197 4 193 113")[[b$cohort]])
d_ca <- approx_discovery(b_coarse)
dp <- d_ca$pairs
red <- dp[, .(ekey, program_a, program_b, entangled, in_family, z_disj, z_part, zs_disj, zs_part, q_disj, q_part, coordinated)]
fwrite(red, file.path(OUT, "frozen_rediscovery_coarse_371.tsv.gz"), sep = "\t")
writeLines(sprintf("%s  frozen_rediscovery_coarse_371.tsv.gz", sha(file.path(OUT, "frozen_rediscovery_coarse_371.tsv.gz"))),
           file.path(OUT, "frozen_rediscovery_coarse_371.tsv.gz.sha256"))
n_fam <- dp[in_family == TRUE, .N]; n_co <- dp[coordinated == TRUE, .N]
same_set <- setequal(dp[coordinated == TRUE, ekey], frozen_3811[coordinated == TRUE, ekey])
ov_frozen <- dp[coordinated == TRUE & ekey %in% frozen_coord$ekey, .N]
diag_add("D1", "family / coordinated / overlap with frozen 2,426", c(n_fam, n_co, ov_frozen), "6195 2328 1918")
diag_add("D1", "ekey set identical to job 3811 file", same_set, "TRUE")
diag_add("D1", "coordinated edges all z_part > 0", all(dp[coordinated == TRUE, z_part] > 0), "TRUE")
tick("rediscovery: family ", n_fam, " coordinated ", n_co, " identical set ", same_set)
if (!same_set || n_co != 2328L || n_fam != 6195L)
  stop_rule(2, "in-run coarse rediscovery on 371 units does not give the job 3811 2,328-edge set",
            c(sprintf("family %d, coordinated %d, overlap with job 3811 set %d, identical %s", n_fam, n_co,
                      dp[coordinated == TRUE & ekey %in% frozen_3811[coordinated == TRUE, ekey], .N], same_set)))

S_co <- edge_sets(dp)
co_ids <- S_co$co; nc_ids <- S_co$nc
mean_r_disc <- mean(tanh(dp$z_part[co_ids]))

# ------------------------------------------------------------- 2. diagnostics
diag_add("D3", "non-coordinated pool / draw overlap / discovery-positive non-coordinated",
         c(length(nc_ids), round(length(co_ids) / length(nc_ids), 4), sum(dp$z_part[nc_ids] > 0)), "3867 0.602 1110")
diag_add("D3", "coordinated cross-label / same-label / technical-involved / cross-label no technical",
         c(sum(pairs0$cross[co_ids]), sum(!pairs0$cross[co_ids]), sum(pairs0$tech_any[co_ids]),
           sum(pairs0$cross[co_ids] & !pairs0$tech_any[co_ids])), "1700 628 495 1313")
diag_add("D3", "non-coordinated cross-label no technical (positive)",
         c(sum(pairs0$cross[nc_ids] & !pairs0$tech_any[nc_ids]), sum(pairs0$cross[nc_ids] & !pairs0$tech_any[nc_ids] & dp$z_part[nc_ids] > 0)),
         "2617 692")
diag_add("D3", "mean discovery partial r coordinated / non-coordinated", c(mean_r_disc, mean(tanh(dp$z_part[nc_ids]))), "0.3538 -0.1357")
c3k <- pairs0$program_a %in% c3 & pairs0$program_b %in% c3
diag_add("D3", "community-3 pairs / entangled / coordinated", c(sum(c3k), sum(dp$entangled[c3k]), sum(dp$coordinated[c3k])), "105 2 89")
rB2 <- fread(P(file.path(rr, "jobs/3811/out/result_B2.tsv")))
diag_add("D4", "coarse vs full (348 NAS-complete): Spearman / mean|dr| / both / coarse only / full only (job 3811 result_B2)",
         c(rB2$spearman_partial_coarse_vs_full, rB2$mean_abs_diff_r, rB2$both, rB2$coarse_only, rB2$full_only), "0.9955 0.021 2088 134 75")
rC <- fread(P(file.path(rr, "jobs/3811/out/result_C.tsv")))
for (k in seq_len(nrow(rC))) diag_add("D5", sprintf("%s residual stage/NAS: n q<0.05 / n p<0.05 / median R2 / chance R2 (job 3811 result_C)", rC$cohort[k]),
                                      c(rC$n_q05[k], rC$n_p05[k], rC$median_r2[k], rC$chance_r2[k]), "")
for (cc in c(DISC, "GSE213621")) {
  x <- if (cc == "GSE213621") g else md[dataset == cc]
  q3 <- function(v) sprintf("%.3f [%.3f-%.3f]", median(v), quantile(v, .25), quantile(v, .75))
  diag_add("D6", sprintf("%s log10 library / entropy / top-100 share: median [IQR]", cc),
           paste(q3(x$log_lib), q3(x$entropy), q3(x$top100), sep = " | "),
           c(GSE162694 = "7.437 6.595 0.520", GSE130970 = "7.504 6.024 0.573", GSE135251 = "7.479 6.318 0.581", GSE213621 = "")[[cc]])
}
# D11 discovery part
for (b in b_coarse) {
  gl <- rep_pass(b)$globals
  diag_add("D11", sprintf("%s alone: positive_fraction_partial_k1 / global_component_share", b$cohort),
           c(gl[["positive_fraction_partial_k1"]], gl[["global_component_share"]]), "")
}
diag_add("D11", "371-unit rediscovery: positive fraction (family) / positive_fraction_partial_k1 (eligible) / global_component_share",
         c(dp[in_family == TRUE, mean(z_part > 0)], d_ca$globals[["positive_fraction_partial_k1"]], d_ca$globals[["global_component_share"]]),
         "0.555 (family)")
# D12 Louvain on the 371-unit coarse coordinated graph
set.seed(SEED)
cm <- t4_modularity(dp, ifelse(dp$coordinated, tanh(dp$z_part), NA_real_), features)
xt <- merge(cm$membership, communities[, .(feature_id, frozen = community)], by = "feature_id")
xtab <- table(new = xt$community, frozen = xt$frozen)
print(xtab)
diag_add("D12", "Louvain (seed 20260923) on 371 coordinated graph: modularity / communities", c(cm$modularity, cm$n_communities), "about 0.358 / 3 at seed 20261015")
for (r in rownames(xtab)) diag_add("D12", sprintf("new community %s: counts of frozen communities 1/2/3/4", r), as.integer(xtab[r, ]), "")
diag_add("D12", "frozen c3 programs in the new community holding most of them", max(table(xt[frozen == 3]$community)), "10 at seed 20261015")
diag_write()

# ------------------------------------------------------------- 3. PC1, PC2
tick("positive controls")
pc_one <- function(h, label) {
  d <- approx_discovery(b_coarse[setdiff(DISC, h)])$pairs
  r <- rep_pass(b_coarse[[h]])
  s <- s1_eval(edge_sets(d), r$z, N_DRAW)
  gt <- graph_tests(r$z, N_DRAW, label)
  list(s1 = cbind(data.table(control = label, held_out = h, n_test = b_coarse[[h]]$n, n_coordinated_disc = d[coordinated == TRUE, .N]), s),
       graph = cbind(data.table(control = label, held_out = h), gt$tab))
}
pc1 <- pc_one("GSE162694", "PC1")
pc2a <- pc_one("GSE130970", "PC2")
pc2b <- pc_one("GSE135251", "PC2")
pc_s1 <- rbind(pc1$s1, pc2a$s1, pc2b$s1)
pc_gr <- rbind(pc1$graph, pc2a$graph, pc2b$graph)
pc1_s1 <- pc1$s1
pc1_gp <- setNames(pc1$graph$p, c("S2a", "S2b", "S3"))
pc1_det <- pc1_s1$n_coordinated_disc == 2033L && abs(pc1_s1$s1 - 0.936055090998524) < 5e-5
pc1_pass <- pc1_det && pc1_s1$excess >= 0.05 && pc1_s1$p <= 0.0125 && all(pc1_gp <= 0.05)
pc_out <- rbind(pc_s1[, .(control, held_out, n_test, statistic = "S1", n_coordinated_disc, observed = s1, null_mean, null_sd, p,
                          excess, signmatched_excess, pool_concordance, pool_pos_concordance, n_draws)],
                pc_gr[, .(control, held_out, n_test = NA_integer_, statistic, n_coordinated_disc = NA_integer_, observed, null_mean, null_sd, p,
                          excess = NA_real_, signmatched_excess = NA_real_, pool_concordance = NA_real_, pool_pos_concordance = NA_real_, n_draws)])
pc_out[, pass_rule_PC1 := ifelse(control == "PC1", pc1_pass, NA)]
fwrite(pc_out, file.path(OUT, "pc_results.tsv"), sep = "\t")
print(pc_out)
diag_add("D7", "PC1 coordinated / S1 / excess / S1 p / S2a p / S2b p / S3 p / pass",
         c(pc1_s1$n_coordinated_disc, f4(pc1_s1$s1), f4(pc1_s1$excess), pc1_s1$p, pc1_gp, pc1_pass), "2033 0.9361 ~0.187 ...")
diag_add("D7", "PC2 S1 (GSE130970 held out / GSE135251 held out)", c(f4(pc2a$s1$s1), f4(pc2b$s1$s1)), "0.764 0.888")
diag_add("D7", "PC2 excess (GSE130970 / GSE135251)", c(f4(pc2a$s1$excess), f4(pc2b$s1$excess)), "0.071 0.151")
diag_write()
if (!pc1_pass)
  stop_rule(3, "positive control PC1 fails its rule",
            c(sprintf("PC1 coordinated %d (expected 2033), S1 %.6f (expected 0.936055), excess %.4f, S1 p %.5f, S2a p %.5f, S2b p %.5f, S3 p %.5f",
                      pc1_s1$n_coordinated_disc, pc1_s1$s1, pc1_s1$excess, pc1_s1$p, pc1_gp[1], pc1_gp[2], pc1_gp[3])))
if (STOP_AFTER == "pc1") { tick("H52_STOP_AFTER=pc1: ending debug run"); quit(save = "no", status = 0) }

# ----------------------------------------------------------- 4. GSE213621 block
tick("GSE213621 block")
bg <- make_block(sm, g, COARSE)
g_rank <- qr(bg$X)$rank; g_test <- sum(bg$coverage$cohort_testable); g_mincov <- min(bg$coverage$cohort_weight_coverage)
diag_add("D8", "GSE213621 n / columns / rank / df / programs testable / min weight coverage",
         c(bg$n, bg$n_par, g_rank, bg$df, g_test, round(g_mincov, 4)), "258 4 4 254 113 0.821")
g1_ok <- bg$n >= 200 && all(table(g$fibrosis_stage) >= 30) && g_rank == ncol(bg$X) && bg$df >= 150 && g_test >= 100
if (!g1_ok) stop_rule(4, "G1 fails at run time", c(sprintf("n %d, bins %s, rank %d of %d, df %d, testable %d", bg$n,
                                                           paste(table(g$fibrosis_stage), collapse = "/"), g_rank, ncol(bg$X), bg$df, g_test)))
RG <- rep_pass(bg)
zg <- RG$z; ops_g <- RG$ops
ncpos_ids <- nc_ids[dp$z_part[nc_ids] > 0]
diag_add("D8", "finite GSE213621 value: coordinated / non-coordinated / discovery-positive non-coordinated",
         c(sum(is.finite(zg[co_ids])), sum(is.finite(zg[nc_ids])), sum(is.finite(zg[ncpos_ids]))), "2328 3867 1110")
g <- add_proxies(g)
for (n in names(MK)) diag_add("D9", sprintf("marker rows found: %s", n), length(mk_rows[[n]]), c(chol = 6, fib = 4, endo = 2, bcell = 2)[[n]])
for (cc in c(DISC, "GSE213621")) {
  x <- if (cc == "GSE213621") g else md[dataset == cc]
  diag_add("D9", sprintf("%s Spearman with fibrosis bin: cholangiocyte proxy / fibroblast proxy", cc),
           c(cor(x$prx_chol, x$fibrosis_stage, method = "spearman"), cor(x$prx_fib, x$fibrosis_stage, method = "spearman")),
           if (cc == "GSE213621") "" else "0.42 to 0.56")
}
diag_add("D11", "GSE213621: positive_fraction_partial_k1 / global_component_share",
         c(RG$globals[["positive_fraction_partial_k1"]], RG$globals[["global_component_share"]]), "")
diag_add("D13", "community-3 pairs eligible (finite) in GSE213621", fin_count(zg, c3), "discovery 103 of 105")
diag_write()

# --------------------------------------------------------------- 5. NC1, NC2
tick("negative controls")
nc_run <- function(label, permute) {
  set.seed(SEED)
  rbindlist(lapply(seq_len(N_NC), function(i) {
    Zp <- permute(bg$Z)
    r <- rep_pass(bg, Z = Zp, ops = ops_g)
    s <- s1_eval(S_co, r$z)
    ge <- set_graph(r$z)
    s3 <- s3_lin_excess(r$z, N_LIN)
    data.table(control = label, perm = i, s1 = s$s1, excess = s$excess, signmatched_excess = s$signmatched_excess,
               s2 = mod_of(ge, part), s3 = s3[["s3"]], s3_null = s3[["s3_null"]], s3_minus_null = s3[["s3"]] - s3[["s3_null"]],
               positive_fraction = r$globals[["positive_fraction_partial_k1"]])
  }))
}
grp_g <- t4_gene_groups(list(bg), weights, features)[[1L]]
nc_draws <- rbind(nc_run("NC1_gene_permutation", t4_permute_genes),
                  nc_run("NC2_program_group_permutation", function(Z) t4_permute_groups(Z, grp_g)))
fwrite(nc_draws, file.path(OUT, "nc_draws.tsv"), sep = "\t")
nc_sum <- nc_draws[, .(n_perm = .N, mean_s1 = mean(s1), mean_excess = mean(excess), sd_excess = sd(excess),
                       q95_excess = unname(quantile(excess, 0.95)), max_excess = max(excess),
                       mean_signmatched_excess = mean(signmatched_excess), mean_s2 = mean(s2, na.rm = TRUE),
                       mean_s3_minus_null = mean(s3_minus_null)), by = control]
nc_sum[, fails := mean_excess > 0.02 | q95_excess >= 0.05 | mean_signmatched_excess > 0.02]
fwrite(nc_sum, file.path(OUT, "nc_results.tsv"), sep = "\t")
print(nc_sum)
for (k in seq_len(nrow(nc_sum))) diag_add("D10", sprintf("%s mean S1 / mean, SD, q95, max excess / mean sign-matched / mean S2 / mean S3-null / fails", nc_sum$control[k]),
                                         nc_sum[k, c(f4(mean_s1), f4(mean_excess), f4(sd_excess), f4(q95_excess), f4(max_excess),
                                                     f4(mean_signmatched_excess), f4(mean_s2), f4(mean_s3_minus_null), fails)],
                                         "S1 0.48-0.50; excess -0.03..+0.01; q95 < 0.04; S2 +/-0.02; S3-null +/-0.10")
diag_write()
if (any(nc_sum$fails))
  stop_rule(5, "negative control NC1 or NC2 fails its rule; result not interpretable",
            nc_sum[, sprintf("%s: mean excess %.4f, q95 %.4f, mean sign-matched %.4f", control, mean_excess, q95_excess, mean_signmatched_excess)])
n_nc_fin <- sum(is.finite(zg[nc_ids]) & is.finite(dp$z_part[nc_ids])); n_co_fin <- sum(is.finite(zg[co_ids]))
if (n_nc_fin < n_co_fin) stop_rule(6, "fewer finite non-coordinated edges than coordinated edges",
                                   sprintf("non-coordinated finite %d, coordinated finite %d", n_nc_fin, n_co_fin))
if (STOP_AFTER == "nc") { tick("H52_STOP_AFTER=nc: ending debug run"); quit(save = "no", status = 0) }

# ------------------------------------------------------------- 6. F1 primary
tick("F1 primary statistics")
f1_s1 <- s1_eval(S_co, zg, N_DRAW)
f1_gr <- graph_tests(zg, N_DRAW, "GSE213621")
diag_add("D13", "mean eligible community-3-size pairs over the S3 null sets", f1_gr$null_eligible_mean, "discovery random sets 103.6")
diag_write()
F1 <- data.table(statistic = c("S1_edge_sign_concordance", f1_gr$tab$statistic),
                 observed = c(f1_s1$s1, f1_gr$tab$observed),
                 null_mean = c(f1_s1$null_mean, f1_gr$tab$null_mean),
                 null_sd = c(f1_s1$null_sd, f1_gr$tab$null_sd),
                 p = c(f1_s1$p, f1_gr$tab$p), n_draws = N_DRAW, p_floor = 1 / (N_DRAW + 1))
F1[, p_holm := p.adjust(p, "holm")]
F1[, holm_rejects := p_holm <= 0.05]
print(F1)

# bootstrap index matrix: within fibrosis bin, drawn once, saved before scoring
set.seed(SEED)
bins <- g$fibrosis_stage
BIDX <- t(vapply(seq_len(N_BOOT), function(b) unlist(lapply(1:3, function(k) {
  pos <- which(bins == k); pos[sample.int(length(pos), length(pos), replace = TRUE)] })), integer(nrow(g))))
fwrite(data.table(b = seq_len(N_BOOT), BIDX), file.path(OUT, "bootstrap_indices.tsv.gz"), sep = "\t")
fwrite(data.table(position = seq_len(nrow(g)), sample_id = g$sample_id, analysis_unit_id = g$analysis_unit_id, fibrosis_bin = bins),
       file.path(OUT, "bootstrap_unit_order.tsv"), sep = "\t")
tick("bootstrap indices written: ", N_BOOT, " x ", ncol(BIDX))
resample_meta <- function(b) {
  mb <- copy(g[BIDX[b, ]]); mb[, orig_id := sample_id]; mb[, sample_id := paste0(sample_id, "_b", seq_len(.N))]
  mb
}
resample_mat <- function(mb) { x <- sm[, mb$orig_id, drop = FALSE]; colnames(x) <- mb$sample_id; x }

ZB <- matrix(NA_real_, N_BOOT, NP)   # coarse GSE213621 z per resample, reused by F2c, F2d
set.seed(SEED)                        # lineage-matched draws inside the bootstrap
boot_rows <- vector("list", N_BOOT)
for (b in seq_len(N_BOOT)) {
  mb <- resample_meta(b)
  blk <- tryCatch(make_block(resample_mat(mb), mb, COARSE), error = function(e) NULL)
  if (is.null(blk)) { boot_rows[[b]] <- data.table(b = b, ok = FALSE); next }
  r <- rep_pass(blk); ZB[b, ] <- r$z
  s <- s1_eval(S_co, r$z)
  ge <- set_graph(r$z); s2 <- mod_of(ge, part)
  s3 <- s3_lin_excess(r$z, N_LIN)
  boot_rows[[b]] <- data.table(b = b, ok = TRUE, s1 = s$s1, excess = s$excess, signmatched_excess = s$signmatched_excess, s2 = s2,
                               s3 = s3[["s3"]], s3_null = s3[["s3_null"]], s3_excess = s3[["s3"]] - s3[["s3_null"]],
                               transport_ratio = mean(tanh(r$z[co_ids])) / mean_r_disc)
  if (b %% 100L == 0L) tick("bootstrap ", b)
}
boot <- rbindlist(boot_rows, fill = TRUE)
fwrite(boot, file.path(OUT, "f1_bootstrap.tsv"), sep = "\t")
tick("bootstrap done: ", sum(boot$ok), " of ", N_BOOT, " resamples scored")
ci <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))
CI <- lapply(boot[ok == TRUE, .(s1, excess, signmatched_excess, s2, s3_excess, transport_ratio)], ci)

obs_s2 <- f1_gr$tab$observed[1]; obs_s3 <- f1_gr$tab$observed[3]
s2a_null <- f1_gr$tab$null_mean[1]; s2b_null <- f1_gr$tab$null_mean[2]
F1[, observed_effect := c(f1_s1$excess, obs_s2, obs_s2, obs_s3 - f1_gr$tab$null_mean[3])]
F1[, effect := c("S1 minus whole non-coordinated pool concordance", "S2", "S2", "S3 minus lineage-null mean")]
F1[, `:=`(boot_lo = c(CI$excess[1], CI$s2[1], CI$s2[1], CI$s3_excess[1]),
          boot_hi = c(CI$excess[2], CI$s2[2], CI$s2[2], CI$s3_excess[2]))]
F1[, unit_condition_threshold := c(0, s2a_null, s2b_null, 0)]
F1[, unit_condition := boot_lo > unit_condition_threshold]
F1[, passes := holm_rejects & unit_condition]
fwrite(F1, file.path(OUT, "f1_results.tsv"), sep = "\t")
print(F1)

pass <- setNames(F1$passes, c("S1", "S2a", "S2b", "S3"))
verdict <- if (all(pass)) "Replicated" else if (pass[["S1"]] && pass[["S2a"]] && pass[["S2b"]]) "Partially replicated" else "Not replicated"
LB <- CI$excess[1]; UB <- CI$excess[2]
qualifier <- if (UB < 0.05) "Excluded" else if (LB >= 0.05) "Clear" else if (LB > 0) "Small" else "Indeterminate"
sm_LB <- CI$signmatched_excess[1]
edge_specific <- sm_LB > 0

# ------------------------------------------------------------- 7. F2 arms
tick("F2 arms")
arm_boot <- function(Sarm, formula, prox) {   # paired bootstrap on the saved indices; no new random draws
  rbindlist(lapply(seq_len(N_BOOT), function(b) {
    mb <- resample_meta(b)
    if (prox) mb <- add_proxies(mb, mb$orig_id)
    blk <- tryCatch(make_block(resample_mat(mb), mb, formula), error = function(e) NULL)
    if (is.null(blk)) return(data.table(b = b, ok = FALSE))
    r <- rep_pass(blk); s <- s1_eval(Sarm, r$z)
    data.table(b = b, ok = TRUE, s1 = s$s1, excess = s$excess, s2 = mod_of(set_graph(r$z), part))
  }), fill = TRUE)
}
# F2a composition proxies: rediscovery on 371 written before the GSE213621 block
b_comp <- lapply(split(md, by = "dataset"), function(mc) make_block(sm, mc, COMP))[DISC]
dA <- approx_discovery(b_comp)$pairs
fwrite(dA[, .(ekey, program_a, program_b, entangled, in_family, z_part, zs_part, q_disj, q_part, coordinated)],
       file.path(OUT, "f2a_rediscovery_371.tsv.gz"), sep = "\t")
S_A <- edge_sets(dA)
bgA <- make_block(sm, g, COMP); rA <- rep_pass(bgA)
f2a <- s1_eval(S_A, rA$z, N_DRAW)
f2a_s2 <- mod_of(set_graph(rA$z), part); f2a_s3 <- dens(c3)
bootA <- arm_boot(S_A, COMP, TRUE)
# F2b library covariates
b_lib <- lapply(split(md, by = "dataset"), function(mc) make_block(sm, mc, LIB))[DISC]
dB <- approx_discovery(b_lib)$pairs
fwrite(dB[, .(ekey, program_a, program_b, entangled, in_family, z_part, zs_part, q_disj, q_part, coordinated)],
       file.path(OUT, "f2b_rediscovery_371.tsv.gz"), sep = "\t")
S_B <- edge_sets(dB)
bgB <- make_block(sm, g, LIB); rB <- rep_pass(bgB)
f2b <- s1_eval(S_B, rB$z, N_DRAW)
bootB <- arm_boot(S_B, LIB, FALSE)
# F2c cross-label edges without technical programs (coarse)
S_C <- edge_sets(dp, keep = pairs0$cross & !pairs0$tech_any)
f2c <- s1_eval(S_C, zg, N_DRAW)
f2c_fallback <- !f2c$null_feasible
# F2d frozen 2,426 edges, frozen saturated-design signs, coarse GSE213621
fz <- frozen_all[match(pairs0$ekey, ekey)]
stopifnot(!anyNA(fz$ekey))
dD <- data.table(in_family = fz$membership_entangled == FALSE | is.na(fz$membership_entangled),
                 coordinated = pairs0$ekey %in% frozen_coord$ekey, z_part = atanh(fz$r_partial))
S_D <- edge_sets(dD)
f2d <- s1_eval(S_D, zg, N_DRAW)
diag_add("F2", "F2a coordinated (371, composition design) / F2b coordinated (371, library design)",
         c(dA[coordinated == TRUE, .N], dB[coordinated == TRUE, .N]), "")
diag_add("F2", "F2c coordinated / pool; F2d coordinated / pool", c(length(S_C$co), length(S_C$nc), length(S_D$co), length(S_D$nc)),
         "1313 2617 2426 3769")
okb <- boot$ok
bootCD <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
  if (!okb[b]) return(data.table(b = b, f2c_excess = NA_real_, f2d_excess = NA_real_, f2c_s1 = NA_real_, f2d_s1 = NA_real_))
  c_ <- s1_eval(S_C, ZB[b, ]); d_ <- s1_eval(S_D, ZB[b, ])
  data.table(b = b, f2c_s1 = c_$s1, f2c_excess = c_$excess, f2d_s1 = d_$s1, f2d_excess = d_$excess)
}))
f2boot <- data.table(b = seq_len(N_BOOT), coarse_excess = boot$excess, coarse_s2 = boot$s2,
                     f2a_s1 = bootA$s1, f2a_excess = bootA$excess, f2a_s2 = bootA$s2,
                     f2b_s1 = bootB$s1, f2b_excess = bootB$excess,
                     f2c_s1 = bootCD$f2c_s1, f2c_excess = bootCD$f2c_excess, f2d_s1 = bootCD$f2d_s1, f2d_excess = bootCD$f2d_excess)
f2boot[, `:=`(ratio_f2a = f2a_excess / coarse_excess, ratio_f2b = f2b_excess / coarse_excess,
              diff_f2a = f2a_excess - coarse_excess, diff_f2b = f2b_excess - coarse_excess, s2_drop_f2a = coarse_s2 - f2a_s2)]
fwrite(f2boot, file.path(OUT, "f2_bootstrap.tsv"), sep = "\t")
F2 <- rbind(cbind(data.table(arm = "F2a_composition_proxies"), f2a), cbind(data.table(arm = "F2b_library_covariates"), f2b),
            cbind(data.table(arm = "F2c_cross_label_no_technical"), f2c), cbind(data.table(arm = "F2d_frozen_2426_edges"), f2d))
F2[, `:=`(excess_lo = c(ci(f2boot$f2a_excess)[1], ci(f2boot$f2b_excess)[1], ci(f2boot$f2c_excess)[1], ci(f2boot$f2d_excess)[1]),
          excess_hi = c(ci(f2boot$f2a_excess)[2], ci(f2boot$f2b_excess)[2], ci(f2boot$f2c_excess)[2], ci(f2boot$f2d_excess)[2]))]
in_holm <- c(TRUE, TRUE, !f2c_fallback, TRUE)
F2[, in_holm_family := in_holm]
F2[in_holm_family == TRUE, p_holm := p.adjust(p, "holm")]
F2[in_holm_family == FALSE, p := NA_real_]
F2[, holm_rejects := p_holm <= 0.05]
r_a <- ci(f2boot$ratio_f2a); r_b <- ci(f2boot$ratio_f2b)
F2[, `:=`(ratio_to_coarse_lo = c(r_a[1], r_b[1], NA, NA), ratio_to_coarse_hi = c(r_a[2], r_b[2], NA, NA))]
F2[, n_coordinated_disc := c(dA[coordinated == TRUE, .N], dB[coordinated == TRUE, .N], length(S_C$co), length(S_D$co))]
fwrite(F2, file.path(OUT, "f2_results.tsv"), sep = "\t")
print(F2)
f2a_notexplained <- F2$excess_lo[1] >= 0.05
f2a_tracks <- r_a[2] < 0.5
f2b_triggered <- r_b[2] < 0.5
f2c_zero_in <- F2$excess_lo[3] <= 0 && F2$excess_hi[3] >= 0

# ------------------------------------------------------- 8. descriptive outputs
tick("descriptive outputs")
DS <- list()
ds_add <- function(item, value, lo = NA_real_, hi = NA_real_, note = "") DS[[length(DS) + 1L]] <<-
  data.table(item = item, value = value, boot_lo = lo, boot_hi = hi, note = note)
# lineage-stratified S2 null
set.seed(SEED)
ge <- set_graph(zg)
nul_ls <- replicate(N_DRAW, { p2 <- part; for (L in unique(lin)) { ix <- which(lin == L); p2[ix] <- p2[ix][sample.int(length(ix))] }; mod_of(ge, p2) })
ds_add("S2 lineage-stratified label-permutation null mean", mean(nul_ls), note = "LOCO null means -0.0027/-0.0012/-0.0027")
ds_add("S2 vs lineage-stratified null: nominal p", (1 + sum(nul_ls >= obs_s2)) / (1 + N_DRAW), note = "nominal")
ds_add("transport ratio (mean GSE213621 r over 2,328 / mean discovery r)", mean(tanh(zg[co_ids])) / mean_r_disc,
       CI$transport_ratio[1], CI$transport_ratio[2], sprintf("denominator %.4f (prespec 0.3538)", mean_r_disc))
for (cl in list(list("same_label", !pairs0$cross), list("cross_label", pairs0$cross), list("technical_involved", pairs0$tech_any),
                list("no_technical", !pairs0$tech_any))) {
  s <- s1_eval(edge_sets(dp, keep = cl[[2]]), zg)
  ds_add(sprintf("S1 by edge class: %s (n coord %d)", cl[[1]], s$n_coord), s$s1, note = sprintf("class pool concordance %.4f; excess %.4f", s$pool_concordance, s$excess))
}
for (k in 1:4) {
  ck <- pairs0$ca == k & pairs0$cb == k
  s <- s1_eval(edge_sets(dp, keep = ck), zg)
  set_graph(zg)
  ds_add(sprintf("community %d: S1 over within-community coordinated edges (n %d)", k, s$n_coord), s$s1)
  ds_add(sprintf("community %d: within-community density", k), dens(names(part)[part == k]))
}
# k = 3 sensitivity
S_k3 <- edge_sets(dp, zcol = "z_part3", fam = "in_family3", co = "coordinated3")
s_k3 <- s1_eval(S_k3, RG$z3, N_DRAW)
ds_add(sprintf("k=3: coordinated edges in 371 rediscovery"), length(S_k3$co))
ds_add("k=3: S1 in GSE213621", s_k3$s1, note = sprintf("null mean %.4f, nominal p %.5f, excess %.4f, sign-matched %.4f",
                                                      s_k3$null_mean, s_k3$p, s_k3$excess, s_k3$signmatched_excess))
# re-derived partition (description only)
set.seed(SEED)
cmg <- t4_modularity(pairs0, zg, features)
mem <- setNames(cmg$membership$community, cmg$membership$feature_id)[features]
ari <- igraph::compare(as.integer(factor(mem)), as.integer(part[features]), method = "adjusted.rand")
tc3 <- table(mem[c3]); best <- as.integer(names(tc3)[which.max(tc3)])
ds_add("re-derived Louvain on GSE213621 positive graph: modularity", cmg$modularity, note = sprintf("%d communities", cmg$n_communities))
ds_add("re-derived partition: adjusted Rand index vs frozen", ari)
ds_add(sprintf("re-derived community %d holds most frozen-c3 programs: count of 15", best), max(tc3),
       note = paste("lineage counts:", paste(names(table(lin[names(mem)[mem == best]])), table(lin[names(mem)[mem == best]]), collapse = ", ")))
# accession-order sensitivity
g[, accnum := as.numeric(sub("^[A-Z]+", "", sample_id))]
stopifnot(!anyNA(g$accnum))
g[, acc_q := factor(cut(accnum, quantile(accnum, 0:4 / 4), include.lowest = TRUE, labels = FALSE))]
bgQ <- make_block(sm, g, ACC); rQ <- rep_pass(bgQ)
s_q <- s1_eval(S_co, rQ$z, N_DRAW)
ds_add("accession-quartile design (GSE213621 only): S1", s_q$s1, note = sprintf("null mean %.4f, nominal p %.5f, excess %.4f", s_q$null_mean, s_q$p, s_q$excess))
ds_add("F2a minus coarse excess (paired)", f2a$excess - f1_s1$excess, ci(f2boot$diff_f2a)[1], ci(f2boot$diff_f2a)[2])
ds_add("F2b minus coarse excess (paired)", f2b$excess - f1_s1$excess, ci(f2boot$diff_f2b)[1], ci(f2boot$diff_f2b)[2])
ds_add("S2 drop under F2a (coarse S2 minus F2a S2, paired)", obs_s2 - f2a_s2, ci(f2boot$s2_drop_f2a)[1], ci(f2boot$s2_drop_f2a)[2])
ds_add("F2a S2 (frozen partition modularity, composition design)", f2a_s2)
ds_add("F2a S3 (community-3 density, composition design)", f2a_s3)
hist193 <- fread(P(file.path(hold, "results/A_coordination_results.tsv")))
for (k in seq_len(nrow(hist193))) ds_add(sprintf("GSE193066 historical (identity-unqualified): %s", hist193$statistic[k]), hist193$observed[k],
                                         note = sprintf("null mean %.4f, p %.5f", hist193$null_mean[k], hist193$p_value[k]))
DESC <- rbindlist(DS)
fwrite(DESC, file.path(OUT, "descriptive.tsv"), sep = "\t")
print(DESC)

# edge-level table
ec <- pairs0[, .(ekey, program_a, program_b, lineage_a = la, lineage_b = lb, community_a = ca, community_b = cb, cross_label = cross, technical_involved = tech_any)]
ec[, `:=`(entangled_disc = dp$entangled, in_family = dp$in_family, coordinated = dp$coordinated, z_partial_disc = dp$z_part,
          zs_partial_disc = dp$zs_part, z_partial_k1_GSE213621 = zg, z_partial_k3_GSE213621 = RG$z3)]
ec[, sign_agree := ifelse(is.finite(z_partial_k1_GSE213621) & is.finite(z_partial_disc), sign(z_partial_k1_GSE213621) == sign(z_partial_disc), NA)]
fwrite(ec, file.path(OUT, "edge_concordance.tsv.gz"), sep = "\t")
diag_write()

# --------------------------------------------------------------- RESULT_H52.md
bench <- fread(P(file.path(rr, "scratch/plan-H52/benchmark_s1.tsv")))[df == 254]
title_support <- verdict == "Replicated" && qualifier == "Clear" && edge_specific && !f2b_triggered
iv <- function(a, b) sprintf("%.4f to %.4f", a, b)
f1row <- function(i) sprintf("| %s | %.4f | %.4f | %.5f | %.5f | %s | %s | %s | %s |", F1$statistic[i], F1$observed[i], F1$null_mean[i],
                             F1$p[i], F1$p_holm[i], F1$holm_rejects[i], iv(F1$boot_lo[i], F1$boot_hi[i]), F1$unit_condition[i], F1$passes[i])
txt <- c(
  "# RESULT H52", "",
  if (DEBUG) "**DEBUG RUN on permuted GSE213621 data and reduced draw counts. Not a result.**" else NULL, "",
  sprintf("Prespecification: PRESPEC_H52.md sha256 %s (verified). Inputs: %d sha256 verified. Seed %d.", EXPECTED_PRESPEC, nrow(INPUTS), SEED),
  sprintf("Units: 371 discovery disease units (GSE130970 65, GSE135251 197, GSE162694 109); 258 GSE213621 disease units (bins 87/94/77). Design ~ sex_final + fibrosis bin (3 levels)."),
  sprintf("Draws: %d per null; unit bootstrap B = %d (within fibrosis bin); NC %d permutations each.", N_DRAW, N_BOOT, N_NC), "",
  "## Stop rules", "",
  "All six stop rules checked in order; none fired.",
  sprintf("Stop rule 2: rediscovery family %d, coordinated %d, identical to the job 3811 set: %s.", n_fam, n_co, same_set),
  sprintf("PC1 (GSE162694 held out): coordinated %d, S1 %.4f, excess %.4f, S1 p %.5f, S2a p %.5f, S2b p %.5f, S3 p %.5f: %s.",
          pc1_s1$n_coordinated_disc, pc1_s1$s1, pc1_s1$excess, pc1_s1$p, pc1_gp[1], pc1_gp[2], pc1_gp[3], if (pc1_pass) "pass" else "fail"),
  sprintf("NC1: mean excess %.4f, q95 %.4f, mean sign-matched %.4f. NC2: mean excess %.4f, q95 %.4f, mean sign-matched %.4f. Both pass.",
          nc_sum$mean_excess[1], nc_sum$q95_excess[1], nc_sum$mean_signmatched_excess[1], nc_sum$mean_excess[2], nc_sum$q95_excess[2], nc_sum$mean_signmatched_excess[2]), "",
  "## Family F1 (GSE213621; Holm over 4; p floor 1/1001)", "",
  "| Statistic | Observed | Null mean | p | Holm p | Holm rejects | 95% unit-bootstrap interval of the effect | Unit condition | Passes |",
  "|---|---:|---:|---:|---:|---|---|---|---|",
  vapply(seq_len(nrow(F1)), f1row, character(1)), "",
  "Effect intervals: S1 row = mixed-pool excess; S2 rows = S2 value (condition: lower bound above the observed null mean); S3 row = density minus per-resample lineage-null mean.", "",
  sprintf("- S1 = %.4f (interval %s); whole-pool concordance %.4f; edge-draw null mean %.4f.", f1_s1$s1, iv(CI$s1[1], CI$s1[2]), f1_s1$pool_concordance, f1_s1$null_mean),
  sprintf("- Mixed-pool excess = %.4f (interval %s).", f1_s1$excess, iv(LB, UB)),
  sprintf("- Sign-matched excess = %.4f (interval %s); concordance of the %d discovery-positive non-coordinated edges %.4f.",
          f1_s1$signmatched_excess, iv(CI$signmatched_excess[1], CI$signmatched_excess[2]), f1_s1$n_noncoord_pos, f1_s1$pool_pos_concordance),
  sprintf("- S2 = %.4f (interval %s); label null mean %.4f; rewiring null mean %.4f.", obs_s2, iv(CI$s2[1], CI$s2[2]), s2a_null, s2b_null),
  sprintf("- S3 = %.4f vs lineage null %.4f; excess interval %s. Community-3 pairs eligible in GSE213621: %d of 105.", obs_s3, f1_gr$tab$null_mean[3],
          iv(CI$s3_excess[1], CI$s3_excess[2]), f1_gr$c3_eligible),
  sprintf("- Benchmarks at df 254 (predicted S1): c = 1 %.3f; c = 0.75 %.3f; c = 0.5 %.3f; c = 0.25 %.3f.", bench$s1_coord[1], bench$s1_coord[2], bench$s1_coord[3], bench$s1_coord[4]),
  sprintf("- Transport ratio %.4f (interval %s).", mean(tanh(zg[co_ids])) / mean_r_disc, iv(CI$transport_ratio[1], CI$transport_ratio[2])), "",
  "## Verdict (computed by the fixed rule)", "",
  sprintf("**Verdict: %s.** Effect qualifier: **%s** (mixed-pool excess interval %s)%s.", verdict, qualifier, iv(LB, UB),
          if (qualifier == "Excluded") sprintf("; LB > 0: %s", LB > 0) else ""),
  sprintf("Edge-specificity check: sign-matched excess LB %.4f; %s.", sm_LB,
          if (edge_specific) "coordinated edges keep their sign more often than weaker positive edges" else
            "LB <= 0, the S1 result is consistent with a cohort-wide shift in edge signs"),
  sprintf("F2b library rule (ratio UB < 0.5): ratio interval %s; triggered: %s.", iv(r_b[1], r_b[2]), f2b_triggered),
  sprintf("Title-supporting reading (Replicated, Clear, sign-matched LB > 0, F2b not triggered): **%s**.", if (title_support) "holds" else "does not hold"),
  "The joint H52/H54 title-clause rule is not applied here.", "",
  "## Family F2 (S1 sensitivity; Holm over the arms in the family)", "",
  "| Arm | Coordinated (disc) | S1 | Null mean | Excess | Excess interval | p | Holm p | Holm rejects |", "|---|---:|---:|---:|---:|---|---:|---:|---|",
  sprintf("| %s | %d | %.4f | %.4f | %.4f | %s | %.5f | %.5f | %s |", F2$arm, F2$n_coordinated_disc, F2$s1, F2$null_mean, F2$excess,
          iv(F2$excess_lo, F2$excess_hi), F2$p, F2$p_holm, F2$holm_rejects), "",
  sprintf("- Composition (F2a): excess LB >= 0.05: %s; F2a/coarse ratio interval %s; ratio UB < 0.5 (tracks proxies): %s; S2 drop %.4f (interval %s).",
          f2a_notexplained, iv(r_a[1], r_a[2]), f2a_tracks, obs_s2 - f2a_s2, iv(ci(f2boot$s2_drop_f2a)[1], ci(f2boot$s2_drop_f2a)[2])),
  sprintf("- Library (F2b): ratio interval %s; rule triggered: %s.", iv(r_b[1], r_b[2]), f2b_triggered),
  sprintf("- Edge class (F2c): excess interval includes 0: %s%s.", f2c_zero_in, if (f2c_fallback) " (fallback: pool too small, no p)" else ""),
  "", "Full tables: diagnostics.tsv, pc_results.tsv, nc_results.tsv, f1_results.tsv, f1_bootstrap.tsv, f2_results.tsv, f2_bootstrap.tsv, descriptive.tsv, edge_concordance.tsv.gz.",
  "", sprintf("Written %s by run_h52.R.", format(Sys.time(), tz = "UTC", usetz = TRUE)))
writeLines(txt, file.path(OUT, "RESULT_H52.md"))
cat(txt, sep = "\n")
tick("DONE")
