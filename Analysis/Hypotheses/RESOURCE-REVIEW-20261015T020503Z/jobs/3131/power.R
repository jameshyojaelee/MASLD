#!/usr/bin/env Rscript
# Gate L07-01 feasibility job. Uses the already-opened GSE193066 holdout only.
# (1) code-path positive control: reproduce holdout S1 = 0.8974 (n = 105).
# (2) observed S1 in the holdout F0-F1 stratum (6 F0 + 37 F1 = 43 participants).
# (3) participant-level sampling SD of S1(coordinated), S1(non-coordinated) and
#     their difference at the candidate's n, by m-out-of-n bootstrap of participants.
# The T4 library is sourced read-only from the repo; nothing in the repo is written.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(Matrix)})
options(digits = 10)
B <- as.integer(Sys.getenv("GATE_B", "100"))
SEED <- 20261015L
script_dir <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/program_observability_map"
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "t4_interaction_lib.R"))
out <- getwd()
t0 <- Sys.time()

T4 <- normalizePath(file.path(sealed_input_root, "..", "systems", "T4_interaction"))
communities <- fread(file.path(T4, "program_communities.tsv"))
coord <- fread(file.path(T4, "coordinated_edges.tsv"))
membership <- fread(program_membership)
allz <- fread(file.path(T4, "program_interaction_edges.tsv.gz"),
              select = c("program_a", "program_b", "z_partial_standardized"))

dge <- readRDS(holdout_dge)
meta <- as.data.table(readRDS(holdout_meta))
cw <- fread(gse193066_crosswalk)
meta <- merge(meta, cw[, .(sample_id = run_id, participant_token, is_first_biopsy)],
              by = "sample_id", all.x = TRUE)
meta <- meta[is_first_biopsy == TRUE & !is.na(fibrosis_stage) & !is.na(nas_score)]
stopifnot(!anyDuplicated(meta$participant_token))
message("holdout participants: ", nrow(meta))

lib <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(lib, rownames(lib), annotation)
rm(lib, dge); invisible(gc())

features <- communities$feature_id
weights <- t4_membership_weights(membership, features)
panel <- intersect(sort(unique(weights$gene_symbol)), rownames(symbol_matrix))
PM <- symbol_matrix[panel, meta$sample_id, drop = FALSE]
rm(symbol_matrix); invisible(gc())
pairs0 <- t4_pair_table(features, weights)
key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
pairs0[, ekey := key(program_a, program_b)]
coord[, ekey := key(program_a, program_b)]
allz[, ekey := key(program_a, program_b)]
pairs0[, z_disc_coord := coord$z_partial_standardized[match(ekey, coord$ekey)]]
pairs0[, z_disc_all := allz$z_partial_standardized[match(ekey, allz$ekey)]]
pairs0[, is_coord := ekey %in% coord$ekey]
message("pairs ", nrow(pairs0), "; coordinated ", sum(pairs0$is_coord),
        "; discovery-positive coordinated ", sum(pairs0$z_disc_coord > 0, na.rm = TRUE))

ops_cache <- new.env()
s1_of <- function(idx) {
  m <- copy(meta[idx])
  m[, orig := sample_id]
  m[, sample_id := paste0(orig, "__", seq_len(.N))]
  X <- PM[, m$orig, drop = FALSE]; colnames(X) <- m$sample_id
  block <- tryCatch(t4_build_block(X, m, weights, features, program_weight_coverage),
                    error = function(e) NULL)
  if (is.null(block)) return(c(s1 = NA, base = NA, diff = NA, n_coord = NA, n_non = NA,
                               n_prog = NA, df = NA))
  gi <- as.numeric(match(block$genes, panel))
  gkey <- paste(length(gi), sum(gi), sum(gi^2), sum(gi * seq_along(gi)),
                paste(as.integer(block$coverage$cohort_testable), collapse = ""))
  ops <- ops_cache[[gkey]]
  if (is.null(ops)) {
    ops <- t4_disjoint_operators(block, pairs0, features, min_retained_weight = 0.8)
    assign(gkey, ops, envir = ops_cache)
  }
  z <- t4_one_pass(list(block), pairs0, list(ops))$z_partial_k1
  okc <- is.finite(z) & pairs0$is_coord
  okn <- is.finite(z) & !pairs0$is_coord & is.finite(pairs0$z_disc_all)
  s1 <- mean(sign(z[okc]) == sign(pairs0$z_disc_coord[okc]))
  base <- mean(sign(z[okn]) == sign(pairs0$z_disc_all[okn]))
  c(s1 = s1, base = base, diff = s1 - base, n_coord = sum(okc), n_non = sum(okn),
    n_prog = sum(block$coverage$cohort_testable), df = block$df)
}

res <- list()
full <- s1_of(seq_len(nrow(meta)))
message(sprintf("FULL n=%d S1=%.10f base=%.4f (published 0.897361912613355; edge-null mean 0.7095)",
                nrow(meta), full[["s1"]], full[["base"]]))
res[["observed_full_105"]] <- data.table(scenario = "observed_all_stage", n = nrow(meta), t(full))
f01 <- which(meta$fibrosis_stage <= 1)
obs01 <- s1_of(f01)
message(sprintf("F0-F1 n=%d S1=%.4f base=%.4f diff=%.4f", length(f01), obs01[["s1"]],
                obs01[["base"]], obs01[["diff"]]))
res[["observed_f01"]] <- data.table(scenario = "observed_F0F1", n = length(f01), t(obs01))
f2p <- which(meta$fibrosis_stage >= 2)
obs2p <- s1_of(f2p)
res[["observed_f2plus"]] <- data.table(scenario = "observed_F2plus", n = length(f2p), t(obs2p))
fwrite(rbindlist(res), file.path(out, "observed.tsv"), sep = "\t")

# Participant-level null: donors permuted independently per gene (T4 native null A).
# Destroys all gene-gene coordination; shows where S1 sits when participants carry no
# coordination at all, for comparison with the non-coordinated-edge baseline.
set.seed(SEED + 1L)
perm_s1 <- function(idx, n_draw) {
  m <- copy(meta[idx]); m[, orig := sample_id]
  X <- PM[, m$orig, drop = FALSE]
  block <- t4_build_block(X, m, weights, features, program_weight_coverage)
  ops <- t4_disjoint_operators(block, pairs0, features, min_retained_weight = 0.8)
  rbindlist(lapply(seq_len(n_draw), function(d) {
    z <- t4_one_pass(list(block), pairs0, list(ops), Zlist = list(t4_permute_genes(block$Z)))$z_partial_k1
    okc <- is.finite(z) & pairs0$is_coord
    okn <- is.finite(z) & !pairs0$is_coord & is.finite(pairs0$z_disc_all)
    data.table(draw = d, s1 = mean(sign(z[okc]) == sign(pairs0$z_disc_coord[okc])),
               base = mean(sign(z[okn]) == sign(pairs0$z_disc_all[okn])))
  }))
}
pn <- rbind(data.table(set = "all_stage_105", perm_s1(seq_len(nrow(meta)), 50L)),
            data.table(set = "F0F1_43", perm_s1(f01, 50L)))
fwrite(pn, file.path(out, "gene_permutation_null.tsv"), sep = "\t")
print(pn[, .(s1_mean = mean(s1), s1_sd = sd(s1), s1_max = max(s1), base_mean = mean(base)), by = set])
message("elapsed after observed: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")

set.seed(SEED)
boot <- list()
scen <- list(list("all_stage_m91", seq_len(nrow(meta)), 91L),
             list("all_stage_m87", seq_len(nrow(meta)), 87L),
             list("all_stage_m66", seq_len(nrow(meta)), 66L),
             list("all_stage_m53", seq_len(nrow(meta)), 53L),
             list("F0F1_m43", f01, 43L),
             list("F0F1_m66", f01, 66L),
             list("F0F1_m53", f01, 53L))
for (sc in scen) {
  for (b in seq_len(B)) {
    idx <- sample(sc[[2]], sc[[3]], replace = TRUE)
    boot[[length(boot) + 1L]] <- data.table(scenario = sc[[1]], m = sc[[3]], b = b, t(s1_of(idx)))
  }
  message(sc[[1]], " done; elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  fwrite(rbindlist(boot), file.path(out, "bootstrap_draws.tsv"), sep = "\t")
}
bd <- rbindlist(boot)
summ <- bd[, .(n_ok = sum(is.finite(diff)),
               s1_mean = mean(s1, na.rm = TRUE), s1_sd = sd(s1, na.rm = TRUE),
               base_mean = mean(base, na.rm = TRUE), base_sd = sd(base, na.rm = TRUE),
               diff_mean = mean(diff, na.rm = TRUE), diff_sd = sd(diff, na.rm = TRUE),
               diff_q025 = quantile(diff, 0.025, na.rm = TRUE),
               frac_diff_le0 = mean(diff <= 0, na.rm = TRUE),
               mde_two_sided_alpha0.025_power0.8 = sd(diff, na.rm = TRUE) * (qnorm(1 - 0.025 / 2) + qnorm(0.8)),
               mde_one_sided_alpha0.025_power0.8 = sd(diff, na.rm = TRUE) * (qnorm(1 - 0.025) + qnorm(0.8))),
           by = .(scenario, m)]
fwrite(summ, file.path(out, "bootstrap_summary.tsv"), sep = "\t")
print(rbindlist(res)); print(summ)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("total elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
