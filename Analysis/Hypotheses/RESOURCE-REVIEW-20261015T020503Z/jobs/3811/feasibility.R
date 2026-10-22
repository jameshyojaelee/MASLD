#!/usr/bin/env Rscript
# Feasibility gate for candidate L10-01 (agent gate-L10-01). Job 3811.
# Does NOT compute any covariance statistic in GSE213621. It only:
#  A. reproduces the frozen 469-donor discovery count with an analytic null
#     (to calibrate the approximation used below);
#  B. reruns approximate discovery on the 774-roster disease units of the three
#     discovery cohorts under the candidate's coarse design and under the full
#     design, and compares the two edge sets;
#  C. measures the stage/NAS structure the coarse design leaves in the residual;
#  D. runs the proposed S1 statistic as an in-data positive control by leaving
#     one discovery cohort out (discover on two, test S1 in the third);
#  E. census of GSE213621 disease units under the coarse design (n, df, programs
#     testable). No edge statistic is computed in GSE213621.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(Matrix)})
options(warn = 1)
SEED <- 20261015L
set.seed(SEED)
QUICK <- identical(Sys.getenv("QUICK"), "1")
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
jobdir <- Sys.getenv("JOBDIR", file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3811"))
out <- file.path(jobdir, if (QUICK) "out_quick" else "out")
dir.create(out, showWarnings = FALSE)
source(file.path(jobdir, "analysis_lib.R"))
source(file.path(jobdir, "t4_interaction_lib.R"))

cand <- file.path(root, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
T4 <- file.path(cand, "workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction")
sealed <- file.path(cand, "workstreams/BULK-PROGRAM-MAP-v9/sealed_inputs")
rel <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
prog_root <- file.path(root, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
annotation <- fread(file.path(cand, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
                    select = c("gene_id", "gene_name"))
communities <- fread(file.path(T4, "program_communities.tsv"))
features <- communities$feature_id
membership <- fread(file.path(prog_root, "program_membership_v2.tsv"))
weights <- t4_membership_weights(membership, features)
cat("frozen features", length(features), " weight symbols", uniqueN(weights$gene_symbol), "\n")
pairs0 <- t4_pair_table(features, weights)
key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
pairs0[, ekey := key(program_a, program_b)]
frozen_coord <- fread(file.path(T4, "coordinated_edges.tsv"))[, ekey := key(program_a, program_b)]
frozen_all <- fread(file.path(T4, "program_interaction_edges.tsv.gz"),
                    select = c("program_a", "program_b", "r_partial", "membership_entangled"))
frozen_all[, ekey := key(program_a, program_b)]

symbols_from <- function(dge) {
  lc <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  collapse_symbols(lc, rownames(lc), annotation)
}

FULL <- ~ sex_final + fibrosis_factor + nas_factor
COARSE <- ~ sex_final + fibrosis_factor
bin3 <- function(f) ifelse(f <= 1, 1L, ifelse(f == 2, 2L, 3L))

make_block <- function(sm, meta_c, formula) {
  T4_RESIDUAL_FORMULA <<- formula
  t4_build_block(sm, meta_c, weights, features, 0.80)
}

# Analytic standardization: per-pair combined Fisher z divided by
# 1/sqrt(sum over cohorts with a finite value of (df - 3)); null mean 0.
approx_discovery <- function(blocks) {
  pairs <- copy(pairs0)
  ops <- lapply(blocks, t4_disjoint_operators, pairs = pairs, features = features,
                min_retained_weight = 0.80)
  elig <- Reduce(`&`, lapply(ops, function(op) { e <- rep(TRUE, nrow(pairs)); e[op$pair_index] <- op$eligible; e }))
  obs <- t4_one_pass(blocks, pairs, ops, want_matrix = TRUE)
  w_c <- vapply(blocks, function(b) b$df - 3L, numeric(1))
  wsum <- function(zc) rowSums(matrix(w_c, nrow(zc), ncol(zc), byrow = TRUE) * is.finite(zc))
  pairs[, `:=`(entangled = !elig,
               z_disj = obs$z_disjoint, z_part = obs$z_partial_k1,
               zs_disj = obs$z_disjoint * sqrt(wsum(obs$zc_disjoint)),
               zs_part = obs$z_partial_k1 * sqrt(wsum(obs$zc_partial[[1L]])))]
  fam <- pairs[, !entangled & is.finite(zs_disj) & is.finite(zs_part)]
  pairs[, `:=`(q_disj = NA_real_, q_part = NA_real_)]
  pairs[fam, q_disj := p.adjust(2 * pnorm(-abs(zs_disj)), "BH")]
  pairs[fam, q_part := p.adjust(2 * pnorm(-abs(zs_part)), "BH")]
  pairs[, coordinated := fam & q_disj < 0.05 & q_part < 0.05 & z_part > 0]
  pairs[, in_family := fam]
  list(pairs = pairs, globals = obs$globals)
}

replicate_z <- function(block) {
  pairs <- copy(pairs0)
  ops <- list(t4_disjoint_operators(block, pairs, features, min_retained_weight = 0.80))
  obs <- t4_one_pass(list(block), pairs, ops)
  pairs[, z_rep := obs$z_partial_k1]
  pairs[, .(ekey, z_rep)]
}

s1_test <- function(disc, rep, n_draws = 1000L) {
  m <- merge(disc[in_family == TRUE, .(ekey, coordinated, z_part)], rep, by = "ekey")
  m <- m[is.finite(z_rep) & is.finite(z_part)]
  co <- m[coordinated == TRUE]; nc <- m[coordinated == FALSE]
  s1 <- mean(sign(co$z_rep) == sign(co$z_part))
  if (nrow(nc) < nrow(co)) return(data.table(n_coord = nrow(co), n_noncoord = nrow(nc), s1 = s1,
                                             null_mean = NA, null_sd = NA, p = NA, null_feasible = FALSE))
  null <- replicate(n_draws, { i <- sample.int(nrow(nc), nrow(co)); mean(sign(nc$z_rep[i]) == sign(nc$z_part[i])) })
  data.table(n_coord = nrow(co), n_noncoord = nrow(nc), draw_overlap = nrow(co) / nrow(nc), s1 = s1,
             null_mean = mean(null), null_sd = sd(null),
             noncoord_all_concordance = mean(sign(nc$z_rep) == sign(nc$z_part)),
             p = (1 + sum(null >= s1)) / (1 + n_draws), null_feasible = TRUE)
}

res <- list()

# ------------------------------------------------------------ A. calibration
if (!QUICK) {
  cat("[A] calibration on sealed 469-donor discovery input, full design\n")
  dgeA <- readRDS(file.path(sealed, "nonholdout_dge.rds"))
  metaA <- as.data.table(readRDS(file.path(sealed, "nonholdout_meta.rds")))
  metaA <- metaA[dataset %in% c("GSE130970", "GSE135251", "GSE162694", "GSE174478") &
                   !is.na(fibrosis_stage) & !is.na(nas_score)]
  smA <- symbols_from(dgeA); rm(dgeA); gc()
  panel <- intersect(sort(unique(weights$gene_symbol)), rownames(smA))
  smA <- smA[panel, metaA$sample_id, drop = FALSE]
  blocksA <- lapply(split(metaA, by = "dataset"), function(mc) make_block(smA, mc, FULL))
  print(rbindlist(lapply(blocksA, function(b) data.table(cohort = b$cohort, n = b$n, n_par = b$n_par, df = b$df,
                                                         testable = sum(b$coverage$cohort_testable)))))
  dA <- approx_discovery(blocksA)
  chk <- merge(dA$pairs[, .(ekey, z_part, coordinated)], frozen_all[, .(ekey, r_partial)], by = "ekey")
  res$A <- data.table(n_units = nrow(metaA), family = dA$pairs[in_family == TRUE, .N],
                      coordinated_approx = dA$pairs[coordinated == TRUE, .N], coordinated_frozen = nrow(frozen_coord),
                      overlap_with_frozen = dA$pairs[coordinated == TRUE & ekey %in% frozen_coord$ekey, .N],
                      cor_r_partial_vs_frozen = cor(tanh(chk$z_part), chk$r_partial, use = "complete.obs"),
                      positive_fraction_partial = dA$pairs[in_family == TRUE, mean(z_part > 0)])
  print(res$A)
  rm(smA); gc()
}

# ------------------------------------------------------- 774 roster inputs
dge <- readRDS(file.path(rel, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
sm <- symbols_from(dge); rm(dge); gc()
panel <- intersect(sort(unique(weights$gene_symbol)), rownames(sm))
cat("panel symbols in 774 universe:", length(panel), "of", uniqueN(weights$gene_symbol), "\n")
sm <- sm[panel, , drop = FALSE]
disc_cohorts <- if (QUICK) c("GSE130970", "GSE162694") else c("GSE130970", "GSE135251", "GSE162694")
md <- meta[dataset %in% disc_cohorts & group_binary == "Disease"]
md[, fib5 := fibrosis_stage]
cat("discovery disease units:", nrow(md), "\n"); print(md[, .N, by = dataset])

blocks_for <- function(m, formula, coarse) {
  m <- copy(m)
  if (coarse) m[, fibrosis_stage := bin3(fib5)] else m[, fibrosis_stage := fib5]
  lapply(split(m, by = "dataset"), function(mc) make_block(sm, mc, formula))
}
bsum <- function(bl) rbindlist(lapply(bl, function(b) data.table(cohort = b$cohort, n = b$n, n_par = b$n_par, df = b$df,
                                                                 panel_genes = length(b$genes),
                                                                 testable = sum(b$coverage$cohort_testable))))

# ------------------------------------------------- E. GSE213621 census only
g <- copy(meta[dataset == "GSE213621" & group_binary == "Disease"])
g[, fib5 := fibrosis_stage]   # already the 3 bins 1/2/3 = F0-F1/F2/F3-F4
bg <- make_block(sm, g, COARSE)
res$E <- data.table(cohort = "GSE213621", n = bg$n, n_par = bg$n_par, df = bg$df, rank = qr(bg$X)$rank,
                    panel_genes = length(bg$genes), testable = sum(bg$coverage$cohort_testable),
                    min_cov = min(bg$coverage$cohort_weight_coverage))
print(res$E); print(table(bin = g$fibrosis_stage, sex = g$sex_final))
rm(bg); gc()

# ------------------------------------------------------------ B. designs
cat("[B] coarse (all 371) vs full and coarse on NAS-complete units\n")
mdn <- md[!is.na(nas_score)]
b_coarse_all <- blocks_for(md, COARSE, TRUE)
b_coarse_nas <- blocks_for(mdn, COARSE, TRUE)
b_full_nas <- blocks_for(mdn, FULL, FALSE)
print(bsum(b_coarse_all)); print(bsum(b_full_nas))
d_ca <- approx_discovery(b_coarse_all)
d_cn <- approx_discovery(b_coarse_nas)
d_fn <- approx_discovery(b_full_nas)
cmp <- merge(d_cn$pairs[, .(ekey, zc = z_part, cc = coordinated, fam = in_family)],
             d_fn$pairs[, .(ekey, zf = z_part, cf = coordinated)], by = "ekey")
cmp <- cmp[fam == TRUE]
res$B <- data.table(
  design = c("coarse_all371", "coarse_nas_complete", "full_nas_complete"),
  n_units = c(nrow(md), nrow(mdn), nrow(mdn)),
  family = c(d_ca$pairs[in_family == TRUE, .N], d_cn$pairs[in_family == TRUE, .N], d_fn$pairs[in_family == TRUE, .N]),
  coordinated = c(d_ca$pairs[coordinated == TRUE, .N], d_cn$pairs[coordinated == TRUE, .N], d_fn$pairs[coordinated == TRUE, .N]),
  positive_fraction_partial = c(d_ca$pairs[in_family == TRUE, mean(z_part > 0)], d_cn$pairs[in_family == TRUE, mean(z_part > 0)],
                                d_fn$pairs[in_family == TRUE, mean(z_part > 0)]),
  overlap_with_frozen_2426 = c(d_ca$pairs[coordinated == TRUE & ekey %in% frozen_coord$ekey, .N],
                               d_cn$pairs[coordinated == TRUE & ekey %in% frozen_coord$ekey, .N],
                               d_fn$pairs[coordinated == TRUE & ekey %in% frozen_coord$ekey, .N]))
print(res$B)
res$B2 <- data.table(
  spearman_partial_coarse_vs_full = cor(cmp$zc, cmp$zf, method = "spearman"),
  mean_abs_diff_r = mean(abs(tanh(cmp$zc) - tanh(cmp$zf))),
  coarse_only = cmp[cc & !cf, .N], full_only = cmp[!cc & cf, .N], both = cmp[cc & cf, .N],
  sign_agree_among_coarse_coord = cmp[cc == TRUE, mean(sign(zc) == sign(zf))],
  mean_r_coarse_coord = cmp[cc == TRUE, mean(tanh(zc))], mean_r_same_edges_full = cmp[cc == TRUE, mean(tanh(zf))])
print(res$B2)
# Community re-derivation on the coarse rediscovery (seeded Louvain), mapped to the frozen partition
set.seed(SEED)
w_new <- ifelse(d_ca$pairs$coordinated, tanh(d_ca$pairs$z_part), NA_real_)
cm <- t4_modularity(d_ca$pairs, w_new, features)
xt <- merge(cm$membership, communities[, .(feature_id, frozen = community)], by = "feature_id")
cat("coarse rediscovery modularity", cm$modularity, " communities", cm$n_communities, "\n")
print(table(new = xt$community, frozen = xt$frozen))
res$B3 <- data.table(modularity = cm$modularity, n_communities = cm$n_communities,
                     frozen_c3_size = sum(xt$frozen == 3),
                     best_new_for_frozen_c3 = as.integer(names(which.max(table(xt[frozen == 3]$community)))),
                     frozen_c3_in_best = max(table(xt[frozen == 3]$community)))
print(res$B3)
fwrite(d_ca$pairs[, .(ekey, program_a, program_b, entangled, in_family, z_disj, z_part, zs_disj, zs_part, q_disj, q_part, coordinated)],
       file.path(out, "approx_rediscovery_coarse_371.tsv.gz"), sep = "\t")

# ------------------------------------------------------------ C. leakage
cat("[C] residual stage and NAS structure left by the coarse design\n")
res$C <- rbindlist(lapply(names(b_coarse_nas), function(cc) {
  b <- b_coarse_nas[[cc]]
  Yt <- t(b$Z); B <- b$xtx_inv %*% crossprod(b$X, Yt)
  U <- Yt - b$X %*% B                           # unscaled residual gene matrix
  P <- U %*% b$W
  P <- P[, colSums(abs(P)) > 0, drop = FALSE]
  mc <- mdn[match(b$samples, sample_id)]
  f5 <- droplevels(factor(mc$fib5)); nb <- droplevels(t4_nas_bin(mc$nas_score))
  Xf <- model.matrix(~ factor(bin3(mc$fib5)) + factor(mc$sex_final) + f5 + nb)
  Xf <- Xf[, qr(Xf)$pivot[seq_len(qr(Xf)$rank)], drop = FALSE]
  Xc <- b$X
  pv <- apply(P, 2L, function(y) { a <- lm.fit(Xc, y); f <- lm.fit(Xf, y)
    rss0 <- sum(a$residuals^2); rss1 <- sum(f$residuals^2); d1 <- ncol(Xf) - ncol(Xc); d2 <- length(y) - ncol(Xf)
    c(r2 = 1 - rss1 / rss0, p = pf(((rss0 - rss1) / d1) / (rss1 / d2), d1, d2, lower.tail = FALSE)) })
  q <- p.adjust(pv["p", ], "BH")
  data.table(cohort = cc, n = nrow(P), extra_df = ncol(Xf) - ncol(Xc), programs = ncol(P),
             chance_r2 = (ncol(Xf) - ncol(Xc)) / (nrow(P) - ncol(Xc)),
             median_r2 = median(pv["r2", ]), max_r2 = max(pv["r2", ]),
             n_q05 = sum(q < 0.05), n_p05 = sum(pv["p", ] < 0.05))
}))
print(res$C)

# ------------------------------------------------- D. in-data positive control
cat("[D] leave-one-cohort-out S1 under the coarse design\n")
res$D <- rbindlist(lapply(disc_cohorts, function(h) {
  d <- approx_discovery(b_coarse_all[setdiff(disc_cohorts, h)])
  r <- replicate_z(b_coarse_all[[h]])
  cbind(data.table(held_out = h, n_test = b_coarse_all[[h]]$n), s1_test(d$pairs, r))
}))
print(res$D)
# Frozen 2,426 edges tested in each discovery cohort under the coarse design (secondary-arm analogue)
frz <- copy(pairs0)[, `:=`(in_family = ekey %in% frozen_all[membership_entangled == FALSE | is.na(membership_entangled), ekey],
                          coordinated = ekey %in% frozen_coord$ekey)]
frz <- merge(frz, frozen_all[, .(ekey, z_part = atanh(r_partial))], by = "ekey")
res$D2 <- rbindlist(lapply(disc_cohorts, function(h) cbind(data.table(test_cohort = h, edges = "frozen_2426_not_independent"),
                                                            s1_test(frz, replicate_z(b_coarse_all[[h]])))))
print(res$D2)

for (nm in names(res)) fwrite(res[[nm]], file.path(out, paste0("result_", nm, ".tsv")), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("DONE\n")
