# H40: alcohol-associated expression signature in MASLD-labeled bulk biopsies.
#
# Usage: Rscript --vanilla h40_alcohol_signature.R <H40_dir> <run_dir> <n_cores>
#   <run_dir> must already contain parse_geo.py outputs (source_design.tsv,
#   counts_GSE142530.tsv.gz, counts_GSE143318.tsv.gz).
#
# Design is fixed in <H40_dir>/PRESPEC_H40.md. Native method: limma-voom on
# bulk RNA-seq counts (source cohorts and the refit cohorts).
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
H40 <- args[1]
RUN <- args[2]
N_CORES <- as.integer(args[3])
stopifnot(dir.exists(H40), dir.exists(RUN), file.exists(file.path(H40, "PRESPEC_H40.md")))

SEED <- 20260923L
B_NULL <- 1000L
B_BOOT <- 2000L
B_EXCL <- 200L
TAIL_RZ <- 2.5
SIG_LFC <- 1
SIG_FDR <- 0.05
CALL_LFC <- 0.5
CALL_FDR <- 0.05
N_CALLS <- 1854L
COHORTS_SRC <- c("GSE142530", "GSE143318")

REFIT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
GENE_META <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"

stop_with <- function(code, ...) {
  msg <- paste0("STOP ", code, ": ", ...)
  writeLines(msg, file.path(RUN, "STOPPED.txt"))
  stop(msg, call. = FALSE)
}
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
out_path <- function(f) file.path(RUN, f)
wtsv <- function(x, f) fwrite(x, out_path(f), sep = "\t")
strip_v <- function(x) sub("[.][0-9]+$", "", x)
summary_list <- list()

# ------------------------------------------------------------------ gene sets
gm <- fread(GENE_META)
gm[, eb := strip_v(gene_id)]
name_counts <- gm[, .N, by = gene_name]
unique_names <- name_counts[N == 1L, gene_name]
sym2ens <- function(sym) gm[gene_name %in% unique_names][match(sym, gene_name), eb]

reg <- fread(file.path(H40, "gene_sets/collectri_HNF4A_regulon.tsv"))
reg[, eb := sym2ens(target)]
hnf_act <- unique(na.omit(reg[mor > 0, eb]))
hnf_all <- unique(na.omit(reg$eb))
col_set <- fread(file.path(H40, "gene_sets/msigdb_reactome_collagen_formation.tsv"))
col_ids <- unique(col_set$ensembl_gene)
sexchr_ids <- unique(c(gm[chromosome == "chrY", eb], gm[gene_name %in% c("XIST", "TSIX"), eb]))
summary_list$gene_sets <- list(
  hnf4a_regulon_targets = nrow(reg), hnf4a_regulon_mapped = sum(!is.na(reg$eb)),
  hnf4a_activating_mapped = length(hnf_act), collagen_ids = length(col_ids),
  chrY_xist_tsix_ids = length(sexchr_ids))

# ------------------------------------------------------------------ refit DGE
dge_refit <- readRDS(file.path(REFIT, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(REFIT, "meta_matched.rds")))
meta <- meta[match(colnames(dge_refit), sample_id)]
if (!identical(meta$sample_id, colnames(dge_refit))) stop_with("S2", "meta/DGE sample mismatch")
if (anyDuplicated(meta$analysis_unit_id)) stop_with("S2", "duplicate analysis units")
refit_ids_all <- rownames(dge_refit)
refit_keep_rows <- !grepl("_PAR_Y$", refit_ids_all)
refit_eb <- strip_v(refit_ids_all[refit_keep_rows])
if (anyDuplicated(refit_eb)) stop_with("S2", "refit version-stripped IDs not unique")
deg_ref <- fread(file.path(REFIT, "deg_results.csv"))
ref_calls <- deg_ref[padj < CALL_FDR & abs(logFC) > CALL_LFC]
if (nrow(ref_calls) != N_CALLS) stop_with("S2", "deg_results.csv has ", nrow(ref_calls), " calls")
sw_full <- fread(file.path(REFIT, "lvqw_sample_weights.tsv"))
sw_full <- sw_full[match(colnames(dge_refit), sample_id), sample_weight]
if (anyNA(sw_full)) stop_with("S2", "missing stored sample weights")

# ------------------------------------------------------------------ source cohorts
design_src <- fread(out_path("source_design.tsv"))
load_counts <- function(coh) {
  x <- fread(out_path(sprintf("counts_%s.tsv.gz", coh)))
  m <- as.matrix(x[, -1, with = FALSE])
  rownames(m) <- x$gene_id
  d <- design_src[cohort == coh][match(colnames(m), column_id)]
  if (anyNA(d$gsm)) stop_with("S1", coh, " column without design row")
  list(counts = m, design = d)
}
src <- setNames(lapply(COHORTS_SRC, load_counts), COHORTS_SRC)

U <- Reduce(intersect, list(rownames(src$GSE142530$counts), rownames(src$GSE143318$counts), refit_eb))
U <- setdiff(U, c(hnf_all, col_ids, sexchr_ids))
summary_list$universe <- list(
  GSE142530_genes = nrow(src$GSE142530$counts), GSE143318_genes = nrow(src$GSE143318$counts),
  refit_genes = length(refit_eb), U = length(U))
log_msg("Universe U = %d genes", length(U))

zscore_rows <- function(L) {
  mu <- rowMeans(L)
  s <- sqrt(rowSums((L - mu)^2) / (ncol(L) - 1))
  z <- (L - mu) / s
  z[!is.finite(s) | s == 0, ] <- NA_real_
  z
}
set_mean <- function(Z, ids) {
  ids <- intersect(ids, rownames(Z))
  if (!length(ids)) return(rep(NA_real_, ncol(Z)))
  colMeans(Z[ids, , drop = FALSE], na.rm = TRUE)
}
n_measured <- function(Z, ids) {
  ids <- intersect(ids, rownames(Z))
  sum(rowSums(is.finite(Z[ids, , drop = FALSE])) > 0)
}
sig_score <- function(Z, up, down) {
  nu <- n_measured(Z, up); nd <- n_measured(Z, down)
  if (nu >= 5 && nd >= 5) return(set_mean(Z, up) - set_mean(Z, down))
  if (nu >= 5) return(set_mean(Z, up))
  if (nd >= 5) return(-set_mean(Z, down))
  rep(NA_real_, ncol(Z))
}

prep_source <- function(coh) {
  s <- src[[coh]]
  dge <- DGEList(counts = s$counts)
  alc <- factor(s$design$alcohol, levels = c(0, 1))
  keep <- filterByExpr(dge, group = alc)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "RLE")
  Z <- zscore_rows(cpm(dge, log = TRUE, prior.count = 2))
  cov <- data.table(column_id = colnames(dge), alcohol = s$design$alcohol, group = s$design$group,
                    jh_explant = s$design$jh_explant,
                    hnf4a = set_mean(Z, hnf_act), collagen = set_mean(Z, col_ids))
  log_msg("%s: %d samples, %d filtered genes, %d in U; HNF4A act measured %d, collagen measured %d",
          coh, ncol(dge), nrow(dge), sum(rownames(dge) %in% U),
          n_measured(Z, hnf_act), n_measured(Z, col_ids))
  list(dge = dge, Z = Z, cov = cov, coh = coh,
       n = list(samples = ncol(dge), filtered = nrow(dge), in_U = sum(rownames(dge) %in% U),
                hnf4a_measured = n_measured(Z, hnf_act), collagen_measured = n_measured(Z, col_ids)))
}
P <- setNames(lapply(COHORTS_SRC, prep_source), COHORTS_SRC)
summary_list$source_prep <- lapply(P, `[[`, "n")
with_cohort <- function(p) copy(p$cov)[, cohort := p$coh][]
wtsv(rbindlist(lapply(P, with_cohort)), "source_covariates.tsv")

fit_alcohol <- function(dge, X) {
  v <- voomWithQualityWeights(dge, X)
  fit <- eBayes(lmFit(v, X))
  tt <- topTable(fit, coef = "alcohol", number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             P.Value = tt$P.Value, padj = tt$adj.P.Val)
}
designs <- list(
  S1 = function(cov) model.matrix(~ alcohol + hnf4a + collagen, data = cov),
  S2 = function(cov) model.matrix(~ alcohol, data = cov))

percohort <- list(S1 = list(), S2 = list())
for (sg in names(designs)) for (coh in COHORTS_SRC) {
  p <- P[[coh]]
  d <- p$dge[rownames(p$dge) %in% U, ]
  X <- designs[[sg]](p$cov)
  log_msg("Per-cohort fit %s %s: %d genes", sg, coh, nrow(d))
  percohort[[sg]][[coh]] <- fit_alcohol(d, X)
}

pooled <- list()
common <- intersect(intersect(rownames(P$GSE142530$dge), rownames(P$GSE143318$dge)), U)
dge_pool <- cbind(P$GSE142530$dge[common, ], P$GSE143318$dge[common, ])
cov_pool <- rbind(with_cohort(P$GSE142530), with_cohort(P$GSE143318))
stopifnot(identical(colnames(dge_pool), cov_pool$column_id))
Xp <- list(S1 = model.matrix(~ cohort + alcohol + hnf4a + collagen, data = cov_pool),
           S2 = model.matrix(~ cohort + alcohol, data = cov_pool))
signatures <- list()
for (sg in names(Xp)) {
  log_msg("Pooled fit %s: %d genes", sg, nrow(dge_pool))
  tt <- fit_alcohol(dge_pool, Xp[[sg]])
  a <- percohort[[sg]]$GSE142530[, .(gene, logFC_GSE142530 = logFC, padj_GSE142530 = padj)]
  b <- percohort[[sg]]$GSE143318[, .(gene, logFC_GSE143318 = logFC, padj_GSE143318 = padj)]
  tt <- merge(merge(tt, a, by = "gene", all.x = TRUE), b, by = "gene", all.x = TRUE)
  tt[, symbol := gm[match(gene, eb), gene_name]]
  tt[, concordant := sign(logFC) == sign(logFC_GSE142530) & sign(logFC) == sign(logFC_GSE143318)]
  tt[, in_signature := padj < SIG_FDR & abs(logFC) > SIG_LFC & concordant %in% TRUE]
  tt[, direction := fifelse(in_signature, fifelse(logFC > 0, "up", "down"), NA_character_)]
  pooled[[sg]] <- tt
  wtsv(tt[order(P.Value)], sprintf("pooled_fit_%s.tsv", sg))
  signatures[[sg]] <- list(up = tt[direction == "up", gene], down = tt[direction == "down", gene])
  log_msg("Signature %s: %d up, %d down", sg, length(signatures[[sg]]$up), length(signatures[[sg]]$down))
}
for (sg in names(percohort)) for (coh in COHORTS_SRC) {
  wtsv(percohort[[sg]][[coh]][order(P.Value)], sprintf("percohort_fit_%s_%s.tsv", sg, coh))
}

# ------------------------------------------------------------------ positive control
auc_mw <- function(x, y) {
  r <- rank(c(x, y))
  (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y))
}
boot_auc <- function(x, y) {
  set.seed(SEED)
  b <- replicate(B_BOOT, auc_mw(sample(x, replace = TRUE), sample(y, replace = TRUE)))
  unname(quantile(b, c(0.025, 0.975)))
}
pc_rows <- list()
for (sg in names(designs)) for (train in COHORTS_SRC) {
  test <- setdiff(COHORTS_SRC, train)
  pt <- percohort[[sg]][[train]]
  up <- pt[padj < SIG_FDR & logFC > SIG_LFC, gene]
  dn <- pt[padj < SIG_FDR & logFC < -SIG_LFC, gene]
  H <- P[[test]]
  sc <- sig_score(H$Z, up, dn)
  cv <- copy(H$cov)
  cv[, score := sc]
  cv[, resid := NA_real_]
  if (!anyNA(sc)) cv[, resid := residuals(lm(score ~ hnf4a + collagen, data = cv))]
  evals <- list(all = cv, no_jh = cv[jh_explant == 0])
  for (ev in names(evals)) {
    e <- evals[[ev]]
    if (ev == "no_jh" && test != "GSE142530") next
    for (sc_type in c("resid", "score")) {
      for (case in c("alcohol", "AH", "AC")) {
        xs <- if (case == "alcohol") e[alcohol == 1, get(sc_type)] else e[group == case, get(sc_type)]
        ys <- e[alcohol == 0, get(sc_type)]
        ok <- length(up) + length(dn) > 0 && !anyNA(c(xs, ys))
        auc <- if (ok) auc_mw(xs, ys) else NA_real_
        ci <- if (ok) boot_auc(xs, ys) else c(NA_real_, NA_real_)
        pc_rows[[length(pc_rows) + 1L]] <- data.table(
          signature = sg, fold = paste0(train, "->", test), evaluation = ev, score_type = sc_type,
          cases = case, n_case = length(xs), n_normal = length(ys),
          fold_sig_up = length(up), fold_sig_down = length(dn),
          up_measured = n_measured(H$Z, up), down_measured = n_measured(H$Z, dn),
          auc = auc, ci_low = ci[1], ci_high = ci[2])
      }
    }
  }
  wtsv(cv, sprintf("pc_scores_%s_%s_to_%s.tsv", sg, train, test))
}
pc <- rbindlist(pc_rows)
wtsv(pc, "positive_control.tsv")
pc_main <- pc[evaluation == "all" & score_type == "resid" & cases == "alcohol"]
pc_pass <- sapply(names(designs), function(sg) {
  r <- pc_main[signature == sg]
  isTRUE(nrow(r) == 2L && !anyNA(r$auc) && !anyNA(r$ci_low) &&
         all(r$fold_sig_up + r$fold_sig_down >= 10) && all(r$auc >= 0.75 & r$ci_low > 0.5))
})
size_ok <- sapply(names(signatures), function(sg) length(unlist(signatures[[sg]])) >= 20L)
summary_list$positive_control <- list(pass = as.list(pc_pass), size_ok = as.list(size_ok),
                                      main = pc_main)
log_msg("Positive control pass: S1=%s S2=%s; size ok: S1=%s S2=%s",
        pc_pass["S1"], pc_pass["S2"], size_ok["S1"], size_ok["S2"])

usable <- names(designs)[pc_pass & size_ok]
primary_sig <- if ("S1" %in% usable) "S1" else if ("S2" %in% usable) "S2" else NA_character_
summary_list$carried <- list(usable = usable, primary_signature = primary_sig,
                             primary_label = if (identical(primary_sig, "S1")) "primary" else if (identical(primary_sig, "S2")) "secondary (S1 failed)" else "none")
if (is.na(primary_sig)) {
  writeLines("No signature passed S3 and S4; refit cohorts not scored.", out_path("STOPPED.txt"))
  write_json(summary_list, out_path("summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  writeLines(capture.output(sessionInfo()), out_path("sessionInfo.txt"))
  log_msg("Stopped after positive control.")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------ refit cohort scoring
L_refit <- cpm(dge_refit, log = TRUE, prior.count = 2)[refit_keep_rows, ]
rownames(L_refit) <- refit_eb
cohort_idx <- split(seq_len(ncol(L_refit)), meta$dataset)
Zr <- matrix(NA_real_, nrow(L_refit), ncol(L_refit), dimnames = dimnames(L_refit))
for (ix in cohort_idx) Zr[, ix] <- zscore_rows(L_refit[, ix, drop = FALSE])
meta[, stratum := fifelse(group_binary == "Control", "Control",
                          fifelse(is.na(fibrosis_stage), "unstaged", paste0("F", fibrosis_stage)))]
meta[, hnf4a := set_mean(Zr, hnf_act)]
meta[, collagen := set_mean(Zr, col_ids)]
meta[, cell := paste(dataset, stratum, sep = "|")]

residualize <- function(score) {
  r <- rep(NA_real_, length(score))
  for (ix in cohort_idx) {
    df <- data.frame(score = score[ix], stratum = factor(meta$stratum[ix]),
                     hnf4a = meta$hnf4a[ix], collagen = meta$collagen[ix],
                     sex = factor(meta$inferred_sex[ix]))
    rhs <- c(if (nlevels(df$stratum) > 1) "stratum", "hnf4a", "collagen", if (nlevels(df$sex) > 1) "sex")
    r[ix] <- residuals(lm(reformulate(rhs, "score"), data = df))
  }
  r
}
robust_z <- function(r) {
  z <- rep(NA_real_, length(r))
  for (ix in cohort_idx) z[ix] <- (r[ix] - median(r[ix])) / mad(r[ix])
  z
}
is_dis <- meta$group_binary == "Disease"
tail_count <- function(rz, thr = TAIL_RZ) sum(rz > thr & is_dis)
tail_by_cohort <- function(rz, thr = TAIL_RZ) {
  vapply(names(cohort_idx), function(c) sum(rz[cohort_idx[[c]]] > thr & is_dis[cohort_idx[[c]]]), numeric(1))
}

measured_in_refit <- function(ids) intersect(ids, rownames(Zr))
mean_logcpm <- rowMeans(L_refit)
U_refit <- intersect(U, rownames(Zr))
bins_all <- cut(rank(mean_logcpm[U_refit], ties.method = "first"),
                breaks = quantile(seq_along(U_refit), probs = seq(0, 1, length.out = 21)),
                include.lowest = TRUE, labels = FALSE)
names(bins_all) <- U_refit

score_one_signature <- function(sg) {
  up <- measured_in_refit(signatures[[sg]]$up)
  dn <- measured_in_refit(signatures[[sg]]$down)
  frac <- (length(up) + length(dn)) / length(unlist(signatures[[sg]]))
  log_msg("%s: %d/%d up and %d/%d down measured in refit (%.1f%%)", sg, length(up),
          length(signatures[[sg]]$up), length(dn), length(signatures[[sg]]$down), 100 * frac)
  sc <- sig_score(Zr, up, dn)
  r <- residualize(sc)
  rz <- robust_z(r)
  tab <- meta[, .(sample_id, analysis_unit_id, dataset, group_binary, stratum, fibrosis_stage,
                  inferred_sex, hnf4a, collagen)]
  tab[, `:=`(score = sc, resid = r, rz = rz, tail = rz > TAIL_RZ & group_binary == "Disease")]
  wtsv(tab, sprintf("refit_scores_%s.tsv", sg))

  # per-cohort description
  desc <- rbindlist(lapply(c(names(cohort_idx), "ALL"), function(c) {
    ix <- if (c == "ALL") seq_len(nrow(tab)) else cohort_idx[[c]]
    rbindlist(lapply(c("Disease", "Control"), function(g) {
      t <- tab[ix][group_binary == g]
      n <- nrow(t)
      rbindlist(lapply(c(2, 2.5, 3), function(thr) {
        k <- sum(t$rz > thr)
        ci <- if (n > 0) binom.test(k, n)$conf.int else c(NA, NA)
        data.table(dataset = c, group = g, threshold = thr, n = n, n_tail = k,
                   frac = if (n > 0) k / n else NA_real_, ci_low = ci[1], ci_high = ci[2],
                   score_median = median(t$score), score_q95 = if (n > 0) unname(quantile(t$score, 0.95)) else NA_real_)
      }))
    }))
  }))
  wtsv(desc, sprintf("tail_summary_%s.tsv", sg))
  wtsv(tab[tail == TRUE, .N, by = .(dataset, stratum, inferred_sex)][order(dataset, stratum)],
       sprintf("tail_composition_%s.tsv", sg))

  T_obs <- tail_count(rz)
  T_obs_c <- tail_by_cohort(rz)

  # Null 1: within cohort x stratum independent permutation of each signature gene
  Zs <- Zr[c(up, dn), , drop = FALSE]
  cells <- split(seq_len(ncol(Zs)), meta$cell)
  set.seed(SEED)
  null1 <- matrix(NA_real_, B_NULL, length(cohort_idx) + 1L,
                  dimnames = list(NULL, c("ALL", names(cohort_idx))))
  for (b in seq_len(B_NULL)) {
    Zp <- Zs
    for (ix in cells) {
      if (length(ix) < 2) next
      Zp[, ix] <- t(apply(Zs[, ix, drop = FALSE], 1, function(v) v[sample.int(length(v))]))
    }
    rzb <- robust_z(residualize(sig_score(Zp, up, dn)))
    null1[b, ] <- c(tail_count(rzb), tail_by_cohort(rzb))
  }

  # Null 2: matched random gene sets
  pool <- setdiff(U_refit, c(up, dn))
  need <- table(factor(bins_all[c(up, dn)], levels = 1:20))
  pool_by_bin <- split(pool, factor(bins_all[pool], levels = 1:20))
  short_bins <- names(need)[need > lengths(pool_by_bin)[names(need)]]
  set.seed(SEED)
  null2 <- matrix(NA_real_, B_NULL, length(cohort_idx) + 1L,
                  dimnames = list(NULL, c("ALL", names(cohort_idx))))
  up_bins <- bins_all[up]; dn_bins <- bins_all[dn]
  for (b in seq_len(B_NULL)) {
    ru <- character(0); rd <- character(0)
    for (k in 1:20) {
      nu <- sum(up_bins == k); nd <- sum(dn_bins == k)
      if (nu + nd == 0) next
      pk <- pool_by_bin[[as.character(k)]]
      draw <- pk[sample.int(length(pk), nu + nd, replace = (nu + nd) > length(pk))]
      ru <- c(ru, draw[seq_len(nu)]); rd <- c(rd, draw[nu + seq_len(nd)])
    }
    rzb <- robust_z(residualize(sig_score(Zr, ru, rd)))
    null2[b, ] <- c(tail_count(rzb), tail_by_cohort(rzb))
  }
  wtsv(as.data.table(null1), sprintf("null1_tail_counts_%s.tsv", sg))
  wtsv(as.data.table(null2), sprintf("null2_tail_counts_%s.tsv", sg))

  obs_all <- c(ALL = T_obs, T_obs_c)
  null_sum <- rbindlist(lapply(names(obs_all), function(c) data.table(
    dataset = c, observed = obs_all[[c]],
    null1_median = median(null1[, c]), null1_q025 = quantile(null1[, c], 0.025), null1_q975 = quantile(null1[, c], 0.975),
    p1 = (1 + sum(null1[, c] >= obs_all[[c]])) / (B_NULL + 1),
    null2_median = median(null2[, c]), null2_q025 = quantile(null2[, c], 0.025), null2_q975 = quantile(null2[, c], 0.975),
    p2 = (1 + sum(null2[, c] >= obs_all[[c]])) / (B_NULL + 1))))
  wtsv(null_sum, sprintf("tail_null_summary_%s.tsv", sg))
  p1 <- null_sum[dataset == "ALL", p1]; p2 <- null_sum[dataset == "ALL", p2]
  log_msg("%s: T_obs=%d; null1 median %.1f p1=%.4f; null2 median %.1f p2=%.4f",
          sg, T_obs, null_sum[dataset == "ALL", null1_median], p1, null_sum[dataset == "ALL", null2_median], p2)
  list(tab = tab, T_obs = T_obs, p1 = p1, p2 = p2, p_tail = max(p1, p2),
       n_up = length(up), n_down = length(dn), frac_measured = frac, short_bins = short_bins)
}

scored <- list()
for (sg in usable) scored[[sg]] <- score_one_signature(sg)
summary_list$tail <- lapply(scored, function(s) list(
  T_obs = s$T_obs, p1 = s$p1, p2 = s$p2, p_tail = s$p_tail, separable = s$p1 < 0.05 && s$p2 < 0.05,
  n_up_measured = s$n_up, n_down_measured = s$n_down, frac_measured = s$frac_measured,
  flag_low_measured = s$frac_measured < 0.8, null2_bins_with_replacement = s$short_bins))
write_json(summary_list, out_path("summary_partial.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

# ------------------------------------------------------------------ refit
mm_refit <- function(keep) {
  info <- data.frame(
    group_binary = factor(dge_refit$samples$group_binary[keep], levels = c("Control", "Disease")),
    dataset = factor(dge_refit$samples$dataset[keep]),
    inferred_sex = factor(meta$inferred_sex[keep]))
  X <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)
  if (qr(X)$rank != ncol(X)) stop_with("S2", "rank-deficient refit design")
  X
}
calls_of <- function(tt) {
  data.table(gene = rownames(tt), logFC = tt$logFC, padj = tt$adj.P.Val)[, call := padj < CALL_FDR & abs(logFC) > CALL_LFC]
}
fit_exact <- function(keep, recompute_rle) {
  d <- dge_refit[, keep]
  if (recompute_rle) d <- calcNormFactors(d, method = "RLE")
  X <- mm_refit(keep)
  v <- voomWithQualityWeights(d, X)
  calls_of(topTable(eBayes(lmFit(v, X)), coef = "group_binaryDisease", number = Inf, sort.by = "none"))
}
fit_fixed <- function(keep, recompute_rle) {
  d <- dge_refit[, keep]
  if (recompute_rle) d <- calcNormFactors(d, method = "RLE")
  X <- mm_refit(keep)
  w <- sw_full[keep]
  v <- voom(d, X, weights = w)
  v$weights <- t(w * t(v$weights))
  calls_of(topTable(eBayes(lmFit(v, X)), coef = "group_binaryDisease", number = Inf, sort.by = "none"))
}
compare_calls <- function(full, ex) {
  m <- merge(full[, .(gene, lf0 = logFC, c0 = call)], ex[, .(gene, lf1 = logFC, c1 = call)], by = "gene")
  stopifnot(nrow(m) == nrow(full))
  lost <- sum(m$c0 & !m$c1); gained <- sum(!m$c0 & m$c1)
  flipped <- sum(m$c0 & m$c1 & sign(m$lf0) != sign(m$lf1))
  list(n_calls = sum(m$c1), lost = lost, gained = gained, flipped = flipped,
       changed = lost + gained + flipped, retained = sum(m$c0 & m$c1 & sign(m$lf0) == sign(m$lf1)),
       spearman_logFC = cor(m$lf0, m$lf1, method = "spearman"))
}

all_units <- seq_len(ncol(dge_refit))
log_msg("Exact reproduction fit (774 units)")
full_exact <- fit_exact(all_units, recompute_rle = FALSE)
same_set <- setequal(full_exact[call == TRUE, gene], ref_calls$gene)
max_dlfc <- max(abs(full_exact$logFC - deg_ref[match(full_exact$gene, gene), logFC]))
summary_list$reproduction <- list(n_calls = sum(full_exact$call), identical_call_set = same_set,
                                  max_abs_logFC_diff = max_dlfc)
if (!same_set || sum(full_exact$call) != N_CALLS) stop_with("S2", "reproduction gives ", sum(full_exact$call), " calls")
log_msg("Reproduction OK: %d calls, max |dlogFC| %.2e", sum(full_exact$call), max_dlfc)
full_fixed <- fit_fixed(all_units, recompute_rle = FALSE)
summary_list$reproduction$fixed_weight_full_calls <- sum(full_fixed$call)
summary_list$reproduction$fixed_vs_exact <- compare_calls(full_exact, full_fixed)

tabP <- scored[[primary_sig]]$tab
rules <- list()
if (scored[[primary_sig]]$T_obs > 0) rules$T <- which(tabP$tail)
rules$F <- unlist(lapply(names(cohort_idx), function(c) {
  ix <- cohort_idx[[c]][tabP$group_binary[cohort_idx[[c]]] == "Disease"]
  k <- ceiling(0.05 * length(ix))
  ix[order(-tabP$rz[ix])][seq_len(k)]
}))

refit_rows <- list(); null_rows <- list(); gene_rows <- list()
for (rn in names(rules)) {
  ex <- rules[[rn]]
  keep <- setdiff(all_units, ex)
  log_msg("Rule %s: excluding %d units; exact fit", rn, length(ex))
  fe <- fit_exact(keep, recompute_rle = TRUE)
  ce <- compare_calls(full_exact, fe)
  ff <- fit_fixed(keep, recompute_rle = TRUE)
  cf <- compare_calls(full_fixed, ff)
  gm_tab <- merge(full_exact[, .(gene, logFC_full = logFC, padj_full = padj, call_full = call)],
                  fe[, .(gene, logFC_excl = logFC, padj_excl = padj, call_excl = call)], by = "gene")
  gm_tab[, `:=`(rule = rn, symbol = gm[match(strip_v(gene), eb), gene_name])]
  gene_rows[[rn]] <- gm_tab[call_full | call_excl]

  # random exclusions matched per cohort x stratum
  need <- table(meta$cell[ex])
  dis_by_cell <- split(which(is_dis), meta$cell[is_dis])
  set.seed(SEED)
  draws <- lapply(seq_len(B_EXCL), function(b) unlist(lapply(names(need), function(cl) {
    pool <- dis_by_cell[[cl]]
    pool[sample.int(length(pool), need[[cl]])]
  })))
  log_msg("Rule %s: %d random-exclusion fixed-weight fits on %d cores", rn, B_EXCL, N_CORES)
  res <- mclapply(draws, function(e) {
    unlist(compare_calls(full_fixed, fit_fixed(setdiff(all_units, e), recompute_rle = TRUE)))
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop("random-exclusion fits failed: ", sum(bad))
  nm <- as.data.table(do.call(rbind, res))[, `:=`(rule = rn, draw = seq_len(B_EXCL))]
  null_rows[[rn]] <- nm
  p_refit <- (1 + sum(nm$changed >= cf$changed)) / (B_EXCL + 1)
  refit_rows[[rn]] <- data.table(
    rule = rn, n_excluded = length(ex), n_units = length(keep),
    exact_n_calls = ce$n_calls, exact_lost = ce$lost, exact_gained = ce$gained, exact_flipped = ce$flipped,
    exact_changed = ce$changed, exact_retained = ce$retained, exact_spearman = ce$spearman_logFC,
    fixed_full_calls = sum(full_fixed$call), fixed_changed = cf$changed, fixed_lost = cf$lost, fixed_gained = cf$gained,
    null_changed_median = median(nm$changed), null_changed_q025 = quantile(nm$changed, 0.025),
    null_changed_q975 = quantile(nm$changed, 0.975), null_lost_median = median(nm$lost),
    null_gained_median = median(nm$gained), p_refit = p_refit)
  log_msg("Rule %s: exact changed %d (lost %d, gained %d, flipped %d); fixed %d vs null median %.1f; p=%.4f",
          rn, ce$changed, ce$lost, ce$gained, ce$flipped, cf$changed, median(nm$changed), p_refit)
}
refit_tab <- rbindlist(refit_rows)
wtsv(refit_tab, "refit_exclusion_summary.tsv")
wtsv(rbindlist(null_rows), "refit_random_exclusion_null.tsv")
wtsv(rbindlist(gene_rows), "refit_exclusion_gene_changes.tsv")
wtsv(data.table(rule = rep(names(rules), lengths(rules)), sample_id = meta$sample_id[unlist(rules)],
                dataset = meta$dataset[unlist(rules)], stratum = meta$stratum[unlist(rules)],
                rz = tabP$rz[unlist(rules)]), "refit_excluded_units.tsv")

# ------------------------------------------------------------------ Holm family
p_tail <- scored[[primary_sig]]$p_tail
fam <- c(tail = p_tail)
if ("T" %in% names(rules)) fam <- c(fam, refit_T = refit_tab[rule == "T", p_refit])
holm <- p.adjust(fam, method = "holm")
summary_list$confirmatory <- list(signature = primary_sig, raw_p = as.list(fam), holm_p = as.list(holm),
                                  family_size = length(fam))
summary_list$refit <- refit_tab
write_json(summary_list, out_path("summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo.txt"))
log_msg("Done.")
