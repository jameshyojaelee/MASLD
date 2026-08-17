#!/usr/bin/env Rscript
# 04: Q1a. Can an axis predict histology in a cohort it has never seen?
#
# Two design decisions that matter more than the model choice.
#
# The axis is ranked WITHIN the held-out cohort before scoring. Cohorts differ in
# mean expression for reasons that have nothing to do with disease, so an unranked
# axis would be graded largely on its ability to recover cohort identity.
#
# The comparator is arm A3, the frozen prespecified NMF axis, not a null. The
# question is not "does the axis carry any information" -- almost anything derived
# from expression will. The question is whether it beats the continuous readout the
# project already owns. The gate therefore requires the lower bound of a donor
# bootstrap CI to sit above A3's point estimate, not merely above zero.

suppressPackageStartupMessages({ library(data.table); library(MASS) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
dir.create(file.path(OUT, "q1"), recursive = TRUE, showWarnings = FALSE)

pre <- read_prespec()
N_BOOT <- pre$q1_gate$an_arm_passes_only_if_all_four_hold$g2_held_out_prediction$n_bootstrap
set.seed(pre$seeds$master)
meta <- load_manifest()

c_index <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  p <- pred[ok]; t <- truth[ok]
  if (length(p) < 10L || length(unique(t)) < 2L) return(NA_real_)
  cc <- 0; dc <- 0
  for (i in seq_along(t)) {
    dt <- t[-i] - t[i]; dp <- p[-i] - p[i]
    cmp <- dt != 0
    cc <- cc + sum(sign(dt[cmp]) == sign(dp[cmp]))
    dc <- dc + sum(sign(dt[cmp]) == -sign(dp[cmp]))
  }
  if ((cc + dc) == 0) return(NA_real_)
  0.5 * (((cc - dc) / (cc + dc)) + 1)
}

qwk <- function(pred_class, truth) {
  ok <- is.finite(pred_class) & is.finite(truth)
  a <- factor(pred_class[ok], levels = sort(unique(c(pred_class[ok], truth[ok]))))
  b <- factor(truth[ok], levels = levels(a))
  O <- table(a, b); n <- sum(O)
  if (n == 0 || nrow(O) < 2) return(NA_real_)
  w <- outer(seq_len(nrow(O)), seq_len(ncol(O)), function(i, j) (i - j)^2) /
    (nrow(O) - 1)^2
  E <- outer(rowSums(O), colSums(O)) / n
  1 - sum(w * O) / sum(w * E)
}

arms <- c("A1", "A2", "A3", "A4", "A5", "A6")
available <- arms[file.exists(file.path(OUT, "arms", paste0(arms, "_axis.tsv")))]
assert_true("A3" %in% available, "Arm A3 (frozen NMF comparator) must exist before Q1a")

load_axis <- function(arm) {
  ax <- fread(file.path(OUT, "arms", paste0(arm, "_axis.tsv")))
  ax[, .(sample_id, axis_raw)]
}
a3 <- setnames(load_axis("A3"), "axis_raw", "nmf_axis")

evaluate <- function(arm, pop_id) {
  d <- merge(population(meta, pop_id), load_axis(arm), by = "sample_id")
  d <- merge(d, a3, by = "sample_id", all.x = TRUE)
  d <- d[is.finite(axis_raw) & !is.na(fibrosis_stage)]
  if (!nrow(d)) return(NULL)
  d[, dataset := droplevels(factor(dataset))]
  d[, axis_rank := rank_within(axis_raw, dataset)]
  d[, nmf_rank := rank_within(nmf_axis, dataset)]
  folds <- levels(d$dataset)
  if (length(folds) < 2L) return(NULL)

  per_fold <- rbindlist(lapply(folds, function(f) {
    tr <- d[dataset != f]; te <- d[dataset == f]
    if (nrow(te) < 10L || uniqueN(te$fibrosis_stage) < 2L) return(NULL)
    ci <- c_index(te$axis_rank, te$fibrosis_stage)
    rho <- suppressWarnings(cor(te$axis_rank, te$fibrosis_stage, method = "spearman"))
    k <- NA_real_
    fit <- try(polr(factor(fibrosis_stage) ~ axis_rank + inferred_sex, data = tr,
                    Hess = TRUE), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      pc <- suppressWarnings(predict(fit, newdata = te, type = "class"))
      k <- qwk(as.numeric(as.character(pc)), te$fibrosis_stage)
    }
    # Nested increment: does the axis add to the frozen NMF axis, in this fold?
    dnmf <- te[is.finite(nmf_rank)]
    d_ci <- NA_real_
    if (nrow(dnmf) >= 20L) {
      d_ci <- c_index(dnmf$axis_rank, dnmf$fibrosis_stage) -
        c_index(dnmf$nmf_rank, dnmf$fibrosis_stage)
    }
    data.table(arm = arm, population = pop_id, fold = f, n = nrow(te),
               c_index = ci, spearman = rho, qwk = k, delta_c_vs_nmf = d_ci)
  }))
  if (!nrow(per_fold)) return(NULL)

  # Pooled C-index with a donor bootstrap. The bootstrap, not the 3-4 folds,
  # carries the inference: four folds cannot support a meta-analysis.
  pooled <- c_index(d$axis_rank, d$fibrosis_stage)
  bs <- replicate(N_BOOT, {
    i <- sample.int(nrow(d), replace = TRUE)
    c_index(d$axis_rank[i], d$fibrosis_stage[i])
  })
  list(per_fold = per_fold,
       overall = data.table(
         arm = arm, population = pop_id, n_donors = nrow(d), n_folds = nrow(per_fold),
         pooled_c_index = pooled,
         boot_lo = quantile(bs, 0.025, na.rm = TRUE),
         boot_hi = quantile(bs, 0.975, na.rm = TRUE),
         mean_fold_c_index = mean(per_fold$c_index, na.rm = TRUE),
         mean_fold_qwk = mean(per_fold$qwk, na.rm = TRUE),
         mean_delta_c_vs_nmf = mean(per_fold$delta_c_vs_nmf, na.rm = TRUE)))
}

folds_all <- list(); overall_all <- list()
for (arm in available) for (pop in c("POP-C", "POP-B")) {
  r <- evaluate(arm, pop)
  if (is.null(r)) next
  folds_all[[paste(arm, pop)]] <- r$per_fold
  overall_all[[paste(arm, pop)]] <- r$overall
  log_step(sprintf("%s %s pooled C=%.4f [%.4f, %.4f]", arm, pop,
                   r$overall$pooled_c_index, r$overall$boot_lo, r$overall$boot_hi))
}
per_fold <- rbindlist(folds_all); overall <- rbindlist(overall_all)

# Gate g2 is judged on POP-C against A3's point estimate.
a3_point <- overall[arm == "A3" & population == "POP-C", pooled_c_index]
overall[, a3_reference_c_index := if (length(a3_point)) a3_point[1] else NA_real_]
overall[, gate_g2_pass := population == "POP-C" & arm != "A3" &
          is.finite(boot_lo) & is.finite(a3_reference_c_index) &
          boot_lo > a3_reference_c_index]

write_tsv_once(per_fold, file.path(OUT, "q1", "loco_per_fold.tsv"))
write_tsv_once(overall, file.path(OUT, "q1", "loco_overall.tsv"))
print(overall[, .(arm, population, n_donors, pooled_c_index, boot_lo, boot_hi,
                  a3_reference_c_index, gate_g2_pass)])
log_step("Q1_LOCO_COMPLETE")
