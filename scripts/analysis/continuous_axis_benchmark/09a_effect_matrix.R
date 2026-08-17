#!/usr/bin/env Rscript
# 09a: build the feature-by-window effect matrix that Q2 adjudicates, with a
# standard error for every entry.
#
# The SEs are not decoration. The whole Q2 null is a parametric bootstrap under
# H_amp using this SE matrix, so an effect matrix without per-entry SEs cannot be
# tested for rank-1 sufficiency at all -- you would have no way to ask whether the
# departure from rank 1 exceeds measurement error.
#
# Window schemes:
#   histology  recorded fibrosis stages F0..F4, each contrasted against F0. This is
#              the scheme the Q1 gate routes to when no derived ordering passes.
#   fixed      equal-donor overlapping windows along a supplied axis, each
#              contrasted against window 1.
#
# Every contrast is cumulative (window w versus window 1) rather than adjacent.
# Adjacent contrasts share a donor group with opposite sign, which induces the
# algebraic negative correlation already characterised in
# stage_axis_geometry/.../d2_shared_group_null.tsv. Cumulative contrasts share only
# the reference, which induces a positive bias that the H_amp bootstrap absorbs
# because it is applied identically to observed and simulated matrices.

suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR)
})
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
SCHEME <- Sys.getenv("CAB_WINDOW_SCHEME", "histology")   # histology | fixed
ARM <- Sys.getenv("CAB_ARM", "A4")
N_WIN <- as.integer(Sys.getenv("CAB_N_WINDOWS", "13"))
OVERLAP <- as.numeric(Sys.getenv("CAB_OVERLAP", "0.5"))
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
assert_true(SCHEME %in% c("histology", "fixed"), "CAB_WINDOW_SCHEME must be histology or fixed")
dir.create(file.path(OUT, "q2"), recursive = TRUE, showWarnings = FALSE)

pre <- read_prespec()
set.seed(pre$seeds$master)
meta <- load_manifest()
dge <- readRDS(substrate_path("dge"))[, meta$sample_id]

popC <- population(meta, "POP-C")
info <- copy(popC)
if (SCHEME == "fixed") {
  ax <- fread(file.path(OUT, "arms", sprintf("%s_axis.tsv", ARM)))
  info <- merge(info, ax[, .(sample_id, axis_raw)], by = "sample_id")
  info <- info[is.finite(axis_raw)]
}
info[, dataset := droplevels(factor(dataset))]
log_step(SCHEME, " scheme, arm ", ARM, ", n=", nrow(info))

# Assign each donor to one or more windows.
if (SCHEME == "histology") {
  info[, win := paste0("F", fibrosis_stage)]
  wins <- paste0("F", sort(unique(info$fibrosis_stage)))
  ref <- wins[1]
  membership <- lapply(wins, function(w) info$sample_id[info$win == w])
  names(membership) <- wins
} else {
  info[, r := rank_within(axis_raw, dataset)]
  setorder(info, r)
  n <- nrow(info)
  width <- ceiling(n / (1 + (N_WIN - 1) * (1 - OVERLAP)))
  step <- max(1L, floor(width * (1 - OVERLAP)))
  starts <- seq(1L, max(1L, n - width + 1L), by = step)[seq_len(N_WIN)]
  starts <- starts[!is.na(starts)]
  membership <- lapply(seq_along(starts), function(i)
    info$sample_id[starts[i]:min(n, starts[i] + width - 1L)])
  names(membership) <- sprintf("W%02d", seq_along(membership))
  wins <- names(membership); ref <- wins[1]
  log_step("window width=", width, " step=", step, " n_windows=", length(wins))
}
assert_true(length(wins) >= 3L, "Need at least three windows")

sizes <- data.table(window = wins, n = vapply(membership, length, integer(1)),
                    is_reference = wins == ref)
write_tsv_once(sizes, file.path(OUT, "q2", sprintf("%s_%s_window_sizes.tsv", ARM, SCHEME)))
print(sizes)

# One voom fit per contrast on the donors in the two windows. Fitting each
# contrast on its own donors keeps the SE honest for that comparison rather than
# borrowing precision from windows not involved.
rows <- list()
for (w in setdiff(wins, ref)) {
  ids <- unique(c(membership[[ref]], membership[[w]]))
  d <- info[sample_id %in% ids]
  d[, grp := factor(fifelse(sample_id %in% membership[[w]], "hi", "lo"), levels = c("lo", "hi"))]
  d <- d[!duplicated(sample_id)]
  if (uniqueN(d$grp) < 2L || min(table(d$grp)) < 5L) {
    log_step("  skipping ", w, ": insufficient donors per group"); next
  }
  X <- model.matrix(~ dataset + inferred_sex + grp, data = d)
  X <- X[, qr(X)$pivot[seq_len(qr(X)$rank)], drop = FALSE]
  cf <- grep("^grphi$", colnames(X), value = TRUE)
  if (!length(cf)) { log_step("  skipping ", w, ": contrast not estimable"); next }
  v <- voomWithQualityWeights(dge[, d$sample_id], X, plot = FALSE)
  fit <- eBayes(lmFit(v$E, X, weights = v$weights))
  j <- match(cf, colnames(X))
  se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  rows[[w]] <- data.table(feature = rownames(v$E), window = w,
                          estimate = fit$coefficients[, j], se = se,
                          n_lo = sum(d$grp == "lo"), n_hi = sum(d$grp == "hi"))
  log_step("  ", w, ": n_lo=", sum(d$grp == "lo"), " n_hi=", sum(d$grp == "hi"))
}
M <- rbindlist(rows)
assert_true(uniqueN(M$window) >= 3L, "Fewer than three usable windows after fitting")
assert_true(all(is.finite(M$se) & M$se > 0), "Non-positive or non-finite SE present")

out <- file.path(OUT, "q2", sprintf("%s_%s_M_gene.tsv.gz", ARM, SCHEME))
write_tsv_once(M, out)
log_step("wrote ", nrow(M), " rows over ", uniqueN(M$window), " windows -> ", out)
log_step("EFFECT_MATRIX_COMPLETE scheme=", SCHEME, " arm=", ARM)
