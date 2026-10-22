# T3 vocabulary invariance: shared library.
#
# The question this system exists to answer. The v9 discovery ran one program
# vocabulary (the frozen 117 Hotspot modules) through one model and reported 27
# of 113 fibrosis-supported programs. That number has no denominator. Two things
# are missing and they are different:
#
#   1. Is 27/113 better than any gene set of the same size and expression level
#      would have achieved on this substrate? Bulk liver expression is heavily
#      correlated with fibrosis, so a random 200-gene set is NOT expected to
#      score zero. Without the matched-random rate the hit rate is
#      uninterpretable.
#   2. Is the recovered biology a property of the liver or a property of one
#      clustering? A competitor framework (Kamzolas, Nat Metab 2026) used WGCNA
#      alone; this one uses Hotspot alone. Neither can tell the two apart.
#
# The design is therefore: identical model, identical scoring, several
# independent vocabularies, plus a size- and expression-matched random null.
#
# SCORING IDENTITY. The canonical Hotspot projection is a WEIGHTED mean of
# within-cohort gene z-scores; Hallmark and every other vocabulary here have no
# weights. Comparing a weighted vocabulary against unweighted ones confounds
# "better genes" with "has weights". Every vocabulary is therefore scored with
# the UNWEIGHTED definition (score_gene_sets), which is the only definition all
# of them can carry. The Hotspot vocabulary is additionally scored with its
# native weighted definition so the canonical 27/113 can be reproduced exactly
# and the weighting contribution read off directly.
#
# LANGUAGE. Nothing here orders donors or stages. Fibrosis stage and NAS are
# recorded cross-sectional histologic scores used as regressors.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(matrixStats)
})

# ---------------------------------------------------------------------------
# Vocabulary assembly
# ---------------------------------------------------------------------------

# Ensembl -> symbol, built from the same GENCODE table collapse_symbols() uses,
# so a vocabulary carrying ENSG ids is placed in the same symbol space as the
# expression matrix rather than silently losing those genes to coverage.
t3_ensembl_symbol_map <- function(annotation) {
  ann <- copy(annotation)
  if (all(c("gene_id", "gene_name") %in% names(ann))) {
    setnames(ann, c("gene_id", "gene_name"), c("gene", "symbol"))
  }
  ann[, gene_key := sub("\\..*$", "", gene)]
  ann[, symbol := toupper(trimws(symbol))]
  ann <- unique(ann[!is.na(symbol) & symbol != "", .(gene_key, symbol)])
  setNames(ann$symbol, ann$gene_key)
}

t3_resolve_symbols <- function(x, ens2sym) {
  x <- toupper(trimws(as.character(x)))
  is_ens <- grepl("^ENSG[0-9]+", x)
  key <- sub("\\..*$", "", x)
  mapped <- unname(ens2sym[key])
  out <- ifelse(is_ens & !is.na(mapped), mapped, x)
  out[!is.na(out) & out != ""]
}

# Returns one long table: vocabulary, feature_id, gene_symbol, direction.
# `direction` is carried only so the published panels can also be scored in
# their native signed form as a secondary view; the primary score for every
# vocabulary ignores it.
t3_build_vocabularies <- function(hotspot_membership_path, hallmark_gmt_path,
                                  dictionary_membership_path, panel_root,
                                  panel_files, ens2sym) {
  out <- list()

  hs <- fread(hotspot_membership_path)
  stopifnot(all(c("program_uid", "mapped_symbol") %in% names(hs)))
  hs <- hs[!is.na(mapped_symbol) & mapped_symbol != ""]
  out$hotspot <- unique(data.table(
    vocabulary = "hotspot_117",
    feature_id = hs$program_uid,
    gene_symbol = toupper(trimws(hs$mapped_symbol)),
    direction = "up"
  ), by = c("feature_id", "gene_symbol"))

  hm <- read_gmt(hallmark_gmt_path)
  out$hallmark <- unique(data.table(
    vocabulary = "hallmark_50",
    feature_id = hm$feature_id,
    gene_symbol = hm$gene_symbol,
    direction = "up"
  ), by = c("feature_id", "gene_symbol"))

  dict <- fread(dictionary_membership_path)
  stopifnot(all(c("program_id", "gene", "source") %in% names(dict)))
  # The dictionary's `hs` arm is a 116-program snapshot on a DIFFERENT id scheme
  # and is NOT the frozen 117; its HALLMARK arm is an msigdbr snapshot rather
  # than the SHA-pinned GMT. Both are taken from their canonical sources above,
  # so only cnmf / bnmf / wgcna are read from here.
  for (src in c("cnmf", "bnmf", "wgcna")) {
    d <- dict[source == src]
    if (!nrow(d)) fail("Dictionary source absent: ", src)
    out[[src]] <- unique(data.table(
      vocabulary = paste0(src, "_", uniqueN(d$program_id)),
      feature_id = d$program_id,
      gene_symbol = t3_resolve_symbols(d$gene, ens2sym),
      direction = "up"
    ), by = c("feature_id", "gene_symbol"))
  }

  panels <- read_published_panels(panel_root, panel_files)
  out$panels <- unique(data.table(
    vocabulary = "published_panels",
    feature_id = panels$feature_id,
    gene_symbol = panels$gene_symbol,
    direction = panels$direction
  ), by = c("feature_id", "gene_symbol"))

  ans <- rbindlist(out, use.names = TRUE)
  ans <- ans[!is.na(gene_symbol) & gene_symbol != ""]
  ans[]
}

# ---------------------------------------------------------------------------
# Expression-matched random gene sets
# ---------------------------------------------------------------------------

# Deciles of mean logCPM over the whole sealed non-holdout matrix. Outcome-blind
# by construction: no histologic variable enters, and the 1,097-sample matrix
# includes samples that carry no fibrosis or NAS score at all.
t3_expression_deciles <- function(symbol_matrix, n_bins = 10L) {
  mu <- rowMeans(symbol_matrix)
  breaks <- stats::quantile(mu, probs = seq(0, 1, length.out = n_bins + 1L),
                            na.rm = TRUE, names = FALSE)
  breaks[[1]] <- -Inf
  breaks[[length(breaks)]] <- Inf
  bin <- as.integer(cut(mu, breaks = breaks, include.lowest = TRUE, labels = FALSE))
  setNames(bin, rownames(symbol_matrix))
}

# One draw matched on BOTH set size and per-decile composition. Matching size
# alone is not enough: highly expressed genes are less noisy, so a set drawn
# uniformly is a weaker competitor than the real program purely on measurement
# grounds and would flatter every vocabulary.
t3_matched_random_set <- function(observed_bins, pool_by_bin) {
  counts <- tabulate(observed_bins, nbins = length(pool_by_bin))
  idx <- lapply(seq_along(counts), function(b) {
    if (!counts[[b]]) return(integer(0))
    pool <- pool_by_bin[[b]]
    if (length(pool) < counts[[b]]) {
      fail("Expression decile ", b, " has ", length(pool),
           " genes but a matched draw needs ", counts[[b]])
    }
    pool[sample.int(length(pool), counts[[b]])]
  })
  sort(unlist(idx, use.names = FALSE))
}

# R matched draws for every feature, as integer indices into `universe`.
# Seeded per feature so a rerun of one feature reproduces exactly.
t3_matched_random_sets <- function(feature_members, universe, bins, n_draws, seed) {
  pool_by_bin <- split(seq_along(universe), bins[universe])
  pool_by_bin <- pool_by_bin[as.character(seq_len(max(bins, na.rm = TRUE)))]
  pool_by_bin[vapply(pool_by_bin, is.null, logical(1))] <- list(integer(0))
  out <- vector("list", length(feature_members))
  names(out) <- names(feature_members)
  for (f in seq_along(feature_members)) {
    idx <- match(feature_members[[f]], universe)
    idx <- idx[!is.na(idx)]
    obs_bins <- bins[universe[idx]]
    set.seed(seed + f)
    out[[f]] <- lapply(seq_len(n_draws), function(r) {
      t3_matched_random_set(obs_bins, pool_by_bin)
    })
  }
  out
}

# ---------------------------------------------------------------------------
# Fast scoring. Algebraically identical to score_gene_sets() in the undirected
# case; run_tests checks the two agree to 1e-10 on synthetic data.
# ---------------------------------------------------------------------------

t3_cohort_z <- function(symbol_matrix, meta, cohorts, universe) {
  out <- lapply(cohorts, function(co) {
    j <- which(meta$dataset == co)
    Z <- zscore_rows(symbol_matrix[universe, j, drop = FALSE])
    if (anyNA(Z)) {
      fail("A universe gene is constant within cohort ", co,
           "; the universe must be non-constant everywhere it is scored")
    }
    Z
  })
  names(out) <- cohorts
  out
}

t3_zscore_rows_complete <- function(M) {
  mu <- rowMeans(M)
  sdv <- matrixStats::rowSds(M)
  out <- (M - mu) / sdv
  bad <- !is.finite(sdv) | sdv <= 0
  if (any(bad)) out[bad, ] <- NA_real_
  out
}

# sets: list of integer vectors indexing rows of each Z (i.e. `universe`).
t3_score_sets_fast <- function(zlist, sets) {
  n_sets <- length(sets)
  sizes <- lengths(sets)
  if (any(sizes == 0L)) fail("An empty gene set reached the scorer")
  n_genes <- nrow(zlist[[1]])
  G <- sparseMatrix(
    i = rep.int(seq_len(n_sets), sizes),
    j = unlist(sets, use.names = FALSE),
    x = rep.int(1 / sizes, sizes),
    dims = c(n_sets, n_genes)
  )
  blocks <- lapply(zlist, function(Z) t3_zscore_rows_complete(as.matrix(G %*% Z)))
  out <- do.call(cbind, blocks)
  rownames(out) <- names(sets)
  out
}

# ---------------------------------------------------------------------------
# Fitting. Per cohort HC3 on the frozen primary design, then REML +
# Knapp-Hartung across cohorts. The estimator is metafor::rma exactly as
# meta_analyze() uses it; only the plumbing around it is faster, and
# run_tests checks this path reproduces fit_feature_models + meta_analyze.
# ---------------------------------------------------------------------------

T3_PRIMARY_FORMULA <- ~ sex_final + fibrosis_centered + nas_centered
T3_AXIS_TERMS <- c(fibrosis = "fibrosis_centered", nas = "nas_centered")

t3_designs <- function(discovery_meta, cohorts) {
  out <- lapply(cohorts, function(co) {
    md <- copy(discovery_meta[dataset == co])
    md[, fibrosis_centered := fibrosis_stage - mean(fibrosis_stage, na.rm = TRUE)]
    md[, nas_centered := nas_score - mean(nas_score, na.rm = TRUE)]
    X <- model.matrix(T3_PRIMARY_FORMULA, data = md)
    if (qr(X)$rank != ncol(X)) fail("Rank-deficient design in cohort ", co)
    list(X = X, sample_id = md$sample_id, n = nrow(md))
  })
  names(out) <- cohorts
  out
}

# Returns list(fibrosis = list(beta, se), nas = list(beta, se)), each a
# features x cohorts matrix.
t3_cohort_effects_fast <- function(S, designs) {
  cohorts <- names(designs)
  features <- rownames(S)
  empty <- matrix(NA_real_, nrow = length(features), ncol = length(cohorts),
                  dimnames = list(features, cohorts))
  beta <- list(fibrosis = empty, nas = empty)
  se <- list(fibrosis = empty, nas = empty)
  for (co in cohorts) {
    d <- designs[[co]]
    Y <- S[, d$sample_id, drop = FALSE]
    fit <- hc3_fit_matrix(Y, d$X)
    if (!isTRUE(fit$estimable)) fail("Matrix HC3 failed in cohort ", co, ": ", fit$failure_reason)
    for (axis in names(T3_AXIS_TERMS)) {
      term <- T3_AXIS_TERMS[[axis]]
      beta[[axis]][, co] <- fit$coefficients[term, features]
      se[[axis]][, co] <- fit$se[term, features]
    }
  }
  list(beta = beta, se = se)
}

t3_meta_block <- function(beta_mat, se_mat, min_cohorts = 3L, mc_cores = 1L) {
  n <- nrow(beta_mat)
  if (!n) return(data.table())
  # Contiguous chunks, so rbind below restores the original feature order. An
  # interleaved split would silently permute every row against beta_mat.
  n_chunks <- max(1L, min(as.integer(mc_cores), n))
  chunks <- split(seq_len(n), rep(seq_len(n_chunks),
                                  each = ceiling(n / n_chunks), length.out = n))
  stopifnot(identical(unlist(chunks, use.names = FALSE), seq_len(n)))
  worker <- function(ii) {
    do.call(rbind, lapply(ii, function(i) {
      b <- beta_mat[i, ]; s <- se_mat[i, ]
      ok <- is.finite(b) & is.finite(s) & s > 0
      if (sum(ok) < min_cohorts) return(c(NA_real_, NA_real_, NA_real_, NA_real_,
                                          NA_real_, NA_real_, sum(ok), NA_real_))
      f <- tryCatch(metafor::rma(yi = b[ok], sei = s[ok], method = "REML", test = "knha"),
                    error = function(e) NULL)
      if (is.null(f)) return(c(NA_real_, NA_real_, NA_real_, NA_real_,
                               NA_real_, NA_real_, sum(ok), NA_real_))
      c(as.numeric(f$b), f$se, f$ci.lb, f$ci.ub, f$pval, f$tau2, sum(ok),
        sum(sign(b[ok]) == sign(as.numeric(f$b))))
    }))
  }
  res <- if (mc_cores > 1L) parallel::mclapply(chunks, worker, mc.cores = mc_cores)
         else lapply(chunks, worker)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) fail("A meta-analysis worker failed")
  M <- do.call(rbind, res)
  data.table(
    feature_id = rownames(beta_mat),
    beta_meta = M[, 1], se_meta = M[, 2], ci_lower = M[, 3], ci_upper = M[, 4],
    p_value_meta = M[, 5], tau2 = M[, 6], n_cohorts = as.integer(M[, 7]),
    n_same_direction = as.integer(M[, 8])
  )
}

# One vocabulary through the whole model, returning both axes with the BH
# structure run_discovery.R uses: one family per model_kind spanning both axes.
t3_fit_vocabulary <- function(S, designs, mc_cores = 1L) {
  eff <- t3_cohort_effects_fast(S, designs)
  out <- rbindlist(lapply(names(T3_AXIS_TERMS), function(axis) {
    dt <- t3_meta_block(eff$beta[[axis]], eff$se[[axis]], mc_cores = mc_cores)
    dt[, axis := axis]
    dt
  }), use.names = TRUE)
  out[, q_value := p.adjust(p_value_meta, method = "BH", n = .N)]
  out[, family_size := .N]
  out[]
}

# Counts only. Used for the matched-random null and the permutation null, where
# the per-feature estimates are never reported.
t3_support_counts <- function(fit) {
  fit[, .(n_supported = sum(q_value < 0.05, na.rm = TRUE),
          n_estimable = sum(is.finite(p_value_meta))), by = axis]
}

# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------

t3_wilson <- function(k, n, conf = 0.95) {
  if (!n) return(list(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lower = max(0, centre - half), upper = min(1, centre + half))
}

t3_empirical_p_ge <- function(observed, null_values) {
  (1 + sum(null_values >= observed)) / (length(null_values) + 1)
}
