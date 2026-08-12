suppressPackageStartupMessages({
  library(data.table)
})

fail <- function(...) stop(paste0(...), call. = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
}

sha256_file <- function(path) {
  if (!file.exists(path)) fail("Missing file for SHA256: ", path)
  out <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(out)) fail("sha256sum returned no output for ", path)
  sub("[[:space:]].*$", "", out[[1]])
}

atomic_dir <- function(destination) {
  parent <- dirname(destination)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination)) fail("Refusing to overwrite: ", destination)
  tmp <- paste0(destination, ".tmp.", Sys.getpid())
  if (file.exists(tmp)) fail("Temporary path already exists: ", tmp)
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  tmp
}

publish_dir <- function(tmp, destination) {
  if (file.exists(destination)) fail("Refusing to overwrite: ", destination)
  if (!file.rename(tmp, destination)) fail("Atomic publish failed: ", destination)
  invisible(destination)
}

write_tsv <- function(x, path) {
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

zscore_vector <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) return(rep(NA_real_, length(x)))
  (x - mu) / sig
}

zscore_rows <- function(mat) {
  means <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1L, stats::sd, na.rm = TRUE)
  out <- sweep(mat, 1L, means, "-")
  ok <- is.finite(sds) & sds > 0
  out[ok, ] <- out[ok, , drop = FALSE] / sds[ok]
  out[!ok, ] <- NA_real_
  out
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  rows <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) return(NULL)
    data.table(
      feature_id = fields[[1]],
      gene_symbol = toupper(fields[3:length(fields)])
    )
  })
  unique(rbindlist(rows, use.names = TRUE), by = c("feature_id", "gene_symbol"))
}

read_published_panels <- function(panel_root, panel_files) {
  out <- lapply(names(panel_files), function(label) {
    path <- file.path(panel_root, panel_files[[label]])
    lines <- readLines(path, warn = FALSE)
    lines <- lines[!grepl("^\\s*#", lines)]
    dt <- fread(text = paste(lines, collapse = "\n"))
    if (!"gene_symbol" %in% names(dt)) fail("Missing gene_symbol in ", path)
    if (!"direction" %in% names(dt)) dt[, direction := "up"]
    dt[, .(
      feature_id = label,
      gene_symbol = toupper(trimws(gene_symbol)),
      direction = tolower(trimws(direction))
    )]
  })
  ans <- rbindlist(out)
  if (ans[!direction %in% c("up", "down"), .N] > 0L) {
    fail("Published signature contains a direction other than up/down")
  }
  if (ans[, uniqueN(direction), by = .(feature_id, gene_symbol)][V1 > 1L, .N] > 0L) {
    fail("Published signature assigns conflicting directions to one gene")
  }
  unique(ans, by = c("feature_id", "gene_symbol", "direction"))
}

collapse_symbols <- function(logcpm, gene_ids, annotation) {
  ann <- copy(annotation)
  if (all(c("gene_id", "gene_name") %in% names(ann))) {
    setnames(ann, c("gene_id", "gene_name"), c("gene", "symbol"))
  }
  if (!all(c("gene", "symbol") %in% names(ann))) {
    fail("Gene annotation must contain gene/symbol or gene_id/gene_name")
  }
  ann[, gene_key := sub("\\..*$", "", gene)]
  ann[, symbol := toupper(trimws(symbol))]
  ann <- ann[!is.na(symbol) & symbol != ""]
  ann <- unique(ann[, .(gene_key, symbol)])
  if (ann[, uniqueN(symbol), by = gene_key][V1 > 1L, .N] > 0L) {
    fail("Gene annotation maps one Ensembl key to multiple symbols")
  }
  keys <- sub("\\..*$", "", gene_ids)
  map <- ann[match(keys, gene_key)]
  ok <- !is.na(map$symbol) & map$symbol != ""
  if (!any(ok)) fail("No expression rows mapped to symbols")
  mapped <- logcpm[ok, , drop = FALSE]
  groups <- map$symbol[ok]
  summed <- rowsum(mapped, groups, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(groups)[rownames(summed)])
  collapsed <- summed / counts
  storage.mode(collapsed) <- "double"
  collapsed
}

make_weight_matrix <- function(symbols, membership, feature_ids, weight_col,
                               direction_col = NULL) {
  m <- copy(membership)
  m[, gene_symbol := toupper(trimws(gene_symbol))]
  m <- m[gene_symbol %in% symbols & feature_id %in% feature_ids]
  if (!is.null(direction_col)) {
    m[, signed_weight := get(weight_col) * fifelse(get(direction_col) == "down", -1, 1)]
    weight_col <- "signed_weight"
  }
  m <- m[, .(weight = sum(get(weight_col))), by = .(feature_id, gene_symbol)]
  W <- matrix(0, nrow = length(symbols), ncol = length(feature_ids),
              dimnames = list(symbols, feature_ids))
  if (nrow(m)) W[cbind(match(m$gene_symbol, symbols), match(m$feature_id, feature_ids))] <- m$weight
  W
}

standardize_scores_by_cohort <- function(score_matrix, cohorts) {
  out <- score_matrix
  for (cohort in unique(cohorts)) {
    j <- which(cohorts == cohort)
    out[, j] <- t(apply(score_matrix[, j, drop = FALSE], 1L, zscore_vector))
  }
  out
}

build_composition_scores <- function(composition, sample_meta, testability,
                                     pseudocount) {
  declared <- testability$celltype
  retained <- testability[testability == "testable", celltype]
  assert_true(length(declared) == 22L && length(retained) == 16L,
              "Composition contract must retain 16 of 22 declared lineages")
  assert_true(all(c("sample_id", "dataset", declared) %in% names(composition)),
              "Composition table lacks contracted columns")
  composition <- composition[match(sample_meta$sample_id, sample_id)]
  assert_true(identical(composition$sample_id, sample_meta$sample_id),
              "Composition and sample metadata order mismatch")
  raw <- as.matrix(composition[, ..retained])
  storage.mode(raw) <- "double"
  observed_per_row <- rowSums(is.finite(raw))
  assert_true(all(observed_per_row %in% c(0L, length(retained))),
              "Composition assay missingness must affect complete sample rows")
  assert_true(all(raw[is.finite(raw)] >= 0),
              "Composition proportions must be nonnegative")

  clr <- matrix(NA_real_, nrow = nrow(raw), ncol = ncol(raw),
                dimnames = dimnames(raw))
  available <- observed_per_row == length(retained)
  P <- raw[available, , drop = FALSE]
  P[P == 0] <- pseudocount
  P <- P / rowSums(P)
  clr[available, ] <- log(P) - rowMeans(log(P))
  scores <- matrix(
    NA_real_, nrow = length(declared), ncol = nrow(sample_meta),
    dimnames = list(declared, sample_meta$sample_id)
  )
  scores[retained, ] <- t(clr)[retained, , drop = FALSE]
  scores <- standardize_scores_by_cohort(scores, sample_meta$dataset)

  coverage <- merge(
    testability[, .(feature_id = celltype, contracted_testability = testability)],
    data.table(
      feature_id = declared,
      n_observed = vapply(declared, function(x) sum(is.finite(composition[[x]])), integer(1L))
    ),
    by = "feature_id", all.x = TRUE
  )
  coverage[, `:=`(
    n_expected = nrow(composition),
    gene_coverage = n_observed / nrow(composition),
    testable = contracted_testability == "testable"
  )]

  ratio_lineages <- c(
    cholangiocytes = "Cholangiocytes",
    fibroblasts = "Fibroblasts",
    macrophages = "Macrophages",
    tcells = "T cells"
  )
  assert_true(all(c("Hepatocytes", unname(ratio_lineages)) %in% names(composition)),
              "Five program-lineage composition columns are unavailable")
  ratios <- data.table(sample_id = sample_meta$sample_id)
  denominator <- composition[["Hepatocytes"]]
  denominator[is.finite(denominator) & denominator == 0] <- pseudocount
  for (label in names(ratio_lineages)) {
    numerator <- composition[[ratio_lineages[[label]]]]
    numerator[is.finite(numerator) & numerator == 0] <- pseudocount
    ratios[[paste0("logratio_", label, "_hepatocytes")]] <-
      log(numerator / denominator)
  }
  ratio_columns <- setdiff(names(ratios), "sample_id")
  ratios[!available, (ratio_columns) := NA_real_]
  list(scores = scores, coverage = coverage, logratios = ratios)
}

score_programs <- function(symbol_matrix, sample_meta, membership, registry,
                           coverage_threshold) {
  features <- registry$program_uid
  membership_all <- membership[, .(
    feature_id = program_uid,
    gene_symbol = mapped_symbol,
    weight = as.numeric(original_l1_weight)
  )]
  total_weight <- membership_all[, .(total_weight = sum(weight)), by = feature_id]
  membership <- membership_all[
    !is.na(gene_symbol) & gene_symbol != "",
    .(weight = sum(weight)), by = .(feature_id, gene_symbol)
  ]
  observed_weight <- membership[gene_symbol %in% rownames(symbol_matrix),
                                .(observed_weight = sum(weight)), by = feature_id]
  coverage <- merge(data.table(feature_id = features), total_weight, by = "feature_id", all.x = TRUE)
  coverage <- merge(coverage, observed_weight, by = "feature_id", all.x = TRUE)
  coverage[is.na(observed_weight), observed_weight := 0]
  coverage[, weight_coverage := observed_weight / total_weight]
  coverage[, testable := is.finite(weight_coverage) & weight_coverage >= coverage_threshold]

  primary <- matrix(NA_real_, nrow = length(features), ncol = ncol(symbol_matrix),
                    dimnames = list(features, colnames(symbol_matrix)))
  unweighted <- primary
  rank_score <- primary
  cohort_coverage_rows <- list()

  ranks <- apply(symbol_matrix, 2L, rank, ties.method = "average", na.last = "keep")
  if (!is.matrix(ranks)) ranks <- matrix(ranks, ncol = ncol(symbol_matrix))
  rownames(ranks) <- rownames(symbol_matrix)
  colnames(ranks) <- colnames(symbol_matrix)
  ranks <- ranks / (nrow(symbol_matrix) + 1) - 0.5

  for (cohort in unique(sample_meta$dataset)) {
    j <- which(sample_meta$dataset == cohort)
    Z <- zscore_rows(symbol_matrix[, j, drop = FALSE])
    variable <- rownames(Z)[rowSums(is.finite(Z)) == ncol(Z)]
    for (feature in features) {
      mm <- membership[feature_id == feature & gene_symbol %in% variable]
      retained <- sum(mm$weight)
      total <- coverage[feature_id == feature, total_weight]
      cohort_coverage_rows[[length(cohort_coverage_rows) + 1L]] <- data.table(
        feature_id = feature, cohort = cohort,
        nonconstant_weight = retained,
        cohort_weight_coverage = retained / total,
        cohort_testable = is.finite(retained) && is.finite(total) &&
          retained / total >= coverage_threshold
      )
      if (!coverage[feature_id == feature, testable] ||
          !is.finite(retained) || !is.finite(total) ||
          retained / total < coverage_threshold) next
      idx <- match(mm$gene_symbol, rownames(Z))
      w <- mm$weight / retained
      primary[feature, j] <- as.numeric(crossprod(w, Z[idx, , drop = FALSE]))
      unweighted[feature, j] <- colMeans(Z[idx, , drop = FALSE])
      rank_score[feature, j] <- as.numeric(crossprod(w, ranks[idx, j, drop = FALSE]))
    }
  }

  list(
    primary = standardize_scores_by_cohort(primary, sample_meta$dataset),
    unweighted = standardize_scores_by_cohort(unweighted, sample_meta$dataset),
    weighted_rank = standardize_scores_by_cohort(rank_score, sample_meta$dataset),
    coverage = coverage,
    cohort_coverage = rbindlist(cohort_coverage_rows)
  )
}

score_gene_sets <- function(symbol_matrix, sample_meta, membership, min_genes,
                            coverage_threshold, directional = FALSE) {
  all_features <- unique(membership$feature_id)
  expected <- membership[, .(n_expected = uniqueN(gene_symbol)), by = feature_id]
  observed <- membership[gene_symbol %in% rownames(symbol_matrix),
                         .(n_observed = uniqueN(gene_symbol)), by = feature_id]
  coverage <- merge(data.table(feature_id = all_features), expected, by = "feature_id")
  coverage <- merge(coverage, observed, by = "feature_id", all.x = TRUE)
  coverage[is.na(n_observed), n_observed := 0L]
  coverage[, gene_coverage := n_observed / n_expected]
  coverage[, testable := n_observed >= min_genes & gene_coverage >= coverage_threshold]

  scores <- matrix(NA_real_, nrow = length(all_features), ncol = ncol(symbol_matrix),
                   dimnames = list(all_features, colnames(symbol_matrix)))
  for (cohort in unique(sample_meta$dataset)) {
    j <- which(sample_meta$dataset == cohort)
    Z <- zscore_rows(symbol_matrix[, j, drop = FALSE])
    for (feature in all_features) {
      if (!coverage[feature_id == feature, testable]) next
      mm <- membership[feature_id == feature & gene_symbol %in% rownames(Z)]
      if (!directional) {
        scores[feature, j] <- colMeans(Z[match(mm$gene_symbol, rownames(Z)), , drop = FALSE], na.rm = TRUE)
      } else {
        up <- unique(mm[direction == "up", gene_symbol])
        down <- unique(mm[direction == "down", gene_symbol])
        up_score <- if (length(up)) colMeans(Z[match(up, rownames(Z)), , drop = FALSE], na.rm = TRUE) else 0
        down_score <- if (length(down)) colMeans(Z[match(down, rownames(Z)), , drop = FALSE], na.rm = TRUE) else 0
        scores[feature, j] <- up_score - down_score
      }
    }
  }
  list(scores = standardize_scores_by_cohort(scores, sample_meta$dataset), coverage = coverage)
}

scores_long <- function(scores, meta, view, definition, coverage) {
  dt <- as.data.table(as.table(scores))
  setnames(dt, c("feature_id", "sample_id", "score"))
  dt[, `:=`(feature_id = as.character(feature_id), sample_id = as.character(sample_id))]
  dt <- merge(dt, meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  dt[, `:=`(view = view, score_definition = definition)]
  dt <- merge(dt, coverage, by = "feature_id", all.x = TRUE, sort = FALSE)
  dt[]
}

hc3_fit <- function(y, data, formula) {
  mf <- tryCatch(
    model.frame(formula, data = data, na.action = na.omit),
    error = function(e) e
  )
  if (inherits(mf, "error")) {
    return(list(estimable = FALSE, failure_reason = conditionMessage(mf)))
  }
  yy <- model.response(mf)
  X <- tryCatch(model.matrix(attr(mf, "terms"), mf), error = function(e) e)
  if (inherits(X, "error")) {
    return(list(estimable = FALSE, failure_reason = conditionMessage(X), n = length(yy)))
  }
  design_rank <- qr(X)$rank
  condition_number <- kappa(X)
  if (length(yy) <= ncol(X)) {
    return(list(
      estimable = FALSE, failure_reason = "insufficient_residual_degrees_of_freedom",
      n = length(yy), rank = design_rank, n_columns = ncol(X),
      condition_number = condition_number
    ))
  }
  if (design_rank != ncol(X)) {
    return(list(
      estimable = FALSE, failure_reason = "rank_deficient_design",
      n = length(yy), rank = design_rank, n_columns = ncol(X),
      condition_number = condition_number
    ))
  }
  xtx_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(xtx_inv)) {
    return(list(
      estimable = FALSE, failure_reason = "singular_crossproduct",
      n = length(yy), rank = design_rank, n_columns = ncol(X),
      condition_number = condition_number
    ))
  }
  beta <- as.numeric(xtx_inv %*% crossprod(X, yy))
  names(beta) <- colnames(X)
  resid <- as.numeric(yy - X %*% beta)
  leverage <- rowSums((X %*% xtx_inv) * X)
  scale <- resid^2 / pmax(1 - leverage, 1e-8)^2
  meat <- crossprod(X, X * scale)
  vcov <- xtx_inv %*% meat %*% xtx_inv
  se <- sqrt(pmax(diag(vcov), 0))
  stat <- beta / se
  df <- nrow(X) - ncol(X)
  p <- 2 * stats::pt(-abs(stat), df = df)
  list(
    estimable = TRUE,
    failure_reason = NA_character_,
    coefficients = beta,
    se = setNames(se, colnames(X)),
    statistic = setNames(stat, colnames(X)),
    p_value = setNames(p, colnames(X)),
    n = nrow(X),
    df = df,
    rank = design_rank,
    n_columns = ncol(X),
    condition_number = condition_number,
    max_leverage = max(leverage),
    vcov = vcov
  )
}

fit_feature_models <- function(scores, meta, cohorts, model_kind = "primary") {
  form <- switch(
    model_kind,
    primary = score ~ sex_final + fibrosis_centered + nas_centered,
    interaction = score ~ sex_final + fibrosis_centered + nas_centered + fibrosis_centered:nas_centered,
    age = score ~ sex_final + age_centered + fibrosis_centered + nas_centered,
    composition = score ~ sex_final + fibrosis_centered + nas_centered +
      logratio_cholangiocytes_hepatocytes + logratio_fibroblasts_hepatocytes +
      logratio_macrophages_hepatocytes + logratio_tcells_hepatocytes,
    marginal_fibrosis = score ~ sex_final + fibrosis_centered,
    disease_control = score ~ sex_final + group_binary,
    nash_nafl = score ~ sex_final + diagnosis_binary,
    fail("Unknown model kind: ", model_kind)
  )
  terms <- switch(
    model_kind,
    primary = c(fibrosis = "fibrosis_centered", nas = "nas_centered"),
    interaction = c(interaction = "fibrosis_centered:nas_centered"),
    age = c(fibrosis = "fibrosis_centered", nas = "nas_centered"),
    composition = c(fibrosis = "fibrosis_centered", nas = "nas_centered"),
    marginal_fibrosis = c(fibrosis_marginal = "fibrosis_centered"),
    disease_control = c(disease_control = "group_binaryDisease"),
    nash_nafl = c(nash_nafl = "diagnosis_binaryNASH"),
    fail("Unknown model kind: ", model_kind)
  )

  rows <- list()
  k <- 0L
  for (cohort in cohorts) {
    md <- meta[dataset == cohort]
    if ("fibrosis_stage" %in% names(md)) {
      md[, fibrosis_centered := fibrosis_stage - mean(fibrosis_stage, na.rm = TRUE)]
    }
    if ("nas_score" %in% names(md)) {
      md[, nas_centered := nas_score - mean(nas_score, na.rm = TRUE)]
    }
    if ("age" %in% names(md)) {
      md[, age_centered := age - mean(age, na.rm = TRUE)]
    }
    for (feature in rownames(scores)) {
      md[, score := as.numeric(scores[feature, match(sample_id, colnames(scores))])]
      fit <- hc3_fit(md$score, md, form)
      for (axis in names(terms)) {
        term <- terms[[axis]]
        k <- k + 1L
        estimable <- isTRUE(fit$estimable) && term %in% names(fit$coefficients)
        failure <- if (estimable) NA_character_ else
          if (!is.null(fit$failure_reason)) fit$failure_reason else "term_not_estimable"
        beta <- if (estimable) fit$coefficients[[term]] else NA_real_
        se <- if (estimable) fit$se[[term]] else NA_real_
        residual_df <- if (!is.null(fit$df)) fit$df else NA_integer_
        critical <- if (is.finite(residual_df)) stats::qt(0.975, df = residual_df) else NA_real_
        rows[[k]] <- data.table(
          feature_id = feature,
          cohort = cohort,
          model_kind = model_kind,
          axis = axis,
          term = term,
          estimable = estimable,
          failure_reason = failure,
          beta = beta,
          se_hc3 = se,
          ci_lower = beta - critical * se,
          ci_upper = beta + critical * se,
          statistic_hc3 = if (estimable) fit$statistic[[term]] else NA_real_,
          p_value = if (estimable) fit$p_value[[term]] else NA_real_,
          n = if (!is.null(fit$n)) fit$n else NA_integer_,
          residual_df = residual_df,
          design_rank = if (!is.null(fit$rank)) fit$rank else NA_integer_,
          n_model_columns = if (!is.null(fit$n_columns)) fit$n_columns else NA_integer_,
          condition_number = if (!is.null(fit$condition_number)) fit$condition_number else NA_real_,
          max_leverage = if (!is.null(fit$max_leverage)) fit$max_leverage else NA_real_
        )
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

meta_analyze <- function(cohort_effects, min_cohorts = 3L) {
  if (!requireNamespace("metafor", quietly = TRUE)) fail("metafor is required")
  groups <- unique(cohort_effects[, .(feature_id, axis, model_kind)])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    d <- cohort_effects[
      feature_id == g$feature_id & axis == g$axis & model_kind == g$model_kind &
        estimable == TRUE & is.finite(beta) & is.finite(se_hc3) & se_hc3 > 0
    ]
    if (nrow(d) < min_cohorts) {
      rows[[i]] <- data.table(
        feature_id = g$feature_id, axis = g$axis, model_kind = g$model_kind,
        estimable = FALSE, failure_reason = "fewer_than_three_cohorts",
        n_cohorts = nrow(d)
      )
      next
    }
    fit <- tryCatch(
      metafor::rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "knha"),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      rows[[i]] <- data.table(
        feature_id = g$feature_id, axis = g$axis, model_kind = g$model_kind,
        estimable = FALSE, failure_reason = conditionMessage(fit), n_cohorts = nrow(d)
      )
      next
    }
    rows[[i]] <- data.table(
      feature_id = g$feature_id,
      axis = g$axis,
      model_kind = g$model_kind,
      estimable = TRUE,
      failure_reason = NA_character_,
      beta_meta = as.numeric(fit$b),
      se_meta = fit$se,
      ci_lower = fit$ci.lb,
      ci_upper = fit$ci.ub,
      statistic_meta = fit$zval,
      p_value_meta = fit$pval,
      tau2 = fit$tau2,
      i2 = fit$I2,
      q_heterogeneity = fit$QE,
      p_heterogeneity = fit$QEp,
      n_cohorts = nrow(d),
      n_same_direction = sum(sign(d$beta) == sign(as.numeric(fit$b)))
    )
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

categorical_marginal_means <- function(scores, meta, cohorts) {
  mean_rows <- list()
  contrast_rows <- list()
  mean_index <- 0L
  contrast_index <- 0L
  categorical_formula <- score ~ sex_final + fibrosis_factor + nas_factor
  for (cohort in cohorts) {
    md <- copy(meta[dataset == cohort])
    md[, `:=`(
      sex_final = factor(sex_final),
      fibrosis_factor = factor(fibrosis_stage, levels = sort(unique(fibrosis_stage))),
      nas_factor = factor(nas_score, levels = sort(unique(nas_score)))
    )]
    keep <- complete.cases(md[, .(sex_final, fibrosis_factor, nas_factor)])
    md <- md[keep]
    for (feature in rownames(scores)) {
      y <- as.numeric(scores[feature, match(md$sample_id, colnames(scores))])
      md[, score := y]
      fit <- hc3_fit(md$score, md, categorical_formula)
      estimable <- isTRUE(fit$estimable)
      failure_reason <- if (estimable) NA_character_ else
        if (!is.null(fit$failure_reason)) fit$failure_reason else "categorical_model_not_estimable"
      critical <- if (estimable) stats::qt(0.975, fit$df) else NA_real_
      for (axis in c("fibrosis", "nas")) {
        target <- if (axis == "fibrosis") "fibrosis_factor" else "nas_factor"
        levels_target <- levels(md[[target]])
        marginal_vectors <- vector("list", length(levels_target))
        for (level_index in seq_along(levels_target)) {
          level <- levels_target[[level_index]]
          adjusted_mean <- NA_real_
          adjusted_se <- NA_real_
          if (estimable) {
            nd <- copy(md)
            nd[[target]] <- factor(level, levels = levels_target)
            X_new <- model.matrix(categorical_formula, data = nd)
            X_new <- X_new[, names(fit$coefficients), drop = FALSE]
            marginal_vectors[[level_index]] <- colMeans(X_new)
            adjusted_mean <- sum(marginal_vectors[[level_index]] * fit$coefficients)
            adjusted_se <- sqrt(max(as.numeric(
              t(marginal_vectors[[level_index]]) %*% fit$vcov %*%
                marginal_vectors[[level_index]]
            ), 0))
          }
          mean_index <- mean_index + 1L
          mean_rows[[mean_index]] <- data.table(
            feature_id = feature, cohort = cohort, axis = axis,
            level = as.numeric(as.character(level)),
            adjusted_mean = adjusted_mean, adjusted_se_hc3 = adjusted_se,
            ci_lower = adjusted_mean - critical * adjusted_se,
            ci_upper = adjusted_mean + critical * adjusted_se,
            estimable = estimable, failure_reason = failure_reason,
            n = if (!is.null(fit$n)) fit$n else nrow(md),
            residual_df = if (!is.null(fit$df)) fit$df else NA_integer_,
            design_rank = if (!is.null(fit$rank)) fit$rank else NA_integer_,
            n_model_columns = if (!is.null(fit$n_columns)) fit$n_columns else NA_integer_,
            condition_number = if (!is.null(fit$condition_number))
              fit$condition_number else NA_real_
          )
        }
        if (length(levels_target) < 2L) next
        for (upper_index in 2:length(levels_target)) {
          contrast_index <- contrast_index + 1L
          estimate <- se <- statistic <- p_value <- NA_real_
          if (estimable) {
            contrast <- marginal_vectors[[upper_index]] -
              marginal_vectors[[upper_index - 1L]]
            estimate <- sum(contrast * fit$coefficients)
            se <- sqrt(max(as.numeric(t(contrast) %*% fit$vcov %*% contrast), 0))
            if (is.finite(se) && se > 0) {
              statistic <- estimate / se
              p_value <- 2 * stats::pt(-abs(statistic), df = fit$df)
            } else if (is.finite(estimate) && estimate == 0) {
              statistic <- 0
              p_value <- 1
            } else if (is.finite(estimate)) {
              statistic <- sign(estimate) * Inf
              p_value <- 0
            }
          }
          contrast_rows[[contrast_index]] <- data.table(
            feature_id = feature, cohort = cohort, axis = axis,
            lower_level = as.numeric(as.character(levels_target[[upper_index - 1L]])),
            upper_level = as.numeric(as.character(levels_target[[upper_index]])),
            estimate = estimate, se_hc3 = se,
            ci_lower = estimate - critical * se,
            ci_upper = estimate + critical * se,
            statistic_hc3 = statistic, p_value = p_value,
            estimable = estimable, failure_reason = failure_reason,
            n = if (!is.null(fit$n)) fit$n else nrow(md),
            residual_df = if (!is.null(fit$df)) fit$df else NA_integer_,
            design_rank = if (!is.null(fit$rank)) fit$rank else NA_integer_,
            n_model_columns = if (!is.null(fit$n_columns)) fit$n_columns else NA_integer_,
            condition_number = if (!is.null(fit$condition_number))
              fit$condition_number else NA_real_
          )
        }
      }
    }
  }
  list(
    means = rbindlist(mean_rows, use.names = TRUE, fill = TRUE),
    contrasts = rbindlist(contrast_rows, use.names = TRUE, fill = TRUE)
  )
}

fit_matrix_coefficients <- function(Y, meta, model_kind = "primary") {
  d <- copy(meta)
  d[, fibrosis_centered := fibrosis_stage - mean(fibrosis_stage, na.rm = TRUE)]
  d[, nas_centered := nas_score - mean(nas_score, na.rm = TRUE)]
  form <- switch(
    model_kind,
    primary = ~ sex_final + fibrosis_centered + nas_centered,
    marginal_fibrosis = ~ sex_final + fibrosis_centered,
    disease_control = ~ sex_final + group_binary,
    nash_nafl = ~ sex_final + diagnosis_binary,
    fail("Unknown matrix model kind: ", model_kind)
  )
  X <- model.matrix(form, data = d)
  if (qr(X)$rank != ncol(X)) return(NULL)
  B <- solve(crossprod(X), crossprod(X, t(Y[, match(d$sample_id, colnames(Y)), drop = FALSE])))
  rownames(B) <- colnames(X)
  B
}

complete_feature_vector_samples <- function(scores, feature_ids, sample_ids) {
  if (!length(feature_ids)) fail("Feature-vector completeness requires at least one feature")
  feature_index <- match(feature_ids, rownames(scores))
  sample_index <- match(sample_ids, colnames(scores))
  if (anyNA(feature_index) || anyNA(sample_index)) {
    fail("Feature-vector completeness received unknown feature or sample identifiers")
  }
  selected <- scores[feature_index, sample_index, drop = FALSE]
  colSums(is.finite(selected)) == nrow(selected)
}

empirical_p <- function(observed, permuted, alternative = "greater") {
  if (!is.finite(observed) || any(!is.finite(permuted))) {
    fail("Empirical p-value requires one finite observed statistic and a complete finite null")
  }
  if (alternative == "greater") {
    (1 + sum(permuted >= observed)) / (length(permuted) + 1)
  } else if (alternative == "less") {
    (1 + sum(permuted <= observed)) / (length(permuted) + 1)
  } else {
    (1 + sum(abs(permuted) >= abs(observed))) / (length(permuted) + 1)
  }
}

cosine_similarity <- function(a, b) {
  denom <- sqrt(sum(a^2, na.rm = TRUE) * sum(b^2, na.rm = TRUE))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(a * b, na.rm = TRUE) / denom
}
