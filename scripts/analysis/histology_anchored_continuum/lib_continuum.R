suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

fail <- function(...) stop(paste0(...), call. = FALSE)

assert_true <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || !condition) fail(message)
  invisible(TRUE)
}

project_root <- function() {
  root <- Sys.getenv("MASLD_PROJECT_ROOT", unset = "")
  if (!nzchar(root)) {
    root <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
  }
  normalizePath(root, mustWork = TRUE)
}

workstream_dir <- function() {
  file.path(project_root(), "scripts/analysis/histology_anchored_continuum")
}

out_root <- function() {
  root <- Sys.getenv("HAC_OUT_ROOT", unset = "")
  assert_true(nzchar(root), "HAC_OUT_ROOT is unset")
  normalizePath(root, mustWork = TRUE)
}

read_prespec <- function() {
  path <- file.path(workstream_dir(), "00_prespecification.json")
  assert_true(file.exists(path), paste0("Prespecification is missing: ", path))
  fromJSON(path, simplifyVector = TRUE)
}

read_env_path <- function(name, default = NULL, must_exist = TRUE) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value) && !is.null(default)) value <- default
  assert_true(nzchar(value), paste0(name, " is unset"))
  if (must_exist) assert_true(file.exists(value), paste0(name, " does not exist: ", value))
  normalizePath(value, mustWork = must_exist)
}

ensure_new_dir <- function(path) {
  assert_true(!dir.exists(path), paste0("Refusing to overwrite directory: ", path))
  ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
  assert_true(ok && dir.exists(path), paste0("Could not create directory: ", path))
  invisible(path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  assert_true(dir.exists(path), paste0("Could not create directory: ", path))
  invisible(path)
}

write_tsv_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite file: ", path))
  ensure_dir(dirname(path))
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

write_json_once <- function(x, path, pretty = TRUE) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite file: ", path))
  ensure_dir(dirname(path))
  write_json(x, path, auto_unbox = TRUE, pretty = pretty, null = "null", na = "string")
  invisible(path)
}

save_rds_once <- function(x, path, compress = "xz") {
  assert_true(!file.exists(path), paste0("Refusing to overwrite file: ", path))
  ensure_dir(dirname(path))
  saveRDS(x, path, compress = compress)
  invisible(path)
}

sha256_file <- function(path) {
  assert_true(file.exists(path), paste0("Cannot hash absent file: ", path))
  line <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(line, "status")
  assert_true(is.null(status) || identical(status, 0L), paste0("sha256sum failed: ", path))
  sub("[[:space:]].*$", "", line[[1L]])
}

base_gene_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))

percent_rank <- function(x) {
  ans <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (n == 1L) ans[ok] <- 0.5
  if (n > 1L) ans[ok] <- (rank(x[ok], ties.method = "average") - 1) / (n - 1)
  ans
}

cosine_similarity <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L) return(NA_real_)
  den <- sqrt(sum(x[ok]^2)) * sqrt(sum(y[ok]^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(x[ok] * y[ok]) / den
}

orient_pca <- function(scores, loadings, reference_loadings, minimum_overlap = 10L) {
  assert_true(!is.null(names(loadings)), "PCA loadings must be named")
  assert_true(!is.null(names(reference_loadings)), "Reference loadings must be named")
  shared <- intersect(names(loadings), names(reference_loadings))
  assert_true(length(shared) >= minimum_overlap,
              paste0("Too few loadings for orientation: ", length(shared)))
  cosine <- cosine_similarity(loadings[shared], reference_loadings[shared])
  assert_true(is.finite(cosine), "Loading cosine is not finite")
  sign <- if (cosine < 0) -1 else 1
  list(
    scores = as.numeric(scores) * sign,
    loadings = loadings * sign,
    orientation_sign = sign,
    orientation_cosine = abs(cosine),
    n_shared_loadings = length(shared)
  )
}

zscore_rows <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- apply(x, 1L, stats::sd, na.rm = TRUE)
  out <- sweep(x, 1L, means, "-")
  variable <- is.finite(sds) & sds > 0
  out[variable, ] <- sweep(out[variable, , drop = FALSE], 1L, sds[variable], "/")
  out[!variable, ] <- NA_real_
  out
}

standardize_vector <- function(x) {
  if (sum(is.finite(x)) < 2L) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

collapse_symbols <- function(expression_matrix, gene_ids, annotation) {
  ann <- copy(annotation)
  assert_true(all(c("gene_id", "gene_name") %in% names(ann)),
              "Annotation requires gene_id and gene_name")
  ann[, gene_key := base_gene_id(gene_id)]
  ann[, symbol := toupper(trimws(gene_name))]
  ann <- unique(ann[!is.na(symbol) & symbol != "", .(gene_key, symbol)])
  assert_true(ann[, uniqueN(symbol), by = gene_key][V1 > 1L, .N] == 0L,
              "One Ensembl gene maps to multiple symbols")
  map <- ann[match(base_gene_id(gene_ids), gene_key)]
  ok <- !is.na(map$symbol) & map$symbol != ""
  assert_true(any(ok), "No expression rows mapped to symbols")
  collapsed <- rowsum(expression_matrix[ok, , drop = FALSE], map$symbol[ok],
                      reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(map$symbol[ok])[rownames(collapsed)])
  collapsed <- collapsed / counts
  storage.mode(collapsed) <- "double"
  collapsed
}

bootstrap_spearman <- function(x, y, strata, replicates, seed, confidence = 0.95) {
  ok <- is.finite(x) & is.finite(y) & !is.na(strata)
  x <- x[ok]
  y <- y[ok]
  strata <- as.character(strata[ok])
  assert_true(length(x) >= 10L, "Too few observations for bootstrap Spearman")
  observed <- suppressWarnings(cor(x, y, method = "spearman"))
  groups <- split(seq_along(x), strata)
  set.seed(seed)
  boot <- replicate(replicates, {
    # Index into g, never sample(g, ...): a singleton stratum makes g a single
    # integer and sample() then draws from 1:g instead of returning that element,
    # silently resampling arbitrary rows into the bootstrap.
    idx <- unlist(lapply(groups, function(g) g[sample.int(length(g), length(g), replace = TRUE)]),
                  use.names = FALSE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  alpha <- (1 - confidence) / 2
  ci <- quantile(boot[is.finite(boot)], c(alpha, 1 - alpha), names = FALSE, type = 8)
  list(rho = observed, ci_low = ci[[1L]], ci_high = ci[[2L]], n = length(x),
       n_finite_bootstrap = sum(is.finite(boot)))
}

permutation_spearman_greater <- function(x, y, replicates, seed) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  assert_true(length(x) >= 8L, "Too few observations for permutation Spearman")
  observed <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(seed)
  null <- replicate(replicates, suppressWarnings(cor(sample(x), y, method = "spearman")))
  p <- (1 + sum(null >= observed, na.rm = TRUE)) / (1 + sum(is.finite(null)))
  list(rho = observed, p_value = p, n = length(x))
}

fixed_effect_meta <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  if (sum(ok) < 2L) {
    return(list(estimable = FALSE, beta = NA_real_, se = NA_real_, z = NA_real_,
                p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                n_cohorts = sum(ok)))
  }
  w <- 1 / se[ok]^2
  b <- sum(w * beta[ok]) / sum(w)
  s <- sqrt(1 / sum(w))
  z <- b / s
  list(estimable = TRUE, beta = b, se = s, z = z,
       p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
       ci_low = b - qnorm(0.975) * s, ci_high = b + qnorm(0.975) * s,
       n_cohorts = sum(ok))
}

p_adjust_complete_family <- function(p, family_size = 117L, method = "BH") {
  ans <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) ans[ok] <- p.adjust(p[ok], method = method, n = family_size)
  ans
}

write_session_info <- function(path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite file: ", path))
  ensure_dir(dirname(path))
  writeLines(capture.output(sessionInfo()), path)
  invisible(path)
}
