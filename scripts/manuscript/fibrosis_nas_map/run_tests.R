#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "analysis_lib.R"))

tol <- 1e-9

# z scoring and zero variance
z <- zscore_vector(c(1, 2, 3))
stopifnot(abs(mean(z)) < tol, abs(stats::sd(z) - 1) < tol)
stopifnot(all(is.na(zscore_vector(c(2, 2, 2)))))

# HC3 checked against a direct matrix implementation on a hand-computable model
d <- data.table(y = c(1, 2, 1, 4, 5, 7), sex = factor(c("F", "M", "F", "M", "F", "M")), x = 0:5)
fit <- hc3_fit(d$y, d, y ~ sex + x)
stopifnot(fit$estimable)
X <- model.matrix(~ sex + x, d)
b <- solve(crossprod(X), crossprod(X, d$y))
e <- d$y - X %*% b
h <- diag(X %*% solve(crossprod(X)) %*% t(X))
v <- solve(crossprod(X)) %*% crossprod(X, X * as.numeric(e)^2 / (1 - h)^2) %*% solve(crossprod(X))
stopifnot(max(abs(fit$coefficients - as.numeric(b))) < tol)
stopifnot(max(abs(fit$se - sqrt(diag(v)))) < tol)
bad_fit <- hc3_fit(c(1, 2, 3), data.table(y = c(1, 2, 3), x = c(1, 1, 1)), y ~ x)
stopifnot(!bad_fit$estimable, bad_fit$failure_reason == "rank_deficient_design")

# Ensembl-version stripping and duplicate-symbol mean collapse
raw <- matrix(c(1, 3, 2, 4, 10, 20), nrow = 3, byrow = TRUE,
              dimnames = list(c("ENSG1.1", "ENSG2.2", "ENSG3.1"), c("s1", "s2")))
annotation <- data.table(
  gene_id = c("ENSG1", "ENSG2", "ENSG3"),
  gene_name = c("DUP", "DUP", "SOLO")
)
collapsed <- collapse_symbols(raw, rownames(raw), annotation)
stopifnot(all.equal(unname(collapsed["DUP", ]), c(1.5, 3.5)))

# Program weighting, coverage, and rank sensitivity
expr <- matrix(c(
  1, 2, 3, 4,
  2, 5, 1, 6,
  5, 4, 3, 2,
  8, 8, 8, 8
), nrow = 4, byrow = TRUE,
dimnames = list(c("A", "B", "C", "CONST"), paste0("s", 1:4)))
meta <- data.table(sample_id = colnames(expr), dataset = rep(c("c1", "c2"), each = 2))
membership <- data.table(
  program_uid = rep(c("p1", "p2"), each = 2),
  mapped_symbol = c("A", "B", "A", "MISSING"),
  original_l1_weight = c(0.25, 0.75, 0.1, 0.9)
)
registry <- data.table(program_uid = c("p1", "p2"))
sc <- score_programs(expr, meta, membership, registry, 0.8)
stopifnot(sc$coverage[feature_id == "p1", testable])
stopifnot(!sc$coverage[feature_id == "p2", testable])
stopifnot(all(is.na(sc$primary["p2", ])))
stopifnot(all(is.finite(sc$primary["p1", ])))
stopifnot(all(is.finite(sc$unweighted["p1", ])))
stopifnot(all(is.finite(sc$weighted_rank["p1", ])))
stopifnot(nrow(sc$cohort_coverage) == 4L)

# Composition CLR, structural untestability, and four lineage log-ratios
declared <- c(
  "Hepatocytes", "Cholangiocytes", "Fibroblasts", "Macrophages", "T cells",
  paste0("Lineage", 6:22)
)
testability <- data.table(
  celltype = declared,
  testability = c(rep("testable", 16), rep("untestable_structurally_unavailable", 6))
)
comp <- data.table(sample_id = paste0("s", 1:4), dataset = rep(c("c1", "c2"), each = 2))
for (j in seq_along(declared)) comp[[declared[[j]]]] <- if (j <= 16) (1:4 + j) else NA_real_
comp_meta <- data.table(sample_id = comp$sample_id, dataset = comp$dataset)
composition_scores <- build_composition_scores(comp, comp_meta, testability, 1e-6)
stopifnot(nrow(composition_scores$scores) == 22L)
stopifnot(sum(composition_scores$coverage$testable) == 16L)
stopifnot(ncol(composition_scores$logratios) == 5L)
composition_manifest <- comp[, .(
  partition = "composition_assay",
  n_samples = sum(rowSums(!is.na(.SD)) > 0L)
), by = dataset, .SDcols = declared[1:16]]
stopifnot(identical(names(composition_manifest), c("dataset", "partition", "n_samples")))
partition_manifest_test <- rbindlist(list(
  data.table(partition = "nonholdout", dataset = "c1", n_samples = 2L),
  composition_manifest
), use.names = TRUE)
stopifnot(
  partition_manifest_test[partition == "composition_assay" & dataset == "c1", n_samples] == 2L
)

# Effect-vector transport excludes whole samples missing any compared feature,
# while ignoring untestable rows that are not part of the reference vector.
transport_scores <- matrix(
  c(1, 2, NA, 4, 5, 6, NA, 8, NA, NA, NA, NA),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("tested_a", "tested_b", "untestable"), paste0("s", 1:4))
)
stopifnot(identical(
  unname(
  complete_feature_vector_samples(
    transport_scores, c("tested_a", "tested_b"), colnames(transport_scores)
  )),
  c(TRUE, TRUE, FALSE, TRUE)
))

axis_regression <- data.table(axis = c("fibrosis", "nas"), estimate = c(0.1, 0.2))
axis_name <- "fibrosis"
stopifnot(identical(axis_regression[axis == axis_name, estimate], 0.1))

# Categorical marginal means carry HC3 uncertainty and adjacent contrasts.
toy <- CJ(
  fibrosis_stage = 0:2, nas_score = 0:2,
  sex_final = c("F", "M"), replicate = 1:4
)
toy[, `:=`(
  sample_id = paste0("toy", .I),
  dataset = "toy_cohort"
)]
toy_pattern <- c(`0` = 0, `1` = 3, `2` = -2)
toy_noise <- c(-0.20, -0.05, 0.08, 0.17)
toy_score <- toy_pattern[as.character(toy$fibrosis_stage)] +
  0.15 * toy$nas_score + 0.10 * (toy$sex_final == "M") +
  toy_noise[toy$replicate]
toy_scores <- matrix(
  toy_score, nrow = 1L,
  dimnames = list("toy_feature", toy$sample_id)
)
toy_linearity <- categorical_marginal_means(
  toy_scores, toy, "toy_cohort"
)
stopifnot(nrow(toy_linearity$means) == 6L)
toy_fibrosis <- toy_linearity$contrasts[axis == "fibrosis"]
stopifnot(
  nrow(toy_fibrosis) == 2L,
  all(toy_fibrosis$estimable),
  all(is.finite(toy_fibrosis$p_value)),
  identical(sign(toy_fibrosis$estimate), c(1, -1))
)

# Empirical p-value uses the locked add-one rule
stopifnot(abs(empirical_p(3, c(0, 1, 2, 3, 4)) - 3 / 6) < tol)
bh <- p.adjust(c(0.01, 0.02, NA, 0.5), method = "BH", n = 4L)
stopifnot(length(bh) == 4L, is.na(bh[[3L]]), abs(bh[[1L]] - 0.04) < tol)

set.seed(20260811L)
labels <- data.table(fibrosis = 1:6, nas = 11:16)
permutation <- sample.int(nrow(labels))
shuffled <- labels[permutation]
stopifnot(all(shuffled$nas - shuffled$fibrosis == 10))

# Cosine and paired delta arithmetic
stopifnot(abs(cosine_similarity(c(1, 0), c(2, 0)) - 1) < tol)
beta_f <- c(a = 1, b = -1)
beta_n <- c(a = 0.5, b = 0.5)
pred <- beta_f * 2 + beta_n * -1
stopifnot(all.equal(unname(pred), c(1.5, -2.5)))

cat("All fibrosis-NAS map unit tests passed\n")
