#!/usr/bin/env Rscript
# KEY MESSAGE: N1 artifact floor. t is a fixed linear projection of the same
# expression matrix, so Cov(x|t) is a rank-one deflation that reads as a
# coherence change. This builds a null that preserves the observed gene
# covariance EXACTLY, recomputes t inside each draw, and therefore contains the
# deflation channel with zero true coherence change.
#
# The null is a sample-space rotation, not an MVN simulation with a shrunk
# covariance. For a centred n x p matrix X, write X = U Y with U an orthonormal
# basis of the complement of the all-ones vector and Y = U'X. For any orthogonal
# R of size (n-1), X_null = U R Y is still centred and satisfies
# X_null' X_null = Y'R'R Y = X'X. The gene covariance is reproduced to machine
# precision, so Cov(x, t) and the deflation survive intact, and no shrinkage
# intensity has to be chosen for the cohorts where p > n.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(splines)
})

options(digits = 17, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 83_focal_program_coherence_n1.R OUTPUT_DIR", call. = FALSE)
out_dir <- args[[1L]]
if (dir.exists(out_dir)) stop("Refusing to overwrite: ", out_dir, call. = FALSE)
dir.create(out_dir, recursive = TRUE)

B <- as.integer(Sys.getenv("HAC_N1_DRAWS", unset = "2000"))
stopifnot(is.finite(B), B >= 1L)
set.seed(20260828)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CAND <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates",
                  "hac-continuum-20260818T024923Z")
REGISTRY <- file.path(CAND, "projection/score_registry.tsv")
LOADINGS <- file.path(CAND, "projection/fixed_projection_loadings.tsv")
DGE <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates",
                 "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION",
                 "arms/F_five/results/integration/merged_dge.rds")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates",
                   "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE",
                   "frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
SIGEXC <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/molecular_layers",
                    "hac-molecular-layers-20260818T173348Z/programs/hotspot",
                    "hotspot_signature_excluded_membership.tsv.gz")

FOCAL <- c(ductular_injury = "hotspot_hepatocytes_48f39dd4d817a10e",
           stromal_ecm     = "hotspot_hepatocytes_f05c535ae5bbc0b9")

reg <- fread(REGISTRY, na.strings = c("", "NA"))[, .(sample_id, dataset, t = fixed_projection_percentile)]
stopifnot(nrow(reg) == 844L)

load_tab <- fread(LOADINGS, na.strings = c("", "NA"))
load_tab <- load_tab[used_for_fixed_projection %in% TRUE & in_resource %in% TRUE]
load_tab[, gene_symbol := toupper(trimws(gene_symbol))]
stopifnot(nrow(load_tab) == 139L)
message("fixed-projection loadings: ", nrow(load_tab), " genes")

dge <- readRDS(DGE)
ann <- unique(fread(ANNOT, select = c("ensembl_base", "gene_name"))[!is.na(gene_name) & nzchar(gene_name)])
sym_of <- setNames(toupper(trimws(ann$gene_name)), ann$ensembl_base)
se <- fread(SIGEXC, na.strings = c("", "NA"))
genes_of <- lapply(FOCAL, function(uid)
  unique(toupper(trimws(se[program_uid == uid & excluded_signature_gene == FALSE, gene_symbol]))))

stat <- function(mat) {
  keep <- apply(mat, 2L, function(v) is.finite(sd(v)) && sd(v) > 0)
  mat <- mat[, keep, drop = FALSE]
  g <- ncol(mat); n <- nrow(mat)
  if (g < 3L || n < 5L) return(c(lev = NA_real_, zbar = NA_real_))
  R <- suppressWarnings(cor(mat)); R[!is.finite(R)] <- 0; diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  off <- pmin(pmax(R[upper.tri(R)], -0.999999), 0.999999)
  c(lev = max(ev) / sum(ev), zbar = mean(atanh(off)) + 1 / (n - 1))
}

# One pass over a matrix and an axis: tercile split, both arms, both statistics.
evaluate <- function(M, tvec, cols_by_program) {
  cut3 <- quantile(tvec, c(1/3, 2/3), names = FALSE, type = 8)
  bot <- which(tvec <= cut3[1]); top <- which(tvec >= cut3[2])
  design <- cbind(1, ns(tvec, df = 4))
  resid <- M - design %*% qr.solve(design, M)
  out <- list()
  for (nm in names(cols_by_program)) {
    cj <- cols_by_program[[nm]]
    for (arm in c("raw", "t_residualized")) {
      X <- if (arm == "raw") M else resid
      a <- stat(X[bot, cj, drop = FALSE]); b <- stat(X[top, cj, drop = FALSE])
      out[[paste0(nm, "|", arm)]] <- c(d_lev = unname(b["lev"] - a["lev"]),
                                       d_zbar = unname(b["zbar"] - a["zbar"]))
    }
  }
  out
}

rows <- list(); draws <- list()
for (coh in sort(unique(reg$dataset))) {
  meta <- reg[dataset == coh]
  sub <- dge[, colnames(dge$counts) %in% meta$sample_id]
  meta <- meta[match(colnames(sub$counts), sample_id)]
  logcpm <- edgeR::cpm(sub, log = TRUE, prior.count = 1)
  syms <- sym_of[sub("\\..*$", "", rownames(logcpm))]
  ok <- !is.na(syms)
  expr <- scale(t(rowsum(logcpm[ok, , drop = FALSE], syms[ok], reorder = FALSE)))
  n <- nrow(expr)

  sig_genes <- intersect(load_tab$gene_symbol, colnames(expr))
  w <- load_tab[match(sig_genes, gene_symbol), discovery_loading_oriented]
  prog_cols <- lapply(genes_of, function(g) intersect(g, colnames(expr)))
  keep_cols <- unique(c(sig_genes, unlist(prog_cols)))
  M <- expr[, keep_cols, drop = FALSE]
  sig_j <- match(sig_genes, keep_cols)
  cols_by_program <- lapply(prog_cols, function(g) match(g, keep_cols))

  # observed
  obs <- evaluate(M, meta$t, cols_by_program)

  # orthonormal basis of the complement of the all-ones vector
  U <- qr.Q(qr(cbind(1, matrix(rnorm(n * (n - 1L)), n, n - 1L))))[, -1L, drop = FALSE]
  Y <- crossprod(U, M)                       # (n-1) x p, M = U Y since M is centred
  stopifnot(max(abs(U %*% Y - M)) < 1e-8)

  null <- vector("list", B)
  for (b in seq_len(B)) {
    R <- qr.Q(qr(matrix(rnorm((n - 1L)^2), n - 1L)))
    Mb <- U %*% (R %*% Y)                    # same gene covariance, exactly
    tb <- rank(as.vector(Mb[, sig_j, drop = FALSE] %*% w), ties.method = "average") / n
    null[[b]] <- evaluate(Mb, tb, cols_by_program)
  }

  for (key in names(obs)) {
    parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
    for (s in c("d_lev", "d_zbar")) {
      nd <- vapply(null, function(x) unname(x[[key]][[s]]), numeric(1))
      nd <- nd[is.finite(nd)]
      o <- unname(obs[[key]][[s]])
      rows[[length(rows) + 1L]] <- data.table(
        program = parts[[1L]], arm = parts[[2L]], statistic = s, dataset = coh,
        n_samples = n, g = length(cols_by_program[[parts[[1L]]]]),
        n_signature_genes = length(sig_genes),
        observed = o, floor_mean = mean(nd), floor_sd = sd(nd),
        floor_p025 = unname(quantile(nd, 0.025)), floor_p975 = unname(quantile(nd, 0.975)),
        floor_p95_abs = unname(quantile(abs(nd - mean(nd)), 0.95)),
        z_vs_floor = (o - mean(nd)) / sd(nd),
        two_sided_p = (1 + sum(abs(nd - mean(nd)) >= abs(o - mean(nd)))) / (1 + length(nd)),
        inside_floor = abs(o - mean(nd)) <= unname(quantile(abs(nd - mean(nd)), 0.95)),
        finite_draws = length(nd), requested_draws = B)
      draws[[length(draws) + 1L]] <- data.table(
        program = parts[[1L]], arm = parts[[2L]], statistic = s, dataset = coh,
        draw = seq_along(nd), null_value = nd)
    }
  }
  message("done ", coh, "  n=", n, "  signature genes=", length(sig_genes))
}

res <- rbindlist(rows)
res[, bh_q := p.adjust(two_sided_p, method = "BH"), by = .(arm, statistic)]
fwrite(res, file.path(out_dir, "focal_coherence_n1_floor.tsv"), sep = "\t")
fwrite(rbindlist(draws), file.path(out_dir, "focal_coherence_n1_draws.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n=== N1 ARTIFACT FLOOR, d_lev ===\n")
print(res[statistic == "d_lev", .(program, dataset, arm, n_samples, g,
        obs = round(observed, 4), floor_mean = round(floor_mean, 4),
        floor_p025 = round(floor_p025, 4), floor_p975 = round(floor_p975, 4),
        z = round(z_vs_floor, 2), p = signif(two_sided_p, 3),
        q = signif(bh_q, 3), inside = inside_floor)], nrows = 100)
