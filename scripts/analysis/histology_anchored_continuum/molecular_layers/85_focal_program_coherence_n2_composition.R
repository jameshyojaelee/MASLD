#!/usr/bin/env Rscript
# KEY MESSAGE: N2 composition-conditional null. Resample t from its conditional
# distribution given the BayesPrism lineage proportions, preserving the
# composition-to-continuum link and breaking only the residual. If the observed
# coherence contrast is a composition artifact, the null reproduces it.
#
# Requested scope is GSE162694, the one cohort whose focal-program contrast
# clears both the split-permutation null and the N1 artifact floor. The other
# four cohorts are computed as labelled context at no extra cost and are NOT the
# requested test.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(splines)
})

options(digits = 17, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 85_focal_program_coherence_n2_composition.R OUTPUT_DIR", call. = FALSE)
out_dir <- args[[1L]]
if (dir.exists(out_dir)) stop("Refusing to overwrite: ", out_dir, call. = FALSE)
dir.create(out_dir, recursive = TRUE)

B <- as.integer(Sys.getenv("HAC_N2_DRAWS", unset = "2000"))
set.seed(20260829)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
REGISTRY <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/candidates",
                      "hac-continuum-20260818T024923Z/projection/score_registry.tsv")
DGE <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates",
                 "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION",
                 "arms/F_five/results/integration/merged_dge.rds")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates",
                   "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE",
                   "frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
SIGEXC <- file.path(ROOT, "RNA-seq/results/histology_anchored_continuum/molecular_layers",
                    "hac-molecular-layers-20260818T173348Z/programs/hotspot",
                    "hotspot_signature_excluded_membership.tsv.gz")
PRISM <- file.path(ROOT, "Analysis/Deconvolution/results")

FOCAL <- c(ductular_injury = "hotspot_hepatocytes_48f39dd4d817a10e",
           stromal_ecm     = "hotspot_hepatocytes_f05c535ae5bbc0b9")
TARGET_COHORT <- "GSE162694"

# Lineages clearing the frozen 1e-4 detection floor in the median sample of
# every cohort, plus Fibroblasts which clears at 1.000 in GSE162694 and is the
# single most fibrosis-relevant lineage. Cholangiocytes is excluded: it clears
# in 0.019 to 0.414 of samples and is unusable in every cohort.
LINEAGES_ALL <- c("Hepatocytes", "Macrophages", "pDCs", "Mono+mono derived cells",
                  "cDC2s", "T cells", "Plasma cells")
LINEAGE_FIBRO <- "Fibroblasts"
FIBRO_COHORTS <- c("GSE130970", "GSE162694", "GSE213621")

# The deposited proportion tables are ragged by one: the header carries 16
# lineage names and each row carries 17 fields, because they were written by R
# write.table with row.names = TRUE. read.table(header = TRUE, row.names = 1) is
# the exact inverse of that writer. The read is then validated by CONTENT, not
# by trusting the column names.
read_prism <- function(cohort) {
  path <- file.path(PRISM, cohort, paste0(cohort, "_bayesprism_proportions.tsv"))
  stopifnot(file.exists(path))
  x <- read.table(path, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(ncol(x) == 16L)
  stopifnot(all(grepl("^[SED]RR[0-9]+$", rownames(x))))       # ids landed in rownames
  stopifnot(all(vapply(x, is.numeric, logical(1))))            # every column numeric
  rs <- rowSums(x)
  stopifnot(all(abs(rs - 1) < 1e-6))                           # rows are a composition
  stopifnot(all(x >= 0))
  as.matrix(x)
}

reg <- fread(REGISTRY, na.strings = c("", "NA"))[, .(sample_id, dataset, t = fixed_projection_percentile)]
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

rows <- list(); draws <- list(); vac <- list()
for (coh in sort(unique(reg$dataset))) {
  meta <- reg[dataset == coh]
  sub <- dge[, colnames(dge$counts) %in% meta$sample_id]
  meta <- meta[match(colnames(sub$counts), sample_id)]
  logcpm <- edgeR::cpm(sub, log = TRUE, prior.count = 1)
  syms <- sym_of[sub("\\..*$", "", rownames(logcpm))]
  ok <- !is.na(syms)
  expr <- scale(t(rowsum(logcpm[ok, , drop = FALSE], syms[ok], reorder = FALSE)))
  n <- nrow(expr)
  prog_cols <- lapply(genes_of, function(g) intersect(g, colnames(expr)))
  keep_cols <- unique(unlist(prog_cols))
  M <- expr[, keep_cols, drop = FALSE]
  cols_by_program <- lapply(prog_cols, function(g) match(g, keep_cols))

  P <- read_prism(coh)
  stopifnot(all(meta$sample_id %in% rownames(P)))
  P <- P[meta$sample_id, , drop = FALSE]
  lin <- LINEAGES_ALL
  if (coh %in% FIBRO_COHORTS) lin <- c(lin, LINEAGE_FIBRO)
  stopifnot(all(lin %in% colnames(P)))
  S <- P[, lin, drop = FALSE]
  # sub-compositional centred log-ratio on the retained lineages
  Slog <- log(pmax(S, 1e-8))
  Z <- Slog - rowMeans(Slog)
  Z <- scale(Z)
  Z <- Z[, apply(Z, 2L, function(v) is.finite(sd(v)) && sd(v) > 0), drop = FALSE]

  tt <- meta$t
  y <- qnorm((rank(tt, ties.method = "average") - 0.5) / n)
  fit <- lm(y ~ Z)
  r2 <- summary(fit)$r.squared
  sigma <- summary(fit)$sigma
  vac[[length(vac) + 1L]] <- data.table(
    dataset = coh, n = n, n_lineages = ncol(Z), fibroblasts_used = LINEAGE_FIBRO %in% lin,
    t_given_Z_r_squared = r2, residual_sd = sigma,
    requested_test = coh == TARGET_COHORT)

  obs <- evaluate(M, tt, cols_by_program)
  fitted_y <- as.vector(fitted(fit))
  null <- vector("list", B)
  for (b in seq_len(B)) {
    ystar <- fitted_y + rnorm(n, 0, sigma)
    tstar <- (rank(ystar, ties.method = "average") - 0.5) / n   # marginal preserved exactly
    null[[b]] <- evaluate(M, tstar, cols_by_program)
  }
  for (key in names(obs)) {
    parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
    for (s in c("d_lev", "d_zbar")) {
      nd <- vapply(null, function(x) unname(x[[key]][[s]]), numeric(1)); nd <- nd[is.finite(nd)]
      o <- unname(obs[[key]][[s]])
      rows[[length(rows) + 1L]] <- data.table(
        dataset = coh, program = parts[[1L]], arm = parts[[2L]], statistic = s,
        requested_test = coh == TARGET_COHORT,
        n_samples = n, g = length(cols_by_program[[parts[[1L]]]]),
        n_lineages = ncol(Z), t_given_Z_r_squared = r2,
        observed = o, null_mean = mean(nd), null_sd = sd(nd),
        null_p025 = unname(quantile(nd, 0.025)), null_p975 = unname(quantile(nd, 0.975)),
        z_vs_null = (o - mean(nd)) / sd(nd),
        two_sided_p = (1 + sum(abs(nd - mean(nd)) >= abs(o - mean(nd)))) / (1 + length(nd)),
        finite_draws = length(nd), requested_draws = B)
      draws[[length(draws) + 1L]] <- data.table(
        dataset = coh, program = parts[[1L]], arm = parts[[2L]], statistic = s,
        draw = seq_along(nd), null_value = nd)
    }
  }
  message("done ", coh, "  n=", n, "  lineages=", ncol(Z), "  R2(t|Z)=", round(r2, 4))
}

res <- rbindlist(rows)
fwrite(res, file.path(out_dir, "focal_coherence_n2_composition.tsv"), sep = "\t")
fwrite(rbindlist(draws), file.path(out_dir, "focal_coherence_n2_draws.tsv.gz"), sep = "\t")
fwrite(rbindlist(vac), file.path(out_dir, "n2_conditioning_vacuity_audit.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n=== VACUITY AUDIT: how much of t does composition explain? ===\n")
print(rbindlist(vac))
cat("\n=== N2, requested test (GSE162694), residualized arm ===\n")
print(res[dataset == TARGET_COHORT & arm == "t_residualized",
          .(program, statistic, g, obs = round(observed, 4),
            null_mean = round(null_mean, 4), null_p025 = round(null_p025, 4),
            null_p975 = round(null_p975, 4), z = round(z_vs_null, 2),
            p = signif(two_sided_p, 3))])
