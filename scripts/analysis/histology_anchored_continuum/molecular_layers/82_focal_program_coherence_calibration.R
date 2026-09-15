#!/usr/bin/env Rscript
# KEY MESSAGE: Calibrate the focal-program tercile coherence contrast against a
# split-permutation null, which is valid for the t-residualized arm because the
# mean channel has already been removed, and is shown but NOT decisive for the
# raw arm because permuting t also destroys the rank-one deflation the raw arm
# is exposed to.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(splines)
})

options(digits = 17, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 82_focal_program_coherence_calibration.R OUTPUT_DIR", call. = FALSE)
out_dir <- args[[1L]]
if (dir.exists(out_dir)) stop("Refusing to overwrite: ", out_dir, call. = FALSE)
dir.create(out_dir, recursive = TRUE)

B <- 2000L
set.seed(20260828)

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

FOCAL <- c(ductular_injury = "hotspot_hepatocytes_48f39dd4d817a10e",
           stromal_ecm     = "hotspot_hepatocytes_f05c535ae5bbc0b9")

reg <- fread(REGISTRY, na.strings = c("", "NA"))[, .(sample_id, dataset, t = fixed_projection_percentile)]
stopifnot(nrow(reg) == 844L)
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

contrast <- function(X, bottom_idx, top_idx) {
  a <- stat(X[bottom_idx, , drop = FALSE]); b <- stat(X[top_idx, , drop = FALSE])
  c(d_lev = unname(b["lev"] - a["lev"]), d_zbar = unname(b["zbar"] - a["zbar"]))
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
  tt <- meta$t
  design <- cbind(1, ns(tt, df = 4))
  resid <- expr - design %*% qr.solve(design, expr)

  cut3 <- quantile(tt, c(1/3, 2/3), names = FALSE, type = 8)
  bot <- which(tt <= cut3[1]); top <- which(tt >= cut3[2])
  n_b <- length(bot); n_t <- length(top); n_all <- length(tt)

  for (nm in names(FOCAL)) {
    gsel <- intersect(genes_of[[nm]], colnames(expr))
    for (arm in c("raw", "t_residualized")) {
      X <- if (arm == "raw") expr else resid
      Xg <- X[, gsel, drop = FALSE]
      obs <- contrast(Xg, bot, top)
      null <- matrix(NA_real_, B, 2L, dimnames = list(NULL, c("d_lev", "d_zbar")))
      for (b in seq_len(B)) {
        perm <- sample.int(n_all)
        null[b, ] <- contrast(Xg, perm[seq_len(n_b)], perm[n_b + seq_len(n_t)])
      }
      for (s in c("d_lev", "d_zbar")) {
        nd <- null[, s]; nd <- nd[is.finite(nd)]
        rows[[length(rows) + 1L]] <- data.table(
          program = nm, dataset = coh, arm = arm, statistic = s,
          n_bottom = n_b, n_top = n_t, g = length(gsel),
          observed = unname(obs[[s]]),
          null_mean = mean(nd), null_sd = sd(nd),
          null_p025 = unname(quantile(nd, 0.025)), null_p975 = unname(quantile(nd, 0.975)),
          z_vs_null = (obs[[s]] - mean(nd)) / sd(nd),
          two_sided_p = (1 + sum(abs(nd - mean(nd)) >= abs(obs[[s]] - mean(nd)))) / (1 + length(nd)),
          finite_draws = length(nd), requested_draws = B)
      }
      draws[[length(draws) + 1L]] <- data.table(
        program = nm, dataset = coh, arm = arm, draw = seq_len(B),
        d_lev = null[, "d_lev"], d_zbar = null[, "d_zbar"])
    }
  }
  message("done ", coh)
}

res <- rbindlist(rows)
res[, bh_q := p.adjust(two_sided_p, method = "BH"), by = .(arm, statistic)]
fwrite(res, file.path(out_dir, "focal_coherence_split_permutation.tsv"), sep = "\t")
fwrite(rbindlist(draws), file.path(out_dir, "focal_coherence_null_draws.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n=== d_lev, observed vs split-permutation null ===\n")
print(res[statistic == "d_lev", .(program, dataset, arm, n_top, g,
        observed = round(observed, 4), null_p025 = round(null_p025, 4),
        null_p975 = round(null_p975, 4), z = round(z_vs_null, 2),
        p = signif(two_sided_p, 3), q = signif(bh_q, 3))], nrows = 100)
