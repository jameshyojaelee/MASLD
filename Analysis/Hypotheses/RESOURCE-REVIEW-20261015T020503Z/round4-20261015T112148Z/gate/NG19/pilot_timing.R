# NG19 gate pilot: time one SuSiE-RSS + coloc.susie simulation cell at one real
# EUR window (FHIP2B, UKBB_ALT; observed PP.H4.susie 0.978) and give a first
# look at recall. Settings copied from GWAS/finemapping/src/06_susie_coloc.R and
# 08_run_eqtl_susie.R (L=10, residual_variance=1 fixed, check_R=FALSE,
# min_abs_corr=0.5 GWAS side, max_iter=500, R + 1e-3 I, eQTL n=1183).
suppressPackageStartupMessages({library(susieR); library(coloc); library(data.table)})
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG19"
COLOC_OVERLAP_MIN <- formals(coloc:::coloc.bf_bf)$overlap.min
source(file.path(FM, "src/coloc_signal_eligibility.R"))
set.seed(20260923)

gene <- "ENSG00000158863"; chr <- 8
t0 <- proc.time()[3]
s_obs <- readRDS(file.path(FM, "results/eqtl_susie/chr8", paste0(gene, "_susie.rds")))
cat("ld_provenance:\n"); print(s_obs$ld_provenance)
pos <- as.integer(sub("^[0-9]+:", "", colnames(s_obs$lbf_variable)))
cat("eQTL window SNPs:", length(pos), " observed eQTL CS:", length(s_obs$sets$cs), "\n")
panel_dir <- file.path(FM, "data/ld_ref", ifelse(is.null(s_obs$ld_provenance$ld_panel) ||
                       s_obs$ld_provenance$ld_panel == "polyfun", "polyfun_eur", "ukbb_eur"))
Sys.setenv(LD_PANEL = ifelse(grepl("polyfun", panel_dir), "polyfun", "ukbb"))
source(file.path(FM, "src/finemapping_functions.R"))
paths <- find_LD_block(chr, min(pos), max(pos), "EUR")
# Memory: read only the window's columns of each block (fixed 16G gate job).
bims <- lapply(paths, function(f) fread(paste0(f, ".bim")))
seen <- integer(0); blocks <- list(); allids <- character(0)
for (b in seq_along(paths)) {
  kb <- which(bims[[b]]$V4 %in% pos & !(bims[[b]]$V4 %in% seen) & !duplicated(bims[[b]]$V4))
  if (!length(kb)) next
  seen <- c(seen, bims[[b]]$V4[kb])
  m <- as.matrix(fread(paste0(paths[b], ".ld"), select = kb))[kb, , drop = FALSE]
  blocks[[length(blocks) + 1]] <- m; allids <- c(allids, paste0(chr, ":", bims[[b]]$V4[kb]))
}
tot <- sum(sapply(blocks, nrow)); R <- matrix(0, tot, tot); off <- 0
for (m in blocks) { ix <- off + seq_len(nrow(m)); R[ix, ix] <- m; off <- off + nrow(m) }
rm(blocks); gc(FALSE); ids <- allids
cat("NA in LD:", sum(is.na(R)), "\n"); R[is.na(R)] <- 0
diag(R) <- diag(R) + 1e-3
dimnames(R) <- list(ids, ids)
cat("matched window p =", length(ids), " LD load s =", round(proc.time()[3] - t0, 1), "\n")
t1 <- proc.time()[3]; U <- tryCatch(chol(R), error = function(e) { ev <- eigen(R, symmetric = TRUE); cat("chol failed; eigen fallback\n"); t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0)))) }); cat("chol s =", round(proc.time()[3] - t1, 1), "\n")

fit <- function(z, n, purity) {
  names(z) <- ids
  args <- list(z = z, R = R, n = n, L = 10L, estimate_residual_variance = FALSE,
               residual_variance = 1, check_R = FALSE, max_iter = 500)
  if (purity) args$min_abs_corr <- 0.5
  set.seed(42); do.call(susie_rss, args)
}
one_run <- function(lam_e, lam_g, shared) {
  j <- sample(length(ids), 1)
  if (shared) k <- j else {
    cand <- which(R[j, ]^2 < 0.1); k <- cand[sample(length(cand), 1)]
  }
  ze <- as.vector(R[, j] * lam_e + crossprod(U, rnorm(length(ids))))
  zg <- as.vector(R[, k] * lam_g + crossprod(U, rnorm(length(ids))))
  te <- proc.time()[3]; se <- fit(ze, 1183, FALSE); te <- proc.time()[3] - te
  tg <- proc.time()[3]; sg <- fit(zg, 337000, TRUE); tg <- proc.time()[3] - tg
  pp4 <- NA; pp3 <- NA; state <- "no_susie_result"
  if (isTRUE(se$converged) && isTRUE(sg$converged) && length(se$sets$cs) > 0 && length(sg$sets$cs) > 0) {
    set.seed(42); res <- tryCatch(coloc.susie(sg, se), error = function(e) NULL)
    if (!is.null(res) && !is.null(res$summary) && nrow(res$summary) > 0 && "PP.H4.abf" %in% names(res$summary)) {
      kept <- as.data.table(res$summary)[!is.na(PP.H4.abf)]
      if (nrow(kept) > 0) { b <- kept[which.max(PP.H4.abf)]; pp4 <- b$PP.H4.abf; pp3 <- b$PP.H3.abf
        state <- ifelse(pp4 > 0.5, "positive", "tested_le_0.5") } else state <- "untestable"
    }
  }
  data.table(lam_e, lam_g, shared, pp4, pp3, state, nE = length(se$sets$cs), nG = length(sg$sets$cs),
             t_eqtl = te, t_gwas = tg)
}
grid <- rbind(data.table(lam = c(4, 4, 4, 8, 8, 8, 20, 20, 20), shared = TRUE),
              data.table(lam = c(8, 8, 8), shared = FALSE))
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) one_run(grid$lam[i], grid$lam[i], grid$shared[i])))
print(res)
fwrite(res, file.path(OUT, "pilot_runs.tsv"), sep = "\t")
cat("median s per run (two fits):", median(res$t_eqtl + res$t_gwas), "\n")
cat("total s:", round(proc.time()[3] - t0, 1), "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
