# Debug only (no simulated outcome): why PC-2 fails for ENSG00000090612. The stored production rds has
# niter 526 and converged TRUE, above the max_iter = 500 cap of 08_run_eqtl_susie.R. This reruns the same
# observed-z eQTL fit (same production loader, same inputs as PC-2) at max_iter 500 and at the refit cap
# 5000 (GWAS/finemapping/src/perf/08_refit.sbatch:33, SUSIE_MAX_ITER default 5000) and compares with the rds.
H80 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H80-20261015T130802Z"
Sys.setenv(H80_SCRIPTS = file.path(H80, "scripts"))
source(file.path(H80, "scripts/h80_lib.R"))
SET <- file.path(H80, "run_setup")
W <- fread(file.path(SET, "windows.tsv"))
w <- W[ensembl == "ENSG00000090612"][1]
prod <- load_production_loader()
gene_dt <- readRDS(file.path(SET, "eqtl_ss", paste0(w$ensembl, ".rds")))
ss <- eqtl_ss_from_rows(gene_dt); se <- eqtl_window(gene_dt)
ld <- prod$get_ld_per_locus(ss, LOCUS = w$ensembl, CHR = w$chr, START = se[1], END = se[2], ancestry = "EUR")
ssf <- ld[[1]]; R <- ld[[2]]; if (!is.matrix(R)) R <- as.matrix(R)
z <- ssf$beta / ssf$standard_error; ids <- paste0(w$chr, ":", ssf$position); names(z) <- ids
diag(R) <- diag(R) + LD_REGULARIZE; colnames(R) <- rownames(R) <- ids
so <- readRDS(w$rds)
cat("stored: niter", so$niter, "converged", so$converged, "n_cs", length(so$sets$cs), "\n")
for (mi in c(500L, 5000L)) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(42)
  f <- susie_rss(z = z, R = R, n = EQTL_N, L = SUSIE_L, estimate_residual_variance = FALSE, residual_variance = 1,
                 check_R = FALSE, max_iter = mi)
  cat("max_iter", mi, ": niter", f$niter, "converged", f$converged, "n_cs", length(f$sets$cs),
      "max|PIP diff| vs stored", signif(max(abs(f$pip - so$pip)), 3), "max|alpha diff|", signif(max(abs(f$alpha - so$alpha)), 3), "\n")
}
# how many sampled Arm A / PC-1 windows use an rds that needed more than 500 iterations
chk <- function(paths) vapply(paths, function(p) readRDS(p)$niter, numeric(1))
nA <- chk(unique(W$rds)); P1 <- fread(file.path(SET, "pc1_windows.tsv")); nP <- chk(unique(P1$rds))
cat("Arm A windows' eQTL rds with niter > 500:", sum(nA > 500), "of", length(nA), "genes;", paste(names(nA)[nA > 500], collapse = " "), "\n")
cat("PC-1 windows' eQTL rds with niter > 500:", sum(nP > 500), "of", length(nP), "genes\n")
