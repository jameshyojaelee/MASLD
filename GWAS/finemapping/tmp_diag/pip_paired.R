#!/usr/bin/env Rscript
# READ-ONLY. Task #15: do the FITS differ between shipped and rebuilt-Gram LD?
#
# lambda_s (task #13) showed susieR's PSD projection does not distort the
# DIAGNOSTIC. It said nothing about the fit: susie_rss does NO PSD check and
# passes R straight to susie_suff_stat, so canonical PIPs were computed on
# genuinely indefinite matrices (min eigenvalue to -0.19). This runs susie_rss
# on both arms over an identical variant set and compares PIPs directly.
SC    <- "/gpfs/commons/home/jameslee/.psd_audit_tmp"
ROOT  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
OUT   <- file.path(ROOT, "tmp_diag", "pip_paired")
RIDGE <- 1e-6
suppressMessages(source(file.path(ROOT, "src", "finemapping_functions.R")))
suppressMessages(library(susieR))

task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
man  <- fread(file.path(SC, "sample_loci.tsv"), colClasses = "character")
if (task > nrow(man)) quit(status = 0)
r <- man[task]
study <- r$study; pop <- r$ld_pop; anc <- r$ancestry; locus <- r$locus
N_tot <- as.numeric(r$N_tot); N_cases <- suppressWarnings(as.numeric(r$N_cases))
if (is.na(N_cases) || N_cases == 0) N_cases <- NA
N_eff <- if (!is.na(N_cases) && (N_tot - N_cases) > 0) 4/(1/N_cases + 1/(N_tot-N_cases)) else N_tot
parts <- strsplit(locus, "\\.")[[1]]; CHR <- as.numeric(parts[1]); BP <- as.numeric(parts[2])
START <- max(1, BP - 0.5e6); END <- BP + 0.5e6
ss_path <- file.path(ROOT, "output", study, paste0(pop, "_0.5Mb"), "ss",
                     paste0(study, "_0.5Mb_", locus, ".txt"))
if (!file.exists(ss_path)) quit(status = 0)
ss <- fread(ss_path)
env_var <- switch(anc, "AFR"="AFR_LD_DIR","AMR"="AMR_LD_DIR","EAS"="EAS_LD_DIR","SAS"="SAS_LD_DIR", NA)
if (is.na(env_var)) quit(status = 0)
ldref <- file.path(ROOT, "data", "ld_ref")
get_arm <- function(d) { do.call(Sys.setenv, setNames(list(file.path(ldref,d)), env_var))
  o <- get_ld_per_locus(copy(ss), locus, CHR, START, END, ancestry = anc); Sys.unsetenv(env_var); o }
A <- get_arm(pop); B <- get_arm(paste0(pop, "_gram"))
if (is.null(A) || is.null(B)) quit(status = 0)
ssA <- A[[1]]; RA <- as.matrix(A[[2]]); ssB <- B[[1]]; RB <- as.matrix(B[[2]])
common <- intersect(ssA$SNP, ssB$SNP)
iA <- match(common, ssA$SNP); iB <- match(common, ssB$SNP)
ssA <- ssA[iA,]; RA <- RA[iA,iA,drop=FALSE]; ssB <- ssB[iB,]; RB <- RB[iB,iB,drop=FALSE]
z <- ssA$beta/ssA$se; p <- length(z)
if (p < 50) quit(status = 0)
stopifnot(max(abs(z - ssB$beta/ssB$se)) == 0)   # identical z in both arms
RA <- RA + RIDGE*diag(p); RB <- RB + RIDGE*diag(p)

fit <- function(M) tryCatch(susie_rss(z = z, R = M, n = N_eff, L = 10,
                                      estimate_residual_variance = FALSE),
                            error = function(e) NULL)
fA <- fit(RA); fB <- fit(RB)
if (is.null(fA) || is.null(fB)) { cat("FIT FAILED\n"); quit(status = 0) }
pA <- fA$pip; pB <- fB$pip
csn <- function(f) if (is.null(f$sets$cs)) 0L else length(f$sets$cs)
lead <- function(pp) common[which.max(pp)]
# credible-set membership overlap (union of all CS variants)
csv <- function(f) if (is.null(f$sets$cs)) character(0) else common[sort(unique(unlist(f$sets$cs)))]
cA <- csv(fA); cB <- csv(fB)
jac <- if (length(union(cA,cB))) length(intersect(cA,cB))/length(union(cA,cB)) else NA_real_
out <- data.table(
  task=task, study=study, ld_pop=pop, ancestry=anc, locus=locus, p=p, N_eff=N_eff,
  min_eig_shipped=min(eigen(RA,symmetric=TRUE,only.values=TRUE)$values),
  min_eig_gram   =min(eigen(RB,symmetric=TRUE,only.values=TRUE)$values),
  n_cs_shipped=csn(fA), n_cs_gram=csn(fB),
  max_pip_shipped=max(pA), max_pip_gram=max(pB),
  lead_shipped=lead(pA), lead_gram=lead(pB), lead_same=identical(lead(pA),lead(pB)),
  pip_cor=suppressWarnings(cor(pA,pB)), pip_max_abs_diff=max(abs(pA-pB)),
  pip_sum_abs_diff=sum(abs(pA-pB)),
  n_pip_gt01_shipped=sum(pA>0.1), n_pip_gt01_gram=sum(pB>0.1),
  n_pip_gt05_shipped=sum(pA>0.5), n_pip_gt05_gram=sum(pB>0.5),
  cs_jaccard=jac)
fwrite(out, file.path(OUT, sprintf("pip_%03d.tsv", task)), sep="\t")
cat(sprintf("task %d %s %s p=%d  pip_cor=%.6f  maxdiff=%.4f  lead_same=%s  cs %d/%d\n",
    task, study, locus, p, out$pip_cor, out$pip_max_abs_diff, out$lead_same,
    out$n_cs_shipped, out$n_cs_gram))
