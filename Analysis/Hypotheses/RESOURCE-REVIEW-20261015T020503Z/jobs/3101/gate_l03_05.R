# Feasibility gate for L03-05 (agent gate-L03-05).
# Computes ONLY: (1) positive control COL1A1/COL1A2 fibrosis slopes per cohort and REML meta,
# (2) call-gene coverage after filterByExpr per cohort,
# (3) null variance of the same-sign rate (within-cohort fibrosis-label permutation) to get the
#     design effect from gene correlation. It does NOT compute the observed same-sign rates.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(metafor)})
set.seed(20261015)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out  <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3101"
p <- function(...) file.path(root, ...)
B <- 1000L

L <- fread(p("RNA-seq/results/manuscript_release/candidates/bulk-loco-summary-20261009T143029Z/loco_gene_dependence.tsv"))
L[, gid := sub("\\.\\d+$", "", gene)]
L[, grp := fifelse(n_cohort_removals_losing_call == 0, "stable",
               fifelse(dependent_on == "GSE213621", "dep213621", "other"))]
print(table(L$grp))
refsign <- setNames(sign(L$full_logFC), L$gid)
grp <- setNames(L$grp, L$gid)

read_mat <- function(path, sep) {
  d <- fread(path, sep = sep, data.table = FALSE, check.names = FALSE)
  rn <- sub("\\.\\d+$", "", as.character(d[[1]]))
  m <- as.matrix(d[, -1, drop = FALSE]); rownames(m) <- rn
  m[!duplicated(rn), , drop = FALSE]
}

coh <- list()
ph <- fread(p("data/GSE281797/metadata/sample_phenotypes_resolved.tsv"))
m <- read_mat(p("data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz"), "\t")
coh$GSE281797 <- list(cnt = m[, as.character(ph$matrix_column)], fib = as.numeric(ph$fibrosis), sex = ph$sex_expressed)

ph <- fread(p("data/GSE225740/metadata/sample_phenotypes_resolved.tsv"))
ph <- ph[!is.na(fibrosis)]
m <- read_mat(p("data/GSE225740/GSE225740_counts_ensembl_stable.tsv.gz"), "\t")
coh$GSE225740 <- list(cnt = m[, as.character(ph$matrix_column)], fib = as.numeric(ph$fibrosis), sex = ph$sex_expressed)

ph <- fread(p("Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"))
m <- read_mat(p("Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE269412_RNA.txt.gz"), "\t")
coh$GSE267145 <- list(cnt = m[, ph$participant_id], fib = as.numeric(ph$fibrosis), sex = ph$sex)

ph <- fread(p("RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T020627Z/external/kamzolas_v1/data/metadata.csv"))
ph <- ph[Dataset == "UCAM"]
m <- read_mat(p("RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T020627Z/external/kamzolas_v1/data/merged_counts.csv"), ",")
coh$EMTAB9815 <- list(cnt = m[, ph[["Sample name"]]], fib = as.numeric(ph$Fibrosis), sex = ph$Sex)

pc_ids <- c(COL1A1 = "ENSG00000108821", COL1A2 = "ENSG00000164692")
pc_rows <- list(); cov_rows <- list(); fits <- list()
for (nm in names(coh)) {
  x <- coh[[nm]]
  stopifnot(ncol(x$cnt) == length(x$fib), !anyNA(x$fib), !anyNA(x$sex))
  cat(nm, "n =", ncol(x$cnt), "fib table:", paste(names(table(x$fib)), table(x$fib), sep = ":", collapse = " "),
      "sex:", paste(names(table(x$sex)), table(x$sex), sep = ":", collapse = " "), "\n")
  y <- DGEList(x$cnt)
  des <- model.matrix(~ fib + sex, data = data.frame(fib = x$fib, sex = factor(x$sex)))
  keep <- filterByExpr(y, des)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  v <- voom(y, des)
  f <- eBayes(lmFit(v, des))
  se <- sqrt(f$s2.post) * f$stdev.unscaled[, "fib"]
  for (g in names(pc_ids)) {
    id <- pc_ids[[g]]
    if (id %in% rownames(f)) {
      pc_rows[[length(pc_rows) + 1]] <- data.frame(cohort = nm, gene = g, n = ncol(y), logFC_per_stage = f$coefficients[id, "fib"],
        se = se[[id]], t = f$t[id, "fib"], p = f$p.value[id, "fib"])
    } else pc_rows[[length(pc_rows) + 1]] <- data.frame(cohort = nm, gene = g, n = ncol(y), logFC_per_stage = NA, se = NA, t = NA, p = NA)
  }
  kept_calls <- intersect(L$gid, rownames(y))
  cov_rows[[nm]] <- data.frame(cohort = nm, n = ncol(y), genes_kept = nrow(y), calls_kept = length(kept_calls),
    stable = sum(grp[kept_calls] == "stable"), dep213621 = sum(grp[kept_calls] == "dep213621"), other = sum(grp[kept_calls] == "other"))
  # label-independent voom weights for the permutation null (design without fibrosis)
  d0 <- model.matrix(~ sex, data = data.frame(sex = factor(x$sex)))
  v0 <- voom(y[kept_calls, ], d0, lib.size = y$samples$lib.size * y$samples$norm.factors)
  fits[[nm]] <- list(E = v0$E, W = v0$weights, fib = x$fib, sex = factor(x$sex), genes = kept_calls)
}
pc <- do.call(rbind, pc_rows)
meta <- do.call(rbind, lapply(names(pc_ids), function(g) {
  s <- pc[pc$gene == g & !is.na(pc$se), ]
  r <- rma(yi = s$logFC_per_stage, sei = s$se, method = "REML")
  r3 <- rma(yi = s$logFC_per_stage[s$cohort != "EMTAB9815"], sei = s$se[s$cohort != "EMTAB9815"], method = "REML")
  data.frame(gene = g, k = r$k, est = as.numeric(r$b), ci_lb = r$ci.lb, ci_ub = r$ci.ub, p_two_sided = r$pval, I2 = r$I2,
             est_noEMTAB = as.numeric(r3$b), p_noEMTAB = r3$pval)
}))
fwrite(pc, file.path(out, "positive_control_per_cohort.tsv"), sep = "\t")
fwrite(meta, file.path(out, "positive_control_meta.tsv"), sep = "\t")
fwrite(do.call(rbind, cov_rows), file.path(out, "call_coverage_filterByExpr.tsv"), sep = "\t")
print(pc); print(meta); print(do.call(rbind, cov_rows))

# gene sets
in_all4 <- Reduce(intersect, lapply(fits, `[[`, "genes"))
in_3 <- Reduce(intersect, lapply(fits[c("GSE281797", "GSE225740", "GSE267145")], `[[`, "genes"))
sets <- list(all4 = list(genes = in_all4, cohorts = names(fits)),
             three_noEMTAB = list(genes = in_3, cohorts = c("GSE281797", "GSE225740", "GSE267145")))
for (s in names(sets)) cat(s, "genes", length(sets[[s]]$genes), "by group:",
  paste(names(table(grp[sets[[s]]$genes])), table(grp[sets[[s]]$genes]), sep = ":", collapse = " "), "\n")

# vectorised weighted least squares of each gene on [1, fib, sex01] with gene-specific voom weights
wslope <- function(E, W, f, s01) {
  one <- rep(1, length(f))
  a11 <- as.vector(W %*% one); a12 <- as.vector(W %*% f); a13 <- as.vector(W %*% s01)
  a22 <- as.vector(W %*% (f^2)); a23 <- as.vector(W %*% (f * s01)); a33 <- as.vector(W %*% s01)
  WY <- W * E
  c1 <- as.vector(WY %*% one); c2 <- as.vector(WY %*% f); c3 <- as.vector(WY %*% s01)
  det <- a11 * (a22 * a33 - a23^2) - a12 * (a12 * a33 - a23 * a13) + a13 * (a12 * a23 - a22 * a13)
  i11 <- (a22 * a33 - a23^2) / det; i12 <- -(a12 * a33 - a13 * a23) / det; i13 <- (a12 * a23 - a13 * a22) / det
  i22 <- (a11 * a33 - a13^2) / det; i23 <- -(a11 * a23 - a12 * a13) / det; i33 <- (a11 * a22 - a12^2) / det
  b1 <- i11 * c1 + i12 * c2 + i13 * c3; b2 <- i12 * c1 + i22 * c2 + i23 * c3; b3 <- i13 * c1 + i23 * c2 + i33 * c3
  res <- E - (b1 + outer(b2, f) + outer(b3, s01))
  s2 <- rowSums(W * res^2) / (length(f) - 3)
  list(b = setNames(b2, rownames(E)), se = setNames(sqrt(s2 * i22), rownames(E)))
}
# check the vectorised fit against lm.wfit on the first 3 genes of each cohort
for (nm in names(fits)) {
  z <- fits[[nm]]; s01 <- as.numeric(z$sex == levels(z$sex)[2]); ws <- wslope(z$E, z$W, z$fib, s01)
  X <- cbind(1, z$fib, s01)
  for (i in 1:3) { ft <- lm.wfit(X, z$E[i, ], z$W[i, ]); stopifnot(abs(ft$coefficients[2] - ws$b[i]) < 1e-8)
    s2 <- sum(z$W[i, ] * ft$residuals^2) / ft$df.residual; R <- chol2inv(ft$qr$qr[1:3, 1:3])
    stopifnot(abs(sqrt(s2 * R[2, 2]) - ws$se[i]) < 1e-8) }
  fits[[nm]]$s01 <- s01
}
cat("vectorised WLS matches lm.wfit\n")

res <- list()
t0 <- Sys.time()
for (bb in seq_len(B)) {
  perm <- lapply(fits, function(z) sample(z$fib))
  per <- lapply(names(fits), function(nm) { z <- fits[[nm]]; wslope(z$E, z$W, perm[[nm]], z$s01) }); names(per) <- names(fits)
  for (s in names(sets)) {
    gs <- sets[[s]]$genes; cs <- sets[[s]]$cohorts
    num <- rowSums(sapply(cs, function(c) per[[c]]$b[gs] / per[[c]]$se[gs]^2))
    ss <- sign(num) == refsign[gs]
    gg <- grp[gs]
    res[[length(res) + 1]] <- data.frame(b = bb, set = s, rate_all = mean(ss), rate_stable = mean(ss[gg == "stable"]),
      rate_dep = mean(ss[gg == "dep213621"]), diff = mean(ss[gg == "stable"]) - mean(ss[gg == "dep213621"]))
  }
  if (bb %% 200 == 0) cat("perm", bb, "elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
res <- do.call(rbind, res)
fwrite(res, file.path(out, "null_same_sign_permutations.tsv"), sep = "\t")
de <- do.call(rbind, lapply(names(sets), function(s) {
  r <- res[res$set == s, ]; gs <- sets[[s]]$genes; gg <- grp[gs]
  n1 <- sum(gg == "stable"); n2 <- sum(gg == "dep213621"); na <- length(gs)
  data.frame(set = s, B = nrow(r), n_all = na, n_stable = n1, n_dep = n2,
    mean_rate_all = mean(r$rate_all), sd_rate_all = sd(r$rate_all), binom_sd_all = sqrt(0.25 / na),
    DE_all = var(r$rate_all) / (0.25 / na),
    mean_diff = mean(r$diff), sd_diff = sd(r$diff), binom_sd_diff = sqrt(0.25 * (1 / n1 + 1 / n2)),
    DE_diff = var(r$diff) / (0.25 * (1 / n1 + 1 / n2)))
}))
fwrite(de, file.path(out, "design_effect.tsv"), sep = "\t")
print(de)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
