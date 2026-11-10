# RV1 step 2: Part A (per-hepatocyte GNMT, MAT1A, CYP2C19) and Part B2 (program coherence).
# Follows PRESPEC_RV1.md. Unit: donor.
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
MEMB <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv")
SEED <- 20261031
FOCAL_GENES <- c("GNMT", "MAT1A", "CYP2C19")
IDENTITY <- c("ALB", "APOA1", "TF", "SERPINA1", "HP")
ECM <- "hotspot_hepatocytes_f05c535ae5bbc0b9"
DUCT <- "hotspot_hepatocytes_48f39dd4d817a10e"
STAGE <- c(Healthy = 0, Steatosis = 1, Steatohepatitis = 2, Cirrhosis = 3)
wr <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")

donors <- fread(file.path(OUT, "donors.tsv"))
read_counts <- function(l) {
  x <- fread(file.path(OUT, sprintf("counts_%s.tsv.gz", l)))
  m <- as.matrix(x[, -1]); rownames(m) <- x$gene; m
}
CNT <- list(Hepatocytes = read_counts("Hepatocytes"), Fibroblasts = read_counts("Fibroblasts"),
            Cholangiocytes = read_counts("Cholangiocytes"))

# Donor table for one cohort and lineage with covariates.
cohort_meta <- function(cohort, lineage, min_nuclei, drop_cirrhosis, labels = "corrected") {
  d <- donors[dataset == cohort]
  d[, diagnosis := if (labels == "corrected") diagnosis_corrected else diagnosis_bf3]
  d <- d[diagnosis %in% names(STAGE)]
  d[, nuclei := get(paste0("nuclei_", lineage))]
  d[, umi := get(paste0("umi_", lineage))]
  d <- d[nuclei >= min_nuclei]
  if (drop_cirrhosis) d <- d[diagnosis != "Cirrhosis"]
  d[, stage := STAGE[diagnosis]]
  d[, male := as.integer(sex == "Male")]
  d[, log_nuc := log10(nuclei)]
  d[, log_upn := log10(umi / nuclei)]
  d
}

logcpm <- function(lineage, d) {
  y <- DGEList(CNT[[lineage]][, d$donor, drop = FALSE])
  keep <- filterByExpr(y, group = factor(d$diagnosis))
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  cpm(y, log = TRUE, prior.count = 1)
}

# OLS for all genes with HC3 SE for one coefficient.
ols_hc3 <- function(Y, X, coef) {
  XtXi <- solve(crossprod(X))
  M <- XtXi %*% t(X)
  B <- Y %*% t(M)
  E <- Y - B %*% t(X)
  h <- rowSums((X %*% XtXi) * X)
  E2 <- (E / matrix(1 - h, nrow(E), ncol(E), byrow = TRUE))^2
  j <- which(colnames(X) == coef)
  se <- sqrt(as.vector(E2 %*% (M[j, ]^2)))
  df <- nrow(X) - ncol(X)
  t <- B[, j] / se
  data.table(gene = rownames(Y), beta = B[, j], se = se, t = t,
             p = 2 * pt(-abs(t), df), df = df, mean_logcpm = rowMeans(Y))
}

fit_cohort <- function(cohort, drop_cirrhosis = TRUE, use_upn = TRUE, term = "stage", labels = "corrected") {
  d <- cohort_meta(cohort, "Hepatocytes", 50, drop_cirrhosis, labels)
  if (term == "disease") d[, disease := as.integer(diagnosis != "Healthy")]
  Y <- logcpm("Hepatocytes", d)
  f <- as.formula(paste("~", term, "+ male + log_nuc", if (use_upn) "+ log_upn" else ""))
  X <- model.matrix(f, d)
  r <- ols_hc3(Y, X, term)
  r[, decile := cut(rank(mean_logcpm), 10, labels = FALSE)]
  r[, pct_in_decile := frank(t) / .N, by = decile]
  list(res = r, d = d, Y = Y, X = X)
}

# ---------------- Part A ----------------
A <- list(); CTRL <- list(); SENS <- list(); PERM <- list()
set.seed(SEED)
for (co in c("GSE202379", "GSE244832")) {
  fc <- fit_cohort(co)
  r <- fc$res
  n <- nrow(fc$d)
  tab <- r[gene %in% FOCAL_GENES]
  missing <- setdiff(FOCAL_GENES, tab$gene)
  if (length(missing)) tab <- rbind(tab, data.table(gene = missing), fill = TRUE)
  tab[, `:=`(cohort = co, n_donors = n, stage_counts = paste(names(table(fc$d$diagnosis)), table(fc$d$diagnosis), sep = "=", collapse = ";"))]
  A[[co]] <- tab
  CTRL[[co]] <- r[gene %in% IDENTITY][, cohort := co]
  # permutation p for focal genes, covariates fixed
  Yf <- fc$Y[intersect(FOCAL_GENES, rownames(fc$Y)), , drop = FALSE]
  obs <- ols_hc3(Yf, fc$X, "stage")$t
  cnt <- numeric(length(obs))
  for (b in 1:10000) {
    Xp <- fc$X; Xp[, "stage"] <- sample(Xp[, "stage"])
    cnt <- cnt + (abs(ols_hc3(Yf, Xp, "stage")$t) >= abs(obs))
  }
  PERM[[co]] <- data.table(cohort = co, gene = rownames(Yf), perm_p = (1 + cnt) / 10001)
  for (s in list(list("disease_vs_healthy", TRUE, TRUE, "disease", "corrected"), list("no_umi_per_nucleus", TRUE, FALSE, "stage", "corrected"),
                 list("bf3_labels", TRUE, TRUE, "stage", "bf3"))) {
    fs <- fit_cohort(co, s[[2]], s[[3]], s[[4]], s[[5]])
    SENS[[paste(co, s[[1]])]] <- fs$res[gene %in% FOCAL_GENES][, `:=`(cohort = co, analysis = s[[1]], n_donors = nrow(fs$d))]
  }
}
fs <- fit_cohort("GSE202379", drop_cirrhosis = FALSE)
SENS[["cirr"]] <- fs$res[gene %in% FOCAL_GENES][, `:=`(cohort = "GSE202379", analysis = "with_cirrhosis", n_donors = nrow(fs$d))]

A <- rbindlist(A, fill = TRUE)
A <- merge(A, rbindlist(PERM), by = c("cohort", "gene"), all.x = TRUE)
A[, holm_p := p.adjust(p, "holm")]
A[, call := fifelse(!is.na(holm_p) & holm_p < 0.05 & beta < 0, "lower_per_hepatocyte", "not_shown")]
wr(A, "A_gene_tests.tsv")
wr(rbindlist(CTRL), "A_controls.tsv")
wr(rbindlist(SENS, fill = TRUE), "A_sensitivity.tsv")

# ---------------- Part B2 ----------------
memb <- fread(MEMB)[cell_type == "hepatocytes"]
programs <- unique(memb[, .(program_uid, module)])

coherence <- function(R, genes) {
  g <- intersect(genes, rownames(R))
  if (length(g) < 5) return(c(NA_real_, length(g)))
  C <- cor(t(R[g, ]), method = "spearman")
  c(mean(C[upper.tri(C)]), length(g))
}

run_b2 <- function(cohort, lineage, min_nuclei, uids) {
  d <- cohort_meta(cohort, lineage, min_nuclei, FALSE)
  if (nrow(d) < 8) return(NULL)
  Y <- logcpm(lineage, d)
  f <- if (length(unique(d$diagnosis)) > 1) ~ factor(diagnosis) + male + log_nuc + log_upn else ~ male + log_nuc + log_upn
  X <- model.matrix(f, d)
  R <- Y - (Y %*% X %*% solve(crossprod(X))) %*% t(X)
  dec <- cut(rank(rowMeans(Y)), 10, labels = FALSE); names(dec) <- rownames(Y)
  rbindlist(lapply(uids, function(u) {
    genes <- memb[program_uid == u, mapped_symbol]
    obs <- coherence(R, genes)
    g <- intersect(genes, rownames(R))
    if (is.na(obs[1])) return(data.table(cohort = cohort, lineage = lineage, program_uid = u, n_donors = nrow(d), n_genes = obs[2]))
    pool <- setdiff(rownames(R), genes)
    by_dec <- split(pool, dec[pool])
    null <- replicate(1000, {
      s <- vapply(g, function(x) { cand <- by_dec[[as.character(dec[x])]]; cand[sample.int(length(cand), 1)] }, "")
      coherence(R, s)[1]
    })
    data.table(cohort = cohort, lineage = lineage, program_uid = u, n_donors = nrow(d), n_genes = obs[2],
               coherence = obs[1], null_median = median(null), null_q95 = quantile(null, 0.95),
               p = (1 + sum(null >= obs[1])) / 1001)
  }))
}

set.seed(SEED)
B2 <- rbindlist(list(
  run_b2("GSE202379", "Hepatocytes", 50, programs$program_uid),
  run_b2("GSE244832", "Hepatocytes", 50, programs$program_uid),
  run_b2("GSE202379", "Fibroblasts", 20, ECM),
  run_b2("GSE244832", "Fibroblasts", 20, ECM),
  run_b2("GSE202379", "Cholangiocytes", 20, DUCT),
  run_b2("GSE244832", "Cholangiocytes", 20, DUCT)), fill = TRUE)
B2 <- merge(B2, programs, by = "program_uid", all.x = TRUE)
B2[, role := fifelse(program_uid == ECM, "Stromal ECM", fifelse(program_uid == DUCT, "Ductular injury", "other"))]
B2[, holm_p_primary := NA_real_]
B2[cohort == "GSE202379" & lineage == "Hepatocytes" & role != "other" & !is.na(p), holm_p_primary := p.adjust(p, "holm")]
wr(B2, "B2_coherence.tsv")

other <- B2[cohort == "GSE202379" & lineage == "Hepatocytes" & role == "other" & !is.na(p)]
sink(file.path(OUT, "summary.txt"))
cat("Part A primary (Holm over 6)\n"); print(A[, .(cohort, gene, n_donors, beta, se, p, holm_p, perm_p, pct_in_decile, call)])
cat("\nIdentity controls\n"); print(rbindlist(CTRL)[, .(cohort, gene, beta, p)])
cat("\nPart B2 primary (GSE202379 hepatocytes, Holm over 2)\n"); print(B2[cohort == "GSE202379" & lineage == "Hepatocytes" & role != "other", .(role, n_donors, n_genes, coherence, null_q95, p, holm_p_primary)])
cat(sprintf("\nCalibration: other GSE202379 hepatocyte programs with p < 0.05: %d of %d\n", sum(other$p < 0.05), nrow(other)))
oth2 <- B2[cohort == "GSE244832" & lineage == "Hepatocytes" & !is.na(p)]
cat(sprintf("GSE244832 hepatocyte programs with p < 0.05: %d of %d\n", sum(oth2$p < 0.05), nrow(oth2)))
cat("\nSecondary lineages\n"); print(B2[lineage != "Hepatocytes", .(cohort, lineage, role, n_donors, n_genes, coherence, null_q95, p)])
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done\n")
