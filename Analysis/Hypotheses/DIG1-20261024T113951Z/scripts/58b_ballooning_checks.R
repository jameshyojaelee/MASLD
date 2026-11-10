# DIG1 job 58b (exploratory, post hoc to job 58). Robustness of the GSE130970 ballooning link and the
# cross-cohort fibrosis link of the inflamed-hepatocyte bulk score. Descriptive; no new test family.
# Scores are rebuilt exactly as in job 58 (within-cohort z of log2 TMM CPM, prior 2, mean over genes).
#   A. GSE130970 (n 71): score ~ steatosis + lobular + ballooning + fibrosis + IMMUNE + HEPMASS for
#      each INFLAMED_HEP gene alone, a hepatocyte-only score (FCAMR, TRIM31; Flex hepatocyte share >= 0.97),
#      and the full score after rank-based inverse normal transform. Leave-one-sample-out range of the
#      ballooning coefficient for the full score.
#   B. Per cohort fibrosis model + IMMUNE + HEPMASS + DUCTULAR (is the stage link the ductular reaction?).
# HC3 errors computed directly, as in job 58.
# Usage: Rscript 58b_ballooning_checks.R <new_out_dir>
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1)
out <- args[1]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
d <- readRDS(file.path(rel, "merged_dge.rds")); m <- readRDS(file.path(rel, "meta_matched.rds"))
m <- m[match(colnames(d$counts), as.character(m$sample_id)), ]
stopifnot(identical(colnames(d$counts), as.character(m$sample_id)))
deg <- read.csv(file.path(rel, "deg_results.csv")); sym <- setNames(deg$symbol, deg$gene)
rn <- rownames(d$counts); g_sym <- ifelse(rn %in% names(sym), sym[rn], rn)
d <- calcNormFactors(d); lc <- cpm(d, log = TRUE, prior.count = 2)
INF <- c("FCAMR", "TRIM31", "CXCL10", "IL32", "CCL20")
SETS <- c(list(INFLAMED_HEP = INF, HEPONLY = c("FCAMR", "TRIM31"),
               IMMUNE = c("PTPRC", "CD3E", "CD2", "CD68", "CD163"), DUCTULAR = c("KRT7", "KRT19", "EPCAM", "SOX9"),
               HEPMASS = c("ALB", "APOA1", "CYP2E1", "HP")), setNames(as.list(INF), paste0("g_", INF)))
sc <- sapply(SETS, function(gs) {
  idx <- which(g_sym %in% gs); v <- rep(NA_real_, ncol(lc))
  for (co in unique(m$dataset)) { j <- which(m$dataset == co); v[j] <- colMeans(t(scale(t(lc[idx, j, drop = FALSE])))) }
  v
})
D <- cbind(m[, c("sample_id", "dataset", "fibrosis_stage", "nas_score", "steatosis_grade", "lobular_grade", "ballooning_grade")], as.data.frame(sc))
fit <- function(f, dat) {
  mod <- lm(f, data = dat); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  b <- coef(mod); se <- sqrt(diag(V))
  data.frame(term = names(b), estimate = b, se = se, p = 2 * pt(-abs(b / se), df.residual(mod)), n = nobs(mod), row.names = NULL)
}
rint <- function(x) qnorm((rank(x) - 0.5) / length(x))
g7 <- subset(D, dataset == "GSE130970"); g7$INFLAMED_RINT <- rint(g7$INFLAMED_HEP)
comp <- "~ steatosis_grade + lobular_grade + ballooning_grade + fibrosis_stage"
rows <- list()
for (s in c("INFLAMED_HEP", "INFLAMED_RINT", "HEPONLY", paste0("g_", INF), "IMMUNE")) for (adj in c("", " + IMMUNE + HEPMASS")) {
  if (s == "IMMUNE" && adj != "") next
  r <- fit(as.formula(paste(s, comp, adj)), g7); r$outcome <- s; r$adjusted <- adj != ""; rows[[length(rows) + 1]] <- r
}
A <- do.call(rbind, rows); A <- subset(A, term != "(Intercept)")
write.table(A, file.path(out, "gse130970_component_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
loo <- sapply(seq_len(nrow(g7)), function(i) {
  r <- fit(as.formula(paste("INFLAMED_HEP", comp, "+ IMMUNE + HEPMASS")), g7[-i, ]); r$estimate[r$term == "ballooning_grade"]
})
mb <- aggregate(cbind(INFLAMED_HEP, HEPONLY, IMMUNE) ~ ballooning_grade, g7, function(x) round(c(mean = mean(x), n = length(x)), 3))
rows <- list()
for (co in c("GSE130970", "GSE135251", "GSE162694", "GSE213621")) {
  dd <- subset(D, dataset == co & !is.na(fibrosis_stage)); has_nas <- sum(!is.na(dd$nas_score)) > 20
  if (has_nas) dd <- subset(dd, !is.na(nas_score))
  for (s in c("INFLAMED_HEP", "HEPONLY")) {
    f <- paste(s, "~ fibrosis_stage", if (has_nas) "+ nas_score" else "", "+ IMMUNE + HEPMASS + DUCTULAR")
    r <- fit(as.formula(f), dd); r$cohort <- co; r$outcome <- s; rows[[length(rows) + 1]] <- r
  }
}
B <- do.call(rbind, rows); B <- subset(B, term %in% c("fibrosis_stage", "nas_score", "DUCTULAR"))
write.table(B, file.path(out, "cohort_fibrosis_ductular_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200)
print(subset(A, term %in% c("ballooning_grade", "lobular_grade", "steatosis_grade", "fibrosis_stage")), digits = 3)
cat("LOO ballooning coefficient (adjusted, full score): min", round(min(loo), 3), "max", round(max(loo), 3), "\n")
print(mb); print(B, digits = 3)
