# DIG1 job 58 (exploratory, post hoc to job 57). Which histological feature does the inflamed-hepatocyte
# state track in bulk biopsies?
# Data: 774-unit bulk release (bulk-corrected-controls-20261008T162308Z), merged DGEList and matched metadata.
# Scores per sample: within-cohort z of log2 CPM (edgeR cpm, TMM, prior 2), mean over set genes.
#   INFLAMED_HEP = FCAMR, TRIM31, UBD, CXCL10, IL32, CCL20 (Flex MASLD hepatocyte share >= 0.75 and top-5%
#                  concentration > 2x null in Flex; job 57). Fixed before this job.
#   Comparators: IMMUNE (PTPRC, CD3E, CD2, CD68, CD163), FIBRO (COL1A1, COL1A2, COL3A1, LUM),
#                DUCTULAR (KRT7, KRT19, EPCAM, SOX9), HEPMASS (ALB, APOA1, CYP2E1, HP).
# Primary family (Holm over 4): GSE130970 (n 71, all four NAS components and fibrosis available)
#   INFLAMED_HEP ~ steatosis + lobular + ballooning + fibrosis, OLS with HC3 errors.
# Descriptive: the same for each comparator; + IMMUNE + HEPMASS as covariates (sensitivity);
#   per cohort score ~ NAS + fibrosis (GSE130970, GSE135251, GSE162694) and ~ fibrosis (GSE213621),
#   fixed-effect inverse-variance summary of the fibrosis and NAS coefficients across cohorts.
# Unit: one sample per analysis unit (the release units). Cross-sectional.
# Usage: Rscript 58_bulk_state_histology.R <new_out_dir>
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1)
out <- args[1]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
set.seed(20261029)
rel <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
d <- readRDS(file.path(rel, "merged_dge.rds")); m <- readRDS(file.path(rel, "meta_matched.rds"))
stopifnot(setequal(colnames(d$counts), as.character(m$sample_id)), !anyDuplicated(m$sample_id))
m <- m[match(colnames(d$counts), as.character(m$sample_id)), ]
stopifnot(identical(colnames(d$counts), as.character(m$sample_id)), identical(as.character(d$samples$dataset), as.character(m$dataset)))
deg <- read.csv(file.path(rel, "deg_results.csv"))
sym <- setNames(deg$symbol, deg$gene)
rn <- rownames(d$counts); g_sym <- ifelse(rn %in% names(sym), sym[rn], rn)
SETS <- list(INFLAMED_HEP = c("FCAMR", "TRIM31", "UBD", "CXCL10", "IL32", "CCL20"),
             IMMUNE = c("PTPRC", "CD3E", "CD2", "CD68", "CD163"), FIBRO = c("COL1A1", "COL1A2", "COL3A1", "LUM"),
             DUCTULAR = c("KRT7", "KRT19", "EPCAM", "SOX9"), HEPMASS = c("ALB", "APOA1", "CYP2E1", "HP"))
d <- calcNormFactors(d)
lc <- cpm(d, log = TRUE, prior.count = 2)
sc <- matrix(NA_real_, ncol(lc), length(SETS), dimnames = list(colnames(lc), names(SETS)))
found <- list()
for (s in names(SETS)) {
  idx <- which(g_sym %in% SETS[[s]]); found[[s]] <- g_sym[idx]
  for (co in unique(m$dataset)) {
    j <- which(m$dataset == co)
    z <- t(scale(t(lc[idx, j, drop = FALSE])))
    sc[j, s] <- colMeans(z, na.rm = TRUE)
  }
}
writeLines(sapply(names(found), function(s) paste(s, paste(found[[s]], collapse = ","))), file.path(out, "set_genes_found.txt"))
D <- cbind(m[, c("sample_id", "dataset", "group_binary", "fibrosis_stage", "nas_score", "steatosis_grade", "lobular_grade", "ballooning_grade")],
           as.data.frame(sc))
write.table(D, file.path(out, "sample_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fit <- function(f, dat) {
  # OLS with HC3 errors computed directly (sandwich is not installed in this env); t on residual df
  mod <- lm(f, data = dat); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  b <- coef(mod); se <- sqrt(diag(V))
  data.frame(term = names(b), estimate = b, se = se, p = 2 * pt(-abs(b / se), df.residual(mod)), n = nobs(mod), row.names = NULL)
}
std <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
rows <- list()
g7 <- subset(D, dataset == "GSE130970")
for (s in names(SETS)) {
  for (adj in c("none", "immune+hepmass")) {
    if (s %in% c("IMMUNE", "HEPMASS") && adj != "none") next
    f <- as.formula(paste(s, "~ steatosis_grade + lobular_grade + ballooning_grade + fibrosis_stage",
                          if (adj == "none") "" else "+ IMMUNE + HEPMASS"))
    r <- fit(f, g7); r$set <- s; r$model <- paste0("GSE130970_components_", adj); rows[[length(rows) + 1]] <- r
  }
}
for (co in c("GSE130970", "GSE135251", "GSE162694", "GSE213621")) {
  dd <- subset(D, dataset == co & !is.na(fibrosis_stage))
  has_nas <- sum(!is.na(dd$nas_score)) > 20
  if (has_nas) dd <- subset(dd, !is.na(nas_score))
  for (s in names(SETS)) {
    f <- as.formula(paste(s, "~ fibrosis_stage", if (has_nas) "+ nas_score" else ""))
    r <- fit(f, dd); r$set <- s; r$model <- paste0(co, "_nas_fibrosis"); rows[[length(rows) + 1]] <- r
    if (s == "INFLAMED_HEP") {
      f2 <- as.formula(paste(s, "~ fibrosis_stage", if (has_nas) "+ nas_score" else "", "+ IMMUNE + HEPMASS"))
      r <- fit(f2, dd); r$set <- s; r$model <- paste0(co, "_nas_fibrosis_adj"); rows[[length(rows) + 1]] <- r
    }
  }
}
R <- do.call(rbind, rows)
write.table(R, file.path(out, "model_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
prim <- subset(R, model == "GSE130970_components_none" & set == "INFLAMED_HEP" & term != "(Intercept)")
prim$p_holm <- p.adjust(prim$p, "holm")
write.table(prim, file.path(out, "primary_holm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
meta <- list()
for (s in names(SETS)) for (tm in c("fibrosis_stage", "nas_score")) for (adj in c("", "_adj")) {
  if (adj == "_adj" && s != "INFLAMED_HEP") next
  x <- subset(R, set == s & term == tm & grepl(paste0("_nas_fibrosis", adj, "$"), model))
  if (nrow(x) < 2) next
  w <- 1 / x$se^2; est <- sum(w * x$estimate) / sum(w); se <- sqrt(1 / sum(w))
  Q <- sum(w * (x$estimate - est)^2)
  meta[[length(meta) + 1]] <- data.frame(set = s, term = tm, adjusted = adj == "_adj", k = nrow(x), estimate = est, se = se,
                                         p = 2 * pnorm(-abs(est / se)), Q = Q, I2 = max(0, (Q - (nrow(x) - 1)) / Q))
}
M <- do.call(rbind, meta)
write.table(M, file.path(out, "fixed_effect_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cr <- cor(sc[, names(SETS)], method = "spearman", use = "pairwise.complete.obs")
write.table(round(cr, 3), file.path(out, "score_correlations.tsv"), sep = "\t", quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200)
print(found); print(prim); print(subset(R, term != "(Intercept)" & grepl("GSE130970_components", model)), digits = 3)
print(M, digits = 3); print(round(cr, 2))
