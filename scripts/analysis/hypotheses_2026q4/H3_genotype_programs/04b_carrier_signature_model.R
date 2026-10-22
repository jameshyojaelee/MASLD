#!/usr/bin/env Rscript
# H3 step 4b: carrier signature for N5, one model per job. Same model as 04_carrier_signature.R
# (PRESPEC_H3.md, "Carrier signature" and Amendment 4). The base fit is voomWithQualityWeights as prespecified.
# Permutation fits reuse the observed voom precision weights, because refitting the weights took longer than 30 minutes per fit.
# Usage: 04b_carrier_signature_model.R CALLS_DIR SCORES_DIR OUT_DIR MODEL   (MODEL = pnpla3 or hsd)
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
calls_dir <- args[[1]]; scores_dir <- args[[2]]; out <- args[[3]]; model <- args[[4]]
stopifnot(model %in% c("pnpla3", "hsd"), !dir.exists(out))
dir.create(out, recursive = TRUE)
SEED <- 20260923L
NPERM <- 10L
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
nm <- c(pnpla3 = "PNPLA3_G_dosage", hsd = "HSD17B13_TA_carrier")[[model]]

pc <- readRDS(file.path(scores_dir, "person_counts.rds"))
ph <- fread(file.path(scores_dir, "person_phenotypes.tsv"))
gt <- fread(file.path(calls_dir, "person_genotypes.tsv"))
ann <- fread(file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))
sym <- ann[!duplicated(ensembl_base), .(gene_id = ensembl_base, symbol = gene_name)]

cm <- pc$counts
stopifnot(identical(colnames(cm), ph$person_id))
dat <- merge(ph[, .(person_id, cohort, stage = as.numeric(stage), sex = fifelse(sex %in% c("F", "M"), sex, NA_character_))],
             gt[, .(person_id, rs738409_call, rs72613567_call)], by = "person_id")
dat[, pnpla3 := unname(c(RR = 0, RA = 1, AA = 2)[rs738409_call])]
dat[, hsd := ifelse(rs72613567_call == "CARRIER", 1, ifelse(rs72613567_call == "NC", 0, NA_real_))]

merge_stage <- function(stage, min_n = 3L) {
  s <- as.numeric(stage)
  repeat {
    tb <- table(s)
    if (length(tb) < 2L || all(tb >= min_n)) break
    lv <- as.numeric(names(tb)); i <- which.min(tb)
    others <- lv[-i]; dd <- abs(others - lv[i])
    s[s == lv[i]] <- others[which(dd == min(dd))[1]]
  }
  s
}

d <- dat[!is.na(get(model)) & !is.na(stage) & !is.na(sex)]
keep <- d[, .(w = sum(get(model) > 0), wo = sum(get(model) == 0)), by = cohort][w >= 5 & wo >= 5, cohort]
d <- d[cohort %in% keep]
d[, st := merge_stage(stage), by = cohort]
d[, cs := paste(cohort, st, sep = ":")]
d[, sexF := as.numeric(sex == "F")]
d[, dose := get(model)]
message(nm, ": persons ", nrow(d), "; cohorts ", paste(unique(d$cohort), collapse = ","),
        "; with allele ", sum(d$dose > 0), "; without ", sum(d$dose == 0))

mk_design <- function(dose) {
  X <- model.matrix(~ factor(cs) + sexF + dose, data = cbind(d[, .(cs, sexF)], dose = dose))
  colnames(X)[ncol(X)] <- "dosage"
  X
}
design <- mk_design(d$dose)
y <- DGEList(cm[, d$person_id, drop = FALSE])
keepg <- filterByExpr(y, design)
y <- calcNormFactors(y[keepg, , keep.lib.sizes = FALSE])
t0 <- Sys.time()
v <- voomWithQualityWeights(y, design, plot = FALSE)
message("voomWithQualityWeights: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
fit <- eBayes(lmFit(v, design))
tt <- topTable(fit, coef = "dosage", number = Inf, sort.by = "none", adjust.method = "BH")
tt$gene_id <- rownames(tt)

# per-cohort fit of the same coefficient (own voom) for the same-sign count
res <- list()
for (co in unique(d$cohort)) {
  dc <- d[cohort == co]
  X <- cbind(model.matrix(~ factor(st), data = dc), sexF = dc$sexF, dosage = dc$dose)
  if (qr(X)$rank < ncol(X)) next
  yc <- calcNormFactors(DGEList(cm[tt$gene_id, dc$person_id, drop = FALSE]))
  res[[co]] <- eBayes(lmFit(voom(yc, X, plot = FALSE), X))$coefficients[, "dosage"]
}
lf <- do.call(cbind, res)
tt$n_cohorts_tested <- ncol(lf)
tt$n_cohorts_same_sign <- rowSums(sign(lf) == sign(tt$logFC))
tt <- merge(tt, sym, by = "gene_id", all.x = TRUE, sort = FALSE)
setDT(tt); setorder(tt, P.Value)
fwrite(tt, file.path(out, paste0("all_genes_", nm, ".tsv.gz")), sep = "\t")
sig <- tt[adj.P.Val < 0.05]
sig[, model := nm]
fwrite(sig, file.path(out, paste0("signature_fdr05_", nm, ".tsv")), sep = "\t")
message("genes at FDR < 0.05: ", nrow(sig))

# permutation null: genotype shuffled within cohort; observed voom weights held fixed
set.seed(SEED)
cnt <- integer(NPERM); minp <- numeric(NPERM)
for (b in seq_len(NPERM)) {
  dp <- copy(d)
  dp[, dose := dose[sample.int(.N)], by = cohort]
  tp <- topTable(eBayes(lmFit(v, mk_design(dp$dose))), coef = "dosage", number = Inf, sort.by = "none", adjust.method = "BH")
  cnt[b] <- sum(tp$adj.P.Val < 0.05); minp[b] <- min(tp$P.Value)
  message("permutation ", b, ": genes at FDR < 0.05 = ", cnt[b], "; smallest p = ", signif(minp[b], 3))
}
fwrite(data.table(model = nm, perm = seq_len(NPERM), genes_fdr05 = cnt, min_p = minp), file.path(out, paste0("permutations_", nm, ".tsv")), sep = "\t")
summ <- data.table(model = nm, persons = nrow(d), cohorts = paste(unique(d$cohort), collapse = ","),
                   n_with_allele = sum(d$dose > 0), n_without_allele = sum(d$dose == 0),
                   genes_tested = nrow(tt), genes_fdr05 = nrow(sig),
                   genes_fdr05_same_sign_in_all_tested_cohorts = sig[n_cohorts_same_sign == n_cohorts_tested, .N],
                   genes_fdr05_same_sign_in_at_least_half = sig[n_cohorts_same_sign >= ceiling(n_cohorts_tested / 2), .N],
                   smallest_p = min(tt$P.Value),
                   perm_genes_fdr05_min = min(cnt), perm_genes_fdr05_median = median(cnt), perm_genes_fdr05_max = max(cnt),
                   lambda_median_t2 = median(tt$t^2) / qchisq(0.5, 1))
fwrite(summ, file.path(out, paste0("summary_", nm, ".tsv")), sep = "\t")
print(summ)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("done")
