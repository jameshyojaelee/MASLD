#!/usr/bin/env Rscript
# H3 step 4: carrier signature for N5. Person-level limma-voom, rule in PRESPEC_H3.md ("Carrier signature").
# Two models: PNPLA3 G dosage (0 to 2) and HSD17B13 TA carrier (0/1). Design: cohort x stage cells + sex + dosage.
# Usage: 04_carrier_signature.R CALLS_DIR SCORES_DIR OUT_DIR
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
calls_dir <- args[[1]]; scores_dir <- args[[2]]; out <- args[[3]]
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
SEED <- 20260923L
NPERM <- 10L
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"

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

# persons and design for one model; cohorts need at least 5 persons with and 5 without the tested allele
build <- function(col) {
  d <- dat[!is.na(get(col)) & !is.na(stage) & !is.na(sex)]
  keep <- d[, .(w = sum(get(col) > 0), wo = sum(get(col) == 0)), by = cohort][w >= 5 & wo >= 5, cohort]
  d <- d[cohort %in% keep]
  d[, st := merge_stage(stage), by = cohort]
  d[, cs := paste(cohort, st, sep = ":")]
  d[, sexF := as.numeric(sex == "F")]
  d[, dose := get(col)]
  d
}

run_model <- function(d, dose = d$dose, return_fit = FALSE) {
  y <- DGEList(cm[, d$person_id, drop = FALSE])
  design <- model.matrix(~ factor(cs) + sexF + dose, data = cbind(d[, .(cs, sexF)], dose = dose))
  colnames(design)[ncol(design)] <- "dosage"
  keepg <- filterByExpr(y, design)
  y <- calcNormFactors(y[keepg, , keep.lib.sizes = FALSE])
  v <- voomWithQualityWeights(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "dosage", number = Inf, sort.by = "none", adjust.method = "BH")
  tt$gene_id <- rownames(tt)
  if (return_fit) list(tt = tt, keepg = keepg) else tt
}

# per-cohort fit of the same dosage coefficient (own voom), for the same-sign count
cohort_logfc <- function(d, genes) {
  res <- list()
  for (co in unique(d$cohort)) {
    dc <- d[cohort == co]
    sf <- factor(dc$st)
    design <- cbind(model.matrix(~ sf), sexF = dc$sexF, dosage = dc$dose)
    if (qr(design)$rank < ncol(design)) next
    y <- DGEList(cm[genes, dc$person_id, drop = FALSE])
    y <- calcNormFactors(y)
    v <- voom(y, design, plot = FALSE)
    res[[co]] <- eBayes(lmFit(v, design))$coefficients[, "dosage"]
  }
  do.call(cbind, res)
}

models <- list(PNPLA3_G_dosage = "pnpla3", HSD17B13_TA_carrier = "hsd")
summ <- list(); sig_all <- list(); perm_rows <- list()
for (nm in names(models)) {
  d <- build(models[[nm]])
  message(nm, ": persons ", nrow(d), "; cohorts ", paste(unique(d$cohort), collapse = ","),
          "; with allele ", sum(d$dose > 0), "; without ", sum(d$dose == 0))
  r <- run_model(d, return_fit = TRUE)
  tt <- r$tt
  lf <- cohort_logfc(d, tt$gene_id)
  tt$n_cohorts_tested <- ncol(lf)
  tt$n_cohorts_same_sign <- rowSums(sign(lf) == sign(tt$logFC))
  tt <- merge(tt, sym, by = "gene_id", all.x = TRUE, sort = FALSE)
  setDT(tt); setorder(tt, P.Value)
  fwrite(tt, file.path(out, paste0("all_genes_", nm, ".tsv.gz")), sep = "\t")
  sig <- tt[adj.P.Val < 0.05]
  sig[, model := nm]
  sig_all[[nm]] <- sig
  fwrite(sig, file.path(out, paste0("signature_fdr05_", nm, ".tsv")), sep = "\t")
  # permutation null: genotype shuffled within cohort
  set.seed(SEED)
  cnt <- integer(NPERM); minp <- numeric(NPERM)
  for (b in seq_len(NPERM)) {
    dp <- copy(d)
    dp[, dose := dose[sample.int(.N)], by = cohort]
    tp <- run_model(dp)
    cnt[b] <- sum(tp$adj.P.Val < 0.05); minp[b] <- min(tp$P.Value)
    message(nm, " permutation ", b, ": genes at FDR < 0.05 = ", cnt[b])
  }
  perm_rows[[nm]] <- data.table(model = nm, perm = seq_len(NPERM), genes_fdr05 = cnt, min_p = minp)
  summ[[nm]] <- data.table(model = nm, persons = nrow(d), cohorts = paste(unique(d$cohort), collapse = ","),
                           n_with_allele = sum(d$dose > 0), n_without_allele = sum(d$dose == 0),
                           genes_tested = nrow(tt), genes_fdr05 = nrow(sig),
                           genes_fdr05_same_sign_in_all_tested_cohorts = sig[n_cohorts_same_sign == n_cohorts_tested, .N],
                           genes_fdr05_same_sign_in_at_least_half = sig[n_cohorts_same_sign >= ceiling(n_cohorts_tested / 2), .N],
                           smallest_p = min(tt$P.Value),
                           perm_genes_fdr05_min = min(cnt), perm_genes_fdr05_median = median(cnt), perm_genes_fdr05_max = max(cnt),
                           lambda_median_t2 = median(tt$t^2) / qchisq(0.5, 1))
  print(summ[[nm]])
}
fwrite(rbindlist(summ), file.path(out, "signature_summary.tsv"), sep = "\t")
fwrite(rbindlist(perm_rows), file.path(out, "signature_permutations.tsv"), sep = "\t")
fwrite(rbindlist(sig_all, fill = TRUE), file.path(out, "carrier_signature_for_N5.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("done")
