#!/usr/bin/env Rscript
# H3 step 4d: HSD17B13 carrier signature, sensitivity to read depth at the site (PRESPEC_H3.md, Amendment 5).
# Base design as in 04b (cohort x stage cells + sex + carrier); adjusted design adds log2(site depth + 1).
# Also repeats the 10 prespecified fixed-weight shuffles on the base design to report the inflation of each permuted fit.
# Usage: 04d_hsd_depth_sensitivity.R CALLS_DIR SCORES_DIR GENO_DIR BASE_SIG_DIR OUT_DIR
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
calls_dir <- args[[1]]; scores_dir <- args[[2]]; geno_dir <- args[[3]]; base_dir <- args[[4]]; out <- args[[5]]
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
SEED <- 20260923L
NPERM <- 10L
nm <- "HSD17B13_TA_carrier"

pc <- readRDS(file.path(scores_dir, "person_counts.rds"))
ph <- fread(file.path(scores_dir, "person_phenotypes.tsv"))
gt <- fread(file.path(calls_dir, "person_genotypes.tsv"))
pvc <- fread(file.path(geno_dir, "person_variant_counts.tsv"))[rsid == "rs72613567", .(person_id, depth)]
stopifnot(!anyDuplicated(pvc$person_id))

cm <- pc$counts
stopifnot(identical(colnames(cm), ph$person_id))
dat <- merge(ph[, .(person_id, cohort, stage = as.numeric(stage), sex = fifelse(sex %in% c("F", "M"), sex, NA_character_))],
             gt[, .(person_id, rs72613567_call)], by = "person_id")
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

d <- dat[!is.na(hsd) & !is.na(stage) & !is.na(sex)]
keep <- d[, .(w = sum(hsd > 0), wo = sum(hsd == 0)), by = cohort][w >= 5 & wo >= 5, cohort]
d <- d[cohort %in% keep]
d[, st := merge_stage(stage), by = cohort]
d[, cs := paste(cohort, st, sep = ":")]
d[, sexF := as.numeric(sex == "F")]
d[, dose := hsd]
d <- merge(d, pvc, by = "person_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(d$depth))
d[, ldepth := log2(depth + 1)]
# the merge above can reorder rows; restore the order used by the base script
d <- d[match(dat[person_id %in% d$person_id, person_id], person_id)]
message("persons ", nrow(d), "; with carrier ", sum(d$dose > 0), "; without ", sum(d$dose == 0))
stopifnot(nrow(d) == 350L)

collin <- d[, .(r_carrier_ldepth = cor(dose, ldepth), median_depth_carrier = as.numeric(median(depth[dose == 1])), median_depth_noncarrier = as.numeric(median(depth[dose == 0]))), by = cohort]
collin <- rbind(collin, d[, .(cohort = "ALL", r_carrier_ldepth = cor(dose, ldepth), median_depth_carrier = as.numeric(median(depth[dose == 1])), median_depth_noncarrier = as.numeric(median(depth[dose == 0])))])
fwrite(collin, file.path(out, "carrier_vs_site_depth.tsv"), sep = "\t")
print(collin)

mk_design <- function(dose, adjust) {
  df <- cbind(d[, .(cs, sexF, ldepth)], dose = dose)
  X <- if (adjust) model.matrix(~ factor(cs) + sexF + ldepth + dose, data = df) else model.matrix(~ factor(cs) + sexF + dose, data = df)
  colnames(X)[ncol(X)] <- "dosage"
  X
}
lam <- function(tt) median(tt$t^2) / qchisq(0.5, 1)

X0 <- mk_design(d$dose, FALSE)
X1 <- mk_design(d$dose, TRUE)
r2 <- summary(lm(d$dose ~ X1[, setdiff(colnames(X1), c("dosage", "(Intercept)"))]))$r.squared
message("R2 of carrier on the other adjusted-design columns: ", round(r2, 3), "; VIF ", round(1 / (1 - r2), 2))

y <- DGEList(cm[, d$person_id, drop = FALSE])
keepg <- filterByExpr(y, X0)
y <- calcNormFactors(y[keepg, , keep.lib.sizes = FALSE])

t0 <- Sys.time()
v0 <- voomWithQualityWeights(y, X0, plot = FALSE)
message("base voomWithQualityWeights: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
tt0 <- topTable(eBayes(lmFit(v0, X0)), coef = "dosage", number = Inf, sort.by = "none", adjust.method = "BH")
tt0$gene_id <- rownames(tt0)

# reproducibility of the base fit against the 04b output
bs <- fread(file.path(base_dir, paste0("all_genes_", nm, ".tsv.gz")))
mm <- merge(tt0[, c("gene_id", "t")], bs[, .(gene_id, t_04b = t)], by = "gene_id")
base_match <- max(abs(mm$t - mm$t_04b))
message("max |t - t_04b| over ", nrow(mm), " genes: ", signif(base_match, 3))

t0 <- Sys.time()
v1 <- voomWithQualityWeights(y, X1, plot = FALSE)
message("adjusted voomWithQualityWeights: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
tt1 <- topTable(eBayes(lmFit(v1, X1)), coef = "dosage", number = Inf, sort.by = "none", adjust.method = "BH")
tt1$gene_id <- rownames(tt1)

ann <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
sym <- ann[!duplicated(ensembl_base), .(gene_id = ensembl_base, symbol = gene_name)]
adj <- merge(as.data.table(tt1), sym, by = "gene_id", all.x = TRUE, sort = FALSE)
setorder(adj, P.Value)
fwrite(adj, file.path(out, paste0("all_genes_", nm, "_depth_adjusted.tsv.gz")), sep = "\t")

both <- merge(as.data.table(tt0)[, .(gene_id, logFC0 = logFC, t0 = t, fdr0 = adj.P.Val)],
              as.data.table(tt1)[, .(gene_id, logFC1 = logFC, t1 = t, p1 = P.Value, fdr1 = adj.P.Val)], by = "gene_id")
sig0 <- both[fdr0 < 0.05]
cmp <- data.table(
  genes_tested = nrow(both),
  base_genes_fdr05 = sum(both$fdr0 < 0.05),
  adjusted_genes_fdr05 = sum(both$fdr1 < 0.05),
  base_signature_genes_still_fdr05_adjusted = sum(sig0$fdr1 < 0.05),
  base_signature_genes_nominal_p05_adjusted = sum(sig0$p1 < 0.05),
  base_signature_genes_same_sign_adjusted = sum(sign(sig0$logFC1) == sign(sig0$logFC0)),
  cor_t_pearson = cor(both$t0, both$t1),
  cor_t_spearman = cor(both$t0, both$t1, method = "spearman"),
  cor_logFC_base_signature = cor(sig0$logFC0, sig0$logFC1),
  median_ratio_logFC_adjusted_over_base_signature = median(sig0$logFC1 / sig0$logFC0),
  lambda_base = lam(tt0), lambda_adjusted = lam(tt1),
  smallest_p_base = min(tt0$P.Value), smallest_p_adjusted = min(tt1$P.Value),
  r2_carrier_on_other_columns = r2, vif_carrier = 1 / (1 - r2),
  max_abs_t_diff_vs_04b_base = base_match)
fwrite(cmp, file.path(out, "compare_to_base.tsv"), sep = "\t")
print(t(cmp))

# the prespecified shuffles, repeated on the base design with the observed voom weights held fixed (same code path and seed as 04b)
set.seed(SEED)
pr <- list()
for (b in seq_len(NPERM)) {
  dp <- copy(d)
  dp[, dose := dose[sample.int(.N)], by = cohort]
  tp <- topTable(eBayes(lmFit(v0, mk_design(dp$dose, FALSE))), coef = "dosage", number = Inf, sort.by = "none", adjust.method = "BH")
  pr[[b]] <- data.table(perm = b, genes_fdr05 = sum(tp$adj.P.Val < 0.05), min_p = min(tp$P.Value), lambda_median_t2 = lam(tp))
}
pr <- rbindlist(pr)
old <- fread(file.path(base_dir, paste0("permutations_", nm, ".tsv")))
pr[, genes_fdr05_04b := old$genes_fdr05]
pr[, min_p_04b := old$min_p]
fwrite(pr, file.path(out, "permutations_replicate_with_lambda.tsv"), sep = "\t")
print(pr)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("done")
