# RV2: replication of the 1,854 adopted bulk calls in GSE281797. Follows PRESPEC_RV2.md.
# Unit: participant (one biopsy each).
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
NPERM <- if (length(args) > 1) as.integer(args[2]) else 1000L
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
DISC <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
CNT <- file.path(ROOT, "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz")
PHE <- file.path(ROOT, "data/GSE281797/metadata/sample_phenotypes_resolved.tsv")
SEED <- 20261102
wr <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")

disc <- fread(DISC)
disc[, gene_id := sub("\\.[0-9]+$", "", gene)]
disc[, call := padj < 0.05 & abs(logFC) > 0.5]
stopifnot(sum(disc$call) == 1854, sum(disc$call & disc$logFC > 0) == 1333, !anyDuplicated(disc$gene_id))

x <- fread(CNT)
counts <- as.matrix(x[, -1]); rownames(counts) <- x$gene_id
phe <- fread(PHE)
stopifnot(all(phe$matrix_column %in% colnames(counts)))
phe[, disease := as.integer(diagnosis %in% c("MASL", "MASH"))]
phe[, sex_used := sex_expressed]

prep <- function(p) {
  y <- DGEList(counts[, p$matrix_column, drop = FALSE])
  keep <- filterByExpr(y, group = factor(p$disease))
  y <- y[keep, , keep.lib.sizes = FALSE]
  calcNormFactors(y, method = "RLE")
}

fit_de <- function(y, p, term = "disease", covars = "sex_used") {
  X <- model.matrix(as.formula(paste("~", paste(c(covars, term), collapse = " + "))), p)
  v <- voomWithQualityWeights(y, X)
  f <- eBayes(lmFit(v, X))
  tt <- topTable(f, coef = term, number = Inf, sort.by = "none")
  data.table(gene_id = rownames(tt), r_logFC = tt$logFC, r_t = tt$t, r_p = tt$P.Value, r_ave = tt$AveExpr)
}

stats <- function(r, calls) {
  m <- merge(calls[, .(gene_id, d = logFC)], r[, .(gene_id, r = r_logFC)], by = "gene_id")
  c(n = nrow(m), sign_agree = mean(sign(m$r) == sign(m$d)), slope = unname(coef(lm(r ~ d, m))[2]))
}

calls <- disc[call == TRUE]

# ---------------- primary ----------------
p91 <- phe[sex_concordant == TRUE]
y91 <- prep(p91)
r91 <- fit_de(y91, p91)
gt <- merge(disc[, .(gene_id, symbol, d_logFC = logFC, d_padj = padj, call)], r91, by = "gene_id", all.x = TRUE)
gt[, testable := !is.na(r_logFC)]
obs <- stats(r91, calls)

set.seed(SEED)
null <- rbindlist(lapply(seq_len(NPERM), function(b) {
  pp <- copy(p91)
  pp[, disease := sample(disease), by = sex_used]
  as.list(stats(fit_de(y91, pp), calls))
}))
null[, perm := seq_len(.N)]
wr(null, "perm_null.tsv.gz")

prim <- data.table(statistic = c("sign_agree", "slope"), n_testable = obs[["n"]],
                   observed = c(obs[["sign_agree"]], obs[["slope"]]),
                   null_median = c(median(null$sign_agree), median(null$slope)),
                   null_q95 = c(quantile(null$sign_agree, 0.95), quantile(null$slope, 0.95)),
                   perm_p = c((1 + sum(null$sign_agree >= obs[["sign_agree"]])) / (NPERM + 1),
                              (1 + sum(null$slope >= obs[["slope"]])) / (NPERM + 1)))
prim[, holm_p := p.adjust(perm_p, "holm")]
prim[statistic == "sign_agree", naive_binom_p := binom.test(round(obs[["sign_agree"]] * obs[["n"]]), obs[["n"]], 0.5, alternative = "greater")$p.value]
wr(prim, "primary.tsv")

# per-gene one-sided p uses the moderated t with the fit's total df
dfr <- {X <- model.matrix(~ sex_used + disease, p91); v <- voomWithQualityWeights(y91, X); f <- eBayes(lmFit(v, X)); f$df.total[1]}
gt[call == TRUE & testable == TRUE, p_onesided := pt(fifelse(d_logFC > 0, r_t, -r_t), df = dfr, lower.tail = FALSE)]
gt[call == TRUE & testable == TRUE, q_onesided := p.adjust(p_onesided, "BH")]
wr(gt, "gene_table.tsv.gz")

# ---------------- descriptive ----------------
desc <- rbindlist(lapply(c("up", "down"), function(dir) {
  s <- stats(r91, calls[if (dir == "up") logFC > 0 else logFC < 0])
  data.table(set = dir, n = s[["n"]], sign_agree = s[["sign_agree"]], slope = s[["slope"]],
             n_rep_q05 = gt[call == TRUE & testable == TRUE & (d_logFC > 0) == (dir == "up"), sum(q_onesided < 0.05)])
}))
allg <- gt[testable == TRUE]
rho_all <- cor(allg$d_logFC, allg$r_logFC, method = "spearman")

# ---------------- sensitivities ----------------
sens <- list()
add <- function(name, p, ...) {
  y <- prep(p); r <- fit_de(y, p, ...); s <- stats(r, calls)
  sens[[name]] <<- data.table(analysis = name, n_participants = nrow(p), n_controls = sum(p$disease == 0),
                              n_testable = s[["n"]], sign_agree = s[["sign_agree"]], slope = s[["slope"]])
}
add("primary_n91", p91)
p94 <- copy(phe); p94[, sex_used := sex_deposited]
add("all94_deposited_sex", p94)
add("controls_F0_only", p91[disease == 1 | fibrosis == 0])
add("plus_age_bmi_diabetes", p91, covars = c("sex_used", "age_at_biopsy", "bmi", "diabetes"))
add("MASH_vs_control", p91[diagnosis %in% c("MASH", "No Pathology")])
sens <- rbindlist(sens)
set.seed(SEED + 1)
null_nosex <- sapply(seq_len(min(NPERM, 200L)), function(b) {
  pp <- copy(p91); pp[, disease := sample(disease)]; stats(fit_de(y91, pp), calls)[["sign_agree"]]
})
sens <- rbind(sens, data.table(analysis = "perm_ignoring_sex_200", n_participants = nrow(p91),
                               n_controls = sum(p91$disease == 0), n_testable = obs[["n"]],
                               sign_agree = (1 + sum(null_nosex >= obs[["sign_agree"]])) / (length(null_nosex) + 1),
                               slope = NA_real_))
wr(sens, "sensitivity.tsv")

sink(file.path(OUT, "summary.txt"))
cat(sprintf("Participants %d (controls %d; F1 controls %d; controls by sex: %s)\n", nrow(p91), sum(p91$disease == 0),
            p91[disease == 0, sum(fibrosis >= 1)], paste(p91[disease == 0, .N, by = sex_used][, paste(sex_used, N)], collapse = ", ")))
cat(sprintf("Genes after filter %d; discovery calls testable %d of 1854\n\n", nrow(y91), obs[["n"]]))
cat("Primary (Holm over 2)\n"); print(prim)
cat("\nBy direction\n"); print(desc)
cat(sprintf("\nSpearman r vs d over %d testable genes: %.3f\n", nrow(allg), rho_all))
cat(sprintf("Calls replicating at one-sided BH q < 0.05: %d of %d\n", gt[call == TRUE & testable == TRUE, sum(q_onesided < 0.05)], obs[["n"]]))
cat("\nSensitivities (last row: permutation p ignoring sex, in the sign_agree column)\n"); print(sens)
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done\n")
