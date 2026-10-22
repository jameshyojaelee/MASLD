# Blind re-check F3cdg (2026-10-15): independent refits for R03-02 (GSE126848 depth) and R03-09
# (pooled 774-unit library-complexity sensitivity with release sample weights fixed).
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
set.seed(20260923)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- commandArgs(TRUE)[1]
A <- file.path(B, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
dge <- readRDS(file.path(A, "merged_dge.rds"))
meta <- as.data.frame(readRDS(file.path(A, "meta_matched.rds")))
meta <- meta[match(colnames(dge), meta$sample_id), ]
stopifnot(identical(meta$sample_id, colnames(dge)))
deg <- fread(file.path(A, "deg_results.csv"))
stopifnot(identical(deg$gene, rownames(dge)))
res <- list(); put <- function(k, v) { res[[k]] <<- paste(format(v, digits = 7), collapse = " | "); cat(k, ":", res[[k]], "\n") }
put("n_genes_units", c(nrow(dge), ncol(dge)))
put("libsize_eq_colsums", isTRUE(all.equal(as.numeric(dge$samples$lib.size), as.numeric(colSums(dge$counts)))))

## ---------------- R03-02 ----------------
ix <- which(meta$dataset == "GSE126848")
d <- dge[, ix]; grp <- factor(meta$group_binary[ix], levels = c("Control", "Disease"))
put("gse126848_units", table(grp))
lib <- d$samples$lib.size; leff <- lib * d$samples$norm.factors
put("median_lib_ctrl_dis", c(median(lib[grp == "Control"]), median(lib[grp == "Disease"])))
put("median_efflib_ctrl_dis", c(median(leff[grp == "Control"]), median(leff[grp == "Disease"])))
pairs_deeper <- sum(outer(lib[grp == "Disease"], lib[grp == "Control"], ">"))
put("pairs_disease_deeper_of", c(pairs_deeper, sum(grp == "Disease") * sum(grp == "Control")))
put("wilcox_exact_p", wilcox.test(lib[grp == "Disease"], lib[grp == "Control"], exact = TRUE)$p.value)
put("cor_group_log2lib", cor(as.integer(grp == "Disease"), log2(lib)))
callf <- function(tt) tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5
fitc <- function(dd, X) { v <- voomWithQualityWeights(dd, X); f <- eBayes(lmFit(v, X)); list(tt = topTable(f, coef = 2, number = Inf, sort.by = "none"), v = v) }
X0 <- model.matrix(~grp)
b <- fitc(d, X0); cb <- callf(b$tt)
put("base_calls_up_down", c(sum(cb), sum(cb & b$tt$logFC > 0), sum(cb & b$tt$logFC < 0)))
# compare with installed per-cohort result
coh <- fread(file.path(B, "figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"))[dataset == "GSE126848"]
coh <- coh[match(rownames(d), gene_id_versioned)]
put("max_abs_logFC_diff_vs_installed", max(abs(coh$logFC - b$tt$logFC)))
put("calls_identical_to_installed", identical(callf(data.frame(adj.P.Val = coh$FDR, logFC = coh$logFC)), cb))
z <- rowSums(d$counts) == 0
put("zero_read_genes_called_up", c(sum(z), sum(cb & z), sum(cb & z & b$tt$logFC > 0)))
put("zero_read_logFC_range", range(b$tt$logFC[z]))
w <- b$v$weights; lc <- log2(leff + 1)
wd <- sapply(which(z), function(i) weighted.mean(lc[grp == "Control"], w[i, grp == "Control"]) - weighted.mean(lc[grp == "Disease"], w[i, grp == "Disease"]))
put("max_abs_zero_logFC_minus_weighted_ctrl_minus_dis_log2Leff1", max(abs(b$tt$logFC[z] - wd)))
put("qw_median_ctrl_dis_percohort", c(median(b$v$targets$sample.weights[grp == "Control"]), median(b$v$targets$sample.weights[grp == "Disease"])))
sw <- fread(file.path(A, "lvqw_sample_weights.tsv"))
sw <- sw[match(colnames(d), sample_id)]
put("qw_median_ctrl_dis_pooled774", c(median(sw$sample_weight[grp == "Control"]), median(sw$sample_weight[grp == "Disease"])))
lr <- log2(lib); le <- log2(leff)
a1 <- fitc(d, model.matrix(~grp + lr)); c1 <- callf(a1$tt)
put("rawlib_calls_kept_zero_kept_spearman", c(sum(c1), sum(c1 & cb), sum(c1 & cb & z), cor(b$tt$logFC, a1$tt$logFC, method = "spearman")))
a2 <- fitc(d, model.matrix(~grp + le)); c2 <- callf(a2$tt)
put("efflib_calls_kept_zero_kept_spearman", c(sum(c2), sum(c2 & cb), sum(c2 & cb & z), cor(b$tt$logFC, a2$tt$logFC, method = "spearman")))
keep <- filterByExpr(d, group = grp)
d3 <- calcNormFactors(d[keep, , keep.lib.sizes = FALSE], method = "RLE")
a3 <- fitc(d3, X0); c3 <- callf(a3$tt)
put("filterByExpr_genes_calls_up_down", c(sum(keep), sum(c3), sum(c3 & a3$tt$logFC > 0), sum(c3 & a3$tt$logFC < 0)))
fwrite(data.table(gene = rownames(d), base_call = cb, zero_all30 = z, rawlib_call = c1, efflib_call = c2),
       file.path(out, "r0302_gse126848_genes_recheck.tsv"), sep = "\t")

## ---------------- R03-09 ----------------
pcall <- deg$padj < 0.05 & abs(deg$logFC) > 0.5
put("pooled_release_calls", sum(pcall))
des <- fread(file.path(A, "model_design.tsv")); stopifnot(identical(des$sample_id, colnames(dge)))
X <- as.matrix(des[, -1]); rownames(X) <- des$sample_id
swa <- fread(file.path(A, "lvqw_sample_weights.tsv")); swa <- swa$sample_weight[match(colnames(dge), swa$sample_id)]
cnt <- dge$counts; cs <- colSums(cnt)
zw <- function(x) ave(x, meta$dataset, FUN = function(y) (y - mean(y)) / sd(y))
p <- t(t(cnt) / cs)
ent <- -colSums(ifelse(p > 0, p * log(p), 0))
covs <- data.frame(
  loglib_raw = zw(log(dge$samples$lib.size)),
  loglib_eff = zw(log(dge$samples$lib.size * dge$samples$norm.factors)),
  top100 = zw(apply(cnt, 2, function(x) sum(sort(x, decreasing = TRUE)[1:100])) / cs),
  entropy = zw(ent),
  detected = zw(colSums(cnt > 0)))
fitp <- function(Xm) {
  v <- voom(dge, Xm, weights = swa)
  v$weights <- sweep(v$weights, 2, swa, "*")
  topTable(eBayes(lmFit(v, Xm)), coef = "group_binaryDisease", number = Inf, sort.by = "none")
}
mods <- list(baseline = NULL, loglib_raw = "loglib_raw", loglib_eff = "loglib_eff", top100 = "top100", entropy = "entropy",
             joint_raw = c("loglib_raw", "top100", "detected"), joint_eff = c("loglib_eff", "top100", "detected"))
rows <- list()
for (nm in names(mods)) {
  Xm <- if (is.null(mods[[nm]])) X else cbind(X, as.matrix(covs[, mods[[nm]], drop = FALSE]))
  tt <- fitp(Xm); cl <- callf(tt)
  rows[[nm]] <- data.table(model = nm, calls = sum(cl), up = sum(cl & tt$logFC > 0), down = sum(cl & tt$logFC < 0),
                           kept = sum(cl & pcall), kept_pct = round(100 * sum(cl & pcall) / sum(pcall), 2),
                           spearman = round(cor(tt$logFC, deg$logFC, method = "spearman"), 4),
                           max_abs_logFC_diff_vs_release = if (nm == "baseline") max(abs(tt$logFC - deg$logFC)) else NA_real_)
  print(rows[[nm]])
}
fwrite(rbindlist(rows), file.path(out, "r0309_refits_recheck.tsv"), sep = "\t")
writeLines(paste(names(res), unlist(res), sep = "\t"), file.path(out, "r0302_recheck_values.tsv"))  # fixed after job 9540 failed here; values in slurm-22256437.out
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
