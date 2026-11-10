# U1 independent refit of the 774-unit pooled disease-versus-control model.
# Written from docs/manuscript/METHODS.md ("Disease-versus-control differential
# expression"), not from the producer. Labels come from the U0 recount
# (my_774_labels.tsv); sex is called here from XIST and Y-linked expression.
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(data.table)
})
set.seed(20261023)
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]; U0 <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
COUNTS <- file.path(R, "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds")
UNTEST <- file.path(R, "figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8/untestable_annotations.tsv")
PROD <- file.path(R, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
logf <- file(file.path(OUT, "u1_log.txt"), "w")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, logf); flush(logf) }
sha <- function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), " +")[[1]][1]
fwrite(data.table(path = c(COUNTS, UNTEST, file.path(PROD, "deg_results.csv"), file.path(U0, "my_774_labels.tsv")),
                  sha256 = sapply(c(COUNTS, UNTEST, file.path(PROD, "deg_results.csv"), file.path(U0, "my_774_labels.tsv")), sha)),
       file.path(OUT, "input_sha256.tsv"), sep = "\t")

counts <- readRDS(COUNTS)
say("counts: ", nrow(counts), " genes x ", ncol(counts), " libraries")
lab <- fread(file.path(U0, "my_774_labels.tsv"))
stopifnot(setequal(lab$sample_id, colnames(counts)), !anyDuplicated(lab$sample_id))
lab <- lab[match(colnames(counts), lab$sample_id)]

# ---- 855 held annotations
ut <- fread(UNTEST)
say("held untestable annotations: ", nrow(ut), "; present in count matrix: ", sum(ut$gene_id %in% rownames(counts)),
    "; present after version strip: ", sum(sub("\\..*", "", ut$gene_id) %in% sub("\\..*", "", rownames(counts))))
say("version-stripped duplicate IDs in counts: ", sum(duplicated(sub("\\..*", "", rownames(counts)))))

# ---- expression sex from counts (independent of any sex field)
base <- sub("\\..*", "", rownames(counts))
lib <- colSums(counts)
xist <- counts[base == "ENSG00000229807", ]
ygenes <- c("ENSG00000129824", "ENSG00000067048", "ENSG00000012817", "ENSG00000183878", "ENSG00000198692", "ENSG00000114374")
ysum <- colSums(counts[base %in% ygenes, , drop = FALSE])
lx <- log2(xist / lib * 1e6 + 1); ly <- log2(ysum / lib * 1e6 + 1)
expr_sex <- ifelse(ly > lx, "M", "F")
margin <- abs(ly - lx)
say("expression sex: ", paste(names(table(expr_sex)), table(expr_sex), collapse = " "), "; min |log2CPM(Y)-log2CPM(XIST)| = ", round(min(margin), 2))
prod_design <- fread(file.path(PROD, "model_design.tsv"))
prod_sex <- ifelse(prod_design$inferred_sexM[match(colnames(counts), prod_design$sample_id)] == 1, "M", "F")
dis <- which(expr_sex != prod_sex)
say("expression sex vs producer design sex disagreements: ", length(dis))
if (length(dis)) print(data.frame(sample = colnames(counts)[dis], dataset = lab$dataset[dis], expr = expr_sex[dis], prod = prod_sex[dis],
                                  reported = lab$sex_reported[dis], lx = round(lx[dis], 2), ly = round(ly[dis], 2)))
say("expression sex vs reported sex disagreements: ", sum(expr_sex != lab$sex_reported, na.rm = TRUE),
    " (reported missing: ", sum(is.na(lab$sex_reported)), ")")
fwrite(data.table(sample_id = colnames(counts), dataset = lab$dataset, expr_sex, prod_sex, reported = lab$sex_reported,
                  log2cpm_xist = lx, log2cpm_y = ly), file.path(OUT, "sex_calls.tsv"), sep = "\t")

# ---- primary fit (Methods)
group <- factor(lab$group_binary, levels = c("Control", "Disease"))
dataset <- factor(lab$dataset)
say("units: ", ncol(counts), " control ", sum(group == "Control"), " disease ", sum(group == "Disease"))
print(table(dataset, group))
fit_one <- function(sex, tag, keep.lib.sizes = FALSE, filter_design = model.matrix(~ 0 + group)) {
  y <- DGEList(counts, group = group)
  keep <- filterByExpr(y, design = filter_design)
  y <- y[keep, , keep.lib.sizes = keep.lib.sizes]
  y <- calcNormFactors(y, method = "RLE")
  design <- model.matrix(~ dataset + sex + group)
  v <- voomWithQualityWeights(y, design)
  fit <- lmFit(v, design)
  eb <- eBayes(fit)
  co <- "groupDisease"
  tt <- topTable(eb, coef = co, number = Inf, sort.by = "none", adjust.method = "BH")
  tr <- topTreat(treat(fit, lfc = 0.25), coef = co, number = Inf, sort.by = "none", adjust.method = "BH")
  se <- eb$stdev.unscaled[, co] * sqrt(eb$s2.post)
  canon <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5
  say(sprintf("[%s] genes=%d canonical=%d up=%d down=%d treat=%d treat_up=%d treat_down=%d df.prior=%.3f",
              tag, nrow(tt), sum(canon), sum(canon & tt$logFC > 0), sum(canon & tt$logFC < 0),
              sum(tr$adj.P.Val < 0.05), sum(tr$adj.P.Val < 0.05 & tr$logFC > 0), sum(tr$adj.P.Val < 0.05 & tr$logFC < 0), eb$df.prior[1]))
  list(y = y, v = v, tt = data.table(gene = rownames(tt), logFC = tt$logFC, SE = se, P = tt$P.Value, padj = tt$adj.P.Val,
                                     treat_p = tr$P.Value, treat_fdr = tr$adj.P.Val, AveExpr = tt$AveExpr, canon = canon), keep = keep)
}
sex_expr <- factor(expr_sex, levels = c("F", "M"))
A <- fit_one(sex_expr, "primary: my labels, expression sex, keep.lib.sizes=FALSE")

# ---- compare with producer deg_results.csv
prod <- fread(file.path(PROD, "deg_results.csv"))
prod[, canon := padj < 0.05 & abs(logFC) > 0.5]
say("producer: genes=", nrow(prod), " canonical=", sum(prod$canon), " up=", sum(prod$canon & prod$logFC > 0),
    " down=", sum(prod$canon & prod$logFC < 0), " treat=", sum(prod$treat_fdr < 0.05))
compare <- function(mine, tag) {
  m <- merge(mine, prod, by = "gene", suffixes = c("_mine", "_prod"))
  say(sprintf("[%s vs producer] gene universe: mine %d, producer %d, shared %d, mine-only %d, producer-only %d",
              tag, nrow(mine), nrow(prod), nrow(m), length(setdiff(mine$gene, prod$gene)), length(setdiff(prod$gene, mine$gene))))
  say(sprintf("  canonical membership differences: mine-only %d, producer-only %d; treat differences: mine-only %d producer-only %d",
              sum(m$canon_mine & !m$canon_prod), sum(!m$canon_mine & m$canon_prod),
              sum(m$treat_fdr_mine < 0.05 & !(m$treat_fdr_prod < 0.05)), sum(!(m$treat_fdr_mine < 0.05) & m$treat_fdr_prod < 0.05)))
  say(sprintf("  max|dlogFC|=%.3g  max|dSE|=%.3g  max|dlog10 padj|=%.3g  max|d treat_fdr|=%.3g  max|dAveExpr|=%.3g",
              max(abs(m$logFC_mine - m$logFC_prod)), max(abs(m$SE_mine - m$SE_prod)),
              max(abs(log10(m$padj_mine) - log10(m$padj_prod))), max(abs(m$treat_fdr_mine - m$treat_fdr_prod)),
              max(abs(m$AveExpr_mine - m$AveExpr_prod))))
  invisible(m)
}
mA <- compare(A$tt, "primary")
fwrite(A$tt, file.path(OUT, "my_refit_774_deg.tsv"), sep = "\t")

# ---- sensitivity of the fit to recipe choices (to locate any mismatch)
sex_prod <- factor(prod_sex, levels = c("F", "M"))
if (length(dis)) { B <- fit_one(sex_prod, "producer design sex"); compare(B$tt, "producer sex") }
C <- fit_one(sex_expr, "keep.lib.sizes=TRUE", keep.lib.sizes = TRUE); compare(C$tt, "keep.lib.sizes=TRUE")
yy <- DGEList(counts, group = group)
say("filterByExpr variants: ~0+group design ", sum(filterByExpr(yy, design = model.matrix(~ 0 + group))),
    "; group= ", sum(filterByExpr(yy, group = group)),
    "; full design ", sum(filterByExpr(yy, design = model.matrix(~ dataset + sex_expr + group))),
    "; no group ", sum(filterByExpr(yy)))

# ---- PCA (Fig S3A): voom log-CPM, regress cohort + sex, top 2000 variable genes
E <- A$v$E
Ec <- removeBatchEffect(E, batch = dataset, batch2 = sex_expr)
for (variant in c("var_after_regress", "var_before_regress")) {
  vv <- if (variant == "var_after_regress") apply(Ec, 1, var) else apply(E, 1, var)
  top <- order(vv, decreasing = TRUE)[1:2000]
  pc <- prcomp(t(Ec[top, ]), center = TRUE, scale. = FALSE)
  pv <- pc$sdev^2 / sum(pc$sdev^2)
  say(sprintf("PCA %s: PC1 %.2f%% PC2 %.2f%%", variant, 100 * pv[1], 100 * pv[2]))
}
# plain voom (no quality weights) E is identical in voomWithQualityWeights; also try logCPM from cpm()
lc <- cpm(A$y, log = TRUE, prior.count = 2)
lcc <- removeBatchEffect(lc, batch = dataset, batch2 = sex_expr)
top <- order(apply(lcc, 1, var), decreasing = TRUE)[1:2000]
pc <- prcomp(t(lcc[top, ])); pv <- pc$sdev^2 / sum(pc$sdev^2)
say(sprintf("PCA edgeR cpm(log) prior 2: PC1 %.2f%% PC2 %.2f%%", 100 * pv[1], 100 * pv[2]))

# ---- endpoint contrast (Fig S3B)
ep <- fread(file.path(U0, "my_endpoint_labels.tsv"))
say("endpoint units: ", nrow(ep), " advanced ", sum(ep$endpoint == "Advanced"), " strict control ", sum(ep$endpoint == "StrictControl"))
j <- match(ep$sample_id, colnames(counts))
for (nf in c("pooled_RLE_factors", "recomputed_RLE_on_subset")) {
  ys <- A$y[, j]
  if (nf == "recomputed_RLE_on_subset") ys <- calcNormFactors(ys, method = "RLE")
  eg <- factor(ep$endpoint, levels = c("StrictControl", "Advanced"))
  ds <- factor(ep$dataset); sx <- sex_expr[j]
  d2 <- model.matrix(~ ds + sx + eg)
  v2 <- voomWithQualityWeights(ys, d2)
  e2 <- eBayes(lmFit(v2, d2))
  t2 <- topTable(e2, coef = "egAdvanced", number = Inf, sort.by = "none")
  c2 <- t2$adj.P.Val < 0.05 & abs(t2$logFC) > 0.5
  rho <- cor(t2$logFC, A$tt$logFC, method = "spearman")
  rec <- sum(c2 & A$tt$canon & sign(t2$logFC) == sign(A$tt$logFC))
  opp <- sum(c2 & A$tt$canon & sign(t2$logFC) != sign(A$tt$logFC))
  say(sprintf("[endpoint %s] genes=%d rho=%.4f endpoint calls=%d recovered=%d of %d; endpoint-only=%d; opposite=%d",
              nf, nrow(t2), rho, sum(c2), rec, sum(A$tt$canon), sum(c2 & !A$tt$canon), opp))
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("done")
close(logf)
