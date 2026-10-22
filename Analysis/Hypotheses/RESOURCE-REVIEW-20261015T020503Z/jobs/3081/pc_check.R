# Gate L03-04 positive-control check. Reproduces the canonical disease-vs-control fit,
# then refits with 15 ILR coordinates from each fraction input. Reports only the
# Stromal ECM program genes among the canonical calls, plus reproduction diagnostics.
# The primary statistic (median ratio over all calls) is deliberately not computed.
suppressMessages({library(data.table); library(edgeR); library(limma)})
set.seed(20261015)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
NEW <- file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz")
OLD <- file.path(ROOT, "Analysis/Deconvolution/results/%s/%s_bayesprism_proportions.tsv")
MEMB <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv")
DS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

d <- readRDS(file.path(P, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(P, "meta_matched.rds")))
m <- m[match(colnames(d), sample_id)]
stopifnot(identical(m$sample_id, colnames(d)))
deg <- fread(file.path(P, "deg_results.csv"))

info <- data.frame(group_binary = factor(m$group_binary, levels = c("Control", "Disease")),
                   dataset = factor(m$dataset), inferred_sex = factor(m$inferred_sex))
base <- model.matrix(~ dataset + inferred_sex + group_binary, data = info)

newf <- fread(NEW); ct <- setdiff(names(newf), c("sample_id", "dataset"))
newf <- as.matrix(newf[match(colnames(d), sample_id), ..ct])
oldl <- lapply(DS, function(g) { x <- read.delim(sprintf(OLD, g, g), check.names = FALSE); x })
oldf <- as.matrix(do.call(rbind, oldl)[colnames(d), ct])
stopifnot(!anyNA(newf), !anyNA(oldf))
H <- contr.helmert(16); H <- sweep(H, 2, sqrt(colSums(H^2)), "/")  # orthonormal Helmert basis
ilr <- function(f) { z <- log(pmax(f, 1e-6)) %*% H; colnames(z) <- sprintf("ilr_%02d", 1:15); z }

fitit <- function(design, label) {
  t0 <- Sys.time()
  v <- voomWithQualityWeights(d, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  ci <- match("group_binaryDisease", colnames(fit$coefficients))
  se <- fit$stdev.unscaled[, ci] * sqrt(fit$s2.post)
  cat(label, "fit seconds", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
  data.table(gene = rownames(tt), logFC = tt$logFC, SE = se, padj = tt$adj.P.Val)
}

u <- fitit(base, "unadjusted")
mm <- merge(u, deg[, .(gene, logFC0 = logFC, SE0 = SE, padj0 = padj, symbol)], by = "gene")
cat("reproduction: genes", nrow(mm), "max |logFC - stored|", signif(max(abs(mm$logFC - mm$logFC0)), 3),
    "max |SE - stored|", signif(max(abs(mm$SE - mm$SE0)), 3),
    "calls refit", mm[padj < 0.05 & abs(logFC) > 0.5, .N], "calls stored", mm[padj0 < 0.05 & abs(logFC0) > 0.5, .N], "\n")

a_new <- fitit(cbind(base, ilr(newf)), "adjusted_rebuilt")
a_old <- fitit(cbind(base, ilr(oldf)), "adjusted_earlier")

memb <- fread(MEMB)[program_uid == "hotspot_hepatocytes_f05c535ae5bbc0b9"]
ecm_sym <- unique(memb$source_gene)
calls <- mm[padj0 < 0.05 & abs(logFC0) > 0.5]
pc <- calls[symbol %in% ecm_sym, .(gene, symbol, logFC_unadj = logFC0, SE_unadj = SE0)]
pc <- merge(pc, a_new[, .(gene, logFC_rebuilt = logFC, SE_rebuilt = SE, padj_rebuilt = padj)], by = "gene")
pc <- merge(pc, a_old[, .(gene, logFC_earlier = logFC, SE_earlier = SE, padj_earlier = padj)], by = "gene")
pc[, ratio_rebuilt := logFC_rebuilt / logFC_unadj]
pc[, ratio_earlier := logFC_earlier / logFC_unadj]
pc[, call_rebuilt := padj_rebuilt < 0.05 & abs(logFC_rebuilt) > 0.5]
pc[, call_earlier := padj_earlier < 0.05 & abs(logFC_earlier) > 0.5]
setorder(pc, symbol)
fwrite(pc, "positive_control_stromal_ecm.tsv", sep = "\t")
print(pc[, .(symbol, logFC_unadj = round(logFC_unadj, 3), logFC_rebuilt = round(logFC_rebuilt, 3), logFC_earlier = round(logFC_earlier, 3),
             ratio_rebuilt = round(ratio_rebuilt, 3), ratio_earlier = round(ratio_earlier, 3), call_rebuilt, call_earlier)])
cat("ECM calls n", nrow(pc), "median ratio rebuilt", round(median(pc$ratio_rebuilt), 3), "earlier", round(median(pc$ratio_earlier), 3),
    "n ratio<1 rebuilt", sum(pc$ratio_rebuilt < 1), "earlier", sum(pc$ratio_earlier < 1),
    "n rebuilt more attenuated than earlier", sum(pc$ratio_rebuilt < pc$ratio_earlier), "\n")
# design-level SE inflation actually realised, median over all genes (not a test statistic)
se_inf <- data.table(rebuilt = median(a_new$SE / u$SE), earlier = median(a_old$SE / u$SE))
cat("median SE inflation all genes: rebuilt", round(se_inf$rebuilt, 3), "earlier", round(se_inf$earlier, 3), "\n")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
