# Light positive-control check (login node): weighted least squares on log-CPM for the
# Stromal ECM genes among the 1,854 calls, using the stored voom quality sample weights
# (no voom observation-level weights). Validated against the stored unadjusted logFC.
suppressMessages({library(data.table); library(edgeR)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
NEW <- file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz")
OLD <- file.path(ROOT, "Analysis/Deconvolution/results/%s/%s_bayesprism_proportions.tsv")
MEMB <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv")
DS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

d <- readRDS(file.path(P, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(P, "meta_matched.rds")))[match(colnames(d), sample_id)]
w <- fread(file.path(P, "lvqw_sample_weights.tsv"))[match(colnames(d), sample_id), sample_weight]
deg <- fread(file.path(P, "deg_results.csv"))
memb <- fread(MEMB)[program_uid == "hotspot_hepatocytes_f05c535ae5bbc0b9"]
calls <- deg[padj < 0.05 & abs(logFC) > 0.5]
pcg <- calls[symbol %in% unique(memb$source_gene)]
cat("Stromal ECM source genes", uniqueN(memb$source_gene), "; in fit", deg[symbol %in% memb$source_gene, .N],
    "; in calls", nrow(pcg), "\n")
# Only the positive-control genes are fitted; other calls are not touched (primary stays blind).
oth <- calls[0]
genes <- pcg$gene
lc <- cpm(d, log = TRUE, prior.count = 0.5)[genes, , drop = FALSE]  # voom uses lib.size*norm.factors and 0.5 offset

info <- data.frame(group = factor(m$group_binary, levels = c("Control", "Disease")), dataset = factor(m$dataset), sex = factor(m$inferred_sex))
X0 <- model.matrix(~ dataset + sex + group, info)
ct <- setdiff(names(fread(NEW, nrows = 1)), c("sample_id", "dataset"))
newf <- as.matrix(fread(NEW)[match(colnames(d), sample_id), ..ct])
oldf <- as.matrix(do.call(rbind, lapply(DS, function(g) read.delim(sprintf(OLD, g, g), check.names = FALSE)))[colnames(d), ct])
H <- contr.helmert(16); H <- sweep(H, 2, sqrt(colSums(H^2)), "/")
ilr <- function(f) log(pmax(f, 1e-6)) %*% H
coef_g <- function(X) { b <- lm.wfit(X, t(lc), w)$coefficients; b["groupDisease", ] }
b0 <- coef_g(X0); bn <- coef_g(cbind(X0, ilr(newf))); bo <- coef_g(cbind(X0, ilr(oldf)))
res <- data.table(gene = genes, set = rep(c("ECM", "other"), c(nrow(pcg), nrow(oth))),
                  symbol = c(pcg$symbol, oth$symbol), stored = c(pcg$logFC, oth$logFC),
                  wls_unadj = b0, wls_rebuilt = bn, wls_earlier = bo)
res[, `:=`(ratio_rebuilt = wls_rebuilt / wls_unadj, ratio_earlier = wls_earlier / wls_unadj)]
cat("approximation check: Pearson(stored, WLS unadjusted) =", round(cor(res$stored, res$wls_unadj), 4),
    "; median |diff| =", round(median(abs(res$stored - res$wls_unadj)), 4), "\n")
print(res[set == "ECM", .(symbol, stored = round(stored, 3), wls_unadj = round(wls_unadj, 3), rebuilt = round(wls_rebuilt, 3),
                          earlier = round(wls_earlier, 3), ratio_rebuilt = round(ratio_rebuilt, 3), ratio_earlier = round(ratio_earlier, 3))])
cat("ECM: n", nrow(pcg), "| ratio<1 rebuilt", res[set == "ECM", sum(ratio_rebuilt < 1)], "earlier", res[set == "ECM", sum(ratio_earlier < 1)],
    "| rebuilt more attenuated than earlier", res[set == "ECM", sum(ratio_rebuilt < ratio_earlier)], "\n")
cat("median ratio (rebuilt, earlier): ECM", round(res[set == "ECM", median(ratio_rebuilt)], 3), round(res[set == "ECM", median(ratio_earlier)], 3), "\n")
# fibroblast shift by group within dataset, each input
for (nm in c("rebuilt", "earlier")) {
  f <- if (nm == "rebuilt") newf else oldf
  y <- log(pmax(f[, "Fibroblasts"], 1e-6))
  b <- lm.wfit(X0, y, w)$coefficients["groupDisease"]
  cat("disease minus control, log fibroblast fraction (dataset+sex adjusted),", nm, ":", round(b, 3), "\n")
}
fwrite(res, "pc_light_results.tsv", sep = "\t")
