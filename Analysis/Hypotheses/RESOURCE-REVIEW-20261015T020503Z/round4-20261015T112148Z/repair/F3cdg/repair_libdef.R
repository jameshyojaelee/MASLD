# F3cdg repair (2026-10-15): library-size definition sensitivity for R03-09 and R03-02.
# R03-09: pooled 774-unit fit, release sample weights fixed (voom with prior weights, then
#   weights multiplied by them, as in the first round and check job 9041 mode C).
#   Library size defined two ways: raw = total counts on the tested genes (lib.size);
#   effective = lib.size x RLE norm factor (the size voom uses for log-CPM).
# R03-02: GSE126848 30-unit model; zero-read log2FC identity and ~group + log2(effective lib).
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
set.seed(20260923)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- commandArgs(TRUE)[1]
arm <- file.path(base, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
dge <- readRDS(file.path(arm, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(arm, "meta_matched.rds")))
meta <- meta[match(colnames(dge), meta$sample_id)]
stopifnot(all(meta$sample_id == colnames(dge)))
log_line <- function(...) { s <- paste0(...); cat(s, "\n"); s }
L <- character()
L <- c(L, log_line("genes in merged_dge: ", nrow(dge), "; units: ", ncol(dge)))
cs <- colSums(dge$counts)
L <- c(L, log_line("lib.size equals colSums over these genes: ", all(abs(dge$samples$lib.size - cs) < 0.5)))

deg <- fread(file.path(arm, "deg_results.csv"))
base_call <- deg$padj < 0.05 & abs(deg$logFC) > 0.5
swr <- fread(file.path(arm, "lvqw_sample_weights.tsv"))
swr <- swr$sample_weight[match(colnames(dge), swr$sample_id)]
stopifnot(!anyNA(swr))
des <- fread(file.path(arm, "model_design.tsv"))
stopifnot(all(des$sample_id == colnames(dge)))
X0 <- as.matrix(des[, -1]); rownames(X0) <- des$sample_id

cnt <- dge$counts
zc <- function(x) ave(x, meta$dataset, FUN = function(y) (y - mean(y)) / sd(y))
lib_raw <- dge$samples$lib.size
lib_eff <- dge$samples$lib.size * dge$samples$norm.factors
top100 <- apply(cnt, 2, function(x) sum(sort(x, decreasing = TRUE)[1:100])) / cs
detected <- colSums(cnt > 0)
cov <- data.table(loglib_raw = zc(log(lib_raw)), loglib_eff = zc(log(lib_eff)),
                  top100 = zc(top100), detected = zc(detected))
L <- c(L, log_line("cor(z log raw lib, z log effective lib) = ", round(cor(cov$loglib_raw, cov$loglib_eff), 4)))

runfit <- function(Xm) {
  v <- voom(dge, Xm, weights = swr)
  v$weights <- t(swr * t(v$weights))
  f <- eBayes(lmFit(v, Xm))
  topTable(f, coef = "group_binaryDisease", number = Inf, sort.by = "none")
}
mods <- list(baseline = character(), loglib_raw = "loglib_raw", joint_raw = c("loglib_raw", "top100", "detected"),
             loglib_eff = "loglib_eff", joint_eff = c("loglib_eff", "top100", "detected"))
rows <- list()
for (nm in names(mods)) {
  Xm <- X0
  if (length(mods[[nm]])) Xm <- cbind(Xm, as.matrix(cov[, mods[[nm]], with = FALSE]))
  t0 <- Sys.time(); tt <- runfit(Xm)
  stopifnot(all(rownames(tt) == deg$gene))
  cl <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5
  r <- data.table(model = nm, covariates = paste(mods[[nm]], collapse = "+"), calls = sum(cl),
                  up = sum(cl & tt$logFC > 0), down = sum(cl & tt$logFC < 0), kept_of_1854 = sum(cl & base_call),
                  kept_fraction = round(sum(cl & base_call) / sum(base_call), 4),
                  up_kept = sum(cl & base_call & deg$logFC > 0), down_kept = sum(cl & base_call & deg$logFC < 0),
                  spearman_vs_release = round(cor(tt$logFC, deg$logFC, method = "spearman"), 4),
                  fit_seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  print(r); rows[[nm]] <- r
}
fwrite(rbindlist(rows), file.path(out, "r0309_libsize_definition.tsv"), sep = "\t")

## R03-02, GSE126848
ix <- which(meta$dataset == "GSE126848"); d <- dge[, ix]
g <- factor(meta$group_binary[ix], levels = c("Control", "Disease"))
Leff <- d$samples$lib.size * d$samples$norm.factors
X <- model.matrix(~g)
v <- voomWithQualityWeights(d, X); f <- eBayes(lmFit(v, X))
tb <- topTable(f, coef = 2, number = Inf, sort.by = "none")
cb <- tb$adj.P.Val < 0.05 & abs(tb$logFC) > 0.5
z <- which(rowSums(d$counts) == 0)
lc <- log2(Leff + 1)
wdiff <- sapply(z, function(i) { w <- v$weights[i, ]
  sum(w[g == "Control"] * lc[g == "Control"]) / sum(w[g == "Control"]) -
  sum(w[g == "Disease"] * lc[g == "Disease"]) / sum(w[g == "Disease"]) })
L <- c(L, log_line("GSE126848 base calls ", sum(cb), " (up ", sum(cb & tb$logFC > 0), ", down ", sum(cb & tb$logFC < 0), ")"))
L <- c(L, log_line("zero-read genes ", length(z), "; called ", sum(cb[z]), "; logFC range ", paste(signif(range(tb$logFC[z]), 7), collapse = " to ")))
L <- c(L, log_line("max |logFC - weighted ctrl-minus-disease mean log2(Leff+1)| over zero-read genes: ", signif(max(abs(tb$logFC[z] - wdiff)), 3)))
L <- c(L, log_line("unweighted ctrl-minus-disease mean log2(Leff+1): ", round(mean(lc[g == "Control"]) - mean(lc[g == "Disease"]), 4)))
L <- c(L, log_line("median raw lib control/disease: ", median(d$samples$lib.size[g == "Control"]), " / ", median(d$samples$lib.size[g == "Disease"])))
L <- c(L, log_line("median effective lib control/disease: ", round(median(Leff[g == "Control"])), " / ", round(median(Leff[g == "Disease"]))))
X2 <- model.matrix(~g + log2(Leff))
v2 <- voomWithQualityWeights(d, X2); f2 <- eBayes(lmFit(v2, X2))
t2 <- topTable(f2, coef = 2, number = Inf, sort.by = "none")
c2 <- t2$adj.P.Val < 0.05 & abs(t2$logFC) > 0.5
L <- c(L, log_line("~group + log2(effective lib): calls ", sum(c2), "; kept of base ", sum(c2 & cb),
                   "; zero-read calls kept ", sum(c2[z] & cb[z]), "; all-gene Spearman ", round(cor(tb$logFC, t2$logFC, method = "spearman"), 4)))
writeLines(L, file.path(out, "repair_libdef_log.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
