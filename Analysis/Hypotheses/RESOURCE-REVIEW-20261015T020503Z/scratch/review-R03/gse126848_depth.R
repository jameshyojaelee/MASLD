# Review R03: GSE126848 per-cohort disease-vs-control (23 control, 7 disease; constant sex) and library depth.
# limma-voomWithQualityWeights on fragment counts; BH over 24,115 genes; unit = genetic-identity cluster.
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
set.seed(20261015)
REL <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
S <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R03"
t0 <- Sys.time()
d <- readRDS(file.path(REL, "merged_dge.rds"))
m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))[match(colnames(d), sample_id)]
cx <- fread(file.path(S, "complexity_774.tsv"))[match(colnames(d), sample_id)]
idx <- which(m$dataset == "GSE126848")
dd <- d[, idx]
grp <- factor(m$group_binary[idx], levels = c("Control", "Disease"))
zl <- as.numeric(scale(log2(cx$lib_size[idx])))
cat("n", length(idx), "control", sum(grp == "Control"), "disease", sum(grp == "Disease"), "\n")
cat("lib size control median", median(cx$lib_size[idx][grp == "Control"]), "disease median", median(cx$lib_size[idx][grp == "Disease"]), "\n")
# 1. baseline (guard: 6,147 calls)
X0 <- model.matrix(~ grp)
f0 <- eBayes(lmFit(voomWithQualityWeights(dd, X0), X0))
b <- topTable(f0, coef = "grpDisease", number = Inf, sort.by = "none")
call0 <- b$adj.P.Val < 0.05 & abs(b$logFC) > 0.5
cat("baseline calls", sum(call0), "time", format(Sys.time() - t0), "\n")
# 2. depth slope learned in controls only (no disease libraries used)
ci <- which(grp == "Control"); Xc <- model.matrix(~ zl[ci])
fc <- eBayes(lmFit(voomWithQualityWeights(dd[, ci], Xc), Xc))
slope <- fc$coefficients[, 2]
dz <- mean(zl[grp == "Disease"]) - mean(zl[grp == "Control"])
pred <- slope * dz
cat("mean z(log2 lib) disease minus control", dz, "\n")
ct_all <- cor.test(pred, b$logFC, method = "spearman", exact = FALSE)
ct_call <- cor.test(pred[call0], b$logFC[call0], method = "spearman", exact = FALSE)
sgn <- mean(sign(pred[call0]) == sign(b$logFC[call0]))
cat("Spearman depth-predicted vs observed disease logFC, all genes:", ct_all$estimate, "p", ct_all$p.value, "\n")
cat("Spearman among baseline calls:", ct_call$estimate, "p", ct_call$p.value, " sign agreement among calls:", sgn, "\n")
sl <- topTable(fc, coef = 2, number = Inf, sort.by = "none")
cat("genes with BH q<0.05 for depth slope within 23 controls:", sum(sl$adj.P.Val < 0.05), "; among baseline calls:", sum(sl$adj.P.Val[call0] < 0.05), "\n")
out <- data.table(gene = rownames(b), disease_logFC = b$logFC, disease_padj = b$adj.P.Val, baseline_call = call0,
                  control_depth_slope = slope, control_depth_slope_padj = sl$adj.P.Val, depth_predicted_logFC = pred)
fwrite(out, file.path(S, "gse126848_depth_check.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(S, "gse126848_sessionInfo.txt"))
cat("total time", format(Sys.time() - t0), "\n")
