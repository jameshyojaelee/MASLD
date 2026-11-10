#!/usr/bin/env Rscript
# Code review U7 (H51): independent b, permutation p and Holm p for the GSE213621 grouped-stage comparison,
# written from METHODS.md "Grouped-stage comparison in GSE213621".
# (1) D0: three discovery cohorts, disease units only, ~ dataset + inferred_sex + bin, voomWithQualityWeights on the stored
#     counts and RLE factors; BH within contrast. (2) G0: GSE213621 disease units, ~ inferred_sex + bin.
# (3) b = OLS slope of G0 logFC on D0 logFC over D0 BH hits kept in GSE213621 (filterByExpr on GSE213621; count compared to 8,449).
# (4) Same-seed check: redraw 1,000 within-stratum label permutations (strata = accession quartile x inferred sex,
#     seed 20260923, Mersenne-Twister) and compare with the stored draws; refit draws 1-3 and compare with stored null coefficients.
# (5) p = (1 + k)/(B + 1) from all 1,000 stored null fits over the same genes; Holm over {F3F4, F2}.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
args <- commandArgs(trailingOnly = TRUE); OUT <- args[1]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
H51 <- file.path(ROOT, "Analysis/Hypotheses/H51-20261015T112617Z")
ACC <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-03/accession_order_bins.tsv")
log <- file(file.path(OUT, "h51_check_log.txt"), "w"); say <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, log) }
inputs <- c(dge = file.path(REL, "merged_dge.rds"), meta = file.path(REL, "meta_matched.rds"), acc = ACC,
            perm = file.path(H51, "stage1/perm_labels.tsv.gz"), nullc = file.path(H51, "stage1/null_coefficients.rds"))
fwrite(data.table(id = names(inputs), path = inputs, sha256 = vapply(inputs, function(f) system2("sha256sum", f, stdout = TRUE), "")),
       file.path(OUT, "inputs_sha256.tsv"), sep = "\t")

dge <- readRDS(inputs["dge"]); meta <- as.data.table(readRDS(inputs["meta"]))
stopifnot(setequal(colnames(dge), meta$sample_id), !anyDuplicated(meta$sample_id))  # same set; order differs, columns selected by name
say("stored DGEList: ", nrow(dge), " genes x ", ncol(dge), " units")
meta[, bin := fifelse(fibrosis_stage %in% 0:1, "F0F1", fifelse(fibrosis_stage == 2, "F2", fifelse(fibrosis_stage %in% 3:4, "F3F4", NA_character_)))]
# GSE213621 records bins as 1/2/3
meta[dataset == "GSE213621", bin := c("F0F1", "F2", "F3F4")[fibrosis_stage]]
fitc <- function(dg, X) { v <- voomWithQualityWeights(dg, X); eBayes(lmFit(v, X)) }

D <- meta[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & group_binary == "Disease"]
say("D0 units: ", nrow(D), " (", paste(D[, .N, by = dataset][order(dataset), paste(dataset, N)], collapse = ", "), "); bins ", paste(D[, .N, by = bin][order(bin), paste(bin, N)], collapse = ", "))
D[, bin := factor(bin, levels = c("F0F1", "F2", "F3F4"))]
XD <- model.matrix(~ dataset + inferred_sex + bin, data = D)
fD <- fitc(dge[, D$sample_id, keep.lib.sizes = TRUE], XD)

G <- meta[dataset == "GSE213621" & group_binary == "Disease"]
acc <- fread(ACC); G <- merge(G, acc[, .(sample_id, acc_q, acc_rank, bin_acc = bin, sex_acc = inferred_sex)], by = "sample_id")
stopifnot(nrow(G) == 258L, all(G$bin == G$bin_acc), all(G$inferred_sex == G$sex_acc))
setorder(G, acc_rank)
say("G0 units: ", nrow(G), "; bins ", paste(G[, .N, by = bin][order(bin), paste(bin, N)], collapse = ", "))
G[, bin := factor(bin, levels = c("F0F1", "F2", "F3F4"))]
XG <- model.matrix(~ inferred_sex + bin, data = G)
dgG <- dge[, G$sample_id, keep.lib.sizes = TRUE]
keepG <- filterByExpr(dgG, XG)
say("GSE213621 filterByExpr (design) keeps ", sum(keepG), " of ", nrow(dgG))
fG <- fitc(dgG[keepG, , keep.lib.sizes = TRUE], XG)

res <- list()
for (cf in c("binF3F4", "binF2")) {
  tD <- topTable(fD, coef = cf, number = Inf, sort.by = "none"); tG <- topTable(fG, coef = cf, number = Inf, sort.by = "none")
  hits <- rownames(tD)[tD$adj.P.Val < 0.05]; gs <- intersect(hits, rownames(tG))
  b <- unname(coef(lm(tG[gs, "logFC"] ~ tD[gs, "logFC"]))[2])
  say(cf, ": D0 BH hits ", length(hits), "; kept in GSE213621 ", length(gs), "; b = ", sprintf("%.4f", b))
  res[[cf]] <- list(gs = gs, b = b, xD = tD[gs, "logFC"])
}

# same-seed permutation draws
set.seed(20260923L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
G[, stratum := paste(acc_q, inferred_sex, sep = "_")]
st <- sort(unique(G$stratum)); P <- matrix(NA_character_, 1000L, nrow(G), dimnames = list(NULL, G$sample_id))
for (k in 1:1000) { lab <- as.character(G$bin); for (s in st) { i <- which(G$stratum == s); lab[i] <- lab[i][sample.int(length(i))] }; P[k, ] <- lab }
stored <- fread(inputs["perm"])
same_order <- identical(names(stored)[-1], G$sample_id)
say("stored permutation columns in acc_rank order: ", same_order)
SP <- as.matrix(stored[, G$sample_id, with = FALSE])
say("same-seed permutations identical to stored: ", sum(rowSums(SP == P) == ncol(P)), " of 1000")

# refit draws 1-3 and compare with stored null coefficients
NULLC <- readRDS(inputs["nullc"])
say("stored null matrices: ", paste(names(NULLC), sapply(NULLC, function(m) paste(dim(m), collapse = "x")), collapse = "; "))
for (k in 1:3) {
  Gi <- copy(G); Gi[, bin := factor(P[k, ], levels = c("F0F1", "F2", "F3F4"))]
  fk <- fitc(dgG[keepG, , keep.lib.sizes = TRUE], model.matrix(~ inferred_sex + bin, data = Gi))
  for (cf in c("binF3F4", "binF2")) {
    sc <- NULLC[[cf]][, k]; com <- intersect(names(sc), rownames(fk))
    say(sprintf("null draw %d %s: genes stored %d mine %d common %d; max |diff| %.3g", k, cf, length(sc), nrow(fk), length(com),
                max(abs(sc[com] - fk$coefficients[com, cf]))))
  }
}
# p and Holm from stored null fits
pv <- c()
for (cf in c("binF3F4", "binF2")) {
  gs <- res[[cf]]$gs; M <- NULLC[[cf]]; miss <- setdiff(gs, rownames(M))
  if (length(miss)) say(cf, ": ", length(miss), " slope genes missing from stored null")
  gs2 <- intersect(gs, rownames(M)); x <- res[[cf]]$xD[match(gs2, gs)]
  xc <- x - mean(x); nb <- as.vector(crossprod(xc, M[gs2, , drop = FALSE])) / sum(xc^2)
  k <- sum(nb >= res[[cf]]$b); p <- (1 + k) / (length(nb) + 1); pv[cf] <- p
  say(sprintf("%s: null slopes mean %.4f sd %.4f max %.4f; k = %d of %d; one-sided p = %.6f (floor %.6f); minimum detectable slope 2.80 x SD = %.3f",
              cf, mean(nb), sd(nb), max(nb), k, length(nb), p, 1 / (length(nb) + 1), 2.80 * sd(nb)))
}
say("Holm over {F3F4, F2}: ", paste(names(pv), sprintf("%.6f", p.adjust(pv, "holm")), collapse = "; "))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt")); close(log)
