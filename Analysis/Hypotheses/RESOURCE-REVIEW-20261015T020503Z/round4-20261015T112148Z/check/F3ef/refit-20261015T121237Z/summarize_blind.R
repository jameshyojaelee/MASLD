#!/usr/bin/env Rscript
# Summarize the blind refit (job 9050) into the numbers each proposed edit prints.
suppressPackageStartupMessages(library(data.table))
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F3ef/refit-20261015T121237Z")
ann <- fread(file.path(root, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))[, .(gene_id, gene_name, ensembl_base)]
rd <- function(fit, con) {
  x <- fread(file.path(out, "fits", paste0(fit, "__", con, ".tsv.gz")))
  merge(x, ann, by = "gene_id", all.x = TRUE)
}
sig <- function(x) x[FDR < 0.05, gene_id]
lines <- character(0)
say <- function(...) { s <- paste0(...); lines <<- c(lines, s); cat(s, "\n") }

# 1. Base fits versus the released table
rel <- fread(file.path(root, "figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"),
             select = c("gene_id_versioned", "contrast", "axis", "logFC", "FDR"))
map_rel <- c(F0_to_F1 = "F0_to_F1", F1_to_F2 = "F1_to_F2", F2_to_F3 = "F2_to_F3", F3_to_F4 = "F3_to_F4",
             F2_vs_F0 = "F2_vs_F0", F3_vs_F0 = "F3_vs_F0", F4_vs_F0 = "F4_vs_F0")
for (con in names(map_rel)) {
  b <- rd("base", con)
  r <- rel[contrast == map_rel[[con]] & axis %in% c("fibrosis_adjacent", "fibrosis")][!duplicated(gene_id_versioned)]
  j <- merge(b, r, by.x = "gene_id", by.y = "gene_id_versioned")
  say("base ", con, ": n_bh05=", sum(b$FDR < 0.05), " released=", sum(r$FDR < 0.05),
      " max|dlogFC|=", signif(max(abs(j$logFC.x - j$logFC.y)), 3), " n_units=", b$n_units[1])
}
base <- list(F1 = rd("base", "F0_to_F1"), F2 = rd("base", "F2_vs_F0"), F3 = rd("base", "F3_vs_F0"), F4 = rd("base", "F4_vs_F0"))
base_sig <- lapply(base, sig)

# 2-4. F1 vs F0: cohort omission, single cohort, complexity covariates
for (dc in c("GSE130970", "GSE135251", "GSE162694")) {
  x <- rd(paste0("drop_", dc), "F1_vs_F0")
  say("F1 drop ", dc, ": n_units=", x$n_units[1], " n_bh05=", sum(x$FDR < 0.05), " of which in released 1338: ",
      length(intersect(sig(x), base_sig$F1)))
  y <- rd(paste0("alone_", dc), "F1_vs_F0")
  say("F1 alone ", dc, ": n_units=", y$n_units[1], " n_bh05=", sum(y$FDR < 0.05))
}
f1cov <- list()
for (cv in c("lib", "gini", "nf", "joint")) {
  x <- rd(cv, "F0_to_F1"); f1cov[[cv]] <- x
  say("F1 cov ", cv, ": n_bh05=", sum(x$FDR < 0.05), " in released 1338: ", length(intersect(sig(x), base_sig$F1)))
}
# 5-6. Adjacent and vs-F0 under covariates
for (cv in c("lib", "gini", "nf", "joint")) for (con in c("F1_to_F2", "F2_to_F3", "F3_to_F4", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")) {
  x <- rd(cv, con)
  say("cov ", cv, " ", con, ": n_bh05=", sum(x$FDR < 0.05))
}
# 7. F2-F4 vs F0 without each cohort
for (dc in c("GSE130970", "GSE135251", "GSE162694")) for (k in c("F2", "F3", "F4")) {
  x <- rd(paste0("drop_", dc), paste0(k, "_vs_F0"))
  say("drop ", dc, " ", k, "_vs_F0: n_units=", x$n_units[1], " n_bh05=", sum(x$FDR < 0.05),
      " in released: ", length(intersect(sig(x), base_sig[[k]])), " of ", length(base_sig[[k]]))
}
# 8. Interaction models
for (k in c("F1", "F4")) {
  x <- rd("interaction", paste0(k, "_vs_F0"))
  b <- base[[k]]
  j <- merge(b[, .(gene_id, logFC, FDR)], x[, .(gene_id, FDR_int, b_GSE130970, b_GSE135251, b_GSE162694)], by = "gene_id")
  calls <- j[FDR < 0.05]
  same3 <- calls[, sign(b_GSE130970) == sign(logFC) & sign(b_GSE135251) == sign(logFC) & sign(b_GSE162694) == sign(logFC)]
  say("interaction ", k, ": genome-wide FDR_int<0.05=", sum(x$FDR_int < 0.05), "; among ", nrow(calls), " pooled calls: ",
      sum(calls$FDR_int < 0.05), "; same sign all three (interaction coefs) = ", round(mean(same3), 4),
      " (", sum(same3), "); per cohort GSE130970/GSE135251/GSE162694 = ",
      paste(round(c(mean(sign(calls$b_GSE130970) == sign(calls$logFC)), mean(sign(calls$b_GSE135251) == sign(calls$logFC)),
                    mean(sign(calls$b_GSE162694) == sign(calls$logFC))), 3), collapse = "/"))
}
# Same-sign check for F1 using separate single-cohort fits instead of interaction coefficients
al <- lapply(c("GSE130970", "GSE135251", "GSE162694"), function(dc) rd(paste0("alone_", dc), "F1_vs_F0")[, .(gene_id, l = logFC)])
jj <- Reduce(function(a, b) merge(a, b, by = "gene_id"), c(list(base$F1[FDR < 0.05, .(gene_id, logFC)]), al))
setnames(jj, c("gene_id", "logFC", "a1", "a2", "a3"))
say("F1 same sign all three (single-cohort fits) = ", round(jj[, mean(sign(a1) == sign(logFC) & sign(a2) == sign(logFC) & sign(a3) == sign(logFC))], 4))
# 9. F0-reference arms
for (arm in c("armB_diseaseF0", "armC_ctrlterm")) for (k in paste0("F", 1:4)) {
  x <- rd(arm, paste0(k, "_vs_F0"))
  say(arm, " ", k, "_vs_F0: n_ref=", x$n_ref[1], " n_cmp=", x$n_cmp[1], " n_bh05=", sum(x$FDR < 0.05))
}
# 10. NAS
for (con in c("NAS1-2_vs_NAS0", "NAS3-4_vs_NAS0", "NAS5-8_vs_NAS0", "NAS3-4_vs_NAS1-2", "NAS5-8_vs_NAS3-4")) {
  x <- rd("nas", con)
  say("nas ", con, ": n_ref=", x$n_ref[1], " n_cmp=", x$n_cmp[1], " n_bh05=", sum(x$FDR < 0.05),
      " up/down=", sum(x$FDR < 0.05 & x$logFC > 0), "/", sum(x$FDR < 0.05 & x$logFC < 0))
}

# 11. Hero genes
hero <- c("GNMT", "MAT1A", "CYP2C19")
fits <- c(list(base = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")),
          setNames(rep(list(c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")), 4), c("lib", "gini", "nf", "joint")),
          setNames(rep(list(paste0("F", 1:4, "_vs_F0")), 5), c("drop_GSE130970", "drop_GSE135251", "drop_GSE162694", "armB_diseaseF0", "armC_ctrlterm")))
H <- rbindlist(lapply(names(fits), function(f) rbindlist(lapply(fits[[f]], function(con) {
  x <- rd(f, con)[gene_name %in% hero, .(gene = gene_name, logFC, q = FDR)]
  x[, `:=`(fit = f, contrast = con)]
}))))
H[contrast == "F0_to_F1", contrast := "F1_vs_F0"]
fwrite(H, file.path(out, "hero_genes_blind.tsv"), sep = "\t")
late <- H[contrast %in% c("F2_vs_F0", "F3_vs_F0", "F4_vs_F0")]
say("hero F2-F4 vs F0, max q over base/arms/covariates: ", signif(late[!grepl("^drop", fit), max(q)], 4),
    "; max q incl. cohort omission: ", signif(late[, max(q)], 4), "; all logFC<0: ", all(late$logFC < 0))
print(dcast(H, gene + contrast ~ fit, value.var = "q"))

# 12. Figure 3F marks
t12 <- fread(file.path(root, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"))
pos <- unique(t12[susie_state_t12 == "positive", ensembl])
adjc <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
getadj <- function(fit) rbindlist(lapply(adjc, function(con) rd(fit, con)[, .(gene_id, ensembl_base, gene_name, contrast = con, logFC, FDR)]))
A <- getadj("base")[ensembl_base %in% pos]
say("3F: positive set ", length(pos), "; in stage table ", uniqueN(A$ensembl_base))
sel <- A[, .(any = any(FDR < 0.05)), by = gene_id][any == TRUE, gene_id]
A54 <- A[gene_id %in% sel]
say("3F base: genes ", length(sel), ", effects ", nrow(A54), ", marks ", sum(A54$FDR < 0.05), " (",
    paste(A54[FDR < 0.05, .N, by = contrast][match(adjc, contrast), N], collapse = "/"), ")")
relm <- fread(file.path(root, "figures/main/fig3_bulk_transcriptomics/source_tables/fig3f_stage_genetics_matrix.tsv"))
say("3F released matrix: genes ", uniqueN(relm$gene_id_versioned), "; identical gene set: ", setequal(relm$gene_id_versioned, sel))
only01 <- A54[, .(m = paste(contrast[FDR < 0.05], collapse = ",")), by = .(gene_id, gene_name)][m == "F0_to_F1"]
say("3F genes whose only mark is F0_to_F1: ", nrow(only01))
J <- getadj("joint"); B1 <- rd("armB_diseaseF0", "F1_vs_F0")[, .(gene_id, FDR_B = FDR, lfc_B = logFC)]
D1 <- rd("drop_GSE130970", "F1_vs_F0")[, .(gene_id, FDR_D = FDR)]
M <- merge(A54[FDR < 0.05], J[, .(gene_id, contrast, FDR_J = FDR, lfc_J = logFC)], by = c("gene_id", "contrast"))
M <- merge(M, B1, by = "gene_id", all.x = TRUE); M <- merge(M, D1, by = "gene_id", all.x = TRUE)
M[, okB := fifelse(contrast == "F0_to_F1", FDR_B < 0.05 & sign(lfc_B) == sign(logFC), TRUE)]
M[, okJ := FDR_J < 0.05 & sign(lfc_J) == sign(logFC)]
M[, filled := okB & okJ]
say("3F marks filled (disease-only F0 and joint, same sign): ", sum(M$filled), " (",
    paste(sapply(adjc, function(c) M[contrast == c, sum(filled)]), collapse = "/"), "); open ", sum(!M$filled),
    "; genes with filled mark ", uniqueN(M[filled == TRUE, gene_id]))
m01 <- M[contrast == "F0_to_F1"]
say("3F F0_to_F1 marks: ", nrow(m01), "; hold disease-only F0: ", sum(m01$okB), "; hold joint: ", sum(m01$okJ),
    "; hold without GSE130970: ", sum(m01$FDR_D < 0.05), " (min q ", signif(min(m01$FDR_D), 3), ")")
lost <- m01[gene_id %in% only01$gene_id & !okB]
say("3F only-F0F1 genes losing the call with disease-only F0: ", nrow(lost), " [", paste(sort(lost$gene_name), collapse = ","), "] q ",
    signif(min(lost$FDR_B), 3), "-", signif(max(lost$FDR_B), 3))
for (cv in c("lib", "gini", "nf")) {
  x <- rd(cv, "F0_to_F1")[, .(gene_id, f = FDR, l = logFC)]
  y <- merge(m01[, .(gene_id, logFC)], x, by = "gene_id")
  say("3F F0_to_F1 marks holding under ", cv, ": ", sum(y$f < 0.05 & sign(y$l) == sign(y$logFC)))
}
Jp <- J[ensembl_base %in% pos]
selJ <- Jp[, .(any = any(FDR < 0.05)), by = gene_id][any == TRUE, gene_id]
J41 <- Jp[gene_id %in% selJ]
say("3F reselection joint: genes ", length(selJ), ", effects ", nrow(J41), ", marks ", sum(J41$FDR < 0.05), " (",
    paste(sapply(adjc, function(c) J41[contrast == c, sum(FDR < 0.05)]), collapse = "/"), "); shared with base 54: ",
    length(intersect(selJ, sel)))
newB <- B1[FDR_B < 0.05 & gene_id %in% A[, unique(gene_id)] & !(gene_id %in% sel)]
say("3F primary genes entering with disease-only F0 at F0_to_F1: ", nrow(newB))
fwrite(M, file.path(out, "fig3f_marks_blind.tsv"), sep = "\t")
writeLines(lines, file.path(out, "blind_summary.txt"))
