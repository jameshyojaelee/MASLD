# DIG1 job 5 (exploratory). Which Hallmark pathways shift with MASLD inside each lineage, and do
# the two independent cohorts agree? limma-voom + camera (competitive test that allows for
# inter-gene correlation) on donor pseudobulk counts, ~ disease + sex. BH over the 50 Hallmarks
# within each lineage x cohort.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
OUT <- file.path(RUN, "out_hallmark")
if (dir.exists(OUT)) stop("refusing to overwrite")
dir.create(OUT)
V5 <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z")
J3 <- file.path(V5, "outputs/job3")
V3 <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
GMT <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt")
BULK_DE <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
annot <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
gl <- strsplit(readLines(GMT), "\t")
sets <- setNames(lapply(gl, function(x) unique(x[-(1:2)])), sapply(gl, `[`, 1))
stopifnot(length(sets) == 50)

dq <- fread(file.path(V5, "outputs/donor_qc.tsv"), colClasses = list(character = "id"))[included == TRUE][order(gsm)]
cellsF <- fread(file.path(J3, "cells_per_donor_set.tsv"))
m3 <- fread(file.path(V3, "donor_metadata.tsv"))
n3 <- fread(file.path(V3, "nuclei_per_donor_lineage.tsv"))
lin <- data.table(cohort = c(rep("GSE298719", 5), rep("GSE256398", 3)),
                  cell_type = c("hepatocytes", "cholangiocytes", "fibroblasts", "macrophages", "tcells",
                                "hepatocytes", "cholangiocytes", "fibroblasts"),
                  set = c("hep", "cholangiocytes", "fibroblasts", "macrophages", "tcells",
                          "hepatocyte", "cholangiocyte", "stellate_fibroblast"))
run_one <- function(i) {
  L <- lin[i]
  if (L$cohort == "GSE298719") {
    pb <- fread(file.path(J3, sprintf("pseudobulk_%s.tsv.gz", L$set)))
    meta <- data.table(donor = dq$gsm, disease = dq$disease, sex_f = dq$sex_f)
    meta[, n := cellsF[set == L$set][match(meta$donor, gsm), n_cells]]
  } else {
    pb <- fread(file.path(V3, sprintf("pseudobulk_%s.tsv.gz", L$set)))
    meta <- data.table(donor = m3$donor, disease = m3$disease, sex_f = as.integer(m3$sex == "female"))
    meta[, n := n3[[L$set]][match(meta$donor, n3$donor)]]
  }
  meta <- meta[n >= 30]
  cnt <- as.matrix(pb[, meta$donor, with = FALSE]); rownames(cnt) <- pb$gene_id
  sym <- annot$gene_name[match(sub("\\..*$", "", rownames(cnt)), annot$ensembl_base)]
  ok <- !is.na(sym) & !duplicated(sym)
  cnt <- cnt[ok, ]; rownames(cnt) <- sym[ok]
  X <- model.matrix(~ disease + sex_f, data = meta)
  y <- DGEList(cnt); y <- y[filterByExpr(y, design = X), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  v <- voom(y, X)
  idx <- ids2indices(sets, rownames(v))
  cam <- as.data.table(camera(v, idx, X, contrast = "disease"), keep.rownames = "hallmark")
  fit <- eBayes(lmFit(v, X))
  tt <- fit$t[, "disease"]
  cam[, `:=`(cohort = L$cohort, cell_type = L$cell_type, n_donors = nrow(meta),
             mean_t = sapply(hallmark, function(h) mean(tt[idx[[h]]])),
             signed = ifelse(Direction == "Up", 1, -1) * -log10(PValue))]
  list(cam = cam, t = data.table(cohort = L$cohort, cell_type = L$cell_type, gene = names(tt), t = tt))
}
res <- lapply(seq_len(nrow(lin)), run_one)
cam <- rbindlist(lapply(res, `[[`, "cam"))
setnames(cam, "FDR", "bh_q_within_lineage")
fwrite(cam, file.path(OUT, "hallmark_camera_by_lineage.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "t")), file.path(OUT, "gene_moderated_t.tsv.gz"), sep = "\t")

# Same Hallmark test on the 774-unit bulk statistic (camera on precomputed t: cameraPR)
bulk <- fread(BULK_DE)
bt <- tapply(bulk$t, bulk$symbol, function(x) x[1])
bt <- bt[!is.na(names(bt)) & names(bt) != ""]
bc <- as.data.table(cameraPR(bt, ids2indices(sets, names(bt))), keep.rownames = "hallmark")
bc[, signed := ifelse(Direction == "Up", 1, -1) * -log10(PValue)]
fwrite(bc, file.path(OUT, "hallmark_cameraPR_bulk774.tsv"), sep = "\t")

w <- dcast(cam, hallmark ~ paste(cohort, cell_type, sep = ":"), value.var = "signed")
w <- merge(w, bc[, .(hallmark, bulk774 = signed)], by = "hallmark")
fwrite(w, file.path(OUT, "hallmark_signed_log10p_wide.tsv"), sep = "\t")
pairs <- list(c("GSE298719:hepatocytes", "GSE256398:hepatocytes"), c("GSE298719:fibroblasts", "GSE256398:fibroblasts"),
              c("GSE298719:cholangiocytes", "GSE256398:cholangiocytes"), c("GSE298719:hepatocytes", "bulk774"),
              c("GSE256398:hepatocytes", "bulk774"), c("GSE298719:fibroblasts", "bulk774"), c("GSE256398:fibroblasts", "bulk774"),
              c("GSE298719:macrophages", "bulk774"))
fwrite(rbindlist(lapply(pairs, function(p) data.table(a = p[1], b = p[2],
  spearman = cor(w[[p[1]]], w[[p[2]]], method = "spearman", use = "complete.obs")))), file.path(OUT, "hallmark_agreement.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
