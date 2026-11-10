# DIG1 job 18a (exploratory). Between-participant partners for the job 17 hepatocyte anchors
# (GNMT + 40 null + 40 changed), so the job 16 agreement can be compared against control anchors.
# Bulk: 636 disease units, partial Spearman of each anchor with every gene, adjusted cohort + sex
# + hepatocyte proxy (units are genetic-identity clusters, not verified donors). Protein:
# PXD051911, 46 MASLD, adjusted batch + age + BMI + sex, anchors quantified only.
# Usage: Rscript 18a_bulk_protein_anchor_partners.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename")
stopifnot(nrow(md) == 58L)
md[, `:=`(sample_id = liver_proteomics_filename, masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
          age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), NAS = as.numeric(nafld_activity_score))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- raw$Genes
x <- x[rowMeans(is.na(x)) <= .5, , drop = FALSE]
x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
gene <- t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(gene) <- md$sample_id
present <- function(gs) intersect(gs, rownames(gene))
anc <- unique(fread(file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_anchor_control/anchor_control.tsv"))[lineage == "Hepatocytes", .(anchor, kind)])
fwrite(anc, file.path(out, "anchors.tsv"), sep = "\t")
part_mat <- function(L, X, anc_genes) {
  # L: features x samples (complete), X: design; returns features x anchors partial Spearman
  Rk <- apply(L, 1, rank)
  Rr <- qr.resid(qr(X), Rk); Rr <- sweep(Rr, 2, sqrt(colSums(Rr^2)) + 1e-12, "/")
  a <- intersect(anc_genes, colnames(Rr))
  crossprod(Rr, Rr[, a, drop = FALSE])
}
M <- md[masld == 1]
Pm <- gene[, M$sample_id]
Pm <- Pm[rowSums(is.na(Pm)) == 0, ]
Xp <- model.matrix(~ batch + age + bmi_num + sex, data = M)
okp <- complete.cases(M[, .(batch, age, bmi_num, sex)])
Cp <- part_mat(Pm[, okp], Xp, anc$anchor)
fwrite(data.table(gene = rownames(Cp), Cp), file.path(out, "protein_partner_rho.tsv.gz"), sep = "\t")

B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
suppressPackageStartupMessages(library(edgeR))
y <- readRDS(file.path(B, "merged_dge.rds")); mb <- readRDS(file.path(B, "meta_matched.rds"))
mb <- mb[match(colnames(y), mb$sample_id)]; stopifnot(identical(mb$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
sym <- an$gene_name[match(sub("\\..*$", "", rownames(y)), an$ensembl_base)]
ok <- !is.na(sym) & !duplicated(sym); y <- y[ok, ]; rownames(y) <- sym[ok]
y <- y[filterByExpr(y, group = mb$group_binary), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
lc <- cpm(y, log = TRUE, prior.count = 1)
zc <- function(v) ave(v, mb$dataset, FUN = function(a) (a - mean(a)) / sd(a))
HEP <- c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB")
mb[, hep_proxy := colMeans(t(apply(lc[HEP, ], 1, zc)))]
D <- which(mb$group_binary == "Disease")
zd <- mb[D]; okb <- complete.cases(zd[, .(dataset, sex_final, hep_proxy)])
Xb <- model.matrix(~ dataset + sex_final + hep_proxy, data = zd[okb])
Cb <- part_mat(lc[, D][, okb], Xb, anc$anchor)
fwrite(data.table(gene = rownames(Cb), Cb), file.path(out, "bulk_partner_rho.tsv.gz"), sep = "\t")
fwrite(data.table(gene = rownames(lc), bulk_mean_logcpm = rowMeans(lc[, D])), file.path(out, "bulk_gene_means.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("anchors in protein:", ncol(Cp), " in bulk:", ncol(Cb), "\n")
