suppressPackageStartupMessages({library(data.table); library(edgeR)})
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
gm <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
us <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/unit_scores.tsv"))
zs <- function(v) (v - mean(v)) / sd(v)
run <- function(cnt, label) {
  sym <- gm$gene_name[match(sub("\\..*$", "", rownames(cnt)), gm$ensembl_base)]
  lib <- colSums(cnt)
  for (nm in c("TMM", "RLE", "none")) {
    nf <- calcNormFactors(DGEList(cnt, lib.size = lib), method = nm)$samples$norm.factors
    for (pr in c(0.5, 1, 2)) {
      L <- cpm(DGEList(cnt, lib.size = lib, norm.factors = nf), log = TRUE, prior.count = pr)
      for (sub in c("all", "staged")) {
        cols <- if (sub == "all") colnames(cnt) else intersect(colnames(cnt), us$sample_alias)
        g <- function(x) { v <- L[which(sym == x)[1], cols]; zs(v) }
        gnmt <- g("GNMT"); p21 <- g("CDKN1A")
        bil <- colMeans(do.call(rbind, lapply(c("SOX9","EPCAM","KRT7","KRT19","CD24","SPP1","KRT23"), g)))
        u <- us[match(cols, sample_alias)]; ok <- !is.na(u$sample_alias)
        cat(label, nm, pr, sub, "GNMT", signif(max(abs(u$j10_GNMT[ok] - gnmt[ok])),3), "p21", signif(max(abs(u$j10_p21[ok] - p21[ok])),3),
            "bil", signif(max(abs(u$j10_biliary[ok] - bil[ok])),3), "\n")
      }
    }
  }
}
d <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds")); run(d$counts, "dge24115")
a <- readRDS(file.path(SS, "out_s2_qc/counts_all_genes.rds")); if (is.list(a) && !is.matrix(a)) a <- a$counts
a <- as.matrix(a); a <- a[, colnames(d$counts)]; cat("all genes dim", dim(a), "\n"); run(a, "allgenes")
