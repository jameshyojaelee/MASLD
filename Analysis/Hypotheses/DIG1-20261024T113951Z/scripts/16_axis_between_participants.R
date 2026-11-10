# DIG1 job 16 (exploratory, descriptive). Does the within-hepatocyte GNMT axis (job 9: per-donor
# partial Spearman of GNMT with each gene across GSE298719 hepatocytes) reappear between
# participants? A. PXD051911 liver protein, 46 MASLD participants: partial Spearman of every
# protein with GNMT protein, adjusted batch + age + BMI + sex (and + fibrosis as a sensitivity).
# B. 774-unit bulk RNA, 636 disease units: partial Spearman of every gene with GNMT, adjusted
# cohort + sex + hepatocyte proxy (units are genetic-identity clusters, not verified donors).
# Agreement: Spearman of the between-participant rho with the within-cell mean z over shared genes,
# for all genes and for hepatocyte-enriched genes (GSE298719 pseudobulk: hepatocyte CPM >= 2x the
# highest of cholangiocytes, fibroblasts, macrophages, T cells). Sign agreement among within-cell
# partners at BH < 0.05. Gene-gene dependence means no p-value is attached to the agreement.
# Usage: Rscript 16_axis_between_participants.R <new_out_dir>
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
md[, GNMTp := gene["GNMT", sample_id]]
base <- c("batch", "age", "bmi_num", "sex")
prank_all <- function(Y, z, a, vars) {
  # partial Spearman of each column of Y with z[[a]]; pairwise-complete rows handled per column
  sapply(seq_len(ncol(Y)), function(k) {
    ok <- complete.cases(z[, c(a, vars), with = FALSE]) & !is.na(Y[, k])
    if (sum(ok) < 20) return(NA_real_)
    xx <- model.matrix(reformulate(vars), data = z[ok])
    rr <- qr.resid(qr(xx), cbind(rank(z[[a]][ok]), rank(Y[ok, k])))
    cor(rr[, 1], rr[, 2]) })
}
M <- md[masld == 1]
Yp <- t(gene[, M$sample_id])
prot <- data.table(gene = colnames(Yp), rho_protein = prank_all(Yp, M, "GNMTp", base),
                   rho_protein_fib = prank_all(Yp, M, "GNMTp", c(base, "Fibrosis")))
prot <- prot[gene != "GNMT"]

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
mb[, GNMTb := lc["GNMT", ]]
D <- which(mb$group_binary == "Disease")
zd <- mb[D]; Xd <- model.matrix(~ dataset + sex_final + hep_proxy, data = zd)
okr <- complete.cases(zd[, .(dataset, sex_final, hep_proxy)])
Rk <- apply(lc[, D][, okr], 1, rank)
Q <- qr(Xd)
Rr <- qr.resid(Q, Rk); g0 <- qr.resid(Q, rank(zd$GNMTb[okr]))
bulk <- data.table(gene = rownames(lc), rho_bulk = as.numeric(crossprod(Rr, g0) / sqrt(colSums(Rr^2) * sum(g0^2))))
bulk <- bulk[gene != "GNMT"]

V5 <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z/outputs/job3")
gm <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
spec_cpm <- function(s) { p <- fread(file.path(V5, sprintf("pseudobulk_%s.tsv.gz", s))); m <- as.matrix(p[, -1]); c <- t(t(m) / colSums(m)) * 1e6
  data.table(gene = gm$gene_name[match(sub("\\..*$", "", p$gene_id), gm$ensembl_base)], cpm = rowMeans(c)) }
sp <- Reduce(function(a, b) merge(a, b, by = "gene"), lapply(c("hep", "cholangiocytes", "fibroblasts", "macrophages", "tcells"), function(s) {
  d <- spec_cpm(s)[!is.na(gene)][, .(cpm = mean(cpm)), by = gene]; setnames(d, "cpm", s); d }))
sp[, hep_enriched := hep >= 2 * pmax(cholangiocytes, fibroblasts, macrophages, tcells)]
wc <- fread(file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_hep_axes/gnmt_within_donor_partners.tsv.gz"))[, .(gene, mean_z, bh_q)]
all <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), list(wc, prot, bulk, sp[, .(gene, hep_enriched)]))
fwrite(all, file.path(out, "gnmt_axis_three_modalities.tsv.gz"), sep = "\t")
ag <- rbindlist(lapply(c("rho_protein", "rho_protein_fib", "rho_bulk"), function(col) rbindlist(lapply(c(FALSE, TRUE), function(he) {
  z <- all[!is.na(mean_z) & !is.na(get(col))]; if (he) z <- z[hep_enriched == TRUE]
  s <- z[bh_q < 0.05]
  data.table(between = col, hep_enriched_only = he, n_genes = nrow(z), spearman_with_within_cell_z = cor(z$mean_z, z[[col]], method = "spearman"),
             n_partners = nrow(s), frac_sign_agree = mean(sign(s$mean_z) == sign(s[[col]]))) }))))
fwrite(ag, file.path(out, "agreement.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(ag)
print(all[gene %in% c("MAT1A", "BHMT", "CTH", "CPS1", "ASS1", "ARG1", "CYP2E1", "CYP3A4", "KRT8", "KRT18", "AKR1B10", "ANXA2", "THBS1", "TIMP1", "ICAM1", "CHI3L1", "APOA1", "AHSG", "LRG1", "TNFRSF12A", "FGFR2", "SOX9")], digits = 2)
