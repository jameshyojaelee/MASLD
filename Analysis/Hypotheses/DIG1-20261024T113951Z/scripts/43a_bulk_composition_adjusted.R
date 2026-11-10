# DIG1 job 43a (exploratory, post hoc; input for 43b). Job 41 found no alignment between the genetic
# direction (liver eQTL x liver-enzyme GWAS) and the bulk disease direction. One reason could be cell
# composition: a liver eQTL mostly reflects hepatocytes, while a bulk disease change can come from
# more stromal or immune cells. Here: per gene, the disease-vs-control effect with and without
# composition proxies.
# Bulk: five-cohort fragment-count DGE (adopted 774-unit release), log CPM > 1 in half of units,
# z within cohort. OLS per gene: z ~ dataset + sex + group [+ proxies]. Proxies (mean z, recalled
# markers): hepatocyte ALB APOA1 APOC3 SERPINA1 TF HP AHSG FGB; stroma COL1A1 COL1A2 COL3A1 LUM DCN
# TAGLN ACTA2 MYL9; immune PTPRC CD74 LYZ CD68 CD14 CD3E; cholangiocyte KRT7 KRT19 EPCAM; endothelial
# PECAM1 VWF CLEC4G. Output only; no test here.
# Usage: Rscript 43a_bulk_composition_adjusted.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
y <- readRDS(file.path(B, "merged_dge.rds")); m <- readRDS(file.path(B, "meta_matched.rds"))
m <- m[match(colnames(y), m$sample_id)]; stopifnot(identical(m$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE, prior.count = 1)
lc <- lc[rowMeans(lc > 1) >= 0.5, ]
ens <- sub("\\..*$", "", rownames(lc)); sym <- an$gene_name[match(ens, an$ensembl_base)]
Z <- t(apply(lc, 1, function(v) ave(v, m$dataset, FUN = function(a) (a - mean(a)) / sd(a))))
PROX <- list(hep = c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"),
             stroma = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "TAGLN", "ACTA2", "MYL9"),
             immune = c("PTPRC", "CD74", "LYZ", "CD68", "CD14", "CD3E"), chol = c("KRT7", "KRT19", "EPCAM"),
             endo = c("PECAM1", "VWF", "CLEC4G"))
used <- list()
for (k in names(PROX)) { gg <- intersect(PROX[[k]], sym); used[[k]] <- gg
  m[, (paste0(k, "_proxy")) := colMeans(Z[match(gg, sym), , drop = FALSE])] }
fwrite(data.table(proxy = names(used), genes = vapply(used, paste, "", collapse = ",")), file.path(out, "proxies_used.tsv"), sep = "\t")
m[, grp := factor(group_binary)]
fit <- function(form) {
  X <- model.matrix(form, data = m); q <- qr(X); co <- qr.coef(q, t(Z)); res <- t(Z) - X %*% co
  df <- nrow(X) - q$rank; s2 <- colSums(res^2) / df; j <- grep("^grp", colnames(X))
  stopifnot(q$rank == ncol(X)); jj <- match(j, q$pivot); se <- sqrt(s2 * chol2inv(qr.R(q))[jj, jj]); b <- co[j, ]; tt <- b / se
  data.table(beta = b, t = tt, p = 2 * pt(-abs(tt), df)) }
u <- fit(~ dataset + sex_final + grp)
a <- fit(~ dataset + sex_final + grp + hep_proxy + stroma_proxy + immune_proxy + chol_proxy + endo_proxy)
res <- data.table(ensembl = ens, symbol = sym, beta_unadj = u$beta, t_unadj = u$t, p_unadj = u$p,
                  beta_adj = a$beta, t_adj = a$t, p_adj = a$p)
res[, `:=`(bh_unadj = p.adjust(p_unadj, "BH"), bh_adj = p.adjust(p_adj, "BH"))]
fwrite(res, file.path(out, "gene_disease_composition_adjusted.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("genes", nrow(res), "units", ncol(Z), "\n")
cat("sign agreement unadj vs adj (bh_unadj<0.05):", res[bh_unadj < 0.05, mean(sign(beta_unadj) == sign(beta_adj))], "\n")
print(res[symbol %in% c("THRB", "HSPA4", "ACTG1", "GAS6", "F2RL1", "PPM1K", "CFHR4", "MTARC2", "RORA", "SPINT2")], digits = 3)
