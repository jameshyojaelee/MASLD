# DIG1 job 10 (exploratory). Do the single-cell hepatocyte couplings appear across participants in
# the 774-unit five-cohort bulk RNA (units are genetic-identity clusters, not verified donors)?
# Fixed scores (z within cohort, mean over genes; sets fixed before reading these data):
#   methyl = GNMT, MAT1A, BHMT; sterol = 19 cholesterol-biosynthesis enzymes;
#   biliary = SOX9, EPCAM, KRT7, KRT19, CD24, SPP1, KRT23; p21 = CDKN1A;
#   hepatocyte content proxy = ALB, APOA1, APOC3, SERPINA1, TF, HP, AHSG, FGB (unchanged in sc).
# Primary family (Holm over 3), disease units only: partial Spearman of GNMT with sterol, biliary
# and CDKN1A, ranks residualised on cohort + sex + hepatocyte proxy. Secondary: + fibrosis stage;
# per-cohort rho; stage association of each score. Usage: Rscript 10_bulk_links.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
SETS <- list(methyl = c("GNMT","MAT1A","BHMT"),
  sterol = c("HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","TM7SF2","MSMO1","NSDHL","HSD17B7","EBP","SC5D","DHCR7","DHCR24"),
  biliary = c("SOX9","EPCAM","KRT7","KRT19","CD24","SPP1","KRT23"), p21 = "CDKN1A",
  hep_proxy = c("ALB","APOA1","APOC3","SERPINA1","TF","HP","AHSG","FGB"))
y <- readRDS(file.path(B, "merged_dge.rds")); m <- readRDS(file.path(B, "meta_matched.rds"))
m <- m[match(colnames(y), m$sample_id)]; stopifnot(identical(m$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
sym <- an$gene_name[match(sub("\\..*$", "", rownames(y)), an$ensembl_base)]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE, prior.count = 1)
need <- unique(unlist(SETS)); have <- need[need %in% sym]
fwrite(data.table(gene = need, present = need %in% sym), file.path(out, "genes_present.tsv"), sep = "\t")
G <- lc[match(have, sym), , drop = FALSE]; rownames(G) <- have
zc <- function(v) ave(v, m$dataset, FUN = function(a) (a - mean(a)) / sd(a))
Z <- t(apply(G, 1, zc))
for (s in names(SETS)) m[, (s) := colMeans(Z[intersect(SETS[[s]], have), , drop = FALSE])]
m[, GNMT := Z["GNMT", ]]
D <- m[group_binary == "Disease"]
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(a = a, b = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df))
}
base <- c("dataset", "sex_final", "hep_proxy")
prim <- rbindlist(lapply(c("sterol", "biliary", "p21"), function(b) pcor(D, "GNMT", b, base)))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_gnmt_pairs_disease_units.tsv"), sep = "\t")
sec <- rbindlist(c(lapply(c("sterol", "biliary", "p21", "methyl"), function(b) pcor(D, "GNMT", b, c(base, "fibrosis_stage"))),
                   lapply(c("sterol", "biliary", "p21"), function(b) pcor(D, "GNMT", b, c("dataset", "sex_final"))),
                   list(pcor(D, "sterol", "biliary", base), pcor(D, "sterol", "p21", base), pcor(D, "biliary", "p21", base))))
fwrite(sec, file.path(out, "secondary_pairs.tsv"), sep = "\t")
per <- rbindlist(lapply(split(D, D$dataset), function(z) {
  v <- if (uniqueN(z$sex_final) > 1) c("sex_final", "hep_proxy") else "hep_proxy"   # single-sex cohorts
  rbindlist(lapply(c("sterol", "biliary", "p21"), function(b) cbind(dataset = z$dataset[1], pcor(z, "GNMT", b, v)))) }))
fwrite(per, file.path(out, "per_cohort_gnmt_pairs.tsv"), sep = "\t")
stage <- rbindlist(lapply(c("GNMT", "methyl", "sterol", "biliary", "p21", "hep_proxy"), function(s) pcor(D, s, "fibrosis_stage", c("dataset", "sex_final"))))
fwrite(stage, file.path(out, "stage_association_disease_units.tsv"), sep = "\t")
grp <- rbindlist(lapply(c("GNMT", "methyl", "sterol", "biliary", "p21", "hep_proxy"), function(s) {
  f <- summary(lm(reformulate(c("group_binary", "dataset", "sex_final"), s), data = m))$coefficients["group_binaryDisease", ]
  data.table(score = s, beta = f[1], se = f[2], p = f[4], n = nrow(m)) }))
fwrite(grp, file.path(out, "disease_vs_control_scores.tsv"), sep = "\t")
fwrite(m[, c("sample_id", "dataset", "group_binary", "fibrosis_stage", names(SETS), "GNMT"), with = FALSE], file.path(out, "unit_scores.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(prim); print(stage); print(grp)
