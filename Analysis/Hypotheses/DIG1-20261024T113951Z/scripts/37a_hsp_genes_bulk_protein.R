# DIG1 job 37a (exploratory; follows job 36). A published NASH biopsy study reports HSPB1 (HSP27)
# mRNA 2.24-fold lower with ballooning and low HSP27 immunostaining in ballooned hepatocytes
# (doi 10.1038/srep22528); another reports liver HSP70 and HSF1 protein falling with NAFLD stage
# (doi 10.1002/oby.20919). Job 36 found HSF1-target RNA rising with ballooning. Gene by gene, which
# heat-shock genes rise and which fall, in bulk RNA and in protein?
# Bulk: five-cohort fragment-count DGE as job 36a (log CPM, z within cohort). Protein: loaders as
# job 31 (quantile-normalised log2, gene median); composition markers from job 33.
# Primary family (Holm over 2), a direct check of the published HSPB1 claim, two-sided:
#   P1  GSE130970 disease units: HSPB1 RNA vs ballooning, partial Spearman | sex + hepatocyte proxy.
#   P2  PXD051911 MASLD (n 46): HSPB1 protein vs ballooning, partial Spearman | batch, age, BMI, sex.
# Descriptive (BH within each table): every gene in HSP_GENES vs ballooning (GSE130970) and NAS
# (three NAS cohorts pooled), each without and with a stroma RNA proxy (COL1A1 COL1A2 COL3A1 LUM DCN
# TAGLN ACTA2 MYL9); protein vs ballooning without and with composition markers; protein disease
# vs control in the three proteomes (difference in means / pooled SD, Welch t).
# Usage: Rscript 37a_hsp_genes_bulk_protein.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; EXT <- file.path(ROOT, "data/external")
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
HSP_GENES <- c("HSPB1", "CRYAB", "HSPB8", "HSPA1A", "HSPA1B", "HSPA6", "HSPA8", "HSPA4", "HSPA4L", "HSPH1", "DNAJB1", "DNAJA1",
               "BAG3", "HSP90AA1", "HSP90AB1", "STIP1", "CHORDC1", "AHSA1", "FKBP4", "HSPD1", "HSPE1", "HSF1", "UCHL1", "CPT1A",
               "PPARGC1A", "HSPA5")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(gene = a, b = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
y <- readRDS(file.path(B, "merged_dge.rds")); m <- readRDS(file.path(B, "meta_matched.rds"))
m <- m[match(colnames(y), m$sample_id)]; stopifnot(identical(m$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
sym <- an$gene_name[match(sub("\\..*$", "", rownames(y)), an$ensembl_base)]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE, prior.count = 1)
ok <- !is.na(sym) & !duplicated(sym); lc <- lc[ok, ]; rownames(lc) <- sym[ok]
lc <- lc[rowMeans(lc > 1) >= 0.5, ]
g <- rownames(lc)
hep <- intersect(c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"), g)
str <- intersect(c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "TAGLN", "ACTA2", "MYL9"), g)
hg <- intersect(HSP_GENES, g)
Z <- t(apply(lc[unique(c(hg, hep, str)), ], 1, function(v) ave(v, m$dataset, FUN = function(a) (a - mean(a)) / sd(a))))
m[, hep_proxy := colMeans(Z[hep, , drop = FALSE])][, stroma_proxy := colMeans(Z[str, , drop = FALSE])]
for (s in hg) m[, (s) := Z[s, ]]
m[, `:=`(Ballooning = as.numeric(ballooning_grade), NAS = as.numeric(nas_score))]
D <- m[group_binary == "Disease"]
G1 <- D[dataset == "GSE130970"]; GN <- D[!is.na(NAS) & dataset %in% c("GSE130970", "GSE135251", "GSE162694")]
bulk <- rbindlist(lapply(hg, function(s) rbind(
  cbind(analysis = "GSE130970 ballooning", pcor(G1, s, "Ballooning", c("sex_final", "hep_proxy"))),
  cbind(analysis = "GSE130970 ballooning", pcor(G1, s, "Ballooning", c("sex_final", "hep_proxy", "stroma_proxy"))),
  cbind(analysis = "NAS pooled", pcor(GN, s, "NAS", c("dataset", "sex_final", "hep_proxy"))),
  cbind(analysis = "NAS pooled", pcor(GN, s, "NAS", c("dataset", "sex_final", "hep_proxy", "stroma_proxy"))))))
bulk[, bh := p.adjust(p, "BH"), by = .(analysis, adjust)]
fwrite(bulk, file.path(out, "bulk_per_gene.tsv"), sep = "\t")
fwrite(data.table(gene = hg, mean_logcpm = rowMeans(lc[hg, ])), file.path(out, "bulk_gene_abundance.tsv"), sep = "\t")
to_gene <- function(x, gsym, maxmiss = 1) {
  ok <- !is.na(gsym) & gsym != "" & !grepl(";", gsym, fixed = TRUE)
  x <- x[ok, , drop = FALSE]; gsym <- gsym[ok]; x[!is.finite(x)] <- NA_real_
  keep <- rowMeans(is.na(x)) <= maxmiss; x <- x[keep, , drop = FALSE]; gsym <- gsym[keep]
  x <- normalizeBetweenArrays(x, method = "quantile")
  idx <- split(seq_len(nrow(x)), gsym)
  t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x)))) }
# PXD051911
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt")); md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
smp <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% smp], by = "liver_proteomics_filename"); setorder(md, liver_proteomics_filename)
x1 <- as.matrix(raw[, md$liver_proteomics_filename, with = FALSE]); storage.mode(x1) <- "double"; x1[x1 <= 0] <- NA; x1 <- log2(x1)
G1 <- to_gene(x1, raw$Genes, 0.5); colnames(G1) <- md$liver_proteomics_filename
a1 <- md[, .(sample = liver_proteomics_filename, grp = fifelse(saf_diagnosis != "No_MASLD", "disease", "control"),
             batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")), age = as.numeric(alder),
             bmi = as.numeric(bmi), sex = factor(gender), Ballooning = as.numeric(hepatocellular_ballooning_score))]
# PXD027722 (cirrhosis kept separate here)
pa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_proteome.csv"))
aa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_annotation_file.csv"), select = 1:3)
setnames(aa, c("patient", "sample", "group")); aa <- aa[sample %in% names(pa)]
x2 <- as.matrix(pa[, aa$sample, with = FALSE]); storage.mode(x2) <- "double"
G2 <- to_gene(x2, pa[["Gene name"]]); colnames(G2) <- aa$sample
a2 <- aa[, .(sample, grp = fcase(group == "control", "control", group == "NASH", "disease", group == "Cirrhosis", "cirrhosis"))]
# PXD065892
pb <- fread(file.path(EXT, "PXD065892/liver_20261017T003606Z/liver_sn_minprob_processed.tsv"))
x3 <- t(as.matrix(pb[, -(1:3)])); storage.mode(x3) <- "double"
G3 <- to_gene(x3, sub("~.*$", "", rownames(x3))); colnames(G3) <- pb$sample
a3 <- data.table(sample = pb$sample, grp = fcase(pb$group == "Healthy", "control", pb$group == "NAFLD", "disease", pb$group == "AIH", "AIH"))
stopifnot(!anyNA(a2$grp), !anyNA(a3$grp))
P <- list(PXD051911 = list(G1, a1), PXD027722 = list(G2, a2), PXD065892 = list(G3, a3))
comp <- fread(file.path(RUN, "out_cct_composition", "scores_PXD051911.tsv"))[, .(sample, MF, immune, macrophage, stroma, ductular, endothelial, hepatocyte)]
pz <- function(G, a, genes) { gg <- intersect(genes, rownames(G)); d <- copy(a)
  for (s in gg) d[, (s) := as.numeric(G[s, sample])]; list(d = d, genes = gg) }
q1 <- pz(G1, a1, HSP_GENES); M1 <- merge(q1$d, comp, by = "sample")[grp == "disease"]
BASE1 <- c("batch", "age", "bmi", "sex"); COMP <- c("MF", "immune", "macrophage", "stroma", "ductular", "endothelial", "hepatocyte")
prot <- rbindlist(lapply(q1$genes, function(s) rbind(
  cbind(analysis = "PXD051911 MASLD ballooning", pcor(M1, s, "Ballooning", BASE1), n_obs = sum(!is.na(M1[[s]]))),
  cbind(analysis = "PXD051911 MASLD ballooning", pcor(M1, s, "Ballooning", c(BASE1, COMP)), n_obs = sum(!is.na(M1[[s]]))))))
prot[, bh := p.adjust(p, "BH"), by = adjust]
fwrite(prot, file.path(out, "protein_ballooning_per_gene.tsv"), sep = "\t")
grp <- rbindlist(lapply(list(c("PXD051911", "disease"), c("PXD027722", "disease"), c("PXD027722", "cirrhosis"), c("PXD065892", "disease"), c("PXD065892", "AIH")), function(k) {
  q <- pz(P[[k[1]]][[1]], P[[k[1]]][[2]], HSP_GENES)
  rbindlist(lapply(q$genes, function(s) { a <- q$d[grp == k[2], get(s)]; b <- q$d[grp == "control", get(s)]; a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 3 || length(b) < 3) return(NULL)
    data.table(cohort = k[1], contrast = paste(k[2], "vs control"), gene = s, n1 = length(a), n0 = length(b),
               d = (mean(a) - mean(b)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2)),
               p = t.test(a, b)$p.value) })) }))
grp[, bh := p.adjust(p, "BH"), by = .(cohort, contrast)]
fwrite(grp, file.path(out, "protein_group_per_gene.tsv"), sep = "\t")
prim <- rbind(cbind(id = "P1", bulk[analysis == "GSE130970 ballooning" & gene == "HSPB1" & adjust == "sex_final+hep_proxy", .(gene, b, adjust, n, rho, p)]),
              cbind(id = "P2", prot[gene == "HSPB1" & adjust == paste(BASE1, collapse = "+"), .(gene, b, adjust, n, rho, p)]))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220); print(prim, digits = 3)
print(dcast(bulk, gene ~ analysis + adjust, value.var = "rho"), digits = 2)
print(dcast(prot, gene ~ adjust, value.var = "rho"), digits = 2)
print(dcast(grp, gene ~ cohort + contrast, value.var = "d"), digits = 2)
