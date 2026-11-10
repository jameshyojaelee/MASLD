# DIG1 job 42b (exploratory; follows job 41). THRB is colocalized with GGT, and the allele that raises
# liver THRB lowers GGT (job 41b). THRB and DIO1 mRNA fall with NASH score in human biopsies
# (doi 10.1530/EC-18-0499, 10.1530/ETJ-24-0080; source-owned). Question: does the liver's thyroid
# hormone target output follow THRB from patient to patient, beyond disease state and hepatocyte
# content? If THRB level limits signalling, the two should co-vary.
# Target set TR4 = CPT1A, ANGPTL4, DIO1, ME1: direct THR targets in human hepatic cells (CPT1A,
# ANGPTL4, DIO1) and rat liver (Dio1, Me1) in doi 10.1371/journal.pone.0240338. Caution: CPT1A and
# ANGPTL4 are also PPAR-alpha targets [recalled].
# Bulk: five-cohort fragment-count DGE (adopted 774-unit release), log CPM, z within cohort.
# Primary family (Holm over 2), two-sided:
#   T1  all units: partial Spearman rho(THRB, TR4 score | dataset, sex, group, hepatocyte proxy),
#       compared with 1,000 random 4-gene sets matched gene by gene on mean log CPM decile
#       (TH genes excluded; seed 20261024). p = two-sided empirical rank in the null.
#   T2  disease units with a Kleiner stage (four cohorts): partial Spearman rho(TR4 score, fibrosis |
#       dataset, sex, hepatocyte proxy); t-approximation p.
# Descriptive (BH within each table): each TH gene vs fibrosis and NAS, without and with the stroma
# proxy; TR4 vs fibrosis adjusted also for THRB; THRB-TR4 coupling by cohort and by group; THRB
# partial rho with every gene (all units); protein in three proteomes.
# Usage: Rscript 42b_thrb_bulk_protein.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; EXT <- file.path(ROOT, "data/external")
TR4_GENES <- c("CPT1A", "ANGPTL4", "DIO1", "ME1")
MACH <- c("THRB", "THRA", "DIO1", "DIO2", "DIO3", "SLC10A1", "SLC16A2", "SLC16A10", "SLCO1C1")
EXTRA <- c("THRSP", "LDLR", "CYP7A1", "SHBG", "PPARGC1A", "PCSK9")  # recalled liver T3 or resmetirom-linked genes [Unverified]
TH <- unique(c(MACH, TR4_GENES, EXTRA))
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  vars <- vars[vapply(vars, function(v) length(unique(z[[v]])) > 1, logical(1))]
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
lc_all <- lc
lc <- lc[rowMeans(lc > 1) >= 0.5, ]
g <- rownames(lc)
fwrite(data.table(gene = TH, in_filtered = TH %in% g, mean_logcpm = ifelse(TH %in% rownames(lc_all), rowMeans(lc_all[intersect(TH, rownames(lc_all)), , drop = FALSE])[TH], NA)),
       file.path(out, "th_gene_abundance.tsv"), sep = "\t")
stopifnot(all(c("THRB", TR4_GENES) %in% g))
hep <- intersect(c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"), g)
str <- intersect(c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "TAGLN", "ACTA2", "MYL9"), g)
zc <- function(v) ave(v, m$dataset, FUN = function(a) (a - mean(a)) / sd(a))
Z <- t(apply(lc, 1, zc))
m[, hep_proxy := colMeans(Z[hep, , drop = FALSE])][, stroma_proxy := colMeans(Z[str, , drop = FALSE])]
thg <- intersect(TH, g)
for (s in thg) m[, (s) := Z[s, ]]
tr4 <- colMeans(Z[TR4_GENES, , drop = FALSE]); tr3 <- colMeans(Z[setdiff(TR4_GENES, "DIO1"), , drop = FALSE])
m[, TR4 := tr4][, TR4_noDIO1 := tr3]
m[, `:=`(NAS = as.numeric(nas_score), Fibrosis = as.numeric(fibrosis_stage), grp = factor(group_binary))]
ADJ1 <- c("dataset", "sex_final", "grp", "hep_proxy")
# T1 with matched null
mu <- rowMeans(lc); dec <- cut(rank(mu), 10, labels = FALSE); names(dec) <- g
pool <- split(setdiff(g, TH), dec[setdiff(g, TH)])
obs1 <- pcor(m, "THRB", "TR4", ADJ1)
NULLR <- vapply(seq_len(1000), function(i) {
  gs <- vapply(TR4_GENES, function(x) sample(pool[[as.character(dec[x])]], 1), "")
  m[, NS := colMeans(Z[gs, , drop = FALSE])]; pcor(m, "THRB", "NS", ADJ1)$rho }, numeric(1))
m[, NS := NULL]
p1 <- min(1, 2 * min((1 + sum(NULLR >= obs1$rho)) / 1001, (1 + sum(NULLR <= obs1$rho)) / 1001))
fwrite(data.table(null_rho = NULLR), file.path(out, "T1_null.tsv"), sep = "\t")
D <- m[grp == "Disease"]
DF <- D[!is.na(Fibrosis) & dataset %in% c("GSE130970", "GSE135251", "GSE162694", "GSE213621")]
DN <- D[!is.na(NAS) & dataset %in% c("GSE130970", "GSE135251", "GSE162694")]
ADJF <- c("dataset", "sex_final", "hep_proxy")
t2 <- pcor(DF, "TR4", "Fibrosis", ADJF)
# T2 matched-null context (same random sets idea, vs fibrosis)
NULLF <- vapply(seq_len(1000), function(i) {
  gs <- vapply(TR4_GENES, function(x) sample(pool[[as.character(dec[x])]], 1), "")
  DF[, NS := colMeans(Z[gs, match(DF$sample_id, m$sample_id), drop = FALSE])]; pcor(DF, "NS", "Fibrosis", ADJF)$rho }, numeric(1))
DF[, NS := NULL]
prim <- rbind(data.table(id = "T1", test = "rho(THRB, TR4 | dataset+sex+group+hep_proxy), all units", n = obs1$n, est = obs1$rho,
                         null_mean = mean(NULLR), null_q025 = quantile(NULLR, 0.025), null_q975 = quantile(NULLR, 0.975), p = p1),
              data.table(id = "T2", test = "rho(TR4, fibrosis | dataset+sex+hep_proxy), disease units", n = t2$n, est = t2$rho,
                         null_mean = mean(NULLF), null_q025 = quantile(NULLF, 0.025), null_q975 = quantile(NULLF, 0.975), p = t2$p))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
# descriptive
desc <- rbindlist(lapply(c(thg, "TR4", "TR4_noDIO1"), function(s) rbind(
  cbind(analysis = "fibrosis", pcor(DF, s, "Fibrosis", ADJF)),
  cbind(analysis = "fibrosis", pcor(DF, s, "Fibrosis", c(ADJF, "stroma_proxy"))),
  cbind(analysis = "fibrosis", pcor(DF, s, "Fibrosis", c("dataset", "sex_final"))),
  cbind(analysis = "NAS", pcor(DN, s, "NAS", ADJF)),
  cbind(analysis = "NAS", pcor(DN, s, "NAS", c(ADJF, "stroma_proxy"))),
  cbind(analysis = "NAS", pcor(DN, s, "NAS", c("dataset", "sex_final"))))))
desc[, bh := p.adjust(p, "BH"), by = .(analysis, adjust)]
fwrite(desc, file.path(out, "th_genes_vs_stage.tsv"), sep = "\t")
med <- rbind(cbind(analysis = "TR4 vs fibrosis", pcor(DF, "TR4", "Fibrosis", ADJF)),
             cbind(analysis = "TR4 vs fibrosis", pcor(DF, "TR4", "Fibrosis", c(ADJF, "THRB"))),
             cbind(analysis = "THRB vs fibrosis", pcor(DF, "THRB", "Fibrosis", ADJF)),
             cbind(analysis = "THRB vs fibrosis", pcor(DF, "THRB", "Fibrosis", c(ADJF, "TR4"))),
             cbind(analysis = "TR4_noDIO1 vs THRB", pcor(m, "THRB", "TR4_noDIO1", ADJ1)))
fwrite(med, file.path(out, "adjusted_for_each_other.tsv"), sep = "\t")
bycoh <- rbindlist(lapply(split(m, m$dataset), function(z) rbind(
  cbind(dataset = z$dataset[1], pcor(z, "THRB", "TR4", c("sex_final", "grp", "hep_proxy"))),
  rbindlist(lapply(TR4_GENES, function(s) cbind(dataset = z$dataset[1], pcor(z, "THRB", s, c("sex_final", "grp", "hep_proxy"))))))))
bycoh <- rbind(bycoh, rbindlist(lapply(levels(m$grp), function(k) cbind(dataset = paste("all", k), pcor(m[grp == k], "THRB", "TR4", c("dataset", "sex_final", "hep_proxy"))))))
fwrite(bycoh, file.path(out, "thrb_tr4_by_cohort.tsv"), sep = "\t")
# THRB partial rho with every gene, all units
X <- model.matrix(~ dataset + sex_final + grp + hep_proxy, data = m); q <- qr(X)
rt <- qr.resid(q, rank(m$THRB)); RX <- qr.resid(q, apply(lc, 1, rank))
r <- as.numeric(crossprod(RX, rt)) / sqrt(colSums(RX^2) * sum(rt^2))
gw <- data.table(gene = g, rho_with_THRB = r)[order(-rho_with_THRB)]
gw[, rank := seq_len(.N)]
fwrite(gw, file.path(out, "thrb_partial_rho_all_genes.tsv.gz"), sep = "\t")
# protein
to_gene <- function(x, gsym, maxmiss = 1) {
  ok <- !is.na(gsym) & gsym != "" & !grepl(";", gsym, fixed = TRUE)
  x <- x[ok, , drop = FALSE]; gsym <- gsym[ok]; x[!is.finite(x)] <- NA_real_
  keep <- rowMeans(is.na(x)) <= maxmiss; x <- x[keep, , drop = FALSE]; gsym <- gsym[keep]
  x <- normalizeBetweenArrays(x, method = "quantile")
  idx <- split(seq_len(nrow(x)), gsym)
  t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x)))) }
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt")); md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
smp <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% smp], by = "liver_proteomics_filename"); setorder(md, liver_proteomics_filename)
x1 <- as.matrix(raw[, md$liver_proteomics_filename, with = FALSE]); storage.mode(x1) <- "double"; x1[x1 <= 0] <- NA; x1 <- log2(x1)
G1 <- to_gene(x1, raw$Genes, 0.5); colnames(G1) <- md$liver_proteomics_filename
a1 <- md[, .(sample = liver_proteomics_filename, grp = fifelse(saf_diagnosis != "No_MASLD", "disease", "control"),
             batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")), age = as.numeric(alder),
             bmi = as.numeric(bmi), sex = factor(gender), Ballooning = as.numeric(hepatocellular_ballooning_score),
             Fibrosis = as.numeric(kleiner_fibrosis_grade), NAS = as.numeric(nafld_activity_score))]
pa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_proteome.csv"))
aa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_annotation_file.csv"), select = 1:3)
setnames(aa, c("patient", "sample", "group")); aa <- aa[sample %in% names(pa)]
x2 <- as.matrix(pa[, aa$sample, with = FALSE]); storage.mode(x2) <- "double"
G2 <- to_gene(x2, pa[["Gene name"]]); colnames(G2) <- aa$sample
a2 <- aa[, .(sample, grp = fcase(group == "control", "control", group == "NASH", "disease", group == "Cirrhosis", "cirrhosis"))]
pb <- fread(file.path(EXT, "PXD065892/liver_20261017T003606Z/liver_sn_minprob_processed.tsv"))
x3 <- t(as.matrix(pb[, -(1:3)])); storage.mode(x3) <- "double"
G3 <- to_gene(x3, sub("~.*$", "", rownames(x3))); colnames(G3) <- pb$sample
a3 <- data.table(sample = pb$sample, grp = fcase(pb$group == "Healthy", "control", pb$group == "NAFLD", "disease", pb$group == "AIH", "AIH"))
stopifnot(!anyNA(a2$grp), !anyNA(a3$grp))
P <- list(PXD051911 = list(G1, a1), PXD027722 = list(G2, a2), PXD065892 = list(G3, a3))
fwrite(rbindlist(lapply(names(P), function(k) data.table(cohort = k, gene = TH, detected = TH %in% rownames(P[[k]][[1]])))), file.path(out, "protein_detection.tsv"), sep = "\t")
pz <- function(G, a, genes) { gg <- intersect(genes, rownames(G)); d <- copy(a)
  for (s in gg) d[, (s) := as.numeric(G[s, sample])]; list(d = d, genes = gg) }
q1 <- pz(G1, a1, TH); M1 <- q1$d[grp == "disease"]
BASE1 <- c("batch", "age", "bmi", "sex")
prot <- rbindlist(lapply(q1$genes, function(s) rbindlist(lapply(c("Fibrosis", "Ballooning", "NAS"), function(b)
  cbind(analysis = paste("PXD051911 MASLD", b), pcor(M1, s, b, BASE1), n_obs = sum(!is.na(M1[[s]])))))))
if (nrow(prot)) prot[, bh := p.adjust(p, "BH"), by = analysis]
fwrite(prot, file.path(out, "protein_stage_per_gene.tsv"), sep = "\t")
grp <- rbindlist(lapply(list(c("PXD051911", "disease"), c("PXD027722", "disease"), c("PXD027722", "cirrhosis"), c("PXD065892", "disease"), c("PXD065892", "AIH")), function(k) {
  q <- pz(P[[k[1]]][[1]], P[[k[1]]][[2]], TH)
  rbindlist(lapply(q$genes, function(s) { a <- q$d[grp == k[2], get(s)]; b <- q$d[grp == "control", get(s)]; a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 3 || length(b) < 3) return(NULL)
    data.table(cohort = k[1], contrast = paste(k[2], "vs control"), gene = s, n1 = length(a), n0 = length(b),
               d = (mean(a) - mean(b)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2)),
               p = t.test(a, b)$p.value) })) }))
if (nrow(grp)) grp[, bh := p.adjust(p, "BH"), by = .(cohort, contrast)]
fwrite(grp, file.path(out, "protein_group_per_gene.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(prim, digits = 3)
print(dcast(desc, gene ~ analysis + adjust, value.var = "rho"), digits = 2)
print(med, digits = 3); print(bycoh[gene == "THRB" & b == "TR4"], digits = 3)
print(gw[gene %in% c(TR4_GENES, MACH, EXTRA)], digits = 3); print(head(gw, 40), digits = 3)
print(prot, digits = 3); print(grp, digits = 3)
