# DIG1 job 44b (exploratory; input for 44c). Protein per colocalized gene with a genetic direction
# (41c table). Loaders and quantile normalisation as jobs 31/37a. Per gene: disease vs control in
# PXD051911, PXD027722 (NASH; cirrhosis separately), PXD065892 (NAFLD; AIH separately) as
# difference in means / pooled SD with Welch p; PXD051911 MASLD partial Spearman with fibrosis,
# ballooning and NAS | batch, age, BMI, sex. BH within each table. No test family here.
# Usage: Rscript 44b_coloc_genes_protein.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; EXT <- file.path(ROOT, "data/external")
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
GL <- unique(fread(file.path(RUN, "out_direction_posthoc_controls/gene_directions_susie_top_dedup.tsv"))$gene)
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  vars <- vars[vapply(vars, function(v) length(unique(z[[v]])) > 1, logical(1))]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(gene = a, b = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
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
             Fibrosis = as.numeric(gsub("[^0-9]", "", kleiner_fibrosis_grade)), NAS = as.numeric(nafld_activity_score))]
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
fwrite(rbindlist(lapply(names(P), function(k) data.table(cohort = k, gene = GL, detected = GL %in% rownames(P[[k]][[1]])))), file.path(out, "protein_detection.tsv"), sep = "\t")
pz <- function(G, a, genes) { gg <- intersect(genes, rownames(G)); d <- copy(a)
  for (s in gg) d[, (s) := as.numeric(G[s, sample])]; list(d = d, genes = gg) }
q1 <- pz(G1, a1, GL); M1 <- q1$d[grp == "disease"]
BASE1 <- c("batch", "age", "bmi", "sex")
prot <- rbindlist(lapply(q1$genes, function(s) rbindlist(lapply(c("Fibrosis", "Ballooning", "NAS"), function(b)
  cbind(analysis = paste("PXD051911 MASLD", b), pcor(M1, s, b, BASE1), n_obs = sum(!is.na(M1[[s]])))))))
if (nrow(prot)) prot[, bh := p.adjust(p, "BH"), by = analysis]
fwrite(prot, file.path(out, "protein_stage_per_gene.tsv"), sep = "\t")
grp <- rbindlist(lapply(list(c("PXD051911", "disease"), c("PXD027722", "disease"), c("PXD027722", "cirrhosis"), c("PXD065892", "disease"), c("PXD065892", "AIH")), function(k) {
  q <- pz(P[[k[1]]][[1]], P[[k[1]]][[2]], GL)
  rbindlist(lapply(q$genes, function(s) { a <- q$d[grp == k[2], get(s)]; b <- q$d[grp == "control", get(s)]; a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 3 || length(b) < 3) return(NULL)
    data.table(cohort = k[1], contrast = paste(k[2], "vs control"), gene = s, n1 = length(a), n0 = length(b),
               d = (mean(a) - mean(b)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2)),
               p = t.test(a, b)$p.value) })) }))
if (nrow(grp)) grp[, bh := p.adjust(p, "BH"), by = .(cohort, contrast)]
fwrite(grp, file.path(out, "protein_group_per_gene.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220); print(grp[order(p)][1:40], digits = 3); print(prot[order(p)][1:30], digits = 3)
