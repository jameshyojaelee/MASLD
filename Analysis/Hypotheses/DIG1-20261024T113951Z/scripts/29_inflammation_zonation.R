# DIG1 job 29 (exploratory; follows job 27). The Complex I loss with Complex II spared appears in
# autoimmune hepatitis (AIH) as strongly as in NAFLD (job 27), so it is not specific to fat. Does it
# track immune infiltration or inflammation, in AIH and in MASLD? And where in the lobule: does it go
# with loss of pericentral (zone 3, where ballooning sits) or periportal protein?
# Cohorts and processing as job 28 (PXD051911 single-gene groups <= 50% missing, quantile, median
# per gene; PXD027722 and PXD065892 deposited imputed matrices). Score = CI - CII (mean z).
# Immune sets (fixed before reading values): MHC_II, ISG (IFN-stimulated), immunoglobulin,
# leukocyte, macrophage; immune composite = mean of the five set scores. Zonation: pericentral
# (GLUL, CYP2E1, CYP1A2, OAT, CYP3A4, RHBG), periportal (HAL, SDS, ASS1, PCK1, HSD17B13, CYP2F1).
# Primary family (Holm over 2), Spearman, two-sided:
#   I1  AIH (PXD065892, n 19): CI - CII vs immune composite.
#   I2  PXD051911 MASLD (n 46): CI - CII vs lobular inflammation, ranks residualised on ballooning,
#       batch, age, BMI and sex.
# Descriptive (BH over the table): every set score vs CI - CII in every disease group; PXD051911
# CI - CII vs each histology grade with and without the other grades; zonation index vs CI - CII.
# Usage: Rscript 29_inflammation_zonation.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; EXT <- file.path(ROOT, "data/external")
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
             bmi = as.numeric(bmi), sex = factor(gender), Ballooning = as.numeric(hepatocellular_ballooning_score),
             Inflammation = as.numeric(lobular_inflammation_score), Steatosis = as.numeric(steatosis_score),
             Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)))]
# PXD027722
pa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_proteome.csv"))
aa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_annotation_file.csv"), select = 1:3)
setnames(aa, c("patient", "sample", "group")); aa <- aa[sample %in% names(pa)]
x2 <- as.matrix(pa[, aa$sample, with = FALSE]); storage.mode(x2) <- "double"
G2 <- to_gene(x2, pa[["Gene name"]]); colnames(G2) <- aa$sample
a2 <- aa[, .(sample, grp = fifelse(group == "control", "control", "disease"))]
# PXD065892
pb <- fread(file.path(EXT, "PXD065892/liver_20261017T003606Z/liver_sn_minprob_processed.tsv"))
x3 <- t(as.matrix(pb[, -(1:3)])); storage.mode(x3) <- "double"
G3 <- to_gene(x3, sub("~.*$", "", rownames(x3))); colnames(G3) <- pb$sample
a3 <- data.table(sample = pb$sample, grp = fcase(pb$group == "Healthy", "control", pb$group == "NAFLD", "disease", pb$group == "AIH", "AIH"))

SETS <- list(
  MHC_II = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74"),
  ISG = c("STAT1", "ISG15", "MX1", "IFIT3", "OAS1", "GBP1", "IFI35", "PSMB8", "PSMB9"),
  immunoglobulin = c("IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGHM", "IGKC", "JCHAIN"),
  leukocyte = c("PTPRC", "LCP1", "CORO1A", "LSP1", "ARHGDIB"),
  macrophage = c("CD163", "CD68", "LYZ", "MSR1", "MARCO", "VSIG4"),
  pericentral = c("GLUL", "CYP2E1", "CYP1A2", "OAT", "CYP3A4", "RHBG"),
  periportal = c("HAL", "SDS", "ASS1", "PCK1", "HSD17B13", "CYP2F1"),
  LD_coat = c("PLIN2", "PLIN4", "ABHD5"), MDB_keratins = c("KRT8", "KRT18", "SQSTM1"))
IMM <- c("MHC_II", "ISG", "immunoglobulin", "leukocyte", "macrophage")
scores <- function(G) { g <- rownames(G); Z <- t(scale(t(G)))
  CI <- setdiff(g[grepl("^NDUF[ABSCV][0-9]", g)], "NDUFA4"); CII <- intersect(c("SDHA", "SDHB", "SDHC", "SDHD"), g)
  sc <- function(m) colMeans(Z[intersect(m, g), , drop = FALSE], na.rm = TRUE)
  d <- data.table(sample = colnames(G), CI_minus_CII = sc(CI) - sc(CII))
  for (s in names(SETS)) d[, (s) := if (length(intersect(SETS[[s]], g))) sc(SETS[[s]]) else NA_real_]
  d[, immune := rowMeans(.SD), .SDcols = IMM][, zonation_pp_minus_pc := periportal - pericentral]
  used <- data.table(set = names(SETS), members = sapply(SETS, function(m) paste(intersect(m, g), collapse = ",")))
  list(d = d, used = used) }
S1 <- scores(G1); S2 <- scores(G2); S3 <- scores(G3)
fwrite(rbind(cbind(cohort = "PXD051911", S1$used), cbind(cohort = "PXD027722", S2$used), cbind(cohort = "PXD065892", S3$used)), file.path(out, "sets_used.tsv"), sep = "\t")
D1 <- merge(a1, S1$d, by = "sample"); D2 <- merge(a2, S2$d, by = "sample"); D3 <- merge(a3, S3$d, by = "sample")
pcor <- function(z, a, b, vars = character()) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- if (length(vars)) model.matrix(reformulate(vars), data = z) else matrix(1, nrow(z), 1)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(z) - ncol(xx) - 1
  data.table(x = a, y = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
base <- c("batch", "age", "bmi", "sex")
M1 <- D1[grp == "disease"]
prim <- rbind(cbind(id = "I1", group = "PXD065892_AIH", pcor(D3[grp == "AIH"], "CI_minus_CII", "immune")),
              cbind(id = "I2", group = "PXD051911_MASLD", pcor(M1, "CI_minus_CII", "Inflammation", c(base, "Ballooning"))))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
Y <- c(names(SETS), "immune", "zonation_pp_minus_pc")
GR <- list(PXD051911_MASLD = list(M1, base), PXD051911_control = list(D1[grp == "control"], character()),
           PXD027722_NASH_cirrhosis = list(D2[grp == "disease"], character()), PXD027722_control = list(D2[grp == "control"], character()),
           PXD065892_NAFLD = list(D3[grp == "disease"], character()), PXD065892_AIH = list(D3[grp == "AIH"], character()),
           PXD065892_healthy = list(D3[grp == "control"], character()))
desc <- rbindlist(lapply(names(GR), function(k) rbindlist(lapply(Y, function(y) {
  z <- GR[[k]][[1]]; if (all(is.na(z[[y]]))) return(NULL); cbind(group = k, pcor(z, "CI_minus_CII", y, GR[[k]][[2]])) }))))
H <- c("Ballooning", "Inflammation", "Steatosis", "Fibrosis")
hist <- rbindlist(lapply(H, function(h) rbind(cbind(group = "PXD051911_MASLD", pcor(M1, "CI_minus_CII", h, base)),
                                              cbind(group = "PXD051911_MASLD", pcor(M1, "CI_minus_CII", h, c(base, setdiff(H, h)))))))
hist2 <- rbindlist(lapply(c(IMM, "immune", "pericentral", "periportal", "zonation_pp_minus_pc"), function(y)
  rbindlist(lapply(H, function(h) cbind(group = "PXD051911_MASLD", pcor(M1, y, h, base))))))
desc <- rbind(desc, hist, hist2); desc[, bh := p.adjust(p, "BH")]
fwrite(desc, file.path(out, "descriptive_bh.tsv"), sep = "\t")
# group contrasts for immune and zonation scores (permutation, 20,000 shuffles)
perm <- function(y, lab, nperm = 20000) { ok <- !is.na(y); y <- y[ok]; lab <- lab[ok]; d <- mean(y[lab]) - mean(y[!lab])
  nd <- replicate(nperm, { l <- sample(lab); mean(y[l]) - mean(y[!l]) }); c(diff_sd = d / sd(y), p = (1 + sum(abs(nd) >= abs(d) - 1e-12)) / (nperm + 1)) }
con <- rbindlist(lapply(list(list("PXD051911", D1, "disease"), list("PXD027722", D2, "disease"), list("PXD065892", D3, "disease"), list("PXD065892", D3, "AIH")),
  function(k) { z <- k[[2]][grp %in% c(k[[3]], "control")]
    rbindlist(lapply(c("CI_minus_CII", Y), function(y) if (all(is.na(z[[y]]))) NULL else data.table(cohort = k[[1]], contrast = paste(k[[3]], "vs control"), score = y, t(perm(z[[y]], z$grp == k[[3]]))))) }))
con[, bh := p.adjust(p, "BH")]
fwrite(con, file.path(out, "group_contrasts.tsv"), sep = "\t")
fwrite(rbind(D1, D2, D3, fill = TRUE), file.path(out, "sample_scores.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(prim, digits = 3); print(desc[x == "CI_minus_CII" & adjust %in% c("", "batch+age+bmi+sex")][order(group, p)], digits = 2); print(hist, digits = 2); print(con, digits = 2)
