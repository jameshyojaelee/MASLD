# DIG1 job 21 (exploratory; follows job 19). Within the 46 MASLD participants of PXD051911, which
# histological feature tracks (a) lipid-droplet coat proteins (PLIN2, PLIN4, ABHD5) and (b) the
# respiratory-chain protein loss, relative to mitochondrial mass proteins? Scores are fixed from
# gene-name patterns used in job 19 (z within protein, mean over members). Partial Spearman adjusted
# for batch, age, BMI and sex. Descriptive; BH over every test in the table.
# Usage: Rscript 21_protein_histology.R <new_out_dir>
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
md[, `:=`(Steatosis = as.numeric(steatosis_score), Ballooning = as.numeric(hepatocellular_ballooning_score),
          Inflammation = as.numeric(lobular_inflammation_score))]
g <- rownames(gene)
SETS <- list(CI = g[grepl("^NDUF[ABSCV][0-9]", g)], CIII = g[grepl("^UQCR|^CYC1$", g)], CIV = g[grepl("^COX[4-8]", g)],
             CV = g[grepl("^ATP5", g)], mtDNA = g[grepl("^MT-", g)],
             mito_mass = intersect(c("VDAC1","VDAC2","VDAC3","SLC25A4","SLC25A5","SLC25A6","TOMM20","TOMM40"), g),
             FAO = intersect(c("CPT2","ACADM","ACADL","ACADVL","ACADS","HADHA","HADHB","ECHS1","HADH","ETFA","ETFB","ETFDH"), g),
             perox_FAO = intersect(c("ACOX1","EHHADH","HSD17B4","ACAA1","SCP2","CAT","PEX11A"), g),
             LD_coat = intersect(c("PLIN2","PLIN4","ABHD5"), g), collagen = g[grepl("^COL[0-9]+A[0-9]$", g)])
score <- function(gs) { z <- t(scale(t(gene[gs, , drop = FALSE]))); colMeans(z, na.rm = TRUE) }
for (s in names(SETS)) md[, (s) := score(SETS[[s]])]
md[, CI_minus_mass := CI - mito_mass]
for (p in c("PLIN2", "PLIN4", "ABHD5", "GNMT")) md[, (p) := gene[p, sample_id]]
fwrite(data.table(set = names(SETS), n = lengths(SETS), genes = sapply(SETS, paste, collapse = ",")), file.path(out, "sets_used.tsv"), sep = "\t")
base <- c("batch", "age", "bmi_num", "sex")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(a = a, b = b, n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
M <- md[masld == 1]
Y <- c(names(SETS), "CI_minus_mass", "PLIN2", "PLIN4", "ABHD5", "GNMT")
H <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
hist <- rbindlist(lapply(Y, function(y) rbindlist(lapply(H, function(h) pcor(M, y, h, base)))))
pairs <- rbindlist(combn(c("CI", "CI_minus_mass", "mito_mass", "perox_FAO", "LD_coat", "GNMT", "collagen"), 2, function(v) pcor(M, v[1], v[2], base), simplify = FALSE))
grp <- rbindlist(lapply(Y, function(y) { z <- md[complete.cases(md[, c(y, "masld", base), with = FALSE])]
  s <- summary(lm(reformulate(c("masld", base), y), data = z))$coefficients["masld", ]; data.table(a = y, b = "MASLD_vs_control", n = nrow(z), rho = s[1], p = s[4]) }))
res <- rbind(cbind(kind = "histology_within_MASLD", hist), cbind(kind = "pairs_within_MASLD", pairs), cbind(kind = "masld_beta_score_or_log2", grp))
res[, bh := p.adjust(p, "BH")]
fwrite(res, file.path(out, "protein_scores_histology.tsv"), sep = "\t")
fwrite(md[, c("sample_id", "masld", "Fibrosis", "Steatosis", "Ballooning", "Inflammation", "NAS", Y), with = FALSE], file.path(out, "participant_scores.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(dcast(hist, a ~ b, value.var = "rho"), digits = 2); print(pairs, digits = 2); print(grp, digits = 2)
