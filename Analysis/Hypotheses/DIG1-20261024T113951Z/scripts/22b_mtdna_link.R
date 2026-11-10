# DIG1 job 22b (exploratory; follows jobs 19 and 21). What kind of mitochondrial protein loss is the
# MASLD respiratory-chain deficit? PXD051911 liver protein, 46 MASLD + 12 controls, same matrix as job
# Follow-up to job 22: is the loss tied to mtDNA nucleoid proteins (TFAM, SSBP1), and is the
# ballooning association spread across every subunit? Same matrix and covariates as job 22.
# Descriptive: per-protein MASLD effect and within-MASLD partial Spearman with ballooning for every
# OXPHOS subunit, nucleoid protein and mitoribosome protein, compared with all other proteins;
# TFAM and SSBP1 vs each complex score within MASLD. BH over each table.
# Usage: Rscript 22b_mtdna_link.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename")
stopifnot(nrow(md) == 58L)
md[, `:=`(sample_id = liver_proteomics_filename, masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
          age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), NAS = as.numeric(nafld_activity_score),
          Steatosis = as.numeric(steatosis_score), Ballooning = as.numeric(hepatocellular_ballooning_score),
          Inflammation = as.numeric(lobular_inflammation_score))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- raw$Genes
x <- x[rowMeans(is.na(x)) <= .5, , drop = FALSE]
x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
gene <- t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(gene) <- md$sample_id
g <- rownames(gene)
CI_all <- setdiff(g[grepl("^NDUF[ABSCV][0-9]", g)], "NDUFA4")
GRP <- c(setNames(rep("CI", length(CI_all)), CI_all),
         setNames(rep("CII", 4), c("SDHA", "SDHB", "SDHC", "SDHD")),
         setNames(rep("CIII", 0), character()), setNames(rep("mtDNA", sum(grepl("^MT-", g))), g[grepl("^MT-", g)]))
GRP <- c(GRP, setNames(rep("CIII", sum(grepl("^UQCR|^CYC1$", g))), g[grepl("^UQCR|^CYC1$", g)]),
         setNames(rep("CIV", sum(grepl("^COX[4-8]", g)) + 1), c(g[grepl("^COX[4-8]", g)], "NDUFA4")),
         setNames(rep("CV", sum(grepl("^ATP5", g))), g[grepl("^ATP5", g)]),
         setNames(rep("mt_ribosome", sum(grepl("^MRP[LS][0-9]", g))), g[grepl("^MRP[LS][0-9]", g)]),
         c(TFAM = "nucleoid", SSBP1 = "nucleoid", POLG = "nucleoid", LRPPRC = "nucleoid"))
GRP <- GRP[names(GRP) %in% g]
base <- c("batch", "age", "bmi_num", "sex")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  c(rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df), n = nrow(z)) }
mm <- model.matrix(~ masld + batch + age + bmi_num + sex, data = md)
fit <- limma::eBayes(limma::lmFit(gene, mm))
eff <- data.table(gene = g, logFC = fit$coefficients[, "masld"], p = fit$p.value[, "masld"])
M <- md[masld == 1]; GM <- gene[, M$sample_id]
bal <- rbindlist(lapply(g, function(p) { z <- copy(M); z[, v := GM[p, ]]; r <- pcor(z, "v", "Ballooning", base)
  data.table(gene = p, rho_ballooning = r[1], p_ballooning = r[2], n = r[3]) }))
tab <- merge(eff, bal, by = "gene")[, group := fifelse(gene %in% names(GRP), GRP[gene], "other")]
tab[, `:=`(bh_masld = p.adjust(p, "BH"), bh_ballooning = p.adjust(p_ballooning, "BH"))]
fwrite(tab[order(group, p)], file.path(out, "per_protein_effects.tsv"), sep = "\t")
summ <- tab[, .(n = .N, mean_logFC = mean(logFC), median_logFC = median(logFC), frac_down = mean(logFC < 0),
                mean_rho_ballooning = mean(rho_ballooning, na.rm = TRUE), frac_rho_neg = mean(rho_ballooning < 0, na.rm = TRUE)), by = group]
oth <- tab[group == "other"]
summ[, wilcox_logFC_vs_other := sapply(group, function(k) if (k == "other") NA else wilcox.test(tab[group == k, logFC], oth$logFC)$p.value)]
summ[, wilcox_rho_vs_other := sapply(group, function(k) if (k == "other") NA else wilcox.test(tab[group == k, rho_ballooning], oth$rho_ballooning)$p.value)]
fwrite(summ, file.path(out, "group_summary.tsv"), sep = "\t")
score <- function(gs) { z <- t(scale(t(gene[gs, , drop = FALSE]))); colMeans(z, na.rm = TRUE) }
for (k in c("CI", "CII", "CIII", "CIV", "CV", "mtDNA", "mt_ribosome")) md[, (k) := score(names(GRP)[GRP == k])]
for (p in c("TFAM", "SSBP1")) md[, (p) := gene[p, sample_id]]
M <- md[masld == 1]
tf <- rbindlist(lapply(c("TFAM", "SSBP1"), function(a) rbindlist(lapply(c("CI", "CII", "CIII", "CIV", "CV", "mtDNA", "mt_ribosome", "Ballooning"), function(b) {
  r <- pcor(M, a, b, base); data.table(a = a, b = b, rho = r[1], p = r[2], n = r[3]) }))))
tf[, bh := p.adjust(p, "BH")]
fwrite(tf, file.path(out, "nucleoid_vs_complexes_within_MASLD.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(summ, digits = 2); print(tf, digits = 2); print(tab[group %in% c("mtDNA", "nucleoid", "CII")], digits = 2)
