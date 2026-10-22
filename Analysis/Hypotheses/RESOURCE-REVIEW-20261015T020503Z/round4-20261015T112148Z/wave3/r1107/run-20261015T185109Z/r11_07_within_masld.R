#!/usr/bin/env Rscript
# R11-07: within-MASLD partial Spearman of GNMT, MAT1A, CYP2C19 protein with NAS and fibrosis.
# Estimator copied from Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R
# (Figure 5C partial_spearman): log2, rows <=50% missing, quantile normalization over all 58
# samples, symbol collapse by median, then within the participant set ranks of protein and
# score residualized on batch + age_z + bmi_z + sex; df = n - rank(design) - 1.
suppressPackageStartupMessages({library(data.table); library(limma)})
set.seed(20260923)  # nothing random
args <- commandArgs(TRUE); OUTD <- args[1]
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
raw <- fread(file.path(BASE, "data/PXD051911/liver_protein_quant.txt"))
meta <- fread(file.path(BASE, "data/PXD051911/meta_data.txt"))
sample_cols <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
meta <- meta[!is.na(liver_proteomics_filename) & liver_proteomics_filename != "" & liver_proteomics_filename %in% sample_cols]
meta <- unique(meta, by = "liver_proteomics_filename")
stopifnot(nrow(meta) == 58L)
meta[, sample_id := liver_proteomics_filename]
meta[, group := factor(fifelse(saf_diagnosis == "No_MASLD", "Control", "MASLD"), levels = c("Control", "MASLD"))]
meta[, acquisition_batch := factor(fifelse(grepl("^2019", sample_id), "2019", "2020"), levels = c("2019", "2020"))]
meta[, `:=`(age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
            Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), NAS = as.numeric(nafld_activity_score))]
meta[, `:=`(age_z = as.numeric(scale(age)), bmi_z = as.numeric(scale(bmi_num)))]
meta[, sample_order := match(sample_id, sample_cols)]; setorder(meta, sample_order)
samples <- meta$sample_id
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
expr <- as.matrix(raw[, ..samples]); storage.mode(expr) <- "double"
expr[!is.finite(expr) | expr <= 0] <- NA_real_; expr <- log2(expr); rownames(expr) <- raw$Genes
expr <- expr[rowMeans(is.na(expr)) <= 0.50, , drop = FALSE]
expr <- normalizeBetweenArrays(expr, method = "quantile")
genes <- c("GNMT", "MAT1A", "CYP2C19")
gexpr <- t(sapply(genes, function(g) { z <- expr[rownames(expr) == g, , drop = FALSE]
  if (nrow(z) == 1L) as.numeric(z[1, ]) else apply(z, 2, median, na.rm = TRUE) }))
colnames(gexpr) <- samples
cat("rows per gene:", sapply(genes, function(g) sum(rownames(expr) == g)), "\n")
nuis <- c("acquisition_batch", "age_z", "bmi_z", "sex")
ps <- function(g, feature, keep_rows) {
  protein <- gexpr[g, samples]; score <- meta[[feature]]
  ok <- keep_rows & is.finite(protein) & is.finite(score) & complete.cases(meta[, ..nuis])
  d <- droplevels(meta[ok]); rhs <- c("age_z", "bmi_z")
  if (nlevels(d$acquisition_batch) > 1L) rhs <- c("acquisition_batch", rhs)
  if (nlevels(d$sex) > 1L) rhs <- c(rhs, "sex")
  x <- model.matrix(reformulate(rhs), data = d); rk <- qr(x)$rank; df <- sum(ok) - rk - 1L
  r <- cor(residuals(lm.fit(x, rank(protein[ok], ties.method = "average"))),
           residuals(lm.fit(x, rank(score[ok], ties.method = "average"))))
  t <- r * sqrt(df / (1 - r^2))
  data.table(gene = g, feature, n = sum(ok), n_masld = sum(d$group == "MASLD"), n_control = sum(d$group == "Control"),
             n_nuisance_terms = rk - 1L, df, rho = r, pvalue = 2 * pt(-abs(t), df))
}
# Step 1: reproduce released 58-participant values
all58 <- rbindlist(lapply(genes, function(g) rbindlist(lapply(c("Fibrosis", "NAS"), ps, g = g, keep_rows = rep(TRUE, 58)))))
rel <- fread(file.path(BASE, "Analysis/Multimodal_Program_Projection/results/proteomics/panel4c_histology_partial.tsv"))
chk <- merge(all58, rel[, .(gene, feature, n_rel = n, rho_rel = rho, p_rel = pvalue)], by = c("gene", "feature"))
chk[, `:=`(abs_diff_rho = abs(rho - rho_rel), abs_diff_p = abs(pvalue - p_rel))]
print(chk); fwrite(chk, file.path(OUTD, "check_all58_vs_released.tsv"), sep = "\t")
stopifnot(all(chk$n == chk$n_rel), max(chk$abs_diff_rho) < 1e-10)
# Step 2: within MASLD
w <- rbindlist(lapply(genes, function(g) rbindlist(lapply(c("NAS", "Fibrosis"), ps, g = g, keep_rows = meta$group == "MASLD"))))
w[, padj_bh := p.adjust(pvalue, method = "BH")]
w[, correction_family := "3_proteins_x_2_scores_within_MASLD"]
w[, estimator := "Fig5C partial Spearman; quantile-normalized log2 abundance; residualized on batch, age, BMI, sex within MASLD"]
print(w); fwrite(w, file.path(OUTD, "r11_07_within_masld_protein_severity.tsv"), sep = "\t")
print(meta[group == "MASLD", .(n = .N, NAS_range = paste(range(NAS, na.rm = TRUE), collapse = "-"),
  Fib_range = paste(range(Fibrosis, na.rm = TRUE), collapse = "-"), nNA_NAS = sum(is.na(NAS)), nNA_Fib = sum(is.na(Fibrosis)))])
writeLines(capture.output(sessionInfo()), file.path(OUTD, "sessionInfo.txt"))
