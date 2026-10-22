# Independent recomputation of the 5C histology column (partial Spearman) from the raw
# PXD051911 liver quant matrix. Own code: single-gene groups, log2, <=50% missing,
# limma quantile normalization, median collapse; rank-residual partial Spearman on
# batch + age + BMI + sex, all participants (58) and within MASLD (46).
suppressPackageStartupMessages({library(data.table); library(limma)})
set.seed(1)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
O <- file.path(B, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute")
raw <- fread(file.path(B, "data/PXD051911/liver_protein_quant.txt"))
meta <- fread(file.path(B, "data/PXD051911/meta_data.txt"))
rows <- fread(file.path(B, "figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/fixed_panel4c_rows.tsv"))
sc <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
meta <- unique(meta[!is.na(liver_proteomics_filename) & liver_proteomics_filename %in% sc], by = "liver_proteomics_filename")
stopifnot(nrow(meta) == 58)
meta[, batch := ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")]
meta[, group := ifelse(saf_diagnosis == "No_MASLD", "Control", "MASLD")]
meta[, `:=`(age = as.numeric(alder), bmi_n = as.numeric(bmi), Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)),
            Steatosis = as.numeric(steatosis_score), Ballooning = as.numeric(hepatocellular_ballooning_score),
            Inflammation = as.numeric(lobular_inflammation_score), NAS = as.numeric(nafld_activity_score))]
smp <- meta$liver_proteomics_filename
r <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
e <- as.matrix(r[, ..smp]); storage.mode(e) <- "double"
e[!is.finite(e) | e <= 0] <- NA; e <- log2(e); rownames(e) <- r$Genes
e <- e[rowMeans(is.na(e)) <= 0.5, , drop = FALSE]
e <- normalizeBetweenArrays(e, method = "quantile")
ge <- t(sapply(split(seq_len(nrow(e)), rownames(e)), function(i) if (length(i) == 1) e[i, ] else apply(e[i, , drop = FALSE], 2, median, na.rm = TRUE)))
colnames(ge) <- smp
cat("proteins after filter (gene level):", nrow(ge), "\n")
feats <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
ps <- function(g, f, keep) {
  x <- ge[g, smp]; y <- meta[[f]]
  ok <- is.finite(x) & is.finite(y) & keep
  d <- meta[ok]
  rhs <- c("scale_age", "scale_bmi")
  d[, scale_age := (age - mean(meta$age)) / sd(meta$age)]
  d[, scale_bmi := (bmi_n - mean(meta$bmi_n)) / sd(meta$bmi_n)]
  if (length(unique(d$batch)) > 1) rhs <- c("batch", rhs)
  if (length(unique(d$gender)) > 1) rhs <- c(rhs, "gender")
  X <- model.matrix(reformulate(rhs), data = d)
  rx <- residuals(lm.fit(X, rank(x[ok]))); ry <- residuals(lm.fit(X, rank(y[ok])))
  rho <- cor(rx, ry); df <- sum(ok) - qr(X)$rank - 1
  t <- rho * sqrt(df / (1 - rho^2)); p <- 2 * pt(-abs(t), df)
  data.table(gene = g, feature = f, n = sum(ok), rho = rho, p = p, df = df)
}
allp <- rbindlist(lapply(rows$gene, function(g) rbindlist(lapply(feats, function(f) ps(g, f, rep(TRUE, 58))))))
allp[, padj_bh125 := p.adjust(p, "BH")]
mas <- rbindlist(lapply(rows$gene, function(g) rbindlist(lapply(feats, function(f) ps(g, f, meta$group == "MASLD")))))
mas[, padj_bh125 := p.adjust(p, "BH")]
fwrite(allp, file.path(O, "5c_recomputed_partial_spearman_raw.tsv"), sep = "\t")
fwrite(mas, file.path(O, "5c_recomputed_partial_spearman_raw_masld_only.tsv"), sep = "\t")
ana <- fread(file.path(B, "figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/panel4c_histology_partial.tsv"))
m <- merge(allp, ana[, .(gene, feature, rho_a = rho, n_a = n, p_a = pvalue)], by = c("gene", "feature"))
cat("max |rho - analysis output|:", max(abs(m$rho - m$rho_a)), " n mismatches:", sum(m$n != m$n_a), " max |p diff|:", max(abs(m$p - m$p_a)), "\n")
mm <- merge(allp, mas, by = c("gene", "feature"), suffixes = c("", "_masld"))
cat("q<0.05 all participants:", sum(allp$padj_bh125 < 0.05), " within MASLD:", sum(mas$padj_bh125 < 0.05), "\n")
cat("sign flips all vs MASLD-only:", sum(sign(mm$rho) != sign(mm$rho_masld)), " median |rho| all:", median(abs(mm$rho)), " MASLD-only:", median(abs(mm$rho_masld)), "\n")
# heatmap matrix (display): removeBatchEffect retaining group
nuis <- cbind(age_z = as.numeric(scale(meta$age)), bmi_z = as.numeric(scale(meta$bmi_n)), sex_male = as.numeric(meta$gender == "Male"))
adj <- removeBatchEffect(ge[, smp], batch = meta$batch, covariates = nuis, design = model.matrix(~ group, data = meta))
adj <- adj[rows$gene, ]
z <- t(scale(t(adj)))
ord <- order(meta$Fibrosis, meta$NAS)
fwrite(data.table(gene = rows$gene, z[, ord]), file.path(O, "5c_recomputed_heatmap_z_from_raw.tsv"), sep = "\t")
fwrite(meta[ord, .(liver_proteomics_filename, group, batch, Fibrosis, NAS, Steatosis, Ballooning, Inflammation)], file.path(O, "5c_recomputed_heatmap_columns.tsv"), sep = "\t")
cat("group by column position:", paste(substr(meta$group[ord], 1, 1), collapse = ""), "\n")
cat("fibrosis by position:", paste(meta$Fibrosis[ord], collapse = ""), "\n")
cat("NAS by position:", paste(meta$NAS[ord], collapse = ","), "\n")
writeLines(capture.output(sessionInfo()), file.path(O, "sessionInfo_5c.txt"))
