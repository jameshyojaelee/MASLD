# DIG1 job 20 (exploratory). Do hepatocyte endocrine changes seen in RNA appear in plasma protein of
# independent participants? Seen before this job, in both single-cell hepatocyte cohorts and bulk:
# lipogenic genes up (FASN, SCD, ACLY, ME1), insulin-induced hepatokines up (ANGPTL8, INHBE),
# insulin/FOXO1-repressed IGFBP1, IGFBP2, IRS2 down, and GH-STAT5 targets down (SOCS2, IGF1) with
# GHR not down. Plasma has not been opened for these proteins. PXD051911 plasma DIA-MS
# (batch-corrected deposit), initial samples with histology: participant is the unit.
# Primary family (Holm over 3), each MASLD vs No_MASLD, lm ~ masld + age + sex + BMI + plasma batch:
#   IGFBP2 (predicted lower), IGF1 (predicted lower), SHBG (predicted lower; hepatocyte-made and
#   repressed by hepatic lipogenesis, Selva et al. 2007 JCI). SHBG RNA was not inspected first.
# Secondary, descriptive: (S1) each protein vs steatosis, ballooning, inflammation, fibrosis, NAS
# (partial Spearman, same covariates); (S2) plasma-wide MASLD effects (BH over mapped proteins) and
# their agreement with hepatocyte RNA logFC for hepatocyte-enriched vs other proteins; (S3) 41
# participants with liver and plasma: liver lipogenic-enzyme protein score and liver GNMT vs plasma
# IGFBP2/SHBG/IGF1 (age, sex, BMI); (S4) 41 initial/follow-up pairs: change in plasma protein vs
# change in steatosis and NAS, adjusted for change in BMI. Usage: Rscript 20_plasma_endocrine.R <out>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PX <- file.path(ROOT, "data/PXD051911"); RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
PRIMARY <- c("IGFBP2", "IGF1", "SHBG")
SHOW <- c(PRIMARY, "IGFALS", "IGFBP3", "IGF2", "AHSG", "FETUB", "RBP4", "SELENOP", "GPLD1", "CLU", "APOA1", "APOB", "APOC3", "APOM", "CRP", "SAA1", "LRG1", "PON1", "SERPINA6", "FCN3", "PGLYRP2")
LIPO <- c("FASN", "ACLY", "ME1", "SCD", "ACACA", "ELOVL6", "THRSP", "PKLR")
# accession -> gene: liver deposit's own Genes column; IGFBP2 and IGF1 are absent there and added
# by UniProt accession (P18065 IGFBP2, P05019 IGF1; checked against rest.uniprot.org 2026-10-24)
lv <- fread(file.path(PX, "liver_protein_quant.txt"))
amap <- unique(lv[!is.na(Genes) & Genes != "", .(acc = unlist(strsplit(ProteinAccessions, ";"))), by = Genes], by = "acc")
amap <- rbind(amap, data.table(Genes = c("IGFBP2", "IGF1"), acc = c("P18065", "P05019")))
pl <- fread(file.path(PX, "plasma_protein_quant_filtered_and_batch_corrected.txt"))
pl[, gene := vapply(strsplit(ProteinAccessions, ";"), function(a) { g <- amap$Genes[match(a, amap$acc)]; g <- g[!is.na(g)]; if (length(g)) g[1] else NA_character_ }, "")]
pl <- pl[!is.na(gene) & !grepl(";", gene, fixed = TRUE)]
pl <- pl[, lapply(.SD, median, na.rm = TRUE), by = gene, .SDcols = setdiff(names(pl), c("ProteinAccessions", "gene"))]
X <- as.matrix(pl[, -1]); rownames(X) <- pl$gene
md <- fread(file.path(PX, "meta_data.txt"))[unique_identifier %in% colnames(X)]
md[, `:=`(masld = as.integer(saf_diagnosis != "No_MASLD"), age = as.numeric(alder), bmi = as.numeric(bmi), sex = factor(gender),
          batch = factor(plasma_batch_effects), Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)),
          steatosis = as.numeric(steatosis_score), ballooning = as.numeric(hepatocellular_ballooning_score),
          inflammation = as.numeric(lobular_inflammation_score), NAS = as.numeric(nafld_activity_score))]
I <- md[sample_group == "initial_sample" & complete.cases(md[, .(age, bmi, sex, batch)])]
fwrite(data.table(proteins_mapped = nrow(X), initial_n = nrow(I), initial_masld = sum(I$masld), initial_control = sum(1 - I$masld),
                  primary_present = paste(PRIMARY[PRIMARY %in% rownames(X)], collapse = ",")), file.path(out, "design.tsv"), sep = "\t")
stopifnot(all(PRIMARY %in% rownames(X)))
Xi <- X[, I$unique_identifier]
mm <- model.matrix(~ masld + age + sex + bmi + batch, data = I)
f <- eBayes(lmFit(Xi, mm))
tt <- data.table(gene = rownames(Xi), logFC = f$coefficients[, "masld"], t = f$t[, "masld"], p = f$p.value[, "masld"],
                 n_obs = rowSums(!is.na(Xi)))[, bh := p.adjust(p, "BH")]
prim <- tt[gene %in% PRIMARY][, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_masld_vs_control.tsv"), sep = "\t")
fwrite(tt[order(p)], file.path(out, "plasma_masld_vs_control_all.tsv"), sep = "\t")
# sensitivity: deposited uncorrected matrix (linear intensities, log2) with batch as covariate, as
# Analysis/Proteomics/scripts/pxd051911_plasma_de.R does
rw <- fread(file.path(PX, "plasma_protein_quant.txt"))[!is.na(Genes) & Genes != ""]
rw[, gene := tstrsplit(Genes, ";", fixed = TRUE)[[1]]]; rw <- rw[!duplicated(gene)]
Ir <- I[plasma_proteomics_filename %in% names(rw)]
R <- log2(as.matrix(rw[, Ir$plasma_proteomics_filename, with = FALSE])); rownames(R) <- rw$gene
R[!is.finite(R)] <- NA_real_
fr <- eBayes(lmFit(R, model.matrix(~ masld + age + sex + bmi + batch, data = Ir)))
fwrite(data.table(gene = rownames(R), logFC = fr$coefficients[, "masld"], t = fr$t[, "masld"], p = fr$p.value[, "masld"], n = nrow(Ir))[gene %in% SHOW],
       file.path(out, "sensitivity_raw_matrix_batch_covariate.tsv"), sep = "\t")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(a = a, b = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
# S1 histology within initial samples
J <- cbind(I, as.data.table(t(Xi[intersect(SHOW, rownames(Xi)), , drop = FALSE])))
s1 <- rbindlist(lapply(intersect(SHOW, rownames(Xi)), function(g) rbindlist(lapply(c("steatosis", "ballooning", "inflammation", "Fibrosis", "NAS"), function(h)
  pcor(J, g, h, c("age", "sex", "bmi", "batch"))))))
s1b <- rbindlist(lapply(PRIMARY, function(g) pcor(J, g, "bmi", c("age", "sex", "batch"))))
fwrite(rbind(s1, s1b), file.path(out, "S1_histology_partial_spearman.tsv"), sep = "\t")
# S2 agreement with hepatocyte RNA
V5 <- file.path(ROOT, "Analysis/Hypotheses/V5-20261021T122723Z/outputs/job3")
gm <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
spec_cpm <- function(s) { p <- fread(file.path(V5, sprintf("pseudobulk_%s.tsv.gz", s))); m <- as.matrix(p[, -1]); c <- t(t(m) / colSums(m)) * 1e6
  data.table(gene = gm$gene_name[match(sub("\\..*$", "", p$gene_id), gm$ensembl_base)], cpm = rowMeans(c)) }
sp <- Reduce(function(a, b) merge(a, b, by = "gene"), lapply(c("hep", "cholangiocytes", "fibroblasts", "macrophages", "tcells"), function(s) {
  d <- spec_cpm(s)[!is.na(gene)][, .(cpm = mean(cpm)), by = gene]; setnames(d, "cpm", s); d }))
sp[, hep_enriched := hep >= 2 * pmax(cholangiocytes, fibroblasts, macrophages, tcells)]
de <- unique(fread(file.path(RUN, "out_pseudobulk_v2/gene_DE_edgeR_all_lineages.tsv.gz")), by = c("lkey", "symbol"))
s2 <- Reduce(function(a, b) merge(a, b, by = "gene", all.x = TRUE), list(tt[, .(gene, plasma_logFC = logFC, plasma_t = t, plasma_bh = bh)],
  de[lkey == "GSE298719:hepatocytes", .(gene = symbol, flex_hep_logFC = logFC)], de[lkey == "GSE256398:hepatocytes", .(gene = symbol, v3_hep_logFC = logFC)],
  sp[, .(gene, hep_enriched)]))
fwrite(s2, file.path(out, "S2_plasma_vs_hepatocyte_rna.tsv"), sep = "\t")
fwrite(rbindlist(lapply(list(TRUE, FALSE), function(he) rbindlist(lapply(c("flex_hep_logFC", "v3_hep_logFC"), function(col) {
  z <- s2[hep_enriched == he & !is.na(get(col))]
  data.table(hep_enriched = he, rna = col, n = nrow(z), spearman = cor(z$plasma_t, z[[col]], method = "spearman"),
             sign_agree_plasma_bh05 = z[plasma_bh < .05, mean(sign(plasma_logFC) == sign(get(col)))], n_plasma_bh05 = z[plasma_bh < .05, .N]) })))),
  file.path(out, "S2_agreement.tsv"), sep = "\t")
# S3 liver protein vs plasma in the same participants
raw <- lv[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
lmd <- fread(file.path(PX, "meta_data.txt"))[liver_proteomics_filename %in% names(raw)]
lmd <- unique(lmd, by = "liver_proteomics_filename")
x <- as.matrix(raw[, lmd$liver_proteomics_filename, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- raw$Genes
x <- x[rowMeans(is.na(x)) <= .5, , drop = FALSE]; x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
L <- t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(L) <- lmd$patient_name
lipo <- intersect(LIPO, rownames(L)); zl <- t(scale(t(L[lipo, , drop = FALSE])))
lv_dt <- data.table(patient_name = colnames(L), liver_lipogenic = colMeans(zl, na.rm = TRUE), liver_GNMT = L["GNMT", ],
                    liver_SHBG = if ("SHBG" %in% rownames(L)) L["SHBG", ] else NA_real_)
P3 <- merge(lv_dt, cbind(I[, .(patient_name, age, sex, bmi, batch, masld)], as.data.table(t(Xi[PRIMARY, , drop = FALSE]))), by = "patient_name")
s3 <- rbindlist(lapply(c("liver_lipogenic", "liver_GNMT", "liver_SHBG"), function(a) rbindlist(lapply(PRIMARY, function(b)
  if (all(is.na(P3[[a]]))) NULL else rbind(pcor(P3, a, b, c("age", "sex", "bmi")), pcor(P3, a, b, c("age", "sex", "bmi", "masld")))))))
fwrite(cbind(s3, liver_lipogenic_proteins = paste(lipo, collapse = ",")), file.path(out, "S3_liver_vs_plasma.tsv"), sep = "\t")
# S4 within-person change
Fu <- md[sample_group == "follow_up_sample"]
pr <- merge(I[, .(patient_name, id0 = unique_identifier, bmi0 = bmi, st0 = steatosis, nas0 = NAS, fib0 = Fibrosis)],
            Fu[, .(patient_name, id1 = unique_identifier, bmi1 = bmi, st1 = steatosis, nas1 = NAS, fib1 = Fibrosis)], by = "patient_name")
D <- pr[, .(patient_name, d_bmi = bmi1 - bmi0, d_steatosis = st1 - st0, d_NAS = nas1 - nas0, d_fibrosis = fib1 - fib0)]
for (g in intersect(SHOW, rownames(X))) D[, (g) := X[g, pr$id1] - X[g, pr$id0]]
s4 <- rbindlist(lapply(intersect(SHOW, rownames(X)), function(g) {
  tt0 <- t.test(D[[g]]); rbind(data.table(a = g, b = "mean_change", adjust = "", n = sum(!is.na(D[[g]])), rho = tt0$estimate, p = tt0$p.value),
    pcor(D, g, "d_steatosis", "d_bmi"), pcor(D, g, "d_NAS", "d_bmi"), pcor(D, g, "d_bmi", "d_steatosis")) }))
fwrite(s4, file.path(out, "S4_within_person_change.tsv"), sep = "\t")
fwrite(D[, .(n = .N, mean_d_bmi = mean(d_bmi), mean_d_steatosis = mean(d_steatosis), mean_d_NAS = mean(d_NAS))], file.path(out, "S4_pairs_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(fread(file.path(out, "design.tsv"))); print(prim)
print(s1[a %in% PRIMARY]); print(fread(file.path(out, "S2_agreement.tsv"))); print(s3); print(s4[a %in% PRIMARY])
