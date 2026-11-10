# DIG1 job 22 (exploratory; follows jobs 19 and 21). What kind of mitochondrial protein loss is the
# MASLD respiratory-chain deficit? PXD051911 liver protein, 46 MASLD + 12 controls, same matrix as job
# 21 (single-gene protein groups, <=50% missing, quantile normalised, median per gene). Scores: z
# within protein across the 58, mean over members; ratio scores are differences of two scores.
# Groups were fixed before this script read any value. NDUFA4 is moved from Complex I to Complex IV.
# Primary family (Holm over 3):
#   H1  MASLD vs control on CI - CII (Complex II is all nuclear-encoded and outside supercomplexes);
#       lm ~ masld + batch + age + BMI + sex.
#   H2  within MASLD, CI - hep_mito (CPS1, OTC, HMGCS2, GLS2: hepatocyte-only matrix enzymes) vs
#       ballooning; partial Spearman, batch + age + BMI + sex.
#   H3  within MASLD, ATF4-target score vs CI - hep_mito; same partial Spearman.
# Descriptive (BH over the whole table): every group's MASLD effect, overall and within each batch;
# within-MASLD partial Spearman of every group with CI - hep_mito, ballooning, steatosis and the
# lipid-droplet coat; INHBE protein vs CI and PLIN2; RNA means of each group from job 19's table
# (bulk logFC, GSE298719 and GSE256398 hepatocyte logFC). Usage: Rscript 22_mito_pattern.R <new_out_dir>
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
CI_N <- intersect(c("NDUFV1", "NDUFV2", "NDUFS1"), g)
SETS <- list(
  CI = CI_all, CI_Nmodule = CI_N, CI_other = setdiff(CI_all, CI_N),
  CII = intersect(c("SDHA", "SDHB", "SDHC", "SDHD"), g), CIII = g[grepl("^UQCR|^CYC1$", g)],
  CIV = union(g[grepl("^COX[4-8]", g)], intersect("NDUFA4", g)), CV = g[grepl("^ATP5", g)],
  mtDNA = g[grepl("^MT-", g)], CYCS = intersect("CYCS", g),
  CI_assembly = intersect(c(paste0("NDUFAF", 1:8), "TIMMDC1", "ACAD9", "ECSIT", "FOXRED1", "TMEM126B", "NUBPL", "DMAC1", "DMAC2"), g),
  mt_proteases = intersect(c("LONP1", "CLPP", "CLPX", "YME1L1", "AFG3L2", "SPG7", "OMA1"), g),
  mt_chaperones = intersect(c("HSPD1", "HSPE1", "HSPA9", "TRAP1"), g),
  mt_nucleoid = intersect(c("TFAM", "POLG", "SSBP1", "LRPPRC", "TWNK"), g),
  mt_ribosome = g[grepl("^MRP[LS][0-9]", g)],
  mito_mass = intersect(c("VDAC1", "VDAC2", "VDAC3", "SLC25A4", "SLC25A5", "SLC25A6", "TOMM20", "TOMM40"), g),
  hep_mito = intersect(c("CPS1", "OTC", "HMGCS2", "GLS2"), g),
  FAO = intersect(c("CPT2", "ACADM", "ACADL", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "HADH", "ETFA", "ETFB", "ETFDH"), g),
  perox_FAO = intersect(c("ACOX1", "EHHADH", "HSD17B4", "ACAA1", "SCP2", "CAT", "PEX11A"), g),
  ROS_oxidases = intersect(c("CYP2E1", "CYP4A11", "ACOX1"), g),
  mt_antioxidant = intersect(c("SOD2", "PRDX3", "TXN2"), g),
  NRF2_targets = intersect(c("NQO1", "HMOX1", "SRXN1", "TXNRD1", "GCLC", "GCLM"), g),
  ATF4_targets = intersect(c("PHGDH", "PSAT1", "PSPH", "SHMT2", "MTHFD2", "MTHFD1L", "ALDH18A1", "PYCR1", "ASNS"), g),
  mitophagy = intersect(c("BNIP3", "BNIP3L", "FUNDC1"), g),
  glycolysis = intersect(c("PKM", "LDHA", "PFKP", "ALDOA", "TPI1", "GAPDH", "ENO1", "PGK1"), g),
  MDB_keratins = intersect(c("KRT8", "KRT18", "SQSTM1"), g),
  LD_coat = intersect(c("PLIN2", "PLIN4", "ABHD5"), g),
  INHBE = intersect("INHBE", g))
SETS <- SETS[lengths(SETS) > 0]
fwrite(data.table(set = names(SETS), n = lengths(SETS), genes = sapply(SETS, paste, collapse = ",")), file.path(out, "sets_used.tsv"), sep = "\t")
score <- function(gs) { z <- t(scale(t(gene[gs, , drop = FALSE]))); colMeans(z, na.rm = TRUE) }
for (s in names(SETS)) md[, (s) := score(SETS[[s]])]
md[, `:=`(CI_minus_CII = CI - CII, CI_minus_hep_mito = CI - hep_mito, CI_minus_mass = CI - mito_mass,
          Nmod_minus_CIother = CI_Nmodule - CI_other, CI_minus_FAO = CI - FAO)]
Y <- c(names(SETS), "CI_minus_CII", "CI_minus_hep_mito", "CI_minus_mass", "Nmod_minus_CIother", "CI_minus_FAO")
base <- c("batch", "age", "bmi_num", "sex")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(a = a, b = b, n = nrow(z), est = r, se = NA_real_, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
beta <- function(z, y, vars) {
  z <- z[complete.cases(z[, c(y, "masld", vars), with = FALSE])]
  s <- summary(lm(reformulate(c("masld", vars), y), data = z))$coefficients["masld", ]
  data.table(a = y, b = "MASLD_vs_control", n = nrow(z), est = s[1], se = s[2], p = s[4]) }
M <- md[masld == 1]
prim <- rbind(cbind(id = "H1", beta(md, "CI_minus_CII", base)),
              cbind(id = "H2", pcor(M, "CI_minus_hep_mito", "Ballooning", base)),
              cbind(id = "H3", pcor(M, "ATF4_targets", "CI_minus_hep_mito", base)))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
grp <- rbindlist(lapply(Y, function(y) beta(md, y, base)))
bt <- rbindlist(lapply(levels(md$batch), function(b) rbindlist(lapply(Y, function(y)
  cbind(beta(md[batch == b], y, c("age", "bmi_num", "sex")), batch = b)))))
H <- c("CI_minus_hep_mito", "Ballooning", "Steatosis", "Fibrosis", "LD_coat")
win <- rbindlist(lapply(setdiff(Y, "CI_minus_hep_mito"), function(y) rbindlist(lapply(setdiff(H, y), function(h) pcor(M, y, h, base)))))
inh <- rbindlist(lapply(c("CI", "CI_minus_hep_mito", "LD_coat", "Steatosis", "Ballooning", "ATF4_targets"), function(h) pcor(M, "INHBE", h, base)))
desc <- rbind(cbind(kind = "masld_effect", grp), cbind(kind = "within_MASLD", win), cbind(kind = "INHBE_within_MASLD", inh))
desc[, bh := p.adjust(p, "BH")]
fwrite(desc, file.path(out, "descriptive_bh.tsv"), sep = "\t")
fwrite(bt, file.path(out, "masld_effect_by_batch.tsv"), sep = "\t")
# RNA means per group from job 19's gene table (bulk is a different set of people)
gd <- fread(file.path(RUN, "out_prot_rna_discord/gene_discordance_disease.tsv.gz"))
rna <- rbindlist(lapply(names(SETS), function(s) { z <- gd[gene %in% SETS[[s]]]
  data.table(set = s, n_in_rna_table = nrow(z), prot = mean(z$logFC_prot), bulk_rna = mean(z$logFC_rna),
             flex_hep = mean(z$flex_hep_logFC, na.rm = TRUE), v3_hep = mean(z$v3_hep_logFC, na.rm = TRUE), d_mean = mean(z$d)) }))
fwrite(rna, file.path(out, "rna_means_by_group.tsv"), sep = "\t")
fwrite(md[, c("sample_id", "batch", "masld", "Fibrosis", "Steatosis", "Ballooning", "Inflammation", "NAS", Y), with = FALSE], file.path(out, "participant_scores.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(prim, digits = 3); print(grp[order(p)], digits = 2); print(rna, digits = 2)
