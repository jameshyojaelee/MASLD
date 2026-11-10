# DIG1 job 31 (exploratory; follows job 30). The ER machinery that inserts and glycosylates membrane
# and secretory proteins (Sec61, TRAP/SSR, OST, signal peptidase, EMC, SRP receptor) loses protein
# without losing RNA with ballooning and in MASLD. "Transmembrane transporter complexes" show the
# same gap, so the pattern may be a general shift for integral membrane proteins (membrane loss in
# injured cells, or extraction). Is it specific to the ER insertion machinery?
# Annotation: TM = HPA protein class "Predicted membrane proteins"; compartments from GO CC (GOALL,
# org.Hs.eg.db as installed): ER membrane GO:0005789, ER lumen GO:0005788, mito inner membrane
# GO:0005743, mito matrix GO:0005759, plasma membrane GO:0005886, cytosol GO:0005829. Classes are
# exclusive, assigned in this order: ER_insertion (fixed list), OXPHOS_mitoribo, ER_TM_other,
# mitoIM_TM_other, PM_TM, other_TM, ER_lumen_nonTM, matrix_nonTM, cytosol_nonTM, other_nonTM.
# Part A, per gene (job 25 ballooning gap = protein rho - RNA rho; job 19 disease d). Contrast A vs B
# = coefficient of A in lm(y ~ A + abundance decile + missing fraction) over genes in A or B;
# abundance and missingness from PXD051911 (mean log2 after quantile normalisation, before the
# per-group NA drop). p by permuting labels within abundance deciles.
# Part B, PXD051911 within 46 MASLD: class scores (mean z) vs ballooning (partial Spearman; batch,
# age, BMI, sex) and the partial correlation matrix between class scores.
# Part C, replication of the specific score ER_insertion - ER_TM_other (and every class score, CCT,
# ER-lumen chaperones) for disease vs control in PXD027722 and PXD065892, plus cirrhosis and AIH.
# Primary family (Holm over 4):
#   M1  ballooning gap, ER_insertion vs ER_TM_other (is the ER effect beyond being an ER TM protein?)
#   M2  ballooning gap, all TM outside OXPHOS/mitoribo/ER_insertion vs all non-TM (general shift?)
#   R1  PXD027722 NASH vs control, ER_insertion - ER_TM_other (permutation, 20,000)
#   R2  PXD065892 NAFLD vs healthy, same
# Descriptive: every class vs cytosol_nonTM for gap and d; all class scores and contrasts (BH).
# Usage: Rscript 31_membrane_topology.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma); library(org.Hs.eg.db); library(AnnotationDbi)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; EXT <- file.path(ROOT, "data/external")
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
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
# annotation
hpa <- fread(file.path(EXT, "hpa/proteinatlas.tsv"), select = c("Gene", "Protein class"))
tm <- unique(hpa[grepl("Predicted membrane proteins", `Protein class`, fixed = TRUE), Gene])
GO <- c(ERM = "GO:0005789", ERL = "GO:0005788", MIM = "GO:0005743", MAT = "GO:0005759", PM = "GO:0005886", CYT = "GO:0005829")
go <- as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unname(GO), keytype = "GOALL", columns = "SYMBOL")))
go <- unique(go[!is.na(SYMBOL), .(go = GOALL, gene = SYMBOL)])
inGO <- function(k) unique(go[go == GO[[k]], gene])
ER_INS <- c("SEC61A1", "SEC61A2", "SEC61B", "SEC61G", "SEC62", "SEC63", "SSR1", "SSR2", "SSR3", "SSR4", "STT3A", "STT3B", "RPN1", "RPN2",
            "DDOST", "DAD1", "OSTC", "TMEM258", "MAGT1", "TUSC3", "KRTCAP2", "OST4", "SPCS1", "SPCS2", "SPCS3", "SEC11A", "SEC11C",
            paste0("EMC", 1:10), "MMGT1", "SRPRA", "SRPRB", "TRAM1", "GET1", "CAMLG")
CCT <- c("TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8")
ER_CHAP <- c("HSPA5", "HSP90B1", "CALR", "PDIA3", "PDIA4", "PDIA6", "P4HB", "DNAJB11", "HYOU1", "MANF", "CRELD2", "SDF2L1")
is_ox <- function(g) (grepl("^NDUF[ABSCV][0-9]", g) | g %in% c("SDHA", "SDHB", "SDHC", "SDHD") | grepl("^UQCR|^CYC1$", g) |
                      grepl("^COX[4-8]", g) | grepl("^ATP5", g) | grepl("^MRP[LS][0-9]", g))
classify <- function(g) {
  t <- g %in% tm
  fcase(g %in% ER_INS, "ER_insertion", is_ox(g), "OXPHOS_mitoribo",
        t & g %in% inGO("ERM"), "ER_TM_other", t & g %in% inGO("MIM"), "mitoIM_TM_other", t & g %in% inGO("PM"), "PM_TM", t, "other_TM",
        g %in% inGO("ERL"), "ER_lumen_nonTM", g %in% inGO("MAT"), "matrix_nonTM", g %in% inGO("CYT"), "cytosol_nonTM", default = "other_nonTM") }
TMCL <- c("ER_TM_other", "mitoIM_TM_other", "PM_TM", "other_TM"); NONTM <- c("ER_lumen_nonTM", "matrix_nonTM", "cytosol_nonTM", "other_nonTM")
# Part A
ab <- data.table(gene = rownames(G1), abundance = rowMeans(G1, na.rm = TRUE), missing = rowMeans(is.na(G1)))
bal <- unique(fread(file.path(RUN, "out_bulk_rna_ballooning/per_gene_rna_vs_protein_ballooning.tsv.gz"))[, .(gene, gap = rho_prot - rho_rna, rho_prot, rho_rna)], by = "gene")
dis <- unique(fread(file.path(RUN, "out_prot_rna_discord/gene_discordance_disease.tsv.gz"))[, .(gene, d, logFC_prot, logFC_rna)], by = "gene")
A <- merge(merge(ab, bal, by = "gene", all.x = TRUE), dis, by = "gene", all.x = TRUE)
A[, class := classify(gene)][, tm := gene %in% tm][, decile := cut(abundance, quantile(abundance, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
fwrite(A, file.path(out, "per_gene_class.tsv.gz"), sep = "\t")
contrast <- function(y, inA, inB, nperm) {
  d <- A[(inA | inB) & !is.na(A[[y]])]; lab <- (inA[(inA | inB) & !is.na(A[[y]])])
  fit1 <- function(l) coef(lm(d[[y]] ~ l + factor(d$decile) + d$missing))[2]
  obs <- fit1(lab); strata <- split(seq_len(nrow(d)), d$decile)
  nd <- replicate(nperm, { l <- lab; for (s in strata) l[s] <- lab[s][sample.int(length(s))]; fit1(l) })
  data.table(stat = y, nA = sum(lab), nB = sum(!lab), meanA = mean(d[[y]][lab]), meanB = mean(d[[y]][!lab]), coef = unname(obs),
             p = (1 + sum(abs(nd) >= abs(obs) - 1e-12)) / (nperm + 1)) }
M1 <- contrast("gap", A$class == "ER_insertion", A$class == "ER_TM_other", 20000)
M2 <- contrast("gap", A$class %in% TMCL, A$class %in% NONTM, 20000)
descA <- rbindlist(lapply(c("gap", "d"), function(y) rbindlist(c(
  lapply(setdiff(unique(A$class), "cytosol_nonTM"), function(k) cbind(A = k, B = "cytosol_nonTM", contrast(y, A$class == k, A$class == "cytosol_nonTM", 2000))),
  list(cbind(A = "ER_insertion", B = "ER_TM_other", contrast(y, A$class == "ER_insertion", A$class == "ER_TM_other", 2000)),
       cbind(A = "all_TM_other", B = "all_nonTM", contrast(y, A$class %in% TMCL, A$class %in% NONTM, 2000)))))))
descA[, bh := p.adjust(p, "BH")]
fwrite(descA, file.path(out, "per_gene_class_contrasts.tsv"), sep = "\t")
fwrite(A[, .(n = .N, n_gap = sum(!is.na(gap)), mean_rho_prot = mean(rho_prot, na.rm = TRUE), mean_rho_rna = mean(rho_rna, na.rm = TRUE),
             mean_gap = mean(gap, na.rm = TRUE), n_d = sum(!is.na(d)), mean_logFC_prot = mean(logFC_prot, na.rm = TRUE),
             mean_logFC_rna = mean(logFC_rna, na.rm = TRUE), mean_d = mean(d, na.rm = TRUE), mean_abundance = mean(abundance)), by = class][order(mean_gap)],
       file.path(out, "per_class_summary.tsv"), sep = "\t")
# sample scores
score <- function(G) { Z <- t(scale(t(G))); cl <- classify(rownames(G))
  S <- c(split(rownames(G), cl), list(CCT = intersect(CCT, rownames(G)), ER_chaperones = intersect(ER_CHAP, rownames(G)),
         CI = rownames(G)[grepl("^NDUF[ABSCV][0-9]", rownames(G)) & rownames(G) != "NDUFA4"], CII = intersect(c("SDHA", "SDHB", "SDHC", "SDHD"), rownames(G))))
  sc <- as.data.table(lapply(S, function(s) colMeans(Z[s, , drop = FALSE], na.rm = TRUE)))
  sc[, ER_ins_minus_ER_TM := ER_insertion - ER_TM_other][, CI_minus_CII := CI - CII]
  list(n = lengths(S), sc = sc) }
S1 <- score(G1); S2 <- score(G2); S3 <- score(G3)
fwrite(rbindlist(lapply(list(PXD051911 = S1, PXD027722 = S2, PXD065892 = S3), function(s) data.table(set = names(s$n), n = s$n)), idcol = "cohort"),
       file.path(out, "class_sizes.tsv"), sep = "\t")
D1 <- cbind(a1, S1$sc); D2 <- cbind(a2, S2$sc); D3 <- cbind(a3, S3$sc)
fwrite(D1, file.path(out, "scores_PXD051911.tsv"), sep = "\t"); fwrite(D2, file.path(out, "scores_PXD027722.tsv"), sep = "\t"); fwrite(D3, file.path(out, "scores_PXD065892.tsv"), sep = "\t")
SC <- Reduce(intersect, list(names(S1$sc), names(S2$sc), names(S3$sc)))
# Part B
M <- D1[grp == "disease"][complete.cases(D1[grp == "disease", .(batch, age, bmi, sex, Ballooning)])]
X <- model.matrix(~ batch + age + bmi + sex, data = M); df <- nrow(X) - ncol(X) - 1
RR <- qr.resid(qr(X), sapply(c(SC, "Ballooning"), function(s) rank(M[[s]])))
pb_bal <- data.table(score = SC, n = nrow(M), rho = as.vector(cor(RR[, SC], RR[, "Ballooning"])))
pb_bal[, p := 2 * pt(-abs(rho * sqrt(df / (1 - rho^2))), df)][, bh := p.adjust(p, "BH")]
fwrite(pb_bal, file.path(out, "PXD051911_scores_vs_ballooning.tsv"), sep = "\t")
cm <- cor(RR[, SC]); fwrite(as.data.table(cm, keep.rownames = "score"), file.path(out, "PXD051911_MASLD_partial_cor_matrix.tsv"), sep = "\t")
# Part C
perm <- function(y, lab, nperm) { ok <- !is.na(y); y <- y[ok]; lab <- lab[ok]; dd <- mean(y[lab]) - mean(y[!lab])
  nd <- replicate(nperm, { l <- sample(lab); mean(y[l]) - mean(y[!l]) })
  c(diff = dd, sd_all = sd(y), p = (1 + sum(abs(nd) >= abs(dd) - 1e-12)) / (nperm + 1), n1 = sum(lab), n0 = sum(!lab)) }
CON <- list(list("PXD051911", D1, "disease"), list("PXD027722", D2, "disease"), list("PXD027722", D2, "cirrhosis"),
            list("PXD065892", D3, "disease"), list("PXD065892", D3, "AIH"))
R12 <- rbindlist(lapply(CON[c(2, 4)], function(k) { d <- k[[2]][grp %in% c(k[[3]], "control")]
  data.table(cohort = k[[1]], contrast = paste(k[[3]], "vs control"), t(perm(d$ER_ins_minus_ER_TM, d$grp == k[[3]], 20000))) }))
prim <- rbind(data.table(id = c("M1", "M2"), what = c("gap: ER_insertion vs ER_TM_other", "gap: TM (non-OXPHOS, non-ER_insertion) vs non-TM"),
                         estimate = c(M1$coef, M2$coef), p = c(M1$p, M2$p), nA = c(M1$nA, M2$nA), nB = c(M1$nB, M2$nB)),
              data.table(id = c("R1", "R2"), what = paste(R12$cohort, R12$contrast, "ER_ins - ER_TM_other"), estimate = R12$diff, p = R12$p,
                         nA = R12$n1, nB = R12$n0))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
sec <- rbindlist(lapply(CON, function(k) { d <- k[[2]][grp %in% c(k[[3]], "control")]
  rbindlist(lapply(SC, function(s) data.table(cohort = k[[1]], contrast = paste(k[[3]], "vs control"), score = s, t(perm(d[[s]], d$grp == k[[3]], 5000))))) }))
sec[, sd_units := diff / sd_all][, bh := p.adjust(p, "BH")]
fwrite(sec, file.path(out, "class_score_contrasts.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200)
print(prim, digits = 3); print(fread(file.path(out, "per_class_summary.tsv")), digits = 2)
print(descA[order(stat, coef)], digits = 2); print(pb_bal[order(rho)], digits = 2)
print(round(cm[c("ER_insertion", "ER_TM_other", "PM_TM", "mitoIM_TM_other", "CI", "cytosol_nonTM"), c("ER_insertion", "ER_TM_other", "PM_TM", "mitoIM_TM_other", "CI", "cytosol_nonTM")], 2))
print(dcast(sec, score ~ cohort + contrast, value.var = "sd_units"), digits = 2)
print(sec[score %in% c("ER_insertion", "ER_TM_other", "ER_ins_minus_ER_TM", "PM_TM", "CCT", "ER_chaperones")], digits = 2)
