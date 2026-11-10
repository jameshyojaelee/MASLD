# DIG1 job 33 (exploratory; follows jobs 30-32). CCT chaperonin (TRiC) protein is higher in MASLD in
# three proteomes and in AIH, rises with ballooning within MASLD, and exceeds its RNA (job 30 d +2.24).
# Its substrates tubulin and actin, keratins 8/18 and vimentin/filamin A/myosin IIA also rise with
# ballooning. CCT is abundant in immune and proliferating cells, so the rise may be a change in cell
# composition rather than in hepatocytes. Does CCT track ballooning (and disease) after adjusting for
# protein markers of the main non-parenchymal compartments?
# Composition scores (mean z of fixed protein markers; none is a CCT substrate set member):
#   immune PTPRC LCP1 CORO1A LSP1 CD74 HLA-DRA; macrophage CD68 CD163 MARCO LYZ;
#   stroma COL1A1 COL1A2 COL3A1 COL6A1 COL6A2 COL6A3 LUM DCN TAGLN ACTA2 MYL9;
#   ductular KRT7 KRT19 EPCAM ANXA4; endothelial PECAM1 VWF STAB1 STAB2 LYVE1 FCGR2B CLEC4G;
#   hepatocyte ALB APOA1 APOC3 SERPINA1 TF HP AHSG FGB. Plus MF and CF from job 31/32.
# Outcome sets: CCT (8), prefoldin, tubulin, actin, KRT8/18, VIM/FLNA/MYH9, proliferation (PCNA, MCM2/3/6).
# Primary family (Holm over 2):
#   C1  PXD051911 MASLD (n 46): CCT vs ballooning, partial Spearman | batch, age, BMI, sex, MF, CF,
#       six composition scores.
#   C2  PXD065892 AIH vs healthy: CCT | MF, CF, six composition scores; Freedman-Lane, 20,000.
# Descriptive (BH): every outcome under four adjustment levels for ballooning and for the five
# group contrasts; CCT vs each composition score within MASLD. Seed 20261024.
# Usage: Rscript 33_cct_composition.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
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
             bmi = as.numeric(bmi), sex = factor(gender), Ballooning = as.numeric(hepatocellular_ballooning_score),
             Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)))]
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
MK <- list(immune = c("PTPRC", "LCP1", "CORO1A", "LSP1", "CD74", "HLA-DRA"), macrophage = c("CD68", "CD163", "MARCO", "LYZ"),
  stroma = c("COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "LUM", "DCN", "TAGLN", "ACTA2", "MYL9"),
  ductular = c("KRT7", "KRT19", "EPCAM", "ANXA4"), endothelial = c("PECAM1", "VWF", "STAB1", "STAB2", "LYVE1", "FCGR2B", "CLEC4G"),
  hepatocyte = c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"))
OUTC <- list(CCT = c("TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8"), prefoldin = c("PFDN1", "PFDN2", "PFDN4", "PFDN5", "PFDN6", "VBP1"),
  tubulin = c("TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBB", "TUBB2A", "TUBB4B", "TUBB6"), actin = c("ACTB", "ACTG1"),
  KRT8_18 = c("KRT8", "KRT18"), VIM_FLNA_MYH9 = c("VIM", "FLNA", "MYH9"), proliferation = c("PCNA", "MCM2", "MCM3", "MCM6"))
COMP <- names(MK)
add_scores <- function(G, an, k) { Z <- t(scale(t(G))); S <- lapply(c(MK, OUTC), intersect, rownames(G))
  sc <- as.data.table(lapply(S, function(s) colMeans(Z[s, , drop = FALSE], na.rm = TRUE)))
  j <- fread(file.path(RUN, "out_membrane_topology", paste0("scores_", k, ".tsv")))
  stopifnot(identical(j$sample, an$sample))
  list(n = data.table(cohort = k, set = names(S), n = lengths(S), genes = sapply(S, paste, collapse = ",")),
       d = cbind(an, sc, j[, .(MF = (ER_TM_other + PM_TM + other_TM) / 3, CF = cytosol_nonTM)])) }
R1 <- add_scores(G1, a1, "PXD051911"); R2 <- add_scores(G2, a2, "PXD027722"); R3 <- add_scores(G3, a3, "PXD065892")
fwrite(rbind(R1$n, R2$n, R3$n), file.path(out, "sets_used.tsv"), sep = "\t")
D1 <- R1$d; D2 <- R2$d; D3 <- R3$d
for (k in list(list("PXD051911", D1), list("PXD027722", D2), list("PXD065892", D3))) fwrite(k[[2]], file.path(out, paste0("scores_", k[[1]], ".tsv")), sep = "\t")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(n = nrow(z), df = df, estimate = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
fl <- function(z, y, g, vars, nperm) {
  z <- z[complete.cases(z[, c(y, vars), with = FALSE])]; lab <- as.numeric(z$grp == g)
  X0 <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z); X1 <- cbind(lab, X0)
  b1 <- function(yy) qr.coef(qr(X1), yy)[1]
  f0 <- qr.fitted(qr(X0), z[[y]]); e0 <- z[[y]] - f0; obs <- b1(z[[y]])
  nd <- replicate(nperm, b1(f0 + e0[sample.int(length(e0))]))
  data.table(n = nrow(z), n1 = sum(lab), df = nrow(z) - ncol(X1), estimate = unname(obs), sd_y = sd(z[[y]]), p = (1 + sum(abs(nd) >= abs(obs) - 1e-12)) / (nperm + 1)) }
BASE1 <- c("batch", "age", "bmi", "sex")
ADJ <- list(none = character(), MF_CF = c("MF", "CF"), composition = c("MF", "CF", COMP), composition_fibrosis = c("MF", "CF", COMP, "Fibrosis"))
M1 <- D1[grp == "disease"]
prim <- rbind(cbind(id = "C1", what = "PXD051911 MASLD: CCT vs ballooning | cov + MF + CF + composition", pcor(M1, "CCT", "Ballooning", c(BASE1, ADJ$composition))),
              cbind(id = "C2", what = "PXD065892 AIH vs healthy: CCT | MF + CF + composition", fl(D3[grp %in% c("AIH", "control")], "CCT", "AIH", ADJ$composition, 20000)), fill = TRUE)
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
Y <- c(names(OUTC), COMP, "MF", "CF")
CON <- list(list("PXD051911", D1, "disease", BASE1), list("PXD027722", D2, "disease", character()), list("PXD027722", D2, "cirrhosis", character()),
            list("PXD065892", D3, "disease", character()), list("PXD065892", D3, "AIH", character()))
desc <- rbindlist(c(
  lapply(names(ADJ), function(a) rbindlist(lapply(setdiff(Y, ADJ[[a]]), function(y)
    cbind(analysis = "PXD051911 MASLD vs ballooning", y = y, adjust = a, pcor(M1, y, "Ballooning", c(BASE1, ADJ[[a]])))))),
  lapply(CON, function(k) rbindlist(lapply(names(ADJ)[1:3], function(a) rbindlist(lapply(setdiff(Y, ADJ[[a]]), function(y)
    cbind(analysis = paste(k[[1]], k[[3]], "vs control"), y = y, adjust = a, fl(k[[2]][grp %in% c(k[[3]], "control")], y, k[[3]], c(k[[4]], ADJ[[a]]), 5000)))))))), fill = TRUE)
desc[, bh := p.adjust(p, "BH")][, sd_units := estimate / sd_y]
fwrite(desc, file.path(out, "adjusted_contrasts.tsv"), sep = "\t")
cc <- rbindlist(lapply(setdiff(Y, "CCT"), function(y) cbind(y = y, pcor(M1, "CCT", y, BASE1))))
cc[, bh := p.adjust(p, "BH")]
fwrite(cc, file.path(out, "PXD051911_MASLD_CCT_partners.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(rbind(R1$n, R2$n, R3$n)[, .(cohort, set, n)][, dcast(.SD, set ~ cohort, value.var = "n")])
print(prim, digits = 3)
print(dcast(desc[analysis == "PXD051911 MASLD vs ballooning"], y ~ adjust, value.var = "estimate"), digits = 2)
print(dcast(desc[analysis == "PXD051911 MASLD vs ballooning"], y ~ adjust, value.var = "p"), digits = 2)
print(dcast(desc[analysis != "PXD051911 MASLD vs ballooning" & y %in% c("CCT", "tubulin", "KRT8_18", "immune", "proliferation", "ductular")], analysis + y ~ adjust, value.var = "sd_units"), digits = 2)
print(dcast(desc[analysis != "PXD051911 MASLD vs ballooning" & y %in% c("CCT", "tubulin", "KRT8_18", "immune", "proliferation", "ductular")], analysis + y ~ adjust, value.var = "p"), digits = 2)
print(cc[order(-estimate)], digits = 2)
