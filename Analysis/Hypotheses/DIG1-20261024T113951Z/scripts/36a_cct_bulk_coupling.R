# DIG1 job 36a (exploratory; follows jobs 33-35). (1) Does bulk RNA show the CCT rise with injury in
# more people than the 46-person proteome? (2) Within diseased liver, is CCT protein coupled to the
# Complex I deficit (one injury axis) or independent of it?
# Bulk: five-cohort fragment-count DGE as job 25 (log CPM, z within cohort, mean over set genes).
# Protein: job 31/33 sample scores (PXD051911, PXD027722, PXD065892).
# Sets (fixed): CCT TCP1 CCT2-8 CCT6A; tubulin; actin; MDB KRT8 KRT18 SQSTM1; YAP targets CCN1 CCN2
# ANKRD1 AMOTL2 AXL NUAK2 LATS2 [Unverified recall of commonly used YAP/TAZ targets]; HSF1 targets
# HSPA1A HSPA1B DNAJB1 HSPH1 BAG3 HSPA4L HSP90AA1; ductular/progenitor KRT7 KRT19 SOX9 EPCAM SPP1;
# CI, CII as job 25.
# Primary family (Holm over 2):
#   B1  NAS cohorts pooled (GSE130970, GSE135251, GSE162694; disease units): CCT RNA vs NAS, partial
#       Spearman | cohort + sex + hepatocyte proxy.
#   B2  PXD051911 MASLD (n 46): CCT protein vs CI - CII, partial Spearman | batch, age, BMI, sex.
# Descriptive (BH): every set vs ballooning (GSE130970) and NAS (pooled); within bulk disease units,
# CCT vs every set | cohort + sex + hepatocyte proxy; protein CCT vs CI - CII with ballooning added
# and in the other disease groups.
# Usage: Rscript 36a_cct_bulk_coupling.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
y <- readRDS(file.path(B, "merged_dge.rds")); m <- readRDS(file.path(B, "meta_matched.rds"))
m <- m[match(colnames(y), m$sample_id)]; stopifnot(identical(m$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
sym <- an$gene_name[match(sub("\\..*$", "", rownames(y)), an$ensembl_base)]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE, prior.count = 1)
ok <- !is.na(sym) & !duplicated(sym); lc <- lc[ok, ]; rownames(lc) <- sym[ok]
lc <- lc[rowMeans(lc > 1) >= 0.5, ]
g <- rownames(lc)
SETS <- list(CCT = c("TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8"),
  tubulin = c("TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBB", "TUBB2A", "TUBB4B", "TUBB6"), actin = c("ACTB", "ACTG1"),
  MDB = c("KRT8", "KRT18", "SQSTM1"), YAP = c("CCN1", "CCN2", "ANKRD1", "AMOTL2", "AXL", "NUAK2", "LATS2"),
  HSF1 = c("HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "BAG3", "HSPA4L", "HSP90AA1"), ductular = c("KRT7", "KRT19", "SOX9", "EPCAM", "SPP1"),
  CI = setdiff(g[grepl("^NDUF[ABSCV][0-9]", g)], "NDUFA4"), CII = c("SDHA", "SDHB", "SDHC", "SDHD"),
  hep_proxy = c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"))
SETS <- lapply(SETS, intersect, g); SETS <- SETS[lengths(SETS) > 0]
fwrite(data.table(set = names(SETS), n = lengths(SETS), genes = sapply(SETS, paste, collapse = ",")), file.path(out, "sets_used.tsv"), sep = "\t")
Z <- t(apply(lc, 1, function(v) ave(v, m$dataset, FUN = function(a) (a - mean(a)) / sd(a))))
for (s in names(SETS)) m[, (s) := colMeans(Z[SETS[[s]], , drop = FALSE], na.rm = TRUE)]
m[, CI_minus_CII := CI - CII][, `:=`(Ballooning = as.numeric(ballooning_grade), NAS = as.numeric(nas_score))]
D <- m[group_binary == "Disease"]
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(a = a, b = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
G1 <- D[dataset == "GSE130970"]; GN <- D[!is.na(NAS) & dataset %in% c("GSE130970", "GSE135251", "GSE162694")]
# protein scores
P <- lapply(c("PXD051911", "PXD027722", "PXD065892"), function(k) { a <- fread(file.path(RUN, "out_cct_composition", paste0("scores_", k, ".tsv")))
  b <- fread(file.path(RUN, "out_membrane_topology", paste0("scores_", k, ".tsv"))); stopifnot(identical(a$sample, b$sample))
  cbind(a, b[, .(CI, CII, CI_minus_CII)]) }); names(P) <- c("PXD051911", "PXD027722", "PXD065892")
P$PXD051911[, sex := factor(sex)][, batch := factor(batch)]
M1 <- P$PXD051911[grp == "disease"]; BASE1 <- c("batch", "age", "bmi", "sex")
prim <- rbind(cbind(id = "B1", pcor(GN, "CCT", "NAS", c("dataset", "sex_final", "hep_proxy"))),
              cbind(id = "B2", pcor(M1, "CCT", "CI_minus_CII", BASE1)))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
Y <- setdiff(c(names(SETS), "CI_minus_CII"), "hep_proxy")
desc <- rbindlist(c(lapply(Y, function(s) rbind(cbind(analysis = "GSE130970 ballooning", pcor(G1, s, "Ballooning", c("sex_final", "hep_proxy"))),
                                                cbind(analysis = "NAS pooled", pcor(GN, s, "NAS", c("dataset", "sex_final", "hep_proxy"))))),
                    lapply(setdiff(Y, "CCT"), function(s) cbind(analysis = "bulk disease units: CCT vs set", pcor(D, "CCT", s, c("dataset", "sex_final", "hep_proxy"))))))
desc[, bh := p.adjust(p, "BH")]
fwrite(desc, file.path(out, "bulk_descriptive.tsv"), sep = "\t")
pc <- rbind(cbind(group = "PXD051911 MASLD", adj = "cov", pcor(M1, "CCT", "CI_minus_CII", BASE1)),
            cbind(group = "PXD051911 MASLD", adj = "cov+ballooning", pcor(M1, "CCT", "CI_minus_CII", c(BASE1, "Ballooning"))),
            cbind(group = "PXD051911 MASLD", adj = "cov+composition", pcor(M1, "CCT", "CI_minus_CII", c(BASE1, "MF", "immune", "macrophage", "stroma", "ductular", "endothelial", "hepatocyte"))),
            cbind(group = "PXD051911 MASLD", adj = "cov: CCT vs CI", pcor(M1, "CCT", "CI", BASE1)),
            cbind(group = "PXD051911 MASLD", adj = "cov: CCT vs CII", pcor(M1, "CCT", "CII", BASE1)),
            cbind(group = "PXD027722 NASH+cirrhosis", adj = "none", pcor(P$PXD027722[grp != "control"], "CCT", "CI_minus_CII", character())),
            cbind(group = "PXD065892 NAFLD", adj = "none", pcor(P$PXD065892[grp == "disease"], "CCT", "CI_minus_CII", character())),
            cbind(group = "PXD065892 AIH", adj = "none", pcor(P$PXD065892[grp == "AIH"], "CCT", "CI_minus_CII", character())),
            cbind(group = "PXD051911 control", adj = "none", pcor(P$PXD051911[grp == "control"], "CCT", "CI_minus_CII", character())))
fwrite(pc, file.path(out, "protein_cct_vs_ci.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200); print(prim, digits = 3); print(desc[order(analysis, p)], digits = 2); print(pc, digits = 2)
