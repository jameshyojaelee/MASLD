# DIG1 job 25 (exploratory; follows jobs 22, 22b). In liver protein (PXD051911) OXPHOS and
# mitoribosome proteins fall with ballooning within MASLD and Complex II is spared. Does bulk RNA show
# the same ballooning link? Five-cohort bulk (job 10 inputs; units are genetic-identity clusters).
# Ballooning is recorded only in GSE130970; NAS in GSE130970, GSE135251, GSE162694.
# Scores: z within cohort of log CPM, mean over genes; sets as job 22 (NDUFA4 in CIV).
# Primary family (Holm over 2), disease units only, partial Spearman, ranks residualised on covariates:
#   R1  GSE130970: CI - CII vs ballooning grade; sex + hepatocyte proxy.
#   R2  NAS cohorts pooled: CI - CII vs NAS; cohort + sex + hepatocyte proxy.
# Descriptive: every set vs ballooning (GSE130970) and NAS (pooled), with and without the hepatocyte
# proxy (BH over the table); per gene RNA rho vs ballooning in GSE130970 set beside the protein rho
# from job 22b for genes in both, with group summaries and the across-gene Spearman of the two.
# Usage: Rscript 25_bulk_rna_ballooning_oxphos.R <new_out_dir>
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
CI_all <- setdiff(g[grepl("^NDUF[ABSCV][0-9]", g)], "NDUFA4")
SETS <- list(
  CI = CI_all, CII = intersect(c("SDHA", "SDHB", "SDHC", "SDHD"), g), CIII = g[grepl("^UQCR|^CYC1$", g)],
  CIV = union(g[grepl("^COX[4-8]", g)], intersect("NDUFA4", g)), CV = g[grepl("^ATP5", g)],
  mtDNA = g[grepl("^MT-", g)], mt_ribosome = g[grepl("^MRP[LS][0-9]", g)],
  mt_nucleoid = intersect(c("TFAM", "POLG", "SSBP1", "LRPPRC", "TWNK"), g),
  mito_mass = intersect(c("VDAC1", "VDAC2", "VDAC3", "SLC25A4", "SLC25A5", "SLC25A6", "TOMM20", "TOMM40"), g),
  hep_mito = intersect(c("CPS1", "OTC", "HMGCS2", "GLS2"), g),
  FAO = intersect(c("CPT2", "ACADM", "ACADL", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "HADH", "ETFA", "ETFB", "ETFDH"), g),
  ATF4_targets = intersect(c("PHGDH", "PSAT1", "PSPH", "SHMT2", "MTHFD2", "MTHFD1L", "ALDH18A1", "PYCR1", "ASNS"), g),
  MDB_keratins = intersect(c("KRT8", "KRT18", "SQSTM1"), g), LD_coat = intersect(c("PLIN2", "PLIN4", "ABHD5"), g),
  hep_proxy = intersect(c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB"), g))
SETS <- SETS[lengths(SETS) > 0]
fwrite(data.table(set = names(SETS), n = lengths(SETS), genes = sapply(SETS, paste, collapse = ",")), file.path(out, "sets_used.tsv"), sep = "\t")
Z <- t(apply(lc, 1, function(v) ave(v, m$dataset, FUN = function(a) (a - mean(a)) / sd(a))))
for (s in names(SETS)) m[, (s) := colMeans(Z[SETS[[s]], , drop = FALSE], na.rm = TRUE)]
m[, CI_minus_CII := CI - CII]
m[, `:=`(Ballooning = as.numeric(ballooning_grade), NAS = as.numeric(nas_score))]
D <- m[group_binary == "Disease"]
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(a = a, b = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
G1 <- D[dataset == "GSE130970"]; GN <- D[!is.na(NAS) & dataset %in% c("GSE130970", "GSE135251", "GSE162694")]
fwrite(rbind(G1[, .N, by = .(dataset, Ballooning)], GN[, .(Ballooning = NA_real_, N = .N), by = dataset]), file.path(out, "design.tsv"), sep = "\t")
prim <- rbind(cbind(id = "R1", pcor(G1, "CI_minus_CII", "Ballooning", c("sex_final", "hep_proxy"))),
              cbind(id = "R2", pcor(GN, "CI_minus_CII", "NAS", c("dataset", "sex_final", "hep_proxy"))))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
Y <- c(setdiff(names(SETS), "hep_proxy"), "CI_minus_CII")
desc <- rbindlist(lapply(Y, function(s) rbind(
  pcor(G1, s, "Ballooning", c("sex_final", "hep_proxy")), pcor(G1, s, "Ballooning", "sex_final"),
  pcor(GN, s, "NAS", c("dataset", "sex_final", "hep_proxy")), pcor(GN, s, "NAS", c("dataset", "sex_final")))))
desc[, bh := p.adjust(p, "BH")]
fwrite(desc, file.path(out, "sets_vs_histology.tsv"), sep = "\t")
# per gene, GSE130970 disease units, RNA vs ballooning (sex + hepatocyte proxy), beside protein rho
G1c <- G1[complete.cases(G1[, .(Ballooning, sex_final, hep_proxy)])]
Zg <- Z[, match(G1c$sample_id, colnames(Z))]; stopifnot(identical(colnames(Zg), G1c$sample_id))
xx <- model.matrix(~ sex_final + hep_proxy, data = G1c)
rb <- qr.resid(qr(xx), rank(G1c$Ballooning))
RG <- qr.resid(qr(xx), apply(Zg, 1, rank))
rr <- as.vector(cor(rb, RG)); df <- nrow(xx) - ncol(xx) - 1
gr <- data.table(gene = rownames(Zg), rho_rna = rr, p_rna = 2 * pt(-abs(rr * sqrt(df / (1 - rr^2))), df))
pp <- fread(file.path(RUN, "out_mtdna_link/per_protein_effects.tsv"))[, .(gene, group, rho_prot = rho_ballooning, logFC_prot = logFC)]
gp <- merge(pp, gr, by = "gene")
fwrite(gp, file.path(out, "per_gene_rna_vs_protein_ballooning.tsv.gz"), sep = "\t")
gs <- gp[, .(n = .N, mean_rho_prot = mean(rho_prot), mean_rho_rna = mean(rho_rna), frac_rna_neg = mean(rho_rna < 0),
             mean_gap = mean(rho_prot - rho_rna)), by = group]
oth <- gp[group == "other"]
gs[, wilcox_rna_vs_other := sapply(group, function(k) if (k == "other") NA_real_ else wilcox.test(gp[group == k, rho_rna], oth$rho_rna)$p.value)]
fwrite(gs[order(mean_rho_prot)], file.path(out, "group_summary.tsv"), sep = "\t")
ct <- cor.test(gp$rho_prot, gp$rho_rna, method = "spearman", exact = FALSE)
fwrite(data.table(n_genes = nrow(gp), spearman = ct$estimate, p = ct$p.value), file.path(out, "across_gene_protein_vs_rna.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(prim, digits = 3); print(desc[adjust != "sex_final" & adjust != "dataset+sex_final"][order(p)], digits = 2); print(gs[order(mean_rho_prot)], digits = 2); print(ct)
