# DIG1 job 19 (exploratory). Which disease changes in liver protein are not explained by the RNA
# change, and the reverse? Two independent sets of participants: PXD051911 DIA-MS (46 MASLD, 12
# controls; batch + age + BMI + sex) and 774-unit five-cohort bulk RNA (genetic-identity clusters;
# cohort + inferred sex). Sets below were named after seeing the published protein pathway table
# (respiratory chain down, lipid droplet and ECM up) and before any RNA-protein comparison.
# Effects: protein limma (08 loading block), RNA limma-voom on the same merged DGE as jobs 10/11.
# Discordance per gene d = (P - a - b R) / sqrt(seP^2 + b^2 seR^2), with a, b from a robust (MASS::rlm)
# fit of protein logFC on RNA logFC over shared genes. Positive d = protein higher than RNA predicts.
# Primary family (Holm over 4): cameraPR on d for GO CC respiratory chain complex,
# collagen-containing ECM, blood microparticle, lipid droplet. Secondary: all GO CC sets with 15-500
# shared genes and all 50 Hallmarks (BH within family). Same for the fibrosis slope within disease
# (protein: 46 MASLD; RNA: disease units with a recorded stage). Single-cell hepatocyte logFC
# (GSE298719, GSE256398) reported per set as a descriptive check on bulk dilution.
# Usage: Rscript 19_protein_rna_discordance.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(MASS)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
GS <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/genesets")
PAPER_PROT <- file.path(ROOT, "figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/protein_de_adjusted.tsv")
PRIMARY <- c("GOCC_RESPIRATORY_CHAIN_COMPLEX", "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", "GOCC_BLOOD_MICROPARTICLE", "GOCC_LIPID_DROPLET")
SHOW <- c("PLIN1","PLIN2","PLIN3","PLIN4","PLIN5","CIDEB","CIDEC","HSD17B13","PNPLA3","COL1A1","COL1A2","COL3A1","COL6A1","COL6A3","LUM","DCN",
          "MT-CO1","MT-CO2","MT-ND1","MT-ND5","NDUFA9","NDUFS1","SDHA","UQCRC1","COX4I1","ATP5F1A","ALB","APOA1","HP","C3","IGKC","IGHG1",
          "GNMT","MAT1A","CYP2C19","AKR1B10","KRT8","KRT18","KRT7","KRT19")
rd <- function(f) { g <- strsplit(readLines(file.path(GS, f)), "\t"); setNames(lapply(g, function(x) unique(x[-(1:2)])), sapply(g, `[`, 1)) }
cc <- rd("go_cc.gmt"); hm <- rd("hallmark.gmt"); stopifnot(all(PRIMARY %in% names(cc)), length(hm) == 50)

# protein (loading block as 08_protein_links.R)
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename"); stopifnot(nrow(md) == 58L)
md[, `:=`(sample_id = liver_proteomics_filename, masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
          age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- raw$Genes
x <- x[rowMeans(is.na(x)) <= .5, , drop = FALSE]; x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
P <- t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(P) <- md$sample_id
pfit <- function(M, mm, coef) { f <- eBayes(lmFit(M, mm)); data.table(gene = rownames(M), logFC = f$coefficients[, coef], se = sqrt(f$s2.post) * f$stdev.unscaled[, coef], t = f$t[, coef], p = f$p.value[, coef]) }
ok <- complete.cases(md[, .(age, bmi_num, sex)])
Pd <- pfit(P[, ok], model.matrix(~ masld + batch + age + bmi_num + sex, data = md[ok]), "masld")
okf <- ok & md$masld == 1 & !is.na(md$Fibrosis)
Pf <- pfit(P[, okf], model.matrix(~ Fibrosis + batch + age + bmi_num + sex, data = md[okf]), "Fibrosis")
pp <- fread(PAPER_PROT)[, .(gene, paper_logFC = logFC)]
chk <- merge(Pd, pp, by = "gene")
fwrite(data.table(n_shared = nrow(chk), spearman_vs_paper = cor(chk$logFC, chk$paper_logFC, method = "spearman"),
                  n_masld = sum(md$masld[ok]), n_control = sum(1 - md$masld[ok]), n_fibrosis_model = sum(okf)), file.path(out, "protein_model_check.tsv"), sep = "\t")

# RNA
y <- readRDS(file.path(B, "merged_dge.rds")); mb <- readRDS(file.path(B, "meta_matched.rds"))
mb <- mb[match(colnames(y), mb$sample_id)]; stopifnot(identical(mb$sample_id, colnames(y)))
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
sb <- an$gene_name[match(sub("\\..*$", "", rownames(y)), an$ensembl_base)]
k <- !is.na(sb) & !duplicated(sb); y <- y[k, ]; rownames(y) <- sb[k]
k <- !is.na(mb$sex_final); y <- y[, k]; mb <- mb[k]
rfit <- function(yy, mm, coef) {
  yy <- yy[filterByExpr(yy, design = mm), , keep.lib.sizes = FALSE]; yy <- calcNormFactors(yy)
  f <- eBayes(lmFit(voom(yy, mm), mm))
  data.table(gene = rownames(f), logFC = f$coefficients[, coef], se = sqrt(f$s2.post) * f$stdev.unscaled[, coef], t = f$t[, coef], p = f$p.value[, coef]) }
mb[, disease := as.integer(group_binary == "Disease")]
Rd <- rfit(y, model.matrix(~ disease + dataset + sex_final, data = mb), "disease")
kf <- mb$disease == 1 & !is.na(mb$fibrosis_stage)
Rf <- rfit(y[, kf], model.matrix(~ fibrosis_stage + dataset + sex_final, data = mb[kf]), "fibrosis_stage")
units <- data.table(contrast = c("disease", "fibrosis_slope"), rna_units = c(ncol(y), sum(kf)), rna_disease_units = c(sum(mb$disease), sum(kf)))

# single-cell hepatocyte logFC
de <- unique(fread(file.path(RUN, "out_pseudobulk_v2/gene_DE_edgeR_all_lineages.tsv.gz")), by = c("lkey", "symbol"))
scF <- de[lkey == "GSE298719:hepatocytes", .(gene = symbol, flex_hep_logFC = logFC)]
scV <- de[lkey == "GSE256398:hepatocytes", .(gene = symbol, v3_hep_logFC = logFC)]

disc <- function(Pt, Rt, lab) {
  m <- merge(Pt, Rt, by = "gene", suffixes = c("_prot", "_rna"))
  fit <- rlm(logFC_prot ~ logFC_rna, data = m, maxit = 100)
  a <- coef(fit)[1]; b <- coef(fit)[2]
  m[, d := (logFC_prot - a - b * logFC_rna) / sqrt(se_prot^2 + b^2 * se_rna^2)]
  m[, `:=`(d_p = 2 * pnorm(-abs(d)))][, d_bh := p.adjust(d_p, "BH")]
  m <- merge(merge(m, scF, by = "gene", all.x = TRUE), scV, by = "gene", all.x = TRUE)
  sets <- function(S, fam) {
    S <- lapply(S, intersect, m$gene); S <- S[lengths(S) >= 15 & lengths(S) <= 500]
    r <- as.data.table(cameraPR(setNames(m$d, m$gene), ids2indices(S, m$gene)), keep.rownames = "set")
    r[, `:=`(family = fam, contrast = lab)]
    r[, c("mean_logFC_prot", "mean_logFC_rna", "mean_flex_hep", "mean_v3_hep", "frac_d_pos") := {
      g <- S[set]; list(sapply(g, function(z) mean(m[gene %in% z, logFC_prot])), sapply(g, function(z) mean(m[gene %in% z, logFC_rna])),
                        sapply(g, function(z) mean(m[gene %in% z, flex_hep_logFC], na.rm = TRUE)), sapply(g, function(z) mean(m[gene %in% z, v3_hep_logFC], na.rm = TRUE)),
                        sapply(g, function(z) mean(m[gene %in% z, d] > 0))) }]
    setnames(r, "FDR", "bh_q_within_family"); r }
  ccr <- sets(cc, "GO_CC"); hmr <- sets(hm, "Hallmark")
  prim <- ccr[set %in% PRIMARY][, p_holm_primary := p.adjust(PValue, "holm")]
  list(genes = m, fit = data.table(contrast = lab, n_genes = nrow(m), intercept = a, slope = b,
                                   spearman_logFC = cor(m$logFC_prot, m$logFC_rna, method = "spearman"),
                                   n_d_bh05_pos = sum(m$d_bh < .05 & m$d > 0), n_d_bh05_neg = sum(m$d_bh < .05 & m$d < 0)),
       sets = rbind(ccr, hmr), prim = prim)
}
res <- list(disease = disc(Pd, Rd, "disease"), fibrosis_slope = disc(Pf, Rf, "fibrosis_slope"))
fwrite(merge(rbindlist(lapply(res, `[[`, "fit")), units, by = "contrast"), file.path(out, "fit_summary.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "prim")), file.path(out, "primary_sets_holm.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "sets")), file.path(out, "all_sets_cameraPR.tsv"), sep = "\t")
for (nm in names(res)) fwrite(res[[nm]]$genes[order(-abs(d))], file.path(out, paste0("gene_discordance_", nm, ".tsv.gz")), sep = "\t")
fwrite(rbindlist(lapply(names(res), function(nm) cbind(contrast = nm, res[[nm]]$genes[gene %in% SHOW]))), file.path(out, "named_genes.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(fread(file.path(out, "protein_model_check.tsv"))); print(fread(file.path(out, "fit_summary.tsv")))
print(rbindlist(lapply(res, `[[`, "prim"))[, .(contrast, set, NGenes, Direction, PValue, p_holm_primary, mean_logFC_prot, mean_logFC_rna, mean_flex_hep)])
