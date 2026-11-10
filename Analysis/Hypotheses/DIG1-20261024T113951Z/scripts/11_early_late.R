# DIG1 job 11 (exploratory). Which hepatocyte changes are present before fibrosis (F0) and which
# appear only with fibrosis? Cross-sectional stage-associated contrasts, not progression.
# A. GSE256398 hepatocyte pseudobulk (17 donors): edgeR QL ~ group + sex, group = Healthy (6),
#    MASLD_F0 (3), F2+ (F2-F3 + cirrhosis, 8). Contrasts F0-H and F2+-H; camera Hallmarks per
#    contrast (BH over 50 within contrast).
# B. 774-unit bulk (units are genetic-identity clusters, not verified donors): limma-voom
#    ~ stage (control, F0..F4) + dataset + sex; units with missing stage dropped. cameraPR on each
#    stage-vs-control moderated t (BH over 50 within contrast).
# C. Agreement: Spearman of Hallmark signed -log10 p between A and B contrasts; gene-level
#    Spearman of V3 F0-H logFC vs bulk F0-control t over genes tested in both.
# Usage: Rscript 11_early_late.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
V3 <- file.path(ROOT, "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk")
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
GMT <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt")
KEY <- c("GNMT", "MAT1A", "BHMT", "CDKN1A", "SOX9", "EPCAM", "KRT8", "AKR1B10", "HMGCS1", "HMGCR", "SQLE", "LDLR", "PCSK9", "INSIG1", "SREBF2")
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
gl <- strsplit(readLines(GMT), "\t")
sets <- setNames(lapply(gl, function(x) unique(x[-(1:2)])), sapply(gl, `[`, 1))
stopifnot(length(sets) == 50)
to_sym <- function(ids) an$gene_name[match(sub("\\..*$", "", ids), an$ensembl_base)]

# A. V3 hepatocytes
m3 <- fread(file.path(V3, "donor_metadata.tsv"))
m3[, grp := factor(fifelse(stage == "Healthy", "H", fifelse(stage == "MASLD_F0", "F0", "F2p")), levels = c("H", "F0", "F2p"))]
pb <- fread(file.path(V3, "pseudobulk_hepatocyte.tsv.gz"))
cnt <- as.matrix(pb[, m3$donor, with = FALSE]); sym <- to_sym(pb$gene_id)
ok <- !is.na(sym) & !duplicated(sym); cnt <- cnt[ok, ]; rownames(cnt) <- sym[ok]
X <- model.matrix(~ 0 + grp + sex, data = m3)
y <- DGEList(cnt); y <- y[filterByExpr(y, design = X), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
y <- estimateDisp(y, X); fit <- glmQLFit(y, X, robust = TRUE)
C <- makeContrasts(F0_vs_H = grpF0 - grpH, F2p_vs_H = grpF2p - grpH, F2p_vs_F0 = grpF2p - grpF0, levels = X)
v <- voom(y, X); idx3 <- ids2indices(sets, rownames(v))
tabA <- rbindlist(lapply(colnames(C), function(k) {
  tt <- topTags(glmQLFTest(fit, contrast = C[, k]), n = Inf)$table
  data.table(contrast = k, gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM, p = tt$PValue, bh = tt$FDR) }))
fwrite(tabA, file.path(out, "v3_hep_de_by_group.tsv.gz"), sep = "\t")
camA <- rbindlist(lapply(colnames(C), function(k) {
  r <- as.data.table(camera(v, idx3, X, contrast = C[, k]), keep.rownames = "hallmark")
  r[, `:=`(contrast = k, signed = ifelse(Direction == "Up", 1, -1) * -log10(PValue))] }))
setnames(camA, "FDR", "bh_q_within_contrast")
fwrite(camA, file.path(out, "v3_hep_camera_by_group.tsv"), sep = "\t")

# B. bulk by stage
y2 <- readRDS(file.path(B, "merged_dge.rds")); mb <- readRDS(file.path(B, "meta_matched.rds"))
mb <- mb[match(colnames(y2), mb$sample_id)]; stopifnot(identical(mb$sample_id, colnames(y2)))
mb[, st := fifelse(group_binary == "Disease", paste0("F", fibrosis_stage), "ctl")]
keep <- mb$st != "FNA" & !is.na(mb$sex_final)
mb <- mb[keep]; y2 <- y2[, keep]
mb[, st := factor(st, levels = c("ctl", "F0", "F1", "F2", "F3", "F4"))]
fwrite(mb[, .N, by = .(dataset, st)][order(dataset, st)], file.path(out, "bulk_units_by_stage.tsv"), sep = "\t")
sb <- to_sym(rownames(y2)); ok <- !is.na(sb) & !duplicated(sb); y2 <- y2[ok, ]; rownames(y2) <- sb[ok]
Xb <- model.matrix(~ st + dataset + sex_final, data = mb)
y2 <- y2[filterByExpr(y2, design = Xb), , keep.lib.sizes = FALSE]; y2 <- calcNormFactors(y2)
vb <- voom(y2, Xb); fb <- eBayes(lmFit(vb, Xb))
stc <- grep("^st", colnames(Xb), value = TRUE)
tabB <- rbindlist(lapply(stc, function(k) data.table(contrast = sub("^st", "", k), gene = rownames(fb), logFC = fb$coefficients[, k], t = fb$t[, k], p = fb$p.value[, k])))
tabB[, bh := p.adjust(p, "BH"), by = contrast]
fwrite(tabB, file.path(out, "bulk_de_by_stage.tsv.gz"), sep = "\t")
camB <- rbindlist(lapply(stc, function(k) {
  r <- as.data.table(cameraPR(fb$t[, k], ids2indices(sets, rownames(fb))), keep.rownames = "hallmark")
  r[, `:=`(contrast = sub("^st", "", k), signed = ifelse(Direction == "Up", 1, -1) * -log10(PValue))] }))
setnames(camB, "FDR", "bh_q_within_contrast")
fwrite(camB, file.path(out, "bulk_cameraPR_by_stage.tsv"), sep = "\t")

# C. agreement
wA <- dcast(camA, hallmark ~ contrast, value.var = "signed"); setnames(wA, -1, paste0("v3_", names(wA)[-1]))
wB <- dcast(camB, hallmark ~ contrast, value.var = "signed"); setnames(wB, -1, paste0("bulk_", names(wB)[-1]))
w <- merge(wA, wB, by = "hallmark"); fwrite(w, file.path(out, "hallmark_signed_wide.tsv"), sep = "\t")
cols <- setdiff(names(w), "hallmark")
fwrite(as.data.table(cor(w[, ..cols], method = "spearman"), keep.rownames = "x"), file.path(out, "hallmark_agreement_spearman.tsv"), sep = "\t")
gA <- tabA[contrast == "F0_vs_H", .(gene, v3_F0_logFC = logFC)]
gAl <- tabA[contrast == "F2p_vs_H", .(gene, v3_F2p_logFC = logFC)]
gB <- dcast(tabB, gene ~ contrast, value.var = "t"); setnames(gB, -1, paste0("bulk_t_", names(gB)[-1]))
g <- Reduce(function(a, b) merge(a, b, by = "gene"), list(gA, gAl, gB))
fwrite(as.data.table(cor(g[, -1], method = "spearman"), keep.rownames = "x"), file.path(out, "gene_agreement_spearman.tsv"), sep = "\t")
fwrite(merge(dcast(tabA[gene %in% KEY], gene ~ contrast, value.var = c("logFC", "p")),
             dcast(tabB[gene %in% KEY], gene ~ contrast, value.var = c("logFC", "p")), by = "gene", all = TRUE,
             suffixes = c("_v3", "_bulk")), file.path(out, "key_genes.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(camA[bh_q_within_contrast < 0.05][order(contrast, PValue)])
print(camB[bh_q_within_contrast < 0.05 & hallmark %in% c("HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_P53_PATHWAY")][order(contrast)])
