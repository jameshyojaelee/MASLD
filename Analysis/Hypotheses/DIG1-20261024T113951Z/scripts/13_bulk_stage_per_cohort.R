# DIG1 job 13 (exploratory). Does the bulk stage ordering (AKR1B10, KRT8, CDKN1A up at F0; GNMT,
# MAT1A down from F1; SOX9, EPCAM up from F2) hold inside each cohort? Per cohort: limma-voom
# ~ stage + sex (sex dropped if single-sex), stage = control, F0..F4, levels with >= 5 units.
# Units are genetic-identity clusters, not verified donors. Nominal p; no family claim.
# Usage: Rscript 13_bulk_stage_per_cohort.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
B <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANNOT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
KEY <- c("AKR1B10", "KRT8", "CDKN1A", "GNMT", "MAT1A", "BHMT", "SOX9", "EPCAM", "KRT7", "HMGCS1", "HMGCR", "SQLE", "FGF7", "FGFR2")
an <- fread(ANNOT, select = c("ensembl_base", "gene_name"))
y0 <- readRDS(file.path(B, "merged_dge.rds")); m0 <- readRDS(file.path(B, "meta_matched.rds"))
m0 <- m0[match(colnames(y0), m0$sample_id)]; stopifnot(identical(m0$sample_id, colnames(y0)))
m0[, st := fifelse(group_binary == "Disease", paste0("F", fibrosis_stage), "ctl")]
sym <- an$gene_name[match(sub("\\..*$", "", rownames(y0)), an$ensembl_base)]
res <- rbindlist(lapply(split(seq_len(nrow(m0)), m0$dataset), function(ix) {
  m <- m0[ix]; y <- y0[, ix]; k <- m$st != "FNA" & !is.na(m$sex_final); m <- m[k]; y <- y[, k]
  lv <- c("ctl", "F0", "F1", "F2", "F3", "F4"); n <- table(factor(m$st, lv)); lv <- lv[n >= 5]
  if (!"ctl" %in% lv || length(lv) < 2) return(NULL)
  k <- m$st %in% lv; m <- m[k]; y <- y[, k]; m[, st := factor(st, lv)]
  X <- if (uniqueN(m$sex_final) > 1) model.matrix(~ st + sex_final, m) else model.matrix(~ st, m)
  y <- DGEList(y$counts); y <- y[filterByExpr(y, design = X), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  f <- eBayes(lmFit(voom(y, X), X)); s <- sym[match(rownames(f), rownames(y0))]
  rbindlist(lapply(grep("^st", colnames(X), value = TRUE), function(cf) {
    w <- which(s %in% KEY)
    data.table(dataset = m$dataset[1], stage = sub("^st", "", cf), n_stage = sum(m$st == sub("^st", "", cf)),
               n_ctl = sum(m$st == "ctl"), gene = s[w], logFC = f$coefficients[w, cf], p = f$p.value[w, cf]) }))
}))
fwrite(res, file.path(out, "per_cohort_stage_vs_control_key_genes.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(dcast(res, dataset + gene ~ stage, value.var = "logFC")[order(gene)], digits = 2)
