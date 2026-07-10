#!/usr/bin/env Rscript
# 05j — DE sensitivity on the exact samples plotted in Fig S3U compact PCA.
#
# This is not a new canonical DEG call. It asks whether the definitive
# Control/Disease points displayed in figs3u_pca_definitive_raw_nas_fib.pdf
# produce a DEG set similar to the canonical limma-voom QW C2 result.

suppressMessages({
  library(edgeR)
  library(limma)
  library(ashr)
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
SDIR <- file.path(RDIR, "sensitivity")
dir.create(SDIR, recursive = TRUE, showWarnings = FALSE)

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
mm <- readRDS(file.path(RDIR, "meta_matched.rds"))
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  select = c("sample_id", "condition", "diagnosis_harmonized", "fibrosis_stage", "nas_score"))

keep <- as.character(dge$samples$dataset) %in% MEGA
dge <- dge[, keep]
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
samp[meta, on = "sample_id",
     `:=`(condition = i.condition, diagnosis_harmonized = i.diagnosis_harmonized,
          fibrosis_stage = i.fibrosis_stage, nas_score = i.nas_score)]

nas <- suppressWarnings(as.numeric(samp$nas_score))
fib <- suppressWarnings(as.numeric(samp$fibrosis_stage))
cond <- samp$condition
dsid <- samp$dataset

def_disease <- (!is.na(nas) & nas >= 5) | (!is.na(fib) & fib >= 3) |
  (dsid == "GSE126848" & cond == "NASH")
def_control <- !def_disease &
  ((!is.na(nas) & nas == 0 & (is.na(fib) | fib == 0)) |
     (dsid == "GSE126848" & cond == "Control") |
     (dsid == "GSE213621" & cond == "Control"))
samp[, def_class := fifelse(def_disease, "Disease",
                            fifelse(def_control, "Control", "Intermediate"))]

# Match figS_pca_definitive.R's compact panel table: controls plus disease
# samples with NAS > 0 and non-missing fibrosis stage.
sel <- samp$def_class == "Control" |
  (samp$def_class == "Disease" & !is.na(nas) & nas > 0 & !is.na(fib))

dge_x <- dge[, sel]
sex <- mm$inferred_sex[match(colnames(dge_x), mm$sample_id)]
info <- data.frame(
  phenotype = factor(samp$def_class[sel], levels = c("Control", "Disease")),
  dataset = droplevels(factor(samp$dataset[sel])),
  inferred_sex = factor(sex),
  nas_score = nas[sel],
  fibrosis_stage = fib[sel]
)
rownames(info) <- colnames(dge_x)
info <- info[!is.na(info$inferred_sex), , drop = FALSE]
dge_x <- dge_x[, rownames(info)]

cat(sprintf("[figS3U plotted] %d genes x %d samples; %d cohorts; %s\n",
            nrow(dge_x), ncol(dge_x), nlevels(info$dataset),
            paste(names(table(info$phenotype)), table(info$phenotype), collapse = " ")))
cat("[figS3U plotted] cohort x phenotype:\n")
print(table(info$dataset, info$phenotype))

dge_x <- calcNormFactors(dge_x)
design <- model.matrix(~ dataset + inferred_sex + phenotype, data = info)
coef_name <- "phenotypeDisease"
stopifnot(coef_name %in% colnames(design))
stopifnot(limma::is.fullrank(design))

v <- voomWithQualityWeights(dge_x, design)
fit <- eBayes(lmFit(v, design))
res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
res$gene <- rownames(res)
dt <- as.data.table(res)
setnames(dt, "adj.P.Val", "padj")
dt[, SE := abs(logFC / t)]

ok <- is.finite(dt$logFC) & is.finite(dt$SE) & dt$SE > 0
ash <- ashr::ash(dt$logFC[ok], dt$SE[ok], mixcompdist = "normal")
dt[, shrunk_logFC := NA_real_]
dt[ok, shrunk_logFC := ash$result$PosteriorMean]
dt[, lfsr := NA_real_]
dt[ok, lfsr := ash$result$lfsr]

gm <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
gm[, eb := sub("[.][0-9]+$", "", gene_id)]
dt[, symbol := gm[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]

out <- dt[, .(gene, symbol, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr)]
out[, method := "limma_voom_qw__C2__figs3u_plotted_points"]
fwrite(out, file.path(SDIR, "figs3u_plotted_points_de.csv"))

meta_out <- data.table(sample_id = colnames(dge_x),
                       dataset = as.character(info$dataset),
                       inferred_sex = as.character(info$inferred_sex),
                       phenotype = as.character(info$phenotype),
                       nas_score = info$nas_score,
                       fibrosis_stage = info$fibrosis_stage)
fwrite(meta_out, file.path(SDIR, "figs3u_plotted_points_de_meta.csv"))

can <- fread(file.path(RDIR, "canonical_deg_results.csv"))
raw1 <- out[padj < 0.05 & abs(logFC) > 0.5]
ash1 <- out[lfsr < 0.05 & abs(shrunk_logFC) > 0.5]
can_raw1 <- can[padj < 0.05 & abs(logFC) > 0.5]
can_ash1 <- can[lfsr < 0.05 & abs(shrunk_logFC) > 0.5]

m <- merge(out[, .(gene, fig_logFC = logFC, fig_padj = padj, fig_slfc = shrunk_logFC, fig_lfsr = lfsr)],
           can[, .(gene, can_logFC = logFC, can_padj = padj, can_slfc = shrunk_logFC, can_lfsr = lfsr)],
           by = "gene")
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
qs <- function(x) as.numeric(quantile(x, c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE))

overlap_raw <- intersect(raw1$gene, can_raw1$gene)
overlap_ash <- intersect(ash1$gene, can_ash1$gene)
fig_only_raw <- setdiff(raw1$gene, can_raw1$gene)
can_only_raw <- setdiff(can_raw1$gene, raw1$gene)

summary_dt <- data.table(
  metric = c(
    "samples_total", "samples_control", "samples_disease", "cohorts",
    "genes_tested_fig", "genes_tested_canonical",
    "fig_raw_Tier1", "fig_raw_Tier1_up", "fig_raw_Tier1_down",
    "canonical_raw_Tier1", "canonical_raw_Tier1_up", "canonical_raw_Tier1_down",
    "fig_ashr_Tier1", "canonical_ashr_Tier1",
    "raw_overlap_n", "raw_overlap_frac_fig", "raw_recall_of_canonical", "raw_jaccard",
    "ashr_overlap_n", "ashr_jaccard",
    "logFC_spearman_all", "logFC_pearson_all", "direction_agree_all",
    "logFC_spearman_canonical_raw_Tier1", "direction_agree_canonical_raw_Tier1",
    "abs_logFC_median_fig_all", "abs_logFC_median_can_all",
    "abs_logFC_median_fig_raw_Tier1", "abs_logFC_median_can_raw_Tier1",
    "fig_only_raw_Tier1", "canonical_only_raw_Tier1"
  ),
  value = c(
    nrow(info), sum(info$phenotype == "Control"), sum(info$phenotype == "Disease"), nlevels(info$dataset),
    nrow(out), nrow(can),
    nrow(raw1), sum(raw1$logFC > 0), sum(raw1$logFC < 0),
    nrow(can_raw1), sum(can_raw1$logFC > 0), sum(can_raw1$logFC < 0),
    nrow(ash1), nrow(can_ash1),
    length(overlap_raw), length(overlap_raw) / nrow(raw1), length(overlap_raw) / nrow(can_raw1),
    jaccard(raw1$gene, can_raw1$gene),
    length(overlap_ash), jaccard(ash1$gene, can_ash1$gene),
    cor(m$fig_logFC, m$can_logFC, method = "spearman", use = "complete.obs"),
    cor(m$fig_logFC, m$can_logFC, method = "pearson", use = "complete.obs"),
    mean(sign(m$fig_logFC) == sign(m$can_logFC), na.rm = TRUE),
    cor(m[gene %in% can_raw1$gene]$fig_logFC, m[gene %in% can_raw1$gene]$can_logFC,
        method = "spearman", use = "complete.obs"),
    mean(sign(m[gene %in% can_raw1$gene]$fig_logFC) == sign(m[gene %in% can_raw1$gene]$can_logFC),
         na.rm = TRUE),
    median(abs(out$logFC), na.rm = TRUE), median(abs(can$logFC), na.rm = TRUE),
    median(abs(raw1$logFC), na.rm = TRUE), median(abs(can_raw1$logFC), na.rm = TRUE),
    length(fig_only_raw), length(can_only_raw)
  )
)
fwrite(summary_dt, file.path(SDIR, "figs3u_plotted_points_de_summary.csv"))

lfc_quant <- rbindlist(list(
  data.table(set = "figS3U_all", stat = c("min", "p10", "p25", "median", "p75", "p90", "max"),
             value = qs(out$logFC)),
  data.table(set = "canonical_all", stat = c("min", "p10", "p25", "median", "p75", "p90", "max"),
             value = qs(can$logFC)),
  data.table(set = "figS3U_raw_Tier1_abs", stat = c("min", "p10", "p25", "median", "p75", "p90", "max"),
             value = qs(abs(raw1$logFC))),
  data.table(set = "canonical_raw_Tier1_abs", stat = c("min", "p10", "p25", "median", "p75", "p90", "max"),
             value = qs(abs(can_raw1$logFC)))
))
fwrite(lfc_quant, file.path(SDIR, "figs3u_plotted_points_lfc_quantiles.csv"))

cat("\n=== Fig S3U plotted-points DE summary ===\n")
print(summary_dt, row.names = FALSE)
cat("\n=== LFC quantiles ===\n")
print(lfc_quant, row.names = FALSE)
cat("\n[done] wrote figs3u_plotted_points_de.csv, _meta.csv, _summary.csv, _lfc_quantiles.csv\n")
