#!/usr/bin/env Rscript

# Grouped NAS activity contrasts using the canonical LVQW engine.
# Groups: NAS0 (reference), NAS1-2, NAS3-4, NAS5-8.
# All QC-passing NAS-annotated samples are eligible regardless of diagnosis.

suppressPackageStartupMessages({
  library(ashr)
  library(data.table)
  library(edgeR)
  library(limma)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUTDIR <- file.path(INT, "results/disease_signatures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

dge_all <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_m <- as.data.table(readRDS(file.path(RDIR, "meta_matched.rds")))
meta <- fread(file.path(INT, "metadata/unified_metadata.csv"),
              select = c("sample_id", "dataset", "nas_score", "group_binary",
                         "diagnosis_harmonized"))
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"),
            select = c("sample_id", "pass_technical"))
meta <- merge(meta, meta_m[, .(sample_id, inferred_sex)], by = "sample_id", all.x = TRUE)
meta <- merge(meta, qc, by = "sample_id", all.x = TRUE)

NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")
meta <- meta[pass_technical == TRUE & dataset %in% NAS_DATASETS &
             sample_id %in% colnames(dge_all)]
meta[, nas_num := suppressWarnings(as.integer(nas_score))]
meta <- meta[!is.na(nas_num) & nas_num %in% 0:8]
meta[, nas_group := fcase(
  nas_num == 0, "NAS0",
  nas_num <= 2, "NAS1-2",
  nas_num <= 4, "NAS3-4",
  default = "NAS5-8"
)]
meta[, nas_group := factor(nas_group,
                           levels = c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"))]

idx <- colnames(dge_all) %in% meta$sample_id
dge <- dge_all[, idx]
meta <- meta[match(colnames(dge), sample_id)]
dge$samples$dataset <- factor(meta$dataset)
dge$samples$inferred_sex <- factor(meta$inferred_sex)
dge$samples$nas_group <- meta$nas_group

info <- data.frame(
  dataset = factor(meta$dataset),
  inferred_sex = factor(meta$inferred_sex),
  nas_group = relevel(meta$nas_group, ref = "NAS0"),
  row.names = meta$sample_id
)
design <- model.matrix(~ dataset + inferred_sex + nas_group, data = info)
stopifnot(qr(design)$rank == ncol(design))

cat(sprintf("Grouped NAS LVQW: %d genes x %d samples\n", nrow(dge), ncol(dge)))
print(meta[, .N, by = .(nas_group, group_binary)][order(nas_group, group_binary)])

# Fit voom quality weights once; extract all three grouped coefficients.
# TREAT (McCarthy & Smyth 2009) is the paper-wide canonical DEG gate (2026-06-29):
# treat() tests H0:|true logFC|<=lfc, folding the effect-size floor into the FDR.
# Mirrors 05h_limma_voom_qw_canonical.R exactly: treat() on the raw lmFit, lfc=0.25.
TREAT_LFC <- as.numeric(Sys.getenv("CANONICAL_TREAT_LFC", "0.25"))
v    <- voomWithQualityWeights(dge, design, plot = FALSE)
fit0 <- lmFit(v, design)
fit  <- eBayes(fit0)
tfit <- treat(fit0, lfc = TREAT_LFC)   # one moderated fit; topTreat() extracts per coef
saveRDS(fit, file.path(OUTDIR, "nas_grouped_vs_nas0_lvqw_fit.rds"))

coef_map <- c(
  "NAS1-2_vs_NAS0" = "nas_groupNAS1-2",
  "NAS3-4_vs_NAS0" = "nas_groupNAS3-4",
  "NAS5-8_vs_NAS0" = "nas_groupNAS5-8"
)

results <- rbindlist(lapply(names(coef_map), function(contrast_name) {
  coef_name <- coef_map[[contrast_name]]
  stopifnot(coef_name %in% colnames(design))
  tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf,
                              sort.by = "none"), keep.rownames = "gene")
  setnames(tt, "adj.P.Val", "padj")
  tt[, SE := abs(logFC / t)]
  ok <- is.finite(tt$logFC) & is.finite(tt$SE) & tt$SE > 0
  afit <- ashr::ash(tt$logFC[ok], tt$SE[ok], mixcompdist = "normal")
  tt[, `:=`(shrunk_logFC = NA_real_, lfsr = NA_real_)]
  tt[ok, `:=`(shrunk_logFC = afit$result$PosteriorMean,
              lfsr = afit$result$lfsr)]
  ttm <- topTreat(tfit, coef = coef_name, number = Inf, sort.by = "none")
  tt[, `:=`(treat_lfc = TREAT_LFC,
            treat_p   = ttm[gene, "P.Value"],
            treat_fdr = ttm[gene, "adj.P.Val"])]
  tt[, `:=`(contrast = contrast_name, method = "limma_voom_qw")]
  tt[, .(gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr,
         treat_lfc, treat_p, treat_fdr, AveExpr, contrast, method)]
}))

sizes <- meta[, .N, by = .(nas_group, group_binary)][order(nas_group, group_binary)]
fwrite(results, file.path(OUTDIR, "nas_grouped_vs_nas0_lvqw.csv"))
fwrite(sizes, file.path(OUTDIR, "nas_grouped_vs_nas0_sample_sizes.csv"))

cat(sprintf("\nCANONICAL DEG counts (TREAT: treat_fdr<0.05 @ lfc=%.2f):\n", TREAT_LFC))
print(results[treat_fdr < 0.05,
              .(up = sum(logFC > 0), down = sum(logFC < 0)),
              by = contrast])
cat("\n[reference] retired ashr gate (lfsr<0.05, |shrunk logFC|>0.3):\n")
print(results[lfsr < 0.05 & abs(shrunk_logFC) > 0.3,
              .(up = sum(shrunk_logFC > 0), down = sum(shrunk_logFC < 0)),
              by = contrast])
