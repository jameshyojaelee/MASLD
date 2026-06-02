#!/usr/bin/env Rscript
# 83_composition_shifts.R
#
# Analysis A2 (v1) — Cell-type composition shifts across MASLD disease stages.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   1. Load MuSiC proportions from all 10 bulk datasets (per-sample per-cell-type).
#   2. Merge with unified metadata (condition, fibrosis_stage, NAS, dataset, sex).
#   3. Propeller-equivalent composition test:
#      - logit-transform (asin sqrt) proportions
#      - limma per cell type with group contrast + dataset/sex as covariates
#      - Contrasts: MASLD vs Healthy; NAS tertile; fibrosis stage F0-F4
#   4. F2-stage boundary test: compare F<=1 vs F>=2 per cell type (no longer central to narrative;
#      "F2 switch" framing retired 2026-05-09; see paper_outline.md).
#   5. Output composition shifts + per-sample proportion matrix for downstream.
#
# v2 (pending install): replace limma with speckle::propeller + scCODA Bayesian CI.
#
# Env: rnaseq
# Outputs: RNA-seq/results/celltype_attribution/

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

BASE      <- Sys.getenv("MASLD_PROJECT_ROOT",
                        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DECONV    <- file.path(BASE, "Analysis/Deconvolution/results")
INT       <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
META      <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
OUTDIR    <- file.path(BASE, "RNA-seq/results/celltype_attribution")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading MuSiC proportions across all datasets...")
datasets <- list.dirs(DECONV, recursive = FALSE, full.names = FALSE)
prop_files <- file.path(DECONV, datasets, paste0(datasets, "_music_prop_weighted.tsv"))
keep <- file.exists(prop_files)
datasets   <- datasets[keep]
prop_files <- prop_files[keep]
message(sprintf("  %d datasets with proportions", length(datasets)))

prop_list <- lapply(seq_along(datasets), function(i) {
  ds <- datasets[i]
  df <- fread(prop_files[i])
  if (ncol(df) < 2) return(NULL)
  # First col is sample; rest are cell types
  sample_col <- names(df)[1]
  setnames(df, sample_col, "sample_id")
  df[, dataset := ds]
  df
})
props <- rbindlist(prop_list, fill = TRUE, use.names = TRUE)
message(sprintf("  Total samples with proportions: %d", nrow(props)))

# Cell types: all columns except sample_id/dataset
celltypes <- setdiff(names(props), c("sample_id", "dataset"))
message(sprintf("  Cell types: %s", paste(celltypes, collapse = ", ")))

message("[2] Joining with metadata...")
meta_dt <- as.data.table(META)
# Identify sample id column in meta
sid_candidates <- c("sample_id", "sample", "sampleID")
sid_col <- intersect(sid_candidates, names(meta_dt))[1]
if (is.na(sid_col)) stop("No sample id column in metadata")
setnames(meta_dt, sid_col, "sample_id")

# Keep relevant columns
keep_cols <- intersect(c("sample_id","dataset","group_binary","condition",
                         "fibrosis_stage","nas","inferred_sex","diagnosis_harmonized"),
                       names(meta_dt))
meta_sub <- unique(meta_dt[, ..keep_cols])

merged <- merge(props, meta_sub, by = c("sample_id","dataset"), all.x = TRUE)
message(sprintf("  Merged rows: %d", nrow(merged)))

# Keep only samples with group_binary assignment
merged <- merged[!is.na(group_binary)]
message(sprintf("  After filtering missing group_binary: %d", nrow(merged)))

message("[3] Composition shift test (propeller-equivalent)...")
# logit transform (arcsin-sqrt)
transf <- function(p) asin(sqrt(pmax(pmin(p, 1 - 1e-6), 1e-6)))

ct_test <- function(ct, contrast_col = "group_binary", levels_pair = c("Control","Disease")) {
  x <- merged[[ct]]
  keep <- !is.na(x)
  if (sum(keep) < 30) return(NULL)
  dat <- merged[keep]
  y <- transf(dat[[ct]])
  grp <- factor(dat[[contrast_col]], levels = levels_pair)
  if (length(unique(grp)) < 2) return(NULL)
  # limma with dataset + sex as covariates
  covars <- c()
  if ("dataset" %in% names(dat) && length(unique(dat$dataset)) > 1) covars <- c(covars, "dataset")
  if ("inferred_sex" %in% names(dat) && length(unique(na.omit(dat$inferred_sex))) > 1) covars <- c(covars, "inferred_sex")
  form <- as.formula(paste("~ grp", if (length(covars)) paste("+", paste(covars, collapse = "+")) else ""))
  design <- model.matrix(form, data = dat)
  # Sanity
  if (nrow(design) != length(y)) return(NULL)
  fit <- lmFit(matrix(y, nrow = 1), design)
  fit <- eBayes(fit)
  coef_name <- paste0("grp", levels_pair[2])
  if (!coef_name %in% colnames(design)) return(NULL)
  tbl <- topTable(fit, coef = coef_name, number = 1, sort.by = "none")
  data.table(
    celltype = ct,
    contrast = paste(levels_pair, collapse = "_vs_"),
    logit_diff = tbl$logFC,
    t = tbl$t,
    pvalue = tbl$P.Value,
    n = sum(keep),
    n_group1 = sum(grp == levels_pair[1], na.rm = TRUE),
    n_group2 = sum(grp == levels_pair[2], na.rm = TRUE)
  )
}

# Primary contrast: MASLD vs Healthy
res_primary <- rbindlist(lapply(celltypes, ct_test,
                                contrast_col = "group_binary",
                                levels_pair = c("Control","Disease")),
                         fill = TRUE)
res_primary[, padj := p.adjust(pvalue, method = "BH")]

message("[4] F2-stage boundary test (F<=1 vs F>=2; no longer central to narrative)...")
merged[, f2_group := fifelse(is.na(fibrosis_stage), NA_character_,
                             fifelse(as.character(fibrosis_stage) %in% c("0","1","F0","F1"),
                                     "F_low", "F_high"))]

res_f2 <- rbindlist(lapply(celltypes, function(ct) {
  sub <- merged[!is.na(f2_group)]
  x <- sub[[ct]]
  keep <- !is.na(x)
  if (sum(keep) < 30) return(NULL)
  sub <- sub[keep]
  y <- transf(sub[[ct]])
  grp <- factor(sub$f2_group, levels = c("F_low","F_high"))
  covars <- c()
  if (length(unique(sub$dataset)) > 1) covars <- c(covars, "dataset")
  if ("inferred_sex" %in% names(sub) && length(unique(na.omit(sub$inferred_sex))) > 1) covars <- c(covars, "inferred_sex")
  form <- as.formula(paste("~ grp", if (length(covars)) paste("+", paste(covars, collapse = "+")) else ""))
  design <- model.matrix(form, data = sub)
  fit <- lmFit(matrix(y, nrow = 1), design)
  fit <- eBayes(fit)
  if (!"grpF_high" %in% colnames(design)) return(NULL)
  tbl <- topTable(fit, coef = "grpF_high", number = 1, sort.by = "none")
  data.table(celltype = ct, contrast = "F_low_vs_F_high",
             logit_diff = tbl$logFC, t = tbl$t, pvalue = tbl$P.Value,
             n = length(y), n_group1 = sum(grp == "F_low"),
             n_group2 = sum(grp == "F_high"))
}), fill = TRUE)
if (nrow(res_f2) > 0) res_f2[, padj := p.adjust(pvalue, method = "BH")]

message("[5] Writing outputs...")
all_shifts <- rbind(res_primary, res_f2, fill = TRUE)
fwrite(all_shifts, file.path(OUTDIR, "composition_shifts.csv"))
fwrite(merged, file.path(OUTDIR, "persample_celltype_proportions.csv"))

summary_lines <- c(
  sprintf("Samples with proportions + group_binary: %d", nrow(merged)),
  sprintf("Cell types: %d", length(celltypes)),
  "",
  "=== MASLD vs Healthy composition shifts (sorted by |t|) ===",
  capture.output(print(res_primary[order(-abs(t))], nrows = 30)),
  "",
  "=== F2-stage boundary (F_low vs F_high) ===",
  capture.output(print(res_f2[order(-abs(t))], nrows = 30)),
  "",
  "NOTE: Uses logit-transformed proportions + limma (propeller-equivalent).",
  "      v2 will use speckle::propeller + scCODA Bayesian CI.")
writeLines(summary_lines, file.path(OUTDIR, "composition_shifts_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
