#!/usr/bin/env Rscript
# ============================================================================
# run_toast_celltype_de.R
#
# TOAST cell-type-specific DE for each F-transition (F0->F1, F1->F2, F2->F3,
# F3->F4) using BayesPrism proportions as mediators. Complements the
# composition-mediation pipeline: gives every gene a per-CT effect size at
# every F-transition (TOAST partial regression) so we can attribute mediation
# IE to a specific CT-of-origin.
#
# Refs:
#   Li 2019 Bioinformatics PMID 31050385 (TOAST)
#   Meng 2023 Brief Bioinformatics PMID 36472568 (CT-DE benchmark)
#   Sun 2022 bioRxiv (CARseq, conceptual cousin)
#
# Inputs:
#   - merged_dge.rds (counts, all cohorts)
#   - bayesprism_proportions.csv
#   - unified_metadata.csv + sample_qc_report.csv
#
# Output: results/mediation/celltype_de_TOAST.csv
#   columns: cell_type, transition, gene, symbol, logFC, p_value, padj
#
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TOAST)
  library(edgeR)
  library(limma)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
META_F <- file.path(INT, "metadata/unified_metadata.csv")
QC_F   <- file.path(INT, "qc/sample_qc_report.csv")
PROP_F <- file.path(INT, "results/progression/cibersortx_celltype_expression",
                    "bayesprism_proportions.csv")
GENC_F <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
OUTDIR <- file.path(INT, "results/mediation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MEDIATORS  <- c("Hepatocyte", "Stellate", "Cholangiocyte", "Macrophage", "Endothelial")
TRANS_PAIRS <- list(
  F1_vs_F0 = c(0L, 1L), F2_vs_F1 = c(1L, 2L),
  F3_vs_F2 = c(2L, 3L), F4_vs_F3 = c(3L, 4L)
)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
message("[1] Loading metadata + QC + counts...")
meta <- fread(META_F)
qc   <- fread(QC_F)
meta <- merge(meta, qc[, .(sample_id, pass_technical)], by = "sample_id",
              all.x = TRUE)
meta_keep <- meta[!is.na(pass_technical) & pass_technical == TRUE &
                    !is.na(fibrosis_stage) & !is.na(sex)]
meta_keep[, fibrosis_stage := as.integer(fibrosis_stage)]

dge  <- readRDS(file.path(RDIR, "merged_dge.rds"))
genc <- fread(GENC_F)
pc_ids <- genc[gene_biotype == "protein_coding", gene_id]
dge <- dge[rownames(dge) %in% pc_ids, , keep.lib.sizes = FALSE]
keep_samp <- colnames(dge) %in% meta_keep$sample_id
dge <- dge[, keep_samp]
meta_aln <- meta_keep[match(colnames(dge), sample_id)]
stopifnot(all(meta_aln$sample_id == colnames(dge)))

prop <- fread(PROP_F)
prop_aln <- prop[match(colnames(dge), sample_id)]
stopifnot(all(prop_aln$sample_id == colnames(dge)))
prop_mat <- as.matrix(prop_aln[, MEDIATORS, with = FALSE])
rownames(prop_mat) <- colnames(dge)
# Normalize to sum-to-1 across mediators (re-compose against missing CTs)
prop_mat <- prop_mat / rowSums(prop_mat)
# Drop smallest mediator to break the simplex collinearity inside TOAST
ct_means <- colMeans(prop_mat)
drop_ct  <- names(ct_means)[which.min(ct_means)]
prop_mat <- prop_mat[, setdiff(MEDIATORS, drop_ct), drop = FALSE]
prop_mat <- prop_mat / rowSums(prop_mat)
ct_kept  <- colnames(prop_mat)
message(sprintf("  Cell types in TOAST: %s (dropped %s for simplex identifiability)",
                paste(ct_kept, collapse = ", "), drop_ct))

# Symbol map
dream <- fread(file.path(RDIR, "dream_results_ashr.csv"))
sym_map <- unique(dream[, .(gene, symbol)])

# ---------------------------------------------------------------------------
# Per-transition TOAST loop
# ---------------------------------------------------------------------------
all_results <- list()

for (tr in names(TRANS_PAIRS)) {
  pair <- TRANS_PAIRS[[tr]]
  pre  <- pair[1]; post <- pair[2]

  message(sprintf("\n[2] === Transition %s (F%d vs F%d) ===", tr, post, pre))
  ix <- which(meta_aln$fibrosis_stage %in% c(pre, post))
  meta_t <- meta_aln[ix]
  if (nrow(meta_t) < 30) {
    message(sprintf("  Skipping (n=%d too small)", nrow(meta_t)))
    next
  }
  message(sprintf("  Donors: pre=%d, post=%d (total=%d)",
                  sum(meta_t$fibrosis_stage == pre),
                  sum(meta_t$fibrosis_stage == post),
                  nrow(meta_t)))

  dge_t <- dge[, ix, keep.lib.sizes = FALSE]
  prop_t <- prop_mat[ix, , drop = FALSE]
  prop_t <- prop_t / rowSums(prop_t)

  grp <- factor(ifelse(meta_t$fibrosis_stage == post, "Post", "Pre"),
                levels = c("Pre", "Post"))

  # Voom transform
  keep_genes <- filterByExpr(dge_t, group = grp)
  dge_t <- dge_t[keep_genes, , keep.lib.sizes = FALSE]
  dge_t <- calcNormFactors(dge_t)

  # Compute voom (use group-only design for variance modelling; dataset
  # adjustment lives in the lm step inside TOAST via design_info)
  ds <- factor(meta_t$dataset)
  if (length(unique(ds)) >= 2) {
    des <- model.matrix(~ grp + ds)
  } else {
    des <- model.matrix(~ grp)
  }
  v <- voom(dge_t, design = des, plot = FALSE)
  Y <- v$E
  message(sprintf("  Genes after filterByExpr: %d", nrow(Y)))

  design_info <- data.frame(group = grp, row.names = colnames(Y))
  design_out  <- makeDesign(design_info, prop_t)
  fit <- fitModel(design_out, Y)

  for (ct in ct_kept) {
    res <- tryCatch({
      csTest(fit, coef = "group", cell_type = ct, contrast_matrix = NULL,
             verbose = FALSE)
    }, error = function(e) {
      message(sprintf("    csTest failed for %s: %s", ct, e$message))
      NULL
    })
    if (is.null(res)) next
    dt <- as.data.table(res, keep.rownames = "gene")
    # Standardize column names: TOAST returns beta, beta_var, p_value, fdr, ...
    if ("p_value" %in% names(dt)) setnames(dt, "p_value", "P.Value")
    if ("fdr" %in% names(dt))     setnames(dt, "fdr",     "padj")
    if ("beta" %in% names(dt))    setnames(dt, "beta",    "logFC")
    dt[, cell_type := ct]
    dt[, transition := tr]
    all_results[[paste(tr, ct, sep = "::")]] <- dt
  }
}

results_all <- rbindlist(all_results, fill = TRUE)
results_all <- merge(results_all, sym_map, by = "gene", all.x = TRUE)
out <- file.path(OUTDIR, "celltype_de_TOAST.csv")
fwrite(results_all, out)
message(sprintf("\n[3] Wrote: %s (%d rows)", out, nrow(results_all)))

# Summary
if ("padj" %in% names(results_all)) {
  smry <- results_all[, .(n_sig_05 = sum(padj < 0.05, na.rm = TRUE),
                           n_sig_10 = sum(padj < 0.10, na.rm = TRUE)),
                       by = .(transition, cell_type)]
  smry <- smry[order(transition, -n_sig_05)]
  print(smry)
  fwrite(smry, file.path(OUTDIR, "celltype_de_TOAST_summary.csv"))
}

message("Done.")
