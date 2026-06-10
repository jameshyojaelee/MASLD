#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Tier 1E — Fibrosis reversal / regression molecular signature (paired biopsies)
#
# GSE193066 (Hoshida/PLS NAFLD cohort) carries paired 1st/2nd biopsies for
# 58 patients (15 progressors / 28 stable / 15 regressors), with fibrosis stage
# and NAS at each timepoint. We run a paired (within-patient) limma-voom
# quality-weighted analysis of the biopsy-2 vs biopsy-1 change, stratified by
# fibrosis trajectory, and test whether regression is the MIRROR IMAGE of
# progression or a distinct program.
#
# Method = repo canonical (limma-voom quality-weighted) + duplicateCorrelation
# blocking on patient (the standard paired-design recipe).
#
# Inputs:
#   merged_counts_raw.rds                          (gene x sample raw counts)
#   results/reversal/gse193066_paired_trajectory.csv (Gate-0 trajectory table)
#   canonical_deg_results.csv                      (gene->symbol map; mirror ref)
# Outputs (RNA-seq/results/reversal/):
#   de_regress_change.csv, de_progress_change.csv, de_stable_change.csv
#   de_reversal_vs_progression.csv
#   mirror_image_summary.txt
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(edgeR); library(limma)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
intg <- file.path(root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
outdir <- file.path(root, "RNA-seq/results/reversal")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- load counts (handle DGEList / list / matrix) --------------------------
obj <- readRDS(file.path(intg, "merged_counts_raw.rds"))
counts <- {
  if (is.matrix(obj)) obj
  else if (inherits(obj, "DGEList")) obj$counts
  else if (is.list(obj) && !is.null(obj$counts)) obj$counts
  else stop("Unrecognized merged_counts_raw.rds structure")
}
cat(sprintf("counts matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# ---- trajectory table ------------------------------------------------------
traj <- read.csv(file.path(outdir, "gse193066_paired_trajectory.csv"),
                 stringsAsFactors = FALSE)
traj <- traj[traj$run %in% colnames(counts), ]
cat(sprintf("paired samples present in matrix: %d (%d patients)\n",
            nrow(traj), length(unique(traj$patient))))
stopifnot(nrow(traj) >= 100)  # expect ~116

# keep only patients with BOTH timepoints present in the matrix
ok <- names(which(table(traj$patient) == 2))
traj <- traj[traj$patient %in% ok, ]
traj$tp    <- factor(ifelse(traj$timepoint == 1, "T1", "T2"), levels = c("T1","T2"))
traj$cls   <- factor(traj$traj_class, levels = c("stable","progress","regress"))
traj$patient <- factor(traj$patient)
cat("class x timepoint table:\n"); print(table(traj$cls, traj$tp))

# ---- DGEList + filter + TMM -------------------------------------------------
cm  <- counts[, traj$run]
grp <- factor(paste(traj$cls, traj$tp, sep = "."))           # 6 groups
dge <- DGEList(counts = cm, group = grp)
keep <- filterByExpr(dge, group = grp)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
cat(sprintf("genes after filterByExpr: %d\n", nrow(dge)))

# ---- design + paired voom (duplicateCorrelation on patient) ----------------
design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
v   <- voomWithQualityWeights(dge, design)
cf  <- duplicateCorrelation(v, design, block = traj$patient)
v   <- voomWithQualityWeights(dge, design, block = traj$patient,
                              correlation = cf$consensus.correlation)
cf  <- duplicateCorrelation(v, design, block = traj$patient)
cat(sprintf("duplicateCorrelation consensus (patient pairing): %.3f\n",
            cf$consensus.correlation))
fit <- lmFit(v, design, block = traj$patient,
             correlation = cf$consensus.correlation)

cm.contr <- makeContrasts(
  regress_change          = regress.T2  - regress.T1,
  progress_change         = progress.T2 - progress.T1,
  stable_change           = stable.T2   - stable.T1,
  reversal_vs_progression = (regress.T2 - regress.T1) - (progress.T2 - progress.T1),
  # stable-referenced (difference out the shared biopsy-to-biopsy / temporal drift,
  # which is large: stable patients with no fibrosis change still move ~1000 genes)
  progression_specific    = (progress.T2 - progress.T1) - (stable.T2 - stable.T1),
  regression_specific     = (regress.T2  - regress.T1)  - (stable.T2 - stable.T1),
  levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm.contr))

# ---- symbol map ------------------------------------------------------------
symmap <- tryCatch({
  d <- read.csv(file.path(intg, "canonical_deg_results.csv"), stringsAsFactors = FALSE)
  setNames(d$symbol, d$gene)
}, error = function(e) NULL)
annotate <- function(tt) {
  tt$gene <- rownames(tt)
  tt$symbol <- if (!is.null(symmap)) symmap[tt$gene] else NA
  tt[, c("gene","symbol", setdiff(colnames(tt), c("gene","symbol")))]
}

for (cn in colnames(cm.contr)) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  write.csv(annotate(tt), file.path(outdir, sprintf("de_%s.csv", cn)), row.names = FALSE)
  nsig <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("  %-24s : %d genes at FDR<0.05\n", cn, nsig))
}

# ---- mirror-image test --------------------------------------------------
tt <- function(cn) topTable(fit2, coef = cn, number = Inf, sort.by = "none")
rg <- tt("regress_change");  pg <- tt("progress_change"); sb <- tt("stable_change")
rs <- tt("regression_specific"); ps <- tt("progression_specific")
g  <- Reduce(intersect, list(rownames(rg), rownames(pg), rownames(sb), rownames(rs), rownames(ps)))
sp <- function(a,b) cor(a[g,"logFC"], b[g,"logFC"], method = "spearman")
# raw (uncorrected) mirror
rho_raw  <- sp(rg, pg)
# magnitude of shared temporal drift: how much does regress/progress track stable?
rho_rg_sb <- sp(rg, sb); rho_pg_sb <- sp(pg, sb)
# stable-corrected mirror (the rigorous test)
rho_corr <- sp(rs, ps)
mv <- g[ rs[g,"P.Value"] < 0.05 | ps[g,"P.Value"] < 0.05 ]
rho_corr_mv <- if (length(mv) > 20) cor(rs[mv,"logFC"], ps[mv,"logFC"], method = "spearman") else NA

sink(file.path(outdir, "mirror_image_summary.txt"))
cat("Tier 1E — fibrosis reversal vs progression (GSE193066 paired biopsies)\n")
cat(sprintf("paired patients: %d  (stable %d / progress %d / regress %d)\n",
            length(unique(traj$patient)),
            sum(traj$cls=="stable")/2, sum(traj$cls=="progress")/2, sum(traj$cls=="regress")/2))
cat(sprintf("duplicateCorrelation (patient): %.3f\n", cf$consensus.correlation))
cat("\nShared biopsy-to-biopsy drift (regress/progress change vs STABLE change):\n")
cat(sprintf("  rho(regress_change, stable_change)  = %.3f\n", rho_rg_sb))
cat(sprintf("  rho(progress_change, stable_change) = %.3f\n", rho_pg_sb))
cat("  (high => a large non-fibrosis longitudinal program dominates the raw change)\n")
cat("\nMIRROR-IMAGE TEST:\n")
cat(sprintf("  RAW       rho(regress_change, progress_change)            = %.3f\n", rho_raw))
cat(sprintf("  CORRECTED rho(regression_specific, progression_specific)  = %.3f (all %d genes)\n", rho_corr, length(g)))
cat(sprintf("            rho on %d stable-corrected moving genes         = %.3f\n", length(mv), rho_corr_mv))
cat("  Interpretation (CORRECTED is the rigorous test, drift removed):\n")
cat("    strongly NEGATIVE => regression is a transcriptional REWIND of progression;\n")
cat("    near-zero/positive => regression is a DISTINCT resolution program.\n")
sink()
cat("\nWrote outputs to", outdir, "\n")
cat(readLines(file.path(outdir, "mirror_image_summary.txt")), sep = "\n")
