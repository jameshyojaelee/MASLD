#!/usr/bin/env Rscript
# 03c_rebuild_corrected_logcpm_c2_nogerhard.R
# ===========================================================================
# Audit P2#20 remediation (2026-06-29) — rebuild the human batch-corrected
# logCPM matrix under the canonical `include_in_mega` gate.
#
# PROBLEM
#   The on-disk canonical artifact
#     results/integration/corrected_logcpm.rds   (31174 x 608, pre-C2)
#   is STALE. Its 608 columns are {GSE126848, GSE130970, GSE135251, GSE167523,
#   PRJNA512027}. This LEAKS 181 Gerhard/PRJNA512027 SRRs (dropped 2026-05-15 for
#   the L0/S0 library-prep confound) and 92 GSE167523 SRRs (no healthy controls),
#   while MISSING the two largest control-bearing mega cohorts GSE213621 and
#   GSE162694. No active script rebuilds it at the canonical path.
#
# FIX
#   Rebuild from the current C2-era 9-cohort merged_dge.rds, subset to exactly the
#   5 cohorts flagged `de.include_in_mega: true` in config/human_datasets.yaml
#   (the 846-sample canonical Disease-vs-Control set), and write to a NEW path.
#   The canonical corrected_logcpm.rds is NOT touched; the 10 downstream consumers
#   are NOT repointed (separate reviewed step).
#
# RECIPE (reconstructed — see provenance below)
#   logcpm = cpm(dge_sub, log = TRUE, prior.count = 1)        # edgeR
#   corrected = limma::removeBatchEffect(
#                 logcpm, batch = dataset,
#                 design = model.matrix(~ group_binary))       # protect disease
#
#   Provenance for this recipe (the original writer is not in git history, so this
#   is reconstructed from the consistent human-side convention across the repo):
#     - limma::removeBatchEffect, batch = dataset, design = ~group_binary is the
#       documented canonical human correction (figure-panel caption: "Fixed-effect =
#       limma::removeBatchEffect (DESeq2 ~dataset)"; git history lines for
#       230_topic_programs.R / 162_loco_nmf_subtypes.py / 208_subtype_coloc).
#     - removeBatchEffect (NOT ComBat) matches the stale file's signature: its
#       per-gene rowMeans are NOT centered at 0 (median ~ -1.4, range -4.9..+16.5),
#       i.e. the per-gene grand mean is preserved, as removeBatchEffect does.
#     - `design = ~group_binary` PROTECTS the disease contrast: cohorts differ in
#       case/control ratio, so an unprotected batch fit would partially absorb the
#       disease signal. The mouse template (M06) uses ComBat with mod=~group_binary
#       for the same protect-disease reason; the human convention is the limma
#       analogue. This is the scientifically correct choice and is documented here
#       as the explicit assumption (the original writer's exact design arg could not
#       be recovered from history; if it was unprotected, re-run with design = NULL).
#
#   logCPM is recomputed on the 846-sample subset and norm factors are
#   recalculated on the subset (calcNormFactors RLE, matching 03_integrate_counts.R)
#   so normalization is calibrated to the mega cohorts rather than inherited from
#   the 9-cohort merge.
#
# Env: rnaseq
# Output: results/integration/corrected_logcpm_c2_nogerhard.rds  (NEW path)
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(yaml)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
YAML <- file.path(BASE, "config/human_datasets.yaml")

OUT_PATH <- file.path(RDIR, "corrected_logcpm_c2_nogerhard.rds")
CANON    <- file.path(RDIR, "corrected_logcpm.rds")  # READ-ONLY here; never written

# --- 1. include_in_mega gate from yaml ---------------------------------------
ds_cfg <- yaml::read_yaml(YAML)$datasets
mega_cohorts <- names(ds_cfg)[vapply(ds_cfg,
  function(d) isTRUE(d$de$include_in_mega), logical(1))]
cat("[1] include_in_mega cohorts (n=", length(mega_cohorts), "): ",
    paste(sort(mega_cohorts), collapse = ", "), "\n", sep = "")

# --- 2. Load the C2-era 9-cohort merged DGE ----------------------------------
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat(sprintf("[2] merged_dge.rds: %d genes x %d samples; cohorts: %s\n",
            nrow(dge), ncol(dge), paste(sort(unique(dge$samples$dataset)), collapse = ", ")))

# Hard guard: Gerhard must already be absent from the C2 merge
stopifnot(!"PRJNA512027" %in% dge$samples$dataset)

# --- 3. Subset to the 5 mega cohorts -----------------------------------------
keep_samp <- dge$samples$dataset %in% mega_cohorts
dge_sub   <- dge[, keep_samp]
dge_sub$samples <- droplevels(dge_sub$samples)
cat(sprintf("[3] Subset to mega cohorts: %d samples\n", ncol(dge_sub)))
print(table(dge_sub$samples$dataset, dge_sub$samples$group_binary))

# --- 4. Re-filter + re-normalize on the subset -------------------------------
# Recompute the expressed-gene universe and RLE norm factors on the 846-sample
# mega set (matches 03_integrate_counts.R: filterByExpr + calcNormFactors RLE).
design_filter <- model.matrix(~ 0 + group_binary, data = dge_sub$samples)
keep_genes <- filterByExpr(dge_sub, design = design_filter)
cat(sprintf("[4] filterByExpr on subset: %d / %d genes retained\n",
            sum(keep_genes), nrow(dge_sub)))
dge_sub <- dge_sub[keep_genes, , keep.lib.sizes = FALSE]
dge_sub <- calcNormFactors(dge_sub, method = "RLE")

# --- 5. logCPM + removeBatchEffect (disease-protected) ------------------------
logcpm <- cpm(dge_sub, log = TRUE, prior.count = 1)
batch  <- factor(dge_sub$samples$dataset)
design <- model.matrix(~ group_binary, data = dge_sub$samples)  # protect disease
corrected <- limma::removeBatchEffect(logcpm, batch = batch, design = design)
cat(sprintf("[5] removeBatchEffect complete: %d genes x %d samples\n",
            nrow(corrected), ncol(corrected)))

# --- 6. Assertions -----------------------------------------------------------
# 6a. No Gerhard SRR in output columns. Build the PRJNA512027 SRR set from the
#     stale canonical file's columns labelled Gerhard (belt-and-suspenders: also
#     check the unified-metadata dataset label).
um <- fread(file.path(INT, "metadata/unified_metadata.csv"))
gerhard_srrs <- um[dataset == "PRJNA512027", sample_id]
cat(sprintf("[6] %d known Gerhard/PRJNA512027 SRRs in unified metadata\n",
            length(gerhard_srrs)))
leaked <- intersect(colnames(corrected), gerhard_srrs)
stopifnot(length(leaked) == 0)
# Also assert no column maps to PRJNA512027 by dataset label
out_ds <- um[match(colnames(corrected), sample_id), dataset]
stopifnot(!any(out_ds == "PRJNA512027", na.rm = TRUE))
stopifnot(all(sort(unique(out_ds)) == sort(mega_cohorts)))
cat("    ASSERT PASS: zero PRJNA512027/Gerhard columns; only mega cohorts present.\n")

# --- 7. Write NEW artifact (never the canonical path) ------------------------
stopifnot(normalizePath(OUT_PATH, mustWork = FALSE) !=
          normalizePath(CANON,    mustWork = FALSE))
saveRDS(corrected, OUT_PATH)
cat(sprintf("\n[7] Wrote NEW artifact: %s\n", OUT_PATH))
cat(sprintf("    dims: %d genes x %d samples\n", nrow(corrected), ncol(corrected)))
cat("    cohort breakdown:\n")
print(table(out_ds))
cat("\nDone. Canonical corrected_logcpm.rds was NOT modified.\n")
