#!/usr/bin/env Rscript
# extract_dialogue_phenoz.R
# Extract DIALOGUE phenotype-Z scores and donor-level MCP scores from the
# rerun_1khvg outputs, produce a long-form table for Fig 2 supplementary
# panel D7 (cNMF k=16 algorithm-independent validation).
#
# Per RESULTS_FINAL_k16.md: at 1K HVG/cell-type, DIALOGUE produced empty
# gene signatures (sig.up/sig.down length 0) but the phenotype-projection
# phase (donor scores, phenoZ) IS populated. Two phenotype columns were
# present in donor_metadata.tsv: condition_binary_num and disease_stage_numeric.
# k = {3, 5, 7} were swept.
#
# Inputs:
#   results_gpu_v2/mcp/dialogue/rerun_1khvg/{phenotype}_k{k}DLG.output_*.rds
#
# Outputs:
#   - dialogue_phenoZ_donor_long.tsv
#       (phenotype, k, mcp_id, cell_type, donor, pheno_z)  -- donor-level MCP scores
#   - dialogue_phenoZ_celltype_long.tsv
#       (phenotype, k, mcp_id, cell_type, pheno_z)         -- per-(cell_type, MCP) phenotype Z stats
#   - dialogue_phenoZ_extraction_summary.tsv
#       (phenotype, k, n_donors, n_celltypes, n_mcps, scores_populated, phenoZ_populated)

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
in_dir  <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/dialogue/rerun_1khvg")
out_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/dialogue")

stopifnot(dir.exists(in_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Discover RDS files at top level of rerun_1khvg/. We use the smaller
# DLG.output_*.rds (NOT DLG.full.output_*.rds) since phenoZ + scores are
# identical between the two and the smaller file loads faster.
rds_files <- list.files(
  in_dir,
  pattern = "^[A-Za-z0-9_]+_k[0-9]+DLG\\.output_[A-Za-z0-9_]+_k[0-9]+\\.rds$",
  full.names = TRUE
)
cat(sprintf("[extract] Found %d DLG.output RDS files\n", length(rds_files)))
if (length(rds_files) == 0) stop("No DIALOGUE RDS files found in rerun_1khvg/")

donor_long_list    <- list()
celltype_long_list <- list()
summary_rows       <- list()

for (f in rds_files) {
  bn <- basename(f)
  # Filename pattern: {phenotype}_k{K}DLG.output_{phenotype}_k{K}.rds
  m <- regmatches(bn, regexec("^(.+)_k([0-9]+)DLG\\.output_.*\\.rds$", bn))[[1]]
  if (length(m) != 3) {
    cat(sprintf("[extract] could not parse %s -- skipping\n", bn))
    next
  }
  phenotype <- m[2]
  k_val     <- as.integer(m[3])
  cat(sprintf("[extract] Loading %s (phenotype=%s, k=%d)\n", bn, phenotype, k_val))

  obj <- tryCatch(readRDS(f), error = function(e) {
    cat(sprintf("  ! readRDS failed: %s\n", conditionMessage(e))); NULL
  })
  if (is.null(obj)) next

  # ---- 1. Per-(cell_type, MCP) phenoZ matrix ----
  pZ <- obj$phenoZ
  phenoZ_populated <- !is.null(pZ) && is.matrix(pZ) && length(pZ) > 0 && any(!is.na(pZ))
  if (phenoZ_populated) {
    pZ_dt <- as.data.table(as.table(pZ))
    setnames(pZ_dt, c("cell_type", "mcp_id", "pheno_z"))
    pZ_dt[, phenotype := phenotype]
    pZ_dt[, k := k_val]
    setcolorder(pZ_dt, c("phenotype", "k", "mcp_id", "cell_type", "pheno_z"))
    celltype_long_list[[length(celltype_long_list) + 1]] <- pZ_dt
  } else {
    cat("  ! phenoZ not populated\n")
  }

  # ---- 2. Per-donor MCP scores (donor x cell_type x MCP) ----
  scores_populated <- !is.null(obj$scores) && length(obj$scores) > 0
  donor_count <- 0L
  ct_count    <- 0L
  mcp_count   <- 0L
  if (scores_populated) {
    per_ct <- list()
    for (ct in names(obj$scores)) {
      sc <- obj$scores[[ct]]
      if (is.null(sc) || nrow(sc) == 0) next
      mcp_cols <- grep("^MCP[0-9]+$", colnames(sc), value = TRUE)
      if (length(mcp_cols) == 0) next
      sc_dt <- as.data.table(sc)
      # Donor identifier lives in `samples` column (each donor has _r1 and _r2 rows)
      if (!"samples" %in% colnames(sc_dt)) {
        cat(sprintf("  ! %s scores missing 'samples' column -- skipping\n", ct))
        next
      }
      melt_dt <- melt(
        sc_dt,
        id.vars      = "samples",
        measure.vars = mcp_cols,
        variable.name = "mcp_id",
        value.name   = "pheno_z"
      )
      # Collapse the duplicated _r1/_r2 cells back to the original donor
      donor_dt <- melt_dt[, .(pheno_z = mean(pheno_z, na.rm = TRUE)),
                         by = .(samples, mcp_id)]
      setnames(donor_dt, "samples", "donor")
      donor_dt[, cell_type := ct]
      per_ct[[ct]] <- donor_dt
    }
    if (length(per_ct) > 0) {
      donor_dt_all <- rbindlist(per_ct, use.names = TRUE)
      donor_dt_all[, phenotype := phenotype]
      donor_dt_all[, k := k_val]
      setcolorder(donor_dt_all,
                  c("phenotype", "k", "mcp_id", "cell_type", "donor", "pheno_z"))
      donor_long_list[[length(donor_long_list) + 1]] <- donor_dt_all
      donor_count <- uniqueN(donor_dt_all$donor)
      ct_count    <- uniqueN(donor_dt_all$cell_type)
      mcp_count   <- uniqueN(donor_dt_all$mcp_id)
    } else {
      scores_populated <- FALSE
    }
  } else {
    cat("  ! scores not populated\n")
  }

  summary_rows[[length(summary_rows) + 1]] <- data.table(
    phenotype        = phenotype,
    k                = k_val,
    n_donors         = donor_count,
    n_celltypes      = ct_count,
    n_mcps           = mcp_count,
    scores_populated = scores_populated,
    phenoZ_populated = phenoZ_populated,
    file             = bn
  )

  cat(sprintf("  donors=%d celltypes=%d MCPs=%d phenoZ=%s scores=%s\n",
              donor_count, ct_count, mcp_count,
              phenoZ_populated, scores_populated))
}

# ---- Write outputs ----
if (length(donor_long_list) > 0) {
  donor_long <- rbindlist(donor_long_list, use.names = TRUE)
  out1 <- file.path(out_dir, "dialogue_phenoZ_donor_long.tsv")
  fwrite(donor_long, out1, sep = "\t")
  cat(sprintf("[extract] WROTE %s (%d rows)\n", out1, nrow(donor_long)))
} else {
  cat("[extract] no donor-level rows extracted\n")
}

if (length(celltype_long_list) > 0) {
  ct_long <- rbindlist(celltype_long_list, use.names = TRUE)
  out2 <- file.path(out_dir, "dialogue_phenoZ_celltype_long.tsv")
  fwrite(ct_long, out2, sep = "\t")
  cat(sprintf("[extract] WROTE %s (%d rows)\n", out2, nrow(ct_long)))
}

if (length(summary_rows) > 0) {
  summary_dt <- rbindlist(summary_rows, use.names = TRUE)
  out3 <- file.path(out_dir, "dialogue_phenoZ_extraction_summary.tsv")
  fwrite(summary_dt, out3, sep = "\t")
  cat(sprintf("[extract] WROTE %s (%d rows)\n", out3, nrow(summary_dt)))
  cat("\n=== Extraction summary ===\n")
  print(summary_dt)
}

cat("[extract] DONE.\n")
