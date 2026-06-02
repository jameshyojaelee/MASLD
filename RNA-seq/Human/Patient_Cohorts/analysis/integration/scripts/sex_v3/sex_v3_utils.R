#!/usr/bin/env Rscript
# sex_v3/sex_v3_utils.R
# ---------------------------------------------------------------------------
# Shared utilities for sex_v3 modules.
#
# Sourced at the top of each module to:
#   - provide atomic file write helpers (tmp + rename) so that a SLURM job
#     killed mid-write does NOT leave a corrupted .rds / .csv that silently
#     propagates wrong numbers to downstream consumers (R5 Issue 2)
#   - provide a sessionInfo dump helper that each module calls at exit so a
#     reviewer can reproduce the exact package + BLAS environment (R5 Issue 2).
#
# Atomic write pattern (R5 Issue 2):
#   1. Write to "<path>.tmp.<pid>" so concurrent module runs do not collide.
#   2. file.rename() to the canonical path — POSIX guarantees rename is atomic
#      on the same filesystem.
#   3. If the script is killed between (1) and (2), only the tmp file exists;
#      downstream consumers see "file not found" rather than a corrupt file.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Contrast routing (added 2026-05-18 for stage-stratified extension)
#
# The sex-modulation pipeline was originally fixed to a single Disease-vs-Ctrl
# contrast on the 5-cohort 847-sample canonical mega-analysis. To extend it
# to MASH-vs-MASL / MASL-vs-Ctrl / MASH-vs-Ctrl without breaking the canonical
# state on disk, each module honors the `CONTRAST_NAME` env var:
#
#   disease_vs_ctrl  (default; existing canonical state, no path changes)
#   mash_vs_masl     -> sex_v3/contrast_mash_vs_masl/...
#   masl_vs_ctrl     -> sex_v3/contrast_masl_vs_ctrl/...
#   mash_vs_ctrl     -> sex_v3/contrast_mash_vs_ctrl/...
#
# Per-contrast sample subsets and Tier-1 universes are documented in
# `contrast_spec()` below.
# ---------------------------------------------------------------------------

# Contrast specifications — sample-subset filter + Tier-1 DEG universe CSV
# + cohort filter rules. Used by Module 01 to compute the per-contrast
# sample mask; used by pillars to pick the right Tier-1 universe.
contrast_spec <- function(contrast_name = Sys.getenv("CONTRAST_NAME", "disease_vs_ctrl")) {
  specs <- list(
    disease_vs_ctrl = list(
      sample_diag       = c("Control", "NAFL", "Borderline", "NASH"),
      arm_ctrl_levels   = "Control",
      arm_dis_levels    = c("NAFL", "Borderline", "NASH"),
      use_yaml_mega     = TRUE,
      tier1_csv         = "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv",
      tier1_logFC_col   = "logFC",
      tier1_padj_col    = "padj",   # custom column name in dream_results_ashr.csv
      pretty_name       = "Disease vs Control"
    ),
    mash_vs_masl = list(
      sample_diag       = c("NAFL", "Borderline", "NASH"),
      arm_ctrl_levels   = "NAFL",
      arm_dis_levels    = c("Borderline", "NASH"),
      use_yaml_mega     = FALSE,   # control-less cohorts unlock here
      tier1_csv         = "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/mash_vs_masl_dream.csv",
      tier1_logFC_col   = "logFC",
      tier1_padj_col    = "adj.P.Val",
      pretty_name       = "MASH vs MASL"
    ),
    masl_vs_ctrl = list(
      sample_diag       = c("Control", "NAFL"),
      arm_ctrl_levels   = "Control",
      arm_dis_levels    = "NAFL",
      use_yaml_mega     = TRUE,
      tier1_csv         = "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/masl_vs_healthy_dream.csv",
      tier1_logFC_col   = "logFC",
      tier1_padj_col    = "adj.P.Val",
      pretty_name       = "MASL vs Control"
    ),
    mash_vs_ctrl = list(
      sample_diag       = c("Control", "Borderline", "NASH"),
      arm_ctrl_levels   = "Control",
      arm_dis_levels    = c("Borderline", "NASH"),
      use_yaml_mega     = TRUE,
      tier1_csv         = "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/mash_vs_healthy_dream.csv",
      tier1_logFC_col   = "logFC",
      tier1_padj_col    = "adj.P.Val",
      pretty_name       = "MASH vs Control"
    )
  )
  if (!contrast_name %in% names(specs)) {
    stop("contrast_spec: unknown CONTRAST_NAME='", contrast_name,
         "'; expected one of: ", paste(names(specs), collapse = ", "))
  }
  spec <- specs[[contrast_name]]
  spec$name <- contrast_name
  spec
}

# Per-contrast SEXV3 + IDIR paths. Default (disease_vs_ctrl) preserves the
# legacy layout `sex_v3/intermediates/`; non-default contrasts get a
# sibling subdir `sex_v3/contrast_<name>/intermediates/`.
contrast_paths <- function(contrast_name = Sys.getenv("CONTRAST_NAME", "disease_vs_ctrl"),
                            base = Sys.getenv("MASLD_PROJECT_ROOT",
                                              "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")) {
  rdir  <- file.path(base, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                     "results/integration")
  if (contrast_name == "disease_vs_ctrl") {
    sexv3 <- file.path(rdir, "sex_v3")
  } else {
    sexv3 <- file.path(rdir, "sex_v3", paste0("contrast_", contrast_name))
  }
  idir <- file.path(sexv3, "intermediates")
  list(base = base, rdir = rdir, sexv3 = sexv3, idir = idir,
       contrast = contrast_name)
}

# Atomic RDS write
write_atomic_rds <- function(obj, path, ...) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(obj, tmp, ...)
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok)) {
    if (file.exists(tmp)) file.remove(tmp)
    stop("write_atomic_rds: failed to rename ", tmp, " -> ", path)
  }
  invisible(path)
}

# Atomic CSV / TSV write (data.table::fwrite passthrough)
write_atomic_csv <- function(dt, path, ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("write_atomic_csv requires data.table.")
  }
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  data.table::fwrite(dt, tmp, ...)
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok)) {
    if (file.exists(tmp)) file.remove(tmp)
    stop("write_atomic_csv: failed to rename ", tmp, " -> ", path)
  }
  invisible(path)
}

# Atomic generic text write (writeLines passthrough)
write_atomic_lines <- function(lines, path, ...) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  con <- file(tmp, "w")
  writeLines(lines, con, ...)
  close(con)
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok)) {
    if (file.exists(tmp)) file.remove(tmp)
    stop("write_atomic_lines: failed to rename ", tmp, " -> ", path)
  }
  invisible(path)
}

# sessionInfo dump — each module calls this once at exit. Writes to
# `<intermediates_dir>/session_info_module_<NN>.txt`. Records R version,
# loaded packages, BLAS / LAPACK identity, OS — enough to diagnose any
# environmental drift between reruns (R5 Issue 2).
dump_session_info <- function(intermediates_dir, module_tag) {
  if (!dir.exists(intermediates_dir)) {
    dir.create(intermediates_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out <- file.path(intermediates_dir,
                   sprintf("session_info_module_%s.txt", module_tag))
  tmp <- paste0(out, ".tmp.", Sys.getpid())
  con <- file(tmp, "w")
  tryCatch({
    writeLines(sprintf("# sessionInfo dump for sex_v3 module %s", module_tag), con)
    writeLines(sprintf("# generated: %s",
                       format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), con)
    writeLines(sprintf("# pid: %s  host: %s",
                       Sys.getpid(), Sys.info()[["nodename"]]), con)
    writeLines("", con)
    capture <- capture.output(print(sessionInfo()))
    writeLines(capture, con)
    writeLines("", con)
    writeLines("# BLAS / LAPACK:", con)
    tryCatch({
      writeLines(sprintf("  BLAS lib:   %s", extSoftVersion()[["BLAS"]]), con)
      # La_library() is available in R >= 3.4
      la_lib <- tryCatch(La_library(), error = function(e) "(unavailable)")
      writeLines(sprintf("  LAPACK lib: %s", la_lib), con)
    }, error = function(e) {
      writeLines(sprintf("  (BLAS/LAPACK info unavailable: %s)",
                         conditionMessage(e)), con)
    })
    writeLines("", con)
    writeLines("# Sys.getenv() — sex_v3 relevant:", con)
    relevant <- c("R_PARALLEL_SEED", "SLURM_CPUS_PER_TASK", "SLURM_JOB_ID",
                  "SLURM_JOB_NAME", "OMP_NUM_THREADS", "MASLD_PROJECT_ROOT",
                  "CHRX_PILOT", "MASHR_N_STRONG", "MASHR_N_RANDOM",
                  "MASHR_N_PCA", "CAL_N_PER_PATTERN", "CAL_N_CARRIER",
                  "CAL_EFFECT_MAG", "DRY_RUN")
    for (var in relevant) {
      val <- Sys.getenv(var, unset = "")
      writeLines(sprintf("  %s = %s", var, val), con)
    }
  }, finally = close(con))
  ok <- file.rename(tmp, out)
  if (!isTRUE(ok)) {
    if (file.exists(tmp)) file.remove(tmp)
    warning("dump_session_info: rename failed ", tmp, " -> ", out)
  }
  invisible(out)
}
