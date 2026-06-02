#!/usr/bin/env Rscript
# sex_v3/09_atlas_refresh.R
# ---------------------------------------------------------------------------
# Module 09 — Atlas refresh driver.
#
# After Module 06 lands sex_deg_classification_v3.csv, this module chains the
# multi-evidence atlas rebuild + downstream causal-evidence refresh:
#
#     27a  ->  75  ->  217
#
# Run as three sbatch submissions linked via --dependency=afterok so that
# atlas evidence assembly (27a) completes before the causal integration (75)
# and stratified-causal projection (217) consume it.
#
# This script is a thin orchestrator — does NOT itself rebuild anything;
# it only emits sbatch commands and logs the submission chain. Run on a
# login node (the orchestration itself is light).
#
# Preconditions (checked before submission):
#   - sex_deg_classification_v3.csv exists (Module 06 output)
#   - 27a backup (27a_assemble_evidence_atlas.R.v2-backup) exists
#   - run_27ab_atlas.sh + 75 + 217 source scripts present
#
# Postconditions (asserted by downstream consumers, not here):
#   - results/multi_evidence/multi_evidence_atlas.csv has a sex_class_v3_*
#     column family (Layer 5 reads from sex_deg_classification_v3.csv via the
#     27a one-line fallback patch added on 2026-05-13).
#
# Log: outputs/team_B/B6_atlas_refresh.log
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
SEXV3  <- file.path(RDIR, "sex_v3")
IDIR   <- file.path(SEXV3, "intermediates")
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)

# Shared utilities (atomic writes + sessionInfo dump) -- R5 Issue 2 fix
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

# Outputs (orchestrator log + a sentinel for downstream agents)
LOG_F  <- file.path(BASE, "outputs/team_B/B6_atlas_refresh.log")
dir.create(dirname(LOG_F), recursive = TRUE, showWarnings = FALSE)

V3_CSV       <- file.path(SEXV3, "sex_deg_classification_v3.csv")
BACKUP_27A   <- file.path(BASE, "RNA-seq/27a_assemble_evidence_atlas.R.v2-backup")
RUN_27AB     <- file.path(BASE, "RNA-seq/run_27ab_atlas.sh")
RUN_75_217   <- file.path(BASE, "RNA-seq/run_75_217_after_atlas.sh")

# Helpers
log_line <- function(msg) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG_F, append = TRUE)
}

# Re-init log
if (file.exists(LOG_F)) file.remove(LOG_F)
cat("=== Module 09 — sex_v3 atlas refresh driver ===\n", file = LOG_F)

log_line("Module 09 (atlas refresh driver) starting.")
log_line(sprintf("BASE = %s", BASE))

# ---------------------------------------------------------------------------
# 1) Preflight checks
# ---------------------------------------------------------------------------
log_line("[1] Preflight checks ...")

ok <- TRUE
for (path in list(V3_CSV = V3_CSV, BACKUP_27A = BACKUP_27A,
                  RUN_27AB = RUN_27AB, RUN_75_217 = RUN_75_217)) {
  exists_flag <- file.exists(path)
  log_line(sprintf("  %s exists = %s : %s",
                   ifelse(exists_flag, "OK", "MISSING"), exists_flag, path))
  if (!exists_flag) ok <- FALSE
}

if (!ok) {
  log_line("Preflight FAILED. Aborting before any sbatch submission.")
  stop("Module 09 preflight failed — see ", LOG_F)
}

# Quick schema sniff on v3 csv: must have assigned_class + sex_class alias
schema_sniff <- names(fread(V3_CSV, nrows = 1))
need <- c("gene", "sex_class", "assigned_class",
          "posterior_P_F_only", "logFC_M", "logFC_F", "interaction_padj")
miss_schema <- setdiff(need, schema_sniff)
if (length(miss_schema) > 0) {
  log_line(sprintf("v3 csv missing schema cols: %s",
                   paste(miss_schema, collapse = ", ")))
  stop("Module 09 schema-sniff failed.")
}
log_line("  v3 csv schema sniff OK.")

# Verify the 27a patch is present (line 942 fallback)
patch_grep <- system2("grep", c("-c", "sex_v3_path",
                                file.path(BASE,
                                          "RNA-seq/27a_assemble_evidence_atlas.R")),
                      stdout = TRUE, stderr = FALSE)
patch_hits <- suppressWarnings(as.integer(patch_grep))
log_line(sprintf("  27a v3 fallback patch hits = %s",
                 ifelse(is.na(patch_hits), "?", patch_hits)))
if (is.na(patch_hits) || patch_hits < 1) {
  log_line("27a v3 fallback patch NOT detected in canonical script.")
  stop("Module 09 abort: 27a patch missing.")
}

# ---------------------------------------------------------------------------
# 2) Build sbatch commands (dependency chain)
# ---------------------------------------------------------------------------
log_line("[2] Building sbatch dependency chain ...")

# Step 1: 27a + 27b (atlas assemble + presets benchmark)
cmd_27 <- sprintf("sbatch --parsable %s", shQuote(RUN_27AB))
log_line(sprintf("  cmd_27 = %s", cmd_27))

# Step 2: 75 + 217 chain, depends on 27a
mk_cmd_75 <- function(jid_27) {
  sprintf("sbatch --parsable --dependency=afterok:%s %s",
          jid_27, shQuote(RUN_75_217))
}

# ---------------------------------------------------------------------------
# 3) Submit (if not DRY_RUN)
# ---------------------------------------------------------------------------
dry_run <- as.logical(Sys.getenv("DRY_RUN", "TRUE"))
log_line(sprintf("[3] DRY_RUN = %s (set DRY_RUN=FALSE to actually submit)",
                 dry_run))

if (dry_run) {
  log_line("DRY_RUN enabled -- emitting commands but NOT submitting.")
  log_line(sprintf("Would submit (1/2): %s", cmd_27))
  log_line(sprintf("Would submit (2/2): %s",
                   mk_cmd_75("<JOBID_FROM_STEP_1>")))
  log_line("Module 09 (dry run) complete.")
  # R5 Issue 2 fix: sessionInfo dump even on dry-run exit
  dump_session_info(IDIR, "09_dryrun")
  invisible(quit(save = "no", status = 0))
}

# Real submission path
log_line("Submitting Step 1: 27a + 27b ...")
jid_27 <- tryCatch(system(cmd_27, intern = TRUE),
                   error = function(e) {
                     log_line(sprintf("sbatch 27a FAILED: %s",
                                      conditionMessage(e)))
                     stop(e)
                   })
log_line(sprintf("  -> 27a job id: %s", jid_27))

cmd_75 <- mk_cmd_75(jid_27)
log_line(sprintf("Submitting Step 2: 75 + 217 (depends on %s) ...", jid_27))
jid_75 <- tryCatch(system(cmd_75, intern = TRUE),
                   error = function(e) {
                     log_line(sprintf("sbatch 75/217 FAILED: %s",
                                      conditionMessage(e)))
                     stop(e)
                   })
log_line(sprintf("  -> 75+217 job id: %s", jid_75))

# ---------------------------------------------------------------------------
# 4) Persist job-id chain for downstream agents
# ---------------------------------------------------------------------------
chain_f <- file.path(SEXV3, "atlas_refresh_chain.tsv")
write_atomic_csv(data.table(step = c("27a_27b", "75_217"),
                            job_id = c(jid_27, jid_75),
                            submitted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                 chain_f, sep = "\t")
log_line(sprintf("[4] Wrote (atomic) sbatch chain manifest: %s", chain_f))

# R5 Issue 2 fix: sessionInfo dump for per-module reproducibility audit
dump_session_info(IDIR, "09")

log_line("Module 09 submission complete. Monitor via squeue.")
