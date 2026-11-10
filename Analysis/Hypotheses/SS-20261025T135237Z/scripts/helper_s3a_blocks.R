# SS S3a parallel block helper (deviation D-S3a-1). Computes named checkpoint blocks of the running S3a job in a
# separate Slurm job, so the main job reuses them instead of fitting them serially. Runs the unchanged
# s3a_replication.R setup (same inputs, seeds and draws), then fits only the requested blocks with the same functions.
# The main job's run_blocks() reuses a block only if ids and draw key are identical, so a mismatched block is refused
# there; this helper checks the draws against the main run before writing anything.
# Usage: Rscript helper_s3a_blocks.R <main S3a out dir> <helper out dir> <tag:block>[,<tag:block>...]
#   tags: null, observed, bootA, bootB
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 3L)
MAIN <- normalizePath(a[1]); HOUT <- a[2]
TARGET <- do.call(rbind, lapply(strsplit(strsplit(a[3], ",")[[1]], ":"), function(x) data.frame(tag = x[1], b = as.integer(x[2]))))
SRC <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_replication.R"
ra <- jsonlite::read_json(file.path(MAIN, "run_args.json"), simplifyVector = TRUE)
commandArgs <- function(trailingOnly = FALSE) c("--input", ra$input, "--stage", ra$stage, "--out", HOUT, "--ilr", ra$ilr,
                                                "--block-size", ra$`block-size`)
ex <- parse(SRC); txt <- vapply(ex, function(e) paste(deparse(e), collapse = " "), "")
stop_at <- grep("^stopifnot\\(length\\(res\\) == length\\(jobs\\)\\)", txt); stopifnot(length(stop_at) == 1L)
skip <- grepl("^RB <- function", txt) | grepl("^NULLS <- ", txt) | grepl("^rm\\(null_res\\)", txt) |
  grepl("^invisible\\(gc\\(\\)\\)", txt) | grepl("^fwrite\\(cbind\\(draw = seq_len\\(B_PERM\\), NULLS", txt) |
  grepl("^log_msg\\(\"Null fits done\"", txt)
skip[seq_along(skip) >= stop_at] <- FALSE; stopifnot(sum(skip) == 6L)
CKM <- file.path(MAIN, "checkpoints"); stopifnot(dir.exists(CKM))
checked <- FALSE
check_draws <- function() {  # helper draws and inputs must equal the main run's
  stopifnot(identical(fread(file.path(MAIN, "input_sha256.tsv"))$sha256, SHA_IN$sha256),
            identical(fread(file.path(MAIN, "script_sha256.tsv"))$sha256, SHA_SC$sha256))
  pm <- as.matrix(fread(file.path(MAIN, "perm_stage_codes.tsv.gz"))[, -1])
  stopifnot(identical(dim(pm), dim(perm)), identical(colnames(pm), colnames(perm)), all(pm == perm))
  z <- readRDS(file.path(CKM, "null_block001.rds")); stopifnot(identical(z$key, perm[z$ids, , drop = FALSE]))
  bi <- fread(file.path(MAIN, "boot_indices.tsv.gz"))
  stopifnot(identical(as.integer(bi[set == "bin_x_sex"]$index), as.integer(t(bootA))),
            identical(as.integer(bi[set == "stage4_x_sex"]$index), as.integer(t(bootB))))
  log_msg("Helper: inputs, scripts and draws match", MAIN); checked <<- TRUE
}
RB <- function(tag, ids, fun, key) {
  if (!checked) check_draws()
  bl <- split(ids, ceiling(seq_along(ids) / BLOCK))
  for (b in TARGET$b[TARGET$tag == tag]) {
    f <- file.path(CKM, sprintf("%s_block%03d.rds", tag, b))
    if (file.exists(f)) { log_msg("Helper: exists, skipped", basename(f)); next }
    r <- mclapply(bl[[b]], fun, mc.cores = NCORES, mc.preschedule = FALSE)
    if (any(vapply(r, inherits, TRUE, "try-error"))) { print(r[vapply(r, inherits, TRUE, "try-error")][1]); stop(tag, " fit failed") }
    tmp <- paste0(f, ".helper", Sys.getenv("SLURM_JOB_ID", "x"), ".tmp")
    saveRDS(list(ids = bl[[b]], key = key(bl[[b]]), res = r), tmp)
    if (file.exists(f)) { file.remove(tmp); log_msg("Helper: main wrote it first", basename(f)) } else {
      file.rename(tmp, f); log_msg("Helper: wrote", basename(f)) }
  }
  NULL
}
for (i in seq_len(stop_at - 1L)) if (!skip[i]) eval(ex[[i]], envir = globalenv())
log_msg("Helper done:", a[3])
