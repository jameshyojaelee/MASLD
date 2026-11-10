# SS shared helpers for the S3 scripts (sourced, not run). Argument parsing (--input <S2 dir> --stage <stage tsv>
# --out <new dir>, plus named extras), refusal of an existing output directory, sha256 checks of frozen/ against
# frozen/SHA256SUMS and of extra inputs against hashes written in the addenda, the stage join (CRLF stripped) and
# sessionInfo. The stage file is read only here; S2 files carry no stage (S2_CONTRACT.md).
suppressPackageStartupMessages({library(data.table); library(digest)})
SS_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SS_RUN <- file.path(SS_ROOT, "Analysis/Hypotheses/SS-20261025T135237Z")
SS_FROZEN <- file.path(SS_RUN, "frozen")

ss_args <- function(extra = character(0)) {
  a <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(a) %% 2 == 0, all(startsWith(a[c(TRUE, FALSE)], "--")))
  v <- setNames(a[c(FALSE, TRUE)], sub("^--", "", a[c(TRUE, FALSE)]))
  need <- c("input", "stage", "out")
  miss <- setdiff(need, names(v)); if (length(miss)) stop("missing argument(s): ", paste(miss, collapse = ", "))
  bad <- setdiff(names(v), c(need, extra)); if (length(bad)) stop("unknown argument(s): ", paste(bad, collapse = ", "))
  if (dir.exists(v[["out"]])) stop("output directory exists; refusing: ", v[["out"]])
  dir.create(v[["out"]], recursive = TRUE)
  as.list(v)
}

sha256_file <- function(p) digest(file = p, algo = "sha256")

# rel: paths relative to frozen/ as written in SHA256SUMS ("./h51/..." or "h51/..."). Returns absolute paths.
ss_frozen <- function(rel) {
  L <- readLines(file.path(SS_FROZEN, "SHA256SUMS"))
  S <- data.table(sha = sub(" .*$", "", L), path = sub("^\\./", "", sub("^[0-9a-f]{64} +\\*?", "", L)))
  vapply(rel, function(r) {
    r <- sub("^\\./", "", r); e <- S$sha[S$path == r]
    if (length(e) != 1) stop("not in frozen/SHA256SUMS: ", r)
    f <- file.path(SS_FROZEN, r); o <- sha256_file(f)
    if (!identical(o, e)) stop("sha256 mismatch for frozen/", r, ": ", o, " != ", e)
    f
  }, character(1), USE.NAMES = FALSE)
}

# Non-frozen inputs whose hashes are fixed in an addendum. h: named character vector path -> sha256.
ss_check_sha <- function(h) {
  for (p in names(h)) { o <- sha256_file(p); if (!identical(o, h[[p]])) stop("sha256 mismatch: ", p, " ", o, " != ", h[[p]]) }
  invisible(names(h))
}

# Stage file: sample_accession, sample_alias, nash_crn_kleiner_fibrosis_stage (CRLF). Codes 1-4 in SteatoSITE
# (1 = F0-F1); reproduction mocks may carry 0-4. NA = unstaged.
ss_read_stage <- function(path) {
  x <- gsub("\r", "", readLines(path)); s <- fread(text = paste(x, collapse = "\n"), na.strings = c("NA", ""))
  stopifnot(all(c("sample_alias", "nash_crn_kleiner_fibrosis_stage") %in% names(s)), !anyDuplicated(s$sample_alias))
  s <- data.table(sample_alias = as.character(s$sample_alias), stage_num = suppressWarnings(as.numeric(s$nash_crn_kleiner_fibrosis_stage)))
  stopifnot(all(is.na(s$stage_num) | s$stage_num %in% 0:4))
  s
}

# Included units joined to stage. Returns the units table (S2 units.tsv, include == TRUE) with stage_num, bin and
# stage4 added, in the column order of dge_included.rds.
ss_units <- function(input, stage_path) {
  u <- fread(file.path(input, "units.tsv")); u <- u[include %in% c(TRUE, "TRUE")]
  st <- ss_read_stage(stage_path)
  u <- merge(u, st, by = "sample_alias", all.x = TRUE, sort = FALSE)
  u[, bin := factor(ifelse(stage_num <= 1, "F0F1", ifelse(stage_num == 2, "F2", "F3F4")), levels = c("F0F1", "F2", "F3F4"))]
  u[, stage4 := factor(ifelse(stage_num <= 1, "S1", paste0("S", stage_num)), levels = paste0("S", 1:4))]
  u[, batch := factor(batch)]
  u
}

ss_finish <- function(out) {
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  cat("DONE ->", out, "\n")
}
