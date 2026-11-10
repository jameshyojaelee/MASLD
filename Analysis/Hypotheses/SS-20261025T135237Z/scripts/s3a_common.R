# SteatoSITE S3 shared definitions (forked from H51 h51_common.R; statistics unchanged). Sourced by s3a_replication.R,
# s3a_build_mock.R and s3b_catalog_sidecar.R. Binding text: PRESPEC_SS.md sections 4-5 and PRESPEC_SS_addendum_S3a.md.
suppressMessages({library(data.table); library(edgeR); library(limma); library(parallel)})

SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FROZEN <- file.path(SS, "frozen")
SEED_PERM <- 20260923L
SEED_BOOT <- 20261007L
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
BINS <- c("F0F1", "F2", "F3F4")
CF <- c(F3F4 = "binF3F4", F2 = "binF2")
SEXCF <- "factor(inferred_sex)M"
STEPS <- c(S2mS1 = "F1_to_F2", S3mS2 = "F2_to_F3", S4mS3 = "F3_to_F4")
STEP_LABEL <- c(S2mS1 = "F0-F1 -> F2 (discovery F1 -> F2)", S3mS2 = "F2 -> F3", S4mS3 = "F3 -> F4 (specimen type may differ)")
PCG <- c(COL1A1 = 1, EPCAM = 1, LUM = 1, GNMT = -1, CYP2C19 = -1)
HERO <- c("GNMT", "MAT1A", "CYP2C19")
QUAL <- c("log10_assigned", "unique_pct", "intron_exon_ratio")
z975 <- qnorm(0.975); z80 <- qnorm(0.80)
MDE_SLOPE_K <- z975 + z80                # slope MDE = 2.80 x null SD
MDE_GENE_K <- qnorm(1 - 0.005 / 2) + z80 # per-gene MDE_g = 3.65 x SE_g

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
sha256 <- function(f) sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE))

# --key value / --flag parser; `flags` are switches without a value
parse_args <- function(args, required, flags = character(), defaults = list()) {
  out <- defaults; i <- 1L
  while (i <= length(args)) {
    k <- sub("^--", "", args[i]); if (!startsWith(args[i], "--")) stop("Unexpected argument: ", args[i])
    if (k %in% flags) { out[[k]] <- TRUE; i <- i + 1L } else { out[[k]] <- args[i + 1L]; i <- i + 2L }
  }
  miss <- setdiff(required, names(out)); if (length(miss)) stop("Missing arguments: ", paste0("--", miss, collapse = " "))
  for (f in flags) if (is.null(out[[f]])) out[[f]] <- FALSE
  out
}

new_outdir <- function(out) {
  if (dir.exists(out) || file.exists(out)) stop("Output dir exists, refusing: ", out)
  dir.create(out, recursive = TRUE); normalizePath(out)
}

# sha256sum-format file: "<hash>  <path>"
read_sums <- function(f) { l <- readLines(f); l <- l[nzchar(l)]
  data.table(expected = sub("^([0-9a-f]{64})\\s+\\*?(.*)$", "\\1", l), file = sub("^([0-9a-f]{64})\\s+\\*?(.*)$", "\\2", l)) }

# Every frozen discovery file against frozen/SHA256SUMS; stops on any mismatch
verify_frozen <- function(outdir) {
  s <- read_sums(file.path(FROZEN, "SHA256SUMS"))
  s[, path := file.path(FROZEN, sub("^\\./", "", file))]
  s[, observed := vapply(path, sha256, "")]
  s[, match := expected == observed]
  fwrite(s, file.path(outdir, "frozen_sha256_check.tsv"), sep = "\t")
  if (!all(s$match)) stop("Frozen sha256 mismatch: ", paste(s[match == FALSE, file], collapse = ","))
  log_msg("FROZEN HASH CHECK PASSED:", nrow(s), "files")
}

# PRESPEC_SS.sha256 (written before S3 submission) must exist and match unless explicitly allowed (mock runs only)
verify_prespec <- function(outdir, allow_unhashed) {
  f <- file.path(SS, "PRESPEC_SS.sha256")
  if (!file.exists(f)) {
    if (!allow_unhashed) stop("PRESPEC_SS.sha256 absent; real S3 runs require it (use --allow-unhashed-prespec only for mock runs)")
    log_msg("NOTE: PRESPEC_SS.sha256 absent; allowed for this run (mock/verification)"); return(invisible(FALSE))
  }
  s <- read_sums(f)
  s[, observed := vapply(file.path(SS, sub("^\\./", "", file)), sha256, "")]
  fwrite(s, file.path(outdir, "prespec_sha256_check.tsv"), sep = "\t")
  if (!all(s$expected == s$observed)) stop("PRESPEC sha256 mismatch: ", paste(s[expected != observed, file], collapse = ","))
  log_msg("PRESPEC HASH CHECK PASSED:", nrow(s), "files"); invisible(TRUE)
}

# Release gene annotation (24,115 genes, release order is checked against the input DGEList)
load_ann <- function() {
  a <- unique(fread(file.path(FROZEN, "stage/stage_extension_all_gene_results.tsv"),
                    select = c("gene_id_versioned", "gene_id_base", "gene_name", "chromosome")))
  stopifnot(nrow(a) == 24115L, !anyDuplicated(a$gene_id_versioned))
  a
}

# Stage file (CRLF): code 1 = F0-F1, 2-4 = F2-F4, NA = unstaged
read_stage <- function(path) {
  l <- sub("\r$", "", readLines(path))
  s <- fread(text = paste(l, collapse = "\n"), colClasses = "character")
  stopifnot(all(c("sample_alias", "nash_crn_kleiner_fibrosis_stage") %in% names(s)), !anyDuplicated(s$sample_alias))
  s[, stage_code := suppressWarnings(as.integer(nash_crn_kleiner_fibrosis_stage))]
  stopifnot(all(is.na(s$stage_code) | s$stage_code %in% 1:4))
  s[, .(sample_alias, stage_code)]
}

bin_of_code <- function(x) factor(fifelse(x <= 1L, "F0F1", fifelse(x == 2L, "F2", "F3F4")), levels = BINS)
stage4_of_code <- function(x) { f <- factor(paste0("S", x), levels = paste0("S", 1:4)); droplevels(f) }
add_labels <- function(info, code) {
  info$code <- code; info$bin <- bin_of_code(code); info$stage4 <- stage4_of_code(code); info$stage_num <- as.numeric(code)
  info
}

# Design: batch and sex enter when they have > 1 level in the units fitted (one-level mock batch reproduces H51 G0)
design_of <- function(info, main, extra = character()) {
  tm <- c(if (length(unique(info$batch)) > 1L) "factor(batch)", if (length(unique(info$inferred_sex)) > 1L) "factor(inferred_sex)", extra, main)
  model.matrix(as.formula(paste("~", paste(tm, collapse = " + "))), as.data.frame(info))
}

# voomWithQualityWeights -> lmFit (-> eBayes), limma defaults (as H51 / job 3831)
fit_qw <- function(d, des, ebayes = TRUE) {
  v <- voomWithQualityWeights(d, des, plot = FALSE)
  f <- lmFit(v, des); if (ebayes) eBayes(f) else f
}

# Adjacent-step contrasts on the stage4 coefficients (S1 reference); absent levels give NA
step_coefs <- function(cm) {
  g <- function(n) if (n %in% colnames(cm)) cm[, n] else rep(NA_real_, nrow(cm))
  cbind(S2mS1 = g("stage4S2"), S3mS2 = g("stage4S3") - g("stage4S2"), S4mS3 = g("stage4S4") - g("stage4S3"))
}
step_matrix <- function(des) {
  cn <- colnames(des); C <- matrix(0, length(cn), 0, dimnames = list(cn, NULL))
  add <- function(nm, plus, minus = NULL) if (all(c(plus, minus) %in% cn)) { v <- setNames(rep(0, length(cn)), cn); v[plus] <- 1; if (!is.null(minus)) v[minus] <- -1
    C <<- cbind(C, matrix(v, ncol = 1, dimnames = list(cn, nm))) }
  add("S2mS1", "stage4S2"); add("S3mS2", "stage4S3", "stage4S2"); add("S4mS3", "stage4S4", "stage4S3")
  C
}

# Per-gene table from an eBayes fit
gene_tab <- function(f, cfs) rbindlist(lapply(cfs, function(cf) data.table(gene = rownames(f$coefficients), coef = cf,
  logFC = unname(f$coefficients[, cf]), SE = unname(f$stdev.unscaled[, cf] * sqrt(f$s2.post)), t = unname(f$t[, cf]),
  df = unname(f$df.total), P = unname(f$p.value[, cf]), FDR = p.adjust(f$p.value[, cf], "BH"))))
p_one <- function(t, df, sign) ifelse(sign > 0, pt(t, df, lower.tail = FALSE), pt(t, df, lower.tail = TRUE))

ols_slope <- function(y, x) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2) }
conc <- function(y, x) { ok <- is.finite(x) & is.finite(y); mean(sign(y[ok]) == sign(x[ok])) }

boot_draw <- function(cells) {
  cl <- sort(unique(cells))
  unlist(lapply(cl, function(cc) { i <- which(cells == cc); i[sample.int(length(i), length(i), replace = TRUE)] }), use.names = FALSE)
}

tbl <- function(dt) { dt <- as.data.frame(dt); c(paste0("| ", paste(names(dt), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
  apply(dt, 1, function(r) paste0("| ", paste(trimws(format(r)), collapse = " | "), " |"))) }
fmt <- function(dt, digits = 4) as.data.table(dt)[, lapply(.SD, function(v) if (is.logical(v)) as.character(v) else
  if (is.numeric(v) && any(abs(v - round(v)) > 0, na.rm = TRUE)) signif(v, digits) else v)]

# Checkpointed runner: `ids` in blocks of `block`; each block file stores its ids, its draw content (`key(ids)`) and the
# results. Existing blocks are reused only if ids and key are identical; new blocks are written via .tmp + rename.
# stop_after > 0 exits with status 3 after that many newly written blocks (restart test only).
.rb <- new.env(); .rb$n_new <- 0L
run_blocks <- function(tag, ids, fun, key, ck, block, stop_after = 0L, cores = 1L) {
  bl <- split(ids, ceiling(seq_along(ids) / block)); out <- vector("list", length(bl))
  for (b in seq_along(bl)) {
    f <- file.path(ck, sprintf("%s_block%03d.rds", tag, b)); kb <- key(bl[[b]])
    if (file.exists(f)) { z <- readRDS(f)
      if (!identical(z$ids, bl[[b]]) || !identical(z$key, kb)) stop("Checkpoint ", basename(f), " does not match the current draws; refusing")
      out[[b]] <- z$res; log_msg("Reused", basename(f)); next }
    r <- mclapply(bl[[b]], fun, mc.cores = cores, mc.preschedule = FALSE)
    if (any(vapply(r, inherits, TRUE, "try-error"))) { print(r[vapply(r, inherits, TRUE, "try-error")][1]); stop(tag, " fit failed") }
    saveRDS(list(ids = bl[[b]], key = kb, res = r), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    out[[b]] <- r; .rb$n_new <- .rb$n_new + 1L; log_msg("Wrote", basename(f))
    if (stop_after > 0L && .rb$n_new >= stop_after) { log_msg("STOP after", .rb$n_new, "new blocks (--stop-after-blocks test)"); quit(save = "no", status = 3) }
  }
  res <- do.call(c, out); stopifnot(length(res) == length(ids)); res
}
