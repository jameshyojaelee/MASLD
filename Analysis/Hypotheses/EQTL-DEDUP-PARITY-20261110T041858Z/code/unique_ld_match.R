# One-to-one allelic matching of summary rows to LD reference records.
# Replaces get_ld_per_locus() only in the EQTL-DEDUP refits (see PRESPEC.md).
# Differences from finemapping_functions.R::get_ld_per_locus():
#   * a (position, unordered allele pair) key must occur in exactly one BIM
#     record and one summary row; all other rows are rejected and logged;
#   * returned rows carry event_id = CHR:POS:A:B (sorted alleles).
# Block loading, orientation, exact-duplicate removal and NA pruning are as before.
suppressPackageStartupMessages({library(data.table); library(Matrix)})

.pair_key <- function(pos, a, b) paste(pos, pmin(a, b), pmax(a, b), sep = ":")

.write_identity_losses <- function(LOCUS, CHR, tab) {
  out <- Sys.getenv("IDENTITY_LOSS_DIR", unset = "")
  if (!nzchar(out) || !nrow(tab)) return(invisible(NULL))
  assay <- get0("gwas_name", envir = globalenv(), ifnotfound = "eqtl")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tab <- cbind(data.table(assay = assay, locus = LOCUS, chr = CHR), tab)
  f <- file.path(out, paste0(assay, "_chr", CHR, "_identity_losses.tsv"))
  fwrite(tab, f, sep = "\t", append = file.exists(f))
}

get_ld_per_locus_unique <- function(ss_per_locus, LOCUS, CHR, START, END, ancestry = "EUR") {
  path_to_LD <- find_LD_block(CHR, START, END, ancestry)
  if (is.null(path_to_LD)) return(NULL)
  bim <- data.table(); mats <- list()
  for (file in path_to_LD) {
    bim_file <- paste0(file, ".bim"); ld_file <- paste0(file, ".ld")
    if (!file.exists(bim_file) || !file.exists(ld_file)) {
      warning(paste("LD files not found:", file)); next
    }
    bim <- rbind(bim, fread(bim_file)); mats[[length(mats) + 1L]] <- as.matrix(fread(ld_file))
  }
  if (nrow(bim) == 0) return(NULL)
  # Dense block-diagonal: same values as bdiag() (src/perf/get_ld_per_locus_cached.R).
  if (length(mats) == 1L) ld <- mats[[1L]] else {
    tot <- sum(vapply(mats, nrow, integer(1))); ld <- matrix(0, tot, tot); off <- 0L
    for (m in mats) { ix <- off + seq_len(nrow(m)); ld[ix, ix] <- m; off <- off + nrow(m) }
  }
  rm(mats)
  setnames(bim, c("chr", "rsid", "dk", "pos", "alt", "ref"))
  bim[, SNP := paste(chr, pos, alt, ref, sep = ":")]
  keep_bim <- !duplicated(bim$SNP)              # exact duplicates, as before
  if (!all(keep_bim)) { ld <- ld[keep_bim, keep_bim]; bim <- bim[keep_bim] }
  bim[, key := .pair_key(pos, alt, ref)]

  ss <- as.data.table(ss_per_locus)
  ss[, row_id := .I]
  ss[, key := .pair_key(position, allele1, allele2)]
  ss[, n_ss := .N, by = key]
  n_bim_tab <- bim[, .(n_bim = .N), by = key]
  ss[n_bim_tab, n_bim := i.n_bim, on = "key"]
  ss[is.na(n_bim), n_bim := 0L]
  rejected <- ss[n_bim > 1L | (n_bim == 1L & n_ss > 1L)]
  if (nrow(rejected)) {
    rejected[, reason := fifelse(n_bim > 1L, "multiple_reference_records", "multiple_summary_rows")]
    rejected[, p := 2 * pnorm(-abs(beta / standard_error))]
    .write_identity_losses(LOCUS, CHR, rejected[, .(position, allele1, allele2, beta,
      standard_error, p, n_bim, n_ss, reason)])
  }
  ok <- ss[n_bim == 1L & n_ss == 1L]
  bim[, ss_row := ok$row_id[match(key, ok$key)]]
  sel <- which(!is.na(bim$ss_row))
  cat(paste0("  Unique allelic matches ", length(sel), " / ", nrow(bim),
             " LD variants; rejected ambiguous summary rows ", nrow(rejected), "\n"))
  if (length(sel) < 10) return(NULL)
  bim_f <- bim[sel]; ld_f <- as.matrix(ld[sel, sel])
  out <- as.data.frame(ss_per_locus)[bim_f$ss_row, , drop = FALSE]
  flip <- out$allele1 == bim_f$ref & out$allele2 == bim_f$alt &
          !(out$allele1 == bim_f$alt & out$allele2 == bim_f$ref)
  out$beta[flip] <- -out$beta[flip]
  out$allele1 <- bim_f$alt; out$allele2 <- bim_f$ref
  out$SNP <- bim_f$SNP
  out$event_id <- paste(CHR, bim_f$pos, pmin(bim_f$alt, bim_f$ref), pmax(bim_f$alt, bim_f$ref), sep = ":")
  stopifnot(!anyDuplicated(out$event_id), nrow(out) == nrow(ld_f))

  if (sum(is.na(ld_f)) != 0) {                   # same pruning as the original
    n_before <- nrow(ld_f)
    for (iter in seq_len(100)) {
      na_per_var <- colSums(is.na(ld_f))
      if (max(na_per_var) == 0) break
      bad <- na_per_var > (nrow(ld_f) * 0.5)
      if (sum(bad) == 0) bad <- na_per_var > 0
      ld_f <- ld_f[!bad, !bad, drop = FALSE]; out <- out[!bad, , drop = FALSE]
      if (nrow(ld_f) < 10) break
    }
    cat("  LD NA pruning removed", n_before - nrow(ld_f), "variants\n")
    if (nrow(ld_f) < 10) return(NULL)
  }
  rownames(out) <- NULL
  list(out, ld_f)
}
