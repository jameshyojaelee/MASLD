# 03_lead_z.R -- H80 job 3: observed lead eQTL |z| per production eQTL fit (readout, gate fix 4).
# lead |z| = max |Beta / SE| in the Broadaway sumstats over exactly the SNPs of the fit
# (colnames of lbf_variable, CHR:POS), matched by ENSG and POS; the 08 row filter is applied first.
# Also reads converged and the credible-set count from every rds (rds-based count for the 13,704).
args0 <- commandArgs(trailingOnly = FALSE)
SCRIPTS <- dirname(normalizePath(sub("^--file=", "", args0[grep("^--file=", args0)])))
Sys.setenv(H80_SCRIPTS = SCRIPTS)
source(file.path(SCRIPTS, "h80_lib.R"))
A <- parse_args(list(setup = NA, out = NA, workers = "4", seed = "20260923", only_chr = "", skip_preamble = "FALSE"))
H80 <- dirname(SCRIPTS); SETUP <- A$setup; OUT <- A$out
WORKERS <- as.integer(A$workers)
stopifnot(as.integer(A$seed) == H80_SEED)
setDTthreads(WORKERS)
log_msg("H80 readout start; setup = ", SETUP, "; out = ", OUT)
st <- if (identical(A$skip_preamble, "TRUE")) "PASS" else job_preamble(H80, SETUP, OUT, cores = WORKERS)
if (!identical(st, "PASS")) {
  writeLines(paste0("SKIPPED: ", st), file.path(OUT, "SKIPPED"))
  log_msg("readout not run: ", st)
  quit(save = "no", status = 0)
}

chrs <- if (nzchar(A$only_chr)) as.integer(strsplit(A$only_chr, ",")[[1]]) else 1:22
gx <- fread(file.path(SETUP, "gene_state_by_fit.tsv"))
out_rows <- list()
for (cc in chrs) {
  rds <- sort(Sys.glob(file.path(EQTL_FITS, paste0("chr", cc), "ENSG*_susie.rds")), method = "radix")
  info <- rbindlist(mclapply(rds, function(p) {
    so <- readRDS(p)
    ids <- colnames(so$lbf_variable)
    data.table(ensembl = sub("_susie\\.rds$", "", basename(p)), chr = cc, rds = p, nE = length(ids), nE_unique_pos = uniqueN(sub("^[0-9]+:", "", ids)),
               rds_converged = isTRUE(so$converged), rds_n_cs = length(so$sets$cs),
               pos = list(as.integer(sub("^[0-9]+:", "", ids))))
  }, mc.cores = WORKERS, mc.preschedule = TRUE))
  f <- file.path(EQTL_DIR, paste0("chr", cc, "_marginal_summary_results.tsv"))
  dt <- fread(f, select = c("POS", "Beta", "SE", "ENSG"), nThread = WORKERS)
  dt <- dt[!is.na(POS) & !is.na(Beta) & !is.na(SE) & SE > 0 & ENSG != ""]   # 08 filter (08:75)
  dt[, absz := abs(Beta / SE)]
  dt <- dt[ENSG %in% info$ensembl, .(ENSG, POS, absz)]
  setkey(dt, ENSG, POS)
  fp <- info[, .(POS = unlist(pos)), by = ensembl]
  fp <- unique(fp)
  m <- dt[fp, on = .(ENSG = ensembl, POS), nomatch = NULL]
  per_pos <- m[, .(absz = max(absz), n_rows = .N), by = .(ENSG, POS)]
  g <- per_pos[, .(lead_absz = max(absz), n_pos_matched = .N, n_dup_pos = sum(n_rows > 1),
                   lead_absz_nodup = if (any(n_rows == 1)) max(absz[n_rows == 1]) else NA_real_), by = ENSG]
  info[, pos := NULL]
  info <- merge(info, g, by.x = "ensembl", by.y = "ENSG", all.x = TRUE)
  out_rows[[length(out_rows) + 1L]] <- info
  log_msg("chr", cc, ": ", nrow(info), " fits; Broadaway rows used ", nrow(m), "; genes with a duplicated position ",
          sum(info$n_dup_pos > 0, na.rm = TRUE), "; fits with a duplicated CHR:POS id ", sum(info$nE_unique_pos < info$nE), "; fits with a position absent from the gene's Broadaway rows ", sum(info$n_pos_matched < info$nE_unique_pos, na.rm = TRUE))
  rm(dt, m, per_pos, fp); gc(FALSE)
}
L <- rbindlist(out_rows)
L <- merge(L, gx[, .(ensembl, susie_state_t12, fit)], by = "ensembl", all.x = TRUE)
fwrite(L, file.path(OUT, "lead_z.tsv"), sep = "\t")
log_msg("lead |z| computed for ", nrow(L), " fits")

# rds-based count for the 13,704 frame genes (prespec n-by-group audit note).
fr <- L[susie_state_t12 == "no_susie_result" & fit == "cs0"]
fr[, rds_cs0 := rds_converged & rds_n_cs == 0]
dis <- fr[rds_cs0 == FALSE]
fwrite(dis, file.path(OUT, "cs0_frame_rds_disagreements.tsv"), sep = "\t")
cnt <- data.table(check = c("frame_cs0_last_summary_row", "frame_cs0_rds_based", "frame_rds_disagreements",
                            "cs_ge1_summary", "fits_with_dup_broadaway_positions", "fits_with_duplicated_fit_ids", "fits_with_unmatched_positions"),
                  value = c(nrow(fr), sum(fr$rds_cs0), nrow(dis), sum(L$fit == "cs_ge1", na.rm = TRUE),
                            sum(L$n_dup_pos > 0, na.rm = TRUE), sum(L$nE_unique_pos < L$nE), sum(L$n_pos_matched < L$nE_unique_pos, na.rm = TRUE)))
fwrite(cnt, file.path(OUT, "readout_counts.tsv"), sep = "\t")
for (i in seq_len(nrow(cnt))) log_msg(cnt$check[i], " = ", cnt$value[i])

used <- c(L$rds, file.path(EQTL_DIR, paste0("chr", chrs, "_marginal_summary_results.tsv")))
log_msg("hashing ", length(used), " rds and Broadaway files opened")
hh <- sha256_files(used, WORKERS)
fwrite(hh, file.path(OUT, "opened_files_sha256.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_R.txt"))
writeLines("DONE", file.path(OUT, "DONE"))
log_msg("readout done")
