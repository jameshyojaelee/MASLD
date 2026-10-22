# 01_arm_a.R -- H80 job 1: Arm A, matched LD, 200 windows x 10 replicates x 19 fits.
args0 <- commandArgs(trailingOnly = FALSE)
SCRIPTS <- dirname(normalizePath(sub("^--file=", "", args0[grep("^--file=", args0)])))
Sys.setenv(H80_SCRIPTS = SCRIPTS)
source(file.path(SCRIPTS, "h80_lib.R"))
A <- parse_args(list(setup = NA, out = NA, workers = "16", seed = "20260923", only_rows = "", n_rep = "10",
                     null_only = "FALSE", skip_preamble = "FALSE"))
H80 <- dirname(SCRIPTS); SETUP <- A$setup; OUT <- A$out
WORKERS <- as.integer(A$workers)
NREP <- as.integer(A$n_rep)
NULL_ONLY <- identical(A$null_only, "TRUE")   # debug mode: lambda = 0 fits only (no outcome computed)
stopifnot(as.integer(A$seed) == H80_SEED)
dir.create(file.path(OUT, "windows"), showWarnings = FALSE)
setDTthreads(1)
log_msg("H80 Arm A start; setup = ", SETUP, "; out = ", OUT, "; workers = ", WORKERS)

st <- if (identical(A$skip_preamble, "TRUE")) "PASS" else job_preamble(H80, SETUP, OUT, cores = 4L)
if (!identical(st, "PASS")) {
  writeLines(paste0("SKIPPED: ", st), file.path(OUT, "SKIPPED"))
  log_msg("Arm A not run: ", st)
  quit(save = "no", status = 0)
}
W <- fread(file.path(SETUP, "windows.tsv"))
streams <- readRDS(file.path(SETUP, "streams.rds"))
setup_hashes <- sha256_files(c(file.path(SETUP, c("windows.tsv", "streams.rds", "STATUS")),
                               file.path(SETUP, "eqtl_ss", paste0(unique(W$ensembl), ".rds")),
                               file.path(SETUP, "gwas_ss", paste0("A", W$row, ".rds"))), 4L)
fwrite(setup_hashes, file.path(OUT, "setup_inputs_sha256.tsv"), sep = "\t")
rows <- if (nzchar(A$only_rows)) as.integer(strsplit(A$only_rows, ",")[[1]]) else W$row
eqtl_lams <- if (NULL_ONLY) 0 else A_EQTL_LAMBDAS

run_window_A <- function(r) {
  w <- W[row == r]
  set_stream(streams, w$stream)
  t0 <- proc.time()[["elapsed"]]
  gene_dt <- readRDS(file.path(SETUP, "eqtl_ss", paste0(w$ensembl, ".rds")))
  ss <- eqtl_ss_from_rows(gene_dt); se <- eqtl_window(gene_dt)
  LE <- ld_load_fast(ss, w$chr, se[1], se[2], PANEL_DIRS[["EUR"]])
  if (is.null(LE)) stop("eQTL-side LD load returned NULL")
  idsE <- snp_ids_of(LE$ss)
  rds_ids <- colnames(readRDS(w$rds)$lbf_variable)
  gss <- readRDS(file.path(SETUP, "gwas_ss", paste0("A", r, ".rds")))
  LG <- ld_load_fast(gss, w$chr, min(gss$position), max(gss$position), w$gwas_panel)
  if (is.null(LG)) stop("GWAS-side LD load returned NULL")
  idsG <- snp_ids_of(LG$ss)
  inE <- match(idsG, idsE)
  EG <- which(!is.na(inE))                       # indices in G of SNPs also in E
  if (length(idsG) < MIN_TRIPLE_SNPS || length(EG) < MIN_TRIPLE_SNPS) stop("fewer than 50 GWAS-side SNPs after LD load")
  t_load <- proc.time()[["elapsed"]] - t0
  eigE <- block_eigen(LE$R, LE$blk); eigG <- block_eigen(LG$R, LG$blk)
  RrE <- reg_matrix(LE$R, idsE); RrG <- reg_matrix(LG$R, idsG)
  t_eig <- proc.time()[["elapsed"]] - t0 - t_load
  pE <- length(idsE); pG <- length(idsG)
  fits <- list(); cells <- list(); reps <- list()
  sum_zz0 <- 0; noise_c <- numeric(0)
  for (rep in seq_len(NREP)) {
    # all draws of the replicate first, in the fixed order c, k, eps_E, eps_G
    cG <- EG[sample.int(length(EG), 1L)]
    cE <- inE[cG]
    r2E <- LE$R[cE, inE[EG]]^2
    r2G <- LG$R[cG, EG]^2
    pool <- EG[r2E < DISTINCT_R2_MAX & r2G < DISTINCT_R2_MAX]
    kG <- if (length(pool) >= DISTINCT_POOL_MIN) pool[sample.int(length(pool), 1L)] else NA_integer_
    epsE <- rnorm(pE); epsG <- rnorm(pG)
    seed_after_draws <- get(".Random.seed", envir = .GlobalEnv)
    nzE <- noise_from(eigE, epsE); nzG <- noise_from(eigG, epsG)
    sum_zz0 <- sum_zz0 + sum(nzE^2)
    noise_c <- c(noise_c, nzE[cE])
    fE <- list(); fGs <- list(); fGd <- list()
    for (lam in eqtl_lams) {
      z <- lam * LE$R[, cE] + nzE; names(z) <- idsE
      fE[[as.character(lam)]] <- fr <- fit_eqtl(z, RrE)
      fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep), fit_row(w, "eqtl", lam, "shared", z, cE, fr))
    }
    if (!NULL_ONLY) {
      for (lam in A_GWAS_SHARED) {
        z <- lam * LG$R[, cG] + nzG; names(z) <- idsG
        fGs[[as.character(lam)]] <- fr <- fit_gwas(z, RrG, w$N_eff)
        fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep), fit_row(w, "gwas", lam, "shared", z, cG, fr))
      }
      if (!is.na(kG)) for (lam in A_GWAS_DISTINCT) {
        z <- lam * LG$R[, kG] + nzG; names(z) <- idsG
        fGd[[as.character(lam)]] <- fr <- fit_gwas(z, RrG, w$N_eff)
        fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep), fit_row(w, "gwas", lam, "distinct", z, kG, fr))
      }
      cl <- rbind(data.table(cell_set = "E", lam_e = A_ESWEEP, lam_g = 10, causal = "shared"),
                  data.table(cell_set = "G", lam_e = 10, lam_g = A_GSWEEP, causal = "shared"),
                  data.table(cell_set = "D", lam_e = sapply(A_DISTINCT_CELLS, `[`, 1), lam_g = sapply(A_DISTINCT_CELLS, `[`, 2),
                             causal = "distinct"))
      for (j in seq_len(nrow(cl))) {
        fe <- fE[[as.character(cl$lam_e[j])]]
        fg <- if (cl$causal[j] == "shared") fGs[[as.character(cl$lam_g[j])]] else fGd[[as.character(cl$lam_g[j])]]
        stt <- if (is.null(fg)) list(state = "distinct_not_applicable", pp4 = NA_real_, pp3 = NA_real_, n_pairs = NA_integer_,
                                     n_eligible = NA_integer_, elig_match = NA, coloc_drew = NA, coloc_restored = NA)
               else coloc_state(fe, fg)
        cells[[length(cells) + 1L]] <- cbind(data.table(row = r, rep = rep), cl[j], as.data.table(stt))
      }
    } else {
      for (lam in eqtl_lams) {
        cells[[length(cells) + 1L]] <- data.table(row = r, rep = rep, cell_set = "E", lam_e = lam, lam_g = NA_real_, causal = "null_only",
                                                  state = NA_character_)
      }
    }
    stream_unchanged <- identical(get(".Random.seed", envir = .GlobalEnv), seed_after_draws)
    reps[[length(reps) + 1L]] <- data.table(row = r, rep = rep, c_id = idsG[cG], k_id = if (is.na(kG)) NA_character_ else idsG[kG],
                                            pool_size = length(pool), epsE_hash = digest(epsE, algo = "xxhash64"),
                                            epsG_hash = digest(epsG, algo = "xxhash64"), noise_at_c = nzE[cE],
                                            stream_unchanged_by_fits = stream_unchanged,
                                            kind_after = RNGkind()[1])
    rm(fE, fGs, fGd); gc(FALSE)
  }
  T_w <- (sum_zz0 - NREP * eigE$trS) / sqrt(2 * NREP * eigE$trS2)
  diag <- data.table(row = r, ensembl = w$ensembl, stratum = w$stratum, ancestry = w$ancestry,
                     pE = pE, pG = pG, pEG = length(EG), ids_E_identical_rds = identical(idsE, rds_ids),
                     n_dup_ids_E = sum(duplicated(idsE)), na_E = LE$n_na, na_G = LG$n_na,
                     min_eig_E = eigE$min_eig, n_clip_E = eigE$n_neg, n_blocks_E = eigE$n_blocks,
                     min_eig_G = eigG$min_eig, n_clip_G = eigG$n_neg, n_blocks_G = eigG$n_blocks,
                     trS_E = eigE$trS, trS2_E = eigE$trS2, eff_df_E = NREP * eigE$trS^2 / eigE$trS2,
                     T_w = T_w, mean_noise_c = mean(noise_c), nc3_T_pass = abs(T_w) <= 4.5,
                     nc3_causal_pass = abs(mean(noise_c)) <= 4 / sqrt(NREP),
                     t_load_s = t_load, t_eig_s = t_eig, t_total_s = proc.time()[["elapsed"]] - t0,
                     peak_rss_gb = as.numeric(system2("bash", c("-c", shQuote(sprintf("grep VmHWM /proc/%d/status | awk '{print $2}'", Sys.getpid()))), stdout = TRUE)) / 1e6,
                     files_used = paste(c(LE$files_used, LG$files_used), collapse = "|"))
  out <- list(diag = diag, fits = rbindlist(fits, fill = TRUE), cells = rbindlist(cells, fill = TRUE), reps = rbindlist(reps))
  saveRDS(out, file.path(OUT, "windows", paste0("A", r, ".rds")))
  data.table(row = r, ok = TRUE, error = NA_character_, t_total_s = diag$t_total_s)
}

ord <- W[row %in% rows][order(-nE), row]
log_msg("Arm A: ", length(ord), " windows, ", NREP, " replicates, ", if (NULL_ONLY) "NULL-ONLY debug" else "full design")
st_rows <- mclapply(ord, function(r) {
  tryCatch(run_window_A(r), error = function(e) data.table(row = r, ok = FALSE, error = conditionMessage(e), t_total_s = NA_real_))
}, mc.cores = WORKERS, mc.preschedule = FALSE)
st_rows <- rbindlist(lapply(seq_along(st_rows), function(i) if (is.data.table(st_rows[[i]])) st_rows[[i]] else
  data.table(row = ord[i], ok = FALSE, error = paste("worker died:", paste(as.character(st_rows[[i]]), collapse = " ")), t_total_s = NA_real_)))
fwrite(st_rows, file.path(OUT, "window_status.tsv"), sep = "\t")
log_msg("Arm A windows finished: ok ", sum(st_rows$ok), " / ", nrow(st_rows))

fs <- file.path(OUT, "windows", paste0("A", ord, ".rds"))
fs <- fs[file.exists(fs)]
res <- lapply(fs, readRDS)
fwrite(rbindlist(lapply(res, `[[`, "diag")), file.path(OUT, "window_diagnostics.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "fits"), fill = TRUE), file.path(OUT, "fits.tsv.gz"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "cells"), fill = TRUE), file.path(OUT, "cells.tsv.gz"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "reps")), file.path(OUT, "reps.tsv.gz"), sep = "\t")
used <- unique(unlist(strsplit(rbindlist(lapply(res, `[[`, "diag"))$files_used, "|", fixed = TRUE)))
used <- c(paste0(used, ".bim"), paste0(used, ".ld"), unique(W[row %in% rows, rds]))
log_msg("hashing ", length(used), " LD / rds files opened by Arm A")
hh <- sha256_files(used, WORKERS)
fwrite(hh, file.path(OUT, "opened_files_sha256.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_R.txt"))
writeLines("DONE", file.path(OUT, "DONE"))
log_msg("Arm A done")
