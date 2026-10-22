# 02_arm_b_pc1.R -- H80 job 2: Arm B (LD mismatch, 120 EUR windows x 10 replicates x 2 generators)
# and PC-1 (observed strength at the 286 S-pos windows, 10 replicates).
args0 <- commandArgs(trailingOnly = FALSE)
SCRIPTS <- dirname(normalizePath(sub("^--file=", "", args0[grep("^--file=", args0)])))
Sys.setenv(H80_SCRIPTS = SCRIPTS)
source(file.path(SCRIPTS, "h80_lib.R"))
A <- parse_args(list(setup = NA, out = NA, workers = "16", seed = "20260923", only_b = "", only_p = "",
                     n_rep = "10", null_only = "FALSE", skip_preamble = "FALSE"))
H80 <- dirname(SCRIPTS); SETUP <- A$setup; OUT <- A$out
WORKERS <- as.integer(A$workers)
NREP <- as.integer(A$n_rep)
NULL_ONLY <- identical(A$null_only, "TRUE")   # debug mode: every lambda set to 0 (no outcome computed)
stopifnot(as.integer(A$seed) == H80_SEED)
dir.create(file.path(OUT, "windows"), showWarnings = FALSE)
setDTthreads(1)
log_msg("H80 Arm B + PC-1 start; setup = ", SETUP, "; out = ", OUT, "; workers = ", WORKERS)

st <- if (identical(A$skip_preamble, "TRUE")) "PASS" else job_preamble(H80, SETUP, OUT, cores = 4L)
if (!identical(st, "PASS")) {
  writeLines(paste0("SKIPPED: ", st), file.path(OUT, "SKIPPED"))
  log_msg("Arm B / PC-1 not run: ", st)
  quit(save = "no", status = 0)
}
W <- fread(file.path(SETUP, "windows.tsv"))
P1 <- fread(file.path(SETUP, "pc1_windows.tsv"))
streams <- readRDS(file.path(SETUP, "streams.rds"))
setup_hashes <- sha256_files(c(file.path(SETUP, c("windows.tsv", "pc1_windows.tsv", "streams.rds", "STATUS")),
                               file.path(SETUP, "eqtl_ss", paste0(unique(c(W[armB == TRUE, ensembl], P1$ensembl)), ".rds")),
                               file.path(SETUP, "gwas_ss", paste0("A", W[armB == TRUE, row], ".rds")),
                               file.path(SETUP, "gwas_ss", paste0("P", P1[!is.na(n_gwas_rows), row], ".rds"))), 4L)
fwrite(setup_hashes, file.path(OUT, "setup_inputs_sha256.tsv"), sep = "\t")
lam <- function(x) if (NULL_ONLY) 0 * x else x

# Variant key that ignores allele order: chr:pos:min(allele):max(allele), from a bim SNP id chr:pos:alt:ref.
key_parts <- function(snp) {
  p <- tstrsplit(snp, ":", fixed = TRUE)
  list(key = paste(p[[1]], p[[2]], pmin(p[[3]], p[[4]]), pmax(p[[3]], p[[4]]), sep = ":"), alt = p[[3]])
}
# Restrict a polyfun load and a 1kg_eur load to the variants present in both; the 1kg matrix is
# sign-aligned to the polyfun allele coding (rows/columns of a variant whose alt/ref are swapped flip sign).
restrict_pair <- function(P, K) {
  kp <- key_parts(colnames(P$R)); kk <- key_parts(colnames(K$R))
  ip <- which(kp$key %in% kk$key)
  jk <- match(kp$key[ip], kk$key)
  sgn <- ifelse(kp$alt[ip] == kk$alt[jk], 1, -1)
  RP <- P$R[ip, ip, drop = FALSE]
  RK <- K$R[jk, jk, drop = FALSE] * outer(sgn, sgn)
  # sign-alignment check: cosine similarity of the two matrices (positive and large when aligned)
  cosine <- sum(RP * RK) / sqrt(sum(RP^2) * sum(RK^2))
  list(RP = RP, RK = RK, blkP = P$blk[ip], blkK = K$blk[jk], ids = snp_ids_of(P$ss)[ip],
       n_poly = ncol(P$R), n_1kg = ncol(K$R), n_both = length(ip), n_flipped = sum(sgn < 0), cosine = cosine,
       files = c(P$files_used, K$files_used))
}

run_window_B <- function(r) {
  w <- W[row == r]
  set_stream(streams, w$armB_stream)
  t0 <- proc.time()[["elapsed"]]
  gene_dt <- readRDS(file.path(SETUP, "eqtl_ss", paste0(w$ensembl, ".rds")))
  ss <- eqtl_ss_from_rows(gene_dt); se <- eqtl_window(gene_dt)
  gss <- readRDS(file.path(SETUP, "gwas_ss", paste0("A", r, ".rds")))
  PE <- ld_load_fast(ss, w$chr, se[1], se[2], PANEL_DIRS[["EUR"]])
  KE <- ld_load_fast(ss, w$chr, se[1], se[2], PANEL_DIRS[["EUR_1KG"]])
  if (is.null(PE) || is.null(KE)) stop("eQTL-side LD load returned NULL")
  E <- restrict_pair(PE, KE); rm(PE, KE); gc(FALSE)
  PG <- ld_load_fast(gss, w$chr, min(gss$position), max(gss$position), w$gwas_panel)
  KG <- ld_load_fast(gss, w$chr, min(gss$position), max(gss$position), PANEL_DIRS[["EUR_1KG"]])
  if (is.null(PG) || is.null(KG)) stop("GWAS-side LD load returned NULL")
  G <- restrict_pair(PG, KG); rm(PG, KG); gc(FALSE)
  inE <- match(G$ids, E$ids)
  EG <- which(!is.na(inE))
  if (length(G$ids) < MIN_TRIPLE_SNPS || length(EG) < MIN_TRIPLE_SNPS) stop("fewer than 50 shared SNPs after restriction to both panels")
  gc(FALSE)
  eig <- list(matched = list(E = block_eigen(E$RP, E$blkP), G = block_eigen(G$RP, G$blkP)),
              mismatched = list(E = block_eigen(E$RK, E$blkK), G = block_eigen(G$RK, G$blkK)))
  RrE <- reg_matrix(E$RP, E$ids); RrG <- reg_matrix(G$RP, G$ids)
  gen_R <- list(matched = list(E = E$RP, G = G$RP), mismatched = list(E = E$RK, G = G$RK))
  pE <- length(E$ids); pG <- length(G$ids)
  fits <- list(); cells <- list(); reps <- list()
  zz0 <- c(matched = 0, mismatched = 0)
  for (rep in seq_len(NREP)) {
    cG <- EG[sample.int(length(EG), 1L)]
    cE <- inE[cG]
    ok <- E$RP[cE, inE[EG]]^2 < DISTINCT_R2_MAX & G$RP[cG, EG]^2 < DISTINCT_R2_MAX &
          E$RK[cE, inE[EG]]^2 < DISTINCT_R2_MAX & G$RK[cG, EG]^2 < DISTINCT_R2_MAX
    pool <- EG[ok]
    kG <- if (length(pool) >= DISTINCT_POOL_MIN) pool[sample.int(length(pool), 1L)] else NA_integer_
    epsE <- rnorm(pE); epsG <- rnorm(pG)
    seed_after_draws <- get(".Random.seed", envir = .GlobalEnv)
    for (gen in c("matched", "mismatched")) {
      nzE <- noise_from(eig[[gen]]$E, epsE); nzG <- noise_from(eig[[gen]]$G, epsG)
      zz0[[gen]] <- zz0[[gen]] + sum(nzE^2)
      fE <- list()
      for (le in B_EQTL_LAMBDAS) {
        z <- lam(le) * gen_R[[gen]]$E[, cE] + nzE; names(z) <- E$ids
        fE[[as.character(le)]] <- fr <- fit_eqtl(z, RrE)
        fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep, generator = gen), fit_row(w, "eqtl", le, "shared", z, cE, fr))
      }
      z <- lam(10) * gen_R[[gen]]$G[, cG] + nzG; names(z) <- G$ids
      fGs <- fr <- fit_gwas(z, RrG, w$N_eff)
      fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep, generator = gen), fit_row(w, "gwas", 10, "shared", z, cG, fr))
      fGd <- NULL
      if (!is.na(kG)) {
        z <- lam(10) * gen_R[[gen]]$G[, kG] + nzG; names(z) <- G$ids
        fGd <- fr <- fit_gwas(z, RrG, w$N_eff)
        fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep, generator = gen), fit_row(w, "gwas", 10, "distinct", z, kG, fr))
      }
      for (le in B_ESWEEP) {
        cells[[length(cells) + 1L]] <- cbind(data.table(row = r, rep = rep, generator = gen, cell_set = "E", lam_e = le, lam_g = 10,
                                                        causal = "shared"), as.data.table(coloc_state(fE[[as.character(le)]], fGs)))
      }
      stt <- if (is.null(fGd)) list(state = "distinct_not_applicable") else coloc_state(fE[["10"]], fGd)
      cells[[length(cells) + 1L]] <- cbind(data.table(row = r, rep = rep, generator = gen, cell_set = "D", lam_e = 10, lam_g = 10,
                                                      causal = "distinct"), as.data.table(stt))
      rm(fE, fGs, fGd)
    }
    reps[[length(reps) + 1L]] <- data.table(row = r, rep = rep, c_id = G$ids[cG], k_id = if (is.na(kG)) NA_character_ else G$ids[kG],
                                            pool_size = length(pool), epsE_hash = digest(epsE, algo = "xxhash64"),
                                            epsG_hash = digest(epsG, algo = "xxhash64"),
                                            stream_unchanged_by_fits = identical(get(".Random.seed", envir = .GlobalEnv), seed_after_draws),
                                            kind_after = RNGkind()[1])
    gc(FALSE)
  }
  T_w <- sapply(c("matched", "mismatched"), function(gen)
    (zz0[[gen]] - NREP * eig[[gen]]$E$trS) / sqrt(2 * NREP * eig[[gen]]$E$trS2))
  diag <- data.table(arm = "B", row = r, ensembl = w$ensembl, stratum = w$stratum, pE = pE, pG = pG, pEG = length(EG),
                     E_n_poly = E$n_poly, E_n_1kg = E$n_1kg, E_n_both = E$n_both, E_n_flipped = E$n_flipped,
                     G_n_poly = G$n_poly, G_n_1kg = G$n_1kg, G_n_both = G$n_both, G_n_flipped = G$n_flipped,
                     ld_cosine_E = E$cosine, ld_cosine_G = G$cosine,
                     min_eig_E_matched = eig$matched$E$min_eig, n_clip_E_matched = eig$matched$E$n_neg,
                     min_eig_E_mismatched = eig$mismatched$E$min_eig, n_clip_E_mismatched = eig$mismatched$E$n_neg,
                     T_w_matched = T_w[["matched"]], T_w_mismatched = T_w[["mismatched"]],
                     eff_df_E_matched = NREP * eig$matched$E$trS^2 / eig$matched$E$trS2,
                     t_total_s = proc.time()[["elapsed"]] - t0,
                     peak_rss_gb = as.numeric(system2("bash", c("-c", shQuote(sprintf("grep VmHWM /proc/%d/status | awk '{print $2}'", Sys.getpid()))), stdout = TRUE)) / 1e6,
                     files_used = paste(unique(c(E$files, G$files)), collapse = "|"))
  out <- list(diag = diag, fits = rbindlist(fits, fill = TRUE), cells = rbindlist(cells, fill = TRUE), reps = rbindlist(reps))
  saveRDS(out, file.path(OUT, "windows", paste0("B", r, ".rds")))
  data.table(task = paste0("B", r), ok = TRUE, error = NA_character_, t_total_s = diag$t_total_s)
}

run_window_P <- function(r) {
  w <- P1[row == r]
  set_stream(streams, w$stream)
  t0 <- proc.time()[["elapsed"]]
  if (is.na(w$n_gwas_rows)) stop("no GWAS rows saved in setup (rds missing)")
  gene_dt <- readRDS(file.path(SETUP, "eqtl_ss", paste0(w$ensembl, ".rds")))
  ss <- eqtl_ss_from_rows(gene_dt); se <- eqtl_window(gene_dt)
  LE <- ld_load_fast(ss, w$chr, se[1], se[2], PANEL_DIRS[["EUR"]])
  gss <- readRDS(file.path(SETUP, "gwas_ss", paste0("P", r, ".rds")))
  if (nrow(gss) < MIN_TRIPLE_SNPS) stop("fewer than 50 GWAS rows")
  LG <- ld_load_fast(gss, w$chr, min(gss$position), max(gss$position), w$gwas_panel)
  if (is.null(LE) || is.null(LG)) stop("LD load returned NULL")
  idsE <- snp_ids_of(LE$ss); idsG <- snp_ids_of(LG$ss)
  rds_ids <- colnames(readRDS(w$rds)$lbf_variable)
  inE <- match(idsG, idsE); EG <- which(!is.na(inE))
  if (length(idsG) < MIN_TRIPLE_SNPS || length(EG) < MIN_TRIPLE_SNPS) stop("fewer than 50 GWAS-side SNPs after LD load")
  zE_obs <- abs(LE$ss$beta / LE$ss$standard_error)
  zG_obs <- abs(LG$ss$beta / LG$ss$standard_error)
  lam_e <- max(zE_obs); lam_g <- max(zG_obs)
  lead <- which.max(zE_obs)
  cG <- match(idsE[lead], idsG)
  lead_in_G <- !is.na(cG)
  if (!lead_in_G) cG <- EG[which.max(zE_obs[inE[EG]])]
  cE <- inE[cG]
  eigE <- block_eigen(LE$R, LE$blk); eigG <- block_eigen(LG$R, LG$blk)
  RrE <- reg_matrix(LE$R, idsE); RrG <- reg_matrix(LG$R, idsG)
  pE <- length(idsE); pG <- length(idsG)
  fits <- list(); cells <- list(); reps <- list()
  for (rep in seq_len(NREP)) {
    epsE <- rnorm(pE); epsG <- rnorm(pG)
    seed_after_draws <- get(".Random.seed", envir = .GlobalEnv)
    zE <- lam(lam_e) * LE$R[, cE] + noise_from(eigE, epsE); names(zE) <- idsE
    zG <- lam(lam_g) * LG$R[, cG] + noise_from(eigG, epsG); names(zG) <- idsG
    fe <- fit_eqtl(zE, RrE); fg <- fit_gwas(zG, RrG, w$N_eff)
    fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep), fit_row(w, "eqtl", lam_e, "shared", zE, cE, fe))
    fits[[length(fits) + 1L]] <- cbind(data.table(row = r, rep = rep), fit_row(w, "gwas", lam_g, "shared", zG, cG, fg))
    cells[[length(cells) + 1L]] <- cbind(data.table(row = r, rep = rep, lam_e = lam_e, lam_g = lam_g), as.data.table(coloc_state(fe, fg)))
    reps[[length(reps) + 1L]] <- data.table(row = r, rep = rep, epsE_hash = digest(epsE, algo = "xxhash64"),
                                            stream_unchanged_by_fits = identical(get(".Random.seed", envir = .GlobalEnv), seed_after_draws),
                                            kind_after = RNGkind()[1])
    rm(fe, fg); gc(FALSE)
  }
  diag <- data.table(arm = "PC1", row = r, ensembl = w$ensembl, ancestry = w$ancestry, study = w$study, pE = pE, pG = pG,
                     pEG = length(EG), ids_E_identical_rds = identical(idsE, rds_ids), lam_e_obs = lam_e, lam_g_obs = lam_g,
                     lead_eqtl_in_G = lead_in_G, causal_id = idsG[cG], t_total_s = proc.time()[["elapsed"]] - t0,
                     peak_rss_gb = as.numeric(system2("bash", c("-c", shQuote(sprintf("grep VmHWM /proc/%d/status | awk '{print $2}'", Sys.getpid()))), stdout = TRUE)) / 1e6,
                     files_used = paste(unique(c(LE$files_used, LG$files_used)), collapse = "|"))
  out <- list(diag = diag, fits = rbindlist(fits, fill = TRUE), cells = rbindlist(cells, fill = TRUE), reps = rbindlist(reps))
  saveRDS(out, file.path(OUT, "windows", paste0("P", r, ".rds")))
  data.table(task = paste0("P", r), ok = TRUE, error = NA_character_, t_total_s = diag$t_total_s)
}

b_rows <- if (nzchar(A$only_b)) as.integer(strsplit(A$only_b, ",")[[1]]) else W[armB == TRUE][order(-nE), row]
p_rows <- if (nzchar(A$only_p)) as.integer(strsplit(A$only_p, ",")[[1]]) else P1$row
if (identical(A$only_b, "none")) b_rows <- integer(0)
if (identical(A$only_p, "none")) p_rows <- integer(0)
tasks <- c(paste0("B", b_rows), paste0("P", p_rows))
log_msg("Arm B windows ", length(b_rows), "; PC-1 windows ", length(p_rows), "; replicates ", NREP,
        if (NULL_ONLY) "; NULL-ONLY debug" else "")
res <- mclapply(tasks, function(tk) {
  r <- as.integer(substring(tk, 2))
  tryCatch(if (startsWith(tk, "B")) run_window_B(r) else run_window_P(r),
           error = function(e) data.table(task = tk, ok = FALSE, error = conditionMessage(e), t_total_s = NA_real_))
}, mc.cores = WORKERS, mc.preschedule = FALSE)
st_rows <- rbindlist(lapply(seq_along(res), function(i) if (is.data.table(res[[i]])) res[[i]] else
  data.table(task = tasks[i], ok = FALSE, error = paste("worker died:", paste(as.character(res[[i]]), collapse = " ")), t_total_s = NA_real_)))
fwrite(st_rows, file.path(OUT, "window_status.tsv"), sep = "\t")
log_msg("tasks finished: ok ", sum(st_rows$ok), " / ", nrow(st_rows))

collect <- function(prefix) {
  fs <- file.path(OUT, "windows", paste0(prefix, sub("^.", "", tasks[startsWith(tasks, prefix)]), ".rds"))
  fs <- fs[file.exists(fs)]
  lapply(fs, readRDS)
}
for (pf in c("B", "P")) {
  rr <- collect(pf)
  if (length(rr) == 0) next
  tag <- if (pf == "B") "armB" else "pc1"
  fwrite(rbindlist(lapply(rr, `[[`, "diag"), fill = TRUE), file.path(OUT, paste0(tag, "_window_diagnostics.tsv")), sep = "\t")
  fwrite(rbindlist(lapply(rr, `[[`, "fits"), fill = TRUE), file.path(OUT, paste0(tag, "_fits.tsv.gz")), sep = "\t")
  fwrite(rbindlist(lapply(rr, `[[`, "cells"), fill = TRUE), file.path(OUT, paste0(tag, "_cells.tsv.gz")), sep = "\t")
  fwrite(rbindlist(lapply(rr, `[[`, "reps"), fill = TRUE), file.path(OUT, paste0(tag, "_reps.tsv.gz")), sep = "\t")
}
dg <- c(if (file.exists(file.path(OUT, "armB_window_diagnostics.tsv"))) fread(file.path(OUT, "armB_window_diagnostics.tsv"))$files_used,
        if (file.exists(file.path(OUT, "pc1_window_diagnostics.tsv"))) fread(file.path(OUT, "pc1_window_diagnostics.tsv"))$files_used)
used <- unique(unlist(strsplit(dg, "|", fixed = TRUE)))
used <- c(paste0(used, ".bim"), paste0(used, ".ld"), unique(c(W[row %in% b_rows, rds], P1[row %in% p_rows, rds])))
log_msg("hashing ", length(used), " LD / rds files opened by Arm B and PC-1")
hh <- sha256_files(used, WORKERS)
fwrite(hh, file.path(OUT, "opened_files_sha256.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_R.txt"))
writeLines("DONE", file.path(OUT, "DONE"))
log_msg("Arm B + PC-1 done")
