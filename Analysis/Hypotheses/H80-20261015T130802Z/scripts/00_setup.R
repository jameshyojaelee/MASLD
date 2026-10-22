# 00_setup.R -- H80 job 0: input hashes (S1), frame recounts (G1/S3), window sampling on stream 0,
# per-window inputs, D2/D4/D5 diagnostics, PC-2 (eQTL fit reproduction), PC-3 (production COLOC
# reproduction and EUR LD panel check). Writes <out>/STATUS: PASS or STOP:<rule>:<reason>.
args0 <- commandArgs(trailingOnly = FALSE)
SCRIPTS <- dirname(normalizePath(sub("^--file=", "", args0[grep("^--file=", args0)])))
Sys.setenv(H80_SCRIPTS = SCRIPTS)
source(file.path(SCRIPTS, "h80_lib.R"))
A <- parse_args(list(out = NA, workers = "4", seed = "20260923", eur_gwas_panel = "polyfun",
                     pc3_parallel = "2", pc2_parallel = "2", debug_n = "0"))
OUT <- A$out
H80 <- dirname(SCRIPTS)
WORKERS <- as.integer(A$workers)
DEBUG_N <- as.integer(A$debug_n)   # >0 only in debug runs: limits PC-2/PC-3 work
stopifnot(as.integer(A$seed) == H80_SEED)
setDTthreads(WORKERS)
dir.create(file.path(OUT, "eqtl_ss"), showWarnings = FALSE)
dir.create(file.path(OUT, "gwas_ss"), showWarnings = FALSE)
diag_lines <- character(0)
D <- function(...) { s <- paste0(...); diag_lines <<- c(diag_lines, s); log_msg(s) }
status_rows <- list()
add_status <- function(check, value, expected, pass) {
  status_rows[[length(status_rows) + 1L]] <<- data.table(check = check, value = as.character(value),
                                                        expected = as.character(expected), pass = pass)
}
finish <- function(status) {
  fwrite(rbindlist(status_rows), file.path(OUT, "setup_status.tsv"), sep = "\t")
  writeLines(diag_lines, file.path(OUT, "diagnostics_setup.txt"))
  writeLines(status, file.path(OUT, "STATUS"))
  log_msg("SETUP STATUS: ", status)
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}

log_msg("H80 setup start; out = ", OUT)
verify_prespec(H80)

# ---- D1 / S1: input hashes ---------------------------------------------------------------------------
h <- verify_inputs(WORKERS)
fwrite(h, file.path(OUT, "input_hash_check.tsv"), sep = "\t")
for (i in seq_len(nrow(h))) D("D1 ", h$role[i], " ", h$kind[i], " observed ", h$observed[i],
                               if (!is.na(h$n_observed[i])) paste0(" n=", h$n_observed[i]) else "",
                               if (h$pass[i]) " OK" else " MISMATCH")
add_status("S1_input_hashes", sum(h$pass), nrow(h), all(h$pass))
if (!all(h$pass)) finish(paste0("STOP:S1:", paste(h[pass == FALSE, role], collapse = ",")))
afsrc <- file.path(FM, "config/gwas_af_sources.tsv")
D("D1 extra input (not in prespec table; read to choose the GWAS file 06 opens): ", afsrc, " sha256 ", sha256_file(afsrc))

# ---- G1 / S3: frame recounts ------------------------------------------------------------------------------
tier <- fread(file.path(FM, "config/gwas_trait_tier.tsv"))
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
t12 <- sort(tier[tier %in% c(1, 2), study_name], method = "radix")
m <- fread(file.path(R2DIR, "susie_coloc_all_gwas.csv"),
           select = c("gwas_name", "ensembl", "chr", "PP.H3.susie", "PP.H4.susie", "method", "eligibility_source"))
m <- m[gwas_name %in% t12]
regc <- reg[, .(gwas_name = study_name, reg_anc = ancestry, trait_type, N_tot, N_cases, sumstats_path)]
s <- merge(m[is.finite(PP.H4.susie)], regc, by = "gwas_name", all.x = TRUE)
win <- unique(s[, .(ensembl, reg_anc)])
spos <- unique(s[PP.H4.susie > 0.9, .(ensembl, reg_anc)])
g <- fread(file.path(R2DIR, "tier12/gene_level_coloc_tier12.csv"), select = c("ensembl", "susie_state_t12"))
sumf <- sort(Sys.glob(file.path(EQTL_FITS, "chr*_summary.tsv")), method = "radix")
e <- rbindlist(lapply(sumf, fread, colClasses = list(character = c("gene", "ensg"))), use.names = TRUE)
e_last <- e[!duplicated(ensg, fromLast = TRUE)]
e_last[, fit := fifelse(status != "ok", paste0("no_fit_", status),
                 fifelse(!(converged %in% TRUE), "not_converged", fifelse(n_cs > 0, "cs_ge1", "cs0")))]
gx <- merge(g, e_last[, .(ensembl = ensg, fit, n_snps, chr_fit = chr)], by = "ensembl", all.x = TRUE)
gx[is.na(fit), fit := "no_summary_row"]
nores_genes <- sort(gx[susie_state_t12 == "no_susie_result" & fit == "cs0", ensembl], method = "radix")
anc_counts <- win[, .N, by = reg_anc]
got <- c(t12 = length(t12), pairs_t12 = nrow(m), susie_pairs = nrow(s), pp4_gt05 = sum(s$PP.H4.susie > 0.5),
         pp4_gt09 = sum(s$PP.H4.susie > 0.9), pp3_gt08 = sum(s$PP.H3.susie > 0.8),
         seval_windows = nrow(win), seval_EUR = anc_counts[reg_anc == "EUR", N], seval_EAS = anc_counts[reg_anc == "EAS", N],
         seval_AFR = anc_counts[reg_anc == "AFR", N], seval_AMR = anc_counts[reg_anc == "AMR", N],
         seval_SAS = anc_counts[reg_anc == "SAS", N], seval_genes = uniqueN(s$ensembl),
         snores_genes = length(nores_genes), spos_windows = nrow(spos), cs_ge1_genes = sum(gx$fit == "cs_ge1"))
expd <- c(t12 = 35, pairs_t12 = 664022, susie_pairs = 12892, pp4_gt05 = 767, pp4_gt09 = 432, pp3_gt08 = 10980,
          seval_windows = 5683, seval_EUR = 3838, seval_EAS = 897, seval_AFR = 447, seval_AMR = 413, seval_SAS = 88,
          seval_genes = 3918, snores_genes = 13704, spos_windows = 286, cs_ge1_genes = 4979)
for (k in names(expd)) {
  D("D2/G1 ", k, " = ", got[[k]], " (expected ", expd[[k]], ")", if (got[[k]] == expd[[k]]) " OK" else " DIFFERS")
  add_status(paste0("G1_", k), got[[k]], expd[[k]], got[[k]] == expd[[k]])
}
if (any(got[names(expd)] != expd)) finish(paste0("STOP:S3:", paste(names(expd)[got[names(expd)] != expd], collapse = ",")))
fwrite(gx, file.path(OUT, "gene_state_by_fit.tsv"), sep = "\t")

# ---- streams --------------------------------------------------------------------------------------------
streams <- make_streams(3000L)
saveRDS(streams, file.path(OUT, "streams.rds"))
D("Seed: RNGkind L'Ecuyer-CMRG, set.seed(", H80_SEED, "); streams 0..3000 by nextRNGStream; stream0[2:4] = ",
  paste(streams[["0"]][2:4], collapse = ","))

# ---- GWAS access (the file 06 would open: 06:165-177) -------------------------------------------------------
af_tab <- fread(afsrc)
gwas_path_of <- function(study) {
  hit <- af_tab[study_name == study & file.exists(af_sumstats_path)]
  if (nrow(hit) == 1L) return(hit$af_sumstats_path[[1L]])
  file.path(FM, reg[study_name == study, sumstats_path])
}
GW <- new.env()
opened_gwas <- character(0)
gwas_rows <- function(study, chr) {
  if (is.null(GW[[study]])) {
    p <- gwas_path_of(study)
    log_msg("reading GWAS ", study, " from ", p)
    x <- fread(p, select = c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval"))
    x[, chromosome := as.character(chromosome)]
    setkey(x, chromosome)
    GW[[study]] <- x
    opened_gwas <<- c(opened_gwas, p)
  }
  GW[[study]][.(as.character(chr))]
}
# GWAS rows at the positions of E, one row per position: the row with the smallest p
# (06:400-412 orders the merge by gwas_pval and keeps the first row per position).
gwas_at <- function(study, chr, pos) {
  x <- gwas_rows(study, chr)
  x <- x[position %in% pos]
  x <- x[order(pval)][!duplicated(position)]
  x[order(position)]
}
gwas_ss_for_ld <- function(x, chr) {
  data.frame(chromosome = as.integer(chr), position = x$position, allele1 = x$allele1, allele2 = x$allele2,
             beta = x$beta, standard_error = x$se, stringsAsFactors = FALSE)
}
eur_gwas_dir <- if (A$eur_gwas_panel == "ukbb") PANEL_DIRS[["EUR_UKBB"]] else PANEL_DIRS[["EUR"]]
gwas_panel_dir <- function(anc) if (anc == "EUR") eur_gwas_dir else PANEL_DIRS[[anc]]
D("EUR GWAS-side LD panel for this run: ", eur_gwas_dir)

rds_path_of <- function(ensg, chr) file.path(EQTL_FITS, paste0("chr", chr), paste0(ensg, "_susie.rds"))

# ---- frames per stratum ----------------------------------------------------------------------------------
gene_chr <- unique(rbind(s[, .(ensembl, chr)], e_last[, .(ensembl = ensg, chr)]))
gene_chr <- gene_chr[!duplicated(ensembl)]
seval_studies <- s[, .(studies = list(sort(unique(gwas_name), method = "radix")), chr = chr[1]), by = .(ensembl, reg_anc)]
mn <- merge(m[ensembl %in% nores_genes], regc[, .(gwas_name, reg_anc)], by = "gwas_name")
nores_studies <- mn[reg_anc == "EUR", .(studies = list(sort(unique(gwas_name), method = "radix")), chr = chr[1]),
                    by = .(ensembl, reg_anc)]
D("S-nores genes with an EUR tier-1/2 master row: ", nrow(nores_studies), " of ", length(nores_genes))
frames <- list(
  EUR_Seval = seval_studies[reg_anc == "EUR"], EUR_Snores = nores_studies,
  EAS = seval_studies[reg_anc == "EAS"], AFR = seval_studies[reg_anc == "AFR"],
  AMR = seval_studies[reg_anc == "AMR"], SAS = seval_studies[reg_anc == "SAS"])
for (k in names(frames)) setorder(frames[[k]], ensembl)
n_plan <- c(EUR_Seval = 60L, EUR_Snores = 60L, EAS = 30L, AFR = 20L, AMR = 15L, SAS = 15L)

check_window <- function(ensg, chr, anc, study) {
  rp <- rds_path_of(ensg, chr)
  if (!file.exists(rp)) return(list(ok = FALSE, reason = "no_rds"))
  so <- readRDS(rp)
  if (!isTRUE(so$converged)) return(list(ok = FALSE, reason = "rds_not_converged"))
  ids <- colnames(so$lbf_variable)
  if (length(ids) < MIN_TRIPLE_SNPS) return(list(ok = FALSE, reason = "eqtl_fit_lt50"))
  pos <- as.integer(sub("^[0-9]+:", "", ids))
  gx_ <- gwas_at(study, chr, pos)
  if (nrow(gx_) < MIN_TRIPLE_SNPS) return(list(ok = FALSE, reason = "gwas_rows_lt50"))
  pl <- ld_plan(gwas_ss_for_ld(gx_, chr), as.integer(chr), min(gx_$position), max(gx_$position), gwas_panel_dir(anc))
  if (is.null(pl)) return(list(ok = FALSE, reason = "gwas_ld_plan_null"))
  gids <- snp_ids_of(pl$ss)
  nEG <- sum(gids %in% ids)
  if (nEG < MIN_TRIPLE_SNPS) return(list(ok = FALSE, reason = "EG_lt50"))
  list(ok = TRUE, reason = "eligible", nE = length(ids), nG = length(gids), nEG = nEG,
       converged = TRUE, n_cs_rds = length(so$sets$cs), rds = rp)
}

# ---- sampling walk on stream 0 ----------------------------------------------------------------------------
set_stream(streams, 0L)
win_rows <- list(); walk_log <- list()
for (k in names(frames)) {
  fr <- frames[[k]]
  perm <- sample.int(nrow(fr))
  kept <- 0L
  for (j in perm) {
    if (kept >= n_plan[[k]]) break
    st <- fr$studies[[j]]
    study <- st[sample.int(length(st), 1L)]
    anc <- fr$reg_anc[j]; ensg <- fr$ensembl[j]; chr <- fr$chr[j]
    r <- check_window(ensg, chr, anc, study)
    walk_log[[length(walk_log) + 1L]] <- data.table(stratum = k, ensembl = ensg, chr = chr, ancestry = anc,
                                                    study = study, n_candidate_studies = length(st), reason = r$reason)
    if (r$ok) {
      kept <- kept + 1L
      rr <- reg[study_name == study]
      win_rows[[length(win_rows) + 1L]] <- data.table(
        stratum = k, ancestry = anc, ensembl = ensg, chr = as.integer(chr), study = study,
        trait_type = rr$trait_type, N_tot = rr$N_tot, N_cases = rr$N_cases,
        N_eff = gwas_n_eff(rr$trait_type, rr$N_tot, rr$N_cases),
        nE = r$nE, nG_plan = r$nG, nEG_plan = r$nEG, n_cs_rds = r$n_cs_rds, rds = r$rds,
        gwas_path = gwas_path_of(study), gwas_panel = gwas_panel_dir(anc))
    }
  }
  D("D2 stratum ", k, ": walked ", sum(sapply(walk_log, function(x) x$stratum == k)), ", kept ", kept, " of ", n_plan[[k]],
    if (kept < n_plan[[k]]) " (FRAME EXHAUSTED: flagged)" else "", if (kept < 10) " (UNTESTABLE, S4)" else "")
  add_status(paste0("D2_kept_", k), kept, n_plan[[k]], kept >= 10)
}
W <- rbindlist(win_rows)
W[, row := .I]
W[, stream := row]
W[, armB := ancestry == "EUR"]
W[armB == TRUE, armB_index := seq_len(.N)]
W[, armB_stream := fifelse(armB, 1000L + armB_index, NA_integer_)]
WL <- rbindlist(walk_log)
fwrite(WL, file.path(OUT, "sampling_walk.tsv"), sep = "\t")
D("D2 ineligible windows skipped by reason: ", if (any(WL$reason != "eligible")) paste(WL[reason != "eligible", .N, by = reason][, paste0(reason, "=", N)], collapse = "; ") else "none")

# PC-3 pairs on stream 0 (after the windows): 5 of the EUR rerun rows with PP.H4 > 0.5, then 5 of those <= 0.5.
rr_eur <- s[reg_anc == "EUR" & eligibility_source == "rerun_shared_posterior_check"]
setorder(rr_eur, gwas_name, ensembl)
hi <- rr_eur[PP.H4.susie > 0.5]; lo <- rr_eur[PP.H4.susie <= 0.5]
D("PC-3 pool: EUR rerun rows ", nrow(rr_eur), " (PP.H4>0.5: ", nrow(hi), "; <=0.5: ", nrow(lo), ")")
add_status("PC3_pool_hi", nrow(hi), 642, nrow(hi) == 642)
add_status("PC3_pool_lo", nrow(lo), 24, nrow(lo) == 24)
pc3 <- rbind(hi[sample.int(nrow(hi), 5L)], lo[sample.int(nrow(lo), 5L)])
pc3 <- pc3[, .(gwas_name, ensembl, chr, PP.H4.susie, PP.H3.susie)]
fwrite(pc3, file.path(OUT, "pc3_pairs.tsv"), sep = "\t")
stream0_end <- get(".Random.seed", envir = .GlobalEnv)
D("stream 0 state after sampling and PC-3 draws: ", digest(stream0_end, algo = "sha256"))

fwrite(W, file.path(OUT, "windows.tsv"), sep = "\t")
D("windows.tsv: ", nrow(W), " windows; by stratum: ", paste(W[, .N, by = stratum][, paste0(stratum, "=", N)], collapse = " "))

# D2: |E| of sampled windows against the frame's.
fq <- function(x) paste(round(quantile(x, c(0, .25, .5, .75, 1), na.rm = TRUE)), collapse = "/")
fe <- merge(rbindlist(lapply(names(frames), function(k) frames[[k]][, .(stratum = k, ensembl)])),
            e_last[, .(ensembl = ensg, n_snps)], by = "ensembl", all.x = TRUE)
for (k in names(frames)) D("D2 |E| min/q1/median/q3/max sampled ", fq(W[stratum == k, nE]), " frame ", fq(fe[stratum == k, n_snps]))

# ---- PC-1 windows -------------------------------------------------------------------------------------------
sp <- s[PP.H4.susie > 0.9]
setorder(sp, ensembl, reg_anc, -PP.H4.susie, gwas_name)
pc1 <- sp[!duplicated(sp[, .(ensembl, reg_anc)])]
setorder(pc1, ensembl, reg_anc)
pc1 <- pc1[, .(ensembl, ancestry = reg_anc, chr = as.integer(chr), study = gwas_name, PP.H4.susie, trait_type, N_tot, N_cases)]
pc1[, N_eff := mapply(gwas_n_eff, trait_type, N_tot, N_cases)]
pc1[, row := .I][, stream := 2000L + row]
pc1[, rds := rds_path_of(ensembl, chr)]
pc1[, rds_converged := vapply(rds, function(p) file.exists(p) && isTRUE(readRDS(p)$converged), logical(1))]
pc1[, gwas_path := vapply(study, gwas_path_of, character(1))]
pc1[, gwas_panel := vapply(ancestry, gwas_panel_dir, character(1))]
fwrite(pc1, file.path(OUT, "pc1_windows.tsv"), sep = "\t")
D("PC-1 windows: ", nrow(pc1), " (", paste(pc1[, .N, by = ancestry][, paste0(ancestry, "=", N)], collapse = " "),
  "); rds converged ", sum(pc1$rds_converged))

# ---- Broadaway rows per gene (08 filter) ----------------------------------------------------------------------
need <- unique(rbind(W[, .(ensembl, chr)], pc1[, .(ensembl, chr)]))
opened_eqtl <- character(0)
for (cc in sort(unique(need$chr))) {
  f <- file.path(EQTL_DIR, paste0("chr", cc, "_marginal_summary_results.tsv"))
  dt <- fread(f, select = c("CHR", "POS", "NEA", "EA", "Beta", "SE", "ENSG", "GeneSymbol"))
  opened_eqtl <- c(opened_eqtl, f)
  dt <- dt[!is.na(POS) & !is.na(Beta) & !is.na(SE) & SE > 0 & ENSG != ""]
  for (gid in need[chr == cc, ensembl]) saveRDS(dt[ENSG == gid], file.path(OUT, "eqtl_ss", paste0(gid, ".rds")))
  rm(dt); gc(FALSE)
  log_msg("Broadaway chr", cc, ": saved ", sum(need$chr == cc), " genes")
}

# ---- GWAS rows per window (positions of E) --------------------------------------------------------------------
save_gwas_ss <- function(key, study, chr, rds) {
  ids <- colnames(readRDS(rds)$lbf_variable)
  pos <- as.integer(sub("^[0-9]+:", "", ids))
  x <- gwas_at(study, chr, pos)
  saveRDS(gwas_ss_for_ld(x, chr), file.path(OUT, "gwas_ss", paste0(key, ".rds")))
  nrow(x)
}
W[, n_gwas_rows := mapply(function(r, st, ch, rd) save_gwas_ss(paste0("A", r), st, ch, rd), row, study, chr, rds)]
pc1[, n_gwas_rows := mapply(function(r, st, ch, rd) if (file.exists(rd)) save_gwas_ss(paste0("P", r), st, ch, rd) else NA_integer_,
                            row, study, chr, rds)]
fwrite(W, file.path(OUT, "windows.tsv"), sep = "\t")
fwrite(pc1, file.path(OUT, "pc1_windows.tsv"), sep = "\t")
for (st in ls(GW)) rm(list = st, envir = GW)
gc(FALSE)

# ---- D4: convergence read from rds ---------------------------------------------------------------------------
D("D4 sampled windows: all ", nrow(W), " rds report converged = TRUE (eligibility); dropped for non-convergence during the walk: ",
  WL[reason == "rds_not_converged", .N])
nc <- gx[fit == "not_converged" & susie_state_t12 %in% c("positive", "tested_le_0.5")]
nc[, rds_converged := vapply(seq_len(.N), function(i) {
  p <- rds_path_of(ensembl[i], chr_fit[i]); file.exists(p) && isTRUE(readRDS(p)$converged) }, logical(1))]
fwrite(nc, file.path(OUT, "D4_notconverged_summary_genes.tsv"), sep = "\t")
D("D4 genes whose last summary row says not converged: positive ", sum(nc$susie_state_t12 == "positive"),
  ", tested ", sum(nc$susie_state_t12 == "tested_le_0.5"), "; rds converged = TRUE in ", sum(nc$rds_converged), " of ", nrow(nc))

# ---- PC-2: eQTL fit reproduction through the 08 code path, plus loader identity (D5) -----------------------------
Sys.unsetenv(c("LD_PANEL", "UKBB_LD_DIR"))
pc2w <- rbind(head(W[stratum == "EUR_Seval"], 10), head(W[stratum == "EUR_Snores"], 10))
if (DEBUG_N > 0) pc2w <- head(pc2w, DEBUG_N)
pc2_one <- function(i) {
  prod <- load_production_loader()
  set_stream(streams, 0L)   # wrapped fits restore whatever state is current; no stream draws here
  w <- pc2w[i]
  gene_dt <- readRDS(file.path(OUT, "eqtl_ss", paste0(w$ensembl, ".rds")))
  ss <- eqtl_ss_from_rows(gene_dt)
  se <- eqtl_window(gene_dt)
  t0 <- proc.time()[["elapsed"]]
  ld <- prod$get_ld_per_locus(ss, LOCUS = w$ensembl, CHR = w$chr, START = se[1], END = se[2], ancestry = "EUR")
  t_prod <- proc.time()[["elapsed"]] - t0
  ssf <- ld[[1]]; R_mat <- ld[[2]]; if (!is.matrix(R_mat)) R_mat <- as.matrix(R_mat)
  z <- ssf$beta / ssf$standard_error
  ids <- paste0(w$chr, ":", ssf$position); names(z) <- ids
  R_reg <- R_mat; diag(R_reg) <- diag(R_reg) + LD_REGULARIZE; colnames(R_reg) <- rownames(R_reg) <- ids
  fr <- fit_eqtl(z, R_reg)
  so <- readRDS(w$rds)
  fit <- fr$res
  ids_ok <- identical(colnames(so$lbf_variable), ids)
  pip_d <- if (ids_ok) max(abs(fit$pip - so$pip)) else NA_real_
  alpha_d <- if (ids_ok && all(dim(fit$alpha) == dim(so$alpha))) max(abs(fit$alpha - so$alpha)) else NA_real_
  Xc <- cov2cor(R_reg)
  cs_st <- susie_get_cs(so, Xcorr = Xc, coverage = 0.95, min_abs_corr = 0.5, n_purity = Inf)
  cs_rr <- susie_get_cs(fit, Xcorr = Xc, coverage = 0.95, min_abs_corr = 0.5, n_purity = Inf)
  n_st <- length(cs_st$cs); n_rr <- length(cs_rr$cs)
  # loader identity: column-select copy against the production cached loader (D5)
  t1 <- proc.time()[["elapsed"]]
  fl <- ld_load_fast(ss, w$chr, se[1], se[2], PANEL_DIRS[["EUR"]])
  t_fast <- proc.time()[["elapsed"]] - t1
  same_R <- identical(unname(fl$R), unname(R_mat)) && identical(dimnames(fl$R), dimnames(R_mat))
  same_ss <- isTRUE(all.equal(as.data.frame(fl$ss), as.data.frame(ssf), check.attributes = FALSE, tolerance = 0))
  data.table(stratum = w$stratum, ensembl = w$ensembl, chr = w$chr, p = length(ids), ids_identical = ids_ok,
             max_pip_diff = pip_d, max_alpha_diff = alpha_d, converged = isTRUE(fit$converged),
             cs_det_stored = n_st, cs_det_rerun = n_rr, cs_raw_stored = length(so$sets$cs), cs_raw_rerun = length(fit$sets$cs),
             cs_raw_sizes_stored = paste(lengths(so$sets$cs), collapse = ";"), cs_raw_sizes_rerun = paste(lengths(fit$sets$cs), collapse = ";"),
             rerun_drew_rng = fr$drew, fast_loader_R_identical = same_R, fast_loader_ss_identical = same_ss,
             fast_loader_n_na = fl$n_na, t_prod_loader_s = t_prod, t_fast_loader_s = t_fast, t_fit_s = fr$time)
}
log_msg("PC-2: ", nrow(pc2w), " windows")
pc2 <- rbindlist(mclapply(seq_len(nrow(pc2w)), pc2_one, mc.cores = as.integer(A$pc2_parallel), mc.preschedule = FALSE))
fwrite(pc2, file.path(OUT, "pc2_results.tsv"), sep = "\t")
pc2[, pass := ids_identical & max_pip_diff < 1e-4 & max_alpha_diff < 1e-4 & cs_det_stored == cs_det_rerun]
for (i in seq_len(nrow(pc2))) D("PC-2 ", pc2$ensembl[i], " p=", pc2$p[i], " ids ", pc2$ids_identical[i], " maxPIPdiff ",
                                 signif(pc2$max_pip_diff[i], 3), " maxAlphaDiff ", signif(pc2$max_alpha_diff[i], 3),
                                 " CS(det) stored/rerun ", pc2$cs_det_stored[i], "/", pc2$cs_det_rerun[i],
                                 " CS(raw) ", pc2$cs_raw_stored[i], "/", pc2$cs_raw_rerun[i],
                                 " fastLoader R/ss identical ", pc2$fast_loader_R_identical[i], "/", pc2$fast_loader_ss_identical[i],
                                 if (pc2$pass[i]) " PASS" else " FAIL")
pc2_pass <- nrow(pc2) == nrow(pc2w) && all(pc2$pass)
add_status("PC2_eqtl_fit_reproduction", sum(pc2$pass), nrow(pc2w), pc2_pass)
loader_ok <- all(pc2$fast_loader_R_identical & pc2$fast_loader_ss_identical)
add_status("D5_fast_loader_identity_eqtl", sum(pc2$fast_loader_R_identical & pc2$fast_loader_ss_identical), nrow(pc2), loader_ok)

# Loader identity on the other panels: GWAS side of the first window of every stratum (EUR polyfun
# and the four Gram panels) and 1kg_eur on the first EUR window (the Arm B generator).
Sys.setenv(AFR_LD_DIR = PANEL_DIRS[["AFR"]], AMR_LD_DIR = PANEL_DIRS[["AMR"]],
           EAS_LD_DIR = PANEL_DIRS[["EAS"]], SAS_LD_DIR = PANEL_DIRS[["SAS"]])
id_cases <- W[!duplicated(stratum)]
if (DEBUG_N > 0) id_cases <- id_cases[stratum %in% c("EUR_Seval", "AFR")]
id_one <- function(i) {
  prod <- load_production_loader()
  w <- id_cases[i]
  gss <- readRDS(file.path(OUT, "gwas_ss", paste0("A", w$row, ".rds")))
  st <- min(gss$position); en <- max(gss$position)
  pr <- prod$get_ld_per_locus(gss, LOCUS = w$ensembl, CHR = w$chr, START = st, END = en, ancestry = w$ancestry)
  fl <- ld_load_fast(gss, w$chr, st, en, gwas_panel_dir(w$ancestry))
  data.table(case = paste0("gwas_side_", w$stratum), ensembl = w$ensembl, panel = basename(gwas_panel_dir(w$ancestry)),
             p = nrow(pr[[1]]), R_identical = identical(unname(fl$R), unname(as.matrix(pr[[2]]))) &&
               identical(dimnames(fl$R), dimnames(pr[[2]])),
             ss_identical = isTRUE(all.equal(as.data.frame(fl$ss), as.data.frame(pr[[1]]), check.attributes = FALSE, tolerance = 0)),
             n_na = fl$n_na)
}
idt <- rbindlist(mclapply(seq_len(nrow(id_cases)), id_one, mc.cores = as.integer(A$pc2_parallel), mc.preschedule = FALSE))
# 1kg_eur (Arm B generator) on the first EUR window, production loader with LD_PANEL=1kg.
w1 <- W[ancestry == "EUR"][1]
gene_dt <- readRDS(file.path(OUT, "eqtl_ss", paste0(w1$ensembl, ".rds")))
ss1 <- eqtl_ss_from_rows(gene_dt); se1 <- eqtl_window(gene_dt)
Sys.setenv(LD_PANEL = "1kg")
prod <- load_production_loader()
pr <- prod$get_ld_per_locus(ss1, LOCUS = w1$ensembl, CHR = w1$chr, START = se1[1], END = se1[2], ancestry = "EUR")
Sys.unsetenv("LD_PANEL")
fl <- ld_load_fast(ss1, w1$chr, se1[1], se1[2], PANEL_DIRS[["EUR_1KG"]])
idt <- rbind(idt, data.table(case = "eqtl_side_1kg_eur", ensembl = w1$ensembl, panel = "1kg_eur", p = nrow(pr[[1]]),
                             R_identical = identical(unname(fl$R), unname(as.matrix(pr[[2]]))) && identical(dimnames(fl$R), dimnames(pr[[2]])),
                             ss_identical = isTRUE(all.equal(as.data.frame(fl$ss), as.data.frame(pr[[1]]), check.attributes = FALSE, tolerance = 0)),
                             n_na = fl$n_na))
rm(prod, pr, fl); gc(FALSE)
fwrite(idt, file.path(OUT, "D5_loader_identity.tsv"), sep = "\t")
for (i in seq_len(nrow(idt))) D("D5 loader identity ", idt$case[i], " ", idt$panel[i], " p=", idt$p[i], " R ", idt$R_identical[i],
                                 " ss ", idt$ss_identical[i], " NA ", idt$n_na[i])
loader_ok <- loader_ok && all(idt$R_identical & idt$ss_identical)
add_status("D5_fast_loader_identity_all", sum(idt$R_identical & idt$ss_identical), nrow(idt), all(idt$R_identical & idt$ss_identical))

# ---- PC-3: unmodified production runner on 10 EUR rerun pairs ----------------------------------------------------
run_pc3 <- function(panel_label) {
  root <- file.path(OUT, paste0("pc3_", panel_label))
  dir.create(file.path(root, "coloc"), recursive = TRUE)
  flt <- file.path(root, "gene_filter.tsv")
  pp <- if (DEBUG_N > 0) head(pc3, DEBUG_N) else pc3
  fwrite(pp[, .(gwas_name, chr, ensembl)], flt, sep = "\t")
  tasks <- unique(pp[, .(gwas_name, chr)])
  one <- function(i) {
    envs <- c(paste0("MASLD_PROJECT_ROOT=", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
              paste0("AFR_LD_DIR=", PANEL_DIRS[["AFR"]]), paste0("AMR_LD_DIR=", PANEL_DIRS[["AMR"]]),
              paste0("EAS_LD_DIR=", PANEL_DIRS[["EAS"]]), paste0("SAS_LD_DIR=", PANEL_DIRS[["SAS"]]),
              paste0("EQTL_SUSIE_DIR=", EQTL_FITS), paste0("COLOC_OUT_ROOT=", file.path(root, "coloc")),
              paste0("COLOC_GENE_FILTER=", flt), paste0("GWAS_NAME=", tasks$gwas_name[i]),
              paste0("CHR_FILTER=", tasks$chr[i]), "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1")
    if (panel_label == "ukbb") envs <- c(envs, "LD_PANEL=ukbb")
    lg <- file.path(root, sprintf("06_%s_chr%s.log", tasks$gwas_name[i], tasks$chr[i]))
    cmd <- sprintf("cd %s && env -u LD_PANEL -u UKBB_LD_DIR %s %s %s > %s 2>&1",
                   shQuote(file.path(FM, "src")), paste(shQuote(envs), collapse = " "),
                   shQuote(file.path(R.home("bin"), "Rscript")), shQuote(file.path(FM, "src/06_susie_coloc.R")), shQuote(lg))
    rc <- system(cmd)
    data.table(gwas_name = tasks$gwas_name[i], chr = tasks$chr[i], rc = rc)
  }
  rcs <- rbindlist(mclapply(seq_len(nrow(tasks)), one, mc.cores = as.integer(A$pc3_parallel), mc.preschedule = FALSE))
  res <- rbindlist(lapply(seq_len(nrow(pp)), function(i) {
    f <- file.path(root, "coloc", pp$gwas_name[i], paste0("susie_coloc_chr", pp$chr[i], ".csv"))
    x <- if (file.exists(f)) fread(f)[ensembl == pp$ensembl[i]] else data.table()
    data.table(gwas_name = pp$gwas_name[i], ensembl = pp$ensembl[i], chr = pp$chr[i],
               master_PP4 = pp$PP.H4.susie[i], master_PP3 = pp$PP.H3.susie[i],
               rerun_PP4 = if (nrow(x) == 1) x$PP.H4.susie else NA_real_,
               rerun_PP3 = if (nrow(x) == 1) x$PP.H3.susie else NA_real_,
               rerun_method = if (nrow(x) == 1) x$method else NA_character_)
  }))
  res[, `:=`(d4 = abs(rerun_PP4 - master_PP4), d3 = abs(rerun_PP3 - master_PP3))]
  res[, pass := is.finite(d4) & is.finite(d3) & d4 < 1e-6 & d3 < 1e-6]
  res[, panel := panel_label]
  list(res = res, rc = rcs)
}
p3 <- run_pc3("ldpanel_unset")
pc3_res <- p3$res
pc3_panel <- "polyfun"
pc3_pass_unset <- all(pc3_res$pass) && nrow(pc3_res) == (if (DEBUG_N > 0) min(DEBUG_N, 10) else 10)
for (i in seq_len(nrow(pc3_res))) D("PC-3 (LD_PANEL unset) ", pc3_res$gwas_name[i], " ", pc3_res$ensembl[i], " master PP4/PP3 ",
                                     signif(pc3_res$master_PP4[i], 8), "/", signif(pc3_res$master_PP3[i], 8), " rerun ",
                                     signif(pc3_res$rerun_PP4[i], 8), "/", signif(pc3_res$rerun_PP3[i], 8),
                                     if (pc3_res$pass[i]) " PASS" else " FAIL")
D("PC-3 06 return codes: ", paste(p3$rc$rc, collapse = ","))
if (!pc3_pass_unset) {
  p3u <- run_pc3("ukbb")
  pc3_res <- rbind(pc3_res, p3u$res)
  for (i in seq_len(nrow(p3u$res))) D("PC-3 (LD_PANEL=ukbb) ", p3u$res$gwas_name[i], " ", p3u$res$ensembl[i],
                                       if (p3u$res$pass[i]) " PASS" else " FAIL")
  pc3_panel <- if (all(p3u$res$pass)) "ukbb" else "none"
}
fwrite(pc3_res, file.path(OUT, "pc3_results.tsv"), sep = "\t")
add_status("PC3_production_reproduction_polyfun", sum(pc3_res[panel == "ldpanel_unset", pass]), nrow(pc3_res[panel == "ldpanel_unset"]), pc3_pass_unset)
writeLines(pc3_panel, file.path(OUT, "production_eur_panel.txt"))

# ---- hashes of every GWAS / Broadaway file opened (setup and PC-3) ----------------------------------------------
opened <- unique(c(opened_gwas, vapply(unique(pc3$gwas_name), gwas_path_of, character(1)),
                   file.path(EQTL_DIR, sprintf("chr%d_marginal_summary_results.tsv", 1:22))))
hh <- sha256_files(opened, WORKERS)
fwrite(hh, file.path(OUT, "opened_files_sha256.tsv"), sep = "\t")
D("D1 hashed ", nrow(hh), " opened GWAS/Broadaway files (", round(sum(hh$bytes) / 1e9, 1), " GB)")

# ---- status -----------------------------------------------------------------------------------------------
if (!pc2_pass) finish("STOP:S2:PC-2_eqtl_fit_reproduction_failed")
if (!loader_ok) finish("STOP:S2:fast_loader_not_identical_to_production_loader")
if (pc3_panel == "none") finish("STOP:S2:PC-3_reproduced_with_neither_panel")
if (pc3_panel == "ukbb" && A$eur_gwas_panel != "ukbb") finish("RESTART:EUR_GWAS_PANEL=ukbb")
if (DEBUG_N > 0) finish("DEBUG_ONLY")
finish("PASS")
