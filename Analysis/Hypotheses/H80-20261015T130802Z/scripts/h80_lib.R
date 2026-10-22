# h80_lib.R -- shared code for H80 (prespecification PRESPEC_H80.md in the run directory).
# Calibration of the production SuSiE-COLOC pipeline at real tier-1/2 LD windows.
# Every production setting below is copied from GWAS/finemapping/src/06_susie_coloc.R and
# 08_run_eqtl_susie.R (hashes verified at run start). Production files are read, never edited.

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(coloc)
  library(parallel)
  library(digest)
})

H80_PRESPEC_SHA <- "29e4a3a6f17d61644b332f8033318c6bc49f711eeb519483122822f3b1b9fe46"
H80_SEED <- 20260923L
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FM <- file.path(ROOT, "GWAS/finemapping")
LDREF <- file.path(FM, "data/ld_ref")
EQTL_DIR <- file.path(ROOT, "data/broadaway_eqtl")
EQTL_FITS <- file.path(FM, "results/eqtl_susie_polyfun")
R2DIR <- file.path(FM, "results/susie_coloc_r2_20260924T053244Z")
ROUND4 <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z")
PLAN_PRESPEC <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H80.md")

# ---- production settings (06:301-312, 08:26-30) --------------------------------------------
EQTL_N <- 1183L
SUSIE_L <- 10L
LD_REGULARIZE <- 1e-3
MIN_TRIPLE_SNPS <- 50L
CIS_WINDOW_BP <- 500000L
PURITY_FLOOR <- 0.5

# ---- simulation design (prespec "Design") -----------------------------------------------------
N_REP <- 10L
A_EQTL_LAMBDAS <- c(0, 3, 4, 5, 6, 7, 8, 10, 15, 20)      # nine E-sweep values plus 20
A_ESWEEP <- c(0, 3, 4, 5, 6, 7, 8, 10, 15)
A_GSWEEP <- c(3, 4, 5, 6, 8, 15)                           # 10 is the E-sweep (10, 10) cell
A_GWAS_SHARED <- c(3, 4, 5, 6, 8, 10, 15)
A_GWAS_DISTINCT <- c(10, 20)
A_DISTINCT_CELLS <- list(c(6, 10), c(10, 10), c(20, 20))
B_EQTL_LAMBDAS <- c(0, 4, 6, 8, 10, 15)
B_ESWEEP <- c(0, 4, 6, 8, 15)
DISTINCT_R2_MAX <- 0.1
DISTINCT_POOL_MIN <- 20L

PANEL_DIRS <- c(
  EUR = file.path(LDREF, "polyfun_eur"),
  AFR = file.path(LDREF, "1kg_afr_gram"),
  AMR = file.path(LDREF, "1kg_amr_gram"),
  EAS = file.path(LDREF, "1kg_eas_gram"),
  SAS = file.path(LDREF, "1kg_sas_gram"),
  EUR_1KG = file.path(LDREF, "1kg_eur"),
  EUR_UKBB = file.path(LDREF, "ukbb_eur")
)

log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%dT%H:%M:%S] "), ..., "\n", sep = "")
  flush(stdout())
}

parse_args <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1L
  while (i <= length(a)) {
    k <- sub("^--", "", a[i])
    if (i + 1L <= length(a) && !startsWith(a[i + 1L], "--")) {
      out[[k]] <- a[i + 1L]; i <- i + 2L
    } else {
      out[[k]] <- "TRUE"; i <- i + 1L
    }
  }
  out
}

# ---- hashing ------------------------------------------------------------------------------------
sha256_file <- function(p) {
  x <- system2("sha256sum", shQuote(p), stdout = TRUE)
  sub(" .*$", "", x[1])
}
sha256_files <- function(paths, cores = 1L) {
  paths <- unique(paths)
  h <- unlist(mclapply(paths, sha256_file, mc.cores = cores, mc.preschedule = FALSE))
  data.table(path = paths, sha256 = h, bytes = file.size(paths))
}
bash_hash <- function(cmd) {
  x <- system2("bash", c("-c", shQuote(paste0("export LC_ALL=C; ", cmd))), stdout = TRUE)
  sub(" .*$", "", x[length(x)])
}

# Single-file input hashes (prespec "Inputs with sha256"; paths relative to ROOT).
INPUT_FILES <- data.table::fread(text = "
role\tpath\tsha256
coloc_r2_master\tGWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv\tacdfa306f1ee45383758bc2522602fc2828c4060080b7ad2d96134fc4a399887
coloc_r2_gene_level_t12\tGWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/tier12/gene_level_coloc_tier12.csv\t6360984750864e7ab56593bce90e8ce665d781ec72a1ba52790d6ba53e83b8d6
tier_table\tGWAS/finemapping/config/gwas_trait_tier.tsv\t865f108de0782a708469b7a4e42943075de79b5ae062dc452f4050034bc15fc7
gwas_registry\tGWAS/finemapping/config/gwas_registry.tsv\t6823e42f5c1373785e2a23e08a4e8657fc3c2363d7d4259f6978e8ecfe8b9b8c
coloc_runner_06\tGWAS/finemapping/src/06_susie_coloc.R\te22b3ffbef48de9b96e220399275edc0113ec751fc765238ccefd076697114cc
eligibility_helper\tGWAS/finemapping/src/coloc_signal_eligibility.R\tfa7a42c621a55e3dee80bdb8bddad5d5049f697e98f6eb1d4e310aef4fb30675
eqtl_fitter_08\tGWAS/finemapping/src/08_run_eqtl_susie.R\tc6f24b8e7716f54e4bbf11d87bb0eb0ba1174d7f2e4fb78742ef4d10df0bf1cd
ld_loader\tGWAS/finemapping/src/finemapping_functions.R\t0654039ffd5cd4017d24b77e8d967ab08a836f3001dbb00586f59054fd997b2f
ld_loader_cached\tGWAS/finemapping/src/perf/get_ld_per_locus_cached.R\t5cbb6d530b2f3d3dbb858a0ca3072311d546caa62c405f5966d910f382175d70
r2_run_env_09b\tGWAS/finemapping/src/perf/09b_coloc_eligibility_batch.sbatch\t70ad56bed1c45aaa0881ccd487a86ad6ee10383ee8f3f05068ad3eda16270ea6
rerun_env_06c\tGWAS/finemapping/src/seqfunc_v2/finemap/06c_coloc_rerun_full.sbatch\t11ac3ca9f1f889218abc774a24f4a1fd6c909c38e94664db2e794bb53b0ebfdc
gate_pilot_script\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG19/pilot_timing.R\t0aa125f736d01026e40c66fdc734478c420263d60cebb5ed603abe7254e8263c
gate_pilot_output\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG19/pilot_runs.tsv\tda2df3ddb6eb3b46f9ce8564050e3e333e21c43e071c023578d1d7b7b8e38fd8
count_script\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80/count_n.py\t0b23c3a2847eeded5bfec17ea1809e6702d6496785656118f5e2b78689f6d2aa
count_output\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80/count_n.out\tfbe0d90abe674f32c7266f18914c431d6f32dc7a90445bd7c90e53be8001a907
frame_script\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80/count_frames.py\t49d9e5aa76dc5ebf09feedc532b1626babdbd3b51b93985a440cb53f17dd3494
frame_output\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80/count_frames.out\t3ba6dd88c7610a2b049bc3a2fe742a0bb47d222a39919182b30317c237d61f5b
gene_state_by_fit\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80/gene_state_by_polyfun_fit.tsv\tbae8b8fa9954f7c045fe63e6f145b50c014d48ab0df47918fdbb69f178af98d6
susie_pairs_t12\tAnalysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80/susie_pairs_tier12.tsv.gz\t9de174fc8e5e6b777bf3d376b79ca6eb40bdbb7fdaf6099d2210275edeaa571f
", sep = "\t")

# Manifest hashes (audit recipe written into the prespec).
INPUT_MANIFESTS <- data.table::fread(text = "
role\tdir\tcmd\tn_expected\tsha256
eqtl_fits_manifest\tGWAS/finemapping/results/eqtl_susie_polyfun\tfind . -type f -name '*_susie.rds' -printf '%P\\t%s\\n'\t18889\t591431ad13092467c0c58cf5a802b1a37913fae379abf43b6614812bbc8b94ac
broadaway_manifest\tdata/broadaway_eqtl\tfind . -maxdepth 1 -name 'chr*_marginal_summary_results.tsv' -printf '%P\\t%s\\t%Ts\\n'\t22\tb439be99b6ec5c65423a83aa61a2c00a1d9e70086434025a71c8d69084c33548
ld_polyfun_eur\tGWAS/finemapping/data/ld_ref/polyfun_eur\tfind . -type f -printf '%P\\t%s\\n'\t12344\t8a1d3a077a27eb375526320c8481acd227af5da4592f0585e54d97ba38311381
ld_1kg_eur\tGWAS/finemapping/data/ld_ref/1kg_eur\tfind . -type f -printf '%P\\t%s\\n'\t10219\te4e466283426f0ff09a980f1e45990ffbce6d1aaeb08026c37e5019ec25dc84e
ld_1kg_afr_gram\tGWAS/finemapping/data/ld_ref/1kg_afr_gram\tfind . -type f -printf '%P\\t%s\\n'\t5110\teb9cf1a634dee88c447db01f37cb7d19b03de562a0374aa63e3a94aa2e0942ee
ld_1kg_amr_gram\tGWAS/finemapping/data/ld_ref/1kg_amr_gram\tfind . -type f -printf '%P\\t%s\\n'\t5110\t675c01a1127d8817349af9835528bbd37683702ff505d7416a529b117354862f
ld_1kg_eas_gram\tGWAS/finemapping/data/ld_ref/1kg_eas_gram\tfind . -type f -printf '%P\\t%s\\n'\t5110\tf2b5cfecc91de44947ae5242455df4734a5462e768f3828ea55d974b7b9ba7f5
ld_1kg_sas_gram\tGWAS/finemapping/data/ld_ref/1kg_sas_gram\tfind . -type f -printf '%P\\t%s\\n'\t5110\td8f4d909b36388b74078f919f371c9b25ed5faba11ac586ffbce4575a09b0a2f
", sep = "\t", quote = "")
SUMMARY_CONCAT_SHA <- "d6a9c1d678377968a1f7a77b6e996557f10e4060dbd7707d57434b29f2f413c0"

verify_prespec <- function(H80) {
  got <- sha256_file(file.path(H80, "PRESPEC_H80.md"))
  got_plan <- sha256_file(PLAN_PRESPEC)
  ok <- identical(got, H80_PRESPEC_SHA) && identical(got_plan, H80_PRESPEC_SHA)
  log_msg("PRESPEC sha256 run copy ", got, " plan ", got_plan, " expected ", H80_PRESPEC_SHA,
          if (ok) " -> OK" else " -> MISMATCH")
  if (!ok) stop("PRESPEC_H80.md hash mismatch; refusing to run")
  invisible(TRUE)
}

# Returns a table of every prespec input hash with pass flags (S1).
verify_inputs <- function(cores = 1L) {
  f <- copy(INPUT_FILES)
  f[, observed := unlist(mclapply(file.path(ROOT, path), sha256_file, mc.cores = cores))]
  f[, kind := "file"]
  f[, n_observed := NA_integer_][, n_expected := NA_integer_]
  m <- copy(INPUT_MANIFESTS)
  m[, observed := vapply(seq_len(.N), function(i)
    bash_hash(sprintf("cd %s && %s | sort | sha256sum", shQuote(file.path(ROOT, dir[i])), cmd[i])),
    character(1))]
  m[, n_observed := vapply(seq_len(.N), function(i)
    as.integer(system2("bash", c("-c", shQuote(sprintf("cd %s && %s | wc -l", shQuote(file.path(ROOT, dir[i])), cmd[i]))),
                       stdout = TRUE)), integer(1))]
  m[, kind := "manifest"]
  s <- data.table(role = "eqtl_fit_summaries_concat", path = "GWAS/finemapping/results/eqtl_susie_polyfun/chr*_summary.tsv",
                  sha256 = SUMMARY_CONCAT_SHA,
                  observed = bash_hash(sprintf("cd %s && cat chr*_summary.tsv | sha256sum", shQuote(EQTL_FITS))),
                  kind = "concat", n_observed = NA_integer_, n_expected = NA_integer_)
  out <- rbindlist(list(
    f[, .(role, path, kind, expected = sha256, observed, n_expected, n_observed)],
    m[, .(role, path = dir, kind, expected = sha256, observed, n_expected, n_observed)],
    s[, .(role, path, kind, expected = sha256, observed, n_expected, n_observed)]), use.names = TRUE)
  out[, pass := expected == observed & (is.na(n_expected) | n_expected == n_observed)]
  out
}

# Common start of jobs 1-4: prespec hash, every input hash, setup status.
job_preamble <- function(H80, SETUP, OUT, cores = 1L) {
  verify_prespec(H80)
  h <- verify_inputs(cores)
  fwrite(h, file.path(OUT, "input_hash_check.tsv"), sep = "\t")
  if (!all(h$pass)) {
    log_msg("INPUT HASH CHECK FAILED: ", paste(h[pass == FALSE, role], collapse = ","))
    return("STOP:S1")
  }
  log_msg("INPUT HASH CHECK PASSED (", nrow(h), " single-file, manifest and concatenation hashes)")
  st_file <- file.path(SETUP, "STATUS")
  if (!file.exists(st_file)) stop("setup STATUS missing: ", st_file)
  st <- readLines(st_file)[1]
  log_msg("SETUP STATUS: ", st)
  if (!identical(st, "PASS")) return(st)
  log_msg("REPRODUCTION CHECKS PASSED in setup (G1 frame recount, PC-2 eQTL fit reproduction, PC-3 production COLOC reproduction)")
  "PASS"
}

# ---- RNG streams (prespec "Seed") ----------------------------------------------------------------
make_streams <- function(max_index = 3000L) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(H80_SEED)
  s <- get(".Random.seed", envir = .GlobalEnv)
  out <- vector("list", max_index + 1L)
  out[[1]] <- s
  for (i in seq_len(max_index)) {
    s <- parallel::nextRNGStream(s)
    out[[i + 1L]] <- s
  }
  names(out) <- as.character(0:max_index)
  out
}
set_stream <- function(streams, i) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", streams[[as.character(i)]], envir = .GlobalEnv)
}

# Production seeding kept, stream protected (prespec audit patch to seeds): save state,
# Mersenne-Twister set.seed(42), call, restore. Records whether the call drew random numbers.
wrapped_call <- function(f) {
  saved <- get(".Random.seed", envir = .GlobalEnv)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(42)
  after_seed <- get(".Random.seed", envir = .GlobalEnv)
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(f(), error = function(e) structure(list(msg = conditionMessage(e)), class = "h80_error"))
  el <- proc.time()[["elapsed"]] - t0
  drew <- !identical(get(".Random.seed", envir = .GlobalEnv), after_seed)
  assign(".Random.seed", saved, envir = .GlobalEnv)
  restored <- identical(get(".Random.seed", envir = .GlobalEnv), saved) &&
    identical(RNGkind()[1], "L'Ecuyer-CMRG")
  list(res = res, drew = drew, restored = restored, time = el)
}

# ---- production fits (08:220-229 eQTL; 06:724-736 GWAS) -------------------------------------------
fit_eqtl <- function(z, R_reg) {
  wrapped_call(function() susie_rss(
    z = z, R = R_reg, n = EQTL_N, L = SUSIE_L,
    estimate_residual_variance = FALSE, residual_variance = 1,
    check_R = FALSE, max_iter = 500))
}
fit_gwas <- function(z, R_reg, n_eff) {
  wrapped_call(function() susie_rss(
    z = z, R = R_reg, n = n_eff, L = SUSIE_L,
    estimate_residual_variance = FALSE, residual_variance = 1,
    check_R = FALSE, min_abs_corr = 0.5, max_iter = 500))
}
gwas_n_eff <- function(trait_type, N_tot, N_cases) {
  N_ctrl <- if (!is.na(N_cases)) N_tot - N_cases else NA_real_
  if (trait_type == "binary" && !is.na(N_cases) && N_cases > 0 && !is.na(N_ctrl) && N_ctrl > 0) {
    4 / (1 / N_cases + 1 / N_ctrl)
  } else {
    N_tot
  }
}
is_err <- function(x) inherits(x, "h80_error")
fit_ncs <- function(s) if (is_err(s) || is.null(s)) NA_integer_ else length(s$sets$cs)

# ---- gene-pair state, copied from 06_susie_coloc.R:740-791 (production rules) ---------------------
# States: fit_error, not_converged, no_eqtl_cs, no_gwas_cs, coloc_error, untestable,
#         tested_h3 (PP.H4 <= 0.5 and PP.H3 > 0.8), tested_other, positive (PP.H4 > 0.5).
coloc_state <- function(fe, fg) {
  na <- list(state = NA_character_, pp4 = NA_real_, pp3 = NA_real_, n_pairs = NA_integer_,
             n_eligible = NA_integer_, elig_match = NA, coloc_drew = NA, coloc_restored = NA)
  if (is_err(fe$res) || is_err(fg$res) || is.null(fe$res) || is.null(fg$res)) {
    na$state <- "fit_error"; return(na)
  }
  se <- fe$res; sg <- fg$res
  if (!isTRUE(se$converged) || !isTRUE(sg$converged)) { na$state <- "not_converged"; return(na) }
  common <- intersect(colnames(sg$lbf_variable), colnames(se$lbf_variable))
  if (length(common) < MIN_TRIPLE_SNPS) { na$state <- "too_few_common"; return(na) }
  w <- wrapped_call(function() coloc.susie(sg, se))
  na$coloc_drew <- w$drew; na$coloc_restored <- w$restored
  sr <- w$res
  if (is_err(sr)) { na$state <- "coloc_error"; return(na) }
  if (is.null(sr$summary) || nrow(sr$summary) == 0 || !("PP.H4.abf" %in% names(sr$summary))) {
    na$state <- if (length(se$sets$cs) == 0) "no_eqtl_cs" else if (length(sg$sets$cs) == 0) "no_gwas_cs" else "coloc_null"
    return(na)
  }
  pt <- coloc_susie_pair_table(sg, se, sr)
  kept <- pt$kept
  na$n_pairs <- nrow(pt$pairs)
  na$n_eligible <- sum(pt$pairs$eligible %in% TRUE)
  na$elig_match <- pt$matches_coloc
  if (coloc_susie_method(kept) == "susie") {
    best <- kept[which.max(kept$PP.H4.abf), ]
    na$pp4 <- best$PP.H4.abf
    na$pp3 <- best$PP.H3.abf
    na$state <- if (na$pp4 > 0.5) "positive" else if (na$pp3 > 0.8) "tested_h3" else "tested_other"
  } else {
    na$state <- "untestable"
  }
  na
}

# ---- LD loading ------------------------------------------------------------------------------------
# Column-select copy of the production cached loader (src/perf/get_ld_per_locus_cached.R), which is
# itself value-identical to finemapping_functions.R::get_ld_per_locus. Changes, and only these:
#   (1) the panel directory is passed explicitly (production resolves it from env vars);
#   (2) only the LD columns/rows of the matched variants are read (fread select), because 16
#       parallel workers reading whole blocks (median 4.5 GB, max 17.7 GB dense) exceed the job's
#       memory; the matched submatrix is the only part production keeps;
#   (3) the LD block of every variant is carried along (needed for the per-block eigendecomposition).
# Value identity with the unmodified production loader is tested in the setup job (D5).
ld_blocks_for <- function(CHR, START, END, panel_dir) {
  blocks_file <- file.path(panel_dir, "approx_LD_blocks.txt")
  if (!file.exists(blocks_file)) stop(paste("LD blocks file not found:", blocks_file))
  ab <- fread(blocks_file)
  m <- ab[chr == CHR & start < END & stop > START]
  if (nrow(m) == 0) return(NULL)
  file.path(panel_dir, paste0("chr", CHR), paste0(m$start, ".", m$stop), paste0(m$start, ".", m$stop))
}

ld_plan <- function(ss_per_locus, CHR, START, END, panel_dir) {
  path_to_LD <- ld_blocks_for(CHR, START, END, panel_dir)
  if (is.null(path_to_LD)) return(NULL)
  bim <- data.frame(); blk <- integer(0); rowin <- integer(0); used <- character(0)
  for (b in seq_along(path_to_LD)) {
    file <- path_to_LD[b]
    bim_file <- paste0(file, ".bim"); ld_file <- paste0(file, ".ld")
    if (!file.exists(bim_file) || !file.exists(ld_file)) { warning(paste("LD files not found:", file)); next }
    BIM <- fread(bim_file)
    bim <- rbind(bim, BIM)
    blk <- c(blk, rep(b, nrow(BIM))); rowin <- c(rowin, seq_len(nrow(BIM)))
    used <- c(used, file)
  }
  if (nrow(bim) == 0) return(NULL)
  colnames(bim) <- c("chr", "rsid", "dk", "pos", "alt", "ref")
  bim$SNP <- paste(bim$chr, bim$pos, bim$alt, bim$ref, sep = ":")
  keep_uniq <- !duplicated(bim$SNP)
  if (!all(keep_uniq)) { bim <- bim[keep_uniq]; blk <- blk[keep_uniq]; rowin <- rowin[keep_uniq] }

  # ---- from here: production matching code, verbatim except index bookkeeping ----
  ss_per_locus$SNP1 <- paste(ss_per_locus$chromosome, ss_per_locus$position,
                             ss_per_locus$allele1, ss_per_locus$allele2, sep = ":")
  ss_per_locus$SNP2 <- paste(ss_per_locus$chromosome, ss_per_locus$position,
                             ss_per_locus$allele2, ss_per_locus$allele1, sep = ":")
  bim$ss_matched <- bim$SNP %in% c(ss_per_locus$SNP2, ss_per_locus$SNP1)
  n_matched <- sum(bim$ss_matched)
  if (n_matched < 10) return(NULL)
  idx <- which(bim$ss_matched)
  bim_filtered <- bim[bim$ss_matched, ]
  ss_filtered <- rbind(
    ss_per_locus %>% dplyr::filter(SNP1 %in% bim_filtered$SNP) %>% dplyr::rename(SNP = SNP1) %>% dplyr::select(!c(SNP2)),
    ss_per_locus %>% dplyr::filter(SNP2 %in% bim_filtered$SNP) %>% dplyr::mutate(beta = -1 * beta) %>%
      dplyr::rename(SNP = SNP2, allele1 = allele2, allele2 = allele1) %>% dplyr::select(!c(SNP1))
  )
  ss_filtered <- ss_filtered[!duplicated(ss_filtered$SNP), ]
  ss_filtered <- ss_filtered[match(bim_filtered$SNP, ss_filtered$SNP), ]
  ss_filtered <- ss_filtered[!is.na(ss_filtered$SNP), ]
  keep <- bim_filtered$SNP %in% ss_filtered$SNP
  if (sum(keep) < 10) return(NULL)
  bim_filtered <- bim_filtered[keep, ]
  idx <- idx[keep]
  ss_filtered <- ss_filtered[match(bim_filtered$SNP, ss_filtered$SNP), ]
  list(ss = ss_filtered, bim = bim_filtered, blk = blk[idx], rowin = rowin[idx],
       paths = path_to_LD, files_used = used)
}

ld_read <- function(plan, nthread = 1L) {
  n <- nrow(plan$bim)
  R <- matrix(0, n, n)
  for (b in sort(unique(plan$blk))) {
    pos_b <- which(plan$blk == b)
    rows_b <- plan$rowin[pos_b]
    cols <- sort(unique(rows_b))
    LD <- fread(paste0(plan$paths[b], ".ld"), select = cols, nThread = nthread)
    M <- as.matrix(LD)
    storage.mode(M) <- "double"
    R[pos_b, pos_b] <- M[rows_b, match(rows_b, cols), drop = FALSE]
    rm(LD, M)
  }
  colnames(R) <- rownames(R) <- plan$bim$SNP
  ss_filtered <- plan$ss
  blk <- plan$blk
  n_na_before <- sum(is.na(R))
  # Production NA handling, verbatim (iteratively remove high-NA variants).
  if (sum(is.na(R)) != 0) {
    ld_mat <- R
    for (iter in seq_len(100)) {
      na_per_var <- colSums(is.na(ld_mat))
      if (max(na_per_var) == 0) break
      bad <- na_per_var > (nrow(ld_mat) * 0.5)
      if (sum(bad) == 0) bad <- na_per_var > 0
      ld_mat <- ld_mat[!bad, !bad, drop = FALSE]
      ss_filtered <- ss_filtered[!bad, ]
      blk <- blk[!bad]
      if (nrow(ld_mat) < 10) break
    }
    R <- ld_mat
    if (nrow(R) < 10) return(NULL)
  }
  list(ss = ss_filtered, R = R, blk = blk, n_na = n_na_before, files_used = plan$files_used)
}

ld_load_fast <- function(ss, CHR, START, END, panel_dir, nthread = 1L) {
  p <- ld_plan(ss, CHR, START, END, panel_dir)
  if (is.null(p)) return(NULL)
  ld_read(p, nthread)
}

`%>%` <- dplyr::`%>%`

# Production loader in an isolated environment (sourced read-only from the repository).
load_production_loader <- function() {
  e <- new.env(parent = globalenv())
  suppressPackageStartupMessages({
    sys.source(file.path(FM, "src/finemapping_functions.R"), envir = e)
    sys.source(file.path(FM, "src/perf/get_ld_per_locus_cached.R"), envir = e)
  })
  e
}

# eQTL ss table exactly as 08 builds it (08:149-158) from a gene's Broadaway rows.
eqtl_ss_from_rows <- function(gene_dt) {
  data.frame(chromosome = gene_dt$CHR, position = gene_dt$POS, allele1 = gene_dt$NEA,
             allele2 = gene_dt$EA, beta = gene_dt$Beta, standard_error = gene_dt$SE,
             stringsAsFactors = FALSE)
}
eqtl_window <- function(gene_dt) {
  START <- min(gene_dt$POS) - CIS_WINDOW_BP
  END <- max(gene_dt$POS) + CIS_WINDOW_BP
  c(max(START, 1L), END)
}

# ---- per-block eigendecomposition and noise ---------------------------------------------------------
# Sigma = clipped(R + 1e-3 I) per LD block; noise = V diag(sqrt(d)) eps (prespec audit clarification).
block_eigen <- function(R, blk) {
  out <- list(); n_neg <- 0L; min_eig <- Inf; trS <- 0; trS2 <- 0
  for (b in sort(unique(blk))) {
    ib <- which(blk == b)
    S <- R[ib, ib, drop = FALSE]
    diag(S) <- diag(S) + LD_REGULARIZE
    e <- eigen(S, symmetric = TRUE)
    d <- e$values
    min_eig <- min(min_eig, min(d))
    n_neg <- n_neg + sum(d < 0)
    d[d < 0] <- 0
    trS <- trS + sum(d); trS2 <- trS2 + sum(d^2)
    out[[length(out) + 1L]] <- list(ib = ib, U = sweep(e$vectors, 2, sqrt(d), "*"))
    rm(S, e)
  }
  list(blocks = out, n_neg = n_neg, min_eig = min_eig, trS = trS, trS2 = trS2,
       n_blocks = length(out))
}
noise_from <- function(be, eps) {
  z <- numeric(length(eps))
  for (bl in be$blocks) z[bl$ib] <- as.vector(bl$U %*% eps[bl$ib])
  z
}
reg_matrix <- function(R, ids) {
  Rr <- R
  diag(Rr) <- diag(Rr) + LD_REGULARIZE
  colnames(Rr) <- rownames(Rr) <- ids
  Rr
}
snp_ids_of <- function(ss) paste0(ss$chromosome, ":", ss$position)

fit_row <- function(w, side, lambda, causal_type, z, causal_idx, fr) {
  s <- fr$res
  data.table(side = side, lambda = lambda, causal_type = causal_type,
             error = if (is_err(s)) s$msg else NA_character_,
             converged = if (is_err(s)) NA else isTRUE(s$converged),
             n_cs = fit_ncs(s),
             cs_sizes = if (is_err(s) || length(s$sets$cs) == 0) "" else paste(lengths(s$sets$cs), collapse = ";"),
             niter = if (is_err(s)) NA_integer_ else as.integer(s$niter),
             lead_absz = max(abs(z)), z_causal = z[causal_idx],
             drew_rng = fr$drew, restored = fr$restored, time_s = fr$time)
}

# Production shared-posterior helper, copied into scripts/ (sha256 checked against the repository copy).
H80_SCRIPTS <- Sys.getenv("H80_SCRIPTS", unset = "")
if (nzchar(H80_SCRIPTS)) {
  stopifnot(identical(sha256_file(file.path(H80_SCRIPTS, "coloc_signal_eligibility.R")),
                      INPUT_FILES[role == "eligibility_helper", sha256]))
  source(file.path(H80_SCRIPTS, "coloc_signal_eligibility.R"))
}
