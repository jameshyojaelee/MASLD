#!/usr/bin/env Rscript
# H23 stage 1: EUR-only SuSiE-RSS refits at the loci that already have a MESuSiE
# multi-ancestry fit, and per-fit summaries (lead, credible set, top gene set).
# Prespecification: Analysis/Hypotheses/H23-20261009T140416Z/PRESPEC_H23.md
#
# Usage: Rscript h23_fit_units.R <new_output_dir> [smoke]
#   smoke = first 3 units per (set, trait) only.
# Env: H23_CORES (default 8).
#
# Method note: susie_rss is the native summary-statistic fine-mapping method for GWAS
# z-scores with an LD reference. Settings mirror src/10b_run_mesusie.R (L=10, LD
# regularisation 1e-3, estimate_residual_variance=FALSE, max_iter=200, min_abs_corr 0.5)
# so that the only difference between the EUR-only and multi-ancestry fits is the
# non-EUR data. Allele alignment and block-diagonal LD follow 10b_run_mesusie.R.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(susieR)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: h23_fit_units.R <new_output_dir> [smoke]")
out_dir <- args[[1]]
smoke   <- length(args) >= 2L && identical(args[[2]], "smoke")
if (file.exists(out_dir)) stop("Refusing to overwrite existing output dir: ", out_dir)
dir.create(out_dir, recursive = TRUE)
n_cores <- as.integer(Sys.getenv("H23_CORES", unset = "8"))

SEED <- 20260923L
set.seed(SEED)
ROOT   <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FM_DIR <- file.path(ROOT, "GWAS/finemapping")
Sys.setenv(MASLD_PROJECT_ROOT = ROOT, LD_PANEL = "polyfun")
Sys.unsetenv("UKBB_LD_DIR")
suppressWarnings(suppressMessages(source(file.path(FM_DIR, "src/finemapping_functions.R"))))

SUSIE_L       <- 10L
LD_REGULARIZE <- 1e-3
MIN_SNPS      <- 50L
MIN_LOCI_FIT  <- 10L   # stop rule: (set, trait) with fewer multi-ancestry loci is not fitted

SETS <- list(
  MVP  = list(dir = "mesusie_mvp", loci = "results/susiex_mvp/shared_loci.csv", pcol = "anchor_p",
              traits = c("ALT", "AST", "Albumin", "Platelet", "NAFLD", "ChronLiver", "Cirrhosis")),
  `4WAY` = list(dir = "mesusie_4way", loci = "results/susiex/shared_loci.csv", pcol = "eur_lead_p",
                traits = c("ALT", "AST", "GGT")))

# ---- registry ---------------------------------------------------------------
reg_lines <- readLines(file.path(FM_DIR, "config/gwas_registry.tsv"))
reg_lines <- reg_lines[!grepl("^\\s*#", reg_lines) & nchar(trimws(reg_lines)) > 0]
REG <- fread(text = paste(reg_lines, collapse = "\n"), sep = "\t")
reg_n    <- function(study) { x <- REG[study_name == study, N_tot]; if (length(x) == 1L) as.numeric(x) else NA_real_ }
reg_path <- function(study) file.path(FM_DIR, REG[study_name == study, sumstats_path][1])

study_B <- function(A, trait) {
  if (!trait %in% c("ALT", "AST")) return(NA_character_)
  if (grepl("^UKBB_", A)) return(paste0("MVP_", trait, "_EUR"))
  if (grepl("^MVP_", A))  return(paste0("UKBB_", trait))
  NA_character_
}

# ---- 1. enumerate multi-ancestry units -------------------------------------
parse_snp <- function(ids) {
  p <- tstrsplit(ids, ":", fixed = TRUE)
  data.table(chr = as.integer(p[[1]]), pos = as.integer(p[[2]]), a1 = p[[3]], a2 = p[[4]])
}
cs_to_idx <- function(cs, snp_ids) {
  if (is.null(cs) || length(cs) == 0L) return(list())
  lapply(cs, function(ix) if (is.numeric(ix)) as.integer(ix) else match(ix, snp_ids))
}

units <- list(); inv <- list()
for (sn in names(SETS)) {
  S <- SETS[[sn]]
  sl <- fread(file.path(FM_DIR, S$loci))
  ptab <- sl[, .(locus_id, anchor_p = as.numeric(get(S$pcol)))]
  for (tr in S$traits) {
    files <- list.files(file.path(FM_DIR, "results", S$dir, tr), pattern = "_mesusie\\.rds$", full.names = TRUE)
    n_files <- length(files)
    if (n_files < MIN_LOCI_FIT) {
      inv[[length(inv) + 1L]] <- data.table(set = sn, trait = tr, n_rds = n_files, n_with_eur_arm = NA_integer_,
                                             n_after_window_collapse = NA_integer_, n_ma_ok = NA_integer_,
                                             note = "not fitted: fewer than 10 multi-ancestry loci (stop rule)")
      next
    }
    rows <- lapply(files, function(f) {
      o <- readRDS(f)
      if (!"EUR" %in% o$ancestries) return(NULL)
      sid <- o$snp_ids
      cs_idx <- cs_to_idx(o$credible_sets, sid)
      non_eur <- setdiff(names(o$studies), "EUR")
      data.table(set = sn, trait = tr, locus_id = o$locus_id, chr = o$chr,
                 win_start = o$window[1], win_end = o$window[2], A = unname(o$studies["EUR"]),
                 n_noneur = sum(vapply(non_eur, function(a) reg_n(unname(o$studies[a])), numeric(1))),
                 ma_converged = isTRUE(o$converged), ma_n_cs = length(cs_idx),
                 ma_arms = paste(o$ancestries, collapse = ","),
                 ma_snp_ids = list(sid), ma_pip = list(as.numeric(o$pip)), ma_cs_idx = list(cs_idx))
    })
    n_eur <- sum(!vapply(rows, is.null, logical(1)))
    dt <- rbindlist(Filter(Negate(is.null), rows))
    dt <- merge(dt, ptab, by = "locus_id", all.x = TRUE, sort = FALSE)
    dt[is.na(anchor_p), anchor_p := 1]
    setorder(dt, anchor_p, locus_id)
    keep <- rep(FALSE, nrow(dt)); kept <- list()
    for (i in seq_len(nrow(dt))) {
      clash <- any(vapply(kept, function(k) k$chr == dt$chr[i] && k$s <= dt$win_end[i] && k$e >= dt$win_start[i], logical(1)))
      if (!clash) { keep[i] <- TRUE; kept[[length(kept) + 1L]] <- list(chr = dt$chr[i], s = dt$win_start[i], e = dt$win_end[i]) }
    }
    dt <- dt[keep]
    dt[, ma_ok := ma_converged & ma_n_cs >= 1L]
    inv[[length(inv) + 1L]] <- data.table(set = sn, trait = tr, n_rds = n_files, n_with_eur_arm = n_eur,
                                           n_after_window_collapse = nrow(dt), n_ma_ok = sum(dt$ma_ok), note = "")
    units[[length(units) + 1L]] <- dt
  }
}
units <- rbindlist(units)
inventory_fit <- rbindlist(inv)
if (smoke) units <- units[, head(.SD, 3L), by = .(set, trait)]
cat(sprintf("Units: %d (smoke=%s)\n", nrow(units), smoke))
units[, uid := .I]
units[, B := unname(mapply(study_B, A, trait))]
units[, nA := vapply(A, reg_n, numeric(1))]
units[, nB := vapply(ifelse(is.na(B), "", B), reg_n, numeric(1))]

# ---- 2. summary statistics windows -----------------------------------------
read_ss_windows <- function(study, win) {
  p <- reg_path(study)
  ss <- fread(p, select = c("chromosome", "position", "allele1", "allele2", "beta", "se"))
  setnames(ss, c("CHR", "POS", "A1", "A2", "BETA", "SE"))
  ss[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS), BETA = as.numeric(BETA), SE = as.numeric(SE))]
  ss <- ss[!is.na(CHR) & !is.na(POS) & !is.na(BETA) & !is.na(SE) & SE > 0]
  out <- rbindlist(lapply(unique(win$chr), function(ch) {
    w <- win[chr == ch]
    ss[CHR == ch][inrange(POS, w$s, w$e)]
  }))
  setkey(out, CHR, POS)
  cat(sprintf("  [%s] %d variants in windows\n", study, nrow(out)))
  out
}
need <- rbindlist(list(units[, .(study = A, chr, s = win_start, e = win_end)],
                       units[!is.na(B), .(study = B, chr, s = win_start, e = win_end)]))
SS <- list()
for (st in unique(need$study)) SS[[st]] <- read_ss_windows(st, need[study == st, .(chr, s, e)])

# ---- 3. gene map ------------------------------------------------------------
# The Broadaway table lists every variant tested in a +-1 Mb cis window of each gene (ALDH2: 111,204,880-113,204,494),
# not only associated pairs. gm_all keeps every tested pair (the original, prespecified map). gm_eqtl keeps pairs with
# PVAL < EQTL_P and is the map for the primary gene endpoint (Addendum B of the prespecification).
EQTL_P <- 1e-5
allw <- units[, .(chr, s = win_start, e = win_end)]
bfiles <- list.files(file.path(ROOT, "data/broadaway_eqtl"), pattern = "_marginal_summary_results\\.tsv$", full.names = TRUE)
tss_list <- list()
gm_raw <- rbindlist(lapply(bfiles, function(f) {
  x <- fread(f, select = c("GeneSymbol", "CHR", "POS", "NEA", "EA", "PVAL"))
  x[, CHR := suppressWarnings(as.integer(CHR))]
  x[, PVAL := suppressWarnings(as.numeric(PVAL))]   # fread reads PVAL as character in chr1, chr7 and chr8; a character PVAL makes "PVAL < 1e-5" a string comparison
  x <- x[!is.na(CHR)]
  # TSS estimate per gene = midpoint of its cis-tested range (the cis window is symmetric around the TSS)
  tss_list[[basename(f)]] <<- x[, .(cis_min = min(POS), cis_max = max(POS), n_tested = .N), by = .(GeneSymbol, chr = CHR)]
  rbindlist(lapply(unique(allw$chr), function(ch) { w <- allw[chr == ch]; x[CHR == ch][inrange(POS, w$s, w$e)] }))
}))
stopifnot(is.numeric(gm_raw$PVAL), !anyNA(gm_raw$PVAL))
gene_tss <- rbindlist(tss_list); gene_tss[, tss := (cis_min + cis_max) / 2]
fwrite(gene_tss, file.path(out_dir, "gene_tss_estimates.tsv.gz"), sep = "\t")
aldh2_rng <- gene_tss[GeneSymbol == "ALDH2", .(chr, cis_min, cis_max, n_tested)]
if (nrow(aldh2_rng) == 1L) fwrite(aldh2_rng, file.path(out_dir, "aldh2_cis_range.tsv"), sep = "\t")
setkey(gene_tss, chr, tss)
nearest_gene <- function(ch, pos) {
  g <- gene_tss[chr == ch]
  if (nrow(g) == 0L || is.na(pos)) return(NA_character_)
  g$GeneSymbol[which.min(abs(g$tss - pos))]
}
gm_raw[, akey := paste(CHR, POS, ifelse(EA < NEA, EA, NEA), ifelse(EA < NEA, NEA, EA), sep = ":")]
gm <- unique(gm_raw[, .(GeneSymbol, akey, CHR, POS)]); setkey(gm, akey)
gm_e <- unique(gm_raw[PVAL < EQTL_P, .(GeneSymbol, akey, CHR, POS)]); setkey(gm_e, akey)
cat(sprintf("Gene-variant pairs in windows: %d tested, %d with PVAL < %g (%.1f%%)\n", nrow(gm), nrow(gm_e), EQTL_P, 100 * nrow(gm_e) / nrow(gm)))
if (nrow(gm_e) / nrow(gm) > 0.2) stop("Share of eQTL pairs below the threshold is implausibly high; check PVAL parsing")
fwrite(gm_e, file.path(out_dir, "gene_links_eqtl_in_windows.tsv.gz"), sep = "\t")
rm(gm_raw)
mk_akey <- function(chr, pos, a1, a2) paste(chr, pos, ifelse(a1 < a2, a1, a2), ifelse(a1 < a2, a2, a1), sep = ":")

# ---- 4. LD, alignment, SuSiE-RSS -------------------------------------------
# Read only the LD rows/columns inside [lo, hi] (the union of the group's windows). Results for the retained
# variants are identical to loading the whole block, as in 10b_run_mesusie.R; blocks are joined block-diagonally.
load_ld <- function(paths, lo, hi) {
  bims <- list(); mats <- list()
  for (file in paths) {
    bf <- paste0(file, ".bim"); lf <- paste0(file, ".ld")
    if (!file.exists(bf) || !file.exists(lf)) next
    BIM <- fread(bf)
    setnames(BIM, c("chr", "rsid", "dk", "pos", "alt", "ref"))
    idx <- which(BIM$pos >= lo & BIM$pos <= hi)
    if (length(idx) == 0L) next
    LD <- fread(lf, select = idx)
    m <- as.matrix(LD)[idx, , drop = FALSE]
    rm(LD)
    bims[[length(bims) + 1L]] <- BIM[idx]; mats[[length(mats) + 1L]] <- m
  }
  if (length(bims) == 0L) return(NULL)
  bim_all <- rbindlist(bims)
  m <- if (length(mats) == 1L) mats[[1]] else as.matrix(bdiag(mats))
  bim_all[, SNP := paste(chr, pos, alt, ref, sep = ":")]
  rownames(m) <- colnames(m) <- bim_all$SNP
  k <- !duplicated(bim_all$SNP)
  list(bim = bim_all[k], ld = m[k, k, drop = FALSE])
}

align_to_ld <- function(ss, ld, keep_pos) {
  s <- ss[POS %in% keep_pos]
  if (nrow(s) == 0L) return(NULL)
  s[, `:=`(fwd = paste(CHR, POS, A1, A2, sep = ":"), rev = paste(CHR, POS, A2, A1, sep = ":"))]
  bim <- ld$bim
  f <- s[fwd %in% bim$SNP]; f[, SNP := fwd]
  r <- s[rev %in% bim$SNP & !(fwd %in% bim$SNP)]
  if (nrow(r) > 0L) { r[, BETA := -BETA]; r[, SNP := rev] }
  a <- rbind(f, r, fill = TRUE)
  a <- a[!duplicated(SNP)][!duplicated(POS)]
  if (nrow(a) == 0L) return(NULL)
  idx <- match(a$SNP, bim$SNP)
  o <- order(idx); a <- a[o]; idx <- idx[o]
  R <- ld$ld[idx, idx, drop = FALSE]
  # alleles in the LD orientation (alt, ref), used for gene-map keys
  a[, `:=`(L1 = bim$alt[idx], L2 = bim$ref[idx])]
  list(ss = a, R = R)
}

fit_susie <- function(al, n) {
  if (is.null(al) || nrow(al$ss) < MIN_SNPS) return(NULL)
  R <- al$R + LD_REGULARIZE * diag(nrow(al$R))
  if (any(is.na(R))) { R[is.na(R)] <- 0; R <- (R + t(R)) / 2 }
  z <- al$ss$BETA / al$ss$SE
  fit <- tryCatch(susie_rss(z = z, R = R, n = n, L = SUSIE_L, estimate_residual_variance = FALSE,
                            max_iter = 200, coverage = 0.95, min_abs_corr = 0.5),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(pos = al$ss$POS, chr = al$ss$CHR, a1 = al$ss$L1, a2 = al$ss$L2, z = z, pip = as.numeric(fit$pip),
       cs = lapply(fit$sets$cs, as.integer), converged = isTRUE(fit$converged))
}

gene_set_for <- function(chr, pos, a1, a2, pip, idx, gmap = gm_e) {
  if (length(idx) == 0L) return(character(0))
  k <- mk_akey(chr[idx], pos[idx], a1[idx], a2[idx])
  hit <- gmap[.(k), nomatch = NULL, on = "akey"]
  if (nrow(hit) == 0L) return(character(0))
  pv <- setNames(pip[idx], k)
  hit[, pip := pv[akey]]
  top <- hit[pip >= max(pip) - 1e-12]
  sort(unique(top$GeneSymbol))
}

# Summarise a fit on a common footing: lead, credible set containing the lead, top gene set.
summarise_fit <- function(chr, pos, a1, a2, pip, z, cs_idx, converged) {
  out <- list(ok = FALSE, n_var = length(pos), n_cs = length(cs_idx), converged = converged,
              lead_pos = NA_integer_, lead_pip = NA_real_, cs_size = NA_integer_, n_ties = NA_integer_, nearest_gene = NA_character_,
              cs_pos = list(integer(0)), genes = list(character(0)), genes_all = list(character(0)))
  if (!isTRUE(converged) || length(cs_idx) == 0L || length(pip) == 0L) return(out)
  mx <- max(pip, na.rm = TRUE)
  cand <- which(pip >= mx - 1e-12)
  out$n_ties <- length(cand)   # variants tied at the maximum PIP; the lead is then the one with the larger |z|
  out$nearest_gene <- NA_character_
  lead <- if (length(cand) == 1L) cand else cand[which.max(abs(z[cand]))]
  in_cs <- which(vapply(cs_idx, function(ix) lead %in% ix, logical(1)))
  if (length(in_cs) == 0L) { out$lead_pos <- pos[lead]; out$lead_pip <- pip[lead]; return(out) }
  j <- in_cs[which.min(vapply(cs_idx[in_cs], length, integer(1)))]
  idx <- cs_idx[[j]]
  out$ok <- TRUE; out$lead_pos <- pos[lead]; out$lead_pip <- pip[lead]
  out$nearest_gene <- nearest_gene(chr[lead], pos[lead])   # sensitivity: expressed gene with the nearest estimated TSS
  out$cs_size <- length(idx); out$cs_pos <- list(pos[idx])
  out$genes <- list(gene_set_for(chr, pos, a1, a2, pip, idx, gm_e))          # primary: eQTL pairs with PVAL < EQTL_P
  out$genes_all <- list(gene_set_for(chr, pos, a1, a2, pip, idx, gm))        # sensitivity: every tested pair
  out
}

# Post hoc nested EUR comparator: z_M = (sqrt(N_A) z_A + (n_add / sqrt(N_B)) z_B) / sqrt(N_A + n_add). With n_add = N_B this is the
# usual sample-size-weighted z meta-analysis (valid when the two cohorts use different effect units). With n_add < N_B the z of B is
# scaled by sqrt(n_add / N_B), the expected z of a subsample of size n_add; its noise variance is then below 1 (about 0.8 here), so
# the matched comparator slightly understates the change that a real subsample would give.
nested_fit <- function(alAb, alB, nA, nB, n_add) {
  m <- merge(alAb$ss[, .(POS, SNP, zA = BETA / SE)], alB$ss[, .(POS, SNP, zB = BETA / SE)], by = c("POS", "SNP"), sort = FALSE)
  if (nrow(m) < MIN_SNPS) return(NULL)
  ix <- match(m$SNP, alAb$ss$SNP)
  ss <- copy(alAb$ss)[ix]
  ss[, `:=`(BETA = (sqrt(nA) * m$zA + (n_add / sqrt(nB)) * m$zB) / sqrt(nA + n_add), SE = 1)]
  fit_susie(list(ss = ss, R = alAb$R[ix, ix, drop = FALSE]), nA + n_add)
}

prefix_list <- function(x, p) setNames(x, paste0(p, names(x)))

process_unit <- function(u, ld) {
  chr <- u$chr; s <- u$win_start; e <- u$win_end
  snp <- u$ma_snp_ids[[1]]; mp <- parse_snp(snp)
  V_ma <- mp$pos
  ssA <- SS[[u$A]][CHR == chr & POS >= s & POS <= e]
  res <- list(uid = u$uid)
  # multi-ancestry fit (stored); z for tie-breaks comes from reference A at the same positions
  ma <- NULL
  alA <- align_to_ld(copy(ssA), ld, V_ma)
  A_ma <- if (!is.null(alA)) fit_susie(alA, u$nA) else NULL
  zA <- rep(0, length(V_ma))
  if (!is.null(A_ma)) { zA <- as.numeric(setNames(A_ma$z, A_ma$pos)[as.character(mp$pos)]); zA[is.na(zA)] <- 0 }
  ma <- summarise_fit(mp$chr, mp$pos, mp$a1, mp$a2, u$ma_pip[[1]], zA, u$ma_cs_idx[[1]], u$ma_converged)
  A1s <- if (!is.null(A_ma)) summarise_fit(A_ma$chr, A_ma$pos, A_ma$a1, A_ma$a2, A_ma$pip, A_ma$z, A_ma$cs, A_ma$converged) else NULL
  res$ma <- ma; res$A_ma <- A1s
  res$n_VMA <- length(V_ma); res$n_A_aligned <- if (is.null(alA)) 0L else nrow(alA$ss)
  # independent EUR baseline (ALT, AST only)
  res$A_b <- NULL; res$B <- NULL; res$n_VB <- NA_integer_
  if (!is.na(u$B) && !is.null(alA)) {
    ssB <- SS[[u$B]][CHR == chr & POS >= s & POS <= e]
    alB <- align_to_ld(copy(ssB), ld, alA$ss$POS)
    if (!is.null(alB)) {
      V_b <- alB$ss$POS
      res$n_VB <- length(V_b)
      B_fit <- fit_susie(alB, u$nB)
      same_set <- identical(sort(V_b), sort(alA$ss$POS))
      alAb <- if (same_set) alA else align_to_ld(copy(ssA), ld, V_b)
      A_b <- if (same_set) A_ma else if (!is.null(alAb)) fit_susie(alAb, u$nA) else NULL
      if (!is.null(B_fit) && !is.null(A_b)) {
        res$B   <- summarise_fit(B_fit$chr, B_fit$pos, B_fit$a1, B_fit$a2, B_fit$pip, B_fit$z, B_fit$cs, B_fit$converged)
        res$A_b <- summarise_fit(A_b$chr, A_b$pos, A_b$a1, A_b$a2, A_b$pip, A_b$z, A_b$cs, A_b$converged)
        # Post hoc nested EUR comparators (Addendum D): study A plus study B by sample-size-weighted z meta-analysis.
        # Mf adds all of B; Mm adds B at the sample size of the non-EUR arms in this unit's multi-ancestry fit.
        if (!is.null(alAb)) {
          for (arm in c("Mf", "Mm")) {
            n_add <- if (arm == "Mf") u$nB else u$n_noneur
            if (!is.finite(n_add) || n_add <= 0) next
            fm <- nested_fit(alAb, alB, u$nA, u$nB, n_add)
            if (!is.null(fm)) res[[arm]] <- summarise_fit(fm$chr, fm$pos, fm$a1, fm$a2, fm$pip, fm$z, fm$cs, fm$converged)
          }
        }
      }
    }
  }
  res
}

# group units by LD block set so each block is read once
get_blocks <- function(chr, s, e) {
  x <- NULL
  invisible(capture.output(x <- suppressWarnings(tryCatch(find_LD_block(chr, s, e, "EUR"), error = function(err) NULL))))
  x
}
units[, ld_key := vapply(seq_len(.N), function(i) {
  x <- get_blocks(chr[i], win_start[i], win_end[i])
  if (is.null(x)) NA_character_ else paste(x, collapse = "|")
}, character(1))]
cat(sprintf("Units without an LD block: %d\n", sum(is.na(units$ld_key))))
groups <- split(units[!is.na(ld_key)], by = "ld_key", keep.by = TRUE)
cat(sprintf("LD block groups: %d, cores: %d\n", length(groups), n_cores))

run_group <- function(g) {
  paths <- strsplit(g$ld_key[1], "|", fixed = TRUE)[[1]]
  ld <- tryCatch(load_ld(paths, min(g$win_start), max(g$win_end)), error = function(e) NULL)
  lapply(seq_len(nrow(g)), function(i) {
    u <- g[i]
    if (is.null(ld)) return(list(uid = u$uid, error = "LD load failed"))
    tryCatch(process_unit(u, ld), error = function(e) list(uid = u$uid, error = conditionMessage(e)))
  })
}
t0 <- Sys.time()
glist <- mclapply(groups, run_group, mc.cores = n_cores, mc.preschedule = FALSE)
cat(sprintf("Fits done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
bad_groups <- which(vapply(glist, function(g) inherits(g, "try-error") || is.null(g), logical(1)))
cat(sprintf("Groups that failed outright (killed or error): %d\n", length(bad_groups)))
fits <- unlist(glist[setdiff(seq_along(glist), bad_groups)], recursive = FALSE)
errs <- Filter(function(x) !is.null(x$error), fits)
cat(sprintf("Unit errors: %d\n", length(errs)))
if (length(errs)) print(rbindlist(lapply(errs, as.data.table), fill = TRUE))

# ---- 5. outputs -------------------------------------------------------------
saveRDS(list(units = units[, !c("ma_snp_ids", "ma_pip", "ma_cs_idx"), with = FALSE], fits = fits,
             inventory_fit = inventory_fit, seed = SEED, smoke = smoke),
        file.path(out_dir, "unit_fits.rds"))
fwrite(inventory_fit, file.path(out_dir, "inventory_fit.tsv"), sep = "\t")
fwrite(units[, !c("ma_snp_ids", "ma_pip", "ma_cs_idx")], file.path(out_dir, "units_enumerated.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Wrote", out_dir, "\n")
