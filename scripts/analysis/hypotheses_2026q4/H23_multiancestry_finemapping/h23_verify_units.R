#!/usr/bin/env Rscript
# H23 check: re-derive per-unit quantities from raw inputs by a separate code path.
#   A. multi-ancestry side: lead (max PIP) and credible set size straight from the MESuSiE RDS.
#   B. EUR-only side: refit a few units from the raw summary statistics and the whole LD block
#      (no row/column subsetting), then compare lead position and credible set size.
# Usage: Rscript h23_verify_units.R <analysis_dir> <fits_dir> <new_output_dir> [n_refit] [comma-separated uids]
suppressPackageStartupMessages({ library(data.table); library(susieR); library(Matrix) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: h23_verify_units.R <analysis_dir> <fits_dir> <new_output_dir> [n_refit]")
an_dir <- args[[1]]; fits_dir <- args[[2]]; out_dir <- args[[3]]
n_refit <- if (length(args) >= 4L) as.integer(args[[4]]) else 3L
if (file.exists(out_dir)) stop("Refusing to overwrite existing output dir: ", out_dir)
dir.create(out_dir, recursive = TRUE)
set.seed(20260923L)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FM <- file.path(ROOT, "GWAS/finemapping")
U <- fread(file.path(an_dir, "units.tsv"))
UE <- fread(file.path(fits_dir, "units_enumerated.tsv"))
dirs <- c(MVP = "mesusie_mvp", `4WAY` = "mesusie_4way")

# ---- A. multi-ancestry side from the raw RDS ------------------------------------
chk_a <- rbindlist(lapply(seq_len(nrow(U)), function(i) {
  u <- U[i]
  f <- file.path(FM, "results", dirs[[u$set]], u$trait, paste0(u$locus_id, "_mesusie.rds"))
  o <- readRDS(f)
  sid <- o$snp_ids; pip <- as.numeric(o$pip)
  pos <- as.integer(vapply(strsplit(sid, ":", fixed = TRUE), `[`, character(1), 2L))
  cs <- lapply(o$credible_sets, function(ix) if (is.numeric(ix)) as.integer(ix) else match(ix, sid))
  mx <- max(pip); ties <- sum(pip >= mx - 1e-12)
  lead <- which.max(pip)
  hit <- which(vapply(cs, function(ix) lead %in% ix, logical(1)))
  size <- if (length(hit)) min(vapply(cs[hit], length, integer(1))) else NA_integer_
  data.table(uid = u$uid, set = u$set, trait = u$trait, locus_id = u$locus_id, n_ties_at_max_pip = ties,
             lead_raw = pos[lead], lead_table = u$ma_lead, cs_size_raw = size, cs_size_table = u$ma_cs_size,
             lead_match = pos[lead] == u$ma_lead, size_match = identical(as.integer(size), as.integer(u$ma_cs_size)))
}), fill = TRUE)
chk_a[n_ties_at_max_pip > 1L, lead_match := NA]
fwrite(chk_a, file.path(out_dir, "verify_A_ma_side.tsv"), sep = "\t")
cat(sprintf("A: %d units; lead mismatches %d; size mismatches %d; units with ties at max PIP %d\n",
            nrow(chk_a), sum(chk_a$lead_match == FALSE, na.rm = TRUE), sum(!chk_a$size_match), sum(chk_a$n_ties_at_max_pip > 1L)))

# ---- B. EUR-only refit from the raw inputs ----------------------------------------
reg_lines <- readLines(file.path(FM, "config/gwas_registry.tsv"))
reg_lines <- reg_lines[!grepl("^\\s*#", reg_lines) & nchar(trimws(reg_lines)) > 0]
REG <- fread(text = paste(reg_lines, collapse = "\n"), sep = "\t")

cand <- U[ma_ok == TRUE & A_ok == TRUE]
cand <- cand[order(set, trait)]
pick <- cand[round(seq(1L, nrow(cand), length.out = min(n_refit, nrow(cand))))]
if (length(args) >= 5L) pick <- cand[uid %in% as.integer(strsplit(args[[5]], ",", fixed = TRUE)[[1]])]   # explicit unit ids
res <- rbindlist(lapply(seq_len(nrow(pick)), function(i) {
  u <- pick[i]; ue <- UE[uid == u$uid]
  # sumstats in the window
  p <- file.path(FM, REG[study_name == u$A, sumstats_path][1])
  ss <- fread(p, select = c("chromosome", "position", "allele1", "allele2", "beta", "se"))
  setnames(ss, c("chr", "pos", "a1", "a2", "beta", "se"))
  ss <- ss[chr == u$chr & pos >= u$win_start & pos <= u$win_end & !is.na(beta) & !is.na(se) & se > 0]
  # whole LD blocks
  stems <- strsplit(ue$ld_key, "|", fixed = TRUE)[[1]]
  bims <- lapply(stems, function(s) { b <- fread(paste0(s, ".bim")); setnames(b, c("chr", "rsid", "dk", "pos", "alt", "ref")); b })
  mats <- lapply(stems, function(s) as.matrix(fread(paste0(s, ".ld"))))
  bim <- rbindlist(bims); M <- if (length(mats) == 1L) mats[[1]] else as.matrix(bdiag(mats))
  bim[, key := paste(chr, pos, alt, ref, sep = ":")]
  # variants used by the multi-ancestry fit
  o <- readRDS(file.path(FM, "results", dirs[[u$set]], u$trait, paste0(u$locus_id, "_mesusie.rds")))
  vpos <- as.integer(vapply(strsplit(o$snp_ids, ":", fixed = TRUE), `[`, character(1), 2L))
  ss <- ss[pos %in% vpos]
  m1 <- merge(ss, bim[, .(pos, alt, ref, row = .I)], by = "pos", allow.cartesian = TRUE)
  fwd <- m1[a1 == alt & a2 == ref]; flp <- m1[a1 == ref & a2 == alt & !(pos %in% fwd$pos)]
  flp[, beta := -beta]
  al <- rbind(fwd, flp)[!duplicated(pos)]
  setorder(al, row)
  R <- M[al$row, al$row]; R <- R + 1e-3 * diag(nrow(R))
  z <- al$beta / al$se
  nA <- as.numeric(REG[study_name == u$A, N_tot])
  fit <- susie_rss(z = z, R = R, n = nA, L = 10L, estimate_residual_variance = FALSE,
                   max_iter = 200, coverage = 0.95, min_abs_corr = 0.5)
  pip <- as.numeric(fit$pip)
  cand_ix <- which(pip >= max(pip) - 1e-12)                       # tie rule from the prespecification: larger |z|
  lead <- if (length(cand_ix) == 1L) cand_ix else cand_ix[which.max(abs(z[cand_ix]))]
  hit <- which(vapply(fit$sets$cs, function(ix) lead %in% ix, logical(1)))
  size <- if (length(hit)) min(vapply(fit$sets$cs[hit], length, integer(1))) else NA_integer_
  data.table(uid = u$uid, set = u$set, trait = u$trait, locus_id = u$locus_id, n_variants_refit = nrow(al), n_ties_refit = length(cand_ix),
             n_variants_table = u$n_A_aligned, converged = isTRUE(fit$converged),
             lead_refit = al$pos[lead], lead_table = u$A_lead, size_refit = size, size_table = u$A_cs_size,
             lead_match = al$pos[lead] == u$A_lead, size_match = identical(as.integer(size), as.integer(u$A_cs_size)))
}), fill = TRUE)
fwrite(res, file.path(out_dir, "verify_B_eur_refit.tsv"), sep = "\t")
print(chk_a[, .(uid, set, trait, lead_match, size_match, n_ties_at_max_pip)])
print(res)

# ---- C. post hoc nested EUR comparators (Addendum D) re-derived from the raw inputs ------------------
# Env H23_VERIFY_NESTED = k: for each (set, trait) in ALT/AST, re-derive the first k units with the lead replaced and the first k
# without it in the Mm arm. The combined z is rebuilt from the two-step description (scale z_B, then add with weight sqrt(n_add)).
k_nest <- as.integer(Sys.getenv("H23_VERIFY_NESTED", unset = "0"))
if (k_nest > 0L) {
  obj <- readRDS(file.path(fits_dir, "unit_fits.rds"))
  stored <- setNames(obj$fits, vapply(obj$fits, function(f) as.character(f$uid), character(1)))
  reg_N <- function(st) as.numeric(REG[study_name == st, N_tot][1])
  Uc <- U[trait %in% c("ALT", "AST") & cmp_MA == TRUE & cmp_Mf == TRUE & cmp_Mm == TRUE]
  pick_c <- rbindlist(lapply(split(Uc, by = c("set", "trait")), function(g) rbind(head(g[rep_Mm == TRUE], k_nest), head(g[rep_Mm == FALSE], k_nest))))
  cat(sprintf("C: re-deriving %d units\n", nrow(pick_c)))

  # noise variance of the Mm combined z relative to a real subsample of the same size, all ALT/AST units with study B
  Ua <- U[trait %in% c("ALT", "AST") & !is.na(B)]
  Ua[, `:=`(nA = vapply(A, reg_N, numeric(1)), nB = vapply(B, reg_N, numeric(1)))]
  Ua[, `:=`(noise_ratio_Mm = (nA + n_noneur^2 / nB) / (nA + n_noneur), nadd_over_nA_Mm = n_noneur / nA, nadd_over_nA_Mf = nB / nA)]
  fwrite(Ua[, .(uid, set, trait, nA, nB, n_noneur, nadd_over_nA_Mm, nadd_over_nA_Mf, noise_ratio_Mm)], file.path(out_dir, "verify_C_noise_ratio.tsv"), sep = "\t")
  print(Ua[, .(n = .N, noise_ratio_min = min(noise_ratio_Mm, na.rm = TRUE), noise_ratio_median = median(noise_ratio_Mm, na.rm = TRUE), noise_ratio_max = max(noise_ratio_Mm, na.rm = TRUE),
               nadd_Mm_min = min(nadd_over_nA_Mm, na.rm = TRUE), nadd_Mm_max = max(nadd_over_nA_Mm, na.rm = TRUE),
               nadd_Mf_min = min(nadd_over_nA_Mf, na.rm = TRUE), nadd_Mf_max = max(nadd_over_nA_Mf, na.rm = TRUE)), by = set])

  # summary statistics inside the windows of the picked units, one read per study
  ss_win <- list()
  for (st in unique(c(pick_c$A, pick_c$B))) {
    ss <- fread(file.path(FM, REG[study_name == st, sumstats_path][1]), select = c("chromosome", "position", "allele1", "allele2", "beta", "se"))
    setnames(ss, c("chr", "pos", "a1", "a2", "beta", "se"))
    ss[, `:=`(chr = as.integer(chr), pos = as.integer(pos), beta = as.numeric(beta), se = as.numeric(se))]
    w <- pick_c[A == st | B == st]
    keep <- rep(FALSE, nrow(ss))
    for (i in seq_len(nrow(w))) keep <- keep | (ss$chr == w$chr[i] & ss$pos >= w$win_start[i] & ss$pos <= w$win_end[i])
    ss_win[[st]] <- ss[keep & !is.na(beta) & !is.na(se) & se > 0]
    rm(ss); invisible(gc())
  }

  chk_c <- rbindlist(lapply(seq_len(nrow(pick_c)), function(i) {
    u <- pick_c[i]; ue <- UE[uid == u$uid]; st0 <- stored[[as.character(u$uid)]]
    o <- readRDS(file.path(FM, "results", dirs[[u$set]], u$trait, paste0(u$locus_id, "_mesusie.rds")))
    vpos <- as.integer(vapply(strsplit(o$snp_ids, ":", fixed = TRUE), `[`, character(1), 2L))
    non_eur <- setdiff(names(o$studies), "EUR")
    n_add_mm <- sum(vapply(non_eur, function(a) reg_N(unname(o$studies[[a]])), numeric(1)))
    stems <- strsplit(ue$ld_key, "|", fixed = TRUE)[[1]]
    bim <- rbindlist(lapply(stems, function(s) { b <- fread(paste0(s, ".bim")); setnames(b, c("chr", "rsid", "dk", "pos", "alt", "ref")); b }))
    mats <- lapply(stems, function(s) as.matrix(fread(paste0(s, ".ld"))))
    M <- if (length(mats) == 1L) mats[[1]] else as.matrix(bdiag(mats))
    z_of <- function(study) {
      ss <- ss_win[[study]][chr == u$chr & pos >= u$win_start & pos <= u$win_end & pos %in% vpos]
      m1 <- merge(ss, bim[, .(pos, alt, ref, row = .I)], by = "pos", allow.cartesian = TRUE)
      fwd <- m1[a1 == alt & a2 == ref]; flp <- m1[a1 == ref & a2 == alt & !(pos %in% fwd$pos)]
      flp[, beta := -beta]
      al <- rbind(fwd, flp)[!duplicated(pos)]
      al[, .(pos, row, z = beta / se)]
    }
    both <- merge(z_of(u$A), z_of(u$B), by = c("pos", "row"), suffixes = c("A", "B"))
    setorder(both, row)
    R <- M[both$row, both$row]; R <- R + 1e-3 * diag(nrow(R))
    nA <- reg_N(u$A); nB <- reg_N(u$B)
    refit <- function(z, n) {
      fit <- susie_rss(z = z, R = R, n = n, L = 10L, estimate_residual_variance = FALSE, max_iter = 200, coverage = 0.95, min_abs_corr = 0.5)
      pip <- as.numeric(fit$pip); cand_ix <- which(pip >= max(pip) - 1e-12)
      lead <- if (length(cand_ix) == 1L) cand_ix else cand_ix[which.max(abs(z[cand_ix]))]
      hit <- which(vapply(fit$sets$cs, function(ix) lead %in% ix, logical(1)))
      cs <- if (length(hit)) fit$sets$cs[[hit[which.min(vapply(fit$sets$cs[hit], length, integer(1)))]]] else integer(0)
      list(lead_pos = both$pos[lead], cs_size = length(cs), cs_pos = both$pos[cs])
    }
    nest_z <- function(n_add) { zBs <- both$zB * sqrt(n_add / nB); (sqrt(nA) * both$zA + sqrt(n_add) * zBs) / sqrt(nA + n_add) }
    set.seed(20260923L)
    ref <- refit(both$zA, nA)
    mf <- refit(nest_z(nB), nA + nB)
    mm <- refit(nest_z(n_add_mm), nA + n_add_mm)
    data.table(uid = u$uid, set = u$set, trait = u$trait, locus_id = u$locus_id, n_var_refit = nrow(both), n_VB_stored = st0$n_VB,
               n_add_mm_raw = n_add_mm, n_add_mm_table = u$n_noneur,
               ref_lead_match = ref$lead_pos == st0$A_b$lead_pos, ref_size_match = ref$cs_size == st0$A_b$cs_size,
               mf_lead_match = mf$lead_pos == st0$Mf$lead_pos, mf_size_match = mf$cs_size == st0$Mf$cs_size,
               mf_rep_raw = !(mf$lead_pos %in% ref$cs_pos), mf_rep_table = u$rep_Mf,
               mm_lead_match = mm$lead_pos == st0$Mm$lead_pos, mm_size_match = mm$cs_size == st0$Mm$cs_size,
               mm_rep_raw = !(mm$lead_pos %in% ref$cs_pos), mm_rep_table = u$rep_Mm)
  }), fill = TRUE)
  fwrite(chk_c, file.path(out_dir, "verify_C_nested_refit.tsv"), sep = "\t")
  print(chk_c)
  cat(sprintf("C: %d units; mismatches in ref lead/size %d/%d, Mf lead/size/rep %d/%d/%d, Mm lead/size/rep %d/%d/%d; variant count mismatches %d\n",
              nrow(chk_c), sum(!chk_c$ref_lead_match), sum(!chk_c$ref_size_match),
              sum(!chk_c$mf_lead_match), sum(!chk_c$mf_size_match), sum(chk_c$mf_rep_raw != chk_c$mf_rep_table),
              sum(!chk_c$mm_lead_match), sum(!chk_c$mm_size_match), sum(chk_c$mm_rep_raw != chk_c$mm_rep_table),
              sum(chk_c$n_var_refit != chk_c$n_VB_stored)))
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
