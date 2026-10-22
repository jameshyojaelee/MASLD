#!/usr/bin/env Rscript
# H23 stage 2: endpoints, tests, controls. Reads the output of h23_fit_units.R.
# Prespecification: Analysis/Hypotheses/H23-20261009T140416Z/PRESPEC_H23.md
# Usage: Rscript h23_analyze.R <fits_dir> <new_output_dir>
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: h23_analyze.R <fits_dir> <new_output_dir>")
fits_dir <- args[[1]]; out_dir <- args[[2]]
if (file.exists(out_dir)) stop("Refusing to overwrite existing output dir: ", out_dir)
dir.create(out_dir, recursive = TRUE)

SEED <- 20260923L
set.seed(SEED)
MIN_UNITS <- as.integer(Sys.getenv("H23_MIN_UNITS", unset = "10"))  # 10 per the prespecification; lowered only for smoke runs
N_BOOT <- 2000L
N_PERM <- 1000L

obj <- readRDS(file.path(fits_dir, "unit_fits.rds"))
units <- obj$units
fits  <- obj$fits
byuid <- setNames(fits, vapply(fits, function(f) as.character(f$uid), character(1)))

# ---- per-unit endpoint table ----------------------------------------------
fld  <- function(f, nm) if (is.null(f) || is.null(f[[nm]])) NA else f[[nm]]
vec  <- function(f, nm) if (is.null(f) || is.null(f[[nm]])) integer(0) else f[[nm]][[1]]
fok  <- function(f) isTRUE(f$ok)

cmp <- function(cmpf, ref) {
  # cmpf = comparator summary, ref = reference summary
  ok <- fok(cmpf) && fok(ref)
  if (!ok) return(list(ok = FALSE, tie = NA, replaced = NA, strict = NA, l2 = NA_real_, gene_ok = NA, gene_changed = NA,
                       gene_ok_all = NA, gene_changed_all = NA, nearest_ok = NA, nearest_changed = NA,
                       lead_c = NA_integer_, lead_r = NA_integer_, size_c = NA_integer_, size_r = NA_integer_))
  gc <- vec(cmpf, "genes"); gr <- vec(ref, "genes")                 # primary: liver eQTL pairs with PVAL < 1e-5
  ga <- vec(cmpf, "genes_all"); gb <- vec(ref, "genes_all")         # sensitivity: every tested pair (prespecified original map)
  list(ok = TRUE,
       tie = isTRUE(cmpf$n_ties > 1L) || isTRUE(ref$n_ties > 1L),
       replaced = !(cmpf$lead_pos %in% vec(ref, "cs_pos")),
       strict = cmpf$lead_pos != ref$lead_pos,
       l2 = log2(cmpf$cs_size / ref$cs_size),
       gene_ok = length(gc) > 0L && length(gr) > 0L,
       gene_changed = if (length(gc) > 0L && length(gr) > 0L) length(intersect(gc, gr)) == 0L else NA,
       nearest_ok = !is.na(cmpf$nearest_gene) && !is.na(ref$nearest_gene),
       nearest_changed = if (!is.na(cmpf$nearest_gene) && !is.na(ref$nearest_gene)) cmpf$nearest_gene != ref$nearest_gene else NA,
       gene_ok_all = length(ga) > 0L && length(gb) > 0L,
       gene_changed_all = if (length(ga) > 0L && length(gb) > 0L) length(intersect(ga, gb)) == 0L else NA,
       lead_c = cmpf$lead_pos, lead_r = ref$lead_pos, size_c = cmpf$cs_size, size_r = ref$cs_size)
}

rows <- lapply(seq_len(nrow(units)), function(i) {
  u <- units[i]; f <- byuid[[as.character(u$uid)]]
  if (is.null(f) || !is.null(f$error)) {
    return(data.table(uid = u$uid, fit_error = if (is.null(f)) "no result" else f$error))
  }
  m <- cmp(f$ma, f$A_ma)
  b <- if (!is.null(f$B) && !is.null(f$A_b)) cmp(f$B, f$A_b) else cmp(NULL, NULL)
  mf <- if (!is.null(f$Mf) && !is.null(f$A_b)) cmp(f$Mf, f$A_b) else cmp(NULL, NULL)   # post hoc nested EUR, all of B added
  mm <- if (!is.null(f$Mm) && !is.null(f$A_b)) cmp(f$Mm, f$A_b) else cmp(NULL, NULL)   # post hoc nested EUR, B added at the non-EUR sample size
  data.table(uid = u$uid, fit_error = NA_character_,
             n_VMA = f$n_VMA, n_A_aligned = f$n_A_aligned, n_VB = f$n_VB,
             ma_ok = fok(f$ma), A_ok = fok(f$A_ma), ma_n_cs = fld(f$ma, "n_cs"), A_n_cs = fld(f$A_ma, "n_cs"),
             A_converged = isTRUE(f$A_ma$converged),
             ma_lead = fld(f$ma, "lead_pos"), A_lead = fld(f$A_ma, "lead_pos"),
             ma_cs_size = fld(f$ma, "cs_size"), A_cs_size = fld(f$A_ma, "cs_size"),
             A_cs_pos = list(vec(f$A_ma, "cs_pos")), ma_genes = paste(vec(f$ma, "genes"), collapse = ";"),
             A_genes = paste(vec(f$A_ma, "genes"), collapse = ";"),
             cmp_MA = m$ok, tie_MA = m$tie, rep_MA = m$replaced, strict_MA = m$strict, l2_MA = m$l2,
             gene_ok_MA = m$gene_ok, gene_ch_MA = m$gene_changed, gene_ok_all_MA = m$gene_ok_all, gene_ch_all_MA = m$gene_changed_all, nearest_ok_MA = m$nearest_ok, nearest_ch_MA = m$nearest_changed,
             cmp_B = b$ok, tie_B = b$tie, rep_B = b$replaced, strict_B = b$strict, l2_B = b$l2,
             gene_ok_B = b$gene_ok, gene_ch_B = b$gene_changed, gene_ok_all_B = b$gene_ok_all, gene_ch_all_B = b$gene_changed_all, nearest_ok_B = b$nearest_ok, nearest_ch_B = b$nearest_changed,
             B_lead = if (b$ok) b$lead_c else NA_integer_,
             cmp_Mf = mf$ok, tie_Mf = mf$tie, rep_Mf = mf$replaced, l2_Mf = mf$l2, gene_ok_Mf = mf$gene_ok, gene_ch_Mf = mf$gene_changed,
             nearest_ok_Mf = mf$nearest_ok, nearest_ch_Mf = mf$nearest_changed,
             cmp_Mm = mm$ok, tie_Mm = mm$tie, rep_Mm = mm$replaced, l2_Mm = mm$l2, gene_ok_Mm = mm$gene_ok, gene_ch_Mm = mm$gene_changed,
             nearest_ok_Mm = mm$nearest_ok, nearest_ch_Mm = mm$nearest_changed)
})
U <- merge(units[, .(uid, set, trait, locus_id, chr, win_start, win_end, A, B, n_noneur, ma_ok_enum = ma_ok)],
           rbindlist(rows, fill = TRUE), by = "uid", all.x = TRUE, sort = FALSE)
for (col in c("cmp_MA", "cmp_B", "cmp_Mf", "cmp_Mm", "ma_ok", "A_ok")) U[is.na(get(col)), (col) := FALSE]
U[, paired := cmp_MA & cmp_B]
# descriptive only: a converged EUR-only fit with zero credible sets is kept apart from a missing or non-converged EUR-only fit
U[, ma_cs_no_A_cs := ma_ok & A_converged & !is.na(A_n_cs) & A_n_cs == 0L]
U[, ma_cs_A_unusable := ma_ok & !A_converged]
U[, A_cs_no_ma_cs := A_ok & !ma_ok]

# ---- inventory ------------------------------------------------------------
inv <- U[, .(n_units = .N, n_fit_error = sum(!is.na(fit_error)), n_ma_ok = sum(ma_ok), n_A_ok = sum(A_ok),
             n_comparable_MA = sum(cmp_MA), n_comparable_B = sum(cmp_B), n_paired = sum(paired),
             n_comparable_Mm = sum(cmp_Mm), n_comparable_Mf = sum(cmp_Mf),
             n_ma_cs_but_no_eur_cs = sum(ma_cs_no_A_cs, na.rm = TRUE),
             n_ma_cs_eur_fit_unusable = sum(ma_cs_A_unusable, na.rm = TRUE)), by = .(set, trait)]
inv <- merge(obj$inventory_fit[, .(set, trait, n_rds, n_with_eur_arm, n_after_window_collapse, note)], inv,
             by = c("set", "trait"), all.x = TRUE, sort = FALSE)
fwrite(inv, file.path(out_dir, "inventory.tsv"), sep = "\t")

# ---- tests ----------------------------------------------------------------
boot_ci <- function(x, f, B = N_BOOT) {
  n <- length(x); if (n < 2L) return(c(NA_real_, NA_real_))
  v <- replicate(B, f(x[sample.int(n, n, replace = TRUE)]))
  as.numeric(quantile(v, c(0.025, 0.975), na.rm = TRUE))
}
boot_ci_idx <- function(n, f, B = N_BOOT) {
  if (n < 2L) return(c(NA_real_, NA_real_))
  v <- replicate(B, f(sample.int(n, n, replace = TRUE)))
  as.numeric(quantile(v, c(0.025, 0.975), na.rm = TRUE))
}

binary_test_primary <- function(d, rep_MA, rep_B, ok_MA, ok_B, label) {
  # paired by locus: units where both MA and baseline comparisons are evaluable for this endpoint
  p <- d[[ok_MA]] & d[[ok_B]]
  p[is.na(p)] <- FALSE
  n <- sum(p)
  out <- list(n_paired = n, n_MA = sum(d[[ok_MA]], na.rm = TRUE), n_B = sum(d[[ok_B]], na.rm = TRUE))
  if (n < MIN_UNITS) return(c(out, list(status = sprintf("untestable: %d paired units < %d", n, MIN_UNITS))))
  x <- d[p]
  a <- as.logical(x[[rep_MA]]); b <- as.logical(x[[rep_B]])
  nb <- sum(a & !b); nc <- sum(!a & b)
  pv <- if (nb + nc == 0L) 1 else binom.test(nb, nb + nc, 0.5, alternative = "greater")$p.value
  ci <- boot_ci_idx(n, function(ix) mean(a[ix]) - mean(b[ix]))
  # unpaired sensitivity: all evaluable units
  A <- as.logical(d[[rep_MA]][d[[ok_MA]] %in% TRUE]); B <- as.logical(d[[rep_B]][d[[ok_B]] %in% TRUE])
  pf <- if (length(A) >= MIN_UNITS && length(B) >= MIN_UNITS)
    fisher.test(matrix(c(sum(A), sum(!A), sum(B), sum(!B)), 2, byrow = TRUE), alternative = "greater")$p.value else NA_real_
  c(out, list(status = "tested", rate_MA = mean(a), rate_B = mean(b), effect = mean(a) - mean(b),
              ci_lo = ci[1], ci_hi = ci[2], disc_MA_only = nb, disc_B_only = nc, p = pv, p_unpaired = pf,
              rate_MA_all = mean(A), rate_B_all = mean(B)))
}

size_test_primary <- function(d) {
  p <- d$paired
  n <- sum(p)
  out <- list(n_paired = n, n_MA = sum(d$cmp_MA), n_B = sum(d$cmp_B))
  if (n < MIN_UNITS) return(c(out, list(status = sprintf("untestable: %d paired units < %d", n, MIN_UNITS))))
  x <- d[p == TRUE]
  diffs <- x$l2_MA - x$l2_B
  pv <- if (all(diffs == 0)) 1 else suppressWarnings(wilcox.test(diffs, mu = 0, alternative = "less", exact = FALSE)$p.value)
  ci <- boot_ci_idx(length(diffs), function(ix) median(diffs[ix]))
  pu <- if (sum(d$cmp_MA) >= MIN_UNITS && sum(d$cmp_B) >= MIN_UNITS)
    suppressWarnings(wilcox.test(d$l2_MA[d$cmp_MA], d$l2_B[d$cmp_B], alternative = "less", exact = FALSE)$p.value) else NA_real_
  c(out, list(status = "tested", rate_MA = median(x$l2_MA), rate_B = median(x$l2_B), effect = median(diffs),
              ci_lo = ci[1], ci_hi = ci[2], disc_MA_only = NA_integer_, disc_B_only = NA_integer_,
              p = pv, p_unpaired = pu, rate_MA_all = median(d$l2_MA[d$cmp_MA]), rate_B_all = median(d$l2_B[d$cmp_B])))
}

# secondary family: no own baseline; pooled ALT+AST baseline within the same set, unpaired
pooled_test <- function(d, base, endpoint) {
  if (endpoint == "lead_replaced") {
    a <- d$rep_MA[d$cmp_MA]; b <- base$rep_B[base$cmp_B]
    out <- list(n_paired = NA_integer_, n_MA = length(a), n_B = length(b))
    if (length(a) < MIN_UNITS || length(b) < MIN_UNITS) return(c(out, list(status = "untestable: units < 10")))
    pv <- fisher.test(matrix(c(sum(a), sum(!a), sum(b), sum(!b)), 2, byrow = TRUE), alternative = "greater")$p.value
    ci <- boot_ci_idx(length(a), function(ix) mean(a[ix]) - mean(b))
    return(c(out, list(status = "tested", rate_MA = mean(a), rate_B = mean(b), effect = mean(a) - mean(b), ci_lo = ci[1], ci_hi = ci[2],
                       disc_MA_only = NA_integer_, disc_B_only = NA_integer_, p = pv, p_unpaired = NA_real_,
                       rate_MA_all = mean(a), rate_B_all = mean(b))))
  }
  if (endpoint == "gene_changed") {
    a <- d$gene_ch_MA[d$cmp_MA & d$gene_ok_MA %in% TRUE]; b <- base$gene_ch_B[base$cmp_B & base$gene_ok_B %in% TRUE]
    out <- list(n_paired = NA_integer_, n_MA = length(a), n_B = length(b))
    if (length(a) < MIN_UNITS || length(b) < MIN_UNITS) return(c(out, list(status = "untestable: units < 10")))
    pv <- fisher.test(matrix(c(sum(a), sum(!a), sum(b), sum(!b)), 2, byrow = TRUE), alternative = "greater")$p.value
    ci <- boot_ci_idx(length(a), function(ix) mean(a[ix]) - mean(b))
    return(c(out, list(status = "tested", rate_MA = mean(a), rate_B = mean(b), effect = mean(a) - mean(b), ci_lo = ci[1], ci_hi = ci[2],
                       disc_MA_only = NA_integer_, disc_B_only = NA_integer_, p = pv, p_unpaired = NA_real_,
                       rate_MA_all = mean(a), rate_B_all = mean(b))))
  }
  a <- d$l2_MA[d$cmp_MA]; b <- base$l2_B[base$cmp_B]
  out <- list(n_paired = NA_integer_, n_MA = length(a), n_B = length(b))
  if (length(a) < MIN_UNITS || length(b) < MIN_UNITS) return(c(out, list(status = "untestable: units < 10")))
  pv <- suppressWarnings(wilcox.test(a, b, alternative = "less", exact = FALSE)$p.value)
  ci <- boot_ci_idx(length(a), function(ix) median(a[ix]) - median(b))
  c(out, list(status = "tested", rate_MA = median(a), rate_B = median(b), effect = median(a) - median(b), ci_lo = ci[1], ci_hi = ci[2],
              disc_MA_only = NA_integer_, disc_B_only = NA_integer_, p = pv, p_unpaired = NA_real_,
              rate_MA_all = median(a), rate_B_all = median(b)))
}

ENDPOINTS <- c("lead_replaced", "cs_size_log2", "gene_changed")
res <- list()
for (sn in unique(U$set)) {
  for (tr in unique(U[set == sn, trait])) {
    d <- U[set == sn & trait == tr]
    fam <- if (tr %in% c("ALT", "AST")) "primary" else "secondary"
    base <- U[set == sn & trait %in% c("ALT", "AST")]
    for (ep in ENDPOINTS) {
      r <- if (fam == "primary") {
        switch(ep,
          lead_replaced = binary_test_primary(d, "rep_MA", "rep_B", "cmp_MA", "cmp_B"),
          cs_size_log2  = size_test_primary(d),
          gene_changed  = { d2 <- copy(d); d2[, ok_MA := cmp_MA & gene_ok_MA %in% TRUE]; d2[, ok_B := cmp_B & gene_ok_B %in% TRUE]
                            binary_test_primary(d2, "gene_ch_MA", "gene_ch_B", "ok_MA", "ok_B") })
      } else pooled_test(d, base, ep)
      res[[length(res) + 1L]] <- c(list(family = fam, set = sn, trait = tr, endpoint = ep), r)
    }
  }
}
tests <- rbindlist(res, fill = TRUE)
for (cn in c("n_paired", "rate_MA", "rate_B", "effect", "p", "q_BH", "ci_lo", "ci_hi")) if (!cn %in% names(tests)) tests[, (cn) := NA_real_]
tests[, q_BH := NA_real_]
# BH over the prespecified family sizes (primary 12, secondary 9). Tests that could not run
# count as p = 1, so n = family size is passed to p.adjust.
FAM_SIZE <- c(primary = 12L, secondary = 9L)
for (fm in c("primary", "secondary")) {
  idx <- which(tests$family == fm & tests$status == "tested" & !is.na(tests$p))
  if (length(idx)) tests[idx, q_BH := p.adjust(p, method = "BH", n = FAM_SIZE[[fm]])]
  tests[family == fm, `:=`(family_size_tested = length(idx), family_size_prespec = FAM_SIZE[[fm]])]
}
# Sensitivity, outside both BH families (added after the smoke run, before the full run): lead_replaced for ALT and AST
# with units dropped when either fit in the comparison has several variants tied at the maximum PIP.
sens <- list()
for (sn in unique(U$set)) for (tr in intersect(c("ALT", "AST"), unique(U[set == sn, trait]))) {
  d2 <- copy(U[set == sn & trait == tr])
  d2[, ok_MA := cmp_MA & !(tie_MA %in% TRUE)]; d2[, ok_B := cmp_B & !(tie_B %in% TRUE)]
  r <- binary_test_primary(d2, "rep_MA", "rep_B", "ok_MA", "ok_B")
  sens[[length(sens) + 1L]] <- c(list(family = "sensitivity_noties", set = sn, trait = tr, endpoint = "lead_replaced"), r)
  # gene endpoint with every tested variant-gene pair as the map (the original prespecified map)
  d3 <- copy(U[set == sn & trait == tr])
  d3[, ok_MA := cmp_MA & gene_ok_all_MA %in% TRUE]; d3[, ok_B := cmp_B & gene_ok_all_B %in% TRUE]
  r3 <- binary_test_primary(d3, "gene_ch_all_MA", "gene_ch_all_B", "ok_MA", "ok_B")
  sens[[length(sens) + 1L]] <- c(list(family = "sensitivity_gene_allpairs", set = sn, trait = tr, endpoint = "gene_changed"), r3)
  # nearest expressed gene (estimated TSS) to the lead variant
  d4 <- copy(U[set == sn & trait == tr])
  d4[, ok_MA := cmp_MA & nearest_ok_MA %in% TRUE]; d4[, ok_B := cmp_B & nearest_ok_B %in% TRUE]
  r4 <- binary_test_primary(d4, "nearest_ch_MA", "nearest_ch_B", "ok_MA", "ok_B")
  sens[[length(sens) + 1L]] <- c(list(family = "sensitivity_nearest_gene", set = sn, trait = tr, endpoint = "nearest_gene_changed"), r4)
}
# Post hoc, exploratory (Addendum D, written after the primary results were read): the same paired tests with a nested EUR comparator in
# place of the independent EUR study. Mm adds study B at the non-EUR sample size, Mf adds all of B. Nominal p, outside both BH families.
for (sn in unique(U$set)) for (tr in intersect(c("ALT", "AST"), unique(U[set == sn, trait]))) for (arm in c("Mm", "Mf")) {
  d <- copy(U[set == sn & trait == tr])
  for (nm in c("cmp", "rep", "l2", "gene_ok", "gene_ch", "nearest_ok", "nearest_ch")) d[, (paste0(nm, "_B")) := get(paste0(nm, "_", arm))]
  d[is.na(cmp_B), cmp_B := FALSE]; d[, paired := cmp_MA & cmp_B]
  d2 <- copy(d); d2[, ok_MA := cmp_MA & gene_ok_MA %in% TRUE]; d2[, ok_B := cmp_B & gene_ok_B %in% TRUE]
  d4 <- copy(d); d4[, ok_MA := cmp_MA & nearest_ok_MA %in% TRUE]; d4[, ok_B := cmp_B & nearest_ok_B %in% TRUE]
  fam <- paste0("exploratory_nested_EUR_", arm)
  sens[[length(sens) + 1L]] <- c(list(family = fam, set = sn, trait = tr, endpoint = "lead_replaced"), binary_test_primary(d, "rep_MA", "rep_B", "cmp_MA", "cmp_B"))
  sens[[length(sens) + 1L]] <- c(list(family = fam, set = sn, trait = tr, endpoint = "cs_size_log2"), size_test_primary(d))
  sens[[length(sens) + 1L]] <- c(list(family = fam, set = sn, trait = tr, endpoint = "gene_changed"), binary_test_primary(d2, "gene_ch_MA", "gene_ch_B", "ok_MA", "ok_B"))
  sens[[length(sens) + 1L]] <- c(list(family = fam, set = sn, trait = tr, endpoint = "nearest_gene_changed"), binary_test_primary(d4, "nearest_ch_MA", "nearest_ch_B", "ok_MA", "ok_B"))
}
if (length(sens)) tests <- rbindlist(list(tests, rbindlist(sens, fill = TRUE)), fill = TRUE)
skipped <- obj$inventory_fit[nchar(note) > 0, .(family = "stop_rule", set, trait, endpoint = "all", status = note)]
if (nrow(skipped)) tests <- rbindlist(list(tests, skipped), fill = TRUE)
fwrite(tests, file.path(out_dir, "tests.tsv"), sep = "\t")

# strict-identity sensitivity (not tested, rates only)
strict <- U[, .(n_cmp_MA = sum(cmp_MA), tie_frac_MA = mean(tie_MA[cmp_MA]), tie_frac_B = mean(tie_B[cmp_B]),
                strict_rate_MA = mean(strict_MA[cmp_MA]), replaced_rate_MA = mean(rep_MA[cmp_MA]),
                n_cmp_B = sum(cmp_B), strict_rate_B = mean(strict_B[cmp_B]), replaced_rate_B = mean(rep_B[cmp_B]),
                median_log2_cs_ratio_MA = median(l2_MA[cmp_MA]), median_log2_cs_ratio_B = median(l2_B[cmp_B])), by = .(set, trait)]
fwrite(strict, file.path(out_dir, "rates_by_set_trait.tsv"), sep = "\t")

# ---- controls -------------------------------------------------------------
ctl <- list()
pc_sites <- data.table(gene = c("PNPLA3", "TM6SF2", "GCKR", "MBOAT7"), chr = c(22L, 19L, 2L, 19L),
                       pos = c(44324727L, 19379549L, 27730940L, 54676763L))
pcu <- rbindlist(lapply(seq_len(nrow(pc_sites)), function(i) {
  s <- pc_sites[i]; x <- U[chr == s$chr & win_start <= s$pos & win_end >= s$pos]
  if (nrow(x) == 0L) return(NULL)
  x[, .(control = "PC_stable_lead", site = s$gene, set, trait, locus_id, comparable = cmp_MA, lead_replaced = rep_MA,
        ma_lead, A_lead)]
}))
n_pc <- if (is.null(pcu)) 0L else sum(pcu$comparable)
stable <- if (n_pc > 0L) mean(!pcu$lead_replaced[pcu$comparable]) else NA_real_
ctl$PC <- list(n_comparable_units = n_pc, stable_fraction = stable,
               pass = !is.na(stable) && n_pc >= 4L && stable >= 0.75)
if (!is.null(pcu)) fwrite(pcu, file.path(out_dir, "control_PC_units.tsv"), sep = "\t")

# ALDH2 expected-change control. The liver eQTL table tests every variant within +-1 Mb of each gene's TSS, so the
# midpoint of the ALDH2 tested range estimates the TSS. The control window is [TSS, TSS + 50 kb] (gene body, approximate).
ec <- NULL
rng <- tryCatch(fread(file.path(fits_dir, "aldh2_cis_range.tsv")), error = function(e) NULL)
if (!is.null(rng) && nrow(rng) == 1L) {
  tss <- (rng$cis_min + rng$cis_max) / 2; lo <- tss; hi <- tss + 50000; ch <- rng$chr
  # evaluable = EUR-only fit converged with no credible set, or comparable; a missing or non-converged EUR-only fit is not evaluable
  ec <- U[set == "4WAY" & chr == ch & win_start <= hi & win_end >= lo,
          .(control = "EC_ALDH2", set, trait, locus_id, ma_ok, A_ok, A_converged, A_n_cs, comparable = cmp_MA, lead_replaced = rep_MA,
            eur_cs_absent = ma_ok & A_converged & !is.na(A_n_cs) & A_n_cs == 0L)]
  ec[, evaluable := eur_cs_absent | comparable]
  ec[, eur_cs_absent_or_lead_changed := eur_cs_absent | (comparable & lead_replaced %in% TRUE)]
}
ctl$EC <- if (is.null(ec) || nrow(ec) == 0L) list(status = "not testable: no 4WAY unit overlaps the ALDH2 gene window", pass = NA) else
  if (!any(ec$evaluable)) list(status = "not testable: no evaluable 4WAY unit at the ALDH2 window", n_units = nrow(ec), pass = NA) else
  list(n_units = nrow(ec), n_evaluable = sum(ec$evaluable), n_pass = sum(ec$eur_cs_absent_or_lead_changed[ec$evaluable]),
       pass = all(ec$eur_cs_absent_or_lead_changed[ec$evaluable]))
if (!is.null(ec)) fwrite(ec, file.path(out_dir, "control_EC_units.tsv"), sep = "\t")

# negative control (a): MA lead of unit i against the reference credible set of a different unit j (same set and trait)
# (a plain function over split() copies: the same code evaluated inside U[, j, by] crashed R on the large Platelet group)
nca_one <- function(g) {
  n <- nrow(g)
  if (n < 2L) return(NULL)
  lead <- g$ma_lead; cs <- g$A_cs_pos
  r <- vapply(seq_len(N_PERM), function(b) {
    p <- sample.int(n); tries <- 0L
    while (any(p == seq_len(n)) && tries < 100L) { p <- sample.int(n); tries <- tries + 1L }
    mean(vapply(seq_len(n), function(i) !(lead[p[i]] %in% cs[[i]]), logical(1)))
  }, numeric(1))
  data.table(set = g$set[1], trait = g$trait[1], n_units = n, mean_replaced_under_shuffle = mean(r))
}
nca <- rbindlist(lapply(split(U[cmp_MA == TRUE], by = c("set", "trait")), nca_one))
ctl$NCa <- list(mean_replaced = if (nrow(nca)) weighted.mean(nca$mean_replaced_under_shuffle, nca$n_units) else NA_real_,
                pass = nrow(nca) > 0 && weighted.mean(nca$mean_replaced_under_shuffle, nca$n_units) >= 0.95)
fwrite(nca, file.path(out_dir, "control_NCa_locus_shuffle.tsv"), sep = "\t")

# negative control (b): cross-trait pairing at overlapping windows in 4WAY
x4 <- U[set == "4WAY" & cmp_MA == TRUE]
cross <- rbindlist(lapply(seq_len(nrow(x4)), function(i) {
  o <- x4[chr == x4$chr[i] & trait != x4$trait[i] & win_start <= x4$win_end[i] & win_end >= x4$win_start[i]]
  if (nrow(o) == 0L) return(NULL)
  data.table(ref_trait = x4$trait[i], ma_trait = o$trait, same_trait_replaced = x4$rep_MA[i],
             cross_replaced = !(o$ma_lead %in% x4$A_cs_pos[[i]]))
}))
ctl$NCb <- if (is.null(cross) || nrow(cross) == 0L) list(n_pairs = 0L) else
  list(n_pairs = nrow(cross), cross_trait_replaced_rate = mean(cross$cross_replaced),
       same_trait_replaced_rate_same_loci = mean(cross$same_trait_replaced))
if (!is.null(cross)) fwrite(cross, file.path(out_dir, "control_NCb_cross_trait.tsv"), sep = "\t")

# ---- outputs --------------------------------------------------------------
U_out <- copy(U)[, A_cs_pos := NULL]
fwrite(U_out, file.path(out_dir, "units.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(ctl), function(k) data.table(control = k, json = jsonlite::toJSON(ctl[[k]], auto_unbox = TRUE, null = "null", na = "string"))),
                 fill = TRUE), file.path(out_dir, "controls.tsv"), sep = "\t")
summ <- list(seed = SEED, smoke = isTRUE(obj$smoke), n_units = nrow(U), n_paired_total = sum(U$paired),
             controls = ctl,
             tests_tested = tests[status == "tested", .(family, set, trait, endpoint, n_paired, effect, p, q_BH)])
writeLines(jsonlite::toJSON(summ, auto_unbox = TRUE, pretty = TRUE, na = "string", digits = 6), file.path(out_dir, "summary.json"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(inv); print(tests[, .(family, set, trait, endpoint, status, n_paired, rate_MA, rate_B, p, q_BH)]); print(ctl)
