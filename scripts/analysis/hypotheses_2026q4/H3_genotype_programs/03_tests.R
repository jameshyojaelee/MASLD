#!/usr/bin/env Rscript
# H3 / H31 step 3: person-level genotype-program tests. Rule: PRESPEC_H3.md (including Amendments 1 to 3).
# Unit: person. Model per cohort: score ~ dosage + factor(stage) + sex, HC3 errors; cohorts pooled by
# DerSimonian-Laird random effects; Holm over 5 variants x 2 programs. No package beyond data.table.
# Usage: 03_tests.R CALLS_DIR SCORES_DIR GENO_DIR OUT_DIR
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
calls_dir <- args[[1]]; scores_dir <- args[[2]]; geno_dir <- args[[3]]; out <- args[[4]]
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
SEED <- 20260923L
B <- 2000L
COHORTS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066", "GSE213621", "GSE240729")
PROG <- c(lipid = "score_lipid", ecm = "score_ecm")
PNAMES <- names(PROG); YCOLS <- unname(PROG)
VARS <- data.table(
  rsid = c("rs72613567", "rs738409", "rs58542926", "rs2642438", "rs1260326"),
  gene = c("HSD17B13", "PNPLA3", "TM6SF2", "MTARC1", "GCKR"),
  tested = c("TA carrier (0/1)", "G count (0-2)", "T count (0-2)", "G count (0-2)", "C count (0-2)"),
  literature = c("TA protective [recalled]", "G risk [recalled]", "T risk [recalled]",
                 "reference A protective [recalled]", "reference T risk [recalled]"))

# ---- inputs ----------------------------------------------------------------------------------
ph <- fread(file.path(scores_dir, "person_phenotypes.tsv"))
gt <- fread(file.path(calls_dir, "person_genotypes.tsv"))
gt[, cohort := NULL]
hs <- fread(file.path(geno_dir, "person_variant_counts.tsv"))[rsid == "rs72613567", .(person_id, depth_hsd = depth)]
dat <- merge(merge(ph, gt, by = "person_id"), hs, by = "person_id", all.x = TRUE)
stopifnot(nrow(dat) == nrow(ph), !anyDuplicated(dat$person_id), all(dat$cohort %in% COHORTS))
dat[, stage := as.numeric(stage)]
dat[, sex2 := fifelse(sex %in% c("F", "M"), sex, NA_character_)]

dose_of <- function(call, rsid) {
  if (rsid == "rs72613567") return(ifelse(call == "CARRIER", 1, ifelse(call == "NC", 0, NA_real_)))
  unname(c(RR = 0, RA = 1, AA = 2)[call])
}
for (rs in VARS$rsid) set(dat, j = paste0("dose_", rs), value = dose_of(dat[[paste0(rs, "_call")]], rs))
dat[, dose_pi := fifelse(pi_status == "MZ", 1, fifelse(pi_status == "MM", 0, NA_real_))]

# ---- statistics ------------------------------------------------------------------------------
# OLS with HC3 errors for coefficient j; Y is a matrix (one column per outcome)
hc3 <- function(X, Y, j) {
  if (qr(X)$rank < ncol(X)) return(NULL)
  XtXi <- solve(crossprod(X))
  Bm <- XtXi %*% crossprod(X, Y)
  E <- Y - X %*% Bm
  H <- rowSums((X %*% XtXi) * X)
  if (max(H) > 1 - 1e-8) return(NULL)
  W <- (E / (1 - H))^2
  se <- vapply(seq_len(ncol(Y)), function(k) sqrt((XtXi %*% crossprod(X * W[, k], X) %*% XtXi)[j, j]), numeric(1))
  list(beta = as.numeric(Bm[j, ]), se = se)
}

dl_pool <- function(b, s) {
  k <- length(b); w <- 1 / s^2
  bfe <- sum(w * b) / sum(w); sefe <- sqrt(1 / sum(w))
  Q <- sum(w * (b - bfe)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (s^2 + tau2)
  bre <- sum(wr * b) / sum(wr); sere <- sqrt(1 / sum(wr)); z <- bre / sere
  list(k = k, beta = bre, se = sere, z = z, p = 2 * pnorm(-abs(z)), tau2 = tau2,
       I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0, Q = Q,
       beta_fe = bfe, se_fe = sefe, p_fe = 2 * pnorm(-abs(bfe / sefe)))
}

# stage levels with fewer than 3 persons are merged into the nearest level (ties to the lower)
merge_stage <- function(stage, min_n = 3L) {
  s <- as.numeric(stage)
  repeat {
    tb <- table(s)
    if (length(tb) < 2L || all(tb >= min_n)) break
    lv <- as.numeric(names(tb)); i <- which.min(tb)
    others <- lv[-i]; dd <- abs(others - lv[i])
    s[s == lv[i]] <- others[which(dd == min(dd))[1]]
  }
  factor(s)
}

design0 <- function(d, sf) {
  sexF <- if (uniqueN(d$sex2) == 2) as.numeric(d$sex2 == "F") else NULL
  cbind(model.matrix(~ sf)[, -1, drop = FALSE], sexF)
}

# one entry per cohort: persons, design without the tested columns, outcomes
prep_dose <- function(col, min_with, min_without, depth_min = NA_real_, ycols = YCOLS) {
  lapply(setNames(COHORTS, COHORTS), function(co) {
    d <- dat[cohort == co & !is.na(get(col)) & !is.na(stage) & !is.na(sex2)]
    if (!is.na(depth_min)) d <- d[!is.na(depth_hsd) & depth_hsd >= depth_min]
    dose <- d[[col]]
    res <- list(cohort = co, n = nrow(d), n_with = sum(dose > 0), n_without = sum(dose == 0),
                n_stage_levels = 0L, eligible = FALSE, reason = "")
    if (nrow(d) == 0) { res$reason <- "no persons"; return(res) }
    sf <- merge_stage(d$stage)
    res$n_stage_levels <- nlevels(sf)
    if (res$n_with < min_with) res$reason <- paste0("fewer than ", min_with, " persons with the allele")
    else if (res$n_without < min_without) res$reason <- paste0("fewer than ", min_without, " persons without the allele")
    else if (nlevels(sf) < 2) res$reason <- "fewer than 2 stage levels"
    else res$eligible <- TRUE
    res$X0 <- design0(d, sf)
    res$Y <- as.matrix(d[, ..ycols])
    res$dose <- dose; res$sf <- sf; res$sex2 <- d$sex2; res$person_id <- d$person_id
    res
  })
}

cohort_est <- function(pr, dose_by = NULL, labels = PNAMES, ycols_idx = seq_along(labels)) {
  rows <- list()
  for (co in names(pr)) {
    p <- pr[[co]]
    if (!isTRUE(p$eligible)) next
    dose <- if (is.null(dose_by)) p$dose else dose_by[[co]]
    f <- hc3(cbind(1, dose, p$X0), p$Y[, ycols_idx, drop = FALSE], 2L)
    if (is.null(f)) next
    rows[[co]] <- data.table(cohort = co, program = labels, n = p$n, n_with = p$n_with, n_without = p$n_without,
                             n_stage_levels = p$n_stage_levels, beta = f$beta, se = f$se)
  }
  rbindlist(rows)
}

stop_rule <- function(pr, min_total) {
  el <- Filter(function(p) isTRUE(p$eligible), pr)
  w <- sum(vapply(el, `[[`, numeric(1), "n_with")); wo <- sum(vapply(el, `[[`, numeric(1), "n_without"))
  list(cohorts = paste(names(el), collapse = ","), k = length(el), n_with = w, n_without = wo,
       n_persons = sum(vapply(el, `[[`, numeric(1), "n")), testable = length(el) >= 2 && w >= min_total && wo >= min_total)
}

pool_table <- function(est, labels = PNAMES) {
  rbindlist(lapply(labels, function(pg) {
    e <- est[program == pg]
    if (nrow(e) < 2) return(data.table(program = pg, k = nrow(e)))
    r <- dl_pool(e$beta, e$se)
    lo <- if (nrow(e) >= 3) {
      l <- lapply(seq_len(nrow(e)), function(i) dl_pool(e$beta[-i], e$se[-i]))
      list(min(vapply(l, `[[`, numeric(1), "beta")), max(vapply(l, `[[`, numeric(1), "beta")), max(vapply(l, `[[`, numeric(1), "p")))
    } else list(NA_real_, NA_real_, NA_real_)
    data.table(program = pg, k = r$k, beta = r$beta, se = r$se, z = r$z, p = r$p, tau2 = r$tau2, I2 = r$I2, Q = r$Q,
               beta_fe = r$beta_fe, p_fe = r$p_fe, loco_beta_min = lo[[1]], loco_beta_max = lo[[2]], loco_p_max = lo[[3]])
  }), fill = TRUE)
}

# pooled model over the eligible cohorts (context only): outcome ~ dose + cohort:stage + sex
pooled_context <- function(pr, labels = PNAMES, ycols_idx = seq_along(labels)) {
  el <- Filter(function(p) isTRUE(p$eligible), pr)
  if (length(el) < 1) return(NULL)
  dd <- rbindlist(lapply(el, function(p) {
    data.table(cohort = p$cohort, dose = p$dose, cs = paste(p$cohort, as.character(p$sf), sep = ":"),
               sexF = as.numeric(p$sex2 == "F"), as.data.table(p$Y[, ycols_idx, drop = FALSE]))
  }))
  X <- cbind(1, dose = dd$dose, model.matrix(~ factor(cs), data = dd)[, -1, drop = FALSE], sexF = dd$sexF)
  f <- hc3(X, as.matrix(dd[, (ncol(dd) - length(ycols_idx) + 1):ncol(dd), with = FALSE]), 2L)
  if (is.null(f)) return(NULL)
  data.table(program = labels, ctx_n = nrow(dd), ctx_beta = f$beta, ctx_se = f$se, ctx_p = 2 * pnorm(-abs(f$beta / f$se)))
}

# ---- analysis-set counts for every variant and cohort ----------------------------------------
count_rows <- list()
add_counts <- function(label, pr) {
  count_rows[[label]] <<- rbindlist(lapply(pr, function(p) data.table(
    variant = label, cohort = p$cohort, persons = p$n, n_with = p$n_with, n_without = p$n_without,
    n_stage_levels = p$n_stage_levels, eligible = isTRUE(p$eligible), reason = p$reason)))
}

# ---- score validity (ECM etc. against stage; no genotype) ------------------------------------
sv <- rbindlist(lapply(COHORTS, function(co) rbindlist(lapply(c(PROG, upr = "score_upr"), function(y) {
  d <- dat[cohort == co & !is.na(stage)]
  ct <- suppressWarnings(cor.test(d$stage, d[[y]], method = "spearman", exact = FALSE))
  f <- hc3(cbind(1, d$stage), matrix(d[[y]], ncol = 1), 2L)
  data.table(cohort = co, score = y, n = nrow(d), spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
             beta_per_stage = f$beta, se_hc3 = f$se)
}))))
fwrite(sv, file.path(out, "score_vs_stage_no_genotype.tsv"), sep = "\t")

# ---- H3 primary family: 5 variants x 2 programs -----------------------------------------------
prep5 <- setNames(lapply(VARS$rsid, function(rs) prep_dose(paste0("dose_", rs), 5, 5)), VARS$rsid)
for (i in seq_len(nrow(VARS))) add_counts(VARS$gene[i], prep5[[i]])
stops <- lapply(prep5, stop_rule, min_total = 20)
stop_tab <- rbindlist(lapply(seq_len(nrow(VARS)), function(i) cbind(VARS[i, .(rsid, gene)], as.data.table(stops[[i]]))))
fwrite(stop_tab, file.path(out, "stop_rule_h3.tsv"), sep = "\t")
message("stop rule:"); print(stop_tab)

obs_est <- list(); obs_pool <- list(); obs_ctx <- list()
for (i in seq_len(nrow(VARS))) {
  rs <- VARS$rsid[i]
  if (!stops[[rs]]$testable) next
  est <- cohort_est(prep5[[rs]])
  obs_est[[rs]] <- cbind(rsid = rs, gene = VARS$gene[i], est)
  pt <- pool_table(est); ctx <- pooled_context(prep5[[rs]])
  obs_pool[[rs]] <- cbind(rsid = rs, gene = VARS$gene[i], merge(pt, ctx, by = "program", all.x = TRUE))
}
est_all <- rbindlist(obs_est); pool_all <- rbindlist(obs_pool, fill = TRUE)
fwrite(est_all, file.path(out, "h3_cohort_estimates.tsv"), sep = "\t")

# permutation null: tested allele shuffled within cohort, whole family recomputed
set.seed(SEED)
test_grid <- CJ(rsid = VARS$rsid, program = PNAMES, sorted = FALSE)
test_grid[, key := paste(rsid, program)]
perm_p <- matrix(NA_real_, B, nrow(test_grid), dimnames = list(NULL, test_grid$key))
perm_z <- perm_p
for (b in seq_len(B)) {
  for (rs in VARS$rsid) {
    if (!stops[[rs]]$testable) next
    pr <- prep5[[rs]]
    db <- lapply(pr, function(p) if (isTRUE(p$eligible)) p$dose[sample.int(length(p$dose))] else NULL)
    est <- cohort_est(pr, dose_by = db)
    for (pg in PNAMES) {
      e <- est[program == pg]
      if (nrow(e) < 2) next
      r <- dl_pool(e$beta, e$se)
      perm_p[b, paste(rs, pg)] <- r$p; perm_z[b, paste(rs, pg)] <- r$z
    }
  }
}
pool_all[, key := paste(rsid, program)]
pool_all[, z_obs := z]
pool_all[, p_perm := NA_real_]
for (i in seq_len(nrow(pool_all))) {
  zz <- perm_z[, pool_all$key[i]]
  set(pool_all, i, "p_perm", (1 + sum(abs(zz) >= abs(pool_all$z[i]), na.rm = TRUE)) / (1 + sum(!is.na(zz))))
}
# Holm over the prespecified family of 10; an untestable variant contributes p = 1
fam <- merge(test_grid, pool_all[, .(key, p)], by = "key", all.x = TRUE, sort = FALSE)
fam[is.na(p), p := 1]
fam[, p_holm := p.adjust(p, "holm")]
pool_all <- merge(pool_all, fam[, .(key, p_holm)], by = "key", all.x = TRUE)
obs_min_p <- min(pool_all$p, na.rm = TRUE)
perm_minp <- apply(perm_p, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
fwer <- (1 + sum(perm_minp <= obs_min_p, na.rm = TRUE)) / (1 + sum(!is.na(perm_minp)))
calib <- mean(perm_p < 0.05, na.rm = TRUE)
pool_all <- merge(pool_all, VARS[, .(rsid, tested, literature)], by = "rsid", all.x = TRUE)
setorder(pool_all, p)
fwrite(pool_all[, !c("key", "z_obs")], file.path(out, "h3_pooled.tsv"), sep = "\t")
neg <- data.table(permutations = B, seed = SEED, tests_per_permutation = sum(!is.na(perm_p[1, ])),
                  share_nominal_p_below_0.05_under_null = calib, expected = 0.05,
                  observed_smallest_p = obs_min_p, permutation_family_wise_p_of_smallest = fwer,
                  null_smallest_p_median = median(perm_minp, na.rm = TRUE),
                  null_smallest_p_5th_pct = quantile(perm_minp, 0.05, na.rm = TRUE))
fwrite(neg, file.path(out, "h3_negative_control.tsv"), sep = "\t")
fwrite(data.table(perm = seq_len(B), min_p = perm_minp), file.path(out, "h3_negative_control_min_p.tsv.gz"), sep = "\t")
message("H3 pooled:"); print(pool_all[, .(gene, program, k, beta = round(beta, 3), se = round(se, 3), p = signif(p, 3), p_holm = signif(p_holm, 3), p_perm, I2 = round(I2, 2))])
message("negative control:"); print(neg)

# ---- HSD17B13 sensitivity: depth >= 15 at the site for both classes ----------------------------
pr15 <- prep_dose("dose_rs72613567", 5, 5, depth_min = 15)
add_counts("HSD17B13 depth>=15", pr15)
st15 <- stop_rule(pr15, 20)
if (st15$testable) {
  e15 <- cohort_est(pr15)
  s15 <- merge(pool_table(e15), pooled_context(pr15), by = "program", all.x = TRUE)
  s15[, `:=`(cohorts = st15$cohorts, n_with = st15$n_with, n_without = st15$n_without)]
  fwrite(e15, file.path(out, "h3_hsd17b13_depth15_cohort_estimates.tsv"), sep = "\t")
  fwrite(s15, file.path(out, "h3_hsd17b13_depth15_pooled.tsv"), sep = "\t")
  message("HSD17B13 depth>=15 sensitivity:"); print(s15[, .(program, k, beta = round(beta, 3), se = round(se, 3), p = signif(p, 3), I2 = round(I2, 2), cohorts)])
} else {
  fwrite(as.data.table(st15), file.path(out, "h3_hsd17b13_depth15_untestable.tsv"), sep = "\t")
}

# ---- Secondary: HSD17B13 carrier x PNPLA3 dosage ----------------------------------------------
inter_cohort <- function() {
  rows <- list()
  for (co in COHORTS) {
    d <- dat[cohort == co & !is.na(dose_rs72613567) & !is.na(dose_rs738409) & !is.na(stage) & !is.na(sex2)]
    cells <- if (nrow(d)) as.vector(table(factor(d$dose_rs72613567, 0:1), factor(d$dose_rs738409 > 0, c(FALSE, TRUE)))) else rep(0L, 4)
    row <- data.table(cohort = co, n = nrow(d), min_cell = min(cells), eligible = FALSE, reason = "")
    sf <- if (nrow(d)) merge_stage(d$stage) else factor(integer(0))
    if (min(cells) < 3) row$reason <- "a carrier x PNPLA3 cell has fewer than 3 persons"
    else if (nlevels(sf) < 2) row$reason <- "fewer than 2 stage levels"
    else {
      X <- cbind(1, car = d$dose_rs72613567, pn = d$dose_rs738409, int = d$dose_rs72613567 * d$dose_rs738409, design0(d, sf))
      f <- hc3(X, as.matrix(d[, ..YCOLS]), 4L)
      if (is.null(f)) row$reason <- "singular design" else {
        row$eligible <- TRUE
        rows[[co]] <- data.table(cohort = co, program = PNAMES, n = nrow(d), min_cell = min(cells), beta = f$beta, se = f$se)
      }
    }
    count_rows[[paste0("inter_", co)]] <<- row[, .(variant = "HSD17B13 x PNPLA3", cohort, persons = n, n_with = NA_integer_, n_without = NA_integer_,
                                                  n_stage_levels = NA_integer_, eligible, reason)]
  }
  rbindlist(rows)
}
if (stops[["rs72613567"]]$testable) {
  ie <- inter_cohort()
  if (nrow(ie)) {
    ip <- pool_table(ie)
    ip[, p_holm := p.adjust(p, "holm")]
    fwrite(ie, file.path(out, "h3_interaction_cohort_estimates.tsv"), sep = "\t")
    fwrite(ip, file.path(out, "h3_interaction_pooled.tsv"), sep = "\t")
    message("interaction:"); print(ip[, .(program, k, beta = round(beta, 3), se = round(se, 3), p = signif(p, 3), p_holm = signif(p_holm, 3), I2 = round(I2, 2))])
  }
}

# ---- Positive control: PNPLA3 G dosage vs fibrosis stage (numeric), per cohort, pooled ----------
pc_rows <- rbindlist(lapply(COHORTS, function(co) {
  d <- dat[cohort == co & !is.na(dose_rs738409) & !is.na(stage) & !is.na(sex2)]
  nw <- sum(d$dose_rs738409 > 0); nwo <- sum(d$dose_rs738409 == 0)
  if (nw < 5 || nwo < 5 || uniqueN(d$stage) < 2) return(data.table(cohort = co, n = nrow(d), n_with = nw, n_without = nwo, eligible = FALSE))
  sexF <- if (uniqueN(d$sex2) == 2) as.numeric(d$sex2 == "F") else NULL
  f <- hc3(cbind(1, d$dose_rs738409, sexF), matrix(d$stage, ncol = 1), 2L)
  data.table(cohort = co, n = nrow(d), n_with = nw, n_without = nwo, eligible = TRUE, beta = f$beta, se = f$se)
}), fill = TRUE)
pcr <- pc_rows[eligible == TRUE]
pcp <- dl_pool(pcr$beta, pcr$se)
pc_sum <- data.table(test = "PNPLA3 G dosage vs fibrosis stage (numeric), adjusted for sex", k = pcp$k, n = sum(pcr$n),
                     beta = pcp$beta, se = pcp$se, z = pcp$z, p_two_sided = pcp$p, p_one_sided_positive = pnorm(-pcp$z),
                     tau2 = pcp$tau2, I2 = pcp$I2, passes = pcp$beta > 0 && pnorm(-pcp$z) < 0.05)
fwrite(pc_rows, file.path(out, "positive_control_cohorts.tsv"), sep = "\t")
fwrite(pc_sum, file.path(out, "positive_control.tsv"), sep = "\t")
message("positive control:"); print(pc_sum)

# ---- H31: Pi*MZ vs Pi*MM, UPR program ----------------------------------------------------------
prpi <- prep_dose("dose_pi", 3, 5, ycols = "score_upr")
add_counts("Pi*MZ vs Pi*MM", prpi)
stpi <- stop_rule(prpi, 20)
h31 <- NULL
if (stpi$testable) {
  epi <- cohort_est(prpi, labels = "upr")
  ppi <- merge(pool_table(epi, "upr"), pooled_context(prpi, "upr"), by = "program", all.x = TRUE)
  set.seed(SEED)
  zp <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    db <- lapply(prpi, function(p) if (isTRUE(p$eligible)) p$dose[sample.int(length(p$dose))] else NULL)
    e <- cohort_est(prpi, dose_by = db, labels = "upr")
    if (nrow(e) >= 2) zp[b] <- dl_pool(e$beta, e$se)$z
  }
  ppi[, p_perm := (1 + sum(abs(zp) >= abs(z), na.rm = TRUE)) / (1 + sum(!is.na(zp)))]
  ppi[, `:=`(cohorts = stpi$cohorts, n_MZ = stpi$n_with, n_MM = stpi$n_without)]
  fwrite(epi, file.path(out, "h31_cohort_estimates.tsv"), sep = "\t")
  fwrite(ppi, file.path(out, "h31_pooled.tsv"), sep = "\t")
  message("H31:"); print(ppi[, .(k, n_MZ, n_MM, beta = round(beta, 3), se = round(se, 3), p = signif(p, 3), p_perm, I2 = round(I2, 2))])
  message("approximate smallest effect with 80% power (alpha 0.05, 2.8 x pooled SE): ", round(2.8 * ppi$se, 2), " SD of the score")
} else {
  fwrite(as.data.table(stpi), file.path(out, "h31_untestable.tsv"), sep = "\t")
  message("H31 untestable by the stop rule"); print(as.data.table(stpi))
}

fwrite(rbindlist(count_rows, fill = TRUE), file.path(out, "analysis_set_counts.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("done")
