#!/usr/bin/env Rscript
# H42 step 2: zonation score versus PNPLA3 I148M dosage in the five corrected bulk cohorts.
# Everything here follows PRESPEC_H42.md. Usage: 02_zonation_models.R GENOTYPE_RUN_DIR OUT_DIR
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(metafor)
})
args <- commandArgs(trailingOnly = TRUE)
geno_dir <- args[[1]]
out <- args[[2]]
stopifnot(dir.exists(geno_dir))
if (dir.exists(out)) stop("Refusing to overwrite ", out)
dir.create(out, recursive = TRUE)
SEED <- 20260923L
set.seed(SEED)
# optional third and fourth arguments lower the iteration counts for a smoke run only
NPERM <- if (length(args) >= 3L) as.integer(args[[3]]) else 1000L
NRAND <- if (length(args) >= 4L) as.integer(args[[4]]) else 1000L
COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

refit <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
dge <- readRDS(file.path(refit, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(refit, "meta_matched.rds")))
comp <- fread("Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz",
              select = c("sample_id", "Hepatocytes"))
degs <- fread(file.path(refit, "deg_results.csv"), select = c("gene", "symbol", "logFC", "padj"))
geno <- fread(file.path(geno_dir, "h42_unit_genotypes.tsv"))
sets <- fread("scripts/analysis/hypotheses_2026q4/H42_zonation/h42_zonation_gene_sets.tsv", fill = TRUE)

# ---------------------------------------------------------------- persons
stopifnot(setequal(colnames(dge$counts), meta$sample_id), nrow(meta) == 774L, !anyNA(meta$sex_final))
stopifnot(identical(as.character(dge$samples[meta$sample_id, "sex"]), meta$sex_final))
pers <- merge(meta[, .(sample_id, dataset, group_binary, sex = sex_final, fibrosis_stage, analysis_unit_id)],
              comp, by = "sample_id", all.x = TRUE)
setnames(pers, "Hepatocytes", "hep")
pers <- merge(pers, geno[, .(analysis_unit_id, dosage, dosage_depth20, depth, n_libraries, libraries_discordant)],
              by = "analysis_unit_id", all.x = TRUE)
pers[, fib := fifelse(is.na(fibrosis_stage) & group_binary == "Control", 0L, as.integer(fibrosis_stage))]
pers[, log10depth := log10(pmax(depth, 1))]
pers <- pers[match(colnames(dge$counts), sample_id)]
stopifnot(identical(pers$sample_id, colnames(dge$counts)), !anyNA(pers$hep), !anyNA(pers$sex))

# ---------------------------------------------------------------- expression
dge <- calcNormFactors(dge, method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 2)
ncpm <- cpm(dge)
symmap <- degs[gene %in% rownames(lcpm)]
pick_gene <- function(s) {
  g <- symmap[symbol == s, gene]
  if (!length(g)) return(NA_character_)
  if (length(g) == 1L) return(g)
  g[which.max(rowMeans(ncpm[g, , drop = FALSE]))]
}
sets[, gene := vapply(symbol, pick_gene, "")]
sets[gene == "", gene := NA_character_]
stopifnot(!anyNA(sets[primary == TRUE, gene]))
cat("primary pericentral:", sets[primary == TRUE & zone == "Pericentral", .N],
    " periportal:", sets[primary == TRUE & zone == "Periportal", .N], "\n")
fwrite(sets, file.path(out, "gene_sets_with_ensembl.tsv"), sep = "\t")

# ---------------------------------------------------------------- scores
build_scores <- function(set_dt, label) {
  res <- list(); ret <- list()
  for (ds in COHORTS) {
    idx <- which(pers$dataset == ds)
    keep <- rowMeans(ncpm[, idx, drop = FALSE] >= 1) >= 0.5
    pc <- set_dt[zone == "Pericentral" & !is.na(gene), gene]; pc <- pc[keep[pc]]
    pp <- set_dt[zone == "Periportal" & !is.na(gene), gene]; pp <- pp[keep[pp]]
    sy <- function(g) paste(symmap[match(g, gene), symbol], collapse = ",")
    ret[[ds]] <- data.table(set = label, dataset = ds, n_pc = length(pc), n_pp = length(pp),
                            pc_genes = sy(pc), pp_genes = sy(pp))
    if (length(pc) < 8L || length(pp) < 8L) next
    zm <- function(g) { x <- lcpm[g, idx, drop = FALSE]; (x - rowMeans(x)) / apply(x, 1, sd) }
    zpc <- colMeans(zm(pc)); zpp <- colMeans(zm(pp))
    sc <- data.table(idx = idx, PC = zpc, PP = zpp, Z = zpc - zpp)
    for (v in c("PC", "PP", "Z")) set(sc, j = v, value = as.numeric(scale(sc[[v]])))
    res[[ds]] <- sc
  }
  list(scores = res, retention = rbindlist(ret))
}

# ---------------------------------------------------------------- per-cohort OLS via QR
make_obj <- function(ds, sc, outcome, xcol = "dosage", kind = "dosage", covs = c("sex", "hep", "fib", "group_binary"),
                     extra = NULL, stage_mode = FALSE) {
  p <- pers[sc$idx]
  y <- sc[[outcome]]
  x <- p[[xcol]]
  ok <- !is.na(x) & !is.na(p$sex) & !is.na(p$hep)
  if (!is.null(extra)) ok <- ok & !is.na(p[[extra]])
  if (kind == "GG") ok <- ok & x %in% c(0, 2)
  if (stage_mode) ok <- ok & !is.na(p$fib)
  d <- p[ok]; yy <- y[ok]; xx <- x[ok]
  if (kind == "carrier") xx <- as.numeric(xx >= 1)
  if (kind == "GG") xx <- as.numeric(xx == 2)
  n1 <- sum(xx > 0)
  n0 <- sum(xx == 0)
  info <- data.table(dataset = ds, n = length(xx), n_exposure0 = n0, n_exposure_pos = n1)
  if (stage_mode) {
    if (length(xx) < 20L || var(xx) == 0) return(list(ok = FALSE, info = info))
  } else if (length(xx) < 20L || n0 < 5L || n1 < 5L) return(list(ok = FALSE, info = info))
  cols <- list()
  if ("sex" %in% covs && nlevels(factor(d$sex)) >= 2L) cols$sex <- factor(d$sex)
  if ("hep" %in% covs) cols$hep <- d$hep
  if ("fib" %in% covs && !anyNA(d$fib) && var(d$fib) > 0) cols$fib <- d$fib
  if ("group_binary" %in% covs && nlevels(factor(d$group_binary)) >= 2L) cols$group_binary <- factor(d$group_binary)
  if (!is.null(extra)) cols$extra <- d[[extra]]
  X0 <- if (length(cols)) model.matrix(~ ., data = as.data.frame(cols)) else matrix(1, length(xx), 1)   # intercept only
  qrX <- qr(X0)
  o <- list(ok = TRUE, info = info, x = xx, y = yy, qrX = qrX, rows = which(ok), df = length(xx) - qrX$rank - 1L,
            terms = paste(names(cols), collapse = "+"))
  o$x_r <- qr.resid(qrX, xx); o$y_r <- qr.resid(qrX, yy)
  o$strata <- paste(d$sex, d$fib, sep = "|")
  o$strata_idx <- split(seq_along(xx), o$strata)
  o
}
est_one <- function(o, x = NULL, y_r = NULL) {
  xr <- if (is.null(x)) o$x_r else qr.resid(o$qrX, x)
  sxx <- sum(xr^2)
  if (sxx < 1e-10) return(c(NA_real_, NA_real_))
  yr <- if (is.null(y_r)) o$y_r else y_r
  b <- sum(xr * yr) / sxx
  rss <- sum(yr^2) - b^2 * sxx
  c(b, sqrt(rss / o$df / sxx))
}
dl <- function(b, se) {
  k <- length(b); v <- se^2; w <- 1 / v; bf <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bf)^2); C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (k - 1)) / C); ws <- 1 / (v + tau2)
  bd <- sum(ws * b) / sum(ws); sed <- sqrt(1 / sum(ws))
  list(beta = bd, se = sed, z = bd / sed, p = 2 * pnorm(-abs(bd / sed)), tau2 = tau2, Q = Q,
       Qp = if (k > 1) pchisq(Q, k - 1, lower.tail = FALSE) else NA_real_,
       I2 = if (k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) else 0, k = k)
}
pool_objs <- function(objs) {
  objs <- Filter(function(o) o$ok, objs)
  if (length(objs) < 1L) return(NULL)
  e <- t(vapply(objs, est_one, numeric(2)))
  keep <- is.finite(e[, 1]) & is.finite(e[, 2])
  e <- e[keep, , drop = FALSE]; objs <- objs[keep]
  if (nrow(e) < 3L) return(list(cohorts = objs, e = e, dl = NULL))   # SR3
  list(cohorts = objs, e = e, dl = dl(e[, 1], e[, 2]))
}
hksj <- function(e) {
  if (is.null(e) || nrow(e) < 3L) return(c(NA, NA, NA))
  f <- rma(yi = e[, 1], sei = e[, 2], method = "DL", test = "knha")
  c(f$ci.lb, f$ci.ub, f$pval)
}
run_variant <- function(label, sc_list, outcome = "Z", xcol = "dosage", kind = "dosage",
                        covs = c("sex", "hep", "fib", "group_binary"), extra = NULL, drop = character()) {
  use <- setdiff(names(sc_list), drop)
  empty_row <- data.table(label = label, outcome = outcome, kind = kind, k_cohorts = NA_integer_, n_persons = NA_integer_,
                          n_exposed = NA_integer_, beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_, tau2 = NA_real_,
                          I2 = NA_real_, Q_p = NA_real_, hk_lo = NA_real_, hk_hi = NA_real_, hk_p = NA_real_, mde80 = NA_real_)
  if (!length(use)) return(list(row = empty_row, per = data.table(), pool = NULL, objs = list()))
  objs <- lapply(use, function(ds) make_obj(ds, sc_list[[ds]], outcome, xcol, kind, covs, extra))
  names(objs) <- use
  pl <- pool_objs(objs)
  per <- rbindlist(lapply(names(objs), function(ds) {
    o <- objs[[ds]]
    if (!o$ok) return(cbind(o$info, data.table(included = FALSE, beta = NA_real_, se = NA_real_, p = NA_real_, terms = NA_character_)))
    e <- est_one(o)
    cbind(o$info, data.table(included = TRUE, beta = e[1], se = e[2],
                             p = 2 * pt(-abs(e[1] / e[2]), o$df), terms = o$terms))
  }))
  per[, `:=`(label = label, outcome = outcome, kind = kind)]
  row <- copy(empty_row)
  if (!is.null(pl) && !is.null(pl$dl)) {
    d <- pl$dl; h <- hksj(pl$e)
    row[, `:=`(k_cohorts = d$k, n_persons = sum(per[included == TRUE, n]), n_exposed = sum(per[included == TRUE, n_exposure_pos]),
               beta = d$beta, se = d$se, z = d$z, p = d$p, tau2 = d$tau2, I2 = d$I2, Q_p = d$Qp,
               hk_lo = h[1], hk_hi = h[2], hk_p = h[3], mde80 = 2.8 * d$se)]
  }
  list(row = row, per = per, pool = pl, objs = objs)
}

# ---------------------------------------------------------------- 0. inputs summary (callable fraction first)
inputs <- pers[, .(persons = .N, callable = sum(!is.na(dosage)), callable_fraction = mean(!is.na(dosage)),
                   dosage0 = sum(dosage == 0, na.rm = TRUE), dosage1 = sum(dosage == 1, na.rm = TRUE),
                   dosage2 = sum(dosage == 2, na.rm = TRUE), carriers = sum(dosage >= 1, na.rm = TRUE),
                   callable_depth20 = sum(!is.na(dosage_depth20)), median_depth = as.numeric(median(depth, na.rm = TRUE))),
               by = dataset][match(COHORTS, dataset)]
fwrite(inputs, file.path(out, "00_inputs_callable_by_cohort.tsv"), sep = "\t")
print(inputs)
total_carriers <- sum(inputs$carriers)
cat("total alt-carrier persons (dosage >= 1):", total_carriers, "\n")
if (total_carriers < 20L) {
  writeLines("SR1 tripped: fewer than 20 alt-carrier persons. No inference run.", file.path(out, "STOPPED_SR1.txt"))
  quit(status = 0)
}

# ---------------------------------------------------------------- 1. primary scores and P1
prim <- build_scores(sets[primary == TRUE], "primary")
fwrite(prim$retention, file.path(out, "01_marker_retention_primary.tsv"), sep = "\t")
sc_prim <- prim$scores
P1 <- run_variant("P1_dosage_on_Z", sc_prim, "Z")
if (is.null(P1$pool) || is.null(P1$pool$dl)) {
  fwrite(P1$per, file.path(out, "03_per_cohort_estimates.tsv"), sep = "\t")
  writeLines("SR3 tripped: fewer than 3 cohorts meet the entry rule. No pooled inference.", file.path(out, "STOPPED_SR3.txt"))
  quit(status = 0)
}
# DL implementation check against metafor
chk <- rma(yi = P1$pool$e[, 1], sei = P1$pool$e[, 2], method = "DL")
stopifnot(abs(chk$beta[1] - P1$row$beta) < 1e-8, abs(chk$tau2 - P1$row$tau2) < 1e-8)

# ---------------------------------------------------------------- 2. secondary family (Holm, m = 4)
S_PC <- run_variant("S1_dosage_on_PC", sc_prim, "PC")
S_PP <- run_variant("S2_dosage_on_PP", sc_prim, "PP")
S_CAR <- run_variant("S3_carrier_on_Z", sc_prim, "Z", kind = "carrier")
S_GG <- run_variant("S4_GG_vs_CC_on_Z", sc_prim, "Z", kind = "GG")
sec <- rbindlist(list(S_PC$row, S_PP$row, S_CAR$row, S_GG$row))
sec[, p_holm := p.adjust(p, method = "holm")]
main <- rbind(P1$row, sec, fill = TRUE)
fwrite(main, file.path(out, "02_meta_results_primary_and_secondary.tsv"), sep = "\t")
per_all <- rbindlist(list(P1$per, S_PC$per, S_PP$per, S_CAR$per, S_GG$per))
fwrite(per_all, file.path(out, "03_per_cohort_estimates.tsv"), sep = "\t")

# ---------------------------------------------------------------- 3. permutation null (NC1)
perm_null <- function(objs, type, n) {
  set.seed(SEED)
  res <- matrix(NA_real_, n, 2)
  for (i in seq_len(n)) {
    e <- matrix(NA_real_, length(objs), 2)
    for (j in seq_along(objs)) {
      o <- objs[[j]]
      xp <- if (type == "within_cohort") sample(o$x) else {
        v <- o$x
        for (ix in o$strata_idx) v[ix] <- o$x[ix][sample.int(length(ix))]
        v
      }
      e[j, ] <- est_one(o, x = xp)
    }
    e <- e[is.finite(e[, 1]) & is.finite(e[, 2]), , drop = FALSE]
    d <- dl(e[, 1], e[, 2])
    res[i, ] <- c(d$beta, d$z)
  }
  data.table(iter = seq_len(n), type = type, beta = res[, 1], z = res[, 2])
}
pobjs <- Filter(function(o) o$ok, P1$objs)
perm <- rbindlist(list(perm_null(pobjs, "within_cohort", NPERM), perm_null(pobjs, "within_cohort_sex_stage", NPERM)))
fwrite(perm, file.path(out, "04_permutation_null.tsv"), sep = "\t")
perm_sum <- perm[, .(n_perm = .N, obs_beta = P1$row$beta, obs_z = P1$row$z,
                     emp_p_beta = (1 + sum(abs(beta) >= abs(P1$row$beta))) / (.N + 1),
                     emp_p_z = (1 + sum(abs(z) >= abs(P1$row$z))) / (.N + 1),
                     null_sd_z = sd(z), frac_abs_z_gt_196 = mean(abs(z) > 1.96)), by = type]
fwrite(perm_sum, file.path(out, "05_permutation_summary.tsv"), sep = "\t")
print(perm_sum)

# ---------------------------------------------------------------- 4. random matched gene-set null (NC2)
pc_int <- Reduce(intersect, lapply(COHORTS, function(ds) {
  idx <- which(pers$dataset == ds); keep <- rowMeans(ncpm[, idx, drop = FALSE] >= 1) >= 0.5
  g <- sets[primary == TRUE & zone == "Pericentral" & !is.na(gene), gene]; g[keep[g]] }))
pp_int <- Reduce(intersect, lapply(COHORTS, function(ds) {
  idx <- which(pers$dataset == ds); keep <- rowMeans(ncpm[, idx, drop = FALSE] >= 1) >= 0.5
  g <- sets[primary == TRUE & zone == "Periportal" & !is.na(gene), gene]; g[keep[g]] }))
cat("intersection marker sets: PC", length(pc_int), " PP", length(pp_int), "\n")
all_marker <- sets[!is.na(gene), gene]
expr_all <- Reduce(`&`, lapply(COHORTS, function(ds) { idx <- which(pers$dataset == ds)
  rowMeans(ncpm[, idx, drop = FALSE] >= 1) >= 0.5 }))
mean_expr <- rowMeans(lcpm)
brks <- quantile(mean_expr[expr_all], probs = seq(0, 1, 0.1))
dec <- pmin(pmax(findInterval(mean_expr, brks, rightmost.closed = TRUE), 1L), 10L)
names(dec) <- names(mean_expr)
pool <- setdiff(names(expr_all)[expr_all], all_marker)
pool_by_dec <- split(pool, dec[pool])
ZM <- lapply(COHORTS, function(ds) { idx <- which(pers$dataset == ds)
  x <- lcpm[pool, idx, drop = FALSE]; (x - rowMeans(x)) / apply(x, 1, sd) })
names(ZM) <- COHORTS
idx_list <- lapply(COHORTS, function(ds) which(pers$dataset == ds)); names(idx_list) <- COHORTS
# observed P1 with intersection sets (apples to apples with the random sets)
int_dt <- sets[primary == TRUE & gene %in% c(pc_int, pp_int)]
sc_int <- build_scores(int_dt, "primary_intersection")$scores
P1_int <- run_variant("P1_intersection_sets", sc_int, "Z")
set.seed(SEED)
draw <- function(markers, used) {
  out_g <- character(length(markers))
  for (i in seq_along(markers)) {
    cand <- setdiff(pool_by_dec[[as.character(dec[markers[i]])]], used)
    out_g[i] <- sample(cand, 1L); used <- c(used, out_g[i])
  }
  out_g
}
rand <- data.table(iter = seq_len(NRAND), beta = NA_real_, z = NA_real_)
for (i in seq_len(NRAND)) {
  g_pc <- draw(pc_int, character()); g_pp <- draw(pp_int, g_pc)
  e <- matrix(NA_real_, length(P1_int$objs), 2); k <- 0L
  for (ds in names(P1_int$objs)) {
    o <- P1_int$objs[[ds]]; if (!o$ok) next
    zr <- colMeans(ZM[[ds]][g_pc, , drop = FALSE]) - colMeans(ZM[[ds]][g_pp, , drop = FALSE])
    zr <- as.numeric(scale(zr))   # SD 1 within cohort, as for the observed score
    # rows of the cohort that entered the fit (same persons as the observed fit)
    yv <- zr[o$rows]
    k <- k + 1L; e[k, ] <- est_one(o, y_r = qr.resid(o$qrX, yv))
  }
  e <- e[seq_len(k), , drop = FALSE]; e <- e[is.finite(e[, 1]) & is.finite(e[, 2]), , drop = FALSE]
  d <- dl(e[, 1], e[, 2]); rand[i, `:=`(beta = d$beta, z = d$z)]
}
fwrite(rand, file.path(out, "06_random_gene_set_null.tsv"), sep = "\t")
rand_sum <- data.table(n_random = NRAND, obs_z_primary = P1$row$z, obs_z_intersection = P1_int$row$z,
                       obs_beta_intersection = P1_int$row$beta,
                       emp_p_primary = (1 + sum(abs(rand$z) >= abs(P1$row$z))) / (NRAND + 1),
                       emp_p_intersection = (1 + sum(abs(rand$z) >= abs(P1_int$row$z))) / (NRAND + 1),
                       q95_abs_z = quantile(abs(rand$z), 0.95), sd_z = sd(rand$z), frac_abs_z_gt_196 = mean(abs(rand$z) > 1.96),
                       n_pc_int = length(pc_int), n_pp_int = length(pp_int))
fwrite(rand_sum, file.path(out, "07_random_gene_set_summary.tsv"), sep = "\t")
print(rand_sum)

# ---------------------------------------------------------------- 5. positive controls
vis <- fread("Analysis/Spatial/results/validation_bulk/deg_zonation_combined.csv", select = c("gene", "dataset", "mean_PP1", "mean_PC1", "zonation_class"))
retained <- unique(unlist(strsplit(prim$retention$pc_genes, ",")))
retained_pp <- unique(unlist(strsplit(prim$retention$pp_genes, ",")))
pca_rows <- merge(rbind(data.table(symbol = retained, zone = "Pericentral"), data.table(symbol = retained_pp, zone = "Periportal")),
                  vis, by.x = "symbol", by.y = "gene", all.x = FALSE)
pca_rows[, agree_row := fifelse(zone == "Pericentral", mean_PC1 > mean_PP1, mean_PP1 > mean_PC1)]
# a gene agrees only if every Visium row where it appears has the expected zone higher (Addendum 1, rule 4)
pca <- pca_rows[, .(zone = zone[1], n_rows = .N, agree = all(agree_row)), by = symbol]
fwrite(pca_rows, file.path(out, "08_positive_control_A_visium_rows.tsv"), sep = "\t")
fwrite(pca, file.path(out, "08_positive_control_A_visium_genes.tsv"), sep = "\t")
pcA <- data.table(control = "PC-A_visium_concordance", n_genes = nrow(pca), n_agree = sum(pca$agree),
                  fraction = mean(pca$agree), pass = nrow(pca) >= 8L && mean(pca$agree) >= 0.70,
                  frac_Guilliams = pca_rows[dataset == "Guilliams et al.", mean(agree_row)],
                  frac_Vu = pca_rows[dataset == "Vu et al.", mean(agree_row)])
visium_discordant <- pca[agree == FALSE, symbol]
stage_objs <- function(outcome) {
  # fibrosis-stage coefficient from the same model as P1 (dosage is a covariate here), Addendum 1 rule 5
  o <- lapply(COHORTS, function(ds) if (ds %in% names(sc_prim)) make_obj(ds, sc_prim[[ds]], outcome, xcol = "fib",
                                      kind = "stage", covs = c("sex", "hep", "group_binary"), extra = "dosage", stage_mode = TRUE) else list(ok = FALSE, info = data.table(dataset = ds)))
  names(o) <- COHORTS; o
}
stage_res <- lapply(c("PP", "PC", "Z"), function(y) {
  ob <- stage_objs(y); pl <- pool_objs(ob)
  per <- rbindlist(lapply(names(ob), function(ds) { o <- ob[[ds]]
    if (!o$ok) return(data.table(dataset = ds, outcome = y, beta = NA_real_, se = NA_real_, n = NA_integer_))
    e <- est_one(o); data.table(dataset = ds, outcome = y, beta = e[1], se = e[2], n = o$info$n) }))
  if (is.null(pl) || is.null(pl$dl))
    return(list(per = per, row = data.table(outcome = y, k = NA_integer_, beta = NA_real_, se = NA_real_, z = NA_real_,
                                            p = NA_real_, tau2 = NA_real_, I2 = NA_real_)))
  list(per = per, row = data.table(outcome = y, k = pl$dl$k, beta = pl$dl$beta, se = pl$dl$se, z = pl$dl$z, p = pl$dl$p,
                                   tau2 = pl$dl$tau2, I2 = pl$dl$I2))
})
stage_tab <- rbindlist(lapply(stage_res, `[[`, "row")); fwrite(stage_tab, file.path(out, "09_positive_control_B_stage.tsv"), sep = "\t")
fwrite(rbindlist(lapply(stage_res, `[[`, "per")), file.path(out, "09b_stage_per_cohort.tsv"), sep = "\t")
pp_row <- stage_tab[outcome == "PP"]; z_row <- stage_tab[outcome == "Z"]
pcB <- data.table(control = "PC-B_stage_direction", PP_stage_beta = pp_row$beta, PP_stage_p = pp_row$p, Z_stage_beta = z_row$beta,
                  Z_stage_p = z_row$p, pass = isTRUE(pp_row$beta < 0 && pp_row$p < 0.05 && z_row$beta > 0))
controls <- rbind(pcA, pcB, fill = TRUE)
fwrite(controls, file.path(out, "10_positive_controls.tsv"), sep = "\t")
print(controls); print(stage_tab)

# ---- post hoc diagnostics of the score itself (added after the smoke run showed PC-B failing; no genotype involved).
# They do not change the reading. They show whether the score responds to disease or stage when group and
# hepatocyte fraction are not both held fixed.
posthoc_row <- function(label, outcome, xcol, covs, kind, stage_mode = FALSE) {
  ob <- lapply(COHORTS, function(ds) if (ds %in% names(sc_prim)) make_obj(ds, sc_prim[[ds]], outcome, xcol = xcol, kind = kind,
                covs = covs, stage_mode = stage_mode) else list(ok = FALSE, info = data.table(dataset = ds)))
  names(ob) <- COHORTS
  pl <- pool_objs(ob)
  per <- rbindlist(lapply(names(ob), function(ds) { o <- ob[[ds]]
    if (!o$ok) return(data.table(dataset = ds, beta = NA_real_, se = NA_real_, n = NA_integer_))
    e <- est_one(o); data.table(dataset = ds, beta = e[1], se = e[2], n = o$info$n) }))
  per[, `:=`(label = label, outcome = outcome)]
  if (is.null(pl) || is.null(pl$dl)) return(list(row = data.table(label = label, outcome = outcome, k = NA_integer_, beta = NA_real_,
                    se = NA_real_, z = NA_real_, p = NA_real_, tau2 = NA_real_, I2 = NA_real_), per = per))
  list(row = data.table(label = label, outcome = outcome, k = pl$dl$k, beta = pl$dl$beta, se = pl$dl$se, z = pl$dl$z,
                        p = pl$dl$p, tau2 = pl$dl$tau2, I2 = pl$dl$I2), per = per)
}
pers[, disease := as.integer(group_binary == "Disease")]
ph <- list()
for (y in c("PP", "PC", "Z")) {
  ph[[length(ph) + 1L]] <- posthoc_row("PH1_disease_vs_control_adj_sex_hep", y, "disease", c("sex", "hep"), "contrast")
  ph[[length(ph) + 1L]] <- posthoc_row("PH2_disease_vs_control_adj_sex", y, "disease", c("sex"), "contrast")
  ph[[length(ph) + 1L]] <- posthoc_row("PH3_stage_adj_sex_only", y, "fib", c("sex"), "stage", stage_mode = TRUE)
  ph[[length(ph) + 1L]] <- posthoc_row("PH4_stage_adj_sex_hep_no_group", y, "fib", c("sex", "hep"), "stage", stage_mode = TRUE)
}
fwrite(rbindlist(lapply(ph, `[[`, "row")), file.path(out, "17_posthoc_score_diagnostics.tsv"), sep = "\t")
fwrite(rbindlist(lapply(ph, `[[`, "per")), file.path(out, "17b_posthoc_score_diagnostics_per_cohort.tsv"), sep = "\t")
print(rbindlist(lapply(ph, `[[`, "row")))

# ---------------------------------------------------------------- 6. sensitivity (descriptive)
deg_set <- degs[!is.na(padj) & padj < 0.05 & abs(logFC) > 0.5, gene]
S_full <- build_scores(sets[full_list == TRUE], "full_list")
S_nondeg <- build_scores(sets[primary == TRUE & !(gene %in% deg_set)], "primary_non_DEG")
S_visium <- build_scores(sets[primary == TRUE & !(symbol %in% visium_discordant)], "primary_visium_concordant")
fwrite(rbind(S_full$retention, S_nondeg$retention, S_visium$retention), file.path(out, "11_marker_retention_sensitivity.tsv"), sep = "\t")
sens <- list(
  run_variant("S_full_repo_list", S_full$scores, "Z"),
  run_variant("S_primary_non_DEG_markers", S_nondeg$scores, "Z"),
  run_variant("S_primary_visium_concordant_markers", S_visium$scores, "Z"),
  run_variant("S_genotype_depth20", sc_prim, "Z", xcol = "dosage_depth20"),
  run_variant("S_plus_log10_genotype_depth", sc_prim, "Z", extra = "log10depth"),
  run_variant("S_no_hepatocyte_fraction_adjustment", sc_prim, "Z", covs = c("sex", "fib", "group_binary")))
for (ds in names(sc_prim)) sens[[length(sens) + 1L]] <- run_variant(paste0("S_leave_out_", ds), sc_prim, "Z", drop = ds)
sens_tab <- rbindlist(lapply(sens, `[[`, "row"))
fwrite(sens_tab, file.path(out, "12_sensitivity.tsv"), sep = "\t")
# pooled regression with cohort fixed effects
pooled_df <- rbindlist(lapply(names(sc_prim), function(ds) { p <- pers[sc_prim[[ds]]$idx]; cbind(p[, .(dataset, dosage, sex, hep, fib, group_binary)], Z = sc_prim[[ds]]$Z) }))
pooled_df[, `:=`(fib_missing = as.integer(is.na(fib)), fib0 = fifelse(is.na(fib), 0L, fib))]
pooled_df <- pooled_df[!is.na(dosage)]
fit <- lm(Z ~ dosage + dataset + sex + hep + fib0 + fib_missing + group_binary, data = pooled_df)
cf <- summary(fit)$coefficients["dosage", ]
fwrite(data.table(label = "S_pooled_fixed_cohort_effects", n = nrow(pooled_df), beta = cf[1], se = cf[2], p = cf[4]),
       file.path(out, "13_pooled_fixed_effects.tsv"), sep = "\t")

# ---------------------------------------------------------------- 7. person-level audit table and reading
person_tab <- rbindlist(lapply(names(sc_prim), function(ds) { p <- pers[sc_prim[[ds]]$idx]
  cbind(p[, .(analysis_unit_id, sample_id, dataset, sex, hep, fib, group_binary, dosage, depth, n_libraries)], sc_prim[[ds]][, .(PC, PP, Z)]) }))
fwrite(person_tab, file.path(out, "14_person_level_table.tsv"), sep = "\t")
by_dose <- person_tab[!is.na(dosage), .(n = .N, mean_Z = mean(Z), sd_Z = sd(Z), mean_PC = mean(PC), mean_PP = mean(PP)), by = .(dataset, dosage)][order(dataset, dosage)]
fwrite(by_dose, file.path(out, "16_unadjusted_score_by_dosage.tsv"), sep = "\t")
print(by_dose)
reading <- {
  ctrl_ok <- isTRUE(pcA$pass) && isTRUE(pcB$pass)
  p1 <- P1$row
  p1_neg_sig <- p1$p < 0.05 && p1$beta < 0
  perm_ok <- perm_sum[type == "within_cohort", emp_p_beta] < 0.05
  nc2_ok <- abs(p1$z) > rand_sum$q95_abs_z
  if (!ctrl_ok) "uninformative_positive_control_failed"
  else if (!p1_neg_sig) "not_supported"
  else if (perm_ok && nc2_ok) "supported"
  else "negative_control_failed"
}
summ <- data.table(item = c("reading", "P1_beta", "P1_se", "P1_p", "P1_k_cohorts", "P1_n_persons", "total_carriers", "controls_pass"),
                   value = c(reading, P1$row$beta, P1$row$se, P1$row$p, P1$row$k_cohorts, P1$row$n_persons, total_carriers,
                             isTRUE(pcA$pass) && isTRUE(pcB$pass)))
fwrite(summ, file.path(out, "15_summary.tsv"), sep = "\t")
print(summ); print(main)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("done:", out, "\n")
