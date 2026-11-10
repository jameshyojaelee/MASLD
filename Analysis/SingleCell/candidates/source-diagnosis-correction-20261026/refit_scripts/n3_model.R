#!/usr/bin/env Rscript
# N3 step 3: donor-level limma-voom models for MASH drug targets in hepatocytes.
# Implements PRESPEC_N3.md. Usage: 03_donor_model.R <RESULT_DIR> <OUT_SUBDIR> <N_PERM> <N_CORES>
suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(parallel); library(jsonlite)
})
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 4L) stop("Usage: 03_donor_model.R RESULT_DIR OUT_SUBDIR N_PERM N_CORES")
res <- argv[[1]]; out <- file.path(res, argv[[2]]); B <- as.integer(argv[[3]]); ncores <- as.integer(argv[[4]])
if (dir.exists(out)) stop("Refusing to overwrite ", out)
dir.create(out, recursive = TRUE)
SEED <- 20260923L
set.seed(SEED)
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
bulk_file <- file.path(repo, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
tick <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

dm <- fread(file.path(res, "donor_meta.tsv"))
# The original source identifies these NAFLD explants as cirrhotic, not simple
# steatosis (Filliol Nature 2022 Methods, doi:10.1038/s41586-022-05289-6).
dm[dataset == "GSE174748" & donor %in% c("SRR14602057", "SRR14602058"),
   `:=`(disease_stage_coarse = "Cirrhosis", stage_ord = 3)]
cnt <- fread(file.path(res, "hep_donor_counts.tsv.gz"))
tl <- fread(file.path(res, "target_list.tsv"))
cts <- fread(file.path(res, "target_celltype_summary.tsv"))
stopifnot(!anyDuplicated(cnt$symbol))

inc <- dm[included == TRUE]
one_donor_ds <- inc[, .N, by = dataset][N < 2L, dataset]
dropped_single <- inc[dataset %in% one_donor_ds, donor]
inc <- inc[!dataset %in% one_donor_ds]
inc[, stage := factor(disease_stage_coarse, levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]
inc[, sex := factor(sex_call, levels = c("Female", "Male"))]
inc[, dataset := factor(dataset)]
inc <- inc[order(dataset, stage, donor)]
mat <- as.matrix(cnt[, inc$donor, with = FALSE]); rownames(mat) <- cnt$symbol
storage.mode(mat) <- "double"
n_by_stage <- inc[, .N, by = stage][order(stage)]
fwrite(inc[, .N, by = .(dataset, stage)][order(dataset, stage)], file.path(out, "donor_counts_by_dataset_stage.tsv"), sep = "\t")
fwrite(inc[, .(donor, dataset, stage, stage_ord, sex, n_hep_cells, hep_umis, F_documented)],
       file.path(out, "donors_in_model.tsv"), sep = "\t")
n_h <- sum(inc$stage == "Healthy"); n_c <- sum(inc$stage == "Cirrhosis")
underpowered <- n_h < 5L || n_c < 5L
tick("donors in model:", nrow(inc), " healthy:", n_h, " cirrhosis:", n_c, " underpowered:", underpowered,
     " single-donor datasets dropped:", paste(dropped_single, collapse = ","))

fit_voom <- function(counts, design, genes = NULL) {
  y <- DGEList(counts)
  keep <- if (is.null(genes)) filterByExpr(y, design, min.count = 10, min.total.count = 15) else rownames(y) %in% genes
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voomWithQualityWeights(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  list(y = y, v = v, fit = fit)
}
coef_table <- function(fit, coef) {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  data.table(symbol = rownames(tt), logFC = tt$logFC, SE = unname(se[rownames(tt)]), t = tt$t,
             P.Value = tt$P.Value, AveExpr = tt$AveExpr)
}
check_rank <- function(design, label) {
  ne <- nonEstimable(design)
  if (!is.null(ne)) stop("Non-estimable coefficients in ", label, ": ", paste(ne, collapse = ", "))
}

# ---------------------------------------------------------------- M1 primary
design1 <- model.matrix(~ dataset + sex + stage_ord, data = inc)
check_rank(design1, "M1")
m1 <- fit_voom(mat, design1)
tested <- intersect(tl$symbol, rownames(m1$y))
untested <- setdiff(tl$symbol, tested)
tick("genes after filter:", nrow(m1$y), " targets testable:", length(tested), " not detected:", length(untested))
m1_all <- coef_table(m1$fit, "stage_ord")
m1_all[, genome_BH := p.adjust(P.Value, "BH")]
fwrite(m1_all, file.path(out, "M1_all_genes.tsv.gz"), sep = "\t")
m1_t <- m1_all[symbol %in% tested]
m1_t[, q_bh := p.adjust(P.Value, "BH")]

# ------------------------------------------------------- M2 categorical stage
design2 <- model.matrix(~ stage + dataset + sex, data = inc)
check_rank(design2, "M2")
m2 <- fit_voom(mat, design2, genes = rownames(m1$y))
m2_list <- lapply(c("stageSteatosis", "stageSteatohepatitis", "stageCirrhosis"), function(cf) {
  x <- coef_table(m2$fit, cf)[symbol %in% tested]; x[, contrast := sub("^stage", "", cf)]; x })
m2_t <- rbindlist(m2_list)
m2_t[, q_bh := p.adjust(P.Value, "BH")]     # one family: testable targets x 3 contrasts
fwrite(m2_t, file.path(out, "M2_categorical_targets.tsv"), sep = "\t")

# -------------------------------------------------------------- S2 documented F
s2_t <- NULL
s2 <- inc[dataset == "GSE202379" & !is.na(F_documented)]
s2[, sex := droplevels(sex)]
if (nrow(s2) >= 20L) {
  d2 <- model.matrix(~ sex + F_documented, data = s2); check_rank(d2, "S2")
  m <- as.matrix(cnt[, s2$donor, with = FALSE]); rownames(m) <- cnt$symbol; storage.mode(m) <- "double"
  s2fit <- fit_voom(m, d2)
  s2_t <- coef_table(s2fit$fit, "F_documented")[symbol %in% tl$symbol]
  s2_t[, q_bh := p.adjust(P.Value, "BH")]
  s2_t[, n_donors := nrow(s2)]
  fwrite(s2_t, file.path(out, "S2_GSE202379_documented_F_targets.tsv"), sep = "\t")
}
tick("S2 donors:", nrow(s2))

# ------------------------------------------------- S3 leave one dataset out
s3_list <- list()
for (ds in c("GSE244832", "GSE202379", "GSE174748")) {
  s <- inc[dataset != ds]; s[, `:=`(dataset = droplevels(dataset), sex = droplevels(sex))]
  d3 <- model.matrix(~ dataset + sex + stage_ord, data = s)
  ne <- nonEstimable(d3)
  if (!is.null(ne)) { s3_list[[ds]] <- data.table(left_out = ds, note = paste("non-estimable:", paste(ne, collapse = ","))); next }
  m <- as.matrix(cnt[, s$donor, with = FALSE]); rownames(m) <- cnt$symbol; storage.mode(m) <- "double"
  f3 <- fit_voom(m, d3)
  x <- coef_table(f3$fit, "stage_ord")[symbol %in% tl$symbol]
  x[, `:=`(left_out = ds, n_donors = nrow(s))]
  s3_list[[ds]] <- x
}
s3 <- rbindlist(s3_list, fill = TRUE)
if ("logFC" %in% names(s3)) {
  s3 <- merge(s3, m1_t[, .(symbol, logFC_M1 = logFC)], by = "symbol", all.x = TRUE)
  s3[, same_sign_as_M1 := sign(logFC) == sign(logFC_M1)]
}
fwrite(s3, file.path(out, "S3_leave_one_dataset_out_targets.tsv"), sep = "\t")

# ------------------------------------------------------------ positive control
ctrl <- c("GNMT", "MAT1A", "CYP2C19")
ctrl_in <- intersect(ctrl, rownames(m1$v$E))
pc <- list(genes_requested = ctrl, genes_in_model = ctrl_in, pass = FALSE, note = NA_character_)
if (length(ctrl_in) >= 2L) {
  Z <- t(scale(t(m1$v$E[ctrl_in, , drop = FALSE])))
  dd <- as.data.frame(inc[match(colnames(m1$v$E), donor), .(dataset, sex, stage_ord)])
  dd$comp <- colMeans(Z)
  lf <- lm(comp ~ dataset + sex + stage_ord, data = dd)
  cf <- summary(lf)$coefficients["stage_ord", ]
  pc$estimate <- unname(cf["Estimate"]); pc$t <- unname(cf["t value"]); pc$p_two_sided <- unname(cf["Pr(>|t|)"])
  pc$p_one_sided_negative <- pt(cf["t value"], lf$df.residual)[[1]]
  pc$pass <- pc$estimate < 0 && pc$p_one_sided_negative < 0.05
  pc$per_gene <- m1_all[symbol %in% ctrl_in, .(symbol, logFC, SE, t, P.Value)]
} else pc$note <- "fewer than 2 control genes pass the expression filter"
tick("positive control pass:", pc$pass)
fwrite(if (!is.null(pc$per_gene)) pc$per_gene else data.table(), file.path(out, "positive_control_per_gene.tsv"), sep = "\t")

# --------------------------------------------------------------- permutations
idx_by_ds <- split(seq_len(nrow(inc)), inc$dataset)
P <- matrix(NA_real_, nrow(inc), B)
for (b in seq_len(B)) {
  s <- inc$stage_ord
  for (ix in idx_by_ds) if (length(ix) > 1L) s[ix] <- s[ix][sample.int(length(ix))]
  P[, b] <- s
}
one_perm <- function(b) {
  d <- copy(inc); d[, stage_ord := P[, b]]
  des <- model.matrix(~ dataset + sex + stage_ord, data = d)
  v <- voomWithQualityWeights(m1$y, des, plot = FALSE)
  fit <- eBayes(lmFit(v, des), robust = TRUE)
  list(t = fit$t[tested, "stage_ord"], p = fit$p.value[tested, "stage_ord"])
}
t0 <- Sys.time()
perm <- mclapply(seq_len(B), one_perm, mc.cores = ncores, mc.set.seed = FALSE)
bad <- vapply(perm, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) stop(sum(bad), " permutation fits failed")
tick("permutations done:", B, " elapsed min:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
Tm <- vapply(perm, function(x) x$t, numeric(length(tested))); Pm <- vapply(perm, function(x) x$p, numeric(length(tested)))
rownames(Tm) <- rownames(Pm) <- tested
if (length(tested) == 1L) { Tm <- matrix(Tm, nrow = 1, dimnames = list(tested, NULL)); Pm <- matrix(Pm, nrow = 1, dimnames = list(tested, NULL)) }
frac_nominal <- colMeans(Pm < 0.05)
n_bh_perm <- apply(Pm, 2, function(p) sum(p.adjust(p, "BH") < 0.05))
t_obs <- setNames(m1_t$t, m1_t$symbol)[tested]
emp_p <- (1 + rowSums(abs(Tm) >= abs(t_obs))) / (B + 1)
maxT <- apply(abs(Tm), 2, max)
p_maxT <- vapply(abs(t_obs), function(x) (1 + sum(maxT >= x)) / (B + 1), numeric(1))
calib_mean <- mean(frac_nominal)
calib_pass <- calib_mean >= 0.03 && calib_mean <= 0.07
n_bh_obs <- sum(m1_t$q_bh < 0.05)
perm_summary <- list(B = B, seed = SEED, mean_fraction_nominal_p_lt_0.05 = calib_mean,
                     quantiles_fraction_nominal = as.list(quantile(frac_nominal, c(0.025, 0.5, 0.975))),
                     calibration_pass = calib_pass,
                     n_BH_lt_0.05_observed = n_bh_obs,
                     n_BH_lt_0.05_perm_distribution = as.list(table(factor(n_bh_perm, levels = 0:max(n_bh_perm, n_bh_obs)))),
                     p_count = (1 + sum(n_bh_perm >= n_bh_obs)) / (B + 1))
fwrite(data.table(perm = seq_len(B), frac_nominal = frac_nominal, n_bh = n_bh_perm), file.path(out, "permutation_summary.tsv"), sep = "\t")
fwrite(data.table(symbol = tested, t_obs = unname(t_obs), emp_p = unname(emp_p), p_maxT = unname(p_maxT)),
       file.path(out, "permutation_target_p.tsv"), sep = "\t")

m1_t <- merge(m1_t, data.table(symbol = tested, emp_p = unname(emp_p), p_maxT = unname(p_maxT)), by = "symbol", all.x = TRUE)
m1_t[, q_emp := p.adjust(emp_p, "BH")]
m1_t[, primary_basis := if (calib_pass) "moderated_t_BH" else "permutation_empirical_BH"]
m1_t[, q_primary := if (calib_pass) q_bh else q_emp]

# ------------------------------------------------------------ bulk replication
bulk <- fread(bulk_file, select = c("gene", "logFC", "P.Value", "padj", "symbol"))
bulk[, ensembl := sub("\\..*$", "", gene)]
stopifnot(!anyDuplicated(bulk$ensembl))
bt <- bulk[ensembl %in% tl$ensembl, .(ensembl, bulk_symbol = symbol, bulk_logFC = logFC, bulk_P = P.Value, bulk_genome_padj = padj)]
bt[, bulk_q_targets := p.adjust(bulk_P, "BH")]

# ------------------------------------------------------- per-target class table
hepc <- cts[cell_type == "Hepatocytes", .(symbol, umis_hep = umis, hep_detect_rate = cells_detected / cells_total)]
allc <- cts[, .(umis_all = sum(umis)), by = symbol]
pcell <- cts[, .(cell_type, symbol, mean_umi_per_cell = umis / cells_total)]
rk <- pcell[, .(hep_per_cell_rank = frank(-mean_umi_per_cell, ties.method = "min")[cell_type == "Hepatocytes"]), by = symbol]
shr <- Reduce(function(a, b) merge(a, b, by = "symbol"), list(hepc, allc, rk))
shr[, hep_share_of_umis := ifelse(umis_all > 0, umis_hep / umis_all, NA_real_)]
shr[, nonhepatocyte_dominant := !is.na(hep_share_of_umis) & hep_share_of_umis < 0.5]
cls <- merge(tl[, .(symbol, ensembl, best_mash_row_stage, n_drugs, drugs, co_targets_listed_for_same_drug)], shr, by = "symbol", all.x = TRUE)
cls <- merge(cls, m1_t[, .(symbol, logFC_per_step = logFC, SE, t, P.Value, q_bh, q_emp, p_maxT, q_primary, primary_basis, AveExpr)], by = "symbol", all.x = TRUE)
cls[, `:=`(testable = symbol %in% tested, MDE80 = 2.8 * SE)]
m2w <- dcast(m2_t, symbol ~ contrast, value.var = c("logFC", "q_bh"))
cls <- merge(cls, m2w, by = "symbol", all.x = TRUE)
cls <- merge(cls, bt, by = "ensembl", all.x = TRUE)
s3sum <- if ("same_sign_as_M1" %in% names(s3)) s3[!is.na(same_sign_as_M1), .(S3_same_sign = sum(same_sign_as_M1), S3_runs = .N), by = symbol] else data.table(symbol = character())
cls <- merge(cls, s3sum, by = "symbol", all.x = TRUE)
if (!is.null(s2_t)) cls <- merge(cls, s2_t[, .(symbol, S2_logFC_per_F = logFC, S2_q_bh = q_bh)], by = "symbol", all.x = TRUE)
cls[, bulk_same_sign := sign(logFC_per_step) == sign(bulk_logFC)]
cls[, class := fifelse(!testable, "not_detected_in_hepatocytes",
               fifelse(q_primary < 0.05,
                       fifelse(!is.na(bulk_q_targets) & bulk_q_targets < 0.05 & bulk_same_sign, "stage_change_replicated_in_bulk",
                               "stage_change_not_replicated_in_bulk"),
                       "no_change_detected"))]
cls[, `:=`(positive_control_pass = pc$pass, calibration_pass = calib_pass, underpowered = underpowered,
           direction_claim_allowed = pc$pass && !underpowered)]
setorder(cls, class, q_primary, na.last = TRUE)
fwrite(cls, file.path(out, "target_classes.tsv"), sep = "\t")

rep_tab <- cls[testable == TRUE & !is.na(bulk_logFC)]
rep_stats <- list(n_targets = nrow(rep_tab), sign_agree = sum(rep_tab$bulk_same_sign))
if (nrow(rep_tab) >= 3L) {
  rep_stats$binom_p <- binom.test(rep_stats$sign_agree, rep_stats$n_targets, 0.5)$p.value
  rep_stats$spearman <- unname(cor.test(rep_tab$logFC_per_step, rep_tab$bulk_logFC, method = "spearman", exact = FALSE)$estimate)
  rep_stats$spearman_p <- cor.test(rep_tab$logFC_per_step, rep_tab$bulk_logFC, method = "spearman", exact = FALSE)$p.value
}

# --------------------------------------------------------------- write outputs
vt <- as.data.table(t(m1$v$E[tested, , drop = FALSE]), keep.rownames = "donor")
fwrite(vt, file.path(out, "target_voom_log2cpm_by_donor.tsv"), sep = "\t")
summ <- list(n_donors = nrow(inc), n_by_stage = setNames(n_by_stage$N, as.character(n_by_stage$stage)),
             underpowered = underpowered, single_donor_datasets_dropped = dropped_single,
             n_genes_after_filter = nrow(m1$y), n_targets = nrow(tl), n_testable = length(tested), not_detected = untested,
             positive_control = pc[intersect(names(pc), c("genes_in_model", "pass", "estimate", "t", "p_one_sided_negative", "note"))],
             permutation = perm_summary, replication = rep_stats,
             class_counts = as.list(table(cls$class)),
             S2_n_donors = nrow(s2), seed = SEED)
write_json(summ, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6, null = "null", na = "null")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(cls[, .(symbol, class, testable, logFC_per_step, SE, q_bh, q_primary, bulk_logFC, bulk_q_targets, hep_share_of_umis)], nrows = 80)
tick("done")
