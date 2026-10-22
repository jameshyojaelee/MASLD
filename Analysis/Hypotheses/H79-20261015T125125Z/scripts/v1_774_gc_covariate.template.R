#!/usr/bin/env Rscript
# V1-774: glucocorticoid-score sensitivity of the adopted 774-unit bulk refit.
# See PRESPEC.md in this directory. Model steps copy fit_disease_model() in
# scripts/manuscript/resource_f_five/refit_corrected_controls.R.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(parallel)
})
set.seed(20260923)
fail <- function(...) stop(..., call. = FALSE)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- Sys.getenv("V1_OUT")
if (!nzchar(OUT) || !dir.exists(OUT)) fail("V1_OUT must be an existing directory")
RES <- Sys.getenv("V1_RES")
if (!nzchar(RES) || !dir.exists(RES)) fail("V1_RES must be an existing directory")
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
N_DRAWS <- as.integer(Sys.getenv("V1_N_DRAWS", "200"))
FIT <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
DEX <- file.path(ROOT, "Analysis/Hypotheses/V1-20261008T191919Z/results/run-20261008T193139Z/liver_dex_did.tsv")
CONC <- file.path(OUT, "inputs/catalog_r8_concordant_genes.tsv")
SURGICAL <- c("GSE162694", "GSE213621")

w <- function(x, name) {
  p <- file.path(RES, name)
  if (file.exists(p)) fail("Refusing to overwrite ", p)
  fwrite(x, p, sep = "\t")
}

# ------------------------------------------------------------------ inputs
dge <- readRDS(file.path(FIT, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(FIT, "meta_matched.rds")))
stored <- fread(file.path(FIT, "deg_results.csv"))
meta <- meta[match(colnames(dge), sample_id)]
stopifnot(identical(meta$sample_id, colnames(dge)), !anyNA(meta$inferred_sex),
          identical(stored$gene, rownames(dge)))
info <- data.frame(
  group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
  dataset = factor(dge$samples$dataset),
  inferred_sex = factor(meta$inferred_sex),
  row.names = colnames(dge)
)
stopifnot(identical(as.character(info$group_binary), meta$group_binary))
COEF <- "group_binaryDisease"

fit_model <- function(score = NULL, keep_voom = FALSE) {
  d <- info
  if (is.null(score)) {
    design <- model.matrix(~ dataset + inferred_sex + group_binary, data = d)
  } else {
    d$gc_score <- score
    design <- model.matrix(~ dataset + inferred_sex + gc_score + group_binary, data = d)
  }
  if (qr(design)$rank != ncol(design)) fail("design is rank deficient")
  v <- voomWithQualityWeights(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = COEF, number = Inf, sort.by = "none")
  out <- list(logFC = tt$logFC, padj = tt$adj.P.Val, P = tt$P.Value)
  if (!is.null(score)) {
    ts <- topTable(fit, coef = "gc_score", number = Inf, sort.by = "none")
    out$score_logFC <- ts$logFC
    out$score_padj <- ts$adj.P.Val
  }
  if (keep_voom) out$E <- v$E
  out
}
is_call <- function(f) f$padj < 0.05 & abs(f$logFC) > 0.5

# ------------------------------------------------------------- step 0
base <- fit_model(keep_voom = TRUE)
E <- base$E
base_call <- is_call(base)
repro <- data.table(
  n_units = ncol(dge), n_control = sum(info$group_binary == "Control"),
  n_disease = sum(info$group_binary == "Disease"), n_genes = nrow(dge),
  calls = sum(base_call), up = sum(base_call & base$logFC > 0),
  down = sum(base_call & base$logFC < 0),
  max_abs_logFC_delta_vs_stored = max(abs(base$logFC - stored$logFC)),
  max_abs_padj_delta_vs_stored = max(abs(base$padj - stored$padj)),
  stored_calls = stored[, sum(padj < 0.05 & abs(logFC) > 0.5)]
)
w(repro, "step0_reproduction.tsv")
if (repro$calls != 1854L || repro$up != 1333L || repro$down != 521L ||
    repro$max_abs_logFC_delta_vs_stored > 1e-6) {
  fail("Step 0 failed: adopted fit not reproduced")
}

genes <- data.table(gene = rownames(dge), symbol = stored$symbol,
                    base_id = sub("[.][0-9]+$", "", rownames(dge)),
                    AveExpr = stored$AveExpr)
genes[, par_y := grepl("_PAR_Y", gene)]
if (anyDuplicated(genes[par_y == FALSE, base_id])) fail("non-unique version-stripped IDs")

# ------------------------------------------------------------- step 1
dex <- fread(DEX)
dex_set <- dex[dex_gene == TRUE]
genes[, dex_tested := !par_y & base_id %in% dex$base]
genes[, in_set := !par_y & base_id %in% dex_set$base]
genes[, dex_logFC := dex$logFC[match(base_id, dex$base)]]
genes[par_y == TRUE, dex_logFC := NA_real_]
genes[, set_sign := ifelse(in_set, sign(dex_logFC), 0)]
core_symbols <- c("FKBP5", "TSC22D3", "PER1", "ZBTB16", "KLF9", "DUSP1", "SGK1", "CDKN1C",
                  "MT1A", "MT1E", "MT1F", "MT1G", "MT1H", "MT1M", "MT1X", "MT2A")
genes[, in_core := !par_y & symbol %in% core_symbols]
genes[, base_logFC := base$logFC][, base_padj := base$padj][, base_call := base_call]
w(genes[in_set | in_core, .(gene, symbol, in_set, set_sign, dex_logFC, in_core, AveExpr,
                            base_logFC, base_padj, base_call)], "step1_gene_sets.tsv")

# ------------------------------------------------------------- step 2
cohort <- as.character(info$dataset)
z_within <- function(idx) {
  m <- E[idx, , drop = FALSE]
  for (co in unique(cohort)) {
    j <- cohort == co
    mu <- rowMeans(m[, j, drop = FALSE])
    s <- apply(m[, j, drop = FALSE], 1, sd)
    m[, j] <- (m[, j, drop = FALSE] - mu) / s
  }
  m
}
score_of <- function(idx, signs) {
  z <- z_within(idx)
  if (any(!is.finite(z))) fail("non-finite z")
  colMeans(z * signs)
}
i_set <- which(genes$in_set)
s_set <- genes$set_sign[i_set]
i_up <- which(genes$in_set & genes$set_sign > 0)
i_core <- which(genes$in_core)
gc_score <- score_of(i_set, s_set)
gc_up <- score_of(i_up, rep(1, length(i_up)))
gc_core <- score_of(i_core, rep(1, length(i_core)))
# Not standardized, for the cohort comparison only.
raw_signed <- colMeans(E[i_set, , drop = FALSE] * s_set)
raw_up <- colMeans(E[i_up, , drop = FALSE])
raw_core <- colMeans(E[i_core, , drop = FALSE])
units <- data.table(sample_id = colnames(dge), dataset = cohort,
                    surgical = cohort %in% SURGICAL,
                    group_binary = as.character(info$group_binary),
                    inferred_sex = as.character(info$inferred_sex),
                    gc_score, gc_up, gc_core, raw_signed, raw_up, raw_core)
units[, raw_signed_centered := raw_signed - mean(raw_signed), by = dataset]
w(units, "step2_unit_scores.tsv")

q <- function(x, p) unname(quantile(x, p))
by_cohort <- units[, .(n = .N, n_control = sum(group_binary == "Control"),
  n_disease = sum(group_binary == "Disease"),
  gc_score_mean = mean(gc_score), gc_score_sd = sd(gc_score),
  gc_score_q05 = q(gc_score, .05), gc_score_median = median(gc_score), gc_score_q95 = q(gc_score, .95),
  gc_score_mean_control = mean(gc_score[group_binary == "Control"]),
  gc_score_mean_disease = mean(gc_score[group_binary == "Disease"]),
  raw_signed_mean = mean(raw_signed), raw_signed_sd = sd(raw_signed),
  raw_signed_centered_q05 = q(raw_signed_centered, .05), raw_signed_centered_q95 = q(raw_signed_centered, .95),
  raw_up_mean = mean(raw_up), raw_up_sd = sd(raw_up),
  raw_core_mean = mean(raw_core), raw_core_sd = sd(raw_core)),
  by = .(dataset, surgical)][order(dataset)]
by_cohort[, gc_score_disease_minus_control := gc_score_mean_disease - gc_score_mean_control]
w(by_cohort, "step2_score_by_cohort.tsv")
by_status <- units[, .(n = .N, gc_score_mean = mean(gc_score), gc_score_sd = sd(gc_score),
  gc_score_median = median(gc_score), gc_up_mean = mean(gc_up), gc_core_mean = mean(gc_core)),
  by = group_binary]
w(by_status, "step2_score_by_status.tsv")
score_lm <- rbindlist(lapply(c("gc_score", "gc_up", "gc_core"), function(s) {
  d <- cbind(info, y = units[[s]])
  m <- lm(y ~ dataset + inferred_sex + group_binary, data = d)
  ci <- confint(m)[COEF, ]
  co <- summary(m)$coefficients[COEF, ]
  data.table(score = s, n_units = nrow(d), estimate_disease_minus_control = co[["Estimate"]],
             ci95_low = ci[[1]], ci95_high = ci[[2]], t = co[["t value"]],
             nominal_p = co[["Pr(>|t|)"]], score_sd = sd(d$y),
             r2_full_model = summary(m)$r.squared)
}))
w(score_lm, "step2_score_disease_lm.tsv")
score_cor <- data.table(pair = c("gc_score~gc_up", "gc_score~gc_core", "gc_up~gc_core"),
  pearson = c(cor(gc_score, gc_up), cor(gc_score, gc_core), cor(gc_up, gc_core)))
w(score_cor, "step2_score_correlations.tsv")

# ------------------------------------------------------------- step 3
fisher_row <- function(universe, direction, label) {
  g <- genes[universe]
  callv <- if (direction == "up") g$base_call & g$base_logFC > 0 else g$base_call & g$base_logFC < 0
  tab <- table(factor(callv, c(TRUE, FALSE)), factor(g$in_set, c(TRUE, FALSE)))
  ft <- fisher.test(tab)
  data.table(background = label, direction = direction, n_background = nrow(g),
    n_calls = sum(callv), n_set = sum(g$in_set), calls_in_set = tab[1, 1],
    frac_calls_in_set = tab[1, 1] / sum(callv), frac_background_in_set = mean(g$in_set),
    calls_in_set_dex_up = sum(callv & g$in_set & g$set_sign > 0),
    calls_in_set_dex_down = sum(callv & g$in_set & g$set_sign < 0),
    odds_ratio = unname(ft$estimate), or_ci95_low = ft$conf.int[1], or_ci95_high = ft$conf.int[2],
    nominal_p_two_sided = ft$p.value)
}
overlap <- rbind(
  fisher_row(genes$dex_tested, "up", "common_background"),
  fisher_row(genes$dex_tested, "down", "common_background"),
  fisher_row(rep(TRUE, nrow(genes)), "up", "all_tested_24115"),
  fisher_row(rep(TRUE, nrow(genes)), "down", "all_tested_24115"))
w(overlap, "step3_overlap_fisher.tsv")

# ------------------------------------------------------------- step 4
compare <- function(f, label) {
  cl <- is_call(f)
  lost <- base_call & !cl
  data.table(model = label, calls = sum(cl), up = sum(cl & f$logFC > 0), down = sum(cl & f$logFC < 0),
    kept = sum(base_call & cl), lost = sum(lost), gained = sum(!base_call & cl),
    kept_frac_of_1854 = sum(base_call & cl) / sum(base_call),
    lost_up = sum(lost & base$logFC > 0), lost_down = sum(lost & base$logFC < 0),
    lost_padj_only = sum(lost & f$padj >= 0.05 & abs(f$logFC) > 0.5),
    lost_lfc_only = sum(lost & f$padj < 0.05 & abs(f$logFC) <= 0.5),
    lost_both = sum(lost & f$padj >= 0.05 & abs(f$logFC) <= 0.5),
    lost_sign_flipped = sum(lost & sign(f$logFC) != sign(base$logFC)),
    lost_in_set = sum(lost & genes$in_set), kept_in_set = sum(base_call & cl & genes$in_set),
    pearson_logFC_all = cor(base$logFC, f$logFC),
    spearman_logFC_all = cor(base$logFC, f$logFC, method = "spearman"),
    pearson_logFC_calls = cor(base$logFC[base_call], f$logFC[base_call]),
    spearman_logFC_calls = cor(base$logFC[base_call], f$logFC[base_call], method = "spearman"),
    median_abs_logFC_ratio_calls = median(abs(f$logFC[base_call]) / abs(base$logFC[base_call])),
    slope_logFC_calls = unname(coef(lm(f$logFC[base_call] ~ 0 + base$logFC[base_call]))),
    score_genes_bh05 = if (is.null(f$score_padj)) NA_integer_ else sum(f$score_padj < 0.05))
}
adj <- fit_model(gc_score)
adj_up <- fit_model(gc_up)
adj_core <- fit_model(gc_core)
w(rbind(compare(base, "adopted_reproduced"), compare(adj, "gc_score_signed_70"),
        compare(adj_up, "gc_up_only"), compare(adj_core, "gc_core_recalled")),
  "step4_refit_summary.tsv")
gt <- copy(genes)[, .(gene, symbol, AveExpr, in_set, set_sign, dex_logFC, in_core,
                      base_logFC, base_padj, base_call)]
gt[, adj_logFC := adj$logFC][, adj_padj := adj$padj][, adj_call := is_call(adj)]
gt[, delta_logFC := adj_logFC - base_logFC]
gt[, score_logFC := adj$score_logFC][, score_padj := adj$score_padj]
gt[, up_logFC := adj_up$logFC][, up_padj := adj_up$padj][, up_call := is_call(adj_up)]
gt[, core_logFC := adj_core$logFC][, core_padj := adj_core$padj][, core_call := is_call(adj_core)]
gt[, status := fifelse(base_call & adj_call, "kept", fifelse(base_call, "lost",
                fifelse(adj_call, "gained", "not_called")))]
fwrite(gt, file.path(RES, "step4_gene_table.tsv.gz"), sep = "\t")
w(gt[order(-abs(delta_logFC))][1:20], "step4_largest_logFC_changes_top20.tsv")
w(gt[base_call == TRUE][order(-abs(delta_logFC))][1:20], "step4_largest_logFC_changes_among_calls_top20.tsv")
w(gt[status %in% c("lost", "gained")][order(status, -abs(base_logFC))], "step4_lost_and_gained.tsv")
conc <- fread(CONC)
c26 <- merge(conc[, .(base_id = ensembl_id, catalog_symbol = symbol)],
             cbind(gt, base_id = genes$base_id, par_y = genes$par_y)[par_y == FALSE],
             by = "base_id", all.x = TRUE)
w(c26, "step4_concordant26.tsv")

# ------------------------------------------------------------- step 5
genes[, decile := cut(AveExpr, quantile(AveExpr, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
pool <- lapply(1:10, function(b) which(genes$decile == b & !genes$in_set))
set_dec <- genes$decile[i_set]
draw_sets <- lapply(seq_len(N_DRAWS), function(k) {
  idx <- integer(length(i_set))
  for (b in unique(set_dec)) {
    pos <- which(set_dec == b)
    idx[pos] <- pool[[b]][sample.int(length(pool[[b]]), length(pos))]
  }
  idx
})
saveRDS(list(draw_sets = draw_sets, signs = s_set, gene = genes$gene),
        file.path(RES, "step5_random_sets.rds"))
DRAWDIR <- file.path(RES, "step5_draws"); dir.create(DRAWDIR)
one_draw <- function(k) {
  sc <- score_of(draw_sets[[k]], s_set)
  f <- fit_model(sc)
  cl <- is_call(f)
  r <- data.table(draw = k, compare(f, "random")[, model := NULL],
             cor_with_gc_score = cor(sc, gc_score),
             score_r2_on_design = summary(lm(sc ~ dataset + inferred_sex + group_binary, data = info))$r.squared,
             lost_among_gc_lost = sum(base_call & !cl & gt$status == "lost"),
             concordant26_called = sum(cl[match(c26$gene, genes$gene)], na.rm = TRUE))
  saveRDS(list(summary = r, call = cl), file.path(DRAWDIR, sprintf("draw_%03d.rds", k)))
  r
}
res <- mclapply(seq_len(N_DRAWS), one_draw, mc.cores = NCORES, mc.preschedule = FALSE)
bad <- !vapply(res, is.data.frame, logical(1))
if (any(bad)) { print(res[bad][[1]]); fail(sum(bad), " random draws failed") }
draws <- rbindlist(res)
w(draws, "step5_random_draws.tsv")
callmat <- vapply(seq_len(N_DRAWS), function(k)
  readRDS(file.path(DRAWDIR, sprintf("draw_%03d.rds", k)))$call, logical(nrow(genes)))
gt[, random_draws_called := rowSums(callmat)]
w(gt[status %in% c("lost", "gained"), .(gene, symbol, status, base_logFC, base_padj, adj_logFC, adj_padj,
     in_set, random_draws_called, n_draws = N_DRAWS)][order(status, random_draws_called)],
  "step5_lost_gained_vs_random.tsv")
obs <- compare(adj, "gc")
sumrow <- function(v) {
  x <- draws[[v]]
  data.table(metric = v, observed_gc_score = obs[[v]], random_median = median(x),
             random_q025 = q(x, .025), random_q975 = q(x, .975), random_min = min(x), random_max = max(x),
             n_draws = length(x),
             empirical_p_ge = (1 + sum(x >= obs[[v]])) / (length(x) + 1),
             empirical_p_le = (1 + sum(x <= obs[[v]])) / (length(x) + 1))
}
w(rbindlist(lapply(c("calls", "kept", "lost", "gained", "pearson_logFC_all", "spearman_logFC_all",
                     "pearson_logFC_calls", "median_abs_logFC_ratio_calls", "score_genes_bh05"), sumrow)),
  "step5_random_summary.tsv")
w(data.table(metric = c("gc_score_r2_on_design", "random_score_r2_on_design_median",
                        "random_abs_cor_with_gc_score_median", "random_abs_cor_with_gc_score_max",
                        "concordant26_called_adopted", "concordant26_called_gc",
                        "concordant26_called_random_min"),
             value = c(score_lm[score == "gc_score", r2_full_model], median(draws$score_r2_on_design),
                       median(abs(draws$cor_with_gc_score)), max(abs(draws$cor_with_gc_score)),
                       sum(c26$base_call, na.rm = TRUE), sum(c26$adj_call, na.rm = TRUE),
                       min(draws$concordant26_called))),
  "step5_extra.tsv")

writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo.txt"))
writeLines("COMPLETE", file.path(RES, "COMPLETE"))
