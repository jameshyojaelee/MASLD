#!/usr/bin/env Rscript
# N3 ambient check, step 4: per-gene ambient readouts, old-versus-new stage association,
# ambient load versus stage, marker correlation, conclusion per PRESPEC.md.
# Usage: 04_ambient_summary.R <RESULT_DIR> <N3_DIR>
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR); library(jsonlite) })
argv <- commandArgs(trailingOnly = TRUE); stopifnot(length(argv) == 2L)
res <- argv[[1]]; n3 <- argv[[2]]
out <- file.path(res, "ambient")
if (dir.exists(out)) stop("Refusing to overwrite ", out)
dir.create(out)
set.seed(20261013L)

genes <- readLines(file.path(res, "genes.txt"))
report_genes <- readLines(file.path(res, "report_genes.txt"))
tl <- fread(file.path(res, "target_list.tsv"))
dm <- fread(file.path(res, "donor_meta.tsv"))
da <- fread(file.path(res, "donor_ambient.tsv"))
qc <- fread(file.path(res, "library_ambient_qc.tsv"))
dl <- fread(file.path(res, "donor_gene_celltype_counts.tsv.gz"))
cts <- fread(file.path(res, "target_celltype_summary.tsv"))
mats <- readRDS(file.path(res, "donor_matrices.rds"))
old <- fread(file.path(n3, "results", "target_classes.tsv"))
new <- fread(file.path(res, "results_corrected", "target_classes.tsv"))
historical <- fread(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N3-20261009T140326Z/results/target_classes.tsv"))
NINE <- historical[q_bh < 0.05][order(q_bh), symbol] # fixed historical diagnostic roster
stopifnot(length(NINE) == 9L)
HEP_GENES <- c("AR", "NR1H4", "NR3C2", "LEPR", "KLB", "GCGR", "ALB", "APOA1")
GROUP3 <- c("FGF1", "PDE4A", "SLC5A1", "VDR", "CNR1")
IMMUNE <- c("PTPRC", "CD3E", "CD2", "CD247", "MS4A1", "CD74")
STROMAL <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM")

# ------------------------------------------------------------- donors as in N3
inc <- dm[included == TRUE]
inc[, stage := factor(disease_stage_coarse, levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]
inc[, sex := factor(sex_call, levels = c("Female", "Male"))]
inc[, dataset := factor(dataset)]
inc <- inc[order(dataset, stage, donor)]
stopifnot(nrow(inc) == 70L)
inc <- merge(inc, da[, .(donor, hep_contamination, frac_umis_outside_called_cells, n_passthrough, n_no_empty_profile)],
             by = "donor", sort = FALSE)
inc[, logit_contam := qlogis(pmin(pmax(hep_contamination, 1e-4), 1 - 1e-4))]

# ---------------------------------------------- M1 refit (same code as 03_donor_model.R) with CI
fit_voom <- function(counts, design, genes = NULL) {
  y <- DGEList(counts)
  keep <- if (is.null(genes)) filterByExpr(y, design, min.count = 10, min.total.count = 15) else rownames(y) %in% genes
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voomWithQualityWeights(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  list(y = y, v = v, fit = fit)
}
ci_table <- function(fit, coef) {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = TRUE)
  data.table(symbol = rownames(tt), logFC = tt$logFC, CI_L = tt$CI.L, CI_R = tt$CI.R, t = tt$t, P.Value = tt$P.Value)
}
load_counts <- function(path) {
  cnt <- fread(path); m <- as.matrix(cnt[, inc$donor, with = FALSE]); rownames(m) <- cnt$symbol
  storage.mode(m) <- "double"; m
}
mat_raw <- load_counts(file.path(n3, "hep_donor_counts.tsv.gz"))
mat_dec <- load_counts(file.path(res, "hep_donor_counts.tsv.gz"))
design1 <- model.matrix(~ dataset + sex + stage_ord, data = inc)
m_raw <- fit_voom(mat_raw, design1); m_dec <- fit_voom(mat_dec, design1)
ci_raw <- ci_table(m_raw$fit, "stage_ord")[symbol %in% tl$symbol]
ci_dec <- ci_table(m_dec$fit, "stage_ord")[symbol %in% tl$symbol]
chk <- function(ci, cls, label) {
  x <- merge(ci, cls[testable == TRUE, .(symbol, logFC_per_step, P.Value)], by = "symbol")
  stopifnot(nrow(x) == sum(cls$testable), nrow(x) == nrow(ci))
  d <- max(abs(x$logFC - x$logFC_per_step)); dp <- max(abs(x$P.Value.x - x$P.Value.y))
  cat(label, "refit vs 03_donor_model.R: max |dlogFC| =", d, " max |dP| =", dp, "\n")
  if (d > 1e-8 || dp > 1e-8) stop("refit does not reproduce 03_donor_model.R for ", label)
  c(max_dlogFC = d, max_dP = dp)
}
repro <- list(raw = chk(ci_raw, old, "raw"), corrected = chk(ci_dec, new, "corrected"))

# --------------------------------------- sensitivity: contamination fraction as covariate
design_s <- model.matrix(~ dataset + sex + stage_ord + logit_contam, data = inc)
stopifnot(is.null(nonEstimable(design_s)))
m_s <- fit_voom(mat_dec, design_s, genes = rownames(m_dec$y))
ci_s <- ci_table(m_s$fit, "stage_ord")[symbol %in% tl$symbol]
ci_s[, q_bh := p.adjust(P.Value, "BH")]

# ---------------------------------------------------- (a) ambient-attributed fraction
hr <- mats$hep_raw[report_genes, inc$donor, drop = FALSE]
hd <- mats$hep_dec[report_genes, inc$donor, drop = FALSE]
frac <- 1 - hd / hr; frac[hr == 0] <- NA
qf <- function(x, p) if (all(is.na(x))) NA_real_ else unname(quantile(x, p, na.rm = TRUE))
amb <- data.table(symbol = report_genes,
                  n_donors_raw_gt0 = rowSums(hr > 0),
                  hep_raw_umis_70 = rowSums(hr), hep_corrected_umis_70 = rowSums(hd),
                  ambient_frac_pooled = 1 - rowSums(hd) / rowSums(hr),
                  ambient_frac_median = apply(frac, 1, qf, 0.5),
                  ambient_frac_q25 = apply(frac, 1, qf, 0.25), ambient_frac_q75 = apply(frac, 1, qf, 0.75))
fr_long <- as.data.table(as.table(frac)); setnames(fr_long, c("symbol", "donor", "ambient_frac"))
fr_long[, raw := as.vector(hr)]; fr_long[, corrected := as.vector(hd)]
fwrite(fr_long[symbol %in% tl$symbol], file.path(out, "ambient_fraction_by_donor_target.tsv.gz"), sep = "\t")

# ---------------------------------------------------- (b) position in the empty-droplet profile
libs <- qc[profile_source == "empty_droplets", library]
ri <- match(report_genes, genes)
cpm <- matrix(NA_real_, length(report_genes), length(libs), dimnames = list(report_genes, libs)); rk <- cpm
for (l in libs) {
  b <- readRDS(file.path(res, "per_library", paste0(l, ".rds")))$background
  cpm[, l] <- b[ri] / sum(b) * 1e6
  rk[, l] <- frank(-b, ties.method = "min")[ri]
}
bgp <- rowSums(mats$background)
rk_pooled <- frank(-bgp, ties.method = "min")[ri]
soup <- data.table(symbol = report_genes, soup_cpm_median = apply(cpm, 1, median), soup_cpm_pooled = bgp[ri] / sum(bgp) * 1e6,
                   soup_rank_median = apply(rk, 1, median), soup_rank_pooled = rk_pooled,
                   n_libraries_zero_in_soup = rowSums(cpm == 0), n_libraries = length(libs))
soup[, `:=`(soup_top100 = soup_rank_median <= 100, soup_top1000 = soup_rank_median <= 1000, soup_top5000 = soup_rank_median <= 5000)]

# ---------------------------------------------------- (c) model-free ratio and expected ambient count
hep_tot <- colSums(mats$hep_raw[, inc$donor]); bg <- mats$background[, inc$donor]; bg_tot <- colSums(bg)
eta <- sweep(bg[report_genes, , drop = FALSE], 2, bg_tot, "/")
share_hep <- sweep(hr, 2, hep_tot, "/")
R <- share_hep / eta; R[eta == 0] <- NA
expct <- sweep(eta, 2, inc$hep_contamination * hep_tot, "*")
oe <- hr / expct; oe[expct == 0] <- NA
mf <- data.table(symbol = report_genes, R_median = apply(R, 1, qf, 0.5), R_q25 = apply(R, 1, qf, 0.25), R_q75 = apply(R, 1, qf, 0.75),
                 n_donors_R_gt1 = rowSums(R > 1, na.rm = TRUE), n_donors_soup_zero = rowSums(eta == 0),
                 obs_over_expected_ambient_median = apply(oe, 1, qf, 0.5),
                 obs_over_expected_ambient_q25 = apply(oe, 1, qf, 0.25), obs_over_expected_ambient_q75 = apply(oe, 1, qf, 0.75),
                 observed_70 = rowSums(hr), expected_ambient_70 = rowSums(expct))

# ---------------------------------------------------- (d) negative-control cell type
rate <- dl[, .(donor, cell_type, symbol, n_cells, raw_rate = raw / ct_raw_umis, dec_rate = fifelse(ct_dec_umis > 0, dec / ct_dec_umis, NA_real_), raw, dec)]
hep <- rate[cell_type == "Hepatocytes", .(donor, symbol, hep_raw_rate = raw_rate, hep_dec_rate = dec_rate, hep_raw = raw, hep_dec = dec)]
tc <- merge(rate[cell_type == "T cells" & n_cells >= 20 & symbol %in% HEP_GENES], hep, by = c("donor", "symbol"))
neg_hepgenes <- tc[, .(negative_control_cell_type = "T cells", n_donors = .N,
                       tcell_raw_umis = sum(raw), tcell_corrected_umis = sum(dec),
                       tcell_ambient_frac_pooled = 1 - sum(dec) / sum(raw),
                       raw_rate_ratio_T_over_hep_median = median(raw_rate / hep_raw_rate, na.rm = TRUE),
                       corrected_rate_ratio_T_over_hep_median = median(dec_rate / hep_dec_rate, na.rm = TRUE)), by = symbol]
fwrite(neg_hepgenes, file.path(out, "negative_control_Tcells_hepatocyte_genes.tsv"), sep = "\t")
pc <- cts[cell_type != "Hepatocytes", .(cell_type, symbol, m = umis / cells_total)]
top <- pc[symbol %in% GROUP3, .SD[which.max(m)], by = symbol][, .(symbol, top_cell_type = cell_type)]
g3 <- merge(merge(rate[symbol %in% GROUP3], top, by = "symbol")[cell_type == top_cell_type & n_cells >= 20], hep, by = c("donor", "symbol"))
neg_group3 <- g3[, .(n_donors = .N, top_type_raw_umis = sum(raw), hep_raw_umis = sum(hep_raw), hep_corrected_umis = sum(hep_dec),
                     raw_rate_ratio_hep_over_top_median = median(hep_raw_rate / raw_rate, na.rm = TRUE),
                     corrected_rate_ratio_hep_over_top_median = median(hep_dec_rate / dec_rate, na.rm = TRUE)), by = .(symbol, top_cell_type)]
fwrite(neg_group3, file.path(out, "group3_hepatocyte_vs_top_cell_type.tsv"), sep = "\t")

# ---------------------------------------------------- ambient load versus stage
ols <- function(y, label) {
  f <- lm(y ~ dataset + sex + stage_ord, data = cbind(inc, y = y))
  cf <- summary(f)$coefficients["stage_ord", ]; ci <- confint(f)["stage_ord", ]
  data.table(outcome = label, n_donors = nrow(inc), estimate_per_step = cf[["Estimate"]], SE = cf[["Std. Error"]],
             CI_L = ci[[1]], CI_R = ci[[2]], t = cf[["t value"]], p = cf[["Pr(>|t|)"]])
}
load_tests <- rbind(
  ols(inc$logit_contam, "primary: logit hepatocyte contamination fraction (decontX)"),
  ols(qlogis(inc$frac_umis_outside_called_cells), "secondary: logit fraction of raw-matrix UMIs outside Cell Ranger-called cells"))
soup_cpm_d <- sweep(bg[NINE, , drop = FALSE], 2, bg_tot, "/") * 1e6
soup_tests <- rbindlist(lapply(NINE, function(g) ols(log2(soup_cpm_d[g, ] + 0.5), paste0("soup log2(CPM+0.5): ", g))))
soup_tests[, q_bh_9 := p.adjust(p, "BH")]
fwrite(rbind(load_tests, soup_tests, fill = TRUE), file.path(out, "ambient_load_vs_stage.tsv"), sep = "\t")
fwrite(inc[, .(donor, dataset, stage, stage_ord, sex, hep_contamination, logit_contam, frac_umis_outside_called_cells,
               n_passthrough, n_no_empty_profile)], file.path(out, "donor_ambient_in_model.tsv"), sep = "\t")
fwrite(inc[, .(n = .N, hep_contamination_median = median(hep_contamination), q25 = quantile(hep_contamination, 0.25),
               q75 = quantile(hep_contamination, 0.75)), by = .(dataset, stage)][order(dataset, stage)],
       file.path(out, "hep_contamination_by_dataset_stage.tsv"), sep = "\t")

# ---------------------------------------------------- marker correlation
marker <- function(m, label) {
  E <- m$v$E; dd <- as.data.frame(inc[, .(dataset, sex, stage_ord)])
  comp <- function(gs) { g <- intersect(gs, rownames(E)); list(score = colMeans(t(scale(t(E[g, , drop = FALSE])))), genes = g) }
  im <- comp(IMMUNE); stc <- comp(STROMAL)
  rz <- function(y) resid(lm(y ~ dataset + sex, data = dd))
  rows <- list()
  for (g in intersect(NINE, rownames(E))) {
    for (cn in c("immune", "stromal")) {
      s <- if (cn == "immune") im$score else stc$score
      ct <- cor.test(rz(E[g, ]), rz(s))
      rows[[length(rows) + 1L]] <- data.table(matrix = label, symbol = g, composite = cn, r = unname(ct$estimate), p = ct$p.value)
    }
  }
  r <- rbindlist(rows); r[, q_bh_18 := p.adjust(p, "BH")]
  adj <- rbindlist(lapply(intersect(NINE, rownames(E)), function(g) {
    d <- cbind(dd, y = E[g, ], immune = im$score, stromal = stc$score)
    a <- summary(lm(y ~ dataset + sex + stage_ord, data = d))$coefficients["stage_ord", ]
    b <- summary(lm(y ~ dataset + sex + stage_ord + immune + stromal, data = d))$coefficients["stage_ord", ]
    data.table(matrix = label, symbol = g, ols_stage = a[["Estimate"]], ols_stage_se = a[["Std. Error"]],
               ols_stage_adj_markers = b[["Estimate"]], ols_stage_adj_markers_se = b[["Std. Error"]], ols_stage_adj_markers_p = b[["Pr(>|t|)"]])
  }))
  list(r = r, adj = adj, immune_genes = im$genes, stromal_genes = stc$genes)
}
mk_raw <- marker(m_raw, "uncorrected"); mk_dec <- marker(m_dec, "corrected")
fwrite(rbind(mk_raw$r, mk_dec$r), file.path(out, "marker_correlation.tsv"), sep = "\t")
fwrite(rbind(mk_raw$adj, mk_dec$adj), file.path(out, "marker_adjusted_stage_ols.tsv"), sep = "\t")

# ---------------------------------------------------- old versus new and conclusion
o <- old[, .(symbol, N3_class = class, hep_share_of_umis, nonhepatocyte_dominant, hep_detect_rate,
             old_testable = testable, old_logFC = logFC_per_step, old_SE = SE, old_q_bh = q_bh, old_q_emp = q_emp, old_p_maxT = p_maxT)]
n <- new[, .(symbol, new_testable = testable, new_logFC = logFC_per_step, new_SE = SE, new_P = P.Value, new_q_bh = q_bh,
             new_q_emp = q_emp, new_p_maxT = p_maxT, new_class = class)]
tab <- merge(o, n, by = "symbol")
tab <- merge(tab, ci_raw[, .(symbol, old_CI_L = CI_L, old_CI_R = CI_R)], by = "symbol", all.x = TRUE)
tab <- merge(tab, ci_dec[, .(symbol, new_CI_L = CI_L, new_CI_R = CI_R)], by = "symbol", all.x = TRUE)
tab <- merge(tab, ci_s[, .(symbol, sens_logFC = logFC, sens_CI_L = CI_L, sens_CI_R = CI_R, sens_q_bh = q_bh)], by = "symbol", all.x = TRUE)
# BH restricted to the targets testable in both fits (reported next to the family of the corrected fit)
both <- tab[old_testable & new_testable, symbol]
tab[symbol %in% both, new_q_bh_common_set := p.adjust(new_P, "BH")]
tab <- Reduce(function(a, b) merge(a, b, by = "symbol", all.x = TRUE), list(tab, amb, soup, mf))
tab[, crit_i := !is.na(ambient_frac_median) & ambient_frac_median >= 0.5]
tab[, crit_ii := old_testable & !new_testable]
tab[, crit_iii := !is.na(obs_over_expected_ambient_median) & obs_over_expected_ambient_median <= 2]
concl <- function(i) with(tab[i], {
  if (crit_ii) return(c("not interpretable because of ambient RNA", "passes the N3 expression filter, fails it after correction"))
  if (nonhepatocyte_dominant && old_testable && (crit_i || crit_iii))
    return(c("not interpretable because of ambient RNA",
             paste(c(if (crit_i) "median ambient-attributed fraction >= 0.5", if (crit_iii) "observed <= 2 x expected ambient count"), collapse = "; ")))
  if (!old_testable) return(if (!new_testable) c("unchanged", "still fails the expression filter") else c("changed", "passes the expression filter after correction"))
  if (old_q_bh < 0.05) {
    ok <- sign(new_logFC) == sign(old_logFC) && new_q_bh < 0.05 && abs(new_logFC) >= 0.75 * abs(old_logFC)
    return(if (ok) c("unchanged", "same sign, corrected BH q < 0.05, |effect| >= 0.75 x N3")
           else c("weakened", paste(c(if (sign(new_logFC) != sign(old_logFC)) "sign differs", if (new_q_bh >= 0.05) "corrected BH q >= 0.05",
                                      if (abs(new_logFC) < 0.75 * abs(old_logFC)) "|effect| < 0.75 x N3"), collapse = "; ")))
  }
  if (new_q_bh >= 0.05) c("unchanged", "corrected BH q >= 0.05") else c("stage trend appears after correction", "corrected BH q < 0.05; N3 q >= 0.05")
})
cc <- t(vapply(seq_len(nrow(tab)), concl, character(2)))
tab[, `:=`(conclusion = cc[, 1], conclusion_reason = cc[, 2])]
setorder(tab, old_q_bh, na.last = TRUE)
fwrite(tab, file.path(out, "old_vs_new_all_targets.tsv"), sep = "\t")
fwrite(Reduce(function(a, b) merge(a, b, by = "symbol"), list(amb, soup, mf)), file.path(out, "ambient_readouts_report_genes.tsv"), sep = "\t")

summ <- list(
  n_donors = nrow(inc), n_libraries = nrow(qc), libraries_corrected = sum(qc$status == "corrected"),
  libraries_empty_droplet_profile = sum(qc$profile_source == "empty_droplets"),
  donors_with_passthrough_library = sum(inc$n_passthrough > 0), donors_with_library_without_empty_profile = sum(inc$n_no_empty_profile > 0),
  hep_contamination_donor_quantiles = as.list(quantile(inc$hep_contamination, c(0, 0.25, 0.5, 0.75, 1))),
  genes_after_filter = list(raw = nrow(m_raw$y), corrected = nrow(m_dec$y)),
  testable = list(old = sum(tab$old_testable), new = sum(tab$new_testable),
                  lost = tab[old_testable & !new_testable, symbol], gained = tab[!old_testable & new_testable, symbol]),
  n_q_lt_0.05 = list(old = sum(tab$old_q_bh < 0.05, na.rm = TRUE), new = sum(tab$new_q_bh < 0.05, na.rm = TRUE),
                     sensitivity_with_contamination_covariate = sum(tab$sens_q_bh < 0.05, na.rm = TRUE)),
  refit_reproduction = repro, conclusion_counts = as.list(table(tab$conclusion)),
  marker_genes = list(immune_uncorrected = mk_raw$immune_genes, stromal_uncorrected = mk_raw$stromal_genes,
                      immune_corrected = mk_dec$immune_genes, stromal_corrected = mk_dec$stromal_genes),
  seed_decontx = 20261013L, seed_permutation = 20260923L)
write_json(summ, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6, null = "null", na = "null")
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_summary.txt"))
print(tab[, .(symbol, N3_class, old_logFC, old_q_bh, new_logFC, new_q_bh, ambient_frac_median, obs_over_expected_ambient_median, conclusion)], nrows = 60)
print(load_tests)
cat("done\n")
