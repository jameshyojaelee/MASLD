#!/usr/bin/env Rscript
# V1: perioperative dexamethasone (DEX) response in GSE320593 versus the
# five-cohort bulk disease signature. Follows PRESPEC_V1.md exactly.
# Usage: Rscript --vanilla v1_dex_overlap.R <V1_dir> <out_dir> <sample_design.tsv>

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
v1_dir <- args[1]
out <- args[2]
design_path <- args[3]
stopifnot(dir.exists(out), file.exists(file.path(v1_dir, "PRESPEC_V1.md")))

repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
dl <- file.path(v1_dir, "downloads_untrusted")
refit_path <- file.path(repo, "RNA-seq/results/manuscript_release/candidates",
  "bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
hot <- file.path(repo, "Analysis/Multimodal_Program_Projection/candidates",
  "program-context-v2-candidate-2026-08-07/hotspot")
gencode_path <- file.path(repo, "RNA-seq/results/manuscript_release/candidates",
  "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs",
  "gencode_v49_gene_metadata.tsv.gz")

SEED <- 20260923
N_NULL <- 10000
N_BINS <- 20
PADJ <- 0.05
LFC <- 0.5
MIN_PROGRAM_GENES <- 8
GR_TARGETS <- c("FKBP5", "TSC22D3", "ZBTB16", "PER1")

summary_rows <- list()
put <- function(key, value) {
  summary_rows[[length(summary_rows) + 1]] <<- data.table(item = key, value = as.character(value))
  cat(key, "=", as.character(value), "\n")
}
cp_ci <- function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  binom.test(x, n)$conf.int[1:2]
}
spearman_ci <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  se <- 1 / sqrt(n - 3)
  c(rho = rho, lo = tanh(atanh(rho) - 1.96 * se), hi = tanh(atanh(rho) + 1.96 * se),
    p = ct$p.value, n = n)
}
fmt <- function(x, d = 4) formatC(x, digits = d, format = "g")

# ---------------------------------------------------------------- gene maps
gc <- fread(gencode_path)
primary_chr <- paste0("chr", c(1:22, "X", "Y", "M"))
sym_primary <- gc[chromosome %in% primary_chr,
  .(n = uniqueN(ensembl_base), base = ensembl_base[1]), by = gene_name][n == 1]
sym_all <- gc[, .(n = uniqueN(ensembl_base), base = ensembl_base[1]), by = gene_name][n == 1]

read_counts <- function(path) {
  x <- fread(path, header = TRUE)
  setnames(x, 1, "symbol")
  m <- as.matrix(x[, -1, with = FALSE])
  rownames(m) <- x$symbol
  storage.mode(m) <- "double"
  stopifnot(!anyDuplicated(rownames(m)), all(m >= 0), all(m == round(m)))
  m
}
map_counts <- function(m, tag) {
  base <- sym_primary$base[match(rownames(m), sym_primary$gene_name)]
  keep <- !is.na(base)
  put(paste0(tag, "_symbols_total"), nrow(m))
  put(paste0(tag, "_symbols_mapped"), sum(keep))
  put(paste0(tag, "_symbols_mapped_pct"), round(100 * mean(keep), 2))
  fwrite(data.table(symbol = rownames(m)[!keep]),
    file.path(out, paste0(tag, "_unmapped_symbols.tsv")), sep = "\t")
  m2 <- m[keep, , drop = FALSE]
  sym <- rownames(m2)
  rownames(m2) <- base[keep]
  stopifnot(!anyDuplicated(rownames(m2)))
  list(m = m2, sym = setNames(sym, rownames(m2)), pct = 100 * mean(keep))
}

# ---------------------------------------------------------------- DEX fit
fit_dex <- function(counts, meta, alt_batch = NULL) {
  meta <- as.data.frame(meta)
  meta$cell <- factor(paste(meta$procedure, meta$arm, sep = "_"),
    levels = c("pre_DEX", "surg_DEX", "pre_NODEX", "surg_NODEX"))
  stopifnot(!anyNA(meta$cell))
  if (!is.null(alt_batch)) {
    meta$batch_alt <- as.integer(meta$batch == alt_batch)
    design <- model.matrix(~ 0 + cell + batch_alt, data = meta)
  } else {
    design <- model.matrix(~ 0 + cell, data = meta)
  }
  colnames(design) <- sub("^cell", "", colnames(design))
  stopifnot(qr(design)$rank == ncol(design))
  y <- DGEList(counts)
  keep <- filterByExpr(y, group = meta$cell)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design)
  dc <- duplicateCorrelation(v, design, block = meta$donor)
  v <- voom(y, design, block = meta$donor, correlation = dc$consensus.correlation)
  dc <- duplicateCorrelation(v, design, block = meta$donor)
  fit <- lmFit(v, design, block = meta$donor, correlation = dc$consensus.correlation)
  cm <- makeContrasts(
    dex_did = (surg_DEX - pre_DEX) - (surg_NODEX - pre_NODEX),
    surg_unpaired = surg_DEX - surg_NODEX,
    pre_arm_balance = pre_DEX - pre_NODEX,
    surgery_nodex = surg_NODEX - pre_NODEX,
    levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  tabs <- lapply(colnames(cm), function(k) {
    tt <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
    data.table(base = rownames(tt), contrast = k, logFC = tt$logFC, AveExpr = tt$AveExpr,
      t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
  })
  tabs <- rbindlist(tabs)
  tabs[, pass := adj.P.Val < PADJ & abs(logFC) > LFC]
  list(tab = tabs, cor = dc$consensus.correlation, n_tested = sum(keep),
    n_samples = ncol(counts), n_donors = length(unique(meta$donor)),
    n_donors_dex = length(unique(meta$donor[meta$arm == "DEX"])),
    n_donors_nodex = length(unique(meta$donor[meta$arm == "NODEX"])))
}

# ---------------------------------------------------------------- design
meta_all <- fread(design_path)
put("design_samples_total", nrow(meta_all))

# ---------------------------------------------------------------- liver
liver_raw <- read_counts(file.path(dl, "GSE320593_liver.counts.csv.gz"))
liver_map <- map_counts(liver_raw, "liver")
meta_l <- meta_all[tissue == "liver"]
stopifnot(setequal(meta_l$library, colnames(liver_map$m)))
meta_l <- meta_l[match(colnames(liver_map$m), library)]

prim <- fit_dex(liver_map$m, meta_l, alt_batch = "D")
put("liver_samples", prim$n_samples)
put("liver_donors", prim$n_donors)
put("liver_donors_DEX_arm", prim$n_donors_dex)
put("liver_donors_NODEX_arm", prim$n_donors_nodex)
put("liver_genes_tested", prim$n_tested)
put("liver_donor_consensus_correlation", fmt(prim$cor))
liver_tab <- prim$tab
liver_tab[, symbol := liver_map$sym[base]]
fwrite(liver_tab, file.path(out, "liver_contrasts_all_genes.tsv.gz"), sep = "\t")

for (k in unique(liver_tab$contrast)) {
  put(paste0("liver_", k, "_n_pass"), liver_tab[contrast == k, sum(pass)])
  put(paste0("liver_", k, "_n_up"), liver_tab[contrast == k, sum(pass & logFC > 0)])
  put(paste0("liver_", k, "_n_down"), liver_tab[contrast == k, sum(pass & logFC < 0)])
}

dex <- liver_tab[contrast == "dex_did"]
dex[, dex_gene := pass]
fwrite(dex[order(P.Value)], file.path(out, "liver_dex_did.tsv"), sep = "\t")
n_dex <- sum(dex$dex_gene)

# Positive control
pc <- dex[match(sym_primary$base[match(GR_TARGETS, sym_primary$gene_name)], base)]
pc[, gr_target := GR_TARGETS]
fwrite(pc[, .(gr_target, base, logFC, adj.P.Val, dex_gene)],
  file.path(out, "positive_control_gr_targets.tsv"), sep = "\t")
n_pc_pos <- sum(pc$logFC > 0, na.rm = TRUE)
put("positive_control_n_positive_of_4", n_pc_pos)
put("positive_control_n_measured", sum(!is.na(pc$logFC)))
s1_pass <- n_pc_pos >= 3
put("S1_positive_control_pass", s1_pass)
put("flag_liver_symbol_map_below_85pct", liver_map$pct < 85)
put("flag_arm_balance_over_50", liver_tab[contrast == "pre_arm_balance", sum(pass)] > 50)
put("S3_dex_genes_below_20", n_dex < 20)

# Sensitivity 1: batch-C-only donors
donors_c <- meta_l[, .(allC = all(batch == "C")), by = donor][allC == TRUE, donor]
meta_c <- meta_l[donor %in% donors_c]
sens_c <- fit_dex(liver_map$m[, meta_c$library], meta_c, alt_batch = NULL)
sc <- sens_c$tab[contrast == "dex_did"]
fwrite(sc, file.path(out, "liver_sens_batchC_dex_did.tsv.gz"), sep = "\t")
put("sensC_donors", sens_c$n_donors)
put("sensC_donors_DEX_arm", sens_c$n_donors_dex)
put("sensC_donors_NODEX_arm", sens_c$n_donors_nodex)
put("sensC_genes_tested", sens_c$n_tested)
put("sensC_n_pass", sum(sc$pass))
cmpC <- merge(dex[, .(base, l1 = logFC, dex_gene)], sc[, .(base, l2 = logFC, pass2 = pass)], by = "base")
put("sensC_primary_dex_genes_tested", cmpC[dex_gene == TRUE, .N])
put("sensC_primary_dex_genes_same_sign", cmpC[dex_gene == TRUE, sum(sign(l1) == sign(l2))])
put("sensC_primary_dex_genes_also_pass", cmpC[dex_gene == TRUE, sum(pass2)])
rC <- spearman_ci(cmpC$l1, cmpC$l2)
put("sensC_spearman_logFC_all_common", fmt(rC["rho"]))

# Sensitivity 2: unpaired surgical contrast (same primary fit)
su <- liver_tab[contrast == "surg_unpaired"]
cmpU <- merge(dex[, .(base, l1 = logFC, dex_gene)], su[, .(base, l2 = logFC, pass2 = pass)], by = "base")
put("unpaired_primary_dex_genes_same_sign", cmpU[dex_gene == TRUE, sum(sign(l1) == sign(l2))])
put("unpaired_primary_dex_genes_also_pass", cmpU[dex_gene == TRUE, sum(pass2)])

# ---------------------------------------------------------------- adipose (secondary)
at_raw <- read_counts(file.path(dl, "GSE320593_AT.counts.csv.gz"))
at_map <- map_counts(at_raw, "adipose")
meta_a <- meta_all[tissue == "adipose"]
stopifnot(setequal(meta_a$library, colnames(at_map$m)))
meta_a <- meta_a[match(colnames(at_map$m), library)]
adip <- fit_dex(at_map$m, meta_a, alt_batch = "B")
at <- adip$tab[contrast == "dex_did"]
at[, symbol := at_map$sym[base]]
fwrite(adip$tab, file.path(out, "adipose_contrasts_all_genes.tsv.gz"), sep = "\t")
put("adipose_samples", adip$n_samples)
put("adipose_donors", adip$n_donors)
put("adipose_donors_DEX_arm", adip$n_donors_dex)
put("adipose_donors_NODEX_arm", adip$n_donors_nodex)
put("adipose_genes_tested", adip$n_tested)
put("adipose_donor_consensus_correlation", fmt(adip$cor))
put("adipose_dex_did_n_pass", sum(at$pass))
put("adipose_pre_arm_balance_n_pass", adip$tab[contrast == "pre_arm_balance", sum(pass)])
cmpA <- merge(dex[, .(base, l1 = logFC, dex_gene)], at[, .(base, l2 = logFC, pass2 = pass)], by = "base")
nA <- cmpA[dex_gene == TRUE, .N]
sA <- cmpA[dex_gene == TRUE, sum(sign(l1) == sign(l2))]
ciA <- cp_ci(sA, nA)
put("adipose_liver_dex_genes_tested", nA)
put("adipose_liver_dex_genes_same_sign", sA)
put("adipose_liver_dex_genes_same_sign_frac", fmt(sA / nA))
put("adipose_liver_dex_genes_same_sign_ci95", paste(fmt(ciA), collapse = ";"))
put("adipose_liver_dex_genes_also_adipose_dex", cmpA[dex_gene == TRUE, sum(pass2)])
rA <- spearman_ci(cmpA$l1, cmpA$l2)
put("adipose_vs_liver_spearman_all_common", fmt(rA["rho"]))
put("adipose_vs_liver_spearman_all_common_n", rA["n"])

if (!s1_pass) {
  put("STOP", "S1 positive control failed; no overlap analysis")
  fwrite(rbindlist(summary_rows), file.path(out, "summary.tsv"), sep = "\t")
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  quit(status = 0)
}

# ---------------------------------------------------------------- refit
ref <- fread(refit_path)
ref <- ref[!grepl("_PAR_Y", gene)]
ref[, base := sub("\\.[0-9]+$", "", gene)]
if (anyDuplicated(ref$base)) stop("S2: version-stripped refit IDs not unique")
ref[, call := !is.na(padj) & padj < PADJ & abs(logFC) > LFC]
n_calls <- sum(ref$call)
put("refit_genes_tested", nrow(ref))
put("refit_calls", n_calls)
if (n_calls != 1854) stop("S2: refit call count is not 1,854")

bg <- merge(
  ref[!is.na(logFC), .(base, symbol, dis_logFC = logFC, dis_padj = padj, dis_AveExpr = AveExpr, call)],
  dex[, .(base, dex_logFC = logFC, dex_padj = adj.P.Val, dex_AveExpr = AveExpr, dex_gene)],
  by = "base")
n_bg <- nrow(bg)
n_calls_bg <- sum(bg$call)
n_dex_bg <- sum(bg$dex_gene)
ov <- bg[call == TRUE & dex_gene == TRUE]
n_ov <- nrow(ov)
put("common_background_genes", n_bg)
put("calls_in_background", n_calls_bg)
put("calls_not_testable_in_dex_data", n_calls - n_calls_bg)
put("dex_genes_in_background", n_dex_bg)
put("overlap_calls_that_are_dex_genes", n_ov)
ci_all <- cp_ci(n_ov, n_calls)
ci_bg <- cp_ci(n_ov, n_calls_bg)
put("overlap_frac_of_1854", fmt(n_ov / n_calls))
put("overlap_frac_of_1854_ci95", paste(fmt(ci_all), collapse = ";"))
put("overlap_frac_of_calls_in_bg", fmt(n_ov / n_calls_bg))
put("overlap_frac_of_calls_in_bg_ci95", paste(fmt(ci_bg), collapse = ";"))
put("overlap_frac_of_dex_genes_that_are_calls", fmt(n_ov / n_dex_bg))
put("background_dex_rate", fmt(n_dex_bg / n_bg))

# Enrichment (two-sided Fisher)
ft <- fisher.test(table(factor(bg$call, c(TRUE, FALSE)), factor(bg$dex_gene, c(TRUE, FALSE))))
put("enrichment_odds_ratio", fmt(ft$estimate))
put("enrichment_or_ci95", paste(fmt(ft$conf.int), collapse = ";"))
put("enrichment_fisher_p_two_sided", fmt(ft$p.value))
p_enrich <- ft$p.value

# Direction among overlapping calls
ov[, agree := sign(dex_logFC) == sign(dis_logFC)]
n_agree <- sum(ov$agree)
ci_ag <- cp_ci(n_agree, n_ov)
put("direction_agree_n", n_agree)
put("direction_agree_frac", fmt(n_agree / n_ov))
put("direction_agree_ci95", paste(fmt(ci_ag), collapse = ";"))
put("overlap_disease_up_DEX_up", ov[dis_logFC > 0 & dex_logFC > 0, .N])
put("overlap_disease_up_DEX_down", ov[dis_logFC > 0 & dex_logFC < 0, .N])
put("overlap_disease_down_DEX_up", ov[dis_logFC < 0 & dex_logFC > 0, .N])
put("overlap_disease_down_DEX_down", ov[dis_logFC < 0 & dex_logFC < 0, .N])
if (n_ov >= 2) {
  pu_d <- mean(ov$dis_logFC > 0)
  pu_x <- mean(ov$dex_logFC > 0)
  put("direction_expected_agree_under_independence", fmt(pu_d * pu_x + (1 - pu_d) * (1 - pu_x)))
  fd <- fisher.test(table(factor(ov$dis_logFC > 0, c(TRUE, FALSE)), factor(ov$dex_logFC > 0, c(TRUE, FALSE))))
  p_dir <- fd$p.value
  put("direction_fisher_odds_ratio", fmt(fd$estimate))
  put("direction_fisher_p_two_sided", fmt(p_dir))
} else {
  p_dir <- NA_real_
  put("direction_fisher_p_two_sided", "NA (overlap < 2)")
}

# Spearman correlations
r_bg <- spearman_ci(bg$dis_logFC, bg$dex_logFC)
r_calls <- spearman_ci(bg[call == TRUE, dis_logFC], bg[call == TRUE, dex_logFC])
put("spearman_bg_rho", fmt(r_bg["rho"]))
put("spearman_bg_ci95_fisher_z", paste(fmt(r_bg[c("lo", "hi")]), collapse = ";"))
put("spearman_bg_gene_level_p_not_calibrated", fmt(r_bg["p"]))
put("spearman_bg_n", r_bg["n"])
put("spearman_calls_rho", fmt(r_calls["rho"]))
put("spearman_calls_ci95_fisher_z", paste(fmt(r_calls[c("lo", "hi")]), collapse = ";"))
put("spearman_calls_gene_level_p_not_calibrated", fmt(r_calls["p"]))
put("spearman_calls_n", r_calls["n"])

# Descriptive bound on attributable calls
calls_bg <- bg[call == TRUE]
calls_bg[, agree := sign(dex_logFC) == sign(dis_logFC)]
for (d in c(1, 0.5, 0.25)) {
  put(paste0("bound_delta_", d, "_all_calls"),
    calls_bg[agree & abs(d * dex_logFC) >= abs(dis_logFC), .N])
  put(paste0("bound_delta_", d, "_dex_gene_calls"),
    calls_bg[agree & dex_gene & abs(d * dex_logFC) >= abs(dis_logFC), .N])
}
fwrite(calls_bg[order(-dex_gene, dex_padj)], file.path(out, "disease_calls_with_dex_stats.tsv"), sep = "\t")

# ---------------------------------------------------------------- matched null
set.seed(SEED)
bg[, bin := ceiling(N_BINS * rank(dis_AveExpr, ties.method = "first") / .N)]
calls_per_bin <- tabulate(bg$bin[bg$call], nbins = N_BINS)
pool_idx <- split(which(!bg$call), factor(bg$bin[!bg$call], levels = seq_len(N_BINS)))
pool_size <- lengths(pool_idx)
n_bins_with_replacement <- sum(pool_size < calls_per_bin)
put("null_bins_sampled_with_replacement", n_bins_with_replacement)
fwrite(data.table(bin = seq_len(N_BINS), calls = calls_per_bin, pool = pool_size,
  AveExpr_min = bg[, min(dis_AveExpr), by = bin][order(bin), V1],
  AveExpr_max = bg[, max(dis_AveExpr), by = bin][order(bin), V1]),
  file.path(out, "null_bins.tsv"), sep = "\t")
dex_vec <- bg$dex_gene
agree_vec <- bg$dex_gene & (sign(bg$dex_logFC) == sign(bg$dis_logFC))
null_ov <- integer(N_NULL)
null_ag <- integer(N_NULL)
for (k in seq_len(N_NULL)) {
  idx <- unlist(lapply(seq_len(N_BINS), function(b) {
    n <- calls_per_bin[b]
    if (n == 0) return(integer(0))
    p <- pool_idx[[b]]
    p[sample.int(length(p), n, replace = length(p) < n)]
  }), use.names = FALSE)
  null_ov[k] <- sum(dex_vec[idx])
  null_ag[k] <- sum(agree_vec[idx])
}
fwrite(data.table(set = seq_len(N_NULL), n_dex = null_ov, n_dex_concordant = null_ag),
  file.path(out, "null_matched_sets.tsv.gz"), sep = "\t")
emp_p2 <- function(obs, null) {
  up <- (1 + sum(null >= obs)) / (length(null) + 1)
  lo <- (1 + sum(null <= obs)) / (length(null) + 1)
  min(1, 2 * min(up, lo))
}
obs_ag <- sum(calls_bg$agree & calls_bg$dex_gene)
p_null <- emp_p2(n_ov, null_ov)
put("null_obs_overlap", n_ov)
put("null_mean_overlap", fmt(mean(null_ov)))
put("null_overlap_q025_q975", paste(quantile(null_ov, c(0.025, 0.975)), collapse = ";"))
put("null_obs_over_expected", fmt(n_ov / mean(null_ov)))
put("null_p_two_sided_overlap", fmt(p_null))
put("null_obs_concordant", obs_ag)
put("null_mean_concordant", fmt(mean(null_ag)))
put("null_concordant_q025_q975", paste(quantile(null_ag, c(0.025, 0.975)), collapse = ";"))
put("null_p_two_sided_concordant", fmt(emp_p2(obs_ag, null_ag)))

# Holm over the confirmatory family
fam <- c(enrichment_fisher = p_enrich, direction_fisher = p_dir, matched_null = p_null)
holm <- p.adjust(fam, method = "holm")
for (nm in names(holm)) put(paste0("holm_", nm), fmt(holm[[nm]]))

# ---------------------------------------------------------------- programs
mem_all <- fread(file.path(hot, "program_membership_v2.tsv"))
reg <- fread(file.path(hot, "program_registry_v2.tsv"))
stopifnot(uniqueN(mem_all$program_uid) == 117)
w_total <- mem_all[, .(n_members = .N, w_total = sum(original_l1_weight)), by = program_uid]
mem <- mem_all[mapped_symbol_status != "unmapped_or_ambiguous" & !is.na(mapped_symbol) & mapped_symbol != "NA"]
mem[, base := ifelse(grepl("^ENSG[0-9]+", mapped_symbol), sub("\\.[0-9]+$", "", mapped_symbol),
  sym_all$base[match(mapped_symbol, sym_all$gene_name)])]
put("program_rows_mapped", nrow(mem))
put("program_rows_mapped_without_v49_id", sum(is.na(mem$base)))
mem <- mem[!is.na(base)]
stopifnot(!anyDuplicated(mem[, .(program_uid, base)]))
mm <- merge(mem[, .(program_uid, base, w = original_l1_weight)],
  bg[, .(base, dex_gene, dex_logFC)], by = "base")
per <- mm[, .(n_measured = .N, w_measured = sum(w), n_dex = sum(dex_gene),
  n_dex_up = sum(dex_gene & dex_logFC > 0), n_dex_down = sum(dex_gene & dex_logFC < 0),
  w_dex_frac = sum(w * dex_gene) / sum(w), w_dex_score = sum(w * dex_logFC) / sum(w),
  dex_members = paste(sort(base[dex_gene]), collapse = ",")), by = program_uid]
per <- merge(w_total, per, by = "program_uid", all.x = TRUE)
per[is.na(n_measured), `:=`(n_measured = 0L, w_measured = 0, n_dex = 0L)]
per[, weight_retained := w_measured / w_total]
per[, frac_dex := ifelse(n_measured > 0, n_dex / n_measured, NA_real_)]
per[, testable := n_measured >= MIN_PROGRAM_GENES]
per[, fisher_p_greater := NA_real_]
per[, fisher_or := NA_real_]
for (i in which(per$testable)) {
  a <- per$n_dex[i]
  b <- per$n_measured[i] - a
  c0 <- n_dex_bg - a
  d0 <- n_bg - per$n_measured[i] - c0
  f <- fisher.test(matrix(c(a, c0, b, d0), 2), alternative = "greater")
  per$fisher_p_greater[i] <- f$p.value
  per$fisher_or[i] <- f$estimate
}
per[testable == TRUE, bh_q := p.adjust(fisher_p_greater, method = "BH")]
per[, flag_enriched := !is.na(bh_q) & bh_q < 0.05]
sym_lookup <- setNames(bg$symbol, bg$base)
per[, dex_member_symbols := vapply(strsplit(dex_members, ","), function(z)
  paste(sym_lookup[z][!is.na(sym_lookup[z])], collapse = ","), character(1))]
per <- merge(reg[, .(program_uid, cell_type, module, module_name)], per, by = "program_uid")
per[, dex_members := NULL]
setorder(per, bh_q, -frac_dex, na.last = TRUE)
fwrite(per, file.path(out, "program_dex_overlap.tsv"), sep = "\t")
put("programs_testable", per[testable == TRUE, .N])
put("programs_flagged_bh05", per[flag_enriched == TRUE, .N])
put("programs_flagged_names", per[flag_enriched == TRUE, paste(module_name, collapse = " | ")])

fwrite(rbindlist(summary_rows), file.path(out, "summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("V1 complete\n")
