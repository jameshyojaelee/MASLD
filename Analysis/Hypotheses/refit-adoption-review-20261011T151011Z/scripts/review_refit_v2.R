#!/usr/bin/env Rscript
# Independent check of the 774-unit five-cohort bulk refit.
# Assay: bulk RNA-seq fragment counts. Model: limma voomWithQualityWeights,
# ~ dataset + inferred_sex + group_binary. Unit: one library per genetic-identity cluster.
# Reads candidate files, refits from the raw count matrices, writes only under the review directory.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
})
set.seed(20260923)
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CAND <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates")
CORR <- file.path(ROOT, "figures/candidates/resource-corrections-20261005-v1")
VC <- file.path(CORR, "verified_counts_r8_20261007")
fit_dir <- function(s) file.path(CAND, paste0("bulk-corrected-controls-", s), "primary_source_controls")
F774 <- fit_dir("20261008T162308Z")
say <- function(...) cat(sprintf(...), "\n")
w <- function(x, name) fwrite(x, file.path(out, name), sep = "\t")

build_dge <- function(counts, meta) {
  stopifnot(all(meta$sample_id %in% colnames(counts)))
  m <- counts[, meta$sample_id, drop = FALSE]
  dge <- DGEList(counts = m)
  dge$samples$dataset <- meta$dataset
  dge$samples$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))
  keep <- filterByExpr(dge, design = model.matrix(~ 0 + group_binary, data = dge$samples))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  calcNormFactors(dge, method = "RLE")
}
fit_model <- function(dge, meta, with_dataset = TRUE) {
  info <- data.frame(
    group_binary = factor(dge$samples$group_binary, levels = c("Control", "Disease")),
    dataset = factor(dge$samples$dataset),
    inferred_sex = factor(meta$inferred_sex[match(colnames(dge), meta$sample_id)]),
    row.names = colnames(dge))
  stopifnot(!anyNA(info))
  terms <- c(if (with_dataset && nlevels(info$dataset) > 1) "dataset",
             if (nlevels(info$inferred_sex) > 1) "inferred_sex", "group_binary")
  design <- model.matrix(reformulate(terms), data = info)
  stopifnot(qr(design)$rank == ncol(design))
  v <- voomWithQualityWeights(dge, design)
  fit0 <- lmFit(v, design)
  fit <- eBayes(fit0)
  tt <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  tr <- topTreat(treat(fit0, lfc = 0.25), coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res <- data.table(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value,
                    padj = tt$adj.P.Val, AveExpr = tt$AveExpr,
                    treat_p = tr$P.Value, treat_fdr = tr$adj.P.Val)
  list(res = res, weights = data.table(sample_id = colnames(dge), sample_weight = v$targets$sample.weights),
       design = design, v = v)
}
call_of <- function(d) d$padj < 0.05 & abs(d$logFC) > 0.5
counts_of <- function(d) data.table(
  n_genes = nrow(d), canonical = sum(call_of(d)), up = sum(call_of(d) & d$logFC > 0),
  down = sum(call_of(d) & d$logFC < 0), treat = sum(d$treat_fdr < 0.05),
  treat_up = sum(d$treat_fdr < 0.05 & d$logFC > 0), treat_down = sum(d$treat_fdr < 0.05 & d$logFC < 0))
compare <- function(a, b, la, lb) {
  m <- merge(a[, .(gene, lfc_a = logFC, padj_a = padj, call_a = call_of(a))],
             b[, .(gene, lfc_b = logFC, padj_b = padj, call_b = call_of(b))], by = "gene", all = TRUE)
  cm <- m[!is.na(lfc_a) & !is.na(lfc_b)]
  data.table(a = la, b = lb, n_a = nrow(a), n_b = nrow(b), n_common = nrow(cm),
    calls_a = sum(m$call_a, na.rm = TRUE), calls_b = sum(m$call_b, na.rm = TRUE),
    both = m[, sum(call_a %in% TRUE & call_b %in% TRUE)],
    only_a = m[, sum(call_a %in% TRUE & !(call_b %in% TRUE))],
    only_a_untested_in_b = m[, sum(call_a %in% TRUE & is.na(lfc_b))],
    only_b = m[, sum(call_b %in% TRUE & !(call_a %in% TRUE))],
    only_b_untested_in_a = m[, sum(call_b %in% TRUE & is.na(lfc_a))],
    both_opposite_sign = cm[call_a & call_b, sum(sign(lfc_a) != sign(lfc_b))],
    calls_changed_among_common = cm[, sum(call_a != call_b)],
    max_abs_dlfc = cm[, max(abs(lfc_a - lfc_b))],
    pearson = cm[, cor(lfc_a, lfc_b)], spearman = cm[, cor(lfc_a, lfc_b, method = "spearman")])
}
read_deg <- function(p) { d <- fread(p); d }

if (mode == "refit774") {
  meta <- as.data.table(readRDS(file.path(F774, "meta_matched.rds")))
  counts <- readRDS(file.path(VC, "assembled_774_primary/counts_raw.rds"))
  stored <- read_deg(file.path(F774, "deg_results.csv"))
  say("counts matrix %d genes x %d libraries", nrow(counts), ncol(counts))

  ## 1. units
  roster <- fread(file.path(CORR, "identity_resource_qualified/bulk_primary_sample_ids.tsv"))
  full <- fread(file.path(CORR, "identity_resource_qualified/bulk_roster.tsv"))
  pairs <- fread(file.path(CORR, "identity_resource_qualified/full_filtered_related_pairs.tsv"))
  setnames(pairs, 1, "IID1")
  in774 <- full[sample_id %in% meta$sample_id]
  unit <- data.table(
    check = c("n_units", "n_control", "n_disease", "sample_id_unique", "analysis_unit_id_unique",
              "meta_equals_roster", "filtered_component_unique_among_774", "screen_component_unique_among_774",
              "source_person_id_nonNA_duplicated", "pairs_kinship_ge_0.354_both_in_774",
              "max_kinship_between_two_774_units", "inferred_sex_missing", "all_identity_qualified",
              "all_primary_selected", "counts_columns_equal_774"),
    value = c(nrow(meta), sum(meta$group_binary == "Control"), sum(meta$group_binary == "Disease"),
              !anyDuplicated(meta$sample_id), !anyDuplicated(meta$analysis_unit_id),
              setequal(meta$sample_id, roster$sample_id) && nrow(roster) == nrow(meta),
              !anyDuplicated(in774$filtered_component), !anyDuplicated(in774$screen_component),
              in774[!is.na(source_person_id) & source_person_id != "", sum(duplicated(source_person_id))],
              pairs[IID1 %in% meta$sample_id & IID2 %in% meta$sample_id & KINSHIP >= 0.354, .N],
              pairs[IID1 %in% meta$sample_id & IID2 %in% meta$sample_id, if (.N) max(KINSHIP) else NA_real_],
              sum(is.na(meta$inferred_sex) | meta$inferred_sex == ""),
              all(in774$identity_qualified %in% c("True", "TRUE", TRUE)),
              all(in774$primary_selected %in% c("True", "TRUE", TRUE)),
              setequal(colnames(counts), meta$sample_id)))
  w(unit, "01_unit_checks.tsv"); print(unit)
  w(pairs[IID1 %in% meta$sample_id & IID2 %in% meta$sample_id][order(-KINSHIP)], "01_related_pairs_within_774.tsv")
  w(meta[, .N, by = .(dataset, group_binary)][order(dataset, group_binary)], "01_cohort_group_counts.tsv")
  w(meta[, .N, by = .(dataset, sex_source, inferred_sex, reported_sex = sex)][order(dataset)], "01_sex_source.tsv")
  w(meta[group_binary == "Control", .N, by = .(dataset, condition, diagnosis_harmonized, nas_score, fibrosis_stage,
      source_control_status, strict_control_nas0)][order(dataset)], "07_control_definition.tsv")
  w(meta[group_binary == "Disease", .N, by = .(dataset, condition, diagnosis_harmonized)][order(dataset)], "07_disease_definition.tsv")

  ## expression-based sex check on the full count matrix
  gm <- fread(file.path(CORR, "annotation_resolution_r8/annotation_correspondence.tsv"))
  cpm_all <- cpm(counts[, meta$sample_id], log = TRUE, prior.count = 1)
  gid <- function(s) { g <- gm[gene_name == s & disposition == "primary_gene", gene_id]; g[g %in% rownames(counts)][1] }
  ysyms <- c("RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY")
  sx <- data.table(sample_id = meta$sample_id, dataset = meta$dataset, inferred_sex = meta$inferred_sex,
                   reported_sex = meta$sex, sex_source = meta$sex_source,
                   xist = cpm_all[gid("XIST"), ], y = colMeans(cpm_all[sapply(ysyms, gid), ]))
  sx[, expr_sex := fifelse(y - xist > 0, "M", "F")]
  w(sx, "01_expression_sex_per_unit.tsv")
  w(sx[, .N, by = .(dataset, inferred_sex, expr_sex)][order(dataset)], "01_expression_sex_vs_inferred.tsv")
  print(sx[, .N, by = .(dataset, inferred_sex, expr_sex)][order(dataset)])
  print(sx[inferred_sex != expr_sex])

  ## 2. rebuild and refit
  dge <- build_dge(counts, meta)
  sd <- readRDS(file.path(F774, "merged_dge.rds"))
  say("my DGE %d genes; stored DGE %d genes x %d; gene sets identical %s; counts identical %s; max norm factor diff %g",
      nrow(dge), nrow(sd), ncol(sd), setequal(rownames(dge), rownames(sd)),
      isTRUE(all.equal(unname(dge$counts[rownames(sd), colnames(sd)]), unname(sd$counts))),
      max(abs(dge$samples[colnames(sd), "norm.factors"] - sd$samples$norm.factors)))
  f <- fit_model(dge, meta)
  res <- f$res
  res[, symbol := stored$symbol[match(gene, stored$gene)]]
  w(res, "02_my_refit_774_deg.tsv")
  m <- merge(res, stored, by = "gene", suffixes = c(".mine", ".stored"))
  rep <- data.table(
    metric = c("genes_mine", "genes_stored", "genes_common", "max_abs_dlogFC", "max_abs_dpadj", "max_abs_dtreat_fdr",
               "canonical_call_sets_identical", "treat_call_sets_identical"),
    value = c(nrow(res), nrow(stored), nrow(m), m[, max(abs(logFC.mine - logFC.stored))],
              m[, max(abs(padj.mine - padj.stored))], m[, max(abs(treat_fdr.mine - treat_fdr.stored))],
              setequal(res[call_of(res), gene], stored[call_of(stored), gene]),
              setequal(res[treat_fdr < 0.05, gene], stored[treat_fdr < 0.05, gene])))
  w(rep, "02_reproduction.tsv"); print(rep)
  cc <- rbind(cbind(src = "mine", counts_of(res)), cbind(src = "stored", counts_of(stored)))
  w(cc, "02_call_counts.tsv"); print(cc)
  w(data.table(coef = colnames(f$design), rank = qr(f$design)$rank), "02_design_columns.tsv")

  ## 7. held annotations and contigs
  held <- fread(file.path(CORR, "annotation_resolution_r8/untestable_annotations.tsv"))
  tested <- merge(data.table(gene_id = res$gene), gm[, .(gene_id, gene_name, chromosome, disposition, gene_biotype)],
                  by = "gene_id", all.x = TRUE)
  primary_chr <- paste0("chr", c(1:22, "X", "Y", "M"))
  leak <- data.table(
    check = c("held_annotations", "held_in_count_matrix", "held_in_tested_set_by_gene_id",
              "held_in_tested_set_by_unversioned_id", "tested_genes_without_annotation_row",
              "tested_genes_not_on_primary_chromosomes", "tested_calls_not_on_primary_chromosomes"),
    value = c(nrow(held), sum(held$gene_id %in% rownames(counts)), sum(held$gene_id %in% res$gene),
              sum(sub("[.].*", "", held$gene_id) %in% sub("[.].*", "", res$gene)),
              sum(is.na(tested$disposition)), tested[!chromosome %in% primary_chr, .N],
              tested[!chromosome %in% primary_chr & gene_id %in% res[call_of(res), gene], .N]))
  w(leak, "07_held_annotation_leak.tsv"); print(leak)
  w(tested[, .N, by = .(disposition)], "07_tested_by_disposition.tsv"); print(tested[, .N, by = disposition])
  w(tested[!chromosome %in% primary_chr], "07_tested_off_primary_chromosomes.tsv")
  # symbols tested that share a name with a held annotation (reads shared with a held copy are not counted)
  w(res[symbol %in% held$gene_name, .(gene, symbol, logFC, padj, call = padj < 0.05 & abs(logFC) > 0.5)],
    "07_tested_genes_sharing_symbol_with_held.tsv")

  ## quality weights
  qw <- merge(f$weights, meta[, .(sample_id, dataset, group_binary)], by = "sample_id")
  sw <- fread(file.path(F774, "lvqw_sample_weights.tsv"))
  say("max abs diff of sample weights versus stored: %g", max(abs(qw$sample_weight - sw$sample_weight[match(qw$sample_id, sw$sample_id)])))
  qw[, lib_size := dge$samples[sample_id, "lib.size"]]
  w(qw, "07_quality_weights.tsv")
  qs <- qw[, .(n = .N, min = min(sample_weight), q05 = quantile(sample_weight, .05), median = median(sample_weight),
               q95 = quantile(sample_weight, .95), max = max(sample_weight), n_below_0.2 = sum(sample_weight < 0.2),
               n_above_3 = sum(sample_weight > 3), share_of_total_weight = sum(sample_weight) / sum(qw$sample_weight)),
           by = .(dataset, group_binary)][order(dataset, group_binary)]
  w(qs, "07_quality_weight_summary.tsv"); print(qs)
  print(qw[, .(share_weight = sum(sample_weight) / sum(qw$sample_weight), share_n = .N / nrow(qw)), by = dataset])
  print(qw[order(sample_weight)][c(1:8, (.N - 7):.N)])

  ## 6. marker genes
  mk <- c("AKR1B10", "THBS2", "COL1A1", "COL1A2", "COL3A1", "FABP4", "TREM2", "SPP1", "LUM", "EPCAM", "GPNMB",
          "CXCL10", "PNPLA3", "IGFBP2", "P4HA1", "FADS2", "ME1", "SCD", "MT1G", "IGF1", "TAT", "CYP7A1", "GNMT",
          "FOSB", "JUN", "MT1X", "CYP2C19", "APOF")
  mkres <- res[symbol %in% mk, .(symbol, gene, logFC, padj, treat_fdr, AveExpr,
                                 canonical = padj < 0.05 & abs(logFC) > 0.5)][order(-logFC)]
  w(mkres, "06_marker_genes.tsv"); print(mkres)
  print(res[call_of(res)][order(-logFC)][1:15, .(symbol, logFC, padj)])
  print(res[call_of(res)][order(logFC)][1:15, .(symbol, logFC, padj)])
  # per-cohort raw disease-minus-control log2 CPM difference for the markers (no model), to show direction by cohort
  lc <- cpm(dge, log = TRUE, prior.count = 1)
  pc <- rbindlist(lapply(mkres$gene, function(g) {
    d <- data.table(dataset = meta$dataset, grp = meta$group_binary, x = lc[g, meta$sample_id])
    d[, .(symbol = mkres$symbol[mkres$gene == g], diff_log2cpm = mean(x[grp == "Disease"]) - mean(x[grp == "Control"])), by = dataset]
  }))
  w(dcast(pc, symbol ~ dataset, value.var = "diff_log2cpm"), "06_marker_genes_raw_diff_by_cohort.tsv")
  print(dcast(pc, symbol ~ dataset, value.var = "diff_log2cpm"))

  ## GSE213621 alone versus the other four (is the pooled effect one cohort?)
  ma <- meta[dataset == "GSE213621"]
  fa <- fit_model(build_dge(counts, ma), ma)$res
  w(fa, "05_gse213621_alone_deg.tsv")
  cmp <- compare(res, fa, "full_774", "GSE213621_alone_326")
  w(cmp, "05_full_vs_gse213621_alone.tsv"); print(cmp)
}

if (mode == "loco") {
  meta <- as.data.table(readRDS(file.path(F774, "meta_matched.rds")))
  counts <- readRDS(file.path(VC, "assembled_774_primary/counts_raw.rds"))
  full <- read_deg(file.path(F774, "deg_results.csv"))
  fc <- full[call_of(full)]
  rows <- list(); detail <- list()
  cohorts <- sort(unique(meta$dataset))
  fits <- parallel::mclapply(cohorts, function(co) {
    mm <- meta[dataset != co]
    fit_model(build_dge(counts, mm), mm)$res
  }, mc.cores = 5, mc.preschedule = FALSE)
  names(fits) <- cohorts
  for (co in cohorts) {
    mm <- meta[dataset != co]
    mine <- fits[[co]]
    stopifnot(is.data.table(mine))
    w(mine, sprintf("05_my_loco_%s_deg.tsv", co))
    stored <- read_deg(file.path(CAND, sprintf("bulk-loco-%s-20261009T141221Z/deg_results.csv", co)))
    ms <- merge(mine, stored, by = "gene")
    x <- merge(fc[, .(gene, symbol, full_lfc = logFC, full_padj = padj)],
               mine[, .(gene, loco_lfc = logFC, loco_padj = padj, loco_p = P.Value)], by = "gene", all.x = TRUE)
    x[, kept := !is.na(loco_lfc) & loco_padj < 0.05 & abs(loco_lfc) > 0.5 & sign(loco_lfc) == sign(full_lfc)]
    x[, tested := !is.na(loco_lfc)]
    x[, same_sign := tested & sign(loco_lfc) == sign(full_lfc)]
    x[, reason := fifelse(kept, "kept", fifelse(!tested, "not_tested",
        fifelse(!same_sign, "opposite_sign",
        fifelse(loco_padj < 0.05, "same_sign_padj_ok_effect_below_0.5",
        fifelse(abs(loco_lfc) > 0.5, "same_sign_effect_ok_padj_above_0.05", "same_sign_both_miss")))))]
    x[, cohort := co]
    detail[[co]] <- x
    lost <- x[kept == FALSE]; lt <- lost[tested == TRUE]
    rows[[co]] <- data.table(left_out = co, n_units = nrow(mm), n_control = sum(mm$group_binary == "Control"),
      n_disease = sum(mm$group_binary == "Disease"), genes_tested = nrow(mine),
      genes_stored = nrow(stored), max_abs_dlogFC_vs_stored = ms[, max(abs(logFC.x - logFC.y))],
      loco_calls = sum(call_of(mine)), loco_calls_stored = sum(call_of(stored)),
      full_calls = nrow(x), kept_same_sign = sum(x$kept), pct_kept = round(100 * mean(x$kept), 1),
      lost = nrow(lost), lost_not_tested = sum(!lost$tested), lost_tested = nrow(lt),
      lost_tested_same_sign = sum(lt$same_sign), lost_opposite_sign = sum(!lt$same_sign),
      lost_same_sign_effect_only = sum(lt$reason == "same_sign_padj_ok_effect_below_0.5"),
      lost_same_sign_padj_only = sum(lt$reason == "same_sign_effect_ok_padj_above_0.05"),
      lost_same_sign_both = sum(lt$reason == "same_sign_both_miss"),
      lost_same_sign_nominal_p_below_0.05 = sum(lt$same_sign & lt$loco_p < 0.05),
      lost_same_sign_abs_lfc_above_0.25 = sum(lt$same_sign & abs(lt$loco_lfc) > 0.25),
      lost_same_sign_padj_below_0.05_and_abs_lfc_above_0.25 = sum(lt$same_sign & lt$loco_padj < 0.05 & abs(lt$loco_lfc) > 0.25),
      median_ratio_loco_over_full_lfc_lost_tested = lt[, median(loco_lfc / full_lfc)],
      pearson_full_calls = x[tested == TRUE, cor(full_lfc, loco_lfc)],
      spearman_full_calls = x[tested == TRUE, cor(full_lfc, loco_lfc, method = "spearman")],
      slope_loco_on_full_calls = x[tested == TRUE, coef(lm(loco_lfc ~ full_lfc))[2]],
      pearson_lost_tested = lt[, cor(full_lfc, loco_lfc)],
      sign_agree_all_full_calls_tested = x[tested == TRUE, mean(same_sign)])
    print(t(rows[[co]]))
  }
  w(rbindlist(rows), "05_loco_recomputed.tsv")
  d <- rbindlist(detail); w(d, "05_loco_per_call_detail.tsv")
  k <- dcast(d, gene + symbol ~ cohort, value.var = "kept")
  k[, n_lost := rowSums(!as.matrix(k[, -(1:2)]))]
  print(k[, .N, by = n_lost][order(n_lost)]); w(k[, .N, by = n_lost][order(n_lost)], "05_loco_losses_per_call.tsv")
  say("lost only when GSE213621 removed: %d", k[n_lost == 1 & GSE213621 == FALSE, .N])
  print(d[cohort == "GSE213621" & reason == "opposite_sign"])
  print(d[cohort == "GSE213621", .N, by = reason])
}

if (mode == "siblings") {
  d774 <- read_deg(file.path(F774, "deg_results.csv"))
  d774s <- read_deg(file.path(fit_dir("20261008T162309Z"), "deg_results.csv"))
  d738 <- read_deg(file.path(fit_dir("20261008T162310Z"), "deg_results.csv"))
  d738s <- read_deg(file.path(fit_dir("20261008T162311Z"), "deg_results.csv"))
  d844 <- read_deg(file.path(fit_dir("20260923T231939Z"), "deg_results.csv"))
  cc <- rbindlist(list(
    cbind(fit = "774_primary_stored", counts_of(d774)), cbind(fit = "774_sensitivity_count_stored", counts_of(d774s)),
    cbind(fit = "738_primary_stored", counts_of(d738)), cbind(fit = "738_sensitivity_count_stored", counts_of(d738s)),
    cbind(fit = "844_prior_stored", counts_of(d844))))
  # reproduction siblings
  for (s in c("20261008T162308Z", "20261008T162309Z", "20261008T162310Z", "20261008T162311Z")) {
    a <- read_deg(file.path(fit_dir(s), "deg_results.csv"))
    b <- read_deg(file.path(fit_dir(paste0(s, "-reproduction")), "deg_results.csv"))
    say("%s versus its reproduction: same genes %s, max abs dlogFC %g, call sets identical %s", s,
        identical(a$gene, b$gene), max(abs(a$logFC - b$logFC)),
        setequal(a[call_of(a), gene], b[call_of(b), gene]))
  }
  cmp <- rbindlist(list(
    compare(d844, d774, "prior_844", "774_primary"),
    compare(d774, d774s, "774_primary", "774_sensitivity_count"),
    compare(d738, d738s, "738_primary", "738_sensitivity_count"),
    compare(d774, d738, "774_primary", "738_primary")))
  w(cmp, "03_04_fit_comparisons.tsv"); print(t(cmp))
  # unversioned-id comparison of 844 vs 774, in case versions differ
  u <- function(d) { d <- copy(d); d[, gene := sub("[.].*", "", gene)]; d[!duplicated(gene)] }
  print(t(compare(u(d844), u(d774), "prior_844_unversioned", "774_unversioned")))
  m <- merge(d844[, .(gene, symbol_844 = symbol, lfc_844 = logFC, padj_844 = padj, call_844 = call_of(d844))],
             d774[, .(gene, symbol, lfc_774 = logFC, padj_774 = padj, call_774 = call_of(d774))], by = "gene", all = TRUE)
  flip <- m[call_844 %in% TRUE & call_774 %in% TRUE & sign(lfc_844) != sign(lfc_774)]
  w(flip, "03_sign_change_genes.tsv"); print(flip)
  lost <- m[call_844 %in% TRUE & !(call_774 %in% TRUE)]
  new <- m[call_774 %in% TRUE & !(call_844 %in% TRUE)]
  w(lost, "03_calls_844_not_774.tsv"); w(new, "03_calls_774_not_844.tsv")
  print(lost[, .(n = .N, untested_in_774 = sum(is.na(lfc_774)), same_sign = sum(sign(lfc_774) == sign(lfc_844), na.rm = TRUE),
                 padj_ok_effect_below = sum(padj_774 < 0.05 & abs(lfc_774) <= 0.5, na.rm = TRUE))])
  print(new[, .(n = .N, untested_in_844 = sum(is.na(lfc_844)), same_sign = sum(sign(lfc_774) == sign(lfc_844), na.rm = TRUE),
                padj_ok_effect_below_in_844 = sum(padj_844 < 0.05 & abs(lfc_844) <= 0.5, na.rm = TRUE))])
  # sensitivity counting arm: genes whose call differs
  ms <- merge(d774[, .(gene, symbol, lfc = logFC, padj, call = call_of(d774))],
              d774s[, .(gene, lfc_s = logFC, padj_s = padj, call_s = call_of(d774s))], by = "gene", all = TRUE)
  w(ms[is.na(call) | is.na(call_s) | call != call_s], "04_774_count_arm_differences.tsv")
  print(ms[(is.na(call) | is.na(call_s)) & (call %in% TRUE | call_s %in% TRUE)])

  # independent refit of the 738 primary and the 774 sensitivity counting arm from raw counts
  for (x in list(list("20261008T162310Z", "assembled_738_primary", d738, "738_primary"),
                 list("20261008T162309Z", "assembled_774_single_candidate_credit", d774s, "774_sensitivity_count"))) {
    meta <- as.data.table(readRDS(file.path(fit_dir(x[[1]]), "meta_matched.rds")))
    counts <- readRDS(file.path(VC, x[[2]], "counts_raw.rds"))
    say("%s: %d units, %d control, %d disease, unit ids unique %s, sex missing %d", x[[4]], nrow(meta),
        sum(meta$group_binary == "Control"), sum(meta$group_binary == "Disease"),
        !anyDuplicated(meta$analysis_unit_id), sum(is.na(meta$inferred_sex)))
    print(meta[, .N, by = .(dataset, group_binary)][order(dataset, group_binary)])
    mine <- fit_model(build_dge(counts, meta), meta)$res
    w(mine, sprintf("04_my_refit_%s_deg.tsv", x[[4]]))
    cc <- rbind(cc, cbind(fit = paste0(x[[4]], "_my_refit"), counts_of(mine)))
    mm <- merge(mine, x[[3]], by = "gene")
    say("%s my refit: genes %d (stored %d), max abs dlogFC %g, call sets identical %s", x[[4]], nrow(mine), nrow(x[[3]]),
        mm[, max(abs(logFC.x - logFC.y))], setequal(mine[call_of(mine), gene], x[[3]][call_of(x[[3]]), gene]))
  }
  w(cc, "04_call_counts_all_fits.tsv"); print(cc)
}
writeLines(capture.output(sessionInfo()), file.path(out, sprintf("sessionInfo_%s.txt", mode)))
say("DONE %s", mode)
