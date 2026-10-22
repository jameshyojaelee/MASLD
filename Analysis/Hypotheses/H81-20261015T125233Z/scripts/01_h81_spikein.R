#!/usr/bin/env Rscript
# H81: binomial-thinning spike-in calibration of the bulk evidence states on the real 774-unit counts.
# Prespecification: PRESPEC_H81.md in this run directory (sha256 checked before anything runs).
# Assay: bulk liver RNA-seq fragment counts (79,095 genes x 774 libraries, integer).
# Method: within-cohort label permutation, then binomial thinning (Gerard 2020, BMC Bioinformatics 21:206)
# of the raw counts, then the shipped limma voomWithQualityWeights recipe: one pooled fit
# (~ dataset + inferred_sex + group_binary), five leave-one-cohort-out fits, five per-cohort fits.
# Unit inside each fit: the genetic-identity cluster (774 units, 138 control / 636 disease, 5 cohorts).
# Unit of the equivalence rule: the cohort (k = 5). Spiked genes and replicates are calibration draws,
# not donors.
# Modes:
#   guard       G1 input hashes, G-A pooled reproduction, G-B LOCO reproduction, G-C per-cohort reproduction
#   replicates  replicates r = 1..n (seed 20260923 + r), each written once
#   debug       reduced data (unit and gene subsets) with non-primary seeds, for code checks only
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(parallel)
})
options(warn = 1, datatable.showProgress = FALSE)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H81 <- file.path(ROOT, "Analysis/Hypotheses/H81-20261015T125233Z")
PRESPEC_SHA <- "ae52c1949dda4f549c6f35c8ce7170815a5a1534e7beaed016aaedc030f46875"
RR <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
AA <- "Analysis/Hypotheses/refit-adoption-review-20261011T151011Z"
FF <- "figures/candidates/bulk-context-774-20261014T140904Z_22236311"
SS <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch"
WW <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H81"

# Input table of the prespec ("Inputs with sha256"), 23 files.
INPUTS <- data.table(
  id = c("counts", "meta", "deg", "merged_dge", "qw", "cohort_group_counts", "producer", "review_refit_v2",
          "repro_02", "loco_05", "loco_detail_05", "stage_ext_R", "restore_R", "cohort_results", "cohort_design_audit",
          "gate_per_gene", "gate_py", "ncount_py", "ncount_out", "count_n_R", "count_n_out", "design_params_py",
          "design_params_out"),
  path = c("figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds",
           file.path(RR, "meta_matched.rds"), file.path(RR, "deg_results.csv"), file.path(RR, "merged_dge.rds"),
           file.path(RR, "lvqw_sample_weights.tsv"), file.path(RR, "cohort_group_counts.tsv"),
           "scripts/manuscript/resource_f_five/refit_corrected_controls.R", file.path(AA, "scripts/review_refit_v2.R"),
           file.path(AA, "out/02_reproduction.tsv"), file.path(AA, "out/05_loco_recomputed.tsv"),
           file.path(AA, "out/05_loco_per_call_detail.tsv"), "scripts/figures/build_pi_stage_extensions.R",
           file.path(FF, "source/restore_constant_sex_cohort.R"), file.path(FF, "cohort_disease_all_gene_results.tsv"),
           file.path(FF, "cohort_disease_design_audit.tsv"), file.path(SS, "gate-L09-01/gate_per_gene.tsv.gz"),
           file.path(SS, "gate-L09-01/gate.py"), file.path(SS, "lens-L09/ncount.py"), file.path(SS, "lens-L09/ncount_output.txt"),
           file.path(WW, "count_n.R"), file.path(WW, "count_n.out"), file.path(WW, "design_params.py"),
           file.path(WW, "design_params.out")),
  bytes = c(52980522, 9082, 5049936, 33106829, 22417, 339, 45376, 22317, 238, 1620, 1457936, 23771, 3086, 25699368,
            453, 1985829, 11836, 5470, 1726, 1283, 691, 1267, 548),
  sha256 = c("f931e340546d679dc69dd0d6f712f4e73b1d805e80c741b22d3ed25498c68675",
             "429e27376ad421f5f5fae19c75f5e760030e727fab65c5768ef562a9939edfc3",
             "0491ed10024388c3497dccb9f99bf669683edd5903c9d21c5308850d7c10bde6",
             "73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b",
             "4e99563bd8d59dcba7efc81ab59ecb7f6f6fd5c4c1e8453226d93e435f0ec615",
             "5fcbcbb36ad151b6a5ee9104957c4379b709f6b19e55e88378efdd923911be8d",
             "5f11b3b580fc4e19772517986dafda9599cd651aa9f1e628db58914985f116d5",
             "7db250050ebf296cdfa025737a162f0325c4e2549136210917dd66f691c2495f",
             "a6cbf5cd28489f9dd5d0705f80d359c0f347d522c011e9ff80daad544984f92d",
             "640ede042f6476e117361b752aa812ea0a0b5522a31ac85018127108b6d746d0",
             "e83ad9f8a93621087ea4416cd5bbca5b57e41b3c972d4a9a1923010523d6a2b2",
             "8f6af913e321086cbc1f5640d411abadcdf345c1dcdd0b6dd785cebe65ce0d13",
             "5b172982a7715c14be26443551c01f4952d12025b7c404a28e4e0870146b5077",
             "327abb8964a605145e632df3cfb683f5a7f1c1834fa424bf26f36ae8acec21b0",
             "ba2554403e4b32730ee00232934afc290b4f4cad24a3226d4e847f7db27f1147",
             "675ae629c5b8edce0f90eb605936857df56047c3f877afd4586fe8f992d84f8a",
             "feb93ca4a1477e5e8f662d05d71a59f4ab50588a4af9af196ae7921f9d8b38fc",
             "6bdcb82821e5ddbf7fea78351173a4ceea6a58da288ef1a41fad4c8196569584",
             "5912cf5c45bb50bae0968e4ca99b81a2a539952f0b11654a411d0271f94d7b30",
             "f979da7d9d1e7c9832b10b06b65ef2dffe3714257d64346d4456b1929bbb168b",
             "4738ea5ef3c65e717119a5f829c1ccfff0c6e65a0f5d78f87b62b7e68c0435f2",
             "09a0c023b65489f97c383e9d596360d670401cf160516f57cfc8e43d2bc3db9a",
             "9d7c88b32af924480abe210cd8454f56375b89b29aa8c0c959b2f14794527fe4"))
inp <- function(k) file.path(ROOT, INPUTS$path[match(k, INPUTS$id)])

# Simulation design (prespec "Simulation design"), fixed.
SETS <- data.table(
  set = c("U0", "U0.5", "U0.75", "U1.0", "U1.5", "H1.0", "H0.5", "S1.0", "S2.0"),
  delta = c(0, 0.5, 0.75, 1.0, 1.5, 1.0, 0.5, 1.0, 2.0),
  type = c("U", "U", "U", "U", "U", "H", "H", "S", "S"),
  per_decile = c(40L, 40L, 40L, 40L, 40L, 40L, 40L, 20L, 20L))
TAU_H <- 0.306
MARGINS <- c(0.5, 0.25)

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[i + 1L]
}
mode <- getarg("mode"); out <- getarg("out")
cores <- as.integer(getarg("cores", "1"))
seed0 <- as.integer(getarg("seed", "20260923"))
n_rep <- as.integer(getarg("n-rep", "50"))
if (is.null(mode) || is.null(out)) stop("need --mode and --out")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(out, sprintf("log_%s.txt", mode))
say <- function(...) {
  s <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(...))
  cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = logf, append = TRUE)
}
write_once <- function(x, path) {
  if (file.exists(path)) stop("refusing to overwrite ", path)
  fwrite(x, path, sep = "\t", na = "NA")
}
sha256_file <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE)[1], "[[:space:]]+")[[1]][1]

check_prespec <- function() {
  h <- sha256_file(file.path(H81, "PRESPEC_H81.md"))
  if (!identical(h, PRESPEC_SHA)) stop("PRESPEC_H81.md sha256 mismatch: ", h)
  say("prespec sha256 OK %s", h)
}
check_inputs <- function() {
  # The copy of the writer's input_sha256.tsv must agree with the table hard-coded from the prespec.
  tsv <- fread(file.path(H81, "input_sha256.tsv"), header = FALSE, col.names = c("bytes", "sha256", "path"))
  agree <- nrow(tsv) == nrow(INPUTS) && all(INPUTS$path %in% tsv$path) &&
    all(tsv$sha256[match(INPUTS$path, tsv$path)] == INPUTS$sha256) &&
    all(tsv$bytes[match(INPUTS$path, tsv$path)] == INPUTS$bytes)
  r <- copy(INPUTS)
  r[, exists := file.exists(file.path(ROOT, path))]
  r[, observed_bytes := ifelse(exists, file.size(file.path(ROOT, path)), NA_real_)]
  r[, observed_sha256 := vapply(seq_len(.N), function(i) if (exists[i]) sha256_file(file.path(ROOT, path[i])) else NA_character_, "")]
  r[, match := exists & observed_bytes == bytes & observed_sha256 == sha256]
  r[, prespec_table_agrees_with_input_sha256_tsv := agree]
  r
}

# ---- Fit functions copied from review_refit_v2.R lines 22-52, TREAT lines removed (prespec step 5) ----
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
  res <- data.table(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value,
                    padj = tt$adj.P.Val, AveExpr = tt$AveExpr)
  res[, SE := logFC / t]
  list(res = res, weights = data.table(sample_id = colnames(dge), sample_weight = v$targets$sample.weights))
}
# Per-cohort fit as build_pi_stage_extensions.R fit_cohort (lines 277-312) and restore_constant_sex_cohort.R:
# subset the pooled filtered DGEList (pooled norm factors kept), units sorted by sample_id,
# ~ inferred_sex + group if sex varies, else ~ group; BH over the cohort fit's genes.
fit_cohort <- function(dge_all, meta, cohort_id) {
  selected <- copy(meta[dataset == cohort_id & group_binary %in% c("Control", "Disease")])
  selected[, group := factor(group_binary, levels = c("Control", "Disease"))]
  varying_sex <- uniqueN(selected$inferred_sex) > 1L
  stopifnot(uniqueN(selected$group) == 2L, !anyNA(selected$inferred_sex))
  setorder(selected, sample_id)
  dge <- dge_all[, match(selected$sample_id, colnames(dge_all))]
  info <- data.frame(inferred_sex = droplevels(factor(selected$inferred_sex)),
                     group = selected$group, row.names = selected$sample_id)
  design <- model.matrix(if (varying_sex) ~ inferred_sex + group else ~ group, info)
  stopifnot("groupDisease" %in% colnames(design), qr(design)$rank == ncol(design))
  transformed <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(transformed, design))
  table <- topTable(fit, coef = "groupDisease", number = Inf, sort.by = "none")
  genes <- rownames(dge)
  r <- data.table(gene = genes, dataset = cohort_id, logFC = table[genes, "logFC"], t = table[genes, "t"],
                  P.Value = table[genes, "P.Value"], FDR = p.adjust(table[genes, "P.Value"], "BH"),
                  covariates = if (varying_sex) "inferred_sex" else "none_constant_sex_omitted")
  r[, SE := logFC / t]
  r
}
# L09-01 candidate rule, coded as ncount.py c1(): DerSimonian-Laird tau2 over the five per-cohort log2FC and SE,
# random-effects mean, Hartung-Knapp SE truncated at 1, TOST p = max of two one-sided t(4) p values at each
# margin, BH over the universe. A gene with a non-finite per-cohort estimate gets p = 1 (kept in the BH family).
eq_rule <- function(L, S, margins = MARGINS) {
  k <- ncol(L)
  ok <- rowSums(!is.finite(L) | !is.finite(S) | S <= 0) == 0
  v <- S^2; w <- 1 / v
  sw <- rowSums(w); fe <- rowSums(w * L) / sw
  Q <- rowSums(w * (L - fe)^2)
  cc <- sw - rowSums(w^2) / sw
  tau2 <- pmax(0, (Q - (k - 1)) / cc)
  ws <- 1 / (v + tau2); sws <- rowSums(ws)
  mu <- rowSums(ws * L) / sws
  qhk <- rowSums(ws * (L - mu)^2) / (k - 1)
  sehk <- sqrt(pmax(qhk, 1) / sws)
  o <- data.table(mu = mu, sehk = sehk, tau2 = tau2, eq_ok = ok)
  for (m in margins) {
    p <- pmax(pt((mu + m) / sehk, df = 4, lower.tail = FALSE), pt((mu - m) / sehk, df = 4))
    p[!ok | is.na(p)] <- 1
    set(o, j = paste0("tost_p_", m), value = p)
    set(o, j = paste0("tost_q_", m), value = p.adjust(p, "BH"))
  }
  o
}
# Expression deciles fixed from the adopted fit: genes sorted by versioned Ensembl ID (C locale),
# AveExpr ranked with ties by that order, 10 equal-count bins (pandas qcut on the ranks).
decile_table <- function(stored) {
  d <- stored[, .(gene, AveExpr)][order(gene)]
  r <- frank(d$AveExpr, ties.method = "first")
  br <- quantile(r, probs = seq(0, 1, 0.1), type = 7, names = FALSE)
  d[, decile := as.integer(cut(r, breaks = br, include.lowest = TRUE, right = TRUE))]
  d[]
}
call_of <- function(lfc, padj) !is.na(lfc) & padj < 0.05 & abs(lfc) > 0.5

# =====================================================================================================
if (mode == "guard") {
  check_prespec()
  gdir <- file.path(out, "guard"); dir.create(gdir, showWarnings = FALSE)
  status <- data.table(rule = c("S1", "S2", "S3", "S4"), triggered = NA, detail = NA_character_)
  checks <- list()
  addc <- function(gate, check, value, target, pass) {
    checks[[length(checks) + 1L]] <<- data.table(gate = gate, check = check, value = as.character(value),
                                                 target = as.character(target), pass = pass)
  }
  # ---- G1 ----
  h <- check_inputs()
  write_once(h, file.path(gdir, "g1_input_hashes.tsv"))
  g1 <- all(h$match) && all(h$prespec_table_agrees_with_input_sha256_tsv)
  addc("G1", "inputs matching prespec sha256 and bytes", sum(h$match), nrow(h), g1)
  addc("G1", "prespec table equals writer input_sha256.tsv", h$prespec_table_agrees_with_input_sha256_tsv[1], TRUE,
       h$prespec_table_agrees_with_input_sha256_tsv[1])
  say("G1 input hashes: %d of %d match", sum(h$match), nrow(h))
  status[rule == "S1", `:=`(triggered = !g1, detail = paste(h[match == FALSE, id], collapse = ","))]
  if (!g1) {
    write_once(status, file.path(gdir, "guard_status.tsv"))
    write_once(rbindlist(checks), file.path(gdir, "guard_checks.tsv"))
    say("STOP S1: input sha256 mismatch: %s", status[rule == "S1", detail])
    quit(save = "no", status = 0)
  }
  say("G1 PASSED: all %d input sha256 match", nrow(h))

  meta <- as.data.table(readRDS(inp("meta")))
  counts <- readRDS(inp("counts"))
  stored <- fread(inp("deg"))
  md <- readRDS(inp("merged_dge"))
  cohorts <- sort(unique(meta$dataset))
  say("counts %d x %d (%s); meta %d units; cohorts %s", nrow(counts), ncol(counts), storage.mode(counts),
      nrow(meta), paste(cohorts, collapse = ","))

  # deciles
  dec <- decile_table(stored)
  fwrite(dec, file.path(gdir, "decile_table.tsv.gz"), sep = "\t")
  ds <- dec[, .(n = .N, min_ave = min(AveExpr), max_ave = max(AveExpr)), by = decile][order(decile)]
  target_n <- c(2412, 2411, 2412, 2411, 2412, 2411, 2411, 2412, 2411, 2412)
  target_edges <- c(-4.738, -2.739, -2.739, -2.02, -2.02, -1.149, -1.149, -0.063, -0.062, 1.284, 1.284, 2.579,
                    2.58, 3.564, 3.565, 4.435, 4.435, 5.496, 5.496, 16.721)
  dec_ok <- identical(as.numeric(ds$n), target_n) &&
    isTRUE(all.equal(as.numeric(t(round(as.matrix(ds[, .(min_ave, max_ave)]), 3))), target_edges))
  addc("deciles", "decile sizes and AveExpr edges equal design_params.out", paste(ds$n, collapse = "/"),
       paste(target_n, collapse = "/"), dec_ok)
  write_once(ds, file.path(gdir, "decile_summary.tsv"))

  # ---- 11 guard fits in parallel ----
  tasks <- c("pooled", paste0("loco:", cohorts), paste0("cohort:", cohorts))
  say("starting %d guard fits on %d cores", length(tasks), cores)
  fits <- mclapply(tasks, function(tk) {
    t0 <- proc.time()[["elapsed"]]
    o <- if (tk == "pooled") {
      d <- build_dge(counts, meta)
      f <- fit_model(d, meta)
      oo <- match(colnames(d), colnames(md))
      list(res = f$res, weights = f$weights,
           dge_check = data.table(
             genes = nrow(d), identical_gene_set_order = identical(rownames(d), rownames(md)),
             max_abs_count_diff = max(abs(d$counts - md$counts[, oo])),
             max_abs_lib_size_diff = max(abs(d$samples$lib.size - md$samples$lib.size[oo])),
             max_abs_norm_factor_diff = max(abs(d$samples$norm.factors - md$samples$norm.factors[oo]))))
    } else if (startsWith(tk, "loco:")) {
      co <- sub("^loco:", "", tk); mm <- meta[dataset != co]
      list(res = fit_model(build_dge(counts, mm), mm)$res)
    } else {
      co <- sub("^cohort:", "", tk)
      list(res = fit_cohort(md, meta, co))
    }
    o$task <- tk; o$elapsed <- proc.time()[["elapsed"]] - t0
    o
  }, mc.cores = cores, mc.preschedule = FALSE)
  names(fits) <- tasks
  bad <- vapply(fits, function(x) inherits(x, "try-error") || is.null(x$res), TRUE)
  if (any(bad)) {
    for (b in which(bad)) say("guard fit %s failed: %s", tasks[b], as.character(fits[[b]]))
    stop("guard fits failed")
  }
  write_once(data.table(task = tasks, elapsed_s = vapply(fits, `[[`, 0, "elapsed")), file.path(gdir, "guard_fit_times.tsv"))

  # ---- G-A ----
  res <- fits$pooled$res
  m <- merge(res, stored[, .(gene, logFC_s = logFC, padj_s = padj)], by = "gene")
  ga_genes <- nrow(res) == 24115L && nrow(stored) == 24115L && nrow(m) == 24115L && setequal(res$gene, stored$gene)
  calls_mine <- res[call_of(logFC, padj), gene]; calls_st <- stored[call_of(logFC, padj), gene]
  ga_calls <- length(calls_mine) == 1854L && setequal(calls_mine, calls_st)
  dl <- m[, max(abs(logFC - logFC_s))]; dp <- m[, max(abs(padj - padj_s))]
  addc("G-A", "genes (mine, stored, common)", sprintf("%d,%d,%d", nrow(res), nrow(stored), nrow(m)), "24115,24115,24115", ga_genes)
  addc("G-A", "canonical calls identical set", length(calls_mine), 1854, ga_calls)
  addc("G-A", "max |dlogFC|", dl, "<1e-8", dl < 1e-8)
  addc("G-A", "max |dpadj|", dp, "<1e-8", dp < 1e-8)
  ga <- ga_genes && ga_calls && dl < 1e-8 && dp < 1e-8
  status[rule == "S2", `:=`(triggered = !ga, detail = sprintf("genes %d calls %d dlogFC %.3g dpadj %.3g", nrow(res), length(calls_mine), dl, dp))]
  say("G-A %s: genes %d, calls %d (identical set %s), max|dlogFC| %.3g, max|dpadj| %.3g",
      if (ga) "PASSED" else "FAILED", nrow(res), length(calls_mine), ga_calls, dl, dp)
  # quality weights against the adopted fit (diagnostic 7 reference)
  sw <- fread(inp("qw"))
  qwd <- max(abs(fits$pooled$weights$sample_weight - sw$sample_weight[match(fits$pooled$weights$sample_id, sw$sample_id)]))
  addc("diag7", "max |quality weight - adopted|", qwd, "report", NA)
  write_once(data.table(stat = c("min", "median", "max"),
                        adopted = c(min(sw$sample_weight), median(sw$sample_weight), max(sw$sample_weight)),
                        guard_refit = c(min(fits$pooled$weights$sample_weight), median(fits$pooled$weights$sample_weight),
                                        max(fits$pooled$weights$sample_weight))),
             file.path(gdir, "quality_weight_range_adopted.tsv"))

  # ---- G-B ----
  tg_tested <- c(24628, 24022, 24071, 24196, 22451); tg_calls <- c(1936, 1860, 1873, 1945, 1960)
  tg_kept <- c(1795, 1750, 1745, 1643, 1236)
  full_calls <- res[call_of(logFC, padj), .(gene, full_lfc = logFC)]
  detail <- rbindlist(lapply(cohorts, function(co) {
    lr <- fits[[paste0("loco:", co)]]$res
    x <- merge(full_calls, lr[, .(gene, loco_lfc = logFC, loco_padj = padj)], by = "gene", all.x = TRUE)
    x[, kept := !is.na(loco_lfc) & loco_padj < 0.05 & abs(loco_lfc) > 0.5 & sign(loco_lfc) == sign(full_lfc)]
    x[, cohort := co]
    x
  }))
  gb <- data.table(left_out = cohorts,
                   genes_tested = vapply(cohorts, function(co) nrow(fits[[paste0("loco:", co)]]$res), 0L),
                   calls = vapply(cohorts, function(co) { r <- fits[[paste0("loco:", co)]]$res; sum(call_of(r$logFC, r$padj)) }, 0L),
                   kept_same_sign = vapply(cohorts, function(co) detail[cohort == co, sum(kept)], 0L))
  gb[, `:=`(target_tested = tg_tested, target_calls = tg_calls, target_kept = tg_kept)]
  kept_all <- detail[, .(all5 = all(kept)), by = gene]
  n_all5 <- kept_all[, sum(all5)]
  stored_detail <- fread(inp("loco_detail_05"))
  cmpk <- merge(detail[, .(gene, cohort, kept)], stored_detail[, .(gene, cohort, kept_s = as.logical(kept))], by = c("gene", "cohort"))
  gb_ok <- all(gb$genes_tested == gb$target_tested) && all(gb$calls == gb$target_calls) &&
    all(gb$kept_same_sign == gb$target_kept) && n_all5 == 1026L && nrow(cmpk) == 9270L && all(cmpk$kept == cmpk$kept_s)
  write_once(gb, file.path(gdir, "gb_loco_counts.tsv"))
  fwrite(detail, file.path(gdir, "gb_loco_per_call.tsv.gz"), sep = "\t")
  addc("G-B", "LOCO genes tested", paste(gb$genes_tested, collapse = "/"), paste(tg_tested, collapse = "/"), all(gb$genes_tested == tg_tested))
  addc("G-B", "LOCO calls", paste(gb$calls, collapse = "/"), paste(tg_calls, collapse = "/"), all(gb$calls == tg_calls))
  addc("G-B", "LOCO kept same sign", paste(gb$kept_same_sign, collapse = "/"), paste(tg_kept, collapse = "/"), all(gb$kept_same_sign == tg_kept))
  addc("G-B", "calls kept in all five (SD-any flags 1854 - this)", sprintf("%d (SD-any %d)", n_all5, 1854L - n_all5), "1026 (828)", n_all5 == 1026L)
  addc("G-B", "per-call kept equals 05_loco_per_call_detail.tsv", sprintf("%d rows, %d equal", nrow(cmpk), sum(cmpk$kept == cmpk$kept_s)), "9270 rows all equal",
       nrow(cmpk) == 9270L && all(cmpk$kept == cmpk$kept_s))
  status[rule == "S3", `:=`(triggered = !gb_ok, detail = sprintf("kept in all five %d", n_all5))]
  say("G-B %s: tested %s calls %s kept %s all-five %d", if (gb_ok) "PASSED" else "FAILED",
      paste(gb$genes_tested, collapse = "/"), paste(gb$calls, collapse = "/"), paste(gb$kept_same_sign, collapse = "/"), n_all5)

  # ---- G-C ----
  ct <- fread(inp("cohort_results"))
  da <- fread(inp("cohort_design_audit"))
  pc <- rbindlist(lapply(cohorts, function(co) fits[[paste0("cohort:", co)]]$res))
  mc <- merge(pc[, .(gene, dataset, lfc = logFC)], ct[, .(gene = gene_id_versioned, dataset, lfc_t = logFC)], by = c("gene", "dataset"))
  gc_rows <- nrow(pc) == 120575L && nrow(ct) == 120575L && nrow(mc) == 120575L
  gc_tol <- mc[, max(abs(lfc - lfc_t))]
  cov_mine <- pc[, .(cov = covariates[1]), by = dataset]
  cov_ok <- all(cov_mine$cov == da$covariates[match(cov_mine$dataset, da$dataset)])
  U <- rownames(md)
  Lm <- sapply(cohorts, function(co) { r <- fits[[paste0("cohort:", co)]]$res; r$logFC[match(U, r$gene)] })
  Sm <- sapply(cohorts, function(co) { r <- fits[[paste0("cohort:", co)]]$res; r$SE[match(U, r$gene)] })
  eq <- eq_rule(Lm, Sm)
  n05 <- sum(eq$tost_q_0.5 < 0.05); n025 <- sum(eq$tost_q_0.25 < 0.05)
  dge_check <- fits$pooled$dge_check
  build_ok <- dge_check$genes == 24115L && dge_check$identical_gene_set_order && dge_check$max_abs_count_diff == 0 &&
    dge_check$max_abs_lib_size_diff == 0 && dge_check$max_abs_norm_factor_diff == 0
  write_once(dge_check, file.path(gdir, "gc_build_dge_vs_merged.tsv"))
  addc("G-C", "per-cohort rows (mine, file, matched)", sprintf("%d,%d,%d", nrow(pc), nrow(ct), nrow(mc)), "120575", gc_rows)
  addc("G-C", "max |dlogFC| per-cohort", gc_tol, "<1e-6", gc_tol < 1e-6)
  addc("G-C", "covariates as design audit", paste(cov_mine$dataset, cov_mine$cov, sep = "=", collapse = ";"), "as cohort_disease_design_audit.tsv", cov_ok)
  addc("G-C", "L09-01 rule calls at margin 0.5", n05, 8548, n05 == 8548L)
  addc("G-C", "L09-01 rule calls at margin 0.25", n025, 0, n025 == 0L)
  addc("G-C", "build_dge reproduces merged_dge.rds (genes, order, counts, lib.size, norm.factors)",
       sprintf("%d genes, order %s, max diffs %g/%g/%g", dge_check$genes, dge_check$identical_gene_set_order,
               dge_check$max_abs_count_diff, dge_check$max_abs_lib_size_diff, dge_check$max_abs_norm_factor_diff), "identical", build_ok)
  gc_core <- n05 == 8548L && n025 == 0L && cov_ok && build_ok && gc_rows
  gc_all <- gc_core && gc_tol < 1e-6
  status[rule == "S4", `:=`(triggered = !gc_core,
                            detail = sprintf("8548 reproduced %s; 0 at 0.25 %s; tolerance %s (%.3g); covariates %s; build_dge %s",
                                             n05 == 8548L, n025 == 0L, gc_tol < 1e-6, gc_tol, cov_ok, build_ok))]
  say("G-C %s: rows %d, max|dlogFC| %.3g, covariates %s, margin 0.5 calls %d, margin 0.25 calls %d, build_dge==merged %s",
      if (gc_all) "PASSED" else if (gc_core) "PASSED except tolerance (arm proceeds)" else "FAILED",
      nrow(mc), gc_tol, cov_ok, n05, n025, build_ok)
  fwrite(cbind(gene = U, eq), file.path(gdir, "gc_equivalence_rule_adopted.tsv.gz"), sep = "\t")
  fwrite(res, file.path(gdir, "ga_pooled_refit.tsv.gz"), sep = "\t")

  ck <- rbindlist(checks)
  write_once(ck, file.path(gdir, "guard_checks.tsv"))
  status[, gate_pass := c(g1, ga, gb_ok, gc_all)]
  write_once(status, file.path(gdir, "guard_status.tsv"))
  print(ck)
  if (status[rule == "S2", triggered]) say("STOP S2: G-A failed; no number from the run is used")
  say("GUARD SUMMARY: G1 %s | G-A %s | G-B %s | G-C %s | deciles %s", g1, ga, gb_ok, gc_all, dec_ok)
  writeLines(capture.output(sessionInfo()), file.path(gdir, "sessionInfo_guard.txt"))
  say("DONE guard")
}

# =====================================================================================================
run_replicate <- function(r, ctx) {
  with(ctx, {
    rdir <- file.path(out, "replicates")
    f_genes <- file.path(rdir, sprintf("rep_%03d_genes.tsv.gz", r))
    f_sum <- file.path(rdir, sprintf("rep_%03d_summary.tsv", r))
    f_done <- file.path(rdir, sprintf("rep_%03d.done", r))
    if (file.exists(f_done) || file.exists(f_genes)) stop("replicate already written: ", r)
    t_start <- proc.time()[["elapsed"]]
    seed_r <- seed_base + r
    set.seed(seed_r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    # 2. permutation within cohort, sorted cohort order, metadata row order
    perm <- meta$group_binary
    for (co in cohorts) {
      idx <- which(meta$dataset == co)
      perm[idx] <- sample(meta$group_binary[idx])
    }
    meta_p <- copy(meta); meta_p[, group_binary := perm]
    lab_tab <- meta_p[, .(n_control = sum(group_binary == "Control"), n_disease = sum(group_binary == "Disease")), by = dataset][order(dataset)]
    chg_frac <- meta[, .(frac_label_changed = mean(group_binary != perm[.I])), by = dataset][order(dataset)]
    # 3. gene assignment within each fixed decile
    rows <- vector("list", 10L * nrow(SETS)); k <- 0L
    for (d in 1:10) {
      g <- dec_genes[[d]]
      s <- sample(g, sum(SETS$per_decile))
      pos <- 0L
      for (i in seq_len(nrow(SETS))) {
        n <- SETS$per_decile[i]
        gg <- s[pos + seq_len(n)]; pos <- pos + n
        k <- k + 1L
        rows[[k]] <- data.table(set = SETS$set[i], set_idx = i, decile = d, pos = seq_len(n), gene = gg,
                                dir = c(rep(1L, n / 2), rep(-1L, n / 2)), delta = SETS$delta[i], type = SETS$type[i],
                                spiked_cohort = if (SETS$type[i] == "S") rep(rep(cohorts, each = 2L), 2L) else NA_character_)
      }
    }
    sp <- rbindlist(rows); setorder(sp, set_idx, decile, pos)
    stopifnot(!anyDuplicated(sp$gene), nrow(sp) == 10L * sum(SETS$per_decile))
    # per-cohort intended log2 change (disease over control): dir * delta_c
    DC <- matrix(0, nrow(sp), length(cohorts), dimnames = list(NULL, cohorts))
    iu <- which(sp$type == "U"); DC[iu, ] <- sp$delta[iu]
    ih <- which(sp$type == "H")  # set order: all H1.0 genes, then all H0.5 genes
    Z <- matrix(rnorm(length(cohorts) * length(ih)), ncol = length(cohorts), byrow = TRUE)
    DC[ih, ] <- sp$delta[ih] + TAU_H * Z
    is_ <- which(sp$type == "S")
    DC[cbind(is_, match(sp$spiked_cohort[is_], cohorts))] <- sp$delta[is_]
    CHG <- DC * sp$dir
    # 4. thinning on the raw matrix, gene by gene in set order, units in metadata order
    sub0 <- counts_m[sp$gene, , drop = FALSE]
    sub <- sub0
    is_ctrl <- perm == "Control"
    ucoh <- match(meta$dataset, cohorts)
    for (i in seq_len(nrow(sp))) {
      ch <- CHG[i, ]
      if (all(ch == 0)) next
      pu <- rep(1, n_units)
      for (cc in which(ch != 0)) {
        side <- if (ch[cc] > 0) is_ctrl else !is_ctrl
        pu[ucoh == cc & side] <- 2^(-abs(ch[cc]))
      }
      idx <- which(pu < 1)
      sub[i, idx] <- rbinom(length(idx), sub[i, idx], pu[idx])
    }
    storage.mode(sub) <- "integer"
    # diagnostic 4: realized change of log2(sum CPM disease / sum CPM control), permuted labels
    lib1 <- lib0 - colSums(sub0 - sub)
    cpm0 <- sweep(sub0, 2, lib0, "/"); cpm1 <- sweep(sub, 2, lib1, "/")
    is_dis <- !is_ctrl
    realized <- vapply(seq_len(nrow(sp)), function(i) {
      u <- if (sp$type[i] == "S") ucoh == match(sp$spiked_cohort[i], cohorts) else rep(TRUE, n_units)
      b <- log2(sum(cpm0[i, u & is_dis]) / sum(cpm0[i, u & is_ctrl]))
      a <- log2(sum(cpm1[i, u & is_dis]) / sum(cpm1[i, u & is_ctrl]))
      a - b
    }, 0)
    rm(cpm0, cpm1)
    counts_r <- counts_m
    counts_r[sp$gene, ] <- sub
    rm(sub0, sub)
    # 5. fits
    tt <- c()
    t0 <- proc.time()[["elapsed"]]
    dge_p <- build_dge(counts_r, meta_p)
    fp <- fit_model(dge_p, meta_p)
    tt["pooled"] <- proc.time()[["elapsed"]] - t0
    fl <- list()
    for (co in cohorts) {
      t0 <- proc.time()[["elapsed"]]
      mm <- meta_p[dataset != co]
      fl[[co]] <- fit_model(build_dge(counts_r, mm), mm)$res
      tt[paste0("loco_", co)] <- proc.time()[["elapsed"]] - t0
    }
    rm(counts_r)
    fc <- list()
    for (co in cohorts) {
      t0 <- proc.time()[["elapsed"]]
      fc[[co]] <- fit_cohort(dge_p, meta_p, co)
      tt[paste0("cohort_", co)] <- proc.time()[["elapsed"]] - t0
    }
    # evidence states over the replicate's pooled universe
    st <- fp$res[, .(gene, logFC, SE, padj, AveExpr_rep = AveExpr)]
    st[, supported := padj < 0.05 & abs(logFC) > 0.5]
    keptm <- matrix(FALSE, nrow(st), length(cohorts)); sigm <- keptm
    for (j in seq_along(cohorts)) {
      co <- cohorts[j]
      x <- fl[[co]][match(st$gene, gene)]
      tested <- !is.na(x$logFC)
      kept <- tested & x$padj < 0.05 & abs(x$logFC) > 0.5 & sign(x$logFC) == sign(st$logFC)
      sk <- tested & x$padj < 0.05 & sign(x$logFC) == sign(st$logFC)
      kept[is.na(kept)] <- FALSE; sk[is.na(sk)] <- FALSE
      keptm[, j] <- kept; sigm[, j] <- sk
      set(st, j = paste0("loco_tested_", co), value = tested)
      set(st, j = paste0("loco_logFC_", co), value = x$logFC)
      set(st, j = paste0("loco_padj_", co), value = x$padj)
      set(st, j = paste0("loco_kept_", co), value = kept)
    }
    st[, sd_any := supported & rowSums(!keptm) > 0]
    st[, sd_sig := supported & rowSums(!sigm) > 0]
    Lm <- matrix(NA_real_, nrow(st), length(cohorts)); Sm <- Lm
    for (j in seq_along(cohorts)) {
      co <- cohorts[j]
      x <- fc[[co]][match(st$gene, gene)]
      Lm[, j] <- x$logFC; Sm[, j] <- x$SE
      set(st, j = paste0("coh_logFC_", co), value = x$logFC)
      set(st, j = paste0("coh_SE_", co), value = x$SE)
      set(st, j = paste0("coh_FDR_", co), value = x$FDR)
    }
    eq <- eq_rule(Lm, Sm)
    st <- cbind(st, eq)
    st[, tested_negative := tost_q_0.5 < 0.05]
    st[, tested_negative_m025 := tost_q_0.25 < 0.05]
    # spiked gene table
    G <- merge(sp, st, by = "gene", all.x = TRUE, sort = FALSE)
    setorder(G, set_idx, decile, pos)
    G[, tested := !is.na(logFC)]
    G[tested == FALSE, `:=`(supported = FALSE, tested_negative = FALSE, tested_negative_m025 = FALSE)]
    G[, realized_change := realized[match(gene, sp$gene)]]
    G[, realized_aligned := realized_change * dir]
    for (j in seq_along(cohorts)) set(G, j = paste0("delta_c_", cohorts[j]), value = DC[match(G$gene, sp$gene), j])
    KM <- as.matrix(G[, paste0("loco_kept_", cohorts), with = FALSE])
    iS <- which(G$type == "S")
    kc <- KM[cbind(iS, match(G$spiked_cohort[iS], cohorts))]
    G[, sd_own := NA]
    G[iS, sd_own := supported & !(kc %in% TRUE)]
    G[, `:=`(rep = r, seed = seed_r)]
    bg <- st[!gene %in% sp$gene]
    smry <- data.table(
      rep = r, seed = seed_r, universe = nrow(st), unspiked_in_universe = nrow(bg),
      spiked_tested = sum(G$tested), supported_unspiked = sum(bg$supported),
      tested_negative_unspiked = sum(bg$tested_negative), tested_negative_all = sum(st$tested_negative),
      tested_negative_all_m025 = sum(st$tested_negative_m025), eq_nonfinite = sum(!st$eq_ok),
      eff_alpha = 0.05 * sum(st$tested_negative) / nrow(st),
      qw_min = min(fp$weights$sample_weight), qw_median = median(fp$weights$sample_weight),
      qw_max = max(fp$weights$sample_weight),
      loco_universe = paste(vapply(cohorts, function(co) nrow(fl[[co]]), 0L), collapse = "/"),
      labels = paste(sprintf("%s:%d/%d", lab_tab$dataset, lab_tab$n_control, lab_tab$n_disease), collapse = ";"),
      frac_label_changed = paste(sprintf("%s:%.4f", chg_frac$dataset, chg_frac$frac_label_changed), collapse = ";"),
      fit_seconds = paste(sprintf("%s=%.0f", names(tt), tt), collapse = ";"),
      total_seconds = proc.time()[["elapsed"]] - t_start)
    fwrite(G, f_genes, sep = "\t", na = "NA")
    fwrite(smry, f_sum, sep = "\t", na = "NA")
    writeLines(format(Sys.time()), f_done)
    sprintf("rep %d done in %.0f s; universe %d; unspiked supported %d", r, smry$total_seconds, smry$universe,
            smry$supported_unspiked)
  })
}

if (mode %in% c("replicates", "debug")) {
  check_prespec()
  debug <- mode == "debug"
  if (!debug) {
    gs <- file.path(out, "guard", "guard_status.tsv")
    if (!file.exists(gs)) stop("guard_status.tsv missing; run --mode guard first")
    gst <- fread(gs)
    if (isTRUE(gst[rule == "S1", triggered]) || isTRUE(gst[rule == "S2", triggered])) {
      say("guard stopped the run (S1 or S2); no replicate is simulated")
      quit(save = "no", status = 0)
    }
    h <- check_inputs()
    if (!all(h$match)) stop("input sha256 mismatch at replicate start: ", paste(h[match == FALSE, id], collapse = ","))
    say("input sha256 re-checked: %d of %d match", sum(h$match), nrow(h))
  }
  meta <- as.data.table(readRDS(inp("meta")))
  counts <- readRDS(inp("counts"))
  stored <- fread(inp("deg"))
  dec <- decile_table(stored)
  cohorts <- sort(unique(meta$dataset))
  seed_base <- seed0
  if (debug) {
    # Debug only: unit and gene subsets, seeds offset by 500000 so no primary replicate is ever computed.
    set.seed(seed0 + 499999L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    keep_u <- unlist(lapply(split(seq_len(nrow(meta)), paste(meta$dataset, meta$group_binary)), function(ix)
      ix[sort(sample.int(length(ix), max(3L, ceiling(length(ix) * as.numeric(getarg("debug-unit-frac", "0.5"))))))]))
    meta <- meta[sort(keep_u)]
    ng <- as.integer(getarg("debug-genes-per-decile", "700"))
    keep_g <- unlist(lapply(1:10, function(d) { g <- dec[decile == d, gene]; g[sort(sample.int(length(g), ng))] }))
    other <- setdiff(rownames(counts), dec$gene)
    keep_o <- other[sort(sample.int(length(other), as.integer(getarg("debug-other-genes", "7000"))))]
    counts <- counts[rownames(counts) %in% c(keep_g, keep_o), meta$sample_id]
    dec <- dec[gene %in% keep_g]
    seed_base <- seed0 + 500000L
    say("DEBUG: %d units, %d genes (%d universe), seed base %d", nrow(meta), nrow(counts), nrow(dec), seed_base)
  }
  stopifnot(all(dec$gene %in% rownames(counts)), all(meta$sample_id %in% colnames(counts)))
  counts_m <- counts[, meta$sample_id, drop = FALSE]
  rm(counts); invisible(gc())
  lib0 <- colSums(counts_m)
  dec_genes <- lapply(1:10, function(d) dec[decile == d][order(gene), gene])
  dir.create(file.path(out, "replicates"), showWarnings = FALSE)
  ctx <- list(meta = meta, counts_m = counts_m, lib0 = lib0, dec_genes = dec_genes, cohorts = cohorts,
              seed_base = seed_base, n_units = nrow(meta), out = out)
  reps <- seq_len(n_rep)
  say("starting %d replicates on %d cores (seed base %d)", length(reps), cores, seed_base)
  res <- mclapply(reps, function(r) {
    o <- tryCatch(run_replicate(r, ctx), error = function(e) {
      msg <- conditionMessage(e)
      writeLines(msg, file.path(out, "replicates", sprintf("rep_%03d.error", r)))
      paste("ERROR", r, msg)
    })
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), o))
    o
  }, mc.cores = cores, mc.preschedule = FALSE)
  for (x in res) say("%s", if (inherits(x, "try-error")) paste("WORKER FAILURE", as.character(x)) else x)
  n_done <- length(list.files(file.path(out, "replicates"), pattern = "\\.done$"))
  say("replicates complete: %d of %d", n_done, length(reps))
  writeLines(capture.output(sessionInfo()), file.path(out, sprintf("sessionInfo_%s.txt", mode)))
  say("DONE %s", mode)
}
