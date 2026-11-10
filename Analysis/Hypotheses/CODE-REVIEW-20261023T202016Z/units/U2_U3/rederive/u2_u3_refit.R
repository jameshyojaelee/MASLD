# U2 (stage and NAS contrasts; Fig 3E, 3F, S3D-F) and U3 (per-cohort fits, leave-cohort-out,
# endpoint PCA; Fig 3C, 3D, S3B, S3C) independent refits.
# Written from docs/manuscript/METHODS.md ("Disease-versus-control differential expression",
# "Per-cohort disease models and GSE126848", "Fibrosis-stage and NAS contrasts") and the
# Figure 3 / S3 legends. Producer code was not copied.
# Unit: one selected library per genetic-identity cluster (774 clusters), not verified donors.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(parallel)
})
set.seed(20261023)  # no random draws are made; set for the record
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]; NCORES <- as.integer(args[2])
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
ARM <- file.path(R, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
IN <- c(
  counts_raw = file.path(R, "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds"),
  meta = file.path(ARM, "meta_matched.rds"),
  prod_dge = file.path(ARM, "merged_dge.rds"),
  pooled_deg = file.path(ARM, "deg_results.csv"),
  coloc_genes = file.path(R, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"),
  coloc_pairs = file.path(R, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz")
)
sha <- function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), " +")[[1]][1]
fwrite(data.table(input = names(IN), path = IN, sha256 = vapply(IN, sha, "")), file.path(OUT, "input_sha256.tsv"), sep = "\t")
logf <- file(file.path(OUT, "log.txt"), "w")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, logf); flush(logf) }
cmp_rows <- list()
record <- function(item, printed, mine, note = "") {
  cmp_rows[[length(cmp_rows) + 1L]] <<- data.table(item = item, printed = as.character(printed), mine = as.character(mine),
                                                  match = identical(as.character(printed), as.character(mine)), note = note)
  say(sprintf("[%s] %s | printed %s | mine %s %s", if (identical(as.character(printed), as.character(mine))) "OK" else "DIFF",
              item, printed, mine, note))
}
base_id <- function(x) sub("\\.[0-9]+$", "", x)

# ------------------------------------------------------------------ 0. pooled DGE from raw counts
raw <- readRDS(IN[["counts_raw"]])
if (is.data.frame(raw)) raw <- as.matrix(raw)
meta <- as.data.table(readRDS(IN[["meta"]]))
stopifnot(setequal(meta$sample_id, colnames(raw)), !anyDuplicated(meta$sample_id), !anyDuplicated(meta$analysis_unit_id))
meta <- meta[match(colnames(raw), sample_id)]
say("raw counts ", nrow(raw), " x ", ncol(raw), "; units ", nrow(meta), " control ", sum(meta$group_binary == "Control"),
    " disease ", sum(meta$group_binary == "Disease"))
gb <- factor(meta$group_binary, levels = c("Control", "Disease"))
y_raw <- DGEList(raw)
keep_pooled <- filterByExpr(y_raw, design = model.matrix(~ 0 + gb))
y <- y_raw[keep_pooled, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "RLE")
say("pooled universe after filterByExpr: ", nrow(y))
prod <- readRDS(IN[["prod_dge"]])
same_genes <- setequal(rownames(prod), rownames(y)) && nrow(prod) == nrow(y)
pm <- match(colnames(y), colnames(prod))
say("producer merged_dge: genes identical set ", same_genes, "; counts identical ",
    if (same_genes) identical(unname(prod$counts[rownames(y), pm]), unname(y$counts)) else NA,
    "; max |lib.size diff| ", max(abs(prod$samples$lib.size[pm] - y$samples$lib.size)),
    "; max |norm.factor diff| ", signif(max(abs(prod$samples$norm.factors[pm] - y$samples$norm.factors)), 3))
record("pooled tested family (genes)", 24115, nrow(y))
pooled <- fread(IN[["pooled_deg"]])
pooled[, call := padj < 0.05 & abs(logFC) > 0.5]
stopifnot(setequal(pooled$gene, rownames(y)))
pooled <- pooled[match(rownames(y), gene)]
say("pooled canonical calls in deg_results.csv: ", sum(pooled$call))

sym <- pooled$symbol

say("version stripping of the 24,115 IDs creates duplicates: ", anyDuplicated(base_id(rownames(y))) > 0)

# ------------------------------------------------------------------ fitting helpers
fit_core <- function(yy, design, coef, nf = "keep") {
  if (nf == "RLE") yy <- calcNormFactors(yy, method = "RLE")
  if (nf == "TMM") yy <- calcNormFactors(yy, method = "TMM")
  stopifnot(qr(design)$rank == ncol(design), coef %in% colnames(design))
  v <- voomWithQualityWeights(yy, design, plot = FALSE)
  f <- eBayes(lmFit(v, design))
  j <- match(coef, colnames(design))
  p <- f$p.value[, j]
  list(tab = data.table(gene = rownames(yy), logFC = f$coefficients[, j], P = p, FDR = p.adjust(p, "BH"),
                        SE = f$stdev.unscaled[, j] * sqrt(f$s2.post)),
       v = v, nf = yy$samples$norm.factors)
}
stage_units <- meta[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & fibrosis_stage %in% 0:4]
stage_units[, fib := paste0("F", fibrosis_stage)]
record("Fig 3B staged units F0/F1/F2/F3/F4", "103/102/88/69/25",
       paste(stage_units[, .N, keyby = fib]$N, collapse = "/"))
stage_units[, ctrl := group_binary == "Control"]
nas_units <- meta[!is.na(nas_score)]
nas_units[, nasg := fcase(nas_score == 0, "NAS0", nas_score <= 2, "NAS1-2", nas_score <= 4, "NAS3-4", nas_score <= 8, "NAS5-8")]
say("NAS units by cohort: ", paste(nas_units[, .N, keyby = dataset][, paste(dataset, N)], collapse = "; "))

# A contrast spec: units table, grouping column, ref, cmp, extra term, nf mode
fit_group <- function(units, gcol, ref, cmp, ctrl_term = FALSE, nf = "keep", interaction = FALSE) {
  u <- units[get(gcol) %in% c(ref, cmp)]
  u[, grp := factor(get(gcol), levels = c(ref, cmp))]
  info <- data.frame(dataset = droplevels(factor(u$dataset)), sex = droplevels(factor(u$inferred_sex)),
                     grp = u$grp, ctrl = factor(u$group_binary, levels = c("Disease", "Control")))
  idx <- match(u$sample_id, colnames(y))
  yy <- y[, idx]
  multi <- nlevels(info$dataset) > 1
  if (interaction) {
    design <- model.matrix(~ sex + dataset + dataset:grp, info)
    v <- voomWithQualityWeights(if (nf == "keep") yy else calcNormFactors(yy, method = nf), design, plot = FALSE)
    f <- lmFit(v, design)
    icols <- grep(":grp", colnames(design), value = TRUE)
    stopifnot(length(icols) == 3, qr(design)$rank == ncol(design))
    L <- matrix(0, ncol(design), 2, dimnames = list(colnames(design), c("d2", "d3")))
    L[icols[2], "d2"] <- 1; L[icols[1], "d2"] <- -1
    L[icols[3], "d3"] <- 1; L[icols[1], "d3"] <- -1
    fc <- eBayes(contrasts.fit(f, L))
    fP <- fc$F.p.value
    coefs <- f$coefficients[, icols]
    return(list(tab = data.table(gene = rownames(yy), F_P = fP, F_FDR = p.adjust(fP, "BH")), coefs = coefs,
                n = nrow(u)))
  }
  form <- if (multi) {
    if (ctrl_term) ~ dataset + sex + ctrl + grp else ~ dataset + sex + grp
  } else if (nlevels(info$sex) > 1) ~ sex + grp else ~ grp
  design <- model.matrix(form, info)
  r <- fit_core(yy, design, paste0("grp", cmp), nf)
  list(tab = r$tab, n_ref = sum(u$grp == ref), n_cmp = sum(u$grp == cmp), n = nrow(u), cohorts = nlevels(info$dataset),
       design = paste(deparse(form), collapse = ""))
}

# ------------------------------------------------------------------ job list
jobs <- list()
add <- function(name, fn) jobs[[name]] <<- fn
fib_pairs <- list(F0_to_F1 = c("F0", "F1"), F1_to_F2 = c("F1", "F2"), F2_to_F3 = c("F2", "F3"), F3_to_F4 = c("F3", "F4"),
                  F2_vs_F0 = c("F0", "F2"), F3_vs_F0 = c("F0", "F3"), F4_vs_F0 = c("F0", "F4"))
cohorts3 <- c("GSE130970", "GSE135251", "GSE162694")
for (k in names(fib_pairs)) local({
  kk <- k; pr <- fib_pairs[[k]]
  add(paste0("base|", kk), function() fit_group(stage_units, "fib", pr[1], pr[2]))
  add(paste0("renorm|", kk), function() fit_group(stage_units, "fib", pr[1], pr[2], nf = "RLE"))
  for (d in cohorts3) local({
    dd <- d
    add(paste0("drop_", dd, "|", kk), function() fit_group(stage_units[dataset != dd], "fib", pr[1], pr[2]))
  })
  if (pr[1] == "F0") {
    add(paste0("armB|", kk), function() fit_group(stage_units[!(fib == "F0" & ctrl)], "fib", pr[1], pr[2]))
    add(paste0("armC|", kk), function() fit_group(stage_units, "fib", pr[1], pr[2], ctrl_term = TRUE))
  }
})
for (d in cohorts3) local({
  dd <- d
  add(paste0("alone_", dd, "|F0_to_F1"), function() fit_group(stage_units[dataset == dd], "fib", "F0", "F1"))
})
add("interaction|F0_to_F1", function() fit_group(stage_units, "fib", "F0", "F1", interaction = TRUE))
add("interaction|F4_vs_F0", function() fit_group(stage_units, "fib", "F0", "F4", interaction = TRUE))
nas_pairs <- list(NAS12_vs_NAS0 = c("NAS0", "NAS1-2"), NAS34_vs_NAS0 = c("NAS0", "NAS3-4"), NAS58_vs_NAS0 = c("NAS0", "NAS5-8"),
                  NAS12_to_NAS34 = c("NAS1-2", "NAS3-4"), NAS34_to_NAS58 = c("NAS3-4", "NAS5-8"))
for (k in names(nas_pairs)) local({
  kk <- k; pr <- nas_pairs[[k]]
  add(paste0("base|", kk), function() fit_group(nas_units, "nasg", pr[1], pr[2]))
  add(paste0("renorm|", kk), function() fit_group(nas_units, "nasg", pr[1], pr[2], nf = "RLE"))
})

# U3 per-cohort fits: pooled universe and pooled RLE factors (Methods), sex term only when sex varies
for (d in sort(unique(meta$dataset))) local({
  dd <- d
  add(paste0("cohort|", dd), function() {
    u <- meta[dataset == dd]
    info <- data.frame(sex = factor(u$inferred_sex), grp = factor(u$group_binary, levels = c("Control", "Disease")))
    design <- if (nlevels(info$sex) > 1) model.matrix(~ sex + grp, info) else model.matrix(~ grp, info)
    r <- fit_core(y[, match(u$sample_id, colnames(y))], design, "grpDisease")
    list(tab = r$tab, n_ctrl = sum(info$grp == "Control"), n_dis = sum(info$grp == "Disease"), sexes = nlevels(info$sex),
         design = paste(colnames(design), collapse = "+"))
  })
})
# GSE126848 filtered refit: (a) filterByExpr on the pooled 24,115-gene subset; (b) on the raw count matrix
add("gse126848_filtered|a_from_24115", function() {
  u <- meta[dataset == "GSE126848"]; g <- factor(u$group_binary, levels = c("Control", "Disease"))
  yy <- y[, match(u$sample_id, colnames(y))]
  k <- filterByExpr(yy, group = g)
  yy <- calcNormFactors(yy[k, , keep.lib.sizes = FALSE], method = "RLE")
  r <- fit_core(yy, model.matrix(~ g), "gDisease")
  list(tab = r$tab, n_genes = nrow(yy))
})
add("gse126848_filtered|b_from_raw", function() {
  u <- meta[dataset == "GSE126848"]; g <- factor(u$group_binary, levels = c("Control", "Disease"))
  yy <- DGEList(raw[, u$sample_id])
  k <- filterByExpr(yy, group = g)
  yy <- calcNormFactors(yy[k, , keep.lib.sizes = FALSE], method = "RLE")
  r <- fit_core(yy, model.matrix(~ g), "gDisease")
  list(tab = r$tab, n_genes = nrow(yy), n_outside_24115 = sum(!rownames(yy) %in% rownames(y)))
})
# Fixed-universe leave-one-cohort-out: three normalization choices
for (d in sort(unique(meta$dataset))) for (nf in c("TMM", "RLE", "keep")) local({
  dd <- d; nn <- nf
  add(paste0("loco_fixed_", nn, "|", dd), function() {
    u <- meta[dataset != dd]
    info <- data.frame(dataset = factor(u$dataset), sex = factor(u$inferred_sex),
                       grp = factor(u$group_binary, levels = c("Control", "Disease")))
    r <- fit_core(y[, match(u$sample_id, colnames(y))], model.matrix(~ dataset + sex + grp, info), "grpDisease", nn)
    list(tab = r$tab, n = nrow(u))
  })
})
# Full leave-one-cohort-out: filter and normalization repeated after removal
for (d in c(sort(unique(meta$dataset)), "NONE")) local({
  dd <- d
  add(paste0("loco_full|", dd), function() {
    u <- meta[dataset != dd]
    g <- factor(u$group_binary, levels = c("Control", "Disease"))
    yy <- DGEList(raw[, u$sample_id])
    k <- filterByExpr(yy, design = model.matrix(~ 0 + g))
    yy <- calcNormFactors(yy[k, , keep.lib.sizes = FALSE], method = "RLE")
    info <- data.frame(dataset = factor(u$dataset), sex = factor(u$inferred_sex), grp = g)
    r <- fit_core(yy, model.matrix(~ dataset + sex + grp, info), "grpDisease")
    list(tab = r$tab, n = nrow(u), n_genes = nrow(yy))
  })
})
# Fig 3D / S3B endpoint
meta[, advanced := (!is.na(nas_score) & nas_score >= 5) | (!is.na(fibrosis_stage) & fibrosis_stage >= 3)]
meta[, src_ctrl := fifelse(!is.na(source_control_status), source_control_status == "Control", group_binary == "Control")]
meta[, strict := src_ctrl & (is.na(nas_score) | nas_score <= 1) & (is.na(fibrosis_stage) | fibrosis_stage <= 0)]
meta[, endpoint := fifelse(advanced, "Advanced", fifelse(strict, "Strict", NA_character_))]
ep_tab <- dcast(meta[!is.na(endpoint), .N, by = .(dataset, endpoint)], dataset ~ endpoint, value.var = "N", fill = 0L)
ep_cohorts <- ep_tab[Advanced >= 5 & Strict >= 5, dataset]
ep <- meta[dataset %in% ep_cohorts & !is.na(endpoint)]
for (nf in c("TMM", "keep")) local({
  nn <- nf
  add(paste0("endpoint|", nn), function() {
    info <- data.frame(dataset = factor(ep$dataset), sex = factor(ep$inferred_sex),
                       grp = factor(ep$endpoint, levels = c("Strict", "Advanced")))
    design <- model.matrix(~ dataset + sex + grp, info)
    r <- fit_core(y[, match(ep$sample_id, colnames(y))], design, "grpAdvanced", nn)
    adj <- removeBatchEffect(r$v$E, batch = info$dataset, covariates = model.matrix(~ sex, info)[, -1, drop = FALSE])
    vv <- apply(adj, 1, var)
    top <- order(vv, decreasing = TRUE)[1:2000]
    pc <- prcomp(t(adj[top, ]), center = TRUE, scale. = FALSE)
    list(tab = r$tab, pve = 100 * pc$sdev^2 / sum(pc$sdev^2))
  })
})

say("fits to run: ", length(jobs), " on ", NCORES, " cores")
t0 <- Sys.time()
res <- mclapply(names(jobs), function(n) tryCatch(jobs[[n]](), error = function(e) list(error = conditionMessage(e))),
                mc.cores = NCORES, mc.preschedule = FALSE)
names(res) <- names(jobs)
say("fits done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
errs <- names(res)[vapply(res, function(r) !is.null(r$error), TRUE)]
if (length(errs)) for (e in errs) say("ERROR ", e, ": ", res[[e]]$error)
saveRDS(res, file.path(OUT, "all_fits.rds"))
nsig <- function(r) sum(r$tab$FDR < 0.05)
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# ------------------------------------------------------------------ U2 results
counts_tab <- rbindlist(lapply(grep("^(base|renorm|drop_|armB|armC|alone_)", names(res), value = TRUE), function(n) {
  r <- res[[n]]; if (!is.null(r$error)) return(NULL)
  p <- strsplit(n, "\\|")[[1]]
  data.table(fit = p[1], contrast = p[2], n_ref = r$n_ref, n_cmp = r$n_cmp, n = r$n, cohorts = r$cohorts,
             family = nrow(r$tab), bh05 = nsig(r), up = sum(r$tab$FDR < 0.05 & r$tab$logFC > 0),
             down = sum(r$tab$FDR < 0.05 & r$tab$logFC < 0), bh05_lfc05 = sum(r$tab$FDR < 0.05 & abs(r$tab$logFC) > 0.5),
             design = r$design)
}))
fwrite(counts_tab, file.path(OUT, "u2_stage_counts.tsv"), sep = "\t")
gc_ <- function(f, con) { x <- counts_tab[counts_tab$fit == f & counts_tab$contrast == con]; if (nrow(x)) x$bh05 else NA }
gn_ <- function(f, con) { x <- counts_tab[counts_tab$fit == f & counts_tab$contrast == con]; if (nrow(x)) x$n else NA }
record("Fig 3E F1/F2/F3/F4 vs F0 (BH<0.05, 24,115 family)", "1,338/4,266/6,190/7,043",
       paste(fmt(c(gc_("base", "F0_to_F1"), gc_("base", "F2_vs_F0"), gc_("base", "F3_vs_F0"), gc_("base", "F4_vs_F0"))), collapse = "/"))
record("Fig 3E unit n F0..F4", "103/102/88/69/25", paste(stage_units[, .N, keyby = fib]$N, collapse = "/"))
record("Fig 3E F0 composition (controls/disease)", "14/89",
       paste0(stage_units[fib == "F0", sum(ctrl)], "/", stage_units[fib == "F0", sum(!ctrl)]))
record("S3E disease-only F0 F1..F4", "783/3,946/6,459/6,929",
       paste(fmt(c(gc_("armB", "F0_to_F1"), gc_("armB", "F2_vs_F0"), gc_("armB", "F3_vs_F0"), gc_("armB", "F4_vs_F0"))), collapse = "/"))
record("S3E control-status term F1..F4", "814/3,679/6,166/6,640",
       paste(fmt(c(gc_("armC", "F0_to_F1"), gc_("armC", "F2_vs_F0"), gc_("armC", "F3_vs_F0"), gc_("armC", "F4_vs_F0"))), collapse = "/"))
record("S3E disease-only F0 units", "89", counts_tab[fit == "armB" & contrast == "F0_to_F1", n_ref])
record("S3F three-cohort F1 vs F0 (units)", "1,338 (205)", sprintf("%s (%s)", fmt(gc_("base", "F0_to_F1")), gn_("base", "F0_to_F1")))
record("S3F F1 vs F0 without GSE130970/GSE135251/GSE162694", "14/929/651",
       paste(c(gc_("drop_GSE130970", "F0_to_F1"), gc_("drop_GSE135251", "F0_to_F1"), gc_("drop_GSE162694", "F0_to_F1")), collapse = "/"))
record("S3F leave-out units", "155/115/140",
       paste(c(gn_("drop_GSE130970", "F0_to_F1"), gn_("drop_GSE135251", "F0_to_F1"), gn_("drop_GSE162694", "F0_to_F1")), collapse = "/"))
record("S3F each cohort alone GSE130970/GSE135251/GSE162694", "252/3/0",
       paste(c(gc_("alone_GSE130970", "F0_to_F1"), gc_("alone_GSE135251", "F0_to_F1"), gc_("alone_GSE162694", "F0_to_F1")), collapse = "/"))
record("S3F alone units", "50/90/65",
       paste(c(gn_("alone_GSE130970", "F0_to_F1"), gn_("alone_GSE135251", "F0_to_F1"), gn_("alone_GSE162694", "F0_to_F1")), collapse = "/"))
record("S3D NAS1-2/3-4/5-8 vs NAS0", "5,338/9,337/9,841",
       paste(fmt(c(gc_("base", "NAS12_vs_NAS0"), gc_("base", "NAS34_vs_NAS0"), gc_("base", "NAS58_vs_NAS0"))), collapse = "/"))
record("S3D NAS0 units (controls + disease)", "46 (45+1)",
       paste0(nas_units[nasg == "NAS0", .N], " (", nas_units[nasg == "NAS0", sum(group_binary == "Control")], "+",
              nas_units[nasg == "NAS0", sum(group_binary == "Disease")], ")"))
record("S3D NAS1-2->3-4 and 3-4->5-8", "3 (65 and 117 units); 4,305 (117 and 167 units)",
       sprintf("%s (%s and %s units); %s (%s and %s units)", fmt(gc_("base", "NAS12_to_NAS34")),
               counts_tab[fit == "base" & contrast == "NAS12_to_NAS34", n_ref], counts_tab[fit == "base" & contrast == "NAS12_to_NAS34", n_cmp],
               fmt(gc_("base", "NAS34_to_NAS58")),
               counts_tab[fit == "base" & contrast == "NAS34_to_NAS58", n_ref], counts_tab[fit == "base" & contrast == "NAS34_to_NAS58", n_cmp]))
record("draft: F2-F4 vs F0 range under single-cohort omission", "1,934 to 9,546",
       paste(fmt(range(counts_tab[grepl("^drop_", fit) & contrast %in% c("F2_vs_F0", "F3_vs_F0", "F4_vs_F0"), bh05])), collapse = " to "))
record("draft: F1->F2 without GSE162694 (units)", "0 (133)",
       sprintf("%s (%s)", gc_("drop_GSE162694", "F1_to_F2"), gn_("drop_GSE162694", "F1_to_F2")))
record("draft: F2->F3 without GSE135251 (units)", "2 (54)",
       sprintf("%s (%s)", gc_("drop_GSE135251", "F2_to_F3"), gn_("drop_GSE135251", "F2_to_F3")))
say("Normalization-reuse sensitivity (RLE recomputed on each contrast subset vs pooled RLE factors):")
nr <- merge(counts_tab[fit == "base", .(contrast, pooled_nf = bh05)], counts_tab[fit == "renorm", .(contrast, subset_RLE = bh05)], by = "contrast")
for (i in seq_len(nrow(nr))) say(sprintf("  %s: pooled factors %d, recomputed RLE %d", nr$contrast[i], nr$pooled_nf[i], nr$subset_RLE[i]))
fwrite(nr, file.path(OUT, "u2_normalization_reuse_sensitivity.tsv"), sep = "\t")

# interaction models
for (k in c("F0_to_F1", "F4_vs_F0")) {
  ir <- res[[paste0("interaction|", k)]]
  if (!is.null(ir$error)) next
  base_tab <- res[[paste0("base|", k)]]$tab
  calls <- base_tab$FDR < 0.05
  same_all <- apply(sign(ir$coefs[calls, , drop = FALSE]) == sign(base_tab$logFC[calls]), 1, all)
  het <- sum(ir$tab$F_FDR[calls] < 0.05)
  say(sprintf("interaction %s: calls %d; pooled sign in all three cohort coefficients %.1f%%; interaction BH<0.05 among calls %d",
              k, sum(calls), 100 * mean(same_all), het))
  if (k == "F0_to_F1") record("draft: F1 calls with pooled sign in all three cohorts", "95.8%", sprintf("%.1f%%", 100 * mean(same_all)))
  if (k == "F4_vs_F0") {
    record("draft: F4 vs F0 calls with cohort-by-stage interaction BH<0.05", "4,047", fmt(het))
    say(sprintf("  F4 pooled sign in all three cohorts %.1f%% (RESULTS.md cross-check 75.4%%)", 100 * mean(same_all)))
  }
}

# hero genes
heroes <- c("GNMT", "MAT1A", "CYP2C19")
hq <- rbindlist(lapply(grep("^(base|armB|armC|drop_GSE[0-9]+)\\|(F2|F3|F4)_vs_F0$", names(res), value = TRUE), function(n) {
  t <- res[[n]]$tab; i <- match(heroes, sym[match(t$gene, rownames(y))])
  data.table(fit = n, gene = heroes, logFC = t$logFC[i], FDR = t$FDR[i])
}))
fwrite(hq, file.path(OUT, "u2_hero_genes.tsv"), sep = "\t")
say(sprintf("hero genes F2-F4 vs F0: all negative %s; max q (base, armB, armC) %.2g; max q without GSE130970 %.2g; max q without GSE135251 %.3f; without GSE162694 %.3f",
            all(hq$logFC < 0), max(hq[grepl("^(base|armB|armC)", fit), FDR]), max(hq[grepl("drop_GSE130970", fit), FDR]),
            max(hq[grepl("drop_GSE135251", fit), FDR]), max(hq[grepl("drop_GSE162694", fit), FDR])))
say("  worst: ", paste(hq[order(-FDR)][1:3, sprintf("%s %s q=%.4f", gene, fit, FDR)], collapse = "; "))

# ------------------------------------------------------------------ Fig 3F
cg <- fread(IN[["coloc_genes"]])
stopifnot(!anyDuplicated(cg$ensembl))
cg[, ensembl := base_id(ensembl)]
susie <- cg[!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > 0.5]
record("primary SuSiE set (tier-1/2 PP.H4>0.5)", "436", nrow(susie))
adj <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
W <- data.table(gene = rownames(y), ens = base_id(rownames(y)), symbol = sym)
for (k in adj) { t <- res[[paste0("base|", k)]]$tab; W[[paste0("q_", k)]] <- t$FDR; W[[paste0("lfc_", k)]] <- t$logFC }
say("SuSiE genes in the 24,115 stage universe: ", sum(susie$ensembl %in% W$ens), " of ", nrow(susie))
S <- W[ens %in% susie$ensembl]
qm <- as.matrix(S[, paste0("q_", adj), with = FALSE]); lm_ <- as.matrix(S[, paste0("lfc_", adj), with = FALSE])
sel <- rowSums(qm < 0.05) > 0
S <- S[sel]; qm <- qm[sel, ]; lm_ <- lm_[sel, ]
record("Fig 3F genes", "54", nrow(S))
record("Fig 3F effects", "216", length(qm))
record("Fig 3F BH marks total (F0F1/F1F2/F2F3/F3F4)", "59 (30/13/11/5)",
       sprintf("%d (%s)", sum(qm < 0.05), paste(colSums(qm < 0.05), collapse = "/")))
record("Fig 3F genes with any q<0.0125 (marks per contrast)", "15 (11/2/2/1)",
       sprintf("%d (%s)", sum(rowSums(qm < 0.0125) > 0), paste(colSums(qm < 0.0125), collapse = "/")))
only01 <- qm[, 1] < 0.05 & rowSums(qm[, 2:4] < 0.05) == 0
record("Fig 3F genes selected only via F0->F1", "29", sum(only01))
record("Fig 3F colour clip (98th pct |log2FC| of displayed effects)", "0.688", sprintf("%.3f", quantile(abs(lm_), 0.98)))
mk <- function(fit, k) { t <- res[[paste0(fit, "|", k)]]$tab; m <- match(S$gene, t$gene); list(q = t$FDR[m], l = t$logFC[m]) }
m01 <- qm[, 1] < 0.05
b <- mk("armB", "F0_to_F1")
record("F0->F1 marks kept (q<0.05, same sign) with disease-only F0", "18 of 30",
       sprintf("%d of %d", sum(m01 & b$q < 0.05 & sign(b$l) == sign(lm_[, 1])), sum(m01)))
d1 <- mk("drop_GSE130970", "F0_to_F1")
record("F0->F1 marks kept without GSE130970", "0 of 30",
       sprintf("%d of %d (min q %.3f)", sum(m01 & d1$q < 0.05 & sign(d1$l) == sign(lm_[, 1])), sum(m01), min(d1$q[m01])))
d2 <- mk("drop_GSE162694", "F1_to_F2")
m12 <- qm[, 2] < 0.05
record("F1->F2 marks kept without GSE162694", "0 of 13", sprintf("%d of %d", sum(m12 & d2$q < 0.05 & sign(d2$l) == sign(lm_[, 2])), sum(m12)))
keep_all <- matrix(TRUE, nrow(qm), 4)
for (j in 1:4) for (d in cohorts3) { x <- mk(paste0("drop_", d), adj[j]); keep_all[, j] <- keep_all[, j] & x$q < 0.05 & sign(x$l) == sign(lm_[, j]) }
ka <- which(qm < 0.05 & keep_all, arr.ind = TRUE)
record("marks kept under every single-cohort omission", "1 (OMG, F2->F3)",
       sprintf("%d (%s)", nrow(ka), paste(S$symbol[ka[, 1]], adj[ka[, 2]], collapse = "; ")))
# PP.H4 lane: gene-level best tier-1/2 SuSiE PP.H4 and its study; check the exact eligible pair carries the same value
pairs <- fread(IN[["coloc_pairs"]])
pairs[, ensembl := base_id(ensembl)]
say("signal-pair eligible column class: ", class(pairs$eligible))
lane <- susie[ensembl %in% S$ens, .(ensembl, pp4 = coloc_best_susie_pp4, study = susie_driving_gwas, trait = driving_trait,
                                    susie_trait = driving_trait, abf_trait = abf_driving_trait)]
best_pair <- pairs[eligible == TRUE][lane, on = .(ensembl, gwas = study), nomatch = 0][, .(pair_max = max(PP.H4)), by = ensembl]
lane <- merge(lane, best_pair, by = "ensembl", all.x = TRUE)
say("Fig 3F lane: genes ", nrow(lane), "; exact-pair max PP.H4 equals gene-level best for all: ",
    isTRUE(all(abs(lane$pair_max - lane$pp4) < 1e-8)), "; missing pair ", sum(is.na(lane$pair_max)))
S[, only_F0F1 := only01]
fwrite(merge(S, lane, by.x = "ens", by.y = "ensembl", all.x = TRUE), file.path(OUT, "u2_fig3f_genes.tsv"), sep = "\t")

# ------------------------------------------------------------------ U3 per-cohort and Fig 3C
pc <- pooled[call == TRUE]
cohort_calls <- rbindlist(lapply(grep("^cohort\\|", names(res), value = TRUE), function(n) {
  r <- res[[n]]; t <- r$tab
  data.table(dataset = sub("cohort\\|", "", n), gene = t$gene, logFC = t$logFC, FDR = t$FDR, call = t$FDR < 0.05 & abs(t$logFC) > 0.5)
}))
cs <- cohort_calls[, .(calls = sum(call), up = sum(call & logFC > 0), down = sum(call & logFC < 0),
                       dir_agree = mean(sign(logFC[match(pc$gene, gene)]) == sign(pc$logFC)),
                       family = .N), by = dataset]
fwrite(cs, file.path(OUT, "u3_cohort_summary.tsv"), sep = "\t")
for (i in seq_len(nrow(cs))) say(sprintf("cohort %s: design %s; n %d/%d; calls %d (%d up, %d down); family %d; direction agreement %.4f",
                                         cs$dataset[i], res[[paste0("cohort|", cs$dataset[i])]]$design, res[[paste0("cohort|", cs$dataset[i])]]$n_ctrl,
                                         res[[paste0("cohort|", cs$dataset[i])]]$n_dis, cs$calls[i], cs$up[i], cs$down[i], cs$family[i], cs$dir_agree[i]))
pr <- function(d) cs[dataset == d]
record("Fig 3C GSE135251 calls", "7,381 (2,578 up, 4,803 down)", sprintf("%s (%s up, %s down)", fmt(pr("GSE135251")$calls), fmt(pr("GSE135251")$up), fmt(pr("GSE135251")$down)))
record("Fig 3C GSE162694 calls", "3,199 (2,518 up, 681 down)", sprintf("%s (%s up, %s down)", fmt(pr("GSE162694")$calls), fmt(pr("GSE162694")$up), fmt(pr("GSE162694")$down)))
record("Fig 3C GSE213621 calls", "2,045 (1,409 up, 636 down)", sprintf("%s (%s up, %s down)", fmt(pr("GSE213621")$calls), fmt(pr("GSE213621")$up), fmt(pr("GSE213621")$down)))
record("Fig 3C GSE130970 calls", "1,282 (349 up, 933 down)", sprintf("%s (%s up, %s down)", fmt(pr("GSE130970")$calls), fmt(pr("GSE130970")$up), fmt(pr("GSE130970")$down)))
record("draft: GSE126848 calls on fixed 24,115 family", "6,147", fmt(pr("GSE126848")$calls))
record("S3C GSE126848 direction agreement", "84.36%", sprintf("%.2f%%", 100 * pr("GSE126848")$dir_agree))
record("S3C direction agreement range (5 cohorts)", "78.3% to 100%", sprintf("%.1f%% to %.1f%%", 100 * min(cs$dir_agree), 100 * max(cs$dir_agree)))
fa <- res[["gse126848_filtered|a_from_24115"]]; fb <- res[["gse126848_filtered|b_from_raw"]]
callsA <- fa$tab[FDR < 0.05 & abs(logFC) > 0.5]
callsB <- fb$tab[FDR < 0.05 & abs(logFC) > 0.5]
record("Fig 3C GSE126848 filtered refit genes and calls", "15,592; 3,440 (1,821 up, 1,619 down)",
       sprintf("%s; %s (%s up, %s down)", fmt(fa$n_genes), fmt(nrow(callsA)), fmt(sum(callsA$logFC > 0)), fmt(sum(callsA$logFC < 0))))
say(sprintf("GSE126848 filterByExpr from the raw matrix instead: %s genes (%d outside the 24,115), %s calls (%s up, %s down)",
            fmt(fb$n_genes), fb$n_outside_24115, fmt(nrow(callsB)), fmt(sum(callsB$logFC > 0)), fmt(sum(callsB$logFC < 0))))
classify <- function(cl) {
  g <- cl[, .(k = .N, nup = sum(sgn > 0), ndn = sum(sgn < 0), has126 = any(dataset == "GSE126848")), by = gene]
  g[, cls := fifelse(k == 1L, "unique", fifelse(nup == 0L | ndn == 0L, paste0("same", k), "mixed"))]
  g
}
four <- cohort_calls[call == TRUE & dataset != "GSE126848", .(gene, dataset, sgn = sign(logFC))]
five_f <- rbind(four, callsA[, .(gene, dataset = "GSE126848", sgn = sign(logFC))])
five_a <- rbind(four, cohort_calls[call == TRUE & dataset == "GSE126848", .(gene, dataset, sgn = sign(logFC))])
g5 <- classify(five_f); g4 <- classify(four); g5a <- classify(five_a)
cls_n <- function(g, c) sum(g$cls == c)
record("Fig 3C memberships", "17,347", fmt(nrow(five_f)))
record("Fig 3C union", "12,180", fmt(nrow(g5)))
record("Fig 3C classes one/same2/same3/same4/same5/mixed", "8,341/1,820/582/186/40/1,211",
       paste(fmt(sapply(c("unique", "same2", "same3", "same4", "same5", "mixed"), function(c) cls_n(g5, c))), collapse = "/"))
record("Fig 3C with GSE126848 call per class", "1,550/493/278/133/40/946",
       paste(fmt(sapply(c("unique", "same2", "same3", "same4", "same5", "mixed"), function(c) sum(g5$cls == c & g5$has126))), collapse = "/"))
record("Fig 3C mixed share with GSE126848", "1,211/3,839 (31.5%)",
       sprintf("%s/%s (%.1f%%)", fmt(cls_n(g5, "mixed")), fmt(sum(g5$k >= 2)), 100 * cls_n(g5, "mixed") / sum(g5$k >= 2)))
record("Fig 3C mixed share without GSE126848", "346/2,565 (13.5%)",
       sprintf("%s/%s (%.1f%%)", fmt(cls_n(g4, "mixed")), fmt(sum(g4$k >= 2)), 100 * cls_n(g4, "mixed") / sum(g4$k >= 2)))
record("draft: mixed share with 6,147 GSE126848 calls", "1,865 of 4,696 (39.7%)",
       sprintf("%s of %s (%.1f%%)", fmt(cls_n(g5a, "mixed")), fmt(sum(g5a$k >= 2)), 100 * cls_n(g5a, "mixed") / sum(g5a$k >= 2)))
# raw-matrix filter variant of the union
g5b <- classify(rbind(four, callsB[, .(gene, dataset = "GSE126848", sgn = sign(logFC))]))
say(sprintf("Union with raw-matrix-filtered GSE126848: union %s, mixed %s/%s (%.1f%%)", fmt(nrow(g5b)), fmt(cls_n(g5b, "mixed")),
            fmt(sum(g5b$k >= 2)), 100 * cls_n(g5b, "mixed") / sum(g5b$k >= 2)))
# single-cohort same-sign calls of the pooled set (draft text)
alone_same <- cohort_calls[call == TRUE][pc, on = "gene", nomatch = 0][sign(logFC) == sign(i.logFC), .N, by = dataset]
say("pooled calls also called alone with same sign: ", paste(alone_same[, paste(dataset, N)], collapse = "; "))
record("draft: GSE213621 alone calls of 1,854 same sign", "1,654", alone_same[dataset == "GSE213621", N])
record("draft: GSE162694 alone calls of 1,854 same sign", "1,014", alone_same[dataset == "GSE162694", N])

# ------------------------------------------------------------------ U3 leave-one-cohort-out
lf <- rbindlist(lapply(grep("^loco_fixed_", names(res), value = TRUE), function(n) {
  p <- strsplit(sub("loco_fixed_", "", n), "\\|")[[1]]; t <- res[[n]]$tab
  m <- match(pooled$gene, t$gene); lc <- t$FDR[m] < 0.05 & abs(t$logFC[m]) > 0.5
  data.table(norm = p[1], left_out = p[2], n = res[[n]]$n, rho = cor(t$logFC[m], pooled$logFC, method = "spearman"),
             retained = sum(lc & pooled$call), retained_same_sign = sum(lc & pooled$call & sign(t$logFC[m]) == sign(pooled$logFC)),
             retained_frac = sum(lc & pooled$call) / sum(pooled$call))
}))
fwrite(lf, file.path(OUT, "u3_loco_fixed.tsv"), sep = "\t")
for (nn in c("TMM", "RLE", "keep")) {
  x <- lf[norm == nn]
  say(sprintf("fixed-universe LOCO, normalization %s: rho %.3f to %.3f; retention %.1f%% to %.1f%% (%s)", nn, min(x$rho), max(x$rho),
              100 * min(x$retained_frac), 100 * max(x$retained_frac),
              paste(x[, sprintf("%s rho %.3f ret %d", left_out, rho, retained)], collapse = "; ")))
}
x <- lf[norm == "TMM"]
record("S3C fixed-universe LOCO rho (producer TMM)", "0.870-0.995", sprintf("%.3f-%.3f", min(x$rho), max(x$rho)))
record("S3C fixed-universe LOCO retention (producer TMM)", "70.7%-96.6%", sprintf("%.1f%%-%.1f%%", 100 * min(x$retained_frac), 100 * max(x$retained_frac)))
x <- lf[norm == "RLE"]
record("S3C fixed-universe LOCO rho (RLE as pooled fit)", "0.870-0.995", sprintf("%.3f-%.3f", min(x$rho), max(x$rho)), "normalization variant")
record("S3C fixed-universe LOCO retention (RLE as pooled fit)", "70.7%-96.6%", sprintf("%.1f%%-%.1f%%", 100 * min(x$retained_frac), 100 * max(x$retained_frac)), "normalization variant")
lfull <- rbindlist(lapply(grep("^loco_full\\|", names(res), value = TRUE), function(n) {
  t <- res[[n]]$tab; m <- match(pooled$gene, t$gene)
  tested <- !is.na(m); lc <- tested & t$FDR[m] < 0.05 & abs(t$logFC[m]) > 0.5
  same <- tested & sign(t$logFC[m]) == sign(pooled$logFC)
  data.table(left_out = sub("loco_full\\|", "", n), n = res[[n]]$n, genes = res[[n]]$n_genes, calls = sum(t$FDR < 0.05 & abs(t$logFC) > 0.5),
             retained_same_sign = sum(pooled$call & lc & same), lost = sum(pooled$call & !(lc & same)),
             lost_untested = sum(pooled$call & !tested),
             tested_sign_kept = sum(pooled$call & tested & same) / sum(pooled$call & tested),
             max_abs_dlogFC = max(abs(t$logFC[m] - pooled$logFC), na.rm = TRUE))
}))
lfull[, retained_pct := round(100 * retained_same_sign / sum(pooled$call), 1)]
fwrite(lfull, file.path(OUT, "u3_loco_full.tsv"), sep = "\t")
print(lfull)
record("pooled refit reproduces 1,854 (leave-none-out)", "1,854", lfull[left_out == "NONE", calls],
       sprintf("max |dlogFC| %.2g", lfull[left_out == "NONE", max_abs_dlogFC]))
x <- lfull[left_out != "NONE"]
record("S3C full LOCO same-sign retention", "66.7%-96.8%", sprintf("%.1f%%-%.1f%%", min(x$retained_pct), max(x$retained_pct)))
record("S3C full LOCO tested calls keeping sign", "99.7%-100%", sprintf("%.1f%%-%.1f%%", 100 * min(x$tested_sign_kept), 100 * max(x$tested_sign_kept)))
record("draft: removing GSE213621 loses (untested)", "618 (108)", sprintf("%d (%d)", x[left_out == "GSE213621", lost], x[left_out == "GSE213621", lost_untested]))

# ------------------------------------------------------------------ U3 Fig 3D / S3B
record("Fig 3D strict controls / advanced / total", "113/279/392",
       sprintf("%d/%d/%d", sum(ep$endpoint == "Strict"), sum(ep$endpoint == "Advanced"), nrow(ep)))
record("Fig 3D cohorts", "4", length(ep_cohorts))
adv <- ep[endpoint == "Advanced"]
record("Fig 3D advanced split F3+ / NAS>=5 with F0-F2", "171/108",
       sprintf("%d/%d", adv[!is.na(fibrosis_stage) & fibrosis_stage >= 3, .N], adv[!(!is.na(fibrosis_stage) & fibrosis_stage >= 3), .N]))
record("Fig 3D GSE213621 F3-F4 bin units", "77", adv[dataset == "GSE213621", .N])
record("Fig 3D exact stages among advanced F0/F1/F2/F3/F4", "7/42/59/69/25",
       paste(sapply(0:4, function(s) adv[dataset != "GSE213621" & fibrosis_stage == s, .N]), collapse = "/"))
say("advanced units with missing fibrosis: ", adv[is.na(fibrosis_stage), .N])
for (nn in c("TMM", "keep")) {
  r <- res[[paste0("endpoint|", nn)]]; t <- r$tab; m <- match(pooled$gene, t$gene)
  ec <- t$FDR[m] < 0.05 & abs(t$logFC[m]) > 0.5
  same <- sign(t$logFC[m]) == sign(pooled$logFC)
  lab <- if (nn == "TMM") "producer TMM" else "pooled RLE factors"
  record(paste0("S3B recovered/subthreshold/endpoint-only (", lab, ")"), "1,755/99/1,195",
         sprintf("%s/%s/%s", fmt(sum(pooled$call & ec & same)), fmt(sum(pooled$call & !(ec & same))), fmt(sum(!pooled$call & ec))),
         if (nn == "keep") "normalization variant" else "")
  record(paste0("S3B all-gene rho (", lab, ")"), "0.953", sprintf("%.3f", cor(t$logFC[m], pooled$logFC, method = "spearman")),
         if (nn == "keep") "normalization variant" else "")
  record(paste0("S3B opposite-direction double calls (", lab, ")"), "0", sum(pooled$call & ec & !same))
  record(paste0("Fig 3D PC1/PC2 variance (", lab, ")"), "9.8%/4.9%", sprintf("%.1f%%/%.1f%%", r$pve[1], r$pve[2]),
         if (nn == "keep") "normalization variant" else "")
}

cmp <- rbindlist(cmp_rows)
fwrite(cmp, file.path(OUT, "printed_vs_mine.tsv"), sep = "\t")
say("items ", nrow(cmp), "; match ", sum(cmp$match), "; differ ", sum(!cmp$match))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
writeLines(c("seed 20261023 (no random draws)", "B: none (no resampling in U2/U3)",
             "BH family: every gene in each fit (24,115 for fixed-universe fits; filtered universe for filtered refits)"),
           file.path(OUT, "run_parameters.txt"))
say("DONE")
close(logf)
