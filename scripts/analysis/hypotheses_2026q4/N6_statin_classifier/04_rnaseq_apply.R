#!/usr/bin/env Rscript
# N6 application step (run only if gates G1 and G2 pass). Cross-modality port: only the
# within-sample rank score of the array-derived gene set is carried to RNA-seq counts.
# Transfer accuracy cannot be measured (no medication data). Reports whether the score is
# associated with histological group and whether it changes the refit's disease calls.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(jsonlite); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 04_rnaseq_apply.R signature_K50_all_arrays.tsv refit_arm_dir out_dir")
sig_f <- args[1]; arm <- args[2]; out <- args[3]
if (dir.exists(out)) stop("Refusing to overwrite output dir")
dir.create(out, recursive = TRUE)
SEED <- 20260923L; NPERM <- 20L
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
set.seed(SEED)

dge <- readRDS(file.path(arm, "merged_dge.rds"))
meta <- readRDS(file.path(arm, "meta_matched.rds"))
deg0 <- fread(file.path(arm, "deg_results.csv"))
meta <- meta[match(colnames(dge), meta$sample_id)]
stopifnot(identical(meta$sample_id, colnames(dge)), !anyNA(meta$inferred_sex), ncol(dge) == 774L)
info <- data.frame(dataset = factor(meta$dataset), inferred_sex = factor(meta$inferred_sex),
                   group_binary = factor(meta$group_binary, levels = c("Control", "Disease")),
                   fibrosis_stage = suppressWarnings(as.numeric(meta$fibrosis_stage)), row.names = meta$sample_id)

# ---- signature coverage and rank score
sig <- fread(sig_f)
sym <- deg0$symbol[match(rownames(dge), deg0$gene)]
lcpm <- cpm(dge, log = TRUE, prior.count = 2)
gene_mean <- rowMeans(lcpm)
pick <- function(s) { ix <- which(sym == s); if (!length(ix)) NA_integer_ else ix[which.max(gene_mean[ix])] }
sig[, row := vapply(symbol, pick, 1L)]
sig[, present := !is.na(row)]
fwrite(sig, file.path(out, "signature_coverage.tsv"), sep = "\t")
frac <- mean(sig$present)
cat("signature genes present:", sum(sig$present), "of", nrow(sig), "\n")
if (frac < 0.8) stop("Fewer than 80% of signature genes present; stopping per PRESPEC")
RK <- apply(lcpm, 2, rank) / nrow(lcpm)
up <- sig[present & sign > 0, row]; dn <- sig[present & sign < 0, row]
mm <- function(ix) if (length(ix)) colMeans(RK[ix, , drop = FALSE]) else rep(0, ncol(RK))
raw <- mm(up) - mm(dn)
score <- as.numeric(scale(raw)); names(score) <- colnames(dge)
fwrite(data.table(sample_id = colnames(dge), dataset = meta$dataset, group_binary = meta$group_binary,
                  fibrosis_stage = info$fibrosis_stage, raw_rank_score = raw, score_z = score),
       file.path(out, "unit_scores.tsv"), sep = "\t")

# ---- association with histological group (units: 774 one-per-person roster)
tid <- function(fit, term) { s <- summary(fit)$coefficients[term, ]; ci <- confint(fit)[term, ]
  data.table(term = term, n = nobs(fit), estimate_sd_units = s[1], se = s[2], lo = ci[1], hi = ci[2], p_nominal = s[4]) }
d <- cbind(info, score = score)
a1 <- tid(lm(score ~ dataset + inferred_sex + group_binary, d), "group_binaryDisease"); a1[, test := "group_binary"]
ds <- d[!is.na(d$fibrosis_stage), ]
a2 <- tid(lm(score ~ dataset + inferred_sex + fibrosis_stage, ds), "fibrosis_stage"); a2[, test := "fibrosis_stage_ordinal"]
assoc <- rbind(a1, a2)
assoc[, p_holm_2 := p.adjust(p_nominal, "holm")]
fwrite(assoc, file.path(out, "association_score_vs_stage.tsv"), sep = "\t")
pc <- rbindlist(lapply(levels(d$dataset), function(ds_) {
  dd <- droplevels(d[d$dataset == ds_, ]); if (nlevels(dd$group_binary) < 2L) return(NULL)
  f <- lm(if (nlevels(dd$inferred_sex) > 1L) score ~ inferred_sex + group_binary else score ~ group_binary, dd)
  data.table(dataset = ds_, n = nrow(dd), n_control = sum(dd$group_binary == "Control"), n_disease = sum(dd$group_binary == "Disease"),
             estimate_sd_units = coef(f)["group_binaryDisease"], p_nominal = summary(f)$coefficients["group_binaryDisease", 4])
}))
fwrite(pc, file.path(out, "association_per_cohort_descriptive.tsv"), sep = "\t")
r2_design <- summary(lm(score ~ dataset + inferred_sex + group_binary, d))$r.squared

# ---- refit: reproduce baseline, then add the score
fit_calls <- function(cov_extra = NULL) {
  inf <- info; if (!is.null(cov_extra)) inf$score <- cov_extra
  fo <- if (is.null(cov_extra)) ~ dataset + inferred_sex + group_binary else ~ dataset + inferred_sex + group_binary + score
  design <- model.matrix(fo, data = inf)
  v <- voomWithQualityWeights(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$score_t <- if (!is.null(cov_extra)) fit$t[tt$gene, "score"] else NA_real_
  tt$score_p <- if (!is.null(cov_extra)) fit$p.value[tt$gene, "score"] else NA_real_
  tt
}
base <- as.data.table(fit_calls())
base[, call := adj.P.Val < 0.05 & abs(logFC) > 0.5]
ref <- deg0[match(base$gene, gene)]
repro <- data.table(n_calls_reference = sum(ref$padj < 0.05 & abs(ref$logFC) > 0.5), n_calls_here = sum(base$call),
                    identical_calls = identical(ref$padj < 0.05 & abs(ref$logFC) > 0.5, base$call),
                    max_abs_dlogFC = max(abs(ref$logFC - base$logFC)), max_abs_dpadj = max(abs(ref$padj - base$adj.P.Val)))
fwrite(repro, file.path(out, "baseline_reproduction.tsv"), sep = "\t")
print(repro)
if (!repro$identical_calls || repro$n_calls_reference != 1854L) stop("Baseline does not reproduce the 1,854 refit calls; stopping per PRESPEC")

summarise_vs_base <- function(tt) {
  tt <- as.data.table(tt); tt[, call := adj.P.Val < 0.05 & abs(logFC) > 0.5]
  data.table(n_calls = sum(tt$call), gained = sum(tt$call & !base$call), lost = sum(!tt$call & base$call),
             sign_flips_among_baseline_calls = sum(base$call & sign(tt$logFC) != sign(base$logFC)),
             logFC_pearson = cor(tt$logFC, base$logFC), max_abs_dlogFC = max(abs(tt$logFC - base$logFC)),
             n_genes_score_padj05 = if (all(is.na(tt$score_p))) NA_integer_ else sum(p.adjust(tt$score_p, "BH") < 0.05))
}
real <- fit_calls(score); sreal <- summarise_vs_base(real); sreal[, covariate := "score_z"]
rt <- as.data.table(real); rt[, call := adj.P.Val < 0.05 & abs(logFC) > 0.5]
chg <- data.table(gene = base$gene, symbol = deg0$symbol[match(base$gene, deg0$gene)], logFC_base = base$logFC, padj_base = base$adj.P.Val,
                  logFC_with_score = rt$logFC, padj_with_score = rt$adj.P.Val, call_base = base$call, call_with_score = rt$call)[call_base != call_with_score]
fwrite(chg, file.path(out, "calls_changed_by_score_covariate.tsv"), sep = "\t")

strata <- interaction(info$dataset, info$group_binary, drop = TRUE)
perm_cov <- function(i) { set.seed(SEED + i); p <- score; for (s in levels(strata)) { ix <- which(strata == s); p[ix] <- score[ix][sample.int(length(ix))] }; p }
perm_res <- mclapply(seq_len(NPERM), function(i) { s <- summarise_vs_base(fit_calls(perm_cov(i))); s[, covariate := paste0("permuted_within_dataset_group_", i)]; s },
                     mc.cores = NCORES, mc.set.seed = FALSE)
bad <- vapply(perm_res, function(z) inherits(z, "try-error"), TRUE)
if (any(bad)) stop("permuted-covariate worker failure")
perm_res <- rbindlist(perm_res)
fwrite(rbind(sreal, perm_res, fill = TRUE), file.path(out, "refit_with_score_covariate_summary.tsv"), sep = "\t")
collin <- data.table(r2_score_on_design = r2_design, n_units = nrow(info))
fwrite(collin, file.path(out, "score_collinearity.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(assoc); print(pc); print(sreal); print(perm_res[, .(gained, lost, logFC_pearson)]); print(collin)
