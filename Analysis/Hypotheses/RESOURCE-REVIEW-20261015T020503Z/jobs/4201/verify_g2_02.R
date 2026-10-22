#!/usr/bin/env Rscript
# verify2-G2-02, job 4201. Independent check of the F1-vs-F0 complexity-covariate claim.
# Assay: bulk liver RNA-seq fragment counts. Model copied from scripts/figures/build_pi_stage_extensions.R
# fit_contrast(): voomWithQualityWeights, design ~ dataset + inferred_sex [+ covariates] + group, eBayes,
# BH over all 24,115 genes, no fold-change filter. Unit: genetic-identity cluster (not verified donor).
# Complexity covariates are recomputed here from the count matrix (not read from the reviewer table).
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table); library(parallel)})
set.seed(4201)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
STG <- file.path(ROOT, "figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions")
COMP <- file.path(ROOT, "figures/candidates/composition-774-20261011T151616Z/composition/accepted_composition_774_plus_unchanged4.tsv.gz")
R03 <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R03/complexity_774.tsv")
OUT <- getwd(); FITS <- file.path(OUT, "fits")
t0 <- Sys.time(); el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

d <- readRDS(file.path(REL, "merged_dge.rds"))
stopifnot(nrow(d) == 24115L, ncol(d) == 774L)
cnt <- d$counts; storage.mode(cnt) <- "double"
lib <- colSums(cnt)
cat("lib.size max abs diff colSums vs d$samples:", max(abs(lib - d$samples$lib.size)), "\n")
nf <- d$samples$norm.factors
nf2 <- calcNormFactors(DGEList(cnt), method = "TMM")$samples$norm.factors
cat("TMM recomputed on 774 vs stored norm.factors: max abs diff", signif(max(abs(nf - nf2)), 4), "\n")
# Gini via the Lorenz curve (different formula from the reviewer's rank-weighted sum)
gini_lorenz <- function(x) { x <- sort(x); n <- length(x); L <- cumsum(x) / sum(x); 1 - (2 * sum(L) - 1) / n }
gini <- apply(cnt, 2, gini_lorenz)
cx <- data.table(sample_id = colnames(d), lib_size = lib, gini = gini, norm_factor = nf)
if (file.exists(R03)) {
  r <- fread(R03)[match(cx$sample_id, sample_id)]
  cat("vs reviewer R03 table: max abs diff lib", max(abs(r$lib_size - cx$lib_size)),
      "gini", signif(max(abs(r$gini - cx$gini)), 4), "nf", signif(max(abs(r$norm_factor - cx$norm_factor)), 4), "\n")
}

s <- fread(file.path(STG, "stage_sample_manifest.tsv"))
stopifnot(nrow(s) == 387L, !anyDuplicated(s$sample_id), all(s$sample_id %in% colnames(d)))
s <- merge(s, cx, by = "sample_id")
zc <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
s[, `:=`(z_loglib = zc(log2(lib_size), dataset), z_gini = zc(gini, dataset), z_nf = zc(norm_factor, dataset),
         r_loglib = log2(lib_size), r_gini = gini, r_nf = norm_factor)]
s01 <- s[fib_group %in% c("F0", "F1")]
s01[, `:=`(y_loglib = zc(log2(lib_size), dataset), y_gini = zc(gini, dataset), y_nf = zc(norm_factor, dataset))]
cat("F0/F1 units:", table(s01$fib_group), "\n")
fwrite(s, file.path(OUT, "roster387_complexity.tsv"), sep = "\t")

# F1 vs F0 difference in each complexity metric, within cohort
bal <- rbindlist(lapply(c("lib_size", "gini", "norm_factor"), function(v) rbindlist(lapply(sort(unique(s01$dataset)), function(ds) {
  a <- s01[dataset == ds & fib_group == "F0"][[v]]; b <- s01[dataset == ds & fib_group == "F1"][[v]]
  data.table(metric = v, dataset = ds, n_F0 = length(a), n_F1 = length(b), median_F0 = median(a), median_F1 = median(b),
             wilcox_p = suppressWarnings(wilcox.test(b, a)$p.value))
}))))
fwrite(bal, file.path(OUT, "f0_f1_complexity_balance.tsv"), sep = "\t")
print(bal, digits = 4)

covsets <- list(base = character(0), lib = "z_loglib", gini = "z_gini", nf = "z_nf",
                joint = c("z_loglib", "z_gini", "z_nf"), joint_raw = c("r_loglib", "r_gini", "r_nf"),
                joint_z205 = c("y_loglib", "y_gini", "y_nf"))
fit_one <- function(lab) {
  ts <- el(); sub <- copy(s01); setorder(sub, sample_id)
  info <- data.frame(dataset = droplevels(factor(sub$dataset)), inferred_sex = droplevels(factor(sub$inferred_sex)),
                     group = factor(sub$fib_group, levels = c("F0", "F1")), row.names = sub$sample_id)
  for (v in covsets[[lab]]) info[[v]] <- sub[[v]]
  X <- model.matrix(as.formula(paste("~", paste(c("dataset", "inferred_sex", covsets[[lab]], "group"), collapse = " + "))), info)
  stopifnot(qr(X)$rank == ncol(X))
  fit <- eBayes(lmFit(voomWithQualityWeights(d[, match(sub$sample_id, colnames(d))], X, plot = FALSE), X))
  tt <- topTable(fit, coef = "groupF1", number = Inf, sort.by = "none")
  g <- data.table(gene_id_versioned = rownames(tt), logFC = tt$logFC, P.Value = tt$P.Value, FDR = p.adjust(tt$P.Value, "BH"))
  fwrite(g, file.path(FITS, paste0(lab, ".genes.tsv.gz")), sep = "\t")
  cat(lab, "n", nrow(sub), "BH<0.05", sum(g$FDR < 0.05), "fit_s", el() - ts, "end", el(), "\n")
  lab
}
done <- mclapply(c("base", "joint", "lib", "gini", "nf", "joint_raw", "joint_z205"),
                 function(l) tryCatch(fit_one(l), error = function(e) paste("ERROR", l, conditionMessage(e))),
                 mc.cores = 4, mc.preschedule = FALSE)
print(unlist(done))

rel <- fread(file.path(STG, "stage_extension_all_gene_results.tsv"), select = c("gene_id_versioned", "logFC", "P.Value", "FDR", "contrast", "axis"))
rel <- rel[axis == "fibrosis" & contrast == "F1_vs_F0"]
stopifnot(nrow(rel) == 24115L)
R <- rel[FDR < 0.05, gene_id_versioned]
cat("released F1_vs_F0 BH<0.05:", length(R), " nominal P<0.01:", sum(rel$P.Value < 0.01), "\n")
summ <- rbindlist(lapply(names(covsets), function(lab) {
  f <- file.path(FITS, paste0(lab, ".genes.tsv.gz")); if (!file.exists(f)) return(NULL)
  g <- fread(f)[match(rel$gene_id_versioned, gene_id_versioned)]
  sig <- g$gene_id_versioned[g$FDR < 0.05]; i <- match(R, g$gene_id_versioned)
  ratio <- g$logFC[i] / rel$logFC[match(R, rel$gene_id_versioned)]
  data.table(fit = lab, bh_sig = length(sig), up = sum(g$FDR < 0.05 & g$logFC > 0), down = sum(g$FDR < 0.05 & g$logFC < 0),
             retained_of_released = length(intersect(sig, R)), retained_frac = round(length(intersect(sig, R)) / length(R), 4),
             bh_sig_frac_of_1338 = round(length(sig) / length(R), 4),
             nominal_p_lt_0.01 = sum(g$P.Value < 0.01),
             median_logFC_ratio_on_released = round(median(ratio), 4),
             max_abs_logFC_diff_vs_released = signif(max(abs(g$logFC - rel$logFC)), 4))
}))
fwrite(summ, file.path(OUT, "f1_vs_f0_refit_summary.tsv"), sep = "\t")
print(summ)

# Within-stage association of expression with each complexity covariate (FWL-equivalent coefficient)
lc <- cpm(d[, s$sample_id], log = TRUE, prior.count = 2)
align <- rbindlist(lapply(c("z_loglib", "z_gini", "z_nf"), function(v) {
  X <- model.matrix(as.formula(paste("~ dataset + inferred_sex + fib_group +", v)), s)
  b <- lmFit(lc, X)$coefficients[, v]
  i <- match(R, names(b)); e <- rel$logFC[match(R, rel$gene_id_versioned)]
  data.table(covariate = v, n_units = nrow(s), n_released_calls = length(R),
             spearman_released_logFC_vs_within_stage_slope = round(cor(e, b[i], method = "spearman"), 4),
             spearman_all_genes = round(cor(rel$logFC, b[match(rel$gene_id_versioned, names(b))], method = "spearman"), 4))
}))
fwrite(align, file.path(OUT, "within_stage_alignment.tsv"), sep = "\t")
print(align)

# Complexity versus fibrosis within cohort (774 meta; same 8-metric family as reviewer R03 for BH)
m <- as.data.table(readRDS(file.path(REL, "meta_matched.rds")))
w <- fread(file.path(REL, "lvqw_sample_weights.tsv"))
p <- sweep(cnt, 2, lib, "/")
cx2 <- data.table(sample_id = colnames(d), lib_size = lib, detected = colSums(cnt > 0),
  entropy = apply(p, 2, function(x) { x <- x[x > 0]; -sum(x * log(x)) }), gini = gini,
  top10 = apply(p, 2, function(x) sum(sort(x, decreasing = TRUE)[1:10])),
  top100 = apply(p, 2, function(x) sum(sort(x, decreasing = TRUE)[1:100])), norm_factor = nf)
cx2 <- merge(merge(cx2, m[, .(sample_id, dataset, fibrosis_stage)], by = "sample_id"), w, by = "sample_id")
cf <- rbindlist(lapply(c("GSE130970", "GSE135251", "GSE162694"), function(ds) rbindlist(lapply(
  c("lib_size", "detected", "entropy", "gini", "top10", "top100", "norm_factor", "sample_weight"), function(v) {
    q <- cx2[dataset == ds & !is.na(fibrosis_stage)]
    ct <- suppressWarnings(cor.test(as.numeric(q$fibrosis_stage), q[[v]], method = "spearman"))
    data.table(dataset = ds, metric = v, n = nrow(q), rho = unname(ct$estimate), p = ct$p.value)
  }))))
cf[, bh_q_24 := p.adjust(p, "BH")]
fwrite(cf, file.path(OUT, "complexity_vs_fibrosis_recomputed.tsv"), sep = "\t")
print(cf[metric %in% c("lib_size", "gini", "norm_factor")], digits = 4)

# Mean deconvolved hepatocyte fraction by fibrosis group on the 387 stage units
comp <- fread(COMP, select = c("sample_id", "Hepatocytes"))
hs <- merge(s[, .(sample_id, fib_group)], comp, by = "sample_id")
print(hs[, .(n = .N, mean_hepatocyte = round(mean(Hepatocytes), 4)), by = fib_group][order(fib_group)])
cat("done", el(), "s\n")
