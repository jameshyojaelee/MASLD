# Feasibility gate L10-03 (agent gate-L10-03). Blinded power check only.
# GSE213621 bin labels are PERMUTED in every fit below; the observed
# GSE213621 stage contrast is never computed.
t0 <- Sys.time()
suppressMessages({library(data.table); library(edgeR); library(limma)})
set.seed(20261015)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rel <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
out <- "out"
dge_all <- readRDS(file.path(rel, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
stopifnot(identical(dim(dge_all), c(24115L, 774L)))
ann <- fread(file.path(root, "figures/candidates/fig3e-s3-774-20261011T151616Z/analysis/stage_extensions/stage_all_gene_results.tsv"),
             select = c("gene_id_versioned", "gene_name"))
ann <- unique(ann)
gname <- ann$gene_name[match(rownames(dge_all), ann$gene_id_versioned)]
pc_genes <- c("COL1A1", "LUM", "EPCAM", "CYP2C19", "GNMT")

bin_of <- function(fs) fifelse(fs <= 1L, "F0F1", fifelse(fs == 2L, "F2", "F3F4"))

## A. Discovery: disease-only grouped refit in the three staged cohorts
disc <- meta[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & group_binary == "Disease" &
             fibrosis_stage %in% 0:4 & !is.na(inferred_sex)]
disc[, bin := factor(bin_of(fibrosis_stage), levels = c("F0F1", "F2", "F3F4"))]
setorder(disc, sample_id)
dd <- dge_all[, match(disc$sample_id, colnames(dge_all))]
des_d <- model.matrix(~ factor(dataset) + factor(inferred_sex) + bin, disc)
v_d <- voomWithQualityWeights(dd, des_d, plot = FALSE)
f_d <- eBayes(lmFit(v_d, des_d))
disc_res <- rbindlist(lapply(c("binF2", "binF3F4"), function(cf) {
  tt <- topTable(f_d, coef = cf, number = Inf, sort.by = "none")
  data.table(gene = rownames(dge_all), gene_name = gname, contrast = cf, logFC = tt$logFC,
             SE = f_d$stdev.unscaled[, cf] * sqrt(f_d$s2.post), P = tt$P.Value,
             FDR = p.adjust(tt$P.Value, "BH"))
}))
fwrite(disc_res[, .(n_units = nrow(disc), n_genes = .N, n_bh05 = sum(FDR < 0.05),
                    n_bh05_up = sum(FDR < 0.05 & logFC > 0), n_bh05_down = sum(FDR < 0.05 & logFC < 0),
                    median_abs_logFC_hits = median(abs(logFC[FDR < 0.05])),
                    sd_logFC_hits = sd(logFC[FDR < 0.05])), by = contrast],
       file.path(out, "discovery_grouped_summary.tsv"), sep = "\t")
fwrite(disc[, .N, by = .(dataset, bin, inferred_sex)][order(dataset, bin, inferred_sex)],
       file.path(out, "discovery_grouped_n.tsv"), sep = "\t")
cat("Discovery fit done", format(Sys.time() - t0), "\n")

## B. GSE213621 disease units, permuted labels only
g <- meta[dataset == "GSE213621" & group_binary == "Disease"]
g[, bin := factor(bin_of(fibrosis_stage), levels = c("F0F1", "F2", "F3F4"))]
setorder(g, sample_id)
dg <- dge_all[, match(g$sample_id, colnames(dge_all))]
lcpm <- cpm(dg, log = TRUE, prior.count = 2)
cpm_raw <- cpm(dg)
expr <- data.table(gene = rownames(dg), gene_name = gname, mean_logCPM = rowMeans(lcpm),
                   frac_units_count_ge10 = rowMeans(dg$counts >= 10), mean_cpm = rowMeans(cpm_raw))
fwrite(expr[, .(n_genes = .N, n_mean_cpm_lt1 = sum(mean_cpm < 1), n_mean_cpm_lt0.5 = sum(mean_cpm < 0.5),
                n_zero_all = sum(rowSums(dg$counts) == 0))], file.path(out, "gse213621_expression_summary.tsv"), sep = "\t")

hits <- lapply(c(binF2 = "binF2", binF3F4 = "binF3F4"), function(cf) disc_res[contrast == cf & FDR < 0.05])
perm_fit <- function(labels, quality) {
  info <- data.frame(sex = factor(g$inferred_sex), bin = labels)
  des <- model.matrix(~ sex + bin, info)
  v <- if (quality) voomWithQualityWeights(dg, des, plot = FALSE) else voom(dg, des, plot = FALSE)
  eBayes(lmFit(v, des))
}
slope_stats <- function(fit, k) {
  rbindlist(lapply(names(hits), function(cf) {
    h <- hits[[cf]]
    y <- fit$coefficients[h$gene, cf]
    x <- h$logFC
    data.table(perm = k, contrast = cf, n_hits = length(x), slope = unname(coef(lm(y ~ x))[2]),
               sign_concordance = mean(sign(y) == sign(x)))
  }))
}
# First permutation: quality weights, as in the candidate; gives realistic per-gene SE.
lab1 <- sample(g$bin)
f1 <- perm_fit(lab1, TRUE)
se1 <- data.table(gene = rownames(dg), gene_name = gname,
                  SE_F3F4 = f1$stdev.unscaled[, "binF3F4"] * sqrt(f1$s2.post),
                  SE_F2 = f1$stdev.unscaled[, "binF2"] * sqrt(f1$s2.post),
                  sigma = f1$sigma, df_total = f1$df.total)
fwrite(se1, file.path(out, "gse213621_blinded_per_gene_se.tsv.gz"), sep = "\t")
null <- list(cbind(slope_stats(f1, 1L), quality_weights = TRUE))
cat("Quality-weight permuted fit done", format(Sys.time() - t0), "\n")
k <- 1L
while (as.numeric(difftime(Sys.time(), t0, units = "secs")) < 330 && k < 400L) {
  k <- k + 1L
  null[[k]] <- cbind(slope_stats(perm_fit(sample(g$bin), FALSE), k), quality_weights = FALSE)
}
null <- rbindlist(null)
fwrite(null, file.path(out, "null_slopes.tsv"), sep = "\t")
cat("Permutations:", k, format(Sys.time() - t0), "\n")

## C. Positive controls: expected grouped effect from discovery, blinded GSE213621 SE
pc <- merge(disc_res[gene_name %in% pc_genes, .(gene, gene_name, contrast, disc_logFC = logFC, disc_FDR = FDR)],
            se1[, .(gene, SE_F3F4, SE_F2)], by = "gene")
pc <- merge(pc, expr[, .(gene, mean_logCPM, frac_units_count_ge10)], by = "gene")
pc[, SE_g := fifelse(contrast == "binF3F4", SE_F3F4, SE_F2)]
z80 <- qnorm(0.8)
pc[, `:=`(MDE_nominal05 = (qnorm(0.975) + z80) * SE_g,
          MDE_alpha005 = (qnorm(1 - 0.005 / 2) + z80) * SE_g,
          MDE_bonf = (qnorm(1 - 0.05 / 24115 / 2) + z80) * SE_g)]
fwrite(pc[, .(gene_name, contrast, disc_logFC, disc_FDR, mean_logCPM, frac_units_count_ge10,
              SE_g, MDE_nominal05, MDE_alpha005, MDE_bonf)][order(gene_name, contrast)],
       file.path(out, "positive_controls.tsv"), sep = "\t")

## D. Per-gene MDE summary among discovery hits
mde <- rbindlist(lapply(names(hits), function(cf) {
  h <- hits[[cf]]
  se <- if (cf == "binF3F4") se1$SE_F3F4 else se1$SE_F2
  names(se) <- se1$gene
  s <- se[h$gene]
  data.table(contrast = cf, n_hits = nrow(h),
             median_SE = median(s), median_MDE_nominal05 = median((qnorm(0.975) + z80) * s),
             median_MDE_alpha005 = median((qnorm(1 - 0.005 / 2) + z80) * s),
             median_MDE_bonf = median((qnorm(1 - 0.05 / 24115 / 2) + z80) * s),
             frac_hits_disc_abs_logFC_gt_MDE_alpha005 = mean(abs(h$logFC) > (qnorm(1 - 0.005 / 2) + z80) * s),
             frac_hits_disc_abs_logFC_gt_MDE_bonf = mean(abs(h$logFC) > (qnorm(1 - 0.05 / 24115 / 2) + z80) * s),
             null_slope_mean = mean(null[contrast == cf, slope]), null_slope_sd = sd(null[contrast == cf, slope]),
             null_conc_mean = mean(null[contrast == cf, sign_concordance]),
             null_conc_sd = sd(null[contrast == cf, sign_concordance]),
             n_perm = null[contrast == cf, .N],
             MDE_slope_onesided_0.025 = (qnorm(0.975) + z80) * sd(null[contrast == cf, slope]),
             MDE_conc_onesided_0.025 = (qnorm(0.975) + z80) * sd(null[contrast == cf, sign_concordance]))
}))
fwrite(mde, file.path(out, "mde_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("DONE", format(Sys.time() - t0), "\n")
