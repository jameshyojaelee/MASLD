#!/usr/bin/env Rscript
# figS_disease_signature_sweep_deg_tables.R
# Companion to figS_disease_signature_sweep.R. Fits each of the 6 DE methods ONCE
# on the FULL 5-cohort MEGA data (Disease vs Control, cohort + sex adjusted) and
# writes the COMPLETE per-gene result table for each — ALL genes, NO lfc/padj
# cutoff. One CSV per method, with the method name in the filename.
#
# Output → figures/supplementary/figS_methods_validation/disease_signature_sweep/deg_tables/
#   deg_limma_voom_full.csv,  deg_dream_full.csv,        deg_deseq2_full.csv
#   deg_edger_qlf_full.csv,   deg_metafor_re_full.csv,   deg_combatseq_deseq2_full.csv
#   deg_all_methods_long.csv  (harmonised long-format concatenation of all six)
#
# Harmonised schema per CSV: method, gene, logFC, stat, pvalue, padj  (+ native extras).

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(BiocParallel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS_METHVAL_DIR, "disease_signature_sweep", "deg_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

N_CPUS      <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
BPPARAM_PAR <- MulticoreParam(N_CPUS)
set.seed(42)

MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

# ── Data loading (identical to the sweep script) ────────────────────────────────
cat("Loading merged DGE ...\n")
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep_s <- samp$dataset %in% MEGA
dge  <- dge[, keep_s]; samp <- samp[keep_s]

# Complete inferred_sex (zero NA) from meta_matched.rds — see sweep script for why.
meta_m <- as.data.frame(readRDS(file.path(INT_RESULTS, "meta_matched.rds")))
samp[, inferred_sex := meta_m$inferred_sex[match(sample_id, meta_m$sample_id)]]
if (anyNA(samp$inferred_sex)) stop("inferred_sex join left NA — sample_id mismatch")

meta <- data.frame(
  group_binary = factor(samp$group_binary, levels = c("Control", "Disease")),
  inferred_sex = factor(samp$inferred_sex),
  dataset      = factor(samp$dataset),
  row.names    = samp$sample_id
)
counts <- dge$counts
cat(sprintf("Full MEGA fit: %d samples, %d genes, %d cohorts (Control=%d, Disease=%d)\n",
            ncol(counts), nrow(counts), nlevels(meta$dataset),
            sum(meta$group_binary == "Control"), sum(meta$group_binary == "Disease")))

# ── Per-method FULL DEG tables ──────────────────────────────────────────────────
# Each returns a data.table with the harmonised core columns
#   gene, logFC, stat, pvalue, padj
# followed by any native extra columns. NO filtering applied.

deg_limma_voom <- function() {
  des <- model.matrix(~ dataset + inferred_sex + group_binary, data = meta)
  v   <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(counts)), des)
  fit <- limma::eBayes(limma::lmFit(v, des))
  tt  <- limma::topTable(fit, coef = "group_binaryDisease", n = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$t,
             pvalue = tt$P.Value, padj = tt$adj.P.Val, AveExpr = tt$AveExpr)
}

deg_dream <- function() {
  suppressPackageStartupMessages(library(variancePartition))
  form <- ~ group_binary + inferred_sex + (1 | dataset)
  v    <- variancePartition::voomWithDreamWeights(
            edgeR::calcNormFactors(edgeR::DGEList(counts)), form, meta, BPPARAM = BPPARAM_PAR)
  fit  <- variancePartition::dream(v, form, meta, BPPARAM = BPPARAM_PAR)
  fit  <- variancePartition::eBayes(fit)
  tt   <- limma::topTable(fit, coef = "group_binaryDisease", n = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, stat = tt$t,
             pvalue = tt$P.Value, padj = tt$adj.P.Val, AveExpr = tt$AveExpr)
}

deg_deseq2 <- function() {
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, meta, ~ dataset + inferred_sex + group_binary)
  suppressWarnings(dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM_PAR, quiet = TRUE))
  # independentFiltering=FALSE so padj is defined for ALL genes (no cutoff/filtering)
  res <- DESeq2::results(dds, contrast = c("group_binary", "Disease", "Control"),
                         independentFiltering = FALSE)
  data.table(gene = rownames(res), logFC = res$log2FoldChange, stat = res$stat,
             pvalue = res$pvalue, padj = res$padj,
             baseMean = res$baseMean, lfcSE = res$lfcSE)
}

deg_edger_qlf <- function() {
  des <- model.matrix(~ dataset + inferred_sex + group_binary, data = meta)
  y   <- edgeR::calcNormFactors(edgeR::DGEList(counts, group = meta$group_binary))
  y   <- edgeR::estimateDisp(y, des)
  fit <- edgeR::glmQLFit(y, des)
  qlf <- edgeR::glmQLFTest(fit, coef = "group_binaryDisease")
  tt  <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  # stat kept signed (sign(logFC)*sqrt(F)) for ranking; native F/PValue/FDR retained
  data.table(gene = rownames(tt), logFC = tt$logFC,
             stat = sign(tt$logFC) * sqrt(pmax(tt$F, 0)),
             pvalue = tt$PValue, padj = tt$FDR, logCPM = tt$logCPM, F = tt$F)
}

deg_combatseq_deseq2 <- function() {
  suppressPackageStartupMessages({ library(sva); library(DESeq2) })
  cat("    ComBat-seq (full data) ...\n")
  cc  <- sva::ComBat_seq(counts, batch = samp$dataset, group = samp$group_binary, full_mod = TRUE)
  cc  <- pmax(cc, 0L)
  dds <- DESeq2::DESeqDataSetFromMatrix(cc, meta, ~ inferred_sex + group_binary)
  suppressWarnings(dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM_PAR, quiet = TRUE))
  res <- DESeq2::results(dds, contrast = c("group_binary", "Disease", "Control"),
                         independentFiltering = FALSE)
  data.table(gene = rownames(res), logFC = res$log2FoldChange, stat = res$stat,
             pvalue = res$pvalue, padj = res$padj,
             baseMean = res$baseMean, lfcSE = res$lfcSE)
}

deg_metafor_re <- function() {
  # Vectorized DerSimonian-Laird RE meta-analysis across all 5 MEGA cohorts.
  per_cohort_dt <- rbindlist(lapply(MEGA, function(coh) {
    idx   <- which(samp$dataset == coh)
    grp_c <- droplevels(factor(samp$group_binary[idx], levels = c("Control", "Disease")))
    if (length(unique(grp_c)) < 2) return(NULL)
    sex_c <- droplevels(factor(samp$inferred_sex[idx]))
    terms <- if (nlevels(sex_c) >= 2) c("sex_c", "grp_c") else "grp_c"
    des_c <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))))
    v_c   <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(counts[, idx])), des_c)
    fit_c <- limma::eBayes(limma::lmFit(v_c, des_c))
    tt_c  <- limma::topTable(fit_c, coef = "grp_cDisease", n = Inf, sort.by = "none")
    data.table(gene = rownames(tt_c), cohort = coh, logFC = tt_c$logFC,
               se = abs(tt_c$logFC / pmax(abs(tt_c$t), 1e-10)))
  }))
  dt_lfc <- dcast(per_cohort_dt, gene ~ cohort, value.var = "logFC")
  dt_se  <- dcast(per_cohort_dt, gene ~ cohort, value.var = "se")
  coh_cols <- setdiff(names(dt_lfc), "gene")
  LFC <- as.matrix(dt_lfc[, ..coh_cols]); SE <- as.matrix(dt_se[, ..coh_cols])
  bad <- !is.finite(SE) | SE <= 0 | !is.finite(LFC); LFC[bad] <- NA; SE[bad] <- NA
  K   <- rowSums(!is.na(SE))
  W   <- 1 / SE^2; W[is.na(W)] <- 0; Ws <- rowSums(W)
  wL  <- W * LFC; wL[is.na(wL)] <- 0
  mu_fe <- rowSums(wL) / pmax(Ws, 1e-10)
  Q   <- rowSums(W * (LFC - mu_fe)^2, na.rm = TRUE)
  C_  <- Ws - rowSums(W^2) / pmax(Ws, 1e-10)
  tau2 <- pmax(0, (Q - pmax(K - 1, 0)) / pmax(C_, 1e-10))
  Wr  <- 1 / (SE^2 + tau2); Wr[is.na(Wr)] <- 0; Wrs <- rowSums(Wr)
  wLr <- Wr * LFC; wLr[is.na(wLr)] <- 0
  mu  <- rowSums(wLr) / pmax(Wrs, 1e-10)
  se  <- sqrt(1 / pmax(Wrs, 1e-10))
  z   <- mu / pmax(se, 1e-10); z[K < 2] <- NA
  pval <- 2 * pnorm(-abs(z))
  data.table(gene = dt_lfc$gene, logFC = mu, stat = z, pvalue = pval,
             padj = p.adjust(pval, "BH"), se_re = se, tau2 = tau2, k_cohorts = K)
}

METHODS <- list(
  limma_voom       = deg_limma_voom,
  dream            = deg_dream,
  deseq2           = deg_deseq2,
  edger_qlf        = deg_edger_qlf,
  metafor_re       = deg_metafor_re,
  combatseq_deseq2 = deg_combatseq_deseq2
)

# ── Run + write ─────────────────────────────────────────────────────────────────
long_all <- list()
for (m in names(METHODS)) {
  cat(sprintf("\n[%s] full DEG table ...\n", m))
  t0 <- proc.time()["elapsed"]
  tab <- METHODS[[m]]()
  tab <- cbind(method = m, tab)
  fn  <- file.path(OUT, sprintf("deg_%s_indep.csv", m))  # method + source IN the filename
  fwrite(tab, fn)
  n_sig <- sum(tab$padj < 0.05, na.rm = TRUE)
  cat(sprintf("[%s] %d genes, %d with padj<0.05 (for reference only; file is UNFILTERED), %.1f s -> %s\n",
              m, nrow(tab), n_sig, proc.time()["elapsed"] - t0, basename(fn)))
  long_all[[m]] <- tab[, .(method, gene, logFC, stat, pvalue, padj)]
}

fwrite(rbindlist(long_all), file.path(OUT, "deg_long_indep.csv"))
cat(sprintf("\nWrote 6 per-method tables + deg_long_indep.csv to:\n  %s\n", OUT))
