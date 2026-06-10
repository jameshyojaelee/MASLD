#!/usr/bin/env Rscript
# Phase 2 — transcript/isoform diversity metrics per gene.
# Metrics (per sample, then summarized): expressed-isoform count N_g,
# dominant-isoform fraction DIF = max_t p_t, Shannon entropy H = -sum p_t log2 p_t,
# effective # isoforms N_eff = 2^H, with p_t = TPM_t / sum_t TPM_t.
# An isoform "counts" only if it is >=5% of the gene AND TPM>=1 (same definition as
# the descriptive panel F); proportions are renormalized over the kept isoforms.
# Shift test: limma on per-sample N_eff (and DIF) ~ cohort + group, controlling cohort.
#
# Usage: Rscript 03_diversity.R <species: human|mouse>
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
species <- if (length(args) >= 1) args[1] else "human"
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RES  <- file.path(PROJ, "RNA-seq/results/isoform_diversity", species)

tpm   <- readRDS(file.path(RES, "tx_tpm.rds"))               # tx x sample
samp  <- fread(file.path(RES, "samples.tsv"))
t2g_f <- if (species == "human") "tx2gene_v49_primary.tsv.gz" else "tx2gene_vM38_primary.tsv.gz"
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index", t2g_f))
COHORT_COL <- "dataset"

# --- sample QC + group filter (+ optional chemistry stratification) ---
setkey(samp, sample_id); samp <- samp[colnames(tpm)]
keep_s <- !is.na(samp$group_binary) & samp$group_binary %in% c("Control","Disease") &
          (is.na(samp$pass_technical) | samp$pass_technical %in% c(TRUE,"TRUE","True"))
LIBTYPE <- Sys.getenv("LIBTYPE", "")          # "", "polyA", or "total-RNA"
sfx <- ""
if (nzchar(LIBTYPE)) { keep_s <- keep_s & !is.na(samp$library_type) & samp$library_type == LIBTYPE
                       sfx <- paste0("_", LIBTYPE) }
cat(sprintf("[diversity] %d/%d samples kept (group+QC%s)\n", sum(keep_s), nrow(samp),
            if (nzchar(LIBTYPE)) paste0("; library_type=", LIBTYPE) else ""))
tpm <- tpm[, keep_s]; samp <- samp[keep_s]
grp <- factor(samp$group_binary, levels = c("Control","Disease"))
coh <- factor(samp[[COHORT_COL]])

# --- expressed transcripts & gene map ---
tx2gene <- tx2gene[match(rownames(tpm), tx2gene$txname)]
stopifnot(nrow(tx2gene) == nrow(tpm))
expr_tx <- rowSums(tpm >= 1) >= 10                # expressed in >=10 samples
cat(sprintf("[diversity] expressed transcripts: %d / %d\n", sum(expr_tx), length(expr_tx)))
M <- tpm[expr_tx, , drop = FALSE]
gene <- tx2gene$geneid[expr_tx]
sym  <- tx2gene$gene_symbol[expr_tx]
gbt  <- tx2gene$gene_biotype[expr_tx]
txbt <- tx2gene$tx_biotype[expr_tx]
ri_gene <- unique(gene[txbt == "retained_intron"])   # chemistry-sensitive isoform class

n_tx_per_gene <- table(gene)
multi <- names(n_tx_per_gene)[n_tx_per_gene >= 2]
cat(sprintf("[diversity] genes expressed: %d ; multi-isoform: %d\n",
            length(n_tx_per_gene), length(multi)))

# --- per-sample metrics for multi-isoform genes (gene-by-gene; memory-safe) ---
idx_by_gene <- split(seq_len(nrow(M)), gene)
idx_by_gene <- idx_by_gene[multi]
ns <- ncol(M)
Neff <- matrix(NA_real_, length(multi), ns, dimnames = list(multi, colnames(M)))
DIF  <- matrix(NA_real_, length(multi), ns, dimnames = list(multi, colnames(M)))
for (i in seq_along(idx_by_gene)) {
  sub  <- M[idx_by_gene[[i]], , drop = FALSE]
  gsum <- colSums(sub)
  p    <- sweep(sub, 2, gsum, "/")
  keep <- (p >= 0.05) & (sub >= 1)        # an isoform counts only if >=5% of gene AND TPM>=1
  sub  <- sub * keep                       # drop trace / EM-spillover isoforms (matches panel F)
  gsum <- colSums(sub)
  ok   <- gsum >= 1                        # gene analyzable in that sample (>=1 kept isoform)
  p    <- sweep(sub, 2, gsum, "/")         # renormalize proportions over the kept isoforms
  pl   <- p * log2(p); pl[!is.finite(pl)] <- 0
  H    <- -colSums(pl)
  Neff[i, ok] <- (2 ^ H)[ok]
  DIF[i, ok]  <- colMaxs(p)[ok]
}
saveRDS(Neff, file.path(RES, sprintf("neff_per_sample%s.rds", sfx)))
saveRDS(DIF,  file.path(RES, sprintf("dif_per_sample%s.rds", sfx)))

# --- group-mean summaries (per gene) for distributions; also single-iso genes ---
mean_grp <- function(mat, g) {
  data.table(gene = rownames(mat),
             ctrl = rowMeans(mat[, grp=="Control", drop=FALSE], na.rm=TRUE),
             dis  = rowMeans(mat[, grp=="Disease", drop=FALSE], na.rm=TRUE))
}
neff_s <- mean_grp(Neff); setnames(neff_s, c("gene","Neff_ctrl","Neff_dis"))
dif_s  <- mean_grp(DIF);  setnames(dif_s,  c("gene","DIF_ctrl","DIF_dis"))

# expressed-isoform count per group (group-mean TPM >=1 & >=5% of gene; matches panel F)
ng_count <- function(grp_lab) {
  cols    <- grp == grp_lab
  meanTPM <- rowMeans(M[, cols, drop=FALSE])
  gsum    <- tapply(meanTPM, gene, sum)
  ex      <- meanTPM >= 1 & (meanTPM / gsum[gene]) >= 0.05
  tapply(ex, gene, sum)
}
ng_c <- ng_count("Control"); ng_d <- ng_count("Disease")
ng_dt <- data.table(gene = names(ng_c), N_g_ctrl = as.integer(ng_c),
                    N_g_dis = as.integer(ng_d[names(ng_c)]))

# --- diversity-shift test (limma, cohort-adjusted) on multi-isoform genes ---
# drop the cohort term for single-cohort strata (e.g. chemistry-restricted runs)
design <- if (nlevels(droplevels(coh)) > 1) model.matrix(~ coh + grp) else model.matrix(~ grp)
# Transform responses to satisfy limma's Gaussian/homoscedastic assumptions (review P1):
# N_eff (=2^entropy) is right-skewed & heteroscedastic -> log2 (the natural entropy scale);
# DIF is a [0,1] proportion -> logit. WARNING: when this run pools chemistries (LIBTYPE
# unset), cohort==chemistry is confounded and the shift test is chemistry-driven (review
# P0) -> interpret the per-LIBTYPE (polyA / total-RNA) runs as primary, pooled as context.
shift_test <- function(mat, transform) {
  ok_genes <- rowSums(!is.na(mat)) >= ncol(mat) * 0.5      # testable in >=50% samples
  m <- mat[ok_genes, , drop = FALSE]
  m <- if (transform == "log2") log2(m)
       else qlogis(pmin(pmax(m, 1e-3), 1 - 1e-3))          # logit, boundary-clamped
  fit <- eBayes(lmFit(m, design))
  cf  <- grep("grpDisease", colnames(design), value = TRUE)
  data.table(gene = rownames(m),
             delta = fit$coefficients[, cf],               # on transformed scale
             t = fit$t[, cf], p = fit$p.value[, cf],
             padj = p.adjust(fit$p.value[, cf], "BH"))
}
neff_test <- shift_test(Neff, "log2");  setnames(neff_test, c("gene","dNeff","Neff_t","Neff_p","Neff_padj"))  # dNeff now log2-scale
dif_test  <- shift_test(DIF,  "logit"); setnames(dif_test,  c("gene","dDIF","DIF_t","DIF_p","DIF_padj"))      # dDIF now logit-scale

# --- assemble per-gene table ---
ann <- unique(data.table(gene = gene, symbol = sym, gene_biotype = gbt))
out <- Reduce(function(a,b) merge(a,b,by="gene",all=TRUE),
              list(ann, ng_dt, neff_s, dif_s, neff_test, dif_test))
out[, n_iso_expressed := as.integer(n_tx_per_gene[gene])]
out[, multi_isoform := gene %in% multi]
out[, has_retained_intron := gene %in% ri_gene]   # flag chemistry-sensitive genes
out[, library_type := if (nzchar(LIBTYPE)) LIBTYPE else "all"]
fwrite(out, file.path(RES, sprintf("isoform_diversity_%s%s.tsv", species, sfx)), sep = "\t")
cat(sprintf("[diversity] wrote isoform_diversity_%s%s.tsv (%d genes; %d Neff-shift padj<0.05; %d retained-intron genes)\n",
            species, sfx, nrow(out), sum(out$Neff_padj < 0.05, na.rm=TRUE), length(ri_gene)))
