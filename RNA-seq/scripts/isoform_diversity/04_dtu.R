#!/usr/bin/env Rscript
# Phase 3 — Differential Transcript Usage (primary: satuRn + stageR).
# Filter (DRIMSeq::dmFilter) -> satuRn fitDTU/testDTU on dtuScaledTPM counts ->
# stageR two-stage (gene screen -> transcript confirmation) at OFDR 0.05.
# Design: ~ cohort + sex + group_binary ; contrast Disease - Control.
# Point-estimate DTU only: the kallisto/0.51.1 build lacks HDF5 (no bootstraps), so a
# swish inferential-replicate arm is not available here (would need salmon/Gibbs requant).
#
# Usage: Rscript 04_dtu.R <species> [contrast: group|stage]
suppressPackageStartupMessages({
  library(data.table); library(DRIMSeq); library(satuRn)
  library(SummarizedExperiment); library(stageR); library(edgeR)
})
args <- commandArgs(trailingOnly = TRUE)
species  <- if (length(args) >= 1) args[1] else "human"
contrast <- if (length(args) >= 2) args[2] else "group"
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RES  <- file.path(PROJ, "RNA-seq/results/isoform_diversity", species)

cts  <- readRDS(file.path(RES, "tx_counts_dtuscaled.rds"))   # tx x sample
samp <- fread(file.path(RES, "samples.tsv")); setkey(samp, sample_id); samp <- samp[colnames(cts)]
t2g_f <- if (species == "human") "tx2gene_v49_primary.tsv.gz" else "tx2gene_vM38_primary.tsv.gz"
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index", t2g_f))
tx2gene <- tx2gene[match(rownames(cts), tx2gene$txname)]

# --- sample filter (+ optional chemistry stratification) ---
keep <- !is.na(samp$group_binary) & samp$group_binary %in% c("Control","Disease") &
        (is.na(samp$pass_technical) | samp$pass_technical %in% c(TRUE,"TRUE","True")) &
        !is.na(samp$sex)
LIBTYPE <- Sys.getenv("LIBTYPE", ""); sfx <- ""
if (nzchar(LIBTYPE)) { keep <- keep & !is.na(samp$library_type) & samp$library_type == LIBTYPE
                       sfx <- paste0(sfx, "_", LIBTYPE) }
# Single-diet-model contrast (mouse): disease = that diet; controls = matched controls
# from the same dataset(s) (handles GSE162876 controls shared by CDAHFD & FPC).
DIET <- Sys.getenv("DIET", "")
if (nzchar(DIET) && "diet_model" %in% names(samp)) {
  dis_ds <- unique(samp$dataset[samp$group_binary == "Disease" & samp$diet_model == DIET])
  keep <- keep & ((samp$group_binary == "Disease" & samp$diet_model == DIET) |
                  (samp$group_binary == "Control" & samp$dataset %in% dis_ds))
  sfx <- paste0(sfx, "_diet-", DIET)
  cat(sprintf("[dtu] DIET=%s disease-datasets: %s\n", DIET, paste(dis_ds, collapse=",")))
}
# Single-cohort DTU (no integration): one cohort = one chemistry, so no confound to
# adjust away; the cohort term drops out and design becomes ~0+group(+sex).
COHORT <- Sys.getenv("COHORT", "")
if (nzchar(COHORT)) { keep <- keep & samp$dataset == COHORT; sfx <- paste0(sfx, "_cohort-", COHORT) }
cts <- cts[, keep]; samp <- samp[keep]
cat(sprintf("[dtu] %d samples, %d Control / %d Disease%s\n", nrow(samp),
            sum(samp$group_binary=="Control"), sum(samp$group_binary=="Disease"),
            if (nzchar(LIBTYPE)) paste0("; library_type=", LIBTYPE) else ""))

# --- DRIMSeq filter (n.small = smallest group) ---
n  <- nrow(samp); n.small <- min(table(samp$group_binary))
counts_df <- data.frame(gene_id = tx2gene$geneid, feature_id = tx2gene$txname,
                        as.data.frame(cts), check.names = FALSE)
sdat <- data.frame(sample_id = samp$sample_id, group = samp$group_binary,
                   cohort = samp$dataset, sex = samp$sex,
                   stage = suppressWarnings(as.numeric(samp$fibrosis_stage)))
d <- dmDSdata(counts = counts_df, samples = sdat)
# gene filter at n.small (smallest group), NOT all n — DTU is a within-gene-proportion
# test; requiring expression in every sample of an 11-dataset pool drops ~45% of mouse
# genes (review P1). Feature-level proportion floor still guards low-usage isoforms.
d <- dmFilter(d, min_samps_feature_expr = n.small, min_feature_expr = 10,
              min_samps_feature_prop = n.small, min_feature_prop = 0.10,
              min_samps_gene_expr = n.small, min_gene_expr = 10)
cat(sprintf("[dtu] after dmFilter: %d transcripts in %d genes\n",
            nrow(counts(d)), length(unique(counts(d)$gene_id))))

# --- satuRn SummarizedExperiment ---
fcts <- as.matrix(counts(d)[, -c(1,2)]); rownames(fcts) <- counts(d)$feature_id
txInfo <- data.frame(isoform_id = counts(d)$feature_id, gene_id = counts(d)$gene_id,
                     row.names = counts(d)$feature_id)
sdat$group  <- factor(sdat$group, levels = c("Control","Disease"))
sdat$cohort <- factor(sdat$cohort)
sdat$sex    <- factor(sdat$sex)
sumExp <- SummarizedExperiment(assays = list(counts = fcts),
                               colData = DataFrame(sdat), rowData = txInfo)

# satuRn: formula references colData columns; build the SAME design to define L.
# Adjust only for covariates that vary in the (possibly chemistry-stratified) subset
# — single-cohort strata (e.g. human polyA = GSE135251 only) drop the cohort term.
adj <- c(if (nlevels(droplevels(sdat$cohort)) > 1) "cohort",
         if (nlevels(droplevels(sdat$sex))    > 1) "sex")
if (contrast == "group") {
  form   <- as.formula(paste(c("~ 0 + group", adj), collapse = " + "))
  design <- model.matrix(form, data = as.data.frame(colData(sumExp)))
  L <- matrix(0, ncol(design), 1, dimnames = list(colnames(design), "Dis_vs_Ctrl"))
  L["groupDisease", 1] <- 1; L["groupControl", 1] <- -1
} else {  # fibrosis stage as continuous trend
  keep2  <- !is.na(sdat$stage)
  sumExp <- sumExp[, keep2]; sdat <- sdat[keep2, ]
  adj <- c(if (nlevels(droplevels(factor(sdat$cohort))) > 1) "cohort",
           if (nlevels(droplevels(factor(sdat$sex)))    > 1) "sex")
  form   <- as.formula(paste(c("~ stage", adj), collapse = " + "))
  design <- model.matrix(form, data = as.data.frame(colData(sumExp)))
  L <- matrix(0, ncol(design), 1, dimnames = list(colnames(design), "stage_trend"))
  L["stage", 1] <- 1
}
metadata(sumExp)$formula <- form
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
sumExp <- satuRn::fitDTU(object = sumExp, formula = form, parallel = TRUE,
                         BPPARAM = BiocParallel::MulticoreParam(ncpu), verbose = TRUE)
sumExp <- satuRn::testDTU(object = sumExp, contrasts = L,
                          diagplot1 = FALSE, diagplot2 = FALSE, sort = FALSE)
res <- rowData(sumExp)[[paste0("fitDTUResult_", colnames(L))]]
res$isoform_id <- rownames(res); res$gene_id <- txInfo[rownames(res), "gene_id"]
setDT(res)

# --- per-transcript Delta-proportion (Disease - Control) on the TESTED dtuScaledTPM ---
fmat <- assay(sumExp, "counts"); gid_f <- txInfo[rownames(fmat), "gene_id"]
gsum <- rowsum(fmat, gid_f); prop <- fmat / gsum[match(gid_f, rownames(gsum)), ]
gv   <- as.character(colData(sumExp)$group)
dprop <- setNames(rowMeans(prop[, gv=="Disease", drop=FALSE], na.rm=TRUE) -
                  rowMeans(prop[, gv=="Control", drop=FALSE], na.rm=TRUE), rownames(fmat))
res[, dprop := dprop[isoform_id]]

# --- stageR two-stage on BOTH empirical and theoretical (raw) p-values ---
simes <- function(p){ p<-sort(p); if(!length(p)) return(NA_real_); min(length(p)*p/seq_along(p)) }  # sort() drops NA; guard all-NA gene
run_stage <- function(pv) {
  names(pv) <- res$isoform_id
  pscr <- tapply(pv, res$gene_id, simes); pscr <- setNames(as.numeric(pscr), names(pscr)); pscr <- pscr[!is.na(pscr)]
  so <- stageWiseAdjustment(stageRTx(pScreen=pscr,
          pConfirmation=matrix(pv, ncol=1, dimnames=list(res$isoform_id,"tx")),
          pScreenAdjusted=FALSE, tx2gene=data.frame(transcript=res$isoform_id, gene=res$gene_id)),
          method="dtu", alpha=0.05, allowNA=TRUE)
  p <- as.data.table(getAdjustedPValues(so, onlySignificantGenes=FALSE, order=FALSE))
  setnames(p, c("txID","gene","transcript"), c("isoform_id","gscreen","tconf")); p[, .(isoform_id, gscreen, tconf)]
}
emp <- run_stage(res$empirical_pval); setnames(emp, c("gscreen","tconf"), c("gene_screen_padj","tx_confirm_padj"))
theo<- run_stage(res$pval);           setnames(theo,c("gscreen","tconf"), c("gene_screen_padj_theo","tx_confirm_padj_theo"))
res <- merge(merge(res, emp, by="isoform_id", all.x=TRUE), theo, by="isoform_id", all.x=TRUE)
res <- merge(res, tx2gene[, .(isoform_id = txname, tx_biotype)], by = "isoform_id", all.x = TRUE)
res[, library_type := if (nzchar(LIBTYPE)) LIBTYPE else "all"]
# effect-size-gated confirmed switch (review P1): stageR-confirmed AND |Δprop| >= 0.10
res[, dtu_confirmed_effect := !is.na(tx_confirm_padj) & tx_confirm_padj < 0.05 & abs(dprop) >= 0.10]

# --- empirical-null calibration diagnostics (review P0): detect null collapse/shift ---
tt <- res$t[is.finite(res$t)]
nullp <- tryCatch({ l <- locfdr::locfdr(tt, plot=0); setNames(as.numeric(l$fp0["mlest", c("delta","sigma","p0")]), c("delta","sigma","p0")) },
                  error = function(e) c(delta=median(tt), sigma=mad(tt), p0=NA_real_))
emp_frac <- mean(res$empirical_pval < 0.05, na.rm=TRUE)
null_status <- if (!is.na(nullp["sigma"]) && (nullp["sigma"] > 1.5 || (!is.na(nullp["p0"]) && nullp["p0"] > 1))) "COLLAPSED-use-theoretical" else
               if (abs(nullp["delta"]) > 0.5) "SHIFTED-design-confound" else "ok"
primary <- if (null_status == "COLLAPSED-use-theoretical") "theoretical" else "empirical"

fwrite(res, file.path(RES, sprintf("dtu_%s_%s%s.tsv", species, contrast, sfx)), sep = "\t")
ndtu  <- uniqueN(res[gene_screen_padj < 0.05, gene_id])
ndtuT <- uniqueN(res[gene_screen_padj_theo < 0.05, gene_id])
neff  <- uniqueN(res[dtu_confirmed_effect == TRUE, gene_id])
cat(sprintf("[dtu] %s/%s%s | empirical DTU genes=%d (effect-gated %d) | theoretical=%d | null delta=%.2f sigma=%.2f p0=%.2f -> %s; PRIMARY=%s\n",
            species, contrast, sfx, ndtu, neff, ndtuT, nullp["delta"], nullp["sigma"], nullp["p0"], null_status, primary))
# append diagnostics row
diag_f <- file.path(RES, "dtu_null_diagnostics.tsv")
fwrite(data.table(stratum=sprintf("%s_%s%s", species, contrast, sfx), n_tx=nrow(res),
                  emp_dtu_genes=ndtu, emp_effect_genes=neff, theo_dtu_genes=ndtuT,
                  null_delta=round(nullp["delta"],3), null_sigma=round(nullp["sigma"],3),
                  null_p0=round(nullp["p0"],3), emp_frac_lt05=round(emp_frac,3),
                  null_status=null_status, primary=primary),
       diag_f, sep="\t", append=file.exists(diag_f))
