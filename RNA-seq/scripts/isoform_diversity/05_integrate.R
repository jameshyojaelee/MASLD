#!/usr/bin/env Rscript
# Phase 3 integration — fuse isoform diversity + DTU + dream DEGs into one per-gene
# table with the DEG x DTU quadrant and the disease/control-dominant isoform.
#
# Usage: Rscript 05_integrate.R <species: human|mouse>
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
species <- if (length(args) >= 1) args[1] else "human"
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RES  <- file.path(PROJ, "RNA-seq/results/isoform_diversity", species)
strip <- function(x) sub("\\..*$", "", x)

div  <- fread(file.path(RES, sprintf("isoform_diversity_%s.tsv", species)))
# --- canonical DTU stratum for integration (review 2026-06-05) -------------------------
# 05 previously read the POOLED dtu_<sp>_group.tsv; for human that pools polyA + total-RNA
# with no library_type covariate, so its empirical null COLLAPSES (chemistry-confounded;
# dtu_null_diagnostics flags null_status=COLLAPSED) and the quadrant table inherited the
# confound. Default to the well-calibrated, chemistry-matched total-RNA stratum (the same
# canonical human DTU used by 09_crossspecies; null_status=ok, primary=empirical). Override
# with the DTU_STRATUM env var (e.g. "group_cohort-GSE213621", "group").
DTU_STRATUM <- Sys.getenv("DTU_STRATUM", if (species == "human") "group_total-RNA" else "group")
dtu  <- fread(file.path(RES, sprintf("dtu_%s_%s.tsv", species, DTU_STRATUM)))
# Calibration guard: integration gates dtu_pos on EMPIRICAL stageR p-values (conservative).
# We do NOT auto-substitute the theoretical null for collapsed strata — it over-calls on
# heterogeneous pools (e.g. mouse_group theoretical = 4,665 vs empirical 0). Instead, warn
# loudly so a non-ok stratum is a visible choice, not a silent one.
diag_dt <- tryCatch(fread(file.path(RES, "dtu_null_diagnostics.tsv")), error = function(e) NULL)
nstat   <- if (!is.null(diag_dt)) diag_dt[stratum == sprintf("%s_%s", species, DTU_STRATUM), null_status][1] else NA
if (length(nstat) && !is.na(nstat) && nstat != "ok")
  cat(sprintf("[integrate] WARNING: DTU stratum %s_%s null_status='%s' (not 'ok'); dtu_pos on empirical p (conservative).\n",
              species, DTU_STRATUM, nstat))
cat(sprintf("[integrate] DTU stratum = %s_%s (null_status=%s)\n", species, DTU_STRATUM,
            if (length(nstat) && !is.na(nstat)) nstat else "unknown"))
stagef <- file.path(RES, sprintf("dtu_%s_stage.tsv", species))
dtu_stage <- if (file.exists(stagef)) fread(stagef) else NULL
t2g_f <- if (species == "human") "tx2gene_v49_primary.tsv.gz" else "tx2gene_vM38_primary.tsv.gz"
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index", t2g_f))
tpm  <- readRDS(file.path(RES, "tx_tpm.rds"))
cts_dtu <- readRDS(file.path(RES, "tx_counts_dtuscaled.rds"))   # tested quantity for isoform naming
samp <- fread(file.path(RES, "samples.tsv")); setkey(samp, sample_id); samp <- samp[colnames(tpm)]

# --- DEG layer (human only: dream; mouse: skip) ---
div[, ensg := strip(gene)]
if (species == "human") {
  dream <- fread(file.path(PROJ,
      "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv"))
  dream[, ensg := strip(gene)]
  dream[, tier1 := !is.na(lfsr) & lfsr < 0.05 & abs(shrunk_logFC) > 0.5]
  dream[, deg_dir := ifelse(shrunk_logFC > 0, "up", "down")]
  deg <- dream[, .(ensg, dream_shrunk_logFC = shrunk_logFC, dream_lfsr = lfsr,
                   tier1, deg_dir)]
  # deg merged onto the full gene UNION after the DTU merge (review P1: keeps DTU+ genes
  # that fail the diversity expression filter, and lets us flag deg_tested vs untested)
} else {
  deg <- NULL; div[, `:=`(tier1 = NA, deg_dir = NA_character_)]
}

# --- DTU gene-level (group) + disease/control-dominant isoform ---
dtu[, ensg := strip(gene_id)]
gene_dtu <- dtu[, .(dtu_gene_padj = suppressWarnings(min(gene_screen_padj, na.rm = TRUE))),
                by = ensg]
gene_dtu[is.infinite(dtu_gene_padj), dtu_gene_padj := NA]
gene_dtu[, dtu_pos := !is.na(dtu_gene_padj) & dtu_gene_padj < 0.05]

# disease/control isoform from the TESTED dtuScaledTPM (review P1: raw-TPM argmax matched
# the satuRn-confirmed transcript only 40% of the time); name the satuRn-confirmed tx when one exists
grp <- factor(samp$group_binary, levels = c("Control","Disease"))
# Restrict isoform-naming samples to the SAME stratum the DTU test used, so the
# disease/control-isoform proportions are chemistry-matched to the test (not pooled).
smask <- rep(TRUE, nrow(samp))
if (grepl("total-RNA", DTU_STRATUM)) {
  smask <- !is.na(samp$library_type) & samp$library_type == "total-RNA"
} else if (grepl("polyA", DTU_STRATUM)) {
  smask <- !is.na(samp$library_type) & samp$library_type == "polyA"
} else if (grepl("cohort-", DTU_STRATUM)) {
  smask <- samp$dataset == sub(".*cohort-", "", DTU_STRATUM)
}
keep_s <- grp %in% c("Control","Disease") & smask &
          (is.na(samp$pass_technical) | samp$pass_technical %in% c(TRUE,"TRUE","True"))
cdtu <- cts_dtu[, keep_s]; grp <- droplevels(grp[keep_s])
tx2gene2 <- tx2gene[match(rownames(cdtu), tx2gene$txname)]
confirmed <- dtu[!is.na(tx_confirm_padj) & tx_confirm_padj < 0.05, isoform_id]
iso_calls <- rbindlist(lapply(gene_dtu[dtu_pos == TRUE, ensg], function(eg) {
  idx <- which(strip(tx2gene2$geneid) == eg & rownames(cdtu) %in% dtu$isoform_id)
  if (length(idx) < 2) return(NULL)
  sub <- cdtu[idx, , drop = FALSE]; ids <- rownames(sub)
  pc  <- rowMeans(sub[, grp=="Control", drop=FALSE]); pc <- pc / sum(pc)
  pd  <- rowMeans(sub[, grp=="Disease", drop=FALSE]); pd <- pd / sum(pd)
  dd  <- pd - pc
  ci  <- which(ids %in% confirmed)                        # satuRn-confirmed transcripts of this gene
  di  <- if (length(ci)) ci[which.max(dd[ci])] else which.max(dd)   # disease iso: confirmed-preferred, max gain
  data.table(ensg = eg,
             disease_iso = ids[di], disease_iso_dprop = dd[di],
             control_iso = ids[which.min(dd)], control_iso_dprop = min(dd))
}), fill = TRUE)
if (nrow(iso_calls) > 0 && "ensg" %in% names(iso_calls)) {
  iso_calls <- merge(iso_calls,
      tx2gene[, .(disease_iso = txname, disease_iso_mane = is_mane_select,
                  disease_iso_canon = is_ensembl_canonical)], by = "disease_iso", all.x = TRUE)
  iso_calls[, switch_away_from_mane := disease_iso_mane == 0]
  gene_dtu <- merge(gene_dtu, iso_calls, by = "ensg", all.x = TRUE)
}  # else: no DTU+ genes (e.g. mouse) -> no disease/control-isoform columns

# --- stage-axis DTU flag ---
if (!is.null(dtu_stage)) {
  dtu_stage[, ensg := strip(gene_id)]
  gs <- dtu_stage[, .(dtu_stage_padj = suppressWarnings(min(gene_screen_padj, na.rm=TRUE))), by=ensg]
  gs[is.infinite(dtu_stage_padj), dtu_stage_padj := NA]
  gs[, dtu_stage_pos := !is.na(dtu_stage_padj) & dtu_stage_padj < 0.05]
  gene_dtu <- merge(gene_dtu, gs, by = "ensg", all.x = TRUE)
}

out <- merge(div, gene_dtu, by = "ensg", all = TRUE)      # full outer: keep DTU+ genes absent from diversity
out[is.na(dtu_pos), dtu_pos := FALSE]

# --- quadrant + focus tags ---
if (species == "human") {
  out <- merge(out, deg, by = "ensg", all.x = TRUE)       # tier1/deg_dir for the full union
  out[, deg_tested := ensg %in% dream$ensg]               # gene present in dream's tested universe
  out[is.na(tier1), tier1 := FALSE]                       # NA tier1 (untested) -> not a Tier1 DEG
  out[, quadrant := fifelse(tier1 & dtu_pos, "DEG+DTU+",
                    fifelse(tier1 & !dtu_pos, "DEG_only",
                    fifelse(!tier1 & dtu_pos &  deg_tested, "pure_switch",
                    fifelse(!tier1 & dtu_pos & !deg_tested, "switch_DEuntested", "neither"))))]
  # effect-size-gated pure switch (review P1: 62% of pure_switch have |Δprop|<0.1)
  out[, pure_switch_strong := quadrant == "pure_switch" &
        !is.na(disease_iso_dprop) & abs(disease_iso_dprop) >= 0.10]
  pc <- tryCatch(fread(file.path(PROJ, "results/library/positive_control.csv")), error=function(e) NULL)
  if (!is.null(pc)) out[, is_positive_control := symbol %in% pc[["Gene symbol"]]]
}

fwrite(out, file.path(RES, sprintf("isoform_master_%s.tsv", species)), sep = "\t")
cat(sprintf("[integrate] %s: %d genes; DTU+ %d; ", species, nrow(out), sum(out$dtu_pos, na.rm=TRUE)))
if (species == "human") cat(sprintf("quadrants: %s\n",
    paste(names(table(out$quadrant)), table(out$quadrant), sep="=", collapse=" ")))
cat("\n[integrate] wrote", sprintf("isoform_master_%s.tsv", species), "\n")
