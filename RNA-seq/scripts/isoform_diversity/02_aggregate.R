#!/usr/bin/env Rscript
# Phase 1 — aggregate per-sample kallisto into transcript-level matrices.
# One tximport gives BOTH: $abundance = TPM (diversity) and $counts = dtuScaledTPM
# (the DTU-correct length scaling, Love/Soneson/Patro 2018). Bootstraps for swish
# are imported later on the FILTERED transcript set (full infReps are too large).
#
# Usage: Rscript 02_aggregate.R <species: human|mouse>
suppressPackageStartupMessages({
  library(tximport); library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
species <- if (length(args) >= 1) args[1] else "human"
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ID   <- file.path(PROJ, "RNA-seq/scripts/isoform_diversity")
RES  <- file.path(PROJ, "RNA-seq/results/isoform_diversity", species)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

if (species == "human") {
  driver  <- fread(file.path(ID, "human_driver.tsv"))
  driver[, subdir := cohort]
  tx2gene <- fread(file.path(PROJ,
      "RNA-seq/results/isoform_diversity/_index/tx2gene_v49_primary.tsv.gz"))
  meta    <- fread(file.path(PROJ,
      "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))
  qc      <- fread(file.path(PROJ,
      "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"))
} else {  # mouse
  driver  <- fread(file.path(ID, "mouse_driver.tsv"))
  driver[, `:=`(subdir = dataset, cohort = dataset)]
  tx2gene <- fread(file.path(PROJ,
      "RNA-seq/results/isoform_diversity/_index/tx2gene_vM38_primary.tsv.gz"))
  meta    <- fread(file.path(PROJ,
      "RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv"))
  qc      <- NULL
}
kroot <- file.path(RES, "kallisto")

# --- file list (one abundance.tsv per sample; kallisto 0.51 build lacks HDF5) ---
files <- file.path(kroot, driver$subdir, driver$sample_id, "abundance.tsv")
names(files) <- driver$sample_id
ok <- file.exists(files)
cat(sprintf("[aggregate] %d/%d abundance.tsv present\n", sum(ok), length(files)))
if (!all(ok)) cat("[aggregate] MISSING:", paste(head(driver$sample_id[!ok], 10), collapse=","), "\n")
files <- files[ok]; driver <- driver[ok]

# tx2gene keyed on bare versioned ENST; kallisto target_id is GENCODE pipe header
# 'ENST|ENSG|...', so ignoreAfterBar=TRUE strips to the ENST to match.
t2g <- tx2gene[, .(txname, geneid)]

cat("[aggregate] tximport (dtuScaledTPM)...\n")
txi <- tximport(files, type = "kallisto", txOut = TRUE,
                countsFromAbundance = "dtuScaledTPM",
                tx2gene = t2g, ignoreAfterBar = TRUE,
                ignoreTxVersion = FALSE, dropInfReps = TRUE)

# Defensive: ensure rownames are bare versioned ENST (strip GENCODE pipe header
# if ignoreAfterBar did not propagate to the txOut matrix rownames).
tx_tpm   <- txi$abundance       # TPM  (diversity)
tx_dtu   <- txi$counts          # dtuScaledTPM counts (DTU)
rownames(tx_tpm) <- sub("\\|.*$", "", rownames(tx_tpm))
rownames(tx_dtu) <- sub("\\|.*$", "", rownames(tx_dtu))
stopifnot(all(rownames(tx_tpm) == rownames(tx_dtu)))
cat(sprintf("[aggregate] matrix: %d transcripts x %d samples\n", nrow(tx_tpm), ncol(tx_tpm)))
matched <- mean(rownames(tx_tpm) %in% tx2gene$txname)
cat(sprintf("[aggregate] %.3f of transcripts found in tx2gene\n", matched))

# --- per-sample library type from non-polyadenylated small-ncRNA fraction ---
# polyA selection strips snoRNA/snRNA/scaRNA; ribo-depletion/total-RNA retains them.
bt    <- tx2gene$tx_biotype[match(rownames(tx_tpm), tx2gene$txname)]
small <- bt %in% c("snoRNA","snRNA","scaRNA")
# library prep is per-cohort, so classify each DATASET by its median small-ncRNA
# fraction (avoids spurious within-cohort splits from biological spread crossing a
# per-sample cut). Threshold 8% cleanly separates polyA (median<3%) from total-RNA.
small_frac <- colSums(tx_tpm[small, , drop = FALSE]) / colSums(tx_tpm) * 100
ltab <- data.table(sample_id = colnames(tx_tpm), small_ncrna_frac = round(small_frac, 3))

# --- sample table: join metadata + QC (robust to species-specific columns) ---
want <- c("sample_id","dataset","group_binary","sex","age","fibrosis_stage","nas_score","diet_model")
for (cc in setdiff(want, names(meta))) meta[[cc]] <- NA      # add missing cols (mouse)
samp <- merge(data.table(sample_id = colnames(tx_tpm)),
              meta[, ..want], by = "sample_id", all.x = TRUE, sort = FALSE)
if (!is.null(qc)) {
  samp <- merge(samp, qc[, .(sample_id, pass_technical, total_counts)],
                by = "sample_id", all.x = TRUE, sort = FALSE)
} else { samp[, `:=`(pass_technical = NA, total_counts = NA)] }
samp <- merge(samp, ltab, by = "sample_id", all.x = TRUE, sort = FALSE)
samp[, lib_median_frac := median(small_ncrna_frac, na.rm = TRUE), by = dataset]
samp[, library_type := ifelse(lib_median_frac < 8, "polyA", "total-RNA")]
setkey(samp, sample_id); samp <- samp[colnames(tx_tpm)]   # preserve matrix col order
cat("[aggregate] group_binary x pass_technical:\n"); print(table(samp$group_binary, samp$pass_technical, useNA="ifany"))
cat("[aggregate] library_type by dataset (cohort-median):\n"); print(table(samp$dataset, samp$library_type, useNA="ifany"))

# --- save ---
saveRDS(tx_tpm, file.path(RES, "tx_tpm.rds"))
saveRDS(tx_dtu, file.path(RES, "tx_counts_dtuscaled.rds"))
fwrite(samp, file.path(RES, "samples.tsv"), sep = "\t")
# lightweight provenance
prov <- data.table(species = species, n_tx = nrow(tx_tpm), n_samp = ncol(tx_tpm),
                   index = "gencode.v49.primary.kidx (507365 tx)", strand = "rf-stranded",
                   bootstraps = 0, scaling = "dtuScaledTPM",
                   note = "kallisto 0.51.1 build lacks HDF5; point estimates only; swish/inferential DTU N/A")
fwrite(prov, file.path(RES, "aggregate_provenance.tsv"), sep = "\t")
cat("[aggregate] DONE ->", RES, "\n")
