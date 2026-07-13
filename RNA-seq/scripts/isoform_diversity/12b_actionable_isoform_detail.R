#!/usr/bin/env Rscript
# Per-transcript isoform breakdown for the isoform-selective Cas13 targets. Targets are read
# DYNAMICALLY from the targetability table (isoform_selective_possible==TRUE) so this tracks the
# current target set automatically. Emits every transcript's mouse-liver TPM / within-gene
# proportion / expressed-flag / canonical-flag, for tiling. Library-design only (not manuscript).
suppressPackageStartupMessages({ library(data.table) })
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse")
strip <- function(x) sub("\\..*$","",x)
TPM_FLOOR <- 0.5     # matches 12_mouse_isoform_targetability.R expressed-isoform floor

tg <- fread(file.path(MRES, "gse213621_targets_mouse_isoform_structure.tsv"))[isoform_selective_possible == TRUE]
cat(sprintf("[12b] %d isoform-selective targets: %s\n", nrow(tg), paste(tg$symbol, collapse=", ")))

mtpm <- readRDS(file.path(MRES, "tx_tpm.rds"))
msamp<- fread(file.path(MRES, "samples.tsv")); setkey(msamp, sample_id); msamp <- msamp[colnames(mtpm)]
mtpm <- mtpm[, msamp$group_binary %in% c("Control","Disease")]
mt2g <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_vM38_primary.tsv.gz"))
mt2g <- mt2g[match(rownames(mtpm), txname)]; mt2g[, mensg := strip(geneid)]
meanTPM <- rowMeans(mtpm)

out <- rbindlist(lapply(seq_len(nrow(tg)), function(i) {
  idx <- which(mt2g$mensg == strip(tg$mouse_ensg[i])); tp <- meanTPM[idx]; prop <- tp / sum(tp)
  data.table(human=tg$symbol[i], mouse=tg$mouse_symbol[i], transcript=rownames(mtpm)[idx],
             mean_TPM=round(tp,2), prop=round(prop,3),
             expressed = prop>=0.05 & tp>=TPM_FLOOR, canonical = mt2g$is_ensembl_canonical[idx]==1)[order(-prop)]
}))
fwrite(out, file.path(MRES, "actionable_isoform_detail.tsv"), sep="\t")
for (g in tg$symbol) {
  d <- out[human==g & (prop>=0.01 | expressed)]
  cat(sprintf("\n=== %s -> %s  (%d tx total, %d shown) ===\n", g, d$mouse[1], nrow(out[human==g]), nrow(d)))
  print(d[, .(transcript, mean_TPM, prop, expressed, canonical)], row.names=FALSE)
}
cat("\n-> actionable_isoform_detail.tsv\n")
