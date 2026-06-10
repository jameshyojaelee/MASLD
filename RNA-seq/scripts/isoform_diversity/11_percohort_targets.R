#!/usr/bin/env Rscript
# GSE213621 per-cohort DTU target list: the effect-gated switches with their
# specific disease/control isoforms + |Δprop| + MANE/canonical/biotype.
# (The per-cohort DTU calibration bar moved to 13_supp_figures.R panel B, 2026-06-05.)
suppressPackageStartupMessages({ library(data.table) })
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/human")
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_v49_primary.tsv.gz"))

# ---- (a) GSE213621 target list ----
dtu <- fread(file.path(RES, "dtu_human_group_cohort-GSE213621.tsv"))
genes <- unique(dtu[dtu_confirmed_effect == TRUE, gene_id])
targets <- rbindlist(lapply(genes, function(g) {
  sub  <- dtu[gene_id == g]
  conf <- sub[!is.na(tx_confirm_padj) & tx_confirm_padj < 0.05]
  di <- conf[which.max(dprop)]                 # disease-gained isoform (confirmed)
  ci <- sub[which.min(dprop)]                  # most disease-lost isoform
  data.table(gene_id = g, disease_iso = di$isoform_id, disease_dprop = round(di$dprop,3),
             disease_iso_biotype = di$tx_biotype,
             control_iso = ci$isoform_id, control_dprop = round(ci$dprop,3),
             gene_screen_padj = signif(di$gene_screen_padj,3))
}))
targets <- merge(targets, tx2gene[, .(disease_iso = txname, symbol = gene_symbol,
                 gene_biotype, disease_iso_mane = is_mane_select,
                 disease_iso_canon = is_ensembl_canonical)], by = "disease_iso", all.x = TRUE)
targets[, switch_away_from_mane := disease_iso_mane == 0]
setcolorder(targets, c("symbol","gene_id","disease_iso","disease_dprop","disease_iso_biotype",
                       "disease_iso_mane","disease_iso_canon","switch_away_from_mane",
                       "control_iso","control_dprop","gene_biotype","gene_screen_padj"))
setorder(targets, -disease_dprop)
fwrite(targets, file.path(RES, "gse213621_dtu_targets.tsv"), sep = "\t")
cat(sprintf("[targets] GSE213621: %d switch genes -> gse213621_dtu_targets.tsv\n", nrow(targets)))
cat(sprintf("  %d switch AWAY from MANE; biotypes of disease isoform: %s\n",
            sum(targets$switch_away_from_mane, na.rm=TRUE),
            paste(names(sort(table(targets$disease_iso_biotype),decreasing=TRUE))[1:3], collapse=", ")))
cat("  top 10 by |Δprop|:\n"); print(head(targets[, .(symbol, disease_dprop, disease_iso_biotype, switch_away_from_mane)], 10))
