#!/usr/bin/env Rscript
# aggregate_pillar_A_lto.R
# ---------------------------------------------------------------------------
# Aggregate dream LTO (leave-2-cohorts-out) results.  Mirrors the LOCO
# aggregator in shape: per-pair recovery of canonical DEGs + per-gene
# fold-recurrence count.  Threshold sweep so the LTO numbers can be
# reported alongside the LFC sweep on Pillar A.
#
# Out: audit_sensitivity/pillar_A_lto_summary_<tag>.csv          (per-pair)
#      audit_sensitivity/pillar_A_lto_per_gene_<tag>.csv          (per-gene)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(yaml) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR <- file.path(RDIR, "loo_cv")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_THR <- 0.05
LFC_THR  <- as.numeric(Sys.getenv("LFC_THR", "0.5"))
THR_TAG  <- if (LFC_THR == 0) "nolfc" else paste0("lfc", sub("^0\\.", "", as.character(LFC_THR)))
cat(sprintf("LTO aggregation at padj<%.2f & |logFC|>%.2f  (tag = %s)\n", PADJ_THR, LFC_THR, THR_TAG))

# Canonical full dream
full <- fread(file.path(RDIR, "dream_results.csv"))
deg_full <- full[padj < PADJ_THR & abs(logFC) > LFC_THR, gene]
N_DEG <- length(deg_full)
cat(sprintf("Canonical DEGs at this threshold: %d\n", N_DEG))

# Read all 10 LTO results
lto_files <- list.files(LOO_DIR, pattern = "^dream_lto_pair\\d{2}_held_.*\\.csv$", full.names = TRUE)
cat(sprintf("LTO files found: %d / 10 expected\n", length(lto_files)))
if (length(lto_files) == 0) stop("No LTO files found — run array first.")

per_pair <- rbindlist(lapply(lto_files, function(f) {
  pair <- as.integer(sub("^dream_lto_pair(\\d{2})_held_.*$", "\\1", basename(f)))
  held <- sub("^dream_lto_pair\\d{2}_held_(.+)\\.csv$", "\\1", basename(f))
  dt <- fread(f, select = c("gene", "logFC", "t", "padj"))
  if (nrow(dt) == 0) {
    return(data.table(pair = pair, held = held, n_DEG_lto = 0L,
                      pct_recovered = NA_real_, rho = NA_real_, jaccard = NA_real_))
  }
  sel <- dt[padj < PADJ_THR & abs(logFC) > LFC_THR, gene]
  m <- merge(full[, .(gene, full_logFC = logFC)],
             dt[, .(gene, lto_logFC = logFC)], by = "gene")
  rho <- suppressWarnings(cor(m$full_logFC, m$lto_logFC, method = "spearman"))
  inter <- length(intersect(deg_full, sel))
  union_n <- length(union(deg_full, sel))
  data.table(pair = pair, held = held,
             n_DEG_lto = length(sel),
             pct_recovered = round(if (N_DEG > 0) inter / N_DEG * 100 else 0, 2),
             rho = round(rho, 4),
             jaccard = round(if (union_n > 0) inter / union_n else 0, 4))
}))[order(pair)]
fwrite(per_pair, file.path(OUT_DIR, sprintf("pillar_A_lto_summary_%s.csv", THR_TAG)))

cat("\n=== LTO per-pair ===\n"); print(per_pair)
cat(sprintf("\nMean recovery: %.1f %%   Mean rho: %.4f   Mean Jaccard: %.4f\n",
            mean(per_pair$pct_recovered, na.rm = TRUE),
            mean(per_pair$rho, na.rm = TRUE),
            mean(per_pair$jaccard, na.rm = TRUE)))
cat(sprintf("Range recovery: %.1f %% - %.1f %%\n",
            min(per_pair$pct_recovered, na.rm = TRUE),
            max(per_pair$pct_recovered, na.rm = TRUE)))

# --- Per-gene LTO recurrence: in how many of 10 LTO folds was the gene a DEG? ---
lto_long <- rbindlist(lapply(lto_files, function(f) {
  dt <- fread(f, select = c("gene", "logFC", "padj"))
  dt[, sel := padj < PADJ_THR & abs(logFC) > LFC_THR]
  dt[, .(gene, sel)]
}))
per_gene_lto <- lto_long[, .(lto_recur = sum(sel, na.rm = TRUE),
                              n_lto_tested = .N), by = gene]
per_gene_lto <- merge(full[, .(gene)], per_gene_lto, by = "gene", all.x = TRUE)
per_gene_lto[is.na(lto_recur), `:=`(lto_recur = 0L, n_lto_tested = 0L)]
per_gene_lto[, is_canonical_DEG := gene %in% deg_full]

fwrite(per_gene_lto, file.path(OUT_DIR, sprintf("pillar_A_lto_per_gene_%s.csv", THR_TAG)))
deg_pg <- per_gene_lto[is_canonical_DEG == TRUE]
if (nrow(deg_pg) > 0) {
  cat(sprintf("\nPer-canonical-DEG LTO recurrence (n=%d):\n", nrow(deg_pg)))
  cat(sprintf("  Mean recurrence: %.2f / 10\n", mean(deg_pg$lto_recur)))
  cat(sprintf("  Recur >= 8/10:    %d (%.1f%%)\n",
              sum(deg_pg$lto_recur >= 8), 100 * mean(deg_pg$lto_recur >= 8)))
  cat(sprintf("  Recur == 10/10:   %d (%.1f%%)\n",
              sum(deg_pg$lto_recur == 10), 100 * mean(deg_pg$lto_recur == 10)))
}
cat("\nDone LTO aggregation.\n")
