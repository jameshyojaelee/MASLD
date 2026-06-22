#!/usr/bin/env Rscript
# KEY MESSAGE: Cohort library-prep choices drive a 100x spread in mt% across
# the 10 cohorts (PRJNA512027 median 51% vs GSE174478 median 0.5%), but the
# spread is stable within cohort and ComBat-seq + dataset random effects in
# dream control for this batch axis.
#
# Output: figS01/panels/figS01_mt_pct.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(edgeR)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- file.path(FIGS01_DIR, "panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dge <- load_merged_dge()
md  <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))

ensg_in_dge <- tstrsplit(rownames(dge$counts), ".", fixed = TRUE, keep = 1L)[[1]]
mt_ensg <- md[chromosome == "chrM", ensembl_base]
mt_idx  <- which(ensg_in_dge %in% mt_ensg)
stopifnot(length(mt_idx) > 0)

mt_counts    <- colSums(dge$counts[mt_idx, , drop = FALSE])
total_counts <- colSums(dge$counts)
qc <- data.table(
  sample_id = rownames(dge$samples),
  dataset   = dge$samples$dataset,
  mt_pct    = 100 * mt_counts / total_counts
)

# Short cohort labels (match fig1_umap.R convention)
cohort_short <- c(
  GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
  GSE162694 = "GSE162694",   GSE167523 = "GSE167523", GSE174478 = "GSE174478",
  GSE193066 = "GSE193066", GSE213621 = "GSE213621", GSE240729 = "GSE240729",
  PRJNA512027 = "PRJNA512027"
)
qc[, cohort := cohort_short[dataset]]
ord <- qc[, .(med = median(mt_pct)), by = cohort][order(med), cohort]
qc[, cohort := factor(cohort, levels = ord)]

# Per-cohort summary table (sidecar CSV)
summ <- qc[, .(n         = .N,
               median    = round(median(mt_pct), 2),
               q25       = round(quantile(mt_pct, 0.25), 2),
               q75       = round(quantile(mt_pct, 0.75), 2),
               max       = round(max(mt_pct), 2)),
           by = cohort][order(median)]
cat("=== Per-cohort mt% summary ===\n"); print(summ)

# 10% reference line is the conventional bulk-RNA-seq QC threshold;
# above 10% suggests RNA degradation / poly-A failure.
p <- ggplot(qc, aes(x = cohort, y = mt_pct, fill = cohort)) +
  geom_violin(scale = "width", width = 0.85, alpha = 0.6, linewidth = 0.2,
              colour = "grey30") +
  geom_jitter(width = 0.12, size = 0.4, alpha = 0.55, colour = "grey20") +
  geom_hline(yintercept = 10, linetype = "dashed", colour = "grey50",
             linewidth = 0.3) +
  annotate("text", x = 0.6, y = 11, label = "10% reference",
           hjust = 0, size = 2.5, colour = "grey40") +
  scale_y_log10(breaks = c(0.1, 1, 5, 10, 30, 60, 100),
                labels = c("0.1", "1", "5", "10", "30", "60", "100")) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(x = NULL, y = "Mitochondrial reads (%, log scale)",
       title = "Per-cohort mitochondrial read fraction",
       subtitle = sprintf(
         "Computed from 34 chrM-encoded genes; %s samples; ordered by median",
         format(nrow(qc), big.mark = ","))) +
  theme_masld() +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1, size = 7),
        plot.subtitle = element_text(size = 8, color = "grey35"))

out_pdf <- file.path(OUT_DIR, "figS01_mt_pct.pdf")
out_csv <- file.path(FIGS01_DIR, "figS01_mt_pct_summary.csv")
ggsave(out_pdf, p, width = 6.5, height = 3.8, device = cairo_pdf)
fwrite(summ, out_csv)
cat("\nWrote:\n  ", out_pdf, "\n  ", out_csv, "\n", sep = "")
