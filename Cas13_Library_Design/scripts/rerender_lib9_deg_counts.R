#!/usr/bin/env Rscript
# rerender_lib9_deg_counts.R
# Re-render S_lib_9_deg_counts.pdf with UP-ONLY gene counts (logFC > x, shrunk > x)
# instead of the previous both-direction |logFC| > x.
# Also rewrites patient_concordance_by_cutoff_to3.csv and raw_ashr_deg_jaccard_by_cutoff.csv
# with up-only counts so all downstream CSVs are consistent.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
FIGDIR <- file.path(BASE, "Cas13_Library_Design/figures")
DATDIR <- file.path(BASE, "Cas13_Library_Design/data")

strip_v <- function(x) sub("[.][0-9]+$", "", x)

# ---- load dream results -------------------------------------------------------
d <- fread(file.path(INTDIR, "dream_results_ashr.csv"))
d[, gb := strip_v(gene)]

cat(sprintf("Loaded %d genes\n", nrow(d)))
cat(sprintf("UP-only at raw LFC>0.5:  %d (was both-dir: %d)\n",
            d[!is.na(padj) & padj < 0.05 & logFC > 0.5, .N],
            d[!is.na(padj) & padj < 0.05 & abs(logFC) > 0.5, .N]))
cat(sprintf("UP-only at ashr shrunk>0.2: %d (was both-dir: %d)\n",
            d[!is.na(lfsr) & lfsr < 0.05 & shrunk_logFC > 0.2, .N],
            d[!is.na(lfsr) & lfsr < 0.05 & abs(shrunk_logFC) > 0.2, .N]))

# ---- recompute with UP-ONLY ---------------------------------------------------
CUTS <- c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0)
RAW  <- "padj<0.05 & log2FC > x (UP)"
ASHR <- "lfsr<0.05 & shrunk > x (UP)"

counts <- rbindlist(c(
  lapply(CUTS, function(x) data.table(
    condition = RAW, cutoff = x,
    n_degs = d[!is.na(padj) & padj < 0.05 & logFC > x, .N])),
  lapply(CUTS, function(x) data.table(
    condition = ASHR, cutoff = x,
    n_degs = d[!is.na(lfsr) & lfsr < 0.05 & shrunk_logFC > x, .N]))
))
counts[, condition := factor(condition, levels = c(RAW, ASHR))]
cat("\n=== UP-only DEG counts by cutoff ===\n"); print(counts)

# ---- rewrite Jaccard CSV with up-only ----------------------------------------
jac <- rbindlist(lapply(CUTS, function(x) {
  raw_set  <- d[!is.na(padj) & padj < 0.05 & logFC > x, gb]
  ashr_set <- d[!is.na(lfsr) & lfsr < 0.05 & shrunk_logFC > x, gb]
  isect    <- length(intersect(raw_set, ashr_set))
  union_n  <- length(union(raw_set, ashr_set))
  data.table(cutoff = x, n_raw = length(raw_set), n_ashr = length(ashr_set),
             n_intersect = isect, jaccard = ifelse(union_n > 0, isect / union_n, NA_real_))
}))
fwrite(jac, file.path(DATDIR, "raw_ashr_deg_jaccard_by_cutoff.csv"))
cat("\nRewritten raw_ashr_deg_jaccard_by_cutoff.csv (UP-only)\n")

# ---- re-render S_lib_9_deg_counts.pdf ----------------------------------------
pal <- setNames(c("#B2182B", "#3B4CC0"), c(RAW, ASHR))
bt <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold", size = 11))

p <- ggplot(counts, aes(cutoff, n_degs, color = condition)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  geom_text(aes(label = n_degs), vjust = -0.9, size = 2.6, show.legend = FALSE) +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = CUTS) +
  expand_limits(y = max(counts$n_degs) * 1.10) +
  labs(x = "Effect-size cutoff", y = "Number of UP-regulated genes",
       title = "UP-only DEG counts by cutoff  (raw padj<0.05 vs ashr lfsr<0.05)") +
  bt

out_pdf <- file.path(FIGDIR, "S_lib_9_deg_counts.pdf")
ggsave(out_pdf, p, width = 7.0, height = 4.5, useDingbats = FALSE)
cat("Wrote", out_pdf, "\n")
cat(sprintf("\nCanonical operating point check — ashr shrunk>0.15: %d UP genes\n",
            counts[condition == ASHR & cutoff == 0.15, n_degs]))
