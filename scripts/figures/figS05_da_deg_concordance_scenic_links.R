#!/usr/bin/env Rscript
# KEY MESSAGE: SCENIC+ peak-gene links give the tightest biological grounding —
# only enhancers whose accessibility is correlated with target gene RNA in
# scRNA. Concordance with bulk dream LFC tests whether single-cell-validated
# regulatory links survive in disease bulk.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_da_deg_concordance_scenic_links.pdf")

MAG  <- "#C9265E"
GRAY <- "#9E9E9E"

# ── SCENIC+ peak-gene links ───────────────────────────────────────────────
links <- fread(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/scenic_plus/enhancer_gene_links.csv"))
links[, peak_id := sprintf("%s:%d-%d", enhancer_chr, enhancer_start, enhancer_end)]

# ── Hep DA (corrected, gene-annotated, contains chrom/start/end + logFC) ──
da <- fread(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv"))
da <- da[cell_type == "Hepatocyte"]
da[, peak_id := sprintf("%s:%d-%d", chrom, start, end)]
da_keep <- da[, .(peak_id, da_logFC = logFC, da_padj = padj)]

# Match links to DA peaks
m <- merge(links, da_keep, by = "peak_id")
cat(sprintf("[SCENIC+ links matched to DA] %d / %d links\n",
            nrow(m), nrow(links)))

# ── Bulk dream LFC (no LFC filter — biological pre-filter via SCENIC+) ───
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("symbol", "logFC", "padj"))
dream <- dream[!is.na(symbol) & symbol != "" & !grepl("^ENSG", symbol) &
               !is.na(logFC)]

d <- merge(m, dream, by.x = "target_gene", by.y = "symbol")
setnames(d, c("logFC", "padj"), c("deg_logFC", "deg_padj"))
cat(sprintf("[pairs] %d SCENIC+ link-DEG pairs\n", nrow(d)))

rho_test <- cor.test(d$da_logFC, d$deg_logFC, method = "spearman", exact = FALSE)
cat(sprintf("[spearman] rho = %.3f, p = %.2e\n",
            rho_test$estimate, rho_test$p.value))

sign_conc <- mean(sign(d$da_logFC) == sign(d$deg_logFC))
cat(sprintf("[sign concord] %d/%d (%.1f%%)\n",
            sum(sign(d$da_logFC) == sign(d$deg_logFC)),
            nrow(d), 100 * sign_conc))

d[, concord := sign(da_logFC) == sign(deg_logFC)]
d[, joint := abs(da_logFC) * abs(deg_logFC)]
setorder(d, -joint)
label_dt <- head(d[concord == TRUE], 10)

x_lim <- max(abs(d$deg_logFC), na.rm = TRUE) * 1.02
y_lim <- max(abs(d$da_logFC),  na.rm = TRUE) * 1.05

p <- ggplot(d, aes(x = deg_logFC, y = da_logFC)) +
  annotate("rect", xmin = 0, xmax = x_lim, ymin = 0, ymax = y_lim,
           fill = MAG, alpha = 0.04) +
  annotate("rect", xmin = -x_lim, xmax = 0, ymin = -y_lim, ymax = 0,
           fill = MAG, alpha = 0.04) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "gray70") +
  geom_point(colour = GRAY, alpha = 0.55, size = 0.7, shape = 16) +
  geom_point(data = d[concord == TRUE], colour = MAG,
             alpha = 0.75, size = 0.9, shape = 16) +
  geom_text_repel(data = label_dt, aes(label = target_gene),
                  size = GEOM_TEXT_6PT, colour = "black",
                  segment.colour = "gray50", segment.size = 0.25,
                  min.segment.length = 0, box.padding = 0.4,
                  max.overlaps = Inf) +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     name = expression("Bulk RNA-seq  "*log[2]*"FC")) +
  scale_y_continuous(limits = c(-y_lim, y_lim),
                     name = expression("Hep ATAC  "*log[2]*"FC")) +
  theme_masld(base_size = 6) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        plot.margin = margin(6, 6, 4, 4))

ggsave(OUT_PDF, p, width = 3.6, height = 3.4, device = cairo_pdf)
cat(sprintf("Saved: %s\n", OUT_PDF))
