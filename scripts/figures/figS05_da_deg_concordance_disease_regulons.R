#!/usr/bin/env Rscript
# KEY MESSAGE: Concordance restricted to the 24 hep disease-regulon target
# genes — the tightest biological framing. Tests whether SCENIC+-identified
# disease-relevant chromatin tracks bulk RNA disease direction.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_da_deg_concordance_disease_regulons.pdf")

MAG  <- "#C9265E"
GRAY <- "#9E9E9E"

# ── Disease regulon target genes (parse ';'-separated target_genes column) ──
dr <- fread(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv"))
target_genes <- unique(unlist(strsplit(paste(dr$target_genes, collapse = ";"), ";")))
target_genes <- target_genes[target_genes != "" & !is.na(target_genes)]
cat(sprintf("[regulons] %d disease regulons, %d unique target genes\n",
            nrow(dr), length(target_genes)))

# ── Hep DA (corrected, gene-annotated, nearest TSS within 50 kb) ────────
da <- fread(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv"))
da <- da[cell_type == "Hepatocyte" &
         !is.na(gene_symbol) & gene_symbol != "" &
         !grepl("^ENSG", gene_symbol)]
# Keep only target-gene peaks; take nearest-TSS peak per gene
da <- da[gene_symbol %in% target_genes]
da_g <- da[, .SD[which.min(abs(distance_to_tss))], by = gene_symbol]
da_g <- da_g[, .(symbol = gene_symbol, da_logFC = logFC, da_padj = padj)]
cat(sprintf("[DA] %d disease-regulon target genes with a Hep DA peak\n",
            nrow(da_g)))

# ── Bulk dream LFC (no filter; biology pre-filtered by disease regulon) ──
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("symbol", "logFC", "padj"))
dream <- dream[!is.na(symbol) & symbol != "" & !grepl("^ENSG", symbol) &
               !is.na(logFC)]

d <- merge(da_g, dream, by = "symbol")
setnames(d, c("logFC", "padj"), c("deg_logFC", "deg_padj"))
cat(sprintf("[pairs] %d disease-regulon target genes\n", nrow(d)))

if (nrow(d) >= 3) {
  rho_test <- cor.test(d$da_logFC, d$deg_logFC, method = "spearman", exact = FALSE)
  cat(sprintf("[spearman] rho = %.3f, p = %.2e\n",
              rho_test$estimate, rho_test$p.value))
}
sign_conc <- mean(sign(d$da_logFC) == sign(d$deg_logFC))
cat(sprintf("[sign concord] %d/%d (%.1f%%)\n",
            sum(sign(d$da_logFC) == sign(d$deg_logFC)),
            nrow(d), 100 * sign_conc))

d[, concord := sign(da_logFC) == sign(deg_logFC)]

x_lim <- max(abs(d$deg_logFC), na.rm = TRUE) * 1.02
y_lim <- max(abs(d$da_logFC),  na.rm = TRUE) * 1.05

p <- ggplot(d, aes(x = deg_logFC, y = da_logFC)) +
  annotate("rect", xmin = 0, xmax = x_lim, ymin = 0, ymax = y_lim,
           fill = MAG, alpha = 0.04) +
  annotate("rect", xmin = -x_lim, xmax = 0, ymin = -y_lim, ymax = 0,
           fill = MAG, alpha = 0.04) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "gray70") +
  geom_point(aes(colour = concord), size = 1.6, shape = 16) +
  geom_text_repel(aes(label = symbol),
                  size = 2.4, colour = "black",
                  segment.colour = "gray50", segment.size = 0.25,
                  min.segment.length = 0, box.padding = 0.4,
                  max.overlaps = Inf) +
  scale_colour_manual(values = c(`TRUE` = MAG, `FALSE` = GRAY),
                      guide = "none") +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     name = expression(bold("Bulk RNA-seq  "*log[2]*"FC"))) +
  scale_y_continuous(limits = c(-y_lim, y_lim),
                     name = expression(bold("Hep ATAC  "*log[2]*"FC"))) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        plot.margin = margin(6, 6, 4, 4))

ggsave(OUT_PDF, p, width = 3.6, height = 3.4, device = cairo_pdf)
cat(sprintf("Saved: %s\n", OUT_PDF))
