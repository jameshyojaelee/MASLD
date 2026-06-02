#!/usr/bin/env Rscript
# KEY MESSAGE: Hep ATAC LFC (MASLD vs Control, covariate-corrected) tracks
# bulk RNA-seq dream LFC across ALL nearest-TSS peak-gene pairs (rank-based).
# n > 1,000, no LFC thresholds.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_da_deg_concordance_nearest_tss.pdf")

MAG  <- "#C9265E"
BLUE <- "#1565C0"
GRAY <- "#9E9E9E"

# ── Hep DA (MASLD vs Control, covariate-corrected, gene-annotated) ───────
da <- fread(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv"))
da <- da[cell_type == "Hepatocyte" &
         !is.na(gene_symbol) & gene_symbol != "" &
         !grepl("^ENSG", gene_symbol) &
         abs(distance_to_tss) < 50000]
# Nearest-TSS peak per gene
da_g <- da[, .SD[which.min(abs(distance_to_tss))], by = gene_symbol]
da_g <- da_g[, .(symbol = gene_symbol, da_logFC = logFC, da_padj = padj)]

# ── Bulk dream LFC (no threshold) ────────────────────────────────────────
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv"),
  select = c("symbol", "logFC", "padj"))
dream <- dream[!is.na(symbol) & symbol != "" & !grepl("^ENSG", symbol) &
               !is.na(logFC)]

d <- merge(da_g, dream, by = "symbol")
setnames(d, c("logFC", "padj"), c("deg_logFC", "deg_padj"))

# Restrict bulk-RNA side to |LFC| >= 0.5 (Tier 1 threshold)
d <- d[abs(deg_logFC) >= 0.5]
cat(sprintf("[pairs] %d peak-gene pairs (bulk |LFC| >= 0.5)\n", nrow(d)))

rho_test <- cor.test(d$da_logFC, d$deg_logFC, method = "spearman",
                     exact = FALSE)
rho   <- unname(rho_test$estimate)
p_rho <- rho_test$p.value
cat(sprintf("[spearman] rho = %.3f, p = %.2e\n", rho, p_rho))

# Sign concordance among non-zero pairs
nonzero <- d[da_logFC != 0 & deg_logFC != 0]
sign_conc <- mean(sign(nonzero$da_logFC) == sign(nonzero$deg_logFC))
cat(sprintf("[sign concord] %d/%d (%.1f%%) across all pairs\n",
            sum(sign(nonzero$da_logFC) == sign(nonzero$deg_logFC)),
            nrow(nonzero), 100 * sign_conc))

# Label genes that are top by |joint signal| AND direction-concordant
d[, joint := abs(da_logFC) * abs(deg_logFC)]
d[, concord := sign(da_logFC) == sign(deg_logFC)]
label_dt <- d[concord == TRUE & da_padj < 0.05 & deg_padj < 0.05]
setorder(label_dt, -joint)
label_dt <- head(label_dt, 8)

# Axis limits
x_lim <- max(abs(d$deg_logFC), na.rm = TRUE) * 1.02
y_lim <- max(abs(d$da_logFC),  na.rm = TRUE) * 1.02

# Annotation: ρ and n
ann <- sprintf("Spearman ρ = %.2f\nP = %.1e\nn = %s pairs",
               rho, p_rho, format(nrow(d), big.mark = ","))

p <- ggplot(d, aes(x = deg_logFC, y = da_logFC)) +
  annotate("rect", xmin = 0, xmax = x_lim, ymin = 0, ymax = y_lim,
           fill = MAG, alpha = 0.04) +
  annotate("rect", xmin = -x_lim, xmax = 0, ymin = -y_lim, ymax = 0,
           fill = MAG, alpha = 0.04) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "gray70") +
  geom_point(colour = GRAY, alpha = 0.35, size = 0.5, shape = 16) +
  geom_point(data = label_dt, colour = MAG, size = 1.5, shape = 16) +
  geom_text_repel(data = label_dt, aes(label = symbol),
                  size = 2.4, colour = "black",
                  segment.colour = "gray50", segment.size = 0.25,
                  min.segment.length = 0, box.padding = 0.4,
                  max.overlaps = Inf) +
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
