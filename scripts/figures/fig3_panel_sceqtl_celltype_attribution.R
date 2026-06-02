#!/usr/bin/env Rscript
# fig3_panel_sceqtl_celltype_attribution.R
#
# Cell-type attribution dot matrix: top sc-eQTL COLOC genes (UKBB ALT GWAS)
# across 4 liver cell types. Dot size = PP4. Color = cell type with best PP4.
# Highlights that genetic risk for MASLD operates through specific cell types —
# information invisible in bulk eQTL.
#
# Output: figures/main/fig3_regulatory_architecture/panels/fig3_panel_sceqtl_celltype_attribution.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")
OUT_PDF   <- file.path(PANEL_DIR, "fig3_panel_sceqtl_celltype_attribution.pdf")
SCEQTL    <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl")

# ── load all cell-type COLOC results (UKBB ALT GWAS) ─────────────────────────
ct_labels <- c(hepatocyte       = "Hepatocyte",
               cholangiocyte    = "Cholangiocyte",
               endothelial_cell = "Endothelial",
               stellate_cell    = "Stellate")

ct_pal <- c(
  Hepatocyte    = masld_colors$fibrosis,  # deep magenta
  Cholangiocyte = "#1565C0",             # deep blue
  Endothelial   = "#00695C",             # teal
  Stellate      = "#F57F17"              # amber
)

all_ct <- rbindlist(lapply(names(ct_labels), function(ct) {
  f <- file.path(SCEQTL, paste0("coloc_results_ukbb_alt_", ct, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  pp4_col <- intersect(c("PP.H4_all", "PP.H4"), names(d))[1]
  d[, .(gene, cell_type = ct_labels[ct], PP4 = get(pp4_col))]
}))

# ── select top genes by best PP4 across cell types ───────────────────────────
best <- all_ct[, .(best_PP4 = max(PP4, na.rm = TRUE)), by = gene]
TOP_N <- 18
EXCLUDE <- c("NT5C", "LCOR", "PRRC1", "PCCB",
             "GTF2B", "DAGLB", "SOX6", "ADAM15", "DANCR", "USP43", "CYFIP2",
             # uninformative lncRNA/AC accession names
             "AL162411.1", "AC020656.2", "AC091114.1", "AC040934.1",
             "AC253536.6", "AC092650.1")
top_genes <- best[best_PP4 > 0.5 & !gene %in% EXCLUDE
                  ][order(-best_PP4)][seq_len(min(TOP_N, .N)), gene]
cat(sprintf("[genes] %d with PP4 > 0.5, showing top %d\n", nrow(best[best_PP4>0.5]), TOP_N))

# full grid: all genes × all cell types
plot_dt <- all_ct[gene %in% top_genes]
plot_dt <- merge(
  CJ(gene = top_genes, cell_type = unname(ct_labels)),
  plot_dt, by = c("gene", "cell_type"), all.x = TRUE)
plot_dt[is.na(PP4), PP4 := 0]

# gene order: by best PP4 descending
gene_order <- best[gene %in% top_genes][order(-best_PP4), gene]
plot_dt[, gene := factor(gene, levels = rev(gene_order))]
plot_dt[, cell_type := factor(cell_type,
  levels = c("Hepatocyte", "Cholangiocyte", "Endothelial", "Stellate"))]

# known MASLD genes to bold
known <- c("HSD17B13", "SPTLC3", "EFHD1", "CHEK2", "PNPLA3", "TM6SF2",
           "MARC1", "GCKR", "FABP1", "HKDC1", "THRB", "SLC12A8", "MLIP")
plot_dt[, gene_face := ifelse(as.character(gene) %in% known, "bold.italic", "plain")]

y_faces <- ifelse(rev(gene_order) %in% known, "bold.italic", "plain")

# ── plot ──────────────────────────────────────────────────────────────────────
n_sig <- nrow(best[best_PP4 > 0.5])

p <- ggplot(plot_dt, aes(x = cell_type, y = gene, fill = PP4)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "gray96", high = masld_colors$up,
                      limits = c(0, 1),
                      breaks = c(0, 0.5, 1),
                      name   = "COLOC\nPP4") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Cell-type-resolved genetic risk for MASLD") +
  theme_masld(base_size = 11) +
  theme(axis.text.x      = element_text(size = 10, face = "bold",
                                         angle = 35, hjust = 1, vjust = 1),
        axis.text.y      = element_text(size = 8.5, face = y_faces),
        panel.grid        = element_blank(),
        axis.line         = element_blank(),
        axis.ticks        = element_blank(),
        legend.position   = "right",
        legend.key.height = unit(1.2, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        legend.text       = element_text(size = 9),
        legend.title      = element_text(size = 9),
        plot.title        = element_text(size = 12, face = "bold",
                                         margin = margin(b = 8)))

ggsave(OUT_PDF, p,
       width  = 5.5,
       height = 2 + TOP_N * 0.28,
       device = cairo_pdf)
cat(sprintf("Saved: %s\n", OUT_PDF))
