#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fig 5e — Incretin/glucagon receptor axis across liver cell types (scRNA).
# The mechanistic keystone of the GLP-1RA section: the approved drug's receptor
# GLP1R is at the single-cell detection floor on EVERY hepatic cell type incl.
# hepatocytes (0.039% of cells), while the glucagon-receptor GCGR — the only
# incretin-axis receptor targeted by co-agonists — is strongly hepatocyte-enriched
# (17.9% of hepatocytes, 6.4x the next cell type). Same-class GPCR GCGR being
# readily detected in the SAME cells rules out dropout: GLP-1 acts off-parenchyma;
# only the glucagon arm engages the disease-executing hepatocyte directly.
#
# Dot SIZE = % cells expressing; COLOUR = mean expression (CP10k). Single-cell
# %-expressing (RevA verified read of the 1.23M-cell atlas), NOT donor pseudobulk.
# NO lollipop; all text BLACK; gene symbols italic; no title/subtitle (caption via
# message()). Data-encoding colour gradient is allowed (it is not text).
#
# Data:   RNA-seq/results/glp1ra/scrna_incretin_axis/scrna_pct_expressing_by_celltype.csv
# Output: figures/main/fig5_convergence/panels/fig5e_glp1ra_incretin_axis.pdf
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
dir.create(file.path(FIG5_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(BASE,
  "RNA-seq/results/glp1ra/scrna_incretin_axis/scrna_pct_expressing_by_celltype.csv"))

# incretin/glucagon axis only (drop lineage markers used by RevA for QC)
gene_order <- c("GLP1R", "GIPR", "GCGR", "GLP2R", "DPP4", "GCG")   # top->bottom
d <- d[gene %in% gene_order]

# cell-type order: parenchyma first, then non-parenchymal lineages
ct_order <- c("Hepatocytes", "Cholangiocytes", "Endothelial cells", "Fibroblasts",
              "Macrophages", "Mono+mono derived cells", "cDC1s", "cDC2s", "pDCs",
              "Basophils", "Neutrophils", "T cells", "Circulating NK/NKT",
              "Resident NK", "B cells", "Plasma cells")
d <- d[cell_type %in% ct_order]
d[, cell_type := factor(cell_type, levels = ct_order)]
d[, gene := factor(gene, levels = rev(gene_order))]   # rev so GLP1R plots at top

p <- ggplot(d, aes(x = cell_type, y = gene)) +
  geom_point(aes(size = pct_cells_expressing, colour = mean_cp10k)) +
  scale_size_continuous(range = c(0.05, 4.6),
                        breaks = c(0.1, 1, 5, 10, 18),
                        name = "% cells\nexpressing") +
  scale_colour_gradient(low = "grey88", high = "#08306b",
                        name = "Mean expr\n(CP10k)") +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.2))

outfile <- file.path(FIG5_DIR, "panels", "fig5e_glp1ra_incretin_axis.pdf")
ggsave(outfile, p, width = 3.5, height = 1.9, device = pdf_device)

message("Fig 5e written: ", outfile)
message("CAPTION: Incretin/glucagon-receptor axis across liver cell types (single-cell ",
        "%-expressing, n=1.23M cells). GLP1R is at the detection floor on all hepatic ",
        "cell types (hepatocytes 0.039%); GCGR is hepatocyte-enriched (17.9%, 6.4x next ",
        "type); same-class GPCR detection rules out dropout. GLP-1 acts off-parenchyma; ",
        "the glucagon (GCGR) arm engages the hepatocyte directly.")
