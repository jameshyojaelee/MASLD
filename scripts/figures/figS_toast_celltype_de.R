#!/usr/bin/env Rscript
# figS_toast_celltype_de.R
# Figure for A1 v2 — TOAST bulk-level cell-type-specific DE.
#
# Panels:
#   A — Bar: TOAST DEG counts per cell type
#   B — Volcano plot per cell type (beta vs -log10 FDR)
#   C — Upset / overlap of TOAST vs A1 v1 per-cell-type attribution
#   D — Known markers: hepatocyte (ALB, APOB), macrophage (CD163), stellate
#       (ACTA2) logFC across TOAST cell types

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "RNA-seq/results/celltype_attribution")
toast <- fread(file.path(IN, "toast_celltype_de_all.csv"))
counts <- fread(file.path(IN, "toast_celltype_deg_counts.csv"))
a1 <- fread(file.path(IN, "celltype_primary_attribution.csv"))

toast[, ensg_base := sub("\\.\\d+$", "", gene)]

# Add symbols via bulk dream
bulk_sym <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
                  select = c("gene","symbol"))
bulk_sym[, ensg_base := sub("\\.\\d+$", "", gene)]
bulk_sym <- unique(bulk_sym[, .(ensg_base, symbol)])
toast <- merge(toast, bulk_sym, by = "ensg_base", all.x = TRUE)

# Panel A
counts[, cell_type := factor(cell_type, levels = rev(cell_type[order(-n_sig_padj05)]))]
pA <- ggplot(counts, aes(n_sig_padj05, cell_type, fill = cell_type)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("%d (fdr<0.1: %d)", n_sig_padj05, n_sig_padj10)),
            hjust = -0.05, size = 2.4) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "TOAST CT-specific DEGs (FDR<0.05)", y = NULL,
       title = "A1 v2: TOAST cell-type-specific DE from bulk") +
  theme_masld()

# Panel B — volcano per cell type
toast[, signif := fdr < 0.05]
vol <- toast[!is.na(fdr) & !is.na(effect_size)]
set.seed(42)
if (nrow(vol) > 30000) vol <- vol[sample(.N, 30000)]
pB <- ggplot(vol, aes(effect_size, -log10(fdr + 1e-50), color = signif)) +
  geom_point(size = 0.35, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  facet_wrap(~ cell_type, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TRUE" = "#C0392B","FALSE" = "grey70"),
                     labels = c("TRUE" = "FDR<0.05","FALSE" = "ns"), name = NULL) +
  labs(x = "TOAST effect size (disease vs control)",
       y = "-log10 FDR",
       title = "TOAST CT-specific DE volcano (per cell type)") +
  theme_masld() +
  theme(strip.text = element_text(size = 7))

# Panel C — Overlap TOAST Hep-intrinsic vs A1 v1 primary=Hepatocytes
hep_toast_deg <- toast[cell_type_original == "Hepatocytes" & fdr < 0.05, unique(symbol)]
hep_a1       <- a1[primary_celltype == "Hepatocytes", unique(symbol)]
bulk_deg     <- a1[attribution_class != "NS_bulk", unique(symbol)]

venn_dt <- data.table(
  label = c("TOAST Hep-intrinsic", "A1 v1 Hep-primary", "Bulk DEG universe",
            "TOAST ∩ A1 Hep", "TOAST ∩ Bulk DEG"),
  N = c(length(hep_toast_deg), length(hep_a1), length(bulk_deg),
        length(intersect(hep_toast_deg, hep_a1)),
        length(intersect(hep_toast_deg, bulk_deg)))
)
venn_dt[, label := factor(label, levels = label)]
pC <- ggplot(venn_dt, aes(label, N)) +
  geom_col(fill = "#27AE60", width = 0.55, color = "white") +
  geom_text(aes(label = N), vjust = -0.3, size = 2.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Gene count",
       title = "TOAST ↔ A1 v1 Hepatocyte-attributed overlap") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5))

# Panel D — canonical markers
canon <- c("ALB","APOB","APOA1","CYP3A4",           # Hepatocyte
           "CD163","CD14","MARCO","VSIG4",          # Macrophage
           "CD3D","CD3E","CD8A","IL7R",             # T cells
           "FCN1","S100A9","VCAN")                  # Mono-derived
mk <- toast[symbol %in% canon]
if (nrow(mk) > 0) {
  mk_w <- dcast(mk[, .(symbol, cell_type_original, effect_size)],
                symbol ~ cell_type_original, value.var = "effect_size")
  mk_long <- melt(mk_w, id.vars = "symbol",
                  variable.name = "celltype", value.name = "effect")
  mk_long[, symbol := factor(symbol, levels = canon)]
  pD <- ggplot(mk_long, aes(celltype, symbol, fill = effect)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", effect)), size = 2.2) +
    scale_fill_gradient2(low = "#2980B9", mid = "grey93", high = "#C0392B",
                         midpoint = 0, name = "TOAST\neffect") +
    labs(x = NULL, y = NULL,
         title = "Canonical markers TOAST effect per CT") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5))
} else {
  pD <- ggplot() + labs(title = "No canonical markers found") + theme_masld()
}

fig <- (pA | pC) / pB / pD + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(file.path(FIGS_CELLTYPE_DIR, "figS_A1v2_toast_celltype_de.pdf"),
       fig, width = 15, height = 13)
message("Saved A1 v2 TOAST figure")
