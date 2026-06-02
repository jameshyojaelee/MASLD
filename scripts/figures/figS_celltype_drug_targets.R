#!/usr/bin/env Rscript
# figS_celltype_drug_targets.R
# Figure for Analysis F3 — Cell-type-specific drug target enrichment.
#
# Panels:
#   A — Curated MASLD drug targets × bulk LFC + primary celltype
#   B — Cell-type drug target load bar chart
#   C — OpenTargets drugs by primary celltype (top 20)
#
# Output: figures/supplementary/figS_celltype_biology/figS_F3_celltype_drug_targets.pdf

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "RNA-seq/results/celltype_drug")
cur <- fread(file.path(IN, "curated_masld_targets_celltype.csv"))
ot  <- fread(file.path(IN, "opentargets_drugs_celltype_annotated.csv"))
drug_top <- fread(file.path(IN, "drug_top_celltype.csv"))

cur[, label := paste0(target_symbol, " (", drug_class, ")")]
cur[, label := factor(label, levels = rev(label[order(abs(bulk_lfc_full))]))]
cur[, pc := as.character(primary_celltype)]
cur[is.na(pc) | pc == "", pc := "no CT attribution"]

pA <- ggplot(cur[!is.na(bulk_lfc_full)],
             aes(bulk_lfc_full, label, fill = pc)) +
  geom_col(width = 0.7, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_text(aes(label = sprintf("p=%.2g", bulk_padj_full),
                x = bulk_lfc_full),
            hjust = ifelse(cur[!is.na(bulk_lfc_full)]$bulk_lfc_full > 0, -0.1, 1.1),
            size = 2) +
  scale_fill_brewer(palette = "Set3", name = "Primary\ncelltype") +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(x = "Bulk dream logFC (ashr)", y = NULL,
       title = "Curated MASLD drug targets: bulk DE + primary cell type") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 5.8),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6))

# Panel B: cell-type drug target load
ct_load <- ot[!is.na(primary_celltype), .(n_distinct_drugs = uniqueN(drug_name),
                                           n_target_hits = .N),
              by = primary_celltype][order(-n_target_hits)]
ct_load[, primary_celltype := factor(primary_celltype,
                                      levels = rev(primary_celltype))]
pB <- ggplot(ct_load, aes(n_target_hits, primary_celltype)) +
  geom_col(fill = "#4472C4", width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("%d (%d drugs)", n_target_hits, n_distinct_drugs)),
            hjust = -0.05, size = 2.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Drug-target gene hits", y = NULL,
       title = "Drug-target gene load per cell type (OpenTargets)") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6.5))

# Panel C: top drugs by primary celltype
top_drugs <- drug_top[1:20]
top_drugs[, drug_name := factor(drug_name, levels = rev(drug_name))]
pC <- ggplot(top_drugs, aes(N, drug_name, fill = primary_celltype)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = sprintf("%d (%s)", N, primary_celltype)),
            hjust = -0.05, size = 2) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Target gene hits in primary celltype", y = NULL,
       title = "Top 20 OpenTargets drugs by primary cell type") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 5.8))

fig <- (pA | pB) / pC +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_F3_celltype_drug_targets.pdf")
ggsave(out_path, fig, width = 14, height = 12)
message("Saved: ", out_path)
