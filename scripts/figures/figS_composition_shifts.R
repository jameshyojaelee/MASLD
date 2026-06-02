#!/usr/bin/env Rscript
# figS_composition_shifts.R
# Figure for Analysis A2 — Cell-type composition shifts in MASLD.
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Panels:
#   A — Forest: MASLD vs Healthy logit shift per cell type, ranked by |t|
#   B — Boxplots of 4 key proportions across disease stages
#   C — F2 inflection: F_low vs F_high composition shift forest
#   D — Alignment with A1 attribution count per cell type
#
# Output: figures/supplementary/figS_celltype_biology/figS_A2_composition_shifts.pdf
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

ATT_DIR <- file.path(BASE, "RNA-seq/results/celltype_attribution")
shifts  <- fread(file.path(ATT_DIR, "composition_shifts.csv"))
props   <- fread(file.path(ATT_DIR, "persample_celltype_proportions.csv"))
attrib  <- fread(file.path(ATT_DIR, "celltype_primary_attribution.csv"))

# ---- Panel A: MASLD vs Healthy forest ---------------------------------------
rep <- shifts[contrast == "Control_vs_Disease"][order(-abs(t))]
rep[, celltype := factor(celltype, levels = rev(celltype))]
rep[, se := abs(logit_diff / t)]
rep[, ci_lo := logit_diff - 1.96 * se]
rep[, ci_hi := logit_diff + 1.96 * se]
rep[, sig := fifelse(padj < 0.05, "sig", "ns")]

pA <- ggplot(rep, aes(logit_diff, celltype, color = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, linewidth = 0.4) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("sig" = "#C0392B", "ns" = "grey60"),
                     labels = c("sig" = "padj<0.05", "ns" = "ns"),
                     name = NULL) +
  labs(x = "Logit shift (Disease - Control)", y = NULL,
       title = "Cell-type composition shift: MASLD vs Healthy",
       subtitle = sprintf("n = %d samples, %d cohorts; propeller-equivalent",
                          rep$n[1], length(unique(props$dataset)))) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6.5))

# ---- Panel B: Proportion boxplots across F-stages ---------------------------
key_cts <- c("Hepatocytes", "Macrophages", "Endothelial cells", "T cells")
key_cts <- intersect(key_cts, names(props))
if (length(key_cts) > 0 && "fibrosis_stage" %in% names(props)) {
  bx <- melt(props[!is.na(fibrosis_stage), c("sample_id","fibrosis_stage", key_cts), with = FALSE],
             id.vars = c("sample_id","fibrosis_stage"),
             variable.name = "celltype", value.name = "proportion")
  bx[, fibrosis_stage := as.character(fibrosis_stage)]
  bx[, fibrosis_stage := sub("^F", "", fibrosis_stage)]
  bx <- bx[fibrosis_stage %in% c("0","1","2","3","4")]
  bx[, fibrosis_stage := factor(fibrosis_stage, levels = c("0","1","2","3","4"))]

  pB <- ggplot(bx, aes(fibrosis_stage, proportion, fill = celltype)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.75) +
    facet_wrap(~ celltype, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = c("Hepatocytes" = "#27AE60",
                                 "Macrophages" = "#C0392B",
                                 "Endothelial cells" = "#4472C4",
                                 "T cells" = "#F39C12"),
                      guide = "none") +
    labs(x = "Fibrosis stage", y = "MuSiC-estimated proportion",
         title = "Composition shift across fibrosis stages") +
    theme_masld() +
    theme(strip.text = element_text(size = 6.5))
} else {
  pB <- ggplot() + labs(title = "No fibrosis_stage available") + theme_masld()
}

# ---- Panel C: F2 inflection -------------------------------------------------
f2 <- shifts[contrast == "F_low_vs_F_high"]
if (nrow(f2) > 0) {
  f2 <- f2[order(-abs(t))]
  f2[, celltype := factor(celltype, levels = rev(celltype))]
  f2[, se := abs(logit_diff / t)]
  f2[, ci_lo := logit_diff - 1.96 * se]
  f2[, ci_hi := logit_diff + 1.96 * se]
  f2[, sig := fifelse(padj < 0.05, "sig", "ns")]

  pC <- ggplot(f2, aes(logit_diff, celltype, color = sig)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, linewidth = 0.4) +
    geom_point(size = 1.8) +
    scale_color_manual(values = c("sig" = "#C0392B", "ns" = "grey60"),
                       labels = c("sig" = "padj<0.05", "ns" = "ns"),
                       name = NULL) +
    labs(x = "Logit shift (F>=2 - F<=1)", y = NULL,
         title = "F2 inflection: composition shift") +
    theme_masld() +
    theme(axis.text.y = element_text(size = 6.5))
} else {
  pC <- ggplot() + labs(title = "No F-stage data") + theme_masld()
}

# ---- Panel D: A1 attribution count vs A2 composition t-stat -----------------
attr_counts <- attrib[!is.na(primary_celltype),
                      .(n_attributed = .N), by = primary_celltype]
# harmonize name mapping (A1 uses underscore; A2 uses space)
attr_counts[, ct_key := gsub("_", " ", primary_celltype)]
attr_counts[ct_key == "Mono+mono derived cells", ct_key := "Mono+mono derived cells"]

merge_dt <- merge(rep[, .(celltype = as.character(celltype), t_shift = t, padj)],
                  attr_counts[, .(celltype = ct_key, n_attributed)],
                  by = "celltype", all = FALSE)

if (nrow(merge_dt) > 0) {
  pD <- ggplot(merge_dt, aes(t_shift, n_attributed, label = celltype)) +
    geom_point(size = 2, color = "#4472C4") +
    geom_text(size = 2, hjust = -0.1, vjust = -0.3) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
    labs(x = "A2 composition shift t-stat (MASLD vs Healthy)",
         y = "A1 DEGs attributed to cell type",
         title = "Composition shift vs. intrinsic DEG load") +
    theme_masld()
} else {
  pD <- ggplot() + labs(title = "No overlap between A1 + A2 CT names") + theme_masld()
}

# ---- Assemble ---------------------------------------------------------------
fig <- (pA + pB) / (pC + pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_A2_composition_shifts.pdf")
ggsave(out_path, fig, width = 13, height = 10)
message("Saved: ", out_path)
