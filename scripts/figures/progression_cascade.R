#!/usr/bin/env Rscript
# KEY MESSAGE: Disease dysregulation expands MONOTONICALLY across fibrosis stages on three orthogonal readouts — stage DEG counts grow, the NMF program mix shifts from quiescent to inflammatory/fibrogenic, and the cellular compartment tips hepatocyte->non-parenchymal — visualizing the multi-step cellular cascade thesis on one shared F0->F4 axis.
# ============================================================================
# progression_cascade.R  — Figure 3F (unified disease-progression cascade)
#
# Three vertically-stacked tracks sharing one discrete F0..F4 x-axis:
#   TRACK 1 — stage-specific DEG counts (up vs down) per stage vs F0
#   TRACK 2 — NMF 6-program dominant-sample composition shift across stages
#   TRACK 3 — cell-type composition shift (hepatocyte decline, NPC rise)
#
# Output: FIG2_DIR/panels/fig3e_progression_cascade.pdf
#   (FIG2_DIR resolves to figures/main/fig3_RNAseq — back-compat constant name.)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)

STAGES <- c("F0", "F1", "F2", "F3", "F4")

# Shared x-axis geometry so F0..F4 columns align vertically across ALL tracks.
# Identical discrete expansion on every panel keeps the data x-positions equal;
# an identical right-side plot margin reserves a label gutter (for p3b's NPC
# end-labels) WITHOUT distorting the x scale, while keeping panel widths equal.
X_EXPAND      <- ggplot2::expansion(mult = c(0.04, 0.04))
RIGHT_GUTTER  <- 64                      # pt of reserved right margin (label gutter)
SHARED_MARGIN <- ggplot2::margin(3, RIGHT_GUTTER, 3, 3)

# ----------------------------------------------------------------------------
# TRACK 1 — stage-specific DEG counts (up vs down) vs F0
# ----------------------------------------------------------------------------
deg <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv"))
# Effect size: prefer ashr-shrunk; fall back to raw logFC if shrunk missing.
deg[, eff := if ("shrunk_logFC" %in% names(deg)) shrunk_logFC else logFC]
deg[is.na(eff), eff := logFC]
deg[, stage := sub("_vs_F0$", "", contrast)]              # F1_vs_F0 -> F1
deg_sig <- deg[padj < 0.05 & abs(eff) > 0.5]
deg_counts <- deg_sig[, .(
    up   = sum(eff > 0, na.rm = TRUE),
    down = sum(eff < 0, na.rm = TRUE)
  ), by = stage]
# F0 vs F0 = baseline, zero by definition (no contrast computed)
deg_counts <- rbind(
  data.table(stage = "F0", up = 0L, down = 0L),
  deg_counts, fill = TRUE)
deg_counts[, stage := factor(stage, levels = STAGES)]
setorder(deg_counts, stage)

deg_long <- melt(deg_counts, id.vars = "stage",
                 measure.vars = c("up", "down"),
                 variable.name = "direction", value.name = "n")
# Down plotted as negative for a back-to-back (diverging) read.
deg_long[, n_signed := ifelse(direction == "down", -n, n)]
deg_long[, direction := factor(direction, levels = c("up", "down"),
                               labels = c("Up vs F0", "Down vs F0"))]

DEG_COLORS <- c("Up vs F0"   = masld_colors$mash,   # #C9265E magenta
                "Down vs F0" = masld_colors$down)   # #1565C0 blue

p1 <- ggplot(deg_long, aes(stage, n_signed, fill = direction)) +
  geom_col(width = 0.66, color = NA) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  geom_text(data = deg_long[n > 0],
            aes(label = n,
                vjust = ifelse(direction == "Up vs F0", -0.35, 1.25)),
            size = 2.4, color = "gray20") +
  scale_fill_manual(values = DEG_COLORS, name = NULL) +
  scale_y_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0.12, 0.18))) +
  scale_x_discrete(expand = X_EXPAND) +
  labs(x = NULL, y = "Stage DEGs\n(vs F0)", fill = "Stage DEG") +
  theme_masld(base_size = 9) +
  theme(legend.key.size = unit(0.28, "cm"),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# TRACK 2 — NMF program composition: fraction of samples dominated by each P1..P6
# ----------------------------------------------------------------------------
nmf  <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
labs <- fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv"))
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  select = c("sample_id", "fibrosis_stage"))

# program_labels uses program_code/biological_label (not program/label).
lab_map <- setNames(labs$biological_label, labs$program_code)

nmf <- merge(nmf[, .(sample_id, dominant_program_code)],
             meta, by = "sample_id")
nmf <- nmf[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
nmf[, stage := factor(paste0("F", fibrosis_stage), levels = STAGES)]
nmf[, prog_label := lab_map[dominant_program_code]]

nmf_comp <- nmf[, .N, by = .(stage, dominant_program_code, prog_label)]
nmf_comp[, frac := N / sum(N), by = stage]
# Order programs from disease/active (top of stack) to quiescent/stable (bottom).
prog_order <- c("P1", "P3", "P2", "P4", "P5", "P6")
prog_labels_ordered <- lab_map[prog_order]
nmf_comp[, prog_label := factor(prog_label, levels = prog_labels_ordered)]

# Colors: active programs get saturated cat hues; the two quiescent-parenchyma
# baselines and the two stable programs anchor in/near control gray.
NMF_COLORS <- c(
  "Progression-Inflammatory" = cat_palette[3],   # red-coral — active disease
  "Fibrogenic"               = cat_palette[6],    # violet     — active fibrosis
  "Quiescent-Parenchyma_1"   = "#BDBDBD",         # near-gray (baseline parenchyma)
  "Quiescent-Parenchyma_2"   = "#D6D6D6",
  "Stable_1"                 = "#9E9E9E",         # control gray
  "Stable_2"                 = "#7D7D7D")

p2 <- ggplot(nmf_comp, aes(stage, frac, fill = prog_label)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = NMF_COLORS, name = "NMF program",
                    breaks = prog_labels_ordered,
                    labels = prog_labels_ordered) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(expand = X_EXPAND) +
  labs(x = NULL, y = "NMF program\n(% samples)") +
  guides(fill = guide_legend(ncol = 1, keyheight = unit(0.26, "cm"))) +
  theme_masld(base_size = 9) +
  theme(legend.position = "right",
        legend.text  = element_text(size = 6),
        legend.title = element_text(size = 6.5, face = "bold"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        plot.margin  = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# TRACK 3 — cell-type composition shift across stages
# ----------------------------------------------------------------------------
ct <- fread(file.path(BASE,
  "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"))
ct <- ct[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
ct[, stage := factor(paste0("F", fibrosis_stage), levels = STAGES)]

CT_KEEP <- c("Hepatocytes", "Fibroblasts", "Macrophages",
             "T cells", "Endothelial cells")
ct_mean <- ct[, lapply(.SD, mean, na.rm = TRUE), by = stage,
              .SDcols = CT_KEEP]
ct_long <- melt(ct_mean, id.vars = "stage",
                variable.name = "celltype", value.name = "frac")
# Sort cell types by direction of change F0->F4: hepatocyte (declining) first,
# then the rising non-parenchymal compartments by magnitude of rise.
delta <- ct_long[, .(d = frac[stage == "F4"] - frac[stage == "F0"]), by = celltype]
ct_order <- c("Hepatocytes",
              as.character(delta[celltype != "Hepatocytes"][order(-d)]$celltype))
ct_long[, celltype := factor(celltype, levels = ct_order)]

# Hepatocyte baseline = control gray (it is the compartment being displaced);
# rising NPC compartments take semantic ct_palette hues.
CT_COLORS <- c(
  "Hepatocytes"       = "#9E9E9E",        # control gray (displaced baseline)
  "Fibroblasts"       = ct_palette[["Fibroblasts"]],
  "Macrophages"       = ct_palette[["Macrophages"]],
  "T cells"           = ct_palette[["T cells"]],
  "Endothelial cells" = ct_palette[["Endothelial cells"]])

# Hepatocyte (78-93%) and the NPC compartments (0-4%) live on incompatible
# scales; a single linear axis crushes the NPC rise. Split Track 3 into two
# sub-panels that still share the F0..F4 x-axis: (3a) hepatocyte decline,
# (3b) the rising non-parenchymal compartments (magnified, NOT % of total ->
# % within compartment is misleading; keep absolute fraction with its own axis).
ct_hep <- ct_long[celltype == "Hepatocytes"]
ct_npc <- ct_long[celltype != "Hepatocytes"]
ct_npc[, celltype := droplevels(celltype)]
lab_npc <- ct_npc[stage == "F4"]

p3a <- ggplot(ct_hep, aes(stage, frac, group = 1)) +
  geom_line(linewidth = 0.8, color = "#9E9E9E") +
  geom_point(size = 1.4, color = "#9E9E9E") +
  geom_text(data = ct_hep[stage %in% c("F0", "F4")],
            aes(label = scales::percent(frac, accuracy = 1)),
            vjust = -0.9, size = 2.2, color = "gray25") +
  annotate("text", x = 3, y = max(ct_hep$frac),
           label = "Hepatocytes", size = 2.4, color = "gray25", vjust = -1.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(min(ct_hep$frac) * 0.97, max(ct_hep$frac) * 1.06)) +
  scale_x_discrete(expand = X_EXPAND) +
  labs(x = NULL, y = "Hepatocyte\nfraction") +
  theme_masld(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = SHARED_MARGIN)

# NPC end-labels render into the reserved right margin (clip="off") instead of
# distorting the x scale — so p3b shares the SAME X_EXPAND as every other track.
p3b <- ggplot(ct_npc, aes(stage, frac, color = celltype, group = celltype)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  ggrepel::geom_text_repel(
    data = lab_npc, aes(label = celltype),
    hjust = 0, direction = "y", nudge_x = 0.6,
    segment.size = 0.2, segment.color = "gray70",
    size = 2.1, min.segment.length = 0, box.padding = 0.15,
    xlim = c(5.35, 6.6), show.legend = FALSE) +
  scale_color_manual(values = CT_COLORS, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.10))) +
  scale_x_discrete(expand = X_EXPAND) +
  coord_cartesian(clip = "off") +
  labs(x = "Fibrosis stage", y = "Non-parenchymal\nfraction") +
  theme_masld(base_size = 9) +
  theme(axis.text.x = element_text(size = 9),
        plot.margin = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# Assemble — vertical stack, shared discrete x; only bottom track shows ticks.
# ----------------------------------------------------------------------------
cascade <- p1 / p2 / p3a / p3b +
  plot_layout(heights = c(1.05, 1.05, 0.7, 0.9), guides = "collect") +
  plot_annotation(
    title = "Disease program expands monotonically with fibrosis stage",
    theme = theme(plot.title = element_text(size = 10, face = "bold")))

save_fig(cascade, file.path(PANEL_DIR, "fig3e_progression_cascade.pdf"),
         width = fig_full_width * 0.92, height = 7.0)

# ----------------------------------------------------------------------------
# Underlying summary tables (caption transparency)
# ----------------------------------------------------------------------------
fwrite(deg_counts, file.path(DATA_DIR, "progression_cascade_deg_counts.csv"))
fwrite(nmf_comp[order(stage, prog_label)],
       file.path(DATA_DIR, "progression_cascade_nmf_composition.csv"))
fwrite(ct_mean, file.path(DATA_DIR, "progression_cascade_celltype_means.csv"))

cat("Wrote fig3e_progression_cascade.pdf\n")
cat("DEG counts per stage (up/down):\n"); print(deg_counts)
cat(sprintf("Hepatocyte fraction F0=%.3f  F4=%.3f  (%.1f%% decline)\n",
            ct_mean[stage == "F0", Hepatocytes],
            ct_mean[stage == "F4", Hepatocytes],
            100 * (1 - ct_mean[stage == "F4", Hepatocytes] /
                       ct_mean[stage == "F0", Hepatocytes])))
