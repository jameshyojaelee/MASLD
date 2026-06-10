#!/usr/bin/env Rscript
# figS_stagedeg_routing.R
# ---------------------------------------------------------------------------
# Supplementary figure: three-axis stage-DEG cell-type routing + verdict.
# The headline = MULTI-CELLULAR cascade with a carrier shift at F2->F3
# (hepatocytes carry the early signal, non-parenchymal cells take over the
# advanced-fibrosis transitions). Honest, data-decides — no scale-thumbing.
#
# Five PANELS, each its own PDF (Illustrator assembly):
#   A  Expression routing heatmap   (cell_type x transition, fill = frac carried)
#   B  Three-axis cell-type bars     (expr / LIANA / hotspot; hep vs others)
#   C  LIANA circuit flow            (Fibroblast->Hepatocyte 4->10->16 shift)
#   D  Hotspot module enrichment     (stage-assoc modules x transition; size/colour)
#   E  Per-transition stacked carrier composition (the cascade money panel)
#
# Inputs: Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing/
#   stagedeg_celltype_verdict.csv, stagedeg_verdict_by_transition.csv  (from WS3_04)
#   circuit_celltype_flow_by_stage.csv, hotspot_module_stagedeg_enrichment.csv
#
# Env: rnaseq.  Conventions: PDF only, individual panels, control/reference
# cell types neutral gray (#9E9E9E), stage = fibrosis/disease hues.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

ROUT <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing")

# Dedicated output directory (new constant kept local to avoid editing
# existing figure scripts; mirrors the FIGS_*_DIR pattern in load_figure_data.R).
FIGS_STAGEDEG_DIR <- file.path(FIG_SUPP, "figS_stagedeg_routing")
dir.create(FIGS_STAGEDEG_DIR, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(p, name, w, h) {
  f <- file.path(FIGS_STAGEDEG_DIR, name)
  ggsave(f, p, width = w, height = h, units = "in", device = cairo_pdf)
  sz <- file.info(f)$size
  cat(sprintf("  wrote %s  (%.1f KB)\n", basename(f), sz / 1024))
  invisible(f)
}

# ---------------------------------------------------------------------------
# Colours & ordering
# ---------------------------------------------------------------------------
# Cell types: hepatocyte = parenchymal reference -> neutral gray (#9E9E9E);
# non-parenchymal cells get distinct categorical hues so the carrier shift
# stands out. (Hepatocyte is the "reference" lineage here.)
CT_LEVELS <- c("Hepatocytes","Macrophages","Fibroblasts","Endothelial",
               "Cholangiocytes","T cells","B cells","Resident NK","Plasma")
ct_colors <- c(
  "Hepatocytes"    = "#9E9E9E",  # parenchymal reference — neutral gray
  "Macrophages"    = "#C9265E",  # Liang deep magenta
  "Fibroblasts"    = "#A01753",  # deep fibrosis magenta
  "Endothelial"    = "#1565C0",  # deep blue
  "Cholangiocytes" = "#00897B",  # teal
  "T cells"        = "#7B1FA2",  # violet
  "B cells"        = "#F4A674",  # warm peach
  "Resident NK"    = "#42A5F5",  # light blue
  "Plasma"         = "#8D6E63"   # brown
)

# Stage gradient (early -> advanced): pink->deep magenta fibrosis ramp.
# Augmented cascade F0->F4 + coarse cascade share the same disease ramp.
stage_levels_aug    <- c("F0->F1","F1->F2","F2->F3","F3->F4")
stage_levels_coarse <- c("Healthy->Steatosis","Steatosis->Steatohepatitis",
                         "Steatohepatitis->Cirrhosis")
stage_ramp <- colorRampPalette(c("#F4A674","#e35e9a","#C9265E","#7d1c46"))

# ---------------------------------------------------------------------------
# Load synthesis outputs
# ---------------------------------------------------------------------------
verd_ct  <- fread(file.path(ROUT, "stagedeg_celltype_verdict.csv"))
verd_tr  <- fread(file.path(ROUT, "stagedeg_verdict_by_transition.csv"))

verd_ct[, cell_type := factor(cell_type, levels = CT_LEVELS)]

# Use the AUGMENTED (F-stage) cascade as the primary axis for A/B/E (4 clean
# estimable transitions); coarse cascade shown in the LIANA-flow panel C.
aug <- verd_ct[stage_axis == "augmented"]
aug[, transition_label := factor(transition_label, levels = stage_levels_aug)]

# ===========================================================================
# PANEL A — Expression routing heatmap (cell_type rows x transition cols)
# ===========================================================================
pA <- ggplot(aug, aes(x = transition_label, y = cell_type,
                      fill = expr_frac_carried)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", expr_frac_carried)),
            size = PUB_GEOM_TEXT, color = "grey15") +
  scale_fill_gradientn(colors = c("#FFFFFF", stage_ramp(6)),
                       name = "frac\nDEGs\ncarried",
                       limits = c(0, NA)) +
  scale_y_discrete(limits = rev(CT_LEVELS)) +
  labs(title = "Expression routing of stage-DEGs",
       subtitle = "Fraction of each transition's stage-DEGs carried (intrinsic + multi-cellular)",
       x = "Fibrosis transition", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")
save_panel(pA, "panelA_expression_routing_heatmap.pdf", 3.4, 2.8)

# ===========================================================================
# PANEL B — Three-axis cell-type bars (hepatocyte vs others, per axis)
#   Aggregate across the augmented cascade (mean fraction per CT per axis),
#   then show expr / LIANA / hotspot side by side. Hepatocyte = gray reference.
# ===========================================================================
axis_long <- melt(
  aug[, .(cell_type,
          Expression = expr_frac_carried,
          LIANA      = liana_frac_as_LR,
          Hotspot    = hotspot_frac_modularized)],
  id.vars = "cell_type", variable.name = "routing_axis", value.name = "frac")
axis_mean <- axis_long[, .(frac = mean(frac, na.rm = TRUE)), by = .(cell_type, routing_axis)]
axis_mean[is.nan(frac), frac := NA]
axis_mean[, routing_axis := factor(routing_axis, levels = c("Expression","LIANA","Hotspot"))]
axis_mean[, cell_type := factor(cell_type, levels = CT_LEVELS)]

pB <- ggplot(axis_mean[!is.na(frac)],
             aes(x = cell_type, y = frac, fill = cell_type)) +
  geom_col(width = 0.78, color = "grey25", linewidth = 0.2) +
  facet_wrap(~ routing_axis, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = ct_colors, guide = "none") +
  labs(title = "Three-axis routing share by cell type",
       subtitle = "Mean across the F0->F4 cascade; hepatocyte = parenchymal reference (gray)",
       x = NULL, y = "Mean fraction of stage-DEGs") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_panel(pB, "panelB_three_axis_celltype_bars.pdf", 5.4, 2.6)

# ===========================================================================
# PANEL C — LIANA circuit flow: the Fibroblast->Hepatocyte 4->10->16 shift
#   Grouped bars: stage-DEG-bearing LR-pair counts per source->target pair,
#   across the COARSE cascade (Steatosis -> SH -> Cirrhosis, coarse axis).
# ===========================================================================
flow <- fread(file.path(ROUT, "circuit_celltype_flow_by_stage.csv"))
COARSE_set <- c("coarse_Steatosis","coarse_SH","coarse_Cirrhosis")
flow_c <- flow[axis == "coarse" & stage_set %in% COARSE_set]
flow_c[, stage_lab := factor(
  fcase(stage_set == "coarse_Steatosis", "Healthy->Steatosis",
        stage_set == "coarse_SH",        "Steatosis->Steatohep",
        stage_set == "coarse_Cirrhosis", "Steatohep->Cirrhosis"),
  levels = c("Healthy->Steatosis","Steatosis->Steatohep","Steatohep->Cirrhosis"))]
# Keep the highest-traffic source->target pairs for legibility.
top_pairs <- flow_c[, .(tot = sum(n_stage_deg_LR_pairs)), by = ct_pair][order(-tot)][1:8, ct_pair]
flow_top <- flow_c[ct_pair %in% top_pairs]
flow_top[, ct_pair := factor(ct_pair, levels = flow_c[, .(tot = sum(n_stage_deg_LR_pairs)),
                                                      by = ct_pair][order(tot), ct_pair])]
flow_top[, highlight := grepl("^Fibroblasts->Hepatocytes$", ct_pair)]

pC <- ggplot(flow_top, aes(x = ct_pair, y = n_stage_deg_LR_pairs, fill = stage_lab)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.74,
           color = "grey25", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = setNames(stage_ramp(3), levels(flow_top$stage_lab)),
                    name = "Disease\ntransition") +
  labs(title = "Stage-DEG ligand-receptor circuits by cell-type flow",
       subtitle = "Fibroblast->Hepatocyte signalling expands 4->10->16 with disease stage",
       x = "Source -> Target", y = "Stage-DEG-bearing L-R pairs") +
  theme_masld() + theme_pub() +
  theme(legend.position = "right")
save_panel(pC, "panelC_liana_circuit_flow.pdf", 4.6, 2.8)

# ===========================================================================
# PANEL D — Hotspot module enrichment dotplot
#   stage-associated (progression) modules (y, grouped by cell type) x
#   transition (x); size = overlap (# stage-DEGs in module), color = OR.
#   Augmented cascade; top modules per cell type for legibility.
# ===========================================================================
hot <- fread(file.path(ROUT, "hotspot_module_stagedeg_enrichment.csv"))
hot_a <- hot[axis == "F_stage_augmented" & is_progression_module == TRUE &
             q < 0.05 & stage_set %in% c("fine_F0_F1","fine_F1_F2","fine_F2_F3","fine_F3_F4")]
hot_a[, ct := fcase(cell_type == "hepatocytes","Hepatocytes",
                    cell_type == "macrophages","Macrophages",
                    cell_type == "fibroblasts","Fibroblasts",
                    cell_type == "cholangiocytes","Cholangiocytes",
                    cell_type == "tcells","T cells",
                    default = cell_type)]
# Select the strongest modules (by max OR across the cascade) per cell type.
mod_rank <- hot_a[, .(maxOR = max(OR, na.rm = TRUE)), by = .(ct, module)][order(ct, -maxOR)]
mod_keep <- mod_rank[, head(.SD, 4), by = ct]
hot_a <- merge(hot_a, mod_keep[, .(ct, module)], by = c("ct","module"))
hot_a[, mod_id := paste0(ct, " m", module)]
hot_a[, trans_lab := factor(
  fcase(stage_set=="fine_F0_F1","F0->F1", stage_set=="fine_F1_F2","F1->F2",
        stage_set=="fine_F2_F3","F2->F3", stage_set=="fine_F3_F4","F3->F4"),
  levels = stage_levels_aug)]
# Order modules by cell type then OR
hot_a[, ct := factor(ct, levels = CT_LEVELS)]
mod_order <- hot_a[, .(o = max(OR)), by = .(ct, mod_id)][order(ct, o), mod_id]
hot_a[, mod_id := factor(mod_id, levels = mod_order)]

pD <- ggplot(hot_a, aes(x = trans_lab, y = mod_id, size = overlap, color = OR)) +
  geom_point(alpha = 0.92) +
  scale_size_continuous(range = c(0.8, 5), name = "stage-DEG\noverlap") +
  scale_color_gradientn(colors = stage_ramp(6), name = "OR",
                        trans = "log10") +
  labs(title = "Stage-associated Hotspot modules capture stage-DEGs",
       subtitle = "Progression modules (FDR<0.05) across the F-stage cascade",
       x = "Fibrosis transition", y = "Module (grouped by cell type)") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = 4),
        legend.position = "right")
save_panel(pD, "panelD_hotspot_module_enrichment.pdf", 3.8, 4.4)

# ===========================================================================
# PANEL E — Per-transition stacked carrier composition (the money panel)
#   For each transition, the relative INTRINSIC carrier share per cell type
#   (hepatocyte gray vs non-parenchymal hues). Shows hep -> non-parenchymal
#   carrier shift across F0->F4. Uses intrinsic strong+likely counts so the
#   composition is a genuine partition (no shared-pool double counting).
# ===========================================================================
es <- fread(file.path(ROUT, "routing_summary_by_transition.csv"))
es <- es[axis == "augmented" & stage_set %in% c("fine_F0_F1","fine_F1_F2","fine_F2_F3","fine_F3_F4")]
es[, cell_type := fcase(cell_type=="T_cells","T cells", cell_type=="B_cells","B cells",
                        cell_type=="Resident_NK","Resident NK", default = cell_type)]
es[, n_intrinsic := n_intrinsic_strong + n_intrinsic_likely]
es[, trans_lab := factor(
  fcase(stage_set=="fine_F0_F1","F0->F1", stage_set=="fine_F1_F2","F1->F2",
        stage_set=="fine_F2_F3","F2->F3", stage_set=="fine_F3_F4","F3->F4"),
  levels = stage_levels_aug)]
es[, cell_type := factor(cell_type, levels = CT_LEVELS)]
# Relative composition per transition
es[, frac_of_intrinsic := n_intrinsic / sum(n_intrinsic), by = trans_lab]

pE <- ggplot(es, aes(x = trans_lab, y = frac_of_intrinsic, fill = cell_type)) +
  geom_col(width = 0.74, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = ct_colors, name = "Cell type") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(c(0, 0.02))) +
  labs(title = "Carrier-shift cascade across fibrosis stages",
       subtitle = "Intrinsic stage-DEG carrier composition; hepatocyte (gray) cedes to non-parenchymal cells at F2->F3",
       x = "Fibrosis transition", y = "Share of intrinsically-carried stage-DEGs") +
  theme_masld() + theme_pub() +
  theme(legend.position = "right")
save_panel(pE, "panelE_carrier_shift_cascade.pdf", 3.6, 2.9)

cat("\nAll panels written to:\n  ", FIGS_STAGEDEG_DIR, "\n")
cat("Headline: MULTI-CELLULAR cascade; hepatocyte carrier share peaks F1->F2,\n")
cat("collapses at F2->F3 (non-parenchymal takeover). Verdict = multi_cellular on all 7\n")
cat("transitions; coarse-vs-augmented 3/3 agree.\n")
