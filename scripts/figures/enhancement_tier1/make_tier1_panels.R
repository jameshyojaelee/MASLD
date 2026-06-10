#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Tier-1 enhancement panels (1A relay, 1B isoform switches, 1C master anchor,
# 1E reversal). Individual PDFs (Illustrator assembly), house style.
# Green-lit 2026-06-09. Writes to figures/supplementary/figS_enhancement_tier1/.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))   # theme_masld/theme_pub/save_fig/masld_colors/palette*
OUT  <- file.path(BASE, "figures/supplementary/figS_enhancement_tier1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CTRL <- "#9E9E9E"
style <- function(p) p + theme_masld() + theme_pub()

# ===================== 1A — cellular relay (carrier handoff) =================
ca <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing/carrier_composition_panel.csv"))
tord <- c("F0->F1","F1->F2","F2->F3","F3->F4")
ca <- ca[transition %in% tord]
ca[, transition := factor(transition, levels = tord)]
# highlight the key carriers; others gray
keep_ct <- c("Hepatocytes","T cells","Fibroblasts","Macrophages","Cholangiocytes","Endothelial")
ca[, ct := ifelse(cell_type %in% keep_ct, cell_type, "Other")]
cap <- ca[, .(frac = sum(expr_frac_carried)), by = .(transition, ct)]
pal <- c(Hepatocytes="#C2185B", `T cells`="#1565C0", Fibroblasts="#2E7D32",
         Macrophages="#EF6C00", Cholangiocytes="#6A1B9A", Endothelial="#00838F", Other=CTRL)
p1a <- style(ggplot(cap, aes(transition, frac, color = ct, group = ct)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = "CRN transition", y = "Stage-DEG carrier share",
       title = "Carrier handoff: hepatocyte -> T-cell/fibroblast at F2->F3"))
save_fig(p1a, file.path(OUT, "panel_1A_carrier_relay.pdf"), width = 4.2, height = 3)

# ===================== 1B — stage-specific isoform switches ==================
sw <- fread(file.path(BASE, "RNA-seq/results/isoform_diversity/curated_stage_dtu_switches.csv"))
sw <- sw[symbol != "?" & !is.na(symbol)][order(-max_abs_dprop)][1:20]
sw[, symbol := factor(symbol, levels = rev(symbol))]
consec <- c(`productive_to_unproductive(LoF)`="#C2185B", unproductive_to_productive="#1565C0",
            protein_isoform_switch="#2E7D32", lncRNA_switch="#6A1B9A", other=CTRL)
p1b <- style(ggplot(sw, aes(max_abs_dprop, symbol, color = consequence)) +
  geom_segment(aes(x = 0, xend = max_abs_dprop, yend = symbol), linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = consec, name = "Switch class") +
  labs(x = "Max |dProportion| (stage DTU)", y = NULL,
       title = "Top stage-specific isoform switches (732 confirmed)"))
save_fig(p1b, file.path(OUT, "panel_1B_isoform_switches.pdf"), width = 4.8, height = 4)

# ===================== 1C — credible master regulators ======================
mr <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/master_anchor/master_regulator_credible.csv"))
mr <- mr[order(-breadth, -abs(mean_activity_diff))][1:18]
mr[, symbol := factor(TF, levels = rev(TF))]
mr[, dir := ifelse(consistent_dir == "up", "Up in disease", "Down in disease")]
p1c <- style(ggplot(mr, aes(mean_activity_diff, symbol, fill = dir)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = CTRL, linewidth = 0.3) +
  scale_fill_manual(values = c(`Up in disease`="#C2185B", `Down in disease`="#1565C0"), name = NULL) +
  labs(x = "Mean TF activity change (disease - control)", y = NULL,
       title = "Cross-cell-type disease master-regulators (RUNX1, THRB anchors)"))
save_fig(p1c, file.path(OUT, "panel_1C_master_regulators.pdf"), width = 4.6, height = 4)

# ===================== 1E — reversal: resolution program ====================
hm <- fread(file.path(BASE, "RNA-seq/results/reversal/fgsea_hallmark_regression_specific.csv"))
hm <- hm[padj < 0.05][order(NES)]
top <- rbind(head(hm, 10), tail(hm, 8))
top[, pathway := gsub("HALLMARK_", "", pathway)]
top[, pathway := factor(pathway, levels = pathway)]
top[, dirn := ifelse(NES < 0, "OFF in regression (disease program)", "ON in regression (metabolic identity)")]
p1e <- style(ggplot(top, aes(NES, pathway, fill = dirn)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, color = CTRL, linewidth = 0.3) +
  scale_fill_manual(values = c("OFF in regression (disease program)"="#C2185B",
                               "ON in regression (metabolic identity)"="#2E7D32"), name = NULL) +
  labs(x = "Regression-specific NES (Hallmark)", y = NULL,
       title = "Fibrosis regression reverses the disease program"))
save_fig(p1e, file.path(OUT, "panel_1E_reversal_hallmark.pdf"), width = 5.2, height = 4)

# 1E rewind scatter: regression-specific vs cross-sectional disease axis
rs <- fread(file.path(BASE, "RNA-seq/results/reversal/de_regression_specific.csv"))
cn <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
m  <- merge(rs[, .(gene, reg = logFC)], cn[, .(gene, dis = logFC)], by = "gene")
m  <- m[is.finite(reg) & is.finite(dis)]
rho <- cor(m$reg, m$dis, method = "spearman")
p1e2 <- style(ggplot(m, aes(dis, reg)) +
  geom_point(size = 0.25, alpha = 0.25, color = CTRL) +
  geom_smooth(method = "lm", se = FALSE, color = "#C2185B", linewidth = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = CTRL) +
  geom_vline(xintercept = 0, linewidth = 0.2, color = CTRL) +
  labs(x = "Disease axis logFC (canonical, disease vs control)",
       y = "Regression-specific logFC",
       title = sprintf("Rewind: regression vs disease axis (rho = %.2f)", rho)))
save_fig(p1e2, file.path(OUT, "panel_1E_rewind_scatter.pdf"), width = 3.6, height = 3.4)

cat("Wrote Tier-1 panels to", OUT, "\n"); print(list.files(OUT))
