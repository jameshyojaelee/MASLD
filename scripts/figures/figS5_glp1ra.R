#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# figS5 — Incretin/glucagon drug-axis supplement to main Fig 5E.
# CONSOLIDATES the five scattered figS_glp1ra panels into three coherent,
# publication-quality composite figures (patchwork; house palette; matches the
# Fig 2/3/4 craft bar). SURVIVORS of the WS5 adversarial gate only
# (RNA-seq/results/glp1ra/WS5_adversarial_verdicts.txt).
#
#   figS5_receptor_localization.pdf  a receptor %-expressing replicates across 7 datasets
#                                    b GLP1R stays at floor healthy AND disease (not a loss);
#                                      GCGR detection-fraction drops with disease
#                                    c GLP1R below the ambient floor even in resolved
#                                      pericentral LSECs (sensitivity-limited, not absent)
#   figS5_axis_dynamics.pdf          a stage-resolved DE (GCGR down @steatosis; DPP4 up @NASH)
#                                    b DPP4 pericentral zonation (Visium)
#   figS5_axis_genetics.pdf            axis is a powered genetic null (0/6 colocalize)
#
# Style: 6pt Helvetica plain, text black, gene symbols italic, NO lollipop, PDF.
# Palette: masld up=#C9265E / down=#1565C0 / control=#9E9E9E; pericentral=#C2185B.
# Output: figures/supplementary/figS_glp1ra/ . Env: rnaseq.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
OUT <- FIGS_GLP1RA_DIR; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
G   <- file.path(BASE, "RNA-seq/results/glp1ra")
GEOM_TXT <- 6 / ggplot2::.pt
UP <- "#C9265E"; DOWN <- "#1565C0"; CTRL <- "#9E9E9E"; PERI <- "#C2185B"
tag_theme <- theme(plot.tag = element_text(size = 6, face = "plain"))

# =============================== FIGURE 1 =====================================
# a — cross-dataset replication of receptor %-expressing in hepatocytes
rep <- fread(file.path(G, "replication/receptor_by_dataset.csv"))
rep[, gene := factor(gene, levels = c("GLP1R", "GIPR", "GCGR", "DPP4"))]
rep <- rep[!is.na(gene)]
rep[, y := pmax(pct_expressing, 0.003)]
per  <- rep[dataset != "POOLED_ALL"]
pool <- rep[dataset == "POOLED_ALL"]
pA <- ggplot(per, aes(gene, y)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey80", linewidth = 0.25) +
  geom_point(position = position_jitter(width = 0.16, height = 0), size = 0.85,
             colour = "grey55", alpha = 0.9) +
  geom_point(data = pool, shape = 23, size = 1.9, fill = UP, colour = "white", stroke = 0.25) +
  annotate("text", x = 0.62, y = 32, hjust = 0, size = GEOM_TXT, colour = "grey40",
           label = "● dataset") +
  annotate("point", x = 0.66, y = 12, shape = 23, size = 1.6, fill = UP, colour = "white", stroke = 0.25) +
  annotate("text", x = 0.75, y = 12, hjust = 0, size = GEOM_TXT, colour = "grey40", label = "pooled") +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 50), labels = c("0.01", "0.1", "1", "10", "50")) +
  labs(x = NULL, y = "% hepatocytes expressing") +
  theme_masld() +
  theme(axis.text.x = element_text(face = "italic"),
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.2))

# b — Healthy vs MASLD hepatocyte %-expressing (per-cell fraction)
cond <- fread(file.path(G, "scrna_incretin_axis/glp1r_pct_by_condition.csv"))
cond[, gene := factor(gene, levels = c("GLP1R", "GCGR", "DPP4"))]
cond[, condition := factor(condition, levels = c("Healthy", "MASLD"))]
pB <- ggplot(cond, aes(condition, pct_cells_expressing, fill = condition)) +
  geom_col(width = 0.68) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(Healthy = CTRL, MASLD = UP), name = NULL) +
  labs(x = NULL, y = "% hepatocytes expressing") +
  theme_masld() +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "bottom", legend.margin = margin(t = -4),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.spacing = unit(0.35, "lines"))

# c — pericentral-LSEC receptor ladder (GLP1R below the ambient floor)
rc <- fread(file.path(G, "lsec/endothelial_subcluster_receptors.csv"))
ladder <- c("STAB2", "ALB", "PECAM1", "VWF", "DPP4", "GCGR", "GIPR", "GLP2R", "GLP1R")
agg <- rc[lsec_subtype == "pericentral_LSEC" & gene %in% ladder,
          .(pct = 100 * sum(n_pos) / sum(n_cells)), by = gene]
agg[, gene := factor(gene, levels = rev(ladder))]
agg[, y := pmax(pct, 0.02)]
agg[, kind := fifelse(gene == "ALB", "ambient (hepatocyte)",
             fifelse(gene %in% c("STAB2", "PECAM1", "VWF"), "endothelial marker",
             "incretin/glucagon receptor"))]
pC <- ggplot(agg, aes(y, gene, colour = kind)) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c("ambient (hepatocyte)" = CTRL,
                                 "endothelial marker" = DOWN,
                                 "incretin/glucagon receptor" = UP), name = NULL) +
  scale_x_log10(breaks = c(0.03, 0.3, 3, 30), labels = c("0.03", "0.3", "3", "30")) +
  labs(x = "% pericentral LSECs expressing", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.2),
        legend.position = "bottom", legend.direction = "vertical",
        legend.margin = margin(t = -4), legend.key.height = unit(0.22, "cm"))

f1 <- (pA | pB | pC) + plot_layout(widths = c(1, 1.15, 1.25)) +
  plot_annotation(tag_levels = "a") & tag_theme
ggsave(file.path(OUT, "figS5_receptor_localization.pdf"), f1,
       width = 7.2, height = 2.5, device = pdf_device)

# =============================== FIGURE 2 =====================================
# a — DONOR-COLLAPSED stage DE (honest null). The run-level heatmap was RETIRED
#     2026-07-11 after the GSE244832 pseudoreplication fix (117 runs -> 18 donors):
#     neither GCGR nor DPP4 survives at the biological-donor level (all 95% CI cross 0).
dc <- fread(file.path(G, "axis_progression/incretin_axis_stage_de_donorcollapsed.csv"))
dc <- dc[subset == "combined_primary" & tested == TRUE & gene %in% c("GCGR", "DPP4", "GIPR", "GLP2R")]
dc[, transition := factor(transition, levels = c("Steatosis_vs_Healthy", "Steatohepatitis_vs_Steatosis"),
                          labels = c("Healthy → Steatosis", "Steatosis → SH"))]
dc[, gene := factor(gene, levels = rev(c("GCGR", "DPP4", "GIPR", "GLP2R")))]
pDA <- ggplot(dc, aes(logFC, gene)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.22, colour = "#4a4a4a", linewidth = 0.4) +
  geom_point(size = 1.3, colour = "#4a4a4a") +
  facet_wrap(~transition, nrow = 1) +
  labs(x = "donor-level log2FC (95% CI)", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(face = "italic"), panel.spacing = unit(0.4, "lines"))

# b — DPP4 pericentral zonation (Visium; descriptive)
z <- fread(file.path(G, "zonation/receptor_zonation.csv")); setnames(z, 1, "zone")
zl <- melt(z[, .(zone, GCGR, DPP4)], id.vars = "zone", variable.name = "gene", value.name = "expr")
zl[, zone := factor(zone, levels = c("PP1", "PP2", "Mid", "PC2", "PC1"))]
pDB <- ggplot(zl, aes(zone, expr, group = gene, colour = gene)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  scale_colour_manual(values = c(GCGR = CTRL, DPP4 = PERI), name = NULL) +
  labs(x = "periportal  →  pericentral", y = "mean expression (Visium)") +
  theme_masld() +
  theme(legend.position = c(0.84, 0.55), legend.text = element_text(face = "italic"),
        legend.key.height = unit(0.28, "cm"))

f2 <- (pDA | pDB) + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "a") & tag_theme
ggsave(file.path(OUT, "figS5_axis_dynamics.pdf"), f2,
       width = 5.6, height = 2.15, device = pdf_device)

# =============================== FIGURE 3 (ambient) ===========================
# decontX ambient correction (657,629 hep, 6 datasets): GCGR/DPP4 hepatocyte signal
# is cell-INTRINSIC (ambient_frac ~0) while cross-lineage counts (PECAM1/STAB2 in hep,
# ALB & GCGR in endo) are ambient -> bulletproofs GCGR hepatocyte-enrichment. GLP1R
# shown for context (decontX addresses false-positives, not dropout).
amb <- fread(file.path(G, "scrna_incretin_axis/decontx_ambient_correction.csv"))
amb[, ct := fifelse(cell_type == "Hepatocytes", "hep", "endo")]
amb[, lab := paste0(gene, " (", ct, ")")]
# GLP1R deliberately EXCLUDED: at 0.04% its ambient_frac is uninformative (decontX
# addresses false positives, not dropout) and a "cell-intrinsic" bar would mislead.
sel <- amb[(gene == "GCGR") | (gene == "DPP4" & ct == "hep") |
           (gene == "ALB" & ct == "endo") |
           (gene %in% c("PECAM1", "STAB2") & ct == "hep")]
sel[, kind := fifelse(ambient_fraction < 0.15, "cell-intrinsic", "ambient (cross-lineage)")]
sel[, lab := factor(lab, levels = sel[order(ambient_fraction)]$lab)]
pAmb <- ggplot(sel, aes(ambient_fraction, lab, fill = kind)) +
  geom_col(width = 0.66) +
  scale_fill_manual(values = c("cell-intrinsic" = UP, "ambient (cross-lineage)" = CTRL), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "decontX ambient fraction", y = NULL) +
  theme_masld() +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom",
        legend.direction = "vertical", legend.margin = margin(t = -4),
        legend.key.height = unit(0.22, "cm"))
ggsave(file.path(OUT, "figS5_ambient_correction.pdf"), pAmb,
       width = 3.0, height = 2.3, device = pdf_device)

# =============================== FIGURE 3 (RETIRED 2026-07-11) =================
# The "powered genetic null" figure was RETRACTED: it tested LIVER eQTL (the wrong
# tissue for these off-parenchyma receptors) and used OTTERS n_snps (prediction-
# model size) as a false "power" proxy. Genetics is now handled as a cited
# paragraph (published drug-target MR — Yan 2024 GIPR-protective) plus a
# tissue-contrast panel built separately. See the CLAUDE.md rule
# "Novelty & source verification" + feedback-verify-literature-and-source-before-claiming.

message("figS5 written to: ", OUT)
message("  figS5_receptor_localization.pdf  figS5_axis_dynamics.pdf  figS5_ambient_correction.pdf")
message("CAPTION figS5_receptor_localization: (a) Receptor %-expressing in hepatocytes replicates ",
        "across all 7 scRNA datasets (grey) and pooled (diamond): GLP1R/GIPR stay <1%, GCGR/DPP4 ",
        "detectable. (b) GLP1R is at the floor in BOTH healthy and MASLD hepatocytes (not a disease ",
        "loss), while GCGR-expressing fraction falls 26.6%->15.6%. (c) Even in 83,388 scVI-resolved ",
        "pericentral LSECs, GLP1R sits three logs below hepatocyte-ambient ALB and below detectable ",
        "receptors -> localization is sensitivity-limited, not demonstrably absent.")
message("CAPTION figS5_axis_dynamics: (a) DONOR-LEVEL stage DE (log2FC + 95% CI; GSE244832 runs collapsed ",
        "to biological donors) — neither GCGR nor DPP4 reaches significance (all CI cross 0; the run-level ",
        "significance was pseudoreplication). GIPR/GLP2R nominal signals are single-dataset (GSE174748, ",
        "2-vs-2 donors) and not reproducible. (b) Descriptive DPP4 pericentral Visium gradient (~3-fold); ",
        "GCGR modestly periportal; single cohort, no inferential test.")
message("CAPTION figS5_ambient_correction: decontX ambient fraction (657,629 hepatocytes, 6 datasets). GCGR ",
        "and DPP4 hepatocyte signal is cell-intrinsic (ambient fraction ~0) while cross-lineage counts ",
        "(PECAM1/STAB2 in hepatocytes; ALB and GCGR in endothelial) are ambient -> the GCGR hepatocyte ",
        "enrichment is not contamination. (GLP1R omitted: decontX addresses ambient/false positives, not ",
        "dropout, so it is uninformative for GLP1R's near-zero signal.)")
