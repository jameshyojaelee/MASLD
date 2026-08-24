#!/usr/bin/env Rscript
# Moved from Fig 4 validation to Fig 5 convergence on 2026-07-08.
# This script name and emitted panel now use the same Fig 5 letter to avoid
# stale Fig 4 rerenders landing in the wrong main-figure directory.
# ─────────────────────────────────────────────────────────────────────────────
# Fig 5c — drug-target calibration keystone, a TWO-EVIDENCE-AXIS scatter. Drug-development
# stage is the DOT COLOUR (not an axis); the two orthogonal evidence channels are
# the axes:
#     x = genetic colocalization (COLOC PP.H4)
#     y = bulk transcriptomic effect (raw canonical logFC, disease vs control)
# This places every famous MASLD drug + every method nomination in the joint
# genetic×transcriptomic evidence plane; colour reads off where each clinical
# fate (Approved / In trials / Failed / Preclinical / Discovery) lands.
#
# Reads the joint plane as gates: x = 0.5 (COLOC significance), y = 0 and ±0.25
# (TREAT interval-null bound, lfc=0.25). The FDA-approved cis-target THRB and the de-novo Phase-1
# discovery RORA sit far-right and below y=0 (genetically colocalizing, suppressed
# in disease); the failed NASH cluster (NR1H4 / MAP3K5 / CCR2 / LOXL2) sits LEFT
# of the COLOC gate (no genetic support) even when DEG-up; and the approved-but-
# non-cis drugs (GLP1R, SLC5A2 / PNPLA3, HSD17B13 ASOs) sit far-left — COLOC is
# SPECIFIC, not high-recall, the honest calibrated behaviour.
#
# Gene set = curated MASLD drug roster (approved/trial/failed) + method nominations
# clearing the paper-wide dual-evidence bar (COLOC>0.5 AND a TREAT bulk DEG, FDR<0.05; see the
# FAIRNESS block below). NO lollipop (memory/feedback-no-lollipop); all text BLACK,
# gene names italic, no titles/subtitles/in-plot annotations (memory/feedback-no-colored-fonts).
#
# Output: figures/main/fig6_gene_catalog/panels/fig6_drug_target_calibration.pdf
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

COLOC_MIN <- 0.5   # reservoir inclusion floor ("significant" colocalization; gate threshold)
FIG5_DATA_DIR <- file.path(FIG5_DIR, "data")
dir.create(file.path(FIG5_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(FIG5_DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Genome-wide drug-development status ──────────────────────────────────────
cls <- fread(file.path(BASE, "data/external/drug_targets/drug_target_classification.tsv"))
setnames(cls, "symbol", "gene")
keep_status <- c("masld_approved", "masld_clinical", "masld_discontinued",
                 "masld_preclinical", "discovery")
cls <- cls[drug_dev_status %in% keep_status & gene != "" & !is.na(gene)]
# NOTE: the 6 manual MASLD fact-check corrections (NR3C2/VDR -> preclinical,
# F2RL1 -> clinical, FADS2/GNMT/KRT8 -> discovery; 2026-06-23, ClinicalTrials.gov +
# PubMed verified) now live AT SOURCE in build_drug_target_classification.R
# (flagged manual_factcheck_override), so drug_dev_status arrives already corrected.

# outcome class (the colour/"stage" channel) from status + best MASLD phase
cls[, outcome := NA_character_]
cls[drug_dev_status == "discovery",         outcome := "Discovery"]
cls[drug_dev_status == "masld_preclinical", outcome := "Preclinical"]
cls[drug_dev_status == "masld_clinical",    outcome := "In trials"]
cls[drug_dev_status == "masld_approved",    outcome := "Approved"]
cls[drug_dev_status == "masld_discontinued", outcome := "Failed / discontinued"]

# ── COLOC PP.H4 best-per-gene (canonical per-GWAS file) ──────────────────────
co <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
            select = c("gene", "PP.H4.susie", "PP.H4.abf"))
co[, pp4 := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
co <- co[is.finite(pp4), .(coloc = max(pp4)), by = gene]
dt <- merge(cls, co, by = "gene", all.x = TRUE)
dt[is.na(coloc), coloc := 0]

# ── Cross-ancestry provenance flag (gene_level_coloc.csv) ────────────────────
# coloc_susie_headline_cross_anc = TRUE when a gene's SuSiE COLOC headline PP.H4 is
# driven by a NON-EUR GWAS. The eQTL panel is EUR, so a non-EUR-GWAS headline holds
# a lower evidentiary bar — open-mark ("*") these anchors in the plot + caption.
xanc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
              select = c("gene", "coloc_susie_headline_cross_anc"))
xanc <- xanc[, .(cross_anc = any(coloc_susie_headline_cross_anc %in% TRUE)), by = gene]
dt <- merge(dt, xanc, by = "gene", all.x = TRUE)
dt[is.na(cross_anc), cross_anc := FALSE]

# ── Bulk DEG (canonical limma-voom-qw C2, raw effect + TREAT FDR) ─────
deg <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("symbol", "logFC", "padj", "treat_lfc", "treat_fdr"))
setnames(deg, "symbol", "gene")
dt <- merge(dt, deg, by = "gene", all.x = TRUE)
dt <- dt[order(-coloc)][!duplicated(gene)]
dt[, abs_lfc := fifelse(is.na(logFC), 0, abs(logFC))]
dt[, plot_lfc := fifelse(is.na(logFC), 0, logFC)]
# Canonical DEG (2026-08-12) = padj < 0.05 AND |log2FC| > 0.5. Used both for reservoir
# INCLUSION below (the paper-wide dual-evidence gate) AND for the glyph outline.
dt[, is_deg := is_canonical_deg(dt)]
# BORDER/COLOUR = bulk-DEG significance = the SAME canonical gate.
# Genes with no DEG data (padj/logFC NA, e.g. GLP1R) are not significant and have no outline.
dt[, is_sig := is_deg]
dt[, sig_class := fifelse(is_sig, "DEG", "n.s.")]
dt[, sig_class := factor(sig_class, levels = c("DEG", "n.s."))]

# ── FAIRNESS / gene set ──────────────────────────────────────────────────────
#    pipeline  = CURATED roster of famous MASLD drugs (approved/trial/failed),
#    reservoir = method NOMINATIONS clearing the SAME paper-wide dual-evidence bar
#                (COLOC > 0.5 AND a TREAT bulk DEG, FDR<0.05 at lfc=0.25).
drug_anchors <- c("THRB", "GLP1R", "SLC5A2", "RORA", "DGAT2", "NR3C2", "VDR",
                  "PNPLA3", "HSD17B13", "NR1H4", "MAP3K5", "CCR2", "CCR5", "LOXL2",
                  "F2RL1")   # F2RL1 promoted to In-trials (NCT05680233); pin it so it isn't dropped from the reservoir set
hero_anchors <- c("HKDC1", "FADS2", "GPAM")
dt[, is_reservoir := drug_dev_status %in% c("discovery", "masld_preclinical")]
dt <- dt[ gene %in% c(drug_anchors, hero_anchors) |
          (is_reservoir & coloc > COLOC_MIN & is_deg) ]

# ── Hepatocyte-EXPRESSED substrate (OUR human scRNA atlas) = SHAPE channel ───────
# Judged directly from hep_mean_cpm + hep_ratio (Cas13_Library_Design/data/
# hep_specificity.csv), NOT the canned hep_substrate label. The old label's pure
# ratio>=1 rule let low-EXPRESSION genes pass on ratio alone (e.g. ASCL1 11 CPM
# ratio 4.7; SPTBN5 4.6 CPM) and demoted real hepatocyte genes whose ambient in
# another cell type exceeds them — so it is renamed "hepatocyte-EXPRESSED" (not
# "enriched") with an honest expression floor:
#   hep_mean_cpm >= HEP_CPM_FLOOR (real hepatocyte expression) AND
#   hep_ratio    >= HEP_RATIO_MIN (not dominated by another cell type's ambient).
# Encoded as the point SHAPE (diamond) so hep-credibility is visible across ALL
# points (relevant: the screen is hepatocyte-autonomous). Genes missing from the
# hep file -> "other".
HEP_CPM_FLOOR <- 10   # mean hepatocyte CPM floor
HEP_RATIO_MIN <- 1    # must not be dominated by another cell type's ambient
hep <- fread(file.path(BASE, "Cas13_Library_Design/data/hep_specificity.csv"),
             select = c("gene_symbol", "hep_mean_cpm", "max_other_cpm",
                        "hep_ratio", "top_other_celltype"))
setnames(hep, "gene_symbol", "gene")
dt <- merge(dt, hep, by = "gene", all.x = TRUE)
dt[, hep_expressed := !is.na(hep_mean_cpm) & hep_mean_cpm >= HEP_CPM_FLOOR &
                      !is.na(hep_ratio)    & hep_ratio    >= HEP_RATIO_MIN]
dt[, hep_class := fifelse(hep_expressed, "hepatocyte-expressed", "other")]
dt[, hep_class := factor(hep_class, levels = c("hepatocyte-expressed", "other"))]

# ── Plot ──────────────────────────────────────────────────────────────────────
# Muted onto the Fig 4c family, identical to Fig 5a clin_status_cols (2026-07-09):
# In trials amber #E0A94F, Failed red #B2182B, Preclinical blue #3B6EA5; Approved
# green and Discovery house-gray unchanged.
outcome_cols <- c("Approved" = "#2E7D32", "In trials" = "#E0A94F",
                  "Failed / discontinued" = "#B2182B",
                  "Preclinical" = "#3B6EA5", "Discovery" = "#9E9E9E")
dt[, outcome := factor(outcome, levels = names(outcome_cols))]
dt[, is_pipeline := outcome %in% c("Approved", "In trials", "Failed / discontinued")]
# Draw order (low -> high = back -> front). Reservoir at the back; famous pipeline
# drugs on top; among co-located pipeline drugs draw the "other" CIRCLES on top of
# the "hepatocyte-expressed" DIAMONDS, and GLP1R dead last — otherwise GLP1R's green
# "other" circle is fully hidden behind the PNPLA3/HSD17B13 diamonds and reads as a
# (hepatocyte-expressed) diamond.
dt[, draw_order := fcase(
  gene == "GLP1R",                          4L,
  is_pipeline & hep_class == "other",       3L,
  is_pipeline,                              2L,
  default =                                 1L)]
setorder(dt, draw_order)

# label set
lab_genes <- c("THRB", "GLP1R", "SLC5A2", "RORA", "NR1H4", "DGAT2", "NR3C2", "VDR",
               "PNPLA3", "HSD17B13", "MAP3K5", "CCR2", "CCR5", "LOXL2", "F2RL1",
               "HKDC1", "FADS2", "GPAM",
               # Discovery-gene labels. Hepatocyte EXPRESSION is judged from OUR human
               # scRNA substrate (Cas13_Library_Design/data/hep_specificity.csv) via the
               # honest CPM-floor + ratio rule (hep_mean_cpm >= 10 CPM AND hep_ratio >= 1),
               # NOT literature focus and NOT the canned hep_substrate label. Genes clearing
               # the floor read as "hepatocyte-expressed" (diamond): e.g. MLIP (507 CPM, 5.4x),
               # GNMT (44, 1.7x). Genes that passed the OLD ratio-only rule but fall UNDER the
               # 10-CPM floor are now "other" (e.g. ASCL1 11 CPM is borderline; SPTBN5 4.6 CPM
               # drops). SPINT2/CDH6/RECQL4/EFHD1/LRRC1 are expressed but not hep-specific.
               # Confirmed non-hepatocyte composition markers (left unlabelled): LYZ/KRT7/IFI30/
               # ULBP2 (immune/cholangiocyte); absent: PKP3/ANGPTL7/LINC01561.
               "SPINT2", "CDH6", "CFHR5", "MLIP", "RECQL4",
               "EFHD1", "GNMT", "LRRC1")

# Compact display set: retain every named target and every clinical pipeline anchor;
# omit only unlabeled discovery/preclinical background points. The full table remains
# in the audit CSV with is_plotted so the presentation filter is explicit.
dt[, is_labeled := gene %in% lab_genes]
dt[, is_plotted := is_labeled | is_pipeline]
plot_dt <- dt[is_plotted == TRUE]

# dual-evidence reservoir counts (COLOC>0.5 & TREAT DEG, FDR<0.05), per row
n_disc <- dt[drug_dev_status == "discovery"         & coloc > COLOC_MIN & is_deg, .N]
n_pre  <- dt[drug_dev_status == "masld_preclinical" & coloc > COLOC_MIN & is_deg, .N]

ymax <- max(abs(plot_dt$plot_lfc), na.rm = TRUE)
ylim <- c(-ymax - 0.15, ymax + 0.15)

# Declutter the COLOC~0 pile-up: several non-cis drugs (GLP1R / SLC5A2 / PNPLA3 /
# HSD17B13) share near-identical (x~0, y~0) coordinates and occlude each other — in
# particular GLP1R's green "other" CIRCLE was hidden behind the PNPLA3/HSD17B13
# "hepatocyte-expressed" DIAMONDS, making GLP1R read as a diamond. Small, seeded
# jitter; the x gate line is 0.5 (COLOC). The y +/-0.25 lines are visual TREAT
# interval-null references (DEG membership = treat_fdr, independent of y position). Draw
# GLP1R-type pipeline drugs last so they sit on top. Noted in the caption.
set.seed(42)
dt[, xj := pmin(pmax(coloc + runif(.N, -0.020, 0.020), -0.02), 1.07)]
dt[, yj := plot_lfc + runif(.N, -0.045, 0.045)]
dt[, label_nudge_x := fcase(
  gene == "GLP1R", -0.035, gene == "MAP3K5", -0.025,
  gene == "PNPLA3", 0.045, gene == "HSD17B13", 0.075,
  gene == "SLC5A2", 0.075, gene == "NR1H4", 0.055,
  default = 0)]
dt[, label_nudge_y := fcase(
  gene == "GLP1R", -0.075, gene == "MAP3K5", -0.105,
  gene == "PNPLA3", 0.105, gene == "HSD17B13", 0.050,
  gene == "SLC5A2", -0.115, gene == "NR1H4", -0.105,
  default = 0)]
plot_dt <- dt[is_plotted == TRUE]
label_dt <- plot_dt[is_labeled == TRUE]
# Trailing "*" on anchors whose SuSiE COLOC headline is non-EUR (cross-ancestry).
label_dt[, lab := fifelse(cross_anc == TRUE, paste0(gene, "*"), gene)]

LABEL_SIZE <- GEOM_TEXT_6PT   # true 6 pt in geom_text_repel size units

p <- ggplot(plot_dt, aes(xj, yj)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = c(-0.25, 0.25), linetype = "22", color = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0.5, linetype = "22", color = "grey75", linewidth = 0.3) +
  geom_point(aes(fill = outcome, shape = hep_class, color = sig_class),
             size = 2.45, stroke = 0.55, alpha = 0.92) +
  # open ring = non-EUR (cross-ancestry) COLOC headline (EUR eQTL panel — lower bar)
  geom_point(data = plot_dt[cross_anc == TRUE], aes(xj, yj),
             shape = 1, size = 3.7, stroke = 0.5, color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = label_dt,
                  aes(label = lab), color = "black", fontface = "italic", size = LABEL_SIZE,
                  box.padding = 0.22, point.padding = 0.14, min.segment.length = 0,
                  segment.size = 0.18, segment.color = "grey65", seed = 42,
                  nudge_x = label_dt$label_nudge_x, nudge_y = label_dt$label_nudge_y,
                  max.overlaps = Inf, max.time = 5, max.iter = 20000) +
  scale_fill_manual(values = outcome_cols, name = NULL, breaks = names(outcome_cols),
                    guide = guide_legend(order = 1,
                              override.aes = list(shape = 21, size = 2.6, alpha = 0.95, color = "grey40"))) +
  # SHAPE = hepatocyte-EXPRESSED in our scRNA (diamond) vs other (circle)
  scale_shape_manual(values = c("hepatocyte-expressed" = 23, "other" = 21), name = "scRNA substrate",
                     guide = guide_legend(order = 2,
                              override.aes = list(fill = "grey60", color = "grey40", size = 2.6))) +
  # BORDER/COLOUR = bulk-DEG significance (TREAT FDR<0.05) -> black ring; n.s. -> transparent
  # (no visible ring). MUST be "transparent", NOT NA: mapping a discrete level to NA makes
  # ggplot treat the whole row as missing and DROP the point (fill and all), not just border.
  scale_color_manual(values = c("DEG" = "black", "n.s." = "transparent"),
                     name = NULL,
                     guide = guide_legend(order = 3,
                              override.aes = list(shape = 21, fill = "grey70", size = 2.6, stroke = 0.7))) +
  scale_x_continuous(limits = c(-0.10, 1.08), breaks = c(0, 0.5, 1.0)) +
  scale_y_continuous(limits = ylim,
                     breaks = c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1)) +
  labs(x = "coloc probabilities (PP.H4)",
       y = "Transcriptomic effect (logFC)") +
  theme_masld_compact() +
  theme(legend.position = "right", legend.key.size = unit(0.2, "cm"),
        legend.box = "vertical", legend.spacing.y = unit(0.01, "cm"),
        text = element_text(family = "Helvetica", size = 6, face = "plain"),
        axis.title = element_text(size = 6, face = "plain"),
        axis.text = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        legend.text = element_text(size = 6, face = "plain"),
        plot.margin = margin(3, 4, 3, 3))

out <- file.path(FIG5_DIR, "panels", "fig6_drug_target_calibration.pdf")
save_fig(p, out, width = 3.31, height = 2.22)
message("Saved: ", out)
message("CAPTION: x = genetic colocalization (COLOC PP.H4); y = canonical raw logFC ",
        "(disease vs control). Fill = drug-development stage; SHAPE = hepatocyte ",
        "substrate (diamond = hepatocyte-expressed: scRNA hepatocyte mean CPM >=10 ",
        "AND hepatocyte:other-lineage ratio >=1; circle = other); black outline = ",
        "TREAT bulk DEG (FDR<0.05 at lfc=0.25), no outline = n.s. Dashed y lines at ",
        "+/-0.25 mark the TREAT interval-null bound. Points jittered slightly (seeded) ",
        "to separate low-COLOC drugs near (0,0). Unlabelled discovery/preclinical ",
        "background points are omitted from the compact display. Open ring + trailing ",
        "'*' mark genes whose SuSiE COLOC headline is non-EUR (cross-ancestry); the ",
        "eQTL panel is EUR, so those hold a lower evidentiary bar.")
.xanc <- label_dt[cross_anc == TRUE, gene]
message("[fig5c] cross-ancestry-flagged anchors: ",
        if (length(.xanc)) paste(.xanc, collapse = ", ") else "none")

fwrite(dt[order(-coloc, -plot_lfc),
          .(gene, drug_dev_status, outcome, coloc = round(coloc, 3),
            logFC = round(logFC, 3), padj = signif(padj, 3),
            treat_lfc = round(treat_lfc, 3), treat_fdr = signif(treat_fdr, 3),
            is_deg, sig_class, cross_anc, is_labeled, is_plotted,
            xj = round(xj, 4), yj = round(yj, 4),
            label_nudge_x, label_nudge_y,
            hep_mean_cpm = round(hep_mean_cpm, 3), hep_ratio = round(hep_ratio, 3),
            top_other_celltype, hep_class,
            max_phase_masld, n_masld_trials)],
       file.path(FIG5_DATA_DIR, "drug_target_calibration.csv"))
cat(sprintf("plotted: %d genes | dual-evidence discovery=%d  preclinical=%d  pipeline=%d\n",
            nrow(plot_dt), n_disc, n_pre, dt[is_pipeline == TRUE, .N]))
print(dt[is_pipeline == TRUE][order(-coloc, -plot_lfc),
         .(gene, outcome, coloc = round(coloc, 2), logFC = round(logFC, 2), is_deg)])

# ── sanity check: hep_class counts + labelled-gene categorization table ──────────
cat("\nhep_class counts:\n")
print(dt[, .N, by = hep_class])
cat("\nlabelled-gene sanity table:\n")
print(dt[gene %in% lab_genes][order(-coloc, -plot_lfc),
         .(gene, outcome, logFC = round(logFC, 3), treat_fdr = signif(treat_fdr, 3),
           hep_mean_cpm = round(hep_mean_cpm, 2), hep_ratio = round(hep_ratio, 2),
           hep_class, sig_class)])
