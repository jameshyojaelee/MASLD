#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fig 4c — the CALIBRATION keystone. The unsupervised prioritizer's genetic
# evidence (cis-eQTL colocalization) tracks real clinical fate across the MASLD
# drug-development gradient, with NO drug labels used to build it: the one
# FDA-APPROVED target (THRB / resmetirom) sits at the top with the strongest
# colocalization, the Phase-3 FAILURE (NR1H4 / obeticholic acid) has null genetic
# support, and the de-novo discovery now in Phase 1 (RORA / TB-840) colocalizes as
# strongly as THRB. Crucially, COLOC is SPECIFIC not high-recall: in-trials drugs
# acting through coding/splice mechanisms (PNPLA3, HSD17B13 ASOs) sit at low COLOC,
# the honest calibrated behaviour.
#
# DESIGN (2026-06-19 v3): driven genome-wide from the new drug-development-status
# classification (data/external/drug_targets/drug_target_classification.tsv),
# NOT the ~30-gene hand-curated table. Adds the two requested lower tiers —
# in-trials AND the high-COLOC "discovery / preclinical" reservoir — so the
# calibration ladder runs Discovery → Preclinical → Ph1 → Ph2 → Ph3 → Approved.
# x = COLOC PP.H4, y = development stage, colour = drug-dev outcome class,
# reservoir alpha emphasises strong (>0.8) colocalization. NO lollipop
# (memory/feedback-no-lollipop); all text BLACK (memory/feedback-no-colored-fonts).
#
# DATA (canonical):
#   drug_dev_status from data/external/drug_targets/drug_target_classification.tsv
#     (genome-wide; encodes curated overrides incl. RORA→Ph1 TB-840 NCT05045534).
#   COLOC PP.H4 best-per-gene from susie_coloc_all_gwas.csv (NEVER atlas conv cols).
#   Categories plotted (per request): masld_approved / masld_clinical /
#     masld_discontinued / masld_preclinical / discovery. drugged_other_indication
#     + undetermined EXCLUDED (not a MASLD development tier).
#
# Output: figures/main/fig4_validation/fig4c_drug_target_calibration_v2.pdf
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

# ── Genome-wide drug-development status ──────────────────────────────────────
cls <- fread(file.path(BASE, "data/external/drug_targets/drug_target_classification.tsv"))
setnames(cls, "symbol", "gene")
keep_status <- c("masld_approved", "masld_clinical", "masld_discontinued",
                 "masld_preclinical", "discovery")
cls <- cls[drug_dev_status %in% keep_status & gene != "" & !is.na(gene)]

# stage ordinal (y) + outcome class (colour) from status + best MASLD phase
# y ladder: 1 Discovery · 2 Preclinical · 3 Ph1 · 4 Ph2 · 5 Ph3 · 6 Approved
phase_to_y <- function(p) ifelse(is.na(p), NA_integer_, as.integer(p) + 2L)
cls[, stage_y := NA_integer_]
cls[, outcome := NA_character_]
cls[drug_dev_status == "discovery",         `:=`(stage_y = 1L, outcome = "Discovery")]
cls[drug_dev_status == "masld_preclinical", `:=`(stage_y = 2L, outcome = "Preclinical")]
# in-trials: place at best MASLD phase; NA-phase-but-has-trials defaults to Ph1
cls[drug_dev_status == "masld_clinical",
    `:=`(stage_y = fifelse(!is.na(max_phase_masld), phase_to_y(max_phase_masld), 3L),
         outcome  = "In trials")]
cls[drug_dev_status == "masld_approved",
    `:=`(stage_y = fifelse(!is.na(max_phase_masld), phase_to_y(max_phase_masld), 6L),
         outcome  = "Approved")]
cls[drug_dev_status == "masld_discontinued",
    `:=`(stage_y = fifelse(!is.na(max_phase_masld), phase_to_y(max_phase_masld), 5L),
         outcome  = "Failed / discontinued")]

# ── COLOC PP.H4 best-per-gene (canonical per-GWAS file) ──────────────────────
co <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
            select = c("gene", "PP.H4.susie", "PP.H4.abf"))
co[, pp4 := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
co <- co[is.finite(pp4), .(coloc = max(pp4)), by = gene]
dt <- merge(cls, co, by = "gene", all.x = TRUE)
dt[is.na(coloc), coloc := 0]

# ── Bulk DEG (canonical limma-voom-qw C2, ashr-shrunk) — the SECOND evidence axis ─
deg <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("symbol", "shrunk_logFC", "lfsr"))
setnames(deg, "symbol", "gene")
dt <- merge(dt, deg, by = "gene", all.x = TRUE)
dt <- dt[order(-coloc)][!duplicated(gene)]
dt[, abs_lfc := fifelse(is.na(shrunk_logFC), 0, abs(shrunk_logFC))]
# Tier-1 DEG = the paper's canonical bar (lfsr < 0.05 & |shrunk_logFC| > 0.5)
dt[, is_deg := !is.na(lfsr) & lfsr < 0.05 & abs_lfc > 0.5]

# ── FAIRNESS: the reservoir is NOT "every gene above PP.H4>0.5" (that cherry-picks
#    colocalizing genes and floods the plot). Both sides are principled sets:
#    · pipeline  = a CURATED roster of the famous MASLD drugs (approved/trial/failed),
#    · reservoir = genes our method NOMINATES — discovery/preclinical that clear the
#      SAME dual-evidence bar used paper-wide: COLOC > 0.5 AND a Tier-1 bulk DEG.
# Drug anchors (always shown, incl. the low-COLOC non-cis contrasts):
#   PNPLA3/HSD17B13 = coding-ASO (in trials, ~0 COLOC); GLP1R/SLC5A2 = approved but
#   incretin/SGLT2 (~0 COLOC, "specific not high-recall" vs THRB); NR1H4/MAP3K5/CCR2/
#   CCR5/LOXL2 = the failed-drug cluster (selonsertib/cenicriviroc/simtuzumab/OCA).
drug_anchors <- c("THRB", "GLP1R", "SLC5A2", "RORA", "DGAT2", "NR3C2", "VDR",
                  "PNPLA3", "HSD17B13", "NR1H4", "MAP3K5", "CCR2", "CCR5", "LOXL2")
# featured preclinical heroes — shown + labelled even if below the Tier-1 DEG bar
# (high COLOC, discussed in text); their true DEG strength is read off point size.
hero_anchors <- c("HKDC1", "FADS2", "GPAM")
dt[, is_reservoir := drug_dev_status %in% c("discovery", "masld_preclinical")]
dt <- dt[ gene %in% c(drug_anchors, hero_anchors) |
          (is_reservoir & coloc > COLOC_MIN & is_deg) ]

# ── Plot ──────────────────────────────────────────────────────────────────────
outcome_cols <- c("Approved" = "#2E7D32", "In trials" = "#F9A825",
                  "Failed / discontinued" = "#C62828",
                  "Preclinical" = "#1565C0", "Discovery" = "#9E9E9E")
dt[, outcome := factor(outcome, levels = names(outcome_cols))]
dt[, is_pipeline := outcome %in% c("Approved", "In trials", "Failed / discontinued")]
# DEG shown by SHAPE (not size, which was confusing): ▲ up-DEG / ▼ down-DEG / ● n.s.
# Shape uses SIGNIFICANCE + direction (lfsr<0.05) so modest-effect-but-significant
# targets (THRB, RORA suppressed in disease) read correctly; the stricter Tier-1 bar
# (|logFC|>0.5, `is_deg`) is used only for reservoir INCLUSION, not the glyph.
dt[, deg_sig := !is.na(lfsr) & lfsr < 0.05]
dt[, deg_dir := fcase(deg_sig & shrunk_logFC > 0, "DEG up",
                      deg_sig & shrunk_logFC < 0, "DEG down",
                      default = "not sig.")]
dt[, deg_dir := factor(deg_dir, levels = c("DEG up", "DEG down", "not sig."))]
deg_shapes <- c("DEG up" = 24, "DEG down" = 25, "not sig." = 21)
set.seed(42)
dt[, yj := stage_y + runif(.N, -0.18, 0.18)]
setorder(dt, is_pipeline)   # draw reservoir first, pipeline drugs on top

# CCR5 left unlabeled — it overlaps CCR2 (same drug, cenicriviroc); CCR2 stands in.
# Discovery labels = top dual-evidence (COLOC+DEG) discovery genes by COLOC, skipping
# cell-composition immune/cholangiocyte markers (LYZ/KRT7/IFI30/ULBP2) and EFHD1 (Fig 2).
# NOVELTY CLAIM IS MASLD-SPECIFIC (the y-axis is the MASLD drug-development ladder):
# these carry MASLD genetic (COLOC) + transcriptomic (DEG) evidence but have NO MASH/MASLD
# therapeutic programme. PubMed-verified 2026-06-19: SPINT2 = HCC tumour suppressor,
# RECQL4 = HCC target, CDH6 = oncology ADC (raludotatug deruxtecan) — i.e. tractable but
# never pursued in MASH. Caption must say "not yet tested in MASH/MASLD", NOT "unstudied".
lab_genes <- c("THRB", "GLP1R", "SLC5A2", "RORA", "NR1H4", "DGAT2", "NR3C2", "VDR",
               "PNPLA3", "HSD17B13", "MAP3K5", "CCR2", "LOXL2",
               "HKDC1", "FADS2", "GPAM",
               "SPINT2", "CDH6", "CFHR5", "MLIP", "RECQL4")

# dual-evidence reservoir counts (COLOC>0.5 & Tier-1 DEG), per row
n_disc <- dt[drug_dev_status == "discovery"        & coloc > COLOC_MIN & is_deg, .N]
n_pre  <- dt[drug_dev_status == "masld_preclinical" & coloc > COLOC_MIN & is_deg, .N]

p <- ggplot(dt, aes(coloc, yj)) +
  geom_vline(xintercept = 0.5, linetype = "22", color = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = 0.8, linetype = "22", color = "grey85", linewidth = 0.25) +
  geom_point(aes(fill = outcome, shape = deg_dir),
             size = 2.6, color = "white", stroke = 0.3, alpha = 0.92) +
  geom_text_repel(data = dt[gene %in% lab_genes],
                  aes(label = gene), color = "black", fontface = "italic", size = 2.2,
                  box.padding = 0.34, point.padding = 0.2, min.segment.length = 0,
                  segment.size = 0.2, segment.color = "grey65", seed = 42, max.overlaps = Inf) +
  scale_fill_manual(values = outcome_cols, name = NULL, breaks = names(outcome_cols),
                    guide = guide_legend(order = 1,
                              override.aes = list(shape = 21, size = 2.4, alpha = 0.95))) +
  scale_shape_manual(values = deg_shapes, name = "bulk DEG",
                     guide = guide_legend(order = 2,
                              override.aes = list(fill = "grey45", size = 2.4))) +
  scale_x_continuous(limits = c(-0.03, 1.08), breaks = c(0, 0.5, 0.8, 1.0)) +
  scale_y_continuous(breaks = 1:6,
                     labels = c("Discovery", "Preclinical", "Phase 1", "Phase 2",
                                "Phase 3", "FDA approved"),
                     limits = c(0.55, 6.45)) +
  labs(x = "Genetic colocalization (COLOC PP.H4)", y = "Drug-development stage") +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(size = PUB_TITLE, face = "bold", color = "black"),
        legend.position = "right", legend.key.size = unit(0.2, "cm"),
        legend.box = "vertical", legend.spacing.y = unit(0.01, "cm"),
        plot.margin = margin(3, 4, 3, 3))

out <- file.path(FIG4_DIR, "fig4c_drug_target_calibration_v2.pdf")
save_fig(p, out, width = fig_half_width + 1.25, height = 2.7)
message("Saved: ", out)

fwrite(dt[order(-stage_y, -coloc),
          .(gene, drug_dev_status, outcome, stage_y, coloc = round(coloc, 3),
            shrunk_logFC = round(shrunk_logFC, 3), lfsr = signif(lfsr, 3), is_deg,
            max_phase_masld, n_masld_trials)],
       file.path(FIG4_DIR, "data", "drug_target_calibration.csv"))
cat(sprintf("plotted: %d genes | dual-evidence discovery=%d  preclinical=%d  pipeline=%d\n",
            nrow(dt), n_disc, n_pre, dt[is_pipeline == TRUE, .N]))
print(dt[is_pipeline == TRUE][order(-stage_y, -coloc),
         .(gene, outcome, coloc = round(coloc, 2), logFC = round(shrunk_logFC, 2), is_deg)])
