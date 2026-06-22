#!/usr/bin/env Rscript
# ============================================================================
# drug_target_orthogonality_landscape.R
# Main Fig 5 (convergence) — PANEL 5c — drug-target orthogonality landscape
# (relocated from Fig 3 RNA-seq 2026-06-18: this is a convergence / clinical-
# calibration story, not an RNA-seq story. Replaces the evidence-layer
# correlation heatmap in the 5c slot; that heatmap is demoted to figS08.)
#
# DRUG-DEVELOPMENT-GRADIENT VALIDATION (not "novel discovery"): the unsupervised
# convergence score, using NO drug/clinical labels, recovers real MASLD targets
# across the entire drug-development gradient, each externally validated AFTER the
# analysis. Approved / trial targets that act by NON-genetic mechanisms cluster in
# the high-disease-expression / no-colocalization quadrant (FASN, SCD, FGF21,
# GLP1R ~ PP.H4 0); the genetically-anchored clinical targets are THRB (resmetirom,
# approved) and RORA (TB-840, RORα agonist, Ph1 — recovered de novo, then entered
# first-in-human testing). HKDC1 and AKR1B10 are PRECLINICALLY VALIDATED anchors
# (hepatocyte Hkdc1-KO protects mice from MASH, Xu 2025 BBA PMID 40020530; AKR1B10
# inhibitors reduce NASH, PMID 34954342) — nominated de novo and independently
# corroborated, NOT novel/untested. HKDC1 stays prominent as the preclinical-stage
# validation hero. The discontinued targets (NR1H4/obeticholic, CCR2-CCR5/
# cenicriviroc, MAP3K5/selonsertib, LOXL2/simtuzumab) cluster off the causal axis.
#
# DATA-DRIVEN (2026-06-18 rebuild; drug_dev_status reframe 2026-06-19): target set
# is the canonical MASLD clinical drug pipeline (data/external/druggability/
# mash_clinical_pipeline.tsv, one point per gene) UNION preclinically-validated
# anchors (HKDC1, AKR1B10) whose drug_dev_status = masld_preclinical in
# data/external/drug_targets/drug_target_classification.tsv. RORA already carries
# clinical status from the pipeline. NOT a hand-curated gene list.
#
# x = bulk_shrunk_logFC (disease direction); y = colocalization PP.H4.
# Faint all-gene background cloud; fill = mechanism class; SHAPE = drug-development
# status (approved / active / discontinued / preclinical); ~14 genes text-labelled.
#
# INTEGRITY: PP.H4 read ONLY from gene_level_coloc.csv (coloc_best_pp4 <=1).
# Atlas coloc_susie_best_pp4 (corrupted, >1) NEVER used. Assert max<=1.
#
# Output: figures/main/fig5_convergence/panels/fig5c_drug_target_orthogonality.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG5_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "fig5c_drug_target_orthogonality.pdf")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
stopifnot(max(coloc$coloc_best_pp4, na.rm = TRUE) <= 1)

m <- merge(atlas[, .(human_symbol, bulk_shrunk_logFC)],
           coloc[, .(human_symbol = gene, coloc_best_pp4)],
           by = "human_symbol")
m <- m[!is.na(bulk_shrunk_logFC) & !is.na(coloc_best_pp4)]

# ---------------------------------------------------------------------------
# Comprehensive, data-driven target set: the canonical MASLD clinical pipeline
# (one point per gene) + preclinically-validated anchors. Mechanism class is
# mapped from the pipeline `mechanism` string; clinical status is collapsed to
# the best outcome per gene (approved > active > discontinued), then OVERRIDDEN
# for any gene the systematic drug_dev_status table flags masld_discontinued
# (that table is the single source of truth — e.g. NR1H4/obeticholic-acid is
# discontinued at the gene level even though one later FXR drug stays active).
# ---------------------------------------------------------------------------
pipe <- fread(file.path(BASE, "data/external/druggability/mash_clinical_pipeline.tsv"),
              fill = TRUE)

# mechanism string -> non-overlapping mechanism class
pipe[, mech_class := fcase(
  grepl("THR-beta", mechanism),                                            "Thyroid (THR-β)",
  grepl("FGF21|GLP-1|GCGR|GIP|klotho", mechanism, ignore.case = TRUE),     "Incretin / FGF21",
  grepl("PPAR|pyruvate carrier",       mechanism, ignore.case = TRUE),     "PPAR",
  grepl("DGAT2|SCD1|ACC inhibitor|FASN|ketohexokinase", mechanism, ignore.case = TRUE), "DNL / lipogenesis",
  grepl("FXR",                         mechanism),                         "FXR",
  grepl("ROR",                         mechanism),                         "RORα",
  grepl("CCR2|CCR5|ASK1|LOXL2",        mechanism, ignore.case = TRUE),     "Fibro-inflammatory",
  default =                                                                "Other metabolic / genetic"
)]

# collapse to one row per gene: best clinical outcome wins
status_rank <- c(approved = 3L, active = 2L, discontinued = 1L)
pipe[, status_r := status_rank[status]]
setorder(pipe, -status_r, -max_phase)
targets <- pipe[, .(status = status[1], max_phase = max(max_phase, na.rm = TRUE),
                    mech_class = mech_class[1]), by = .(human_symbol = gene_symbol)]

# Authority: the systematic drug_dev_status table is the single source of truth
# for drug-development status (spec docs/manuscript/working/drug_dev_status_reframe_spec.md).
# Override the pipeline best-outcome collapse for any gene it flags
# masld_discontinued, so e.g. NR1H4 (obeticholic acid, discontinued at the gene
# level) renders as 'discontinued' rather than picking up a later active FXR
# drug. This realises the NR1H4/OCA "discontinued = calibration win" narrative.
ddc <- fread(file.path(BASE,
  "data/external/drug_targets/drug_target_classification.tsv"))
ddc_disc <- ddc[drug_dev_status == "masld_discontinued", symbol]
targets[human_symbol %in% ddc_disc, status := "discontinued"]

# preclinically-validated anchors: nominated de novo by the convergence score,
# then independently corroborated in animal MASH studies (drug_dev_status =
# masld_preclinical). HKDC1 = preclinical-stage validation hero (hepatocyte
# Hkdc1-KO protects mice; Xu 2025 BBA PMID 40020530). AKR1B10 = AKR1B10
# inhibitors reduce NASH (PMID 34954342). Status is read from the systematic
# drug_dev_status table, not hardcoded by hand.
disc_genes <- c("HKDC1", "AKR1B10")
stopifnot(all(ddc[symbol %in% disc_genes, drug_dev_status] == "masld_preclinical"))
disc <- data.table(human_symbol = disc_genes,
                   status = "preclinical", max_phase = NA_real_,
                   mech_class = "Preclinical (validated)")
targets <- rbind(targets, disc[!human_symbol %in% targets$human_symbol])

# attach disease logFC (atlas) + colocalization (gene_level_coloc.csv ONLY)
atlas_u <- atlas[!duplicated(human_symbol), .(human_symbol, bulk_shrunk_logFC)]
coloc_u <- coloc[gene != "" & !is.na(gene) & !duplicated(gene),
                 .(human_symbol = gene, coloc_best_pp4)]
named <- merge(targets, atlas_u, by = "human_symbol", all.x = TRUE)
named <- merge(named,   coloc_u, by = "human_symbol", all.x = TRUE)
# absent from the tested eGene set -> no colocalization signal (plot at PP.H4 = 0)
named[is.na(coloc_best_pp4), coloc_best_pp4 := 0]
miss_x <- named[is.na(bulk_shrunk_logFC), human_symbol]
if (length(miss_x))
  cat(sprintf("[warn] %d target(s) lack bulk logFC, dropped: %s\n",
              length(miss_x), paste(miss_x, collapse = ", ")))
named <- named[!is.na(bulk_shrunk_logFC)]

mech_levels <- c("Thyroid (THR-β)", "Incretin / FGF21", "PPAR", "DNL / lipogenesis",
                 "FXR", "RORα", "Fibro-inflammatory", "Other metabolic / genetic", "Preclinical (validated)")
named[, mech_class := factor(mech_class, levels = mech_levels)]
named[, status     := factor(status, levels = c("approved", "active", "discontinued", "preclinical"))]

# Selective labels (~13) spanning the calibration story.
# NB: GLP1R (semaglutide) is intentionally absent — it is not hepatically
# expressed (systemic incretin mechanism), so it is not in the hepatic atlas
# and cannot be placed on a disease-logFC axis (see [warn] above).
label_set <- c("THRB", "RORA", "HKDC1", "DGAT2",            # anchored / top
               "FASN", "SCD", "FGF21", "PNPLA3",            # un-anchored active pipeline
               "NR1H4", "CCR2", "MAP3K5", "LOXL2",          # failed / discontinued
               "AKR1B10")                                   # preclinical (validated)
named_lab <- named[human_symbol %in% label_set]

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
gp <- function(g) { v <- named[human_symbol == g, coloc_best_pp4]
                    if (length(v)) v[1] else NA_real_ }
cat(sprintf("[hero] %d targets plotted (%d pipeline + %d preclinical-validated anchor)\n",
            nrow(named), sum(named$status != "preclinical"), sum(named$status == "preclinical")))
cat(sprintf("[hero] status counts: %s\n",
            paste(names(table(named$status)), as.integer(table(named$status)),
                  sep = "=", collapse = " | ")))
cat(sprintf("[hero] anchors: THRB=%.3f RORA=%.3f HKDC1=%.3f | un-anchored: GLP1R=%.3f FASN=%.3f FGF21=%.3f\n",
            gp("THRB"), gp("RORA"), gp("HKDC1"), gp("GLP1R"), gp("FASN"), gp("FGF21")))
cat(sprintf("[hero] discontinued PP.H4: NR1H4=%.3f CCR2=%.3f MAP3K5=%.3f LOXL2=%.3f\n",
            gp("NR1H4"), gp("CCR2"), gp("MAP3K5"), gp("LOXL2")))

fwrite(named[order(-coloc_best_pp4)], file.path(DATA_DIR, "drug_target_orthogonality_landscape.csv"))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
mech_cols <- c(
  `Thyroid (THR-β)`           = "#00695C",
  `Incretin / FGF21`          = "#F39C12",
  `PPAR`                      = "#1565C0",
  `DNL / lipogenesis`         = "#C9265E",
  `FXR`                       = "#7B1FA2",
  `RORα`                      = "#AD1457",
  `Fibro-inflammatory`        = "#795548",
  `Other metabolic / genetic` = "#9E9E9E",
  `Preclinical (validated)`   = "#212121"
)
# Shape = drug-development status (all fillable shapes so the mechanism fill always reads)
status_shapes <- c(approved = 22, active = 21, discontinued = 25, preclinical = 23)

p <- ggplot() +
  geom_point(data = m, aes(x = bulk_shrunk_logFC, y = coloc_best_pp4),
             color = "grey85", size = 0.4, alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.25, color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, color = "grey55") +
  geom_point(data = named, aes(x = bulk_shrunk_logFC, y = coloc_best_pp4,
                               fill = mech_class, shape = status),
             size = 2.3, color = "white", stroke = 0.3) +
  ggrepel::geom_text_repel(
    data = named_lab, aes(x = bulk_shrunk_logFC, y = coloc_best_pp4,
                      label = human_symbol, color = mech_class),
    size = 2.1, fontface = "italic",
    segment.size = 0.25, segment.color = "grey55", min.segment.length = 0,
    max.overlaps = Inf, box.padding = 0.8, point.padding = 0.45,
    force = 12, force_pull = 0.2, max.iter = 1e5, max.time = 1, seed = 7,
    show.legend = FALSE) +
  scale_fill_manual(values = mech_cols, name = "Mechanism", drop = FALSE) +
  scale_color_manual(values = mech_cols, guide = "none") +
  scale_shape_manual(values = status_shapes, name = "Drug-development status", drop = FALSE) +
  scale_y_continuous(name = "Colocalization PP.H4", limits = c(0, 1.02),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_x_continuous(name = expression("Disease bulk log"[2]*"FC (shrunk)"),
                     expand = expansion(mult = c(0.10, 0.02))) +
  labs(title = "c", subtitle = "Drug targets: expression vs genetic anchoring") +
  theme_masld(base_size = 7) +
  theme(
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 6.3, color = "grey30", margin = margin(b = 4)),
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.01, "cm"),
    legend.margin    = margin(t = 0, b = 0),
    legend.title     = element_text(size = 5.8, face = "bold"),
    legend.text      = element_text(size = 5),
    legend.key.size  = unit(0.16, "cm")
  ) +
  guides(
    fill  = guide_legend(order = 1, nrow = 3, byrow = TRUE,
                         override.aes = list(shape = 21, size = 2)),
    shape = guide_legend(order = 2, nrow = 1,
                         override.aes = list(fill = "grey45", size = 2))
  )

ggsave(OUT_PDF, p,
       width  = 88 / 25.4,
       height = 92 / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
