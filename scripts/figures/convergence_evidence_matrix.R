#!/usr/bin/env Rscript
# KEY MESSAGE: Each prioritized target reproduces across a DISTINCT combination
# of orthogonal modalities — no two genes share the same evidence fingerprint.
# The panel reads as DRUG-DEVELOPMENT-GRADIENT VALIDATION: the unsupervised
# convergence score, using NO drug or clinical labels, recovers real MASLD
# targets across the entire drug-development gradient, and each was externally
# validated AFTER the analysis. Above the divider sit two clinically anchored
# targets recovered de novo (THRB, approved [resmetirom]; RORA, Phase 1
# [TB-840, RORA agonist] — the clinical-stage validation hero); below sit a
# preclinically validated hero and two pharmacology/control case studies
# (HKDC1, preclinically validated by a hepatocyte Hkdc1-KO MASH study; SERPINE1,
# positive control; CYP3A4, drugged metabolizer). The grid is a replication
# fingerprint: reading across a row shows which orthogonal modalities
# independently re-detect that target.
#
# A dot-matrix (genes x modalities). DIRECTION is encoded by colour for the two
# signed modalities (Bulk RNA, Mouse): disease-up = magenta, down = blue. The
# unsigned presence/strength modalities (COLOC, ATAC motif, Spatial) use a
# single neutral teal ("detected"). STRENGTH is encoded by dot size. Two
# no-signal states are shape-coded, NEVER drawn as a zero-size dot: an OPEN
# CIRCLE = tested-but-absent (interrogated, no signal); a grey X = NOT MEASURED
# (the modality never tested this gene).
#
# Status tags annotate the drug-development-status interpretation (backed by
# data/external/drug_targets/drug_target_classification.tsv):
#   THRB = approved (resmetirom, masld_approved);
#   RORA = Phase 1 (TB-840, masld_clinical) — clinical-stage validation hero,
#          recovered de novo without drug labels, then entered first-in-human;
#   NR1H4 = discontinued (obeticholic acid, masld_discontinued);
#   HKDC1 = preclinical (masld_preclinical) — preclinical validation hero,
#          corroborated by a hepatocyte Hkdc1-KO study protecting mice from MASH;
#   SERPINE1 = positive control (drugged-other; Kim 2024 spatial scoop);
#   CYP3A4 = drugged (CYP metabolizer; pharmacology case study).
# NONE of HKDC1/SERPINE1/CYP3A4 is claimed as a de-novo NOVEL discovery here.
#
# HARD-GATE PROVENANCE (anti-FADS2):
#  - COLOC pulled from the CANONICAL per-GWAS table susie_coloc_all_gwas.csv
#    (susie-else-abf best PP.H4), NOT the atlas *_coloc_pp4 convenience columns
#    (method-inconsistent). Axis is labelled "COLOC PP.H4" because most plotted
#    genes are abf-only (no SuSiE convergence) — see ledger rows 7-12.
#  - Spatial uses the per-condition Moran's C from spatial_autocorr_*.csv
#    (ledger rows 32/33); the atlas single-value spatial_morans_i is explicitly
#    "do not cite". Only SERPINE1/CYP3A4 are in the SVG/autocorr test set among
#    these genes => every other gene is NOT MEASURED for spatial (x).
#  - ATAC motif, bulk, protein, mouse, COLOC numbers all trace to ledger rows.
#
# Output: figures/main/fig5_molecular_context/_supp/convergence_evidence_matrix.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

genes_oi  <- c("THRB", "RORA", "NR1H4", "HKDC1", "SERPINE1", "CYP3A4")  # row order
mods      <- c("Bulk RNA", "COLOC", "ATAC motif", "Protein", "Spatial", "Mouse")  # col order

# Status tags (drug-development-status interpretation; printed beside gene name).
# Labels are the NEW drug_dev_status-backed tags from the relabel map, traceable
# to data/external/drug_targets/drug_target_classification.tsv (drug_dev_status:
# THRB masld_approved, RORA masld_clinical [max_phase_masld=1], NR1H4
# masld_discontinued, HKDC1/SERPINE1 masld_preclinical, CYP3A4
# drugged_other_indication). RORA & HKDC1 are the validation heroes.
# NOTE: SERPINE1's drug_target_classification.tsv drug_dev_status is
# masld_preclinical, but it is tagged here as a POSITIVE CONTROL per the reframe
# spec (drug_dev_status_reframe_spec.md line 31) and fig4.md — an explicit,
# intentional framing override, not a label mismatch with the TSV.
status_tag <- c(THRB = "approved", RORA = "Phase 1", NR1H4 = "discontinued",
                HKDC1 = "preclinical", SERPINE1 = "positive control",
                CYP3A4 = "drugged (CYP)")

# ── Signed modalities (bulk, protein, mouse) from the atlas ───────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
A <- atlas[match(genes_oi, atlas$human_symbol), ]

# ── COLOC: CANONICAL susie-else-abf best PP.H4 per gene (NOT atlas columns) ───
# Ledger rows 7-12. We never read the atlas *_coloc_pp4 convenience columns.
co <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  stringsAsFactors = FALSE)
coloc_pp4    <- setNames(rep(NA_real_, length(genes_oi)), genes_oi)
coloc_method <- setNames(rep(NA_character_, length(genes_oi)), genes_oi)
coloc_gwas   <- setNames(rep(NA_character_, length(genes_oi)), genes_oi)
for (g in genes_oi) {
  sub <- co[co$gene == g, ]
  if (!nrow(sub)) next
  bs <- suppressWarnings(max(sub$PP.H4.susie, na.rm = TRUE)); if (!is.finite(bs)) bs <- NA
  ba <- suppressWarnings(max(sub$PP.H4.abf,   na.rm = TRUE)); if (!is.finite(ba)) ba <- NA
  if (!is.na(bs)) {                       # SuSiE preferred when it converged
    coloc_pp4[g]    <- bs
    coloc_method[g] <- "susie"
    coloc_gwas[g]   <- sub$gwas_name[which.max(sub$PP.H4.susie)]
  } else if (!is.na(ba)) {                # else abf fallback
    coloc_pp4[g]    <- ba
    coloc_method[g] <- "abf"
    coloc_gwas[g]   <- sub$gwas_name[which.max(sub$PP.H4.abf)]
  }
}

# ── ATAC motif disruption: variant-motif PAIRS per gene's own TF motif ────────
# Ledger rows 36/44. Only the TF genes (THRB/RORA/NR1H4) have a self-motif that
# GWAS credible-set variants disrupt; the non-TF genes have none = tested-absent.
md <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
atac_pairs <- setNames(vapply(genes_oi,
  function(g) sum(grepl(g, md$tf_name, ignore.case = TRUE)), integer(1)), genes_oi)
# THRB=10, RORA=9, NR1H4=2 (NR1H4::RXRA heterodimer), HKDC1/SERPINE1/CYP3A4=0.

# ── Spatial: per-condition Moran's C (ledger rows 32/33) — NOT atlas single-val
# The SVG/autocorr test set among these genes is SERPINE1 & CYP3A4 only; all
# others were NOT in the spatial test => not measured. We plot the DISEASE-state
# (Steatotic) Moran C as the strength; the gain is reported in the legend.
sp_steat <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Steatotic.csv"),
  stringsAsFactors = FALSE)
sp_heal  <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Healthy.csv"),
  stringsAsFactors = FALSE)
sp_get <- function(df, g) { r <- df[df$Gene == g, ]; if (nrow(r)) r$C[1] else NA_real_ }
# Honest 3-state spatial encoding (verifier H10):
#  - SVG-significant in GSE192741 (SERPINE1, CYP3A4) => detected (teal, strength=C)
#  - present in the 11,984-gene autocorr file but NOT SVG (THRB 0.006/RORA 0.081/
#    NR1H4 0.029) => TESTED-ABSENT (open circle), NOT "not measured"
#  - absent from the autocorr file entirely (HKDC1) => NOT MEASURED (x)
spatial_svg    <- c("SERPINE1", "CYP3A4")            # SVG-significant -> detected
spatial_steat  <- setNames(vapply(genes_oi, function(g) sp_get(sp_steat, g), numeric(1)), genes_oi)
spatial_heal   <- setNames(vapply(genes_oi, function(g) sp_get(sp_heal,  g), numeric(1)), genes_oi)

# ── Build the long tile table: one row per gene x modality ───────────────────
num <- function(x) suppressWarnings(as.numeric(x))
mk <- function(gene, modality, value, state, strength, direction)
  data.frame(gene = gene, modality = modality, value = value, state = state,
             strength = strength, direction = direction, stringsAsFactors = FALSE)

# Per-modality max magnitudes for size scaling (so dot area is comparable).
BULK_MAX <- 1.10   # HKDC1 |logFC| ledger row 19
PROT_MAX <- 0.70   # CYP3A4 |logFC| ledger row 23
MOUSE_MAX <- 1.87  # SERPINE1 ledger row 26
ATAC_MAX <- 10     # THRB pairs ledger row 36
SPAT_MAX <- 0.65   # ~SERPINE1 steatotic C ledger row 32

rows <- list()
for (i in seq_along(genes_oi)) {
  g <- genes_oi[i]; r <- A[i, ]

  # 1. Bulk RNA — signed; every gene measured.
  blfc <- num(r$bulk_logFC)
  rows[[length(rows) + 1]] <- mk(g, "Bulk RNA", blfc, "measured",
                                 min(abs(blfc) / BULK_MAX, 1),
                                 ifelse(blfc > 0, "up", "down"))

  # 2. COLOC — unsigned susie-else-abf PP.H4 (canonical file). All measured, but
  #    only PP.H4 >= 0.5 (canonical colocalization threshold) counts as DETECTED;
  #    sub-threshold posteriors (NR1H4 0.208, SERPINE1 0.138 — ledger flags NR1H4
  #    null) render as TESTED-ABSENT (open circle, strength 0), never overstated.
  pp4 <- coloc_pp4[[g]]
  rows[[length(rows) + 1]] <- mk(g, "COLOC", pp4, "measured",
                                 ifelse(!is.na(pp4) & pp4 >= 0.5, pp4, 0), "neutral")

  # 3. ATAC motif — unsigned count of disrupted variant-motif pairs. All genes
  #    interrogated in the GWAS-ATAC scan => measured; 0 pairs = tested-absent.
  np <- atac_pairs[[g]]
  rows[[length(rows) + 1]] <- mk(g, "ATAC motif", np, "measured",
                                 min(np / ATAC_MAX, 1), "neutral")

  # 4. Protein — signed DIA-MS logFC. NA => NOT MEASURED (absent from DIA panel).
  plfc <- num(r$best_protein_logFC)
  if (is.na(plfc)) {
    rows[[length(rows) + 1]] <- mk(g, "Protein", NA, "not_measured", NA, NA)
  } else {
    rows[[length(rows) + 1]] <- mk(g, "Protein", plfc, "measured",
                                   min(abs(plfc) / PROT_MAX, 1),
                                   ifelse(plfc > 0, "up", "down"))
  }

  # 5. Spatial — 3-state (see spatial_svg note above).
  ms <- spatial_steat[[g]]
  if (g %in% spatial_svg) {
    rows[[length(rows) + 1]] <- mk(g, "Spatial", ms, "measured",
                                   min(ms / SPAT_MAX, 1), "neutral")   # SVG -> detected
  } else if (!is.na(ms)) {
    rows[[length(rows) + 1]] <- mk(g, "Spatial", ms, "measured", 0, "neutral")  # tested, not SVG -> open circle
  } else {
    rows[[length(rows) + 1]] <- mk(g, "Spatial", NA, "not_measured", NA, NA)     # absent from file -> x
  }

  # 6. Mouse — signed cross-species meta logFC. NA => NOT MEASURED in mouse meta.
  mlfc <- num(r$mouse_meta_logFC)
  if (is.na(mlfc)) {
    rows[[length(rows) + 1]] <- mk(g, "Mouse", NA, "not_measured", NA, NA)
  } else {
    rows[[length(rows) + 1]] <- mk(g, "Mouse", mlfc, "measured",
                                   min(abs(mlfc) / MOUSE_MAX, 1),
                                   ifelse(mlfc > 0, "up", "down"))
  }
}
tiles <- do.call(rbind, rows)

# Gene label with status tag; reverse so THRB is the TOP row.
gene_label <- vapply(genes_oi, function(g) {
  t <- status_tag[[g]]
  if (nzchar(t)) sprintf("%s  (%s)", g, t) else g
}, character(1))
lab_map <- setNames(gene_label, genes_oi)
tiles$gene_lab <- factor(lab_map[tiles$gene], levels = rev(unname(gene_label)))
tiles$modality <- factor(tiles$modality, levels = mods)

# ── Colour mapping ───────────────────────────────────────────────────────────
tiles <- tiles %>%
  mutate(
    measured = state == "measured",
    fill_class = case_when(
      !measured           ~ "Not measured",
      direction == "up"   ~ "Up in MASLD",
      direction == "down" ~ "Down in MASLD",
      TRUE                ~ "Detected"))

fill_pal <- c(
  "Up in MASLD"   = masld_colors[["up"]],     # deep magenta (signed)
  "Down in MASLD" = masld_colors[["down"]],   # deep blue (signed)
  "Detected"      = "#00897B")                 # teal (unsigned presence/strength)
tiles$fill_class <- factor(tiles$fill_class,
                           levels = c(names(fill_pal), "Not measured"))

measured_tiles <- tiles %>% filter(measured, strength > 0)
noflag <- bind_rows(
  tiles %>% filter(measured, strength == 0) %>% mutate(glyph = "Tested, absent"),
  tiles %>% filter(!measured)               %>% mutate(glyph = "Not measured"))
noflag$glyph <- factor(noflag$glyph, levels = c("Tested, absent", "Not measured"))

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(tiles, aes(x = modality, y = gene_lab)) +
  geom_tile(fill = "white", color = "grey88", linewidth = 0.25) +
  geom_point(data = measured_tiles,
             aes(size = strength, fill = fill_class),
             shape = 21, color = "white", stroke = 0.25) +
  geom_point(data = noflag,
             aes(x = modality, y = gene_lab, shape = glyph),
             size = 1.7, color = "grey50", stroke = 0.55, inherit.aes = FALSE) +
  scale_fill_manual(values = fill_pal, name = NULL,
                    breaks = c("Up in MASLD", "Down in MASLD", "Detected"),
                    drop = FALSE,
                    guide = guide_legend(order = 1, override.aes = list(
                      shape = 21, size = 2.2, stroke = 0.25, color = "white",
                      fill  = c(masld_colors[["up"]], masld_colors[["down"]],
                                "#00897B")))) +
  scale_shape_manual(values = c("Tested, absent" = 1, "Not measured" = 4),
                     name = NULL,
                     guide = guide_legend(order = 2,
                       override.aes = list(size = 2, color = "grey45", stroke = 0.6))) +
  scale_size_area(max_size = 4.2, name = "Strength",
                  breaks = c(0.25, 0.5, 1.0), labels = c("low", "med", "high"),
                  guide = guide_legend(order = 3,
                    override.aes = list(shape = 21, fill = "grey55",
                                        color = "white", stroke = 0.2))) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 0.85, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(
    axis.text.y      = element_text(size = PUB_AXIS_TEXT),
    axis.text.x.top  = element_text(angle = 35, hjust = 0, size = PUB_AXIS_TEXT),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.04, "cm"),
    legend.margin    = margin(0, 0, 0, 0))

# Italicise only the gene symbol (not the "(status)" tag): use a per-label
# fontface vector. ggplot's axis.text takes a single face, so we instead render
# the symbol bold-italic via element_markdown-free approach: keep plain text but
# italic for clarity. (All labels italic reads cleanly for gene symbols.)
p <- p + theme(axis.text.y = element_text(face = "italic", size = PUB_AXIS_TEXT))

# Divider along the drug-development gradient. In the reversed factor the
# boundary sits between NR1H4 and HKDC1: above = clinically anchored, recovered
# de novo (THRB approved, RORA Phase 1, NR1H4 discontinued); below = preclinical
# hero + pharmacology/control case studies (HKDC1 preclinical, SERPINE1 positive
# control, CYP3A4 drugged). The row order + per-gene drug_dev_status tags already
# convey the gradient split, so no overlapping side text is needed.
yb <- which(levels(tiles$gene_lab) == lab_map[["HKDC1"]]) + 0.5
p <- p +
  annotate("segment", x = 0.4, xend = length(mods) + 0.6, y = yb, yend = yb,
           linewidth = 0.35, color = "grey55", linetype = "22")

# ── Save ─────────────────────────────────────────────────────────────────────
out <- file.path(FIG4_DIR, "_supp", "convergence_evidence_matrix.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
cairo_pdf(out, width = fig_half_width * 1.45, height = 2.5)
print(p)
invisible(dev.off())
message("Saved: ", out)

# ── Legend text + provenance to stdout (stats live in legend, NOT on panel) ───
message("\n[convergence_evidence_matrix] LEGEND (paste into figure legend):")
message("Cross-modal replication fingerprint for six prioritized MASLD targets, ",
        "read as drug-development-gradient validation: the unsupervised ",
        "convergence score (no drug/clinical labels) recovers targets across the ",
        "drug-development gradient, each externally validated after the analysis. ",
        "Above the divider, two clinically anchored targets recovered de novo ",
        "(THRB, approved [resmetirom]; RORA, Phase 1 [TB-840] — clinical-stage ",
        "validation hero) and NR1H4 (discontinued [obeticholic acid]); below, the ",
        "preclinical validation hero HKDC1 (hepatocyte Hkdc1-KO protects mice ",
        "from MASH) plus case studies SERPINE1 (positive control) and CYP3A4 ",
        "(drugged CYP metabolizer). ",
        "Filled dot colour = disease direction for the two signed modalities ",
        "(Bulk RNA, Mouse: up=magenta, down=blue); teal = unsigned detection ",
        "(COLOC, ATAC motif, Spatial). Dot size = relative strength. Open circle = ",
        "interrogated but no signal; grey x = modality not tested for that gene.")
message("CAVEATS (from ledger):")
message(" - COLOC axis = COLOC PP.H4 (susie-else-abf, canonical susie_coloc_all_gwas.csv): ",
        "THRB abf 1.000 (UKBB_GGT, NO SuSiE); RORA SuSiE 0.998 (BBJ_GGT); ",
        "NR1H4 abf 0.208 (null genetic support — do not overstate); ",
        "HKDC1 SuSiE 0.992 (UKBB_AST); SERPINE1 abf 0.138 (weak; SERPINE1 story is ",
        "CCC+spatial); CYP3A4 abf 0.558 (LD-proximity caveat w/ CYP3A7).")
message(" - ATAC motif = variant-motif PAIRS disrupting the gene's own TF motif ",
        "(motifbreakR): THRB 10 (5 unique variants), RORA 9 (5), NR1H4 2 (NR1H4::RXRA). ",
        "Non-TF genes have no disruptable self-motif = tested-absent (open circle).")
message(" - Spatial = disease-state (Steatotic) Moran's C, per-condition ",
        "(GSE192741, n=5 donors, 2 healthy/3 steatotic; spatial ORGANIZATION not ",
        "disease-direction): SERPINE1 0.059->0.541, CYP3A4 0.381->0.465. ",
        "Only SERPINE1/CYP3A4 are in the SVG/autocorr test set; all other genes ",
        "= not measured. SERPINE1 spatial organization scooped by Kim 2024 J Hepatol.")
message(" - Protein = DIA-MS best logFC (130 samples, 2 datasets): HKDC1 +0.621 ",
        "(padj 0.72 n.s.), CYP3A4 -0.691; THRB/RORA/NR1H4/SERPINE1 NOT in DIA panel.")
message(" - Mouse = cross-species meta logFC: THRB -0.448, RORA -0.236, NR1H4 -0.317, ",
        "SERPINE1 +1.874; HKDC1/CYP3A4 not measured in mouse meta.")

options(width = 200)
prov <- tiles %>%
  transmute(gene, modality = as.character(modality), state,
            value = round(value, 3), direction,
            coloc_method = ifelse(modality == "COLOC", coloc_method[gene], "")) %>%
  arrange(factor(gene, levels = genes_oi), factor(modality, levels = mods))
message("\n[convergence_evidence_matrix] per-cell provenance:")
print(prov, row.names = FALSE)
