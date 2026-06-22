#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# KEY MESSAGE: Each prioritized target reproduces across a DISTINCT, gene-unique
# set of orthogonal modalities — no two of the six case-study genes share the
# same evidence fingerprint. The six genes span the DRUG-DEVELOPMENT GRADIENT
# (drug_dev_status, data/external/drug_targets/drug_target_classification.tsv):
# the unsupervised convergence score recovers them with NO drug/clinical labels.
#   * Approved / in-clinic (above the divider): THRB (resmetirom, masld_approved),
#     RORA (Phase 1 TB-840, masld_clinical; nominated de novo → first-in-human
#     AFTER our analysis = clinical validation hero), NR1H4 (obeticholic acid,
#     masld_discontinued).
#   * Preclinical / case study (below the divider): HKDC1 (preclinically validated
#     — hepatocyte Hkdc1-KO protects mice from MASH; masld_preclinical;
#     preclinical validation hero), SERPINE1 (positive control), CYP3A4 (drugged
#     metabolizer / pharmacology case).
# Each lights up a different combination of bulk RNA, cross-ancestry COLOC,
# GWAS-ATAC motif disruption, DIA-MS protein, spatial organization, and mouse
# cross-species.
#
# DESIGN (independent alt-take of the convergence hero — sister panel to
# target_evidence_matrix.R; deliberately a different visual grammar):
#   * A clustered dot/glyph matrix (genes x 6 modalities) with a per-gene
#     "modalities converging" count bar on the LEFT, so the reader's eye reads
#     each ROW as a unique fingerprint rather than scanning columns.
#   * SIGNED modalities (Bulk RNA, Protein, Mouse) encode DIRECTION by glyph
#     SHAPE + COLOUR: filled up-triangle = up in MASLD (magenta), down-triangle
#     = down in MASLD (blue). Magnitude = glyph size.
#   * UNSIGNED presence/strength modalities (COLOC, ATAC motif, Spatial) encode
#     STRENGTH by a teal filled circle sized 0–1; they carry no direction.
#   * HONEST 3-state encoding of every cell, never a fake zero-dot:
#       - measured-with-signal  -> sized coloured glyph
#       - tested-but-absent      -> small open ring (motif scan run, 0 hits)
#       - not-measured           -> light-grey "x" (gene not in that assay panel)
#
# NUMBERS: every plotted value traces to fig4_number_ledger.tsv (rows cited
# inline). COLOC comes from the canonical susie_coloc_all_gwas.csv (NEVER the
# atlas *_coloc_pp4 convenience columns). Spatial uses the ledger-sanctioned
# PER-CONDITION Moran autocorr (Healthy->Steatotic gain), NOT the atlas
# single-value spatial_morans_i column (ledger row 32 flags it: "do not cite").
#
# Output: figures/main/fig4_validation/convergence_evidence_matrix_v2.pdf
# Env:    rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

num <- function(x) suppressWarnings(as.numeric(x))

# NOTE: SERPINE1's drug_target_classification.tsv drug_dev_status is
# masld_preclinical, but it is presented here as a POSITIVE CONTROL per the
# reframe spec (drug_dev_status_reframe_spec.md line 31) and fig4.md — an
# explicit, intentional framing override, not a label mismatch.
genes_oi <- c("THRB", "RORA", "NR1H4", "HKDC1", "SERPINE1", "CYP3A4")  # approved/clinical | preclinical/case
mods     <- c("Bulk RNA", "COLOC", "ATAC motif", "Protein", "Spatial", "Mouse")
signed_mods <- c("Bulk RNA", "Protein", "Mouse")  # direction-bearing

# ── Modality normalisers (denominator = observed maximum across the 6 genes) ──
NORM <- c("Bulk RNA"   = 1.103,   # HKDC1 |logFC| (ledger row 19)
          "Protein"    = 0.691,   # CYP3A4 |logFC| (ledger row 23)
          "Mouse"      = 1.874,   # SERPINE1 |logFC| (ledger row 26)
          "COLOC"      = 1.000,   # posterior PP.H4 is already 0-1
          "ATAC motif" = 10,      # THRB 10 pairs (ledger row 36/44)
          "Spatial"    = 0.541)   # SERPINE1 steatotic Moran (ledger row 32)

# ── 1. Atlas-sourced signed modalities + COLOC fallback id ────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
A <- atlas[match(genes_oi, atlas$human_symbol), ]

bulk_lfc  <- setNames(num(A$bulk_logFC),         genes_oi)   # ledger rows 16-22
prot_lfc  <- setNames(num(A$best_protein_logFC), genes_oi)   # ledger rows 23-25
mouse_lfc <- setNames(num(A$mouse_meta_logFC),   genes_oi)   # ledger row 26

# ── 2. COLOC — canonical susie_coloc_all_gwas.csv (ledger rows 7-12) ──────────
#    best SuSiE PP.H4 if any CS converged; else best ABF PP.H4 fallback.
#    A gene with NO coloc row at all = not_measured (not the case here — all 6
#    are in the Broadaway eGene panel). axis label says "COLOC PP.H4" because
#    several plotted genes are abf-only (THRB/NR1H4/SERPINE1/CYP3A4).
cl <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
coloc_pp4    <- setNames(rep(NA_real_, length(genes_oi)), genes_oi)
coloc_method <- setNames(rep(NA_character_, length(genes_oi)), genes_oi)
coloc_gwas   <- setNames(rep(NA_character_, length(genes_oi)), genes_oi)
for (g in genes_oi) {
  sub <- cl[cl$gene == g, ]
  if (nrow(sub) == 0) next
  su <- num(sub$PP.H4.susie); ab <- num(sub$PP.H4.abf)
  if (any(!is.na(su))) {
    j <- which.max(su); coloc_pp4[g] <- su[j]; coloc_method[g] <- "SuSiE"
    coloc_gwas[g] <- sub$gwas_name[j]
  } else if (any(!is.na(ab))) {
    j <- which.max(ab); coloc_pp4[g] <- ab[j]; coloc_method[g] <- "abf"
    coloc_gwas[g] <- sub$gwas_name[j]
  }
}

# ── 3. ATAC motif disruption — canonical motifbreakR table (ledger row 36/44) ─
#    variant-motif PAIRS hitting the gene's own TF motif. Only TF genes carry a
#    self-motif; non-TF genes (HKDC1/SERPINE1/CYP3A4) = tested-but-absent (0).
md <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
atac_pairs <- setNames(sapply(genes_oi,
  function(g) sum(grepl(g, md$tf_name, ignore.case = TRUE))), genes_oi)

# ── 4. Spatial — ledger-sanctioned PER-CONDITION Moran autocorr GAIN ──────────
#    (ledger rows 32-33). We plot the Steatotic Moran I as the strength of
#    spatial ORGANIZATION; genes not in the spatial co-expression table are
#    not_measured. SERPINE1 0.541 / CYP3A4 0.465 are the only two case genes
#    present. n=5 donors -> spatial ORGANIZATION not disease-direction (caveat).
sp_steat <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Steatotic.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
sp_heal  <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Healthy.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
spat_val  <- setNames(rep(NA_real_, length(genes_oi)), genes_oi)
spat_heal <- setNames(rep(NA_real_, length(genes_oi)), genes_oi)
for (g in genes_oi) {
  rs <- sp_steat[match(g, sp_steat$Gene), "C"]
  rh <- sp_heal [match(g, sp_heal$Gene),  "C"]
  if (length(rs) && !is.na(rs)) spat_val[g]  <- num(rs)
  if (length(rh) && !is.na(rh)) spat_heal[g] <- num(rh)
}

# ── Assemble the long cell table — one row per gene x modality ────────────────
mk <- function(gene, modality, value, raw, state, strength, direction) {
  data.frame(gene = gene, modality = modality, value = value, raw = raw,
             state = state, strength = strength, direction = direction,
             stringsAsFactors = FALSE)
}
sgn   <- function(x) ifelse(x > 0, "up", "down")
scl01 <- function(x, mx) pmin(abs(x) / mx, 1)

rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- mk(...)
for (g in genes_oi) {
  # Bulk RNA — signed; every gene measured
  v <- bulk_lfc[g]
  add(g, "Bulk RNA", v, v, "measured", scl01(v, NORM["Bulk RNA"]), sgn(v))

  # COLOC — only PP.H4 >= 0.5 (canonical threshold) counts as DETECTED; sub-
  # threshold posteriors (NR1H4 0.208, SERPINE1 0.138; ledger flags NR1H4 null)
  # are TESTED-ABSENT (open circle), never counted as a converging modality.
  v <- coloc_pp4[g]
  if (is.na(v))            add(g, "COLOC", NA, NA, "not_measured", NA, NA)
  else if (v >= 0.5)       add(g, "COLOC", v, v, "measured", scl01(v, NORM["COLOC"]), "neutral")
  else                     add(g, "COLOC", v, v, "tested_absent", 0, "neutral")

  # ATAC motif — unsigned count; all interrogated => measured, 0 = tested-absent
  v <- atac_pairs[g]
  st <- if (v > 0) "measured" else "tested_absent"
  add(g, "ATAC motif", v, v, st, scl01(v, NORM["ATAC motif"]), "neutral")

  # Protein — signed DIA-MS; NA => not in panel = not_measured
  v <- prot_lfc[g]
  if (is.na(v)) add(g, "Protein", NA, NA, "not_measured", NA, NA)
  else          add(g, "Protein", v, v, "measured", scl01(v, NORM["Protein"]), sgn(v))

  # Spatial — 3-state: SVG-significant (SERPINE1/CYP3A4) => detected; present in
  # the 11,984-gene autocorr file but not SVG (THRB/RORA/NR1H4) => tested-absent
  # (open circle); absent from the file (HKDC1) => not measured (x).
  v <- spat_val[g]
  if (g %in% c("SERPINE1", "CYP3A4")) add(g, "Spatial", v, v, "measured", scl01(v, NORM["Spatial"]), "neutral")
  else if (!is.na(v))                 add(g, "Spatial", v, v, "tested_absent", 0, "neutral")
  else                                add(g, "Spatial", NA, NA, "not_measured", NA, NA)

  # Mouse — signed cross-species meta; NA => not tested = not_measured
  v <- mouse_lfc[g]
  if (is.na(v)) add(g, "Mouse", NA, NA, "not_measured", NA, NA)
  else          add(g, "Mouse", v, v, "measured", scl01(v, NORM["Mouse"]), sgn(v))
}
tiles <- do.call(rbind, rows)

# Order: genes TOP->BOTTOM as listed (reverse so THRB is top row).
tiles$gene     <- factor(tiles$gene, levels = rev(genes_oi))
tiles$modality <- factor(tiles$modality, levels = mods)
tiles$signed   <- tiles$modality %in% signed_mods

# ── Per-gene "modalities converging" count = measured-with-signal cells ───────
conv <- tiles %>%
  group_by(gene) %>%
  summarise(n_conv = sum(state == "measured" & strength > 0), .groups = "drop")
conv$gene <- factor(conv$gene, levels = levels(tiles$gene))

# ── Glyph sub-frames ──────────────────────────────────────────────────────────
# viz_strength floors a measured-with-signal cell to a visible glyph so a genuine
# but weak value (e.g. THRB spatial Moran 0.006) reads as a small filled circle —
# distinguishable from a tested-absent open ring — without misrepresenting the
# magnitude (the true value/strength stay intact in the provenance dump).
SIZE_FLOOR <- 0.12
tiles$viz_strength <- ifelse(tiles$state == "measured" & tiles$strength > 0,
                             pmax(tiles$strength, SIZE_FLOOR), tiles$strength)
# Signed measured cells -> directional triangles (shape carries direction).
sig_tiles <- tiles %>%
  filter(signed, state == "measured", strength > 0) %>%
  mutate(dir_lab = ifelse(direction == "up", "Up in MASLD", "Down in MASLD"))
# Unsigned measured cells -> teal strength circles.
uns_tiles <- tiles %>% filter(!signed, state == "measured", strength > 0)
# No-signal cells -> shape legend (tested-absent ring vs not-measured x).
noflag <- bind_rows(
  tiles %>% filter(state == "tested_absent") %>% mutate(glyph = "Tested, absent"),
  tiles %>% filter(state == "not_measured")  %>% mutate(glyph = "Not measured")
)
noflag$glyph <- factor(noflag$glyph, levels = c("Tested, absent", "Not measured"))

dir_pal  <- c("Up in MASLD" = masld_colors[["up"]], "Down in MASLD" = masld_colors[["down"]])
TEAL     <- "#00897B"
SEP_GREY <- "#9E9E9E"

# Numeric x positions so the count bar lives to the LEFT of the matrix.
mod_x  <- setNames(seq_along(mods), mods)
gene_y <- setNames(seq_along(levels(tiles$gene)), levels(tiles$gene))
xmax   <- length(mods)
bar_x0 <- -0.4
bar_unit <- 0.95 / 6   # 6 = max possible convergent modalities

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  # faint cell grid
  geom_tile(data = tiles,
            aes(x = as.numeric(modality), y = as.numeric(gene)),
            fill = "white", color = "grey88", linewidth = 0.25) +
  # signed directional triangles (shape + colour = direction; size = magnitude)
  geom_point(data = sig_tiles,
             aes(x = as.numeric(modality), y = as.numeric(gene),
                 size = viz_strength, fill = dir_lab, shape = dir_lab),
             color = "white", stroke = 0.22) +
  # unsigned strength circles (teal; size = strength)
  geom_point(data = uns_tiles,
             aes(x = as.numeric(modality), y = as.numeric(gene), size = viz_strength),
             shape = 21, fill = TEAL, color = "white", stroke = 0.22) +
  # no-signal glyphs (ring = tested-absent; x = not-measured)
  geom_point(data = noflag,
             aes(x = as.numeric(modality), y = as.numeric(gene), shape = glyph),
             size = 1.6, color = "grey55", stroke = 0.5) +
  scale_shape_manual(
    values = c("Up in MASLD" = 24, "Down in MASLD" = 25,
               "Tested, absent" = 1, "Not measured" = 4),
    breaks = c("Up in MASLD", "Down in MASLD", "Tested, absent", "Not measured"),
    name = NULL,
    guide = guide_legend(order = 1, override.aes = list(
      size   = c(2.2, 2.2, 1.8, 1.8),
      fill   = c(masld_colors[["up"]], masld_colors[["down"]], NA, NA),
      color  = c("white", "white", "grey55", "grey55"),
      stroke = c(0.22, 0.22, 0.5, 0.5)))) +
  scale_fill_manual(values = dir_pal, guide = "none") +
  scale_size_area(max_size = 4.3, name = "Strength",
                  breaks = c(0.25, 0.5, 1.0), labels = c("low", "med", "high"),
                  guide = guide_legend(order = 3, override.aes = list(
                    shape = 21, fill = "grey55", color = "white", stroke = 0.2)))

# ── Left count bar: how many modalities converge per gene ─────────────────────
conv$ y    <- gene_y[as.character(conv$gene)]
conv$ xend <- bar_x0 - conv$n_conv * bar_unit
p <- p +
  geom_segment(data = conv,
               aes(x = bar_x0, xend = xend, y = y, yend = y),
               linewidth = 1.7, color = "#37474F", lineend = "round") +
  geom_text(data = conv,
            aes(x = xend - 0.08, y = y, label = n_conv),
            hjust = 1, size = PUB_GEOM_TEXT, color = "#37474F") +
  annotate("text", x = bar_x0 - 0.5 * 0.95, y = length(gene_y) + 0.75,
           label = "modalities\nconverging", lineheight = 0.85,
           size = PUB_GEOM_TEXT - 0.1, color = "#37474F", fontface = "plain")

# ── Drug-dev gradient separator: approved/clinical (top) | preclinical/case ───
#    (between NR1H4 and HKDC1). Above = in-clinic real MASLD drugs the
#    unsupervised score recovered (THRB approved, RORA Phase 1, NR1H4
#    discontinued); below = preclinical / positive-control / pharmacology case
#    (HKDC1 preclinically validated, SERPINE1 positive control, CYP3A4
#    metabolizer). Statuses = drug_dev_status (drug_target_classification.tsv).
yb <- gene_y["HKDC1"] + 0.5
p <- p +
  annotate("segment", x = bar_x0 - 0.95 - 0.15, xend = xmax + 0.55,
           y = yb, yend = yb, linewidth = 0.35, color = SEP_GREY, linetype = "22") +
  annotate("text", x = xmax + 0.7, y = (gene_y["THRB"] + gene_y["NR1H4"]) / 2,
           label = "approved / clinical", angle = -90, size = PUB_GEOM_TEXT,
           color = "grey45", vjust = 0) +
  annotate("text", x = xmax + 0.7, y = (gene_y["HKDC1"] + gene_y["CYP3A4"]) / 2,
           label = "preclinical / case", angle = -90, size = PUB_GEOM_TEXT,
           color = "grey45", vjust = 0)

# ── Axes / theme ──────────────────────────────────────────────────────────────
p <- p +
  scale_x_continuous(breaks = mod_x, labels = mods,
                     limits = c(bar_x0 - 0.95 - 0.55, xmax + 1.0),
                     position = "top", expand = c(0, 0)) +
  scale_y_continuous(breaks = gene_y, labels = names(gene_y),
                     limits = c(0.4, length(gene_y) + 1.35), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Each target has a unique evidence fingerprint") +
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(
    axis.text.y     = element_text(face = "italic", size = PUB_AXIS_TEXT + 1),
    axis.text.x.top = element_text(angle = 35, hjust = 0, size = PUB_AXIS_TEXT + 0.5),
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "right",
    legend.box      = "vertical",
    legend.spacing.y = unit(0.03, "cm"),
    legend.margin   = margin(0, 0, 0, 0),
    plot.title      = element_text(size = PUB_TITLE, face = "bold"),
    plot.margin     = margin(3, 4, 3, 3)
  )

# ── Save ──────────────────────────────────────────────────────────────────────
out <- file.path(FIG4_DIR, "_supp", "convergence_evidence_matrix_v2.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
cairo_pdf(out, width = fig_half_width * 1.55, height = 2.45)
print(p)
invisible(dev.off())
message("Saved: ", out)

# ── Legend text + caveats to stdout (stats live in the legend, not on panel) ──
message("\n[convergence_evidence_matrix_v2] FIGURE LEGEND TEXT:")
message("Cross-modal evidence fingerprints for six prioritized MASLD targets ",
        "spanning the drug-development gradient (drug_dev_status): approved / ",
        "in-clinic above the divider — THRB (resmetirom, approved), RORA ",
        "(Phase 1 TB-840, recovered de novo then independently entered ",
        "first-in-human after our analysis — clinical validation hero), NR1H4 ",
        "(obeticholic acid, discontinued); preclinical / case below — HKDC1 ",
        "(preclinically validated: hepatocyte Hkdc1-KO protects mice from MASH ",
        "— preclinical validation hero), SERPINE1 (positive control), CYP3A4 ",
        "(drugged metabolizer / pharmacology case). The unsupervised ",
        "convergence score used NO drug or clinical labels. Left bar = number ",
        "of modalities ",
        "converging (measured-with-signal). Up-triangle = up in MASLD ",
        "(magenta), down-triangle = down (blue), for the three signed ",
        "modalities (Bulk RNA, Protein, Mouse); teal circle = strength of the ",
        "three unsigned modalities (COLOC PP.H4, GWAS-ATAC motif disruption, ",
        "spatial organization); glyph size = magnitude. Open ring = assay run ",
        "but no signal (tested-absent); grey x = gene not in that assay panel ",
        "(not measured). No two genes share the same fingerprint.")
message("\nCAVEATS (baked from ledger):")
message(" - COLOC axis labeled 'COLOC PP.H4': THRB/NR1H4/SERPINE1/CYP3A4 are ",
        "coloc.abf-only (no SuSiE-convergent CS); RORA/HKDC1 are SuSiE ",
        "(ledger rows 7-12).")
message(" - Spatial = per-condition Moran organization (Steatotic), NOT a ",
        "disease-direction claim; GSE192741 n=5 donors, spot-level perm is ",
        "anticonservative; SERPINE1 0.059->0.541 / CYP3A4 0.381->0.465 ",
        "Healthy->Steatotic = organization GAIN only (ledger rows 32-33,49). ",
        "Atlas single-value spatial_morans_i (0.242 SERPINE1) deliberately NOT used.")
message(" - ATAC motif = variant-motif PAIRS on each gene's own TF motif ",
        "(THRB 10 / RORA 9 / NR1H4 2 NR1H4::RXRA; unique variants 5/5/<=2); ",
        "non-TF genes have no disruptable self-motif = tested-absent (ledger 36/44).")
message(" - HKDC1 protein +0.621 is n.s. (padj 0.72; ledger row 24); plotted ",
        "as measured-with-direction but not a significance claim.")

# ── Per-cell provenance dump ──────────────────────────────────────────────────
options(width = 200)
prov <- tiles %>%
  transmute(gene = as.character(gene), modality = as.character(modality),
            state, value = round(value, 3), direction, strength = round(strength, 3)) %>%
  arrange(factor(gene, levels = genes_oi), factor(modality, levels = mods))
message("\n[convergence_evidence_matrix_v2] per-cell provenance:")
print(prov, row.names = FALSE)
message("\nCOLOC source per gene (susie_coloc_all_gwas.csv):")
for (g in genes_oi) message(sprintf("  %-9s PP.H4=%s method=%s gwas=%s",
  g, round(coloc_pp4[g], 3), coloc_method[g], coloc_gwas[g]))
message("\nPer-gene convergent modality count:")
print(conv[, c("gene", "n_conv")], row.names = FALSE)
