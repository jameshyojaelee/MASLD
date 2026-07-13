#!/usr/bin/env Rscript
# KEY MESSAGE: Each of the six prioritized case-study genes — spanning the
# drug-development gradient that the unsupervised convergence score recovered
# WITHOUT any drug/clinical labels — reproduces across a DISTINCT set of
# orthogonal modalities. No two genes share the same evidence fingerprint.
# Drug-development status (data/external/drug_targets/drug_target_classification.tsv):
#   APPROVED / CLINICAL / DISCONTINUED (reached-clinic tier) above the divider:
#     THRB    — Approved (masld_approved; resmetirom)
#     RORA    — Clinical, Phase 1 (masld_clinical; TB-840, RORalpha agonist) — validation hero
#     NR1H4   — Discontinued (masld_discontinued; OCA/obeticholic acid withdrawn Sep-2025)
#   PRECLINICAL / CASE below the divider:
#     HKDC1   — Preclinically validated (masld_preclinical; hepatocyte Hkdc1-KO protects mice from MASH) — validation hero
#     SERPINE1 — Positive control (raw status masld_preclinical; framed as positive
#                control per spec line 31 — SERPINE1 inhibitors reduce NASH, Kim 2024)
#     CYP3A4  — Drugged metabolizer (drugged_other_indication); pharmacology case study
# The point: the score reached real targets across this entire gradient
# independently, and RORA/HKDC1 were externally validated AFTER the analysis —
# drug-development-gradient validation, not a bare novelty claim.
#
# A dot-matrix (genes x modalities). DIRECTION is encoded by colour (disease-up
# magenta / down blue) for the two signed modalities (Bulk RNA, Mouse); the
# unsigned presence-or-strength modalities (COLOC, ATAC motif, Spatial) use a
# neutral->strong teal ramp. STRENGTH is encoded by dot size. Cells that were
# NOT MEASURED (e.g. protein for the four genes absent from the DIA-MS panel)
# are rendered as a light grey "x", NEVER as a zero/minimum-strength dot —
# that would falsely imply measured-but-absent.
#
# Output: figures/main/fig4_validation/target_evidence_matrix.pdf
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

# ── Data ──────────────────────────────────────────────────────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)

genes_oi <- c("THRB", "RORA", "NR1H4", "HKDC1", "SERPINE1", "CYP3A4")  # row order
mods     <- c("Bulk RNA", "COLOC", "ATAC motif", "Protein", "Spatial", "Mouse")  # col order

A <- atlas[match(genes_oi, atlas$human_symbol), ]

# ── ATAC motif disruption (variant-motif pairs per gene's own TF motif) ───────
# Sourced honestly from the GWAS-ATAC motifbreakR output (Scripts 55-57). For
# these case-study genes the spatial-atlas gwas_motif_disrupted column is empty
# (the gwas-ATAC layer was not joined onto these rows), so we read the canonical
# per-TF disruption table directly — the same source as fig4c_nr_triad.R. Only
# the TF genes (THRB/RORA/NR1H4) have a self-motif that GWAS variants disrupt;
# the non-TF genes (HKDC1/SERPINE1/CYP3A4) have no disruptable motif = absent.
motif_f <- file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv")
atac_pairs <- setNames(rep(0L, length(genes_oi)), genes_oi)
if (file.exists(motif_f)) {
  md <- read.csv(motif_f, stringsAsFactors = FALSE, check.names = FALSE)
  for (g in genes_oi) {
    atac_pairs[g] <- sum(grepl(g, md$tf_name, ignore.case = TRUE))
  }
}
# THRB=10, RORA=9, NR1H4=2 (NR1H4::RXRA heterodimer), HKDC1/SERPINE1/CYP3A4=0.

# ── Build the long tile table: one row per gene x modality ───────────────────
# For each cell: state in {measured, not_measured}; value = signed/strength;
# strength = scaled 0-1 magnitude used for dot size; direction = up/down/neutral.

num <- function(x) suppressWarnings(as.numeric(x))
is_true <- function(x) !is.na(x) & (x %in% c(TRUE, "True", "TRUE", "true"))

mk <- function(gene, modality, value, state, strength, direction) {
  data.frame(gene = gene, modality = modality, value = value,
             state = state, strength = strength, direction = direction,
             stringsAsFactors = FALSE)
}

rows <- list()
for (i in seq_along(genes_oi)) {
  g  <- genes_oi[i]
  r  <- A[i, ]

  # 1. Bulk RNA — signed; every gene measured. Strength = |logFC| scaled.
  blfc <- num(r$bulk_logFC)
  rows[[length(rows)+1]] <- mk(g, "Bulk RNA", blfc, "measured",
                               min(abs(blfc) / 1.10, 1),                # 1.10 = HKDC1 max
                               ifelse(blfc > 0, "up", "down"))

  # 2. COLOC — unsigned posterior PP.H4 (SuSiE; ABF fallback). All measured.
  pp4 <- num(r$coloc_susie_best_pp4)
  if (is.na(pp4)) pp4 <- num(r$coloc_abf_best_pp4)
  rows[[length(rows)+1]] <- mk(g, "COLOC", pp4, "measured", pp4, "neutral")

  # 3. ATAC motif — unsigned count of disrupted variant-motif pairs.
  #    All genes interrogated (in the GWAS-ATAC scan) => measured; 0 = absent.
  np <- atac_pairs[[g]]
  rows[[length(rows)+1]] <- mk(g, "ATAC motif", np, "measured",
                               min(np / 10, 1), "neutral")             # 10 = THRB max

  # 4. Protein — signed DIA-MS logFC. NA => NOT MEASURED (gene absent from panel).
  plfc <- num(r$best_protein_logFC)
  if (is.na(plfc)) {
    rows[[length(rows)+1]] <- mk(g, "Protein", NA, "not_measured", NA, NA)
  } else {
    rows[[length(rows)+1]] <- mk(g, "Protein", plfc, "measured",
                                 min(abs(plfc) / 0.70, 1),             # ~0.69 max
                                 ifelse(plfc > 0, "up", "down"))
  }

  # 5. Spatial — unsigned Moran's I. The spatial-variability test was run on the
  #    top-3000 HVG set only; a gene with NA Moran's I was NOT in that test set,
  #    so it is NOT MEASURED for spatial variability (x), NOT a tested-absent (o).
  #    Only THRB/RORA/NR1H4/HKDC1 here (the SVG-tested case genes carry a value).
  mi  <- num(r$spatial_morans_i)
  if (is.na(mi)) {
    rows[[length(rows)+1]] <- mk(g, "Spatial", NA, "not_measured", NA, NA)
  } else {
    rows[[length(rows)+1]] <- mk(g, "Spatial", mi, "measured",
                                 min(mi / 0.65, 1), "neutral")         # ~0.644 max
  }

  # 6. Mouse — signed cross-species meta logFC. NA logFC => not tested in the
  #    mouse meta-analysis = NOT MEASURED (distinct from a measured null).
  mlfc <- num(r$mouse_meta_logFC)
  if (is.na(mlfc)) {
    rows[[length(rows)+1]] <- mk(g, "Mouse", NA, "not_measured", NA, NA)
  } else {
    rows[[length(rows)+1]] <- mk(g, "Mouse", mlfc, "measured",
                                 min(abs(mlfc) / 1.87, 1),             # 1.87 = SERPINE1 max
                                 ifelse(mlfc > 0, "up", "down"))
  }
}
tiles <- do.call(rbind, rows)

# Factor ordering: genes top->bottom as given (reverse so THRB is the TOP row),
# modalities left->right as given.
tiles$gene     <- factor(tiles$gene, levels = rev(genes_oi))
tiles$modality <- factor(tiles$modality, levels = mods)

# ── Colour mapping ───────────────────────────────────────────────────────────
# Signed modalities use semantic up/down; unsigned use a neutral->strong teal.
# We map a single fill aesthetic per row by computing a colour string up front
# so the legend can carry just the qualitative direction key.
tiles <- tiles %>%
  mutate(
    measured = state == "measured",
    fill_class = case_when(
      !measured                ~ "Not measured",
      direction == "up"        ~ "Up in MASLD",
      direction == "down"      ~ "Down in MASLD",
      TRUE                     ~ "Detected"
    )
  )

# Colour key = the three meanings a FILLED dot can carry. The two no-signal
# states (tested-but-absent, not-measured) are encoded by SHAPE, not colour,
# and get their own legend so every glyph in the matrix is explained.
fill_pal <- c(
  "Up in MASLD"   = masld_colors[["up"]],     # deep magenta (signed)
  "Down in MASLD" = masld_colors[["down"]],   # deep blue (signed)
  "Detected"      = "#00897B"                 # teal (unsigned presence/strength)
)
tiles$fill_class <- factor(tiles$fill_class,
                           levels = c(names(fill_pal), "Not measured"))

measured_tiles <- tiles %>% filter(measured, strength > 0)
# Two no-signal states share one shape legend: tested-but-absent (open circle)
# vs not-measured (x). Keep them in one frame so the legend is complete.
noflag <- bind_rows(
  tiles %>% filter(measured, strength == 0) %>% mutate(glyph = "Tested, absent"),
  tiles %>% filter(!measured)               %>% mutate(glyph = "Not measured")
)
noflag$glyph <- factor(noflag$glyph, levels = c("Tested, absent", "Not measured"))

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(tiles, aes(x = modality, y = gene)) +
  # faint cell grid so absent cells read as "tested but empty"
  geom_tile(fill = "white", color = "grey88", linewidth = 0.25) +
  # measured-with-signal: sized + coloured filled dots
  geom_point(data = measured_tiles,
             aes(size = strength, fill = fill_class),
             shape = 21, color = "white", stroke = 0.25) +
  # the two no-signal states, shape-coded (open circle = tested-absent; x = not measured)
  geom_point(data = noflag,
             aes(x = modality, y = gene, shape = glyph),
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
                       override.aes = list(size = 2, color = "grey45",
                                           stroke = 0.6))) +
  scale_size_area(max_size = 4.2, name = "Strength",
                  breaks = c(0.25, 0.5, 1.0),
                  labels = c("low", "med", "high"),
                  guide = guide_legend(order = 3,
                    override.aes = list(shape = 21, fill = "grey55",
                                        color = "white", stroke = 0.2))) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 0.85, clip = "off") +
  labs(x = NULL, y = NULL, title = "Cross-modal evidence") +
  theme_masld() + theme_pub() +
  theme(
    axis.text.y      = element_text(face = "italic", size = PUB_AXIS_TEXT + 1),
    axis.text.x.top  = element_text(angle = 35, hjust = 0, size = PUB_AXIS_TEXT + 0.5),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.04, "cm"),
    legend.margin    = margin(0, 0, 0, 0),
    plot.title       = element_text(size = PUB_TITLE, face = "plain")
  )

# Separator line along the drug-development gradient: APPROVED/CLINICAL above
# (THRB approved, RORA Phase 1, NR1H4 discontinued) | PRECLINICAL/CASE below
# (HKDC1 preclinically validated, SERPINE1 positive control, CYP3A4 pharmacology case).
# y positions: rev order => NR1H4 (row 4 from top) | HKDC1 below it.
# In the reversed factor, the boundary sits between "NR1H4" and "HKDC1".
yb <- which(levels(tiles$gene) == "HKDC1") + 0.5
p <- p +
  annotate("segment", x = 0.4, xend = length(mods) + 0.6, y = yb, yend = yb,
           linewidth = 0.35, color = "grey55", linetype = "22")

# ── Save ─────────────────────────────────────────────────────────────────────
out <- file.path(FIG4_DIR, "_supp", "target_evidence_matrix.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
cairo_pdf(out, width = fig_half_width * 1.35, height = 2.4)
print(p)
invisible(dev.off())
message("Saved: ", out)

# Emit the exact per-cell values to stdout for the figure caption / provenance.
options(width = 200)
prov <- tiles %>%
  transmute(gene = as.character(gene), modality = as.character(modality),
            state, value = round(value, 3), direction) %>%
  arrange(factor(gene, levels = genes_oi), factor(modality, levels = mods))
message("\n[target_evidence_matrix] per-cell provenance:")
print(prov, row.names = FALSE)
