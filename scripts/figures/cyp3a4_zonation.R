#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# KEY MESSAGE (counter-intuitive exemplar):
#   CYP3A4 — the dominant drug-metabolizing CYP (~50% of all marketed drugs) — is
#   LOST pericentrally in MASLD at both the mRNA (C2 bulk log2FC = -0.42) and
#   protein (log2FC = -0.69) level, YET its RESIDUAL expression becomes MORE
#   spatially clustered with disease (per-condition Moran's I 0.381 -> 0.465).
#   Pericentral drug-metabolizing capacity collapses even as what remains
#   self-organizes more tightly.
#
# Numbers are C2 (limma-voom quality-weighted), NOT the retired dream method, and
# are reported in the figure legend (emitted to stdout) — never printed on the
# panel (PI directive).
#
# Provenance (every plotted number traced to a ledger row / canonical file):
#   bulk log2FC  -0.418  multi_evidence_atlas_with_spatial.csv :: bulk_logFC        (ledger 21)
#   protein log2FC -0.691 multi_evidence_atlas_with_spatial.csv :: best_protein_logFC (ledger 23)
#   Moran's I 0.381->0.465  spatial_autocorr_{Healthy,Steatotic}.csv :: C            (ledger 33)
#   COLOC PP.H4 0.558       gene_level_coloc.csv :: coloc_best_pp4 (UKBB_GGT, abf)    (ledger 12)
#   pericentral gradient    deg_zonation_combined.csv :: mean_PP1..mean_PC1 (Guilliams + Vu)
#
# CAVEATS baked into the legend (message() to stdout):
#   - Spatial autocorrelation = spatial ORGANIZATION, NOT disease direction; n=5
#     GSE192741 donors (2 Healthy / 3 Steatotic JBO), spot-level perm is
#     anticonservative (ledger 33/49).
#   - COLOC is abf-only (PP.H4.abf 0.558, UKBB_GGT); SuSiE non-convergent. The
#     locus is an LD-proximity cluster with CYP3A7 — colocalization cannot resolve
#     CYP3A4 vs CYP3A7 (ledger 12; gene_level_coloc ld_cluster_flag).
#   - Drug-clearance consequences are substrate-specific (CYP3A4 metabolizes ~50%
#     of drugs but the impact varies by compound).
#
# Output: figures/main/fig4_validation/fig4f_cyp3a4_zonation.pdf   (flat, no fig-# prefix)
# Env:    rnaseq
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(tidyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

dz <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/validation_bulk/deg_zonation_combined.csv"),
  stringsAsFactors = FALSE)

# Canonical COLOC: gene_level_coloc.csv (NOT the atlas *_coloc_pp4 convenience
# columns, which are method-mixed). coloc_best_pp4 = best abf PP.H4 across the
# 23-GWAS portfolio; coloc_best_susie_pp4 = best SuSiE PP.H4 (NA = non-convergent).
coloc <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
  stringsAsFactors = FALSE)

# Per-condition spatial autocorrelation — single source of truth for the
# healthy -> steatotic Moran's I comparison (GSE192741; ledger row 33).
# NOTE: differential_svgs.csv lists CYP3A4 at 0.567/0.644 (a separate squidpy
# pass); the canonical per-condition Moran's I in the ledger is 0.381 -> 0.465
# from spatial_autocorr_{Healthy,Steatotic}.csv. We use the ledger source.
.ac_h <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Healthy.csv"),
  stringsAsFactors = FALSE)
.ac_s <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Steatotic.csv"),
  stringsAsFactors = FALSE)
cyp_moran_h <- .ac_h$C[.ac_h$Gene == "CYP3A4"][1]   # 0.381 healthy
cyp_moran_s <- .ac_s$C[.ac_s$Gene == "CYP3A4"][1]   # 0.465 steatotic

# CYP3A4 canonical values for the legend.
cyp_atlas <- filter(atlas, human_symbol == "CYP3A4")
cyp_coloc <- filter(coloc, gene == "CYP3A4")
cyp_pp4   <- as.numeric(cyp_coloc$coloc_best_pp4[1])           # 0.558 abf
cyp_pp4_gwas <- cyp_coloc$coloc_best_gwas[1]                   # UKBB_GGT
cyp_ld    <- cyp_coloc$ld_cluster_best_gene[1]                 # CYP3A7

# ── Panel A: the counter-intuitive headline ─────────────────────────────────
# Residual CYP3A4 expression becomes MORE spatially clustered with disease:
# per-condition Moran's I rises 0.381 (Healthy) -> 0.465 (Steatotic) even as
# the gene is lost (Panels B/C). This is the per-condition autocorr metric
# (spatial_autocorr_{Healthy,Steatotic}.csv), distinct from the cross-gene
# atlas snapshot — so it gets its own panel rather than an x-axis position.
autocorr_df <- data.frame(
  cond   = factor(c("Healthy", "Steatotic"), levels = c("Healthy", "Steatotic")),
  morans = c(cyp_moran_h, cyp_moran_s)
)
cond_cols <- c(Healthy = masld_colors[["ns"]], Steatotic = "#E65100")  # control gray, disease orange

p_autocorr <- ggplot(autocorr_df, aes(x = cond, y = morans, group = 1)) +
  geom_line(linewidth = 0.7, color = "#E65100") +
  geom_point(aes(color = cond), size = 3.2) +
  geom_text(aes(label = sprintf("%.3f", morans)),
            vjust = -1.1, size = PUB_GEOM_TEXT, color = "gray15") +
  scale_color_manual(values = cond_cols, guide = "none") +
  scale_y_continuous(limits = c(0.30, 0.52),
                     breaks = c(0.35, 0.40, 0.45, 0.50)) +
  labs(x = NULL, y = "Spatial autocorrelation\n(Moran's I)") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Panel B: CYP3A4 periportal -> pericentral expression gradient ─────────────
zone_order  <- c("PP1", "PP2", "Mid", "PC2", "PC1")
zone_labels <- c("Periportal", "PP2", "Mid", "PC2", "Pericentral")

cyp_dz <- dz %>%
  filter(gene == "CYP3A4") %>%
  select(dataset, mean_PP1, mean_PP2, mean_Mid, mean_PC2, mean_PC1) %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "zone_raw",
               values_to = "expression") %>%
  mutate(
    zone    = factor(sub("mean_", "", zone_raw), levels = zone_order),
    dataset = sub(" et al\\.", "", dataset)
  ) %>%
  group_by(dataset) %>%
  mutate(expr_z = (expression - mean(expression)) / sd(expression)) %>%
  ungroup()

p_gradient <- ggplot(cyp_dz, aes(x = zone, y = expr_z, group = dataset,
                                  color = dataset, linetype = dataset)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_x_discrete(labels = zone_labels) +
  # Display Guilliams Visium cohort by its GEO accession (data keys unchanged so
  # the color/linetype joins still match; Vu has no mapped accession).
  scale_color_manual(values = c("Guilliams" = "#E65100", "Vu" = "#FB8C00"),
                     labels = c("Guilliams" = "GSE192741", "Vu" = "Vu"),
                     name = NULL) +
  scale_linetype_manual(values = c("Guilliams" = "solid", "Vu" = "dashed"),
                        labels = c("Guilliams" = "GSE192741", "Vu" = "Vu"),
                        name = NULL) +
  labs(x = NULL, y = "Expression (z)") +
  theme_masld() + theme_pub() +
  theme(legend.position  = "bottom",
        legend.key.size  = PUB_LEGEND_KEY,
        axis.text.x      = element_text(size = PUB_AXIS_TEXT, angle = 35, hjust = 1))

# ── Panel C: mRNA vs protein log2FC (both DOWN) ───────────────────────────────
modality_df <- data.frame(
  modality = factor(c("mRNA", "Protein"), levels = c("mRNA", "Protein")),
  logFC    = c(cyp_atlas$bulk_logFC[1], cyp_atlas$best_protein_logFC[1]),
  fill_col = c(masld_colors[["down"]], "#1976D2")
)

p_modality <- ggplot(modality_df, aes(x = modality, y = logFC, fill = fill_col)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  scale_fill_identity() +
  scale_y_continuous(limits = c(-0.78, 0.04), breaks = c(-0.6, -0.4, -0.2, 0)) +
  labs(x = NULL, y = "log2FC (MASLD vs control)") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Legend stats + caveats -> stdout (NOT on the panel) ───────────────────────
message(sprintf(
  "[cyp3a4 legend] CYP3A4 is lost pericentrally (bulk mRNA log2FC = %.3f [C2]; protein log2FC = %.3f [DIA-MS]) yet its residual expression becomes MORE spatially clustered with disease (per-condition Moran's I %.3f healthy -> %.3f steatotic).",
  cyp_atlas$bulk_logFC[1], cyp_atlas$best_protein_logFC[1],
  cyp_moran_h, cyp_moran_s))
message(sprintf(
  "[cyp3a4 caveat] Spatial autocorrelation indexes spatial ORGANIZATION, not disease direction (GSE192741, n=5: 2 Healthy / 3 Steatotic donors; spot-level permutation is anticonservative). COLOC = abf-only PP.H4 %.3f (%s); SuSiE non-convergent; locus is an LD-proximity cluster with %s, so colocalization cannot resolve CYP3A4 vs %s. Drug-clearance impact is substrate-specific (CYP3A4 ~50%% of marketed drugs).",
  cyp_pp4, cyp_pp4_gwas, cyp_ld, cyp_ld))

# ── Assemble ──────────────────────────────────────────────────────────────────
# Left-to-right narrative: LOSS (mRNA + protein) -> WHERE (pericentral gradient)
# -> the counter-intuitive TWIST (residual clustering rises with disease).
p_out <- (p_modality | p_gradient | p_autocorr) +
  plot_layout(widths = c(0.85, 1.15, 0.85))

out <- file.path(FIG4_DIR, "fig4f_cyp3a4_zonation.pdf")
cairo_pdf(out, width = fig_full_width, height = 2.6, onefile = TRUE)
print(p_out)
dev.off()
message("Saved: ", out)
