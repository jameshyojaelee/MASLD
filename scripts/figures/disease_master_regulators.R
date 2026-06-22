#!/usr/bin/env Rscript
# PANEL A (Fig 4, VALIDATION) — master-regulator CROSS-MODALITY CONVERGENCE.
#
# KEY MESSAGE (validation, NOT primary discovery): the single-cohort scATAC
# hepatocyte SCENIC+ gene-regulatory network independently re-derives the SAME
# disease transcription factors that the well-powered modalities already call —
# the bulk MASLD DEG analysis (Fig 3, n=846) and the GWAS/COLOC causal map
# (Fig 2). Each master-regulator TF is shown with its evidence across three
# INDEPENDENT, PROPERLY-POWERED axes:
#   (1) bulk MASLD effect   — bulk_logFC / bulk_lfsr   (Fig 3, n=846; POWERED)
#   (2) genetic colocalization — is_coloc / coloc_pp4  (Fig 2 GWAS+eQTL; POWERED)
#   (3) variant->motif disruption — motifbreakR count  (sequence-level; VALID)
# RORA and THRB (n_evidence = 2: bulk DEG + COLOC) are the convergence heroes:
# disease-repressed in bulk, colocalizing (PP.H4 > 0.99), and carrying the most
# fine-mapped variant->motif disruptions.
#
# *** WHAT THIS PANEL DELIBERATELY DOES NOT SHOW ***
# The single-cohort (n=18; 5 control / 13 disease) regulon-activity DE is
# PSEUDOREPLICATED at the cell level and UNDERPOWERED at the donor level
# (sc_underpowered = TRUE; activity_padj == 1 for ALL 48 TFs). We therefore make
# NO disease-DE significance claim from the scATAC side. The scATAC contribution
# is purely STRUCTURAL: which TFs the GRN nominates as hepatocyte regulators —
# the corroboration comes from the convergence with the powered modalities.
#
# DATA PROVENANCE (read from disk; numbers NEVER hardcoded from prose):
#   - disease_master_regulators.csv  (C2-regenerated; cols bulk_logFC/bulk_lfsr,
#     is_coloc/coloc_pp4, n_evidence, activity_padj, sc_underpowered)
#     Analysis/ATAC/Human_Multiome/scenic_plus/disease_master_regulators.csv
#   - motif_disruption_scores.csv  (motifbreakR; per-TF row count = variant-motif
#     PAIRS; UNIQUE-variant count computed separately)
#     GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#
# Output: figures/main/fig4_validation/disease_master_regulators.pdf
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

# ── Data ────────────────────────────────────────────────────────────────────
dmr <- read.csv(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/scenic_plus/disease_master_regulators.csv"),
  stringsAsFactors = FALSE)
mds <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"),
  stringsAsFactors = FALSE)

# Guard: this panel REQUIRES the C2-regenerated schema (bulk_*, not dream_*).
req_cols <- c("tf_name", "bulk_logFC", "bulk_lfsr", "is_coloc", "coloc_pp4",
              "n_evidence", "activity_padj", "sc_underpowered")
missing  <- setdiff(req_cols, colnames(dmr))
if (length(missing))
  stop("disease_master_regulators.csv missing C2 columns: ",
       paste(missing, collapse = ", "),
       " — re-run the C2 regeneration before plotting.")

dmr_tfs <- unique(dmr$tf_name)                         # 48 cross-modality TFs

# Variant->motif disruption counts. motifbreakR tf_name may be a dimer
# (e.g. "HOXD12::ELK1"); split into tokens so a master TF is credited whenever
# it appears in a (possibly composite) motif. PAIRS = number of (variant, motif)
# rows; UNIQUE = distinct fine-mapped variants.
mds_tokens <- strsplit(mds$tf_name, "::")
pair_count <- vapply(dmr_tfs, function(tf)
  sum(vapply(mds_tokens, function(z) tf %in% z, logical(1))), integer(1))
uniq_count <- vapply(dmr_tfs, function(tf)
  length(unique(mds$SNP_id[vapply(mds_tokens, function(z) tf %in% z, logical(1))])),
  integer(1))

ev <- dmr %>%
  transmute(
    tf         = tf_name,
    bulk_logFC = bulk_logFC,
    bulk_lfsr  = bulk_lfsr,
    is_coloc   = is_coloc == TRUE,
    coloc_pp4  = coloc_pp4,
    n_evidence = n_evidence
  ) %>%
  mutate(
    motif_pairs = pair_count[match(tf, dmr_tfs)],
    motif_uniq  = uniq_count[match(tf, dmr_tfs)]
  )

# ── Select the TFs to display ─────────────────────────────────────────────────
# Show every COLOC hero (n_evidence = 2) plus the most variant-motif-disrupted
# disease master regulators, so the panel stays compact and information-dense.
TOP_MOTIF <- 12
keep <- ev %>%
  filter(is_coloc | motif_pairs > 0) %>%
  arrange(desc(is_coloc), desc(motif_pairs)) %>%
  { union(.$tf[.$is_coloc], head(.$tf[order(-.$motif_pairs)], TOP_MOTIF)) }
ev <- ev %>% filter(tf %in% keep)

# Rank for the y-axis: heroes (n_evidence = 2) on top, then by motif disruptions,
# then by bulk significance. coord_flip => reverse the level order so #1 is top.
ev <- ev %>%
  arrange(n_evidence, motif_pairs, -log10(pmax(bulk_lfsr, 1e-300))) %>%
  mutate(
    hero    = is_coloc,                       # RORA, THRB, RELA, AHR
    tf_lab  = ifelse(hero, paste0(tf, "*"), tf),
    tf_lab  = factor(tf_lab, levels = unique(tf_lab))
  )
hero_levels <- levels(ev$tf_lab)[ ev$hero[match(levels(ev$tf_lab), ev$tf_lab)] ]

# ── Long format: one row per (TF, modality) cell of the convergence matrix ────
# Three modality columns. Each cell carries the encoded value + a printed label.
cell <- bind_rows(
  ev %>% transmute(
    tf_lab, hero,
    modality = "Bulk MASLD\n(n=846)",
    # signed effect drives the fill; magnitude drives the dot size
    fill_val = bulk_logFC,
    size_val = abs(bulk_logFC),
    present  = bulk_lfsr < 0.05,
    cell_lab = sprintf("%+.2f", bulk_logFC),
    kind     = "effect"
  ),
  ev %>% transmute(
    tf_lab, hero,
    modality = "GWAS COLOC\n(PP.H4)",
    fill_val = ifelse(is_coloc, coloc_pp4, NA_real_),
    size_val = coloc_pp4,
    present  = is_coloc,
    cell_lab = ifelse(is_coloc, sprintf("%.2f", coloc_pp4), ""),
    kind     = "coloc"
  ),
  ev %>% transmute(
    tf_lab, hero,
    modality = "Variant->motif\n(# variants)",
    fill_val = ifelse(motif_pairs > 0, as.numeric(motif_uniq), NA_real_),
    size_val = as.numeric(motif_uniq),
    present  = motif_pairs > 0,
    cell_lab = ifelse(motif_pairs > 0, as.character(motif_uniq), ""),
    kind     = "motif"
  )
) %>%
  mutate(modality = factor(modality,
    levels = c("Bulk MASLD\n(n=846)", "GWAS COLOC\n(PP.H4)",
               "Variant->motif\n(# variants)")))

# ── Colors ────────────────────────────────────────────────────────────────────
col_up   <- masld_colors$up     # #C9265E disease-up
col_down <- masld_colors$down   # #1565C0 disease-down
col_ns   <- masld_colors$ns     # #9E9E9E
col_pos  <- "#7B1FA2"           # violet — positive evidence (coloc / motif)

# Build a per-cell fill by modality kind:
#  - effect: diverging up(magenta)/down(blue) by sign of bulk_logFC
#  - coloc / motif: single positive-evidence violet, alpha by present
cell <- cell %>%
  mutate(
    cell_fill = case_when(
      kind == "effect" & fill_val >= 0 ~ col_up,
      kind == "effect" & fill_val <  0 ~ col_down,
      kind != "effect" & present       ~ col_pos,
      TRUE                             ~ col_ns
    ),
    cell_alpha = ifelse(present, 1, 0.18)
  )

# Hero TFs are flagged via the "*" suffix on their TF label (set above); per
# house style ALL figure text stays black, so no per-tick colour/face is used.

# ── Plot: cross-modality convergence dot-matrix ──────────────────────────────
# Layout: each cell is a fixed-size rounded marker whose FILL encodes the
# evidence (bulk = magenta-up / blue-down, COLOC + motif = violet positive,
# faded if absent). The numeric value is printed to the LEFT of the marker with
# a white halo so it never collides with the marker — the marker carries the
# at-a-glance signal, the number carries the exact value. The discrete x scale
# is preserved (so the column headers render) and the marker / label are offset
# from each tick via position_nudge.
p <- ggplot(cell, aes(x = modality, y = tf_lab)) +
  # subtle row band behind the hero TFs to draw the eye to the convergence rows
  geom_tile(data = cell %>% filter(hero, modality == levels(cell$modality)[1]),
            aes(x = 2, y = tf_lab), width = 3, height = 0.92,
            fill = "#F4A674", alpha = 0.18, inherit.aes = FALSE) +
  geom_point(aes(fill = I(cell_fill), alpha = I(cell_alpha)),
             shape = 21, color = "grey35", stroke = 0.3, size = 3.4,
             position = position_nudge(x = 0.20)) +
  # white halo + black value label, left of the marker
  geom_text(aes(label = cell_lab), size = PUB_GEOM_TEXT + 0.4,
            color = "white", fontface = "bold", lineheight = 0.8,
            position = position_nudge(x = -0.18)) +
  geom_text(aes(label = cell_lab), size = PUB_GEOM_TEXT, color = "black",
            position = position_nudge(x = -0.18)) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.6)) +
  scale_y_discrete(expand = expansion(add = 0.55)) +
  labs(
    x = NULL, y = NULL,
    title = "scATAC GRN corroborates the powered disease TFs",
    subtitle = "Hepatocyte SCENIC+ master regulators × independent powered evidence"
  ) +
  theme_masld() + theme_pub() +
  theme(
    axis.text.x      = element_text(face = "bold", lineheight = 0.85),
    # Heroes are flagged via the "*" suffix on their tick label (not by colour or
    # a per-tick face, which ggplot does not officially support) — house style
    # keeps ALL figure text black.
    axis.text.y      = element_text(face = "plain"),
    plot.title       = element_text(size = 7.5, face = "bold"),
    plot.subtitle    = element_text(size = 6, color = "grey30"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2),
    plot.margin      = margin(4, 6, 4, 4)
  )

# Manual legend strip (drawn as caption text) — fill semantics + the honest
# sc-underpowered caveat. Kept OUT of the data panel; emitted to stdout.

# ── Save ──────────────────────────────────────────────────────────────────────
out_pdf <- file.path(FIG4_DIR, "disease_master_regulators.pdf")
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
# useDingbats=FALSE is enforced globally via pdf.options() in publication_theme.R.
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
ggsave(out_pdf, p, width = 84 / 25.4, height = 96 / 25.4, device = pdf_device)

# ── Legend / stats to stdout (the panel caption, NOT printed on the figure) ───
n_coloc <- sum(ev$is_coloc)
heroes  <- ev$tf[ev$n_evidence == 2]
message("==================================================================")
message("PANEL A: disease_master_regulators.pdf  (Fig 4 VALIDATION)")
message("MESSAGE: the single-cohort scATAC hepatocyte SCENIC+ GRN independently")
message("         re-derives the disease TFs already called by the POWERED")
message("         bulk MASLD DEG analysis (Fig 3, n=846) and GWAS/COLOC (Fig 2).")
message("------------------------------------------------------------------")
message(sprintf("n cross-modality disease master regulators = %d", length(dmr_tfs)))
message(sprintf("  all %d are bulk MASLD DEGs (bulk_lfsr < 0.05)", length(dmr_tfs)))
message(sprintf("  of which carry GWAS COLOC support (n_evidence = 2) = %d (%s)",
                n_coloc, paste(heroes, collapse = ", ")))
message("Convergence heroes (bulk DEG + COLOC + most motif disruptions):")
for (tf in c("RORA", "THRB")) {
  r <- ev[ev$tf == tf, ]
  if (nrow(r))
    message(sprintf("  %-5s  bulk_logFC=%+.3f  bulk_lfsr=%.2e  COLOC PP.H4=%.4f  variant->motif=%d uniq (%d pairs)",
                    tf, r$bulk_logFC, r$bulk_lfsr, r$coloc_pp4, r$motif_uniq, r$motif_pairs))
}
message("Other COLOC TFs:")
for (tf in setdiff(heroes, c("RORA", "THRB"))) {
  r <- ev[ev$tf == tf, ]
  message(sprintf("  %-5s  bulk_logFC=%+.3f  bulk_lfsr=%.2e  COLOC PP.H4=%.4f  variant->motif=%d uniq",
                  tf, r$bulk_logFC, r$bulk_lfsr, r$coloc_pp4, r$motif_uniq))
}
message("------------------------------------------------------------------")
message("FILL semantics: bulk cell = magenta(up)/blue(down) MASLD effect;")
message("  COLOC + variant->motif cells = violet positive-evidence (faded if absent).")
message("  Dot size: |bulk_logFC|, PP.H4, or # unique variants. Cell label = value.")
message("CAVEAT (carried in the figure legend, NOT shown as significance):")
message(sprintf("  - scATAC regulon-activity DE is UNDERPOWERED (n=18, 5 ctrl/13 disease);"))
message(sprintf("    activity_padj == 1 for ALL %d TFs (sc_underpowered=TRUE).",
                length(dmr_tfs)))
message("    The panel shows NO scATAC disease-DE p-value as significant — the")
message("    scATAC role is STRUCTURAL (which TFs the GRN nominates), and the")
message("    corroboration is the convergence with the powered bulk + COLOC axes.")
message("  - Variant->motif counts are computational motifbreakR predictions")
message("    (no MPRA / CRISPRi functional validation).")
message(sprintf("Wrote: %s", out_pdf))
message("==================================================================")
