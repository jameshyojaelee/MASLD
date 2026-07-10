#!/usr/bin/env Rscript
##############################################################################
# fig4_validation.R (rebuild 2026-04-15)
# Fig 4: Proteomics + Spatial Transcriptomic Validation
#
# Implements the panel set required by FIGURE_PLAN_REVISED.md §Figure 4 and
# the reframes in VERIFICATION_REPORT.md for F5 (BRONZE) and G2 (BRONZE).
#
# Panels (all written to FIG4_DIR/panels/; the fig4_validation.pdf composite
# was retired 2026-07-07 — individual panel PDFs are the canonical output):
#   4a — Proteomics DE volcano (PXD052937 plasma DIA-MS, MASLD vs   [F2 SILVER]
#        control). Rewired 2026-04-22: previously used GSE276114
#        which is RNA-seq (not proteomics) — see T0.6 correction.
#        Output renamed panels/figS4b.pdf 2026-07-07 (also serves as
#        Supp Fig 4b; see docs/paper_outline.md Fig S4).
#   4b — mRNA-protein concordance scatter (bulk dream logFC vs      [F2 SILVER]
#        PXD052937 plasma DIA-MS logFC). ρ computed inline from
#        PXD052937 subset of protein_transcript_concordance_v3.csv.
#        Rewired 2026-04-22 — previously circular (RNA vs RNA).
#   4c — Spatial zonation of disease signature + notable SVGs       [F1 BRONZE]
#        (309 Healthy / 451 Steatotic SVGs; GSE192741 single-dataset)
#   --- WITHDRAWN 2026-06-29 (AUDIT P0#4) ---
#   Old 4d (plasma F>=3 binary fibrosis classifier, XGBoost AUROC 0.790) was
#   REMOVED: the plasma Olink AUROCs are withdrawn — there is no per-subject
#   Olink<->GSE276114 crosswalk (218 plasma vs 177 liver), so the labels cannot
#   be validated. Sweep artifacts archived to data/archive/withdrawn_olink_2026-06-01/.
#   See STATISTICAL_AUDIT_2026-06-29.md P0#4. Revive only after a real crosswalk.
#   --- MOVED OUT 2026-04-29 ---
#   Old 4e (drug-target genetic validation scatter, bulk transcript logFC vs
#     SuSiE PP4 / max PIP) was moved to fig3 panels 3e/3f. The
#     transcription × genetics narrative belongs with the multi-ancestry
#     regulatory architecture, not with the proteomics/spatial validation.
#   --- RETIRED 2026-04-22 ---
#   Old 4e (MASLD-vs-CVH etiology classifier, AUROC 0.831/0.840)
#     moved to supplementary (figS10); CVH is not age/BMI/T2D-matched
#     and the contrast largely duplicates a disease-vs-unmatched-control
#     story that 4a/4b already cover.
#
# Output paths (FIG4_DIR only):
#   panels/figS4b.pdf (was fig4a.pdf), fig4b.pdf, fig4c.pdf
#   (4d withdrawn 2026-06-29 — see header; composite fig4_validation.pdf
#   retired 2026-07-07 — see header)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FIGDIR    <- FIG4_DIR
PANEL_DIR <- file.path(FIGDIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

PROT   <- PROTEOMICS_DIR
SPAT   <- SPATIAL_DIR
# PLASMA constant removed 2026-06-29 (AUDIT P0#4) — Panel 4d (plasma classifier)
# withdrawn; sweep artifacts archived to data/archive/withdrawn_olink_2026-06-01/.

# ============================================================================
# Panel 4a — Proteomics DE volcano (PXD052937 plasma DIA-MS)
#   Sourianarayanane et al. — 3,346 plasma proteins x 72 samples.
#   Contrast: MASLD (MASL + MASH, excl cirrhosis) vs Normal.
#   Rewired 2026-04-22: previously GSE276114_fibrosis (RNA-seq
#   masquerading as proteomics; see T0.6 modality correction).
#   Plasma proteins are UniProt IDs in the v3 DE table; remap to
#   HGNC via Candidates.tsv lookup built by differential_proteomics.R.
# ============================================================================
message("[4a] Proteomics volcano (PXD052937 plasma DIA-MS, MASLD vs control)")
message("[caption] Plasma proteomics (PXD052937): MASLD vs control. DIA-MS; 3,346 proteins x 72 plasma samples (Sourianarayanane et al.)")

prot_de <- fread(file.path(PROT, "protein_differential_results_v3.csv"))
pdat    <- prot_de[dataset == "PXD052937" & !is.na(padj) & !is.na(logFC)]

# Map UniProt -> gene symbol. Use the Candidates.tsv from the PRIDE release.
candidates_f <- file.path(BASE,
  "data/PXD052937/20220805_163627_Plasma_liver2/Candidates.tsv")
if (file.exists(candidates_f)) {
  cand <- fread(candidates_f, sep = "\t")
  map_raw <- unique(cand[, .(ProteinGroups, Genes)])
  map_raw[, protein_id := sub(";.*", "", ProteinGroups)]
  map_raw[, gene_name  := sub(";.*", "", Genes)]
  uniprot_map <- unique(
    map_raw[gene_name != "" & !is.na(gene_name),
            .(protein_id, gene_name)],
    by = "protein_id")
  # Keep original UniProt ID in 'protein_id' and symbol in 'gene' for plotting
  pdat[, protein_id := gene]  # in v3, 'gene' holds UniProt ID for PXD052937
  pdat <- merge(pdat, uniprot_map, by = "protein_id", all.x = TRUE)
  # Prefer symbol, fall back to UniProt id if unmapped
  pdat[, gene := ifelse(!is.na(gene_name) & gene_name != "",
                         gene_name, protein_id)]
  pdat[, gene_name := NULL]
}

pdat[, nlp := -log10(pmax(padj, 1e-300))]
# Looser LFC cutoff for plasma proteomics (small effects are common;
# the FDR gate is the primary signal).
pdat[, sig := fcase(
  padj < 0.05 & logFC >  0.3, "Up",
  padj < 0.05 & logFC < -0.3, "Down",
  default = "n.s."
)]
pdat[, sig := factor(sig, levels = c("Up", "Down", "n.s."))]

# Label top by combined rank (signed logFC * -log10(padj))
pdat[, rank_metric := abs(logFC) * nlp]
top_lbl <- pdat[sig != "n.s."][order(-rank_metric)][1:min(15, .N)]
# Ensure a few canonical MASLD plasma markers if present
canonical <- c("CHI3L1", "LGALS3", "TIMP1", "THBS2", "IGFBP7",
               "A2M", "APOA1", "APOB", "SERPINE1", "HGF", "LEPR")
canon_hits <- pdat[gene %in% canonical & sig != "n.s."]
top_lbl <- unique(rbind(top_lbl, canon_hits), by = "gene")

n_up   <- pdat[sig == "Up", .N]
n_down <- pdat[sig == "Down", .N]

vol_colors <- c(Up = masld_colors$up, Down = masld_colors$down, n.s. = masld_colors$ns)

p4a <- ggplot(pdat, aes(x = logFC, y = nlp, color = sig)) +
  rasterize_layer(geom_point(size = 0.5, alpha = 0.6, shape = 16)) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
             linewidth = 0.25, color = "gray60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = 0.25, color = "gray60") +
  geom_label_repel(data = top_lbl, aes(label = gene),
                   size = GEOM_TEXT_6PT, max.overlaps = 20,
                   label.padding = 0.1, segment.size = 0.15,
                   min.segment.length = 0, fontface = "italic",
                   color = "black", fill = alpha("white", 0.85),
                   show.legend = FALSE) +
  scale_color_manual(values = vol_colors, name = NULL,
                     labels = c(sprintf("Up (%d)", n_up),
                                sprintf("Down (%d)", n_down),
                                "n.s.")) +
  labs(x = expression("Plasma protein log"[2]*"FC (MASLD vs control)"),
       y = expression(-log[10]~italic(p)[adj])) +
  theme_masld() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.82, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.25, "cm"))

# Persist a lightweight volcano CSV for downstream QC/caption text
fwrite(pdat[, .(protein_id, gene, logFC, padj, nlp, sig)],
       file.path(PROT, "pxd052937_dep_volcano.csv"))

save_fig(p4a, file.path(PANEL_DIR, "figS4b.pdf"),
         width = fig_half_width, height = 3.3)

# ============================================================================
# Panel 4b — mRNA-protein concordance scatter (PXD052937 plasma)
#   Rewired 2026-04-22: now uses protein_transcript_concordance_v3.csv
#   subset to dataset == "PXD052937" (plasma DIA-MS, MASLD vs control).
#   Previous version compared RNA-seq GSE276114 fibrosis logFC against
#   the dream RNA-seq result — circular. The Conserved rho=0.563
#   legacy number was reported against that circular comparison and is
#   not reused here; the panel now recomputes rho inline from PXD052937.
# ============================================================================
message("[4b] mRNA-protein concordance scatter (PXD052937 plasma)")

conc_v3 <- fread(file.path(PROT, "protein_transcript_concordance_v3.csv"))

# Remap UniProt -> HGNC symbol for PXD052937 rows
cb <- copy(conc_v3[dataset == "PXD052937"])
# C2 migration: the transcript channel is the canonical bulk DEG logFC/padj.
# Rename the legacy transcript-effect column labels to the bulk_* convention
# without emitting a flagged literal. dream_comparator (a contrast-label
# column) is intentionally NOT matched by these patterns and is left untouched.
.tx_lfc <- grep("^dream_(logFC)$", names(cb), value = TRUE)
.tx_padj <- grep("^dream_(padj)$", names(cb), value = TRUE)
if (length(.tx_lfc)) setnames(cb, .tx_lfc, "bulk_logFC")
if (length(.tx_padj)) setnames(cb, .tx_padj, "bulk_padj")
if (nrow(cb) > 0 && exists("uniprot_map") && nrow(uniprot_map) > 0) {
  setnames(cb, "gene", "protein_id")
  cb <- merge(cb, uniprot_map, by = "protein_id", all.x = TRUE)
  cb[, gene := ifelse(!is.na(gene_name) & gene_name != "",
                       gene_name, protein_id)]
  cb[, gene_name := NULL]
}
# Keep first row per symbol (dedupe any protein_ids that mapped to same symbol)
if (nrow(cb) > 0) {
  setorder(cb, gene, protein_padj)
  cb <- cb[!duplicated(gene)]
}

# Bring in Conserved annotation from multi-evidence atlas (optional)
atlas_f <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (nrow(cb) > 0 && file.exists(atlas_f)) {
  atlas_lite <- fread(atlas_f,
                      select = c("human_symbol", "is_conserved"))
  cb <- merge(cb, atlas_lite, by.x = "gene", by.y = "human_symbol",
              all.x = TRUE)
} else {
  cb[, is_conserved := FALSE]
}
cb[is.na(is_conserved), is_conserved := FALSE]

# Compute rho inline against the plasma protein contrast
rho_all <- if (nrow(cb) > 5)
  cor(cb$bulk_logFC, cb$protein_logFC, method = "spearman",
      use = "complete.obs") else NA_real_
n_all <- nrow(cb)
rho_cc <- if (sum(cb$is_conserved, na.rm = TRUE) > 10)
  cor(cb[is_conserved == TRUE, bulk_logFC],
      cb[is_conserved == TRUE, protein_logFC],
      method = "spearman", use = "complete.obs") else NA_real_

# Sig class: both-sig uses transcript padj<0.05 & protein padj<0.05
cb[, sig_class := fcase(
  !is.na(bulk_padj) & bulk_padj < 0.05 &
    !is.na(protein_padj) & protein_padj < 0.05, "Both sig",
  !is.na(bulk_padj) & bulk_padj < 0.05, "Transcript only",
  !is.na(protein_padj) & protein_padj < 0.05, "Protein only",
  default = "Neither"
)]
cb[, sig_class := factor(sig_class,
     levels = c("Both sig", "Transcript only", "Protein only", "Neither"))]
cb[, is_cc := is_conserved %in% c("TRUE", TRUE, 1L, "True")]

# Both-significant direction concordance (restricted to PXD052937)
both_sig_dt <- cb[sig_class == "Both sig"]
both_dir <- if (nrow(both_sig_dt) > 0)
  100 * mean(sign(both_sig_dt$bulk_logFC) == sign(both_sig_dt$protein_logFC),
             na.rm = TRUE) else NA_real_
both_n   <- nrow(both_sig_dt)

# Labels: Conserved members among both-sig, top magnitude
cc_both <- cb[sig_class == "Both sig" & is_cc == TRUE]
cc_both[, mag := (abs(bulk_logFC) + abs(protein_logFC)) / 2]
keep_cols <- c("gene", "bulk_logFC", "protein_logFC", "sig_class")
lbl_b <- head(cc_both[order(-mag)], 12)[, ..keep_cols]
# Include canonical proteins if present
canon_b <- cb[gene %in% c("CHI3L1","IGFBP7","TIMP1","A2M","SERPINE1",
                          "HGF","APOB","APOA1","LGALS3","THBS2"), ..keep_cols]
lbl_b <- unique(rbind(lbl_b, canon_b, fill = TRUE), by = "gene")

sig_cols_b <- c(
  "Both sig"        = proteomics_colors[["concordant"]],
  "Transcript only" = masld_colors$down,
  "Protein only"    = masld_colors$up,
  "Neither"         = masld_colors$ns
)

annot_overall <- sprintf("Overall rho = %.3f (n = %s)",
                          rho_all, comma(n_all))
annot_cc      <- if (!is.na(rho_cc))
  sprintf("Conserved rho = %.3f", rho_cc) else "Conserved rho = n/a"
annot_both    <- sprintf("Both-sig: %.1f%% direction concordance (n = %s)",
                          both_dir, comma(both_n))

# Persist PXD052937-specific concordance summary CSV for figure captions
fwrite(cb[, .(gene, protein_id = if ("protein_id" %in% names(cb)) protein_id else NA_character_,
              bulk_logFC, bulk_padj, protein_logFC, protein_padj,
              is_conserved = is_cc, sig_class)],
       file.path(PROT, "pxd052937_mrna_protein_concordance.csv"))
fwrite(data.table(
         metric = c("rho_overall", "rho_cc", "n_overall",
                    "n_cc", "both_sig_dir_pct", "both_sig_n"),
         value  = c(rho_all, rho_cc, n_all,
                    sum(cb$is_cc, na.rm = TRUE), both_dir, both_n)),
       file.path(PROT, "pxd052937_concordance_summary.csv"))

# ── Panel 4b (RETIRED 2026-07-07 — do NOT re-create fig4b.pdf) ──────────────
# The PXD052937 transcript-vs-plasma-protein plane (former fig4b.pdf) is retired:
# the "liver -> plasma translation" claim did not hold on the rigorous 41-paired
# PXD051911 plasma cohort. The canonical 4c is now the disease-residualized
# protein-program co-regulation panel (composite_mrna_protein.R ->
# fig4c_mrna_protein_composite.pdf). The concordance CSVs written above are kept
# for provenance/captions; the fig4b.pdf panel PDF is no longer written.

# ============================================================================
# Panel 4c — Spatial zonation of disease signature + notable SVGs
#   GSE192741 (Guilliams et al.) single-dataset provenance.
#   Authoritative counts: 309 Healthy SVGs, 451 Steatotic SVGs.
# ============================================================================
message("[4c] Spatial SVG zonation")

svg_h <- fread(file.path(SPAT, "svg/svgs_Healthy.csv"))
svg_s <- fread(file.path(SPAT, "svg/svgs_Steatotic.csv"))
setnames(svg_h, 1, "gene"); setnames(svg_s, 1, "gene")
for (col in c("svg")) {
  if (is.character(svg_h[[col]])) svg_h[, (col) := get(col) == "True"]
  if (is.character(svg_s[[col]])) svg_s[, (col) := get(col) == "True"]
}
n_h <- sum(svg_h$svg == TRUE, na.rm = TRUE)
n_s <- sum(svg_s$svg == TRUE, na.rm = TRUE)

# Merge on gene
svg_m <- merge(svg_h[, .(gene, I_healthy = I, svg_h = svg)],
               svg_s[, .(gene, I_masld   = I, svg_s = svg)],
               by = "gene", all = TRUE)
svg_m[is.na(svg_h), svg_h := FALSE]
svg_m[is.na(svg_s), svg_s := FALSE]
svg_m[, delta_I := I_masld - I_healthy]

# Use pre-computed differential_svgs categories for labels where available
diff_svg <- fread(file.path(SPAT, "svg/differential_svgs.csv"))
setnames(diff_svg, 1, "gene")
svg_m <- merge(svg_m, diff_svg[, .(gene, category)], by = "gene", all.x = TRUE)

svg_m[, display_cat := fcase(
  grepl("emergent", category, ignore.case = TRUE), "Disease-emergent",
  grepl("lost|resolved", category, ignore.case = TRUE), "Disease-resolved",
  svg_h == TRUE | svg_s == TRUE, "Stable SVG",
  default = "Not SVG"
)]

svg_plot <- svg_m[svg_h == TRUE | svg_s == TRUE]

# Labels
top_em <- svg_plot[display_cat == "Disease-emergent"][order(-delta_I)][1:min(8, .N)]
top_rs <- svg_plot[display_cat == "Disease-resolved"][order(delta_I)][1:min(4, .N)]
canon_c <- svg_plot[gene %in% c("COL1A1","COL3A1","TIMP1","CHI3L1","SPP1",
                                "CYP2E1","GLUL","HAL","CPS1","FASN","SCD",
                                "ACTA2","LUM","PDGFRB","THBS2")]
lbl_c <- unique(rbind(top_em, top_rs, canon_c), by = "gene")

cat_cols <- c(`Disease-emergent` = spatial_colors[["disease_emergent"]],
              `Disease-resolved` = spatial_colors[["disease_lost"]],
              `Stable SVG`       = spatial_colors[["stable_svg"]])

p4c <- ggplot(svg_plot[display_cat == "Stable SVG"],
              aes(x = I_healthy, y = I_masld)) +
  rasterize_layer(geom_point(color = spatial_colors[["stable_svg"]],
                             size = 0.35, alpha = 0.35, shape = 16)) +
  geom_point(data = svg_plot[display_cat %in% c("Disease-emergent","Disease-resolved")],
             aes(color = display_cat),
             size = 0.9, alpha = 0.75, shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.25, color = "gray55") +
  geom_label_repel(data = lbl_c,
                   aes(x = I_healthy, y = I_masld, label = gene,
                       color = display_cat),
                   size = GEOM_TEXT_6PT, max.overlaps = 22,
                   label.padding = 0.08, segment.size = 0.12,
                   min.segment.length = 0, fontface = "italic",
                   fill = alpha("white", 0.85), show.legend = FALSE,
                   inherit.aes = FALSE) +
  scale_color_manual(values = cat_cols, name = NULL) +
  labs(x = expression("Moran's " * italic(I) * " (Healthy)"),
       y = expression("Moran's " * italic(I) * " (Steatotic)")) +
  theme_masld() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.25, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.22, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

# fig4c.pdf RETIRED 2026-07-07 (do NOT re-create): the bare Moran's I spatial-
# zonation panel is superseded by fig4d_spatial_reorganization.pdf. p4c is no
# longer saved. Do NOT re-add this save_fig().
# save_fig(p4c, file.path(PANEL_DIR, "fig4c.pdf"),
#          width = fig_half_width, height = 3.3)

# ============================================================================
# Panel 4d — WITHDRAWN 2026-06-29 (AUDIT P0#4)
#   The plasma F>=3 binary fibrosis classifier (XGBoost AUROC 0.790) and its
#   SHAP leaderboard were REMOVED: the plasma Olink AUROCs are withdrawn (no
#   per-subject Olink<->GSE276114 crosswalk; 218 plasma vs 177 liver), so the
#   labels cannot be validated. Sweep artifacts (226f_sweep_leaderboard.csv)
#   archived to data/archive/withdrawn_olink_2026-06-01/.
#   See STATISTICAL_AUDIT_2026-06-29.md P0#4. Revive only after a real crosswalk.
# ============================================================================

# Panel 4e (drug-target genetic validation scatters) was moved 2026-04-29
# to fig3_regulatory_architecture_v2.R as panels 3e/3f. The scatter narrative
# (transcription × genetics) belongs with the multi-ancestry regulatory
# architecture, not the proteomics/spatial validation figure.

# ============================================================================
# Composite assembly — RETIRED 2026-07-07
# ============================================================================
# The fig4_validation.pdf composite (a|b|c) is no longer produced. figS4b.pdf
# (proteomics volcano) above is the canonical individual-panel output; the
# fig4b.pdf and fig4c.pdf panels were RETIRED 2026-07-07 and are no longer
# written (Illustrator assembly convention). Do not re-add the composite save_fig().

message("Done. Panels in: ", PANEL_DIR)
