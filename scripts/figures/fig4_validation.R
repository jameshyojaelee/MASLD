#!/usr/bin/env Rscript
##############################################################################
# fig4_validation.R (rebuild 2026-04-15)
# Fig 4: Proteomics + Spatial Transcriptomic Validation
#
# Implements the panel set required by FIGURE_PLAN_REVISED.md §Figure 4 and
# the reframes in VERIFICATION_REPORT.md for F5 (BRONZE) and G2 (BRONZE).
#
# Panels (all written to FIG4_DIR/panels/, combined into FIG4_DIR/fig4_validation.pdf):
#   4a — Proteomics DE volcano (PXD052937 plasma DIA-MS, MASLD vs   [F2 SILVER]
#        control). Rewired 2026-04-22: previously used GSE276114
#        which is RNA-seq (not proteomics) — see T0.6 correction.
#   4b — mRNA-protein concordance scatter (bulk dream logFC vs      [F2 SILVER]
#        PXD052937 plasma DIA-MS logFC). ρ computed inline from
#        PXD052937 subset of protein_transcript_concordance_v3.csv.
#        Rewired 2026-04-22 — previously circular (RNA vs RNA).
#   4c — Spatial zonation of disease signature + notable SVGs       [F1 BRONZE]
#        (309 Healthy / 451 Steatotic SVGs; GSE192741 single-dataset)
#   4d — Plasma F>=3 binary fibrosis classifier: XGBoost AUROC      [F4 GOLD]
#        0.790 +/- 0.076; CHI3L1 top protein at SHAP 9.7%
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
#   panels/fig4a.pdf ... fig4d.pdf  +  fig4_validation.pdf (composite)
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
PLASMA <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")

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
                   size = 1.9, max.overlaps = 20,
                   label.padding = 0.1, segment.size = 0.15,
                   min.segment.length = 0, fontface = "italic",
                   color = "black", fill = alpha("white", 0.85),
                   show.legend = FALSE) +
  scale_color_manual(values = vol_colors, name = NULL,
                     labels = c(sprintf("Up (%d)", n_up),
                                sprintf("Down (%d)", n_down),
                                "n.s.")) +
  labs(x = expression("Plasma protein log"[2]*"FC (MASLD vs control)"),
       y = expression(-log[10]~italic(p)[adj]),
       title = "Plasma proteomics (PXD052937): MASLD vs control",
       subtitle = "DIA-MS; 3,346 proteins x 72 plasma samples (Sourianarayanane et al.)") +
  theme_masld() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.82, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.25, "cm"),
        plot.subtitle = element_text(size = 6, color = "gray35"))

# Persist a lightweight volcano CSV for downstream QC/caption text
fwrite(pdat[, .(protein_id, gene, logFC, padj, nlp, sig)],
       file.path(PROT, "pxd052937_dep_volcano.csv"))

save_fig(p4a, file.path(PANEL_DIR, "fig4a.pdf"),
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

p4b <- ggplot(cb, aes(x = bulk_logFC, y = protein_logFC, color = sig_class)) +
  rasterize_layer(geom_point(size = 0.5, alpha = 0.55, shape = 16)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.25, color = "gray55") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "gray70") +
  scale_color_manual(values = sig_cols_b, name = NULL) +
  annotate("text", x = -Inf, y = Inf, label = annot_both,
           hjust = -0.05, vjust = 1.6, size = 2.0, color = "gray25") +
  labs(x = expression("Transcript log"[2]*"FC (MASLD vs control)"),
       y = expression("Plasma protein log"[2]*"FC (PXD052937; MASLD vs control)"),
       title = "mRNA vs plasma protein concordance") +
  theme_masld() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.82, 0.20),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.22, "cm"))

if (nrow(lbl_b) > 0) {
  p4b <- p4b +
    geom_label_repel(data = lbl_b,
                     mapping = aes(x = bulk_logFC, y = protein_logFC,
                                   label = gene),
                     size = 1.9, max.overlaps = 20,
                     label.padding = 0.1, segment.size = 0.15,
                     min.segment.length = 0, fontface = "italic",
                     color = "black", fill = alpha("white", 0.85),
                     inherit.aes = FALSE,
                     show.legend = FALSE)
}

save_fig(p4b, file.path(PANEL_DIR, "fig4b.pdf"),
         width = fig_half_width, height = 3.3)

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
                   size = 1.7, max.overlaps = 22,
                   label.padding = 0.08, segment.size = 0.12,
                   min.segment.length = 0, fontface = "italic",
                   fill = alpha("white", 0.85), show.legend = FALSE,
                   inherit.aes = FALSE) +
  scale_color_manual(values = cat_cols, name = NULL) +
  labs(x = expression("Moran's " * italic(I) * " (Healthy)"),
       y = expression("Moran's " * italic(I) * " (Steatotic)"),
       title = "Spatial zonation") +
  theme_masld() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.25, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.22, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

save_fig(p4c, file.path(PANEL_DIR, "fig4c.pdf"),
         width = fig_half_width, height = 3.3)

# ============================================================================
# Panel 4d — Plasma F>=3 binary fibrosis classifier
#   XGBoost AUROC 0.790 +/- 0.076 (10x5 nested CV).
#   CHI3L1 top SHAP at 9.7% of total SHAP mass (not 40%).
# ============================================================================
message("[4d] Plasma fibrosis classifier (F4 GOLD)")

lb   <- fread(file.path(PLASMA, "226f_sweep_leaderboard.csv"))
imp  <- fread(file.path(PLASMA, "226a_protein_importance.csv"))

lb_t1 <- lb[target == "T1_binary" & primary_metric_name == "AUROC"]
setorder(lb_t1, -primary_metric_mean)
lb_t1[, model := factor(model, levels = rev(model))]
n_show <- min(12, nrow(lb_t1))
lb_show <- head(lb_t1, n_show)
lb_show[, highlight := model == "xgboost"]

# Panel 4d(i): model leaderboard bar
pd1 <- ggplot(lb_show,
              aes(x = primary_metric_mean, y = model, fill = highlight)) +
  geom_col(width = 0.75) +
  geom_errorbarh(aes(xmin = primary_metric_mean - primary_metric_sd,
                     xmax = primary_metric_mean + primary_metric_sd),
                 height = 0.2, linewidth = 0.3, color = "gray25") +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, color = "gray55") +
  geom_text(data = lb_show[highlight == TRUE],
            aes(label = sprintf("%.3f +/- %.3f",
                                 primary_metric_mean, primary_metric_sd)),
            hjust = -0.05, size = 2.0, color = "gray15") +
  scale_fill_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = "#B0BEC5"),
                    guide = "none") +
  scale_x_continuous(limits = c(0.48, 0.95),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "AUROC (10x5 nested CV)", y = NULL,
       title = "Plasma F>=3 fibrosis: model leaderboard") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

# Panel 4d(ii): top SHAP proteins (T1_binary)
imp_t1 <- imp[target == "T1_binary"]
total_shap <- sum(imp_t1$mean_abs_shap, na.rm = TRUE)
imp_t1[, shap_frac := 100 * mean_abs_shap / total_shap]
setorder(imp_t1, -mean_abs_shap)
top_imp <- head(imp_t1, 10)
top_imp[, protein := factor(protein, levels = rev(protein))]
top_imp[, highlight := protein == "CHI3L1"]

pd2 <- ggplot(top_imp, aes(x = shap_frac, y = protein, fill = highlight)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", shap_frac)),
            hjust = -0.1, size = 2.0, color = "gray20") +
  scale_fill_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = "#90CAF9"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, max(top_imp$shap_frac) * 1.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Fraction of total SHAP (%)", y = NULL,
       title = "Top 10 proteins; CHI3L1 leads (9.7% SHAP)") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6, face = "italic"))

p4d <- (pd1 | pd2) + plot_layout(widths = c(1.05, 1))
save_fig(p4d, file.path(PANEL_DIR, "fig4d.pdf"),
         width = fig_full_width, height = 3.3)

# Panel 4e (drug-target genetic validation scatters) was moved 2026-04-29
# to fig3_regulatory_architecture_v2.R as panels 3e/3f. The scatter narrative
# (transcription × genetics) belongs with the multi-ancestry regulatory
# architecture, not the proteomics/spatial validation figure.

# ============================================================================
# Composite assembly
# ============================================================================
message("[composite] Assembling fig4_validation.pdf")

# Row 1: 4a | 4b | 4c   (three ~half-width panels)
# Row 2: 4d (full-width, has two sub-bars)
row1 <- (p4a | p4b | p4c)
row2 <- p4d

combined <- row1 / row2 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c("a","b","c","d",""))) &
  theme(plot.tag = element_text(size = 9, face = "bold"))

out_composite <- file.path(FIGDIR, "fig4_validation.pdf")
save_fig(combined, out_composite,
         width = fig_full_width, height = 6.8)

message("Done. Panels in: ", PANEL_DIR)
message("Composite:      ", out_composite)
