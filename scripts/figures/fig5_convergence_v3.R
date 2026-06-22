##############################################################################
# fig5_convergence_v3.R — Fig 5 REDESIGN (2026-04-15, updated 2026-04-17)
#
# Binding spec: docs/manuscript/FIGURE_PLAN_REVISED.md §"Figure 5"
# Retractions enforced: D5 numbers banned (OR=1.67, 3.6%, Jaccard=0.012,
# "0/5,605 at 6/6 layers"). No "true NASH core" or "discrete F2 switch"
# language. Continuous 1.25x / Wilcoxon p=4.3e-43 only permitted in text.
#
# Panels:
#   5a  multi-evidence matrix heatmap (written directly by fig5_convergence.R → panels/fig5a.pdf)
#   5b  per-source DEG recovery + atlas sources-active histogram
#
# Output: $FIG5_DIR/panels/fig5{a,b}.pdf + $FIG5_DIR/fig5_convergence.pdf
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUTDIR  <- FIG5_DIR
PANDIR  <- file.path(OUTDIR, "panels")
dir.create(PANDIR, recursive = TRUE, showWarnings = FALSE)

# panels/fig5a.pdf is written directly by fig5_convergence.R — no copy needed here.

# ═══════════════════════════════════════════════════════════════════════════
# Load atlas
# ═══════════════════════════════════════════════════════════════════════════
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat(sprintf("Atlas: %d genes x %d cols\n", nrow(atlas), ncol(atlas)))

# ═══════════════════════════════════════════════════════════════════════════
# Compute 6-modality convergence flags — mirrors fig5_convergence.R exactly.
#   M1: Human DEG     — sig in ≥1 of 6 RNA contrasts (padj < 0.05)
#   M2: Mouse DEG     — mouse_meta_padj < 0.05
#   M3: Genetic       — max COLOC PP.H4 > 0.8
#   M4: Regulatory    — GWAS variant in hepatocyte peak with motif disruption
#   M5: Proteomics    — protein_padj < 0.05 (any dataset; no direction filter)
#   M6: Spatial       — spatial_max_I > 0.05
# ═══════════════════════════════════════════════════════════════════════════

# --- M3: genetic — canonical SuSiE-COLOC best PP.H4 (post gate-swap 2026-06-19) ---
# `coloc_susie_best_pp4` is the canonical per-gene COLOC (SuSiE with abf fallback);
# do NOT pmax the method-inconsistent per-trait `*_coloc_pp4` atlas columns.
if ("coloc_susie_best_pp4" %in% names(atlas)) {
  atlas[, max_coloc_pp4 := as.numeric(coloc_susie_best_pp4)]
} else {
  coloc_pp4_cols <- intersect(c("broadaway_coloc_pp4", "ukbb_alt_coloc_pp4",
    "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
    "bbj_alt_coloc_pp4"), names(atlas))
  if (length(coloc_pp4_cols) > 0) {
    atlas[, max_coloc_pp4 := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = coloc_pp4_cols]
  } else {
    atlas[, max_coloc_pp4 := NA_real_]
  }
}

# --- M4: GWAS-ATAC hepatocyte motif disruption (always derived fresh) ---
atlas[, n_motifs_disrupted := NA_integer_]
ann_file   <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv")
motif_file <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv")
if (file.exists(ann_file) && file.exists(motif_file)) {
  ann_atac   <- fread(ann_file)
  motif_atac <- fread(motif_file)
  hep_ann    <- ann_atac[cell_type == "Hepatocytes" & !is.na(nearest_gene) & nearest_gene != ""]
  var_gene   <- unique(hep_ann[, .(variant_id, nearest_gene, max_pip)])
  mg <- merge(motif_atac[, .(SNP_id, tf_name, motif_in_disease_regulon, alleleDiff)],
              var_gene, by.x = "SNP_id", by.y = "variant_id", allow.cartesian = TRUE)
  atac_gene_dt <- mg[, .(n_motifs_disrupted = uniqueN(tf_name)), by = nearest_gene]
  atlas[, n_motifs_disrupted := NULL]
  atlas <- merge(atlas, atac_gene_dt, by.x = "human_symbol", by.y = "nearest_gene", all.x = TRUE)
  cat(sprintf("  GWAS-ATAC: %d genes with hepatocyte motif disruption\n",
              sum(!is.na(atlas$n_motifs_disrupted))))
}

# --- M5: Proteomics (no direction concordance filter) ---
if (!"best_protein_logFC" %in% names(atlas)) {
  prot_file <- file.path(BASE,
    "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv")
  if (file.exists(prot_file)) {
    prot <- fread(prot_file)
    gene_prot <- prot[protein_padj < 0.05,
      .(best_protein_logFC = protein_logFC[which.min(protein_padj)],
        best_protein_padj  = min(protein_padj)),
      by = gene][!duplicated(gene)]
    atlas <- merge(atlas, gene_prot, by.x = "human_symbol", by.y = "gene", all.x = TRUE)
  } else {
    atlas[, best_protein_logFC := NA_real_]
  }
}
cat(sprintf("  Proteomics: %d genes with sig protein evidence\n",
            sum(!is.na(atlas$best_protein_logFC))))

# --- M6: Spatial (composite Moran's I from both Visium datasets) ---
if (!"spatial_max_I" %in% names(atlas)) {
  if ("spatial_morans_i" %in% names(atlas)) {
    # Use atlas column as proxy
    atlas[, spatial_max_I := spatial_morans_i]
  } else {
    sp_g <- file.path(BASE, "Analysis/Spatial/results/validation_bulk/spatial_bulk_merged_Guilliams_et_al.csv")
    sp_v <- file.path(BASE, "Analysis/Spatial/results/validation_bulk/spatial_bulk_merged_Vu_et_al.csv")
    mi_cols <- character(0)
    if (file.exists(sp_g)) {
      sg <- fread(sp_g)[, .(symbol, morans_I_g = I)][!duplicated(symbol)]
      atlas <- merge(atlas, sg, by.x = "human_symbol", by.y = "symbol", all.x = TRUE)
      mi_cols <- c(mi_cols, "morans_I_g")
    }
    if (file.exists(sp_v)) {
      sv <- fread(sp_v)[, .(symbol, morans_I_v = I)][!duplicated(symbol)]
      atlas <- merge(atlas, sv, by.x = "human_symbol", by.y = "symbol", all.x = TRUE)
      mi_cols <- c(mi_cols, "morans_I_v")
    }
    if (length(mi_cols) > 0) {
      atlas[, spatial_max_I := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = mi_cols]
    } else {
      atlas[, spatial_max_I := NA_real_]
    }
  }
}

# --- M1: Human DEG — sig in ≥1 of 6 RNA contrasts ---
deg_contrast_cols_v3 <- list(
  list(lfc = "bulk_logFC",        padj = "bulk_padj",        lfc_thresh = 0.5),
  list(lfc = "nafl_vs_ctrl_logFC", padj = "nafl_vs_ctrl_padj", lfc_thresh = 0.0),
  list(lfc = "nash_vs_ctrl_logFC", padj = "nash_vs_ctrl_padj", lfc_thresh = 0.0),
  list(lfc = "nafl_vs_nash_logFC", padj = "nafl_vs_nash_padj", lfc_thresh = 0.0),
  list(lfc = "adv_fib_logFC",      padj = "adv_fib_padj",      lfc_thresh = 0.0),
  list(lfc = "cirrhosis_logFC",    padj = "cirrhosis_padj",    lfc_thresh = 0.0)
)
atlas[, m1_human_deg := Reduce(`|`, lapply(deg_contrast_cols_v3, function(d) {
  if (!d$lfc %in% names(atlas)) return(rep(FALSE, nrow(atlas)))
  lfc  <- as.numeric(atlas[[d$lfc]])
  padj <- if (d$padj %in% names(atlas)) as.numeric(atlas[[d$padj]]) else rep(1, nrow(atlas))
  !is.na(padj) & padj < 0.05 & abs(lfc) > d$lfc_thresh
}))]
atlas[, m2_mouse_deg   := !is.na(mouse_meta_padj)   & mouse_meta_padj   < 0.05]
atlas[, m3_genetic     := !is.na(max_coloc_pp4)      & max_coloc_pp4     > 0.8]
atlas[, m4_regulatory  := !is.na(n_motifs_disrupted) & n_motifs_disrupted > 0]
atlas[, m5_proteomics  := !is.na(best_protein_logFC)]
atlas[, m6_spatial     := !is.na(spatial_max_I)      & spatial_max_I     > 0.05]

mod_cols <- c("m1_human_deg", "m2_mouse_deg", "m3_genetic",
              "m4_regulatory", "m5_proteomics", "m6_spatial")
atlas[, n_convergence := rowSums(as.matrix(.SD)), .SDcols = mod_cols]

# Primary DEG denominator: MASLD vs ctrl (bulk_padj < 0.05, |logFC| > 0.5)
is_deg <- !is.na(atlas$bulk_padj) & atlas$bulk_padj < 0.05 & abs(atlas$bulk_logFC) > 0.5
n_deg  <- sum(is_deg, na.rm = TRUE)
n_tot  <- nrow(atlas)

# ── Top panel: THRESHOLD-FREE expression concordance (Spearman ρ + GSEA) ─────
# Hard-threshold DEG overlap is arbitrary, discards effect-size/significance
# information, and double-thresholding underpowers low-power assays (proteomics
# n=130 vs bulk n=846). Instead, for each DE-comparable modality we report the
# Spearman correlation of the human disease effect (bulk_logFC) vs the modality's
# disease effect over ALL shared genes (threshold-free), with a bootstrap 95% CI,
# plus an fgsea test that the human DEG-up / DEG-down signature is concordantly
# enriched in the modality's ranked statistic (rank = sign(logFC)·−log10 padj;
# the make_ranks idiom from Cross_Species_Concordance/02a_fgsea_concordance.R).
# EXCLUDED from this panel (kept only in the convergence count below):
#   · Spatial — its contrasts are zonal/niche (sh-vs-pt etc.), NOT disease-vs-
#     control, so they don't rank-correlate with bulk (ρ≈0 / negative). Spatial
#     stays as a spatially-structured-evidence flag in the histogram.
#   · Genetic (COLOC posterior) and Regulatory/ATAC (GWAS-locus motif count) are
#     not DE effect sizes → set/count evidence, convergence-count only.
suppressPackageStartupMessages(library(fgsea))

deg_up   <- unique(atlas[is_deg & bulk_logFC > 0, human_symbol])
deg_down <- unique(atlas[is_deg & bulk_logFC < 0, human_symbol])
deg_sets <- list(DEG_up = deg_up, DEG_down = deg_down)

p_to_stars <- function(p) ifelse(is.na(p), "",
  ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", ""))))

de_mods <- list(
  list(s = "Mouse DEG",         lfc = "mouse_meta_logFC",    padj = "mouse_meta_padj"),
  list(s = "Proteomics",        lfc = "best_protein_logFC",  padj = "best_protein_padj"),
  list(s = "Single-cell (hep)", lfc = "sc_hepatocyte_logFC", padj = "sc_hepatocyte_padj")
)

set.seed(42)
conc <- rbindlist(lapply(de_mods, function(m) {
  d <- atlas[, .(g = human_symbol, hl = bulk_logFC, ml = get(m$lfc), mp = get(m$padj))]
  d <- d[!is.na(g) & !is.na(hl) & !is.na(ml) & !is.na(mp)][!duplicated(g)]
  n <- nrow(d)
  rho <- if (n >= 10) cor(d$hl, d$ml, method = "spearman") else NA_real_
  bs  <- if (n >= 10) replicate(1000, { i <- sample.int(n, replace = TRUE);
              cor(d$hl[i], d$ml[i], method = "spearman") }) else rep(NA_real_, 2)
  ci  <- if (n >= 10) unname(quantile(bs, c(0.025, 0.975))) else c(NA_real_, NA_real_)
  rk  <- setNames(sign(d$ml) * -log10(pmax(d$mp, 1e-300)), d$g); rk <- rk[is.finite(rk)]
  fg  <- suppressWarnings(fgsea(deg_sets, rk, minSize = 10, maxSize = 6000, eps = 0))
  nes_up <- fg[pathway == "DEG_up",   NES]; p_up <- fg[pathway == "DEG_up",   padj]
  nes_dn <- fg[pathway == "DEG_down", NES]; p_dn <- fg[pathway == "DEG_down", padj]
  cand <- c()
  if (length(nes_up) == 1 && isTRUE(nes_up > 0)) cand <- c(cand, p_up)
  if (length(nes_dn) == 1 && isTRUE(nes_dn < 0)) cand <- c(cand, p_dn)
  conc_p <- if (length(cand)) suppressWarnings(min(cand, na.rm = TRUE)) else NA_real_
  if (!is.finite(conc_p)) conc_p <- NA_real_
  data.table(source = m$s, n_shared = n, rho = rho, ci_lo = ci[1], ci_hi = ci[2],
             nes_up   = if (length(nes_up) == 1) nes_up else NA_real_,
             nes_down = if (length(nes_dn) == 1) nes_dn else NA_real_,
             gsea_padj = conc_p)
}))
conc[, star := p_to_stars(gsea_padj)]
conc[, source := factor(source, levels = rev(source))]
cat("\nExpression concordance (Spearman ρ + DEG-signature GSEA):\n"); print(conc)

sa_dist <- atlas[, .N, by = n_convergence][order(n_convergence)]
setnames(sa_dist, "n_convergence", "sa")
sa_dist[, sa := as.integer(sa)]
cat("\nConvergence distribution (0–6 modalities):\n"); print(sa_dist)

max_sa         <- max(sa_dist$sa, na.rm = TRUE)
top_sa_n       <- sa_dist[sa == max_sa, N]
convergent_ge3 <- sum(sa_dist[sa >= 3, N])
convergent_ge4 <- sum(sa_dist[sa >= 4, N])
convergent_ge5 <- sum(sa_dist[sa >= 5, N])
cat(sprintf("\nHeadline: %d genes at max %d/6; %d at >=3/6; %d at >=4/6; %d at >=5/6\n",
            top_sa_n, max_sa, convergent_ge3, convergent_ge4, convergent_ge5))

# ═══════════════════════════════════════════════════════════════════════════
# PANEL 5b — per-source DEG recovery + atlas sources-active histogram
# ═══════════════════════════════════════════════════════════════════════════

# 5b-top: threshold-free expression concordance — Spearman ρ bars (bootstrap 95%
# CI), bar-tip label = ρ + DEG-signature GSEA significance stars.
rho_hi <- max(conc$ci_hi, conc$rho, na.rm = TRUE)
p5a_top <- ggplot(conc, aes(x = rho, y = source)) +
  geom_col(fill = masld_colors$mash, width = 0.6, alpha = 0.9) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.16,
                 linewidth = 0.3, color = "grey25") +
  geom_text(aes(x = ci_hi, label = sprintf("ρ=%.2f%s", rho, star)),
            hjust = -0.12, size = 2.1, family = "Helvetica") +
  scale_x_continuous(limits = c(0, rho_hi * 1.38),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Spearman ρ  (human disease logFC vs modality;  ✱✱✱ DEG GSEA padj<0.001)",
       y = NULL,
       title = "b",
       subtitle = "Expression concordance of human DEGs") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7, color = "grey30"),
        axis.title.x = element_text(size = 6.2),
        axis.text.y = element_text(size = 7))

# 5b-bottom: atlas-wide convergence histogram (0–6 modalities)
p5a_bot <- ggplot(sa_dist, aes(x = factor(sa), y = N)) +
  geom_col(aes(fill = sa >= 3), width = 0.78, alpha = 0.95) +
  geom_text(aes(label = comma(N)),
            vjust = -0.3, size = 2.1, family = "Helvetica") +
  scale_fill_manual(values = c(`FALSE` = "grey65",
                               `TRUE`  = masld_colors$conserved),
                    guide = "none") +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.15))) +
  # PI directive (2026-06-11): short single title line only. Atlas-wide
  # convergence counts (>=3/6, >=4/6, >=5/6, max-modality genes) are printed to
  # stdout above and belong in the figure caption, not on the panel.
  labs(x = "Number of independent modalities active (0 to 6)",
       y = "Genes",
       subtitle = "Atlas-wide convergence") +
  theme_masld(base_size = 7) +
  theme(plot.subtitle = element_text(size = 7, color = "grey30"))

p5b <- p5a_top / p5a_bot + plot_layout(heights = c(1, 1))

ggsave(file.path(PANDIR, "fig5b.pdf"), p5b,
       width = 4.2, height = 4.6, device = cairo_pdf)
cat("Saved: panels/fig5b.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════
# Sidecar CSVs: concordance stats (ρ/CI/NES/padj) + convergence histogram
# ═══════════════════════════════════════════════════════════════════════════
fwrite(conc, file.path(OUTDIR, "fig5b_concordance.csv"))
fwrite(sa_dist, file.path(OUTDIR, "fig5b_convergence_histogram.csv"))

cat("\nDONE.\n")
