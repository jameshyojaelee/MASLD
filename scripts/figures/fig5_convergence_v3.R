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

# --- M3: max COLOC PP.H4 across GWAS sources ---
coloc_pp4_cols <- intersect(c("broadaway_coloc_pp4", "ukbb_alt_coloc_pp4",
  "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
  "bbj_alt_coloc_pp4"), names(atlas))
if (length(coloc_pp4_cols) > 0) {
  atlas[, max_coloc_pp4 := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = coloc_pp4_cols]
} else {
  atlas[, max_coloc_pp4 := NA_real_]
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

cov_tbl <- data.table(
  source  = c("Human DEG", "Mouse DEG", "Genetic causal",
               "Regulatory/ATAC", "Proteomics", "Spatial"),
  deg_n   = sapply(mod_cols, function(s) sum(atlas[[s]] & is_deg, na.rm = TRUE)),
  atlas_n = sapply(mod_cols, function(s) sum(atlas[[s]], na.rm = TRUE))
)
cov_tbl[, deg_pct   := 100 * deg_n / n_deg]
cov_tbl[, atlas_pct := 100 * atlas_n / n_tot]
cov_tbl[, source := factor(source, levels = source)]

cat("\nPer-modality recovery:\n"); print(cov_tbl)

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

# 5b-top: per-modality coverage of primary human DEGs (horizontal bars)
p5a_top <- ggplot(cov_tbl, aes(x = deg_pct, y = source)) +
  geom_col(fill = masld_colors$mash, width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", deg_n, deg_pct)),
            hjust = -0.05, size = 2.1, family = "Helvetica") +
  scale_y_discrete(limits = rev(levels(cov_tbl$source))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = sprintf("%% of %s primary DEGs with evidence", comma(n_deg)),
       y = NULL,
       title = "b",
       subtitle = "Per-modality coverage of primary human DEGs") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7, color = "grey30"),
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
  labs(x = "Number of independent modalities active (0 to 6)",
       y = "Genes",
       subtitle = sprintf(
         "Atlas-wide: %s genes >=3/6 modalities; %s >=4/6; %s >=5/6; %s at max %d/6",
         comma(convergent_ge3),
         comma(convergent_ge4),
         comma(convergent_ge5),
         comma(top_sa_n),
         max_sa)) +
  theme_masld(base_size = 7) +
  theme(plot.subtitle = element_text(size = 7, color = "grey30"))

p5b <- p5a_top / p5a_bot + plot_layout(heights = c(1, 1))

ggsave(file.path(PANDIR, "fig5b.pdf"), p5b,
       width = 4.2, height = 4.6, device = cairo_pdf)
cat("Saved: panels/fig5b.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════
# Sidecar CSV: per-source DEG recovery numbers used in 5b
# ═══════════════════════════════════════════════════════════════════════════
fwrite(cov_tbl, file.path(OUTDIR, "fig5b_modality_recovery.csv"))
fwrite(sa_dist, file.path(OUTDIR, "fig5b_convergence_histogram.csv"))

cat("\nDONE.\n")
