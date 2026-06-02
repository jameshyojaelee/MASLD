#!/usr/bin/env Rscript
# 02b_ssgsea_concordance.R
# ---------------------------------------------------------------------------
# Per-sample pathway activity via ssGSEA (GSVA), then correlate activity
# between human and mouse conditions
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(msigdbr)
  library(ggplot2)
})

pdf.options(useDingbats = FALSE)

cat("=== Phase 2b: ssGSEA Concordance ===\n\n")

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT   <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT   <- file.path(H_INT, "results/gene_annotation")
INT_DIR <- file.path(H_INT, "results/integration")
MOUSE_UI <- file.path(BASE, "Mouse/Unified_Integration")
WD      <- file.path(BASE, "Analysis/Cross_Species_Concordance")
RES     <- file.path(WD, "results")

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")  # 2026-05-29: LIDPAD archived/dropped

# ============================================================
#  Gene sets (Hallmark only for ssGSEA — smaller, more interpretable)
# ============================================================
cat("Loading Hallmark gene sets...\n")
h_hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
h_pathways <- split(h_hallmark$gene_symbol, h_hallmark$gs_name)

m_hallmark <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
m_pathways <- split(m_hallmark$gene_symbol, m_hallmark$gs_name)
cat(sprintf("  Hallmark sets: %d\n", length(h_pathways)))

# ============================================================
#  Human ssGSEA
# ============================================================
cat("\n=== Human ssGSEA ===\n")
h_dge    <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
h_logcpm <- edgeR::cpm(h_dge, log = TRUE, prior.count = 1)
h_meta   <- as.data.table(readRDS(file.path(INT_DIR, "meta_matched.rds")))

# Map ENSG to symbol
h_annot <- fread(file.path(ANNOT, "human_ensg_to_symbol.tsv"))
h_annot[, gene_base := gsub("\\..*", "", gene_id)]

rownames_base <- gsub("\\..*", "", rownames(h_logcpm))
sym_map <- h_annot[match(rownames_base, gene_base), symbol]
valid <- !is.na(sym_map) & sym_map != ""

# Aggregate duplicate symbols by averaging (not dropping)
h_mat <- h_logcpm[valid, ]
rownames(h_mat) <- sym_map[valid]
if (any(duplicated(rownames(h_mat)))) {
  dup_syms <- unique(rownames(h_mat)[duplicated(rownames(h_mat))])
  uniq_rows <- which(!rownames(h_mat) %in% dup_syms)
  avg_list <- lapply(dup_syms, function(sym) {
    colMeans(h_mat[rownames(h_mat) == sym, , drop = FALSE])
  })
  avg_mat <- do.call(rbind, avg_list)
  rownames(avg_mat) <- dup_syms
  h_mat <- rbind(h_mat[uniq_rows, , drop = FALSE], avg_mat)
  cat(sprintf("  Averaged %d duplicate symbols\n", length(dup_syms)))
}

cat(sprintf("  Expression matrix: %d genes × %d samples\n", nrow(h_mat), ncol(h_mat)))

# Run ssGSEA
cat("  Running ssGSEA (this may take a few minutes)...\n")
h_param <- ssgseaParam(h_mat, h_pathways, minSize = 10, maxSize = 500)
h_ssgsea <- gsva(h_param, verbose = FALSE)
cat(sprintf("  ssGSEA result: %d pathways × %d samples\n", nrow(h_ssgsea), ncol(h_ssgsea)))

saveRDS(h_ssgsea, file.path(RES, "ssgsea_scores_human.rds"))

# ============================================================
#  Mouse ssGSEA
# ============================================================
cat("\n=== Mouse ssGSEA ===\n")
m_logcpm <- readRDS(file.path(MOUSE_UI, "results/corrected_logcpm.rds"))
m_meta   <- as.data.table(readRDS(file.path(MOUSE_UI, "results/meta_matched.rds")))

# Map ENSMUSG to symbol
m_annot <- fread(file.path(ANNOT, "mouse_ensmusg_to_symbol.tsv"))
m_annot_base <- gsub("\\..*", "", m_annot$gene_id)

rownames_mbase <- gsub("\\..*", "", rownames(m_logcpm))
sym_map_m <- m_annot$symbol[match(rownames_mbase, m_annot_base)]
valid_m <- !is.na(sym_map_m) & sym_map_m != ""

# Aggregate duplicate symbols by averaging
m_mat <- m_logcpm[valid_m, ]
rownames(m_mat) <- sym_map_m[valid_m]
if (any(duplicated(rownames(m_mat)))) {
  dup_syms <- unique(rownames(m_mat)[duplicated(rownames(m_mat))])
  uniq_rows <- which(!rownames(m_mat) %in% dup_syms)
  avg_list <- lapply(dup_syms, function(sym) {
    colMeans(m_mat[rownames(m_mat) == sym, , drop = FALSE])
  })
  avg_mat <- do.call(rbind, avg_list)
  rownames(avg_mat) <- dup_syms
  m_mat <- rbind(m_mat[uniq_rows, , drop = FALSE], avg_mat)
  cat(sprintf("  Averaged %d duplicate symbols\n", length(dup_syms)))
}

cat(sprintf("  Expression matrix: %d genes × %d samples\n", nrow(m_mat), ncol(m_mat)))

cat("  Running ssGSEA...\n")
m_param <- ssgseaParam(m_mat, m_pathways, minSize = 10, maxSize = 500)
m_ssgsea <- gsva(m_param, verbose = FALSE)
cat(sprintf("  ssGSEA result: %d pathways × %d samples\n", nrow(m_ssgsea), ncol(m_ssgsea)))

saveRDS(m_ssgsea, file.path(RES, "ssgsea_scores_mouse.rds"))

# ============================================================
#  Concordance: compare mean pathway activity disease-vs-ctrl
# ============================================================
cat("\n=== Computing ssGSEA Concordance ===\n")

# Human: disease vs control mean scores
h_meta[, is_disease := condition %in% c("NAFL", "NASH", "NASH_Fibrosis")]
h_disease_ids <- h_meta[is_disease == TRUE, sample_id]
h_control_ids <- h_meta[is_disease == FALSE, sample_id]

h_disease_ids <- intersect(h_disease_ids, colnames(h_ssgsea))
h_control_ids <- intersect(h_control_ids, colnames(h_ssgsea))

h_diff <- rowMeans(h_ssgsea[, h_disease_ids]) - rowMeans(h_ssgsea[, h_control_ids])

# Mouse: per-diet disease vs control
ssgsea_conc <- data.table()
for (diet in DIETS) {
  # Get disease samples for this diet
  m_disease_ids <- intersect(
    m_meta[diet_model == diet & group_binary == "Disease", sample_id],
    colnames(m_ssgsea))
  # Get matching control samples (from datasets that contain this diet)
  diet_datasets <- m_meta[diet_model == diet, unique(dataset)]
  m_control_ids <- intersect(
    m_meta[dataset %in% diet_datasets & group_binary == "Control", sample_id],
    colnames(m_ssgsea))

  if (length(m_disease_ids) < 2 || length(m_control_ids) < 2) {
    cat(sprintf("  %s: insufficient samples (%d disease, %d control), skipping\n",
      diet, length(m_disease_ids), length(m_control_ids)))
    next
  }

  m_diff <- rowMeans(m_ssgsea[, m_disease_ids, drop = FALSE]) -
            rowMeans(m_ssgsea[, m_control_ids, drop = FALSE])

  # Align pathways
  common <- intersect(names(h_diff), names(m_diff))
  rho <- cor(h_diff[common], m_diff[common], method = "spearman")

  row <- data.table(
    diet = diet,
    n_pathways = length(common),
    rho_ssgsea = round(rho, 4),
    n_disease = length(m_disease_ids),
    n_control = length(m_control_ids)
  )
  ssgsea_conc <- rbindlist(list(ssgsea_conc, row))

  cat(sprintf("  %s: ρ_ssGSEA=%.3f (%d pathways, %d vs %d samples)\n",
    diet, rho, length(common), length(m_disease_ids), length(m_control_ids)))
}

fwrite(ssgsea_conc, file.path(RES, "ssgsea_concordance.csv"))

# ============================================================
#  Severity-matched ssGSEA: per-signature contrasts (Component E)
# ============================================================
cat("\n=== Severity-Matched ssGSEA Concordance ===\n")

human_contrasts <- list(
  nafl_vs_nash = list(
    disease = c("NASH", "NASH_Fibrosis"),
    control = c("NAFL"),
    label   = "NAFL-vs-NASH"
  ),
  nafl_specific = list(
    disease = c("NAFL"),
    control = c("Control", "Control_Obese"),
    label   = "NAFL-specific"
  )
  # fibrosis: continuous → skip for ssGSEA binary contrast
  # disease_vs_ctrl: already computed above
)

ssgsea_severity <- data.table()
for (contrast_name in names(human_contrasts)) {
  cdef <- human_contrasts[[contrast_name]]
  h_dis_ids <- intersect(h_meta[condition %in% cdef$disease, sample_id], colnames(h_ssgsea))
  h_ctl_ids <- intersect(h_meta[condition %in% cdef$control, sample_id], colnames(h_ssgsea))

  if (length(h_dis_ids) < 2 || length(h_ctl_ids) < 2) {
    cat(sprintf("  %s: insufficient human samples (%d vs %d), skipping\n",
      cdef$label, length(h_dis_ids), length(h_ctl_ids)))
    next
  }

  h_diff <- rowMeans(h_ssgsea[, h_dis_ids, drop = FALSE]) -
            rowMeans(h_ssgsea[, h_ctl_ids, drop = FALSE])

  for (diet in DIETS) {
    m_disease_ids <- intersect(
      m_meta[diet_model == diet & group_binary == "Disease", sample_id],
      colnames(m_ssgsea))
    diet_datasets <- m_meta[diet_model == diet, unique(dataset)]
    m_control_ids <- intersect(
      m_meta[dataset %in% diet_datasets & group_binary == "Control", sample_id],
      colnames(m_ssgsea))

    if (length(m_disease_ids) < 2 || length(m_control_ids) < 2) next

    m_diff <- rowMeans(m_ssgsea[, m_disease_ids, drop = FALSE]) -
              rowMeans(m_ssgsea[, m_control_ids, drop = FALSE])

    common <- intersect(names(h_diff), names(m_diff))
    rho <- cor(h_diff[common], m_diff[common], method = "spearman")

    row <- data.table(
      human_signature = contrast_name,
      diet = diet,
      n_pathways = length(common),
      rho_ssgsea = round(rho, 4),
      n_human_disease = length(h_dis_ids),
      n_human_control = length(h_ctl_ids),
      n_mouse_disease = length(m_disease_ids),
      n_mouse_control = length(m_control_ids)
    )
    ssgsea_severity <- rbindlist(list(ssgsea_severity, row))

    cat(sprintf("  %s × %s: ρ=%.3f (%d pathways)\n",
      cdef$label, diet, rho, length(common)))
  }
}

# Merge pooled + severity-matched
ssgsea_conc[, human_signature := "disease_vs_ctrl"]
ssgsea_combined <- rbindlist(list(
  ssgsea_conc[, .(human_signature, diet, n_pathways, rho_ssgsea)],
  ssgsea_severity[, .(human_signature, diet, n_pathways, rho_ssgsea)]
), fill = TRUE)
fwrite(ssgsea_combined, file.path(RES, "ssgsea_concordance.csv"))

cat("\nSaved: ssgsea_concordance.csv (with severity-matched contrasts)\n")
cat("=== Phase 2b complete ===\n")
