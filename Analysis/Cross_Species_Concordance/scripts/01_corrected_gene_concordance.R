#!/usr/bin/env Rscript
# 01_corrected_gene_concordance.R
# ---------------------------------------------------------------------------
# Corrected cross-species concordance: 4 human signatures × 5 mouse diets
# Fixes severity-mismatch bias by using severity-matched comparisons
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

pdf.options(useDingbats = FALSE)

cat("=== Phase 1: Corrected Gene-Level Concordance ===\n\n")

# ============================================================
#  Paths
# ============================================================
BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT   <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT   <- file.path(H_INT, "results/gene_annotation")
DS_DIR  <- file.path(H_INT, "results/disease_signatures")
INT_DIR <- file.path(H_INT, "results/integration")
MOUSE_PD <- file.path(BASE, "Mouse/Unified_Integration/results/per_diet")
WD      <- file.path(BASE, "Analysis/Cross_Species_Concordance")
RES     <- file.path(WD, "results")
PLOTS   <- file.path(WD, "plots")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)

# 2026-05-29: LIDPAD dropped — its per-diet DE was archived (results/per_diet/archive_dropped_diets/);
# concordance now uses the 4 surviving originally-intended diets. (NASH_diet/Western_diet exist in
# the active set but were never part of this concordance; left out to preserve the original analysis.)
DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")

# Concordance significance threshold
# NOTE: Primary DEGs defined by padj < 0.05, |logFC| > 0.3.
# For cross-species concordance we use padj < 0.05 for both human and mouse.
PADJ_THRESH <- 0.05

# ============================================================
#  Load ortholog mapping
# ============================================================
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
cat("One-to-one orthologs:", nrow(ortho), "\n\n")

# ============================================================
#  Load 4 human signatures
# ============================================================
cat("=== Loading Human Signatures ===\n")

load_human <- function(path, label, padj_col = "adj.P.Val") {
  dt <- fread(path)
  dt[, gene_base := gsub("\\..*", "", gene)]
  # Always use raw logFC for symmetric comparison with mouse (raw logFC).
  # If ashr columns are present, keep them as auxiliary but do NOT substitute.
  if ("shrunk_logFC" %in% names(dt)) {
    cat(sprintf("  NOTE: ashr columns present in %s — using raw logFC (not shrunk_logFC) for symmetry with mouse\n", label))
  }
  # Normalize padj column
  if (padj_col %in% names(dt) && padj_col != "padj") {
    setnames(dt, padj_col, "padj", skip_absent = TRUE)
  }
  if (!"padj" %in% names(dt) && "adj.P.Val" %in% names(dt)) {
    setnames(dt, "adj.P.Val", "padj")
  }
  cat(sprintf("  %s: %d genes, %d DEGs (padj<0.05)\n",
    label, nrow(dt), sum(dt$padj < 0.05, na.rm = TRUE)))
  return(dt)
}

# Use RAW (unshrunk) logFC for disease-vs-control, symmetric with mouse raw logFC.
# Canonical hard cutover (2026-06-08): canonical_deg_results.csv carries the raw
# `logFC` column (plus `shrunk_logFC` as auxiliary). load_human() deliberately keeps
# raw logFC — shrinking only the human side would bias the Spearman rho.
disease_dream_path <- file.path(INT_DIR, "canonical_deg_results.csv")
if (!file.exists(disease_dream_path)) {
  stop("canonical_deg_results.csv not found at: ", disease_dream_path)
}
human_sigs <- list(
  disease_vs_ctrl = load_human(disease_dream_path,
                                "Disease-vs-Control", "padj"),
  nafl_specific   = load_human(file.path(RES, "nafl_vs_ctrl_dream.csv"),
                                "NAFL-vs-Control", "adj.P.Val"),
  nafl_vs_nash    = load_human(file.path(DS_DIR, "nafl_vs_nash_dream.csv"),
                                "NAFL-vs-NASH", "adj.P.Val"),
  fibrosis        = load_human(file.path(DS_DIR, "fibrosis_dream.csv"),
                                "Fibrosis slope", "adj.P.Val")
)

# ============================================================
#  Load 5 mouse per-diet DE results
# ============================================================
cat("\n=== Loading Mouse Per-Diet DE ===\n")

mouse_diets <- list()
for (diet in DIETS) {
  f <- file.path(MOUSE_PD, paste0(diet, "_de_results.csv"))
  dt <- fread(f)
  dt[, mouse_base := gsub("\\..*", "", gene)]
  mouse_diets[[diet]] <- dt
  cat(sprintf("  %s: %d genes, %d DEGs (padj<0.05)\n",
    diet, nrow(dt), sum(dt$adj.P.Val < 0.05)))
}

# ============================================================
#  Compute concordance for all 4×5 = 20 comparisons
# ============================================================
cat("\n=== Computing 4×5 Concordance Matrix ===\n")

results_all <- data.table()
per_gene_all <- data.table()

for (sig_name in names(human_sigs)) {
  hsig <- human_sigs[[sig_name]]

  # Map human → mouse via orthologs
  h_mapped <- merge(
    hsig[, .(gene_base, h_lfc = logFC, h_padj = padj, h_t = t)],
    ortho[, .(human_gene_id, mouse_gene_id, human_symbol, mouse_symbol)],
    by.x = "gene_base", by.y = "human_gene_id"
  )

  for (diet in DIETS) {
    dt_mouse <- mouse_diets[[diet]]

    paired <- merge(
      h_mapped,
      dt_mouse[, .(mouse_base, m_lfc = logFC, m_padj = adj.P.Val, m_t = t)],
      by.x = "mouse_gene_id", by.y = "mouse_base"
    )

    n_paired <- nrow(paired)
    rho_all <- cor(paired$h_lfc, paired$m_lfc, method = "spearman", use = "complete.obs")

    both_sig <- paired[h_padj < PADJ_THRESH & m_padj < PADJ_THRESH]
    concordant <- both_sig[sign(h_lfc) == sign(m_lfc)]
    discordant <- both_sig[sign(h_lfc) != sign(m_lfc)]
    conc_rate <- if (nrow(both_sig) > 0) nrow(concordant) / nrow(both_sig) * 100 else NA

    rho_sig <- if (nrow(both_sig) > 10) {
      cor(both_sig$h_lfc, both_sig$m_lfc, method = "spearman", use = "complete.obs")
    } else NA

    # Jaccard
    h_sig_ids <- paired[h_padj < PADJ_THRESH, mouse_gene_id]
    m_sig_ids <- paired[m_padj < PADJ_THRESH, mouse_gene_id]
    jaccard <- length(intersect(h_sig_ids, m_sig_ids)) / length(union(h_sig_ids, m_sig_ids))

    # Top-N Jaccard
    for (topN in c(500, 1000)) {
      h_topN <- paired[order(-abs(h_t))][1:min(topN, n_paired), mouse_gene_id]
      m_topN <- paired[order(-abs(m_t))][1:min(topN, n_paired), mouse_gene_id]
      assign(paste0("j", topN), length(intersect(h_topN, m_topN)) / length(union(h_topN, m_topN)))
    }

    row <- data.table(
      human_signature = sig_name,
      diet = diet,
      n_paired = n_paired,
      rho_all = round(rho_all, 4),
      rho_sig = round(rho_sig, 4),
      n_both_sig = nrow(both_sig),
      n_concordant = nrow(concordant),
      n_discordant = nrow(discordant),
      concordance_pct = round(conc_rate, 1),
      jaccard = round(jaccard, 4),
      jaccard_top500 = round(j500, 4),
      jaccard_top1000 = round(j1000, 4)
    )
    results_all <- rbindlist(list(results_all, row))

    # Per-gene concordance for this comparison
    paired[, `:=`(
      sig_name = sig_name,
      diet_name = diet,
      both_sig = h_padj < PADJ_THRESH & m_padj < PADJ_THRESH,
      concordant = h_padj < PADJ_THRESH & m_padj < PADJ_THRESH & sign(h_lfc) == sign(m_lfc),
      discordant = h_padj < PADJ_THRESH & m_padj < PADJ_THRESH & sign(h_lfc) != sign(m_lfc)
    )]
    per_gene_all <- rbindlist(list(per_gene_all,
      paired[, .(mouse_gene_id, human_symbol, sig_name, diet_name,
                  h_lfc, m_lfc, h_padj, m_padj, both_sig, concordant, discordant)]),
      fill = TRUE)

    cat(sprintf("  %s × %s: ρ=%.3f, ρ_sig=%.3f, conc=%d/%d (%.1f%%)\n",
      sig_name, diet, rho_all, ifelse(is.na(rho_sig), 0, rho_sig),
      nrow(concordant), nrow(both_sig), ifelse(is.na(conc_rate), 0, conc_rate)))
  }
  cat("\n")
}

fwrite(results_all, file.path(RES, "gene_concordance_matrix_20x.csv"))
cat("Saved: gene_concordance_matrix_20x.csv\n")

# ============================================================
#  Per-gene concordance: ALL signatures (Components A+B)
# ============================================================
cat("\n=== Per-Gene Concordance Summary (All Signatures) ===\n")

# Helper: classify genes for one signature
classify_genes <- function(pg, sig) {
  sg <- pg[sig_name == sig]
  gs <- sg[, .(
    h_significant = any(h_padj < PADJ_THRESH),
    n_mouse_sig   = sum(m_padj < PADJ_THRESH),
    n_diets_sig   = sum(both_sig),
    n_concordant  = sum(concordant),
    n_discordant  = sum(discordant),
    mean_h_lfc    = mean(h_lfc),
    diets_concordant  = paste(diet_name[concordant == TRUE], collapse = ";"),
    diets_discordant  = paste(diet_name[discordant == TRUE], collapse = ";"),
    diets_mouse_sig   = paste(diet_name[m_padj < PADJ_THRESH], collapse = ";")
  ), by = .(mouse_gene_id, human_symbol)]

  # Mutually exclusive classification via fcase (first match wins)
  gs[, category := fcase(
    n_concordant >= 3,                                       "Conserved",
    n_discordant > n_concordant & n_diets_sig >= 2,          "Species_Discordant",
    n_concordant %in% 1:2 & n_discordant <= 1,               "Moderate_Concordance",
    n_mouse_sig == 1 & n_concordant <= 1,                    "Diet_Selective",
    h_significant & n_mouse_sig == 0,                        "Human_Enriched",
    !h_significant & n_mouse_sig >= 2,                       "Mouse_Specific",
    !h_significant & n_mouse_sig == 0,                       "Not_Significant",
    default = "Unclassified"
  )]
  gs[, signature := sig]
  return(gs)
}

# Run classification for each signature
all_sigs <- unique(per_gene_all$sig_name)
all_summaries <- rbindlist(lapply(all_sigs, function(s) classify_genes(per_gene_all, s)))

# Save full per-signature classification
fwrite(all_summaries, file.path(RES, "gene_concordance_all_signatures.csv"))
cat("Per-signature category counts:\n")
print(all_summaries[, .N, by = .(signature, category)][order(signature, -N)])

# -----------------------------------------------------------------------
# Primary classification — NAFL-vs-NASH anchor
# Rationale: mouse diet models specifically induce steatohepatitis, so the
# NAFL→NASH progression signal in humans has the highest empirical concordance
# with mouse models (ρ=0.323 for FPC vs ρ=0.178 for disease_vs_ctrl).
# These genes are "progression drivers conserved across species."
# -----------------------------------------------------------------------
primary <- all_summaries[signature == "nafl_vs_nash"]
setnames(primary, "category", "primary_category")
primary[, signature := NULL]

# -----------------------------------------------------------------------
# Secondary classification — disease_vs_ctrl anchor
# Rationale: captures genes broadly elevated in MASLD vs healthy regardless
# of disease stage. A gene can be stably elevated in NAFL AND NASH (not
# change in NAFL→NASH comparison) yet still be a valid Cas13 target.
# These genes are "disease maintenance genes conserved across species."
# -----------------------------------------------------------------------
primary_dvc <- all_summaries[signature == "disease_vs_ctrl", .(
  mouse_gene_id,
  human_symbol,
  dvc_category       = category,
  dvc_n_concordant   = n_concordant,
  dvc_n_discordant   = n_discordant,
  dvc_diets_concordant = diets_concordant,
  dvc_h_significant  = h_significant,
  dvc_mean_h_lfc     = mean_h_lfc
)]

# -----------------------------------------------------------------------
# Best-case across ALL signatures (highest concordant diet count)
# -----------------------------------------------------------------------
best <- all_summaries[, .(
  best_n_concordant = max(n_concordant),
  best_category = category[which.max(n_concordant)],
  best_signature = signature[which.max(n_concordant)],
  n_signatures_concordant = sum(n_concordant > 0),
  all_categories = paste(unique(category), collapse = ";")
), by = .(mouse_gene_id, human_symbol)]

# -----------------------------------------------------------------------
# Merge all three into one summary
# -----------------------------------------------------------------------
gene_summary <- merge(primary, best, by = c("mouse_gene_id", "human_symbol"), all.x = TRUE)
gene_summary <- merge(gene_summary, primary_dvc, by = c("mouse_gene_id", "human_symbol"), all.x = TRUE)

cat("\nPrimary (NAFL-vs-NASH) categories:\n")
print(gene_summary[, .N, by = primary_category][order(-N)])
cat("\nDisease-vs-Control anchor categories:\n")
print(gene_summary[, .N, by = dvc_category][order(-N)])
cat("\nBest-case (across all 4 signatures) categories:\n")
print(gene_summary[, .N, by = best_category][order(-N)])

# Dual-anchor summary
cat("\nDual-anchor Conserved (NAFL-vs-NASH AND disease_vs_ctrl):",
    gene_summary[primary_category == "Conserved" & dvc_category == "Conserved", .N], "\n")

fwrite(gene_summary, file.path(RES, "gene_concordance_per_gene.csv"))
cat("Saved: gene_concordance_per_gene.csv\n")

# ============================================================
#  Summary: Translatability ranking by signature
# ============================================================
cat("\n=== TRANSLATABILITY RANKINGS ===\n")
for (sig in unique(results_all$human_signature)) {
  cat(sprintf("\n  %s:\n", sig))
  sub <- results_all[human_signature == sig][order(-rho_all)]
  print(sub[, .(diet, rho_all, rho_sig, concordance_pct, jaccard_top500)])
}

cat("\n=== Phase 1 complete ===\n")
