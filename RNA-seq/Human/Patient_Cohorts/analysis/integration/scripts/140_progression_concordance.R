#!/usr/bin/env Rscript
# 140_progression_concordance.R
# ---------------------------------------------------------------------------
# Cross-species concordance for progression contrasts C3/C5/C6/C8/C9.
#
# Extends the concordance analysis from Script 01 (which covers C1/C2/C4)
# to the new progression contrasts produced by Script 130:
#   C3: F3-F4 vs F0-F2 (advanced vs early fibrosis)
#   C5: NAS >= 5 vs NAS < 5 (clinical NASH threshold)
#   C6: NASH + F3-F4 vs NAFL + F0-F1 (extreme endpoints)
#   C8: F4 vs F0-F3 (cirrhosis binary)
#   C9: F2-F4 vs F0-F1 (F2 inflection point)
#
# For each human contrast x mouse diet pair, computes:
#   - Spearman rho of logFC (all shared orthologs)
#   - Spearman rho of logFC (both-significant genes only)
#   - Direction concordance (% of both-sig genes with same logFC sign)
#   - Jaccard overlap at padj < 0.05
#   - Jaccard top-500 and top-1000 by |t-statistic|
#   - Per-gene concordance flags
#
# Mouse diets: MCD, HFD, CDAHFD, FPC, LIDPAD
#
# Ortholog mapping: same one-to-one ortholog table used by Script 01
# (gene_annotation/ortholog_mapping.tsv, Ensembl biomaRt).
#
# Concordance threshold: padj < 0.05 (same as Script 01, stricter than
# library design padj < 0.1 to reduce false-positive concordance inflation).
#
# Output (to results/progression/):
#   progression_concordance_all.csv     — per-gene concordance for all pairs
#   progression_concordance_summary.csv — 5x5 concordance metrics matrix
#
# Prerequisites: Script 130 must have run (produces c*_dream.csv files).
# Dependencies: data.table
#
# Usage: Rscript 140_progression_concordance.R
# SLURM: cpu, 4 CPU, 16GB RAM, ~5min
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Script 140: Progression Cross-Species Concordance ===\n\n")

# ============================================================
#  Paths
# ============================================================
BASE     <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H_INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PROG_DIR <- file.path(H_INT, "results/progression")
ANNOT    <- file.path(H_INT, "results/gene_annotation")
MOUSE_PD <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")

ODIR <- PROG_DIR
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC", "LIDPAD")

# Match Script 01 threshold
PADJ_THRESH <- 0.05

# ============================================================
#  Load ortholog mapping (same file as Script 01)
# ============================================================
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
cat(sprintf("One-to-one orthologs: %d\n\n", nrow(ortho)))

# ============================================================
#  Load human progression signatures (C3, C5, C6, C8, C9)
# ============================================================
cat("=== Loading Human Progression Signatures ===\n")

# Contrast registry: ID -> filename, description
contrast_registry <- list(
  c3_adv_vs_early_fib = list(
    file = "c3_adv_vs_early_fib_dream.csv",
    desc = "C3: F3-F4 vs F0-F2 (Advanced vs Early Fibrosis)"
  ),
  c5_nas_ge5_vs_lt5 = list(
    file = "c5_nas_ge5_vs_lt5_dream.csv",
    desc = "C5: NAS >= 5 vs NAS < 5 (Clinical NASH Threshold)"
  ),
  c6_extreme_endpoints = list(
    file = "c6_extreme_endpoints_dream.csv",
    desc = "C6: NASH+F3-F4 vs NAFL+F0-F1 (Extreme Endpoints)"
  ),
  c8_cirrhosis = list(
    file = "c8_cirrhosis_dream.csv",
    desc = "C8: F4 vs F0-F3 (Cirrhosis Binary)"
  ),
  c9_f2_inflection = list(
    file = "c9_f2_inflection_dream.csv",
    desc = "C9: F2-F4 vs F0-F1 (F2 Inflection Point)"
  )
)

# Load each signature with consistent column names
load_human_sig <- function(path, label) {
  if (!file.exists(path)) {
    cat(sprintf("  WARNING: %s not found — skipping %s\n", basename(path), label))
    return(NULL)
  }
  dt <- fread(path)
  # Strip version suffix from gene IDs for ortholog matching
  dt[, gene_base := gsub("\\..*", "", gene)]

  # Script 130 outputs padj (renamed from adj.P.Val). Handle both.
  if (!"padj" %in% names(dt) && "adj.P.Val" %in% names(dt)) {
    setnames(dt, "adj.P.Val", "padj")
  }

  n_deg <- sum(dt$padj < PADJ_THRESH, na.rm = TRUE)
  cat(sprintf("  %s: %d genes, %d DEGs (padj<%.2f)\n",
    label, nrow(dt), n_deg, PADJ_THRESH))
  return(dt)
}

human_sigs <- list()
for (cid in names(contrast_registry)) {
  info <- contrast_registry[[cid]]
  path <- file.path(PROG_DIR, info$file)
  sig  <- load_human_sig(path, info$desc)
  if (!is.null(sig)) human_sigs[[cid]] <- sig
}

if (length(human_sigs) == 0) {
  stop("No human progression signatures found. Run Script 130 first.")
}
cat(sprintf("\nLoaded %d / %d progression signatures.\n\n", length(human_sigs),
  length(contrast_registry)))

# ============================================================
#  Load 5 mouse per-diet DE results
# ============================================================
cat("=== Loading Mouse Per-Diet DE ===\n")

mouse_diets <- list()
for (diet in DIETS) {
  f <- file.path(MOUSE_PD, paste0(diet, "_de_results.csv"))
  if (!file.exists(f)) {
    cat(sprintf("  WARNING: %s not found — skipping\n", basename(f)))
    next
  }
  dt <- fread(f)
  dt[, mouse_base := gsub("\\..*", "", gene)]
  mouse_diets[[diet]] <- dt
  cat(sprintf("  %s: %d genes, %d DEGs (padj<%.2f)\n",
    diet, nrow(dt), sum(dt$adj.P.Val < PADJ_THRESH), PADJ_THRESH))
}
cat("\n")

# ============================================================
#  Compute concordance for all N_sig x 5 comparisons
# ============================================================
cat(sprintf("=== Computing %d x %d Concordance Matrix ===\n",
  length(human_sigs), length(mouse_diets)))

results_all  <- data.table()
per_gene_all <- data.table()

for (sig_name in names(human_sigs)) {
  hsig <- human_sigs[[sig_name]]

  # Map human genes to mouse orthologs
  h_mapped <- merge(
    hsig[, .(gene_base, h_lfc = logFC, h_padj = padj, h_t = t)],
    ortho[, .(human_gene_id, mouse_gene_id, human_symbol, mouse_symbol)],
    by.x = "gene_base", by.y = "human_gene_id"
  )

  for (diet in names(mouse_diets)) {
    dt_mouse <- mouse_diets[[diet]]

    # Merge human-mapped orthologs with mouse DE results
    paired <- merge(
      h_mapped,
      dt_mouse[, .(mouse_base, m_lfc = logFC, m_padj = adj.P.Val, m_t = t)],
      by.x = "mouse_gene_id", by.y = "mouse_base"
    )

    n_paired <- nrow(paired)
    if (n_paired == 0) {
      cat(sprintf("  %s x %s: 0 paired genes — skipping\n", sig_name, diet))
      next
    }

    # --- Spearman rho across all shared orthologs ---
    rho_all <- cor(paired$h_lfc, paired$m_lfc, method = "spearman",
      use = "complete.obs")

    # --- Both-significant subset ---
    both_sig   <- paired[h_padj < PADJ_THRESH & m_padj < PADJ_THRESH]
    concordant <- both_sig[sign(h_lfc) == sign(m_lfc)]
    discordant <- both_sig[sign(h_lfc) != sign(m_lfc)]
    conc_rate  <- if (nrow(both_sig) > 0) {
      nrow(concordant) / nrow(both_sig) * 100
    } else NA_real_

    # --- Spearman rho on significant genes only ---
    rho_sig <- if (nrow(both_sig) > 10) {
      cor(both_sig$h_lfc, both_sig$m_lfc, method = "spearman",
        use = "complete.obs")
    } else NA_real_

    # --- Jaccard index at significance threshold ---
    h_sig_ids <- paired[h_padj < PADJ_THRESH, mouse_gene_id]
    m_sig_ids <- paired[m_padj < PADJ_THRESH, mouse_gene_id]
    n_union   <- length(union(h_sig_ids, m_sig_ids))
    jaccard   <- if (n_union > 0) {
      length(intersect(h_sig_ids, m_sig_ids)) / n_union
    } else NA_real_

    # --- Top-N Jaccard (by |t-statistic|) ---
    j500  <- NA_real_
    j1000 <- NA_real_
    for (topN in c(500, 1000)) {
      n_use  <- min(topN, n_paired)
      h_topN <- paired[order(-abs(h_t))][seq_len(n_use), mouse_gene_id]
      m_topN <- paired[order(-abs(m_t))][seq_len(n_use), mouse_gene_id]
      n_u    <- length(union(h_topN, m_topN))
      jval   <- if (n_u > 0) length(intersect(h_topN, m_topN)) / n_u else NA_real_
      if (topN == 500)  j500  <- jval
      if (topN == 1000) j1000 <- jval
    }

    # --- Assemble summary row ---
    row <- data.table(
      human_signature = sig_name,
      diet            = diet,
      n_paired        = n_paired,
      rho_all         = round(rho_all, 4),
      rho_sig         = round(rho_sig, 4),
      n_both_sig      = nrow(both_sig),
      n_concordant    = nrow(concordant),
      n_discordant    = nrow(discordant),
      concordance_pct = round(conc_rate, 1),
      jaccard         = round(jaccard, 4),
      jaccard_top500  = round(j500, 4),
      jaccard_top1000 = round(j1000, 4)
    )
    results_all <- rbindlist(list(results_all, row))

    # --- Per-gene concordance table ---
    paired[, `:=`(
      sig_name   = sig_name,
      diet_name  = diet,
      both_sig   = h_padj < PADJ_THRESH & m_padj < PADJ_THRESH,
      concordant = h_padj < PADJ_THRESH & m_padj < PADJ_THRESH &
                   sign(h_lfc) == sign(m_lfc),
      discordant = h_padj < PADJ_THRESH & m_padj < PADJ_THRESH &
                   sign(h_lfc) != sign(m_lfc)
    )]
    per_gene_all <- rbindlist(list(per_gene_all,
      paired[, .(mouse_gene_id, human_symbol, mouse_symbol, sig_name, diet_name,
                  h_lfc, m_lfc, h_padj, m_padj, both_sig, concordant, discordant)]),
      fill = TRUE)

    cat(sprintf("  %s x %s: n=%d, rho=%.3f, rho_sig=%s, conc=%d/%d (%.1f%%)\n",
      sig_name, diet, n_paired, rho_all,
      ifelse(is.na(rho_sig), "NA", sprintf("%.3f", rho_sig)),
      nrow(concordant), nrow(both_sig),
      ifelse(is.na(conc_rate), 0, conc_rate)))
  }
  cat("\n")
}

# ============================================================
#  Per-gene concordance classification
#  (Same classification scheme as Script 01)
# ============================================================
cat("=== Per-Gene Concordance Classification ===\n")

classify_genes_progression <- function(pg, sig) {
  sg <- pg[sig_name == sig]
  gs <- sg[, .(
    h_significant    = any(h_padj < PADJ_THRESH),
    n_mouse_sig      = sum(m_padj < PADJ_THRESH),
    n_diets_sig      = sum(both_sig),
    n_concordant     = sum(concordant),
    n_discordant     = sum(discordant),
    mean_h_lfc       = mean(h_lfc),
    diets_concordant = paste(diet_name[concordant == TRUE], collapse = ";"),
    diets_discordant = paste(diet_name[discordant == TRUE], collapse = ";"),
    diets_mouse_sig  = paste(diet_name[m_padj < PADJ_THRESH], collapse = ";")
  ), by = .(mouse_gene_id, human_symbol)]

  # Mutually exclusive classification (same logic as Script 01)
  gs[, category := fcase(
    n_concordant >= 3,                                  "Conserved",
    n_discordant > n_concordant & n_diets_sig >= 2,     "Species_Discordant",
    n_concordant %in% 1:2 & n_discordant <= 1,          "Moderate_Concordance",
    n_mouse_sig == 1 & n_concordant <= 1,               "Diet_Selective",
    h_significant & n_mouse_sig == 0,                   "Human_Enriched",
    !h_significant & n_mouse_sig >= 2,                  "Mouse_Specific",
    !h_significant & n_mouse_sig == 0,                  "Not_Significant",
    default = "Unclassified"
  )]
  gs[, signature := sig]
  return(gs)
}

all_sigs <- unique(per_gene_all$sig_name)
all_classified <- rbindlist(lapply(all_sigs, function(s) {
  classify_genes_progression(per_gene_all, s)
}))

cat("Per-signature category counts:\n")
print(all_classified[, .N, by = .(signature, category)][order(signature, -N)])

# ============================================================
#  Save outputs
# ============================================================

# 1. Per-gene concordance (all pairs)
fwrite(per_gene_all, file.path(ODIR, "progression_concordance_all.csv"))
cat(sprintf("\nSaved: progression_concordance_all.csv (%d rows)\n", nrow(per_gene_all)))

# 1b. Per-gene classified concordance (with category labels)
fwrite(all_classified, file.path(ODIR, "progression_concordance_classified.csv"))
cat(sprintf("Saved: progression_concordance_classified.csv (%d rows)\n", nrow(all_classified)))

# 2. Summary concordance matrix
fwrite(results_all, file.path(ODIR, "progression_concordance_summary.csv"))
cat(sprintf("Saved: progression_concordance_summary.csv (%d rows)\n", nrow(results_all)))

# ============================================================
#  Translatability rankings: which mouse model best matches
#  human PROGRESSION signatures?
# ============================================================
cat("\n=== TRANSLATABILITY RANKINGS (Progression Contrasts) ===\n")

for (sig in unique(results_all$human_signature)) {
  desc <- contrast_registry[[sig]]$desc
  cat(sprintf("\n  %s\n", desc))
  sub <- results_all[human_signature == sig][order(-rho_all)]
  print(sub[, .(diet, rho_all, rho_sig, concordance_pct, jaccard, jaccard_top500)])
}

# ============================================================
#  Best mouse model per contrast (highest rho_all)
# ============================================================
cat("\n=== Best Mouse Model per Progression Contrast ===\n")
best_per_contrast <- results_all[, .SD[which.max(rho_all)],
  by = human_signature][, .(human_signature, diet, rho_all, concordance_pct, jaccard)]
print(best_per_contrast)

# ============================================================
#  Comparison with existing contrasts (C1, C2, C4 from Script 01)
# ============================================================
cat("\n=== Comparison with Existing Concordance (Script 01) ===\n")
existing_file <- file.path(BASE, "Analysis/Cross_Species_Concordance/results",
  "gene_concordance_matrix_20x.csv")
if (file.exists(existing_file)) {
  existing <- fread(existing_file)
  cat("Script 01 concordance (existing contrasts):\n")
  existing_best <- existing[, .SD[which.max(rho_all)],
    by = human_signature][, .(human_signature, diet, rho_all, concordance_pct, jaccard)]
  print(existing_best)

  # Combined ranking across all contrasts
  cat("\nCombined ranking (all contrasts, best mouse model per contrast):\n")
  combined <- rbindlist(list(
    existing_best[, .(source = "Script01_original", human_signature, diet,
                      rho_all, concordance_pct, jaccard)],
    best_per_contrast[, .(source = "Script140_progression", human_signature, diet,
                          rho_all, concordance_pct, jaccard)]
  ))
  print(combined[order(-rho_all)])

  # Overall best mouse model (mean rho across all progression contrasts)
  cat("\nOverall mouse model ranking (mean rho across progression contrasts):\n")
  diet_ranking <- results_all[, .(
    mean_rho = round(mean(rho_all, na.rm = TRUE), 4),
    mean_concordance = round(mean(concordance_pct, na.rm = TRUE), 1),
    mean_jaccard = round(mean(jaccard, na.rm = TRUE), 4),
    n_contrasts = .N
  ), by = diet][order(-mean_rho)]
  print(diet_ranking)
} else {
  cat("  Script 01 output not found — skipping comparison.\n")

  # Still rank mouse models based on progression contrasts alone
  cat("\nOverall mouse model ranking (mean rho across progression contrasts):\n")
  diet_ranking <- results_all[, .(
    mean_rho = round(mean(rho_all, na.rm = TRUE), 4),
    mean_concordance = round(mean(concordance_pct, na.rm = TRUE), 1),
    mean_jaccard = round(mean(jaccard, na.rm = TRUE), 4),
    n_contrasts = .N
  ), by = diet][order(-mean_rho)]
  print(diet_ranking)
}

# ============================================================
#  Conserved gene counts per progression contrast
# ============================================================
cat("\n=== Conserved Genes per Progression Contrast ===\n")
cc_counts <- all_classified[category == "Conserved", .N, by = signature][order(-N)]
print(cc_counts)

# Cross-reference with Script 01 Conserved
existing_cc_file <- file.path(BASE, "Analysis/Cross_Species_Concordance/results",
  "gene_concordance_per_gene.csv")
if (file.exists(existing_cc_file)) {
  existing_cc <- fread(existing_cc_file)
  orig_cc_genes <- existing_cc[primary_category == "Conserved" |
                                dvc_category == "Conserved", mouse_gene_id]

  for (sig in unique(all_classified$signature)) {
    prog_cc <- all_classified[signature == sig & category == "Conserved",
      mouse_gene_id]
    n_overlap <- length(intersect(prog_cc, orig_cc_genes))
    n_novel   <- length(setdiff(prog_cc, orig_cc_genes))
    cat(sprintf("  %s: %d CC genes (%d overlap with Script 01, %d novel)\n",
      sig, length(prog_cc), n_overlap, n_novel))
  }
}

cat("\n=== Script 140 complete ===\n")
