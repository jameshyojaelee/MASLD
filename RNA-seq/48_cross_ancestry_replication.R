#!/usr/bin/env Rscript
# 48_cross_ancestry_replication.R
# ---------------------------------------------------------------------------
# Cross-Ancestry Colocalization Replication Analysis
#
# Compares COLOC results across 4 ancestry groups:
#   - European: UKBB ALT/AST/GGT, PDFF, FinnGen NAFLD/NASH/HCC, Broadaway, GTEx
#   - East Asian: BBJ ALT/AST/GGT
#   - African: Pan-UKBB AFR ALT/AST/GGT
#   - Central/South Asian: Pan-UKBB CSA ALT/AST/GGT
# Plus progression-specific GWAS:
#   - Ghouse Cirrhosis/HCC, deCODE NAFL/Cirrhosis/HCC
#
# Analyses:
#   1. Per-ancestry best PP.H4 per gene (EUR, EAS, AFR, CSA)
#   2. PP.H4 correlation: pairwise across ancestries
#   3. Ancestry-specific loci: PP.H4 > 0.5 in one, < 0.1 in others
#   4. Known anchor genes: PNPLA3, TM6SF2, HSD17B13, MBOAT7
#   5. Cross-ancestry validated targets → n_ancestry_sig counting
#   6. Progression-specific analysis: liver enzyme vs cirrhosis/HCC overlap
#
# Input: COLOC results from Scripts 35b/35g/35h, 46, 47, 49, 50
# Output: RNA-seq/results/causal_inference/cross_ancestry/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUTDIR <- file.path(BASE, "RNA-seq/results/causal_inference/cross_ancestry")
FIGDIR <- file.path(OUTDIR, "figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Script 48: Cross-Ancestry COLOC Replication ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Load all COLOC results
# ==============================================================================
cat("--- Step 1: Loading COLOC results ---\n")

load_coloc <- function(path, source_name) {
  if (!file.exists(path)) {
    cat("  ", source_name, ": not found\n")
    return(data.table(gene = character(), PP.H4 = numeric(), source = character()))
  }
  dt <- fread(path)
  if (nrow(dt) == 0) {
    cat("  ", source_name, ": empty\n")
    return(data.table(gene = character(), PP.H4 = numeric(), source = character()))
  }
  # Handle both legacy (PP.H4) and new SuSiE-COLOC (PP.H4.abf, PP.H4.susie) schemas
  # Prefer PP.H4.susie when available, fall back to PP.H4.abf, then PP.H4
  if ("PP.H4.susie" %in% names(dt) && any(!is.na(dt$PP.H4.susie))) {
    # Use SuSiE PP.H4 with ABF fallback for genes where SuSiE failed
    dt[, PP.H4 := fifelse(!is.na(PP.H4.susie), PP.H4.susie, PP.H4.abf)]
    cat("  ", source_name, ": using SuSiE PP.H4 (with ABF fallback)\n")
  } else if ("PP.H4.abf" %in% names(dt)) {
    setnames(dt, "PP.H4.abf", "PP.H4")
  }
  if (!"PP.H4" %in% names(dt)) {
    cat("  ", source_name, ": missing PP.H4 column\n")
    return(data.table(gene = character(), PP.H4 = numeric(), source = character()))
  }
  # Normalize gene column name
  gene_col <- intersect(c("gene", "gene_symbol", "GeneSymbol"), names(dt))[1]
  if (is.na(gene_col)) {
    cat("  ", source_name, ": missing gene column\n")
    return(data.table(gene = character(), PP.H4 = numeric(), source = character()))
  }
  if (gene_col != "gene") setnames(dt, gene_col, "gene")
  dt <- dt[, .(gene, PP.H4)]
  dt[, source := source_name]
  cat("  ", source_name, ":", nrow(dt), "genes,",
      sum(dt$PP.H4 > 0.5, na.rm = TRUE), "with PP.H4 > 0.5\n")
  return(dt)
}

# European COLOC sources
eur_sources <- list(
  "UKBB_ALT"     = file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv"),
  "UKBB_AST"     = file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb_ast/coloc_results.csv"),
  "UKBB_GGT"     = file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb_ggt/coloc_results.csv"),
  "UKBB_PDFF"    = file.path(BASE, "RNA-seq/results/causal_inference/broadaway_pdff/coloc_results.csv"),
  "FinnGen_NAFLD" = file.path(BASE, "RNA-seq/results/causal_inference/finngen_nafld/coloc_results.csv"),
  "FinnGen_NASH"  = file.path(BASE, "RNA-seq/results/causal_inference/finngen_nash/coloc_results.csv"),
  "FinnGen_HCC"   = file.path(BASE, "RNA-seq/results/causal_inference/finngen_hcc/coloc_results.csv"),
  "Broadaway_Ghodsian" = file.path(BASE, "RNA-seq/results/causal_inference/broadaway/coloc_results.csv"),
  "GTEx_Ghodsian" = file.path(BASE, "RNA-seq/results/causal_inference/ghodsian/coloc_results.csv")
)

# European progression COLOC sources
eur_sources[["Ghouse_Cirrhosis"]] <- file.path(BASE, "RNA-seq/results/causal_inference/ghouse_cirrhosis/coloc_results.csv")
eur_sources[["Ghouse_HCC"]]       <- file.path(BASE, "RNA-seq/results/causal_inference/ghouse_hcc/coloc_results.csv")
# deCODE COLOC not available — Sveinbjornsson studies were run through the SuSiE-COLOC
# pipeline (GWAS/finemapping/results/susie_coloc/) but not the legacy ABF pipeline.
# They'll be picked up when Script 48 is updated to read from gene_level_coloc.csv.

# East Asian COLOC sources
eas_sources <- list(
  "BBJ_ALT" = file.path(BASE, "RNA-seq/results/causal_inference/bbj_alt/coloc_results.csv"),
  "BBJ_AST" = file.path(BASE, "RNA-seq/results/causal_inference/bbj_ast/coloc_results.csv"),
  "BBJ_GGT" = file.path(BASE, "RNA-seq/results/causal_inference/bbj_ggt/coloc_results.csv")
)

# African COLOC sources (Pan-UKBB)
afr_sources <- list(
  "PanUKBB_AFR_ALT" = file.path(BASE, "RNA-seq/results/causal_inference/panukbb_afr_alt/coloc_results.csv"),
  "PanUKBB_AFR_AST" = file.path(BASE, "RNA-seq/results/causal_inference/panukbb_afr_ast/coloc_results.csv"),
  "PanUKBB_AFR_GGT" = file.path(BASE, "RNA-seq/results/causal_inference/panukbb_afr_ggt/coloc_results.csv")
)

# Central/South Asian COLOC sources (Pan-UKBB)
csa_sources <- list(
  "PanUKBB_CSA_ALT" = file.path(BASE, "RNA-seq/results/causal_inference/panukbb_csa_alt/coloc_results.csv"),
  "PanUKBB_CSA_AST" = file.path(BASE, "RNA-seq/results/causal_inference/panukbb_csa_ast/coloc_results.csv"),
  "PanUKBB_CSA_GGT" = file.path(BASE, "RNA-seq/results/causal_inference/panukbb_csa_ggt/coloc_results.csv")
)

eur_results <- rbindlist(lapply(names(eur_sources), function(n) load_coloc(eur_sources[[n]], n)), fill = TRUE)
eas_results <- rbindlist(lapply(names(eas_sources), function(n) load_coloc(eas_sources[[n]], n)), fill = TRUE)
afr_results <- rbindlist(lapply(names(afr_sources), function(n) load_coloc(afr_sources[[n]], n)), fill = TRUE)
csa_results <- rbindlist(lapply(names(csa_sources), function(n) load_coloc(csa_sources[[n]], n)), fill = TRUE)

cat("\n  Total EUR results:", nrow(eur_results), "gene-source pairs\n")
cat("  Total EAS results:", nrow(eas_results), "gene-source pairs\n")
cat("  Total AFR results:", nrow(afr_results), "gene-source pairs\n")
cat("  Total CSA results:", nrow(csa_results), "gene-source pairs\n")

# ==============================================================================
# 2. Per-ancestry best PP.H4 per gene
# ==============================================================================
cat("\n--- Step 2: Computing per-ancestry best PP.H4 ---\n")

# European: best PP.H4 across all European COLOC sources
if (nrow(eur_results) > 0) {
  eur_best <- eur_results[, .(
    best_pp4_eur = max(PP.H4, na.rm = TRUE),
    n_eur_sources_sig = sum(PP.H4 > 0.5),
    best_eur_source = source[which.max(PP.H4)]
  ), by = gene]
} else {
  eur_best <- data.table(gene = character(), best_pp4_eur = numeric(),
                          n_eur_sources_sig = integer(), best_eur_source = character())
}

# East Asian: best PP.H4 across all BBJ sources
if (nrow(eas_results) > 0) {
  eas_best <- eas_results[, .(
    best_pp4_eas = max(PP.H4, na.rm = TRUE),
    n_eas_sources_sig = sum(PP.H4 > 0.5),
    best_eas_source = source[which.max(PP.H4)]
  ), by = gene]
} else {
  eas_best <- data.table(gene = character(), best_pp4_eas = numeric(),
                          n_eas_sources_sig = integer(), best_eas_source = character())
}

# African: best PP.H4 across Pan-UKBB AFR sources
if (nrow(afr_results) > 0) {
  afr_best <- afr_results[, .(
    best_pp4_afr = max(PP.H4, na.rm = TRUE),
    n_afr_sources_sig = sum(PP.H4 > 0.5),
    best_afr_source = source[which.max(PP.H4)]
  ), by = gene]
} else {
  afr_best <- data.table(gene = character(), best_pp4_afr = numeric(),
                          n_afr_sources_sig = integer(), best_afr_source = character())
}

# Central/South Asian: best PP.H4 across Pan-UKBB CSA sources
if (nrow(csa_results) > 0) {
  csa_best <- csa_results[, .(
    best_pp4_csa = max(PP.H4, na.rm = TRUE),
    n_csa_sources_sig = sum(PP.H4 > 0.5),
    best_csa_source = source[which.max(PP.H4)]
  ), by = gene]
} else {
  csa_best <- data.table(gene = character(), best_pp4_csa = numeric(),
                          n_csa_sources_sig = integer(), best_csa_source = character())
}

cat("  EUR genes tested:", nrow(eur_best), "\n")
cat("  EUR genes with PP.H4 > 0.5:", sum(eur_best$best_pp4_eur > 0.5), "\n")
cat("  EAS genes tested:", nrow(eas_best), "\n")
cat("  EAS genes with PP.H4 > 0.5:", sum(eas_best$best_pp4_eas > 0.5), "\n")
cat("  AFR genes tested:", nrow(afr_best), "\n")
cat("  AFR genes with PP.H4 > 0.5:", sum(afr_best$best_pp4_afr > 0.5), "\n")
cat("  CSA genes tested:", nrow(csa_best), "\n")
cat("  CSA genes with PP.H4 > 0.5:", sum(csa_best$best_pp4_csa > 0.5), "\n")

# ==============================================================================
# 3. Cross-ancestry comparison
# ==============================================================================
cat("\n--- Step 3: Cross-ancestry comparison (4-way) ---\n")

# Merge all 4 ancestry best tables
cross <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE),
                list(eur_best, eas_best, afr_best, csa_best))
# Track which ancestries were actually tested (NA = untested, 0 = tested but no signal)
for (col in c("best_pp4_eur", "best_pp4_eas", "best_pp4_afr", "best_pp4_csa")) {
  if (!col %in% names(cross)) cross[, (col) := NA_real_]
}
cross[, n_ancestry_tested := (!is.na(best_pp4_eur)) + (!is.na(best_pp4_eas)) +
                              (!is.na(best_pp4_afr)) + (!is.na(best_pp4_csa))]

# Count how many ancestry groups have PP.H4 > 0.5 (NA-safe: untested counts as FALSE)
cross[, n_ancestry_sig := (fifelse(is.na(best_pp4_eur), FALSE, best_pp4_eur > 0.5)) +
                           (fifelse(is.na(best_pp4_eas), FALSE, best_pp4_eas > 0.5)) +
                           (fifelse(is.na(best_pp4_afr), FALSE, best_pp4_afr > 0.5)) +
                           (fifelse(is.na(best_pp4_csa), FALSE, best_pp4_csa > 0.5))]
# Validated ancestry count: EUR+EAS only (AFR/CSA are exploratory due to
# low power and prior-driven posteriors — see 50c prior sensitivity analysis)
cross[, n_validated_ancestry := (best_pp4_eur > 0.5) + (best_pp4_eas > 0.5)]

# Classification — distinguish validated (EUR/EAS) from exploratory (AFR/CSA)
cross[, class := "Neither"]
cross[n_validated_ancestry >= 2, class := "Multi_ancestry_validated"]
cross[n_validated_ancestry == 1 & best_pp4_eur > 0.5, class := "EUR_only"]
cross[n_validated_ancestry == 1 & best_pp4_eas > 0.5, class := "EAS_only"]
# AFR/CSA-only genes (no EUR/EAS signal) are exploratory
cross[n_validated_ancestry == 0 & best_pp4_afr > 0.5, class := "AFR_exploratory"]
cross[n_validated_ancestry == 0 & best_pp4_csa > 0.5, class := "CSA_exploratory"]

cat("\n  Cross-ancestry classification:\n")
class_counts <- cross[, .N, by = class][order(-N)]
for (i in seq_len(nrow(class_counts))) {
  cat("    ", class_counts$class[i], ":", class_counts$N[i], "genes\n")
}

cat("\n  n_ancestry_sig distribution (all 4 ancestries):\n")
for (k in 0:4) {
  cat("    ", k, "ancestries:", sum(cross$n_ancestry_sig == k), "genes\n")
}
cat("\n  n_validated_ancestry distribution (EUR+EAS only):\n")
for (k in 0:2) {
  cat("    ", k, "validated:", sum(cross$n_validated_ancestry == k), "genes\n")
}

# Pairwise correlations among shared tested genes
ancestry_cols <- c(EUR = "best_pp4_eur", EAS = "best_pp4_eas",
                   AFR = "best_pp4_afr", CSA = "best_pp4_csa")
cat("\n  Pairwise Spearman rho (among shared genes):\n")
for (i in seq_along(ancestry_cols)) {
  for (j in seq_along(ancestry_cols)) {
    if (j <= i) next
    col_a <- ancestry_cols[i]; col_b <- ancestry_cols[j]
    if (!col_a %in% names(cross) || !col_b %in% names(cross)) next
    shared_pair <- cross[get(col_a) > 0 & get(col_b) > 0]
    if (nrow(shared_pair) > 5) {
      rho <- cor(shared_pair[[col_a]], shared_pair[[col_b]], method = "spearman")
      cat("    ", names(ancestry_cols)[i], "vs", names(ancestry_cols)[j],
          ": rho =", round(rho, 3), "(N =", nrow(shared_pair), ")\n")
    }
  }
}

# Legacy: shared EUR vs EAS for backwards compatibility
shared <- cross[best_pp4_eur > 0 & best_pp4_eas > 0]

# ==============================================================================
# 4. Anchor gene validation
# ==============================================================================
cat("\n--- Step 4: Anchor gene validation ---\n")
cat("  Known ancestry-variable MASLD risk genes:\n")

# Expected patterns:
# PNPLA3: MAF ~42% EAS vs ~22% EUR → potentially stronger in BBJ
# TM6SF2: similar MAF → similar signal
# HSD17B13: lower in EAS → may lose signal
# MBOAT7: similar MAF
anchor <- c("PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7", "GCKR",
            "SLC39A8", "SORT1", "RORA", "THRB", "DGAT2")
for (g in anchor) {
  row <- cross[gene == g]
  if (nrow(row) > 0) {
    eur_val <- ifelse("best_pp4_eur" %in% names(row), round(row$best_pp4_eur[1], 3), NA)
    eas_val <- ifelse("best_pp4_eas" %in% names(row), round(row$best_pp4_eas[1], 3), NA)
    afr_val <- ifelse("best_pp4_afr" %in% names(row), round(row$best_pp4_afr[1], 3), NA)
    csa_val <- ifelse("best_pp4_csa" %in% names(row), round(row$best_pp4_csa[1], 3), NA)
    cat(sprintf("    %s: EUR=%.3f, EAS=%.3f, AFR=%.3f, CSA=%.3f → %s (n_ancestry=%d)\n",
                g, eur_val, eas_val, afr_val, csa_val,
                row$class[1], row$n_ancestry_sig[1]))
  } else {
    cat(sprintf("    %s: not tested in any COLOC\n", g))
  }
}

# ==============================================================================
# 5. FinnGen replication of European COLOC
# ==============================================================================
cat("\n--- Step 5: FinnGen replication of UKBB COLOC ---\n")

# Compare FinnGen NAFLD with UKBB liver enzyme COLOC
fg_nafld <- eur_results[source == "FinnGen_NAFLD"]
ukbb_alt <- eur_results[source == "UKBB_ALT"]

if (nrow(fg_nafld) > 0 && nrow(ukbb_alt) > 0) {
  fg_comp <- merge(
    fg_nafld[, .(gene, PP.H4_fg = PP.H4)],
    ukbb_alt[, .(gene, PP.H4_ukbb = PP.H4)],
    by = "gene"
  )
  if (nrow(fg_comp) > 5) {
    cat("  FinnGen NAFLD vs UKBB ALT:\n")
    cat("    Shared genes:", nrow(fg_comp), "\n")
    cat("    Pearson r:", round(cor(fg_comp$PP.H4_fg, fg_comp$PP.H4_ukbb), 3), "\n")
    cat("    Both PP.H4>0.5:", sum(fg_comp$PP.H4_fg > 0.5 & fg_comp$PP.H4_ukbb > 0.5), "\n")
    cat("    FinnGen-only:", sum(fg_comp$PP.H4_fg > 0.5 & fg_comp$PP.H4_ukbb <= 0.5), "\n")
    cat("    UKBB-only:", sum(fg_comp$PP.H4_fg <= 0.5 & fg_comp$PP.H4_ukbb > 0.5), "\n")
  }
}

# ==============================================================================
# 6. Liver enzyme trait concordance (EUR vs EAS)
# ==============================================================================
cat("\n--- Step 6: Liver enzyme COLOC concordance (EUR vs EAS) ---\n")

for (trait in c("ALT", "AST", "GGT")) {
  eur_trait <- eur_results[source == paste0("UKBB_", trait)]
  eas_trait <- eas_results[source == paste0("BBJ_", trait)]

  if (nrow(eur_trait) > 0 && nrow(eas_trait) > 0) {
    trait_comp <- merge(
      eur_trait[, .(gene, PP.H4_eur = PP.H4)],
      eas_trait[, .(gene, PP.H4_eas = PP.H4)],
      by = "gene"
    )
    if (nrow(trait_comp) > 5) {
      cat(sprintf("  %s: %d shared genes, r=%.3f, both>0.5=%d, EUR-only=%d, EAS-only=%d\n",
                  trait, nrow(trait_comp),
                  cor(trait_comp$PP.H4_eur, trait_comp$PP.H4_eas),
                  sum(trait_comp$PP.H4_eur > 0.5 & trait_comp$PP.H4_eas > 0.5),
                  sum(trait_comp$PP.H4_eur > 0.5 & trait_comp$PP.H4_eas <= 0.5),
                  sum(trait_comp$PP.H4_eur <= 0.5 & trait_comp$PP.H4_eas > 0.5)))
    }
  } else {
    cat(sprintf("  %s: insufficient data for comparison\n", trait))
  }
}

# ==============================================================================
# 6b. AFR/CSA liver enzyme concordance
# ==============================================================================
cat("\n--- Step 6b: AFR/CSA liver enzyme COLOC concordance ---\n")

for (pop in c("AFR", "CSA")) {
  pop_results <- if (pop == "AFR") afr_results else csa_results
  pop_label <- if (pop == "AFR") "African" else "Central/South Asian"

  for (trait in c("ALT", "AST", "GGT")) {
    eur_trait <- eur_results[source == paste0("UKBB_", trait)]
    pop_trait <- pop_results[source == paste0("PanUKBB_", pop, "_", trait)]

    if (nrow(eur_trait) > 0 && nrow(pop_trait) > 0) {
      tc <- merge(
        eur_trait[, .(gene, PP.H4_eur = PP.H4)],
        pop_trait[, .(gene, PP.H4_pop = PP.H4)],
        by = "gene"
      )
      if (nrow(tc) > 5) {
        cat(sprintf("  %s %s: %d shared, r=%.3f, both>0.5=%d, EUR-only=%d, %s-only=%d\n",
                    pop, trait, nrow(tc),
                    cor(tc$PP.H4_eur, tc$PP.H4_pop),
                    sum(tc$PP.H4_eur > 0.5 & tc$PP.H4_pop > 0.5),
                    sum(tc$PP.H4_eur > 0.5 & tc$PP.H4_pop <= 0.5),
                    pop, sum(tc$PP.H4_eur <= 0.5 & tc$PP.H4_pop > 0.5)))
      }
    } else {
      cat(sprintf("  %s %s: insufficient data\n", pop, trait))
    }
  }
}

# ==============================================================================
# 6c. Progression-specific analysis
# ==============================================================================
cat("\n--- Step 6c: Progression-specific gene analysis ---\n")

# Identify progression sources
progression_sources <- c("Ghouse_Cirrhosis", "Ghouse_HCC",
                          "deCODE_NAFL", "deCODE_Cirrhosis", "deCODE_HCC",
                          "FinnGen_NAFLD", "FinnGen_NASH", "FinnGen_HCC")

# Liver enzyme sources
enzyme_sources <- c("UKBB_ALT", "UKBB_AST", "UKBB_GGT")

prog_results <- eur_results[source %in% progression_sources]
enz_results  <- eur_results[source %in% enzyme_sources]

if (nrow(prog_results) > 0 && nrow(enz_results) > 0) {
  prog_genes <- prog_results[PP.H4 > 0.5, unique(gene)]
  enz_genes  <- enz_results[PP.H4 > 0.5, unique(gene)]

  both_prog_enz <- intersect(prog_genes, enz_genes)
  prog_only     <- setdiff(prog_genes, enz_genes)
  enz_only      <- setdiff(enz_genes, prog_genes)

  cat("  Liver enzyme genes (PP.H4>0.5):", length(enz_genes), "\n")
  cat("  Progression genes (PP.H4>0.5):", length(prog_genes), "\n")
  cat("  Overlap:", length(both_prog_enz), "\n")
  cat("  Progression-specific:", length(prog_only), "\n")
  if (length(prog_only) > 0) {
    cat("    ", paste(head(prog_only, 20), collapse = ", "), "\n")
  }
  cat("  Enzyme-specific:", length(enz_only), "\n")
}

# ==============================================================================
# 7. Save outputs
# ==============================================================================
cat("\n--- Step 7: Saving outputs ---\n")

# Full cross-ancestry comparison table
setorder(cross, -best_pp4_eur, -best_pp4_eas)
fwrite(cross, file.path(OUTDIR, "cross_ancestry_comparison.csv"))
cat("  Saved: cross_ancestry_comparison.csv\n")

# Multi-ancestry validated targets (EUR+EAS both >0.5)
validated <- cross[n_validated_ancestry >= 2][order(-n_validated_ancestry, -best_pp4_eur)]
if (nrow(validated) > 0) {
  fwrite(validated, file.path(OUTDIR, "cross_ancestry_validated_targets.csv"))
  cat("  Saved: cross_ancestry_validated_targets.csv (", nrow(validated), "genes, EUR+EAS validated)\n")
}

# Exploratory targets (AFR/CSA-only, no EUR/EAS signal)
exploratory <- cross[n_validated_ancestry == 0 & n_ancestry_sig >= 1]
if (nrow(exploratory) > 0) {
  setorder(exploratory, class, -best_pp4_afr, -best_pp4_csa)
  fwrite(exploratory, file.path(OUTDIR, "cross_ancestry_exploratory_targets.csv"))
  cat("  Saved: cross_ancestry_exploratory_targets.csv (", nrow(exploratory),
      "genes, AFR/CSA exploratory — prior-sensitive, underpowered)\n")
}

# Ancestry-specific loci (PP.H4 > 0.5 in exactly one validated ancestry)
ancestry_specific <- cross[n_ancestry_sig >= 1]
if (nrow(ancestry_specific) > 0) {
  ancestry_specific[, max_pp4 := pmax(best_pp4_eur, best_pp4_eas,
                                       best_pp4_afr, best_pp4_csa, na.rm = TRUE)]
  setorder(ancestry_specific, class, -max_pp4)
  fwrite(ancestry_specific, file.path(OUTDIR, "ancestry_specific_loci.csv"))
  cat("  Saved: ancestry_specific_loci.csv (", nrow(ancestry_specific), "genes)\n")
  for (cls in c("Multi_ancestry_validated", "EUR_only", "EAS_only",
                "AFR_exploratory", "CSA_exploratory")) {
    n <- sum(ancestry_specific$class == cls)
    if (n > 0) cat("    ", cls, ":", n, "\n")
  }
}

# Summary statistics
eur_eas_rho <- ifelse(nrow(shared) > 5,
                      round(cor(shared$best_pp4_eur, shared$best_pp4_eas,
                                method = "spearman"), 3), NA_real_)

summary_dt <- data.table(
  metric = c(
    "total_genes_tested",
    "eur_genes_sig", "eas_genes_sig", "afr_genes_sig", "csa_genes_sig",
    "multi_ancestry_validated_eur_eas",
    "multi_ancestry_any_2plus",
    "eur_only", "eas_only",
    "afr_exploratory", "csa_exploratory",
    "n_ancestry_sig_2", "n_ancestry_sig_3", "n_ancestry_sig_4",
    "n_validated_ancestry_2",
    "pp4_eur_eas_spearman"
  ),
  value = c(
    nrow(cross),
    sum(cross$best_pp4_eur > 0.5),
    sum(cross$best_pp4_eas > 0.5),
    sum(cross$best_pp4_afr > 0.5),
    sum(cross$best_pp4_csa > 0.5),
    sum(cross$n_validated_ancestry >= 2),
    sum(cross$n_ancestry_sig >= 2),
    sum(cross$class == "EUR_only"),
    sum(cross$class == "EAS_only"),
    sum(cross$class == "AFR_exploratory"),
    sum(cross$class == "CSA_exploratory"),
    sum(cross$n_ancestry_sig == 2),
    sum(cross$n_ancestry_sig == 3),
    sum(cross$n_ancestry_sig == 4),
    sum(cross$n_validated_ancestry == 2),
    eur_eas_rho
  )
)
fwrite(summary_dt, file.path(OUTDIR, "cross_ancestry_summary.csv"))
cat("  Saved: cross_ancestry_summary.csv\n")

# ==============================================================================
# 8. Figures
# ==============================================================================
cat("\n--- Step 8: Generating figures ---\n")

if (nrow(shared) > 5) {
  # PP.H4 scatter: EUR vs EAS
  pdf(file.path(FIGDIR, "cross_ancestry_pp4_scatter.pdf"), width = 7, height = 7)
  p1 <- ggplot(shared, aes(x = best_pp4_eur, y = best_pp4_eas)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "red", alpha = 0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", alpha = 0.5) +
    labs(title = "Cross-Ancestry COLOC Replication",
         subtitle = paste0("Best PP.H4 per gene: European (UKBB/FinnGen) vs East Asian (BBJ)\n",
                          "Spearman rho = ", round(cor(shared$best_pp4_eur, shared$best_pp4_eas,
                                                        method = "spearman"), 3),
                          "; N = ", nrow(shared), " shared genes"),
         x = "Best PP.H4 (European GWAS)", y = "Best PP.H4 (East Asian GWAS)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # Label anchor genes if present
  anchor_in_data <- shared[gene %in% anchor]
  if (nrow(anchor_in_data) > 0) {
    p1 <- p1 + ggrepel::geom_text_repel(
      data = anchor_in_data,
      aes(label = gene),
      size = 3, max.overlaps = 20
    )
  }
  print(p1)
  dev.off()
  cat("  Saved: cross_ancestry_pp4_scatter.pdf\n")
}

# Classification barplot (4-way)
if (nrow(cross) > 0) {
  class_dt <- cross[class != "Neither", .N, by = class]
  if (nrow(class_dt) > 0) {
    class_levels <- c("Multi_ancestry_validated", "EUR_only", "EAS_only",
                       "AFR_exploratory", "CSA_exploratory")
    class_dt[, class := factor(class, levels = class_levels)]
    class_dt <- class_dt[!is.na(class)]
    pdf(file.path(FIGDIR, "cross_ancestry_classification.pdf"), width = 8, height = 4)
    p2 <- ggplot(class_dt, aes(x = class, y = N, fill = class)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c("Multi_ancestry_validated" = "#2C7BB6",
                                    "EUR_only" = "#ABD9E9", "EAS_only" = "#5E3C99",
                                    "AFR_exploratory" = "#E66101",
                                    "CSA_exploratory" = "#FDB863")) +
      labs(title = "COLOC Signal by Ancestry Group",
           subtitle = "EUR/EAS = validated; AFR/CSA = exploratory (underpowered, prior-sensitive)",
           x = NULL, y = "Number of genes") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "none",
            axis.text.x = element_text(angle = 20, hjust = 1))
    print(p2)
    dev.off()
    cat("  Saved: cross_ancestry_classification.pdf\n")
  }

  # n_ancestry_sig distribution barplot
  nsig_dt <- cross[, .N, by = n_ancestry_sig]
  pdf(file.path(FIGDIR, "n_ancestry_sig_distribution.pdf"), width = 6, height = 4)
  p3 <- ggplot(nsig_dt, aes(x = factor(n_ancestry_sig), y = N)) +
    geom_col(fill = "#2C7BB6", width = 0.6) +
    labs(title = "Cross-Ancestry Replication Depth",
         subtitle = "Number of ancestry groups with PP.H4 > 0.5 per gene",
         x = "Number of ancestry groups with PP.H4 > 0.5",
         y = "Number of genes") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  print(p3)
  dev.off()
  cat("  Saved: n_ancestry_sig_distribution.pdf\n")
}

# Per-trait EUR vs EAS scatter
for (trait in c("ALT", "AST", "GGT")) {
  eur_trait <- eur_results[source == paste0("UKBB_", trait)]
  eas_trait <- eas_results[source == paste0("BBJ_", trait)]
  if (nrow(eur_trait) > 0 && nrow(eas_trait) > 0) {
    tc <- merge(eur_trait[, .(gene, PP.H4_eur = PP.H4)],
                eas_trait[, .(gene, PP.H4_eas = PP.H4)], by = "gene")
    if (nrow(tc) > 10) {
      pdf(file.path(FIGDIR, paste0("cross_ancestry_", tolower(trait), "_scatter.pdf")),
          width = 6, height = 6)
      p <- ggplot(tc, aes(x = PP.H4_eur, y = PP.H4_eas)) +
        geom_point(alpha = 0.3, size = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = paste0(trait, " COLOC: UKBB (EUR) vs BBJ (EAS)"),
             subtitle = paste0("r = ", round(cor(tc$PP.H4_eur, tc$PP.H4_eas), 3),
                              "; N = ", nrow(tc)),
             x = paste0("PP.H4 UKBB ", trait, " (EUR)"),
             y = paste0("PP.H4 BBJ ", trait, " (EAS)")) +
        theme_bw(base_size = 11)
      print(p)
      dev.off()
      cat(paste0("  Saved: cross_ancestry_", tolower(trait), "_scatter.pdf\n"))
    }
  }
}

cat("\n=== Script 48: Complete ===\n")
cat("Ancestry groups: EUR, EAS (BBJ), AFR (Pan-UKBB), CSA (Pan-UKBB)\n")
cat("Progression: Ghouse (Cirrhosis, HCC), deCODE (NAFL, Cirrhosis, HCC)\n")
cat("End time:", format(Sys.time()), "\n")
cat("\nKey outputs:\n")
cat("  ", file.path(OUTDIR, "cross_ancestry_comparison.csv"), "\n")
cat("  ", file.path(OUTDIR, "cross_ancestry_validated_targets.csv"), "\n")
cat("  ", file.path(OUTDIR, "ancestry_specific_loci.csv"), "\n")
cat("  ", file.path(OUTDIR, "cross_ancestry_summary.csv"), "\n")
