#!/usr/bin/env Rscript
# 20_drug_repurposing_cmap.R
#
# Strategy 11 v2: Drug Repurposing via CMap-style Analysis
#
# Key changes from v1:
#   - Uses t-statistic (not logFC) as ranking metric for better sensitivity
#   - Primary: signatureSearch LINCS L1000 query (HepG2 cell line) if installed
#   - Secondary: fgsea with MSigDB C2:CGP using t-stat ranking + increased nPermSimple
#   - C2:CGP results reframed as "disease signature concordance" (supplementary)
#   - (MR convergence check for drug targets — REMOVED 2026-04-22; MR ditched from paper)
#
# Inputs:
#   - Dream mega-analysis results (dream_results.csv)
#   - MSigDB C2:CGP gene sets via msigdbr
#   - Optionally: signatureSearch + signatureSearchData (Bioconductor)
#
# Outputs (results/drug_repurposing/):
#   - lincs_reversal_hits.csv            — LINCS L1000 reversal hits (if available)
#   - lincs_top50_reversals.csv          — top 50 LINCS reversals
#   - cgp_disease_concordance.csv        — C2:CGP concordance (supplementary)
#   - cgp_reversal_hits.csv              — C2:CGP reversal hits (fgsea, t-stat)
#   - top50_reversal_compounds.csv       — top 50 fgsea reversal compounds
#   - (multi_layer_drug_targets.csv — retired 2026-04-22 with MR ditch)

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(biomaRt)
  library(patchwork)
})

# ============================================================
#  Theme & palette (Sanjana lab standard)
# ============================================================
pdf.options(useDingbats = FALSE)

theme_pub <- theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title       = element_text(size = 8, face = "bold"),
    plot.subtitle    = element_text(size = 7, color = "grey40"),
    axis.title       = element_text(size = 8),
    axis.text        = element_text(size = 6),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 7),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 7, face = "bold")
  )

pal <- list(
  magenta = "#e14b9d", pink = "#e35070", purple = "#d358c7",
  blue = "#4baeef", orange = "#e1b172", green = "#30d796",
  teal = "#2bbfbd", grey = "#808080"
)

cat("=== Strategy 11 v2: Drug Repurposing via CMap-style Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# Check for signatureSearch availability
HAS_SIGSEARCH <- requireNamespace("signatureSearch", quietly = TRUE) &&
                 requireNamespace("signatureSearchData", quietly = TRUE)
if (HAS_SIGSEARCH) {
  cat("  signatureSearch: AVAILABLE (will run LINCS L1000 query)\n")
} else {
  cat("  signatureSearch: NOT INSTALLED (skipping LINCS L1000; using fgsea only)\n")
}

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_DIR  <- file.path(BASE_DIR, "RNA-seq")
RESULTS_DIR <- file.path(RNASEQ_DIR, "results/drug_repurposing")
FIG_DIR     <- file.path(BASE_DIR, "figures")

DREAM_FILE  <- file.path(RNASEQ_DIR,
  "Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")
# MR_FILE removed 2026-04-22 — MR ditched from paper.

FGSEA_PADJ_THRESHOLD <- 0.05
FGSEA_MIN_SIZE       <- 15
FGSEA_MAX_SIZE       <- 500

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Part A: Build Disease Signature from Dream Results (t-statistic ranking)
# ==============================================================================
cat("--- Part A: Building disease signature from dream results ---\n")

dream <- fread(DREAM_FILE)
cat("  Loaded dream results:", nrow(dream), "genes\n")

# Strip Ensembl version suffixes
dream[, ensembl_id := sub("\\.\\d+$", "", gene)]

# Remove duplicates: keep the one with smallest P.Value (most significant)
setorder(dream, P.Value)
dream_dedup <- dream[!duplicated(ensembl_id)]
cat("  After deduplication:", nrow(dream_dedup), "unique Ensembl IDs\n")

# Use t-statistic as ranking metric (better sensitivity than logFC alone)
ranked_stats <- dream_dedup$t
names(ranked_stats) <- dream_dedup$ensembl_id
ranked_stats <- ranked_stats[!is.na(ranked_stats)]
ranked_stats <- sort(ranked_stats, decreasing = TRUE)

cat("  Ranked gene list:", length(ranked_stats), "genes (ranked by t-statistic)\n")
cat("  t-stat range:", round(min(ranked_stats), 3), "to",
    round(max(ranked_stats), 3), "\n\n")

# Build gene symbol mapping: local annotation cache (faster, no network) with biomaRt fallback
ann_cache_path <- file.path(BASE_DIR, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
cat("  Building Ensembl -> HGNC symbol mapping...\n")
if (file.exists(ann_cache_path)) {
  cat("  Using local annotation cache...\n")
  ann_cache <- fread(ann_cache_path)
  ensembl_to_symbol <- ann_cache[symbol != "" & !is.na(symbol),
                                  .(ensembl_id = gene_base, symbol)]
  ensembl_to_symbol <- ensembl_to_symbol[!duplicated(ensembl_id)]
  cat("  Mapped", nrow(ensembl_to_symbol), "Ensembl IDs to HGNC symbols\n")
} else {
  cat("  Querying biomaRt for Ensembl -> HGNC symbol mapping...\n")
  ensembl_to_symbol <- tryCatch({
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                       mirror = "useast")
    bm <- as.data.table(getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters    = "ensembl_gene_id",
      values     = names(ranked_stats),
      mart       = mart
    ))
    setnames(bm, c("ensembl_id", "symbol"))
    bm <- bm[symbol != ""]
    bm <- bm[!duplicated(ensembl_id)]
    cat("  Mapped", nrow(bm), "Ensembl IDs to HGNC symbols\n")
    bm
  }, error = function(e) {
    cat("  WARNING: biomaRt failed:", conditionMessage(e), "\n")
    data.table(ensembl_id = character(), symbol = character())
  })
}

# ==============================================================================
# Part B: Primary Analysis — signatureSearch LINCS L1000 (if available)
# ==============================================================================
lincs_results <- NULL

if (HAS_SIGSEARCH) {
  cat("\n--- Part B: LINCS L1000 Query via signatureSearch ---\n")
  suppressPackageStartupMessages({
    library(signatureSearch)
    library(signatureSearchData)
    library(ExperimentHub)
  })

  tryCatch({
    # Load LINCS L1000 database
    cat("  Loading LINCS L1000 reference database...\n")
    eh <- ExperimentHub()

    # Get the LINCS L1000 Level 5 data (gene expression signatures)
    lincs_db_path <- eh[["EH3226"]]  # LINCS L1000 compound signatures

    # Build query from top up/down-regulated genes
    # signatureSearch expects a list with "upset" and "downset" gene symbols
    n_query <- 150  # top 150 up and down

    # Map Ensembl IDs to gene symbols for the query
    top_up_ensembl   <- names(head(ranked_stats, n_query))
    top_down_ensembl <- names(tail(ranked_stats, n_query))

    top_up_symbols   <- ensembl_to_symbol[ensembl_id %in% top_up_ensembl, symbol]
    top_down_symbols <- ensembl_to_symbol[ensembl_id %in% top_down_ensembl, symbol]

    cat(sprintf("  Query: %d up-regulated genes, %d down-regulated genes\n",
                length(top_up_symbols), length(top_down_symbols)))

    if (length(top_up_symbols) >= 10 && length(top_down_symbols) >= 10) {
      # Run LINCS query
      qsig <- qSig(query = list(upset = top_up_symbols,
                                  downset = top_down_symbols),
                    gess_method = "LINCS",
                    refdb = lincs_db_path)

      cat("  Running GESS LINCS query (this may take several minutes)...\n")
      lincs_res <- gess_lincs(qSig = qsig, sortby = "NCS", tau = TRUE,
                               workers = 4)
      lincs_dt <- as.data.table(result(lincs_res))

      cat("  LINCS query complete:", nrow(lincs_dt), "perturbations tested\n")

      # Filter for HepG2 cell line (liver-relevant)
      if ("cell" %in% names(lincs_dt)) {
        hepg2_hits <- lincs_dt[grepl("HepG2", cell, ignore.case = TRUE)]
        cat("  HepG2 cell line hits:", nrow(hepg2_hits), "\n")
      } else {
        hepg2_hits <- lincs_dt
        cat("  Cell line column not found; using all results\n")
      }

      # Reversal hits: negative connectivity score (NCS < 0 = reversing disease signature)
      # tau < -90 is considered a strong reversal in CMap convention
      lincs_reversals <- hepg2_hits[NCS < 0]
      lincs_reversals <- lincs_reversals[order(NCS)]
      cat("  LINCS reversal hits (NCS < 0):", nrow(lincs_reversals), "\n")

      strong_reversals <- hepg2_hits[tau < -90]
      cat("  Strong reversals (tau < -90):", nrow(strong_reversals), "\n")

      # Save LINCS results
      fwrite(lincs_reversals, file.path(RESULTS_DIR, "lincs_reversal_hits.csv"))
      cat("  Saved: lincs_reversal_hits.csv\n")

      lincs_top50 <- head(lincs_reversals, 50)
      fwrite(lincs_top50, file.path(RESULTS_DIR, "lincs_top50_reversals.csv"))
      cat("  Saved: lincs_top50_reversals.csv\n")

      lincs_results <- lincs_reversals
    } else {
      cat("  WARNING: Not enough mapped gene symbols for LINCS query.\n")
    }
  }, error = function(e) {
    cat("  ERROR in signatureSearch:", conditionMessage(e), "\n")
    cat("  Falling back to fgsea-only analysis.\n")
  })
} else {
  cat("\n--- Part B: LINCS L1000 SKIPPED (signatureSearch not installed) ---\n")
  cat("  Install with: BiocManager::install(c('signatureSearch', 'signatureSearchData'))\n\n")
}

# ==============================================================================
# Part C: Secondary Analysis — fgsea with C2:CGP (t-stat ranked)
# ==============================================================================
cat("\n--- Part C: fgsea with MSigDB C2:CGP (t-statistic ranking) ---\n")

# Get C2:CGP (Chemical and Genetic Perturbations) gene sets
cgp_df <- as.data.table(msigdbr(species = "Homo sapiens",
                                 collection = "C2",
                                 subcollection = "CGP"))
cat("  C2:CGP gene sets loaded:", uniqueN(cgp_df$gs_name), "sets\n")

# Build pathway list using Ensembl IDs
cgp_df <- cgp_df[ensembl_gene != ""]
pathways <- split(cgp_df$ensembl_gene, cgp_df$gs_name)

# Filter to pathways that have genes in our ranked list
pathway_overlap <- vapply(pathways, function(genes) {
  sum(genes %in% names(ranked_stats))
}, integer(1))
cat("  Pathways with >=", FGSEA_MIN_SIZE, "overlapping genes:",
    sum(pathway_overlap >= FGSEA_MIN_SIZE), "\n")

# Run fgsea with increased nPermSimple to avoid NA p-values
cat("  Running fgsea (nPermSimple=10000 for stability)...\n")
set.seed(42)
fgsea_res <- fgsea(pathways     = pathways,
                   stats        = ranked_stats,
                   minSize      = FGSEA_MIN_SIZE,
                   maxSize      = FGSEA_MAX_SIZE,
                   nPermSimple  = 10000)

fgsea_dt <- as.data.table(fgsea_res)
n_tested <- nrow(fgsea_dt)
n_valid  <- sum(!is.na(fgsea_dt$padj))
n_na     <- sum(is.na(fgsea_dt$padj))
cat(sprintf("  fgsea completed: %d gene sets tested (%d valid, %d NA)\n",
            n_tested, n_valid, n_na))
cat("  Significant (padj < 0.05):",
    sum(fgsea_dt$padj < FGSEA_PADJ_THRESHOLD, na.rm = TRUE), "\n")

# ==============================================================================
# Part D: Filter Reversal and Concordant Hits
# ==============================================================================
cat("\n--- Part D: Filtering reversal and concordant hits ---\n")

# Reversal hits: NES < 0 = signature opposes disease transcriptome
reversal_hits <- fgsea_dt[!is.na(NES) & NES < 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD]
setorder(reversal_hits, NES)
cat("  Reversal hits (NES < 0, padj < 0.05):", nrow(reversal_hits), "\n")

# Disease-concordant hits: NES > 0 = signature parallels disease
concordant_hits <- fgsea_dt[!is.na(NES) & NES > 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD]
cat("  Concordant hits (NES > 0, padj < 0.05):", nrow(concordant_hits), "\n")

# Classify gene sets by name patterns for interpretability
classify_cgp <- function(pathway_name) {
  pn <- toupper(pathway_name)
  if (grepl("LIVER|HEPAT|HCC|NAFLD|NASH|STEATO|CIRR", pn)) return("Liver/MASLD")
  if (grepl("DRUG|TREAT|COMPOUND|DOXO|ETOPO|TAMOX|METFORM", pn)) return("Drug_Treatment")
  if (grepl("OBESE|OBESITY|BMI|ADIPOSE|FAT|LIPID|CHOLEST", pn)) return("Metabolic")
  if (grepl("INFLAM|TNF|IL6|NFKB|IMMUNE|MACROPHAGE", pn)) return("Inflammatory")
  if (grepl("FIBROSIS|STELLATE|COLLAGEN|TGF", pn)) return("Fibrosis")
  return("Other_CGP")
}

fgsea_dt[, cgp_class := vapply(pathway, classify_cgp, character(1))]
reversal_hits[, cgp_class := vapply(pathway, classify_cgp, character(1))]

cat("\n  Reversal hits by class:\n")
if (nrow(reversal_hits) > 0) {
  print(reversal_hits[, .N, by = cgp_class][order(-N)])
}

# ==============================================================================
# Part E: Extract top 50 and save results
# ==============================================================================
cat("\n--- Part E: Saving results ---\n")

top50 <- head(reversal_hits, 50)

# Clean up gene set names for display
clean_name <- function(x) {
  x <- gsub("^[A-Z]+_", "", x)
  x <- gsub("_DN$|_UP$", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- gsub("(^|\\s)(\\w)", "\\1\\U\\2", x, perl = TRUE)
  ifelse(nchar(x) > 60, paste0(substr(x, 1, 57), "..."), x)
}
top50[, display_name := clean_name(pathway)]

# Helper to save fgsea results
save_fgsea <- function(dt, file) {
  out <- copy(dt)
  if ("leadingEdge" %in% names(out)) {
    out[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  }
  fwrite(out, file)
  cat("  Saved:", file, "(", nrow(out), "rows)\n")
}

# Disease concordance (C2:CGP as supplementary)
save_fgsea(fgsea_dt[order(padj)],
           file.path(RESULTS_DIR, "cgp_disease_concordance.csv"))

# Reversal hits
save_fgsea(reversal_hits,
           file.path(RESULTS_DIR, "cgp_reversal_hits.csv"))

# Top 50
save_fgsea(top50,
           file.path(RESULTS_DIR, "top50_reversal_compounds.csv"))

# ==============================================================================
# Part F: MR Convergence Analysis — REMOVED 2026-04-22
# ==============================================================================
# MR permanently ditched from the paper; TWAS + COLOC + INTACT is the causal
# framework. The prior MR convergence block read causal_inference_summary.csv
# and wrote multi_layer_drug_targets.csv. TWAS-only drug-target convergence
# happens downstream via Script 42 (regulon_drug_gwas_convergence) and Script
# 45a (integrate_all_sources).
cat("\n--- Part F: MR convergence analysis — SKIPPED (MR ditched 2026-04-22) ---\n")

# ==============================================================================
# Part G: Figures
# ==============================================================================
cat("\n--- Part G: Generating figures ---\n")

# ---- Figure 1: NES distribution histogram ----
cat("  Generating NES histogram...\n")

plot_dt <- fgsea_dt[!is.na(NES)]
p1 <- ggplot(plot_dt, aes(x = NES)) +
  geom_histogram(bins = 60, fill = "grey80", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_histogram(data = plot_dt[NES < 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD],
                 aes(x = NES), bins = 60,
                 fill = pal$teal, color = pal$teal, alpha = 0.7, linewidth = 0.3) +
  geom_histogram(data = plot_dt[NES > 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD],
                 aes(x = NES), bins = 60,
                 fill = pal$magenta, color = pal$magenta, alpha = 0.7, linewidth = 0.3) +
  annotate("text", x = min(plot_dt$NES, na.rm = TRUE) * 0.7, y = Inf,
           label = paste0("Reversal: ", nrow(reversal_hits)),
           hjust = 0, vjust = 1.5, color = pal$teal, fontface = "bold", size = 2) +
  annotate("text", x = max(plot_dt$NES, na.rm = TRUE) * 0.7, y = Inf,
           label = paste0("Concordant: ", nrow(concordant_hits)),
           hjust = 1, vjust = 1.5, color = pal$magenta, fontface = "bold", size = 2) +
  labs(title = "Drug Repurposing: NES Distribution (t-statistic ranked)",
       subtitle = paste0("MSigDB C2:CGP (", n_valid,
                         " gene sets tested); ranked by dream t-statistic"),
       x = "Normalized Enrichment Score (NES)",
       y = "Count") +
  theme_pub

ggsave(file.path(FIG_DIR, "drug_1_nes_histogram.pdf"),
       p1, width = 4.5, height = 3, device = cairo_pdf)
cat("  Saved: figures/drug_1_nes_histogram.pdf\n")

# ---- Figure 2: Top 20 reversal compounds bar chart ----
cat("  Generating top 20 reversal bar chart...\n")

if (nrow(reversal_hits) > 0) {
  top20 <- head(reversal_hits, min(20, nrow(reversal_hits)))
  top20[, display_name := clean_name(pathway)]
  top20[, display_name := factor(display_name, levels = rev(unique(display_name)))]

  p2 <- ggplot(top20, aes(x = NES, y = display_name, fill = -log10(padj))) +
    geom_col(width = 0.7) +
    scale_fill_gradient(low = pal$blue, high = pal$teal,
                        name = expression(-log[10](p[adj]))) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    labs(title = "Top 20 Reversal Perturbation Signatures",
         subtitle = "Gene sets anticorrelated with MASLD transcriptome",
         x = "Normalized Enrichment Score (NES)",
         y = NULL) +
    theme_pub +
    theme(panel.grid.major.y = element_blank())

  class_summary <- reversal_hits[, .(n_sets = .N,
                                      mean_NES = mean(NES, na.rm = TRUE)),
                                  by = cgp_class][order(-n_sets)]

  p3 <- ggplot(class_summary, aes(x = reorder(cgp_class, n_sets),
                                   y = n_sets, fill = mean_NES)) +
    geom_col(width = 0.7) +
    scale_fill_gradient2(low = pal$teal, mid = "grey90", high = pal$magenta,
                         midpoint = 0, name = "Mean NES") +
    coord_flip() +
    labs(title = "Reversal Hits by Gene Set Class",
         subtitle = "Classification of C2:CGP sets",
         x = NULL, y = "Number of reversal gene sets") +
    theme_pub
    
  p_combined_pharma <- p2 + p3 + plot_layout(ncol = 2, widths = c(1.5, 1))
  ggsave(file.path(FIG_DIR, "drug_reversal_combined_v3.pdf"),
         p_combined_pharma, width = 7.5, height = 4.5, device = cairo_pdf)
  cat("  Saved: figures/drug_reversal_combined_v3.pdf\n")
} else {
  cat("  No reversal hits to plot.\n")
}

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Strategy 11 v2: Drug Repurposing Complete ===\n")
cat("  Ranking metric:         t-statistic (dream mega-analysis)\n")
cat("  LINCS L1000:           ", ifelse(HAS_SIGSEARCH, "RAN", "SKIPPED (not installed)"), "\n")
if (!is.null(lincs_results)) {
  cat("  LINCS reversal hits:   ", nrow(lincs_results), "\n")
}
cat("  C2:CGP sets tested:    ", n_valid, "\n")
cat("  fgsea reversal hits:   ", nrow(reversal_hits), "\n")
cat("  fgsea concordant hits: ", nrow(concordant_hits), "\n")
cat("  Output dir:            ", RESULTS_DIR, "\n")
cat("  Figures dir:           ", FIG_DIR, "\n")
cat("End time:", format(Sys.time()), "\n")
