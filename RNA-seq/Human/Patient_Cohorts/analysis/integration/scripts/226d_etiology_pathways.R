#!/usr/bin/env Rscript
# 226d_etiology_pathways.R
# ---------------------------------------------------------------------------
# Etiology-specific protein programs from per-class SHAP values + fgsea.
#
# Loads the T4 etiology 3-class SHAP .npz from Script 226a
# (shap_values: 177 subjects x 1461 proteins x 3 classes) via reticulate,
# computes per-protein signed mean SHAP per class, builds ranked gene lists,
# and runs fgsea pathway enrichment (Hallmark, KEGG MEDICUS, Reactome).
#
# Inputs:
#   - 226a_shap_T4_etiology_3.npz    (per-class SHAP from 226a)
#   - olink.qc.finished.mendeley.data.txt (protein gene symbols)
#
# Outputs (to results/multiprogram/):
#   226d_etiology_pathways.csv     — pathway, NES, padj, pval, size,
#                                     etiology_class, collection
#   226d_etiology_protein_lists.csv — protein, masld_shap, cvh_shap,
#                                      arld_shap, top_class
#   226d_pathway_matrix.csv        — wide NES: pathways x etiology classes
#
# Usage: Rscript 226d_etiology_pathways.R
# SLURM: --partition=cpu --cpus-per-task=4 --mem=32G --time=02:00:00
# Env:   micromamba activate rnaseq
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(reticulate)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 226d: Etiology-Specific Pathway Enrichment ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# 1. Load per-class SHAP values via reticulate
# ============================================================
cat("Loading per-class SHAP values...\n")

npz_path <- file.path(OUTDIR, "226a_shap_T4_etiology_3.npz")
if (!file.exists(npz_path)) {
  stop(sprintf(
    "SHAP file not found: %s\n  Run 226a_plasma_shap_panels.py first.",
    npz_path
  ))
}

np <- import("numpy")
data <- np$load(npz_path, allow_pickle = TRUE)

# shap_values: (177, 1461, 3) — subjects x proteins x classes
shap_vals <- data["shap_values"]
protein_names_np <- data["protein_names"]
protein_names <- as.character(protein_names_np)

n_subjects <- dim(shap_vals)[1]
n_proteins <- dim(shap_vals)[2]
n_classes  <- dim(shap_vals)[3]
cat(sprintf("  SHAP array: %d subjects x %d proteins x %d classes\n",
            n_subjects, n_proteins, n_classes))
cat(sprintf("  Protein names: %d (first: %s, last: %s)\n",
            length(protein_names), protein_names[1],
            protein_names[length(protein_names)]))

# Validate dimensions
stopifnot(n_classes == 3)
stopifnot(length(protein_names) == n_proteins)

# Class labels matching 226a encoding: MASLD=0, CVH=1, ARLD=2
class_labels <- c("MASLD", "CVH", "ARLD")

# ============================================================
# 2. Compute per-protein signed mean SHAP per class
# ============================================================
cat("Computing per-protein SHAP statistics...\n")

# Signed mean SHAP: colMeans(shap_vals[,,k]) for each class k
# This preserves the direction — positive SHAP means the protein pushes
# the prediction towards that class
signed_shap <- matrix(NA_real_, nrow = n_proteins, ncol = n_classes,
                      dimnames = list(protein_names, class_labels))

for (k in seq_len(n_classes)) {
  # R uses 1-based indexing; class k in R is k-1 in Python
  class_slice <- shap_vals[, , k]  # (n_subjects, n_proteins)
  signed_shap[, k] <- colMeans(class_slice, na.rm = TRUE)
}

# Also compute mean |SHAP| per protein per class (for ranking)
abs_shap <- matrix(NA_real_, nrow = n_proteins, ncol = n_classes,
                   dimnames = list(protein_names, class_labels))
for (k in seq_len(n_classes)) {
  class_slice <- shap_vals[, , k]
  abs_shap[, k] <- colMeans(abs(class_slice), na.rm = TRUE)
}

# Determine top_class: which class each protein most distinguishes
# (highest mean |SHAP|)
top_class_idx <- apply(abs_shap, 1, which.max)
top_class <- class_labels[top_class_idx]

# Build protein-level output
protein_dt <- data.table(
  protein    = protein_names,
  masld_shap = signed_shap[, "MASLD"],
  cvh_shap   = signed_shap[, "CVH"],
  arld_shap  = signed_shap[, "ARLD"],
  masld_abs_shap = abs_shap[, "MASLD"],
  cvh_abs_shap   = abs_shap[, "CVH"],
  arld_abs_shap  = abs_shap[, "ARLD"],
  top_class  = top_class
)

# Sort by max absolute SHAP across classes
protein_dt[, max_abs_shap := pmax(masld_abs_shap, cvh_abs_shap, arld_abs_shap)]
setorder(protein_dt, -max_abs_shap)

cat(sprintf("  Top 10 proteins by max |SHAP|:\n"))
for (i in 1:min(10, nrow(protein_dt))) {
  row <- protein_dt[i]
  cat(sprintf("    %2d. %-12s  top_class=%-5s  |SHAP|=%.4f  "
              , i, row$protein, row$top_class, row$max_abs_shap))
  cat(sprintf("(M=%.4f, C=%.4f, A=%.4f)\n",
              row$masld_shap, row$cvh_shap, row$arld_shap))
}

# ============================================================
# 3. Create ranked gene lists for fgsea
# ============================================================
cat("\nBuilding ranked gene lists for fgsea...\n")

# Use signed mean SHAP as the ranking statistic.
# Positive SHAP = protein pushes prediction towards that class.
ranked_lists <- list()
for (k in seq_len(n_classes)) {
  label <- class_labels[k]
  stats <- signed_shap[, k]
  names(stats) <- protein_names

  # Remove NaN/Inf
  stats <- stats[is.finite(stats)]

  # Remove duplicates (should not happen with Olink, but safety)
  stats <- stats[!duplicated(names(stats))]

  # Sort descending
  stats <- sort(stats, decreasing = TRUE)

  ranked_lists[[label]] <- stats
  cat(sprintf("  %s: %d proteins ranked (range %.4f to %.4f)\n",
              label, length(stats), max(stats), min(stats)))
}

# ============================================================
# 4. Load MSigDB pathway gene sets
# ============================================================
cat("\nLoading MSigDB gene sets (gene_symbol for proteomics)...\n")

# Hallmark
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
hallmark <- lapply(hallmark, function(x) unique(x[!is.na(x) & x != ""]))
cat(sprintf("  Hallmark: %d pathways\n", length(hallmark)))

# KEGG MEDICUS (newer curated KEGG)
kegg_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CP:KEGG_MEDICUS")
kegg <- split(kegg_df$gene_symbol, kegg_df$gs_name)
kegg <- lapply(kegg, function(x) unique(x[!is.na(x) & x != ""]))
cat(sprintf("  KEGG MEDICUS: %d pathways\n", length(kegg)))

# Reactome
reactome_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                       subcollection = "CP:REACTOME")
reactome <- split(reactome_df$gene_symbol, reactome_df$gs_name)
reactome <- lapply(reactome, function(x) unique(x[!is.na(x) & x != ""]))
cat(sprintf("  Reactome: %d pathways\n", length(reactome)))

all_sets <- c(hallmark, kegg, reactome)
cat(sprintf("  Total: %d pathways\n", length(all_sets)))

# Track collection membership for each pathway
pathway_collection <- c(
  setNames(rep("Hallmark", length(hallmark)), names(hallmark)),
  setNames(rep("KEGG_MEDICUS", length(kegg)), names(kegg)),
  setNames(rep("Reactome", length(reactome)), names(reactome))
)

# ============================================================
# 5. Run fgsea per etiology class
# ============================================================
cat("\nRunning fgsea...\n")

all_results <- list()

for (label in class_labels) {
  cat(sprintf("  %s:\n", label))

  stats <- ranked_lists[[label]]

  # Check overlap with gene sets
  all_pathway_genes <- unique(unlist(all_sets))
  n_overlap <- sum(names(stats) %in% all_pathway_genes)
  cat(sprintf("    %d / %d proteins overlap with pathway gene sets (%.1f%%)\n",
              n_overlap, length(stats), 100 * n_overlap / length(stats)))

  res <- fgsea(
    pathways    = all_sets,
    stats       = stats,
    minSize     = 5,
    maxSize     = 500,
    nPermSimple = 10000
  )

  res[, etiology_class := label]
  res[, collection := pathway_collection[pathway]]

  # Convert leadingEdge list to semicolon-separated string
  res[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]

  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  cat(sprintf("    %d pathways tested, %d significant (padj < 0.05)\n",
              nrow(res), n_sig))

  all_results[[label]] <- res
}

gsea_dt <- rbindlist(all_results)
cat(sprintf("\nTotal fgsea results: %d rows\n", nrow(gsea_dt)))

# ============================================================
# 6. Save outputs
# ============================================================

# Output 1: Full pathway results
out1 <- file.path(OUTDIR, "226d_etiology_pathways.csv")
# Use write.csv to avoid data.table fwrite conflict with fgsea's nMoreExtreme column
out_df <- as.data.frame(gsea_dt[, .(pathway, NES, padj, pval, size,
                                      etiology_class, collection, ES)])
out_df$leadingEdge <- sapply(gsea_dt$leadingEdge, paste, collapse = ";")
write.csv(out_df, out1, row.names = FALSE)
cat(sprintf("Saved: %s (%d rows)\n", out1, nrow(gsea_dt)))

# Output 2: Protein-level SHAP lists
out2 <- file.path(OUTDIR, "226d_etiology_protein_lists.csv")
fwrite(protein_dt[, .(protein, masld_shap, cvh_shap, arld_shap,
                       masld_abs_shap, cvh_abs_shap, arld_abs_shap,
                       top_class, max_abs_shap)],
       out2)
cat(sprintf("Saved: %s (%d proteins)\n", out2, nrow(protein_dt)))

# Output 3: Wide NES matrix (pathways x etiology classes)
# Pivot: rows = pathway, columns = etiology_class, values = NES
nes_wide <- dcast(gsea_dt, pathway + collection ~ etiology_class,
                  value.var = "NES", fill = NA_real_)
setorder(nes_wide, collection, pathway)
out3 <- file.path(OUTDIR, "226d_pathway_matrix.csv")
fwrite(nes_wide, out3)
cat(sprintf("Saved: %s (%d pathways)\n", out3, nrow(nes_wide)))

# ============================================================
# 7. Summary statistics
# ============================================================
cat("\n=== Summary ===\n")

# Per-class significance counts
for (label in class_labels) {
  sub <- gsea_dt[etiology_class == label]
  n_up   <- sum(sub$padj < 0.05 & sub$NES > 0, na.rm = TRUE)
  n_down <- sum(sub$padj < 0.05 & sub$NES < 0, na.rm = TRUE)
  cat(sprintf("  %s: %d significant (padj<0.05): %d up, %d down\n",
              label, n_up + n_down, n_up, n_down))
}

# Per-collection significance counts
cat("\n  Per collection:\n")
for (coll in c("Hallmark", "KEGG_MEDICUS", "Reactome")) {
  sub <- gsea_dt[collection == coll & padj < 0.05]
  if (nrow(sub) > 0) {
    counts <- sub[, .N, by = etiology_class]
    cat(sprintf("    %s: %s\n", coll,
                paste(sprintf("%s=%d", counts$etiology_class, counts$N),
                      collapse = ", ")))
  } else {
    cat(sprintf("    %s: 0 significant\n", coll))
  }
}

# Top 5 pathways per class (by |NES|)
cat("\n  Top 5 pathways per class (by |NES|, padj<0.25):\n")
for (label in class_labels) {
  sub <- gsea_dt[etiology_class == label & padj < 0.25]
  sub <- sub[order(-abs(NES))][1:min(5, nrow(sub))]
  cat(sprintf("\n  %s:\n", label))
  if (nrow(sub) == 0) {
    cat("    (none at padj<0.25)\n")
  } else {
    for (i in seq_len(nrow(sub))) {
      row <- sub[i]
      cat(sprintf("    %d. %s  NES=%.2f  padj=%.3g  [%s]\n",
                  i, row$pathway, row$NES, row$padj, row$collection))
    }
  }
}

# Class-specific pathway divergence: pathways significant in one class but
# not others (|NES| > 1.5 in one, |NES| < 0.5 in others)
cat("\n  Class-specific pathways (|NES|>1.5 in one, <0.5 in others, padj<0.1):\n")
sig_pathways <- gsea_dt[padj < 0.1 & abs(NES) > 1.5, unique(pathway)]
for (pw in sig_pathways) {
  pw_sub <- nes_wide[pathway == pw]
  if (nrow(pw_sub) == 0) next
  vals <- unlist(pw_sub[, class_labels, with = FALSE])
  if (any(is.na(vals))) next
  # Check if exactly one class has |NES| > 1.5 and others < 0.5
  high <- abs(vals) > 1.5
  low  <- abs(vals) < 0.5
  if (sum(high) == 1 && sum(low) >= 1) {
    which_high <- class_labels[which(high)]
    cat(sprintf("    %s → %s-specific (NES: M=%.2f, C=%.2f, A=%.2f)\n",
                pw, which_high, vals[1], vals[2], vals[3]))
  }
}

cat(sprintf("\nCompleted: %s\n", as.character(Sys.time())))
cat("=== 226d done ===\n")
