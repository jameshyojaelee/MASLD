#!/usr/bin/env Rscript
# 08b_unified_pathway_gsea.R
# ---------------------------------------------------------------------------
# Unified GSEA across ALL contrasts and ALL gene set collections.
# Supersedes Scripts 08, 14c, 14f, 133, 115b.
#
# Contrasts:
#   C1  Disease vs Control          (dream_results.csv)
#   C2  NASH vs NAFL                (nafl_vs_nash_dream.csv)
#   Fibrosis cumulative  F1-F4 vs F0 (fibrosis_stage_dream.csv)
#   Fibrosis consecutive F0->F4     (fibrosis_consecutive_dream.csv)
#   Progression binary              (c3-c13 dream CSVs)
#   Ordinal                         (c17_fibrosis_ordinal_dream.csv)
#
# Collections: H, C2:CP:KEGG_MEDICUS, C2:CP:REACTOME,
#              C5:GO:BP, C5:GO:MF, C5:GO:CC
#
# GO semantic clustering via rrvgo (Wang similarity, threshold 0.7)
# to address reviewer R1 section 6.A (pathway non-independence).
#
# Output: results/integration/gsea_kallisto/
#
# Usage: Rscript 08b_unified_pathway_gsea.R
# SLURM: cpu, 16 CPU, 64G, 48h
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

cat("=== 08b: Unified Pathway GSEA ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("fgsea version:", as.character(packageVersion("fgsea")), "\n")
cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n\n")

# Check rrvgo availability upfront
has_rrvgo <- requireNamespace("rrvgo", quietly = TRUE)
if (has_rrvgo) {
  library(rrvgo)
  cat("rrvgo version:", as.character(packageVersion("rrvgo")), "-- GO clustering enabled\n")
} else {
  warning("rrvgo not installed. GO semantic clustering will be SKIPPED. ",
          "Install with: BiocManager::install('rrvgo')")
}

# Also need org.Hs.eg.db for rrvgo GO-to-term mapping
has_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (has_rrvgo && !has_orgdb) {
  warning("org.Hs.eg.db not installed. rrvgo clustering requires it. Clustering SKIPPED.")
  has_rrvgo <- FALSE
}

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Using", ncores, "threads for fgsea (BPPARAM not needed; fgsea uses internal parallelism)\n\n")

# --- Paths ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
SIGS <- file.path(INT, "results/disease_signatures")
PROG <- file.path(INT, "results/progression")
OUTDIR <- file.path(RDIR, "gsea_kallisto")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", OUTDIR, "\n\n")

# ============================================================================
# 1. Load all contrasts
# ============================================================================
cat("=== Loading contrasts ===\n")

load_contrast <- function(file, contrast_name, t_col = "t", gene_col = "gene",
                          contrast_col = NULL) {
  if (!file.exists(file)) {
    cat("  WARNING: file not found, skipping:", file, "\n")
    return(NULL)
  }
  dt <- fread(file)
  # Standardise column names
  if (!t_col %in% names(dt)) {
    stop("t-statistic column '", t_col, "' not found in ", basename(file))
  }
  if (!gene_col %in% names(dt)) {
    stop("gene column '", gene_col, "' not found in ", basename(file))
  }

  # If file has a contrast column, split into multiple contrasts
  if (!is.null(contrast_col) && contrast_col %in% names(dt)) {
    contrasts <- unique(dt[[contrast_col]])
    out <- lapply(contrasts, function(ctr) {
      sub <- dt[get(contrast_col) == ctr]
      sub[, gene_base := gsub("\\..*", "", get(gene_col))]
      sub <- sub[order(-get(t_col))]
      # Remove duplicates (keep first = highest |t|)
      sub <- sub[!duplicated(gene_base)]
      ranks <- sub[[t_col]]
      names(ranks) <- sub$gene_base
      # Remove NAs and Inf
      ranks <- ranks[is.finite(ranks)]
      label <- paste0(contrast_name, "::", ctr)
      cat("  ", label, ":", length(ranks), "genes\n")
      list(name = label, ranks = ranks)
    })
    return(out)
  }

  # Single contrast
  dt[, gene_base := gsub("\\..*", "", get(gene_col))]
  dt <- dt[order(-get(t_col))]
  dt <- dt[!duplicated(gene_base)]
  ranks <- dt[[t_col]]
  names(ranks) <- dt$gene_base
  ranks <- ranks[is.finite(ranks)]
  cat("  ", contrast_name, ":", length(ranks), "genes\n")
  list(list(name = contrast_name, ranks = ranks))
}

# Accumulate all contrasts
all_contrasts <- list()

# C1: Disease vs Control
c1 <- load_contrast(file.path(RDIR, "dream_results.csv"),
                     "DiseaseVsControl")
all_contrasts <- c(all_contrasts, c1)

# C2: NASH vs NAFL (uses adj.P.Val not padj — but we only need t)
c2 <- load_contrast(file.path(SIGS, "nafl_vs_nash_dream.csv"),
                     "NASHvsNAFL")
all_contrasts <- c(all_contrasts, c2)

# Fibrosis cumulative (F1-F4 vs F0)
fib_cum <- load_contrast(file.path(SIGS, "fibrosis_stage_dream.csv"),
                          "FibCumul", contrast_col = "contrast")
if (!is.null(fib_cum)) all_contrasts <- c(all_contrasts, fib_cum)

# Fibrosis consecutive (F0->F1, F1->F2, ...)
fib_con <- load_contrast(file.path(SIGS, "fibrosis_consecutive_dream.csv"),
                          "FibConsec", contrast_col = "contrast")
if (!is.null(fib_con)) all_contrasts <- c(all_contrasts, fib_con)

# Progression contrasts (c3-c13)
progression_files <- list(
  c("c3_adv_vs_early_fib_dream.csv",    "C3_AdvVsEarlyFib"),
  c("c4_nafl_vs_ctrl_dream.csv",        "C4_NAFLvsCtrl"),
  c("c5_nas_ge5_vs_lt5_dream.csv",      "C5_NAS_ge5_vs_lt5"),
  c("c6_extreme_endpoints_dream.csv",   "C6_ExtremeEndpoints"),
  c("c8_cirrhosis_dream.csv",           "C8_Cirrhosis"),
  c("c9_f2_inflection_dream.csv",       "C9_F2Inflection"),
  c("c11_nash_vs_ctrl_dream.csv",       "C11_NASHvsCtrl"),
  c("c12_early_vs_late_nash_dream.csv", "C12_EarlyVsLateNASH"),
  c("c13_nash_vs_nafl_fib_adj_dream.csv", "C13_NASHvsNAFL_FibAdj")
)

for (pf in progression_files) {
  fpath <- file.path(PROG, pf[1])
  pc <- load_contrast(fpath, pf[2])
  if (!is.null(pc)) all_contrasts <- c(all_contrasts, pc)
}

# Ordinal fibrosis
ord <- load_contrast(file.path(PROG, "c17_fibrosis_ordinal_dream.csv"),
                      "C17_FibOrdinal")
if (!is.null(ord)) all_contrasts <- c(all_contrasts, ord)

cat("\nTotal contrasts loaded:", length(all_contrasts), "\n\n")

# ============================================================================
# 2. Load gene set collections (Ensembl IDs)
# ============================================================================
cat("=== Loading gene set collections ===\n")

load_geneset <- function(collection, subcollection = NULL, label) {
  cat("Loading", label, "...\n")
  args <- list(species = "Homo sapiens", collection = collection)
  if (!is.null(subcollection)) args$subcollection <- subcollection

  df <- tryCatch(do.call(msigdbr, args), error = function(e) {
    warning("Failed to load ", label, ": ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(df) || nrow(df) == 0) {
    cat("  WARNING:", label, "returned 0 rows. Skipping.\n")
    return(NULL)
  }

  # Use ensembl_gene to match dream gene IDs
  n_mapped <- sum(!is.na(df$ensembl_gene) & df$ensembl_gene != "")
  cat("  ", label, ":", length(unique(df$gs_name)), "sets,",
      n_mapped, "/", nrow(df), "members with Ensembl IDs (",
      round(100 * n_mapped / nrow(df), 1), "%)\n")

  sets <- split(df$ensembl_gene, df$gs_name)
  sets <- lapply(sets, function(x) unique(x[!is.na(x) & x != ""]))
  # Drop empty sets
  sets <- sets[sapply(sets, length) >= 1]
  list(sets = sets, label = label)
}

collections <- list(
  load_geneset("H", NULL, "Hallmark"),
  load_geneset("C2", "CP:KEGG_MEDICUS", "KEGG_MEDICUS"),
  load_geneset("C2", "CP:REACTOME", "Reactome"),
  load_geneset("C5", "GO:BP", "GO_BP"),
  load_geneset("C5", "GO:MF", "GO_MF"),
  load_geneset("C5", "GO:CC", "GO_CC")
)
# Remove NULLs
collections <- Filter(Negate(is.null), collections)

# If KEGG_MEDICUS failed, try KEGG legacy
kegg_loaded <- any(sapply(collections, function(x) x$label == "KEGG_MEDICUS"))
if (!kegg_loaded) {
  cat("KEGG_MEDICUS not available, trying legacy CP:KEGG...\n")
  kegg_legacy <- load_geneset("C2", "CP:KEGG_LEGACY", "KEGG_LEGACY")
  if (is.null(kegg_legacy)) {
    kegg_legacy <- load_geneset("C2", "CP:KEGG", "KEGG")
  }
  if (!is.null(kegg_legacy)) collections <- c(collections, list(kegg_legacy))
}

cat("\nCollections loaded:", length(collections), "\n\n")

# ============================================================================
# 3. Run fgsea: all contrasts x all collections
# ============================================================================
cat("=== Running fgsea ===\n")
cat("Parameters: nPermSimple=10000, minSize=15, maxSize=500, eps=0\n\n")

all_results <- list()
result_idx <- 0L

for (coll in collections) {
  cat("--- Collection:", coll$label, "(", length(coll$sets), "sets) ---\n")

  for (ctr in all_contrasts) {
    tryCatch({
      res <- fgsea(
        pathways    = coll$sets,
        stats       = ctr$ranks,
        minSize     = 15,
        maxSize     = 500,
        nPermSimple = 10000,
        eps         = 0
      )
      res$contrast   <- ctr$name
      res$collection <- coll$label

      n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
      if (n_sig > 0) {
        cat("  ", ctr$name, ":", n_sig, "sig /", nrow(res), "\n")
      }

      # Convert leadingEdge list to semicolon-separated string
      res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]

      result_idx <- result_idx + 1L
      all_results[[result_idx]] <- res
    }, error = function(e) {
      cat("  ERROR in", ctr$name, "/", coll$label, ":", conditionMessage(e), "\n")
    })
  }
}

cat("\nCombining results...\n")
gsea_all <- rbindlist(all_results, fill = TRUE)
cat("Total rows:", nrow(gsea_all), "\n")
cat("Significant (padj < 0.05):", sum(gsea_all$padj < 0.05, na.rm = TRUE), "\n\n")

# ============================================================================
# 4. Save main output
# ============================================================================
fwrite(gsea_all, file.path(OUTDIR, "gsea_all_contrasts.csv"))
cat("Saved: gsea_all_contrasts.csv\n")

# ============================================================================
# 5. NES heatmap matrix for Hallmark
# ============================================================================
cat("\n=== Building Hallmark NES heatmap matrix ===\n")

hallmark_rows <- gsea_all[collection == "Hallmark"]
if (nrow(hallmark_rows) > 0) {
  nes_wide <- dcast(hallmark_rows, pathway ~ contrast, value.var = "NES")
  fwrite(nes_wide, file.path(OUTDIR, "gsea_nes_heatmap_hallmark.csv"))
  cat("Saved: gsea_nes_heatmap_hallmark.csv (",
      nrow(nes_wide), "pathways x", ncol(nes_wide) - 1, "contrasts)\n")
} else {
  cat("WARNING: No Hallmark results found. NES heatmap skipped.\n")
}

# ============================================================================
# 6. GO semantic clustering via rrvgo
# ============================================================================
cluster_go <- function(gsea_dt, ont, out_file) {
  cat("\n--- Clustering", ont, "---\n")

  sig <- gsea_dt[collection == paste0("GO_", ont) & padj < 0.05]
  cat("Significant terms across all contrasts:", nrow(sig), "\n")

  if (nrow(sig) == 0) {
    cat("No significant terms; skipping clustering.\n")
    return(NULL)
  }

  if (!has_rrvgo) {
    cat("rrvgo not available; writing unclustered significant terms.\n")
    fwrite(sig, out_file)
    return(sig)
  }

  # Cluster per contrast (semantic similarity is GO-structure-based, not contrast-specific,
  # but NES/padj differ per contrast so we cluster within each)
  clustered_list <- list()
  for (ctr_name in unique(sig$contrast)) {
    sig_ctr <- sig[contrast == ctr_name]
    if (nrow(sig_ctr) < 2) {
      sig_ctr$cluster <- seq_len(nrow(sig_ctr))
      sig_ctr$is_representative <- TRUE
      sig_ctr$cluster_size <- 1L
      clustered_list[[ctr_name]] <- sig_ctr
      next
    }

    # Extract GO IDs from pathway names (GOBP_xxx -> GO:xxx mapping via rrvgo)
    # msigdbr pathway names are like GOBP_REGULATION_OF_xxx; we need GO IDs.
    # rrvgo::calculateSimMatrix needs GO IDs. Get the ID mapping from msigdbr.
    go_ont_map <- c("BP" = "BP", "MF" = "MF", "CC" = "CC")
    ont_full <- go_ont_map[ont]

    # Re-query msigdbr for GO IDs
    go_df <- msigdbr(species = "Homo sapiens", collection = "C5",
                     subcollection = paste0("GO:", ont))
    go_id_map <- unique(go_df[, c("gs_name", "gs_exact_source")])
    go_id_map <- setDT(go_id_map)

    sig_ctr <- merge(sig_ctr, go_id_map, by.x = "pathway", by.y = "gs_name", all.x = TRUE)
    sig_ctr <- sig_ctr[!is.na(gs_exact_source) & gs_exact_source != ""]

    if (nrow(sig_ctr) < 2) {
      sig_ctr$cluster <- seq_len(nrow(sig_ctr))
      sig_ctr$is_representative <- TRUE
      sig_ctr$cluster_size <- 1L
      clustered_list[[ctr_name]] <- sig_ctr
      next
    }

    # Deduplicate GO IDs (some pathways may share an ID after merge)
    sig_ctr <- sig_ctr[!duplicated(gs_exact_source)]

    # Calculate similarity matrix
    sim_mat <- tryCatch({
      calculateSimMatrix(sig_ctr$gs_exact_source,
                         orgdb    = "org.Hs.eg.db",
                         ont      = ont_full,
                         method   = "Wang")
    }, error = function(e) {
      cat("  rrvgo similarity failed for", ctr_name, ":", conditionMessage(e), "\n")
      return(NULL)
    })

    if (is.null(sim_mat) || nrow(sim_mat) < 2) {
      sig_ctr$cluster <- seq_len(nrow(sig_ctr))
      sig_ctr$is_representative <- TRUE
      sig_ctr$cluster_size <- 1L
      clustered_list[[ctr_name]] <- sig_ctr
      next
    }

    # Use -log10(padj) as scores (higher = more significant)
    scores <- setNames(-log10(pmax(sig_ctr$padj, 1e-300)), sig_ctr$gs_exact_source)
    # Keep only terms present in sim_mat
    common <- intersect(names(scores), rownames(sim_mat))
    scores <- scores[common]
    sig_ctr <- sig_ctr[gs_exact_source %in% common]

    reduced <- tryCatch({
      reduceSimMatrix(sim_mat, scores = scores, threshold = 0.7,
                      orgdb = "org.Hs.eg.db")
    }, error = function(e) {
      cat("  rrvgo reduceSimMatrix failed for", ctr_name, ":", conditionMessage(e), "\n")
      return(NULL)
    })

    if (!is.null(reduced)) {
      reduced_dt <- as.data.table(reduced)
      # Map cluster info back
      merge_cols <- intersect(c("go", "cluster", "parentTerm"), names(reduced_dt))
      if ("go" %in% merge_cols && "cluster" %in% merge_cols) {
        sig_ctr <- merge(sig_ctr, reduced_dt[, ..merge_cols],
                         by.x = "gs_exact_source", by.y = "go", all.x = TRUE)
      }
      # Identify representatives: rows where go == parent (self-referencing = cluster head)
      if ("parent" %in% names(reduced_dt)) {
        rep_ids <- reduced_dt[go == parent, go]
        sig_ctr[, is_representative := (gs_exact_source %in% rep_ids)]
      } else {
        # Fallback: pick highest-scoring term per cluster as representative
        sig_ctr[, is_representative := FALSE]
        if ("cluster" %in% names(sig_ctr)) {
          sig_ctr[, .is_best := padj == min(padj, na.rm = TRUE), by = cluster]
          sig_ctr[(.is_best), is_representative := TRUE]
          sig_ctr[, .is_best := NULL]
        } else {
          sig_ctr[, is_representative := TRUE]
        }
      }
      # Cluster sizes
      if ("cluster" %in% names(sig_ctr)) {
        sig_ctr[, cluster_size := .N, by = cluster]
      } else {
        sig_ctr[, cluster := seq_len(.N)]
        sig_ctr[, cluster_size := 1L]
        sig_ctr[, is_representative := TRUE]
      }
    } else {
      sig_ctr$cluster <- seq_len(nrow(sig_ctr))
      sig_ctr$is_representative <- TRUE
      sig_ctr$cluster_size <- 1L
    }

    clustered_list[[ctr_name]] <- sig_ctr
  }

  clustered <- rbindlist(clustered_list, fill = TRUE)
  fwrite(clustered, out_file)
  n_rep <- sum(clustered$is_representative, na.rm = TRUE)
  cat("Saved:", basename(out_file), "(", nrow(clustered), "rows,",
      n_rep, "cluster representatives)\n")
  return(clustered)
}

cluster_go(gsea_all, "BP", file.path(OUTDIR, "gsea_go_bp_clustered.csv"))
cluster_go(gsea_all, "MF", file.path(OUTDIR, "gsea_go_mf_clustered.csv"))

# ============================================================================
# 7. Summary counts
# ============================================================================
cat("\n=== Summary counts ===\n")

summary_dt <- gsea_all[, .(
  n_sig_up   = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
  n_sig_down = sum(padj < 0.05 & NES < 0, na.rm = TRUE),
  n_tested   = .N,
  top_pathway_up   = fifelse(any(padj < 0.05 & NES > 0, na.rm = TRUE),
                              pathway[which.max(fifelse(padj < 0.05 & NES > 0, NES, -Inf))],
                              NA_character_),
  top_pathway_down = fifelse(any(padj < 0.05 & NES < 0, na.rm = TRUE),
                              pathway[which.min(fifelse(padj < 0.05 & NES < 0, NES, Inf))],
                              NA_character_)
), by = .(contrast, collection)]

summary_dt[, n_sig_total := n_sig_up + n_sig_down]
setorder(summary_dt, contrast, collection)
fwrite(summary_dt, file.path(OUTDIR, "gsea_summary_counts.csv"))
cat("Saved: gsea_summary_counts.csv\n")

cat("\nSummary (sig up / sig down / tested):\n")
print(summary_dt[, .(contrast, collection, n_sig_up, n_sig_down, n_tested)], nrows = 100)

# ============================================================================
# Done
# ============================================================================
cat("\n=== 08b completed:", as.character(Sys.time()), "===\n")
cat("Output directory:", OUTDIR, "\n")
cat("Files:\n")
cat("  gsea_all_contrasts.csv\n")
cat("  gsea_nes_heatmap_hallmark.csv\n")
cat("  gsea_go_bp_clustered.csv\n")
cat("  gsea_go_mf_clustered.csv\n")
cat("  gsea_summary_counts.csv\n")
