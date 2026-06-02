#!/usr/bin/env Rscript
# 31_lincs_annotation.R
# ---------------------------------------------------------------------------
# LINCS L1000 Compound Annotation & Integration
#
# Post-processes Script 20 LINCS results to:
#   1. Map BRD IDs to compound names via LINCS RESTful API / local metadata
#   2. Annotate MOA using LINCS perturbation metadata
#   3. Cross-reference with DGIdb drug-gene interactions
#   4. Integrate with network proximity scores (Script 32, if available)
#   5. Produce final ranked compound table for manuscript
#
# Inputs:
#   - results/drug_repurposing/lincs_reversal_compounds.csv
#   - results/drug_repurposing/dgidb_drug_gene_interactions.csv
#   - results/drug_repurposing/pharmacotranscriptomics_summary.csv
#   - results/drug_repurposing/network_proximity/network_proximity_scores.csv (optional)
#   - causal_inference/causal_inference_summary.csv
#
# Outputs:
#   - results/drug_repurposing/lincs_annotated_compounds.csv
#   - results/drug_repurposing/lincs_final_ranked.csv
#   - figures/lincs_annotated_barplot.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)  # for rescale()
})

cat("=== Script 31: LINCS L1000 Annotation & Integration ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_DIR  <- file.path(BASE_DIR, "RNA-seq")
RESULTS_DIR <- file.path(RNASEQ_DIR, "results/drug_repurposing")
FIG_DIR     <- file.path(BASE_DIR, "figures")

LINCS_FILE  <- file.path(RESULTS_DIR, "lincs_reversal_compounds.csv")
LINCS_TOP50 <- file.path(RESULTS_DIR, "lincs_top50_reversals.csv")
DGIDB_FILE  <- file.path(RESULTS_DIR, "dgidb_drug_gene_interactions.csv")
PHARMA_FILE <- file.path(RESULTS_DIR, "pharmacotranscriptomics_summary.csv")
NETPROX_FILE <- file.path(RESULTS_DIR, "network_proximity/network_proximity_scores.csv")
MR_FILE     <- file.path(RNASEQ_DIR, "results/causal_inference/causal_inference_summary.csv")

# LINCS metadata (local cache)
LINCS_META_DIR <- file.path(BASE_DIR, "data/lincs_metadata")
dir.create(LINCS_META_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load LINCS reversal compounds
# ==============================================================================
cat("--- Step 1: Loading LINCS results ---\n")

if (!file.exists(LINCS_FILE)) {
  cat("  LINCS results not found:", LINCS_FILE, "\n")
  cat("  Run Script 20 first. Exiting.\n")
  quit(save = "no", status = 0)
}

lincs <- fread(LINCS_FILE)
cat("  LINCS reversal compounds:", nrow(lincs), "\n")
cat("  Columns:", paste(names(lincs), collapse = ", "), "\n")

# ==============================================================================
# 2. Map BRD IDs to compound names via LINCS API
# ==============================================================================
cat("\n--- Step 2: Mapping BRD IDs to compound names ---\n")

# Try local metadata cache first
PERT_CACHE <- file.path(LINCS_META_DIR, "lincs_pert_info.csv")

if (file.exists(PERT_CACHE)) {
  cat("  Loading cached LINCS perturbation metadata...\n")
  pert_info <- fread(PERT_CACHE)
} else {
  cat("  Querying LINCS API for BRD compound metadata...\n")
  cat("  (This queries the CLUE Connectivity Map API)\n")

  # LINCS CLUE API endpoint for perturbation info
  brd_ids <- unique(lincs$pert)
  cat("  Unique BRD IDs to query:", length(brd_ids), "\n")

  pert_info <- data.table(
    pert_id = character(),
    pert_iname = character(),
    moa = character(),
    target = character(),
    pubchem_cid = character()
  )

  # --- Strategy: Use local LINCS metadata (faster + more reliable than API) ---
  # Primary: GEO supplementary file (GSE92742) — already downloaded (4.7MB, 51K rows)
  # The uncompressed .txt file is the correct source; .txt.gz also available as backup
  PERT_META_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz"
  PERT_META_LOCAL <- file.path(BASE_DIR, "data/lincs_metadata", "GSE92742_Broad_LINCS_pert_info.txt")
  PERT_META_GZ   <- paste0(PERT_META_LOCAL, ".gz")

  if (!file.exists(PERT_META_LOCAL) || file.size(PERT_META_LOCAL) < 1e6) {
    # Try decompressing the .gz if available
    if (file.exists(PERT_META_GZ) && file.size(PERT_META_GZ) > 500e3) {
      cat("  Decompressing existing .gz file...\n")
      system2("gunzip", args = c("-k", PERT_META_GZ))
    }
    # If still missing, download
    if (!file.exists(PERT_META_LOCAL) || file.size(PERT_META_LOCAL) < 1e6) {
      cat("  Downloading LINCS perturbation metadata...\n")
      dir.create(dirname(PERT_META_LOCAL), recursive = TRUE, showWarnings = FALSE)
      tryCatch({
        download.file(PERT_META_URL, PERT_META_GZ, mode = "wb", quiet = FALSE)
        system2("gunzip", args = c("-f", PERT_META_GZ))
        fsize <- file.size(PERT_META_LOCAL)
        cat("  Downloaded and decompressed:", PERT_META_LOCAL, "(", round(fsize / 1e6, 1), "MB)\n")
      }, error = function(e) {
        cat("  WARNING: Could not download LINCS metadata:", conditionMessage(e), "\n")
        cat("  Download manually on login node:\n")
        cat("    wget -P", dirname(PERT_META_LOCAL), PERT_META_URL, "\n")
      })
    }
  } else {
    cat("  Using existing LINCS metadata:", PERT_META_LOCAL, "\n")
    cat("  File size:", round(file.size(PERT_META_LOCAL) / 1e6, 1), "MB\n")
  }

  if (file.exists(PERT_META_LOCAL)) {
    cat("  Loading local LINCS metadata...\n")
    pert_meta <- fread(PERT_META_LOCAL, sep = "\t", fill = TRUE)
    cat("  LINCS metadata:", nrow(pert_meta), "perturbations,",
        ncol(pert_meta), "columns\n")
    cat("  Columns:", paste(head(names(pert_meta), 10), collapse = ", "), "\n")

    # Map BRD IDs to compound info
    if ("pert_id" %in% names(pert_meta)) {
      matched <- pert_meta[pert_id %in% brd_ids]
      if (nrow(matched) > 0) {
        pert_info <- data.table(
          pert_id = matched$pert_id,
          pert_iname = if ("pert_iname" %in% names(matched)) matched$pert_iname else NA_character_,
          moa = if ("moa" %in% names(matched)) matched$moa else NA_character_,
          target = if ("target" %in% names(matched)) matched$target else NA_character_,
          pubchem_cid = if ("pubchem_cid" %in% names(matched)) as.character(matched$pubchem_cid) else NA_character_
        )
        cat("  Matched", nrow(pert_info), "/", length(brd_ids), "BRD IDs from local metadata\n")

        # If local file lacks moa/target columns (GEO format), try to fill from
        # signatureSearch output columns (MOAss, t_gn_sym) if available in lincs data
        if (all(is.na(pert_info$moa) | pert_info$moa == "") &&
            "MOAss" %in% names(lincs)) {
          cat("  Local metadata lacks MOA/target; backfilling from signatureSearch output...\n")
          ss_moa <- unique(lincs[!is.na(MOAss) & MOAss != "", .(pert, MOAss)])
          ss_tgt <- unique(lincs[!is.na(t_gn_sym) & t_gn_sym != "", .(pert, t_gn_sym)])
          if (nrow(ss_moa) > 0) {
            pert_info <- merge(pert_info, ss_moa, by.x = "pert_id", by.y = "pert", all.x = TRUE)
            pert_info[is.na(moa) & !is.na(MOAss), moa := MOAss]
            pert_info[, MOAss := NULL]
          }
          if (nrow(ss_tgt) > 0) {
            pert_info <- merge(pert_info, ss_tgt, by.x = "pert_id", by.y = "pert", all.x = TRUE)
            pert_info[is.na(target) & !is.na(t_gn_sym), target := t_gn_sym]
            pert_info[, t_gn_sym := NULL]
          }
          cat("  Compounds with MOA after backfill:",
              sum(!is.na(pert_info$moa) & pert_info$moa != ""), "\n")
        }
      }
    }
  }

  # Fallback: CLUE API for any remaining unmapped BRDs
  mapped_ids <- if (nrow(pert_info) > 0) pert_info$pert_id else character(0)
  unmapped_ids <- setdiff(brd_ids, mapped_ids)

  if (length(unmapped_ids) > 0) {
    cat("  Attempting CLUE API for", length(unmapped_ids), "unmapped BRDs...\n")
    HAS_JSONLITE <- requireNamespace("jsonlite", quietly = TRUE)
    if (HAS_JSONLITE) {
      library(jsonlite)
      batch_size <- 100
      n_batches <- ceiling(length(unmapped_ids) / batch_size)

      for (b in seq_len(n_batches)) {
        start_idx <- (b - 1) * batch_size + 1
        end_idx <- min(b * batch_size, length(unmapped_ids))
        batch_ids <- unmapped_ids[start_idx:end_idx]

        url <- paste0(
          "https://api.clue.io/api/perts?filter=",
          URLencode(paste0('{"where":{"pert_id":{"inq":["',
                            paste(batch_ids, collapse = '","'),
                            '"]}}}'), reserved = TRUE)
        )

        tryCatch({
          resp <- fromJSON(url, flatten = TRUE)
          if (is.data.frame(resp) && nrow(resp) > 0) {
            batch_dt <- data.table(
              pert_id = resp$pert_id,
              pert_iname = if ("pert_iname" %in% names(resp)) resp$pert_iname else NA_character_,
              moa = if ("moa" %in% names(resp)) sapply(resp$moa, function(x) paste(x, collapse = "; ")) else NA_character_,
              target = if ("target" %in% names(resp)) sapply(resp$target, function(x) paste(x, collapse = "; ")) else NA_character_,
              pubchem_cid = if ("pubchem_cid" %in% names(resp)) as.character(resp$pubchem_cid) else NA_character_
            )
            pert_info <- rbind(pert_info, batch_dt, fill = TRUE)
          }
          if (b %% 5 == 0) cat("  API Batch", b, "/", n_batches, "complete\n")
        }, error = function(e) {
          cat("  API error for batch", b, ":", conditionMessage(e), "\n")
        })

        Sys.sleep(0.5)
      }
    } else {
      cat("  jsonlite not available. Skipping API fallback.\n")
    }
  }

  if (nrow(pert_info) > 0) {
    fwrite(pert_info, PERT_CACHE)
    cat("  Cached", nrow(pert_info), "perturbation records\n")
  }
}

# Merge annotations into LINCS results
if (nrow(pert_info) > 0) {
  lincs <- merge(lincs, pert_info, by.x = "pert", by.y = "pert_id", all.x = TRUE)
  n_mapped <- sum(!is.na(lincs$pert_iname) & lincs$pert_iname != "")
  cat("  Mapped", n_mapped, "/", nrow(lincs), "compounds to names\n")
} else {
  lincs[, pert_iname := NA_character_]
  lincs[, moa := NA_character_]
  lincs[, target := NA_character_]
  cat("  No CLUE API results. BRD IDs retained.\n")
}

# ==============================================================================
# 3. Cross-reference with DGIdb
# ==============================================================================
cat("\n--- Step 3: Cross-referencing with DGIdb ---\n")

if (file.exists(DGIDB_FILE)) {
  dgidb <- fread(DGIDB_FILE)
  cat("  DGIdb entries:", nrow(dgidb), "\n")

  # Build drug -> target gene lookup from DGIdb
  # Handle both raw DGIdb format (drug_name/gene) and pharma summary (dgidb_drugs/symbol)
  dgidb_drugs <- list()
  drug_col <- if ("dgidb_drugs" %in% names(dgidb)) "dgidb_drugs" else
              if ("drug_name" %in% names(dgidb)) "drug_name" else NULL
  gene_col <- if ("symbol" %in% names(dgidb)) "symbol" else
              if ("gene" %in% names(dgidb)) "gene" else NULL

  if (!is.null(drug_col) && !is.null(gene_col)) {
    cat("  Using columns:", drug_col, "->", gene_col, "\n")
    for (i in seq_len(nrow(dgidb))) {
      drugs <- trimws(unlist(strsplit(as.character(dgidb[[drug_col]][i]), ";")))
      sym <- as.character(dgidb[[gene_col]][i])
      if (is.na(sym) || sym == "") next
      for (d in drugs) {
        if (is.na(d) || d == "") next
        dgidb_drugs[[toupper(d)]] <- unique(c(dgidb_drugs[[toupper(d)]], sym))
      }
    }
  } else {
    cat("  WARNING: Could not identify drug/gene columns in DGIdb file\n")
    cat("  Available columns:", paste(names(dgidb), collapse = ", "), "\n")
  }
  cat("  DGIdb drugs mapped:", length(dgidb_drugs), "\n")

  # Match LINCS compounds to DGIdb targets
  lincs[, dgidb_targets := NA_character_]
  lincs[, dgidb_n_targets := 0L]
  for (i in seq_len(nrow(lincs))) {
    compound <- toupper(if (!is.na(lincs$pert_iname[i])) lincs$pert_iname[i] else lincs$pert[i])
    if (compound %in% names(dgidb_drugs)) {
      targets <- dgidb_drugs[[compound]]
      lincs$dgidb_targets[i] <- paste(targets, collapse = ";")
      lincs$dgidb_n_targets[i] <- length(targets)
    }
  }
  cat("  LINCS compounds with DGIdb targets:", sum(lincs$dgidb_n_targets > 0), "\n")
} else {
  cat("  DGIdb file not found. Skipping DGIdb cross-reference.\n")
}

# ==============================================================================
# 4. Integrate with MR causal evidence
# ==============================================================================
cat("\n--- Step 4: Integrating MR causal evidence ---\n")

if (file.exists(MR_FILE)) {
  mr <- fread(MR_FILE)
  mr_causal <- mr[(!is.na(twas_p) & twas_p < 0.05) |
                   (!is.na(mr_p) & mr_p < 0.05), gene]
  cat("  MR/TWAS-nominal genes:", length(mr_causal), "\n")

  # Check if any LINCS targets overlap with MR-causal genes
  lincs[, mr_causal_overlap := FALSE]
  lincs[, mr_overlap_genes := NA_character_]
  for (i in seq_len(nrow(lincs))) {
    targets <- character(0)
    tgt_val <- as.character(lincs$target[i])
    if (length(tgt_val) == 1 && !is.na(tgt_val) && nchar(tgt_val) > 0) {
      targets <- c(targets, trimws(unlist(strsplit(tgt_val, ";"))))
    }
    dgi_val <- as.character(lincs$dgidb_targets[i])
    if (length(dgi_val) == 1 && !is.na(dgi_val) && nchar(dgi_val) > 0) {
      targets <- c(targets, trimws(unlist(strsplit(dgi_val, ";"))))
    }
    targets <- unique(targets[targets != "" & targets != "NA"])
    overlap <- intersect(targets, mr_causal)
    if (length(overlap) > 0) {
      lincs$mr_causal_overlap[i] <- TRUE
      lincs$mr_overlap_genes[i] <- paste(overlap, collapse = ";")
    }
  }
  cat("  LINCS compounds with MR-causal target overlap:",
      sum(lincs$mr_causal_overlap), "\n")
} else {
  cat("  MR results not found. Skipping causal integration.\n")
}

# ==============================================================================
# 5. Integrate network proximity (if available)
# ==============================================================================
cat("\n--- Step 5: Integrating network proximity ---\n")

if (file.exists(NETPROX_FILE)) {
  netprox <- fread(NETPROX_FILE)
  cat("  Network proximity results:", nrow(netprox), "drugs\n")

  # Match by drug name (LINCS: compounds may appear as "LINCS:BRD-..." in netprox)
  lincs[, net_z_combined := NA_real_]
  lincs[, net_p_closest := NA_real_]
  for (i in seq_len(nrow(lincs))) {
    # Try matching by BRD ID or compound name
    brd_key <- paste0("LINCS:", lincs$pert[i])
    name_key <- if (!is.na(lincs$pert_iname[i])) toupper(lincs$pert_iname[i]) else ""

    match_idx <- which(netprox$drug == brd_key |
                       toupper(netprox$drug) == name_key)
    if (length(match_idx) > 0) {
      lincs$net_z_combined[i] <- netprox$z_combined[match_idx[1]]
      lincs$net_p_closest[i] <- netprox$p_closest[match_idx[1]]
    }
  }
  n_matched <- sum(!is.na(lincs$net_z_combined))
  cat("  LINCS compounds with network proximity:", n_matched, "\n")
} else {
  cat("  Network proximity not yet computed. Run Script 32 first.\n")
}

# ==============================================================================
# 6. Compute composite drug score and rank
# ==============================================================================
cat("\n--- Step 6: Computing composite drug score ---\n")

# Ensure all columns exist (defensive — some may be missing if steps were skipped)
if (!"mr_causal_overlap" %in% names(lincs)) lincs[, mr_causal_overlap := FALSE]
if (!"dgidb_n_targets" %in% names(lincs)) lincs[, dgidb_n_targets := 0L]
if (!"net_z_combined" %in% names(lincs)) lincs[, net_z_combined := NA_real_]

# Composite score components:
#   1. LINCS reversal strength (|NCS| or |WTCS|, higher = better)
#   2. Statistical significance (WTCS_FDR)
#   3. MR-causal target overlap (binary bonus)
#   4. DGIdb target count (more targets = more actionable)
#   5. Network proximity z-score (more negative = closer)

lincs[, score_reversal := scales::rescale(-NCS, to = c(0, 1))]
lincs[, score_sig := scales::rescale(-log10(pmax(WTCS_FDR, 1e-10)), to = c(0, 1))]
lincs[, score_mr := ifelse(mr_causal_overlap, 1, 0)]
lincs[, score_dgidb := scales::rescale(pmin(dgidb_n_targets, 10), to = c(0, 1))]
lincs[, score_network := ifelse(!is.na(net_z_combined),
                                 scales::rescale(-net_z_combined, to = c(0, 1)),
                                 0)]

# Weighted composite (reversal 35%, significance 20%, MR 20%, DGIdb 15%, network 10%)
lincs[, composite_score := 0.35 * score_reversal +
                            0.20 * score_sig +
                            0.20 * score_mr +
                            0.15 * score_dgidb +
                            0.10 * score_network]

setorder(lincs, -composite_score)

# Display name: prefer pert_iname, fall back to BRD ID
lincs[, display_name := ifelse(!is.na(pert_iname) & pert_iname != "",
                                pert_iname, pert)]

# ==============================================================================
# 7. Save annotated and ranked results
# ==============================================================================
cat("\n--- Step 7: Saving results ---\n")

fwrite(lincs, file.path(RESULTS_DIR, "lincs_annotated_compounds.csv"))
cat("  Saved: lincs_annotated_compounds.csv (", nrow(lincs), "rows)\n")

# Top 50 ranked
top50 <- head(lincs, 50)
fwrite(top50, file.path(RESULTS_DIR, "lincs_final_ranked.csv"))
cat("  Saved: lincs_final_ranked.csv (", nrow(top50), "rows)\n")

# ==============================================================================
# 8. Visualization
# ==============================================================================
cat("\n--- Step 8: Generating figures ---\n")

# Top 20 annotated barplot
top20 <- head(lincs, 20)
top20[, display_name := factor(display_name, levels = rev(display_name))]

if (nrow(top20) > 0) {
  pdf(file.path(FIG_DIR, "lincs_annotated_top20.pdf"), width = 10, height = 7)

  p <- ggplot(top20, aes(x = composite_score, y = display_name,
                          fill = mr_causal_overlap)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#D7191C", "FALSE" = "#2C7BB6"),
                      name = "MR-causal\noverlap",
                      labels = c("TRUE" = "Yes", "FALSE" = "No")) +
    labs(title = "Top 20 LINCS L1000 Reversal Compounds",
         subtitle = "Composite score: reversal (35%) + significance (20%) + MR (20%) + DGIdb (15%) + network (10%)",
         x = "Composite Drug Score",
         y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          panel.grid.major.y = element_blank())

  print(p)
  dev.off()
  cat("  Saved: figures/lincs_annotated_top20.pdf\n")
}

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Script 31: LINCS Annotation Complete ===\n")
cat("  Total compounds annotated:", nrow(lincs), "\n")
cat("  Compounds with names:", sum(!is.na(lincs$pert_iname) & lincs$pert_iname != ""), "\n")
cat("  Compounds with MOA:", sum(!is.na(lincs$moa) & lincs$moa != ""), "\n")
cat("  Compounds with MR-causal overlap:", sum(lincs$mr_causal_overlap, na.rm = TRUE), "\n")
cat("  Compounds with DGIdb targets:", sum(lincs$dgidb_n_targets > 0, na.rm = TRUE), "\n")
cat("  Compounds with network proximity:", sum(!is.na(lincs$net_z_combined)), "\n")
cat("  Top compound:", top20$display_name[1], "(score:", round(top20$composite_score[1], 3), ")\n")
cat("End time:", format(Sys.time()), "\n")
