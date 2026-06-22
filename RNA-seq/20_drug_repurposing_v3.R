#!/usr/bin/env Rscript
# 20_drug_repurposing_v3.R
#
# Strategy 11 v3: Comprehensive Pharmacotranscriptomics for MASLD Drug Repurposing
#
# Two-arm analysis:
#   Arm 1:  LINCS L1000 signature reversal (HepG2 cell line, signatureSearch)
#   Arm 2a: DGIdb drug-gene interactions (GraphQL API, no auth)
#   Arm 2b: Open Targets known drugs for causal genes (GraphQL API, no auth)
#   Plus:   C2:CGP disease concordance (retained from v2)
#   Integration: convergence scoring, MR overlap, composite supplementary figure
#
# Inputs:
#   - canonical_deg_results.csv  (limma-voom-qw C2 mega-analysis)
#   - consensus_degs.csv         (tiered consensus DEGs)
#   - causal_inference_summary.csv (MR + TWAS results)
#   - human_ensg_to_symbol.tsv   (annotation cache)
#
# Outputs (results/drug_repurposing/):
#   - lincs_reversal_compounds.csv      Arm 1: all LINCS HepG2 reversal hits
#   - lincs_top50_reversals.csv         Arm 1: top 50 compounds
#   - dgidb_drug_gene_interactions.csv  Arm 2a: DGIdb drug-gene pairs
#   - opentargets_known_drugs.csv       Arm 2b: Open Targets drug annotations
#   - cgp_disease_concordance.csv       C2:CGP fgsea full results
#   - cgp_reversal_hits.csv             C2:CGP reversal hits
#   - convergent_drug_targets.csv       Integration: LINCS + DGIdb intersection
#   - pharmacotranscriptomics_summary.csv  Master gene table, all layers
#   - mr_convergent_drug_targets.csv    MR-causal genes with drug evidence
#   - figures/figS_pharmacotranscriptomics.pdf  Composite supplementary figure

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(httr2)
  library(jsonlite)
})

cat("=== Strategy 11 v3: Comprehensive Pharmacotranscriptomics ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# Check optional packages
HAS_SIGSEARCH <- requireNamespace("signatureSearch", quietly = TRUE) &&
                 requireNamespace("signatureSearchData", quietly = TRUE)
HAS_COMPLEXHM <- requireNamespace("ComplexHeatmap", quietly = TRUE)
HAS_ORGDB     <- requireNamespace("org.Hs.eg.db", quietly = TRUE)

cat("Package availability:\n")
cat("  signatureSearch:", ifelse(HAS_SIGSEARCH, "YES", "NO (Arm 1 will be skipped)"), "\n")
cat("  ComplexHeatmap: ", ifelse(HAS_COMPLEXHM, "YES", "NO (UpSet panel will use bar chart fallback)"), "\n")
cat("  org.Hs.eg.db:   ", ifelse(HAS_ORGDB, "YES", "NO (Entrez mapping limited)"), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- Sys.getenv("MASLD_PROJECT_ROOT",
                 "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RNASEQ_DIR  <- file.path(BASE_DIR, "RNA-seq")
RESULTS_DIR <- file.path(RNASEQ_DIR, "results/drug_repurposing")
FIG_DIR     <- file.path(BASE_DIR, "figures")
INT_DIR     <- file.path(RNASEQ_DIR,
  "Human/Patient_Cohorts/analysis/integration/results/integration")

# NOTE: Uses canonical_deg_results.csv (limma-voom-qw C2) for t-statistic ranking
# (fgsea + LINCS). t-statistics drive the disease query signature. DEG significance
# comes from consensus_degs.csv bulk_sig column (padj < 0.05, |logFC| > 0.5).
DREAM_FILE     <- file.path(INT_DIR, "canonical_deg_results.csv")
CONSENSUS_FILE <- file.path(INT_DIR, "consensus_degs.csv")
MR_FILE        <- file.path(RNASEQ_DIR, "results/causal_inference/causal_inference_summary.csv")
ANN_CACHE_PATH <- file.path(RNASEQ_DIR,
  "Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")

FGSEA_PADJ_THRESHOLD <- 0.05
FGSEA_MIN_SIZE       <- 15
FGSEA_MAX_SIZE       <- 500
LINCS_N_QUERY        <- 150
SLURM_CPUS           <- as.integer(Sys.getenv("SLURM_CPUS", "8"))

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Load publication theme
theme_file <- file.path(BASE_DIR, "scripts/figures/publication_theme.R")
if (file.exists(theme_file)) {
  source(theme_file)
  cat("  Loaded publication theme\n")
} else {
  cat("  WARNING: publication_theme.R not found; using theme_bw fallback\n")
  theme_masld <- function(...) theme_bw(base_size = 7)
  save_fig <- function(plot, filename, width = 7.1, height = 4, dpi = 300) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    ggsave(filename, plot, width = width, height = height, dpi = dpi)
    ggsave(sub("\\.pdf$", ".png", filename), plot, width = width, height = height, dpi = dpi)
  }
  masld_colors <- list(up = "#B2182B", down = "#2166AC", ns = "#CCCCCC",
                       deg = "#2166AC", mr = "#B2182B", twas = "#D6604D",
                       gwas = "#4393C3", deconv = "#92C5DE")
  fig_full_width <- 7.1  # 180mm / 25.4
}

# ==============================================================================
# Part A: Build Disease Signature & Gene ID Maps
# ==============================================================================
cat("\n--- Part A: Building disease signature & gene ID maps ---\n")

# A1. Load dream results
dream <- fread(DREAM_FILE)
cat("  Loaded dream results:", nrow(dream), "genes\n")

dream[, ensembl_id := sub("\\.\\d+$", "", gene)]
setorder(dream, P.Value)
dream_dedup <- dream[!duplicated(ensembl_id)]
cat("  After deduplication:", nrow(dream_dedup), "unique Ensembl IDs\n")

# A2. Build t-statistic ranking (for fgsea + LINCS query)
ranked_stats <- dream_dedup$t
names(ranked_stats) <- dream_dedup$ensembl_id
ranked_stats <- ranked_stats[!is.na(ranked_stats)]
ranked_stats <- sort(ranked_stats, decreasing = TRUE)
cat("  Ranked gene list:", length(ranked_stats), "genes (by t-statistic)\n")
cat("  t-stat range:", round(min(ranked_stats), 3), "to", round(max(ranked_stats), 3), "\n")

# A3. Build ID mappings: Ensembl → symbol, Ensembl → Entrez
cat("  Building ID mappings...\n")

# Ensembl → symbol (local annotation cache)
if (file.exists(ANN_CACHE_PATH)) {
  cat("    Using local annotation cache for Ensembl → symbol...\n")
  ann_cache <- fread(ANN_CACHE_PATH)
  ensembl_to_symbol <- ann_cache[symbol != "" & !is.na(symbol),
                                  .(ensembl_id = gene_base, symbol)]
  ensembl_to_symbol <- ensembl_to_symbol[!duplicated(ensembl_id)]
  cat("    Mapped", nrow(ensembl_to_symbol), "Ensembl IDs to HGNC symbols\n")
} else {
  cat("    WARNING: Annotation cache not found at:", ANN_CACHE_PATH, "\n")
  cat("    Falling back to biomaRt...\n")
  ensembl_to_symbol <- tryCatch({
    library(biomaRt)
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                       mirror = "useast")
    bm <- as.data.table(getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters    = "ensembl_gene_id",
      values     = names(ranked_stats),
      mart       = mart
    ))
    setnames(bm, c("ensembl_id", "symbol"))
    bm[symbol != "" & !duplicated(ensembl_id)]
  }, error = function(e) {
    cat("    ERROR: biomaRt failed:", conditionMessage(e), "\n")
    data.table(ensembl_id = character(), symbol = character())
  })
  cat("    Mapped", nrow(ensembl_to_symbol), "Ensembl IDs to HGNC symbols\n")
}

# Ensembl → Entrez (for LINCS signatureSearch which uses Entrez)
if (HAS_ORGDB) {
  cat("    Building Ensembl → Entrez mapping via org.Hs.eg.db...\n")
  ensembl_to_entrez <- tryCatch({
    eg <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                 keys = names(ranked_stats),
                                 keytype = "ENSEMBL",
                                 column = "ENTREZID",
                                 multiVals = "first")
    dt <- data.table(ensembl_id = names(eg), entrez_id = unname(eg))
    dt[!is.na(entrez_id)]
  }, error = function(e) {
    cat("    WARNING: org.Hs.eg.db mapping failed:", conditionMessage(e), "\n")
    data.table(ensembl_id = character(), entrez_id = character())
  })
  cat("    Mapped", nrow(ensembl_to_entrez), "Ensembl IDs to Entrez IDs\n")
} else {
  cat("    org.Hs.eg.db not available; Entrez mapping skipped\n")
  ensembl_to_entrez <- data.table(ensembl_id = character(), entrez_id = character())
}

# A4. Load consensus DEGs for DGIdb query
cat("  Loading consensus DEGs...\n")
consensus <- fread(CONSENSUS_FILE)
cat("  Consensus DEGs:", nrow(consensus), "total\n")
cat("    Significant DEGs:", nrow(consensus[bulk_sig == TRUE]), "\n")

# Upregulated highly significant DEGs (bulk_dir == 1)
up_consensus <- consensus[bulk_sig == TRUE & bulk_dir == 1]
cat("  Upregulated Significant DEGs:", nrow(up_consensus), "genes\n")

# Map to symbols
up_consensus[, ensembl_id := sub("\\.\\d+$", "", gene)]
# Use symbol from consensus_degs.csv if available; otherwise merge from annotation
if ("symbol" %in% names(up_consensus) && sum(!is.na(up_consensus$symbol) & up_consensus$symbol != "") > 0) {
  up_with_sym <- copy(up_consensus)
} else {
  up_with_sym <- merge(up_consensus[, setdiff(names(up_consensus), "symbol"), with = FALSE],
                       ensembl_to_symbol, by = "ensembl_id", all.x = TRUE)
}
up_symbols  <- unique(up_with_sym[!is.na(symbol) & symbol != "", symbol])
cat("  With HGNC symbols:", length(up_symbols), "genes → DGIdb input\n\n")


# ==============================================================================
# Part B: Arm 1 — LINCS L1000 Signature Reversal
# ==============================================================================
lincs_results     <- NULL
lincs_results_all <- NULL
sig_reversals     <- NULL

if (HAS_SIGSEARCH) {
  cat("--- Part B: Arm 1 — LINCS L1000 Signature Reversal ---\n")
  suppressPackageStartupMessages({
    library(signatureSearch)
    library(signatureSearchData)
    library(ExperimentHub)
    library(rhdf5)
  })

  tryCatch({
    # B1. Load LINCS2 database
    cat("  Loading LINCS L1000 reference database (EH7297)...\n")
    cat("  (First run downloads ~10GB HDF5; subsequent runs use cache)\n")
    eh <- ExperimentHub()
    lincs_db_path <- eh[["EH7297"]]
    cat("  LINCS2 database loaded:", lincs_db_path, "\n")

    # B2. Identify HepG2 signatures for ref_trts filter (if HDF5 structure allows)
    hepg2_trts <- NULL
    tryCatch({
      cat("  Scanning HDF5 for HepG2 signatures...\n")
      h5_contents <- h5ls(lincs_db_path, recursive = FALSE)
      if ("colnames" %in% h5_contents$name) {
        sig_names <- h5read(lincs_db_path, "colnames")
        hepg2_idx <- grep("HEPG2", sig_names, ignore.case = TRUE, value = TRUE)
        cat("  Total signatures:", length(sig_names), "\n")
        cat("  HepG2 signatures:", length(hepg2_idx), "\n")
        if (length(hepg2_idx) >= 10) {
          hepg2_trts <- hepg2_idx
        } else {
          # Expand to all liver cell lines
          liver_idx <- grep("HEPG2|HUH7|HEP3B", sig_names, ignore.case = TRUE, value = TRUE)
          if (length(liver_idx) >= 10) hepg2_trts <- liver_idx
          cat("  Liver cell line signatures:", length(liver_idx), "\n")
        }
      } else {
        cat("  HDF5 has no 'colnames' dataset; will filter post-hoc\n")
      }
    }, error = function(e) {
      cat("  Could not scan HDF5 columns:", conditionMessage(e), "\n")
      cat("  Will run full query and filter post-hoc by cell line\n")
    })

    # B3. Build query: top 150 up + 150 down by t-stat → Entrez IDs
    top_up_ens   <- names(head(ranked_stats, LINCS_N_QUERY))
    top_down_ens <- names(tail(ranked_stats, LINCS_N_QUERY))

    # Map to Entrez for signatureSearch
    up_entrez <- ensembl_to_entrez[ensembl_id %in% top_up_ens, entrez_id]
    dn_entrez <- ensembl_to_entrez[ensembl_id %in% top_down_ens, entrez_id]
    cat("  Query genes: ", length(up_entrez), "up /", length(dn_entrez), "down (Entrez)\n")

    if (length(up_entrez) >= 10 && length(dn_entrez) >= 10) {
      # B4. Create qSig and run GESS
      qsig <- qSig(query = list(upset = up_entrez, downset = dn_entrez),
                    gess_method = "LINCS",
                    refdb = lincs_db_path)

      cat("  Running GESS LINCS query (may take 10-30 min)...\n")
      cat("  Using", SLURM_CPUS, "workers\n")
      if (!is.null(hepg2_trts)) {
        cat("  Restricting to", length(hepg2_trts), "liver cell line signatures via ref_trts\n")
        lincs_res <- gess_lincs(qSig = qsig, sortby = "NCS", tau = TRUE,
                                 ref_trts = hepg2_trts, workers = SLURM_CPUS)
      } else {
        lincs_res <- gess_lincs(qSig = qsig, sortby = "NCS", tau = TRUE,
                                 workers = SLURM_CPUS)
      }
      lincs_dt <- as.data.table(result(lincs_res))
      cat("  LINCS query complete:", nrow(lincs_dt), "perturbations tested\n")

      # B4b. Map BRD IDs → drug names via CLUE Repurposing Hub compound info
      clue_cache <- file.path(RESULTS_DIR, "clue_compoundinfo_beta.txt")
      if (!file.exists(clue_cache)) {
        cat("  Downloading CLUE Repurposing Hub compound info...\n")
        clue_url <- "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt"
        tryCatch({
          download.file(clue_url, clue_cache, quiet = TRUE, mode = "w")
          cat("  Downloaded:", clue_cache, "\n")
        }, error = function(e) {
          cat("  WARNING: CLUE download failed:", conditionMessage(e), "\n")
          cat("  BRD IDs will be retained as-is.\n")
        })
      }
      if (file.exists(clue_cache)) {
        clue_info <- fread(clue_cache, sep = "\t", quote = "")
        cat("  CLUE compound info:", nrow(clue_info), "entries\n")
        # Key columns: pert_id (BRD), cmap_name, moa, target
        clue_cols <- intersect(names(clue_info), c("pert_id", "cmap_name", "moa", "target"))
        if ("pert_id" %in% clue_cols && "cmap_name" %in% clue_cols) {
          clue_map <- clue_info[, .SD, .SDcols = clue_cols]
          clue_map <- clue_map[!duplicated(pert_id)]
          lincs_dt <- merge(lincs_dt, clue_map, by.x = "pert", by.y = "pert_id", all.x = TRUE)
          n_mapped <- sum(!is.na(lincs_dt$cmap_name) & lincs_dt$cmap_name != "")
          cat("  BRD→drug name mapped:", n_mapped, "/", nrow(lincs_dt), "perturbations\n")
        } else {
          cat("  WARNING: CLUE file missing pert_id/cmap_name columns; BRD IDs retained\n")
          lincs_dt[, cmap_name := NA_character_]
        }
      } else {
        lincs_dt[, cmap_name := NA_character_]
      }

      # Filter to HepG2 / liver cell lines
      if ("cell" %in% names(lincs_dt)) {
        liver_hits <- lincs_dt[grepl("HEPG2|HUH7|HEP3B", cell, ignore.case = TRUE)]
        cat("  Liver cell line hits:", nrow(liver_hits), "\n")
      } else {
        liver_hits <- lincs_dt
      }

      # B5. Classify reversals using WTCS_FDR (Tau is empty/unreliable)
      lincs_reversals_all <- liver_hits[NCS < 0][order(NCS)]
      lincs_reversals <- liver_hits[NCS < 0 & !is.na(WTCS_FDR) & WTCS_FDR < 0.25][order(NCS)]
      sig_reversals   <- liver_hits[NCS < 0 & !is.na(WTCS_FDR) & WTCS_FDR < 0.05][order(NCS)]

      cat("  All reversal hits (NCS < 0):", nrow(lincs_reversals_all), "\n")
      cat("  Significant reversals (WTCS_FDR < 0.25):", nrow(lincs_reversals), "\n")
      cat("  Stringent reversals (WTCS_FDR < 0.05):", nrow(sig_reversals), "\n")

      # B6. Save
      fwrite(lincs_reversals, file.path(RESULTS_DIR, "lincs_reversal_compounds.csv"))
      cat("  Saved: lincs_reversal_compounds.csv\n")

      lincs_top50 <- head(lincs_reversals, 50)
      fwrite(lincs_top50, file.path(RESULTS_DIR, "lincs_top50_reversals.csv"))
      cat("  Saved: lincs_top50_reversals.csv\n")

      lincs_results     <- lincs_reversals
      lincs_results_all <- liver_hits  # keep full set for figure panel
    } else {
      cat("  WARNING: Not enough mapped Entrez IDs for LINCS query.\n")
      cat("  Skipping LINCS analysis.\n")
    }
  }, error = function(e) {
    cat("  ERROR in signatureSearch:", conditionMessage(e), "\n")
    cat("  Continuing without LINCS results.\n")
  })
} else {
  cat("--- Part B: Arm 1 — LINCS SKIPPED (signatureSearch not installed) ---\n")
  cat("  Install with: BiocManager::install(c('signatureSearch', 'signatureSearchData'))\n\n")
}


# ==============================================================================
# Part C: Arm 2a — DGIdb Drug-Gene Interactions
# ==============================================================================
cat("\n--- Part C: Arm 2a — DGIdb Drug-Gene Interactions ---\n")

dgidb_results <- tryCatch({
  DGIDB_URL <- "https://dgidb.org/api/graphql"
  BATCH_SIZE <- 50
  SLEEP_SEC  <- 0.5

  # DGIdb GraphQL query
  dgidb_query_template <- '
  query($genes: [String!]!) {
    genes(names: $genes) {
      nodes {
        name
        conceptId
        interactions {
          interactionScore
          interactionTypes {
            type
            directionality
          }
          drug {
            name
            conceptId
            approved
          }
          interactionAttributes {
            name
            value
          }
          publications {
            pmid
          }
          sources {
            fullName
          }
        }
      }
    }
  }'

  cat("  Querying DGIdb for", length(up_symbols), "genes in batches of", BATCH_SIZE, "...\n")

  all_interactions <- list()
  batches <- split(up_symbols, ceiling(seq_along(up_symbols) / BATCH_SIZE))

  for (i in seq_along(batches)) {
    batch <- batches[[i]]
    cat("    Batch", i, "/", length(batches), "(", length(batch), "genes)...")

    resp <- tryCatch({
      req <- request(DGIDB_URL) |>
        req_headers("Content-Type" = "application/json") |>
        req_body_json(list(query = dgidb_query_template,
                           variables = list(genes = batch))) |>
        req_timeout(60) |>
        req_retry(max_tries = 3, backoff = ~ 2)
      resp_body_json(req_perform(req))
    }, error = function(e) {
      cat(" ERROR:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(resp) && !is.null(resp$data$genes$nodes)) {
      nodes <- resp$data$genes$nodes
      for (node in nodes) {
        gene_name <- node$name
        for (ix in node$interactions) {
          drug_name    <- ix$drug$name %||% NA_character_
          drug_id      <- ix$drug$conceptId %||% NA_character_
          approved     <- ix$drug$approved %||% FALSE
          score        <- ix$interactionScore %||% NA_real_

          # Extract interaction type and directionality
          itype <- NA_character_
          idir  <- NA_character_
          if (length(ix$interactionTypes) > 0) {
            itype <- paste(vapply(ix$interactionTypes,
                                   function(x) x$type %||% "", character(1)),
                           collapse = ";")
            idir  <- paste(vapply(ix$interactionTypes,
                                   function(x) x$directionality %||% "", character(1)),
                           collapse = ";")
          }

          # Source databases
          sources <- if (length(ix$sources) > 0) {
            paste(vapply(ix$sources, function(x) x$fullName %||% "", character(1)),
                  collapse = ";")
          } else NA_character_

          # Number of publications
          n_pubs <- length(ix$publications)

          all_interactions[[length(all_interactions) + 1L]] <- data.table(
            gene          = gene_name,
            drug_name     = drug_name,
            drug_id       = drug_id,
            approved      = approved,
            interaction_score = score,
            interaction_type  = itype,
            directionality    = idir,
            n_publications    = n_pubs,
            sources           = sources
          )
        }
      }
      cat(" OK (", length(nodes), "genes returned)\n")
    } else {
      cat(" no data\n")
    }

    if (i < length(batches)) Sys.sleep(SLEEP_SEC)
  }

  if (length(all_interactions) > 0) {
    dgidb_dt <- rbindlist(all_interactions)
    cat("  DGIdb results:", nrow(dgidb_dt), "drug-gene pairs\n")
    cat("  Unique genes with drugs:", uniqueN(dgidb_dt$gene), "\n")
    cat("  Unique drugs:", uniqueN(dgidb_dt$drug_name), "\n")
    cat("  Approved drugs:", sum(dgidb_dt$approved, na.rm = TRUE), "interactions\n")

    # Flag therapeutically relevant interactions for upregulated genes
    # (inhibitor/antagonist = can suppress overexpressed target)
    dgidb_dt[, therapeutic := grepl("inhibitor|antagonist|blocker|suppressor|negative",
                                    interaction_type, ignore.case = TRUE)]

    fwrite(dgidb_dt, file.path(RESULTS_DIR, "dgidb_drug_gene_interactions.csv"))
    cat("  Saved: dgidb_drug_gene_interactions.csv\n")
    dgidb_dt
  } else {
    cat("  No DGIdb interactions found.\n")
    data.table(gene = character(), drug_name = character(),
               approved = logical(), interaction_type = character(),
               therapeutic = logical())
  }
}, error = function(e) {
  cat("  ERROR in DGIdb query:", conditionMessage(e), "\n")
  cat("  Continuing without DGIdb results.\n")
  data.table(gene = character(), drug_name = character(),
             approved = logical(), interaction_type = character(),
             therapeutic = logical())
})


# ==============================================================================
# Part D: Arm 2b — Open Targets Known Drugs
# ==============================================================================
cat("\n--- Part D: Arm 2b — Open Targets Known Drugs ---\n")

opentargets_results <- tryCatch({
  OT_URL <- "https://api.platform.opentargets.org/api/v4/graphql"

  # Select genes for OT query: MR/TWAS-causal + DGIdb-druggable Tier 1+2
  causal_genes <- data.table()
  if (file.exists(MR_FILE)) {
    mr_dt <- fread(MR_FILE)

    # Primary: MR/TWAS-causal genes
    causal_ens <- mr_dt[
      (!is.na(twas_fdr) & twas_fdr < 0.05) |
      (!is.na(mr_p) & mr_p < 0.05) |
      (!is.na(twas_p) & twas_p < 0.001),
      ensembl_id
    ]
    cat("  MR/TWAS-causal genes:", length(causal_ens), "\n")

    # Secondary: DGIdb-druggable Significant genes (from Part C results)
    druggable_symbols <- character(0)
    if (nrow(dgidb_results) > 0) {
      druggable_symbols <- unique(dgidb_results$gene)
    }
    druggable_ens <- ensembl_to_symbol[symbol %in% druggable_symbols, ensembl_id]
    cat("  DGIdb-druggable genes:", length(druggable_ens), "\n")

    # Combine and deduplicate
    ot_query_ens <- unique(c(causal_ens, druggable_ens))
    causal_genes <- mr_dt[ensembl_id %in% ot_query_ens]
    # Also add DGIdb-only genes that may not be in mr_dt
    dgidb_only_ens <- setdiff(druggable_ens, mr_dt$ensembl_id)
    if (length(dgidb_only_ens) > 0) {
      dgidb_only_dt <- ensembl_to_symbol[ensembl_id %in% dgidb_only_ens,
                                          .(ensembl_id, gene = symbol)]
      causal_genes <- rbindlist(list(causal_genes, dgidb_only_dt), fill = TRUE)
    }
    causal_genes <- causal_genes[!duplicated(ensembl_id)]
    cat("  Total genes for OT query:", nrow(causal_genes), "\n")
  } else {
    cat("  MR results not found; skipping Open Targets\n")
  }

  if (nrow(causal_genes) > 0) {
    ot_query <- '
    query($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        knownDrugs {
          uniqueDrugs
          rows {
            drug {
              name
              drugType
              maximumClinicalTrialPhase
              isApproved
              mechanismsOfAction {
                rows {
                  mechanismOfAction
                  actionType
                }
              }
            }
            disease {
              name
              id
            }
            phase
            status
            urls {
              url
              name
            }
          }
        }
      }
    }'

    cat("  Querying Open Targets for", nrow(causal_genes), "targets...\n")
    all_ot <- list()

    for (i in seq_len(nrow(causal_genes))) {
      ens_id <- causal_genes$ensembl_id[i]
      gene_sym <- causal_genes$gene[i]

      if (i %% 20 == 0) cat("    Processed", i, "/", nrow(causal_genes), "\n")

      resp <- tryCatch({
        req <- request(OT_URL) |>
          req_headers("Content-Type" = "application/json") |>
          req_body_json(list(query = ot_query,
                             variables = list(ensemblId = ens_id))) |>
          req_timeout(30) |>
          req_retry(max_tries = 3, backoff = ~ 2)
        resp_body_json(req_perform(req))
      }, error = function(e) NULL)

      if (!is.null(resp) && !is.null(resp$data$target$knownDrugs$rows)) {
        rows <- resp$data$target$knownDrugs$rows
        for (r in rows) {
          drug_name <- r$drug$name %||% NA_character_
          drug_type <- r$drug$drugType %||% NA_character_
          max_phase <- r$drug$maximumClinicalTrialPhase %||% NA_integer_
          is_approved <- r$drug$isApproved %||% FALSE
          disease_name <- r$disease$name %||% NA_character_
          disease_id   <- r$disease$id %||% NA_character_
          trial_phase  <- r$phase %||% NA_integer_
          trial_status <- r$status %||% NA_character_

          # Mechanism of action
          moa <- NA_character_
          action_type <- NA_character_
          if (!is.null(r$drug$mechanismsOfAction$rows) &&
              length(r$drug$mechanismsOfAction$rows) > 0) {
            moa <- paste(vapply(r$drug$mechanismsOfAction$rows,
                                 function(x) x$mechanismOfAction %||% "", character(1)),
                         collapse = "; ")
            action_type <- paste(vapply(r$drug$mechanismsOfAction$rows,
                                         function(x) x$actionType %||% "", character(1)),
                                 collapse = "; ")
          }

          all_ot[[length(all_ot) + 1L]] <- data.table(
            gene           = gene_sym,
            ensembl_id     = ens_id,
            drug_name      = drug_name,
            drug_type      = drug_type,
            mechanism       = moa,
            action_type     = action_type,
            max_clinical_phase = max_phase,
            is_approved     = is_approved,
            disease_name    = disease_name,
            disease_id      = disease_id,
            trial_phase     = trial_phase,
            trial_status    = trial_status
          )
        }
      }

      Sys.sleep(0.2)
    }

    if (length(all_ot) > 0) {
      ot_dt <- rbindlist(all_ot)
      cat("  Open Targets results:", nrow(ot_dt), "drug-target-disease rows\n")
      cat("  Unique targets:", uniqueN(ot_dt$gene), "\n")
      cat("  Unique drugs:", uniqueN(ot_dt$drug_name), "\n")
      cat("  Approved drugs:", sum(ot_dt$is_approved, na.rm = TRUE), "rows\n")

      # Count liver/metabolic disease hits
      liver_hits <- ot_dt[grepl("liver|hepat|steato|fibrosis|cirrh|NAFLD|NASH|MASLD",
                                 disease_name, ignore.case = TRUE)]
      cat("  Liver/MASLD disease rows:", nrow(liver_hits), "\n")

      fwrite(ot_dt, file.path(RESULTS_DIR, "opentargets_known_drugs.csv"))
      cat("  Saved: opentargets_known_drugs.csv\n")
      ot_dt
    } else {
      cat("  No Open Targets drug annotations found.\n")
      data.table(gene = character(), drug_name = character(),
                 max_clinical_phase = integer(), is_approved = logical())
    }
  } else {
    data.table(gene = character(), drug_name = character(),
               max_clinical_phase = integer(), is_approved = logical())
  }
}, error = function(e) {
  cat("  ERROR in Open Targets query:", conditionMessage(e), "\n")
  data.table(gene = character(), drug_name = character(),
             max_clinical_phase = integer(), is_approved = logical())
})


# ==============================================================================
# Part E: C2:CGP Disease Concordance (retained from v2)
# ==============================================================================
cat("\n--- Part E: C2:CGP Disease Concordance (fgsea) ---\n")

cgp_df <- as.data.table(msigdbr(species = "Homo sapiens",
                                 collection = "C2",
                                 subcollection = "CGP"))
cat("  C2:CGP gene sets loaded:", uniqueN(cgp_df$gs_name), "sets\n")

cgp_df <- cgp_df[ensembl_gene != ""]
pathways <- split(cgp_df$ensembl_gene, cgp_df$gs_name)

cat("  Running fgsea (nPermSimple=10000)...\n")
set.seed(42)
fgsea_res <- fgsea(pathways     = pathways,
                   stats        = ranked_stats,
                   minSize      = FGSEA_MIN_SIZE,
                   maxSize      = FGSEA_MAX_SIZE,
                   nPermSimple  = 10000)
fgsea_dt <- as.data.table(fgsea_res)

n_tested <- nrow(fgsea_dt)
n_valid  <- sum(!is.na(fgsea_dt$padj))
cat("  fgsea completed:", n_tested, "gene sets tested (", n_valid, "valid)\n")
cat("  Significant (padj < 0.05):", sum(fgsea_dt$padj < FGSEA_PADJ_THRESHOLD, na.rm = TRUE), "\n")

# Reversal and concordant hits
reversal_hits   <- fgsea_dt[NES < 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD][order(NES)]
concordant_hits <- fgsea_dt[NES > 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD][order(-NES)]
cat("  Reversal hits (NES < 0, padj < 0.05):", nrow(reversal_hits), "\n")
cat("  Concordant hits (NES > 0, padj < 0.05):", nrow(concordant_hits), "\n")

# Classify CGP sets
classify_cgp <- function(pn) {
  pn <- toupper(pn)
  if (grepl("LIVER|HEPAT|HCC|NAFLD|NASH|STEATO|CIRR", pn)) return("Liver/MASLD")
  if (grepl("DRUG|TREAT|COMPOUND|DOXO|ETOPO|TAMOX|METFORM", pn)) return("Drug_Treatment")
  if (grepl("OBESE|OBESITY|BMI|ADIPOSE|FAT|LIPID|CHOLEST", pn)) return("Metabolic")
  if (grepl("INFLAM|TNF|IL6|NFKB|IMMUNE|MACROPHAGE", pn)) return("Inflammatory")
  if (grepl("FIBROSIS|STELLATE|COLLAGEN|TGF", pn)) return("Fibrosis")
  return("Other_CGP")
}
fgsea_dt[, cgp_class := vapply(pathway, classify_cgp, character(1))]
reversal_hits[, cgp_class := vapply(pathway, classify_cgp, character(1))]

# Save fgsea helper
save_fgsea <- function(dt, file) {
  out <- copy(dt)
  if ("leadingEdge" %in% names(out)) {
    out[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  }
  fwrite(out, file)
  cat("  Saved:", basename(file), "(", nrow(out), "rows)\n")
}

save_fgsea(fgsea_dt[order(padj)], file.path(RESULTS_DIR, "cgp_disease_concordance.csv"))
save_fgsea(reversal_hits, file.path(RESULTS_DIR, "cgp_reversal_hits.csv"))


# ==============================================================================
# Part F: Integration Layer
# ==============================================================================
cat("\n--- Part F: Integration Layer ---\n")

# F1. Build master gene table from upregulated Significant DEGs
master <- data.table(
  ensembl_id = up_with_sym$ensembl_id,
  symbol     = up_with_sym$symbol,
  bulk_sig   = up_with_sym$bulk_sig,
  bulk_logFC = up_with_sym$bulk_logFC,
  bulk_padj  = up_with_sym$bulk_padj
)
master <- master[!is.na(symbol) & symbol != ""]
master <- master[!duplicated(ensembl_id)]
cat("  Master table: ", nrow(master), "upregulated Significant genes\n")

# F2. Annotate with DGIdb
if (nrow(dgidb_results) > 0) {
  dgidb_summary <- dgidb_results[, .(
    dgidb_n_drugs    = uniqueN(drug_name),
    dgidb_n_approved = sum(approved, na.rm = TRUE),
    dgidb_inhibitor  = any(therapeutic, na.rm = TRUE),
    dgidb_drugs      = paste(unique(drug_name)[1:min(5, uniqueN(drug_name))], collapse = "; ")
  ), by = .(gene)]

  master <- merge(master, dgidb_summary, by.x = "symbol", by.y = "gene", all.x = TRUE)
  master[is.na(dgidb_n_drugs), dgidb_n_drugs := 0L]
  master[, dgidb_druggable := dgidb_n_drugs > 0]
  cat("  DGIdb druggable genes:", sum(master$dgidb_druggable), "\n")
} else {
  master[, `:=`(dgidb_n_drugs = 0L, dgidb_druggable = FALSE,
                dgidb_inhibitor = FALSE, dgidb_drugs = NA_character_,
                dgidb_n_approved = 0L)]
}

# F3. Annotate with causal evidence (MR/TWAS)
if (file.exists(MR_FILE)) {
  mr_dt <- fread(MR_FILE)
  mr_annot <- mr_dt[, .(ensembl_id, twas_z, twas_p, twas_fdr, mr_p, mr_fdr, causal_score)]
  mr_annot[, twas_fdr_sig := !is.na(twas_fdr) & twas_fdr < 0.05]
  mr_annot[, twas_nominal := !is.na(twas_p) & twas_p < 0.001]
  mr_annot[, mr_nominal   := !is.na(mr_p) & mr_p < 0.05]

  master <- merge(master, mr_annot, by = "ensembl_id", all.x = TRUE)
  master[is.na(twas_fdr_sig), twas_fdr_sig := FALSE]
  master[is.na(twas_nominal), twas_nominal := FALSE]
  master[is.na(mr_nominal), mr_nominal := FALSE]
  cat("  TWAS FDR sig:", sum(master$twas_fdr_sig), "\n")
  cat("  TWAS nominal:", sum(master$twas_nominal), "\n")
  cat("  MR nominal:  ", sum(master$mr_nominal), "\n")
} else {
  master[, `:=`(twas_fdr_sig = FALSE, twas_nominal = FALSE, mr_nominal = FALSE,
                twas_z = NA_real_, twas_p = NA_real_, twas_fdr = NA_real_,
                mr_p = NA_real_, mr_fdr = NA_real_, causal_score = NA_real_)]
}

# F4. Annotate with Open Targets
if (nrow(opentargets_results) > 0) {
  ot_summary <- opentargets_results[, .(
    ot_has_drug   = TRUE,
    ot_n_drugs    = uniqueN(drug_name),
    ot_max_phase  = max(max_clinical_phase, na.rm = TRUE),
    ot_n_approved = sum(is_approved, na.rm = TRUE),
    ot_drugs      = paste(unique(drug_name)[1:min(5, uniqueN(drug_name))], collapse = "; ")
  ), by = .(gene)]

  master <- merge(master, ot_summary, by.x = "symbol", by.y = "gene", all.x = TRUE)
  master[is.na(ot_has_drug), ot_has_drug := FALSE]
  master[is.na(ot_max_phase), ot_max_phase := 0L]
  cat("  Open Targets annotated:", sum(master$ot_has_drug), "genes\n")
} else {
  master[, `:=`(ot_has_drug = FALSE, ot_n_drugs = 0L, ot_max_phase = 0L,
                ot_n_approved = 0L, ot_drugs = NA_character_)]
}

# F5. Annotate with CGP leading-edge membership (reversal + concordant)
cgp_reversal_le_genes  <- unique(unlist(reversal_hits$leadingEdge))
cgp_concordant_le_genes <- unique(unlist(concordant_hits$leadingEdge))
master[, cgp_reversal_le   := ensembl_id %in% cgp_reversal_le_genes]
master[, cgp_concordant_le := ensembl_id %in% cgp_concordant_le_genes]
cat("  CGP reversal leading-edge:", sum(master$cgp_reversal_le), "genes\n")
cat("  CGP concordant leading-edge:", sum(master$cgp_concordant_le), "genes\n")

# F6. Annotate with LINCS reversal membership
# LINCS results are compound-level; mark genes as "LINCS reversal" if they were
# in the query gene set (top up/down by t-stat) AND LINCS found reversal compounds.
master[, lincs_reversal := FALSE]
if (!is.null(lincs_results) && nrow(lincs_results) > 0) {
  # The query used top 150 up-regulated Entrez IDs; map back to Ensembl
  lincs_query_up_ens <- names(head(ranked_stats, LINCS_N_QUERY))
  master[, lincs_reversal := ensembl_id %in% lincs_query_up_ens]
  cat("  LINCS reversal query genes in master:", sum(master$lincs_reversal), "genes\n")
  cat("  LINCS reversal compounds found:", nrow(lincs_results), "\n")
}

# F7. Compute convergence score (0-10)
master[, convergence_score := 0]
master[dgidb_druggable == TRUE,  convergence_score := convergence_score + 2]
master[dgidb_inhibitor == TRUE,  convergence_score := convergence_score + 1]
master[twas_fdr_sig == TRUE,     convergence_score := convergence_score + 2]
master[twas_nominal == TRUE & twas_fdr_sig == FALSE,
                                 convergence_score := convergence_score + 1]
master[mr_nominal == TRUE,       convergence_score := convergence_score + 1]
master[ot_has_drug == TRUE,      convergence_score := convergence_score + 1]
master[ot_max_phase >= 3,        convergence_score := convergence_score + 1]
master[cgp_concordant_le == TRUE, convergence_score := convergence_score + 1]
master[lincs_reversal == TRUE,   convergence_score := convergence_score + 1]

cat("  Convergence score distribution:\n")
print(master[, .N, by = convergence_score][order(-convergence_score)])

# F8. Match LINCS compounds to DGIdb drugs via cmap_name (not BRD ID)
if (!is.null(lincs_results) && nrow(lincs_results) > 0 && nrow(dgidb_results) > 0) {
  cat("  Matching LINCS compounds to DGIdb drugs...\n")
  if ("cmap_name" %in% names(lincs_results)) {
    lincs_drugs  <- unique(toupper(gsub("[^A-Za-z0-9]", "", lincs_results$cmap_name)))
    lincs_drugs  <- lincs_drugs[lincs_drugs != "" & !is.na(lincs_drugs)]
    dgidb_drugs  <- unique(toupper(gsub("[^A-Za-z0-9]", "", dgidb_results$drug_name)))
    dgidb_drugs  <- dgidb_drugs[dgidb_drugs != "" & !is.na(dgidb_drugs)]
    overlap      <- intersect(lincs_drugs, dgidb_drugs)
    cat("  LINCS drug names (mapped):", length(lincs_drugs), "\n")
    cat("  DGIdb drug names:", length(dgidb_drugs), "\n")
    cat("  LINCS-DGIdb drug name overlap:", length(overlap), "compounds\n")
    if (length(overlap) > 0) cat("  Overlap drugs:", paste(head(overlap, 20), collapse = ", "), "\n")
  } else if ("pert" %in% names(lincs_results)) {
    # Fallback to BRD IDs (won't match DGIdb names)
    cat("  WARNING: cmap_name not available; using BRD IDs (expect 0 overlap)\n")
    lincs_drugs  <- unique(toupper(gsub("[^A-Za-z0-9]", "", lincs_results$pert)))
    dgidb_drugs  <- unique(toupper(gsub("[^A-Za-z0-9]", "", dgidb_results$drug_name)))
    overlap      <- intersect(lincs_drugs, dgidb_drugs)
    cat("  LINCS-DGIdb drug name overlap:", length(overlap), "compounds\n")
  }
}

# Save master summary
setorder(master, -convergence_score)
fwrite(master, file.path(RESULTS_DIR, "pharmacotranscriptomics_summary.csv"))
cat("  Saved: pharmacotranscriptomics_summary.csv (", nrow(master), "genes)\n")

# Save convergent drug targets (convergence_score >= 3)
convergent <- master[convergence_score >= 3]
fwrite(convergent, file.path(RESULTS_DIR, "convergent_drug_targets.csv"))
cat("  Saved: convergent_drug_targets.csv (", nrow(convergent), "genes with score >= 3)\n")


# ==============================================================================
# Part G: MR Convergence
# ==============================================================================
cat("\n--- Part G: MR Convergence ---\n")

if (file.exists(MR_FILE)) {
  mr_dt <- fread(MR_FILE)
  mr_causal <- mr_dt[(!is.na(mr_p) & mr_p < 0.05) | (!is.na(twas_p) & twas_p < 0.05)]
  cat("  MR/TWAS-implicated genes (nominal p < 0.05):", nrow(mr_causal), "\n")

  if (nrow(mr_causal) > 0) {
    mr_ensembl <- mr_causal$ensembl_id
    mr_symbols <- mr_causal$gene

    # Intersect with each evidence layer
    # LINCS: query up/down genes that were used to find reversal compounds
    set_lincs <- if (!is.null(lincs_results) && nrow(lincs_results) > 0) {
      lincs_query_ens <- c(names(head(ranked_stats, LINCS_N_QUERY)),
                           names(tail(ranked_stats, LINCS_N_QUERY)))
      ensembl_to_symbol[ensembl_id %in% lincs_query_ens, symbol]
    } else character(0)
    set_dgidb   <- if (nrow(dgidb_results) > 0) unique(dgidb_results$gene) else character(0)
    set_cgp_le  <- ensembl_to_symbol[ensembl_id %in% c(cgp_reversal_le_genes, cgp_concordant_le_genes), symbol]
    set_ot      <- if (nrow(opentargets_results) > 0) unique(opentargets_results$gene) else character(0)

    conv_lincs  <- intersect(mr_symbols, set_lincs)
    conv_dgidb  <- intersect(mr_symbols, set_dgidb)
    conv_cgp    <- intersect(mr_symbols, set_cgp_le)
    conv_ot     <- intersect(mr_symbols, set_ot)

    all_conv_symbols <- unique(c(conv_lincs, conv_dgidb, conv_cgp, conv_ot))
    cat("  MR-convergent with LINCS:", length(conv_lincs), "\n")
    cat("  MR-convergent with DGIdb:", length(conv_dgidb), "\n")
    cat("  MR-convergent with CGP:  ", length(conv_cgp), "\n")
    cat("  MR-convergent with OT:   ", length(conv_ot), "\n")
    cat("  Total unique MR-convergent:", length(all_conv_symbols), "\n")

    if (length(all_conv_symbols) > 0) {
      conv_dt <- mr_causal[gene %in% all_conv_symbols]
      conv_dt[, lincs_convergent := gene %in% conv_lincs]
      conv_dt[, dgidb_convergent := gene %in% conv_dgidb]
      conv_dt[, cgp_convergent   := gene %in% conv_cgp]
      conv_dt[, ot_convergent    := gene %in% conv_ot]
      conv_dt[, n_drug_layers    := lincs_convergent + dgidb_convergent +
                                     cgp_convergent + ot_convergent]

      setorder(conv_dt, -n_drug_layers)
      fwrite(conv_dt, file.path(RESULTS_DIR, "mr_convergent_drug_targets.csv"))
      cat("  Saved: mr_convergent_drug_targets.csv (", nrow(conv_dt), "rows)\n")
    } else {
      fwrite(data.table(note = "No convergent genes between MR-causal and drug evidence"),
             file.path(RESULTS_DIR, "mr_convergent_drug_targets.csv"))
      cat("  No convergent genes found.\n")
    }
  } else {
    fwrite(data.table(note = "No MR-causal genes at nominal threshold"),
           file.path(RESULTS_DIR, "mr_convergent_drug_targets.csv"))
  }
} else {
  cat("  MR results not found; skipping convergence.\n")
  fwrite(data.table(note = "MR results not available"),
         file.path(RESULTS_DIR, "mr_convergent_drug_targets.csv"))
}


# ==============================================================================
# Part H: Composite Supplementary Figure
# ==============================================================================
cat("\n--- Part H: Generating composite supplementary figure ---\n")

# Panel a: LINCS NCS distribution (or CGP NES if LINCS unavailable)
cat("  Panel a: NCS/NES distribution...\n")
if (!is.null(lincs_results_all) &&
    nrow(lincs_results_all) > 0 && "NCS" %in% names(lincs_results_all)) {
  # Full LINCS results for histogram; highlight FDR<0.05 hits, label top 10 by drug name
  lincs_plot_dt <- copy(lincs_results_all)
  lincs_plot_dt[, fdr_sig := !is.na(WTCS_FDR) & WTCS_FDR < 0.05 & NCS < 0]
  n_fdr <- sum(lincs_plot_dt$fdr_sig, na.rm = TRUE)
  n_disc <- sum(!is.na(lincs_plot_dt$WTCS_FDR) & lincs_plot_dt$WTCS_FDR < 0.25 & lincs_plot_dt$NCS < 0, na.rm = TRUE)

  # Top 10 compounds by NCS for labeling (use cmap_name if available)
  top10 <- head(lincs_plot_dt[fdr_sig == TRUE][order(NCS)], 10)
  if ("cmap_name" %in% names(top10)) {
    top10[, display_name := fifelse(!is.na(cmap_name) & cmap_name != "", cmap_name, pert)]
  } else {
    top10[, display_name := pert]
  }

  pa <- ggplot(lincs_plot_dt, aes(x = NCS)) +
    geom_histogram(bins = 50, fill = "grey70", color = "grey50", linewidth = 0.2) +
    geom_histogram(data = lincs_plot_dt[fdr_sig == TRUE], aes(x = NCS), bins = 50,
                   fill = masld_colors$down, alpha = 0.8, linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = min(lincs_plot_dt$NCS, na.rm = TRUE) * 0.6, y = Inf,
             label = paste0("FDR<0.05: ", n_fdr, "\nFDR<0.25: ", n_disc),
             hjust = 0, vjust = 1.5, color = masld_colors$down, size = 2.2, fontface = "bold") +
    theme_masld() +
    labs(title = "LINCS L1000 signature reversal (HepG2)",
         subtitle = paste0(nrow(lincs_plot_dt), " compound signatures; ", n_fdr, " FDR<0.05 reversals"),
         x = "Normalized Connectivity Score (NCS)", y = "Count")
} else {
  # Fallback: CGP NES distribution
  plot_dt <- fgsea_dt[!is.na(NES)]
  pa <- ggplot(plot_dt, aes(x = NES)) +
    geom_histogram(bins = 60, fill = "grey70", color = "grey50", linewidth = 0.2) +
    geom_histogram(data = plot_dt[NES < 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD],
                   aes(x = NES), bins = 60,
                   fill = masld_colors$down, alpha = 0.7, linewidth = 0.2) +
    geom_histogram(data = plot_dt[NES > 0 & !is.na(padj) & padj < FGSEA_PADJ_THRESHOLD],
                   aes(x = NES), bins = 60,
                   fill = masld_colors$up, alpha = 0.7, linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = min(plot_dt$NES, na.rm = TRUE) * 0.6, y = Inf,
             label = paste0("Reversal: ", nrow(reversal_hits)),
             hjust = 0, vjust = 1.5, color = masld_colors$down, size = 2.2, fontface = "bold") +
    annotate("text", x = max(plot_dt$NES, na.rm = TRUE) * 0.6, y = Inf,
             label = paste0("Concordant: ", nrow(concordant_hits)),
             hjust = 1, vjust = 1.5, color = masld_colors$up, size = 2.2, fontface = "bold") +
    theme_masld() +
    labs(title = "C2:CGP disease signature concordance",
         subtitle = paste0(n_valid, " gene sets; t-statistic ranked"),
         x = "Normalized Enrichment Score (NES)", y = "Count")
}

# Panel b: DGIdb druggability stacked bar
cat("  Panel b: DGIdb druggability...\n")
if (nrow(dgidb_results) > 0) {
  # Summarize interaction types by DEG direction
  dgidb_tier <- merge(
    dgidb_results[, .(gene, interaction_type)],
    master[, .(symbol, bulk_dir = fifelse(bulk_logFC > 0, "Up", "Down"))],
    by.x = "gene", by.y = "symbol", all.x = TRUE
  )
  dgidb_tier[is.na(bulk_dir), bulk_dir := "Unknown"]

  # Clean interaction types
  dgidb_tier[, itype_clean := fifelse(
    grepl("inhibitor", interaction_type, ignore.case = TRUE), "Inhibitor",
    fifelse(grepl("agonist|activator", interaction_type, ignore.case = TRUE), "Activator",
    fifelse(grepl("antagonist|blocker", interaction_type, ignore.case = TRUE), "Antagonist",
    fifelse(grepl("modulator", interaction_type, ignore.case = TRUE), "Modulator",
    fifelse(grepl("antibody", interaction_type, ignore.case = TRUE), "Antibody",
    "Other")))))]

  itype_counts <- dgidb_tier[, .N, by = .(bulk_dir, itype_clean)]

  pb <- ggplot(itype_counts, aes(x = bulk_dir, y = N, fill = itype_clean)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_brewer(palette = "Set2", name = "Interaction") +
    theme_masld() +
    labs(title = "DGIdb drug-gene interactions",
         x = "DEG direction", y = "Drug-gene interactions")
} else {
  pb <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "DGIdb: no results") +
    theme_void()
}

# Panel c: Multi-layer evidence stacked bar (how many layers per gene)
cat("  Panel c: Evidence layer distribution...\n")

# Count evidence layers per gene
master[, n_evidence_layers := 0L]
master[dgidb_druggable == TRUE,               n_evidence_layers := n_evidence_layers + 1L]
master[dgidb_inhibitor == TRUE,               n_evidence_layers := n_evidence_layers + 1L]
master[cgp_concordant_le == TRUE,             n_evidence_layers := n_evidence_layers + 1L]
master[lincs_reversal == TRUE,                n_evidence_layers := n_evidence_layers + 1L]
master[twas_fdr_sig == TRUE | twas_nominal == TRUE, n_evidence_layers := n_evidence_layers + 1L]
master[mr_nominal == TRUE,                    n_evidence_layers := n_evidence_layers + 1L]
master[ot_has_drug == TRUE,                   n_evidence_layers := n_evidence_layers + 1L]

layer_dist <- master[, .N, by = n_evidence_layers][order(n_evidence_layers)]
layer_dist[, label := fifelse(n_evidence_layers >= 4,
                               paste0(n_evidence_layers, "+"),
                               as.character(n_evidence_layers))]

pc <- ggplot(layer_dist, aes(x = factor(n_evidence_layers), y = N)) +
  geom_col(fill = masld_colors$deg, width = 0.7) +
  geom_text(aes(label = N), vjust = -0.3, size = 2.2) +
  theme_masld() +
  labs(title = "Evidence layers per gene",
       subtitle = paste0(nrow(master), " upregulated Tier 1+2 genes"),
       x = "Number of evidence layers", y = "Genes")

# Also save per-layer coverage as a separate horizontal bar
layer_coverage <- data.table(
  layer   = c("DGIdb druggable", "DGIdb inhibitor", "CGP concordant LE",
              "LINCS query", "TWAS/MR causal", "Open Targets"),
  n_genes = c(sum(master$dgidb_druggable),
              sum(master$dgidb_inhibitor, na.rm = TRUE),
              sum(master$cgp_concordant_le),
              sum(master$lincs_reversal),
              sum(master$twas_fdr_sig | master$twas_nominal | master$mr_nominal),
              sum(master$ot_has_drug))
)
layer_coverage <- layer_coverage[n_genes > 0]  # Only show non-empty layers
cat("  Per-layer coverage:\n")
print(layer_coverage)

# Panel d: Top convergent drug-target dot plot (DGIdb + LINCS sources)
cat("  Panel d: Convergent drug-target dot plot...\n")
{
  dot_rows <- list()

  # DGIdb drug-target pairs
  if (nrow(dgidb_results) > 0) {
    top_conv <- master[dgidb_druggable == TRUE & convergence_score >= 2][
      order(-convergence_score)][1:min(15, .N)]
    if (nrow(top_conv) > 0) {
      top_drug_per_gene <- dgidb_results[gene %in% top_conv$symbol,
                                          .SD[which.max(interaction_score)],
                                          by = gene]
      dgidb_dot <- merge(
        top_conv[, .(symbol, convergence_score)],
        top_drug_per_gene[, .(gene, drug_name, approved)],
        by.x = "symbol", by.y = "gene", all.x = TRUE
      )
      dgidb_dot <- dgidb_dot[!is.na(drug_name)]
      if (nrow(dgidb_dot) > 0) {
        dgidb_dot[, source := "DGIdb"]
        dot_rows[[length(dot_rows) + 1L]] <- dgidb_dot
      }
    }
  }

  # LINCS top reversal drugs (if cmap_name available)
  if (!is.null(lincs_results) && nrow(lincs_results) > 0 && "cmap_name" %in% names(lincs_results)) {
    lincs_top <- lincs_results[!is.na(cmap_name) & cmap_name != ""][order(NCS)][1:min(5, .N)]
    if (nrow(lincs_top) > 0) {
      # These are compound-level, not gene-level; attach to top convergent gene
      lincs_dot <- data.table(
        symbol = "MASLD signature",
        convergence_score = NA_real_,
        tier = NA_character_,
        drug_name = lincs_top$cmap_name,
        approved = NA,
        source = "LINCS"
      )
      dot_rows[[length(dot_rows) + 1L]] <- lincs_dot
    }
  }

  if (length(dot_rows) > 0) {
    dot_dt <- rbindlist(dot_rows, fill = TRUE)
    dot_dt[, label := paste0(symbol, " — ", drug_name)]
    dot_dt <- dot_dt[!duplicated(label)]
    dot_dt[, label := factor(label, levels = rev(label))]
    dot_dt[is.na(convergence_score), convergence_score := 0]

    pd <- ggplot(dot_dt, aes(x = convergence_score, y = label,
                               color = source)) +
      geom_point(size = 3) +
      scale_color_manual(values = c("DGIdb" = masld_colors$deg,
                                     "LINCS" = masld_colors$twas),
                         name = "Source") +
      theme_masld() +
      labs(title = "Top drug-target pairs",
           x = "Convergence score", y = NULL)
  } else {
    pd <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No drug-target pairs") +
      theme_void()
  }
}

# Panel e: Open Targets clinical phase bar chart
cat("  Panel e: Open Targets clinical phases...\n")
if (nrow(opentargets_results) > 0) {
  phase_counts <- opentargets_results[, .(n_drugs = uniqueN(drug_name)),
                                        by = .(max_clinical_phase)]
  phase_counts[, phase_label := fifelse(max_clinical_phase == 4, "Approved",
                                  fifelse(max_clinical_phase == 3, "Phase III",
                                  fifelse(max_clinical_phase == 2, "Phase II",
                                  fifelse(max_clinical_phase == 1, "Phase I",
                                  "Preclinical"))))]
  phase_counts[, phase_label := factor(phase_label,
                                        levels = c("Preclinical", "Phase I", "Phase II",
                                                   "Phase III", "Approved"))]

  pe <- ggplot(phase_counts, aes(x = phase_label, y = n_drugs)) +
    geom_col(fill = masld_colors$twas, width = 0.7) +
    theme_masld() +
    labs(title = "Open Targets: clinical phases",
         subtitle = "For MR/TWAS-implicated targets",
         x = "Max clinical phase", y = "Unique drugs") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  # Summary statistics table panel as fallback
  n_lincs <- if (!is.null(lincs_results)) nrow(lincs_results) else 0L
  n_lincs_fdr <- if (!is.null(sig_reversals)) nrow(sig_reversals) else 0L
  n_dgidb_genes <- sum(master$dgidb_druggable, na.rm = TRUE)
  n_dgidb_drugs <- if (nrow(dgidb_results) > 0) uniqueN(dgidb_results$drug_name) else 0L
  n_cgp <- sum(master$cgp_concordant_le, na.rm = TRUE)
  n_conv3 <- sum(master$convergence_score >= 3, na.rm = TRUE)

  summary_dt <- data.table(
    Metric = c("LINCS reversals (FDR<0.05)",
               "LINCS reversals (FDR<0.25)",
               "DGIdb druggable genes",
               "DGIdb unique drugs",
               "CGP concordant LE genes",
               "Convergence score >= 3"),
    Value  = c(n_lincs_fdr, n_lincs, n_dgidb_genes, n_dgidb_drugs, n_cgp, n_conv3)
  )
  summary_dt[, Metric := factor(Metric, levels = rev(Metric))]

  pe <- ggplot(summary_dt, aes(x = Value, y = Metric)) +
    geom_col(fill = masld_colors$twas, width = 0.6) +
    geom_text(aes(label = Value), hjust = -0.2, size = 2.5) +
    theme_masld() +
    labs(title = "Pharmacotranscriptomics summary",
         x = "Count", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3)))
}

# Composite figure using patchwork
cat("  Compositing figure...\n")
composite <- (pa | pb) / (pc | pd | pe) +
  plot_annotation(
    title = "Supplementary Figure: Pharmacotranscriptomics",
    subtitle = "Comprehensive drug repurposing analysis for MASLD consensus DEGs",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7, color = "grey40")
    )
  ) +
  plot_layout(heights = c(1, 1))

fig_path <- file.path(FIG_DIR, "figS_pharmacotranscriptomics.pdf")
save_fig(composite, fig_path, width = fig_full_width, height = 10)
cat("  Saved:", fig_path, "\n")
cat("  Saved:", sub("\\.pdf$", ".png", fig_path), "\n")


# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Strategy 11 v3: Comprehensive Pharmacotranscriptomics Complete ===\n")
cat("  Arm 1 (LINCS):     ", ifelse(!is.null(lincs_results),
                                      paste(nrow(lincs_results), "reversal compounds (WTCS_FDR<0.25)"),
                                      "SKIPPED"), "\n")
if (!is.null(sig_reversals)) {
  cat("  LINCS FDR<0.05:   ", nrow(sig_reversals), "compounds\n")
}
cat("  Arm 2a (DGIdb):    ", nrow(dgidb_results), "drug-gene pairs\n")
cat("  Arm 2b (OT):       ", nrow(opentargets_results), "drug-target-disease rows\n")
cat("  C2:CGP reversal:   ", nrow(reversal_hits), "hits\n")
cat("  C2:CGP concordant: ", nrow(concordant_hits), "hits\n")
cat("  CGP concordant LE: ", sum(master$cgp_concordant_le), "genes in master\n")
cat("  Master table:      ", nrow(master), "genes scored\n")
cat("  Convergent (>=3):  ", sum(master$convergence_score >= 3), "genes\n")
cat("  Output dir:        ", RESULTS_DIR, "\n")
cat("  Figures dir:       ", FIG_DIR, "\n")
cat("End time:", format(Sys.time()), "\n")
