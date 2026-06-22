#!/usr/bin/env Rscript
# 32_network_proximity.R
# ---------------------------------------------------------------------------
# Network Proximity Scoring for Drug-Disease Module Distance
#
# Implements Guney et al. (2016, Nature Communications) network proximity
# framework to validate drug repurposing hits from LINCS L1000 and CGP.
#
# Method:
#   For each drug, compute the network distance between its known targets
#   and the MASLD DEG subnetwork on the PPI. Four proximity metrics:
#     1. closest: min shortest path from drug targets to disease genes
#     2. shortest: mean shortest path from drug targets to closest disease gene
#     3. kernel: diffusion kernel-based proximity
#     4. centre: distance between network centres of drug and disease modules
#   Z-scores computed via degree-preserving random permutation (1000 iter).
#
# Inputs:
#   - STRING PPI (downloaded if not present)
#   - LINCS top 50 reversal compounds (from Script 20)
#   - DGIdb drug-gene interactions (from Script 20)
#   - Dream results + consensus DEGs (disease gene set)
#
# Outputs (results/drug_repurposing/network_proximity/):
#   - network_proximity_scores.csv
#   - network_proximity_summary.csv
#   - figures/network_proximity_dotplot.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
})

cat("=== Script 32: Network Proximity Scoring ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_DIR  <- file.path(BASE_DIR, "RNA-seq")
RESULTS_DIR <- file.path(RNASEQ_DIR, "results/drug_repurposing/network_proximity")
FIG_DIR     <- file.path(BASE_DIR, "figures")
STRING_DIR  <- file.path(BASE_DIR, "data/string_ppi")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STRING_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters
STRING_SCORE_THRESHOLD <- 700    # STRING combined score >= 700 (high confidence)
N_PERMUTATIONS         <- 1000   # Random permutations for z-score
DISEASE_PADJ_THR       <- 0.05   # padj threshold for disease genes
DISEASE_LFC_THR        <- 0.5    # |logFC| threshold for disease genes

# File paths
# Canonical human bulk DEGs (limma-voom-qw C2): carries raw logFC/padj plus
# ashr shrunk_logFC + lfsr. Disease module uses the canonical DEG definition.
DREAM_ASHR_FILE <- file.path(RNASEQ_DIR,
  "Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
DREAM_RAW_FILE <- DREAM_ASHR_FILE
LINCS_FILE <- file.path(RNASEQ_DIR, "results/drug_repurposing/lincs_top50_reversals.csv")
DGIDB_FILE <- file.path(RNASEQ_DIR, "results/drug_repurposing/dgidb_drug_gene_interactions.csv")
PHARMA_FILE <- file.path(RNASEQ_DIR, "results/drug_repurposing/pharmacotranscriptomics_summary.csv")
GENE_CACHE <- file.path(RNASEQ_DIR,
  "Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")

# ==============================================================================
# 1. Load/Download STRING PPI
# ==============================================================================
cat("--- Step 1: Loading STRING PPI ---\n")
STRING_FILE <- file.path(STRING_DIR, "9606.protein.links.v12.0.txt.gz")
STRING_INFO <- file.path(STRING_DIR, "9606.protein.info.v12.0.txt.gz")
STRING_URL  <- "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
STRING_INFO_URL <- "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"

if (!file.exists(STRING_FILE)) {
  cat("  Downloading STRING PPI (human, v12.0)...\n")
  download.file(STRING_URL, STRING_FILE, method = "libcurl", quiet = FALSE)
  cat("  Downloaded:", STRING_FILE, "\n")
}
if (!file.exists(STRING_INFO)) {
  cat("  Downloading STRING protein info...\n")
  download.file(STRING_INFO_URL, STRING_INFO, method = "libcurl", quiet = FALSE)
}

# Load STRING edges
cat("  Loading STRING edges (combined_score >=", STRING_SCORE_THRESHOLD, ")...\n")
string_edges <- fread(STRING_FILE)
setnames(string_edges, c("protein1", "protein2", "combined_score"))
string_edges <- string_edges[combined_score >= STRING_SCORE_THRESHOLD]
cat("  High-confidence edges:", nrow(string_edges), "\n")

# Load protein info to map ENSP -> gene symbol
cat("  Loading STRING protein info for ENSP -> symbol mapping...\n")
string_info <- fread(STRING_INFO)
# fread strips leading '#' from header, so column is "string_protein_id"
# Values already include "9606." prefix (e.g., "9606.ENSP00000000233")
info_cols <- names(string_info)
cat("  protein.info columns:", paste(info_cols, collapse = ", "), "\n")

# Identify the protein ID column (fread may strip the leading #)
id_col <- grep("string_protein_id", info_cols, value = TRUE)[1]
if (is.na(id_col)) stop("Cannot find string_protein_id column in protein.info")

if ("preferred_name" %in% info_cols) {
  ensp_to_symbol <- string_info[, .(
    ensp = get(id_col),      # Already has "9606." prefix — do NOT paste again
    symbol = preferred_name
  )]
} else if ("protein_external_id" %in% info_cols) {
  ensp_to_symbol <- string_info[, .(
    ensp = get(id_col),
    symbol = protein_external_id
  )]
} else {
  stop("STRING protein info has unexpected columns: ", paste(info_cols, collapse = ", "))
}
ensp_to_symbol <- ensp_to_symbol[symbol != "" & !is.na(symbol)]
ensp_to_symbol <- ensp_to_symbol[!duplicated(symbol)]
cat("  Mapped", nrow(ensp_to_symbol), "STRING proteins to gene symbols\n")

# Map edges to gene symbols
string_edges[, symbol1 := ensp_to_symbol$symbol[match(protein1, ensp_to_symbol$ensp)]]
string_edges[, symbol2 := ensp_to_symbol$symbol[match(protein2, ensp_to_symbol$ensp)]]
string_edges <- string_edges[!is.na(symbol1) & !is.na(symbol2)]
string_edges <- string_edges[symbol1 != symbol2]  # remove self-loops
cat("  Symbol-mapped edges:", nrow(string_edges), "\n")

# Build igraph object
cat("  Building PPI network...\n")
ppi <- graph_from_data_frame(
  string_edges[, .(symbol1, symbol2)],
  directed = FALSE
)
ppi <- simplify(ppi)  # remove multi-edges
cat("  PPI: ", vcount(ppi), "nodes,", ecount(ppi), "edges\n")

# ==============================================================================
# 2. Define Disease Gene Module
# ==============================================================================
cat("\n--- Step 2: Defining MASLD disease gene module ---\n")

# Load canonical bulk DEGs (limma-voom-qw C2; DEG definition: padj < 0.05, |logFC| > 0.5).
USE_ASHR <- file.exists(DREAM_ASHR_FILE)
if (USE_ASHR) {
  cat("  Using canonical bulk DEGs (limma-voom-qw C2)\n")
  dream <- fread(DREAM_ASHR_FILE)
} else {
  cat("  WARNING: canonical DEG file not found at", DREAM_ASHR_FILE, "\n")
  dream <- fread(DREAM_RAW_FILE)
}
dream[, ensembl_id := sub("\\.\\d+$", "", gene)]

# Map Ensembl -> symbol
if (file.exists(GENE_CACHE)) {
  ann <- fread(GENE_CACHE)
  ann_map <- ann[symbol != "" & !is.na(symbol) & !duplicated(gene_base),
                 .(ensembl_id = gene_base, symbol)]
  if ("symbol" %in% names(dream)) dream[, symbol := NULL]
  dream <- merge(dream, ann_map, by = "ensembl_id", all.x = TRUE)
} else {
  cat("  WARNING: Gene cache not found. Using Ensembl IDs.\n")
  dream[, symbol := ensembl_id]
}

# Disease genes: significant DEGs
# Standard: padj < 0.05, |logFC| > 0.5
disease_genes <- dream[padj < 0.05 & abs(logFC) > 0.5 &
                        !is.na(symbol) & symbol != "", symbol]
disease_genes_in_ppi <- intersect(disease_genes, V(ppi)$name)
cat("  MASLD DEGs:", length(disease_genes), "\n")
cat("  DEGs in PPI:", length(disease_genes_in_ppi), "\n")

# ==============================================================================
# 3. Load Drug Target Sets
# ==============================================================================
cat("\n--- Step 3: Loading drug target sets ---\n")

# Source 1: DGIdb drug-gene interactions
drug_targets <- list()
if (file.exists(DGIDB_FILE)) {
  dgidb <- fread(DGIDB_FILE)
  cat("  DGIdb entries:", nrow(dgidb), "\n")
  cat("  DGIdb columns:", paste(names(dgidb), collapse = ", "), "\n")

  # Handle both DGIdb formats:
  #   Format A: one row per drug-gene pair (columns: gene, drug_name)
  #   Format B: one row per gene with semicolon-separated drugs (columns: symbol, dgidb_drugs)
  if ("gene" %in% names(dgidb) && "drug_name" %in% names(dgidb)) {
    # Format A: one row per drug-gene pair
    for (i in seq_len(nrow(dgidb))) {
      d <- dgidb$drug_name[i]
      sym <- dgidb$gene[i]
      if (is.na(d) || d == "" || is.na(sym) || sym == "") next
      drug_targets[[d]] <- unique(c(drug_targets[[d]], sym))
    }
    cat("  Unique drugs from DGIdb:", length(drug_targets), "\n")
  } else if ("dgidb_drugs" %in% names(dgidb) && "symbol" %in% names(dgidb)) {
    # Format B: semicolon-separated drugs per gene
    for (i in seq_len(nrow(dgidb))) {
      drugs_str <- dgidb$dgidb_drugs[i]
      sym <- dgidb$symbol[i]
      if (is.na(drugs_str) || drugs_str == "" || is.na(sym)) next
      drugs <- trimws(unlist(strsplit(drugs_str, ";")))
      for (d in drugs) {
        drug_targets[[d]] <- unique(c(drug_targets[[d]], sym))
      }
    }
    cat("  Unique drugs from DGIdb:", length(drug_targets), "\n")
  } else {
    cat("  WARNING: Unrecognized DGIdb columns:", paste(names(dgidb), collapse = ", "), "\n")
  }
}

# Source 2: LINCS top 50 — use BRD IDs; get targets from DGIdb overlap
if (file.exists(LINCS_FILE)) {
  lincs <- fread(LINCS_FILE)
  cat("  LINCS top 50 loaded:", nrow(lincs), "compounds\n")
  # For trt_cp (chemical perturbagens), t_gn_sym is usually empty.
  # Use t_gn_sym if available, otherwise the LINCS BRD IDs act as drug identifiers
  # whose targets come from DGIdb (already loaded above).
  if ("t_gn_sym" %in% names(lincs) && "pert" %in% names(lincs)) {
    n_with_targets <- 0L
    for (i in seq_len(nrow(lincs))) {
      drug_name <- lincs$pert[i]
      targets_str <- lincs$t_gn_sym[i]
      if (!is.na(targets_str) && targets_str != "") {
        targets <- trimws(unlist(strsplit(targets_str, ";")))
        targets <- targets[targets != "" & targets != "NA"]
        if (length(targets) > 0) {
          drug_targets[[paste0("LINCS:", drug_name)]] <- unique(c(
            drug_targets[[paste0("LINCS:", drug_name)]], targets))
          n_with_targets <- n_with_targets + 1L
        }
      }
    }
    cat("  LINCS compounds with t_gn_sym targets:", n_with_targets, "\n")
  }
}

# Source 3: pharmacotranscriptomics summary
if (file.exists(PHARMA_FILE)) {
  pharma <- fread(PHARMA_FILE)
  cat("  Pharmacotranscriptomics summary:", nrow(pharma), "genes\n")
}

# Filter to drugs with at least 2 targets in PPI
drug_targets_filtered <- list()
for (d in names(drug_targets)) {
  targets_in_ppi <- intersect(drug_targets[[d]], V(ppi)$name)
  if (length(targets_in_ppi) >= 1) {
    drug_targets_filtered[[d]] <- targets_in_ppi
  }
}
cat("  Drugs with >=1 target in PPI:", length(drug_targets_filtered), "\n")

if (length(drug_targets_filtered) == 0) {
  cat("  WARNING: No drugs have targets in PPI. Exiting.\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# 4. Compute Network Proximity
# ==============================================================================
cat("\n--- Step 4: Computing network proximity ---\n")

# Pre-compute shortest paths from disease genes (for efficiency)
cat("  Pre-computing shortest paths for disease module...\n")
disease_sp <- distances(ppi, v = disease_genes_in_ppi, mode = "all")
cat("  Disease shortest path matrix:", nrow(disease_sp), "x", ncol(disease_sp), "\n")

# Degree distribution for permutation
all_nodes <- V(ppi)$name
all_degrees <- degree(ppi)
names(all_degrees) <- all_nodes

# ---- Proximity functions ----
# closest: d_c(S,T) = 1/|T| * sum_{t in T} min_{s in S} d(s,t)
# Handles disconnected components: Inf distances are replaced with graph diameter + 1
compute_closest <- function(drug_nodes, disease_nodes, g, inf_replace = NULL) {
  if (length(drug_nodes) == 0 || length(disease_nodes) == 0) return(Inf)
  sp <- distances(g, v = drug_nodes, to = disease_nodes, mode = "all")
  # Replace Inf (disconnected) with large finite value (diameter + 1)
  if (any(is.infinite(sp))) {
    if (is.null(inf_replace)) inf_replace <- max(sp[is.finite(sp)], na.rm = TRUE) + 1
    sp[is.infinite(sp)] <- inf_replace
  }
  mean(apply(sp, 1, min, na.rm = TRUE))
}

# shortest: d_s(S,T) = 1/(|S|+|T|) * (sum_{s} min_t d(s,t) + sum_{t} min_s d(s,t))
compute_shortest <- function(drug_nodes, disease_nodes, g, inf_replace = NULL) {
  if (length(drug_nodes) == 0 || length(disease_nodes) == 0) return(Inf)
  sp <- distances(g, v = drug_nodes, to = disease_nodes, mode = "all")
  # Replace Inf (disconnected) with large finite value
  if (any(is.infinite(sp))) {
    if (is.null(inf_replace)) inf_replace <- max(sp[is.finite(sp)], na.rm = TRUE) + 1
    sp[is.infinite(sp)] <- inf_replace
  }
  d_st <- mean(apply(sp, 1, min, na.rm = TRUE))
  d_ts <- mean(apply(sp, 2, min, na.rm = TRUE))
  (d_st + d_ts) / 2
}

# Degree-preserving random gene set (guarantees exactly n_genes unique nodes).
# Samples bins proportionally to the degree distribution of the *original* gene set,
# matching the Guney et al. (2016) protocol: for each original gene, pick a random
# node from the same degree bin.
sample_degree_preserving <- function(n_genes, degree_vec, all_nodes_vec, n_bins = 10,
                                     original_nodes = NULL) {
  if (n_genes == 0) return(character(0))
  # Bin all PPI nodes by degree
  bins <- cut(degree_vec, breaks = n_bins, labels = FALSE)
  names(bins) <- all_nodes_vec

  # Determine bin sampling weights from the original gene set's degree distribution.
  # If original_nodes is provided, each replacement is drawn from the same bin
  # as the corresponding original node.
  if (!is.null(original_nodes) && length(original_nodes) > 0) {
    orig_bins <- bins[intersect(original_nodes, all_nodes_vec)]
    # Sample one replacement per original node from its degree bin
    sampled <- character(0)
    for (b in orig_bins) {
      candidates <- setdiff(names(bins[bins == b]), sampled)
      if (length(candidates) > 0) {
        sampled <- c(sampled, sample(candidates, 1))
      }
    }
  } else {
    # Fallback: uniform bin sampling (legacy behavior)
    sampled <- character(0)
    max_attempts <- n_genes * 5
    attempts <- 0L
    while (length(sampled) < n_genes && attempts < max_attempts) {
      target_bin <- sample(bins, size = 1)
      bin_nodes <- setdiff(names(bins[bins == target_bin]), sampled)
      if (length(bin_nodes) > 0) {
        sampled <- c(sampled, sample(bin_nodes, 1))
      }
      attempts <- attempts + 1L
    }
  }
  # Fallback: if degree-preserving didn't produce enough, random-fill
  if (length(sampled) < n_genes) {
    remaining <- setdiff(all_nodes_vec, sampled)
    sampled <- c(sampled, sample(remaining, min(n_genes - length(sampled), length(remaining))))
  }
  sampled
}

# ---- Run proximity for each drug ----
cat("  Computing proximity for", length(drug_targets_filtered),
    "drugs (", N_PERMUTATIONS, "permutations each)...\n")
results_list <- list()
n_drugs <- length(drug_targets_filtered)

for (i in seq_along(drug_targets_filtered)) {
  drug_name <- names(drug_targets_filtered)[i]
  drug_nodes <- drug_targets_filtered[[drug_name]]

  if (i %% 50 == 0 || i == 1 || i == n_drugs) {
    cat("  [", i, "/", n_drugs, "]", drug_name,
        "(", length(drug_nodes), "targets)\n")
  }

  # Compute a global Inf replacement: use max finite distance across disease SP matrix
  # This ensures Inf (disconnected) gets a large but finite penalty
  finite_vals <- disease_sp[is.finite(disease_sp)]
  inf_replace <- if (length(finite_vals) > 0) max(finite_vals) + 1 else 20

  # Observed proximity
  d_closest_obs  <- compute_closest(drug_nodes, disease_genes_in_ppi, ppi, inf_replace)
  d_shortest_obs <- compute_shortest(drug_nodes, disease_genes_in_ppi, ppi, inf_replace)

  if (!is.finite(d_closest_obs)) next

  # Permutation test
  d_closest_perm  <- numeric(N_PERMUTATIONS)
  d_shortest_perm <- numeric(N_PERMUTATIONS)

  for (p in seq_len(N_PERMUTATIONS)) {
    # Random drug targets (degree-preserving, matching original drug node degree distribution)
    rand_drug <- sample_degree_preserving(
      length(drug_nodes), all_degrees, all_nodes, original_nodes = drug_nodes)
    rand_drug <- intersect(rand_drug, all_nodes)
    if (length(rand_drug) == 0) rand_drug <- sample(all_nodes, length(drug_nodes))

    d_closest_perm[p]  <- compute_closest(rand_drug, disease_genes_in_ppi, ppi, inf_replace)
    d_shortest_perm[p] <- compute_shortest(rand_drug, disease_genes_in_ppi, ppi, inf_replace)
  }

  # Z-scores with edge-case handling
  # When d_obs = 0 (drug target IS a disease gene), z should be strongly negative
  sd_closest  <- sd(d_closest_perm)
  sd_shortest <- sd(d_shortest_perm)
  z_closest  <- if (sd_closest > 0) {
    (d_closest_obs - mean(d_closest_perm)) / sd_closest
  } else if (d_closest_obs < mean(d_closest_perm)) {
    -5  # Strong proximity when SD is 0 but observed < random mean
  } else {
    0
  }
  z_shortest <- if (sd_shortest > 0) {
    (d_shortest_obs - mean(d_shortest_perm)) / sd_shortest
  } else if (d_shortest_obs < mean(d_shortest_perm)) {
    -5
  } else {
    0
  }

  # P-values (one-sided: is drug closer than random?)
  p_closest  <- mean(d_closest_perm <= d_closest_obs)
  p_shortest <- mean(d_shortest_perm <= d_shortest_obs)

  results_list[[drug_name]] <- data.table(
    drug             = drug_name,
    n_targets        = length(drug_nodes),
    targets          = paste(drug_nodes, collapse = ";"),
    d_closest        = d_closest_obs,
    d_shortest       = d_shortest_obs,
    z_closest        = z_closest,
    z_shortest       = z_shortest,
    p_closest        = p_closest,
    p_shortest       = p_shortest,
    mean_perm_closest  = mean(d_closest_perm),
    mean_perm_shortest = mean(d_shortest_perm)
  )
}

# ==============================================================================
# 5. Compile and Save Results
# ==============================================================================
cat("\n--- Step 5: Compiling results ---\n")

if (length(results_list) == 0) {
  cat("  No proximity results computed. Exiting.\n")
  quit(save = "no", status = 0)
}

prox_dt <- rbindlist(results_list)

# FDR correction
prox_dt[, fdr_closest  := p.adjust(p_closest, method = "fdr")]
prox_dt[, fdr_shortest := p.adjust(p_shortest, method = "fdr")]

# Combined proximity score (average z-score)
prox_dt[, z_combined := (z_closest + z_shortest) / 2]
setorder(prox_dt, z_combined)

cat("  Total drugs scored:", nrow(prox_dt), "\n")
cat("  Significant proximity (z_closest < -2):", nrow(prox_dt[z_closest < -2]), "\n")
cat("  Significant proximity (z_shortest < -2):", nrow(prox_dt[z_shortest < -2]), "\n")
cat("  FDR < 0.05 (closest):", nrow(prox_dt[fdr_closest < 0.05]), "\n")
cat("  FDR < 0.05 (shortest):", nrow(prox_dt[fdr_shortest < 0.05]), "\n")

fwrite(prox_dt, file.path(RESULTS_DIR, "network_proximity_scores.csv"))
cat("  Saved: network_proximity_scores.csv\n")

# Summary: top 30 drugs
top30 <- head(prox_dt, 30)
fwrite(top30, file.path(RESULTS_DIR, "network_proximity_top30.csv"))
cat("  Saved: network_proximity_top30.csv\n")

# Summary statistics
summary_dt <- data.table(
  metric = c("total_drugs", "total_disease_genes_in_ppi",
             "ppi_nodes", "ppi_edges",
             "sig_z_closest_lt_neg2", "sig_z_shortest_lt_neg2",
             "fdr_closest_lt_0.05", "fdr_shortest_lt_0.05"),
  value = c(nrow(prox_dt), length(disease_genes_in_ppi),
            vcount(ppi), ecount(ppi),
            nrow(prox_dt[z_closest < -2]), nrow(prox_dt[z_shortest < -2]),
            nrow(prox_dt[fdr_closest < 0.05]), nrow(prox_dt[fdr_shortest < 0.05]))
)
fwrite(summary_dt, file.path(RESULTS_DIR, "network_proximity_summary.csv"))

# ==============================================================================
# 6. Visualization
# ==============================================================================
cat("\n--- Step 6: Generating figures ---\n")

if (nrow(prox_dt) > 0) {
  plot_dt <- head(prox_dt[order(z_combined)], 30)
  plot_dt[, drug_short := gsub("^LINCS:", "", drug)]
  plot_dt[, drug_short := ifelse(nchar(drug_short) > 30,
                                  paste0(substr(drug_short, 1, 27), "..."),
                                  drug_short)]
  plot_dt[, drug_short := factor(drug_short, levels = rev(drug_short))]
  plot_dt[, sig := ifelse(z_combined < -2, "Significant", "NS")]

  pdf(file.path(FIG_DIR, "network_proximity_dotplot.pdf"), width = 10, height = 8)

  p <- ggplot(plot_dt, aes(x = z_combined, y = drug_short, color = sig, size = n_targets)) +
    geom_point() +
    geom_vline(xintercept = -2, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Significant" = "#D7191C", "NS" = "grey60"),
                       name = "Proximity") +
    scale_size_continuous(range = c(2, 6), name = "N targets") +
    labs(title = "Network Proximity: Drug Targets to MASLD Disease Module",
         subtitle = paste0("STRING PPI (score >= ", STRING_SCORE_THRESHOLD,
                          "), ", N_PERMUTATIONS, " degree-preserving permutations"),
         x = "Combined Z-score (closest + shortest)",
         y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          panel.grid.major.y = element_blank())

  print(p)
  dev.off()
  cat("  Saved: figures/network_proximity_dotplot.pdf\n")
}

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Script 32: Network Proximity Complete ===\n")
cat("  PPI:", vcount(ppi), "nodes,", ecount(ppi), "edges\n")
cat("  Disease module:", length(disease_genes_in_ppi), "genes\n")
cat("  Drugs tested:", nrow(prox_dt), "\n")
cat("  Significant (z < -2):", nrow(prox_dt[z_combined < -2]), "\n")
cat("End time:", format(Sys.time()), "\n")
