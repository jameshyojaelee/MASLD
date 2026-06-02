#!/usr/bin/env Rscript
# 103b_pathway_features_decoupler.R
# ---------------------------------------------------------------------------
# Size-corrected per-sample pathway scores via decoupleR ULM/MLM.
#
# Addresses reviewer concern R1 §6.B: ssGSEA scores correlate with gene set
# size — larger sets produce higher absolute scores.  decoupleR's ULM and MLM
# fit a linear model of gene-set membership on expression and produce
# t-statistics that are NOT biased by set size.
#
# Input:  results/integration/merged_dge.rds  (39,276 genes × 847 samples)
# Output (to results/integration/):
#   pathway_scores_ulm.rds                — sample × pathway matrix (ULM t-stat)
#   pathway_scores_mlm.rds                — sample × pathway matrix (MLM t-stat)
#   pathway_scores_ssgsea.rds             — sample × pathway matrix (legacy)
#   pathway_size_correction_diagnostic.csv — size-independence evidence
#
# Usage:  Rscript 103b_pathway_features_decoupler.R
# SLURM:  cpu, 8 CPUs, 64G, 48h, --qos=nslab, --job-name=ssGSEA-ULM
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(msigdbr)
  library(decoupleR)
  library(GSVA)
  library(BiocParallel)
})

set.seed(42)
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
ME     <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUTDIR <- RDIR
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 103b: Pathway Features — decoupleR ULM/MLM + ssGSEA ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Cores:", NCORES, "\n\n")

# ============================================================
# STEP 1: Load expression data and compute TMM logCPM
# ============================================================
cat("--- Step 1: Load DGE and compute logCPM ---\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("DGE loaded:", nrow(dge), "genes x", ncol(dge), "samples\n")

dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("logCPM matrix:", nrow(logcpm), "x", ncol(logcpm), "\n\n")

# ============================================================
# STEP 2: Map Ensembl IDs to gene symbols
# ============================================================
cat("--- Step 2: Map Ensembl IDs to gene symbols ---\n")

# Load symbol mapping from multi-evidence atlas
atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
atlas <- atlas[!is.na(human_symbol) & human_symbol != ""]
atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
atlas <- atlas[!duplicated(ensembl_clean)]

# Clean rownames of logCPM
ensembl_clean <- sub("\\..*", "", rownames(logcpm))
idx <- match(ensembl_clean, atlas$ensembl_clean)
mapped <- !is.na(idx)

cat("Genes mapped to symbols:", sum(mapped), "of", nrow(logcpm),
    sprintf("(%.1f%%)\n", 100 * sum(mapped) / nrow(logcpm)))

# Subset to mapped genes
logcpm_sym <- logcpm[mapped, ]
symbols <- atlas$human_symbol[idx[mapped]]

# Handle duplicate symbols: keep the one with highest mean expression
mean_expr <- rowMeans(logcpm_sym)
dup_syms <- symbols[duplicated(symbols)]
if (length(dup_syms) > 0) {
  cat("Resolving", length(unique(dup_syms)), "duplicate symbols",
      "(keeping highest mean expression)\n")
  keep <- rep(TRUE, length(symbols))
  for (s in unique(dup_syms)) {
    which_dup <- which(symbols == s)
    best <- which_dup[which.max(mean_expr[which_dup])]
    keep[setdiff(which_dup, best)] <- FALSE
  }
  logcpm_sym <- logcpm_sym[keep, ]
  symbols <- symbols[keep]
}

rownames(logcpm_sym) <- symbols
cat("Final symbol-mapped matrix:", nrow(logcpm_sym), "genes x",
    ncol(logcpm_sym), "samples\n\n")

# ============================================================
# STEP 3: Load MSigDB gene sets
# ============================================================
cat("--- Step 3: Load MSigDB gene sets ---\n")

# Hallmark
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat("Hallmark:", length(hallmark_sets), "gene sets\n")

# KEGG_MEDICUS (fall back to KEGG if unavailable)
kegg_df <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2",
          subcollection = "CP:KEGG_MEDICUS"),
  error = function(e) {
    cat("  CP:KEGG_MEDICUS not found, falling back to CP:KEGG\n")
    msigdbr(species = "Homo sapiens", collection = "C2",
            subcollection = "CP:KEGG")
  }
)
kegg_sets <- if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
  split(kegg_df$gene_symbol, kegg_df$gs_name)
} else {
  cat("WARNING: No KEGG gene sets found.\n")
  list()
}
cat("KEGG:", length(kegg_sets), "gene sets\n")

# Reactome
reactome_df <- msigdbr(species = "Homo sapiens", collection = "C2",
                        subcollection = "CP:REACTOME")
reactome_sets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
cat("Reactome:", length(reactome_sets), "gene sets\n")

# Combine all
all_sets <- c(hallmark_sets, kegg_sets, reactome_sets)
cat("\nTotal gene sets:", length(all_sets), "\n")

# Filter: at least 10 and at most 500 genes present in our data
set_sizes_raw <- sapply(all_sets, function(gs) sum(gs %in% rownames(logcpm_sym)))
all_sets <- all_sets[set_sizes_raw >= 10 & set_sizes_raw <= 500]
cat("Gene sets after size filter (10-500 mapped genes):", length(all_sets), "\n")

# Record effective set sizes (genes actually present in expression matrix)
set_sizes <- sapply(all_sets, function(gs) sum(gs %in% rownames(logcpm_sym)))

# Tag each set with its collection
set_collection <- fifelse(grepl("^HALLMARK_", names(all_sets)), "Hallmark",
                  fifelse(grepl("^KEGG_|^KEGG_MEDICUS_", names(all_sets)), "KEGG",
                  fifelse(grepl("^REACTOME_", names(all_sets)), "Reactome", "Other")))
cat("  Hallmark:", sum(set_collection == "Hallmark"),
    " KEGG:", sum(set_collection == "KEGG"),
    " Reactome:", sum(set_collection == "Reactome"), "\n\n")

# ============================================================
# STEP 4: Build decoupleR network data frame
# ============================================================
cat("--- Step 4: Build decoupleR network ---\n")

# decoupleR expects: source (pathway), target (gene), mor (mode of regulation)
net_list <- lapply(names(all_sets), function(pw) {
  genes <- intersect(all_sets[[pw]], rownames(logcpm_sym))
  data.frame(source = pw, target = genes, mor = 1, stringsAsFactors = FALSE)
})
net <- do.call(rbind, net_list)
net <- as.data.frame(net)
cat("Network:", nrow(net), "edges,", length(unique(net$source)), "pathways,",
    length(unique(net$target)), "genes\n\n")

# ============================================================
# STEP 5: Run decoupleR ULM
# ============================================================
cat("--- Step 5: Run ULM (univariate linear model) ---\n")
t0 <- Sys.time()

ulm_res <- run_ulm(
  mat      = logcpm_sym,
  network  = net,
  .source  = "source",
  .target  = "target",
  .mor     = "mor",
  minsize  = 10
)

t1 <- Sys.time()
cat("ULM completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")
cat("ULM result rows:", nrow(ulm_res), "\n")

# Pivot to pathway × sample matrix
ulm_wide <- dcast(as.data.table(ulm_res), source ~ condition, value.var = "score")
ulm_mat <- as.matrix(ulm_wide[, -1, with = FALSE])
rownames(ulm_mat) <- ulm_wide$source
cat("ULM matrix:", nrow(ulm_mat), "pathways x", ncol(ulm_mat), "samples\n")
cat("ULM score range:", round(range(ulm_mat, na.rm = TRUE), 3), "\n\n")

# ============================================================
# STEP 6: Run decoupleR MLM
# ============================================================
cat("--- Step 6: Run MLM (multivariate linear model) ---\n")
t0 <- Sys.time()

mlm_res <- tryCatch(
  run_mlm(
    mat      = logcpm_sym,
    network  = net,
    .source  = "source",
    .target  = "target",
    .mor     = "mor",
    minsize  = 10
  ),
  error = function(e) {
    cat("WARNING: MLM failed (likely colinear pathways):", conditionMessage(e), "\n")
    cat("Skipping MLM — ULM is the primary size-corrected method.\n")
    NULL
  }
)

t1 <- Sys.time()
if (!is.null(mlm_res)) {
  cat("MLM completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")
  cat("MLM result rows:", nrow(mlm_res), "\n")
  mlm_wide <- dcast(as.data.table(mlm_res), source ~ condition, value.var = "score")
  mlm_mat <- as.matrix(mlm_wide[, -1, with = FALSE])
  rownames(mlm_mat) <- mlm_wide$source
  cat("MLM matrix:", nrow(mlm_mat), "pathways x", ncol(mlm_mat), "samples\n")
  cat("MLM score range:", round(range(mlm_mat, na.rm = TRUE), 3), "\n\n")
} else {
  mlm_mat <- NULL
  cat("MLM skipped — will not be included in outputs.\n\n")
}

# ============================================================
# STEP 7: Run ssGSEA (legacy, for comparison)
# ============================================================
cat("--- Step 7: Run ssGSEA (legacy comparison) ---\n")
t0 <- Sys.time()

gsva_version <- packageVersion("GSVA")
cat("GSVA version:", as.character(gsva_version), "\n")

if (gsva_version >= "1.50") {
  ssgsea_param <- ssgseaParam(exprData = logcpm_sym, geneSets = all_sets,
                               normalize = TRUE)
  ssgsea_scores <- gsva(ssgsea_param, verbose = TRUE,
                         BPPARAM = MulticoreParam(NCORES))
} else {
  ssgsea_scores <- gsva(logcpm_sym, all_sets, method = "ssgsea",
                         kcdf = "Gaussian", parallel.sz = NCORES, verbose = TRUE)
}

t1 <- Sys.time()
cat("\nssGSEA completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")
cat("ssGSEA matrix:", nrow(ssgsea_scores), "pathways x", ncol(ssgsea_scores),
    "samples\n")
cat("ssGSEA score range:", round(range(ssgsea_scores), 3), "\n\n")

# ============================================================
# STEP 8: Size-bias diagnostic
# ============================================================
cat("--- Step 8: Size-bias diagnostic ---\n")

# For each method × collection, compute Spearman rho(|mean_score|, set_size)
methods <- list(
  ULM    = ulm_mat,
  MLM    = mlm_mat,
  ssGSEA = ssgsea_scores
)

diag_rows <- list()
for (mname in names(methods)) {
  mat <- methods[[mname]]
  if (is.null(mat)) { cat("  Skipping", mname, "(NULL)\n"); next }
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  # Mean absolute score per pathway
  mean_abs <- rowMeans(abs(mat), na.rm = TRUE)

  for (coll in c("Hallmark", "KEGG", "Reactome", "All")) {
    if (coll == "All") {
      idx_coll <- seq_along(set_collection)
    } else {
      idx_coll <- which(set_collection == coll)
    }

    # Intersect with pathways present in this method's matrix
    pw_in_mat <- intersect(names(all_sets)[idx_coll], rownames(mat))
    if (length(pw_in_mat) < 5) next

    sizes_sub  <- set_sizes[pw_in_mat]
    scores_sub <- mean_abs[pw_in_mat]

    rho_test <- cor.test(scores_sub, sizes_sub, method = "spearman",
                         exact = FALSE)

    diag_rows[[length(diag_rows) + 1]] <- data.table(
      method     = mname,
      collection = coll,
      n_sets     = length(pw_in_mat),
      spearman_rho       = round(rho_test$estimate, 4),
      spearman_p         = signif(rho_test$p.value, 4),
      mean_set_size      = round(mean(sizes_sub), 1),
      median_set_size    = median(sizes_sub),
      mean_abs_score     = round(mean(scores_sub), 4)
    )
  }
}

diag_dt <- rbindlist(diag_rows)
cat("\nSize-bias diagnostic:\n")
print(diag_dt)

# ============================================================
# STEP 9: Save all outputs
# ============================================================
cat("\n--- Step 9: Save outputs ---\n")

saveRDS(ulm_mat, file.path(OUTDIR, "pathway_scores_ulm.rds"))
cat("Saved:", file.path(OUTDIR, "pathway_scores_ulm.rds"), "\n")

if (!is.null(mlm_mat)) {
  saveRDS(mlm_mat, file.path(OUTDIR, "pathway_scores_mlm.rds"))
  cat("Saved:", file.path(OUTDIR, "pathway_scores_mlm.rds"), "\n")
} else {
  cat("MLM skipped — no file saved.\n")
}

if (!is.matrix(ssgsea_scores)) ssgsea_scores <- as.matrix(ssgsea_scores)
saveRDS(ssgsea_scores, file.path(OUTDIR, "pathway_scores_ssgsea.rds"))
cat("Saved:", file.path(OUTDIR, "pathway_scores_ssgsea.rds"), "\n")

fwrite(diag_dt, file.path(OUTDIR, "pathway_size_correction_diagnostic.csv"))
cat("Saved:", file.path(OUTDIR, "pathway_size_correction_diagnostic.csv"), "\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== Summary ===\n")
cat("ULM:    ", nrow(ulm_mat), "pathways x", ncol(ulm_mat), "samples\n")
if (!is.null(mlm_mat)) cat("MLM:    ", nrow(mlm_mat), "pathways x", ncol(mlm_mat), "samples\n") else cat("MLM:     skipped (colinear pathways)\n")
cat("ssGSEA: ", nrow(ssgsea_scores), "pathways x", ncol(ssgsea_scores), "samples\n")

# Highlight the key finding
cat("\nSize-bias summary (All collections, spearman rho of |mean_score| vs set_size):\n")
for (m in c("ULM", "MLM", "ssGSEA")) {
  row <- diag_dt[method == m & collection == "All"]
  if (nrow(row) == 1) {
    flag <- if (abs(row$spearman_rho) < 0.15) "OK (size-independent)" else
            "BIASED (size-dependent)"
    cat(sprintf("  %-7s: rho = %+.3f  (p = %.2e)  %s\n",
                m, row$spearman_rho, row$spearman_p, flag))
  }
}

cat("\nFinished:", as.character(Sys.time()), "\n")
