#!/usr/bin/env Rscript
# 147_transition_subtyping.R
# ---------------------------------------------------------------------------
# Transition Activity Scores (TAS): ssGSEA-based model that assigns 7
# continuous transition activity scores to each of 1,444 MASLD patients.
#
# Pipeline:
#   1. Load merged DGE (1,444 samples) + transition programs + metadata
#   2. Build 14 directional gene sets (7 transitions x UP/DOWN)
#   3. Compute ssGSEA scores using GSVA package
#   4. Compute TAS per transition: TAS_T = ssGSEA_UP - ssGSEA_DOWN, z-norm
#   5. Assign dominant transition per patient (stage-adjacency constrained)
#   6. Run 6-test validation battery
#   7. Save outputs
#
# Input:
#   - results/integration/merged_dge.rds (34K genes x 1,444 samples)
#   - results/progression/transition_programs.csv (180K rows)
#   - results/staging_classifier/modeling_metadata.csv
#   - results/progression/consensus_pseudotime.csv
#   - results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv
#   - results/progression/fate_probabilities.csv
#
# Output (all to results/progression/):
#   - transition_genesets.csv
#   - transition_activity_scores.csv
#   - transition_subtype_assignments.csv
#   - transition_subtype_validation.csv
#
# SLURM: io, 8 CPUs, 32G, 48h
# Env:   micromamba activate rnaseq
#
# Usage:
#   sbatch --job-name=stg147_tas \
#          --partition=io --cpus-per-task=8 --mem=32G --time=48:00:00 \
#          --output=logs/147_transition_subtyping_%j.out \
#          --error=logs/147_transition_subtyping_%j.err \
#          --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && \
#                  micromamba activate rnaseq && \
#                  cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \
#                  Rscript 147_transition_subtyping.R'"
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(msigdbr)
})

# Install GSVA if not present
if (!requireNamespace("GSVA", quietly = TRUE)) {
  cat("GSVA not found. Installing from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}
library(GSVA)

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "progression")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(INT, "scripts/logs"), showWarnings = FALSE, recursive = TRUE)

NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

cat("=== 147: Transition Activity Scores (TAS) ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Cores:", NCORES, "\n\n")

# ============================================================
# STEP 1: Load data
# ============================================================
cat("--- Step 1: Load data ---\n")

# 1a. Merged DGE -> TMM logCPM
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
cat("  DGE:", nrow(dge), "genes x", ncol(dge), "samples\n")
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("  logCPM matrix:", nrow(logcpm), "x", ncol(logcpm), "\n")

# 1b. Gene symbol mapping
cat("Loading gene annotation...\n")
annot <- fread(file.path(RDIR, "gene_annotation/human_ensg_to_symbol.tsv"))
# annot has: gene_id, symbol, gene_type, gene_base
annot <- annot[!is.na(symbol) & symbol != ""]
# Create clean Ensembl ID (strip version) for matching to logcpm rownames
annot[, ensg_clean := sub("\\..*", "", gene_id)]
annot <- annot[!duplicated(ensg_clean)]

# Map logCPM rownames (Ensembl with version) to symbols
ensg_clean <- sub("\\..*", "", rownames(logcpm))
idx <- match(ensg_clean, annot$ensg_clean)
mapped <- !is.na(idx)
cat("  Genes mapped to symbols:", sum(mapped), "of", nrow(logcpm),
    sprintf("(%.1f%%)\n", 100 * sum(mapped) / nrow(logcpm)))

logcpm_sym <- logcpm[mapped, ]
symbols <- annot$symbol[idx[mapped]]

# Handle duplicate symbols: keep the one with highest variance
var_expr <- apply(logcpm_sym, 1, var)
dup_syms <- symbols[duplicated(symbols)]
if (length(dup_syms) > 0) {
  cat("  Resolving", length(unique(dup_syms)), "duplicate symbols (keeping highest variance)\n")
  keep <- rep(TRUE, length(symbols))
  for (s in unique(dup_syms)) {
    which_dup <- which(symbols == s)
    best <- which_dup[which.max(var_expr[which_dup])]
    keep[setdiff(which_dup, best)] <- FALSE
  }
  logcpm_sym <- logcpm_sym[keep, ]
  symbols <- symbols[keep]
}
rownames(logcpm_sym) <- symbols
cat("  Final symbol-mapped matrix:", nrow(logcpm_sym), "genes x", ncol(logcpm_sym), "samples\n")

# 1c. Transition programs
cat("Loading transition programs...\n")
tp <- fread(file.path(OUTDIR, "transition_programs.csv"))
cat("  Transition programs:", nrow(tp), "rows,", uniqueN(tp$transition), "transitions\n")
cat("  Transitions:", paste(unique(tp$transition), collapse = ", "), "\n")

# 1d. Modeling metadata
cat("Loading modeling metadata...\n")
meta <- fread(file.path(RDIR, "staging_classifier/modeling_metadata.csv"))
meta <- meta[match(colnames(logcpm_sym), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(logcpm_sym)))
cat("  Metadata matched:", nrow(meta), "samples\n")

# 1e. Pseudotime
cat("Loading consensus pseudotime...\n")
ptime <- fread(file.path(OUTDIR, "consensus_pseudotime.csv"))
cat("  Pseudotime:", nrow(ptime), "samples\n")

# 1f. BayesPrism cell-type proportions
cat("Loading BayesPrism proportions...\n")
cellprop <- fread(file.path(OUTDIR, "cibersortx_celltype_expression/bayesprism_proportions.csv"))
cat("  Cell-type proportions:", nrow(cellprop), "samples x", ncol(cellprop) - 1, "cell types\n")

# 1g. Fate probabilities
cat("Loading fate probabilities...\n")
fate <- fread(file.path(OUTDIR, "fate_probabilities.csv"))
cat("  Fate probabilities:", nrow(fate), "samples\n")

cat("\n")

# ============================================================
# STEP 2: Build 14 directional gene sets (7 transitions x UP/DOWN)
# ============================================================
cat("--- Step 2: Build directional gene sets ---\n")

transitions <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4",
                  "NAS01_to_NAS24", "NAS24_to_NAS5", "NAS5_to_NAS68")

gene_sets <- list()
geneset_info <- list()

for (tr in transitions) {
  sub <- tp[transition == tr & !is.na(gene_symbol) & gene_symbol != ""]

  # Primary filter: tau > 0.8, adj.P.Val < 0.05
  sig <- sub[tau > 0.8 & adj.P.Val < 0.05]

  # If fewer than 50 genes, relax to tau > 0.6
  tau_threshold <- 0.8
  if (nrow(sig) < 50) {
    sig <- sub[tau > 0.6 & adj.P.Val < 0.05]
    tau_threshold <- 0.6
    cat(sprintf("  %s: relaxed to tau > 0.6 (%d genes)\n", tr, nrow(sig)))
  }

  # Split UP/DOWN
  sig_up   <- sig[logFC > 0][order(-abs(t))]
  sig_down <- sig[logFC < 0][order(-abs(t))]

  # Cap at 200 genes each by |t-statistic|
  if (nrow(sig_up) > 200)   sig_up   <- sig_up[1:200]
  if (nrow(sig_down) > 200) sig_down <- sig_down[1:200]

  # Filter to genes present in our expression matrix
  sig_up   <- sig_up[gene_symbol %in% rownames(logcpm_sym)]
  sig_down <- sig_down[gene_symbol %in% rownames(logcpm_sym)]

  cat(sprintf("  %s: tau>%.1f => UP=%d, DOWN=%d genes\n",
              tr, tau_threshold, nrow(sig_up), nrow(sig_down)))

  # Store gene sets for ssGSEA (need at least 5 genes to run)
  up_name   <- paste0(tr, "_UP")
  down_name <- paste0(tr, "_DOWN")

  if (nrow(sig_up) >= 5) {
    gene_sets[[up_name]] <- unique(sig_up$gene_symbol)
  } else {
    cat(sprintf("    WARNING: %s has only %d UP genes, skipping\n", tr, nrow(sig_up)))
  }
  if (nrow(sig_down) >= 5) {
    gene_sets[[down_name]] <- unique(sig_down$gene_symbol)
  } else {
    cat(sprintf("    WARNING: %s has only %d DOWN genes, skipping\n", tr, nrow(sig_down)))
  }

  # Collect info for output
  if (nrow(sig_up) > 0) {
    geneset_info[[length(geneset_info) + 1]] <- sig_up[, .(
      gene_symbol, transition = tr, direction = "UP",
      tau, padj = adj.P.Val, t_stat = t
    )]
  }
  if (nrow(sig_down) > 0) {
    geneset_info[[length(geneset_info) + 1]] <- sig_down[, .(
      gene_symbol, transition = tr, direction = "DOWN",
      tau, padj = adj.P.Val, t_stat = t
    )]
  }
}

geneset_dt <- rbindlist(geneset_info)
cat("\nTotal gene set genes:", nrow(geneset_dt), "\n")
cat("Gene sets created:", length(gene_sets), "of 14 expected\n")
cat("Gene set sizes:\n")
for (nm in names(gene_sets)) {
  cat(sprintf("  %-20s: %d genes\n", nm, length(gene_sets[[nm]])))
}
cat("\n")

# ============================================================
# STEP 3: Compute ssGSEA scores using GSVA
# ============================================================
cat("--- Step 3: Compute ssGSEA scores ---\n")
cat("This may take 10-30 minutes for", ncol(logcpm_sym), "samples and",
    length(gene_sets), "gene sets...\n")
t0 <- Sys.time()

# Detect GSVA API version (>=1.50 uses new param objects)
gsva_version <- packageVersion("GSVA")
cat("GSVA version:", as.character(gsva_version), "\n")

if (gsva_version >= "1.50") {
  ssgsea_param <- ssgseaParam(exprData = logcpm_sym, geneSets = gene_sets,
                               normalize = TRUE)
  ssgsea_scores <- gsva(ssgsea_param, verbose = TRUE,
                         BPPARAM = BiocParallel::MulticoreParam(NCORES))
} else {
  ssgsea_scores <- gsva(logcpm_sym, gene_sets, method = "ssgsea",
                         kcdf = "Gaussian", parallel.sz = NCORES, verbose = TRUE)
}

t1 <- Sys.time()
cat("\nssGSEA completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")
cat("ssGSEA score matrix:", nrow(ssgsea_scores), "gene sets x", ncol(ssgsea_scores), "samples\n")
cat("Score range:", round(range(ssgsea_scores), 3), "\n\n")

# ============================================================
# STEP 4: Compute TAS per transition (UP - DOWN, then z-normalize)
# ============================================================
cat("--- Step 4: Compute Transition Activity Scores ---\n")

tas_mat <- matrix(NA_real_, nrow = ncol(logcpm_sym), ncol = length(transitions),
                  dimnames = list(colnames(logcpm_sym), transitions))

for (tr in transitions) {
  up_name   <- paste0(tr, "_UP")
  down_name <- paste0(tr, "_DOWN")

  has_up   <- up_name %in% rownames(ssgsea_scores)
  has_down <- down_name %in% rownames(ssgsea_scores)

  if (has_up && has_down) {
    raw <- ssgsea_scores[up_name, ] - ssgsea_scores[down_name, ]
  } else if (has_up) {
    raw <- ssgsea_scores[up_name, ]
    cat(sprintf("  %s: only UP available (no DOWN set)\n", tr))
  } else if (has_down) {
    raw <- -ssgsea_scores[down_name, ]
    cat(sprintf("  %s: only DOWN available (no UP set)\n", tr))
  } else {
    cat(sprintf("  WARNING: %s has neither UP nor DOWN gene set, setting to 0\n", tr))
    raw <- rep(0, ncol(logcpm_sym))
  }

  # Z-normalize across samples
  mu <- mean(raw, na.rm = TRUE)
  sd_val <- sd(raw, na.rm = TRUE)
  if (sd_val > 0) {
    tas_mat[, tr] <- (raw - mu) / sd_val
  } else {
    tas_mat[, tr] <- 0
  }
}

cat("TAS matrix:", nrow(tas_mat), "samples x", ncol(tas_mat), "transitions\n")
cat("TAS ranges per transition:\n")
for (tr in transitions) {
  cat(sprintf("  %-20s: [%.2f, %.2f], mean=%.3f, sd=%.3f\n",
              tr, min(tas_mat[, tr], na.rm = TRUE), max(tas_mat[, tr], na.rm = TRUE),
              mean(tas_mat[, tr], na.rm = TRUE), sd(tas_mat[, tr], na.rm = TRUE)))
}
cat("\n")

# ============================================================
# STEP 5: Assign dominant transition per patient
# ============================================================
cat("--- Step 5: Assign dominant transition ---\n")

# Define stage adjacency: which transitions are allowed per fibrosis stage
fib_transitions <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
nas_transitions_list <- c("NAS01_to_NAS24", "NAS24_to_NAS5", "NAS5_to_NAS68")

# Fibrosis adjacency: each stage can be in the "from" or "to" of adjacent transitions
# F0 -> can be F0_to_F1
# F1 -> can be F0_to_F1 or F1_to_F2
# F2 -> can be F1_to_F2 or F2_to_F3
# F3 -> can be F2_to_F3 or F3_to_F4
# F4 -> can be F3_to_F4
fib_adjacency <- list(
  "0" = c("F0_to_F1"),
  "1" = c("F0_to_F1", "F1_to_F2"),
  "2" = c("F1_to_F2", "F2_to_F3"),
  "3" = c("F2_to_F3", "F3_to_F4"),
  "4" = c("F3_to_F4")
)

# NAS adjacency (using NAS groups from transition programs)
# NAS 0-1 -> NAS01_to_NAS24
# NAS 2-4 -> NAS01_to_NAS24 or NAS24_to_NAS5
# NAS 5   -> NAS24_to_NAS5 or NAS5_to_NAS68
# NAS 6-8 -> NAS5_to_NAS68
nas_adjacency <- list(
  "NAS_0_1" = c("NAS01_to_NAS24"),
  "NAS_2_4" = c("NAS01_to_NAS24", "NAS24_to_NAS5"),
  "NAS_5"   = c("NAS24_to_NAS5", "NAS5_to_NAS68"),
  "NAS_6_8" = c("NAS5_to_NAS68")
)

# Map NAS score to NAS group
nas_score_to_group <- function(ns) {
  ifelse(is.na(ns), NA_character_,
    ifelse(ns <= 1, "NAS_0_1",
      ifelse(ns <= 4, "NAS_2_4",
        ifelse(ns == 5, "NAS_5", "NAS_6_8"))))
}

# Assign dominant fibrosis transition
assign_dominant <- function(tas_row, allowed_transitions, all_transitions, quiescent_threshold = 1.0) {
  # Check if all allowed TAS scores are below threshold -> quiescent
  allowed_vals <- tas_row[all_transitions %in% allowed_transitions]
  if (length(allowed_vals) == 0 || all(is.na(allowed_vals))) return("unassigned")
  if (all(abs(allowed_vals) < quiescent_threshold, na.rm = TRUE)) return("quiescent")
  # Pick the allowed transition with highest TAS
  names(allowed_vals) <- all_transitions[all_transitions %in% allowed_transitions]
  return(names(which.max(allowed_vals)))
}

# Build assignment table
assign_dt <- data.table(
  sample_id = colnames(logcpm_sym),
  fibrosis_stage = meta$fibrosis_stage,
  nas_score = meta$nas_score
)

# Fibrosis dominant
assign_dt$fib_dominant <- NA_character_
for (i in seq_len(nrow(assign_dt))) {
  fib_s <- as.character(assign_dt$fibrosis_stage[i])
  if (is.na(fib_s) || fib_s == "" || !fib_s %in% names(fib_adjacency)) {
    assign_dt$fib_dominant[i] <- "unassigned"
  } else {
    allowed <- fib_adjacency[[fib_s]]
    assign_dt$fib_dominant[i] <- assign_dominant(
      tas_mat[i, ], allowed, fib_transitions)
  }
}

# NAS dominant
assign_dt$nas_group <- nas_score_to_group(assign_dt$nas_score)
assign_dt$nas_dominant <- NA_character_
for (i in seq_len(nrow(assign_dt))) {
  ng <- assign_dt$nas_group[i]
  if (is.na(ng) || !ng %in% names(nas_adjacency)) {
    assign_dt$nas_dominant[i] <- "unassigned"
  } else {
    allowed <- nas_adjacency[[ng]]
    assign_dt$nas_dominant[i] <- assign_dominant(
      tas_mat[i, ], allowed, nas_transitions_list)
  }
}

# Quiescent flag (all |TAS_z| < 1.0 for both fibrosis and NAS)
assign_dt$is_quiescent <- apply(tas_mat, 1, function(x) all(abs(x) < 1.0, na.rm = TRUE))

cat("Fibrosis dominant assignments:\n")
print(table(assign_dt$fib_dominant, useNA = "ifany"))
cat("\nNAS dominant assignments:\n")
print(table(assign_dt$nas_dominant, useNA = "ifany"))
cat("\nQuiescent samples:", sum(assign_dt$is_quiescent), "of", nrow(assign_dt),
    sprintf("(%.1f%%)\n", 100 * sum(assign_dt$is_quiescent) / nrow(assign_dt)))
cat("\n")

# ============================================================
# STEP 6: Validation battery (6 tests)
# ============================================================
cat("--- Step 6: Validation battery ---\n")

validation_results <- list()

# Helper: safe boolean mask that replaces NA with FALSE
safe_mask <- function(...) {
  m <- Reduce(`&`, list(...))
  m[is.na(m)] <- FALSE
  m
}

# ---------- V1: Cell-type concordance ----------
# F1->F2 dominant patients should have higher stellate fraction than others at same fibrosis stage
cat("V1: Cell-type concordance (stellate in F1->F2 dominant)...\n")
v1_result <- tryCatch({
  # Merge cell-type proportions
  cp <- cellprop[match(assign_dt$sample_id, cellprop$sample_id)]
  # Focus on fibrosis stage 1 and 2 (where F1_to_F2 is allowed)
  mask_f12 <- safe_mask(
    as.character(assign_dt$fibrosis_stage) %in% c("1", "2"),
    assign_dt$fib_dominant != "unassigned"
  )
  if (sum(mask_f12) < 10) stop("Too few samples for V1")
  is_f1f2_dom <- safe_mask(assign_dt$fib_dominant == "F1_to_F2", mask_f12)
  is_other    <- safe_mask(assign_dt$fib_dominant != "F1_to_F2", mask_f12)
  if (sum(is_f1f2_dom) < 3 || sum(is_other) < 3) stop("Too few per group for V1")
  wt <- wilcox.test(cp$Stellate[is_f1f2_dom], cp$Stellate[is_other], alternative = "greater")
  cat(sprintf("  F1->F2 dominant Stellate: median=%.4f vs others=%.4f, p=%.3e\n",
              median(cp$Stellate[is_f1f2_dom], na.rm = TRUE),
              median(cp$Stellate[is_other], na.rm = TRUE),
              wt$p.value))
  data.table(test = "V1_celltype_concordance",
             description = "Stellate higher in F1->F2 dominant vs others at F1/F2",
             statistic = wt$statistic, p_value = wt$p.value,
             effect_size = median(cp$Stellate[is_f1f2_dom], na.rm = TRUE) -
                           median(cp$Stellate[is_other], na.rm = TRUE),
             pass = wt$p.value < 0.05)
}, error = function(e) {
  cat("  V1 FAILED:", conditionMessage(e), "\n")
  data.table(test = "V1_celltype_concordance",
             description = conditionMessage(e),
             statistic = NA_real_, p_value = NA_real_,
             effect_size = NA_real_, pass = FALSE)
})
validation_results[[1]] <- v1_result

# ---------- V2: Pseudotime ordering ----------
# Within-stage: higher TAS for "forward" transition correlates with higher pseudotime
cat("V2: Pseudotime ordering (within-stage TAS vs pseudotime)...\n")
v2_result <- tryCatch({
  pt <- ptime[match(assign_dt$sample_id, ptime$sample_id)]
  # Use consensus pseudotime
  pt_vals <- pt$pseudotime_consensus
  # Within fibrosis stage 2 (can be F1_to_F2 or F2_to_F3):
  # Forward transition = F2_to_F3, so higher TAS_F2_to_F3 -> higher pseudotime
  mask_f2 <- safe_mask(as.character(assign_dt$fibrosis_stage) == "2", !is.na(pt_vals))
  if (sum(mask_f2) < 10) stop("Too few F2 samples with pseudotime")
  rho <- cor(tas_mat[mask_f2, "F2_to_F3"], pt_vals[mask_f2],
             method = "spearman", use = "complete.obs")
  ct <- cor.test(tas_mat[mask_f2, "F2_to_F3"], pt_vals[mask_f2],
                 method = "spearman")
  cat(sprintf("  F2 stage: TAS_F2_to_F3 vs pseudotime rho=%.3f, p=%.3e (n=%d)\n",
              rho, ct$p.value, sum(mask_f2)))
  data.table(test = "V2_pseudotime_ordering",
             description = "Spearman TAS_F2_to_F3 vs pseudotime within F2",
             statistic = rho, p_value = ct$p.value,
             effect_size = rho, pass = ct$p.value < 0.05 & rho > 0)
}, error = function(e) {
  cat("  V2 FAILED:", conditionMessage(e), "\n")
  data.table(test = "V2_pseudotime_ordering",
             description = conditionMessage(e),
             statistic = NA_real_, p_value = NA_real_,
             effect_size = NA_real_, pass = FALSE)
})
validation_results[[2]] <- v2_result

# ---------- V3: Fate probability concordance ----------
# F2->F3 and F3->F4 dominant patients should have higher P(F4)
cat("V3: Fate probability concordance (advanced dominant vs P(F4))...\n")
v3_result <- tryCatch({
  ft <- fate[match(assign_dt$sample_id, fate$sample_id)]
  fate_f4 <- ft$fate_prob_F4
  # Focus on stages 2 and 3 where both forward and backward transitions are possible
  mask_f23 <- safe_mask(
    as.character(assign_dt$fibrosis_stage) %in% c("2", "3"),
    !is.na(fate_f4),
    assign_dt$fib_dominant != "unassigned"
  )
  if (sum(mask_f23) < 10) stop("Too few F2/F3 samples with fate probabilities")
  # "Advanced" = F2_to_F3 or F3_to_F4 dominant; "Earlier" = others
  is_adv <- safe_mask(assign_dt$fib_dominant %in% c("F2_to_F3", "F3_to_F4"), mask_f23)
  is_early <- safe_mask(!(assign_dt$fib_dominant %in% c("F2_to_F3", "F3_to_F4")), mask_f23)
  if (sum(is_adv) < 3 || sum(is_early) < 3) stop("Too few per group for V3")
  wt <- wilcox.test(fate_f4[is_adv], fate_f4[is_early], alternative = "greater")
  cat(sprintf("  Advanced dominant P(F4): median=%.3f vs early=%.3f, p=%.3e\n",
              median(fate_f4[is_adv], na.rm = TRUE),
              median(fate_f4[is_early], na.rm = TRUE),
              wt$p.value))
  data.table(test = "V3_fate_probability",
             description = "P(F4) higher in F2->F3/F3->F4 dominant vs others at F2/F3",
             statistic = wt$statistic, p_value = wt$p.value,
             effect_size = median(fate_f4[is_adv], na.rm = TRUE) -
                           median(fate_f4[is_early], na.rm = TRUE),
             pass = wt$p.value < 0.05)
}, error = function(e) {
  cat("  V3 FAILED:", conditionMessage(e), "\n")
  data.table(test = "V3_fate_probability",
             description = conditionMessage(e),
             statistic = NA_real_, p_value = NA_real_,
             effect_size = NA_real_, pass = FALSE)
})
validation_results[[3]] <- v3_result

# ---------- V4: Within-stage heterogeneity ----------
# Shannon entropy of dominant assignments per fibrosis stage
cat("V4: Within-stage heterogeneity (Shannon entropy)...\n")
v4_result <- tryCatch({
  entropy_vals <- list()
  for (fstage in c("0", "1", "2", "3", "4")) {
    mask <- safe_mask(
      as.character(assign_dt$fibrosis_stage) == fstage,
      assign_dt$fib_dominant != "unassigned"
    )
    if (sum(mask) < 5) next
    tab <- table(assign_dt$fib_dominant[mask])
    probs <- tab / sum(tab)
    probs <- probs[probs > 0]
    H <- -sum(probs * log2(probs))
    H_max <- log2(length(probs))
    entropy_vals[[fstage]] <- data.table(
      fibrosis_stage = fstage, H = H, H_max = H_max,
      H_norm = ifelse(H_max > 0, H / H_max, 0),
      n = sum(mask), n_categories = length(probs)
    )
    cat(sprintf("  F%s: H=%.3f (H_max=%.3f, H_norm=%.3f, n=%d, k=%d)\n",
                fstage, H, H_max, H / max(H_max, 1e-10), sum(mask), length(probs)))
  }
  ent_dt <- rbindlist(entropy_vals)
  mean_hnorm <- mean(ent_dt$H_norm, na.rm = TRUE)
  cat(sprintf("  Mean normalized entropy: %.3f\n", mean_hnorm))
  # Pass if mean H_norm > 0.3 (indicating real heterogeneity, not trivial assignment)
  data.table(test = "V4_within_stage_heterogeneity",
             description = "Shannon entropy of fib_dominant per stage",
             statistic = mean_hnorm, p_value = NA_real_,
             effect_size = mean_hnorm,
             pass = mean_hnorm > 0.3)
}, error = function(e) {
  cat("  V4 FAILED:", conditionMessage(e), "\n")
  data.table(test = "V4_within_stage_heterogeneity",
             description = conditionMessage(e),
             statistic = NA_real_, p_value = NA_real_,
             effect_size = NA_real_, pass = FALSE)
})
validation_results[[4]] <- v4_result

# ---------- V5: LOCO cross-dataset ----------
# KS-test: are TAS distributions similar across datasets?
cat("V5: LOCO cross-dataset consistency (KS-test)...\n")
v5_result <- tryCatch({
  datasets <- unique(meta$dataset)
  ks_pvals <- c()
  for (tr in fib_transitions[1:2]) {  # Test first 2 fibrosis transitions
    for (i in seq_along(datasets)) {
      for (j in seq(i + 1, length(datasets))) {
        if (j > length(datasets)) break
        d1_mask <- meta$dataset == datasets[i]
        d2_mask <- meta$dataset == datasets[j]
        if (sum(d1_mask) < 10 || sum(d2_mask) < 10) next
        ks <- ks.test(tas_mat[d1_mask, tr], tas_mat[d2_mask, tr])
        ks_pvals <- c(ks_pvals, ks$p.value)
      }
    }
  }
  # Fraction of dataset pairs that are NOT significantly different
  frac_consistent <- mean(ks_pvals > 0.05, na.rm = TRUE)
  cat(sprintf("  %d pairwise KS tests, %.1f%% consistent (p > 0.05)\n",
              length(ks_pvals), 100 * frac_consistent))
  # Pass if > 50% of pairs are consistent
  data.table(test = "V5_loco_cross_dataset",
             description = "Fraction of dataset pairs with KS p>0.05",
             statistic = frac_consistent, p_value = NA_real_,
             effect_size = frac_consistent,
             pass = frac_consistent > 0.5)
}, error = function(e) {
  cat("  V5 FAILED:", conditionMessage(e), "\n")
  data.table(test = "V5_loco_cross_dataset",
             description = conditionMessage(e),
             statistic = NA_real_, p_value = NA_real_,
             effect_size = NA_real_, pass = FALSE)
})
validation_results[[5]] <- v5_result

# ---------- V6: Pathway coherence ----------
# F0->F1 dominant should show immune-related enrichment (ALLOGRAFT_REJECTION / INFLAMMATORY_RESPONSE)
cat("V6: Pathway coherence (Hallmark enrichment of F0->F1 dominant)...\n")
v6_result <- tryCatch({
  hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

  # Get F0->F1 UP gene set
  f0f1_genes <- gene_sets[["F0_to_F1_UP"]]
  if (is.null(f0f1_genes) || length(f0f1_genes) < 10) stop("F0->F1 UP gene set too small")

  # Test overlap with immune/inflammatory hallmarks
  immune_hallmarks <- c("HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INFLAMMATORY_RESPONSE",
                         "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                         "HALLMARK_COMPLEMENT", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  all_genes_in_matrix <- rownames(logcpm_sym)

  best_pval <- 1
  best_pathway <- ""
  best_or <- NA_real_
  for (pw in immune_hallmarks) {
    if (!pw %in% names(hallmark_sets)) next
    pw_genes <- hallmark_sets[[pw]]
    # Fisher's exact test
    a <- sum(f0f1_genes %in% pw_genes)  # in both
    b <- length(f0f1_genes) - a  # in F0F1 only
    c <- sum(pw_genes %in% all_genes_in_matrix) - a  # in pathway only
    d <- length(all_genes_in_matrix) - a - b - c  # in neither
    ft <- fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater")
    if (ft$p.value < best_pval) {
      best_pval <- ft$p.value
      best_pathway <- pw
      best_or <- ft$estimate
    }
  }
  if (is.na(best_or)) stop("No immune hallmarks found in MSigDB")
  cat(sprintf("  Best immune overlap: %s (OR=%.2f, p=%.3e)\n",
              best_pathway, best_or, best_pval))
  data.table(test = "V6_pathway_coherence",
             description = paste0("F0->F1 UP overlap with ", best_pathway),
             statistic = best_or, p_value = best_pval,
             effect_size = best_or,
             pass = best_pval < 0.05)
}, error = function(e) {
  cat("  V6 FAILED:", conditionMessage(e), "\n")
  data.table(test = "V6_pathway_coherence",
             description = conditionMessage(e),
             statistic = NA_real_, p_value = NA_real_,
             effect_size = NA_real_, pass = FALSE)
})
validation_results[[6]] <- v6_result

validation_dt <- rbindlist(validation_results)
cat("\n=== Validation Summary ===\n")
cat(sprintf("  Passed: %d / %d tests\n", sum(validation_dt$pass, na.rm = TRUE), nrow(validation_dt)))
print(validation_dt[, .(test, pass, statistic = round(statistic, 4),
                         p_value = ifelse(is.na(p_value), NA_character_,
                                          formatC(p_value, format = "e", digits = 2)))])
cat("\n")

# ============================================================
# STEP 7: Save outputs
# ============================================================
cat("--- Step 7: Save outputs ---\n")

# 7a. Gene sets
fwrite(geneset_dt, file.path(OUTDIR, "transition_genesets.csv"))
cat("  Saved transition_genesets.csv:", nrow(geneset_dt), "rows\n")

# 7b. Transition activity scores
tas_dt <- data.table(sample_id = rownames(tas_mat))
for (tr in transitions) {
  # Use cleaner column names
  col_name <- paste0("TAS_", tr)
  tas_dt[[col_name]] <- tas_mat[, tr]
}
# Add dominant assignments
tas_dt$fib_dominant <- assign_dt$fib_dominant
tas_dt$nas_dominant <- assign_dt$nas_dominant
fwrite(tas_dt, file.path(OUTDIR, "transition_activity_scores.csv"))
cat("  Saved transition_activity_scores.csv:", nrow(tas_dt), "samples x",
    ncol(tas_dt), "columns\n")

# 7c. Subtype assignments (with stage context)
subtype_dt <- data.table(
  sample_id = assign_dt$sample_id,
  fibrosis_stage = assign_dt$fibrosis_stage,
  nas_score = assign_dt$nas_score,
  fib_dominant = assign_dt$fib_dominant,
  nas_dominant = assign_dt$nas_dominant,
  is_quiescent = assign_dt$is_quiescent
)
fwrite(subtype_dt, file.path(OUTDIR, "transition_subtype_assignments.csv"))
cat("  Saved transition_subtype_assignments.csv:", nrow(subtype_dt), "samples\n")

# 7d. Validation results
fwrite(validation_dt, file.path(OUTDIR, "transition_subtype_validation.csv"))
cat("  Saved transition_subtype_validation.csv:", nrow(validation_dt), "tests\n")

# Final summary
cat("\n=== 147: COMPLETE ===\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("Key outputs in", OUTDIR, ":\n")
cat("  transition_genesets.csv         -", nrow(geneset_dt), "gene set members\n")
cat("  transition_activity_scores.csv  -", nrow(tas_dt), "samples x 7 TAS\n")
cat("  transition_subtype_assignments.csv -", nrow(subtype_dt), "samples\n")
cat("  transition_subtype_validation.csv  -", nrow(validation_dt), "validation tests\n")
