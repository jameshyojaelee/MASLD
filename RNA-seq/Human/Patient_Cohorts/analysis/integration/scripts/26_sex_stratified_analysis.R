#!/usr/bin/env Rscript
# 26_sex_stratified_analysis.R
# ---------------------------------------------------------------------------
# Strategy 4 — Sex-stratified MASLD analysis
#
#   1. Dream interaction model: group_binary * inferred_sex + (1|dataset)
#      (NO deconvolution covariates — those belong in script 25 only)
#   2. Male-only and female-only dream runs (within-stratum DE)
#   3. Sex-DEG classification: Male_specific / Female_specific / Shared / Divergent
#   4. Pathway enrichment by sex (fgsea / Hallmark)
#   5. Cross-species sex concordance (human <-> mouse MCD)
#   6. Figures: M-vs-F LFC scatter, pathway NES dot plot
#
# Output: results/integration/sex_interaction_dream.csv
#         results/integration/sex_stratified_results_male.csv
#         results/integration/sex_stratified_results_female.csv
#         results/integration/sex_deg_classification.csv
#         results/integration/sex_fgsea_male.csv
#         results/integration/sex_fgsea_female.csv
#         results/integration/sex_cross_species_concordance.csv
#         results/integration/sex_lfc_scatter.pdf
#         results/integration/sex_pathway_nes_dotplot.pdf
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE     <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT      <- file.path(BASE, "analysis/integration")
RDIR     <- file.path(INT, "results/integration")
MOUSE_SEX <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/InHouse_MCD/meta_analysis/sex_bias_mcd_specific_genes.csv"
ORTHO_MAP <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"

# ---------------------------------------------------------------------------
# Parallel setup — pass BPPARAM explicitly to dream(); do NOT use register()
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# ---------------------------------------------------------------------------
# Model selection (DREAM_MODEL env var)
# ---------------------------------------------------------------------------
model_choice <- Sys.getenv("DREAM_MODEL", "base")
cat("Dream model variant:", model_choice, "\n")

# ============================================================================
# Load data (shared across all dream runs)
# ============================================================================
cat("\n===== Loading data =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Cohort selection: enforce config/human_datasets.yaml `include_in_mega` (same
# as Script 05). Without this, sex-stratified DE would be on 8 cohorts while
# the primary mega-analysis is on 5, producing mismatched downstream tables
# (atlas reads sex_deg_classification.csv + sex_interaction_dream.csv).
ycfg <- yaml::read_yaml(file.path(
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
  "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("yaml-declared mega cohorts (k =", length(mega_cohorts), "):",
    paste(mega_cohorts, collapse = ", "), "\n")

keep_samples <- dge$samples$dataset %in% mega_cohorts
excl_extra <- Sys.getenv("EXCL_EXTRA", "")
if (nchar(excl_extra) > 0) {
  excl_extra_vec <- trimws(strsplit(excl_extra, ",")[[1]])
  keep_samples <- keep_samples & !dge$samples$dataset %in% excl_extra_vec
  cat("  Excluding extra cohorts:", paste(excl_extra_vec, collapse=", "), "\n")
}
dge_mega <- dge[, keep_samples]
cat("Samples for mega-analysis:", ncol(dge_mega), "\n")

# Filter to samples that passed the sex check (pass_sex == TRUE in QC report).
# Since script 03 now uses pass_technical (not pass_all), the DGE may contain
# samples that failed only the sex check. Remove them for sex-stratified analysis.
qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) {
  cat("Removing", sum(sex_fail), "samples that failed sex check (pass_sex=FALSE)\n")
  dge_mega <- dge_mega[, !sex_fail]
}
cat("Samples after sex check filter:", ncol(dge_mega), "\n")

# Load metadata with inferred sex
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# NOTE: Deconvolution covariates (Hepatocytes, Macrophages) are intentionally
# EXCLUDED from this primary model. They belong exclusively in script 25
# (25_deconv_attribution.R). Including them here absorbs composition-driven
# disease signal and collapses DEG counts to near-zero.
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

# Drop samples with NA sex
na_sex <- is.na(info$inferred_sex)
cat("Samples with NA inferred_sex:", sum(na_sex), "\n")
if (any(na_sex)) {
  dge_mega <- dge_mega[, !na_sex]
  info     <- info[!na_sex, , drop = FALSE]
  cat("Samples remaining after dropping NA sex:", ncol(dge_mega), "\n")
}

cat("Group x sex distribution:\n")
print(table(info$group_binary, info$inferred_sex))
cat("Group x dataset distribution:\n")
print(table(info$group_binary, info$dataset))

# ============================================================================
# PART 1 — Dream interaction model (group_binary * inferred_sex)
# ============================================================================
cat("\n===== PART 1: Dream interaction model =====\n")

form_int <- switch(model_choice,
  "base"      = ~ group_binary * inferred_sex + (1|dataset),
  "randslope" = ~ group_binary * inferred_sex + (1 + group_binary | dataset),
  stop("Unknown DREAM_MODEL: ", model_choice)
)
cat("Interaction formula:", deparse(form_int), "\n")

cat("Running voomWithDreamWeights (interaction)...\n")
v_int <- suppressWarnings(voomWithDreamWeights(dge_mega, form_int, info, BPPARAM = param))

cat("Running dream (interaction)...\n")
fit_int <- suppressWarnings(dream(v_int, form_int, info, BPPARAM = param))

# Find the interaction coefficient (handles varying factor level names)
all_coefs <- colnames(fit_int$coefficients)
cat("Available coefficients:", paste(all_coefs, collapse = ", "), "\n")
int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary", all_coefs, value = TRUE)
if (length(int_coef) == 0) stop("No interaction coefficient found among: ", paste(all_coefs, collapse = ", "))
if (length(int_coef) > 1) int_coef <- int_coef[1]
cat("Using interaction coefficient:", int_coef, "\n")

# Extract the interaction coefficient
res_int <- topTable(fit_int, coef = int_coef,
                    number = Inf, sort.by = "none")
res_int$gene <- rownames(res_int)

res_int_dt <- as.data.table(res_int)
setnames(res_int_dt, "adj.P.Val", "padj")

sig_int <- res_int_dt[padj < 0.1]
cat("\n===== DREAM INTERACTION RESULTS =====\n")
cat("Total genes tested:", nrow(res_int_dt), "\n")
cat("Sex-differential DEGs (interaction padj < 0.1):", nrow(sig_int), "\n")
cat("  Stronger in males (interaction logFC > 0):", sum(sig_int$logFC > 0), "\n")
cat("  Stronger in females (interaction logFC < 0):", sum(sig_int$logFC < 0), "\n")

fwrite(res_int_dt, file.path(RDIR, "sex_interaction_dream.csv"))
cat("Saved: sex_interaction_dream.csv\n")

# ============================================================================
# PART 2 — Male-only and female-only dream runs (within-stratum DE)
# ============================================================================
cat("\n===== PART 2: Sex-stratified dream runs =====\n")

run_stratum_dream <- function(sex_label, dge_full, info_full, param, rdir) {
  cat("\n--- Stratum:", sex_label, "---\n")

  info_s <- info_full[info_full$inferred_sex == sex_label, , drop = FALSE]
  sample_ids_s <- rownames(info_s)
  dge_s <- dge_full[, colnames(dge_full) %in% sample_ids_s]

  cat("  Samples:", ncol(dge_s), "\n")
  cat("  Group distribution:\n")
  print(table(info_s$group_binary, info_s$dataset))

  # Require both groups present in at least one dataset
  n_ctrl <- sum(info_s$group_binary == "Control")
  n_dis  <- sum(info_s$group_binary == "Disease")
  if (n_ctrl == 0 || n_dis == 0) {
    warning("  Stratum ", sex_label, " missing controls or disease samples — skipping.")
    return(NULL)
  }

  # Drop datasets with only one group (would cause singular fit)
  ds_counts <- table(info_s$dataset, info_s$group_binary)
  keep_ds <- rownames(ds_counts)[apply(ds_counts, 1, function(x) all(x > 0))]
  if (length(keep_ds) < nrow(ds_counts)) {
    dropped <- setdiff(rownames(ds_counts), keep_ds)
    cat("  Dropping single-group datasets:", paste(dropped, collapse = ", "), "\n")
    info_s <- info_s[info_s$dataset %in% keep_ds, , drop = FALSE]
    info_s$dataset <- droplevels(info_s$dataset)
    dge_s  <- dge_s[, colnames(dge_s) %in% rownames(info_s)]
  }

  # Require at least 2 datasets for random effect
  n_ds <- nlevels(info_s$dataset)
  if (n_ds < 2) {
    warning("  Only ", n_ds, " dataset(s) in stratum ", sex_label,
            " — random effect not estimable; skipping.")
    return(NULL)
  }

  form_s <- ~ group_binary + (1|dataset)
  cat("  Formula:", deparse(form_s), "\n")

  v_s <- suppressWarnings(
    voomWithDreamWeights(dge_s, form_s, info_s, BPPARAM = param)
  )
  fit_s <- suppressWarnings(
    dream(v_s, form_s, info_s, BPPARAM = param)
  )

  res_s <- topTable(fit_s, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res_s$gene <- rownames(res_s)

  res_s_dt <- as.data.table(res_s)
  setnames(res_s_dt, "adj.P.Val", "padj")

  n_sig <- nrow(res_s_dt[padj < 0.05 & abs(logFC) > 0.25])
  cat("  DEGs padj < 0.05 & |LFC| > 0.25:", n_sig, "\n")

  out_label <- if (grepl("^M", sex_label)) "male" else "female"
  out_file  <- file.path(rdir, paste0("sex_stratified_results_", out_label, ".csv"))
  fwrite(res_s_dt, out_file)
  cat("  Saved:", basename(out_file), "\n")

  return(res_s_dt)
}

# Detect sex labels from data (could be "M"/"F" or "Male"/"Female")
sex_levels <- levels(info$inferred_sex)
cat("Sex factor levels:", paste(sex_levels, collapse = ", "), "\n")
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]
cat("Using male label:", male_label, ", female label:", female_label, "\n")

# ---- Matched-dataset restriction for fair sex comparison ----
# Only use datasets that have controls in BOTH sexes
ds_m_ctrl <- unique(as.character(info$dataset[info$group_binary == "Control" & info$inferred_sex == male_label]))
ds_f_ctrl <- unique(as.character(info$dataset[info$group_binary == "Control" & info$inferred_sex == female_label]))
matched_ds <- intersect(ds_m_ctrl, ds_f_ctrl)
dropped_ds <- setdiff(levels(info$dataset), matched_ds)

cat("\nMatched-dataset restriction (controls in both sexes):\n")
cat("  Datasets with male controls:", paste(ds_m_ctrl, collapse = ", "), "\n")
cat("  Datasets with female controls:", paste(ds_f_ctrl, collapse = ", "), "\n")
cat("  Matched datasets:", paste(matched_ds, collapse = ", "), "\n")
if (length(dropped_ds) > 0)
  cat("  Dropped (missing controls in one sex):", paste(dropped_ds, collapse = ", "), "\n")

keep_matched <- info$dataset %in% matched_ds
dge_matched <- dge_mega[, keep_matched]
info_matched <- info[keep_matched, , drop = FALSE]
info_matched$dataset <- droplevels(info_matched$dataset)

cat("  Samples after matching:", ncol(dge_matched), "\n")
cat("  Group x sex (matched):\n")
print(table(info_matched$group_binary, info_matched$inferred_sex))

res_male   <- run_stratum_dream(male_label,   dge_matched, info_matched, param, RDIR)
res_female <- run_stratum_dream(female_label, dge_matched, info_matched, param, RDIR)

# ============================================================================
# PART 3 — Interaction-based sex-DEG classification
# ============================================================================
# RATIONALE: The previous classification used significance in sex-stratified
# models (sig in F but not M → "Female_specific"). This is confounded by
# power: with ~101 female vs ~52 male controls, the female stratum has ~2x
# the statistical power. Result: 7,586 "Female_specific" vs 300 "Male_specific"
# reflects sample size, not biology.
#
# NEW APPROACH: The interaction term (group_binary:inferred_sex) from the
# joint dream model directly tests whether the disease effect differs by sex.
# This is the correct statistical test for sex-specificity, unaffected by
# differential power in the stratified runs.
#
# Categories:
#   Concordant     — interaction padj >= 0.05 (disease effect same in both sexes)
#   Female_biased  — interaction padj < 0.05, |disease effect| stronger in females
#   Male_biased    — interaction padj < 0.05, |disease effect| stronger in males
#   Divergent      — interaction padj < 0.05, opposite-direction effects
# ============================================================================
cat("\n===== PART 3: Interaction-based sex-DEG classification =====\n")
cat("Using interaction model as primary arbiter of sex-specificity.\n")
cat("Previous significance-based classification confounded by power asymmetry.\n\n")

PADJ_THRESH      <- 0.05
LFC_THRESH       <- 0.25   # minimum |disease logFC| for stratified significance
LFC_DIFF_THRESH  <- 0.5    # for sex_differential flag
INT_PADJ_THRESH  <- 0.05   # interaction model FDR threshold
DIVERGENT_MIN_LFC <- 0.1   # minimum |logFC| in both sexes to call divergent

if (!is.null(res_male) && !is.null(res_female)) {
  # Join stratified results on gene
  cls_dt <- merge(
    res_male[,   .(gene, logFC_M = logFC, padj_M = padj)],
    res_female[, .(gene, logFC_F = logFC, padj_F = padj)],
    by = "gene", all = TRUE
  )

  # Add interaction model results from PART 1
  cls_dt <- merge(cls_dt,
    res_int_dt[, .(gene, interaction_logFC = logFC, interaction_padj = padj)],
    by = "gene", all.x = TRUE
  )

  # --- Stratified significance (for legacy and diagnostics) ---
  cls_dt[, sig_M := (!is.na(padj_M) & padj_M < PADJ_THRESH &
                     !is.na(logFC_M) & abs(logFC_M) > LFC_THRESH)]
  cls_dt[, sig_F := (!is.na(padj_F) & padj_F < PADJ_THRESH &
                     !is.na(logFC_F) & abs(logFC_F) > LFC_THRESH)]
  cls_dt[, same_sign := !is.na(logFC_M) & !is.na(logFC_F) & sign(logFC_M) == sign(logFC_F)]
  cls_dt[, lfc_diff  := abs(logFC_M - logFC_F)]

  # --- Legacy classification (preserved for comparison) ---
  cls_dt[, sex_class_stratified := fcase(
    sig_M & !sig_F,                          "Male_specific",
    !sig_M & sig_F,                          "Female_specific",
    sig_M & sig_F & same_sign,               "Shared",
    sig_M & sig_F & !same_sign,              "Divergent",
    default = "Not_significant"
  )]

  # --- PRIMARY: Interaction-based classification ---
  cls_dt[, sex_dimorphic := !is.na(interaction_padj) & interaction_padj < INT_PADJ_THRESH]

  cls_dt[, sex_class := fcase(
    # Not sex-dimorphic: no evidence disease effect differs by sex
    !sex_dimorphic, "Concordant",

    # Opposite direction effects (qualitative interaction)
    # Require both effects to be non-trivial to avoid noise-driven sign flips
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      sign(logFC_F) != sign(logFC_M) &
      abs(logFC_F) > DIVERGENT_MIN_LFC & abs(logFC_M) > DIVERGENT_MIN_LFC,
      "Divergent",

    # Same direction (or one near zero): female effect stronger in magnitude
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      abs(logFC_F) >= abs(logFC_M), "Female_biased",

    # Same direction (or one near zero): male effect stronger in magnitude
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      abs(logFC_M) > abs(logFC_F), "Male_biased",

    default = "Unclassified"
  )]

  # --- Sex differential flag (legacy compat: interaction + large LFC diff) ---
  cls_dt[, sex_differential := sex_dimorphic & !is.na(lfc_diff) & (lfc_diff > LFC_DIFF_THRESH)]

  # --- Report new classification ---
  cat("===== INTERACTION-BASED CLASSIFICATION (padj <", INT_PADJ_THRESH, ") =====\n")
  print(cls_dt[, .N, by = sex_class][order(-N)])

  cat("\nSex-dimorphic genes (interaction padj <", INT_PADJ_THRESH, "):",
      sum(cls_dt$sex_dimorphic), "\n")
  cat("  Female_biased:", sum(cls_dt$sex_class == "Female_biased"), "\n")
  cat("  Male_biased:", sum(cls_dt$sex_class == "Male_biased"), "\n")
  cat("  Divergent:", sum(cls_dt$sex_class == "Divergent"), "\n")
  cat("  Unclassified:", sum(cls_dt$sex_class == "Unclassified"), "\n")

  # Diagnostic: median effect sizes by class
  cat("\nMedian |logFC| by new classification:\n")
  for (sc in c("Female_biased", "Male_biased", "Divergent", "Concordant")) {
    sub <- cls_dt[sex_class == sc]
    if (nrow(sub) > 0) {
      cat(sprintf("  %-16s  N=%5d  |logFC_F|=%.3f  |logFC_M|=%.3f  |diff|=%.3f\n",
                  sc, nrow(sub),
                  median(abs(sub$logFC_F), na.rm = TRUE),
                  median(abs(sub$logFC_M), na.rm = TRUE),
                  median(sub$lfc_diff, na.rm = TRUE)))
    }
  }

  # Legacy comparison
  cat("\nLegacy stratified classification (for reference):\n")
  print(cls_dt[, .N, by = sex_class_stratified][order(-N)])

  cat("\nCross-tabulation (rows=legacy, cols=new):\n")
  print(table(cls_dt$sex_class_stratified, cls_dt$sex_class))

  cat("\nSex_differential (interaction + |LFC diff| > 0.5):",
      sum(cls_dt$sex_differential), "\n")

  fwrite(cls_dt, file.path(RDIR, "sex_deg_classification.csv"))
  cat("Saved: sex_deg_classification.csv\n")

  # Use for downstream fgsea and figures
  sex_meta <- cls_dt
} else {
  cat("WARNING: one or both strata failed — classification skipped.\n")
  sex_meta <- NULL
}

# ============================================================================
# PART 4 — Pathway enrichment by sex (fgsea + Hallmark)
# ============================================================================
cat("\n===== PART 4: Pathway enrichment by sex =====\n")

if (!is.null(sex_meta) && !is.null(res_male) && !is.null(res_female)) {
  hallmark    <- msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_dt <- as.data.table(hallmark)
  pathways_list <- split(hallmark_dt$ensembl_gene, hallmark_dt$gs_name)

  # Strip version suffixes from gene IDs for matching
  strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

  # Rank by t-statistic (matching stage GSEA scripts 14c/14f)
  make_ranks <- function(res_dt) {
    res_dt <- res_dt[!is.na(t)]
    res_dt[, gene_base := strip_ver(gene)]
    res_dt <- res_dt[order(-t)]
    ranks <- res_dt$t
    names(ranks) <- res_dt$gene_base
    ranks[!duplicated(names(ranks))]
  }

  ranks_M <- make_ranks(res_male)
  ranks_F <- make_ranks(res_female)

  cat("Male ranked genes:", length(ranks_M), "\n")
  cat("Female ranked genes:", length(ranks_F), "\n")

  set.seed(42)
  fgsea_M <- fgsea(pathways = pathways_list, stats = ranks_M, minSize = 15, maxSize = 500, nPermSimple = 10000)
  fgsea_F <- fgsea(pathways = pathways_list, stats = ranks_F, minSize = 15, maxSize = 500, nPermSimple = 10000)

  fgsea_M_dt <- as.data.table(fgsea_M)
  fgsea_F_dt <- as.data.table(fgsea_F)
  fgsea_M_dt[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
  fgsea_F_dt[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]

  cat("Male fgsea   — significant pathways (padj < 0.1):", nrow(fgsea_M_dt[padj < 0.1]), "\n")
  cat("Female fgsea — significant pathways (padj < 0.1):", nrow(fgsea_F_dt[padj < 0.1]), "\n")

  fwrite(fgsea_M_dt, file.path(RDIR, "sex_fgsea_male.csv"))
  fwrite(fgsea_F_dt, file.path(RDIR, "sex_fgsea_female.csv"))
  cat("Saved: sex_fgsea_male.csv, sex_fgsea_female.csv\n")
} else {
  cat("Skipping fgsea — stratum results unavailable.\n")
  fgsea_M_dt <- NULL
  fgsea_F_dt <- NULL
}

# ============================================================================
# PART 5 — Cross-species sex concordance
# ============================================================================
cat("\n===== PART 5: Cross-species sex concordance =====\n")

if (file.exists(ORTHO_MAP) && file.exists(MOUSE_SEX) && !is.null(sex_meta)) {
  ortho     <- fread(ORTHO_MAP)
  mouse_sex <- fread(MOUSE_SEX)
  cat("Ortholog map:", nrow(ortho), "rows\n")
  cat("Mouse MCD sex-biased genes:", nrow(mouse_sex), "\n")

  mouse_sex[, mouse_gene_base := sub("\\.[0-9]+$", "", gene)]

  mouse_sex_ortho <- merge(
    mouse_sex,
    ortho,
    by.x  = "mouse_gene_base",
    by.y  = "mouse_ensembl_gene_id",
    all.x = TRUE
  )

  human_sex_sub <- sex_meta[, .(
    gene,
    sex_class,
    logFC_M, padj_M,
    logFC_F, padj_F,
    lfc_diff,
    sex_differential
  )]
  human_sex_sub[, gene_base := sub("\\.[0-9]+$", "", gene)]

  xspecies <- merge(
    mouse_sex_ortho,
    human_sex_sub,
    by.x  = "human_ensembl_gene_id",
    by.y  = "gene_base",
    all.x = TRUE
  )

  cat("Mouse sex-biased genes with human sex data:", sum(!is.na(xspecies$sex_class)), "\n")
  cat("  of which classified in human:\n")
  print(xspecies[!is.na(sex_class), .N, by = sex_class][order(-N)])

  fwrite(xspecies, file.path(RDIR, "sex_cross_species_concordance.csv"))
  cat("Saved: sex_cross_species_concordance.csv\n")
} else {
  cat("Skipping cross-species — ortholog map, mouse sex file, or classification unavailable.\n")
}

# ============================================================================
# PART 6 — Figures
# ============================================================================
cat("\n===== PART 6: Figures =====\n")

# --- Figure (a): M vs F logFC scatter colored by interaction-based sex_class ---
if (!is.null(sex_meta)) {
  class_colors <- c(
    "Male_biased"   = "#3B82F6",
    "Female_biased" = "#EF4444",
    "Divergent"     = "#F59E0B",
    "Concordant"    = "#D1D5DB",
    "Unclassified"  = "#D1D5DB"
  )

  p_scatter <- ggplot() +
    geom_point(
      data = sex_meta[sex_class == "Concordant"],
      aes(x = logFC_F, y = logFC_M),
      color = class_colors["Concordant"], alpha = 0.2, size = 0.5
    ) +
    geom_point(
      data = sex_meta[!sex_class %in% c("Concordant", "Unclassified")],
      aes(x = logFC_F, y = logFC_M, color = sex_class),
      alpha = 0.7, size = 1.5
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = class_colors, name = "Sex class\n(interaction FDR<0.05)") +
    labs(
      x     = "Female dream logFC",
      y     = "Male dream logFC",
      title = "Interaction-based sex-dimorphic MASLD DEGs"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 13)
    )

  ggsave(file.path(RDIR, "sex_lfc_scatter.pdf"), p_scatter, width = 8, height = 6)
  cat("Saved: sex_lfc_scatter.pdf\n")
}

# --- Figure (b): Pathway NES dot plot by sex ---
if (!is.null(fgsea_M_dt) && !is.null(fgsea_F_dt)) {
  fgsea_M_plot <- fgsea_M_dt[, .(pathway = gsub("HALLMARK_", "", pathway),
                                  NES, padj, sex = "Male")]
  fgsea_F_plot <- fgsea_F_dt[, .(pathway = gsub("HALLMARK_", "", pathway),
                                  NES, padj, sex = "Female")]
  fgsea_combined <- rbind(fgsea_M_plot, fgsea_F_plot)

  sig_pathways <- fgsea_combined[padj < 0.25, unique(pathway)]
  fgsea_plot   <- fgsea_combined[pathway %in% sig_pathways]

  if (nrow(fgsea_plot) > 0) {
    pathway_order <- fgsea_plot[, .(mean_abs_NES = mean(abs(NES))), by = pathway][
      order(mean_abs_NES)]$pathway
    fgsea_plot[, pathway := factor(pathway, levels = pathway_order)]

    p_pathway <- ggplot(fgsea_plot, aes(x = NES, y = pathway,
                                        size = -log10(padj), color = sex)) +
      geom_point(alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("Male" = "#3B82F6", "Female" = "#EF4444")) +
      scale_size_continuous(name = "-log10(padj)", range = c(1, 6)) +
      labs(
        x     = "Normalized Enrichment Score (NES)",
        y     = NULL,
        title = "Hallmark Pathway Enrichment by Sex"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title  = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 8)
      )

    fig_height <- max(5, length(sig_pathways) * 0.35 + 2)
    ggsave(file.path(RDIR, "sex_pathway_nes_dotplot.pdf"), p_pathway,
           width = 10, height = fig_height)
    cat("Saved: sex_pathway_nes_dotplot.pdf\n")
  } else {
    cat("No pathways significant at padj < 0.25 in either sex; skipping NES dot plot.\n")
  }
}

# ============================================================================
# Summary
# ============================================================================
cat("\n===== STRATEGY 4 COMPLETE =====\n")
cat("Outputs in", RDIR, ":\n")
cat("  sex_interaction_dream.csv           — Dream group x sex interaction\n")
cat("  sex_stratified_results_male.csv     — Male-only dream DE results\n")
cat("  sex_stratified_results_female.csv   — Female-only dream DE results\n")
cat("  sex_deg_classification.csv          — Sex-class labels for all genes\n")
cat("  sex_fgsea_male.csv                  — fgsea Hallmark (male ranked list)\n")
cat("  sex_fgsea_female.csv                — fgsea Hallmark (female ranked list)\n")
cat("  sex_cross_species_concordance.csv   — Mouse MCD sex bias <-> human orthologs\n")
cat("  sex_lfc_scatter.pdf                 — M vs F logFC scatter\n")
cat("  sex_pathway_nes_dotplot.pdf         — Hallmark NES dot plot by sex\n")
cat("Done:", as.character(Sys.time()), "\n")
