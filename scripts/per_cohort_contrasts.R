#!/usr/bin/env Rscript
# per_cohort_contrasts.R
# ---------------------------------------------------------------------------
# Per-cohort DE across multiple contrast families (Phase 7C).
#
# Contrasts:
#   disease_vs_control  : (NAFL + NASH + Borderline) vs Control
#   MASH_vs_control     : NASH vs Control
#   MASL_vs_control     : NAFL vs Control
#   MASH_vs_MASL        : NASH vs NAFL
#   F{1..4}_vs_F0       : fibrosis stage vs F0 baseline
#
# Engine: limma-voom with sex covariate when annotated / inferred; age when
# present and non-constant. Same pattern as 02_per_study_de.R but simplified
# (no config-driven formula; we build contrasts per combo here).
#
# Input:
#   - RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_counts_raw.rds
#   - RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/meta_matched.rds
#   - RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv
#
# Output:
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study_contrasts/{CONTRAST}/{COHORT}_de.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
BASE <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts")
INT  <- file.path(BASE, "analysis/integration")
OUT_ROOT <- file.path(INT, "results/per_study_contrasts")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# ---- Validity table (plan §3) ----
# Rows = cohort, cols = contrast (T/F)
validity <- fread(text = "
cohort,disease_vs_control,MASH_vs_control,MASL_vs_control,MASH_vs_MASL,F_stage
GSE126848,T,T,T,T,F
GSE130970,T,T,T,T,T
GSE135251,T,T,T,T,T
GSE162694,T,T,T,T,T
GSE167523,F,F,F,F,F
GSE174478,F,F,F,T,T
GSE193066,F,F,F,T,T
GSE213621,T,F,F,F,T
GSE240729,F,F,F,F,T
PRJNA512027,T,T,T,T,T
")
for (c in setdiff(names(validity), "cohort")) validity[[c]] <- validity[[c]] == "T"

F_STAGE_CONTRASTS <- c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")
DX_CONTRASTS <- c("disease_vs_control", "MASH_vs_control", "MASL_vs_control", "MASH_vs_MASL")
ALL_CONTRASTS <- c(DX_CONTRASTS, F_STAGE_CONTRASTS)

for (ct in ALL_CONTRASTS) {
  dir.create(file.path(OUT_ROOT, ct), showWarnings = FALSE, recursive = TRUE)
}

# ---- Load data ----
merged <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
for (col in names(meta)) if (is.character(meta[[col]])) meta[get(col) == "", (col) := NA]
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_ids <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing:", length(pass_ids), "/", nrow(meta), "\n")
merged <- merged[, pass_ids]
meta   <- meta[sample_id %in% pass_ids]

# ---- Helper to run one DE comparison ----
# group: factor with 2 levels — level 1 = reference/control, level 2 = case
run_two_group_de <- function(counts_sub, meta_sub, group, cohort, contrast_name) {
  # Drop NA group
  keep_s <- !is.na(group)
  if (sum(keep_s) < 6 || length(unique(group[keep_s])) < 2) {
    cat(sprintf("  SKIP %s / %s (insufficient samples: n=%d, groups=%d)\n",
                cohort, contrast_name, sum(keep_s), length(unique(group[keep_s]))))
    return(invisible(NULL))
  }
  group <- droplevels(factor(group[keep_s]))
  if (length(levels(group)) != 2) {
    cat(sprintf("  SKIP %s / %s (not 2 groups)\n", cohort, contrast_name))
    return(invisible(NULL))
  }
  # Per-group counts
  tab <- table(group)
  if (any(tab < 2)) {
    cat(sprintf("  SKIP %s / %s (group size <2: %s)\n", cohort, contrast_name,
                paste(names(tab), tab, sep="=", collapse=",")))
    return(invisible(NULL))
  }
  counts_sub <- counts_sub[, meta_sub$sample_id[keep_s]]
  meta_sub <- meta_sub[keep_s, ]
  meta_sub$.grp <- group

  # Covariates: sex (annotated or inferred), age (only if non-constant w/o NA)
  covars <- character(0)
  if ("sex" %in% names(meta_sub) && sum(!is.na(meta_sub$sex)) == nrow(meta_sub) &&
      length(unique(meta_sub$sex)) > 1) {
    meta_sub$sex <- factor(meta_sub$sex); covars <- c(covars, "sex")
  } else if ("inferred_sex" %in% names(meta_sub) &&
             sum(!is.na(meta_sub$inferred_sex)) == nrow(meta_sub) &&
             length(unique(meta_sub$inferred_sex)) > 1) {
    meta_sub$inferred_sex <- factor(meta_sub$inferred_sex); covars <- c(covars, "inferred_sex")
  }
  if ("age" %in% names(meta_sub)) {
    a <- suppressWarnings(as.numeric(meta_sub$age))
    if (!all(is.na(a)) && sum(is.na(a)) == 0 && length(unique(a)) > 2) {
      meta_sub$age <- a; covars <- c(covars, "age")
    }
  }

  form_rhs <- paste(c("0 + .grp", covars), collapse = " + ")
  design <- model.matrix(as.formula(paste("~", form_rhs)), data = meta_sub)
  # Contrast: level2 - level1
  lv <- levels(meta_sub$.grp)
  cname <- paste0(".grp", lv[2], " - ", ".grp", lv[1])
  contr <- tryCatch(makeContrasts(contrasts = cname, levels = design),
                    error = function(e) NULL)
  if (is.null(contr)) {
    cat(sprintf("  SKIP %s / %s (contrast fail)\n", cohort, contrast_name))
    return(invisible(NULL))
  }

  dge <- DGEList(counts = counts_sub)
  keep_g <- filterByExpr(dge, design = design)
  dge <- dge[keep_g, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  dt <- data.table(
    gene    = res$gene,
    logFC   = res$logFC,
    t       = res$t,
    P.Value = res$P.Value,
    padj    = res$adj.P.Val,
    AveExpr = res$AveExpr
  )
  outfile <- file.path(OUT_ROOT, contrast_name, paste0(cohort, "_de.csv"))
  fwrite(dt, outfile)
  n_sig <- sum(dt$padj < 0.1, na.rm = TRUE)
  cat(sprintf("  OK   %s / %s  n=%d  genes=%d  padj<0.1=%d  -> %s\n",
              cohort, contrast_name, nrow(meta_sub), nrow(dt), n_sig, outfile))
}

# ---- Iterate cohorts ----
cohorts <- intersect(unique(meta$dataset), validity$cohort)

for (co in cohorts) {
  cat("\n============================\n  ", co, "\n============================\n")
  msk <- meta$dataset == co
  meta_sub <- as.data.frame(meta[msk])
  counts_sub <- merged[, meta[msk, sample_id]]
  v_row <- validity[cohort == co]

  dx <- meta_sub$diagnosis_harmonized
  # Treat Borderline as disease (NASH-ish, NAS 3-4)
  disease_any <- ifelse(is.na(dx), NA, ifelse(dx %in% c("NAFL","NASH","Borderline"), "Disease",
                                              ifelse(dx == "Control", "Control", NA)))

  if (isTRUE(v_row$disease_vs_control)) {
    g <- factor(disease_any, levels = c("Control", "Disease"))
    run_two_group_de(counts_sub, meta_sub, g, co, "disease_vs_control")
  } else {
    cat(sprintf("  SKIP %s / disease_vs_control (gated)\n", co))
  }

  if (isTRUE(v_row$MASH_vs_control)) {
    g <- factor(ifelse(dx == "NASH", "NASH", ifelse(dx == "Control", "Control", NA)),
                levels = c("Control", "NASH"))
    run_two_group_de(counts_sub, meta_sub, g, co, "MASH_vs_control")
  } else {
    cat(sprintf("  SKIP %s / MASH_vs_control (gated)\n", co))
  }

  if (isTRUE(v_row$MASL_vs_control)) {
    g <- factor(ifelse(dx == "NAFL", "NAFL", ifelse(dx == "Control", "Control", NA)),
                levels = c("Control", "NAFL"))
    run_two_group_de(counts_sub, meta_sub, g, co, "MASL_vs_control")
  } else {
    cat(sprintf("  SKIP %s / MASL_vs_control (gated)\n", co))
  }

  if (isTRUE(v_row$MASH_vs_MASL)) {
    g <- factor(ifelse(dx == "NASH", "NASH", ifelse(dx == "NAFL", "NAFL", NA)),
                levels = c("NAFL", "NASH"))
    run_two_group_de(counts_sub, meta_sub, g, co, "MASH_vs_MASL")
  } else {
    cat(sprintf("  SKIP %s / MASH_vs_MASL (gated)\n", co))
  }

  if (isTRUE(v_row$F_stage)) {
    fs <- suppressWarnings(as.integer(meta_sub$fibrosis_stage))
    for (k in 1:4) {
      contrast_name <- paste0("F", k, "_vs_F0")
      g <- factor(ifelse(fs == 0, "F0", ifelse(fs == k, paste0("F", k), NA)),
                  levels = c("F0", paste0("F", k)))
      run_two_group_de(counts_sub, meta_sub, g, co, contrast_name)
    }
  } else {
    for (k in 1:4) cat(sprintf("  SKIP %s / F%d_vs_F0 (gated)\n", co, k))
  }
}

cat("\nDone.\n")
