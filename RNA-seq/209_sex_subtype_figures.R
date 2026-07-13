#!/usr/bin/env Rscript
# 209_sex_subtype_figures.R — Sex-subtype synthesis figures for stratified causal architecture
#
# Generates 4 supplementary panels:
#   Panel A: Enrichment barplot — COLOC OR for female vs male vs divergent vs S2-specific DEGs (95% CI)
#   Panel B: S2 genetic risk score ROC curve (if s2_genetic_risk.csv available)
#   Panel C: Scatter — sex-causal score (x) vs S2-specificity (y)
#   Panel D: Top 15 sex/subtype-specific causal genes (dot plot: PIP x sex bias)
#
# Depends on outputs from Scripts 207 + 208:
#   - RNA-seq/results/stratified_causal/sex_coloc_enrichment.csv
#   - RNA-seq/results/stratified_causal/subtype_coloc_enrichment.csv
#   - RNA-seq/results/stratified_causal/s2_genetic_risk.csv
#   - RNA-seq/results/stratified_causal/sex_causal_scores.csv
#   - RNA-seq/results/stratified_causal/subtype_causal_scores.csv
#
# If upstream outputs missing, computes panels directly from primary data:
#   - GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   - RNA-seq/Human/.../results/integration/sex_deg_classification.csv
#   - RNA-seq/results/subtypes/subtype_markers.csv
#   - RNA-seq/results/subtypes/nmf_assignments.csv
#
# Output: figures/supplementary/figS_sex_dimorphism/figS_sex_subtype_*.pdf
# SLURM: cpu partition, 4 CPUs, 32GB, 48h

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
STRAT_DIR   <- file.path(BASE, "RNA-seq/results/stratified_causal")
COLOC_FILE  <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
SEX_FILE    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                         "results/integration/sex_deg_classification.csv")
MARKERS_FILE <- file.path(BASE, "RNA-seq/results/subtypes/subtype_markers.csv")
NMF_FILE    <- file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv")
ATLAS_FILE  <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

FIG_DIR <- file.path(BASE, "figures/supplementary/figS_sex_dimorphism")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Check for pre-computed Script 207/208 outputs
# ---------------------------------------------------------------------------
precomputed_files <- c(
  sex_enrich    = file.path(STRAT_DIR, "sex_coloc_enrichment.csv"),
  sub_enrich    = file.path(STRAT_DIR, "subtype_coloc_enrichment.csv"),
  s2_risk       = file.path(STRAT_DIR, "s2_genetic_risk.csv"),
  sex_causal    = file.path(STRAT_DIR, "sex_causal_scores.csv"),
  sub_causal    = file.path(STRAT_DIR, "subtype_causal_scores.csv")
)

have_precomputed <- all(file.exists(precomputed_files))

if (have_precomputed) {
  cat("Found all pre-computed 207/208 outputs. Loading directly.\n")
} else {
  missing <- names(precomputed_files)[!file.exists(precomputed_files)]
  cat("Missing pre-computed files:", paste(missing, collapse = ", "), "\n")
  cat("Will compute panels from primary data.\n")
}

# ---------------------------------------------------------------------------
# Load primary data (always needed for fallback or supplementation)
# ---------------------------------------------------------------------------
# COLOC gene-level
if (!file.exists(COLOC_FILE)) {
  cat("CRITICAL: COLOC gene-level file not found:", COLOC_FILE, "\n")
  cat("Waiting for Script 207/208 or COLOC pipeline...\n")
  quit(save = "no", status = 0)
}
coloc <- fread(COLOC_FILE)
# Clean: drop rows with empty gene name
coloc <- coloc[gene != "" & !is.na(gene)]
cat("COLOC genes loaded:", nrow(coloc), "\n")

# Sex DEG classification
if (!file.exists(SEX_FILE)) {
  cat("CRITICAL: Sex DEG classification not found:", SEX_FILE, "\n")
  cat("Waiting for Script 207/208...\n")
  quit(save = "no", status = 0)
}
sex_dt <- fread(SEX_FILE)
cat("Sex classification loaded:", nrow(sex_dt), "genes\n")
cat("  Classes:", paste(names(table(sex_dt$sex_class)), collapse = ", "), "\n")

# Subtype markers
if (!file.exists(MARKERS_FILE)) {
  cat("WARNING: Subtype markers not found:", MARKERS_FILE, "\n")
  markers_dt <- NULL
} else {
  markers_dt <- fread(MARKERS_FILE)
  cat("Subtype markers loaded:", nrow(markers_dt), "genes\n")
}

# Multi-evidence atlas (for symbol mapping)
gene_map <- NULL
if (file.exists(ATLAS_FILE)) {
  atlas <- fread(ATLAS_FILE, select = c("ensembl_id", "human_symbol"))
  atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
  gene_map <- unique(atlas[human_symbol != "" & !is.na(human_symbol),
                           .(ensembl_clean, symbol = human_symbol)])
  cat("Gene map loaded:", nrow(gene_map), "mappings\n")
}

# Helper: map ensembl -> symbol
add_symbol <- function(dt, gene_col = "gene") {
  if ("symbol" %in% names(dt)) return(dt)
  if (is.null(gene_map)) {
    dt[, symbol := get(gene_col)]
    return(dt)
  }
  dt[, ensembl_clean := sub("\\..*", "", get(gene_col))]
  dt <- merge(dt, gene_map, by = "ensembl_clean", all.x = TRUE)
  dt[is.na(symbol), symbol := ensembl_clean]
  dt[, ensembl_clean := NULL]
  dt
}

# ===========================================================================
# PANEL A: Enrichment barplot — COLOC OR per sex/subtype category
# ===========================================================================
cat("\n=== Panel A: COLOC enrichment by sex/subtype class ===\n")

if (have_precomputed && file.exists(precomputed_files["sex_enrich"]) &&
    file.exists(precomputed_files["sub_enrich"])) {
  # Load pre-computed enrichment
  sex_enrich <- fread(precomputed_files["sex_enrich"])
  sub_enrich <- fread(precomputed_files["sub_enrich"])
  enrich_dt <- rbind(sex_enrich, sub_enrich, fill = TRUE)
} else {
  # Compute Fisher's exact enrichment from primary data
  cat("Computing enrichment from primary data...\n")

  # Build universe: genes present in both COLOC and sex classification
  sex_dt[, ensembl_clean := sub("\\..*", "", gene)]
  coloc[, ensembl_clean := sub("\\..*", "", ensembl)]

  # COLOC positive: PP.H4 > 0.5
  coloc_genes <- coloc[coloc_best_pp4 > 0.5, ensembl_clean]
  cat("  COLOC PP.H4>0.5 genes:", length(coloc_genes), "\n")

  universe <- unique(c(sex_dt$ensembl_clean, coloc$ensembl_clean))
  cat("  Universe size:", length(universe), "\n")

  # Function for Fisher's exact
  fisher_enrich <- function(test_genes, coloc_pos, univ, label) {
    test_set <- unique(test_genes)
    a <- length(intersect(test_set, coloc_pos))
    b <- length(setdiff(coloc_pos, test_set))
    c_val <- length(setdiff(test_set, coloc_pos))
    d <- length(setdiff(univ, union(test_set, coloc_pos)))
    mat <- matrix(c(a, b, c_val, d), nrow = 2)
    ft <- fisher.test(mat)
    data.table(
      category = label,
      n_genes = length(test_set),
      n_coloc_overlap = a,
      odds_ratio = ft$estimate,
      ci_low = ft$conf.int[1],
      ci_high = ft$conf.int[2],
      pvalue = ft$p.value
    )
  }

  # Sex categories
  enrich_list <- list()
  # Auto-detect v2 (interaction-based) vs v1 (stratified) class names
  sex_classes_209 <- unique(sex_dt$sex_class)
  use_v2_209 <- "Female_biased" %in% sex_classes_209
  female_cls <- if (use_v2_209) "Female_biased" else "Female_specific"
  male_cls   <- if (use_v2_209) "Male_biased"   else "Male_specific"
  for (cls in c(female_cls, male_cls, "Divergent")) {
    test_g <- sex_dt[sex_class == cls, ensembl_clean]
    if (length(test_g) > 0) {
      enrich_list[[cls]] <- fisher_enrich(test_g, coloc_genes, universe, cls)
    }
  }

  # Subtype: k-program refactor (2026-04-20).
  # subtype_markers.csv now has `program_code` + `program` columns from the
  # k-program NMF refactor. Emit per-program enrichment plus legacy
  # S1_upregulated / S2_upregulated aliases where S2 = Fibrotic program
  # and S1 = dominant program of samples where Fibrotic is NOT dominant.
  if (!is.null(markers_dt)) {
    markers_dt[, ensembl_clean := sub("\\..*", "", gene)]
    prog_lbls_path <- file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv")
    prog_lbls <- if (file.exists(prog_lbls_path)) fread(prog_lbls_path) else NULL
    if ("program_code" %in% names(markers_dt) && !is.null(prog_lbls)) {
      for (j in seq_len(nrow(prog_lbls))) {
        pc <- prog_lbls$program_code[j]
        lbl <- prog_lbls$biological_label[j]
        pg <- markers_dt[program_code == pc & direction == "up" & padj < 0.05,
                         unique(ensembl_clean)]
        if (length(pg) > 0) {
          cat_name <- paste0(gsub("[^A-Za-z0-9]+", "_", lbl), "_upregulated")
          enrich_list[[cat_name]] <- fisher_enrich(pg, coloc_genes, universe, cat_name)
        }
      }
      # Legacy aliases: S2 = Fibrotic program, S1 = all other programs collapsed
      fib_pc <- prog_lbls$program_code[prog_lbls$program_code == "P3"]  # was: label_category == "Fibrogenic"
      if (length(fib_pc) >= 1) {
        s2_up <- markers_dt[program_code == fib_pc[1] & direction == "up" & padj < 0.05,
                            unique(ensembl_clean)]
        s1_up <- markers_dt[program_code != fib_pc[1] & direction == "up" & padj < 0.05,
                            unique(ensembl_clean)]
      } else {
        s2_up <- character(0); s1_up <- character(0)
      }
    } else {
      # Legacy path (pre-refactor marker file)
      s2_up <- if ("subtype" %in% names(markers_dt))
        markers_dt[subtype == "S2" & direction == "up", unique(ensembl_clean)] else character(0)
      s1_up <- if ("subtype" %in% names(markers_dt))
        markers_dt[subtype == "S1" & direction == "up", unique(ensembl_clean)] else character(0)
    }
    if (length(s2_up) > 0) {
      enrich_list[["S2_upregulated"]] <- fisher_enrich(s2_up, coloc_genes, universe, "S2_upregulated")
    }
    if (length(s1_up) > 0) {
      enrich_list[["S1_upregulated"]] <- fisher_enrich(s1_up, coloc_genes, universe, "S1_upregulated")
    }
  }

  enrich_dt <- rbindlist(enrich_list)
  enrich_dt[, padj := p.adjust(pvalue, method = "BH")]

  cat("Enrichment results:\n")
  print(enrich_dt[, .(category, n_genes, n_coloc_overlap, odds_ratio = round(odds_ratio, 2),
                       pvalue = signif(pvalue, 3))])

  # Save intermediate
  dir.create(STRAT_DIR, recursive = TRUE, showWarnings = FALSE)
  fwrite(enrich_dt, file.path(STRAT_DIR, "sex_subtype_coloc_enrichment_209.csv"))
}

# Build Panel A
if (nrow(enrich_dt) == 0) {
  cat("WARNING: No enrichment results — all sex classes had 0 overlap with COLOC genes.\n")
  cat("  Creating placeholder Panel A.\n")
  enrich_dt <- data.table(category = "No data", n_genes = 0, n_coloc_overlap = 0,
                           odds_ratio = 1, ci_low = 1, ci_high = 1, pvalue = 1, padj = 1)
}
enrich_plot <- copy(enrich_dt)
# Harmonize CI column names (207 uses ci_lower/ci_upper, local Fisher uses ci_low/ci_high)
if ("ci_upper" %in% names(enrich_plot) && !"ci_high" %in% names(enrich_plot)) {
  setnames(enrich_plot, c("ci_lower", "ci_upper"), c("ci_low", "ci_high"))
}
# Harmonize category column (207 uses sex_class, local Fisher uses category)
if ("sex_class" %in% names(enrich_plot) && !"category" %in% names(enrich_plot)) {
  setnames(enrich_plot, "sex_class", "category")
}
# Cap extreme CIs for display
enrich_plot[, ci_high_cap := pmin(ci_high, 10)]
# Significance stars
enrich_plot[, sig_label := fifelse(pvalue < 0.001, "***",
                            fifelse(pvalue < 0.01, "**",
                              fifelse(pvalue < 0.05, "*", "ns")))]

# Category display names and colors
# Build labels/colors dynamically for v1 or v2 class names
cat_labels <- c(
  Female_biased   = "Female-biased",
  Female_specific = "Female-specific",
  Male_biased     = "Male-biased",
  Male_specific   = "Male-specific",
  Divergent       = "Sex-divergent",
  S2_upregulated  = "S2-upregulated",
  S1_upregulated  = "S1-upregulated"
)
cat_colors <- c(
  Female_biased   = masld_colors$female,
  Female_specific = masld_colors$female,
  Male_biased     = masld_colors$male,
  Male_specific   = masld_colors$male,
  Divergent       = "#E91E63",
  S2_upregulated  = "#C2185B",
  S1_upregulated  = "#42A5F5"
)

enrich_plot[, category_label := cat_labels[category]]
enrich_plot[, category_label := factor(category_label,
  levels = rev(cat_labels[cat_labels %in% enrich_plot$category_label]))]

p_a <- ggplot(enrich_plot, aes(x = odds_ratio, y = category_label, color = category)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high_cap), height = 0.2, linewidth = 0.4) +
  geom_text(aes(x = ci_high_cap + 0.15, label = sig_label), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = cat_colors, guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  labs(
    title = "COLOC enrichment by sex & subtype class",
    x = "Odds ratio (95% CI)",
    y = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold"))

cat("Panel A: done\n")

# ===========================================================================
# PANEL B: S2 genetic risk score ROC curve
# ===========================================================================
cat("\n=== Panel B: S2 genetic risk ROC ===\n")

p_b <- placeholder("Panel b: S2 genetic risk ROC\n(requires Script 208 output)")

if (have_precomputed && file.exists(precomputed_files["s2_risk"])) {
  s2_risk <- fread(precomputed_files["s2_risk"])
  # Expect columns: model, fpr, tpr, auroc (or similar)
  if (all(c("fpr", "tpr") %in% names(s2_risk))) {
    # If multiple models, use 'model' column
    if ("model" %in% names(s2_risk)) {
      auroc_labels <- s2_risk[, .(auroc = max(auroc, na.rm = TRUE)), by = model]
      auroc_labels[, label := paste0(model, " (AUC=", sprintf("%.3f", auroc), ")")]
      p_b <- ggplot(s2_risk, aes(x = fpr, y = tpr, color = model)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
        geom_line(linewidth = 0.6) +
        scale_color_manual(
          values = c("COLOC-only" = masld_colors$twas,
                     "COLOC+TWAS" = masld_colors$up,
                     "Full transcriptomic" = masld_colors$conserved),
          labels = auroc_labels$label
        ) +
        labs(
          title = "S2 subtype prediction from causal features",
          x = "False positive rate",
          y = "True positive rate",
          color = NULL
        ) +
        theme_masld(base_size = 7) +
        theme(
          legend.position = c(0.65, 0.25),
          legend.background = element_blank(),
          plot.title = element_text(size = 8, face = "bold")
        )
    } else {
      # Single model
      auroc_val <- if ("auroc" %in% names(s2_risk)) s2_risk$auroc[1] else NA
      p_b <- ggplot(s2_risk, aes(x = fpr, y = tpr)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
        geom_line(linewidth = 0.6, color = masld_colors$up) +
        annotate("text", x = 0.6, y = 0.3,
                 label = paste0("AUROC = ", sprintf("%.3f", auroc_val)),
                 size = 2.5, color = masld_colors$up) +
        labs(
          title = "S2 subtype prediction from causal features",
          x = "False positive rate",
          y = "True positive rate"
        ) +
        theme_masld(base_size = 7) +
        theme(plot.title = element_text(size = 8, face = "bold"))
    }
    cat("Panel B: loaded from s2_genetic_risk.csv\n")
  } else {
    cat("Panel B: s2_genetic_risk.csv lacks expected columns; using placeholder\n")
  }
} else {
  # Fallback: compute a simple genetic risk ROC from COLOC + NMF
  cat("No pre-computed S2 risk data. Attempting fallback from primary data...\n")

  if (file.exists(NMF_FILE) && nrow(coloc) > 0) {
    nmf <- fread(NMF_FILE)
    nmf[, is_s2 := as.integer(nmf_subtype == "S2")]

    # Load atlas for per-sample expression of COLOC genes (not available here)
    # Instead, show gene-level S2 enrichment as a proxy
    cat("  Cannot compute per-sample ROC without expression matrix. Using placeholder.\n")
  }
}

# ===========================================================================
# PANEL C: Scatter — sex-causal score vs S2-specificity
# ===========================================================================
cat("\n=== Panel C: Sex-causal score vs S2-specificity scatter ===\n")

if (have_precomputed && file.exists(precomputed_files["sex_causal"]) &&
    file.exists(precomputed_files["sub_causal"]) &&
    file.info(precomputed_files["sub_causal"])$size > 0 &&
    nrow(fread(precomputed_files["sub_causal"], nrows = 1)) > 0) {
  sex_scores <- fread(precomputed_files["sex_causal"])
  sub_scores <- fread(precomputed_files["sub_causal"])
  # Merge by gene
  # Find a column present in BOTH sex_scores and sub_scores
  merge_col <- intersect(
    intersect(c("human_symbol", "gene", "ensembl_id", "ens_base", "symbol"), names(sex_scores)),
    names(sub_scores)
  )[1]
  if (is.na(merge_col)) {
    # Fallback: use human_symbol in both (rename if needed)
    if ("human_symbol" %in% names(sex_scores) && "human_symbol" %in% names(sub_scores)) {
      merge_col <- "human_symbol"
    } else {
      cat("WARNING: No common merge column between sex_scores and sub_scores. Skipping Panel C.\n")
      merge_col <- NULL
    }
  }
  scatter_dt <- if (!is.null(merge_col)) {
    cat("  Merging sex_scores and sub_scores by:", merge_col, "\n")
    # Dedup: subtype_causal_scores.csv (k-program refactor) has one row per
    # (gene, program_code) — need to collapse to one row per gene before
    # the cross join with sex_scores. Also coerce merge-key to character
    # (all-NA columns get typed as logical by fread and break bmerge).
    sex_scores_u <- unique(as.data.table(sex_scores), by = merge_col)
    sub_scores_u <- unique(as.data.table(sub_scores), by = merge_col)
    sex_scores_u[[merge_col]] <- as.character(sex_scores_u[[merge_col]])
    sub_scores_u[[merge_col]] <- as.character(sub_scores_u[[merge_col]])
    dt <- merge(sex_scores_u, sub_scores_u, by = merge_col, suffixes = c(".sex", ".sub"))
    # Reconcile suffixed columns after merge
    # For columns duplicated in both: prefer .sex version, then .sub
    for (col in c("coloc_best_pp4", "coloc_best_gwas", "sex_class", "lfc_diff",
                  "sex_causal_score", "sex_causal_signed")) {
      sex_col <- paste0(col, ".sex")
      if (sex_col %in% names(dt) && !col %in% names(dt)) {
        setnames(dt, sex_col, col)
      }
    }
    for (col in c("s2_specificity", "cohens_d", "mean_S1", "mean_S2", "mean_diff",
                  "div_pval", "div_padj", "s2_causal_score")) {
      sub_col <- paste0(col, ".sub")
      if (sub_col %in% names(dt) && !col %in% names(dt)) {
        setnames(dt, sub_col, col)
      }
    }
    # Derive s2_specificity from mean_diff if not present (pre-computed path)
    if (!"s2_specificity" %in% names(dt)) {
      if ("mean_diff" %in% names(dt)) {
        dt[, s2_specificity := mean_diff]
      } else if ("cohens_d" %in% names(dt)) {
        dt[, s2_specificity := cohens_d]
      } else {
        dt[, s2_specificity := 0]
      }
    }
    dt
  } else {
    data.table()
  }
} else {
  # Compute from primary data
  cat("Computing sex-causal and S2-causal scores from primary data...\n")

  # Sex-causal score: PP4 x |female_logFC - male_logFC|
  sex_dt[, ensembl_clean := sub("\\..*", "", gene)]
  coloc[, ensembl_clean := sub("\\..*", "", ensembl)]

  sc <- merge(sex_dt[, .(ensembl_clean, lfc_diff, sex_class, logFC_F, logFC_M)],
              coloc[, .(ensembl_clean, coloc_best_pp4)],
              by = "ensembl_clean")
  sc[, sex_causal_score := coloc_best_pp4 * abs(lfc_diff)]

  # S2-specificity: use marker logFC for S2 (or |S2_logFC| from markers file)
  if (!is.null(markers_dt)) {
    markers_dt[, ensembl_clean := sub("\\..*", "", gene)]
    s2_markers <- markers_dt[subtype == "S2", .(ensembl_clean, s2_logFC = logFC, s2_padj = adj.P.Val)]
    s1_markers <- markers_dt[subtype == "S1", .(ensembl_clean, s1_logFC = logFC, s1_padj = adj.P.Val)]

    sc <- merge(sc, s2_markers, by = "ensembl_clean", all.x = TRUE)
    sc <- merge(sc, s1_markers, by = "ensembl_clean", all.x = TRUE)
    # S2-specificity: positive = S2-up, negative = S1-up
    sc[, s2_specificity := fifelse(!is.na(s2_logFC), s2_logFC,
                            fifelse(!is.na(s1_logFC), -s1_logFC, 0))]
  } else {
    sc[, s2_specificity := 0]
  }

  sc <- add_symbol(sc, "ensembl_clean")
  scatter_dt <- sc

  # Save intermediate
  fwrite(scatter_dt, file.path(STRAT_DIR, "sex_s2_scatter_209.csv"))
}

# Ensure symbol column exists for labeling
if (!"symbol" %in% names(scatter_dt) && "human_symbol" %in% names(scatter_dt)) {
  setnames(scatter_dt, "human_symbol", "symbol")
} else if (!"symbol" %in% names(scatter_dt)) {
  scatter_dt <- add_symbol(scatter_dt, names(scatter_dt)[1])
}

# Determine x/y columns
if ("sex_causal_score" %in% names(scatter_dt)) {
  x_col <- "sex_causal_score"
} else {
  x_col <- grep("sex.*score|sex_causal", names(scatter_dt), value = TRUE, ignore.case = TRUE)[1]
}
if ("s2_specificity" %in% names(scatter_dt)) {
  y_col <- "s2_specificity"
} else {
  y_col <- grep("s2.*score|s2.*spec|subtype.*score", names(scatter_dt), value = TRUE, ignore.case = TRUE)[1]
}

if (!is.null(x_col) && !is.null(y_col) && !is.na(x_col) && !is.na(y_col)) {
  # Identify genes to label: top by combined score
  scatter_dt[, combined := abs(get(x_col)) + abs(get(y_col))]
  # Highlight: both female-biased AND S2-specific AND genetically causal
  pp4_col <- grep("pp4|coloc", names(scatter_dt), value = TRUE, ignore.case = TRUE)[1]
  if (!is.null(pp4_col) && !is.na(pp4_col)) {
    scatter_dt[, highlight := get(x_col) > quantile(get(x_col), 0.9, na.rm = TRUE) &
                              get(y_col) > quantile(get(y_col), 0.8, na.rm = TRUE) &
                              get(pp4_col) > 0.5]
  } else {
    scatter_dt[, highlight := combined > quantile(combined, 0.95, na.rm = TRUE)]
  }

  # Top 15 for labeling
  top_label <- scatter_dt[order(-combined)][1:min(15, .N)]

  p_c <- ggplot(scatter_dt, aes(x = get(x_col), y = get(y_col))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    geom_point(aes(color = highlight), size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns), guide = "none") +
    geom_text_repel(
      data = top_label,
      aes(label = symbol),
      size = 2, max.overlaps = 15, segment.size = 0.2, min.segment.length = 0,
      color = "gray20"
    ) +
    labs(
      title = "Sex-causal score vs S2-specificity",
      x = "Sex-causal score (PP.H4 x |female - male| LFC)",
      y = "S2-specificity (logFC S2 vs S1)"
    ) +
    theme_masld(base_size = 7) +
    theme(plot.title = element_text(size = 8, face = "bold"))

  cat("Panel C: done\n")
} else {
  cat("Panel C: could not identify score columns. Using placeholder.\n")
  p_c <- placeholder("Panel c: Sex-causal vs S2-specificity\n(missing score columns)")
}

# ===========================================================================
# PANEL D: Top 15 sex/subtype-specific causal genes (dot plot)
# ===========================================================================
cat("\n=== Panel D: Top 15 causal genes dot plot ===\n")

if (exists("scatter_dt") && nrow(scatter_dt) > 0 &&
    any(c("sex_causal_score", "s2_causal_score") %in% names(scatter_dt))) {
  # Use scatter_dt which has COLOC PP4, sex class, S2 specificity
  dot_dt <- copy(scatter_dt)

  # Ensure required columns exist (may be missing in pre-computed merge path)
  if (!"sex_causal_score" %in% names(dot_dt) && "coloc_best_pp4" %in% names(dot_dt) &&
      "lfc_diff" %in% names(dot_dt)) {
    dot_dt[, sex_causal_score := coloc_best_pp4 * abs(lfc_diff)]
  }
  if (!"combined" %in% names(dot_dt)) {
    sc1 <- if ("sex_causal_score" %in% names(dot_dt)) dot_dt$sex_causal_score else 0
    sc2 <- if ("s2_causal_score" %in% names(dot_dt)) dot_dt$s2_causal_score else 0
    dot_dt[, combined := nafill(sc1, fill = 0) + nafill(sc2, fill = 0)]
  }
  if (!"symbol" %in% names(dot_dt) && "human_symbol" %in% names(dot_dt)) {
    setnames(dot_dt, "human_symbol", "symbol")
  } else if (!"symbol" %in% names(dot_dt)) {
    dot_dt <- add_symbol(dot_dt, names(dot_dt)[1])
  }

  dot_dt <- dot_dt[!is.na(coloc_best_pp4) & coloc_best_pp4 > 0.3]

  # Classification for coloring
  if ("sex_class" %in% names(dot_dt)) {
    # Support both v1 (Female_specific) and v2 (Female_biased) class names
    fem_vals <- c("Female_specific", "Female_biased")
    mal_vals <- c("Male_specific", "Male_biased")
    dot_dt[, class := fifelse(sex_class %in% c(fem_vals, "Divergent") & s2_specificity > 0,
                              "Female + S2",
                       fifelse(sex_class %in% fem_vals, "Female-biased",
                        fifelse(sex_class %in% mal_vals, "Male-biased",
                         fifelse(s2_specificity > 0.2, "S2-specific", "Other"))))]
  } else {
    dot_dt[, class := "Other"]
  }

  # Select top 15 by combined sex-causal + S2 score, requiring COLOC PP4 > 0.5
  dot_top <- dot_dt[coloc_best_pp4 > 0.5][order(-combined)][1:min(15, .N)]

  if (nrow(dot_top) > 0) {
    dot_top <- dot_top[!duplicated(symbol)]  # guard against dup symbols
    dot_top[, symbol := factor(symbol, levels = rev(unique(symbol)))]

    class_colors <- c(
      "Female + S2"     = "#880E4F",
      "Female-specific" = masld_colors$female,
      "Male-specific"   = masld_colors$male,
      "S2-specific"     = masld_colors$up,
      "Other"           = "gray50"
    )

    p_d <- ggplot(dot_top, aes(x = coloc_best_pp4, y = symbol)) +
      geom_point(aes(size = abs(lfc_diff), color = class)) +
      scale_color_manual(values = class_colors, name = "Category") +
      scale_size_continuous(range = c(1.5, 5), name = "|Sex LFC diff|") +
      scale_x_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)) +
      labs(
        title = "Top causal genes by sex/subtype specificity",
        x = "COLOC PP.H4",
        y = NULL
      ) +
      theme_masld(base_size = 7) +
      theme(
        plot.title = element_text(size = 8, face = "bold"),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm")
      )

    cat("Panel D: done (", nrow(dot_top), "genes plotted)\n")
  } else {
    cat("Panel D: No genes with PP.H4>0.5 and combined score. Using placeholder.\n")
    p_d <- placeholder("Panel d: Top causal genes\n(insufficient data)")
  }
} else {
  cat("Panel D: scatter_dt unavailable. Using placeholder.\n")
  p_d <- placeholder("Panel d: Top 15 sex/subtype-specific\ncausal genes")
}

# ===========================================================================
# Assemble and save
# ===========================================================================
cat("\n=== Assembling combined figure ===\n")

combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 9, face = "bold")
    )
  )

# Save combined figure
save_fig(combined,
         file.path(FIG_DIR, "figS_sex_subtype_combined.pdf"),
         width = fig_full_width, height = 5.5)
cat("Saved:", file.path(FIG_DIR, "figS_sex_subtype_combined.pdf"), "\n")

# Save individual panels
save_fig(p_a, file.path(FIG_DIR, "figS_sex_subtype_enrichment.pdf"),
         width = fig_half_width, height = 2.5)
save_fig(p_b, file.path(FIG_DIR, "figS_sex_subtype_risk.pdf"),
         width = fig_half_width, height = 2.5)
save_fig(p_c, file.path(FIG_DIR, "figS_sex_subtype_scatter.pdf"),
         width = fig_half_width, height = 2.5)
save_fig(p_d, file.path(FIG_DIR, "figS_sex_subtype_dotplot.pdf"),
         width = fig_half_width, height = 3)

cat("\nAll panels saved to:", FIG_DIR, "\n")
cat("Done.\n")
