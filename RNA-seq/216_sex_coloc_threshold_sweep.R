#!/usr/bin/env Rscript
# 216_sex_coloc_threshold_sweep.R
# ---------------------------------------------------------------------------
# Threshold sweep of sex × COLOC enrichment (Team B1, task A2).
#
# For each PP.H4 threshold in {0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9}, per GWAS
# family (enzyme / imaging / steatosis-disease / cirrhosis / HCC) × per
# ancestry, run Fisher's exact test of:
#   - sex class membership (Female_biased / Male_biased / Divergent /
#     Concordant) among DEGs
#   - whether the gene reaches the COLOC threshold in that GWAS family
#
# Uses canonical sex-class DEG file + per-GWAS SuSiE-COLOC table; this is the
# pooled-sex enrichment (the sex-stratified A1 analysis is in 216a + the
# downstream COLOC SLURM array).
#
# Output:
#   RNA-seq/results/stratified_causal/sex_coloc_threshold_sweep.csv
# Columns: gwas_family, ancestry, threshold, sex_class, n_class, n_coloc,
#          n_total_class, n_coloc_in_class, odds_ratio, ci_lower, ci_upper,
#          pvalue, padj
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

outdir <- file.path(BASE, "RNA-seq/results/stratified_causal")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

OUT_FILE <- file.path(outdir, "sex_coloc_threshold_sweep.csv")

# ===========================================================================
# 1. Load per-GWAS SuSiE-COLOC table + ancestry/family classification
# ===========================================================================
cat("=== Loading per-GWAS SuSiE-COLOC ===\n")
coloc_all <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))

# Canonical posterior: SuSiE first, ABF fallback when SuSiE is missing.
coloc_all[, pp4 := fifelse(!is.na(PP.H4.susie), PP.H4.susie, PP.H4.abf)]

# GWAS family
coloc_all[, gwas_family := fcase(
  grepl("ALT|AST|GGT", gwas_name, ignore.case = TRUE),     "enzyme",
  grepl("PDFF", gwas_name, ignore.case = TRUE),            "imaging",
  grepl("NAFLD|NASH|Anstee|Steatosis", gwas_name, ignore.case = TRUE), "steatosis_disease",
  grepl("Cirrhosis|CHIRHEP", gwas_name, ignore.case = TRUE), "cirrhosis",
  grepl("HCC|HEPATOCELLU|EXALLC", gwas_name, ignore.case = TRUE), "hcc",
  default = "other"
)]

# Ancestry (parsed from GWAS name; canonical mapping based on registry)
coloc_all[, ancestry := fcase(
  grepl("BBJ_", gwas_name),                       "EAS",
  grepl("_EAS$", gwas_name),                      "EAS",
  grepl("PanUKBB_AFR_", gwas_name),               "AFR",
  grepl("PanUKBB_CSA_", gwas_name),               "SAS",
  default = "EUR"
)]

cat("  Per-GWAS rows:", nrow(coloc_all), "\n")
cat("  GWAS families:\n");  print(coloc_all[, .N, by = gwas_family][order(-N)])
cat("  Ancestries:\n");     print(coloc_all[, .N, by = ancestry][order(-N)])

# ===========================================================================
# 2. Load sex-stratified DEG classification
# ===========================================================================
cat("\n=== Loading sex DEG classification ===\n")
# v3 mashr Bayesian preferred; v2 fallback. v3 CSV provides v2-compatible `sex_class` alias.
sex_v3_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/sex_deg_classification_v3.csv")
sex_v2_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
sex_path <- if (file.exists(sex_v3_path)) sex_v3_path else sex_v2_path
cat("  Source:", basename(dirname(sex_path)), "/", basename(sex_path), "\n")
if (!file.exists(sex_path)) {
  stop("sex_deg_classification.csv not found at ", sex_path)
}
sex <- fread(sex_path)
cat("  Rows:", nrow(sex), "\n")
sex_col_candidates <- c("sex_class", "sex_classification", "sex_dimorphic_class")
sex_col <- intersect(sex_col_candidates, names(sex))[1]
if (is.na(sex_col)) {
  stop("No sex-class column found in ", sex_path,
       " — looked for: ", paste(sex_col_candidates, collapse = ", "))
}
cat("  Using sex-class column:", sex_col, "\n")
setnames(sex, sex_col, "sex_class")
cat("  Sex classes:\n"); print(sex[, .N, by = sex_class][order(-N)])

# Choose gene identifier column for join
id_candidates <- c("ensembl_id", "ensembl", "gene_id", "Geneid")
sex_id <- intersect(id_candidates, names(sex))[1]
if (is.na(sex_id)) {
  symbol_candidates <- c("gene", "symbol", "human_symbol", "gene_symbol")
  sex_id <- intersect(symbol_candidates, names(sex))[1]
  setnames(sex, sex_id, "gene")
  join_key <- "gene"
} else {
  setnames(sex, sex_id, "ensembl")
  join_key <- "ensembl"
}
cat("  Joining on:", join_key, "\n")

# ===========================================================================
# 3. Threshold sweep
# ===========================================================================
THRESHOLDS <- c(0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9)
sex_classes <- setdiff(unique(sex$sex_class), c(NA, ""))

run_one <- function(family_in, ancestry_in, threshold_in) {
  sub_coloc <- coloc_all[gwas_family == family_in &
                         ancestry == ancestry_in &
                         !is.na(pp4)]
  if (nrow(sub_coloc) == 0) return(NULL)

  # Best PP4 per gene within this family/ancestry
  if (join_key == "ensembl") {
    coloc_gene <- sub_coloc[, .(pp4 = max(pp4, na.rm = TRUE)), by = ensembl]
    setnames(coloc_gene, "ensembl", join_key)
  } else {
    coloc_gene <- sub_coloc[, .(pp4 = max(pp4, na.rm = TRUE)), by = gene]
  }

  merged <- merge(sex[, c(join_key, "sex_class"), with = FALSE],
                  coloc_gene, by = join_key, all.x = TRUE)
  merged[, hit := !is.na(pp4) & pp4 >= threshold_in]

  rows <- list()
  for (sc in sex_classes) {
    in_class    <- merged$sex_class == sc
    is_concord  <- merged$sex_class == "Concordant"
    if (sum(in_class) == 0) next

    n_class           <- sum(in_class)
    n_coloc_in_class  <- sum(in_class & merged$hit)
    n_other           <- sum(is_concord)
    n_coloc_in_other  <- sum(is_concord & merged$hit)
    if (n_other == 0) next

    mat <- matrix(c(n_coloc_in_class,           n_class - n_coloc_in_class,
                    n_coloc_in_other,           n_other - n_coloc_in_other),
                  nrow = 2, byrow = TRUE)
    ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
    if (is.null(ft)) next

    rows[[length(rows) + 1L]] <- data.table(
      gwas_family       = family_in,
      ancestry          = ancestry_in,
      threshold         = threshold_in,
      sex_class         = sc,
      n_class           = n_class,
      n_coloc_in_class  = n_coloc_in_class,
      n_concordant      = n_other,
      n_coloc_concordant = n_coloc_in_other,
      odds_ratio        = unname(ft$estimate),
      ci_lower          = ft$conf.int[1],
      ci_upper          = ft$conf.int[2],
      pvalue            = ft$p.value
    )
  }
  if (length(rows) == 0) return(NULL)
  rbindlist(rows)
}

families <- setdiff(unique(coloc_all$gwas_family), c("other", ""))
ancestries <- setdiff(unique(coloc_all$ancestry), c("", NA))

cat("\n=== Sweeping", length(families), "families ×",
    length(ancestries), "ancestries ×",
    length(THRESHOLDS), "thresholds ×",
    length(sex_classes), "sex classes ===\n")

results <- list()
for (fam in families) {
  for (anc in ancestries) {
    for (thr in THRESHOLDS) {
      r <- run_one(fam, anc, thr)
      if (!is.null(r)) results[[length(results) + 1L]] <- r
    }
  }
}

if (length(results) == 0) {
  cat("WARN — no sweep rows produced. Check input tables.\n")
  quit(status = 0)
}

out <- rbindlist(results, fill = TRUE)
out[, padj := p.adjust(pvalue, method = "BH")]

fwrite(out, OUT_FILE)
cat("\nWrote:", OUT_FILE, "(", nrow(out), "rows )\n")

cat("\n--- Significant (padj < 0.1) ---\n")
sig <- out[padj < 0.1][order(padj)]
if (nrow(sig) > 0) print(sig) else cat("  (none)\n")

cat("\nEnd:", format(Sys.time()), "\n")
