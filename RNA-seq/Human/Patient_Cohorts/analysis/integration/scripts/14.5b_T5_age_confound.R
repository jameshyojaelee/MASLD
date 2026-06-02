#!/usr/bin/env Rscript
# 14.5b_T5_age_confound.R
# ---------------------------------------------------------------------------
# A7 T5: Sex x age confound check
#
# The canonical age:sex 3-way interaction file (Script 14e Part 2) is empty
# (insufficient samples in 8 cells of disease x age x sex). We use TWO
# proxies for age-confounding instead:
#   (a) age_stratified_dream.csv: disease:age_group interaction (proxy for
#       age-dependent disease effects).
#   (b) age_continuous_dream.csv: continuous age main effect (genes where
#       expression tracks age across all subjects).
#
# We then overlap Female_biased / Male_biased / Divergent gene sets from
# sex_deg_classification.csv with the "age-significant" gene set (padj<0.05).
# If >50% of a sex_class is in the age-significant set, we flag the class as
# "age-confounded" overall.
#
# Outputs:
#   audit_sensitivity/sex_age_confound/female_biased_age_filtered.csv
#   audit_sensitivity/sex_age_confound/male_biased_age_filtered.csv
#   audit_sensitivity/sex_age_confound/divergent_age_filtered.csv
#   audit_sensitivity/sex_age_confound/summary.csv
# ---------------------------------------------------------------------------

set.seed(42)

suppressPackageStartupMessages({ library(data.table) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
SIGS <- file.path(INT, "results/disease_signatures")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_age_confound")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("=== T5: sex x age confound check ===\n")

# --- Load canonical sex classification ---
sex_cls <- fread(file.path(RDIR, "sex_deg_classification.csv"))
cat("Sex classification rows:", nrow(sex_cls), "\n")
cat("sex_class counts:\n"); print(table(sex_cls$sex_class, useNA = "ifany"))

# --- Load age-related dream tables ---
strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

read_age_table <- function(fp, label) {
  if (!file.exists(fp) || file.info(fp)$size < 100) {
    cat("WARN: missing/empty:", fp, "\n"); return(NULL)
  }
  dt <- fread(fp)
  if (!"gene" %in% names(dt) || !"padj" %in% names(dt)) return(NULL)
  dt[, gene_base := strip_ver(gene)]
  dt
}

age_strat <- read_age_table(file.path(SIGS, "age_stratified_dream.csv"),       "age_stratified")
age_cont  <- read_age_table(file.path(SIGS, "age_continuous_dream.csv"),       "age_continuous")
age_3way  <- read_age_table(file.path(SIGS, "age_sex_interaction_dream.csv"),  "age_sex_3way")

age_proxy_used <- character(0)
if (!is.null(age_3way) && nrow(age_3way) > 0) age_proxy_used <- c(age_proxy_used, "age_sex_3way")
if (!is.null(age_strat))                       age_proxy_used <- c(age_proxy_used, "age_stratified")
if (!is.null(age_cont))                        age_proxy_used <- c(age_proxy_used, "age_continuous")
cat("Age proxies available:", paste(age_proxy_used, collapse = ", "), "\n")

# --- Define age-significant set (union across available proxies) ---
age_sig_set <- character(0)
proxy_records <- list()
if (!is.null(age_3way) && nrow(age_3way) > 0) {
  s <- age_3way[!is.na(padj) & padj < 0.05, gene_base]
  proxy_records[["age_sex_3way"]] <- s
  age_sig_set <- union(age_sig_set, s)
}
if (!is.null(age_strat)) {
  s <- age_strat[!is.na(padj) & padj < 0.05, gene_base]
  proxy_records[["age_stratified"]] <- s
  age_sig_set <- union(age_sig_set, s)
}
if (!is.null(age_cont)) {
  s <- age_cont[!is.na(padj) & padj < 0.05, gene_base]
  proxy_records[["age_continuous"]] <- s
  age_sig_set <- union(age_sig_set, s)
}
cat("Total age-significant genes (union, padj<0.05):", length(age_sig_set), "\n")

# --- Annotate sex_cls with age-confound flag ---
sex_cls[, gene_base := strip_ver(gene)]
sex_cls[, age_significant := gene_base %in% age_sig_set]
# Per-proxy flags for transparency
for (nm in names(proxy_records)) {
  sex_cls[[paste0("age_sig_", nm)]] <- sex_cls$gene_base %in% proxy_records[[nm]]
}

# --- Overlap by sex_class ---
classes_of_interest <- c("Female_biased", "Male_biased", "Divergent")
summary_rows <- list()
for (cl in classes_of_interest) {
  sub <- sex_cls[sex_class == cl]
  n_total <- nrow(sub)
  n_age_sig <- sum(sub$age_significant)
  prop <- if (n_total > 0) n_age_sig / n_total else NA_real_
  age_confounded_class <- isTRUE(prop > 0.5)
  cat(sprintf("  %-14s N=%5d age_sig=%5d (%.1f%%) class_confounded=%s\n",
              cl, n_total, n_age_sig, 100 * prop, age_confounded_class))
  # Save filtered (genes NOT in age-significant set) table
  fwrite(sub, file.path(OUT, sprintf("%s_age_filtered.csv", tolower(cl))))
  summary_rows[[cl]] <- data.table(
    sex_class = cl,
    n_total   = n_total,
    n_age_significant = n_age_sig,
    prop_age_significant = prop,
    class_flagged_age_confounded = age_confounded_class
  )
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE)
summary_dt[, age_proxies_used := paste(age_proxy_used, collapse = ";")]
fwrite(summary_dt, file.path(OUT, "summary.csv"))
cat("\nSummary:\n"); print(summary_dt)

# --- Also write the union age-significant set for traceability ---
fwrite(data.table(gene_base = age_sig_set), file.path(OUT, "age_significant_union.csv"))
cat("Done:", as.character(Sys.time()), "\n")
