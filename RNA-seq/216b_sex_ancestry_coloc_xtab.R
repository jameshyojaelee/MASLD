#!/usr/bin/env Rscript
# 216b_sex_ancestry_coloc_xtab.R
# ---------------------------------------------------------------------------
# Sex x Ancestry x COLOC three-way cross-tabulation.
#
# Combines:
#   (1) Pooled 28-GWAS SuSiE-COLOC long table
#       (GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv)
#       -> sex_class = "pooled", ancestry from registry, gwas_family parsed
#   (2) Pan-UKBB sex-stratified COLOC long table from B1
#       (GWAS/finemapping/results/susie_coloc/PanUKBB_{F,M}_{ALT,AST,GGT}/...)
#       -> sex_class in {F,M}, ancestry = EUR, gwas_family = liver_enzyme
#   (3) FinnGen sex-stratified COLOC from this agent (B3)
#       -> currently UNAVAILABLE (see outputs/team_B/B3_blocked.md). Script
#       transparently picks up the rows if and when they appear.
#
# Output:
#   RNA-seq/results/stratified_causal/sex_ancestry_coloc.csv
#
# Columns:
#   gene, ensembl, sex_class, ancestry, gwas_family, gwas_name,
#   pp4_max, pp4_susie_max, pp4_abf_max, method_best,
#   replicated_in_2plus_ancestries (per gene x sex_class x gwas_family),
#   n_gwas_in_cell
#
# Also writes:
#   sex_ancestry_coloc_anova.txt — 3-way linear-model summary
#   sex_ancestry_coloc_summary.csv — N(PP4>0.5) per cell
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
COLOC_DIR <- file.path(BASE_DIR, "GWAS/finemapping/results/susie_coloc")
REGISTRY  <- file.path(BASE_DIR, "GWAS/finemapping/config/gwas_registry.tsv")
OUT_DIR   <- file.path(BASE_DIR, "RNA-seq/results/stratified_causal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_LONG    <- file.path(OUT_DIR, "sex_ancestry_coloc.csv")
OUT_ANOVA   <- file.path(OUT_DIR, "sex_ancestry_coloc_anova.txt")
OUT_SUMMARY <- file.path(OUT_DIR, "sex_ancestry_coloc_summary.csv")

cat("=== Sex x Ancestry x COLOC cross-tab ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Map a GWAS name to its phenotype family.
gwas_family_of <- function(name) {
  fcase(
    grepl("ALT$",          name),                "ALT",
    grepl("AST$",          name),                "AST",
    grepl("GGT$",          name),                "GGT",
    grepl("PDFF",          name),                "PDFF",
    grepl("HCC|HEPATOCELL",name, ignore.case = TRUE),  "HCC",
    grepl("Cirrhosis|CHIRHEP|CIRR", name, ignore.case = TRUE), "Cirrhosis",
    grepl("NASH",          name, ignore.case = TRUE), "NASH",
    grepl("NAFLD",         name, ignore.case = TRUE), "NAFLD",
    default = "Other"
  )
}

best_pp4 <- function(susie, abf) {
  pmax(susie, abf, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Load registry for ancestry annotation
# ---------------------------------------------------------------------------
reg <- as.data.table(read.delim(REGISTRY, stringsAsFactors = FALSE))
setkey(reg, study_name)

# ---------------------------------------------------------------------------
# Source (1): pooled 28-GWAS
# ---------------------------------------------------------------------------
pooled_file <- file.path(COLOC_DIR, "susie_coloc_all_gwas.csv")
if (!file.exists(pooled_file)) {
  stop("Missing pooled COLOC table: ", pooled_file)
}
cat("Loading pooled 28-GWAS table:", pooled_file, "\n")
pooled <- fread(pooled_file)
cat("  rows:", nrow(pooled), " | unique GWAS:",
    length(unique(pooled$gwas_name)), "\n")

# Annotate ancestry + gwas_family + sex_class
pooled[, ancestry    := reg[.(pooled$gwas_name), ancestry]]
pooled[is.na(ancestry), ancestry := "EUR"]  # registry default
pooled[, gwas_family := gwas_family_of(gwas_name)]
pooled[, sex_class   := "pooled"]
pooled[, pp4_max     := best_pp4(PP.H4.susie, PP.H4.abf)]
pooled[, pp4_susie_max := PP.H4.susie]
pooled[, pp4_abf_max   := PP.H4.abf]
pooled[, method_best := fifelse(!is.na(PP.H4.susie) &
                                  PP.H4.susie >= ifelse(is.na(PP.H4.abf), -Inf,
                                                        PP.H4.abf),
                                "susie", "abf")]

keep_cols <- c("gene", "ensembl", "sex_class", "ancestry", "gwas_family",
               "gwas_name", "pp4_max", "pp4_susie_max", "pp4_abf_max",
               "method_best")
pooled_long <- pooled[, ..keep_cols]

# ---------------------------------------------------------------------------
# Source (2): Pan-UKBB sex-stratified (B1) — optional
# ---------------------------------------------------------------------------
panukbb_dirs <- list.files(COLOC_DIR,
                           pattern = "^PanUKBB_[FM]_(ALT|AST|GGT)$",
                           full.names = TRUE)
sex_rows <- list()

if (length(panukbb_dirs) > 0) {
  cat("\nFound", length(panukbb_dirs), "Pan-UKBB sex-strat dirs (B1):\n")
  for (d in panukbb_dirs) {
    tag <- basename(d)
    files <- list.files(d, pattern = "^susie_coloc_chr[0-9]+\\.csv$",
                        full.names = TRUE)
    if (length(files) == 0) { cat("  ", tag, ": no per-chr files\n"); next }
    dt <- rbindlist(lapply(files, fread), fill = TRUE)
    parts <- strsplit(tag, "_", fixed = TRUE)[[1]]
    sex_class <- parts[2]
    trait     <- parts[3]
    dt[, sex_class    := sex_class]
    dt[, ancestry     := "EUR"]
    dt[, gwas_family  := trait]
    dt[, gwas_name    := tag]
    dt[, pp4_max         := best_pp4(PP.H4.susie, PP.H4.abf)]
    dt[, pp4_susie_max   := PP.H4.susie]
    dt[, pp4_abf_max     := PP.H4.abf]
    dt[, method_best := fifelse(!is.na(PP.H4.susie) &
                                  PP.H4.susie >= ifelse(is.na(PP.H4.abf),
                                                        -Inf, PP.H4.abf),
                                "susie", "abf")]
    sex_rows[[tag]] <- dt[, ..keep_cols]
    cat("  ", tag, ":", nrow(dt), "rows\n")
  }
} else {
  cat("\n[INFO] No Pan-UKBB sex-strat COLOC dirs yet (B1 still running or pending).\n")
}

# ---------------------------------------------------------------------------
# Source (3): FinnGen sex-stratified (B3) — optional, currently blocked
# ---------------------------------------------------------------------------
fg_sex_dirs <- list.files(COLOC_DIR,
                          pattern = "^FINNGEN_[FM]_(NAFLD|NASH|HCC)$",
                          full.names = TRUE)
if (length(fg_sex_dirs) > 0) {
  cat("\nFound", length(fg_sex_dirs), "FinnGen sex-strat dirs (B3):\n")
  for (d in fg_sex_dirs) {
    tag <- basename(d)
    files <- list.files(d, pattern = "^susie_coloc_chr[0-9]+\\.csv$",
                        full.names = TRUE)
    if (length(files) == 0) { cat("  ", tag, ": no per-chr files\n"); next }
    dt <- rbindlist(lapply(files, fread), fill = TRUE)
    parts <- strsplit(tag, "_", fixed = TRUE)[[1]]
    sex_class <- parts[2]
    pheno     <- parts[3]
    dt[, sex_class    := sex_class]
    dt[, ancestry     := "EUR"]
    dt[, gwas_family  := pheno]
    dt[, gwas_name    := tag]
    dt[, pp4_max         := best_pp4(PP.H4.susie, PP.H4.abf)]
    dt[, pp4_susie_max   := PP.H4.susie]
    dt[, pp4_abf_max     := PP.H4.abf]
    dt[, method_best := fifelse(!is.na(PP.H4.susie) &
                                  PP.H4.susie >= ifelse(is.na(PP.H4.abf),
                                                        -Inf, PP.H4.abf),
                                "susie", "abf")]
    sex_rows[[tag]] <- dt[, ..keep_cols]
    cat("  ", tag, ":", nrow(dt), "rows\n")
  }
} else {
  cat("\n[INFO] No FinnGen sex-strat COLOC dirs (B3 blocked: see B3_blocked.md).\n")
}

# ---------------------------------------------------------------------------
# Combine
# ---------------------------------------------------------------------------
all_long <- rbindlist(c(list(pooled_long), sex_rows), fill = TRUE)

# Per-gene x sex_class x gwas_family roll-up
cell <- all_long[, .(
  pp4_max         = max(pp4_max,         na.rm = TRUE),
  pp4_susie_max   = max(pp4_susie_max,   na.rm = TRUE),
  pp4_abf_max     = max(pp4_abf_max,     na.rm = TRUE),
  n_gwas_in_cell  = .N,
  n_ancestries    = uniqueN(ancestry),
  ancestries      = paste(sort(unique(ancestry)), collapse = ",")
), by = .(ensembl, gene, sex_class, gwas_family)]

# Replication: same (gene, sex_class, gwas_family) with PP4>0.5 in 2+ ancestries
ancestry_replication <- all_long[
  pp4_max > 0.5,
  .(n_ancestries_replicated = uniqueN(ancestry)),
  by = .(ensembl, sex_class, gwas_family)
]
cell <- merge(cell, ancestry_replication,
              by = c("ensembl", "sex_class", "gwas_family"), all.x = TRUE)
cell[is.na(n_ancestries_replicated), n_ancestries_replicated := 0L]
cell[, replicated_in_2plus_ancestries := n_ancestries_replicated >= 2L]

# Replace -Inf (no observations) with NA
for (col in c("pp4_max", "pp4_susie_max", "pp4_abf_max")) {
  cell[is.infinite(get(col)), (col) := NA_real_]
}

fwrite(cell, OUT_LONG)
cat("\nWrote", nrow(cell), "rows to", OUT_LONG, "\n")

# ---------------------------------------------------------------------------
# Cell-level summary (counts of PP4>0.5 per sex_class x ancestry x family)
# ---------------------------------------------------------------------------
# Pre-cell summary needs ancestry as a column, so summarise the long table
# directly (not the rolled-up cell)
summary_dt <- all_long[, .(
  n_genes_tested = uniqueN(ensembl),
  n_pp4_gt_05    = uniqueN(ensembl[pp4_max > 0.5]),
  n_pp4_gt_08    = uniqueN(ensembl[pp4_max > 0.8]),
  n_pp4_gt_09    = uniqueN(ensembl[pp4_max > 0.9])
), by = .(sex_class, ancestry, gwas_family)]
fwrite(summary_dt, OUT_SUMMARY)
cat("Wrote cell summary to", OUT_SUMMARY, "\n")

# ---------------------------------------------------------------------------
# 3-way model: pp4_max ~ sex_class * ancestry * gwas_family  (per gene fixed)
# (random-effect by gene via lme4 if available; else fixed-effects ANOVA)
# ---------------------------------------------------------------------------
cat("\n--- 3-way ANOVA / mixed model ---\n")

# Drop NA pp4 rows; keep only genes seen in 2+ cells (for meaningful contrast)
mod_dt <- all_long[!is.na(pp4_max)]
# Only declare factor levels that are actually observed — otherwise empty
# levels make aov / lmer choke ("contrasts ... 2 or more levels").
mod_dt[, sex_class   := factor(sex_class,
                               levels = intersect(c("pooled", "F", "M"),
                                                   unique(sex_class)))]
mod_dt[, ancestry    := factor(ancestry,
                               levels = intersect(c("EUR", "EAS", "AFR", "SAS"),
                                                   unique(ancestry)))]
mod_dt[, gwas_family := factor(gwas_family)]

# Build the model formula adaptively so we never include a single-level term.
predictors <- c()
if (nlevels(mod_dt$sex_class)   > 1) predictors <- c(predictors, "sex_class")
if (nlevels(mod_dt$ancestry)    > 1) predictors <- c(predictors, "ancestry")
if (nlevels(mod_dt$gwas_family) > 1) predictors <- c(predictors, "gwas_family")

# Logit-like transform so PP4 bounded in (0,1) is suitable for linear model
mod_dt[, pp4_logit := log((pmin(pmax(pp4_max, 1e-4), 1 - 1e-4)) /
                          (1 - pmin(pmax(pp4_max, 1e-4), 1 - 1e-4)))]

sink(OUT_ANOVA)
cat("Sex x Ancestry x GWAS family ANOVA on logit(PP4_max)\n")
cat("=======================================================\n\n")
cat("Total observations:", nrow(mod_dt), "\n")
cat("Unique genes:", uniqueN(mod_dt$ensembl), "\n")
cat("Unique cells (sex_class x ancestry x gwas_family):",
    nrow(unique(mod_dt[, .(sex_class, ancestry, gwas_family)])), "\n\n")

cat("--- Levels ---\n")
print(table(mod_dt$sex_class))
print(table(mod_dt$ancestry))
print(table(mod_dt$gwas_family))

cat("\n--- Mean logit(PP4) per cell ---\n")
print(mod_dt[, .(mean_logit = mean(pp4_logit, na.rm = TRUE),
                 n          = .N),
             by = .(sex_class, ancestry, gwas_family)
             ][order(sex_class, ancestry, gwas_family)])

# Mixed model with random intercept per gene (lme4 — optional)
have_lme4 <- requireNamespace("lme4", quietly = TRUE)
have_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)

if (length(predictors) == 0) {
  cat("\nWARNING: no predictor has >1 level (only pooled data, single ancestry,\n")
  cat("         and single gwas family). Skipping ANOVA.\n")
  cat("\n--- Sex_class / ancestry / gwas_family levels seen ---\n")
  cat("sex_class:  ", paste(levels(mod_dt$sex_class), collapse = ", "), "\n")
  cat("ancestry:   ", paste(levels(mod_dt$ancestry), collapse = ", "), "\n")
  cat("gwas_family:", paste(levels(mod_dt$gwas_family), collapse = ", "), "\n")
} else {
  rhs <- paste(predictors, collapse = " * ")
  fmla_fixed <- as.formula(paste("pp4_logit ~", rhs))
  fmla_mixed <- as.formula(paste("pp4_logit ~", rhs, "+ (1 | ensembl)"))

  cat("\n--- Predictors with >=2 levels:", paste(predictors, collapse = ", "), "---\n")
  cat("Fixed effects formula:", deparse(fmla_fixed), "\n")
  if (have_lme4) cat("Mixed model formula:  ", deparse(fmla_mixed), "\n")

  ctrl <- if (have_lme4) {
    lme4::lmerControl(check.conv.singular = lme4::.makeCC("ignore",
                                                          tol = 1e-4))
  } else NULL

  fit <- NULL
  if (have_lme4) {
    cat("\n--- lme4 mixed model: ", deparse(fmla_mixed), " ---\n", sep = "")
    fit_call <- function(lmer_fn) {
      lmer_fn(fmla_mixed, data = mod_dt, control = ctrl)
    }
    if (have_lmerTest) {
      fit <- tryCatch(fit_call(lmerTest::lmer),
                      error = function(e) {
                        cat("  lmerTest failed (", conditionMessage(e),
                            "), falling back to lme4::lmer\n", sep = "")
                        NULL
                      })
    }
    if (is.null(fit)) {
      fit <- tryCatch(fit_call(lme4::lmer),
                      error = function(e) {
                        cat("  lme4::lmer failed (", conditionMessage(e), ")\n",
                            sep = "")
                        NULL
                      })
    }
  }

  if (!is.null(fit)) {
    print(summary(fit))
    cat("\n--- ANOVA on lmer fit ---\n")
    print(tryCatch(anova(fit, type = 2),
                   error = function(e) anova(fit)))
  } else {
    cat("\n--- Fixed-effects ANOVA fallback ---\n")
    fit <- tryCatch(aov(fmla_fixed, data = mod_dt),
                    error = function(e) {
                      cat("  aov failed (", conditionMessage(e), ")\n", sep = "")
                      NULL
                    })
    if (!is.null(fit)) print(summary(fit))
  }
}
sink()

cat("Wrote ANOVA to", OUT_ANOVA, "\n")

cat("\n=== Done ===\n")
cat("End:", format(Sys.time()), "\n")

# Print a console preview
cat("\n--- Cell summary preview ---\n")
print(summary_dt[order(-n_pp4_gt_05)][1:min(20, nrow(summary_dt))])
