#!/usr/bin/env Rscript
# 216a_aggregate_sex_strat_coloc.R
# ---------------------------------------------------------------------------
# Aggregate per-chromosome SuSiE-COLOC outputs for the 6 sex-stratified
# Pan-UKBB / Neale Round 2 strata (F+M × ALT/AST/GGT) into one long table,
# then compute per-gene sex-specificity scores.
#
# Output:
#   RNA-seq/results/stratified_causal/sex_stratified_coloc_panukbb.csv
#
# Columns:
#   gene, ensembl, gwas, trait, sex, pp4_susie, pp4_abf, method, n_snps,
#   top_snp, top_snp_PP, sex_specificity_log10ratio_F_M, sex_concordance
#
# `sex_specificity_log10ratio_F_M` uses paired (per-trait) F vs M PP.H4 with
# +/- pseudocount 1e-3 to avoid log(0):
#     log10((PP4_F + eps) / (PP4_M + eps))
# Positive => female-biased causal effect; negative => male-biased.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
COLOC_DIR <- file.path(BASE_DIR, "GWAS/finemapping/results/susie_coloc")
OUT_DIR   <- file.path(BASE_DIR, "RNA-seq/results/stratified_causal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_FILE <- file.path(OUT_DIR, "sex_stratified_coloc_panukbb.csv")

STRATA <- c("PanUKBB_F_ALT", "PanUKBB_F_AST", "PanUKBB_F_GGT",
            "PanUKBB_M_ALT", "PanUKBB_M_AST", "PanUKBB_M_GGT")

cat("=== Aggregating sex-stratified SuSiE-COLOC outputs ===\n")
cat("Start:", format(Sys.time()), "\n\n")

all_rows <- list()

for (g in STRATA) {
  d <- file.path(COLOC_DIR, g)
  if (!dir.exists(d)) {
    cat("  WARNING — missing dir:", d, "\n")
    next
  }
  files <- list.files(d, pattern = "^susie_coloc_chr[0-9]+\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    cat("  WARNING — no per-chr outputs for", g, "\n")
    next
  }
  dt_list <- lapply(files, function(f) tryCatch(fread(f), error = function(e) NULL))
  dt_list <- dt_list[!sapply(dt_list, is.null)]
  if (length(dt_list) == 0) next
  dt <- rbindlist(dt_list, fill = TRUE)

  # Parse stratum tag into trait + sex
  parts  <- strsplit(g, "_", fixed = TRUE)[[1]]
  sex    <- parts[2]
  trait  <- parts[3]

  dt[, gwas := g]
  dt[, trait := trait]
  dt[, sex := sex]

  cat("  ", g, ":", nrow(dt), "genes |",
      "susie>0.5 =", sum(dt$PP.H4.susie > 0.5, na.rm = TRUE), "|",
      "abf>0.5 =", sum(dt$PP.H4.abf > 0.5, na.rm = TRUE), "\n")

  all_rows[[g]] <- dt
}

if (length(all_rows) == 0) {
  cat("\nERROR — no sex-stratified COLOC outputs found. Did the SLURM jobs finish?\n")
  quit(status = 0)  # non-fatal so the threshold-sweep step still runs
}

long <- rbindlist(all_rows, fill = TRUE)

# Drop rows with empty/NA ensembl ID — these are abf_only artifacts from upstream
# `06_susie_coloc.R` where an eQTL block lacked a valid ENSG; they break dcast by
# colliding on (ensembl="", gene="", trait, sex) and silently demote the entire
# wide table to fun.aggregate=length() (i.e. counts not PP4 values).
n_pre <- nrow(long)
long <- long[!is.na(ensembl) & ensembl != ""]
if (nrow(long) < n_pre) {
  cat("  Dropped", n_pre - nrow(long),
      "rows with empty/NA ensembl (abf_only upstream artifacts)\n")
}

# Take best canonical posterior per gene per stratum.
# Prefer SuSiE; fall back to ABF when SuSiE didn't converge.
long[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE), by = seq_len(nrow(long))]

# Wide form: per gene x trait -> F vs M
wide <- dcast(long, ensembl + gene + trait ~ sex,
              value.var = c("PP.H4.susie", "PP.H4.abf", "pp4_best",
                            "method", "n_snps", "top_snp", "top_snp_PP"))

eps <- 1e-3
wide[, sex_specificity_log10ratio_F_M :=
       log10((pmax(pp4_best_F, 0, na.rm = TRUE) + eps) /
             (pmax(pp4_best_M, 0, na.rm = TRUE) + eps))]

# Concordance class (canonical PP4 threshold 0.5)
wide[, sex_concordance := fcase(
  is.na(pp4_best_F) | is.na(pp4_best_M),                "incomplete",
  pp4_best_F > 0.5 & pp4_best_M > 0.5,                  "both_sexes",
  pp4_best_F > 0.5 & pp4_best_M <= 0.5,                 "female_only",
  pp4_best_M > 0.5 & pp4_best_F <= 0.5,                 "male_only",
  default = "neither"
)]

fwrite(long, file.path(OUT_DIR, "sex_stratified_coloc_panukbb_long.csv"))
fwrite(wide, OUT_FILE)

cat("\n--- Wide table summary (per-trait F vs M, PP4 > 0.5) ---\n")
print(wide[, .N, by = .(trait, sex_concordance)][order(trait, -N)])

cat("\nWrote:\n  ", OUT_FILE, "\n  ",
    file.path(OUT_DIR, "sex_stratified_coloc_panukbb_long.csv"), "\n")
cat("End:", format(Sys.time()), "\n")
