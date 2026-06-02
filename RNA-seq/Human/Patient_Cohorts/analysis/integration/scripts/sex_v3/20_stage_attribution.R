#!/usr/bin/env Rscript
# sex_v3/20_stage_attribution.R
# ---------------------------------------------------------------------------
# Stage-attribution tagger — joins the 4 sex × disease interaction contrasts
# into a single per-gene tag answering "where does sex modulation act in the
# disease trajectory?"
#
# Inputs (read from each contrast's per-contrast SEXV3 dir):
#   - Disease vs Control (D, anchor): sex_v3/sex_deg_classification_v6.csv
#       (canonical cross-pillar output with `class_v6_consensus` ∈
#        {Strong, Moderate, Suggestive, Power_limited, Uncertain},
#        and `class_v6_combined` ∈ {Suggestive_F_biased, Suggestive_M_biased,
#        Suggestive_Divergent, Suggestive_NA, Power_limited, Uncertain})
#   - MASH vs MASL  (P): sex_v3/contrast_mash_vs_masl/sex_deg_classification_v6.csv
#   - MASL vs Ctrl  (E): sex_v3/contrast_masl_vs_ctrl/interaction_fixed.csv  (Tier B)
#   - MASH vs Ctrl  (L): sex_v3/contrast_mash_vs_ctrl/interaction_fixed.csv  (Tier B)
#
# Tagging rules (asymmetric — well-powered contrasts use cross-pillar
# Suggestive; underpowered contrasts use magnitude + direction):
#   pan-stage              D & P ≥ Suggestive, same sign, E + L |β_int|>0.5 concordant
#   progression-only       P ≥ Suggestive, D below Suggestive
#   D-Suggestive-only      D ≥ Suggestive, P below Suggestive
#   early-onset-direction  D + P below Suggestive, |β_int_E|>0.5, |β_int_L|<0.3
#   late-onset-direction   D + P below Suggestive, |β_int_L|>0.5, |β_int_E|<0.3
#   discordant             Sign flips between any of {D, P, E, L} with significant evidence
#   no-signal              All four |β_int| < 0.3 AND no Suggestive call (the vast majority)
#
# Output: stage_attribution_table.csv at the canonical sex_v3 root.
#         Row scope = all ~4,370 Tier-1 Disease-vs-Ctrl DEGs (kallisto canonical; the spec's
#         "all Tier-1 disease DEGs" tagging scope).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

INT_RESULTS <- file.path(BASE,
                         "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                         "results/integration")
SEXV3_ROOT  <- file.path(INT_RESULTS, "sex_v3")

# ---------------------------------------------------------------------------
# 0) Anchor universe: Tier-1 Disease-vs-Ctrl DEGs (padj<0.05 & |logFC|>0.3)
# ---------------------------------------------------------------------------
tier1_csv <- file.path(INT_RESULTS, "dream_results_ashr.csv")
stopifnot(file.exists(tier1_csv))
tier1 <- fread(tier1_csv)
tier1_padj_col  <- intersect(c("adj.P.Val", "padj"), names(tier1))[1]
tier1_logFC_col <- intersect(c("logFC", "log2FoldChange"), names(tier1))[1]
tier1_gene_col  <- intersect(c("gene", "Gene"), names(tier1))[1]
stopifnot(!is.na(c(tier1_padj_col, tier1_logFC_col, tier1_gene_col)))
tier1_genes <- tier1[!is.na(get(tier1_padj_col)) & get(tier1_padj_col) < 0.05 &
                      abs(get(tier1_logFC_col)) > 0.3, get(tier1_gene_col)]
cat("Anchor (Tier-1 Disease-vs-Ctrl) DEGs:", length(tier1_genes), "\n")

# Symbol lookup (multi_evidence atlas has ensembl_id ↔ human_symbol)
sym_csv <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (file.exists(sym_csv)) {
  sym_dt <- fread(sym_csv, select = c("ensembl_id", "human_symbol"))
  setnames(sym_dt, "ensembl_id", "gene")
  setnames(sym_dt, "human_symbol", "symbol")
} else {
  sym_dt <- data.table(gene = tier1_genes, symbol = NA_character_)
}

# ---------------------------------------------------------------------------
# 1) Load Disease-vs-Ctrl consensus (D)
# ---------------------------------------------------------------------------
d_csv <- file.path(SEXV3_ROOT, "sex_deg_classification_v6.csv")
stopifnot(file.exists(d_csv))
D <- fread(d_csv)
# Standardize column names:
beta_col_D <- intersect(c("beta_int", "beta_interaction"), names(D))[1]
padj_col_D <- intersect(c("padj_int_5k_random", "padj_int_5k"), names(D))[1]
class_col_D <- intersect(c("class_v6_consensus", "class_v6"), names(D))[1]
combined_col_D <- intersect(c("class_v6_combined"), names(D))[1]
stopifnot(!is.na(c(beta_col_D, padj_col_D, class_col_D)))
D_sub <- D[, .(
  gene = gene,
  beta_int_D    = get(beta_col_D),
  padj_int_D    = get(padj_col_D),
  class_D       = get(class_col_D),
  class_dir_D   = if (!is.na(combined_col_D)) get(combined_col_D) else NA_character_
)]
cat("D (Disease vs Ctrl):", nrow(D_sub), " rows;",
    " class_v6_consensus tally:\n")
print(table(D_sub$class_D, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 2) Load MASH-vs-MASL consensus (P, Tier A)
# ---------------------------------------------------------------------------
p_csv <- file.path(SEXV3_ROOT, "contrast_mash_vs_masl", "sex_deg_classification_v6.csv")
if (file.exists(p_csv)) {
  P <- fread(p_csv)
  beta_col_P <- intersect(c("beta_int", "beta_interaction"), names(P))[1]
  padj_col_P <- intersect(c("padj_int_5k_random", "padj_int_5k"), names(P))[1]
  class_col_P <- intersect(c("class_v6_consensus"), names(P))[1]
  combined_col_P <- intersect(c("class_v6_combined"), names(P))[1]
  P_sub <- P[, .(
    gene = gene,
    beta_int_P  = get(beta_col_P),
    padj_int_P  = get(padj_col_P),
    class_P     = get(class_col_P),
    class_dir_P = if (!is.na(combined_col_P)) get(combined_col_P) else NA_character_
  )]
  cat("P (MASH vs MASL):", nrow(P_sub), " rows;",
      " class_v6_consensus tally:\n")
  print(table(P_sub$class_P, useNA = "ifany"))
} else {
  warning("MASH-vs-MASL consensus CSV not found at ", p_csv,
          " — emitting NA for P columns")
  P_sub <- data.table(gene = tier1_genes,
                       beta_int_P = NA_real_, padj_int_P = NA_real_,
                       class_P = NA_character_, class_dir_P = NA_character_)
}

# ---------------------------------------------------------------------------
# 3) Load Tier B contrasts (E, L)
# ---------------------------------------------------------------------------
load_tierB <- function(contrast_name, suffix) {
  f <- file.path(SEXV3_ROOT, paste0("contrast_", contrast_name), "interaction_fixed.csv")
  if (!file.exists(f)) {
    warning(contrast_name, " interaction_fixed.csv not found at ", f,
            " — emitting NA for ", suffix, " columns")
    return(data.table(gene = tier1_genes))
  }
  dt <- fread(f)
  setnames(dt, c("beta_int", "se_int", "padj_int_5k", "cohort_Q_pval", "I2_pct", "n_cohorts_identifiable"),
                c(paste0("beta_int_",  suffix),
                  paste0("se_int_",    suffix),
                  paste0("padj_int_",  suffix),
                  paste0("cohort_Q_pval_", suffix),
                  paste0("I2_pct_",    suffix),
                  paste0("n_cohorts_id_", suffix)),
           skip_absent = TRUE)
  cols_to_keep <- intersect(
    c("gene",
      paste0("beta_int_",  suffix),
      paste0("padj_int_",  suffix),
      paste0("cohort_Q_pval_", suffix),
      paste0("I2_pct_",    suffix),
      paste0("n_cohorts_id_", suffix)),
    names(dt))
  dt[, ..cols_to_keep]
}
E_sub <- load_tierB("masl_vs_ctrl", "E")
L_sub <- load_tierB("mash_vs_ctrl", "L")

# ---------------------------------------------------------------------------
# 4) Master join — anchor on Tier-1 Disease-vs-Ctrl genes
# ---------------------------------------------------------------------------
master <- data.table(gene = tier1_genes)
master <- merge(master, sym_dt,  by = "gene", all.x = TRUE)
master <- merge(master, D_sub,   by = "gene", all.x = TRUE)
master <- merge(master, P_sub,   by = "gene", all.x = TRUE)
master <- merge(master, E_sub,   by = "gene", all.x = TRUE)
master <- merge(master, L_sub,   by = "gene", all.x = TRUE)
cat("Master table after join:", nrow(master), " genes x ", ncol(master), " cols\n")

# Guarantee beta_int_E/L + padj_int_E/L exist (Tier B contrasts may be absent)
for (col in c("beta_int_E", "beta_int_L", "padj_int_E", "padj_int_L",
              "beta_int_D", "beta_int_P", "padj_int_D", "padj_int_P",
              "class_D", "class_P", "class_dir_D", "class_dir_P")) {
  if (!col %in% names(master)) {
    if (grepl("^beta_int_|^padj_int_", col)) {
      master[, (col) := NA_real_]
    } else {
      master[, (col) := NA_character_]
    }
  }
}

# ---------------------------------------------------------------------------
# 5) Apply stage-attribution tag rules
# ---------------------------------------------------------------------------
is_suggestive <- function(x) {
  !is.na(x) & x %in% c("Strong", "Moderate", "Suggestive")
}
sign_safe <- function(x) {
  s <- sign(x)
  s[is.na(s)] <- 0L
  s
}

# Direction-of-effect helpers
dir_D <- sign_safe(master$beta_int_D)
dir_P <- sign_safe(master$beta_int_P)
dir_E <- sign_safe(master$beta_int_E)
dir_L <- sign_safe(master$beta_int_L)

# Magnitude flags (for E + L: descriptive Tier B threshold)
strong_E <- !is.na(master$beta_int_E) & abs(master$beta_int_E) > 0.5
strong_L <- !is.na(master$beta_int_L) & abs(master$beta_int_L) > 0.5
weak_E   <- !is.na(master$beta_int_E) & abs(master$beta_int_E) < 0.3
weak_L   <- !is.na(master$beta_int_L) & abs(master$beta_int_L) < 0.3
weak_D   <- !is.na(master$beta_int_D) & abs(master$beta_int_D) < 0.3
weak_P   <- !is.na(master$beta_int_P) & abs(master$beta_int_P) < 0.3

sug_D <- is_suggestive(master$class_D)
sug_P <- is_suggestive(master$class_P)

# Discordant detector: signs disagree among contrasts that carry evidence
has_evid <- function(b, p, sug) (!is.na(b) & abs(b) > 0.3) | (!is.na(p) & p < 0.20) | sug
evid_D <- has_evid(master$beta_int_D, master$padj_int_D, sug_D)
evid_P <- has_evid(master$beta_int_P, master$padj_int_P, sug_P)
evid_E <- strong_E
evid_L <- strong_L

signs_mat <- cbind(
  D = ifelse(evid_D, dir_D, NA),
  P = ifelse(evid_P, dir_P, NA),
  E = ifelse(evid_E, dir_E, NA),
  L = ifelse(evid_L, dir_L, NA)
)
n_evidence_contrasts <- rowSums(!is.na(signs_mat))
# A row is "discordant" if it has ≥ 2 evidence-bearing contrasts AND they
# disagree on sign.
sign_disagree <- apply(signs_mat, 1, function(s) {
  s <- s[!is.na(s)]
  if (length(s) < 2) return(FALSE)
  any(s > 0) && any(s < 0)
})

# Build tags in priority order. First match wins.
tag <- rep(NA_character_, nrow(master))

# 1. discordant first (overrides everything if evidence-bearing contrasts disagree)
tag[is.na(tag) & sign_disagree] <- "discordant"

# 2. pan-stage: D and P both Suggestive same sign AND E + L magnitude-concordant
pan_stage <- is.na(tag) & sug_D & sug_P & dir_D == dir_P & dir_D != 0 &
             strong_E & strong_L & dir_E == dir_D & dir_L == dir_D
tag[pan_stage] <- "pan-stage"

# 3. progression-only: P Suggestive but D below Suggestive
prog_only <- is.na(tag) & sug_P & !sug_D
tag[prog_only] <- "progression-only"

# 4. D-Suggestive-only: D Suggestive but P below Suggestive
d_only <- is.na(tag) & sug_D & !sug_P
tag[d_only] <- "D-Suggestive-only"

# 5. Both Suggestive but different signs or no E/L concordance: tag as
# both-stage (a softer alternative to pan-stage when E/L is undefined / weak)
both_sug <- is.na(tag) & sug_D & sug_P
tag[both_sug] <- "both-stage"

# 6. early-onset-direction: D + P below Suggestive, |β_E|>0.5, |β_L|<0.3
early_dir <- is.na(tag) & !sug_D & !sug_P & strong_E & weak_L
tag[early_dir] <- "early-onset-direction"

# 7. late-onset-direction: D + P below Suggestive, |β_L|>0.5, |β_E|<0.3
late_dir <- is.na(tag) & !sug_D & !sug_P & strong_L & weak_E
tag[late_dir] <- "late-onset-direction"

# 8. no-signal: everything else (no Suggestive, no large magnitude in E or L)
tag[is.na(tag)] <- "no-signal"

master[, stage_attribution_tag := tag]

# Direction tag — derive from D's sub-class when available, else from P
dir_tag <- master[, fifelse(
  !is.na(class_dir_D) & class_dir_D == "Suggestive_F_biased",  "Female",
  fifelse(!is.na(class_dir_D) & class_dir_D == "Suggestive_M_biased",  "Male",
  fifelse(!is.na(class_dir_D) & class_dir_D == "Suggestive_Divergent", "Divergent",
  fifelse(!is.na(class_dir_P) & class_dir_P == "Suggestive_F_biased",  "Female",
  fifelse(!is.na(class_dir_P) & class_dir_P == "Suggestive_M_biased",  "Male",
  fifelse(!is.na(class_dir_P) & class_dir_P == "Suggestive_Divergent", "Divergent",
                                                                       NA_character_))))))]
master[, sex_bias_direction := dir_tag]

# ---------------------------------------------------------------------------
# 6) Write out
# ---------------------------------------------------------------------------
OUT <- file.path(SEXV3_ROOT, "stage_attribution_table.csv")
write_atomic_csv(master, OUT)
cat("\nWrote:", OUT, "  rows:", nrow(master), "  cols:", ncol(master), "\n\n")

cat("=== stage_attribution_tag distribution ===\n")
print(table(master$stage_attribution_tag, useNA = "ifany"))

cat("\n=== Direction × tag cross-tab ===\n")
print(table(master$sex_bias_direction, master$stage_attribution_tag, useNA = "ifany"))

dump_session_info(file.path(SEXV3_ROOT, "intermediates"), "20")
cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
