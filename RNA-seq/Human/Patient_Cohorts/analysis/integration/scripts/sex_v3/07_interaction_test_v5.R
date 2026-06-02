#!/usr/bin/env Rscript
# sex_v3/07_interaction_test_v5.R
# ---------------------------------------------------------------------------
# Layer 1 — Interaction-test classifier (sex_v5).
#
# Replaces the mashr-based v4 consensus, which produced biologically
# implausible 2,704 F_only vs 3 M_only (900:1). The interaction-test is
# symmetric in F/M by construction — any F:M asymmetry in the output reflects
# biology, not asymmetric power.
#
# Reads:
#   sex_interaction_dream_v3.csv  — cached dream M2 interaction term (β_int, t_int, p_int, padj genome-wide)
#   stratified_ashr_v4.csv         — per-arm β_F, β_M, lfsr_F, lfsr_M, power_M_at_F
#   dream_results_ashr.csv         — Tier-1 disease DEG screen (padj_main<0.05, |LFC|>0.3 → ~4,370, kallisto canonical)
#   sex_deg_classification_v3.csv  — for chr / gene_symbol / gene_biotype joins
#
# Writes:
#   interaction_classifier_v5.csv  — per-gene class_v5_interaction + evidence_tier_A_B_C
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(ashr) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")
INT_CSV   <- file.path(SEXV3, "sex_interaction_dream_v3.csv")
STRAT_CSV <- file.path(SEXV3, "stratified_ashr_v4.csv")
TIER1_CSV <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                       "results/integration/dream_results_ashr.csv")
V3_CSV    <- file.path(SEXV3, "sex_deg_classification_v3.csv")
OUT_CSV   <- file.path(SEXV3, "interaction_classifier_v5.csv")

cat("== Layer 1 interaction-test v5 ==\n")

# Load inputs
v3int  <- fread(INT_CSV)
strat  <- fread(STRAT_CSV)
tier1  <- fread(TIER1_CSV)
v3meta <- fread(V3_CSV, select = c("gene", "gene_symbol", "ensembl_base", "chr",
                                    "chr_category", "gene_biotype"))

# Tier-1 universe: padj_main<0.05 AND |logFC|>0.3 → ~4370 (kallisto canonical; was 5378 under STAR)
tier1[, tier1_DEG := !is.na(padj) & padj < 0.05 & abs(logFC) > 0.3]
tier1_genes <- tier1[tier1_DEG == TRUE, gene]
cat("Tier-1 DEG universe (padj_main<0.05 + |LFC|>0.3):", length(tier1_genes), "\n")

# Build per-gene table
dt <- data.table(
  gene     = v3int$gene,
  beta_int = v3int$logFC,    # logFC IS the interaction coefficient β_M - β_F
  se_int   = v3int$se_int,
  t_int    = v3int$t,
  p_int    = v3int$P.Value,
  padj_int_full = v3int$padj  # genome-wide BH from dream
)

# Recompute padj over Tier-1 universe (independent-filtering BH)
dt[, padj_int_5k := NA_real_]
dt[gene %in% tier1_genes, padj_int_5k := p.adjust(p_int, method = "BH")]
# Also compute Storey q-value (sensitivity)
dt[, qvalue_int := NA_real_]
tryCatch({
  if (requireNamespace("qvalue", quietly = TRUE)) {
    q <- qvalue::qvalue(dt$p_int[!is.na(dt$p_int)], pi0 = NULL)
    dt[!is.na(p_int), qvalue_int := q$qvalues]
  }
}, error = function(e) cat("  qvalue not available; skipping qvalue_int\n"))

# ashr-lfsr on interaction effect (sensitivity)
keep <- !is.na(dt$beta_int) & !is.na(dt$se_int) & dt$se_int > 0
ash_int <- ashr::ash(dt$beta_int[keep], dt$se_int[keep],
                     mixcompdist = "normal", mode = "estimate")
dt[, lfsr_int_ashr := NA_real_]
dt[keep, lfsr_int_ashr := ash_int$result$lfsr]

# Merge per-arm + meta
dt <- merge(dt, strat[, .(gene, beta_F = beta_F_strat, se_F = se_F_strat,
                          lfsr_F = lfsr_F_strat, beta_M = beta_M_strat,
                          se_M = se_M_strat, lfsr_M = lfsr_M_strat,
                          power_M_at_F)], by = "gene", all.x = TRUE)
dt <- merge(dt, v3meta, by = "gene", all.x = TRUE)

cat("\nEmpirical anchors:\n")
cat("  mean |t_int|:", round(mean(abs(dt$t_int), na.rm=TRUE), 3), "\n")
cat("  padj_int_full<0.05:", sum(dt$padj_int_full < 0.05, na.rm=TRUE), "\n")
cat("  padj_int_5k<0.05:",   sum(dt$padj_int_5k   < 0.05, na.rm=TRUE), "\n")
cat("  padj_int_5k<0.10:",   sum(dt$padj_int_5k   < 0.10, na.rm=TRUE), "\n")
cat("  padj_int_5k<0.20:",   sum(dt$padj_int_5k   < 0.20, na.rm=TRUE), "\n")
cat("  lfsr_int_ashr<0.05:", sum(dt$lfsr_int_ashr < 0.05, na.rm=TRUE), "\n")

# Decision tree (plan §1.4)
classify_v5 <- function(bF, bM, bI, lF, lM, pI, pow, ch, chcat) {
  if (is.na(pI) || is.na(bF) || is.na(bM)) {
    # No padj_int_5k (not in Tier-1) — fall back to per-arm class
    if (!is.na(lF) && !is.na(lM)) {
      if (lF < 0.05 && lM < 0.05 && sign(bF) == sign(bM)) return("Concordant")
      if (lF < 0.05 && lM < 0.05 && sign(bF) != sign(bM)) return("Divergent")
    }
    return("Not_DEG_or_NotInUniverse")
  }
  same_sign <- (sign(bF) == sign(bM)) && bF != 0 && bM != 0
  opp_sign  <- (sign(bF) != sign(bM)) && bF != 0 && bM != 0

  if (pI >= 0.20) {
    if (!is.na(lF) && !is.na(lM)) {
      if (lF < 0.05 && lM < 0.05 && same_sign) return("Concordant")
      if (lF < 0.05 || lM < 0.05)              return("Concordant_single_arm")
    }
    return("Not_DEG")
  }

  # padj_int < 0.20 — passes some interaction tier
  if (opp_sign) {
    if (!is.na(lF) && !is.na(lM) && lF < 0.05 && lM < 0.05) return("Divergent")
    if ((!is.na(lF) && lF < 0.05) || (!is.na(lM) && lM < 0.05))
      return("Divergent_one_sided")
  }
  if (same_sign) {
    rF <- abs(bF); rM <- abs(bM)
    if (rF >= 2 * rM) {
      # F dominant
      if (!is.na(lF) && lF < 0.05 && !is.na(lM)) {
        if (lM > 0.20) return("Female_biased")
        if (lM >= 0.05 && lM <= 0.20 && !is.na(pow) && pow < 0.5)
          return("Female_biased_M_underpowered")
      }
    }
    if (rM >= 2 * rF) {
      if (!is.na(lM) && lM < 0.05 && !is.na(lF)) {
        if (lF > 0.20) return("Male_biased")
        if (lF >= 0.05 && lF <= 0.20)
          return("Male_biased_F_underpowered")
      }
    }
    return("Sex_modifier")
  }
  return("Uncertain")
}

dt[, class_v5_interaction := mapply(classify_v5,
                                    beta_F, beta_M, beta_int,
                                    lfsr_F, lfsr_M, padj_int_5k,
                                    power_M_at_F, chr, chr_category)]

# Edge-case post-processing
# chrY F-only is biologically impossible
dt[chr == "chrY" & class_v5_interaction == "Female_biased",
   class_v5_interaction := "Male_biased"]
dt[, suspicious_flag := NA_character_]
dt[chr == "chrY" & class_v5_interaction == "Male_biased" &
   beta_F != 0 & sign(beta_F) == sign(beta_M),
   suspicious_flag := "chrY_reassigned_from_F_only"]

# chrX_aware annotation
dt[, chrX_aware_flag := fcase(
  chr_category == "chrX_escape", "chrX_escape_expected",
  chr == "chrX" & class_v5_interaction %in% c("Female_biased", "Male_biased"),
    "chrX_review",
  default = NA_character_)]

# Evidence tier
dt[, evidence_tier_A_B_C := fcase(
  !is.na(padj_int_5k) & padj_int_5k < 0.05, "A",
  !is.na(padj_int_5k) & padj_int_5k < 0.10, "B",
  !is.na(padj_int_5k) & padj_int_5k < 0.20 & !is.na(beta_int) & abs(beta_int) > 0.5, "C",
  default = NA_character_)]

cat("\n== class_v5_interaction distribution ==\n")
print(table(dt$class_v5_interaction, useNA = "ifany"))
cat("\n== Tier breakdown ==\n")
print(table(dt$class_v5_interaction, dt$evidence_tier_A_B_C, useNA = "ifany"))

# Sanity checks
`%||%` <- function(a, b) if (is.null(a) || (length(a)==1 && is.na(a))) b else a

cat("\n== Sanity ==\n")
sanity_genes <- c("XIST", "DDX3Y", "KDM5D", "RPS4Y1", "ASCL1")
for (g in sanity_genes) {
  row <- dt[gene_symbol == g][1]
  if (nrow(row) > 0) {
    cat(sprintf("  %s: class=%s  tier=%s  bF=%.3f bM=%.3f padj_int_5k=%.3g\n",
                g, row$class_v5_interaction, row$evidence_tier_A_B_C,
                row$beta_F %||% NA, row$beta_M %||% NA, row$padj_int_5k %||% NA))
  } else {
    cat(sprintf("  %s: not in atlas\n", g))
  }
}

# Reorder columns + write
setcolorder(dt, c("gene", "gene_symbol", "ensembl_base", "chr", "chr_category",
                  "gene_biotype",
                  "beta_F", "se_F", "lfsr_F",
                  "beta_M", "se_M", "lfsr_M",
                  "power_M_at_F",
                  "beta_int", "se_int", "t_int", "p_int",
                  "padj_int_5k", "padj_int_full", "lfsr_int_ashr", "qvalue_int",
                  "class_v5_interaction", "evidence_tier_A_B_C",
                  "chrX_aware_flag", "suspicious_flag"))

tmp <- paste0(OUT_CSV, ".tmp"); fwrite(dt, tmp); file.rename(tmp, OUT_CSV)
cat("\nWrote:", OUT_CSV, "  rows:", nrow(dt), "  cols:", ncol(dt), "\n")
