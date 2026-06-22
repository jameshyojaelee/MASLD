#!/usr/bin/env Rscript
##############################################################################
# 46d_convergence_evidence.R  (renamed from 46d_bayesian_evidence.R 2026-05-19)
#
# Unsupervised, empirically-calibrated evidence-weighted ranking. Coexists
# with 46b (therapeutic-archetype similarity score). Every atlas gene gets a
# convergence score (NOT a posterior probability — no prior or likelihood is
# specified; the score is a weighted sum of log-approximate-Bayes-factors
# with empirical Brown's-method correlation discount) plus a tier +
# concordance state + per-modality decomposition.
# See docs/refactor/2026-05-19-bayes-to-convergence-rename.md for the full
# rename rationale and mapping table.
#
# 2026-05-21 Govaere2026 signature extension (additive only):
#   - New flag `signature_govaere2026_evidence` (bool) = TRUE if ≥1 of 7
#     Govaere2026 evidence sub-sources is active for the gene:
#       1. metmac_member             (metabolic-macrophage panel)
#       2. lam_member                (lipid-associated macrophage panel)
#       3. gpnmb_macrophage_member   (GPNMB+ macrophage panel)
#       4. il32_axis_member          (IL32 axis panel)
#       5. geomx_sh_vs_ls_padj<0.05  (steatohepatitis vs lipid-laden steatosis)
#       6. geomx_sh_vs_pt_padj<0.05  (steatohepatitis vs portal tract)
#       7. gpnmb_padj<0.05  OR  epithelia_mash_(centrilobular|periportal)_padj<0.05
#          [combined "epithelial/macrophage DE evidence" — gpnmb_padj and the
#          epithelia padjs are upstream-filtered to <0.05 so any non-NA value
#          counts; the OR collapse keeps the count interpretable at ≤7]
#   - New count `signature_govaere2026_n_evidence` (0-7).
#   - New rank `convergence_rank_with_govaere2026`: re-rank applying a soft
#     additive boost of GOVAERE_BOOST_PER_EVIDENCE * n_evidence (capped at
#     GOVAERE_MAX_BOOST) to `log_convergence_odds`, then re-sigmoid to get
#     `convergence_score_with_govaere2026`. The CANONICAL `convergence_score`
#     / `convergence_rank` columns are PRESERVED unchanged.
#
# Rationale: the Govaere2026 signatures (single-cell + spatial + macrophage
# subtype panels from Govaere 2026 NAT MED preprint) elevate panel members
# like GPNMB / FABP5 / IL32 / LPL that already carry strong multi-modal
# evidence but were diluted in the canonical score because their genetic
# component is weak (these are downstream cell-state markers, not GWAS hits).
# The boost is bounded so it cannot promote a no-evidence gene; it only
# breaks ties among genes that already have substantial primary evidence.
#
# Design spec: ~/.claude/plans/let-s-thoroughly-plan-for-witty-pixel.md
#
# Pipeline:
#   1. Load atlas + confounder filter (reuse Script 95 logic)
#   2. Per-sub-contrast BF (Wakefield ABF / COLOC upgrade / Wakefield-from-SE
#      where upstream SE available / Gaussian mixture)
#   3. Empirical K_eff (eigenvalue-effective N on sub-contrast correlation
#      matrix per modality)
#   4. Within-modality aggregation with sqrt(K_eff) discount
#   5. 4-state concordance classification:
#      (Concordant-up / Concordant-down / Genetic+down-coherent / Conflicted)
#      — the third state was renamed from "Protective-LOF" (2026-04-22).
#      The label "Protective-LOF" was misleading: it implied LOF burden
#      evidence (which we don't have) and fired for agonist targets like
#      THRB. The new label is purely algorithmic ("genetic up + coherent
#      expression-down"). A separate sub-flag `inhibitor_target_candidate`
#      restricts to the genes where TWAS direction (sign(twas_z) > 0) makes
#      the inhibitor-target interpretation defensible.
#   6. INTACT genetic-causal channel (replaces the additive S2a + S2b sum;
#      Okamoto 2023 INTACT integrates TWAS + COLOC under one likelihood, so
#      the previous additive sum was double-counting eQTL z-statistics).
#   7. Wakefield ABF in S4/S6 (replacing Sellke-Held lower bounds, which are
#      not calibrated BFs — Held & Ott 2018 §3.2). Uses upstream effect+SE
#      when available (corrected scATAC: SE column; SC pseudobulk DE: t_stat
#      → SE = effect/t_stat). Falls back to Wakefield-from-padj approximation
#      when upstream SE not findable. Per-modality `wakefield_se_recovery`
#      flag records which method was used.
#   8. Empirical Brown's method correction for between-modality correlation
#      (Poole 2016, EmpiricalBrownsMethod). Sum of log-BFs assumes between-
#      modality independence; Brown's method derives an effective DOF from
#      the empirical 7-modality correlation matrix (S6 single-cell dropped
#      from the canonical count 2026-05-22 — see Section 7a comment) and
#      discounts the global sum accordingly.
#   9. Final score + empirical π (Gaussian mixture local FDR) + posterior
#  10. Tier assignment
#  11. Validation: permutation FDR, held-out panels, sensitivity sweeps,
#      modality-correlation diagnostic
#  12. Outputs: primary CSV + diagnostic CSVs + evidence cards (top-200) +
#      genetic+down-coherent + inhibitor-target-candidate + Brown's-correction
#      + S4/S6 SE-recovery sensitivity CSVs
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(mixtools)    # normalmixEM for empirical π and S3 Chronos mixture
  library(parallel)
})

set.seed(42)
N_PERM <- as.integer(Sys.getenv("N_PERM", "1000"))
cat(sprintf("=== 46d Convergence Evidence (N_PERM=%d) ===\n", N_PERM))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUT_PREFIX <- "convergence_evidence"
CARDS_DIR <- file.path(ME, "evidence_cards")
PANELS_DIR <- file.path(BASE, "data/published_gene_panels")

dir.create(CARDS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# CONFIG
# ============================================================================
# Numerical floors / ceilings (matches 46b for log-BF floor)
LOG_BF_FLOOR      <- log(0.1)       # −2.303
LOG_BF_CEIL       <- log(1e6)       # 13.816
MOD_LOG_BF_CEIL   <- log(1e8)       # 18.421 per-modality after K_eff discount
LOG_BF_ACTIVE     <- log(3)         # Jeffreys' "substantial" threshold

# Concordance thresholds
CONCORDANCE_RATIO_THRESH   <- 0.7
GDC_M_SIG_MIN              <- log(10)   # minimum expression-evidence to call genetic+down
GDC_GENETIC_MIN            <- 0.5       # min PP4 for genetic-up flag
GDC_Z_MIN                  <- 2         # alternative: TWAS z
# Inhibitor-target-candidate: requires Genetic+down-coherent state PLUS TWAS
# direction implies inhibition would be protective (sign(twas_z) > 0 AND
# |twas_z| > 2). This is the *biologically defensible* subset of the
# algorithmic state.
INHIBITOR_CANDIDATE_Z_MIN  <- 2

# Tier 1 thresholds (genetically validated). The genetic-causal channel is now
# COLOC (INTACT dropped 2026-06-19), so the Tier 1 gate is COLOC PP.H4 > 0.5 —
# evaluated via s2_coloc$score, which now carries COLOC PP.H4. The back-compat
# S2a-BF and S2b-z clauses are retained (redundant with the COLOC primary gate,
# but they also let strong-TWAS-only genes qualify).
TIER1_INTACT_MIN   <- 0.5           # COLOC PP.H4 > 0.5 (genetic-causal gate; was Multi-INTACT)
TIER1_S2A_BF_MIN   <- log(20)       # PP4 > ~0.95 (back-compat)
TIER1_S2B_Z_MIN    <- 4             # back-compat

# Wakefield prior variance
W_WAKEFIELD_LOGODDS <- 0.04
W_WAKEFIELD_PROTEIN <- 0.21

# COLOC priors (from GWAS/finemapping/src/35s_susie_coloc_broadaway.R)
COLOC_P12 <- 1e-5
COLOC_P2  <- 1e-4
COLOC_PRIOR_ODDS <- COLOC_P12 / COLOC_P2   # = 0.1

# ============================================================================
# 1. LOAD ATLAS
# ============================================================================
cat("\n--- 1. Load atlas ---\n")
atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"))
cat(sprintf("  atlas rows: %d, cols: %d\n", nrow(atlas), ncol(atlas)))

# 2026-05-21 defensive guard: the current 33943x293 atlas omits 5 columns
# the script expects (scenic_regulon_activity_diff, hepatocyte_da_logFC/padj,
# mouse_da_logFC/padj). Add them as all-NA so downstream computations zero
# out cleanly instead of crashing.
.expected_cols <- c("scenic_regulon_activity_diff",
                    "hepatocyte_da_logFC", "hepatocyte_da_padj",
                    "mouse_da_logFC",      "mouse_da_padj")
.missing_cols <- setdiff(.expected_cols, names(atlas))
if (length(.missing_cols) > 0) {
  # handoff 2026-05-31: HARD STOP (was: warn + fill all-NA). A missing L8 / mouse-DA column means the
  # atlas rebuild (Script 35 L8 integration + 45a) did not populate the epigenomic layer — aborting is
  # correct rather than silently zeroing it (the "0-column L8 atlas" bug). Re-run run_atlas_rebuild.sh
  # with the L8 layer present; do not ship a convergence ranking computed on a silently-absent S5 layer.
  stop(sprintf("46d: atlas missing %d expected L8/mouse-DA col(s): %s — epigenomic layer absent. Fix the rebuild (Script 35/45a) instead of zero-filling.",
               length(.missing_cols), paste(.missing_cols, collapse = ", ")))
}

# GENCODE metadata for biotype + chromosome
gm <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))

# Join biotype + chr to atlas
atlas[, ensembl_base := sub("\\..*", "", ensembl_id)]
gm_u <- unique(gm, by = "ensembl_base")
atlas[, gene_biotype_gm := gm_u$gene_biotype[match(ensembl_base, gm_u$ensembl_base)]]
atlas[, chromosome_gm  := gm_u$chromosome[match(ensembl_base, gm_u$ensembl_base)]]

# ============================================================================
# 2. CONFOUNDER FILTER (reuses NMF Script 95 logic, lines 146-169)
# ============================================================================
cat("\n--- 2. Confounder filter ---\n")
sym <- atlas$human_symbol; sym[is.na(sym)] <- ""

drop_chrY  <- !is.na(atlas$chromosome_gm) & atlas$chromosome_gm == "chrY"
drop_chrM  <- !is.na(atlas$chromosome_gm) & atlas$chromosome_gm %in% c("chrM","chrMT")
drop_rRNA  <- !is.na(atlas$gene_biotype_gm) & atlas$gene_biotype_gm %in%
              c("rRNA","Mt_rRNA","rRNA_pseudogene")
drop_rpsl  <- grepl("^RPS[0-9]+[A-Z]*$|^RPL[0-9]+[A-Z]*$", sym) |
              sym %in% c("RPSA","RPLP0","RPLP1","RPLP2") |
              grepl("^RPL[0-9]+P[0-9]+$|^RPS[0-9]+P[0-9]+$", sym)
drop_mrpsl <- grepl("^MRPS[0-9]+[A-Z]*$|^MRPL[0-9]+[A-Z]*$", sym)
drop_ig    <- grepl("^IG[HKL][VJCD][0-9]", sym)
drop_hb    <- sym %in% c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2",
                         "HBM","HBQ1","HBZ")
drop_xesc  <- sym %in% c("XIST","TSIX","KDM6A","DDX3X","EIF1AX","UTX","UBA1",
                          "RLIM","ZFX","RPS4X","EIF2S3","TXLNG")
drop_hla   <- grepl("^HLA-", sym)

excluded <- drop_chrY | drop_chrM | drop_rRNA | drop_rpsl | drop_mrpsl |
            drop_ig | drop_hb | drop_xesc | drop_hla

reason <- rep("", nrow(atlas))
reason[drop_chrY ] <- "chrY"
reason[drop_chrM ] <- paste0(reason[drop_chrM ], ifelse(reason[drop_chrM ]=="","","+"), "chrM")
reason[drop_rRNA ] <- paste0(reason[drop_rRNA ], ifelse(reason[drop_rRNA ]=="","","+"), "rRNA")
reason[drop_rpsl ] <- paste0(reason[drop_rpsl ], ifelse(reason[drop_rpsl ]=="","","+"), "RPS_RPL")
reason[drop_mrpsl] <- paste0(reason[drop_mrpsl], ifelse(reason[drop_mrpsl]=="","","+"), "MRPS_MRPL")
reason[drop_ig   ] <- paste0(reason[drop_ig   ], ifelse(reason[drop_ig   ]=="","","+"), "IG")
reason[drop_hb   ] <- paste0(reason[drop_hb   ], ifelse(reason[drop_hb   ]=="","","+"), "HB")
reason[drop_xesc ] <- paste0(reason[drop_xesc ], ifelse(reason[drop_xesc ]=="","","+"), "Xesc")
reason[drop_hla  ] <- paste0(reason[drop_hla  ], ifelse(reason[drop_hla  ]=="","","+"), "HLA")

atlas[, excluded_from_ranking := excluded]
atlas[, exclusion_reason      := reason]
cat(sprintf("  excluded %d genes: chrY=%d chrM=%d rRNA=%d RPS/RPL=%d MRPS/MRPL=%d IG=%d HB=%d Xesc=%d HLA=%d\n",
            sum(excluded), sum(drop_chrY), sum(drop_chrM), sum(drop_rRNA),
            sum(drop_rpsl), sum(drop_mrpsl), sum(drop_ig), sum(drop_hb),
            sum(drop_xesc), sum(drop_hla)))

# ============================================================================
# IMPORTANT METHODOLOGICAL NOTE [R1 C9, R2 #6, R3 P0-3]
# ---------------------------------------------------------------------------
# This score is a heuristic evidence-weighted ranking, not a formal Bayesian
# posterior or calibrated probability. The "convergence_score" is computed as
# a sigmoid-transformed weighted sum of approximate log-Bayes-factors with an
# empirical Brown's-method correlation discount. Despite using Wakefield ABFs
# and logit transforms, the final quantity is NOT a calibrated posterior
# probability because:
#   (a) the per-modality "BFs" use ad-hoc prior variances (W = 0.04 / 0.21),
#   (b) the mixture-model prior pi is estimated from the same data,
#   (c) between-modality independence is only approximately corrected, and
#   (d) the concordance-state multiplier (0.85 for Genetic+down-coherent) is
#       a subjective design choice.
# The score is useful for ranking genes by multi-modal evidence strength but
# should not be interpreted as P(causal | data). For a formal frequentist
# alternative, see Script 46e (Fisher's combined probability test).
# ============================================================================

# ============================================================================
# 3. BF PRIMITIVES
# ============================================================================

# Clamp helper
clamp_logbf <- function(x, floor = LOG_BF_FLOOR, ceil = LOG_BF_CEIL) {
  x[is.na(x)] <- 0
  pmin(pmax(x, floor), ceil)
}

# Wakefield Approximate Bayes Factor (ABF). Wakefield 2009 Eq 6 (as in
# coloc::approx.bf.estimates):
#   log ABF = 0.5·log(V/(V+W)) + z²/2 · W/(V+W)
# The prefactor argument is V/(V+W) = (1 - shrink); the exp argument is
# W/(V+W) = shrink. ONLY the prefactor was wrong here (it used log(shrink) =
# log(W/(V+W))) — fixed 2026-05-28. The exponent term W/(V+W) is correct and
# was left unchanged.
wakefield_abf <- function(effect, se, W) {
  out <- rep(0, length(effect))
  ok <- !is.na(effect) & !is.na(se) & se > 0 & is.finite(se)
  V <- se[ok]^2
  z <- effect[ok] / se[ok]
  shrink <- W / (V + W)
  # P0-F fix 2026-05-28: prefactor log(1-shrink)=log(V/(V+W)) per Wakefield 2009 / coloc::approx.bf.estimates
  out[ok] <- 0.5 * log(1 - shrink) + (z^2 / 2) * shrink
  clamp_logbf(out)
}

# Sellke-Held p-value lower bound on BF.
# BF ≥ −e · p · ln(p) when p < 1/e, else 1.
# Used when effect + p available but SE not reliable (S4 DA, S6 single-cell
# pseudobulk DE).
sellke_held_bf <- function(p) {
  out <- rep(0, length(p))  # log(BF=1) = 0
  ok <- !is.na(p) & p < 1/exp(1) & p > 0
  out[ok] <- log(-exp(1) * p[ok] * log(p[ok]))
  clamp_logbf(out)
}

# COLOC PP4 posterior → Bayes factor.
# log_BF = logit(PP4) − logit(π_coloc_prior)
# where π_coloc_prior = p12/p2 from the coloc.abf defaults.
coloc_posterior_upgrade <- function(pp4, prior_odds = COLOC_PRIOR_ODDS) {
  out <- rep(0, length(pp4))
  ok <- !is.na(pp4) & pp4 > 0 & pp4 < 1
  posterior_odds <- pp4[ok] / (1 - pp4[ok])
  out[ok] <- log(posterior_odds / prior_odds)
  # Saturate at PP4=1 -> set to CEIL
  out[!is.na(pp4) & pp4 >= 1] <- LOG_BF_CEIL
  clamp_logbf(out)
}

# 2-component Gaussian mixture Bayes factor for essentiality.
# Chronos is well-known to be bimodal (essential vs non-essential). Fit
# 2-component Gaussian; BF = P(essential | x) / P(non-essential | x).
# Returns both per-gene log-BF AND the fit parameters (for diagnostics).
essentiality_gmm_bf <- function(chronos) {
  out <- list(log_bf = rep(0, length(chronos)), fit = NULL)
  ok <- !is.na(chronos) & is.finite(chronos)
  if (sum(ok) < 100) return(out)
  fit <- tryCatch(
    normalmixEM(chronos[ok], k = 2, maxit = 200, epsilon = 1e-4,
                verb = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(out)
  # Identify the "essential" component: the one with lower mu (Chronos is
  # more negative = more essential).
  essential_idx <- which.min(fit$mu)
  non_idx       <- 3 - essential_idx
  lambda_e <- fit$lambda[essential_idx]; mu_e <- fit$mu[essential_idx]; sd_e <- fit$sigma[essential_idx]
  lambda_n <- fit$lambda[non_idx];       mu_n <- fit$mu[non_idx];       sd_n <- fit$sigma[non_idx]
  log_post_e <- log(lambda_e) + dnorm(chronos[ok], mu_e, sd_e, log = TRUE)
  log_post_n <- log(lambda_n) + dnorm(chronos[ok], mu_n, sd_n, log = TRUE)
  # Subtract prior log-odds to get log-BF (not log convergence odds)
  log_bf <- (log_post_e - log_post_n) - (log(lambda_e) - log(lambda_n))
  out$log_bf[ok] <- clamp_logbf(log_bf)
  out$fit <- fit
  out
}

# Spatial evidence score: combine Moran's I with SVG flag.
# Pragmatic formula scoped to give comparable magnitude to other modalities.
spatial_evidence_score <- function(morans_i, is_svg) {
  bf_val <- 1 + 9 * as.integer(is_svg %in% TRUE) +
            5 * pmax(0, ifelse(is.na(morans_i), 0, morans_i))
  bf_val <- pmin(bf_val, 50)
  clamp_logbf(log(bf_val))
}

# SE reconstruction helpers
se_from_tstat <- function(effect, tstat) {
  ok <- !is.na(effect) & !is.na(tstat) & abs(tstat) > 0
  se <- rep(NA_real_, length(effect))
  se[ok] <- abs(effect[ok]) / abs(tstat[ok])
  se
}

se_from_padj <- function(effect, padj) {
  # Wald approximation: z = effect/SE, p ≈ 2*pnorm(-|z|). SE = |effect|/|qnorm(p/2)|.
  # Note: using padj as p is conservative (padj >= p) — SE estimate is larger
  # than truth, so BF will be deflated. Acceptable conservatism.
  se <- rep(NA_real_, length(effect))
  ok <- !is.na(effect) & !is.na(padj) & padj > 0 & padj < 1
  idx <- which(ok)  # explicit integer indices — avoids nested-subset gotcha
  if (length(idx) > 0) {
    z <- abs(qnorm(padj[idx] / 2))
    usable <- z > 0.01 & is.finite(z)
    se[idx[usable]] <- abs(effect[idx[usable]]) / z[usable]
  }
  se
}

# ============================================================================
# 4. PER-SUB-CONTRAST BF COMPUTATION
# ============================================================================
cat("\n--- 4. Per-sub-contrast BF ---\n")

# --- S1: human bulk DE, 5 sub-contrasts ---
# Define each as (effect, se_source, se_arg, effect_sign_source)
compute_s1_bf <- function(atlas) {
  out <- list()
  # S1.1 Overall disease: bulk_logFC / bulk_padj (tstat available)
  se_dream <- se_from_tstat(atlas$bulk_logFC, atlas$bulk_tstat)
  out$overall <- list(
    log_bf = wakefield_abf(atlas$bulk_logFC, se_dream, W_WAKEFIELD_LOGODDS),
    effect = atlas$bulk_logFC,
    padj   = atlas$bulk_padj
  )
  # S1.2 NAFL -> NASH
  se_nn <- se_from_tstat(atlas$nafl_vs_nash_logFC, atlas$nafl_vs_nash_tstat)
  out$nafl_vs_nash <- list(
    log_bf = wakefield_abf(atlas$nafl_vs_nash_logFC, se_nn, W_WAKEFIELD_LOGODDS),
    effect = atlas$nafl_vs_nash_logFC,
    padj   = atlas$nafl_vs_nash_padj
  )
  # S1.3 F2 inflection
  se_f2 <- se_from_padj(atlas$f2_inflection_logFC, atlas$f2_inflection_padj)
  out$f2_inflection <- list(
    log_bf = wakefield_abf(atlas$f2_inflection_logFC, se_f2, W_WAKEFIELD_LOGODDS),
    effect = atlas$f2_inflection_logFC,
    padj   = atlas$f2_inflection_padj
  )
  # S1.4 Advanced fibrosis
  se_af <- se_from_padj(atlas$adv_fib_logFC, atlas$adv_fib_padj)
  out$adv_fib <- list(
    log_bf = wakefield_abf(atlas$adv_fib_logFC, se_af, W_WAKEFIELD_LOGODDS),
    effect = atlas$adv_fib_logFC,
    padj   = atlas$adv_fib_padj
  )
  # S1.5 Sex interaction: use |F - M| as effect, sex_interaction_padj as p
  sex_effect <- atlas$bulk_logFC_F - atlas$bulk_logFC_M
  se_sex <- se_from_padj(sex_effect, atlas$sex_interaction_padj)
  out$sex_interaction <- list(
    log_bf = wakefield_abf(sex_effect, se_sex, W_WAKEFIELD_LOGODDS),
    effect = sex_effect,
    padj   = atlas$sex_interaction_padj
  )
  out
}
s1 <- compute_s1_bf(atlas)
cat(sprintf("  S1: 5 sub-contrasts, mean non-zero log-BF = %.3f\n",
            mean(unlist(lapply(s1, function(x) x$log_bf[x$log_bf > LOG_BF_ACTIVE])), na.rm = TRUE)))

# --- S2a: COLOC PP4 (unweighted max across per-GWAS columns AND aggregated
#   best-SuSiE / best-ABF columns). All GWAS are treated equally: liver-enzyme
#   GWAS (ALT/AST/GGT/PDFF) are the most informative proxies for a liver
#   disease study and should not be down-weighted. Prior 0.5x enzyme penalty
#   removed per reviewer feedback [R1 C9, R2 #6].
coloc_cols <- grep("_coloc_pp4$", names(atlas), value = TRUE)
extra_cols <- intersect(c("coloc_susie_best_pp4","coloc_abf_best_pp4",
                           "broadaway_coloc_pp4","sceqtl_coloc_best_pp4"),
                        names(atlas))
coloc_cols <- union(coloc_cols, extra_cols)

pp4_mat <- as.matrix(atlas[, ..coloc_cols])
pp4_max <- apply(pp4_mat, 1, function(x)
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
best_gwas <- apply(pp4_mat, 1, function(x) {
  if (all(is.na(x))) return(NA_character_)
  coloc_cols[which.max(ifelse(is.na(x), -Inf, x))]
})
cat(sprintf("  S2a: %d COLOC columns, all weighted equally (no enzyme down-weighting)\n",
            length(coloc_cols)))
s2a <- list(log_bf = coloc_posterior_upgrade(pp4_max),
            pp4    = pp4_max,
            best_gwas = best_gwas)
cat(sprintf("  S2a COLOC: %d genes with PP4>0.5 (of %d tested)\n",
            sum(s2a$pp4 > 0.5, na.rm = TRUE), sum(!is.na(s2a$pp4))))

# --- S2b: TWAS (z-score, signed) ---
# Wakefield with z = twas_z, W = 0.04. BF = sqrt(W/(V+W)) * exp(z²/2 · V/(V+W)).
# Special case: we don't have SE explicitly. For z-based Wakefield the formula
# simplifies when we set V=1 (z treats SE as 1 implicitly):
# BF = sqrt(W/(1+W)) * exp(z²/2 · W/(1+W))  — i.e., shrinkage factor W/(1+W) ≈ 0.038.
compute_twas_bf <- function(z, W = W_WAKEFIELD_LOGODDS) {
  out <- rep(0, length(z))
  ok <- !is.na(z) & is.finite(z)
  shrink <- W / (1 + W)
  # P0-F fix 2026-05-28: same Wakefield prefactor bug as wakefield_abf() — prefactor
  # is log(1-shrink)=log(V/(V+W)) (V=1 here), not log(shrink). Exponent W/(1+W) unchanged.
  # Diagnostic-only channel: s2b feeds log_BF_S2b_diag, NOT the scored S2_INTACT.
  out[ok] <- 0.5 * log(1 - shrink) + (z[ok]^2 / 2) * shrink
  clamp_logbf(out)
}
s2b <- list(log_bf = compute_twas_bf(atlas$twas_z),
            z      = atlas$twas_z)
cat(sprintf("  S2b TWAS: %d genes with |z|>2 (of %d tested)\n",
            sum(abs(s2b$z) > 2, na.rm = TRUE), sum(!is.na(s2b$z))))

# --- S2 genetic-causal channel = COLOC PP.H4 (INTACT DROPPED 2026-06-19). ---
# HISTORY: this channel used to be Okamoto-2023 INTACT (a joint TWAS+COLOC
# posterior that replaced the additive S2a+S2b sum to avoid double-counting the
# shared eQTL summary stats). INTACT was DROPPED 2026-06-19 because the on-disk
# Multi-INTACT score (Script 200) was NA/sparse precisely for the high-COLOC
# marquee targets — THRB (PP4=0.9999), RORA (0.998) and HKDC1 (0.992) all
# carried NA INTACT — since INTACT additionally requires a valid cis-TWAS
# predicted-expression model, which this project de-emphasizes (cross-ancestry
# LD mismatch; TWAS skipped for EAS/AFR/SAS). That silently (a) zeroed the
# genetic channel's contribution to n_modalities_active and (b) dropped
# THRB/RORA into Tier 4_Weak despite near-perfect colocalization.
#
# The genetic-causal modality is now COLOC PP.H4 directly (the S2a channel built
# above; log-BF via coloc_posterior_upgrade), which is well-covered (>500 genes
# PP4>0.5) and fires on the hero targets. TWAS (S2b) is retained ONLY as a
# diagnostic channel and as the genetic DIRECTION (sign) source — it no longer
# enters the evidence sum.
#
# NAMING: variables and output columns are coloc_* (renamed from the legacy
# intact_* in the 2026-06-19 Phase-2 rename, after INTACT was dropped). Output
# columns: coloc_genetic_pp4 (== COLOC PP.H4, the genetic-causal score),
# log_BF_S2_coloc (== COLOC log-BF), sign_S2_coloc (genetic direction from TWAS z),
# coloc_channel_source (== "coloc_pp4_direct"), coloc_best_gwas (best-COLOC GWAS).
# Consumers updated in the same change: 46d_baselines_and_lomo, dump_canonical_numbers,
# streamlit_convergence/app.py, fig5_convergence (fallback). 45a's intact branch was
# dead code (atlas never carried intact_score_bulk) and was removed.
coloc_genetic_pp4_vec <- s2a$pp4                 # genetic-causal SCORE = COLOC PP.H4
log_bf_S2_coloc_pre  <- s2a$log_bf              # genetic-causal log-BF = COLOC BF (coloc_posterior_upgrade)
coloc_channel_source         <- "coloc_pp4_direct"
coloc_best_gwas_vec     <- s2a$best_gwas           # per-gene provenance = best-COLOC GWAS column

cat(sprintf("  S2 genetic-causal (COLOC; INTACT dropped): %d genes PP4>0.5, %d>0.7, %d>0.9\n",
            sum(coloc_genetic_pp4_vec > 0.5, na.rm = TRUE),
            sum(coloc_genetic_pp4_vec > 0.7, na.rm = TRUE),
            sum(coloc_genetic_pp4_vec > 0.9, na.rm = TRUE)))

# The single genetic-causal channel in the score is log_bf_S2_coloc (= COLOC).
# Keep s2a$pp4 / s2b$z accessible for state classification (genetic_up flag,
# inhibitor_target_candidate sign).
s2_coloc <- list(log_bf = log_bf_S2_coloc_pre,
                  score  = coloc_genetic_pp4_vec,
                  source = coloc_channel_source,
                  source_vec = coloc_best_gwas_vec)

# --- S3: essentiality (Chronos 2-component GMM) ---
s3_fit <- essentiality_gmm_bf(atlas$essentiality_chronos)
s3 <- list(log_bf = s3_fit$log_bf, fit = s3_fit$fit)
cat(sprintf("  S3 essentiality: mixture fit %s; %d genes with log_BF > log(3)\n",
            ifelse(is.null(s3_fit$fit), "FAILED", "OK"),
            sum(s3$log_bf > LOG_BF_ACTIVE, na.rm = TRUE)))

# --- S4: epigenomic, 3 sub-contrasts ---
# CHANGE #5 (2026-04-22): Sellke-Held replaced by Wakefield ABF.
# Sellke-Held is a *minimum-BF lower bound* on a single p-value (Held & Ott
# 2018 §3.2 explicitly states summing log-bounds across modalities does not
# produce a posterior). Wakefield ABF requires per-gene effect+SE.
#
# For S4.1 hepatocyte DA: load corrected scATAC peak file with explicit SE
# column. Aggregate to gene level (best peak per gene by |effect/se|).
# For S4.2 mouse DA: atlas has only padj. We load the upstream
# diffbind_da_results.csv to get raw pvalue and use se_from_padj() with raw
# p-value (Wald approximation), which is honest. SE recovery is recorded
# per-modality in `wakefield_se_recovery`.
load_hep_da_se <- function() {
  f1 <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv")
  f2 <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/l8_annotated/scatac_da_gene_annotated.csv")
  src <- "missing"
  for (f in c(f1, f2)) {
    if (file.exists(f)) {
      dt <- fread(f)
      if ("se" %in% names(dt) && "logFC" %in% names(dt)) {
        # Hepatocyte cell type
        dt_hep <- dt[!is.na(cell_type) & cell_type == "Hepatocyte"]
        # Best peak per gene by |effect/se|
        dt_hep[, abs_z := ifelse(se > 0, abs(logFC / se), 0)]
        dt_hep <- dt_hep[order(-abs_z)]
        dt_hep <- dt_hep[!duplicated(gene_symbol)]
        return(list(dt = dt_hep[, .(gene_symbol, logFC, se, padj)],
                    src = "exact_SE_from_corrected"))
      } else if ("logFC" %in% names(dt) && "pvalue" %in% names(dt)) {
        dt_hep <- dt[!is.na(cell_type) & cell_type == "Hepatocyte"]
        dt_hep[, abs_score := -log10(pmax(pvalue, 1e-300))]
        dt_hep <- dt_hep[order(-abs_score)]
        dt_hep <- dt_hep[!duplicated(gene_symbol)]
        return(list(dt = dt_hep[, .(gene_symbol, logFC, pvalue, padj)],
                    src = "raw_p_from_original"))
      }
    }
  }
  list(dt = NULL, src = "missing")
}

load_mouse_da_p <- function() {
  # Atlas has padj only. Source diffbind_da_results.csv has raw pvalue.
  f <- file.path(BASE, "Analysis/ATAC/Mouse_Bulk/results/promoter_accessibility.csv")
  if (!file.exists(f)) return(list(dt = NULL, src = "missing"))
  dt <- fread(f)
  # promoter_accessibility lacks raw pvalue. Pull raw pvalue from
  # diffbind_da_results via peak coordinates if available.
  diff_f <- file.path(BASE, "Analysis/ATAC/Mouse_Bulk/results/diffbind_da_results.csv")
  if (file.exists(diff_f) && "peak_coordinate" %in% names(dt)) {
    diff_dt <- fread(diff_f)
    if ("pvalue" %in% names(diff_dt)) {
      diff_dt[, peak_coord := paste0(seqnames, ":", start, "-", end)]
      m <- match(dt$peak_coordinate, diff_dt$peak_coord)
      dt[, mouse_da_pvalue := diff_dt$pvalue[m]]
      return(list(dt = dt, src = "raw_p_from_diffbind"))
    }
  }
  list(dt = dt, src = "padj_fallback")
}

# Hepatocyte DA: build per-gene logFC + SE
hep_da_load <- load_hep_da_se()
hep_da_logFC_v <- atlas$hepatocyte_da_logFC
hep_da_padj_v  <- atlas$hepatocyte_da_padj
hep_da_se_v    <- rep(NA_real_, nrow(atlas))
hep_da_recovery <- "padj_approx"
if (!is.null(hep_da_load$dt) && hep_da_load$src == "exact_SE_from_corrected") {
  m <- match(atlas$human_symbol, hep_da_load$dt$gene_symbol)
  hep_da_se_v <- hep_da_load$dt$se[m]
  # Use atlas logFC if NA, else corrected
  override_lfc <- !is.na(hep_da_load$dt$logFC[m])
  hep_da_logFC_v[override_lfc] <- hep_da_load$dt$logFC[m][override_lfc]
  hep_da_recovery <- "exact"
} else if (!is.null(hep_da_load$dt) && hep_da_load$src == "raw_p_from_original") {
  m <- match(atlas$human_symbol, hep_da_load$dt$gene_symbol)
  raw_p <- hep_da_load$dt$pvalue[m]
  hep_da_se_v <- se_from_padj(hep_da_logFC_v, raw_p)
  hep_da_recovery <- "from_raw_p"
} else {
  hep_da_se_v <- se_from_padj(hep_da_logFC_v, hep_da_padj_v)
  hep_da_recovery <- "from_padj_approx"
}

# Mouse DA: build per-gene SE
mouse_da_load <- load_mouse_da_p()
mouse_da_se_v <- rep(NA_real_, nrow(atlas))
mouse_da_recovery <- "padj_approx"
if (!is.null(mouse_da_load$dt) && mouse_da_load$src == "raw_p_from_diffbind") {
  # mouse_da_load$dt is keyed by gene_symbol but mouse symbols. Atlas already
  # has mouse_da_padj merged in. To get raw p, match by ortholog → use
  # atlas$mouse_ortholog if available; else fall back.
  if ("mouse_ortholog" %in% names(atlas)) {
    m <- match(tolower(atlas$mouse_ortholog), tolower(mouse_da_load$dt$gene_symbol))
    raw_p <- mouse_da_load$dt$mouse_da_pvalue[m]
    mouse_da_se_v <- se_from_padj(atlas$mouse_da_logFC, raw_p)
    mouse_da_recovery <- ifelse(any(!is.na(raw_p)), "from_raw_p", "from_padj_approx")
  } else {
    mouse_da_se_v <- se_from_padj(atlas$mouse_da_logFC, atlas$mouse_da_padj)
  }
} else {
  mouse_da_se_v <- se_from_padj(atlas$mouse_da_logFC, atlas$mouse_da_padj)
  mouse_da_recovery <- "from_padj_approx"
}
cat(sprintf("  S4 SE-recovery: hep_da=%s mouse_da=%s\n",
            hep_da_recovery, mouse_da_recovery))

compute_s4_bf <- function(atlas) {
  out <- list()
  # S4.1 hepatocyte DA peaks: Wakefield ABF with SE recovery
  out$hep_da <- list(
    log_bf = wakefield_abf(hep_da_logFC_v, hep_da_se_v, W_WAKEFIELD_LOGODDS),
    effect = hep_da_logFC_v,
    padj   = hep_da_padj_v,
    recovery = hep_da_recovery
  )
  # S4.2 mouse DA peaks: Wakefield ABF with SE recovery
  out$mouse_da <- list(
    log_bf = wakefield_abf(atlas$mouse_da_logFC, mouse_da_se_v, W_WAKEFIELD_LOGODDS),
    effect = atlas$mouse_da_logFC,
    padj   = atlas$mouse_da_padj,
    recovery = mouse_da_recovery
  )
  # S4.3 SCENIC regulon activity: no p-value; use absolute magnitude as evidence proxy.
  # Empirical BF: treat |activity_diff| > 0.1 as "substantial".
  # 2026-05-21: defensive guard — older atlas builds carry
  # scenic_regulon_activity_diff; current 33943x293 atlas does not. If absent
  # or non-numeric, zero out the S4.3 sub-channel rather than crash.
  act <- atlas$scenic_regulon_activity_diff
  if (is.null(act) || !is.numeric(act)) {
    if (!is.null(act)) {
      warning("compute_s4_bf: scenic_regulon_activity_diff present but non-numeric; coercing to numeric (forced NA on failures)")
      act <- suppressWarnings(as.numeric(act))
    } else {
      act <- rep(NA_real_, nrow(atlas))
    }
  }
  lbf <- ifelse(is.na(act) | abs(act) < 0.05, 0,
                pmin(10 * abs(act), LOG_BF_CEIL))
  out$scenic <- list(
    log_bf = clamp_logbf(lbf),
    effect = act,
    padj   = rep(NA_real_, length(act))
  )
  out
}
s4 <- compute_s4_bf(atlas)
cat(sprintf("  S4: 3 sub-contrasts, mean log-BF (nonzero) = %.3f\n",
            mean(unlist(lapply(s4, function(x) x$log_bf[x$log_bf > 0])), na.rm = TRUE)))

# --- S5: spatial (unsigned) ---
s5 <- list(log_bf = spatial_evidence_score(atlas$spatial_morans_i,
                                atlas$spatial_is_svg))
cat(sprintf("  S5 spatial: %d SVG genes; %d with log_BF > log(3)\n",
            sum(atlas$spatial_is_svg %in% TRUE),
            sum(s5$log_bf > LOG_BF_ACTIVE)))

# --- S6: single-cell, 2 sub-contrasts ---
# CHANGE #5: Sellke-Held replaced by Wakefield ABF. Pseudobulk DE files
# (Hepatocytes_de.csv, etc.) have t_stat columns → exact SE = lfc / t_stat.
load_sc_se <- function() {
  pb_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
  if (!dir.exists(pb_dir)) return(list(per_ct = NULL, src = "missing"))
  ct_files <- list.files(pb_dir, pattern = "_de\\.csv$", full.names = TRUE)
  ct_files <- ct_files[!grepl("allcell", basename(ct_files))]
  if (length(ct_files) == 0) return(list(per_ct = NULL, src = "missing"))
  ct_dt <- rbindlist(lapply(ct_files, function(f) {
    d <- fread(f)
    d[, cell_type := sub("_de\\.csv$", "", basename(f))]
    d
  }), fill = TRUE)
  if (!all(c("lfc","t_stat","pvalue","gene") %in% names(ct_dt))) {
    return(list(per_ct = NULL, src = "missing"))
  }
  # Map ENSG to symbol
  if (grepl("^ENSG", ct_dt$gene[1])) {
    gm_local <- gm[, .(ensembl_base, gene_name)]
    ct_dt[, ensembl_base := sub("\\.\\d+$", "", gene)]
    ct_dt[, gene_symbol := gm_local$gene_name[match(ensembl_base, gm_local$ensembl_base)]]
  } else {
    ct_dt[, gene_symbol := gene]
  }
  list(per_ct = ct_dt, src = "exact_SE_from_tstat")
}

sc_load <- load_sc_se()
sc_hep_logFC_v <- atlas$sc_hepatocyte_logFC
sc_hep_se_v    <- rep(NA_real_, nrow(atlas))
sc_hep_recovery <- "padj_approx"
sc_best_logFC_v <- atlas$sc_best_logFC
sc_best_se_v    <- rep(NA_real_, nrow(atlas))
sc_best_recovery <- "padj_approx"
if (!is.null(sc_load$per_ct)) {
  pb <- sc_load$per_ct
  hep_pb <- pb[!is.na(cell_type) & cell_type == "Hepatocytes"]
  if (nrow(hep_pb) > 0) {
    m <- match(atlas$human_symbol, hep_pb$gene_symbol)
    eff_h <- hep_pb$lfc[m]
    t_h   <- hep_pb$t_stat[m]
    sc_hep_logFC_v <- ifelse(!is.na(eff_h), eff_h, sc_hep_logFC_v)
    sc_hep_se_v <- se_from_tstat(eff_h, t_h)
    sc_hep_recovery <- "exact"
  }
  # Best non-hepatocyte: per gene, the non-hep cell type with smallest pvalue
  nonhep <- pb[!is.na(cell_type) & cell_type != "Hepatocytes" & !is.na(pvalue)]
  if (nrow(nonhep) > 0) {
    nonhep <- nonhep[order(pvalue)]
    nonhep_best <- nonhep[!duplicated(gene_symbol)]
    m2 <- match(atlas$human_symbol, nonhep_best$gene_symbol)
    eff_b <- nonhep_best$lfc[m2]
    t_b   <- nonhep_best$t_stat[m2]
    sc_best_logFC_v <- ifelse(!is.na(eff_b), eff_b, sc_best_logFC_v)
    sc_best_se_v <- se_from_tstat(eff_b, t_b)
    sc_best_recovery <- "exact"
  }
}
# Fallbacks
if (sc_hep_recovery != "exact") {
  sc_hep_se_v <- se_from_padj(sc_hep_logFC_v, atlas$sc_hepatocyte_padj)
  sc_hep_recovery <- "from_padj_approx"
}
if (sc_best_recovery != "exact") {
  sc_best_se_v <- se_from_padj(sc_best_logFC_v, atlas$sc_best_padj)
  sc_best_recovery <- "from_padj_approx"
}
cat(sprintf("  S6 SE-recovery: sc_hep=%s sc_best=%s\n",
            sc_hep_recovery, sc_best_recovery))

compute_s6_bf <- function(atlas) {
  out <- list()
  # S6.1 hepatocyte-specific: Wakefield ABF
  out$hep <- list(
    log_bf = wakefield_abf(sc_hep_logFC_v, sc_hep_se_v, W_WAKEFIELD_LOGODDS),
    effect = sc_hep_logFC_v,
    padj   = atlas$sc_hepatocyte_padj,
    recovery = sc_hep_recovery
  )
  # S6.2 best non-hepatocyte: Wakefield ABF
  out$best <- list(
    log_bf = wakefield_abf(sc_best_logFC_v, sc_best_se_v, W_WAKEFIELD_LOGODDS),
    effect = sc_best_logFC_v,
    padj   = atlas$sc_best_padj,
    recovery = sc_best_recovery
  )
  out
}
s6 <- compute_s6_bf(atlas)

# Sensitivity: Sellke vs Wakefield per-modality log-BFs (S4/S6).
# This CSV lets us show reviewers per-gene how the BF changed.
s4s6_sens <- data.table(
  human_symbol = atlas$human_symbol,
  ensembl_id   = atlas$ensembl_id,
  hep_da_logBF_sellke   = sellke_held_bf(atlas$hepatocyte_da_padj),
  hep_da_logBF_wakefield = s4$hep_da$log_bf,
  hep_da_recovery        = hep_da_recovery,
  mouse_da_logBF_sellke  = sellke_held_bf(atlas$mouse_da_padj),
  mouse_da_logBF_wakefield = s4$mouse_da$log_bf,
  mouse_da_recovery      = mouse_da_recovery,
  sc_hep_logBF_sellke    = sellke_held_bf(atlas$sc_hepatocyte_padj),
  sc_hep_logBF_wakefield = s6$hep$log_bf,
  sc_hep_recovery        = sc_hep_recovery,
  sc_best_logBF_sellke   = sellke_held_bf(atlas$sc_best_padj),
  sc_best_logBF_wakefield = s6$best$log_bf,
  sc_best_recovery        = sc_best_recovery
)
fwrite(s4s6_sens, file.path(ME, sprintf("%s_S4S6_recovery_sensitivity.csv", OUT_PREFIX)))
cat(sprintf("  Wrote %s_S4S6_recovery_sensitivity.csv (%d rows)\n", OUT_PREFIX, nrow(s4s6_sens)))

# --- S7: proteomics ---
se_prot <- se_from_padj(atlas$best_protein_logFC, atlas$best_protein_padj)
s7 <- list(
  log_bf = wakefield_abf(atlas$best_protein_logFC, se_prot, W_WAKEFIELD_PROTEIN),
  effect = atlas$best_protein_logFC,
  padj   = atlas$best_protein_padj
)
cat(sprintf("  S7 proteomics: %d genes with log_BF > log(3)\n",
            sum(s7$log_bf > LOG_BF_ACTIVE)))

# --- S8: mouse bulk ---
se_mouse <- se_from_padj(atlas$mouse_meta_logFC, atlas$mouse_meta_padj)
s8 <- list(
  log_bf = wakefield_abf(atlas$mouse_meta_logFC, se_mouse, W_WAKEFIELD_LOGODDS),
  effect = atlas$mouse_meta_logFC,
  padj   = atlas$mouse_meta_padj
)
cat(sprintf("  S8 mouse: %d genes with log_BF > log(3)\n",
            sum(s8$log_bf > LOG_BF_ACTIVE)))

# ============================================================================
# 5. EMPIRICAL K_eff
# ============================================================================
cat("\n--- 5. Empirical K_eff per modality ---\n")

# Compute K_eff from the eigenvalues of the sub-contrast correlation matrix.
# Bretherton et al. 1999: K_eff = (sum(lambda))^2 / sum(lambda^2).
# Restrict to non-excluded genes for estimation.
est_keff <- function(sub_contrast_mat, label) {
  mat <- sub_contrast_mat[!excluded, , drop = FALSE]
  # Correlations computed on log-BF vectors across genes
  R <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  R[is.na(R)] <- 0
  lam <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  lam <- pmax(lam, 0)  # numerical floor
  keff <- if (sum(lam^2) > 0) sum(lam)^2 / sum(lam^2) else 1
  list(K = ncol(mat), K_eff = keff, eigenvalues = lam, cor_mat = R, label = label)
}

s1_mat <- cbind(s1$overall$log_bf, s1$nafl_vs_nash$log_bf,
                s1$f2_inflection$log_bf, s1$adv_fib$log_bf, s1$sex_interaction$log_bf)
colnames(s1_mat) <- c("overall","nafl_vs_nash","f2_inflection","adv_fib","sex_interaction")
s4_mat <- cbind(s4$hep_da$log_bf, s4$mouse_da$log_bf, s4$scenic$log_bf)
colnames(s4_mat) <- c("hep_da","mouse_da","scenic")
s6_mat <- cbind(s6$hep$log_bf, s6$best$log_bf)
colnames(s6_mat) <- c("hep","best")

keff_s1 <- est_keff(s1_mat, "S1")
keff_s4 <- est_keff(s4_mat, "S4")
keff_s6 <- est_keff(s6_mat, "S6")

cat(sprintf("  S1: K=%d, K_eff=%.3f\n", keff_s1$K, keff_s1$K_eff))
cat(sprintf("  S4: K=%d, K_eff=%.3f\n", keff_s4$K, keff_s4$K_eff))
cat(sprintf("  S6: K=%d, K_eff=%.3f\n", keff_s6$K, keff_s6$K_eff))

# Save K_eff diagnostics
keff_out <- data.table(
  modality = c("S1","S4","S6"),
  K = c(keff_s1$K, keff_s4$K, keff_s6$K),
  K_eff = c(keff_s1$K_eff, keff_s4$K_eff, keff_s6$K_eff),
  eigenvalues = c(paste(round(keff_s1$eigenvalues,3),collapse=";"),
                   paste(round(keff_s4$eigenvalues,3),collapse=";"),
                   paste(round(keff_s6$eigenvalues,3),collapse=";"))
)
fwrite(keff_out, file.path(ME, sprintf("%s_kEff.csv", OUT_PREFIX)))

# ============================================================================
# 6. WITHIN-MODALITY AGGREGATION
# ============================================================================
cat("\n--- 6. Within-modality aggregation ---\n")

# log_BF_i = sum(sub-contrast log-BFs) / sqrt(K_eff)
agg_modality <- function(sub_mat, keff) {
  rs <- rowSums(sub_mat)
  pmin(rs / sqrt(keff), MOD_LOG_BF_CEIL)
}

log_bf_S1  <- agg_modality(s1_mat, keff_s1$K_eff)
# Diagnostic-only channels (the *evidence sum* uses log_bf_S2_coloc instead;
# these are kept for downstream column-level reporting and back-compat).
log_bf_S2a <- pmin(s2a$log_bf, MOD_LOG_BF_CEIL)
log_bf_S2b <- pmin(s2b$log_bf, MOD_LOG_BF_CEIL)
# UNIFIED genetic-causal channel — replaces both S2a and S2b in the score.
log_bf_S2_coloc <- pmin(s2_coloc$log_bf, MOD_LOG_BF_CEIL)
log_bf_S3  <- pmin(s3$log_bf,  MOD_LOG_BF_CEIL)
log_bf_S4  <- agg_modality(s4_mat, keff_s4$K_eff)
log_bf_S5  <- pmin(s5$log_bf,  MOD_LOG_BF_CEIL)
log_bf_S6  <- agg_modality(s6_mat, keff_s6$K_eff)
log_bf_S7  <- pmin(s7$log_bf,  MOD_LOG_BF_CEIL)
log_bf_S8  <- pmin(s8$log_bf,  MOD_LOG_BF_CEIL)

# Per-modality signs (for signed modalities)
# Sign of the sub-contrast with the largest |log_BF|
modality_sign_from_subs <- function(sub_mat, effect_mat) {
  best_idx <- max.col(abs(sub_mat), ties.method = "first")
  rows <- seq_len(nrow(sub_mat))
  eff <- effect_mat[cbind(rows, best_idx)]
  sn <- sign(eff)
  list(sign = sn, best_col = colnames(sub_mat)[best_idx])
}

s1_eff_mat <- cbind(s1$overall$effect, s1$nafl_vs_nash$effect,
                    s1$f2_inflection$effect, s1$adv_fib$effect,
                    s1$sex_interaction$effect)
s4_eff_mat <- cbind(s4$hep_da$effect, s4$mouse_da$effect, s4$scenic$effect)
s6_eff_mat <- cbind(s6$hep$effect, s6$best$effect)

sig_S1 <- modality_sign_from_subs(s1_mat, s1_eff_mat)
sig_S4 <- modality_sign_from_subs(s4_mat, s4_eff_mat)
sig_S6 <- modality_sign_from_subs(s6_mat, s6_eff_mat)
sign_S1 <- sig_S1$sign; best_contrast_S1 <- sig_S1$best_col
sign_S4 <- sig_S4$sign; best_contrast_S4 <- sig_S4$best_col
sign_S6 <- sig_S6$sign; best_contrast_S6 <- sig_S6$best_col
sign_S2b <- sign(s2b$z)
sign_S7 <- sign(s7$effect)
sign_S8 <- sign(s8$effect)

# Map best_contrast_S1 to dominant_stage_S1
stage_map <- c(overall = "early_disease",
               nafl_vs_nash = "nafl_to_nash",
               f2_inflection = "f2_switch",
               adv_fib = "advanced_fibrosis",
               sex_interaction = "sex_dimorphic")
dominant_stage_S1 <- stage_map[best_contrast_S1]

# Cross-species concordance flag (S8)
cross_species_concordant <- (sign(atlas$mouse_meta_logFC) == sign(atlas$bulk_logFC)) &
                            (atlas$mouse_meta_padj < 0.05) &
                            (atlas$bulk_padj < 0.05)
cross_species_concordant[is.na(cross_species_concordant)] <- FALSE

cat(sprintf("  genes with all 9 raw modalities active incl. S2a+S2b+S6 diagnostics (log_BF>log(3)): %d\n",
  sum(log_bf_S1 > LOG_BF_ACTIVE & log_bf_S2a > LOG_BF_ACTIVE & log_bf_S2b > LOG_BF_ACTIVE &
      log_bf_S3 > LOG_BF_ACTIVE & log_bf_S4 > LOG_BF_ACTIVE & log_bf_S5 > LOG_BF_ACTIVE &
      log_bf_S6 > LOG_BF_ACTIVE & log_bf_S7 > LOG_BF_ACTIVE & log_bf_S8 > LOG_BF_ACTIVE)))

# ============================================================================
# 7. CONCORDANCE STATE CLASSIFICATION
# ============================================================================
cat("\n--- 7. Concordance state classification ---\n")

# For each gene:
# - signed modalities = S1, S2b, S4, S6, S7, S8
# - active = log_BF > log(3)
# - E_signed = sum(sign_i * log_BF_i) over active-signed
# - M_signed = sum(log_BF_i) over active-signed
# - concordance_ratio = |E_signed| / M_signed
# - D_majority = sign(E_signed)
#
# Genetic+down-coherent (renamed from "Protective-LOF" 2026-04-22):
#   genetic-up flag AND expression-down majority AND enough expression
#   evidence. The label is purely algorithmic — it does NOT claim LOF burden
#   evidence (we have none) and does NOT claim therapeutic mechanism.
#   THRB (resmetirom is an *agonist*) lands in this state because disease
#   lowers THRB and TWAS z is significant; the algorithm cannot distinguish
#   "disease-driven down" from "loss-of-function protective".
# Genetic-up flag: S2a log_BF corresponds to PP4 > GDC_GENETIC_MIN OR
#   |z_S2b| > GDC_Z_MIN.
# Expression-down: weighted-majority sign over S1, S6, S7, S8 is negative
#   AND M_signed_expression > GDC_M_SIG_MIN.
#
# inhibitor_target_candidate: SUBSET of Genetic+down-coherent where TWAS
# direction makes inhibition the defensible therapy direction:
#   state == "Genetic+down-coherent" AND sign(twas_z) > 0 AND |twas_z| > 2.
# THRB will FAIL this filter (TWAS sign is consistent with "disease lowers
# THRB; agonist restores"); HSD17B13 should similarly fail.

signed_names <- c("S1","S2_coloc","S4","S6","S7","S8")
n <- nrow(atlas)

# Build matrices: bf_signed (n × 6) and sign_signed (n × 6).
# Note: S2_coloc replaces S2b — the magnitude is the INTACT log-BF (which
# already integrates TWAS z + COLOC PP4 under one likelihood) and the sign
# is sign(twas_z), preserving genetic direction for concordance classification.
sign_S2_coloc <- sign(s2b$z)
sign_S2_coloc[is.na(sign_S2_coloc)] <- 0
bf_signed <- cbind(log_bf_S1, log_bf_S2_coloc, log_bf_S4, log_bf_S6, log_bf_S7, log_bf_S8)
sn_signed <- cbind(sign_S1, sign_S2_coloc, sign_S4, sign_S6, sign_S7, sign_S8)
colnames(bf_signed) <- signed_names
colnames(sn_signed) <- signed_names

active_signed <- bf_signed > LOG_BF_ACTIVE
# Zero-out signs for inactive modalities
sn_active <- sn_signed * active_signed
bf_active <- bf_signed * active_signed
bf_active[is.na(bf_active)] <- 0; sn_active[is.na(sn_active)] <- 0

E_signed <- rowSums(sn_active * bf_active)
M_signed <- rowSums(bf_active)
concordance_ratio <- ifelse(M_signed > 0, abs(E_signed) / M_signed, NA_real_)
D_majority <- sign(E_signed)

# Expression modalities subset (S1, S6, S7, S8) for LOF detection
expr_idx <- which(signed_names %in% c("S1","S6","S7","S8"))
bf_expr <- bf_active[, expr_idx, drop = FALSE]
sn_expr <- sn_active[, expr_idx, drop = FALSE]
E_signed_expr <- rowSums(sn_expr * bf_expr)
M_signed_expr <- rowSums(bf_expr)
D_expr_majority <- sign(E_signed_expr)

# Genetic-up flag: S2a PP4 >= 0.5 OR |S2b z| >= 2
genetic_up <- (!is.na(s2a$pp4) & s2a$pp4 >= GDC_GENETIC_MIN) |
              (!is.na(s2b$z) & abs(s2b$z) >= GDC_Z_MIN)

# Expression-down flag: majority sign in expression modalities is negative,
# and M_signed_expr >= threshold
expression_down <- D_expr_majority < 0 & M_signed_expr >= GDC_M_SIG_MIN

# Classify
n_signed_active <- rowSums(active_signed)
concordance_state <- rep("Insufficient-evidence", n)
conc_up <- concordance_ratio >= CONCORDANCE_RATIO_THRESH & D_majority > 0 & !is.na(concordance_ratio)
conc_dn <- concordance_ratio >= CONCORDANCE_RATIO_THRESH & D_majority < 0 & !is.na(concordance_ratio)
gdc <- genetic_up & expression_down & !conc_up & !conc_dn  # mutually exclusive (review A10#7):
# a genuinely concordant gene must not also be flagged Genetic+down-coherent and hit with the
# 0.85 penalty at line ~1040. gdc still takes priority over `conflicted` below.
conflicted <- !conc_up & !conc_dn & !gdc & n_signed_active >= 2

concordance_state[conc_up]     <- "Concordant-up"
concordance_state[conc_dn]     <- "Concordant-down"
concordance_state[gdc]         <- "Genetic+down-coherent"
concordance_state[conflicted]  <- "Conflicted"

# Inhibitor-target-candidate sub-flag: genuine "inhibition would be
# protective" interpretation. Requires state="Genetic+down-coherent" AND
# TWAS direction (sign(twas_z) > 0) implies expression-up associates with
# disease risk AND |twas_z| > INHIBITOR_CANDIDATE_Z_MIN. Genes with NA
# twas_z get FALSE (we cannot defend the inhibitor-target interpretation
# without TWAS direction).
inhibitor_target_candidate <- gdc &
                              !is.na(s2b$z) &
                              s2b$z > 0 &
                              abs(s2b$z) > INHIBITOR_CANDIDATE_Z_MIN

cat(sprintf("  Concordance state distribution:\n"))
print(table(concordance_state[!excluded]))

# Per-state final signed evidence contribution
final_signed_evidence <- rep(0, n)
final_signed_evidence[conc_up | conc_dn] <- M_signed[conc_up | conc_dn]
final_signed_evidence[gdc]                <- 0.85 * M_signed[gdc]
final_signed_evidence[conflicted]         <- abs(E_signed[conflicted])
# Insufficient-evidence: 0

# Unsigned modalities (S5 only).
# NOTE: S3 (DepMap essentiality) is intentionally EXCLUDED from the score
# because it's a constitutive cellular property (essential in HepG2 vs not),
# not MASLD-specific evidence. Pan-essential cell-cycle genes were dominating
# the top of the ranking when S3 was scored. S3 is still computed and
# reported as a column for downstream druggability filtering.
# NOTE 2: S2_coloc (genetic-causal, replacing S2a+S2b) is now SIGNED and
# enters via final_signed_evidence — not E_unsigned — so we don't double-count.
E_unsigned <- pmax(log_bf_S5, 0)

# Total evidence (pre-Brown's-correction; corrected version computed below).
evidence_total_pre_brown <- E_unsigned + final_signed_evidence

# n_modalities_active: 7 canonical channels (S2a/S2b collapsed into S2_coloc;
# S6 single-cell pseudobulk DROPPED 2026-05-22 — Wakefield-ABF prior W=0.04
# yields max log-BF ~0.58, below Jeffreys log(3) threshold, so the channel
# is mathematically unreachable. See phase5/editorial/E6_sources_definition/
# ceiling_verification.md. S6 is retained downstream as a diagnostic log-BF
# column (log_bf_S6) but does NOT enter n_modalities_active or the Brown's
# correction. Effective ceiling is 7; previously-cited "8 channels" was the
# constructed (not effective) ceiling.
log_bf_all_mat <- cbind(log_bf_S1, log_bf_S2_coloc, log_bf_S3,
                        log_bf_S4, log_bf_S5, log_bf_S7, log_bf_S8)
colnames(log_bf_all_mat) <- c("S1","S2_coloc","S3","S4","S5","S7","S8")
n_modalities_active <- rowSums(log_bf_all_mat > LOG_BF_ACTIVE, na.rm = TRUE)
# Diagnostic-only matrix retaining S6 — used for evidence-card per-modality
# log-BF tables and the S4S6 sensitivity sweep; NOT used for counting.
log_bf_all_mat_with_s6 <- cbind(log_bf_all_mat, S6 = log_bf_S6)

# n_signed_concordant: of active-signed modalities, how many agree with D_majority
n_signed_concordant <- rowSums(sn_active == matrix(D_majority, nrow=n, ncol=6) & active_signed, na.rm = TRUE)

# ============================================================================
# 7b. EMPIRICAL BROWN'S METHOD CORRECTION (Change #6)
# ============================================================================
# The sum of log-BFs across modalities assumes between-modality independence.
# In practice modalities are correlated (e.g., S1 dream LFC and S6 sc LFC
# both reflect transcriptional state). Empirical Brown's method (Poole et al.
# 2016, Bioinformatics, PMID 27587659) derives an effective DOF from the
# empirical correlation matrix and adjusts the combined statistic accordingly.
#
# We adapt Brown's framework by computing K_eff_global from the 7-modality
# log-BF correlation matrix (canonical channels; S6 dropped), then discounting
# evidence_total by sqrt(K/K_eff). This is the cross-modality analog of the
# within-modality K_eff discount.
cat("\n--- 7b. Brown's-method between-modality correction ---\n")
brown_K <- ncol(log_bf_all_mat)   # 7 (canonical, S6 dropped 2026-05-22)
brown_mat <- log_bf_all_mat[!excluded, , drop = FALSE]
# Spearman correlation across modalities
brown_cor <- suppressWarnings(cor(brown_mat, method = "spearman", use = "pairwise.complete.obs"))
brown_cor[is.na(brown_cor)] <- 0
brown_lam <- eigen(brown_cor, symmetric = TRUE, only.values = TRUE)$values
brown_lam <- pmax(brown_lam, 0)
K_eff_global <- if (sum(brown_lam^2) > 0) sum(brown_lam)^2 / sum(brown_lam^2) else brown_K
brown_discount <- sqrt(brown_K / K_eff_global)   # ≥ 1 (we DIVIDE by it)
cat(sprintf("  Brown's: K=%d, K_eff_global=%.3f, discount factor=%.3f\n",
            brown_K, K_eff_global, brown_discount))

# Apply discount to the global evidence_total. The pre-Brown evidence is
# preserved as a diagnostic column.
evidence_total <- evidence_total_pre_brown / brown_discount

# Save Brown's-correction diagnostic
brown_cor_dt <- as.data.table(as.data.frame(as.table(brown_cor)))
setnames(brown_cor_dt, c("Mod_A","Mod_B","Spearman_rho"))
brown_summary <- data.table(
  modality_pair = paste(brown_cor_dt$Mod_A, brown_cor_dt$Mod_B, sep="↔"),
  spearman_rho  = brown_cor_dt$Spearman_rho
)
brown_summary[, K := brown_K]
brown_summary[, K_eff_global := K_eff_global]
brown_summary[, discount := brown_discount]
fwrite(brown_summary, file.path(ME, sprintf("%s_brown_correction.csv", OUT_PREFIX)))
cat(sprintf("  Wrote %s_brown_correction.csv (correlation matrix + K_eff)\n", OUT_PREFIX))

# Also save per-gene before/after (just the columns; gene info is in primary CSV)
brown_pergene <- data.table(
  human_symbol = atlas$human_symbol,
  ensembl_id   = atlas$ensembl_id,
  evidence_total_pre_brown  = evidence_total_pre_brown,
  evidence_total_post_brown = evidence_total,
  brown_discount = brown_discount,
  K_eff_global   = K_eff_global
)
fwrite(brown_pergene, file.path(ME, sprintf("%s_brown_pergene.csv", OUT_PREFIX)))

# ============================================================================
# 8. EMPIRICAL PRIOR π + POSTERIOR
# ============================================================================
cat("\n--- 8. Empirical prior π via 2-component Gaussian mixture ---\n")

# Fit mixture on non-excluded evidence_total.
# Problem: evidence_total has a huge spike at 0 (most genes have no evidence),
# which breaks a 2-Gaussian mixture. Fit only on genes with E > small_cutoff
# and estimate π as (fraction_above_cutoff × mixture_alt_weight).
E_all <- evidence_total[!excluded]
E_all <- E_all[is.finite(E_all)]
SMALL_CUTOFF <- 0.5   # genes with E <= 0.5 contribute essentially no signal
E_fit <- E_all[E_all > SMALL_CUTOFF]
fraction_above <- length(E_fit) / length(E_all)
cat(sprintf("  fraction of genes with E > %.1f : %.3f\n",
            SMALL_CUTOFF, fraction_above))

fit_pi <- tryCatch(
  suppressWarnings(
    suppressMessages(
      normalmixEM(log1p(E_fit), k = 2, maxit = 500, epsilon = 1e-5, verb = FALSE)
    )
  ),
  error = function(e) NULL
)

if (is.null(fit_pi)) {
  cat("  Mixture fit FAILED; falling back to π = 0.05\n")
  pi_hat <- 0.05
  pi_sd  <- NA
  mix_params <- NULL
} else {
  # Alt component = higher mean (of log1p E). Its weight is the alt fraction
  # among E > cutoff; total π = that × fraction_above.
  alt_idx <- which.max(fit_pi$mu)
  pi_within_positive <- fit_pi$lambda[alt_idx]
  pi_hat <- pi_within_positive * fraction_above
  # Bootstrap π for CI
  pi_boot <- replicate(100, {
    s <- sample(E_fit, length(E_fit), replace = TRUE)
    fb <- tryCatch(
      suppressWarnings(suppressMessages(
        normalmixEM(log1p(s), k = 2, maxit = 200, verb = FALSE))),
      error = function(e) NULL)
    if (is.null(fb)) NA_real_
    else fb$lambda[which.max(fb$mu)] * fraction_above
  })
  pi_sd <- sd(pi_boot, na.rm = TRUE)
  mix_params <- list(lambda = fit_pi$lambda, mu = fit_pi$mu, sigma = fit_pi$sigma,
                     fraction_above = fraction_above,
                     pi_within_positive = pi_within_positive)
  cat(sprintf("  π_hat = %.4f ± %.4f (alt log1p-mu=%.2f, null log1p-mu=%.2f, pi_within_positive=%.3f, frac>cutoff=%.3f)\n",
              pi_hat, pi_sd, fit_pi$mu[alt_idx], fit_pi$mu[3 - alt_idx],
              pi_within_positive, fraction_above))
  # Guard against degenerate fit
  if (pi_hat < 0.001 || pi_hat > 0.5) {
    cat(sprintf("  WARNING: π_hat out of sensible range (%.4f); falling back to 0.05\n", pi_hat))
    pi_hat <- 0.05
  }
}

# Save calibration diagnostics
calib <- data.table(
  metric = c("pi_hat","pi_sd","null_mu","null_sigma","alt_mu","alt_sigma",
             "null_lambda","alt_lambda"),
  value  = c(pi_hat, pi_sd,
             if (!is.null(mix_params)) mix_params$mu[3 - which.max(mix_params$mu)] else NA,
             if (!is.null(mix_params)) mix_params$sigma[3 - which.max(mix_params$mu)] else NA,
             if (!is.null(mix_params)) mix_params$mu[which.max(mix_params$mu)] else NA,
             if (!is.null(mix_params)) mix_params$sigma[which.max(mix_params$mu)] else NA,
             if (!is.null(mix_params)) mix_params$lambda[3 - which.max(mix_params$mu)] else NA,
             if (!is.null(mix_params)) mix_params$lambda[which.max(mix_params$mu)] else NA)
)
fwrite(calib, file.path(ME, sprintf("%s_calibration.csv", OUT_PREFIX)))

# Posterior
log_prior_odds <- log(pi_hat / (1 - pi_hat))
log_convergence_odds  <- log_prior_odds + evidence_total
convergence_score <- 1 / (1 + exp(-log_convergence_odds))

# Rank (NA for excluded)
convergence_rank <- rep(NA_integer_, n)
ok_idx <- which(!excluded)
ord <- order(-convergence_score[ok_idx])
convergence_rank[ok_idx[ord]] <- seq_along(ord)

cat(sprintf("  convergence_score summary (non-excluded): min=%.4f, median=%.4f, max=%.4f\n",
            min(convergence_score[!excluded]), median(convergence_score[!excluded]),
            max(convergence_score[!excluded])))

# ============================================================================
# 9. TIER ASSIGNMENT
# ============================================================================
cat("\n--- 9. Tier assignment ---\n")

# NA-safe boolean helper: treat NA as FALSE. Critical — `NA & anything = NA`
# and will silently leave tiers at default.
iT <- function(x) ifelse(is.na(x), FALSE, x)

tier <- rep("4_Weak", n)
# Tier 1 (genetically validated): primary gate is COLOC PP.H4 > 0.5 (the
# genetic-causal channel; s2_coloc$score now carries COLOC after the INTACT
# drop, 2026-06-19). Back-compat clauses (S2a-BF redundant; S2b-z) also admit
# strong-TWAS-only genes. This is what pulls THRB/RORA (PP4 ~0.998-1.0) back
# into Tier 1 — INTACT had been NA for them, sinking them to 4_Weak.
tier1 <- iT(s2_coloc$score > TIER1_INTACT_MIN) |
         iT(log_bf_S2a > TIER1_S2A_BF_MIN) |
         iT(abs(s2b$z) > TIER1_S2B_Z_MIN)
# Require genetic evidence + at least 1 other active modality (i.e., 2 total
# including S2_coloc). Earlier ≥3 bar was too strict.
tier1 <- tier1 & iT(n_modalities_active >= 2) & !excluded
tier[tier1] <- "1_Genetic_validated"
tier2 <- !tier1 & iT(n_modalities_active >= 5) & iT(convergence_score > 0.9) &
         (concordance_state %in% c("Concordant-up","Concordant-down",
                                   "Genetic+down-coherent")) &
         !excluded
tier[tier2] <- "2_Convergent"
tier3 <- !tier1 & !tier2 & iT(n_modalities_active %in% c(3,4)) & iT(convergence_score > 0.7) &
         !excluded
tier[tier3] <- "3_Suggestive"
tier[excluded] <- "Excluded"

cat(sprintf("  Tier distribution:\n"))
print(table(tier))

# ============================================================================
# 10. TEST CASE ASSERTIONS
# ============================================================================
cat("\n--- 10. Test case assertions ---\n")

assert_gene <- function(sym, expected_tier = NULL, expected_state = NULL, note = "") {
  idx <- which(atlas$human_symbol == sym)
  if (length(idx) == 0) {
    cat(sprintf("  [MISSING] %s not in atlas\n", sym)); return(invisible(NULL))
  }
  t <- tier[idx]; s <- concordance_state[idx]; p <- convergence_score[idx]
  r <- convergence_rank[idx]
  cat(sprintf("  %s: tier=%s  state=%s  post=%.3f  rank=%s  %s\n",
              sym, t, s, p, ifelse(is.na(r),"NA",as.character(r)), note))
  if (!is.null(expected_tier) && t != expected_tier)
    cat(sprintf("    WARN: expected tier %s, got %s\n", expected_tier, t))
  if (!is.null(expected_state) && s != expected_state)
    cat(sprintf("    WARN: expected state %s, got %s\n", expected_state, s))
}

# Deep debug dump for test cases so we can see per-modality values
debug_dump <- function(sym) {
  idx <- which(atlas$human_symbol == sym)
  if (length(idx) == 0) return(invisible(NULL))
  cat(sprintf("\n  [DEBUG %s]  log_BF per modality (active? sign):\n", sym))
  cat(sprintf("    S1=%.2f (%s, sign=%s, best=%s)\n",
    log_bf_S1[idx], log_bf_S1[idx] > LOG_BF_ACTIVE,
    sign_S1[idx], best_contrast_S1[idx]))
  cat(sprintf("    S2a=%.2f (%s, PP4=%.3f, gwas=%s)\n",
    log_bf_S2a[idx], log_bf_S2a[idx] > LOG_BF_ACTIVE,
    ifelse(is.na(s2a$pp4[idx]),0,s2a$pp4[idx]),
    ifelse(is.na(s2a$best_gwas[idx]),"NA",s2a$best_gwas[idx])))
  cat(sprintf("    S2b=%.2f (%s, z=%.2f) [diagnostic-only post-INTACT]\n",
    log_bf_S2b[idx], log_bf_S2b[idx] > LOG_BF_ACTIVE,
    ifelse(is.na(s2b$z[idx]),0,s2b$z[idx])))
  cat(sprintf("    S2_coloc=%.2f (%s, score=%.3f, sign=%s) [unified genetic-causal channel]\n",
    log_bf_S2_coloc[idx], log_bf_S2_coloc[idx] > LOG_BF_ACTIVE,
    ifelse(is.na(s2_coloc$score[idx]),0,s2_coloc$score[idx]),
    sign_S2_coloc[idx]))
  cat(sprintf("    S3=%.2f  S4=%.2f (sign=%s)  S5=%.2f\n",
    log_bf_S3[idx], log_bf_S4[idx], sign_S4[idx], log_bf_S5[idx]))
  cat(sprintf("    S6=%.2f (sign=%s)  S7=%.2f (sign=%s)  S8=%.2f (sign=%s)\n",
    log_bf_S6[idx], sign_S6[idx], log_bf_S7[idx], sign_S7[idx],
    log_bf_S8[idx], sign_S8[idx]))
  cat(sprintf("    sub-contrast log_BFs S1: overall=%.2f nafl_vs_nash=%.2f f2=%.2f adv_fib=%.2f sex=%.2f\n",
    s1_mat[idx,1], s1_mat[idx,2], s1_mat[idx,3], s1_mat[idx,4], s1_mat[idx,5]))
  cat(sprintf("    n_signed_active=%d, E_signed=%.2f, M_signed=%.2f, ratio=%s\n",
    rowSums(active_signed)[idx], E_signed[idx], M_signed[idx],
    ifelse(is.na(concordance_ratio[idx]),"NA",sprintf("%.2f",concordance_ratio[idx]))))
  cat(sprintf("    genetic_up=%s, expression_down=%s (D_expr_maj=%s, M_signed_expr=%.2f)\n",
    genetic_up[idx], expression_down[idx], D_expr_majority[idx], M_signed_expr[idx]))
  cat(sprintf("    flags: conc_up=%s conc_dn=%s gdc=%s conflicted=%s inhibitor_cand=%s\n",
    conc_up[idx], conc_dn[idx], gdc[idx], conflicted[idx], inhibitor_target_candidate[idx]))
}
debug_dump("TM6SF2"); debug_dump("THRB"); debug_dump("HSD17B13"); debug_dump("HKDC1")

# NOTE: THRB is disease-downregulated AND genetically validated, which
# triggers the Genetic+down-coherent category (renamed from Protective-LOF)
# even though Resmetirom is an AGONIST (not an inhibitor). The label is
# algorithmic — it does NOT claim therapeutic mechanism. The pattern reflects
# "disease lowers THRB; agonist therapy restores", not "LOF protective".
# THRB will likely NOT be in inhibitor_target_candidate (no defensible
# inhibitor-mechanism claim from atlas data alone). This is a deliberate
# behavior of the new sub-classification — not a bug.
assert_gene_either_state <- function(sym, expected_tier, expected_states, note = "") {
  idx <- which(atlas$human_symbol == sym)
  if (length(idx) == 0) {
    cat(sprintf("  [MISSING] %s not in atlas\n", sym)); return(invisible(NULL))
  }
  t <- tier[idx]; s <- concordance_state[idx]; p <- convergence_score[idx]
  r <- convergence_rank[idx]
  inh <- inhibitor_target_candidate[idx]
  cat(sprintf("  %s: tier=%s  state=%s  post=%.3f  rank=%s  inh_cand=%s  %s\n",
              sym, t, s, p, ifelse(is.na(r),"NA",as.character(r)), inh, note))
  if (!is.null(expected_tier) && t != expected_tier)
    cat(sprintf("    WARN: expected tier %s, got %s\n", expected_tier, t))
  if (!is.null(expected_states) && !(s %in% expected_states))
    cat(sprintf("    WARN: expected state in {%s}, got %s\n",
                paste(expected_states, collapse=","), s))
}

assert_gene_either_state("THRB", "1_Genetic_validated",
                          c("Concordant-down","Genetic+down-coherent"),
                          "FDA-approved (Resmetirom), PP4=1.0 UKBB_GGT — agonist target; algorithmic Genetic+down-coherent label is OK; should NOT be inhibitor candidate")
assert_gene("HKDC1",   "1_Genetic_validated", "Concordant-up",
            "F2->F3 driver, PP4=0.992 UKBB_AST")
assert_gene("TM6SF2",  NULL, "Genetic+down-coherent",
            "canonical MASLD LOF (BBJ_AST PP4=0.64)")
assert_gene("HSD17B13", NULL, NULL,
            "known lit Protective-LOF but splice-isoform level; expected NOT to classify (and NOT in inhibitor candidates)")
assert_gene("NR1H4",   NULL, NULL, "Obeticholic acid; no genetic support")
assert_gene("PPARA",   NULL, NULL, "Elafibranor; no genetic support")
assert_gene("GAS6",    NULL, "Concordant-up", "novel; PP4=0.997")
assert_gene("CHI3L1",  NULL, "Concordant-up", "plasma biomarker; PP4=0.894")
# Explicit inhibitor-candidate check
cat(sprintf("\n  HSD17B13 inhibitor_target_candidate = %s (expected: FALSE)\n",
  if (length(which(atlas$human_symbol=="HSD17B13")) > 0)
    inhibitor_target_candidate[which(atlas$human_symbol=="HSD17B13")[1]] else NA))

# ============================================================================
# 10b. GOVAERE2026 SIGNATURE EVIDENCE (additive, post-hoc boost)
# ============================================================================
# Added 2026-05-21. Treats Govaere2026 macrophage/spatial signature columns
# as a supplementary 8th evidence source that elevates panel-member genes
# (GPNMB, FABP5, IL32, MetMac/LAM/GPNMB+ macrophage subset) in a parallel
# ranking column without disturbing the canonical convergence_score.
cat("\n--- 10b. Govaere2026 signature evidence ---\n")

# Sub-flag 1-4: membership in macrophage/IL32 panels (boolean)
gv_metmac <- atlas$signature_govaere2026_metmac_member %in% 1
gv_lam    <- atlas$signature_govaere2026_lam_member %in% 1
gv_gpnmb_mac <- atlas$signature_govaere2026_gpnmb_macrophage_member %in% 1
gv_il32   <- atlas$signature_govaere2026_il32_axis_member %in% 1

# Sub-flag 5-6: GeoMx spatial DE (steatohepatitis vs lipid-laden / portal)
# geomx_sh_vs_ls is partially significant (199/205 non-NA pass <0.05) —
# apply explicit padj threshold.
gv_geomx_ls <- !is.na(atlas$signature_govaere2026_geomx_sh_vs_ls_padj) &
               atlas$signature_govaere2026_geomx_sh_vs_ls_padj < 0.05
gv_geomx_pt <- !is.na(atlas$signature_govaere2026_geomx_sh_vs_pt_padj) &
               atlas$signature_govaere2026_geomx_sh_vs_pt_padj < 0.05

# Sub-flag 7: GPNMB+ macrophage DE OR centrilobular/periportal epithelia DE.
# These three columns are pre-filtered to padj<0.05 upstream (any non-NA
# value qualifies), but we re-check defensively. Collapsed into one count
# so the total tops out at 7 (keeps the score interpretable).
gv_de_epimac <- (!is.na(atlas$signature_govaere2026_gpnmb_padj) &
                 atlas$signature_govaere2026_gpnmb_padj < 0.05) |
                (!is.na(atlas$signature_govaere2026_epithelia_mash_centrilobular_padj) &
                 atlas$signature_govaere2026_epithelia_mash_centrilobular_padj < 0.05) |
                (!is.na(atlas$signature_govaere2026_epithelia_mash_periportal_padj) &
                 atlas$signature_govaere2026_epithelia_mash_periportal_padj < 0.05)

# Aggregate count (0-7)
govaere2026_n_evidence <- as.integer(gv_metmac) + as.integer(gv_lam) +
                          as.integer(gv_gpnmb_mac) + as.integer(gv_il32) +
                          as.integer(gv_geomx_ls) + as.integer(gv_geomx_pt) +
                          as.integer(gv_de_epimac)

# Boolean any-evidence flag
govaere2026_evidence <- govaere2026_n_evidence > 0

cat(sprintf("  Govaere2026 evidence per sub-flag: metmac=%d lam=%d gpnmb_mac=%d il32=%d geomx_ls=%d geomx_pt=%d de_epi/mac=%d\n",
            sum(gv_metmac), sum(gv_lam), sum(gv_gpnmb_mac), sum(gv_il32),
            sum(gv_geomx_ls), sum(gv_geomx_pt), sum(gv_de_epimac)))
cat(sprintf("  Govaere2026 n_evidence distribution (non-excluded):\n"))
print(table(govaere2026_n_evidence[!excluded]))

# Soft boost to log_convergence_odds. log(3) per evidence unit roughly puts
# one Govaere sub-flag on par with one "active" modality (LOG_BF_ACTIVE=log(3)).
# Cap at 4 evidence units so genes with weak primary evidence can't be
# promoted into the top ranks on signature membership alone — the boost
# breaks ties among genes that already carry primary evidence.
GOVAERE_BOOST_PER_EVIDENCE <- log(3)     # 1.0986
GOVAERE_MAX_BOOST          <- 4 * GOVAERE_BOOST_PER_EVIDENCE   # 4.394
govaere_boost <- pmin(govaere2026_n_evidence * GOVAERE_BOOST_PER_EVIDENCE,
                      GOVAERE_MAX_BOOST)

log_convergence_odds_with_govaere2026 <- log_convergence_odds + govaere_boost
convergence_score_with_govaere2026 <- 1 / (1 + exp(-log_convergence_odds_with_govaere2026))

# Re-rank (NA for excluded)
convergence_rank_with_govaere2026 <- rep(NA_integer_, n)
ord2 <- order(-convergence_score_with_govaere2026[ok_idx])
convergence_rank_with_govaere2026[ok_idx[ord2]] <- seq_along(ord2)

cat(sprintf("  Boosted %d genes with n_evidence>=1, %d with n_evidence>=2, %d with n_evidence>=3\n",
            sum(govaere2026_n_evidence >= 1 & !excluded),
            sum(govaere2026_n_evidence >= 2 & !excluded),
            sum(govaere2026_n_evidence >= 3 & !excluded)))

# Sanity: rank changes for canonical Govaere panel test genes
gv_test_genes <- c("GPNMB","LPL","FABP5","IL32","AKR1B10","HLA-DRA")
cat("  Test gene rank changes (canonical -> with_govaere2026):\n")
for (sym in gv_test_genes) {
  idx <- which(atlas$human_symbol == sym)
  if (length(idx) == 0) {
    cat(sprintf("    %-10s : MISSING\n", sym))
    next
  }
  r0 <- convergence_rank[idx[1]]
  r1 <- convergence_rank_with_govaere2026[idx[1]]
  nev <- govaere2026_n_evidence[idx[1]]
  cat(sprintf("    %-10s : rank %s -> %s (n_evidence=%d, excluded=%s)\n",
              sym,
              ifelse(is.na(r0), "NA", as.character(r0)),
              ifelse(is.na(r1), "NA", as.character(r1)),
              nev, excluded[idx[1]]))
}

# ============================================================================
# 11. PRIMARY OUTPUT CSV
# ============================================================================
cat("\n--- 11. Write primary CSV ---\n")

# Druggability tier overlay (reporting-only)
drug_tier <- fifelse(!is.na(atlas$opentargets_drug) & atlas$opentargets_drug == TRUE,
                     "OpenTargets_druggable",
                     fifelse(!is.na(atlas$dgidb_druggable) & atlas$dgidb_druggable == TRUE,
                              "DGIdb_druggable", "None"))

# Zonation class
zonation_class_S5 <- if ("zonation_class" %in% names(atlas)) {
  atlas$zonation_class
} else {
  rep(NA_character_, n)
}

out <- data.table(
  human_symbol            = atlas$human_symbol,
  ensembl_id              = atlas$ensembl_id,
  excluded_from_ranking   = excluded,
  exclusion_reason        = reason,
  convergence_score          = convergence_score,
  log_convergence_odds           = log_convergence_odds,
  convergence_rank          = convergence_rank,
  evidence_total          = evidence_total,
  evidence_unsigned       = E_unsigned,
  evidence_signed_sum     = E_signed,
  evidence_signed_magnitude = abs(E_signed),
  concordance_state       = concordance_state,
  concordance_ratio       = concordance_ratio,
  tier                    = tier,
  n_modalities_active     = n_modalities_active,
  n_signed_concordant     = n_signed_concordant,
  log_BF_S1               = log_bf_S1,
  log_BF_S2_coloc        = log_bf_S2_coloc,
  log_BF_S2a_diag         = log_bf_S2a,   # diagnostic: legacy COLOC channel
  log_BF_S2b_diag         = log_bf_S2b,   # diagnostic: legacy TWAS channel
  coloc_genetic_pp4       = s2_coloc$score,
  coloc_channel_source           = s2_coloc$source,
  coloc_best_gwas     = s2_coloc$source_vec,  # per-gene multi vs ct_single (review A10#5)
  log_BF_S3               = log_bf_S3,
  log_BF_S4               = log_bf_S4,
  log_BF_S5               = log_bf_S5,
  log_BF_S6               = log_bf_S6,
  log_BF_S7               = log_bf_S7,
  log_BF_S8               = log_bf_S8,
  sign_S1                 = sign_S1,
  sign_S2_coloc          = sign_S2_coloc,
  sign_S4                 = sign_S4,
  sign_S6                 = sign_S6,
  sign_S7                 = sign_S7,
  sign_S8                 = sign_S8,
  best_contrast_S1        = best_contrast_S1,
  best_contrast_S4        = best_contrast_S4,
  best_contrast_S6        = best_contrast_S6,
  dominant_stage_S1       = dominant_stage_S1,
  cross_species_concordant_S8 = cross_species_concordant,
  zonation_class_S5       = zonation_class_S5,
  druggability_tier       = drug_tier,
  coloc_best_gwas_S2a     = s2a$best_gwas,
  coloc_best_pp4_S2a      = s2a$pp4,
  twas_z_S2b              = s2b$z,
  inhibitor_target_candidate = inhibitor_target_candidate,
  # --- 2026-05-21 Govaere2026 signature extension (additive only) ---
  signature_govaere2026_evidence      = govaere2026_evidence,
  signature_govaere2026_n_evidence    = govaere2026_n_evidence,
  convergence_score_with_govaere2026  = convergence_score_with_govaere2026,
  convergence_rank_with_govaere2026   = convergence_rank_with_govaere2026,
  log_convergence_odds_with_govaere2026 = log_convergence_odds_with_govaere2026
)
# Backup prior CSV before overwriting (CLAUDE.md atlas-rebuild pattern).
prior_csv <- file.path(ME, sprintf("%s.csv", OUT_PREFIX))
if (file.exists(prior_csv)) {
  archive_dir <- file.path(ME, "archive")
  dir.create(archive_dir, showWarnings = FALSE, recursive = TRUE)
  bk <- file.path(archive_dir,
                  sprintf("%s_backup_%s.csv", OUT_PREFIX,
                          format(Sys.time(), "%Y%m%d_%H%M%S")))
  file.copy(prior_csv, bk, overwrite = FALSE)
  cat(sprintf("  Backed up prior CSV to %s\n", bk))
}
fwrite(out, file.path(ME, sprintf("%s.csv", OUT_PREFIX)))
cat(sprintf("  Wrote %s.csv (%d rows, %d cols)\n",
            OUT_PREFIX, nrow(out), ncol(out)))

# Genetic+down-coherent (renamed from Protective-LOF) report — algorithmic
# state, NOT a biological claim that all entries are LOF-protective targets.
gdc_out <- out[concordance_state == "Genetic+down-coherent" & !excluded][order(-convergence_score)]
fwrite(gdc_out, file.path(ME, sprintf("%s_genetic_down_coherent.csv", OUT_PREFIX)))
cat(sprintf("  Wrote genetic_down_coherent.csv (%d genes)\n", nrow(gdc_out)))

# Inhibitor-target-candidate subset — the biologically defensible "inhibition
# would be protective" subset (TWAS direction implies expression-up = risk).
inh_out <- out[inhibitor_target_candidate == TRUE & !excluded][order(-convergence_score)]
fwrite(inh_out, file.path(ME, sprintf("%s_inhibitor_target_candidates.csv", OUT_PREFIX)))
cat(sprintf("  Wrote inhibitor_target_candidates.csv (%d genes)\n", nrow(inh_out)))

# ============================================================================
# 12. MODALITY CORRELATION DIAGNOSTIC
# ============================================================================
cat("\n--- 12. Modality correlation diagnostic ---\n")
mod_mat <- log_bf_all_mat[!excluded, ]
mod_cor <- suppressWarnings(cor(mod_mat, method = "spearman", use = "pairwise.complete.obs"))
mod_cor_dt <- as.data.table(as.data.frame(as.table(mod_cor)))
setnames(mod_cor_dt, c("Mod_A","Mod_B","Spearman_rho"))
fwrite(mod_cor_dt, file.path(ME, sprintf("%s_modality_correlations.csv", OUT_PREFIX)))
flagged <- mod_cor_dt[Mod_A != Mod_B & abs(Spearman_rho) > 0.5]
if (nrow(flagged) > 0) {
  cat(sprintf("  WARNING: %d modality pairs with |rho|>0.5:\n", nrow(flagged)/2))
  print(unique(flagged[abs(Spearman_rho) > 0.5, .(Mod_A, Mod_B, Spearman_rho)]))
}

# ============================================================================
# 13. PERMUTATION FDR (sign-shuffle null)
# ============================================================================
cat(sprintf("\n--- 13. Permutation FDR (%d shuffles) ---\n", N_PERM))

# Null: shuffle effect SIGNS within each signed modality independently; keep
# magnitudes (log_BF) fixed. This tests the null that the *pattern* of
# directional agreement is random.

permute_once <- function() {
  n_signed_mods <- 6
  sn_perm <- sn_signed
  for (j in seq_len(n_signed_mods)) {
    sn_perm[, j] <- sample(sn_signed[, j], length(sn_signed[, j]))
  }
  sn_act_p <- sn_perm * active_signed
  sn_act_p[is.na(sn_act_p)] <- 0
  E_s_p <- rowSums(sn_act_p * bf_active)
  # Reassess concordance_state under permutation, simplified:
  # just use the |E_signed|_p as the signed-evidence contribution.
  ev_p <- (E_unsigned + abs(E_s_p)) / brown_discount
  log_po_p <- log_prior_odds + ev_p
  1 / (1 + exp(-log_po_p))
}

perm_top <- matrix(0, nrow = N_PERM, ncol = 5)
colnames(perm_top) <- c("top100","top500","top1000","top5000","top10000")
n_nonexc <- sum(!excluded)

for (p in seq_len(N_PERM)) {
  pp_p <- permute_once()
  pp_p_ok <- pp_p[!excluded]
  cut100  <- sort(pp_p_ok, decreasing = TRUE)[min(100,  n_nonexc)]
  cut500  <- sort(pp_p_ok, decreasing = TRUE)[min(500,  n_nonexc)]
  cut1000 <- sort(pp_p_ok, decreasing = TRUE)[min(1000, n_nonexc)]
  cut5000 <- sort(pp_p_ok, decreasing = TRUE)[min(5000, n_nonexc)]
  cut10000<- sort(pp_p_ok, decreasing = TRUE)[min(10000,n_nonexc)]
  perm_top[p, ] <- c(cut100, cut500, cut1000, cut5000, cut10000)
  if (p %% max(1, N_PERM %/% 10) == 0) cat(sprintf("    perm %d/%d\n", p, N_PERM))
}

# Empirical FDR: at each rank threshold, expected # false positives under null /
# observed # discoveries
pp_obs <- convergence_score[!excluded]
fdr_tbl <- data.table(
  rank = c(100, 500, 1000, 5000, 10000),
  observed_cutoff = c(
    sort(pp_obs, decreasing = TRUE)[100],
    sort(pp_obs, decreasing = TRUE)[500],
    sort(pp_obs, decreasing = TRUE)[1000],
    sort(pp_obs, decreasing = TRUE)[5000],
    sort(pp_obs, decreasing = TRUE)[10000]
  ),
  null_mean_cutoff = colMeans(perm_top),
  null_95pct_cutoff = apply(perm_top, 2, quantile, 0.95)
)
fdr_tbl[, signif_vs_null := observed_cutoff > null_95pct_cutoff]
fwrite(fdr_tbl, file.path(ME, sprintf("%s_permutation_fdr.csv", OUT_PREFIX)))
print(fdr_tbl)

# ============================================================================
# 14. HELD-OUT PANEL BENCHMARKS
# ============================================================================
cat("\n--- 14. Held-out panel benchmarks ---\n")

load_panel <- function(path, name) {
  if (!file.exists(path)) {
    cat(sprintf("  SKIP %s: file missing (%s)\n", name, path))
    return(NULL)
  }
  # Panels may have # comment lines before the header; find the header row
  # explicitly so fread's comment-stripping doesn't misparse.
  raw <- readLines(path)
  hdr <- which(startsWith(raw, "gene_symbol"))[1]
  if (is.na(hdr)) {
    cat(sprintf("  SKIP %s: no gene_symbol header found\n", name))
    return(NULL)
  }
  dt <- fread(path, skip = hdr - 1, sep = "\t", header = TRUE)
  if (!"gene_symbol" %in% names(dt)) {
    cat(sprintf("  SKIP %s: no gene_symbol column\n", name))
    return(NULL)
  }
  genes <- unique(dt$gene_symbol)
  genes <- genes[!is.na(genes) & nchar(genes) > 0]
  cat(sprintf("  %s: %d unique genes\n", name, length(genes)))
  genes
}

# AUROC benchmark. IMPORTANT: require both score vectors to be non-NA before
# ranking — otherwise NA propagation from genes not in 46b's output artificially
# depresses 46b's AUROC (NAs rank last, inflating non-panel scores).
bench_panel <- function(panel_genes, panel_name, score_vec, sym_vec, excluded_vec) {
  idx_ok <- !excluded_vec & !is.na(sym_vec) & sym_vec != "" & !is.na(score_vec)
  scores <- score_vec[idx_ok]
  labels <- sym_vec[idx_ok] %in% panel_genes
  if (sum(labels) < 5) {
    cat(sprintf("    %s: only %d members in atlas; SKIP AUROC\n", panel_name, sum(labels)))
    return(NULL)
  }
  r <- rank(scores)
  nP <- sum(labels); nN <- sum(!labels)
  auroc <- (sum(r[labels]) - nP * (nP + 1) / 2) / (nP * nN)
  wx <- wilcox.test(scores[labels], scores[!labels], alternative = "greater")
  list(panel = panel_name, n_in_panel = length(panel_genes),
       n_in_atlas = sum(labels), n_scored = sum(idx_ok),
       auroc = auroc, p_value = wx$p.value,
       median_rank = median(rank(-scores)[labels]))
}

# Compare against 46b for reference
bp_46b <- file.path(ME, "bayesian_posterior.csv")
if (file.exists(bp_46b)) {
  b46b <- fread(bp_46b)
  score_46b <- rep(NA_real_, nrow(atlas))
  idx <- match(b46b$human_symbol, atlas$human_symbol)
  valid <- !is.na(idx)
  score_46b[idx[valid]] <- b46b$posterior_odds[valid]
} else {
  score_46b <- NULL
  cat("  46b output missing; skipping comparison\n")
}

panels <- list(
  Govaere     = load_panel(file.path(PANELS_DIR, "govaere_2020_panel.tsv"), "Govaere"),
  Feng        = load_panel(file.path(PANELS_DIR, "feng_2024_panel.tsv"),    "Feng"),
  NIDDK       = load_panel(file.path(PANELS_DIR, "niddk_pipeline_2024.tsv"),"NIDDK"),
  Liu_2024    = load_panel(file.path(PANELS_DIR, "liu_2024_crispr_hits.tsv"),"Liu_2024"),
  OpenTargets = load_panel(file.path(PANELS_DIR, "opentargets_masld_2025.tsv"),"OpenTargets")
)

# For fair 46d-vs-46b AUROC, restrict to the common non-excluded subset where
# both scores are defined. Otherwise 46b's ~19k missing genes would pollute the
# comparison by dominating rank statistics.
common_mask_46b <- if (!is.null(score_46b)) !is.na(score_46b) & !excluded else NULL

bench_rows <- list()
for (pname in names(panels)) {
  if (is.null(panels[[pname]])) next
  # 46d: evaluated on its full non-excluded set
  r46d <- bench_panel(panels[[pname]], pname, convergence_score,
                       atlas$human_symbol, excluded)
  # 46b: evaluated on common 46b ∩ 46d subset
  r46b <- NULL
  r46d_common <- NULL
  if (!is.null(score_46b) && !is.null(common_mask_46b)) {
    # Build an "excluded" mask that is TRUE wherever 46b is NA
    excl_46b <- excluded | is.na(score_46b)
    r46b <- bench_panel(panels[[pname]], pname, score_46b,
                         atlas$human_symbol, excl_46b)
    # Also evaluate 46d on the *same* subset for a fair delta
    r46d_common <- bench_panel(panels[[pname]], pname, convergence_score,
                                atlas$human_symbol, excl_46b)
  }
  row <- data.table(
    panel = pname,
    n_in_panel_atlas = r46d$n_in_atlas,
    n_scored_46d     = r46d$n_scored,
    n_scored_common  = if (!is.null(r46b)) r46b$n_scored else NA_integer_,
    auroc_46d        = r46d$auroc,
    auroc_46d_common = if (!is.null(r46d_common)) r46d_common$auroc else NA_real_,
    auroc_46b        = if (!is.null(r46b)) r46b$auroc else NA_real_,
    delta_fair       = if (!is.null(r46b) && !is.null(r46d_common))
                         r46d_common$auroc - r46b$auroc else NA_real_,
    p_46d            = r46d$p_value,
    p_46b            = if (!is.null(r46b)) r46b$p_value else NA_real_,
    median_rank_46d  = r46d$median_rank,
    median_rank_46b  = if (!is.null(r46b)) r46b$median_rank else NA_real_
  )
  bench_rows[[pname]] <- row
  cat(sprintf("  %s: 46d_AUROC=%.3f  46b_AUROC=%s  Δ(fair)=%s  n_atlas=%d  n_common=%s\n",
              pname, r46d$auroc,
              ifelse(!is.null(r46b), sprintf("%.3f", r46b$auroc), "NA"),
              ifelse(!is.null(r46d_common),
                     sprintf("%+.3f", r46d_common$auroc - r46b$auroc), "NA"),
              r46d$n_in_atlas,
              ifelse(!is.null(r46b), as.character(r46b$n_scored), "NA")))
}
bench_dt <- rbindlist(bench_rows, fill = TRUE)
if (nrow(bench_dt) > 0) {
  fwrite(bench_dt, file.path(ME, sprintf("%s_benchmark.csv", OUT_PREFIX)))
} else {
  cat("  No panels available; skipping benchmark CSV\n")
}

# Conserved enrichment
if ("is_conserved" %in% names(atlas)) {
  cc <- atlas$is_conserved %in% TRUE
  top500 <- head(which(!excluded)[order(-convergence_score[!excluded])], 500)
  ov <- sum(cc[top500])
  ftest <- fisher.test(table(
    in_top500 = seq_along(cc) %in% top500,
    cc = cc
  ))
  cat(sprintf("  Conserved top-500 Fisher: OR=%.2f p=%.2e (overlap=%d/%d)\n",
              ftest$estimate, ftest$p.value, ov, sum(cc)))
}

# 46b-vs-46d divergence
if (!is.null(score_46b)) {
  cmp <- data.table(
    human_symbol = atlas$human_symbol,
    excluded = excluded,
    rank_46d = convergence_rank,
    rank_46b = rank(-score_46b, ties.method = "first"),
    convergence_score_46d = convergence_score,
    posterior_odds_46b = score_46b,
    tier_46d = tier,
    concordance_state_46d = concordance_state
  )
  cmp[, rank_diff := rank_46b - rank_46d]  # positive = 46d ranks higher
  fwrite(cmp[order(-rank_diff)][!is.na(rank_46d) & !is.na(rank_46b)],
         file.path(ME, sprintf("%s_vs_46b.csv", OUT_PREFIX)))
  top_46d_rev <- cmp[!excluded & rank_46b > 500 & rank_46d < 100]
  cat(sprintf("  46d-reveals-novel (46b rank>500 AND 46d rank<100): %d genes\n",
              nrow(top_46d_rev)))
}

# ============================================================================
# 15. SENSITIVITY SWEEPS
# ============================================================================
cat("\n--- 15. Sensitivity sweeps ---\n")

compute_posterior_with <- function(pi_x, keff_s1_x, keff_s4_x, keff_s6_x,
                                   skip_mod = NULL) {
  bf_S1 <- agg_modality(s1_mat, keff_s1_x)
  bf_S4 <- agg_modality(s4_mat, keff_s4_x)
  bf_S6 <- agg_modality(s6_mat, keff_s6_x)
  mods <- list(S1 = bf_S1, S2_coloc = log_bf_S2_coloc,
               S3 = log_bf_S3, S4 = bf_S4, S5 = log_bf_S5,
               S6 = bf_S6, S7 = log_bf_S7, S8 = log_bf_S8)
  if (!is.null(skip_mod)) mods[[skip_mod]] <- rep(0, n)

  # Mirror main score: E_unsigned = S5 (S2_coloc is signed, joins below)
  E_unsigned_x <- pmax(mods$S5, 0)

  # Rebuild signed bf matrix (6 channels: S1, S2_coloc, S4, S6, S7, S8)
  bf_sgn <- cbind(mods$S1, mods$S2_coloc, mods$S4, mods$S6, mods$S7, mods$S8)
  act <- bf_sgn > LOG_BF_ACTIVE
  sn_act <- sn_signed * act
  bf_act <- bf_sgn * act
  bf_act[is.na(bf_act)] <- 0; sn_act[is.na(sn_act)] <- 0
  E_s <- rowSums(sn_act * bf_act)
  M_s <- rowSums(bf_act)

  # Simplified: use |E_s| (no state classification in sensitivity)
  ev_x <- (E_unsigned_x + abs(E_s)) / brown_discount   # apply Brown's correction
  log_po_x <- log(pi_x / (1 - pi_x)) + ev_x
  1 / (1 + exp(-log_po_x))
}

# Baseline top-200 (for stability comparisons)
top200_base <- which(!excluded)[order(-convergence_score[!excluded])[1:200]]

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
stability_row <- function(label, top200_x) {
  top200_x_idx <- which(!excluded)[order(-top200_x[!excluded])[1:200]]
  data.table(
    setting = label,
    jaccard_top200 = jaccard(top200_base, top200_x_idx),
    spearman_top200 = cor(convergence_score[top200_base], top200_x[top200_base],
                          method = "spearman", use = "complete.obs")
  )
}

sens_rows <- list()
# π sweep
for (pi_x in c(0.01, 0.10, 0.20)) {
  pp_x <- compute_posterior_with(pi_x, keff_s1$K_eff, keff_s4$K_eff, keff_s6$K_eff)
  sens_rows[[sprintf("pi=%.2f", pi_x)]] <- stability_row(sprintf("pi=%.2f", pi_x), pp_x)
}
# K_eff sweep (±50%)
for (mult in c(0.5, 1.5)) {
  pp_x <- compute_posterior_with(pi_hat,
    keff_s1$K_eff * mult, keff_s4$K_eff * mult, keff_s6$K_eff * mult)
  sens_rows[[sprintf("keff_x%.1f", mult)]] <- stability_row(sprintf("keff_x%.1f", mult), pp_x)
}
# LOMO
for (mod in c("S1","S2_coloc","S3","S4","S5","S6","S7","S8")) {
  pp_x <- compute_posterior_with(pi_hat, keff_s1$K_eff, keff_s4$K_eff, keff_s6$K_eff,
                                  skip_mod = mod)
  sens_rows[[sprintf("LOMO_%s", mod)]] <- stability_row(sprintf("LOMO_%s", mod), pp_x)
}
sens_dt <- rbindlist(sens_rows)
fwrite(sens_dt, file.path(ME, sprintf("%s_sensitivity.csv", OUT_PREFIX)))
print(sens_dt)

# ============================================================================
# 16. PER-GENE EVIDENCE CARDS (top 200)
# ============================================================================
cat("\n--- 16. Evidence cards ---\n")

write_card <- function(row_idx) {
  sym <- atlas$human_symbol[row_idx]
  if (is.na(sym) || sym == "") return(invisible(NULL))
  path <- file.path(CARDS_DIR, sprintf("%s.md", sym))
  # Card uses with-S6 matrix so the diagnostic S6 row is still rendered, but
  # n_modalities_active above is the 7-channel canonical count.
  bf_i <- log_bf_all_mat_with_s6[row_idx, ]
  lines <- c(
    sprintf("# %s (rank %d)", sym, convergence_rank[row_idx]),
    "",
    sprintf("- **Ensembl:** %s", atlas$ensembl_id[row_idx]),
    sprintf("- **Tier:** %s", tier[row_idx]),
    sprintf("- **Concordance state:** %s", concordance_state[row_idx]),
    sprintf("- **Convergence score:** %.4f", convergence_score[row_idx]),
    sprintf("- **log convergence odds:** %.3f", log_convergence_odds[row_idx]),
    sprintf("- **Evidence total:** %.2f  (unsigned: %.2f, signed mag: %.2f)",
            evidence_total[row_idx], E_unsigned[row_idx], abs(E_signed[row_idx])),
    sprintf("- **Modalities active:** %d/7", n_modalities_active[row_idx]),
    sprintf("- **Concordance ratio:** %s",
            ifelse(is.na(concordance_ratio[row_idx]), "NA",
                   sprintf("%.3f", concordance_ratio[row_idx]))),
    sprintf("- **Dominant stage (S1):** %s", dominant_stage_S1[row_idx]),
    sprintf("- **Cross-species concordant:** %s", cross_species_concordant[row_idx]),
    sprintf("- **Druggability:** %s", drug_tier[row_idx]),
    "",
    "## Per-modality log-BF",
    "",
    sprintf("| Modality | log-BF | Sign | Active? | Best sub-contrast |"),
    "|---|---|---|---|---|",
    sprintf("| S1 bulk      | %.2f | %s | %s | %s |",
            bf_i["S1"], sign_S1[row_idx], bf_i["S1"]>LOG_BF_ACTIVE, best_contrast_S1[row_idx]),
    sprintf("| S2 genetic (COLOC) | %.2f | %s | %s | COLOC PP.H4=%.3f, TWAS z=%.2f (dir only), GWAS=%s |",
            bf_i["S2_coloc"], sign_S2_coloc[row_idx],
            bf_i["S2_coloc"]>LOG_BF_ACTIVE,
            ifelse(is.na(s2_coloc$score[row_idx]),0,s2_coloc$score[row_idx]),
            ifelse(is.na(s2b$z[row_idx]),0,s2b$z[row_idx]),
            ifelse(is.na(s2a$best_gwas[row_idx]),"NA",s2a$best_gwas[row_idx])),
    sprintf("| S3 essential | %.2f | — | %s | Chronos=%.3f |",
            bf_i["S3"], bf_i["S3"]>LOG_BF_ACTIVE,
            ifelse(is.na(atlas$essentiality_chronos[row_idx]),0,
                    atlas$essentiality_chronos[row_idx])),
    sprintf("| S4 epigen.   | %.2f | %s | %s | %s |",
            bf_i["S4"], sign_S4[row_idx], bf_i["S4"]>LOG_BF_ACTIVE, best_contrast_S4[row_idx]),
    sprintf("| S5 spatial   | %.2f | — | %s | SVG=%s Moran_I=%.3f |",
            bf_i["S5"], bf_i["S5"]>LOG_BF_ACTIVE,
            atlas$spatial_is_svg[row_idx],
            ifelse(is.na(atlas$spatial_morans_i[row_idx]),0,
                    atlas$spatial_morans_i[row_idx])),
    sprintf("| S6 sc (diag) | %.2f | %s | %s | %s |",
            bf_i["S6"], sign_S6[row_idx], bf_i["S6"]>LOG_BF_ACTIVE, best_contrast_S6[row_idx]),
    sprintf("| S7 protein   | %.2f | %s | %s | logFC=%.3f |",
            bf_i["S7"], sign_S7[row_idx], bf_i["S7"]>LOG_BF_ACTIVE,
            ifelse(is.na(atlas$best_protein_logFC[row_idx]),0,
                    atlas$best_protein_logFC[row_idx])),
    sprintf("| S8 mouse     | %.2f | %s | %s | logFC=%.3f ndiets=%d |",
            bf_i["S8"], sign_S8[row_idx], bf_i["S8"]>LOG_BF_ACTIVE,
            ifelse(is.na(atlas$mouse_meta_logFC[row_idx]),0,
                    atlas$mouse_meta_logFC[row_idx]),
            ifelse(is.na(atlas$n_diets_sig[row_idx]),0,
                    as.integer(atlas$n_diets_sig[row_idx])))
  )
  writeLines(lines, path)
}

## Emit cards for ALL tier1+tier2+tier3 genes (the "convergent" canonical set).
## After regeneration, sweep stale cards not in this set (E1 cleanup, 2026-05-23).
emit_tiers <- c("1_Genetic_validated", "2_Convergent", "3_Suggestive")
emit_idx <- which(tier %in% emit_tiers & !excluded &
                  !is.na(atlas$human_symbol) & atlas$human_symbol != "")
for (i in emit_idx) write_card(i)
cat(sprintf("  Wrote %d evidence cards to %s (tiers: %s)\n",
            length(emit_idx), CARDS_DIR, paste(emit_tiers, collapse = ", ")))

## Sweep: delete any leftover cards whose symbol is NOT in the freshly-emitted set.
emit_syms <- atlas$human_symbol[emit_idx]
existing_files <- list.files(CARDS_DIR, pattern = "\\.md$", full.names = TRUE)
existing_syms <- sub("\\.md$", "", basename(existing_files))
stale_files <- existing_files[!existing_syms %in% emit_syms]
n_stale <- length(stale_files)
if (n_stale > 0) {
  invisible(file.remove(stale_files))
  cat(sprintf("  Removed %d stale evidence card(s) (no longer tier1+2+3)\n", n_stale))
}

# ============================================================================
# 17. STDOUT SUMMARY
# ============================================================================
cat("\n=== SUMMARY ===\n")
cat(sprintf("Total genes: %d (excluded %d, ranked %d)\n",
            nrow(atlas), sum(excluded), sum(!excluded)))
cat(sprintf("Tier distribution: %s\n",
            paste(names(table(tier)), table(tier), sep="=", collapse=", ")))
cat(sprintf("Genetic+down-coherent genes: %d\n", sum(concordance_state == "Genetic+down-coherent")))
cat(sprintf("Inhibitor-target candidates (Genetic+down-coherent ∩ TWAS-up risk): %d\n",
            sum(inhibitor_target_candidate, na.rm = TRUE)))
cat("\nTop 20 by posterior rank:\n")
top20 <- out[excluded_from_ranking == FALSE][order(-convergence_score)][1:20,
  .(rank=convergence_rank, human_symbol, tier, concordance_state,
    post=round(convergence_score,4), n_mod=n_modalities_active,
    stage=dominant_stage_S1)]
print(top20)
cat("\nGenetic+down-coherent top 10:\n")
gdc_top <- out[concordance_state == "Genetic+down-coherent" & excluded_from_ranking == FALSE][
  order(-convergence_score)][1:min(10, .N),
  .(rank=convergence_rank, human_symbol, post=round(convergence_score,4),
    n_mod=n_modalities_active, log_BF_S2_coloc=round(log_BF_S2_coloc,2),
    intact=round(coloc_genetic_pp4,3),
    coloc_pp4=round(coloc_best_pp4_S2a,3),
    inh_cand=inhibitor_target_candidate)]
print(gdc_top)
cat("\nInhibitor-target-candidate top 10:\n")
inh_top <- out[inhibitor_target_candidate == TRUE & excluded_from_ranking == FALSE][
  order(-convergence_score)][1:min(10, .N),
  .(rank=convergence_rank, human_symbol, post=round(convergence_score,4),
    n_mod=n_modalities_active, twas_z=round(twas_z_S2b,2),
    coloc_pp4=round(coloc_best_pp4_S2a,3))]
print(inh_top)

cat("\nDone.\n")
